import os, sys
import time, glob
import random
from math import log, sqrt, floor, exp
import numpy as np
import matplotlib.pyplot as plt

import multiprocessing
from functools import partial 
from tqdm import tqdm
# import pickle
# import bisect

sys.path.append('../Protocols')
sys.path.append('../Lib')
import protocol_lib
import mod_trace
        
        
class Membrane():
    def __init__(self):
        self.C = 1  # [uF/cm^2] : The membrane capacitance
        self.V = -84.622        
    
    def dot_V(self, Iions):
        '''
        in [mV]
        Membrane potential
        '''
        return -(1.0/self.C) * Iions
    
class Stimulus():
    def __init__(self):
        self.amplitude = 25 # [uA/cm^2]
        self.I = 0
        
    def cal_stimulation(self, pace):
        self.I = self.amplitude * pace     
        return self.I
                
class INa():
    def __init__(self):
        self.m = 0.01
        self.h = 0.99
        self.j = 0.98
        
        self.gNaBar = 4   # [mS/cm^2]
        self.gNaC = 0.003 # [mS/cm^2]
        self.ENa = 50     # [mV]
                
    def dot_m(self, m, V):
        '''
        The activation parameter       
        '''
        alpha = (V + 47) / (1 - exp(-0.1 * (V + 47)))
        beta  = 40 * exp(-0.056 * (V + 72))        
        return alpha * (1 - m) - beta * m  # 
    
    def dot_h(self, h, V):
        '''
        An inactivation parameter   
        '''
        alpha = 0.126 * exp(-0.25 * (V + 77))
        beta  = 1.7 / (1 + exp(-0.082 * (V + 22.5)))
        return alpha * (1 - h) - beta * h  # 
    
    def dot_j(self, j, V):
        '''
        An inactivation parameter
        '''
        alpha = 0.055 * exp(-0.25 * (V + 78)) / (1 + exp(-0.2 * (V + 78)))
        beta  = 0.3 / (1 + exp(-0.1 * (V + 32)))
        return alpha * (1 - j) - beta * j  # 
            
    def I(self, m, h, j, V):        
        """
        in [uA/cm^2]
        The excitatory inward sodium current
        """                
        return (self.gNaBar * m**3 * h * j + self.gNaC) * (V - self.ENa)
        
        
class Isi():
    def __init__(self):
        self.d = 0.003
        self.f = 0.99
        self.Cai = 2e-7
        
        self.gsBar = 0.09    

    def dot_d(self, d, V):
        alpha = 0.095 * exp(-0.01 * (V + -5)) / (exp(-0.072 * (V + -5)) + 1)
        beta  = 0.07 * exp(-0.017 * (V + 44)) / (exp(0.05 * (V + 44)) + 1)
        return alpha * (1 - d) - beta * d
        
    def dot_f(self, f, V):
        alpha = 0.012 * exp(-0.008 * (V + 28)) / (exp(0.15 * (V + 28)) + 1)
        beta  = 0.0065 * exp(-0.02 * (V + 30)) / (exp(-0.2 * (V + 30)) + 1)
        return alpha * (1 - f) - beta * f
        
    def dot_Cai(self, Isi, Cai):
        '''
        desc: The intracellular Calcium concentration
        in [mol/L]
        '''
        return -1e-7 * Isi + 0.07 * (1e-7 - Cai)
   
    def I(self, d, f, Cai, V):        
        """
        in [uA/cm^2]
        The slow inward current, primarily carried by calcium ions. Called either
        "iCa" or "is" in the paper.
        """        
        Es = - 82.3 - 13.0287 * log(Cai)        # in [mV]
        return self.gsBar * d * f * (V - Es)         # in [uA/cm^2]
    

class IK1():
    def __init__(self):
        self.xx = 0
    
    def I(self, V):        
        """
        in [uA/cm^2]
        A time-independent outward potassium current exhibiting inward-going rectification
        """
        return 0.35 * (4 * (exp(0.04 * (V + 85)) - 1) / (exp(0.08 * (V + 53)) + exp(0.04 * (V + 53))) + 0.2 * (V + 23) / (1 - exp(-0.04 * (V + 23)))
    )
        
        
class Ix1():
    def __init__(self):
        self.x1 = 0.0004
        
    def dot_x1(self, x1, V):
        alpha = 0.0005 * exp(0.083 * (V + 50)) / (exp(0.057 * (V + 50)) + 1)
        beta  = 0.0013 * exp(-0.06 * (V + 20)) / (exp(-0.04 * (V + 333)) + 1)
        return alpha * (1 - x1) - beta * x1
    
    def I(self, x1, V):        
        """
        in [uA/cm^2]
        A voltage- and time-dependent outward current, primarily carried by potassium ions 
        """        
        return x1 * 0.8 * (exp(0.04 * (V + 77)) - 1) / exp(0.04 * (V + 35))
    
        
        
class BR1977():
    """    
    Beeler and Reuter 1977
    """
    def __init__(self, protocol=None):

        self.name = "BR1977"
        
        self.current_response_info = mod_trace.CurrentResponseInfo(protocol)
        
        self.protocol = protocol
        
        self.membrane = Membrane()
        self.stimulus = Stimulus()        
        self.ina = INa()       
        self.isi = Isi()       
        self.ik1 = IK1()       
        self.ix1 = Ix1()               

        self.y0 = [self.membrane.V, self.ina.m, self.ina.h, self.ina.j, self.isi.d, self.isi.f, self.isi.Cai, self.ix1.x1]
        self.params = []

    def set_result(self, t, y, log=None):
        self.times =  t
        self.V = y[0]    
        self.m = y[1]
        self.h = y[2]
        self.j = y[3]
        self.d = y[4]
        self.f = y[5]
        self.Cai = y[6]
        self.x1 = y[7]               
        
    def differential_eq(self, t, y):    
        V, m, h, j, d, f, Cai, x1 = y
        
        # INa        
        INa = self.ina.I(m, h, j, V)
        dot_m = self.ina.dot_m(m, V)
        dot_h = self.ina.dot_h(h, V)
        dot_j = self.ina.dot_j(j, V)        
                
        # Isi        
        Isi = self.isi.I(d, f, Cai, V)
        dot_d = self.isi.dot_d(d, V)       
        dot_f = self.isi.dot_f(f, V)       
        dot_Cai = self.isi.dot_Cai(Isi, Cai)       
            
        # IK1
        IK1 = self.ik1.I(V)        
    
        # Ix1           
        Ix1 = self.ix1.I(x1, V)
        dot_x1 = self.ix1.dot_x1(x1, V)
                   
        # Membrane potential        
        dot_V = self.membrane.dot_V(IK1 + Ix1 + INa + Isi - self.stimulus.I)
        
        if self.current_response_info:
            current_timestep = [
                mod_trace.Current(name='I_Na', value=INa),
                mod_trace.Current(name='I_si', value=Isi),
                mod_trace.Current(name='I_K1', value=IK1),
                mod_trace.Current(name='I_x1', value=Ix1),               
            ]
            self.current_response_info.currents.append(current_timestep)
            
        return [dot_V, dot_m, dot_h, dot_j, dot_d, dot_f, dot_Cai, dot_x1]
    
    
    def response_diff_eq(self, t, y):
        
        if isinstance(self.protocol, protocol_lib.PacingProtocol)  :                      
            face = self.protocol.pacing(t)
            self.stimulus.cal_stimulation(face) # Stimulus    
        else:                         
            y[0] = self.protocol.get_voltage_at_time(t)
                    
        return self.differential_eq(t, y)


    def diff_eq_solve_ivp(self, t, y):
        return self.response_diff_eq(t, y)
        
    def diff_eq_odeint(self, y, t, *p):
        return self.response_diff_eq(t, y)


def main():
    start_time = time.time()

    print("--- %s seconds ---"%(time.time()-start_time))


if __name__ == '__main__':
    main()