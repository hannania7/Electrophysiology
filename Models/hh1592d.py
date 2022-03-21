import os, sys
import time, glob

import random
import numpy as np
import matplotlib.pyplot as plt

import multiprocessing
from functools import partial 
from tqdm import tqdm

# import pickle
# import bisect

class Membrane():
    def __init__(self):
        
        self.C = 1        # [uF/cm^2] : The membrane capacity per unit area (See table 3 in HH1952d)        
    
    def dot_V(self, Iions):      # INa  
        """        
        """                                
        return -(1/self.C) * Iions
    
    
class Potassium():
    def __init__(self):
                
        self.Ek = -87 # [mV]
        self.GK_max = 36 # [mS/cm^2]

    def dot_n(self, n, V):        
        """
        
        """        
        a1 = 0.01 * (-V - 65.0) / (np.exp((-V - 65.0) / 10.0) - 1)
        b1 = 0.125 * np.exp((-V - 75.0) / 80.0)        
        return a1*(1-n) - b1*n
    
    def I(self, n, V):        
        """        
        """                      
        return self.GK_max * n**4 * (V - self.Ek)  # [uA/cm^2] Current carried by potassium ions
    
    
class Sodium():
    def __init__(self):
        
        self.ENa = 40.0  # [mV]
        self.GNa_max = 120.0  # [mS/cm^2]

    def dot_m(self, m, V):        
        """
        
        """        
        a2 = 0.1 * (-V - 50.0) / (np.exp((-V - 50.0) / 10.0) - 1)
        b2 = 4.0 * np.exp((-V - 75.0) / 18.0)     
        return a2 * (1.0 - m) - b2 * m     
    
    def dot_h(self, h, V):        
        """
        
        """        
        a3 = 0.07 * np.exp((-V - 75.0) / 20.0)
        b3 = 1.0 / (np.exp((-V - 45.0) / 10.0) + 1.0)        
        return a3 * (1.0 - h) - b3 * h      
    
    def I(self, m, h, V):      # INa  
        """        
        """                      
        return self.GNa_max * m**3 * h * (V - self.ENa)  # desc: Current carried by Sodium ions in [uA/cm^2] 

    
class Leak():
    def __init__(self):
        
        self.Eleak = -64.387
        self.GLeak_max = 0.3 # [mS/cm^2]
    
    def I(self, V):      # INa  
        """        
        """                    
        return self.GLeak_max * (V - self.Eleak)
    
    
class HH1592d():
    """
    Hodgkin, Huxley, 1952d, Journal of Physiology
    """
    def __init__(self, protocol):

        self.name = "HH1592d"
        
        self.protocol = protocol
                
        self.membrane_V0  = -75.0
        self.potassium_n0 = 0.317
        self.sodium_m0    = 0.05
        self.sodium_h0    = 0.595
        
        self.Membrane = Membrane()
        self.K = Potassium()
        self.Na = Sodium()
        self.Leak = Leak()
        
        self.y0 = [self.membrane_V0, self.potassium_n0, self.sodium_m0, self.sodium_h0]
        self.params = []

    def set_result(self, t, y, log=[]):
        
        self.times =  t
        self.V = y[0]      
        self.n = y[1]
        self.m = y[2]
        self.h = y[3]
#         if len(log)>0:
#             if 'IK' in log:
#                 self.IK = self.K.I(self.n, self.V)
#             if 'INa' in log:
#                 self.INa = self.Na.I(self.m, self.h, self.V)
#             if 'ILeak' in log:
#                 self.ILeak = self.Leak.I( self.V )        
              
    def differential_eq(self, t, y):    
        V, n, m, h = y
       
        # Potassium current
        dn = self.K.dot_n(n, V)
        IK = self.K.I(n,V)        

        # Sodium current   
        INa = self.Na.I(m, h, V)
        dm = self.Na.dot_m(m, V)        
        dh = self.Na.dot_h(h, V)

        # Leak current        
        ILeak = self.Leak.I(V)

        # Stimulus
        Vhold = -60  # [mV] : A temporary holding potential
        A = 100 # An amplification factor for the holding current        
        i_stim = (V - Vhold) * A * self.protocol.pacing(t)

        # Membrane        
        dV =self.Membrane.dot_V(INa + IK + ILeak + i_stim)
                    
        return [dV, dn, dm, dh]

    def response_diff_eq(self, t, y):
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