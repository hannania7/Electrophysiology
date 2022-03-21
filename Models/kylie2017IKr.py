import os, sys
import time, glob

import random
import math
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

# Parameter range setting
parameter_ranges = []
dataset_dir = 'Kylie'
if 'Kylie' in dataset_dir:    
    if 'rmax600' in dataset_dir:
        # Kylie: rmax = 600 
        parameter_ranges = []
        parameter_ranges.append( [100, 500000] )
        parameter_ranges.append( [0.0001, 598] )
        parameter_ranges.append( [0.0001, 260] )
        parameter_ranges.append( [0.0001, 598] )
        parameter_ranges.append( [0.0001, 130] )
        parameter_ranges.append( [0.0001, 598] )
        parameter_ranges.append( [0.0001, 260] )
        parameter_ranges.append( [0.0001, 598] )
        parameter_ranges.append( [0.0001, 130] )
        print("Kylie-rmax600 dataset has been selected.")
    else :
        # Kylie
        parameter_ranges.append([100, 500000])
        parameter_ranges.append( [0.0001, 1000000])
        parameter_ranges.append( [0.0001, 384])
        parameter_ranges.append( [0.0001, 1000000] )
        parameter_ranges.append( [0.0001, 192] )
        parameter_ranges.append( [0.0001, 1000000] )
        parameter_ranges.append( [0.0001, 384] )
        parameter_ranges.append( [0.0001, 1000000] )
        parameter_ranges.append( [0.0001, 192] )
        print("Kylie dataset has been selected.")

elif 'RealRange' in dataset_dir:
        parameter_ranges.append([3134, 500000])                 # g
        parameter_ranges.append( [0.0001, 2.6152843264828003])  # p1
        parameter_ranges.append( [43.33271226094526, 259])      # p2
        parameter_ranges.append( [0.001, 0.5] )                 # p3
        parameter_ranges.append( [15, 75] )                     # p4
        parameter_ranges.append( [0.8, 410] )                   # p5
        parameter_ranges.append( [0.0001, 138.] )               # p6
        parameter_ranges.append( [1.0, 59] )                    # p7
        parameter_ranges.append( [1.6, 90] )                    # p8
        print("RealRange dataset has been selected.")

parameter_ranges = np.array(parameter_ranges)
print(parameter_ranges.shape)

class Kylie2017IKr():
    def __init__(self, protocol):
        self.protocol = protocol
                
        self.open0 = 0
        self.active0 = 1
        R = 8.314472 # [J/mol/K]
        T = 310     # [K]  # 36-37oC (BT)
        F = 9.64853415e4 #[C/mol]
        self.RTF = R * T / F        
        self.Ki = 110  # [mM]
        #Ki = 125 [mM]  # for iPSC solution
        self.Ko = 4    # [mM]
        #Ko = 3.75 [mM]  # for iPSC solution
        self.EK = self.RTF * log(self.Ko / self.Ki)  # in [V]
        
        self.g = 0.1524 * 1e3 # [pA/V]
        self.p1 = 2.26e-4 * 1e3 # [1/s]
        self.p2 = 0.0699 * 1e3  # [1/V]
        self.p3 = 3.45e-5 * 1e3 # [1/s]
        self.p4 = 0.05462 * 1e3 # [1/V]
        self.p5 = 0.0873 * 1e3  # [1/s]
        self.p6 = 8.91e-3 * 1e3 # [1/V]
        self.p7 = 5.15e-3 * 1e3 # [1/s]
        self.p8 = 0.03158 * 1e3 # [1/V]

        self.y0 = [self.open0, self.active0]
        self.params = [self.p1, self.p2, self.p3, self.p4, self.p5, self.p6, self.p7, self.p8]

    def set_params(self, g, p1, p2, p3, p4, p5, p6, p7, p8):
        self.g = g # [pA/V]
        self.p1 = p1 # [1/s]
        self.p2 = p2  # [1/V]
        self.p3 = p3 # [1/s]
        self.p4 = p4 # [1/V]
        self.p5 = p5  # [1/s]
        self.p6 = p6 # [1/V]
        self.p7 = p7 # [1/s]
        self.p8 = p8 # [1/V]
        self.params = [g, p1, p2, p3, p4, p5, p6, p7, p8]

    def set_result(self, times, y, log=None):
        self.times =  times
        self.V = np.array(self.protocol.get_voltage_clamp_protocol(times))
        self.open = y[0]  
        self.active = y[1]                  
        self.IKr = self.g * self.open * self.active * (self.V - self.EK)
    
    def differential_eq(self, t, y):    
        a, r = y
        V = self.protocol.get_voltage_at_time(t)
        k1 = self.p1*exp(self.p2*V)
        k2 = self.p3*exp(-self.p4*V)
        k3 = self.p5*exp(self.p6*V)
        k4 = self.p7*exp(-self.p8*V)
        tau_a = 1.0/(k1+k2)
        tau_r = 1.0/(k3+k4)
        a_inf = k1/(k1+k2)
        r_inf = k4/(k3+k4) 
        da = (a_inf-a)/tau_a
        dr = (r_inf-r)/tau_r          
        return [da, dr]

    def response_diff_eq(self, t, y):
        
        if isinstance(self.protocol, protocol_lib.PacingProtocol)  :                      
            face = self.protocol.pacing(t)
            self.stimulus.cal_stimulation(face) # Stimulus    
                            
        return self.differential_eq(t, y)


    def diff_eq_solve_ivp(self, t, y):
        return self.response_diff_eq(t, y)
        
    def diff_eq_odeint(self, y, t, *p):
        return self.response_diff_eq(t, y)
        
        
        
    def simulate_odeint(self, t, g,p1,p2,p3,p4,p5,p6,p7,p8):       
    
        def myode(ar, t):           
            a, r = ar        
            V = self.protocol.voltage_at_time(t)    
            k1 = p1*exp(p2*V)
            k2 = p3*exp(-p4*V)
            k3 = p5*exp(p6*V)
            k4 = p7*exp(-p8*V)
            tau_a = 1/(k1+k2)
            tau_r = 1/(k3+k4)
            a_inf = k1/(k1+k2)
            r_inf = k4/(k3+k4)    
            dot_a = (a_inf-a)/tau_a
            dot_r = (r_inf-r)/tau_r
            return [dot_a, dot_r]

        ar = odeint(myode, [self.open0, self.active0], t)
        a = ar[:, 0]
        r = ar[:, 1]            
        V = self.voltage(t)        
        IKr = g * a * r * (V - self.EK)
        return IKr
    
    def curve_fitting(self, times, data, p0, bounds=None, method=None):
        fit_p, pcov = curve_fit(self.simulate_odeint, times, data, p0=p0, bounds=bounds ,method=method)
        return fit_p
        


def main():
    start_time = time.time()

    print("--- %s seconds ---"%(time.time()-start_time))


if __name__ == '__main__':
    main()