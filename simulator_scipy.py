import os, sys
import copy
import time, glob

import random
from numbalsoda import lsoda_sig, lsoda
from numba import njit, cfunc
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import ode, solve_ivp, odeint
from scipy.optimize import curve_fit, least_squares

sys.path.append('./Protocols')
sys.path.append('./Lib')
import protocol_lib
import mod_trace

class Simulator:
    """
    
    """
    def __init__(self, model):
        '''
        '''        
        # Get model
        self.model = model
                        
        self.name = None        
        self.bcl = 1000
        self.vhold = 0  # -80e-3      -88.0145 
                        
#         self.times = np.linspace(0, self.bcl, 1000)  # np.arange(self._bcl)        
        
     
    # def set_times(self, times):
    #     self.times = times
    #     print("Times has been set.")

    def pre_simulate(self, pre_step=5000, protocol=None, 
                           method='BDF', max_step=np.inf, atol=1E-6, rtol=1E-3, t_eval=None):
        '''
        ''' 
        protocol_temp = copy.copy(self.model.protocol)
        self.model.protocol = protocol           
        if protocol == 'constant':            
            self.model.protocol = protocol_lib.VoltageClampProtocol( [protocol_lib.VoltageClampStep(voltage=self.model.y0[0], duration=pre_step)] )
        elif protocol =='pacing':
            self.model.protocol = protocol_lib.PacingProtocol(level=-1, start=-10, length=-1, period=-1, multiplier=0, default_time_unit='ms')
        elif protocol==None:            
            self.model.protocol = protocol_temp        
        
        sol = solve_ivp(self.model.response_diff_eq,  (0, pre_step), y0=self.model.y0, t_eval=t_eval,
                                dense_output=False, 
                                method=method, # RK45 | LSODA | DOP853 | Radau | BDF | RK23
                                max_step=max_step, atol=atol, rtol=rtol )        
            
        self.model.y0 = copy.copy(sol.y[:,-1])

        self.model.protocol = protocol_temp

        return sol.y[:,-1]
                
        
    def simulate(self, t_span : list, 
                       t_eval=None, log=None, method='BDF', max_step=np.inf, atol=1E-6, rtol=1E-3):
        '''
        method : RK45 | LSODA | DOP853 | Radau | BDF | RK23
        '''                            
        self.solution = solve_ivp(self.model.diff_eq_solve_ivp, t_span, y0=self.model.y0, t_eval=t_eval,
                                dense_output=False, 
                                method=method, # RK45 | LSODA | DOP853 | Radau | BDF | RK23
                                max_step=max_step, atol=atol, rtol=rtol )    
        
        self.model.current_response_info = mod_trace.CurrentResponseInfo()
        list(map(self.model.differential_eq, self.solution.t, self.solution.y.transpose()))        
        
        self.model.set_result(self.solution.t, self.solution.y, log)
        return self.solution


    def pre_simulate2(self, pre_step=5000, protocol=None, 
                           method='BDF', max_step=np.inf, atol=1E-6, rtol=1E-3, t_eval=None):
        '''
        ''' 
        protocol_temp = copy.copy(self.model.protocol)
        self.model.protocol = protocol           
        if protocol == 'constant':            
            self.model.protocol = protocol_lib.VoltageClampProtocol( [protocol_lib.VoltageClampStep(voltage=self.model.y0[0], duration=pre_step)] )
        elif protocol =='pacing':
            self.model.protocol = protocol_lib.PacingProtocol(level=-1, start=-10, length=-1, period=-1, multiplier=0, default_time_unit='ms')
        elif protocol==None:            
            self.model.protocol = protocol_temp        
        
        times = np.linspace(0, pre_step, pre_step+10) 
        sol = odeint(self.model.diff_eq_odeint, t=times, y0=self.model.y0, hmax=max_step, rtol=rtol, atol=atol )   

        self.model.y0 = copy.copy(sol[-1, :])

        self.model.protocol = protocol_temp

        return sol[-1, :]


    def simulate2(self, times, log=None, max_step=np.inf, atol=1E-6, rtol=1E-3):
        '''
        
        '''            
        self.solution = odeint(self.model.diff_eq_odeint, t=times, y0=self.model.y0, hmax=max_step, rtol=rtol, atol=atol )      
        self.solution = self.solution.transpose(1,0)   
        self.model.set_result(times, self.solution, log) 
        return self.solution  
        


        
if __name__=='__main__':
    
    start_time = time.time()  
   

    print("--- %s seconds ---"%(time.time()-start_time))