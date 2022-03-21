import os, sys
import time, glob

import random
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import ode, solve_ivp, odeint
from scipy.optimize import curve_fit, least_squares

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
        self.V = -80.0
        self.dt = 0.01
        self.record_time_step = 1

    # def set_times(self, times):
    #     self.times = times
    #     print("Times has been set.")

    def cal_dt(self, max_step : float) -> float:        
        if self.dt>max_step:
            self.dt = max_step
        return self.dt    
    
    def simulate(self, end_time : float,
                 log=[], max_step=float('inf'), default_time_unit='ms'):
        '''
        '''                               
        # self.nt = end_time/dt
                   
        # initial values
        nIter = 0
        current_time = 0          
        current_y = np.array(self.model.y0)

        self.times = [current_time]
        self.y_li = [self.model.y0]
        
        while current_time<=end_time:    
                        
            # integration
            df = self.model.response_diff_eq(current_time, current_y)
            
            # calculate time step
            dt = self.cal_dt(max_step)
            
            next_y = current_y + dt*np.array(df)
            
            # update values            
            current_y = next_y
            current_time += dt 
            nIter += 1

        #     if current_time
            self.times.append(current_time)
            self.y_li.append(current_y)  

        self.model.set_result(np.array(self.times), np.array(self.y_li).T, log)


        
if __name__=='__main__':
    
    start_time = time.time()
   
   

    print("--- %s seconds ---"%(time.time()-start_time))