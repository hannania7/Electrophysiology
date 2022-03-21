"""Contains classes to store the result of a genetic algorithm run.

Additionally, the classes in this module allow for figure generation.
"""

from abc import ABC
import copy
import enum
import math
import random
from typing import Dict, List, Union
from os import listdir, mkdir

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

import ga_configs
import mod_protocols as protocols
import mod_trace as trace
import mod_kernik as kernik


#############################################
from scipy.integrate import ode, solve_ivp
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import pickle
import bisect

sys.path.append('../')
sys.path.append('../Protocols')
import protocol_lib
import simulator_scipy
import simulator_myokit
import myokit
from Models.br1977 import BR1977
from Models.ord2011 import ORD2011
import mod_trace


def get_model_response_JK( model, protocol, prestep=None):    
    
    model.cell.mode = 1
    simulator = simulator_scipy.Simulator(model)     
        
    if prestep == None:
        print("There is no pre-step simulation.")
    elif prestep == 5000:        
        y0 = [-8.00000003e+01,  6.94549002e+00,  6.94553614e+00,  1.44766826e+02,
                1.44766919e+02,  5.46283800e-05,  5.38550879e-05,  1.25377970e+00,
                1.25388392e+00,  1.63694063e-02,  3.83078124e-01,  3.83078124e-01,
                3.83078124e-01,  1.83137288e-01,  3.83078124e-01,  8.60298196e-04,
                2.65750243e-01,  1.36775744e-01,  1.71654793e-03,  9.98192733e-01,
                9.98192733e-01,  8.74934836e-04,  9.98192733e-01,  9.98192733e-01,
                1.55207580e-08,  9.99999920e-01,  9.99999921e-01,  9.99999920e-01,
                9.99999920e-01,  9.99999920e-01,  4.72523502e-04,  9.99999920e-01,
                9.99999920e-01,  2.60425715e-05,  2.54957029e-05,  4.27866636e-04,
                4.72094402e-04,  9.98307893e-01,  6.06464770e-07,  7.58083578e-07,
                2.45432407e-04]
        simulator.model.y0 = y0        
    else:        
        simulator.pre_simulate( pre_step=prestep, protocol='constant')    
    solution = simulator.simulate( [0, protocol.get_voltage_change_endpoints()[-1]], method='BDF', max_step=1, atol=1e-06, rtol=1e-6)     
    command_voltages = [protocol.get_voltage_at_time(t) for t in solution.t]    

    tr = trace.Trace(protocol,
                     cell_params=None,
                     t=solution.t,
                     y=command_voltages,  # simulator.model.V,
                     command_voltages=command_voltages,
                     current_response_info=simulator.model.current_response_info,
                     default_unit=None)        
    # print(solution)
    return tr


def get_model_response_with_myokit( simulator, protocol, prestep=None):    
    
    model, p, s = myokit.load( "../mmt-model-files/ohara-cipa-v1-2017_VC.mmt" )    
                      
    simulator = simulator_myokit.Simulator(model, protocol, max_step=1.0, abs_tol=1e-8, rel_tol=1e-8, vhold=-80) # 1e-12, 1e-14 # 1e-08, 1e-10  # max_step=1, atol=1E-2, rtol=1E-4 # defalt: abs_tol=1e-06, rel_tol=0.0001    
    # simulator.reset_simulation_with_new_protocol( protocol )
    simulator.simulation.set_constant('cell.mode', 1)  
    
    if prestep == None:
        print("There is no pre-step simulation.")
    elif prestep == 15000:        
        y0 = [-8.69999996e+01,  6.94732336e+00,  6.94736848e+00,  1.44992431e+02,
               1.44992434e+02,  5.48328391e-05,  5.40431668e-05,  1.25617506e+00,
               1.25618638e+00,  8.12231733e-03,  6.62326077e-01,  6.62326077e-01,
               6.62326077e-01,  4.14582271e-01,  6.62326077e-01,  2.27721811e-04,
               4.79645030e-01,  2.87189165e-01,  1.07103663e-03,  9.99468797e-01,
               9.99468797e-01,  5.45740810e-04,  9.99468797e-01,  9.99468797e-01,
               2.96634937e-09,  9.99999988e-01,  9.99999988e-01,  9.99999988e-01,
               9.99999988e-01,  9.99999988e-01,  4.78979614e-04,  9.99999988e-01,
               9.99999988e-01,  9.28750206e-06,  9.23466020e-06,  1.96054631e-04,
               2.15667189e-04,  9.97012407e-01,  1.27419629e-07,  1.59274616e-07,
               2.47073549e-04]
        simulator.set_initial_values(y0)      
    else:
        simulator.pre_simulate(pre_step=prestep, sim_type=1)
    
    d = simulator.simulate(protocol.get_voltage_change_endpoints()[-1], log_times=None, extra_log=['ina.INa', 'inal.INaL', 'ito.Ito', 'ical.ICaL', 'ikr.IKr', 'iks.IKs', 'ik1.IK1'])     
    times = d['engine.time']    
    command_voltages = [protocol.get_voltage_at_time(t) for t in times]    
    tr = trace.Trace(protocol,
                     cell_params=None,
                     t=times,
                     y=command_voltages,  # simulator.model.V,
                     command_voltages=command_voltages,
                     current_response_info=simulator.current_response_info,
                     default_unit=None)         
    return tr