"""Contains protocols to act in silico to probe cellular mechanics."""

import bisect
from typing import List, Union
import random
from math import floor
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from protocol_lib import VoltageClampProtocol
from protocol_lib import VoltageClampStep
from protocol_lib import VoltageClampRamp

def hERG_CiPA():    
    steps = []
    steps.append( VoltageClampStep(voltage=-80, duration=100) )
    steps.append( VoltageClampStep(voltage=-90, duration=100) )
    steps.append( VoltageClampStep(voltage=-80, duration=100) )
    steps.append( VoltageClampStep(voltage=40, duration=500) )
    steps.append( VoltageClampRamp(voltage_start=40, voltage_end=-80, duration=100)) # ramp step
    steps.append( VoltageClampStep(voltage=-80, duration=3000) )
    return VoltageClampProtocol(steps)  # steps=steps


def cav12_CiPA():    
    steps = []
    steps.append( VoltageClampStep(voltage=-80, duration=100) )
    steps.append( VoltageClampStep(voltage=-90, duration=100) )
    steps.append( VoltageClampStep(voltage=-80, duration=100) )
    steps.append( VoltageClampStep(voltage=0, duration=40) )
    steps.append( VoltageClampStep(voltage=30, duration=200) )
    steps.append( VoltageClampRamp(voltage_start=30, voltage_end=-80, duration=100)) # ramp step
    steps.append( VoltageClampStep(voltage=-80, duration=3000) )
    return VoltageClampProtocol(steps)  # steps=steps


def lateNav15_CiPA():    
    steps = []
    steps.append( VoltageClampStep(voltage=-95, duration=50) )
    steps.append( VoltageClampStep(voltage=-120, duration=200) )
    steps.append( VoltageClampStep(voltage=-15, duration=40) )
    steps.append( VoltageClampStep(voltage=40, duration=200) )    
    steps.append( VoltageClampRamp(voltage_start=40, voltage_end=-95, duration=100)) # ramp step
    steps.append( VoltageClampStep(voltage=-95, duration=3000) )
    return VoltageClampProtocol(steps)  # steps=steps


def leak_staircase():    
    tpre  = 0.2           # Time before step to variable V
    tstep = 0.5           # Time at variable V
    tpost = 0.1           # Time after step to variable V
    vhold = -80e-3
    vmin = -60e-3#-100e-3
    vmax = 40e-3
    vres = 20e-3        # Difference in V between steps
    v = np.arange(vmin, vmax + vres, vres)

    VC_protocol = VoltageClampProtocol()  # steps=steps

    # Leak estimate
    VC_protocol.add( VoltageClampStep(voltage=vhold, duration=0.25) )
    VC_protocol.add( VoltageClampStep(voltage=-120e-3, duration=0.05) )
    VC_protocol.add( VoltageClampRamp(voltage_start=-0.12, voltage_end=-0.08, duration=0.4)) # ramp step
    # Staircase
    VC_protocol.add( VoltageClampStep(voltage=vhold, duration=0.2) )
    VC_protocol.add( VoltageClampStep(voltage=40e-3, duration=1.0) )
    VC_protocol.add( VoltageClampStep(voltage=-120e-3, duration=0.5) )
    VC_protocol.add( VoltageClampStep(voltage=vhold, duration=1.0) )
    for vstep in v[1::]:
        VC_protocol.add( VoltageClampStep(voltage=vstep, duration=tstep) )
        VC_protocol.add( VoltageClampStep(voltage=vstep-vres, duration=tstep) )    
    for vstep in v[::-1][:-1]:
        VC_protocol.add( VoltageClampStep(voltage=vstep, duration=tstep) )
        VC_protocol.add( VoltageClampStep(voltage=vstep-2*vres, duration=tstep) )    
    VC_protocol.add( VoltageClampStep(voltage=vhold, duration=1.0 - tstep) )    # extend a bit the ending...
    # EK estimate
    VC_protocol.add( VoltageClampStep(voltage=40e-3, duration=tstep) )
    VC_protocol.add( VoltageClampStep(voltage=-70e-3, duration=10e-3) )  # Michael's suggestion
    VC_protocol.add( VoltageClampRamp(voltage_start=-70e-3, voltage_end=-110e-3, duration=0.1))  # second ramp step  
    VC_protocol.add( VoltageClampStep(voltage=-120e-3, duration=tstep-110e-3) )  # 
    VC_protocol.add( VoltageClampStep(voltage=vhold, duration=100) )
    
    return VC_protocol


