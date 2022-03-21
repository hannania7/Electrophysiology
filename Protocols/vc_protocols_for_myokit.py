"""Contains protocols to act in silico to probe cellular mechanics."""

import bisect
from typing import List, Union
import random
from math import floor
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import myokit
import protocol_lib

def leak_staircase(model, return_capmask=False):
    # My 'test6_v3'/staircase-ramp protocol
    # model: myokit model
    # return_capmask: if True, return an extra function that takes time series
    #                 as argument and return a mask to filter off capacitance
    #                 effect.
    tpre  = 0.2           # Time before step to variable V
    tstep = 0.5           # Time at variable V
    tpost = 0.1           # Time after step to variable V
    vhold = -80e-3
    vmin = -60e-3#-100e-3
    vmax = 40e-3
    vres = 20e-3        # Difference in V between steps
    v = np.arange(vmin, vmax + vres, vres)

    steps = []
    # Leak estimate
    steps += [(vhold, 0.25)]
    steps += [(-120e-3, 0.05)]
    steps += [(-30e-3, 400e-3)]  # ramp step
    # Staircase
    steps += [(vhold, 0.2)]
    steps += [(40e-3, 1.0)]
    steps += [(-120e-3, 0.5)]
    steps += [(vhold, 1.0)]
    for vstep in v[1::]:
        steps += [(vstep, tstep)]
        steps += [(vstep-vres, tstep)]
    for vstep in v[::-1][:-1]:
        steps += [(vstep, tstep)]
        steps += [(vstep-2*vres, tstep)]
    steps += [(vhold, 1.0 - tstep)]  # extend a bit the ending...
    # EK estimate
    steps += [(40e-3, tstep)]
    steps += [(-70e-3, 10e-3)]  # Michael's suggestion
    steps += [(-120e-3, tstep - 10e-3)]  # second ramp step
    steps += [(vhold, 100)]
    # Set ramp bit
    
    model.get('membrane.V').set_rhs(
                'piecewise('
                +
                'engine.time >= 0.300 and engine.time < 0.700001,'
                + '-150e-3 + 0.1 * engine.time'
                +
                ', engine.time >= 14.410 and engine.time < 14.510001,'
                + ' + 5.694 - 0.4 * engine.time'
                +
                ', engine.pace)')

    capmask = None
    if return_capmask:

        def capmask(times, capmaskdt=capmaskdt):
            fcap = np.ones(times.shape)
            currentt = 0
            for v, dur in steps:
                idxi = np.where(times > currentt)[0][0] - 1  # inclusive
                idxf = np.where(times > currentt + capmaskdt)[0][0]
                fcap[idxi:idxf] = 0
                currentt += dur
            return fcap
        
    protocol = myokit.Protocol()
    for f, t in steps:
        protocol.add_step(f, t)

    return model, protocol, capmask




