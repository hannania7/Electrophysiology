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
from numba import jit
# import pickle
# import bisect
sys.path.append('../')
sys.path.append('../Lib')
sys.path.append('../Protocols')
from pacing_protocol import PacingProtocol
import mod_trace as trace

from numbalsoda import lsoda_sig, lsoda
from numba import njit, cfunc
        
'''
'Endocardial' : 0,
'Epicardial' : 1,
'Mid-myocardial' : 2,


This version have conductance adjusment values such as GNa_adj, GNaL_adj, Gto_adj
and does not have cell_mode.
'''        


        
class ORD2011():
    """    
    O'Hara-Rudy CiPA v1.0 (2011)
    """
    def __init__(self):
        
        self.name = "ORD2011"    
        self.y0 = [1,2,2,2,2]

def ORD2011_initial_values():
    """    
    O'Hara-Rudy CiPA v1.0 (2011)
    """
    V = -87
    return [V]
   
# def ORD2011_AP(t, y):     odeint
@cfunc(lsoda_sig)
def ORD2011_AP(t, y, dy, p):    
# @njit
# def ORD2011_fc(y, t, *p):    
    # dy = np.zeros(41)
    # Protocol ==============================================================================================
    level=1
    start=20
    length=0.5
    period=1000
    multiplier=0
  
    default_time_unit='ms'
    time_conversion = 1000.0
    if default_time_unit == 's':
        time_conversion = 1.0
        default_unit = 'standard'
    else:
        time_conversion = 1000.0
        default_unit = 'milli'
        
    # stim_amplitude = protocol.stim_amplitude * 1E-3 * self._time_conversion
    stim_start = start * 1E-3 * time_conversion
    stim_duration = length * 1E-3 * time_conversion
#         stim_end = self._end * 1E-3 * self._time_conversion
    i_stim_period = period * 1E-3 *time_conversion 

    if time_conversion == 1:
        denom = 1E9
    else:
        denom = 1

    pace = (1.0 if t - stim_start - i_stim_period*floor((t - stim_start)/i_stim_period) <= stim_duration \
#                 and time <= stim_end \
            and t >= stim_start else 0) / denom
    
    # Model ==============================================================================================    
    V = y[0]
    Nai = y[1]
    Na_ss = y[2]
    Ki = y[3]
    K_ss = y[4]
    Cai = y[5]
    cass = y[6]
    Ca_nsr = y[7]
    Ca_jsr = y[8]
    m = y[9]
    hf = y[10]
    hs = y[11]
    j = y[12]
    hsp = y[13]
    jp = y[14]
    mL = y[15]
    hL = y[16]
    hLp = y[17]
    a = y[18]
    iF = y[19]
    iS = y[20]
    ap = y[21]
    iFp = y[22]
    iSp = y[23]
    d = y[24]
    ff = y[25]
    fs = y[26]
    fcaf = y[27]
    fcas = y[28]
    jca = y[29]
    nca = y[30]
    ffp = y[31]
    fcafp = y[32]
    xf = y[33]
    xs = y[34]
    xs1 = y[35]
    xs2 = y[36]
    xk1 = y[37]
    Jrelnp = y[38]
    Jrelp = y[39]
    CaMKt = y[40]

    # Stimulus
    amplitude = -80 # [uA/cm^2]
    I_stim = amplitude * pace     
    
    # Phys
    R = 8314           # [J/kmol/K] : Gas constant
    T = 310            # [K] : Temperature
    F = 96485          # [C/mol] : Faraday's constant
    RTF  = R*T/F
    FRT  = F/(R*T)
    FFRT = F*F/(R*T)  
    zna = 1
    zca = 2
    zk = 1
                 
    # Cell
    mode = 1  # The type of cell. Endo=0, Epi=1, Mid=2
    L = 0.01  # [cm] Cell length
    rad = 0.0011  # [cm] cell radius
    vcell = 1000 * 3.14 * rad * rad * L # [uL] Cell volume
    Ageo = 2*3.14 * rad * rad + 2*3.14 * rad * L   # [cm^2] Geometric cell area
    Acap = 2 * Ageo            # [cm^2] Capacitative membrane area
    vmyo = 0.68 * vcell        # [uL] Volume of the cytosolic compartment
    vnsr = 0.0552 * vcell      # [uL] Volume of the NSR compartment
    vjsr = 0.0048 * vcell      # [uL] Volume of the JSR compartment
    vss = 0.02 * vcell         # [uL] Volume of the Submembrane space near the T-tubules
    AF = Acap / F         # F : Faraday's constant

    # Extra
    Nao = 140 # [mmol/L] : Extracellular Na+ concentration
    Cao = 1.8 # [mmol/L] : Extracellular Ca2+ concentration
    Ko  = 5.4 # [mmol/L] : Extracellular K+ concentration

    # Nernst
    PKNa = 0.01833          # desc: Permeability ratio K+ to Na+
    ENa = RTF * log(Nao / Nai)      # in [mV]  desc: Reversal potential for Sodium currents
    EK = RTF * log(Ko / Ki)      # in [mV]  desc: Reversal potential for Potassium currents                
    EKs = RTF * log((Ko + PKNa * Nao) / (Ki + PKNa * Nai)) # desc: Reversal potential for IKs  in [mV]

    # CaMKt
    KmCaMK = 0.15
    aCaMK  = 0.05
    bCaMK  = 0.00068
    CaMKo  = 0.05
    KmCaM  = 0.0015
    CaMKb  = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass)
    CaMKa  = CaMKb + CaMKt
    d_CaMKt = aCaMK * CaMKb * CaMKa - bCaMK * CaMKt
    camk_f = 1 / (1 + KmCaMK / CaMKa) # Fraction of phosphorylated channels        

    dy[40] = d_CaMKt


    # INa =========================================================================
    #: Maximum conductance of INa channels
    GNa_max = 75.0         # 2011 : 75
    GNa_adj = 1.0   
    mtD1 = 6.765
    mtD2 = 8.552
    mtV1 = 11.64
    mtV2 = 34.77
    mtV3 = 77.42
    mtV4 = 5.955            
    tm  = 1.0 / (mtD1 * exp((V + mtV1) / mtV2) + mtD2 * exp(-(V + mtV3) / mtV4)) # desc: Time constant for m-gate   in [ms]
    mssV1 = 39.57
    mssV2 = 9.871
    mss  = 1.0 / (1.0 + exp(-(V + mssV1)/mssV2))  # desc: Steady state value for m-gate         
    d_m = (mss - m) / tm           

    hssV1 = 82.9 
    hssV2 = 6.086 
    shift_INa_inact = 0.0
    hss = 1.0 / (1.0 + exp((V + hssV1-shift_INa_inact) / hssV2))   # desc: Steady-state value for h-gate
    thf = 1.0 / (1.432e-5 * exp(-(V + 1.196 - shift_INa_inact) / 6.285) + 6.1490 * exp((V + 0.5096 - shift_INa_inact) / 20.27)) # desc: Time constant for fast development of inactivation in INa   in [ms]
    ths = 1.0 / (0.009794 * exp(-(V + 17.95-shift_INa_inact) / 28.05) + 0.3343 * exp((V + 5.7300 - shift_INa_inact) / 56.66))  # desc: Time constant for slow development of inactivation in INa  in [ms]
    Ahf = 0.99 # : Fraction of INa channels with fast inactivation
    Ahs = 1.0 - Ahf # : Fraction of INa channels with slow inactivation        
    d_hf = (hss - hf) / thf   # desc: Fast componennt of the inactivation gate for INa channels
    d_hs = (hss - hs) / ths   # desc: Slow componennt of the inactivation gate for non-phosphorylated INa channels
    
    h = Ahf * hf + Ahs * hs   # desc: Inactivation gate for INa
    tj = 2.038 + 1 / (0.02136 * exp(-(V + 100.6 - shift_INa_inact) / 8.281) + 0.3052 * exp((V + 0.9941 - shift_INa_inact) / 38.45)) # desc: Time constant for j-gate in INa in [ms]
    jss = hss  # desc: Steady-state value for j-gate in INa
    d_j = (jss - j) / tj # desc: Recovery from inactivation gate for non-phosphorylated INa channels

    # Phosphorylated channels
    thsp = 3 * ths # desc: Time constant for h-gate of phosphorylated INa channels  in [ms]
    hssp = 1 / (1 + exp((V + 89.1 - shift_INa_inact) / 6.086)) # desc: Steady-state value for h-gate of phosphorylated INa channels        
    d_hsp = (hssp - hsp) / thsp # desc: Slow componennt of the inactivation gate for phosphorylated INa channels

    hp = Ahf * hf + Ahs * hsp # desc: Inactivation gate for phosphorylated INa channels
    tjp = 1.46 * tj # desc: Time constant for the j-gate of phosphorylated INa channels     in [ms]        
    d_jp = (jss - jp) / tjp # desc: Recovery from inactivation gate for phosphorylated INa channels
    
    # Current       
    GNa = GNa_max * GNa_adj         
    INa = GNa * (V - ENa) * m**3 * ((1.0 - camk_f) * h * j + camk_f * hp * jp) # in [uA/uF]  desc: Fast sodium current

    dy[9] = d_m
    dy[10] = d_hf 
    dy[11] = d_hs
    dy[12] = d_j
    dy[13] = d_hsp
    dy[14] = d_jp


    # INaL =========================================================================
    #: Maximum conductance           Endo = 0, Epi = 1, Mid = 2
    GNaL_max = 0.0075*0.6   # 'Endocardial' : 0.0075,  'Epicardial' : 0.0075*0.6,   'Mid-myocardial' : 0.0075*1.0,
    GNaL_adj = 1.0              

    mLss = 1.0 / (1.0 + exp(-(V + 42.85) / 5.264)) # desc: Steady state value of m-gate for INaL
    tmL = tm
    d_mL = (mLss - mL) / tmL # desc: Activation gate for INaL

    thL = 200.0 # [ms] : Time constant for inactivation of non-phosphorylated INaL channels
    hLss = 1.0 / (1.0 + exp((V + 87.61) / 7.488))  # desc: Steady-state value for inactivation of non-phosphorylated INaL channels
    
    d_hL = (hLss - hL) / thL # desc: Inactivation gate for non-phosphorylated INaL channels

    hLssp = 1.0 / (1.0 + exp((V + 93.81) / 7.488)) # desc: Steady state value for inactivation of phosphorylated INaL channels
    thLp = 3 * thL # in [ms] desc: Time constant for inactivation of phosphorylated INaL channels

    d_hLp = (hLssp - hLp) / thLp  # desc: Inactivation gate for phosphorylated INaL channels
            
    # Current        
    GNaL = GNaL_max * GNaL_adj
    INaL = GNaL * (V - ENa) * mL * ((1 - camk_f) * hL + camk_f * hLp)
 
    dy[15] = d_mL
    dy[16] = d_hL 
    dy[17] = d_hLp
 

    # Ito ============================================================================================
    #: Maximum conductance
    Gto_max = 0.08   # 'Endocardial' : 0.02,  'Epicardial' : 0.08,   'Mid-myocardial' : 0.08
    Gto_adj = 1.0   
    
    ass = 1.0 / (1.0 + exp(-(V - 14.34) / 14.82))  # desc: Steady-state value for Ito activation
    one = 1.0 / (1.2089 * (1 + exp(-(V - 18.4099) / 29.3814)))
    two = 3.5 / (1 + exp((V + 100) / 29.3814))
    ta = 1.0515 / (one + two)  # desc: Time constant for Ito activation  in [ms]        
    d_a = (ass - a) / ta   # desc: Ito activation gate
    
    iss = 1.0 / (1.0 + exp((V + 43.94) / 5.711))   # desc: Steady-state value for Ito inactivation
    delta_epi = 1.0
    if mode==1:  delta_epi = 1.0 - (0.95 / (1 + exp((V + 70.0) / 5.0)))   # desc: Adjustment for different cell types
    tiF_b = (4.562 + 1 / (0.3933 * exp(-(V+100) / 100) + 0.08004 * exp((V + 50) / 16.59)))  # desc: Time constant for fast component of Ito inactivation     in [ms]
    tiS_b = (23.62 + 1 / (0.001416 * exp(-(V + 96.52) / 59.05) + 1.780e-8 * exp((V + 114.1) / 8.079)))   # desc: Time constant for slow component of Ito inactivation  in [ms]
    tiF = tiF_b * delta_epi
    tiS = tiS_b * delta_epi
    AiF = 1.0 / (1.0 + exp((V - 213.6) / 151.2))  # desc: Fraction of fast inactivating Ito channels
    AiS = 1.0 - AiF        
    d_iF = (iss - iF) / tiF  # desc: Fast component of Ito activation        
    d_iS = (iss - iS) / tiS  # desc: Slow component of Ito activation
    
    i = AiF * iF + AiS * iS     # desc: Inactivation gate for non-phosphorylated Ito
    assp=1.0/(1.0+exp(-(V-24.34)/14.82))        
    d_ap = (assp - ap) / ta
        
    dti_develop = 1.354 + 1e-4 / (exp((V - 167.4) / 15.89) + exp(-(V - 12.23) / 0.2154))
    dti_recover = 1 - 0.5 / (1 + exp((V+70) / 20))
    tiFp = dti_develop * dti_recover * tiF # desc: Time constant for fast component of inactivation of phosphorylated Ito channels   in [ms]
    tiSp = dti_develop * dti_recover * tiS  # desc: Time constant for slot component of inactivation of phosphorylated Ito channels  in [ms]        
    d_iFp = (iss - iFp) / tiFp # desc: Fast component of inactivation of phosphorylated Ito channels        
    d_iSp = (iss - iSp) / tiSp # desc: Slow component of inactivation of phosphorylated Ito channels
    
    ip = AiF * iFp + AiS * iSp  # desc: Inactivation gate for phosphorylated Ito channels
    
    # Current        
    Gto = Gto_max * Gto_adj                
    Ito = Gto * (V - EK) * ((1 - camk_f) * a * i + camk_f * ap * ip) # desc: Transient outward Potassium current
    
    dy[18] = d_a
    dy[19] = d_iF 
    dy[20] = d_iS
    dy[21] = d_ap
    dy[22] = d_iFp 
    dy[23] = d_iSp


    # ICaL ========================================================================
    #: Maximum conductance
    PCa_max = 1.2*0.0001  # 'Endocardial' : 0.0001  'Epicardial' : 0.0001*1.2,   'Mid-myocardial' : 0.0001*2.5
    PCa_adj = 1.0
    
    vfrt = V * FRT
    vffrt = V * FFRT
    
    # Activation
    dss = 1.0 / (1.0 + exp(-(V + 3.94) / 4.23)) # Steady-state value for activation gate of ICaL channel
    td = 0.6 + 1.0 / (exp(-0.05 * (V + 6)) + exp(0.09 * (V + 14)))  # Time constant for activation gate of ICaL channel   in [ms]
    d_d = (dss - d) / td   # Activation gate of ICaL channel
    
    # Inactivation
    fss = 1.0 / (1.0 + exp((V + 19.58) / 3.696)) # Steady-state value for inactivation gate of ICaL channel
    tff = 7.0 + 1.0 / (0.0045 * exp(-(V + 20) / 10) + 0.0045 * exp((V + 20) / 10))  # Time constant for fast inactivation of ICaL channels  in [ms]
    tfs = 1000 + 1.0 / (0.000035 * exp(-(V + 5) / 4) + 0.000035 * exp((V + 5) / 6)) # Time constant for fast inactivation of ICaL channels  in [ms]
    Aff = 0.6  # Fraction of ICaL channels with fast inactivation
    Afs = 1.0 - Aff # Fraction of ICaL channels with slow inactivation
    d_ff = (fss - ff) / tff   # Fast inactivation of ICaL channels
    d_fs = (fss - fs) / tfs   # Slow inactivation of ICaL channels        
    f = Aff * ff + Afs * fs  # Inactivation of ICaL channels
    
    # Ca-dependent inactivation
    fcass = fss  # Steady-state value for Ca-dependent inactivation of ICaL channels
    tfcaf = 7.0 + 1.0 / (0.04 * exp(-(V - 4.0) / 7.0) + 0.04 * exp((V - 4.0) / 7.0))  # Time constant for fast Ca-dependent inactivation of ICaL channels in [ms]
    tfcas = 100.0 + 1 / (0.00012 * exp(-V / 3) + 0.00012 * exp(V / 7))  # Time constant for slow Ca-dependent inactivation of ICaL channels  in [ms]
    Afcaf = 0.3 + 0.6 / (1 + exp((V - 10) / 10))  # Fraction of ICaL channels with fast Ca-dependent inactivation
    Afcas = 1.0 - Afcaf    # Fraction of ICaL channels with slow Ca-dependent inactivation
    d_fcaf = (fcass - fcaf) / tfcaf   # Fast Ca-dependent inactivation of ICaL channels
    d_fcas = (fcass - fcas) / tfcas   # Slow Ca-dependent inactivation of ICaL channels
    fca = Afcaf * fcaf + Afcas * fcas    # Ca-dependent inactivation of ICaL channels
    
    # Recovery from Ca-dependent inactivation
    tjca = 75.0  # [ms] : Time constant of recovery from Ca-dependent inactivation
    d_jca = (fcass - jca) / tjca  # Recovery from Ca-dependent inactivation
    
    # Inactivation of phosphorylated channels
    tffp = 2.5 * tff   # in [ms]  desc: Time constant for fast inactivation of phosphorylated ICaL channels
    d_ffp = (fss - ffp) / tffp   # Fast inactivation of phosphorylated ICaL channels
    fp = Aff * ffp + Afs * fs    # Inactivation of phosphorylated ICaL channels
    
    # Ca-dependent inactivation of phosphorylated channels
    tfcafp = 2.5 * tfcaf   # in [ms]  desc: Time constant for fast Ca-dependent inactivation of phosphorylated ICaL channels
    d_fcafp = (fcass - fcafp) / tfcafp  # Fast Ca-dependent inactivation of phosphorylated ICaL channels
    fcap = Afcaf * fcafp + Afcas * fcas   # Ca-dependent inactivation of phosphorylated ICaL channels
            
    Kmn = 0.002
    k2n = 1000
    km2n = jca * 1.0
    anca = 1.0 / (k2n / km2n + (1 + Kmn / cass)**4.0)  # Fraction of channels in Ca-depdent inactivation mode
    d_nca = anca * k2n - nca*km2n   # Fraction of channels in Ca-depdent inactivation mode
    
    # Total currents through ICaL channel
    # v0 = 0
    # B_1 = 2.0 * FRT
    # A_1 = ( 4.0 * FFRT * (cass * exp(2 * vfrt) - 0.341 * Cao) ) / B_1
    # U_1 = B_1 * (V-v0)
    # PhiCaL = (A_1*U_1)/(exp(U_1)-1.0)
    # if -1e-7<=U_1 and U_1<=1e-7 :  PhiCaL = A_1 * (1.0-0.5*U_1)        
    # B_2 = FRT
    # A_2 = ( 0.75 * FFRT * (Na_ss * exp(vfrt) - 0.341 * Nao) ) / B_2
    # U_2 = B_2 * (V-v0)
    # PhiCaNa = (A_2*U_2)/(exp(U_2)-1.0)
    # if (-1e-7<=U_2 and U_2<=1e-7) :  PhiCaNa = A_2 * (1.0-0.5*U_2) 
    # B_3 = FRT
    # A_3 = ( 0.75 * FFRT * (K_ss * exp(vfrt) - 0.341 * Ko) ) / B_3
    # U_3 = B_3 * (V-v0)
    # PhiCaK = (A_3*U_3)/(exp(U_3)-1.0)
    # if (-1e-7<=U_3 and U_3<=1e-7) :  PhiCaK = A_3 * (1.0-0.5*U_3)

    PhiCaL  = 4 * vffrt *(       cass  * exp(2 * vfrt) - 0.341 * Cao) / (exp(2 * vfrt) - 1)
    PhiCaNa = 1 * vffrt *(0.75 * Na_ss   * exp(1 * vfrt) - 0.75  * Nao) / (exp(1 * vfrt) - 1)
    PhiCaK  = 1 * vffrt *(0.75 * K_ss * exp(1 * vfrt) - 0.75  * Ko ) / (exp(1 * vfrt) - 1)
            
    PCa = PCa_max * PCa_adj
                
    PCap   = 1.1      * PCa
    PCaNa  = 0.00125  * PCa
    PCaK   = 3.574e-4 * PCa
    PCaNap = 0.00125  * PCap
    PCaKp  = 3.574e-4 * PCap
    
    g  = d * (f  * (1.0 - nca) + jca * fca  * nca)   # Conductivity of non-phosphorylated ICaL channels
    gp = d * (fp * (1.0 - nca) + jca * fcap * nca)   # Conductivity of phosphorylated ICaL channels
    
    ICaL   = (1.0 - camk_f) * PCa   * PhiCaL  * g + camk_f * PCap   * PhiCaL  * gp  # L-type Calcium current   in [uA/uF]
    ICaNa  = (1.0 - camk_f) * PCaNa * PhiCaNa * g + camk_f * PCaNap * PhiCaNa * gp   # Sodium current through ICaL channels  in [uA/uF]
    ICaK   = (1.0 - camk_f) * PCaK  * PhiCaK  * g + camk_f * PCaKp  * PhiCaK  * gp   # Potassium current through ICaL channels  in [uA/uF]
    
    dy[24] = d_d
    dy[25] = d_ff 
    dy[26] = d_fs
    dy[27] = d_fcaf
    dy[28] = d_fcas 
    dy[29] = d_jca
    dy[30] = d_nca 
    dy[31] = d_ffp
    dy[32] = d_fcafp


    # IKr
    #: Maximum conductance
    GKr_max = 0.046*1.3   # 'Endocardial' : 0.046  'Epicardial' : 0.046*1.3,   'Mid-myocardial' : 0.046*0.8
    GKr_adj = 1.0             
    # Activation
    sx = 1.0 / (1.0 + exp((V + 8.337) / -6.789))  # desc: Steady-state value for IKr activation
    txf = 12.98 + 1.0 / (0.36520 * exp((V - 31.66) / 3.869) + 4.123e-5 * exp((V - 47.78) / -20.38))  # Time constant for fast IKr activation
    txs = 1.865 + 1.0 / (0.06629 * exp((V - 34.70) / 7.355) + 1.128e-5 * exp((V - 29.74) / -25.94))  # Time constant for slow IKr activation
    d_xf = (sx - xf) / txf   # Fast activation of IKr channels
    d_xs = (sx - xs) / txs  # Slow activation of IKr channels
    Axf = 1.0 / (1.0 + exp((V + 54.81) / 38.21)) # Fraction of IKr channels with fast activation
    Axs = 1.0 - Axf   # Fraction of IKr channels with slow activation
    x = Axf * xf + Axs * xs   # Activation of IKr channels
    # Inactivation
    r = 1.0 / (1.0 + exp((V + 55.0) / 75.0)) * 1.0 / (1.0 + exp((V - 10.0) / 30.0))  # Inactivation of IKr channels
    # Current        
    GKr = GKr_max * GKr_adj
    IKr = GKr * sqrt(Ko / 5.4) * x * r * (V - EK) # Rapid delayed Potassium current  in [uA/uF]

    dy[33] = d_xf 
    dy[34] = d_xs    
    

    # Iks
    #: Maximum conductance
    GKs_max = 0.0034*1.4  # 'Endocardial' : 0.0034  'Epicardial' : 0.0034*1.4,   'Mid-myocardial' : 0.0034
    GKs_adj = 1.0
    
    xs1ss  = 1.0 / (1.0 + exp(-(V + 11.60) / 8.932)) # desc: Steady-state value for activation of IKs channels
    txs1_max = 817.3
    txs1 = txs1_max + 1.0 / (2.326e-4 * exp((V + 48.28) / 17.80) + 0.001292 * exp(-(V + 210) / 230)) # desc: Time constant for slow, low voltage IKs activation
                    
    d_xs1 = (xs1ss - xs1) / txs1  # desc: Slow, low voltage IKs activation
    
    xs2ss = xs1ss
    txs2 = 1.0 / (0.01 * exp((V - 50) / 20) + 0.0193 * exp(-(V + 66.54) / 31.0)) # desc: Time constant for fast, high voltage IKs activation
    
    d_xs2 = (xs2ss - xs2) / txs2   # desc: Fast, high voltage IKs activation
    
    KsCa = 1.0 + 0.6 / (1.0 + (3.8e-5 / Cai)**1.4) # desc: Maximum conductance for IKs
    
    GKs = GKs_max * GKs_adj        
    IKs = GKs * KsCa * xs1 * xs2 * (V - EKs)  # Slow delayed rectifier Potassium current
        
    dy[35] = d_xs1 
    dy[36] = d_xs2   


    # IK1
    #: Maximum conductance
    GK1_max = 0.1908*1.2 # 'Endocardial' : 0.1908  'Epicardial' : 0.1908*1.2,   'Mid-myocardial' : 0.1908*1.3
    GK1_adj = 1.0                    
    
    xk1ss = 1 / (1 + exp(-(V + 2.5538 * Ko + 144.59) / (1.5692 * Ko + 3.8115))) # Steady-state value for activation of IK1 channels  : sx in 2011
    txk1 = 122.2 / (exp(-(V + 127.2) / 20.36) + exp((V + 236.8) / 69.33))  # Time constant for activation of IK1 channels  : tx in 2011        
    d_xk1 = (xk1ss - xk1) / txk1  # Activation of IK1 channels        
    rk1 = 1.0 / (1.0 + exp((V + 105.8 - 2.6 * Ko) / 9.493))   # Inactivation of IK1 channels    : r in 2011            
    # desc: Conductivity of IK1 channels, cell-type dependent        
    GK1 = GK1_max * GK1_adj
    IK1 = GK1 * sqrt(Ko) * rk1 * xk1 * (V - EK)  # Inward rectifier Potassium current
    
    dy[37] = d_xk1 
    

    # INaCa =======================================================================
    kna1   = 15.0        
    kna2   = 5.0        
    kna3   = 88.12        
    kasymm = 12.5        
    wna    = 6.0e4
    wca    = 6.0e4
    wnaca  = 5.0e3 
    kcaon  = 1.5e6        
    kcaoff = 5.0e3        
    qna = 0.5224
    qca = 0.1670
    hca    = exp(qca * V * FRT)
    hna    = exp(qna * V * FRT)
    
    # Parameters h
    h1  = 1.0 + Nai / kna3 * (1 + hna)
    h2  = (Nai * hna) / (kna3 * h1)
    h3  = 1.0 / h1
    h4  = 1.0 + Nai / kna1 * (1 + Nai / kna2)
    h5  = Nai * Nai / (h4 * kna1 * kna2)
    h6  = 1.0 / h4
    h7  = 1.0 + Nao / kna3 * (1 + 1 / hna)
    h8  = Nao / (kna3 * hna * h7)
    h9  = 1.0 / h7
    h10 = kasymm + 1 + Nao / kna1 * (1 + Nao / kna2)
    h11 = Nao * Nao / (h10 * kna1 * kna2)
    h12 = 1.0 / h10
    
    # Parameters k
    k1   = h12 * Cao * kcaon
    k2   = kcaoff
    k3p  = h9 * wca
    k3pp = h8 * wnaca
    k3   = k3p + k3pp
    k4p  = h3 * wca / hca
    k4pp = h2 * wnaca
    k4   = k4p + k4pp
    k5   = kcaoff
    k6   = h6 * Cai * kcaon
    k7   = h5 * h2 * wna
    k8   = h8 * h11 * wna
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    KmCaAct = 150.0e-6
    allo    = 1 / (1 + (KmCaAct / Cai)**2.0)
    JncxNa  = 3 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
    JncxCa  = E2 * k2 - E1 * k1  
    Gncx_b = 0.0008        
    Gncx = Gncx_b*1.4
    if mode==0: Gncx = Gncx_b
    elif mode==1: Gncx = Gncx_b*1.1
                
    INaCa = 0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa)  # Sodium/Calcium exchange current  in [uA/uF]
    

    # INaCa_ss =================================================================================
    h1  = 1.0 + Na_ss / kna3 * (1 + hna)
    h2  = (Na_ss * hna)/(kna3 * h1)
    h3  = 1 / h1
    h4  = 1 + Na_ss / kna1 * (1 + Na_ss / kna2)
    h5  = Na_ss * Na_ss /(h4 * kna1 * kna2)
    h6  = 1 / h4
    h7  = 1 + Nao / kna3 * (1 + 1 / hna)
    h8  = Nao / (kna3 * hna * h7)
    h9  = 1/h7
    h10 = kasymm + 1 + Nao / kna1 * (1 + Nao / kna2)
    h11 = Nao * Nao / (h10 * kna1 * kna2)
    h12 = 1/h10
    k1   = h12 * Cao * kcaon
    k2   = kcaoff
    k3p  = h9 * wca
    k3pp = h8 * wnaca
    k3   = k3p + k3pp
    k4p  = h3 * wca / hca
    k4pp = h2 * wnaca
    k4   = k4p + k4pp
    k5   = kcaoff
    k6   = h6 * cass * kcaon
    k7   = h5 * h2 * wna
    k8   = h8 * h11 * wna
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    KmCaAct = 150.0e-6
    allo    = 1 / (1 + (KmCaAct / cass)**2)
    JncxNa  = 3 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
    JncxCa  = E2 * k2 - E1 * k1
    
    INaCa_ss = 0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa)  # Sodium/Calcium exchange current into the T-Tubule subspace  in [uA/uF]
    

    # INaK ========================================================
    k1p = 949.5
    k1m = 182.4
    k2p = 687.2
    k2m = 39.4
    k3p = 1899.0
    k3m = 79300.0
    k4p = 639.0
    k4m = 40.0
    Knai0 = 9.073
    Knao0 = 27.78
    delta = -0.1550
    Knai = Knai0 * exp(delta * V * FRT / 3.0)
    Knao = Knao0 * exp((1.0-delta) * V * FRT / 3.0)
    Kki    = 0.5
    Kko    = 0.3582
    MgADP  = 0.05
    MgATP  = 9.8
    Kmgatp = 1.698e-7
    H      = 1.0e-7
    eP     = 4.2
    Khp    = 1.698e-7
    Knap   = 224.0
    Kxkur  = 292.0
    P = eP / (1 + H / Khp + Nai / Knap + Ki / Kxkur)
    a1 = (k1p * (Nai / Knai)**3) / ((1 + Nai / Knai)**3 + (1 + Ki / Kki)**2 - 1)
    b1 = k1m * MgADP
    a2 = k2p
    b2 = (k2m * (Nao / Knao)**3) / ((1 + Nao / Knao)**3 + (1 + Ko / Kko)**2 - 1)
    a3 = (k3p * (Ko / Kko)**2  ) / ((1 + Nao / Knao)**3 + (1 + Ko / Kko)**2 - 1)
    b3 = (k3m * P * H)/(1 + MgATP / Kmgatp)
    a4 = (k4p * MgATP / Kmgatp) / (1 + MgATP / Kmgatp)
    b4 = (k4m * (Ki / Kki)**2) / ((1 + Nai / Knai)**3 + (1 + Ki / Kki)**2 - 1)
    x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2
    x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4
    x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1
    x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    JnakNa = 3 * (E1 * a3 - E2 * b3)
    JnakK  = 2 * (E4 * b1 - E3 * a1)
    Pnak_b = 30        
    Pnak = Pnak_b
    if mode==1: Pnak = Pnak_b*0.9
    elif mode==2: Pnak = Pnak_b*0.7
                
    INaK = Pnak * (zna * JnakNa + zk * JnakK)  # Sodium/Potassium ATPase current    in [uA/uF]
   

   # IKb ================================================    
    xkb = 1.0 / (1.0 + exp(-(V - 14.48) / 18.34))
    GKb_b = 0.003
    GKb = GKb_b
    if mode==1: GKb = GKb_b*0.6
                
    IKb = GKb * xkb * (V - EK)  # Background Potassium current   in [uA/uF]
    

    # INab =========================================
    B = FRT
    v0 = 0
    PNab = 3.75e-10
    A = PNab * FFRT * (Nai * exp(V * FRT) - Nao) / B            
    U = B * (V - v0)
    INab = (A*U)/(exp(U)-1.0)
    # if -1e-7<=U and U<=1e-7: INab = A*(1.0-0.5*U)   # Background Sodium current   in [uA/uF] <- 2017 version
               

    # ICab ====================================================
    B = 2 * FRT
    v0 = 0
    PCab = 2.5e-8
    A = PCab * 4.0 * FFRT * (Cai * exp( 2.0 * V * FRT) - 0.341 * Cao) / B
    U = B * (V - v0)
    ICab = (A*U)/(exp(U)-1.0)
    # if -1e-7<=U and U<=1e-7:  ICab = A*(1.0-0.5*U) # Background Calcium current  in [uA/uF] <- 2017 version
     
    # IpCa ==============================================
    GpCa = 0.0005
    KmCap = 0.0005
    IpCa = GpCa * Cai / (KmCap + Cai)  # Sarcolemmal Calcium pump current   in [uA/uF]

    # Ryr ================================================    
    bt=4.75
    a_rel=0.5*bt        
    Jrel_inf_temp = a_rel * -ICaL / (1 + (1.5 / Ca_jsr)**8)
    Jrel_inf = Jrel_inf_temp
    if (mode==2):   Jrel_inf = Jrel_inf_temp * 1.7      
    tau_rel_temp = bt / (1.0 + 0.0123 / Ca_jsr)
    tau_rel = tau_rel_temp
    if (tau_rel_temp < 0.001):   tau_rel = 0.001                        
    d_Jrelnp = (Jrel_inf - Jrelnp) / tau_rel   
    btp = 1.25*bt
    a_relp = 0.5*btp
    Jrel_temp = a_relp * -ICaL / (1 + (1.5 / Ca_jsr)**8)
    Jrel_infp = Jrel_temp
    if mode==2:   Jrel_infp = Jrel_temp * 1.7                
    tau_relp_temp = btp / (1 + 0.0123 / Ca_jsr)
    tau_relp = tau_relp_temp
    if tau_relp_temp < 0.001:   tau_relp = 0.001        
    d_Jrelp = (Jrel_infp - Jrelp) / tau_relp     
    Jrel_scaling_factor = 1.0
    Jrel = Jrel_scaling_factor * (1.0 - camk_f) * Jrelnp + camk_f * Jrelp # desc: SR Calcium release flux via Ryanodine receptor  in [mmol/L/ms]

    dy[38] = d_Jrelnp 
    dy[39] = d_Jrelp 
    
    # Diff ==========================================================
    JdiffNa = (Na_ss - Nai) / 2.0   # (sodium.Na_ss - sodium.Nai) / 2
    JdiffK  = (K_ss  - Ki)  / 2.0  # (potassium.K_ss  - potassium.Ki)  / 2
    Jdiff   = (cass - Cai) / 0.2     # (calcium.cass - calcium.Cai) / 0.2


    # Serca =====================================================
    upScale = 1.0
    if (mode == 1):   upScale = 1.3
    Jupnp = upScale * (0.004375 * Cai / (Cai + 0.00092))
    Jupp  = upScale * (2.75 * 0.004375 * Cai / (Cai + 0.00092 - 0.00017))        
    Jleak = 0.0039375 * Ca_nsr / 15.0 # in [mmol/L/ms]
    Jup_b = 1.0
    Jup = Jup_b * ((1.0 - camk_f) * Jupnp + camk_f * Jupp - Jleak) # desc: Total Ca2+ uptake, via SERCA pump, from myoplasm to nsr in [mmol/L/ms]        
    Jtr = (Ca_nsr - Ca_jsr) / 100.0   # desc: Ca2+ translocation from nsr to jsr    in [mmol/L/ms]
        
        
    # Sodium ========================================
    cm = 1.0
    INa_tot    = INa + INaL + INab + 3*INaCa + 3*INaK
    d_Nai   = -INa_tot * AF * cm / vmyo + JdiffNa * vss / vmyo # Intracellular Potassium concentration
    
    INa_ss_tot = ICaNa + 3*INaCa_ss
    d_Na_ss = -INa_ss_tot * AF * cm / vss - JdiffNa

    dy[1] = d_Nai 
    dy[2] = d_Na_ss 

        
    # Potassium ============================================
    cm = 1.0
    IK_tot = Ito + IKr + IKs + IK1 + IKb - 2 * INaK        
    d_Ki  = -(IK_tot + I_stim) * cm * AF / vmyo + JdiffK * vss / vmyo # Intracellular Potassium concentration
    
    IK_ss_tot = ICaK
    d_K_ss = -IK_ss_tot * cm * AF / vss - JdiffK  # Potassium concentration in the T-Tubule subspace
    
    dy[3] = d_Ki 
    dy[4] = d_K_ss 


    # Calcium ============================================
    cm = 1.0
    cmdnmax_b = 0.05
    cmdnmax = cmdnmax_b  
    if mode == 1:   cmdnmax = 1.3*cmdnmax_b        
    kmcmdn  = 0.00238
    trpnmax = 0.07
    kmtrpn  = 0.0005
    BSRmax  = 0.047
    KmBSR   = 0.00087
    BSLmax  = 1.124
    KmBSL   = 0.0087
    csqnmax = 10.0
    kmcsqn  = 0.8
    '''
    desc: Intracellular Calcium concentratium
    in [mmol/L]
    '''
    ICa_tot = IpCa + ICab - 2*INaCa 
    a = kmcmdn + Cai
    b = kmtrpn + Cai
    Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / (a*a) + trpnmax * kmtrpn / (b*b))        
    d_Cai = Bcai * (-ICa_tot * cm * Acap / (2*F*vmyo) - Jup  * vnsr / vmyo + Jdiff * vss / vmyo )
            
    '''
    desc: Calcium concentratium in the T-Tubule subspace
    in [mmol/L]
    ''' 
    ICa_ss_tot = ICaL - 2 * INaCa_ss
    a = KmBSR + cass
    b = KmBSL + cass
    Bcass = 1.0 / (1.0 + BSRmax * KmBSR / (a*a) + BSLmax * KmBSL / (b*b))
    d_cass = Bcass * (-ICa_ss_tot * cm * Acap / (2*F*vss) + Jrel * vjsr / vss - Jdiff )
        
    '''
    desc: Calcium concentration in the NSR subspace
    in [mmol/L]
    '''         
    d_Ca_nsr = Jup - Jtr * vjsr / vnsr
        
    '''
    desc: Calcium concentration in the JSR subspace
    in [mmol/L]
    '''        
    a = kmcsqn + Ca_jsr    
    Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / (a*a) )                
    d_Ca_jsr = Bcajsr * (Jtr - Jrel)

    dy[5] = d_Cai 
    dy[6] = d_cass 
    dy[7] = d_Ca_nsr 
    dy[8] = d_Ca_jsr 


    # Membrane potential     
    I_ion = INa + INaL + INaCa + INaK + INab + ICaNa + INaCa_ss + IpCa + ICab + ICaL + Ito + IKr + IKs + IK1 + IKb + ICaK
    d_V_li = -( I_ion + I_stim)

    dy[0] = d_V_li

    # return dy



    

def main():
    start_time = time.time()
    ord2017 = ORD2017(None)    
    print("--- %s seconds ---"%(time.time()-start_time))


if __name__ == '__main__':
    main()


