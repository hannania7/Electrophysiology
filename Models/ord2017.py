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

sys.path.append('../')
sys.path.append('../Lib')
from Protocols import pacing_protocol
import mod_trace as trace
        
        
class Membrane():
    def __init__(self):        
        self.V = -87.0
        self.y0 = [self.V]

    def d_V(self, Iions):
        '''
        in [mV]
        Membrane potential
        '''
        return [-Iions]
    
class Stimulus():
    def __init__(self):
        self.amplitude = -80 # [uA/cm^2]
        self.I = 0
        
    def cal_stimulation(self, pace=0):
        self.I = self.amplitude * pace     
        return self.I

class Phys():
    def __init__(self):
        '''
        The type of cell. Endo=0, Epi=1, Mid=2
        '''
        self.R = 8314           # [J/kmol/K] : Gas constant
        self.T = 310            # [K] : Temperature
        self.F = 96485          # [C/mol] : Faraday's constant
        self.RTF  = self.R*self.T/self.F
        self.FRT  = self.F/(self.R*self.T)
        self.FFRT = self.F*self.F/(self.R*self.T)  
        self.zna = 1
        self.zca = 2
        self.zk = 1

class Cell():
    def __init__(self, phys):
        '''
        Cell geometry
        Page 6        
        '''                        
        self.mode = 0  # The type of cell. Endo=0, Epi=1, Mid=2
        self.L = 0.01  # [cm] Cell length
        self.rad = 0.0011  # [cm] cell radius
        self.vcell = 1000 * 3.14 * self.rad * self.rad * self.L # [uL] Cell volume
        self.Ageo = 2*3.14 * self.rad * self.rad + 2*3.14 * self.rad * self.L   # [cm^2] Geometric cell area
        self.Acap = 2 * self.Ageo            # [cm^2] Capacitative membrane area
        self.vmyo = 0.68 * self.vcell        # [uL] Volume of the cytosolic compartment
        self.vnsr = 0.0552 * self.vcell      # [uL] Volume of the NSR compartment
        self.vjsr = 0.0048 * self.vcell      # [uL] Volume of the JSR compartment
        self.vss = 0.02 * self.vcell         # [uL] Volume of the Submembrane space near the T-tubules
        self.AF = self.Acap / phys.F         # F : Faraday's constant

class Extra():
    '''
  
    '''
    def __init__(self):        
        self.Nao = 140 # [mmol/L] : Extracellular Na+ concentration
        self.Cao = 1.8 # [mmol/L] : Extracellular Ca2+ concentration
        self.Ko  = 5.4 # [mmol/L] : Extracellular K+ concentration




class CaMK(): 
    def __init__(self):
        '''
        CaMKII signalling
        '''
        # initial value
        self.CaMKt = 0.0125840447  # CaMKt=0

        self.y0 = [self.CaMKt]
                
    def d_CaMKt(self, CaMKt, cass):
        '''          
        '''
        KmCaMK = 0.15
        aCaMK  = 0.05
        bCaMK  = 0.00068
        CaMKo  = 0.05
        KmCaM  = 0.0015
        
        CaMKb  = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass)

        d_CaMKt = aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt        

        CaMKa  = CaMKb + CaMKt
        self.f = 1 / (1 + KmCaMK / CaMKa) # Fraction of phosphorylated channels        
        
        return [d_CaMKt]

class Nernst():
    '''
    
    '''
    def __init__(self, phys, cell, extra):               
        self.phys = phys
        self.cell = cell
        self.extra = extra  
            
    def calculate(self, Nai, Ki):
        PKNa = 0.01833          # desc: Permeability ratio K+ to Na+
        self.ENa = self.phys.RTF * np.log(self.extra.Nao / Nai)      # in [mV]  desc: Reversal potential for Sodium currents
        self.EK = self.phys.RTF * np.log(self.extra.Ko / Ki)      # in [mV]  desc: Reversal potential for Potassium currents                
        self.EKs = self.phys.RTF * np.log((self.extra.Ko + PKNa * self.extra.Nao) / (Ki + PKNa * Nai)) # desc: Reversal potential for IKs  in [mV]


class INa():
    '''
    INa :: Fast Sodium current
    Page 6
    
    The fast sodium current is modelled using a Hodgkin-Huxley type formulation
    including activation (m), slow and fast components of inactivation (h) and
    recovery from inactivation (j). The slow component of inactivation and
    recovery from inactivation have an alternative formulation for CaMKII-
    phosphorylated channels.
    '''
    def __init__(self, phys, cell, extra):       
        self.phys = phys
        self.cell = cell
        self.extra = extra   
        
        # initial values    
        self.m = 0.007344121102   # m=0                        10
        self.hf = 0.6981071913     # hf=1                       11
        self.hs = 0.6980895801     # hs=1                       12
        self.j = 0.6979908432     # j=1                        13
        self.hsp = 0.4549485525     # hsp=1                      14
        self.jp = 0.6979245865     # jp=1                       

        self.y0 = [self.m, self.hf, self.hs, self.j, self.hsp, self.jp]

        #: Maximum conductance of INa channels
        self.GNa = 75.0        
            
    def diff_eq(self, V, m, hf, hs, j, hsp, jp, camk, nernst):
        '''
        Activation gate for INa channels
        '''
        mtD1 = 6.765
        mtD2 = 8.552
        mtV1 = 11.64
        mtV2 = 34.77
        mtV3 = 77.42
        mtV4 = 5.955            
        self.tm  = 1.0 / (mtD1 * np.exp((V + mtV1) / mtV2) + mtD2 * np.exp(-(V + mtV3) / mtV4)) # desc: Time constant for m-gate   in [ms]
        mssV1 = 39.57
        mssV2 = 9.871
        mss  = 1.0 / (1.0 + np.exp(-(V + mssV1)/mssV2))  # desc: Steady state value for m-gate         
        d_m = (mss - m) / self.tm           

        hssV1 = 82.9 
        hssV2 = 6.086 
        shift_INa_inact = 0.0
        hss = 1.0 / (1.0 + np.exp((V + hssV1-shift_INa_inact) / hssV2))   # desc: Steady-state value for h-gate
        thf = 1.0 / (1.432e-5 * np.exp(-(V + 1.196 - shift_INa_inact) / 6.285) + 6.1490 * np.exp((V + 0.5096 - shift_INa_inact) / 20.27)) # desc: Time constant for fast development of inactivation in INa   in [ms]
        ths = 1.0 / (0.009794 * np.exp(-(V + 17.95-shift_INa_inact) / 28.05) + 0.3343 * np.exp((V + 5.7300 - shift_INa_inact) / 56.66))  # desc: Time constant for slow development of inactivation in INa  in [ms]
        Ahf = 0.99 # : Fraction of INa channels with fast inactivation
        Ahs = 1.0 - Ahf # : Fraction of INa channels with slow inactivation        
        d_hf = (hss - hf) / thf   # desc: Fast componennt of the inactivation gate for INa channels
        d_hs = (hss - hs) / ths   # desc: Slow componennt of the inactivation gate for non-phosphorylated INa channels
        
        h = Ahf * hf + Ahs * hs   # desc: Inactivation gate for INa
        tj = 2.038 + 1 / (0.02136 * np.exp(-(V + 100.6 - shift_INa_inact) / 8.281) + 0.3052 * np.exp((V + 0.9941 - shift_INa_inact) / 38.45)) # desc: Time constant for j-gate in INa in [ms]
        jss = hss  # desc: Steady-state value for j-gate in INa
        d_j = (jss - j) / tj # desc: Recovery from inactivation gate for non-phosphorylated INa channels

        # Phosphorylated channels
        thsp = 3 * ths # desc: Time constant for h-gate of phosphorylated INa channels  in [ms]
        hssp = 1 / (1 + np.exp((V + 89.1 - shift_INa_inact) / 6.086)) # desc: Steady-state value for h-gate of phosphorylated INa channels        
        d_hsp = (hssp - hsp) / thsp # desc: Slow componennt of the inactivation gate for phosphorylated INa channels

        hp = Ahf * hf + Ahs * hsp # desc: Inactivation gate for phosphorylated INa channels
        tjp = 1.46 * tj # desc: Time constant for the j-gate of phosphorylated INa channels     in [ms]        
        d_jp = (jss - jp) / tjp # desc: Recovery from inactivation gate for phosphorylated INa channels
        
        # Current        
        fINap = camk.f
        INa = self.GNa * (V - nernst.ENa) * m**3 * ((1.0 - fINap) * h * j + fINap * hp * jp) # in [uA/uF]  desc: Fast sodium current
    
        return [d_m, d_hf, d_hs, d_j, d_hsp, d_jp], INa

class INaL():
    '''
    INaL :: Late component of the Sodium current
    Page 7
    '''
    def __init__(self, phys, cell, extra):        
        self.phys = phys
        self.cell = cell
        self.extra = extra  
        
        # initial values            
        self.mL         = 0.0001882617273  # mL=0                       16
        self.hL         = 0.5008548855     # hL=1                       17
        self.hLp        = 0.2693065357     # hLp=1                      

        self.y0 = [self.mL, self.hL, self.hLp]

        #: Maximum conductance          
            
    def diff_eq(self, V, mL, hL, hLp, camk, nernst, ina):
        mLss = 1.0 / (1.0 + np.exp(-(V + 42.85) / 5.264)) # desc: Steady state value of m-gate for INaL
        tmL = ina.tm
        d_mL = (mLss - mL) / tmL # desc: Activation gate for INaL

        thL = 200.0 # [ms] : Time constant for inactivation of non-phosphorylated INaL channels
        hLss = 1.0 / (1.0 + np.exp((V + 87.61) / 7.488))  # desc: Steady-state value for inactivation of non-phosphorylated INaL channels
        
        d_hL = (hLss - hL) / thL # desc: Inactivation gate for non-phosphorylated INaL channels

        hLssp = 1.0 / (1.0 + np.exp((V + 93.81) / 7.488)) # desc: Steady state value for inactivation of phosphorylated INaL channels
        thLp = 3 * thL # in [ms] desc: Time constant for inactivation of phosphorylated INaL channels

        d_hLp = (hLssp - hLp) / thLp  # desc: Inactivation gate for phosphorylated INaL channels
                
        # Current
        GNaL_b = 0.019957499999999975      
        GNaL = GNaL_b
        if self.cell.mode == 1:  GNaL = GNaL_b*0.6  # desc: Adjustment for different cell types            
        fINaLp = camk.f
        INaL = GNaL * (V - nernst.ENa) * mL * ((1 - fINaLp) * hL + fINaLp * hLp)
    
        return [d_mL, d_hL, d_hLp], INaL


class Ito():
    '''
    Ito :: Transient outward Potassium current
    page 8
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
        
        # initial values            
        self.a    = 0.001001097687   # a=0                        19
        self.iF   = 0.9995541745     # iF=1                       20
        self.iS   = 0.5865061736     # iS=1
        self.ap   = 0.0005100862934  # ap=0                       22
        self.iFp  = 0.9995541823     # iFp=1                      23
        self.iSp  = 0.6393399482     # iSp=1  

        self.y0 = [self.a, self.iF, self.iS, self.ap, self.iFp, self.iSp]

        #: Maximum conductance
         
            
    def diff_eq(self, V, a, iF, iS, ap, iFp, iSp, camk, nernst):
        
        ass = 1.0 / (1.0 + np.exp(-(V - 14.34) / 14.82))  # desc: Steady-state value for Ito activation
        one = 1.0 / (1.2089 * (1 + np.exp(-(V - 18.4099) / 29.3814)))
        two = 3.5 / (1 + np.exp((V + 100) / 29.3814))
        ta = 1.0515 / (one + two)  # desc: Time constant for Ito activation  in [ms]        
        d_a = (ass - a) / ta   # desc: Ito activation gate
        
        iss = 1.0 / (1.0 + np.exp((V + 43.94) / 5.711))   # desc: Steady-state value for Ito inactivation
        delta_epi = 1.0
        if self.cell.mode==1:  delta_epi = 1.0 - (0.95 / (1 + np.exp((V + 70.0) / 5.0)))   # desc: Adjustment for different cell types
        tiF_b = (4.562 + 1 / (0.3933 * np.exp(-(V+100) / 100) + 0.08004 * np.exp((V + 50) / 16.59)))  # desc: Time constant for fast component of Ito inactivation     in [ms]
        tiS_b = (23.62 + 1 / (0.001416 * np.exp(-(V + 96.52) / 59.05) + 1.780e-8 * np.exp((V + 114.1) / 8.079)))   # desc: Time constant for slow component of Ito inactivation  in [ms]
        tiF = tiF_b * delta_epi
        tiS = tiS_b * delta_epi
        AiF = 1.0 / (1.0 + np.exp((V - 213.6) / 151.2))  # desc: Fraction of fast inactivating Ito channels
        AiS = 1.0 - AiF        
        d_iF = (iss - iF) / tiF  # desc: Fast component of Ito activation        
        d_iS = (iss - iS) / tiS  # desc: Slow component of Ito activation
        
        i = AiF * iF + AiS * iS     # desc: Inactivation gate for non-phosphorylated Ito
        assp=1.0/(1.0+np.exp(-(V-24.34)/14.82))        
        d_ap = (assp - ap) / ta
            
        dti_develop = 1.354 + 1e-4 / (np.exp((V - 167.4) / 15.89) + np.exp(-(V - 12.23) / 0.2154))
        dti_recover = 1 - 0.5 / (1 + np.exp((V+70) / 20))
        tiFp = dti_develop * dti_recover * tiF # desc: Time constant for fast component of inactivation of phosphorylated Ito channels   in [ms]
        tiSp = dti_develop * dti_recover * tiS  # desc: Time constant for slot component of inactivation of phosphorylated Ito channels  in [ms]        
        d_iFp = (iss - iFp) / tiFp # desc: Fast component of inactivation of phosphorylated Ito channels        
        d_iSp = (iss - iSp) / tiSp # desc: Slow component of inactivation of phosphorylated Ito channels
        
        ip = AiF * iFp + AiS * iSp  # desc: Inactivation gate for phosphorylated Ito channels
        
        # Current
        Gto_b = 0.02
        Gto = Gto_b
        if self.cell.mode==1 :  Gto = Gto_b*4.0
        elif self.cell.mode==2 :  Gto = Gto_b*4.0        
        fItop = camk.f        
        Ito = Gto * (V - nernst.EK) * ((1 - fItop) * a * i + fItop * ap * ip) # desc: Transient outward Potassium current
    
        return [d_a, d_iF, d_iS, d_ap, d_iFp, d_iSp], Ito
    

class ICaL():
    '''
    ICaL  :: L-type Calcium current
    ICaNa :: Sodium current through the L-type Calcium channel
    ICaK  :: Potassium current through the L-type Calcium channel
    Page 9
    The ICaL channel is modeled using activation, inactivation (fast and slow),
    Ca-dependent inactivation (fast and slow) and recovery from Ca-dependent
    inactivation.
    Inactivation and Ca-dependent inactivation have an alternative formulation
    for CaMKII phosphorylated channels.
    '''
    def __init__(self, phys, cell, extra):    
        self.phys = phys
        self.cell = cell
        self.extra = extra      
        
        # initial values            
        self.d          = 2.34e-9          # d=0                        25
        self.ff         = 0.9999999909     # ff=1                       26
        self.fs         = 0.9102412777     # fs=1
        self.fcaf       = 0.9999999909     # fcaf=1                     28
        self.fcas       = 0.9998046777     # fcas=1
        self.jca        = 0.9999738312     # jca=1             30
        self.nca        = 0.002749414044   # nca=0
        self.ffp        = 0.9999999909     # ffp=1             32
        self.fcafp      = 0.9999999909     # fcafp=1           33

        self.y0 = [self.d, self.ff, self.fs, self.fcaf, self.fcas, self.jca, self.nca, self.ffp, self.fcafp]

        #: Maximum conductance
                                
    def diff_eq(self, V, 
                d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp, 
                calcium_cass, sodium_Na_ss, potassium_K_ss,
                camk, nernst):
        
        vfrt = V * self.phys.FRT
        vffrt = V * self.phys.FFRT
        
        # Activation
        dss = 1.0 / (1.0 + np.exp(-(V + 3.94) / 4.23)) # Steady-state value for activation gate of ICaL channel
        td = 0.6 + 1.0 / (np.exp(-0.05 * (V + 6)) + np.exp(0.09 * (V + 14)))  # Time constant for activation gate of ICaL channel   in [ms]
        d_d = (dss - d) / td   # Activation gate of ICaL channel
        
        # Inactivation
        fss = 1.0 / (1.0 + np.exp((V + 19.58) / 3.696)) # Steady-state value for inactivation gate of ICaL channel
        tff = 7.0 + 1.0 / (0.0045 * np.exp(-(V + 20) / 10) + 0.0045 * np.exp((V + 20) / 10))  # Time constant for fast inactivation of ICaL channels  in [ms]
        tfs = 1000 + 1.0 / (0.000035 * np.exp(-(V + 5) / 4) + 0.000035 * np.exp((V + 5) / 6)) # Time constant for fast inactivation of ICaL channels  in [ms]
        Aff = 0.6  # Fraction of ICaL channels with fast inactivation
        Afs = 1.0 - Aff # Fraction of ICaL channels with slow inactivation
        d_ff = (fss - ff) / tff   # Fast inactivation of ICaL channels
        d_fs = (fss - fs) / tfs   # Slow inactivation of ICaL channels        
        f = Aff * ff + Afs * fs  # Inactivation of ICaL channels
        
        # Ca-dependent inactivation
        fcass = fss  # Steady-state value for Ca-dependent inactivation of ICaL channels
        tfcaf = 7.0 + 1.0 / (0.04 * np.exp(-(V - 4.0) / 7.0) + 0.04 * np.exp((V - 4.0) / 7.0))  # Time constant for fast Ca-dependent inactivation of ICaL channels in [ms]
        tfcas = 100.0 + 1 / (0.00012 * np.exp(-V / 3) + 0.00012 * np.exp(V / 7))  # Time constant for slow Ca-dependent inactivation of ICaL channels  in [ms]
        Afcaf = 0.3 + 0.6 / (1 + np.exp((V - 10) / 10))  # Fraction of ICaL channels with fast Ca-dependent inactivation
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
        anca = 1.0 / (k2n / km2n + (1 + Kmn / calcium_cass)**4.0)  # Fraction of channels in Ca-depdent inactivation mode
        d_nca = anca * k2n - nca*km2n   # Fraction of channels in Ca-depdent inactivation mode
        
        # Total currents through ICaL channel
        v0 = 0
        B_1 = 2.0 * self.phys.FRT
        A_1 = ( 4.0 * self.phys.FFRT * (calcium_cass * np.exp(2 * vfrt) - 0.341 * self.extra.Cao) ) / B_1
        U_1 = B_1 * (V-v0)
        PhiCaL = (A_1*U_1)/(np.exp(U_1)-1.0)
        if -1e-7<=U_1 and U_1<=1e-7 :  PhiCaL = A_1 * (1.0-0.5*U_1)        
        B_2 = self.phys.FRT
        A_2 = ( 0.75 * self.phys.FFRT * (sodium_Na_ss * np.exp(vfrt) - 0.341 * self.extra.Nao) ) / B_2
        U_2 = B_2 * (V-v0)
        PhiCaNa = (A_2*U_2)/(np.exp(U_2)-1.0)
        if (-1e-7<=U_2 and U_2<=1e-7) :  PhiCaNa = A_2 * (1.0-0.5*U_2) 
        B_3 = self.phys.FRT
        A_3 = ( 0.75 * self.phys.FFRT * (potassium_K_ss * np.exp(vfrt) - 0.341 * self.extra.Ko) ) / B_3
        U_3 = B_3 * (V-v0)
        PhiCaK = (A_3*U_3)/(np.exp(U_3)-1.0)
        if (-1e-7<=U_3 and U_3<=1e-7) :  PhiCaK = A_3 * (1.0-0.5*U_3)
        
        PCa_b = 0.0001007
        PCa = PCa_b
        if self.cell.mode==1: PCa = PCa_b*1.2
        elif self.cell.mode==2: PCa = PCa_b*2.5
                    
        PCap   = 1.1      * PCa
        PCaNa  = 0.00125  * PCa
        PCaK   = 3.574e-4 * PCa
        PCaNap = 0.00125  * PCap
        PCaKp  = 3.574e-4 * PCap
        flCaLp = camk.f
        g  = d * (f  * (1.0 - nca) + jca * fca  * nca)   # Conductivity of non-phosphorylated ICaL channels
        gp = d * (fp * (1.0 - nca) + jca * fcap * nca)   # Conductivity of phosphorylated ICaL channels
        
        ICaL   = (1.0 - flCaLp) * PCa   * PhiCaL  * g + flCaLp * PCap   * PhiCaL  * gp  # L-type Calcium current   in [uA/uF]
        ICaNa  = (1.0 - flCaLp) * PCaNa * PhiCaNa * g + flCaLp * PCaNap * PhiCaNa * gp   # Sodium current through ICaL channels  in [uA/uF]
        ICaK   = (1.0 - flCaLp) * PCaK  * PhiCaK  * g + flCaLp * PCaKp  * PhiCaK  * gp   # Potassium current through ICaL channels  in [uA/uF]
    
        return [d_d, d_ff, d_fs, d_fcaf, d_fcas, d_jca, d_nca, d_ffp, d_fcafp,], ICaL, ICaNa, ICaK


class IKr():
    '''
    IKr :: Rapid delayed rectifier Potassium current
    Page 11
    Modelled with activation (fast and slow) and an instantaneous inactivation.
    '''
    def __init__(self, phys, cell, extra):   
        self.phys = phys
        self.cell = cell
        self.extra = extra       
        
        # initial values            
        self.IC1         = 0.999637         #
        self.IC2         = 6.83208e-05      #
        self.C1          = 1.80145e-08      #
        self.C2          = 8.26619e-05      #
        self.O           = 0.00015551       #
        self.IO          = 5.67623e-05      #
        self.IObound     = 0                #
        self.Obound      = 0                #
        self.Cbound      = 0                #
        self.D           = 0                #

        self.y0 = [self.IC1, self.IC2, self.C1, self.C2, self.O, self.IO, self.IObound, self.Obound, self.Cbound, self.D]

        #: Maximum conductance
         
            
    def diff_eq(self, V, IC1, IC2, C1, C2, O, IO, IObound, Obound, Cbound, D, camk, nernst):
        
        A1 = 0.0264
        B1 = 4.631E-05
        q1 = 4.843
        A2 = 4.986E-06 
        B2 = -0.004226
        q2 = 4.23  
        A3 = 0.001214
        B3 = 0.008516
        q3 = 4.962
        A4 = 1.854E-05
        B4 = -0.04641 
        q4 = 3.769
        A11 = 0.0007868
        B11 = 1.535E-08
        q11 = 4.942
        A21 = 5.455E-06
        B21 = -0.1688
        q21 = 4.156
        A31 = 0.005509
        B31 = 7.771E-09 
        q31 = 4.22
        A41 = 0.001416
        B41 = -0.02877
        q41 = 1.459
        A51 = 0.4492
        B51 = 0.008595
        q51 = 5
        A52 = 0.3181
        B52 = 3.613E-08
        q52 = 4.663
        A53 = 0.149
        B53 = 0.004668
        q53 = 2.412
        A61 = 0.01241
        B61 = 0.1725
        q61 = 5.568
        A62 = 0.3226
        B62 = -0.0006575
        q62 = 5
        A63 = 0.008978
        B63 = -0.02215
        q63 = 5.682 
        Kmax = 0
        Ku = 0
        n = 1
        halfmax = 1
        Kt = 0
        Vhalf = 1
        Temp = 37
        d_IC1 = (-(A11 * np.exp(B11*V) * IC1 * np.exp((Temp-20.0)*np.log(q11)/10.0) - (A21 * np.exp(B21*V) * IC2 * np.exp((Temp-20.0)*np.log(q21)/10.0)))) \
                + A51 * np.exp(B51*V) * C1 * np.exp((Temp-20.0)*np.log(q51)/10.0) - (A61 * np.exp(B61*V) * IC1 * np.exp((Temp-20.0)*np.log(q61)/10.0))
        d_IC2 = (A11 * np.exp(B11*V) * IC1 * np.exp((Temp-20.0)*np.log(q11)/10.0) - (A21 * np.exp(B21*V) * IC2 * np.exp((Temp-20.0)*np.log(q21)/10.0)) \
                - (A3 * np.exp(B3*V) * IC2 * np.exp((Temp-20.0)*np.log(q3)/10.0) - (A4 * np.exp(B4*V) * IO * np.exp((Temp-20.0)*np.log(q4)/10.0))))  \
                    + A52 * np.exp(B52*V) * C2 * np.exp((Temp-20.0)*np.log(q52)/10.0) - (A62 * np.exp(B62*V) * IC2 * np.exp((Temp-20.0)*np.log(q62)/10.0))
        d_C1 = -(A1 * np.exp(B1*V) * C1 * np.exp((Temp-20.0)*np.log(q1)/10.0) - (A2 * np.exp(B2*V) * C2 * np.exp((Temp-20.0)*np.log(q2)/10.0))) \
                -(A51 * np.exp(B51*V) * C1 * np.exp((Temp-20.0)*np.log(q51)/10.0) - (A61 * np.exp(B61*V) * IC1 * np.exp((Temp-20.0)*np.log(q61)/10.0)))
        d_C2 = A1 * np.exp(B1*V) * C1 * np.exp((Temp-20.0)*np.log(q1)/10.0) - (A2 * np.exp(B2*V) * C2 * np.exp((Temp-20.0)*np.log(q2)/10.0))  \
                - (A31 * np.exp(B31*V) * C2 * np.exp((Temp-20.0)*np.log(q31)/10.0) - (A41 * np.exp(B41*V) * O   * np.exp((Temp-20.0)*np.log(q41)/10.0))) \
                    - (A52 * np.exp(B52*V) * C2 * np.exp((Temp-20.0)*np.log(q52)/10.0) - (A62 * np.exp(B62*V) * IC2 * np.exp((Temp-20.0)*np.log(q62)/10.0)))
        d_O = A31 * np.exp(B31*V) * C2 * np.exp((Temp-20.0)*np.log(q31)/10.0) - (A41 * np.exp(B41*V) * O * np.exp((Temp-20.0)*np.log(q41)/10.0)) \
                - (A53 * np.exp(B53*V) * O * np.exp((Temp-20.0)*np.log(q53)/10.0) - (A63 * np.exp(B63*V) * IO * np.exp((Temp-20.0)*np.log(q63)/10.0))) \
                    - ((Kmax*Ku*np.exp(n*np.log(D)))/(np.exp(n*np.log(D))+halfmax) * O  - Ku * Obound)
        d_IO = (A3 * np.exp(B3*V) * IC2 * np.exp((Temp-20.0)*np.log(q3)/10.0) - (A4 * np.exp(B4*V) * IO * np.exp((Temp-20.0)*np.log(q4)/10.0))) \
                + A53 * np.exp(B53*V) * O * np.exp((Temp-20.0)*np.log(q53)/10.0) - (A63 * np.exp(B63*V) * IO * np.exp((Temp-20.0)*np.log(q63)/10.0)) \
                    - ((Kmax*Ku*np.exp(n*np.log(D)))/(np.exp(n*np.log(D))+halfmax) * IO - (Ku*A53*np.exp(B53*V)*np.exp((Temp-20.0)*np.log(q53)/10.0))/(A63*np.exp(B63*V)*np.exp((Temp-20)*np.log(q63)/10)) * IObound)
        d_IObound = ( (Kmax * Ku * np.exp(n*np.log(D))) / (np.exp(n*np.log(D)) + halfmax) * IO \
                    - ( (Ku*A53*np.exp(B53*V)*np.exp((Temp-20.0)*np.log(q53)/10.0))/(A63*np.exp(B63*V)*np.exp((Temp-20.0)*np.log(q63)/10.0)) * IObound ) ) \
                        + Kt / (1 + np.exp(-(V-Vhalf)/6.789)) * Cbound - Kt * IObound
        d_Obound = ( (Kmax * Ku * np.exp(n*np.log(D))) / (np.exp(n*np.log(D)) + halfmax) *  O - Ku * Obound ) + Kt / (1 + np.exp(-(V-Vhalf)/6.789)) * Cbound - Kt * Obound
        d_Cbound = -( Kt / (1 + np.exp(-(V-Vhalf)/6.789)) * Cbound - Kt * Obound ) - ( Kt / (1 + np.exp(-(V-Vhalf)/6.789)) * Cbound - Kt * IObound )
        d_D = 0
        
        # Current
        GKr_b = 0.04658545454545456  
        GKr = GKr_b   # desc: Maximum conductance of IKr channels
        if self.cell.mode ==1:  GKr = GKr_b*1.3
        elif self.cell.mode==2:  GKr = GKr_b*0.8
            
        IKr = GKr * np.sqrt(self.extra.Ko / 5.4) * O * (V - nernst.EK) # desc: Rapid delayed Potassium current  in [uA/uF]
    
        return [d_IC1, d_IC2, d_C1, d_C2, d_O, d_IO, d_IObound, d_Obound, d_Cbound, d_D], IKr
    
class IKs():
    '''
    IKs :: Slow delayed rectifier Potassium current
    Page 11
    Modelled with two activation channels
    '''
    def __init__(self, phys, cell, extra):        
        self.phys = phys
        self.cell = cell
        self.extra = extra  
        
        # initial values            
        self.xs1 = 0.2707758025     # xs1=0                      36
        self.xs2 = 0.0001928503426  # xs2=0

        self.y0 = [self.xs1, self.xs2]

        #: Maximum conductance
                    
    def diff_eq(self, V, xs1, xs2, Cai, camk, nernst):
        
        xs1ss  = 1.0 / (1.0 + np.exp(-(V + 11.60) / 8.932)) # desc: Steady-state value for activation of IKs channels
        txs1_max = 817.3
        txs1 = txs1_max + 1.0 / (2.326e-4 * np.exp((V + 48.28) / 17.80) + 0.001292 * np.exp(-(V + 210) / 230)) # desc: Time constant for slow, low voltage IKs activation
                        
        d_xs1 = (xs1ss - xs1) / txs1  # desc: Slow, low voltage IKs activation
        
        xs2ss = xs1ss
        txs2 = 1.0 / (0.01 * np.exp((V - 50) / 20) + 0.0193 * np.exp(-(V + 66.54) / 31.0)) # desc: Time constant for fast, high voltage IKs activation
        
        d_xs2 = (xs2ss - xs2) / txs2   # desc: Fast, high voltage IKs activation
        
        KsCa = 1.0 + 0.6 / (1.0 + (3.8e-5 / Cai)**1.4) # desc: Maximum conductance for IKs
        GKs_b = 0.006358000000000001
        GKs = GKs_b
        if self.cell.mode==1 : GKs = GKs_b * 1.4  # desc: Conductivity adjustment for cell type        
        IKs = GKs * KsCa * xs1 * xs2 * (V - nernst.EKs)  # Slow delayed rectifier Potassium current
            
        return [d_xs1, d_xs2], IKs
    
class IK1():
    '''
    IK1 :: Inward rectifier Potassium current
    Page 12
    Modelled with an activation channel and an instantaneous inactivation channel
    '''
    def __init__(self, phys, cell, extra):    
        self.phys = phys
        self.cell = cell
        self.extra = extra      
        
        # initial values            
        self.xk1 = 0.9967597594    
        
        self.y0 = [self.xk1]

        #: Maximum conductance
                    
    def diff_eq(self, V, xk1, camk, nernst):
        
        xk1ss = 1 / (1 + np.exp(-(V + 2.5538 * self.extra.Ko + 144.59) / (1.5692 * self.extra.Ko + 3.8115))) # Steady-state value for activation of IK1 channels
        txk1 = 122.2 / (np.exp(-(V + 127.2) / 20.36) + np.exp((V + 236.8) / 69.33))  # Time constant for activation of IK1 channels
        
        d_xk1 = (xk1ss - xk1) / txk1  # Activation of IK1 channels
        
        rk1 = 1.0 / (1.0 + np.exp((V + 105.8 - 2.6 * self.extra.Ko) / 9.493))   # Inactivation of IK1 channels        
        
        # desc: Conductivity of IK1 channels, cell-type dependent
        GK1_b = 0.3239783999999998
        GK1 = GK1_b  
        if self.cell.mode==1 :  GK1 = GK1_b * 1.2
        elif self.cell.mode==2 :  GK1 = GK1_b * 1.3        
        IK1 = GK1 * np.sqrt(self.extra.Ko) * rk1 * xk1 * (V - nernst.EK)  # Inward rectifier Potassium current
            
        return [d_xk1], IK1
    
    
class INaCa():
    '''
    INaCa :: Sodium/Calcium exchange current
    page 12
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    

        #: Maximum conductance
                     
    def calculate(self, V, Nai, Cai):
        
        self.kna1   = 15.0        
        self.kna2   = 5.0        
        self.kna3   = 88.12        
        self.kasymm = 12.5        
        self.wna    = 6.0e4
        self.wca    = 6.0e4
        self.wnaca  = 5.0e3 
        self.kcaon  = 1.5e6        
        self.kcaoff = 5.0e3        
        qna = 0.5224
        qca = 0.1670
        self.hca    = np.exp(qca * V * self.phys.FRT)
        self.hna    = np.exp(qna * V * self.phys.FRT)
        
        # Parameters h
        h1  = 1.0 + Nai / self.kna3 * (1 + self.hna)
        h2  = (Nai * self.hna) / (self.kna3 * h1)
        h3  = 1.0 / h1
        h4  = 1.0 + Nai / self.kna1 * (1 + Nai / self.kna2)
        h5  = Nai * Nai / (h4 * self.kna1 * self.kna2)
        h6  = 1.0 / h4
        h7  = 1.0 + self.extra.Nao / self.kna3 * (1 + 1 / self.hna)
        h8  = self.extra.Nao / (self.kna3 * self.hna * h7)
        h9  = 1.0 / h7
        h10 = self.kasymm + 1 + self.extra.Nao / self.kna1 * (1 + self.extra.Nao / self.kna2)
        h11 = self.extra.Nao * self.extra.Nao / (h10 * self.kna1 * self.kna2)
        h12 = 1.0 / h10
        
        # Parameters k
        k1   = h12 * self.extra.Cao * self.kcaon
        k2   = self.kcaoff
        k3p  = h9 * self.wca
        k3pp = h8 * self.wnaca
        k3   = k3p + k3pp
        k4p  = h3 * self.wca / self.hca
        k4pp = h2 * self.wnaca
        k4   = k4p + k4pp
        k5   = self.kcaoff
        k6   = h6 * Cai * self.kcaon
        k7   = h5 * h2 * self.wna
        k8   = h8 * h11 * self.wna
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
        self.Gncx = Gncx_b
        if self.cell.mode==1: self.Gncx = Gncx_b*1.1
        elif self.cell.mode==2: self.Gncx = Gncx_b*1.4
                    
        return 0.8 * self.Gncx * allo * (self.phys.zna * JncxNa + self.phys.zca * JncxCa)  # Sodium/Calcium exchange current  in [uA/uF]
    
class INaCa_ss():
    '''
    INaCa_ss :: Sodium/Calcium exchanger current into the L-type subspace
    Page 12
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, Na_ss, cass, inaca):
        
        h1  = 1.0 + Na_ss / inaca.kna3 * (1 + inaca.hna)
        h2  = (Na_ss * inaca.hna)/(inaca.kna3 * h1)
        h3  = 1 / h1
        h4  = 1 + Na_ss / inaca.kna1 * (1 + Na_ss / inaca.kna2)
        h5  = Na_ss * Na_ss /(h4 * inaca.kna1 * inaca.kna2)
        h6  = 1 / h4
        h7  = 1 + self.extra.Nao / inaca.kna3 * (1 + 1 / inaca.hna)
        h8  = self.extra.Nao / (inaca.kna3 * inaca.hna * h7)
        h9  = 1/h7
        h10 = inaca.kasymm + 1 + self.extra.Nao / inaca.kna1 * (1 + self.extra.Nao / inaca.kna2)
        h11 = self.extra.Nao * self.extra.Nao / (h10 * inaca.kna1 * inaca.kna2)
        h12 = 1/h10
        k1   = h12 * self.extra.Cao * inaca.kcaon
        k2   = inaca.kcaoff
        k3p  = h9 * inaca.wca
        k3pp = h8 * inaca.wnaca
        k3   = k3p + k3pp
        k4p  = h3 * inaca.wca / inaca.hca
        k4pp = h2 * inaca.wnaca
        k4   = k4p + k4pp
        k5   = inaca.kcaoff
        k6   = h6 * cass * inaca.kcaon
        k7   = h5 * h2 * inaca.wna
        k8   = h8 * h11 * inaca.wna
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
        
        return 0.2 * inaca.Gncx * allo * (self.phys.zna * JncxNa + self.phys.zca * JncxCa)  # Sodium/Calcium exchange current into the T-Tubule subspace  in [uA/uF]
    
class INaK():
    '''
    INaK :: Sodium/Potassium ATPase current
    Page 14
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, V, Nai, Na_ss, Ki, K_ss):
        
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
        Knai = Knai0 * np.exp(delta * V * self.phys.FRT / 3.0)
        Knao = Knao0 * np.exp((1.0-delta) * V * self.phys.FRT / 3.0)
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
        b2 = (k2m * (self.extra.Nao / Knao)**3) / ((1 + self.extra.Nao / Knao)**3 + (1 + self.extra.Ko / Kko)**2 - 1)
        a3 = (k3p * (self.extra.Ko / Kko)**2  ) / ((1 + self.extra.Nao / Knao)**3 + (1 + self.extra.Ko / Kko)**2 - 1)
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
        if self.cell.mode==1: Pnak = Pnak_b*0.9
        elif self.cell.mode==2: Pnak = Pnak_b*0.7
                    
        return Pnak * (self.phys.zna * JnakNa + self.phys.zk * JnakK)  # Sodium/Potassium ATPase current    in [uA/uF]
   
class IKb():
    '''
    IKb :: Background Potassium current
    Page 15
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, V, nernst):
        
        xkb = 1.0 / (1.0 + np.exp(-(V - 14.48) / 18.34))
        GKb_b = 0.003
        GKb = GKb_b
        if self.cell.mode==1: GKb = GKb_b*0.6
                    
        return GKb * xkb * (V - nernst.EK)  # Background Potassium current   in [uA/uF]
    
class INab():
    '''
    INab :: Background Sodium current
    Page 15
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, V, Nai):
        
        B = self.phys.FRT
        v0 = 0
        PNab = 3.75e-10
        A = PNab * self.phys.FFRT * (Nai * np.exp(V * self.phys.FRT) - self.extra.Nao) / B            
        U = B * (V - v0)
        INab = (A*U)/(np.exp(U)-1.0)
        if -1e-7<=U and U<=1e-7: INab = A*(1.0-0.5*U)   # Background Sodium current   in [uA/uF]
                    
        return INab
    
class ICab():
    '''
    ICab :: Background Calcium current
    Page 15
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, V, Cai):
        
        B = 2 * self.phys.FRT
        v0 = 0
        PCab = 2.5e-8
        A = PCab * 4.0 * self.phys.FFRT * (Cai * np.exp( 2.0 * V * self.phys.FRT) - 0.341 * self.extra.Cao) / B
        U = B * (V - v0)
        ICab = (A*U)/(np.exp(U)-1.0)
        if -1e-7<=U and U<=1e-7:  ICab = A*(1.0-0.5*U) # Background Calcium current  in [uA/uF]
                    
        return ICab
    
class IpCa():
    '''
    IpCa :: Sarcolemmal Calcium pump current
    Page 15
    '''
    def __init__(self, phys, cell, extra):      
        self.phys = phys
        self.cell = cell
        self.extra = extra    
            
    def calculate(self, Cai):
        
        GpCa = 0.0005
        KmCap = 0.0005
        IpCa = GpCa * Cai / (KmCap + Cai)  # Sarcolemmal Calcium pump current   in [uA/uF]
                    
        return IpCa
    
class Ryr():
    '''
    '''
    def __init__(self, phys, cell, extra):
        self.phys = phys
        self.cell = cell
        self.extra = extra  

        # initial values        
        self.Jrelnp      = 2.5e-7           # Jrelnp=0
        self.Jrelp       = 3.12e-7          # Jrelp=0                    

        self.y0 = [self.Jrelnp, self.Jrelp]
        
    def diff_eq(self, V, Jrelnp, Jrelp, Ca_jsr, ICaL, camk):

        bt=4.75
        a_rel=0.5*bt
        Jrel_inf_temp = a_rel * -ICaL / (1 + (1.5 / Ca_jsr)**8)
        Jrel_inf = Jrel_inf_temp
        if (self.cell.mode == 2):   Jrel_inf = Jrel_inf_temp * 1.7            
        tau_rel_temp = bt / (1.0 + 0.0123 / Ca_jsr)
        tau_rel = tau_rel_temp
        if (tau_rel_temp < 0.001):   tau_rel = 0.001                        
        d_Jrelnp = (Jrel_inf - Jrelnp) / tau_rel                       

        btp = 1.25*bt
        a_relp = 0.5*btp
        Jrel_temp = a_relp * -ICaL / (1 + (1.5 / Ca_jsr)**8)
        Jrel_infp = Jrel_temp
        if self.cell.mode == 2:   Jrel_infp = Jrel_temp * 1.7
        tau_relp_temp = btp / (1 + 0.0123 / Ca_jsr)
        tau_relp = tau_relp_temp
        if tau_relp_temp < 0.001:   tau_relp = 0.001        
        d_Jrelp = (Jrel_infp - Jrelp) / tau_relp            
            
        fJrelp = camk.f
        Jrel_scaling_factor = 1.0
        self.Jrel = Jrel_scaling_factor * (1 - fJrelp) * Jrelnp + fJrelp * Jrelp # desc: SR Calcium release flux via Ryanodine receptor  in [mmol/L/ms]

        return [d_Jrelnp, d_Jrelp]

class Diff():
    '''
    Diffusion fluxes
    Page 16
    '''
    def __init__(self):
        self.JdiffNa = None
        self.JdiffK  = None
        self.Jdiff   = None

    def calculate(self, Nai, Na_ss, Ki, K_ss, Cai, cass):
        self.JdiffNa = (Na_ss - Nai) / 2   # (sodium.Na_ss - sodium.Nai) / 2
        self.JdiffK  = (K_ss  - Ki)  / 2  # (potassium.K_ss  - potassium.Ki)  / 2
        self.Jdiff   = (cass - Cai) / 0.2     # (calcium.cass - calcium.Cai) / 0.2

class Serca():
    '''
    Jup :: Calcium uptake via SERCA pump
    Page 17
    '''
    def __init__(self, phys, cell, extra):
        self.phys = phys
        self.cell = cell
        self.extra = extra  
        
    def calculate(self, Cai, Ca_nsr, Ca_jsr, camk):
        upScale = 1.0
        if (self.cell.mode == 1):   upScale = 1.3
        Jupnp = upScale * (0.004375 * Cai / (Cai + 0.00092))
        Jupp  = upScale * (2.75 * 0.004375 * Cai / (Cai + 0.00092 - 0.00017))
        fJupp = camk.f
        Jleak = 0.0039375 * Ca_nsr / 15.0 # in [mmol/L/ms]
        Jup_b = 1.0
        self.Jup = Jup_b * ((1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak) # desc: Total Ca2+ uptake, via SERCA pump, from myoplasm to nsr in [mmol/L/ms]

class Trans_flux():
    '''
    '''
    def __init__(self):
        self.Jtr = None
        
    def calculate(self, calcium_Ca_nsr, calcium_Ca_jsr):
        self.Jtr = (calcium_Ca_nsr - calcium_Ca_jsr) / 100.0   # desc: Ca2+ translocation from nsr to jsr    in [mmol/L/ms]
        

class Sodium():
    '''
    Intracellular Sodium concentrations
    Page 18
    '''
    def __init__(self, phys, cell, extra): 
        self.phys = phys
        self.cell = cell
        self.extra = extra       
        
        # initial values
        self.Nai      = 7.268004498      # nai=7.268004498            2
        self.Na_ss    = 7.268089977      # nass=nai                   3

        self.y0 = [self.Nai, self.Na_ss]

    def diff_eq(self, Nai : float, Na_ss : float, 
                INa, INaL, INaCa, INaK, INab, ical_ICaNa, INaCa_ss, 
                diff): 
        cm = 1.0
        INa_tot = INa + INaL + 3*INaCa + 3*INaK + INab 
        d_Nai = (-INa_tot * self.cell.Acap * cm) / (self.phys.F * self.cell.vmyo) + diff.JdiffNa * self.cell.vss / self.cell.vmyo   
    
        INa_ss_tot = ical_ICaNa + 3*INaCa_ss
        d_Na_ss = -INa_ss_tot * cm * self.cell.Acap / (self.phys.F * self.cell.vss) - diff.JdiffNa

        return [d_Nai, d_Na_ss]
        
class Potassium():
    '''
    Intracellular Potassium concentrations
    Page 18
    '''
    def __init__(self, phys, cell, extra):        
        self.phys = phys
        self.cell = cell
        self.extra = extra  
                             
        # initial values
        self.Ki    = 144.6555918      # ki=144.6555918             4
        self.K_ss  = 144.6555651      # kss=ki                     5

        self.y0 = [self.Ki, self.K_ss]

    def diff_eq(self, Ki, K_ss, 
                Ito, IKr, IKs, IK1, IKb, INaK, i_stim, ical_ICaK,
                diff):
        cm = 1.0
        
        IK_tot = Ito + IKr + IKs + IK1 + IKb - 2 * INaK 
        d_Ki = -(IK_tot + i_stim) * cm * self.cell.Acap / (self.phys.F * self.cell.vmyo) + diff.JdiffK * self.cell.vss / self.cell.vmyo

        IK_ss_tot = ical_ICaK
        d_K_ss = -IK_ss_tot * cm * self.cell.Acap / (self.phys.F * self.cell.vss) - diff.JdiffK

        return [d_Ki, d_K_ss]
    
class Calcium():
    '''
    Intracellular Calcium concentrations and buffers
    Page 18
    '''
    def __init__(self, phys, cell, extra):        
        self.phys = phys
        self.cell = cell
        self.extra = extra  
        
        # initial values
        self.Cai = 8.6e-05 
        self.cass = 8.49e-05  
        self.Ca_nsr = 1.619574538
        self.Ca_jsr = 1.571234014 
        
        self.y0 = [self.Cai, self.cass, self.Ca_nsr, self.Ca_jsr]

    def diff_eq(self, Cai, cass, Ca_nsr, Ca_jsr, 
                IpCa, ICab, INaCa, ICaL, INaCa_ss, 
                ryr, serca, diff, trans_flux):
        cm = 1.0
        cmdnmax_b = 0.05
        cmdnmax = cmdnmax_b  
        if self.cell.mode == 1:
            cmdnmax = 1.3*cmdnmax_b        
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
        d_Cai = Bcai * (-ICa_tot * cm * self.cell.Acap / (2*self.phys.F*self.cell.vmyo) - serca.Jup  * self.cell.vnsr / self.cell.vmyo + diff.Jdiff * self.cell.vss / self.cell.vmyo )
                
        '''
        desc: Calcium concentratium in the T-Tubule subspace
        in [mmol/L]
        ''' 
        ICa_ss_tot = ICaL - 2 * INaCa_ss
        a = KmBSR + cass
        b = KmBSL + cass
        Bcass = 1 / (1 + BSRmax * KmBSR / (a*a) + BSLmax * KmBSL / (b*b))
        d_cass = Bcass * (-ICa_ss_tot * cm * self.cell.Acap / (2*self.phys.F*self.cell.vss) + ryr.Jrel * self.cell.vjsr / self.cell.vss - diff.Jdiff )
            
        '''
        desc: Calcium concentration in the NSR subspace
        in [mmol/L]
        '''         
        d_Ca_nsr = serca.Jup - trans_flux.Jtr * self.cell.vjsr / self.cell.vnsr
            
        '''
        desc: Calcium concentration in the JSR subspace
        in [mmol/L]
        '''         
        Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / (kmcsqn + Ca_jsr)**2)                
        d_Ca_jsr = Bcajsr * (trans_flux.Jtr - ryr.Jrel)

        return [d_Cai, d_cass, d_Ca_nsr, d_Ca_jsr]
        
class ORD2017():
    """    
    O'Hara-Rudy CiPA v1.0 (2017)
    """
    def __init__(self, protocol=None):
        
        self.name = "ORD2017"
        
        self.phys = Phys()
        self.cell = Cell(self.phys)        
        self.extra = Extra()
        
        self.current_response_info = trace.CurrentResponseInfo(protocol)
        
        self.protocol = protocol
        
        self.membrane = Membrane()
        self.stimulus = Stimulus()      

        self.nernst = Nernst(self.phys, self.cell, self.extra)
        self.camk = CaMK()  
        
        self.ina = INa(self.phys, self.cell, self.extra)
        self.inal = INaL(self.phys, self.cell, self.extra)
        self.ito = Ito(self.phys, self.cell, self.extra)        
        self.ical = ICaL(self.phys, self.cell, self.extra)
        self.ikr = IKr(self.phys, self.cell, self.extra)
        self.iks = IKs(self.phys, self.cell, self.extra)
        self.ik1 = IK1(self.phys, self.cell, self.extra)
        
        self.inaca = INaCa(self.phys, self.cell, self.extra)
        self.inacass = INaCa_ss(self.phys, self.cell, self.extra)
        self.inak = INaK(self.phys, self.cell, self.extra)
        self.ikb = IKb(self.phys, self.cell, self.extra)
        self.inab = INab(self.phys, self.cell, self.extra)
        self.icab = ICab(self.phys, self.cell, self.extra)
        self.ipca = IpCa(self.phys, self.cell, self.extra)
                
        self.ryr = Ryr(self.phys, self.cell, self.extra)
        self.serca = Serca(self.phys, self.cell, self.extra)
        self.diff = Diff()
        self.trans_flux = Trans_flux()
        
        self.sodium = Sodium(self.phys, self.cell, self.extra)
        self.potassium = Potassium(self.phys, self.cell, self.extra)      
        self.calcium = Calcium(self.phys, self.cell, self.extra)

        self.y0 = self.membrane.y0 + self.camk.y0 + \
                    self.ina.y0 + self.inal.y0 + self.ito.y0 + self.ical.y0 + self.ikr.y0 + self.iks.y0 + self.ik1.y0 +\
                        self.ryr.y0 + \
                            self.sodium.y0 + self.potassium.y0 + self.calcium.y0
        self.params = []

    def set_result(self, t, y, log=None):
        self.times =  t
        self.V = y[0]    
                         
    def differential_eq(self, t, y):    
        V, CaMKt, \
            m, hf, hs, j, hsp, jp, \
                mL, hL, hLp, \
                    a, iF, iS, ap, iFp, iSp, \
                        d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp,\
                            IC1, IC2, C1, C2, O, IO, IObound, Obound, Cbound, D,\
                                xs1, xs2, \
                                    xk1, \
                                        Jrelnp, Jrelp, \
                                            Nai, Na_ss, Ki, K_ss, Cai, cass, Ca_nsr, Ca_jsr = y

        # Calculate Nernst  
        self.nernst.calculate(Nai, Ki)        
        
        # CaMKt
        d_CaMKt_li = self.camk.d_CaMKt(CaMKt, cass)
        
        # currents  
        d_INa_li, INa = self.ina.diff_eq(V, m, hf, hs, j, hsp, jp, self.camk, self.nernst)
        d_INaL_li, INaL = self.inal.diff_eq(V, mL, hL, hLp, self.camk, self.nernst, self.ina)
        d_Ito_li, Ito = self.ito.diff_eq(V, a, iF, iS, ap, iFp, iSp, self.camk, self.nernst)
        d_ICaL_li, ICaL, ICaNa, ICaK = self.ical.diff_eq(V, d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp, cass, Na_ss, K_ss, self.camk, self.nernst)
        d_IKr_li, IKr = self.ikr.diff_eq(V, IC1, IC2, C1, C2, O, IO, IObound, Obound, Cbound, D, self.camk, self.nernst)
        d_IKs_li, IKs = self.iks.diff_eq(V, xs1, xs2, Cai, self.camk, self.nernst) 
        d_IK1_li, IK1 = self.ik1.diff_eq(V, xk1, self.camk, self.nernst) 
        
        INaCa = self.inaca.calculate(V, Nai, Cai)
        INaCa_ss = self.inacass.calculate(Na_ss, cass, self.inaca)
        INaK = self.inak.calculate(V, Nai, Na_ss, Ki, K_ss)
        IKb = self.ikb.calculate(V, self.nernst)
        INab = self.inab.calculate(V, Nai)
        ICab = self.icab.calculate(V, Cai)
        IpCa = self.ipca.calculate(Cai)


        d_ryr_li = self.ryr.diff_eq(V, Jrelnp, Jrelp, Ca_jsr, ICaL, self.camk)
        self.diff.calculate(Nai, Na_ss, Ki, K_ss, Cai, cass)        
        self.serca.calculate(Cai, Ca_nsr, Ca_jsr, self.camk)
        self.trans_flux.calculate(Ca_nsr, Ca_jsr)

        d_sodium_li = self.sodium.diff_eq( Nai, Na_ss, 
                                           INa, INaL, INaCa, INaK, INab, ICaNa, INaCa_ss,
                                           self.diff)

        d_potassium_li = self.potassium.diff_eq(Ki, K_ss, 
                                                Ito, IKr, IKs, IK1, IKb, INaK, self.stimulus.I, ICaK,
                                                self.diff)

        d_calcium_li = self.calcium.diff_eq(Cai, cass, Ca_nsr, Ca_jsr, 
                                            IpCa, ICab, INaCa, ICaL, INaCa_ss, 
                                            self.ryr, self.serca, self.diff, self.trans_flux)   
                   
        # Membrane potential     
        I_ion = INa + INaL + INaCa + INaK + INab + ICaNa + INaCa_ss + IpCa + ICab + ICaL + Ito + IKr + IKs + IK1 + IKb + ICaK
        d_V_li = self.membrane.d_V( I_ion + self.stimulus.I )
        
        if self.current_response_info:
            current_timestep = [
                trace.Current(name='INa', value=INa),
                trace.Current(name='INaL', value=INaL),                
            ]
            self.current_response_info.currents.append(current_timestep)
            
        return d_V_li + d_CaMKt_li + \
                    d_INa_li + d_INaL_li + d_Ito_li + d_ICaL_li + d_IKr_li + d_IKs_li + d_IK1_li + \
                        d_ryr_li + \
                            d_sodium_li + d_potassium_li + d_calcium_li
    
    
    def response_diff_eq(self, t, y):
        
        if type(self.protocol) == pacing_protocol.PacingProtocol :
            if self.protocol.type=='AP':            
                face = self.protocol.pacing(t)
                self.stimulus.cal_stimulation(face) # Stimulus    
            
            elif self.protocol.type=='VC':
                y[0] = self.protocol.voltage_at_time(t)

        else:                         
            y[0] = self.protocol.get_voltage_at_time(t)
                    
        return self.differential_eq(t, y)
   



def main():
    start_time = time.time()
    ord2017 = ORD2017(None)    
    print("--- %s seconds ---"%(time.time()-start_time))


if __name__ == '__main__':
    main()