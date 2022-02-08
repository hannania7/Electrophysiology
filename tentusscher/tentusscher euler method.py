import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('font', family='malgun gothic')
plt.rcParams['axes.unicode_minus'] = False

R = 8314.472  
T = 310   
F = 96485.3415 
RTF = R*T/F
FRT = F/(R*T)
FFRT = F*F/(R*T)

current_time = 0
next_time = 1000
dt = 0.001
V = -86.2
m = 0
h = 0.75
j = 0.75
r = 0
s = 1
Xs = 0
d = 0
f = 1
fca = 1
Nai = 11.6
Cai = 0.0002
Casr = 0.2
Ki = 138.3
xr1 = 0
xr2 = 1
g = 1

a12 = [-86.2]
times = [current_time]
Cai2 = [0.0002]
INa2 = [0]
Ito2 = [0]
IKr2 = [0]
IKs2 = [0]
ICaL2 = [0]
iK12 = [0]

i = 0
while i <= next_time - current_time:
    GNa = 14.838

    m1 = 1 / (np.power(1 + np.exp((-56.86 - V) / 9.03), 2))
    # print('m1 : ', m1)
    h1 = 1 / (np.power(1 + np.exp((V + 71.55) / 7.43), 2))
    # print('h1 : ', h1)
    j1 = 1 / (np.power(1 + np.exp((V + 71.55) / 7.43), 2))
    am = 1 / (1 + np.exp((-60 - V) / 5))
    bm = 0.1 / (1 + np.exp((V + 35) / 5)) + 0.1 / (1 + np.exp((V - 50) / 200))
    Tm = am * bm

    if V < -40:
        ah = 0.057 * np.exp(-(V + 80) / 6.8)
    else:
        ah = 0
    if V < -40:
        Bh = 2.7 * np.exp(0.079 * V) + 310000 * np.exp(0.3485 * V)
    else:
        Bh = 0.77 / (0.13 * (1 + np.exp((V + 10.66) / -11.1)))   
        
    Th = 1 / (ah + Bh)

    if V < -40:
        aj = ((-25428 * np.exp(0.2444 * V) - 6.948e-6 * np.exp(-0.04391 * V)) * (V + 37.78)) / (1 + np.exp(0.311 * (V + 79.23)))
    else:
        aj = 0
    if V < -40:
        Bj = (0.02424 * np.exp(-0.01052 * V)) / (1 + np.exp(-0.1378 * (V + 40.14)))
    else:
        Bj = (0.6 * np.exp(0.057 * V)) / (1 + np.exp(-0.1 * (V + 32)))            

    Tj = 1 / (aj + Bj)

    dm = (m1 - m) / Tm
    dh = (h1 - h) / Th
    dj = (j1 - j) / Tj  
    Nao = 140
    ENa = RTF * np.log(Nao / Nai)

    INa = GNa * np.power(m,3) * h * j * (V - ENa)

    GK1 = 5.405
    Ko = 5.4
    Ek = RTF * np.log(Ko / Ki)
    aK1 = 0.1 / (1 + np.exp(0.06 * (V - Ek - 200)))
    bK1 = (3 * np.exp(0.0002 * (V - Ek + 100)) + np.exp(0.1 * (V - Ek - 10))) / ((1 + np.exp(-0.5 * (V - Ek))))
    xK1 = aK1 / (aK1 + bK1)
    iK1 = GK1 * xK1 * (V - Ek)

    Gtoendo = 0.073
    r1 = 1 / (1 + (np.exp((20 - V) / (6))))
    Tr = 9.5 * np.exp(np.power((V + 40),2) / -1800) + 0.8
    dr = (r1 - r) / Tr

    s1 = 1 / ((1 + (np.exp((V + 28) / 5))))
    Ts = (1000 * np.exp(np.power((V + 67),2) / -1000)) + 8
    ds = (s1 - s) / Ts
    Ko = 5.4
    Ek = RTF * np.log(Ko / Ki)
    Ito = Gtoendo * r * s * (V - Ek)

    GKr = 0.096
    Ko = 5.4
    xr11 = 1 / (1 + (np.exp((-26 - V) / 7)))
    xr21 = 1 / (1 + (np.exp((V + 88) / 24)))
    axr1 = 450 / (1 + (np.exp((-45 - V) / 10)))
    Bxr1 = 6 / (1 + (np.exp((V + 30) / 11.5)))
    Txr1 = axr1 * Bxr1
    axr2 = 3 / (1 + (np.exp((-60 - V) / 20)))
    Bxr2 = 1.12 / (1 + (np.exp((V - 60) / 20)))
    Txr2 = axr2 * Bxr2
    dxr1 = (xr11 - xr1) / Txr1
    dxr2 = (xr21 - xr2) / Txr2
    Ek = RTF * np.log(Ko / Ki)
    IKr = GKr * np.sqrt(Ko / 5.4) * xr1 * xr2 * (V - Ek)

    Gksepiendo = 0.245
    Xs1 = 1 / (1 + (np.exp((-5 - V) / (14))))
    axs = 1100 / (np.sqrt(1 + (np.exp((-10 - V) / (6)))))
    bxs = 1 / (1 + (np.exp((V - 60) / (20))))
    Txs = axs * bxs
    dXs = (Xs1 - Xs) / Txs
    
    Ko = 5.4
    pKNa = 0.03
    Nao = 140
    EKs = RTF * np.log((Ko + pKNa * Nao) / (Ki + pKNa * Nai))
    IKs = Gksepiendo * np.power(Xs,2) * (V - EKs)

    GCaL = 0.000175
    Tfca = 2
    d1 = 1 / (1 + np.exp((-5 - V) / 7.5))
    f1 = 1 / (1 + np.exp((V + 20) / 7))

    afca = 1 / (1 + np.power(Cai / 0.000325,8))
    bfca = 0.1 / (1 + np.exp((Cai - 0.0005) / 0.0001))
    rfca = 0.2 / (1 + np.exp((Cai - 0.00075) / 0.0008))
    fca1 = (afca + bfca + rfca + 0.23) / 1.46
    d_fca = (fca1 - fca) / Tfca
    ad = 1.4 / (1 + np.exp((-35 - V) / 13)) + 0.25
    Bd = 1.4 / (1 + np.exp((V + 5) / 5))
    rd = 1 / (1 + np.exp((50 - V) / 20))
    Td = ad * Bd + rd
    Tf = 1125 * np.exp(-(np.power(V + 27,2)) / 240) +  80 + 165 / (1 + np.exp((25 - V) / 10)) 
        
    dd = (d1 - d) / Td
    df = (f1 - f) / Tf
    
    if 0.01 * d_fca > 0 and V > -60:
        dfca = 0
    else:
        dfca = d_fca
        
    Cao = 2

    ICaL = GCaL * d * f * fca * 4 * V * FFRT * (Cai * np.exp(2 * V * FRT) - 0.341 * Cao) / (np.exp(2 * V * FRT) - 1)

    kNaCa = 1000
    z = 0.35 
    y = 2.5 
    KmNai = 87.5 
    KmCa = 1.38 
    ksat = 0.1
    
    Cao = 2
    Nao = 140
    up = np.exp(z * V * FRT) * (np.power(Nai,3)) * Cao - np.exp((z - 1) * V * FRT) * (np.power(Nao,3)) * Cai * y
    down = (np.power(KmNai,3) + np.power(Nao,3)) * (KmCa + Cao) * (1 + ksat * np.exp((z - 1) * V * FRT))
    INaCa = kNaCa * (up / down)

    PNaK = 1.362
    KmK = 1
    KmNa = 40
    
    Ko = 5.4
    INaK = ((((PNaK * Ko)  / (Ko + KmK)) * Nai) / (Nai + KmNa)) / (1 + 0.1245 * np.exp(-0.1 * V * FRT) + 0.0353 * np.exp(-V * FRT))
    
    KpCa = 0.0005
    GpCa = 0.825
    IpCa = (GpCa * Cai) / (KpCa + Cai)

    GpK = 0.0146
    Ko = 5.4
    
    Ek = RTF * np.log(Ko / Ki)
    IpK = GpK * ((V - Ek) / (1 + (np.exp((25 - V) / (5.98)))))

    GbCa = 0.000592
    
    Cao = 2
    ECa = 0.5 * RTF * np.log(Cao / Cai)
    IbCa = GbCa * (V - ECa)
    GbNa = 0.00029
    
    Nao = 140
    ENa = RTF * np.log(Nao / Nai)
    IbNa = GbNa * (V - ENa)

    Bufc = 0.15
    Kbufc = 0.001
    Bufsr = 10 
    Kbufsr = 0.3
    Vc = 0.016404
    Kup = 0.00025
    brel = 0.25
    Vsr = 0.001094
    Vleak = 8e-5
    Vmaxup = 0.000425
    arel = 0.016464
    crel = 0.008232
    Cm = 0.185

    if Cai < 0.00035:
        gm = 1 / (1 + np.power(Cai / 0.00035,6)) 
    else:
        gm = 1 / (1 + np.power(Cai / 0.00035,16))
    Tg = 2 
    d_g = (gm-g) / Tg
    if 0.01*d_g > 0 and V > -60:
        dg = 0 
    else: 
        dg = d_g       
        
    Ileak = Vleak * (Casr - Cai)
    Iup = Vmaxup / (1 + (np.power(Kup,2)) / (np.power(Cai,2)))
    Irel = ((arel * (np.power(Casr,2))) / (np.power(brel,2) + np.power(Casr,2)) + crel) * d * g
    Caibufc = 1 / (1 + (Bufc * Kbufc) / (np.power(Cai + Kbufc,2)))
    Casrbufsr = 1 / (1 + (Bufsr * Kbufsr) / (np.power(Casr + Kbufsr,2)))
    dCai = Caibufc * (Ileak - Iup + Irel - (ICaL + IbCa + IpCa - (2 * INaCa)) / (2 * Vc * F) * Cm)   
    dCasr = Casrbufsr * Vc / Vsr * (Iup - (Irel + Ileak))

    Iion = (INa
        + iK1
        + Ito
        + IKr
        + IKs
        + ICaL
        + INaCa
        + INaK
        + IpCa
        + IpK
        + IbCa
        + IbNa)
    

    if current_time >= 0 and current_time <= 1:
        i_stim = -52
        dV = -(Iion + i_stim)
    else:
        i_stim = 0
        dV = -(Iion + i_stim)
    # All : dNai = (-(INa + IbNa + 3 * INaK + 3 * INaCa) / (Vc * F)) * Cm 
    # All : dKi = (-((iK1 + Ito + IKr + IKs + IpK + i_stim) - (2 * INaK)) /(Vc * F)) * Cm
    # IKr : dKi = (-((IKr + i_stim) ) /(Vc * F)) * Cm     
    
    
    Cm = 0.185
    Vc = 0.016404
    
    dNai = (INa + IbNa + 3 * INaK + 3 * INaCa) / -(Vc * F) * Cm 
    dKi = ((iK1 + Ito + IKr + IKs + IpK + i_stim) - (2 * INaK)) /-(Vc * F) * Cm

    # update time
    current_time += dt

    # integrate
    V_Next = V + dV*dt
    m_Next = m + dm*dt 
    h_Next = h + dh*dt
    j_Next = j + dj*dt
    r_Next = r + dr*dt 
    s_Next = s + ds*dt
    Xs_Next = Xs + dXs*dt
    d_Next = d + dd*dt
    f_Next = f + df*dt
    fca_Next = fca + dfca*dt
    Nai_Next = Nai + dNai*dt
    Cai_Next = Cai + dCai*dt
    Casr_Next = Casr + dCasr*dt
    g_Next = g + dg*dt
    Ki_Next = Ki + dKi*dt
    xr1_Next = xr1 + dxr1*dt
    xr2_Next = xr2 + dxr2*dt

    # update values
    V = V_Next
    m = m_Next
    h = h_Next
    j = j_Next
    r = r_Next
    s = s_Next
    Xs = Xs_Next
    d = d_Next
    f = f_Next
    fca = fca_Next
    Nai = Nai_Next
    Cai = Cai_Next
    Casr = Casr_Next
    g = g_Next
    Ki = Ki_Next
    xr1 = xr1_Next
    xr2 = xr2_Next

    times.append(current_time)
    a12.append(V_Next)
    Cai2.append(Cai_Next)
    INa2.append(INa)    
    Ito2.append(Ito)
    IKr2.append(IKr)
    IKs2.append(IKs)
    ICaL2.append(ICaL)
    iK12.append(iK1)

np.shape(a12)
# Euler method call
plt.title("IKr + iK1 + Ito + IKs + ICaL + INaCa + INaK + IpCa + IpK + IbCa + IbNa + INa")
plt.subplot(2, 4, 1)
plt.plot(times, a12)
plt.title('Action potential')
plt.subplot(2, 4, 2)
plt.plot(times, Cai2)
plt.title('Cai')
plt.subplot(2, 4, 3)
plt.plot(times, INa2)
plt.title('INa')
plt.subplot(2, 4, 4)
plt.plot(times, Ito2)
plt.title('Ito')
plt.subplot(2, 4, 5)
plt.plot(times, IKr2)
plt.title('IKr')
plt.subplot(2, 4, 6)
plt.plot(times, IKs2)
plt.title('IKs')
plt.subplot(2, 4, 7)
plt.plot(times, ICaL2)
plt.title('ICaL')
plt.subplot(2, 4, 8)
plt.plot(times, iK12)
plt.title('iK1')
plt.show()