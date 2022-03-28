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
print(next_time/dt)
Vm = -91.364111716184041
Na_i = 6.8151373016677805
K_i = 139.0561998319366
Ca_2_tot_jnc = 0.00039795809096128409
Ca_2_tot_iz = 0.00025268299350756843
Ca_2_tot_blk = 0.0000667801021116907
Ca_jnc = 0.00039795809096128409
Ca_iz = 0.00025268299350756843
Ca_blk = 0.0000667801021116907
Ca_2_SRup = 0.56250295256576965
Ca_2_tot_SRrl = 1.9128092931058895
TnChCa2 = 0.107546349
CaMCa2 = 0.000116997
SRCa2 = 0.00092639
TnChCa3 = 0.107546349
CaMCa3 = 0.000116997
SRCa3 = 0.00092639
# Lb_jnc = 0.032868783807565416
# Lb_iz = 0.011563775702743251
# Hb_jnc = 0.226665087396588
# Hb_iz = 0.099458321251370288
# internal_ion_concentrations 
Ca_Total2 = 4.0180173572968586e-4
ATPi = 4.657102729020499
p_RP_Na = 0.3556412697995689
p_AP_Na = 1.779648367445368e-5
p_AI_Na = 0.40285968661346977
# sodium_current_ultra_slow_gate
y2 = 0.5861887862983165
p_RP_CaL = 0.9968480629364956
p_AP_CaL = 1.5445004166497696e-6
p_AI_CaL = 8.77325391245903e-4
p_U = 0.17246483915629204
p_UCa = 6.098246017787626e-5
p_C = 0.4250747299372254
# L_type_Ca_channel_ultra_slow_gate
y3 = 0.9985266538252986
# T_type_Ca_channel_y1_gate
y4 = 1.6882718240109127e-5
y5 = 0.8585352091865849
# time_independent_potassium_current_y_gate
#TIPCY.y = 0.6080573900752752
# rapid_time_dependent_potassium_current_y1_gate
#RTDPC1.y1 = 0.0018339931180983765
#RTDPC2.y2 = 0.20443083454225305
#RTDPC3.y3 = 0.967887666264921
# slow_time_dependent_potassium_current_y1_gate
#STDPC1.y1 = 0.09738789658609195
#STDPC2.y2 = 0.09745345578743213
# transient_outward_current_y1_gate
#TOC1.y1 = 7.956883250874798e-4
#TOC2.y2 = 0.9999125083105881
# sodium_calcium_exchanger_y_gate
y6 = 0.9891789193465331
# sodium_potassium_pump_y_gate
y7 = 0.5910747147428818
# SR_calcium_pump_y_gate
y8 = 0.46108441538480216
p_open_RyR = 3.4314360001543243e-4
p_close_RyR = 0.19135178123107768
# Ca_concentrations_in_SR
Ca_Total3 = 9.455741736977666
Caup = 2.611712901567567
X = 0.9573749975411884
pCa = 0.02490898775497523
pCaCB = 0.001990153835322864
pCB = 4.2941813853474524e-4
p_O_NaT2 = 0.000000706725155695262
p_O_NaT3 = 0.000000706725155695262
p_I_2_NaT = 0.0117704053067285
p_I_s_NaT = 0.304002781414015
p_O_NaL2 = 0.00000295214591324261
p_O_NaL3 = 0.00000295214591324261
p_I_1_NaL = 0.00254273877063925
p_I_2_NaL = 0.0118261382165599
p_I_s_NaL = 0.303220346353844
chi_r_fast = 0.00000486210633393005
chi_r_slow = 0.437041249050081
P_7 = 0.0831770174499825
P_8_13 = 0.281082409575779
P_1_6 = 0.435289193632868
p_E_1_NCX_blk = 0.111872123711613
p_I_1_NCX_blk = 0.203023555446362
p_I_2_NCX_blk = 0.684869019924837
p_E_1_NCX_iz = 0.238718640001014
p_I_1_NCX_iz = 0.13771129457898
p_I_2_NCX_iz = 0.622892868847556
Y_ooo = 0.00000172489315884865
Y_ooc = 0.00000142034754677507
Y_coo = 0.0000138422676498755
Y_coc = 0.992110534408681
Y_cco = 0.0000000953816272498217
Y_oco = 0.00000000000156949238162028
Y_occ = 0.0000000249594301562175
Y_co_iz = 0.992251726297519
Y_oo_iz = 0.00000314564543512061
Y_oc_iz = 0.000000024556270151713
Y_co_blk = 0.992424981547859
Y_oo_blk = 0.00000314619469048683
Y_oc_blk = 0.0000000240070147854924
TSCa_3 = 0.0086895187306879739
TSCa_3W = 0.0003391299963111013    
TSCa_3S = 0.00013952053460170792    
TS_S = 0.000876347322180234
TS_W = 0.000492054058977473
hw = 0.000100147615113241
hp = 0.00600014761511324
Pb_spm2 = 0.594875991179992
Pb_spm3 = 0.594875991179992
# random 
p_C_NaT = 0
p_C_NaL = 0
y9 = 0
y10 = 0
y11 = 0
p_O_v = 0
p_blk_O_c = 0
p_blk_C_2 = 0
p_iz_O_c = 0
p_iz_C_2 = 0
y_1Kto = 0
y_2Kto = 0


times = []
Vm2 = []
# Na_i2 = []
# K_i2 = []
# Ca_2_tot_jnc2 = []
# Ca_2_tot_iz2 = []
# Ca_2_tot_blk2 = []
# Ca_jnc2 = []
# Ca_iz2 = []
# Ca_blk2 = []
# Ca_2_SRup2 = []
# Ca_2_tot_SRrl2 = []
# TnChCa2 = []
# CaMCa2 = []
# SRCa2 = []
# TnChCa3 = []
# CaMCa3 = []
# SRCa3 = []
# # Lb_jnc = 0.032868783807565416
# # Lb_iz = 0.011563775702743251
# # Hb_jnc = 0.226665087396588
# # Hb_iz = 0.099458321251370288
# # internal_ion_concentrations 
# Ca_Total2 = []
# ATPi2 = []
# p_RP_Na2 = []
# p_AP_Na2 = []
# p_AI_Na2 = []
# # sodium_current_ultra_slow_gate
# SCUSG.y2 = []
# p_RP_CaL2 = []
# p_AP_CaL2 = []
# p_AI_CaL2 = []
# p_U2 = []
# p_UCa2 = []
# p_C2 = []
# # L_type_Ca_channel_ultra_slow_gate
# y3 = []
# # T_type_Ca_channel_y1_gate
# y4 = []
# y5 = []
# # time_independent_potassium_current_y_gate
# #TIPCY.y = 0.6080573900752752
# # rapid_time_dependent_potassium_current_y1_gate
# #RTDPC1.y1 = 0.0018339931180983765
# #RTDPC2.y2 = 0.20443083454225305
# #RTDPC3.y3 = 0.967887666264921
# # slow_time_dependent_potassium_current_y1_gate
# #STDPC1.y1 = 0.09738789658609195
# #STDPC2.y2 = 0.09745345578743213
# # transient_outward_current_y1_gate
# #TOC1.y1 = 7.956883250874798e-4
# #TOC2.y2 = 0.9999125083105881
# # sodium_calcium_exchanger_y_gate
# y6 = []
# # sodium_potassium_pump_y_gate
# y7 = []
# # SR_calcium_pump_y_gate
# y8 = []
# p_open_RyR2 = []
# p_close_RyR2 = []
# # Ca_concentrations_in_SR
# Ca_Total3 = []
# Caup2 = []
# X2 = []
# pCa2 = []
# pCaCB2 = []
# pCB2 = []
# p_O_NaT2 = []
# p_O_NaT3 = []
# p_I_2_NaT2 = []
# p_I_s_NaT2 = []
# p_O_NaL2 = []
# p_O_NaL3 = []
# p_I_1_NaL2 = []
# p_I_2_NaL2 = []
# p_I_s_NaL2 = []
# chi_r_fast2 = []
# chi_r_slow2 = []
# P_72 = []
# P_8_132 = []
# P_1_62 = []
# p_E_1_NCX_blk2 = []
# p_I_1_NCX_blk2 = []
# p_I_2_NCX_blk2 = []
# p_E_1_NCX_iz2 = []
# p_I_1_NCX_iz2 = []
# p_I_2_NCX_iz2 = []
# Y_ooo2 = []
# Y_ooc2 = []
# Y_coo2 = []
# Y_coc2 = []
# Y_cco2 = []
# Y_oco2 = []
# Y_occ2 = []
# Y_co_iz2 = []
# Y_oo_iz2 = []
# Y_oc_iz2 = []
# Y_co_blk2 = []
# Y_oo_blk2 = []
# Y_oc_blk2 = []
# TSCa_32 = []
# TSCa_3W2 = []
# TSCa_3S2 = []
# TS_S2 = []
# TS_W2 = []
# hw2 = []
# hp2 = []
# Pb_spm2 = []
# Pb_spm3 = []
# # random 
# p_C_NaT2 = []
# p_C_NaL2 = []
# TGP.y9 = []
# TGP.y10 = []
# TGP.y11 = []
# p_O_v2 = []
# p_blk_O_c2 = []
# p_blk_C_22 = []
# p_iz_O_c2 = []
# p_iz_C_22 = []
# y_1Kto2 = []
# y_2Kto2 = []

while current_time < next_time:
    R = 8.3143
    T = 310
    F = 96.4867
    k_off_CaM = 0.238
    k_on_CaM = 34
    k_off_TnCh = 0.000032
    k_on_TnCh = 2.37
    k_off_SR = 0.06
    k_on_SR = 100
    k_off_L_iz = 1.3
    k_on_L_iz = 100
    k_off_H_iz = 0.03
    k_on_H_iz = 100
    k_off_L_jnc = 1.3
    k_on_L_jnc = 100
    k_off_H_jnc = 0.03
    k_on_H_jnc = 100
    k_off_CSQN = 65
    k_on_CSQN = 100
    k_on_CaM = 34
    B_tot_CaM = 0.024
    k_off_CaM = 0.238
    k_on_TnCh = 2.37
    B_tot_TnCh = 0.12
    k_off_TnCh = 0.000032
    k_on_SR = 100
    B_tot_SR = 0.0171
    k_off_SR = 0.06
    L_bound_iz = 0.0075621764602356
    B_tot_L_iz = 0.6078
    k_off_L_iz = 1.3
    k_on_L_iz = 100
    H_bound_iz = 0.0769149150028914
    B_tot_H_iz = 0.2178
    k_off_H_iz = 0.03
    k_on_H_iz = 100
    B_tot_L_jnc = 1.1095
    Ca_2_jnc = 0
    k_off_L_jnc = 1.3
    k_on_L_jnc = 100
    B_tot_H_jnc = 0.398
    k_off_H_jnc = 0.03
    k_on_H_jnc = 100
    B_tot_CSQN = 3
    k_off_CSQN = 65
    k_on_CSQN = 100
    G_dCa_jnciz = 3395.88
    Sc_Cell = 1
    G_dCa_izblk = 3507.78
    P_trans = 4.8037
    T_L = 147.51
    K_L = 0.0044
    ATP = 6
    P_CaL_Ca = 14.21
    Ca_o = 1.8
    Na_o = 140
    K_o = 4.5
    f_CaL_blk = 0.1
    f_CaL_iz = 0.15
    f_CaL_jnc = 0.75
    k_I2O = 0.0001312
    k_I1I2 = 0.00534
    k_I1O = 0.01
    f_LSM = 0.13125
    P_Na = 8.1072
    Mg_2_cyt = 0.8
    SPM = 5.0
    G_K1 = 1.353
    f_mode1 = 0.9
    G_Kr = 0.0166
    f_Ks_iz = 0.1
    f_Ks_blk = 0.9
    P_Ks_K = 0.002782
    G_Kto = 0.0312
    P_Kpl = 0.0000172
    f_Cab_blk = 0.9
    P_Cab = 0.00006822
    f_Cab_iz = 0.1
    P_bNSC_K = 0.00014
    P_bNSC_Na = 0.00035
    f_l_Ca_blk = 0.9
    P_l_Ca_Na = 0.00273
    f_l_Ca_iz = 0.1
    ATP_cyt = 6.67701543987464
    G_KATP = 17.674
    delta_Nai = -0.14
    K_d_Nai_0 = 5
    delta_Nao = 0.44
    K_d_Nao_0 = 26.8
    delta_Ki = -0.14
    K_d_Ki_0 = 18.8
    delta_Ko = 0.23
    K_d_Ko_0 = 0.8
    K_d_MgATP = 0.6
    MgATP_cyt = 6.631643709767415
    k_1_plus = 0.72
    k_1_minus = 0.08
    k_2_plus = 0.08
    k_2_minus = 0.008
    k_3_plus = 4
    k_3_minus = 8000
    k_4_plus = 0.3
    k_4_minus = 0.2
    Pi = 0.50872066859173026
    H = 0.0001
    Amp_NaK = 25.178
    Stoi_NaK_Na = 3
    Stoi_NaK_K = -2
    K_m_act = 0.004
    alpha_1_on = 0.002
    alpha_1_off = 0.0015
    alpha_2_on = 0.00006
    alpha_2_off = 0.02
    beta_1_on = 0.0012
    beta_1_off = 0.0000005
    beta_2_on = 0.18
    beta_2_off = 0.0002
    K_m_Nai = 20.74854
    K_m_Nao = 87.5
    K_m_Cai = 0.0184
    K_m_Cao = 1.38
    k_3 = 1.0
    k_4 = 1.0
    Amp_NCX = 30.53
    f_NCX_blk = 0.9
    f_NCX_iz = 0.1
    K_m = 0.0005
    Amp_PMCA = 0.19
    f_PMCA_blk = 0.9
    f_PMCA_iz = 0.1
    J_L = 0.000913
    J_R = 0.02
    g_D = 0.065
    Q_10 = 3
    sloc0 = 0.1
    f_n = 7
    N_RyR = 10
    p_O_RyR_base = 0.000075
    P_RyR = 5191
    K_dCai = 0.0027
    K_dCasr = 1.378
    Amp_SERCA = 106.4448
    I_Kto_Na = 0
    I_app = 0
    halfSL = 0.91
    TS_tot = 23
    propFh = 28000
    Za = 0.0023
    Yv = 1.5
    Yd = 0.0333
    Yc = 1
    Lc = 1.2
    Zb = 0.1397
    Yb = 0.1816
    rate_f = 0.0023
    convertF = 15
    eqvhalfSL = 1.15
    Zp = 0.2095
    Yp = 0.1397
    Zr = 7.2626
    Yr = 0.1397
    Zq = 0.3724
    Yq = 0.2328
    hwr = 0.0001
    rate_B = 0.5
    hpr = 0.006
    dt = 0.01
    MgADP = 36.3
    Cm = 177
    # from huvec2
    Ca_2_blk= Ca_2_tot_blk - (CaMCa2 + TnChCa2 + SRCa2 + 3 * (TSCa_3 + TSCa_3W + TSCa_3S) / 1000)
    dCaMCa2 = k_on_CaM * Ca_2_blk * (B_tot_CaM - CaMCa2) - k_off_CaM * CaMCa2
    dTnChCa2 = k_on_TnCh * Ca_2_blk * (B_tot_TnCh - TnChCa2) - k_off_TnCh * TnChCa2
    dSRCa2 = k_on_SR * Ca_2_blk * (B_tot_SR - SRCa2) - k_off_SR * SRCa2
 
    # from huvec2
    K_dL_iz = k_off_L_iz / k_on_L_iz
    K_dH_iz = k_off_H_iz / k_on_H_iz
    K_dL_jnc = k_off_L_jnc / k_on_L_jnc
    K_dH_jnc = k_off_H_jnc / k_on_H_jnc
    K_d_CSQN_Ca = k_off_CSQN / k_on_CSQN
    P_CaL_Na = 0.0000185 * P_CaL_Ca
    P_CaL_K = 0.000367 * P_CaL_Ca
    P_Ks_Na = 0.04 * P_Ks_K
    P_l_Ca_K = P_l_Ca_Na
    alpha_2_plus = k_2_plus
    delta_RTF = 2 * F / (R * T)
    f_L = J_L / g_D
    f_R = J_R / g_D
    k_oc = Q_10 * 0.5664
    Cm = 192.46 * Sc_Cell
    V_cell = 120 * (37.62 * Sc_Cell) * 8.4
    V_jnc = 0.008 * V_cell
    V_iz = 0.035 * V_cell
    V_blk = 0.68 * V_cell
    V_SRt = 0.06 * V_cell
    V_SRrl = 0.2 * V_SRt
    V_SRup = 0.8 * V_SRt
    V_cyt = V_jnc + V_iz + V_blk
    Yvd = Yv

    # from huvec2
    TS = TS_tot - TSCa_3 - TSCa_3W - TSCa_3S - TS_S - TS_W
    rate_g = Za + Yv * (1 - np.exp(-propFh * np.power((hw - hwr),2)))
    rate_gd = Yd + Yc * np.power((halfSL - Lc),2) + Yvd * (1 - np.exp(-propFh * np.power((hw - hwr),2)))
    dTSCa_3 = Yb * TS * np.power((Ca_2_blk * 1000),3) - Zb * TSCa_3 + rate_g * TSCa_3W - rate_f * np.exp(-convertF * np.power((halfSL - eqvhalfSL),2)) * TSCa_3
    dTSCa_3W = rate_f * np.exp(-convertF * np.power((halfSL - eqvhalfSL),2)) * TSCa_3 - rate_g * TSCa_3W + Zp * TSCa_3S - Yp * TSCa_3W
    dTSCa_3S = Yp * TSCa_3W - Zp * TSCa_3S + Zr * TS_S * np.power((Ca_2_blk * 1000),3) - Yr * TSCa_3S
    dTS_S = Yr * TSCa_3S - Zr * TS_S * np.power((Ca_2_blk * 1000),3) + Zq * TS_W - Yq * TS_S
    dTS_W = Yq * TS_S - Zq * TS_W - rate_gd * TS_W
    dhw = -rate_B * (hw - hwr)
    dhp = -rate_B * (hp - hpr)

    # from huvec2
    L_free_iz = B_tot_L_iz - L_bound_iz
    H_free_iz = B_tot_H_iz - H_bound_iz
    for n in range(0,10):
        Ca_2_iz = Ca_2_tot_iz / (1 + L_free_iz / K_dL_iz + H_free_iz / K_dH_iz)
        L_free_iz = B_tot_L_iz / (1 + Ca_2_iz / K_dL_iz)
        H_free_iz = B_tot_H_iz / (1 + Ca_2_iz / K_dH_iz)
    L_bound_iz = B_tot_L_iz - L_free_iz
    H_bound_iz = B_tot_H_iz - H_free_iz

    # from huvec2
    for n in range(0,10):
        L_free_jnc = B_tot_L_jnc / (1 + Ca_2_jnc / K_dL_jnc)
        H_free_jnc = B_tot_H_jnc / (1 + Ca_2_jnc / K_dH_jnc)
        Ca_2_jnc = Ca_2_tot_jnc / (1 + L_free_jnc / K_dL_jnc + H_free_jnc / K_dH_jnc)
    # print(Ca_2_jnc)
    # from huvec2
    #[IC]

    #dCa_2_tot_jnc = -I_tot_Ca_jnc * Cm / (V_jnc * 2 * F) + J_Ca_rel / V_jnc - J_Ca_jnciz / V_jnc
    #dCa_2_tot_iz = -I_tot_Ca_iz * Cm / (V_iz * 2 * F) + J_Ca_jnciz / V_iz - J_Ca_izblk / V_iz
    #dCa_2_tot_blk = -I_tot_Ca_blk * Cm / (V_blk * 2 * F) - J_SERCA / V_blk + J_Ca_izblk / V_blk
    #dCa_2_SRup = J_SERCA / V_SRup - J_trans_SR / V_SRup
    #dCa_2_tot_SRrl = J_trans_SR / V_SRrl - J_Ca_rel / V_SRrl
    #dNai = -I_tot_Na * Cm / (V_cyt * F)
    #dKi = -(I_tot_K + I_app) * Cm / (V_cyt * F)





    # from huvec2
    a = 1
    b = B_tot_CSQN - Ca_2_tot_SRrl + K_d_CSQN_Ca
    c = -K_d_CSQN_Ca * Ca_2_tot_SRrl
    Ca_2_SRrl = (-b + np.sqrt(np.power((b),2) - 4 * a * c)) / (2 * a)
    # from huvec2
    alpha_plus = 1 / (3.734 * np.exp(-Vm / 8.5) + 0.35 * np.exp(-Vm / 3500))
    # print(alpha_plus)
    alpha_minus = 1 / (4.65 * np.exp(Vm / 15) + 1.363 * np.exp(Vm / 100))
    epsilon_minus = 1 / (8084 * np.exp(Vm / 10) + 158 * np.exp(Vm / 1000)) + 1 / (134736 * np.exp(-Vm / 5) + 337 * np.exp(-Vm / 2000))

    # Cell volume
    Vcell = 115 * 20 * 8
    Vblk = 0.65 * Vcell
    Viz = 0.035 * Vcell
    Vjnc = 0.011 * Vcell
    Vcyt = Vjnc + Viz + Vblk
    VSR = 0.035 * Vcell
    VSRup = 0.6 * VSR
    VSRrl = 0.4 * VSR

    # external_ion_concentrations
    Na_o = 140
    Ca_o = 1.8
    K_o = 4.5

    # Substrates (fixed)
    MgATPcyt = 6.677
    MgADPcyt = 0.0228
    Picyt = 0.443

    # internal_ion_concentrations
    CMDN_max = 0.05
    K_mCMDN = 0.00238
    Vi = 8000
    Cm = 177
    
    # Boundary ion diffusion
    G_dCa_Jnciz = 3215.8
    G_dCa_izblk = 2076.1
    J_Ca_jnciz = G_dCa_Jnciz * (Ca_jnc - Ca_iz)
    J_Ca_izblk = G_dCa_izblk * (Ca_iz - Ca_blk)


    # Ca2+ buffering
    B_Total_CaM = 0.024
    B_Total_TnCh = 0.12
    B_Total_SR = 0.0171
    B_Total_H_iz = 0.2178
    B_Total_L_iz = 0.6078
    B_Total_H_jnc = 0.398
    B_Total_L_jnc = 1.1095
    B_Total_CSQN = 3
    b1 = (CMDN_max - Ca_Total2) + K_mCMDN
    c1 = K_mCMDN * Ca_Total2
    Ca_i = (np.sqrt(np.power(b1,2) + 4 * c1) - b1) / 2

    # L_type_Ca_channel
    P_CaL = 8712
    p_open_CaL = (p_AP_CaL * (p_U + p_UCa) * y3) / (np.power(1 + (1.4 / ATPi),3))

    # constant_field_equations(F * Vm) / (R * T)) * (Na_i - (Na_o * np.exp((-F * Vm) / (R * T)))) / (1 - np.exp((-F * Vm) / (R * T))))
    if Vm == 0 :
        CF_Na = Na_o
        CF_Ca = -Ca_o
        CF_K = K_i
    else:
        CF_Na = ((F * Vm) / (R * T)) * (Na_i - (Na_o * np.exp((-F * Vm) / (R * T)))) / (1 - np.exp((-F * Vm) / (R * T)))
        CF_Ca = ((2 * F * Vm) / (R * T)) * (Ca_i - (Ca_o * np.exp((-2 * F * Vm) / (R * T)))) / (1 - np.exp((-2 * F * Vm) / (R * T)))
        CF_K = ((F * Vm) / (R * T)) * (K_i - (K_o * np.exp((-F * Vm) / (R * T)))) / (1 - np.exp((-F * Vm) / (R * T)))
    E_K = ((R * T) / F) * np.log(K_o / K_i)
    i_CaL_Ca = P_CaL * CF_Ca * p_open_CaL
    i_CaL_Na = 0.0000185 * P_CaL * CF_Na * p_open_CaL
    i_CaL_K = 0.000365 * P_CaL * CF_K * p_open_CaL

    #i_net_Na = i_Na_Na + i_Ks_Na + TOC.i_to_Na + LTi_CaL_Na + BNC.i_bNSC_Na + BIC.i_lCa_Na + 3 * i_NaK + 3 * i_NaCa
    #i_net_K = i_K1 + RTDPC.i_Kr + TOC.i_to_K + BKATP.i_KATP + i_Ks_K + i_Na_K + LTi_CaL_K + BNC.i_bNSC_K + BIC.i_lCa_K + BKPI.i_Kpl - (2 * i_NaK)

    # T_type_Ca_channel_y1_gate
    alpha_y1 = 1 / ((0.019 * np.exp(-Vm / 5.6)) + (0.82 * np.exp(-Vm / 250)))
    beta_y1 = 1 / ((40 * np.exp(Vm / 6.3)) + (1.5 * np.exp(Vm / 10000)))
    dy4 = alpha_y1 * (1 - y4) - (beta_y1 * y4)

    # T_type_Ca_channel_y2_gate
    alpha_y2 = 1 / ((62000 * np.exp(Vm / 10.1)) + (30 * np.exp(Vm / 3000)))
    beta_y2 = 1 / (0.0006 * np.exp(-Vm / 6.7) + (1.2 * np.exp(-Vm / 25)))
    dy5 = alpha_y2 * (1 - y5) - (beta_y2 * y5)

    # T_type_Ca_channel
    P_CaT = 612
    i_Ca_T = P_CaT * CF_Ca * y4 * y5
    # from huvec2
    p_O_RyR = Y_ooo + Y_coo + Y_cco + Y_oco
    p_O_RyR_t = p_O_RyR + p_O_RyR_base

    J_Ca_rel = P_RyR * p_O_RyR_t * (Ca_2_SRrl - Ca_2_jnc) * Sc_Cell
    Ca_2_nd_00 = Ca_2_jnc
    Ca_2_nd_L0 = (Ca_2_nd_00 + f_L * delta_RTF * Vm * np.exp(-delta_RTF * Vm) / (1 - np.exp(-delta_RTF * Vm)) * Ca_o) / (1 + f_L * delta_RTF * Vm / (1 - np.exp(-delta_RTF * Vm)))
    Ca_2_nd_0R = (Ca_2_nd_00 + f_R * Ca_2_SRrl) / (1 + f_R)
    Ca_2_nd_LR = (Ca_2_nd_00 + f_R * Ca_2_SRrl + f_L * delta_RTF * Vm * np.exp(-delta_RTF * Vm) / (1 - np.exp(-delta_RTF * Vm)) * Ca_o) / (1 + f_R + f_L * delta_RTF * Vm / (1 - np.exp(-delta_RTF * Vm)))
    epsilon_plus_00 = (Ca_2_nd_00 * alpha_plus) / (T_L * K_L)
    epsilon_plus_L0 = (Ca_2_nd_L0 * alpha_plus) / (T_L * K_L)
    epsilon_plus_0R = (Ca_2_nd_0R * alpha_plus) / (T_L * K_L)
    epsilon_plus_LR = (Ca_2_nd_LR * alpha_plus) / (T_L * K_L)

    k_co_00 = Q_10 * 0.4 / (1 + np.power((0.025 / Ca_2_nd_00),2.7))
    k_co_L0 = Q_10 * 0.4 / (1 + np.power((0.025 / Ca_2_nd_L0),2.7))
    k_co_0R = Q_10 * 0.4 / (1 + np.power((0.025 / Ca_2_nd_0R),2.7))
    k_co_LR = Q_10 * 0.4 / (1 + np.power((0.025 / Ca_2_nd_LR),2.7))
    f_t_00 = k_co_00 / (k_co_00 + k_oc)
    f_t_L0 = k_co_L0 / (k_co_L0 + k_oc)
    k_rco_0 = f_n * f_t_00 * k_co_0R * (sloc0 + Ca_2_SRrl)
    k_rco_L = f_n * f_t_L0 * k_co_LR * (sloc0 + Ca_2_SRrl)
    p_C_0 = k_oc / (k_oc + f_t_00 * (k_rco_0 / (f_n * f_t_00)))
    p_C_L = k_oc / (k_oc + f_t_00 * (k_rco_L / (f_n * f_t_L0)))
    k_roc_0 = k_oc * np.power((p_C_0),(N_RyR - 1) * 0.74)
    k_roc_L = k_oc * np.power((p_C_L),(N_RyR - 1) * 0.74)    
    Y_ccc = 1 - (Y_ooo + Y_ooc + Y_coo + Y_coc + Y_cco + Y_oco + Y_occ)
    dY_ooo = k_rco_L * Y_ooc + alpha_plus * Y_coo + epsilon_minus * Y_oco - (k_roc_L + alpha_minus + epsilon_plus_LR) * Y_ooo
    dY_ooc = alpha_plus * Y_coc + k_roc_L * Y_ooo + epsilon_minus * Y_occ - (alpha_minus + k_rco_L + epsilon_plus_L0) * Y_ooc
    dY_coo = k_rco_0 * Y_coc + alpha_minus * Y_ooo + epsilon_minus * Y_cco - (k_roc_0 + alpha_plus + epsilon_plus_0R) * Y_coo
    dY_coc = k_roc_0 * Y_coo + alpha_minus * Y_ooc + epsilon_minus * Y_ccc - (k_rco_0 + alpha_plus + epsilon_plus_00) * Y_coc
    dY_cco = k_rco_0 * Y_ccc + alpha_minus * Y_oco + epsilon_plus_0R * Y_coo - (k_roc_0 + alpha_plus + epsilon_minus) * Y_cco
    dY_oco = k_rco_0 * Y_occ + alpha_plus * Y_cco + epsilon_plus_LR * Y_ooo - (k_roc_0 + alpha_minus + epsilon_minus) * Y_oco
    dY_occ = alpha_plus * Y_ccc + k_roc_0 * Y_oco + epsilon_plus_L0 * Y_ooc - (alpha_minus + k_rco_0 + epsilon_minus) * Y_occ



    epsilon_plus_iz = (Ca_2_iz * alpha_plus) / (T_L * K_L)
    epsilon_plus_blk = (Ca_2_blk * alpha_plus) / (T_L * K_L)
    Ca_2_iz_loc = (Ca_2_iz + f_L * delta_RTF * Vm * np.exp(-delta_RTF * Vm) / (1 - np.exp(-delta_RTF * Vm)) * Ca_o) / (1 + f_L * delta_RTF * Vm / (1 - np.exp(-delta_RTF * Vm)))
    Ca_2_blk_loc = (Ca_2_blk + f_L * (delta_RTF * Vm * np.exp(-delta_RTF * Vm)) / (1 - np.exp(-delta_RTF * Vm)) * Ca_o) / (1 + f_L * delta_RTF * Vm / (1 - np.exp(-delta_RTF * Vm)))
    epsilon_plus_iz_loc = (Ca_2_iz_loc * alpha_plus) / (T_L * K_L)
    epsilon_plus_blk_loc = (Ca_2_blk_loc * alpha_plus) / (T_L * K_L)



    p_O_LCC = Y_ooo + Y_ooc
    ATP = 6
    ATPfactor = 1 / (np.power(1 + (1.4 / ATP),3))

    GHK_Ca_LR = 2 * Vm / (R * T / F) * (Ca_2_nd_LR - Ca_o * np.exp(-2 * Vm / (R * T / F))) / (1 - np.exp(-2 * Vm / (R * T / F)))
    GHK_Ca_L0 = 2 * Vm / (R * T / F) * (Ca_2_nd_L0 - Ca_o * np.exp(-2 * Vm / (R * T / F))) / (1 - np.exp(-2 * Vm / (R * T / F)))
    GHK_Ca_iz = 2 * Vm / (R * T / F) * (Ca_2_iz - Ca_o * np.exp(-2 * Vm / (R * T / F))) / (1 - np.exp(-2 * Vm / (R * T / F)))
    GHK_Ca_blk = 2 * Vm / (R * T / F) * (Ca_2_blk - Ca_o * np.exp(-2 * Vm / (R * T / F))) / (1 - np.exp(-2 * Vm / (R * T / F)))
    GHK_Na = ((1 * F * Vm) / (R * T)) * (Na_i - (Na_o * np.exp((-1 * F * Vm) / R * T) / (1 - np.exp(-1 * F * Vm / R * T)))) 
    GHK_K = ((1 * F * Vm) / (R * T)) * (K_i - (K_o * np.exp((-1 * F * Vm) / R * T) / (1 - np.exp(-1 * F * Vm / R * T)))) 
    GHK_Ca = ((1 * F * Vm) / (R * T)) * (Ca_i - (Ca_o * np.exp((-1 * F * Vm) / R * T) / (1 - np.exp(-1 * F * Vm / R * T))))  

    I_CaL_Ca_blk = f_CaL_blk * P_CaL_Ca * GHK_Ca_blk * Y_oo_blk * ATPfactor
    I_CaL_Ca_iz = f_CaL_iz * P_CaL_Ca * GHK_Ca_iz * Y_oo_iz * ATPfactor
    print(GHK_Ca_iz)
    I_CaL_Ca_LR = f_CaL_jnc * P_CaL_Ca * GHK_Ca_LR * Y_ooo * ATPfactor
    I_CaL_Ca_L0 = f_CaL_jnc * P_CaL_Ca * GHK_Ca_L0 * Y_ooc * ATPfactor
    I_CaL_Ca_jnc = I_CaL_Ca_LR + I_CaL_Ca_L0    
    I_CaL_Na_blk = f_CaL_blk * P_CaL_Na * GHK_Na * Y_oo_blk * ATPfactor
    I_CaL_Na_iz = f_CaL_iz * P_CaL_Na * GHK_Na * Y_oo_iz * ATPfactor
    I_CaL_Na_jnc = f_CaL_jnc * P_CaL_Na * GHK_Na * p_O_LCC * ATPfactor
    I_CaL_K_blk = f_CaL_blk * P_CaL_K * GHK_K * Y_oo_blk * ATPfactor
    I_CaL_K_iz = f_CaL_iz * P_CaL_K * GHK_K * Y_oo_iz * ATPfactor
    I_CaL_K_jnc = f_CaL_jnc * P_CaL_K * GHK_K * p_O_LCC * ATPfactor

    I_CaL = (I_CaL_Ca_jnc + I_CaL_Na_jnc + I_CaL_K_jnc) + (I_CaL_Ca_iz + I_CaL_Na_iz + I_CaL_K_iz) + (I_CaL_Ca_blk + I_CaL_Na_blk + I_CaL_K_blk) 

    Y_cc_iz = 1 - (Y_co_iz + Y_oo_iz + Y_oc_iz)
    dY_co_iz = epsilon_minus * Y_cc_iz + alpha_minus * Y_oo_iz - (epsilon_plus_iz + alpha_plus) * Y_co_iz
    dY_oo_iz = alpha_plus * Y_co_iz + epsilon_minus * Y_oc_iz - (alpha_minus + epsilon_plus_iz_loc) * Y_oo_iz
    dY_oc_iz = epsilon_plus_iz_loc * Y_oo_iz + alpha_plus * Y_cc_iz - (epsilon_minus + alpha_minus) * Y_oc_iz

    Y_cc_blk = 1 - (Y_co_blk + Y_oo_blk + Y_oc_blk)
    dY_co_blk = epsilon_minus * Y_cc_blk + alpha_minus * Y_oo_blk - (epsilon_plus_blk + alpha_plus) * Y_co_blk
    dY_oo_blk = alpha_plus * Y_co_blk + epsilon_minus * Y_oc_blk - (alpha_minus + epsilon_plus_blk_loc) * Y_oo_blk
    dY_oc_blk = epsilon_plus_blk_loc * Y_oo_blk + alpha_plus * Y_cc_blk - (epsilon_minus + alpha_minus) * Y_oc_blk

    # Background calcium current
    f_Cab_blk = 0.9
    f_Cab_iz = 0.1
    I_Cab_blk = 0.0006822 * f_Cab_blk * GHK_Ca
    I_Cab_iz = 0.0006822 * f_Cab_iz * GHK_Ca
    I_Cab = I_Cab_iz + I_Cab_blk

    # background_Cab_current
    P_Cab = 0.04
    i_Cab = P_Cab * CF_Ca

    # sodium_calcium_exchanger
    P_NaCa = 6.81
    k3 = 1
    k4 = 1
    Km_Nai = 8.75
    Km_Nao = 87.5
    Km_Cai = 0.00138
    Km_Cao = 1.38
    Partition = 0.32
    p_E1Na = 1 / (1 + (np.power((Km_Nai / Na_i),3) * (1 + (Ca_i / Km_Cai))))
    p_E2Na = 1 / (1 + (np.power((Km_Nao / Na_o),3) * (1 + (Ca_o / Km_Cao))))
    p_E1Ca = 1 / (1 + ((Km_Cai / Ca_i) * (1 + np.power((Na_i / Km_Nai),3))))
    p_E2Ca = 1 / (1 + ((Km_Cao / Ca_o) * (1 + np.power((Na_o / Km_Nao),3))))
    k1 = np.exp((Partition * F * Vm) / (R * T))
    k2 = np.exp(((Partition - 1) * F * Vm) / (R * T))
    i_NaCa = P_NaCa * Cm * (k1 * p_E1Na * y6 - (k2 * p_E2Na * (1 - y6)))

    I_tot_Ca_jnc = I_CaL_Ca_jnc
    # Plasma membrane Ca2+-ATPase current
    f_PMCA_iz = 0.1
    f_PMCA_blk = 0.9
    Amp_PMCA = 0.1615
    K_m = 0.0005
    I_PMCA_blk = f_PMCA_blk * Amp_PMCA * (np.power((Ca_blk),1.6)) / (np.power((K_m),1.6) + np.power((Ca_blk),1.6))
    I_PMCA_iz = f_PMCA_iz * Amp_PMCA * (np.power((Ca_iz),1.6)) / (np.power((K_m),1.6) + np.power((Ca_iz),1.6))
    I_PMCA = I_PMCA_iz + I_PMCA_blk

# Na+/Ca2+ exchange current
    f_NCX_iz = 0.1
    f_NCX_blk = 0.9
    AMP_NCX = 61.06   
    k_1 = np.exp((0.32 * F * Vm) / (R * T))
    k_2 = np.exp(((0.32 - 1) * F * Vm) / (R * T))
    k_3 = 1
    k_4 = 1
    a_1_on = 0.002
    a_1_off = 0.0015
    B_1_on = 0.0012
    B_1_off = 0.0000005
    a_2_on = 0.00006
    a_2_off = 0.02
    B_2_on = 0.18
    B_2_off = 0.0002
    K_m_Nao = 87.5
    K_m_Nai = 20.74854
    K_m_Cao = 1.38
    K_m_Cai = 0.0184
    K_m_act = 0.004 
    
    p_E_2_NCX_blk = 1 - p_E_1_NCX_blk - p_I_1_NCX_blk - p_I_2_NCX_blk
    p_E_2_NCX_iz = 1 - p_E_1_NCX_iz - p_I_1_NCX_iz - p_I_2_NCX_iz
    q_blk_E_1Na = 1 / (1 + np.power((K_m_Nai / Na_i),3) * (1 + (Ca_blk / K_m_Cai)))
    q_iz_E_1Na = 1 / (1 + np.power((K_m_Nai / Na_i),3) * (1 + (Ca_iz / K_m_Cai)))
    q_blk_E_1Ca = 1 / (1 + (K_m_Cai / Ca_blk) * (1 + np.power((Na_i / K_m_Nai),3)))
    q_iz_E_1Ca = 1 / (1 + (K_m_Cai / Ca_iz) * (1 + np.power((Na_i / K_m_Nai),3)))
    q_E_2Na = 1 / (1 + np.power((K_m_Nao / Na_o),3) * (1 + (Ca_o / K_m_Cao)))
    q_E_2Ca = 1 / (1 + (K_m_Cao / Ca_o) * (1 + np.power((Na_o / K_m_Nao),3)))
    f_Caina_blk = Ca_blk / (Ca_blk + K_m_act) 
    f_Caina_iz = Ca_iz / (Ca_iz + K_m_act)
    a_1_blk = q_blk_E_1Na * (f_Caina_blk * a_1_on + (1 - f_Caina_blk) * a_1_off)
    a_1_iz = q_iz_E_1Na * (f_Caina_iz * a_1_on + (1 - f_Caina_iz) * a_1_off)
    B_1_blk = f_Caina_blk * B_1_on + (1 - f_Caina_blk) * B_1_off
    B_1_iz = f_Caina_iz * B_1_on + (1 - f_Caina_iz) * B_1_off
    a_2_blk = f_Caina_blk * a_2_on + (1 - f_Caina_blk) * a_2_off
    a_2_iz = f_Caina_iz * a_2_on + (1 - f_Caina_blk) * a_2_off
    B_2_blk = f_Caina_blk * B_2_on + (1 - f_Caina_blk) * B_2_off
    B_2_iz = f_Caina_iz * B_2_on + (1 - f_Caina_iz) * B_2_off
    a_E = k_2 * q_E_2Na + k_4 * q_E_2Ca
    B_E_blk = k_1 * q_blk_E_1Na + k_3 * q_blk_E_1Ca
    B_E_iz = k_1 * q_iz_E_1Na + k_3 * q_iz_E_1Ca

    dp_E_1_NCX_blk = p_E_2_NCX_blk * a_E + p_I_1_NCX_blk * B_1_blk + p_I_2_NCX_blk * B_2_blk + p_E_1_NCX_blk * (B_E_blk + a_1_blk + a_2_blk)
    dp_E_1_NCX_iz = p_E_2_NCX_iz * a_E + p_I_1_NCX_iz * B_1_iz + p_I_2_NCX_iz * B_2_iz + p_E_1_NCX_iz * (B_E_iz + a_1_iz + a_2_iz)
    dp_I_1_NCX_blk = p_E_1_NCX_blk * a_1_blk - p_I_1_NCX_blk * B_1_blk
    dp_I_1_NCX_iz = p_E_1_NCX_iz * a_1_iz - p_I_1_NCX_iz * B_1_iz
    dp_I_2_NCX_blk = p_E_1_NCX_blk * a_2_blk - p_I_2_NCX_blk * B_2_blk
    dp_I_2_NCX_iz = p_E_1_NCX_iz * a_2_iz - p_I_2_NCX_iz * B_2_iz
 
    v_cyc_NCX_blk = k_1 * q_blk_E_1Na * p_E_1_NCX_blk - k_2 * q_E_2Na * p_E_2_NCX_blk
    v_cyc_NCX_iz = k_1 * q_iz_E_1Na * p_E_1_NCX_iz - k_2 * q_E_2Na * p_E_2_NCX_iz
    I_NCX_blk = f_NCX_blk * AMP_NCX * v_cyc_NCX_blk
    I_NCX_iz = f_NCX_iz * AMP_NCX * v_cyc_NCX_iz
    I_NCX_Na_blk = 3 * I_NCX_blk
    I_NCX_Na_iz = 3 * I_NCX_iz
    I_NCX_Ca_blk = -2 * I_NCX_blk
    I_NCX_Ca_iz = -2 * I_NCX_iz

    I_tot_Ca_iz = I_CaL_Ca_iz + I_PMCA_iz + I_NCX_Ca_iz + I_Cab_iz

    i_net_Ca = i_CaL_Ca + i_Ca_T + i_Cab - (2 * i_NaCa)
    dCa_2_tot_jnc = ((I_tot_Ca_jnc * Cm) / (Vjnc * 2 * F)) + (J_Ca_rel / Vjnc) - (J_Ca_jnciz / Vjnc)    
    dCa_2_tot_iz = ((I_tot_Ca_iz * Cm) / (Viz * 2 * F)) + (J_Ca_jnciz / Viz) - (J_Ca_izblk / Viz)

    I_tot_Ca_blk = I_CaL_Ca_blk + I_PMCA_blk + I_NCX_Ca_blk + I_Cab_blk
    # from huvec2
    chi_r_infinity = 1 / (1 + np.exp(-(Vm + 8.337) / 6.789))
    tau_chi_r_fast = 12.98 + 1 / (0.3652 * np.exp((Vm - 31.66) / 3.869) + 0.00004123 * np.exp(-(Vm - 47.78) / 20.38))
    tau_chi_r_slow = 1.865 + 1 / (0.06629 * np.exp((Vm - 34.70) / 7.355) + 0.00001128 * np.exp(-(Vm - 19.74) / 25.94))
    A_chi_r_fast = 1 / (1 + np.exp((Vm + 4.81) / 38.21))
    A_chi_r_slow = 1 - A_chi_r_fast
    dchi_r_fast = (chi_r_infinity - chi_r_fast) / tau_chi_r_fast
    dchi_r_slow = (chi_r_infinity - chi_r_slow) / tau_chi_r_slow
    chi_r = A_chi_r_fast * chi_r_fast + A_chi_r_slow * chi_r_slow
    R_Kr = 1 / ((1 + np.exp((Vm + 55) / 75)) * (1 + np.exp((Vm - 10) / 30)))
    p_O_Kr = chi_r * R_Kr
    chi_Kr = np.sqrt(K_o / 4.5)
    I_Kr_K = G_Kr * chi_Kr * (Vm - E_K) * p_O_Kr

    MgATP = 5e3
    alpha_1 = 25900 * MgATP
    alpha_2 = 2540 / (1 + np.power((K_dCai / Ca_2_blk),1.7))
    alpha_3 = 5.35 / (1 + np.power((Ca_2_SRup / K_dCasr),1.7))
    beta_1 = 0.1972 / (1 + np.power((Ca_2_blk / K_dCai),1.7))
    beta_2 = (25435 * MgADP) / (1 + np.power((K_dCasr / Ca_2_SRup),1.7))
    beta_3 = 149 * Pi
    v_cyc = (6.86 * (alpha_1 * alpha_2 * alpha_3 - beta_1 * beta_2 * beta_3)) / (alpha_2 * alpha_3 + beta_1 * alpha_3 + beta_1 * beta_2 + alpha_1 * alpha_3 + beta_2 * alpha_1 + beta_2 * beta_3 + alpha_1 * alpha_2 + beta_3 * beta_1 + beta_3 * alpha_2)
    J_SERCA = (Amp_SERCA * v_cyc / (2 * F)) * Sc_Cell

    # Release site of the SR (SRrl)
    K_d_CSQN_Ca = k_off_CSQN / k_on_CSQN
    a = 1
    b = B_Total_CSQN - Ca_2_tot_SRrl + K_d_CSQN_Ca
    c = K_d_CSQN_Ca * Ca_2_tot_SRrl
    Ca_SRrl = (-b + np.sqrt(np.power(b,2) - 4 * a * c)) / (2 * a)

    # Ca2+ transfer from SR uptake site to release site
    P_trans = 22.745
    J_trans_SR = P_trans * (Ca_2_SRup - Ca_SRrl)
    J_CaTnC_3Ca = (3 * (TSCa_3 + TSCa_3W + TSCa_3S)) / 1000


    dCa_2_tot_blk = ((I_tot_Ca_blk * Cm) / (Vblk * 2 * F)) + (J_SERCA / Vblk) - (J_Ca_izblk / Vblk) - J_CaTnC_3Ca
    dCa_2_SRup = (J_SERCA / V_SRup) - (J_trans_SR / V_SRup)

    # LCC open; RyR open
    g_D = 0.065
    J_R = 0.02
    J_L = 0.000913
    Z_Ca = 1
    pV = (Z_Ca * F * Vm) / (R * T)
    Ca_nd = (Ca_jnc + (J_L / g_D) * Ca_SRrl + (J_L / g_D) * (pV * np.exp(-pV)) / (1 - np.exp(-pV)) * Ca_o) / (1 + (J_L / g_D) + (J_L / g_D) * (pV / (1 - np.exp(-pV))))
    dCa_blk = Ca_nd

    # RyR channel and the Ca2+ leak from SR
    K_RyR = 20.5
    t_R = 0.1
    up_R = 0.5
    T_R = 2.43
    seta_R = 0.00048
    c = 0.01
    d = 100
    p_O = Y_ooo + Y_coo + Y_cco + Y_oco
    p_O_basal = 0.0000504
    p_RyR = 3532
    B0 = np.power((Ca_nd),2) / (t_R * np.power((Ca_nd),2) + np.power(K_RyR,2))
    B1 = up_R / t_R
    u0 = (np.power((Ca_nd),2) + c * np.power((K_RyR),2)) / (T_R * np.power((Ca_nd),2) + np.power((K_RyR),2))
    u1 = seta_R * d * (np.power((Ca_nd),2) + c * np.power((K_RyR),2)) / (T_R * (d * np.power((Ca_nd),2) + c * np.power((K_RyR),2))) * 1 / (1 + np.power((0.8 / Ca_SRrl),2))
    J_rel_SR = p_RyR * (p_O + p_O_basal) * (Ca_SRrl - Ca_jnc)

    dCa_2_tot_SRrl = (J_trans_SR / V_SRrl) - (J_rel_SR / V_SRrl)

    # slow_time_dependent_potassium_current
    P_Ks_K = 0.6955
    P_Na_IKs = 0.04 * P_Ks_K
    f_Ks_iz = 0.1
    f_Ks_blk = 0.9
    p_O_Ks_blk = np.power((p_O_v),2) * (0.99 * p_blk_O_c + 0.01)
    p_O_Ks_iz = np.power((p_O_v),2) * (0.99 * p_iz_O_c + 0.01)
    I_Ks_K_blk = f_Ks_blk * P_Ks_K * GHK_K * p_O_Ks_blk
    I_Ks_Na_blk = f_Ks_blk * P_Ks_Na * GHK_Na * p_O_Ks_blk 
    I_Ks_K_iz = f_Ks_iz * P_Ks_K * GHK_K * p_O_Ks_iz
    I_Ks_Na_iz = f_Ks_iz * P_Ks_Na * GHK_Na * p_O_Ks_iz
    #I_Ks = (I_Ks_K_iz + I_Ks_Na_iz) + (I_Ks_K_blk + I_Ks_Na_blk)

    # from Huvec2
    f_C_Na = 1 / (1 + np.exp(-(Vm + 48) / 7))
    k_C2O = 1 / (0.0025 * np.exp(Vm / -8.0) + 0.15 * np.exp(Vm / -100.0))
    k_OC = 1 / (30 * np.exp(Vm / 12.0) + 0.53 * np.exp(Vm / 50.0))
    k_OI2 = 1 / (0.0433 * np.exp(Vm / -27.0) + 0.34 * np.exp(Vm / -2000.0))
    k_C2I2 = 0.5 / (1.0 + (k_I2O * k_OC) / (k_OI2 * k_C2O))
    k_I2C = 0.5 - k_C2I2
    k_Isb = 1 / (300000.0 * np.exp(Vm / 10.0) + 50000 * np.exp(Vm / 16.0))
    k_Isf = 1 / (0.016 * np.exp(Vm / -9.9) + 8.0 * np.exp(Vm / -45.0))
    k_OI1 = k_OI2
    k_I1C = k_I2C
    k_C2I1 = k_C2I2
    p_C_NaT = 1.0 - p_I_2_NaT - p_I_s_NaT - p_O_NaT2
    dp_O_NaT2 = k_I2O * p_I_2_NaT + f_C_Na * k_C2O * p_C_NaT - (k_OC + k_OI2) * p_O_NaT2
    p_C_NaL = 1.0 - p_I_2_NaL - p_I_s_NaL - p_I_1_NaL - p_O_NaL2
    dp_O_NaL2 = k_I1O * p_I_1_NaL + f_C_Na * k_C2O * p_C_NaL - (k_OC + k_OI1) * p_O_NaL2

    I_NaT_Na = (1 - f_LSM) * P_Na * GHK_Na * p_O_NaT2
    I_NaT_K = (1 - f_LSM) * P_Na * 0.1 * GHK_K * p_O_NaT2
    I_NaL_Na = f_LSM * P_Na * GHK_Na * p_O_NaL2
    I_NaL_K = f_LSM * P_Na * 0.1 * GHK_K * p_O_NaL2

 # Na+/K+ pump current
    k_10 = 0.72
    k_11 = 0.08
    k_20 = 0.08
    k_21 = 0.008
    k_30 = 4
    k_31 = 8000
    k_40 = 0.3
    k_41 = 0.2
    Pi = 0.50872066859173026
    H = 0.0001

    K_d_Naio = 5
    K_d_Naoo = 26.8
    K_d_Kio = 18.8
    K_d_Koo = 0.8
    K_d_MgATP = 0.6
    semo_Nai = -0.14
    semo_Nao = 0.44
    semo_Ki = -0.14
    semo_Ko = 0.23
    K_d_Nao = K_d_Naoo * np.exp((semo_Nao * F * Vm) / (R * T))
    K_d_Nai = K_d_Naio * np.exp((semo_Nai * F * Vm) / (R * T))
    K_d_Ko = K_d_Koo * np.exp((semo_Ko * F * Vm) / (R * T))
    K_d_Ki = K_d_Kio * np.exp((semo_Ki * F * Vm) / (R * T))
    Na_il = Na_i / K_d_Nai
    Na_ol = Na_o / K_d_Nao
    K_il = K_i / K_d_Ki
    K_ol = K_o / K_d_Ko
    MgATPl = MgATPcyt / K_d_MgATP


 
    Amp_NaK = 24.35
    Stoi_Nak_Na = 3
    Stoi_Nak_K = -2
    a_10 = (k_10 * np.power((Na_il),3)) / (np.power((1 + Na_il),3) + np.power((1 + K_il),2) - 1)
    a_20 = k_20 
    a_30 = (k_30 * np.power((K_ol),2)) / (np.power((1 + Na_ol),3) + np.power((1 + K_ol),2) - 1)
    a_40 = (k_40 * MgATPl) / (1 + MgATPl)

    a_11 = k_11 * MgADP
    a_21 = (k_21 * np.power((Na_ol),3)) / (np.power((1 + Na_ol),3) + np.power((1 + K_ol),2) - 1)
    a_31 = (k_31 * Pi * H) / (1 + MgATPl)
    a_41 = (k_41 * np.power((K_il),2)) / (np.power((1 + Na_il),3) + np.power((1 + K_il),2) - 1)


    P_14_15 = 1 - P_1_6 - P_7 - P_8_13

    dP_1_6 = -a_10 * P_1_6 + a_11 * P_7 + a_40 * P_14_15 - a_41 * P_1_6
    dP_7 = -a_20 * P_7 + a_21 * P_8_13 + a_10 * P_1_6 - a_11 * P_7
    dP_8_13 = -a_30 * P_8_13 + a_31 * P_14_15 + a_20 * P_7 - a_21 * P_8_13

    v_cyc_NaK = a_20 * P_7 - a_21 * P_8_13
    I_NaK = Amp_NaK * v_cyc_NaK
    I_NaK_Na = Stoi_Nak_Na * I_NaK
    I_NaK_K = Stoi_Nak_K * I_NaK

    # Background non-selective cation current
    P_bNSC_K = 0.00014
    P_bNSC_Na = 0.00035
    I_bNSC_K = P_bNSC_K * GHK_K
    I_bNSC_Na = P_bNSC_Na * GHK_Na

    I_bNSC = I_bNSC_K + I_bNSC_Na 

    # from huvec2
    P_l_Ca_K = P_l_Ca_Na
    p_O_blk = 1.0 / (1.0 + np.power((0.0012 / Ca_2_blk),3))
    p_O_iz = 1.0 / (1.0 + np.power((0.0012 / Ca_2_iz),3))
    I_LCCa_Na_blk = P_l_Ca_Na * f_l_Ca_blk * GHK_Na * p_O_blk
    I_LCCa_Na_iz = P_l_Ca_Na * f_l_Ca_iz * GHK_Na * p_O_iz
    I_LCCa_K_blk = P_l_Ca_K * f_l_Ca_blk * GHK_K * p_O_blk
    I_LCCa_K_iz = P_l_Ca_K * f_l_Ca_iz * GHK_K * p_O_iz
    #I_LCCa = I_LCCa_Na_iz + I_LCCa_K_iz + I_LCCa_Na_blk + I_LCCa_K_blk

    I_tot_Na = (I_CaL_Na_jnc + I_CaL_Na_iz + I_CaL_Na_blk) + (I_NCX_Na_iz + I_NCX_Na_blk) + (I_Ks_Na_iz + I_Ks_Na_blk) + I_NaT_Na + I_NaL_Na + I_NaK_Na + I_Kto_Na + I_bNSC_Na + (I_LCCa_Na_iz + I_LCCa_Na_blk)
    dNa_i = (I_tot_Na * Cm) / (V_cyt * F)

    # SR_calcium_pump
    k1 = 0.01
    k3 = 1
    k4 = 0.01
    Km_CaSR = 0.08
    Km_CaCyto = 0.0008
    Km_ATP = 0.1
    i_max = 162500
    p_E1Ca = 1 / (1 + (Km_CaSR / Caup))
    p_E2Ca = 1 / (1 + (Km_CaCyto / Ca_i))
    p_E1 = 1 - p_E1Ca
    p_E2 = 1 - p_E2Ca
    k2 = 1 / (1 + (Km_ATP / ATPi))
    # SR_calcium_pump_y_gate
    alpha_y = k2 * p_E2Ca + k4 * p_E2
    beta_y = k1 * p_E1Ca + k3 * p_E1
    dy8 = alpha_y * (1 - y8) - (beta_y * y8)
    
    i_SR_U = i_max * (k1 * p_E1Ca * y8 - (k2 * p_E2Ca * (1 - y8)))



    # Ca_concentrations_in_SR
    CSQN_max = 10
    K_mCSQN = 0.8
    V_rel = 160
    V_up = 400
    b1 = (CSQN_max - Ca_Total3) + K_mCSQN 
    c1 = K_mCSQN * Ca_Total3
    Carel = (np.sqrt(np.power(b1,2) + 4 * c1) - b1) / 2
    # SR_T_current
    P_SR_T = 386
    i_SR_T = P_SR_T * (Caup - Carel)

    # SR_L_current
    P_SR_L = 459
    i_SR_L = P_SR_L * (Caup - Ca_i)

    # L_type_Ca_channel_Ca_dependent_gate
    k_U_UCa = 6.954
    k_C_CCa = 6.954
    k_C_U = 0.143
    k_U_C = 0.35
    k_CCa_UCa = 0.0003
    k_UCa_CCa = 0.35
    k_CCa_C = 0.0042
    I_CaL = 0.0676 * CF_Ca
    CaDiadic = I_CaL * p_open_CaL 
    Cacm = Ca_i - (0.3 * I_CaL)
    CaEffC = Cacm * p_AP_CaL
    CaEffU = CaEffC + Ca_i * (1 - p_AP_CaL)
    k_UUCa_Ca = k_U_UCa * CaEffU
    k_CCCa_Ca = k_C_CCa * CaEffC
    p_CCa = 1 - p_C - p_U - p_UCa
    k_UCa_U = (k_CCa_C * k_C_U * k_U_UCa * k_UCa_CCa) / (k_U_C * k_C_CCa * k_CCa_UCa)
    dp_U = p_C * k_C_U + p_UCa * k_UCa_U - (p_U * (k_UUCa_Ca + k_U_C))
    dp_UCa = p_U * k_UUCa_Ca + p_CCa * k_CCa_UCa - (p_UCa * (k_UCa_CCa + k_UCa_U))
    dp_C = p_CCa * k_CCa_C + p_U * k_U_C - (p_C * (k_C_U + k_C_CCa * Cacm * p_AP_CaL))

    # RyR_channel
    P_RyR = 62000
    k4 = 0.000849
    Diadid_Factor = -150
    i_RyR = P_RyR * (Carel - Ca_i) * p_open_RyR
    dp_open_RyR = p_close_RyR * k1 - (p_open_RyR * k2)
    dp_close_RyR = k3 * (1 - (p_open_RyR + p_close_RyR))-((k1 + k4) * p_close_RyR)
    k1 = 280000 * np.power((Ca_i),2) + Diadid_Factor * CaDiadic
    k2 = 0.08 / (1 + (0.36 / Carel))
    k3 = 0.000377 * np.power((Carel),2)

    dCa_Total3 = (i_SR_T - i_RyR) / (2 * F * V_rel)
    dCaup = (-i_SR_U - i_SR_T - i_SR_L) / (2 * F * V_up)

    # NL_model
    T_t = 0.07
    Y_1 = 39
    Y_2 = 0.0039
    Y_3 = 0.03
    Y_4 = 0.12
    Y_d = 0.027
    Z_1 = 0.03
    Z_2 = 0.0039
    Z_3 = 1560
    L_a = 1.17
    L = 0.9623799975411884
    KForceEC = 140000
    ZeroForceEL = 0.97
    KForceLinearEc = 200
    ForceFactor = 1800000
    B = 1.2
    h_c = 0.005
    h = L - X 
    p = 1 - pCa - pCaCB - pCB
    Q_b = Y_1 * Ca_i * p - (Z_1 * pCa)
    EffFraction = np.exp(-20 * np.power((L - L_a),2))
    Q_a = Y_2 * pCa * EffFraction - (Z_2 * pCaCB) 
    Q_r = Y_3 * pCaCB - (Z_3 * pCB * Ca_i)
    Q_d = Y_4 * pCB
    dX = B * (h - h_c)    
    Q_d1 = Y_d * np.power(dX,2) * pCB
    Q_d2 = Y_d * np.power(dX,2) * pCaCB
    dpCa = Q_b - Q_a
    dpCaCB = Q_a - Q_r - Q_d2
    dpCB = Q_r - Q_d - Q_d1
    dCaidt = T_t * (Q_d2 + Q_r - Q_b)
    dATPdt = (-0.4) * pCaCB * T_t 
    CBBound = T_t * (pCaCB + pCB)
    NewCBF = ForceFactor * CBBound
    ForceEcomp = KForceEC * np.power((ZeroForceEL - L),5) + KForceLinearEc * (ZeroForceEL - L)
    ForceCB = NewCBF * h
    ForceExt=-(ForceEcomp) + ForceCB


    dCa_Total2 = -((i_net_Ca - i_SR_U - i_RyR - i_SR_L) / (2 * F * Vi)) + dCaidt


    dCaMCa3 = (k_on_CaM * Ca_blk * (B_Total_CaM - CaMCa3)) - (k_off_CaM * CaMCa3)
    dTnChCa3 = (k_on_TnCh * Ca_blk * (B_Total_TnCh - TnChCa3)) - (k_off_TnCh * TnChCa3)
    dSRCa3 = (k_on_SR * Ca_blk * (B_Total_SR - SRCa3)) - (k_off_SR * SRCa3)

    # Intermediate zone (iz)
    #L_free_iz = B_Total_L_iz / (1 + (Ca_iz / K_dL_iz))
    #K_dL_iz = k_off_L_iz / k_on_L_iz
    #H_free_iz = B_Total_H_iz / (1 + (Ca_iz / K_dH_iz))
    #K_dH_iz = k_off_H_iz / k_on_H_iz
    #dCa_iz = Ca_tot_iz / (1 + (Lf_iz / K_dL_iz) + (Hf_iz / K_dH_iz))

    Lf_iz = B_Total_L_iz / (1 + (Ca_iz / K_dL_iz))
    K_dL_iz = k_off_L_iz / k_on_L_iz
    Hf_iz = B_Total_H_iz / (1 + (Ca_iz / K_dH_iz))
    K_dH_iz = k_off_H_iz / k_on_H_iz
    dCa_iz = Ca_2_tot_iz / (1 + (Lf_iz / K_dL_iz) + (Hf_iz / K_dH_iz))

    # Junctional space (jnc)
    #L_free_jnc = B_Total_L_jnc / (1 + (Ca_jnc / K_dL_jnc))
    #K_dL_jnc = k_off_L_jnc / k_on_L_jnc
    #H_free_jnc = B_Total_H_jnc / (1 + (Ca_jnc / K_dH_jnc))
    #K_dH_jnc = k_off_H_jnc / k_on_H_jnc
    #dCa_jnc = Ca_tot_jnc / (1 + (Lf_jnc / K_dL_jnc) + (Hf_jnc / K_dH_jnc))

    Lf_jnc = B_Total_L_jnc / (1 + (Ca_jnc / K_dL_jnc))
    K_dL_jnc = k_off_L_jnc / k_on_L_jnc
    Hf_jnc = B_Total_H_jnc / (1 + (Ca_jnc / K_dH_jnc))
    K_dH_jnc = k_off_H_jnc / k_on_H_jnc
    dCa_jnc = Ca_2_tot_jnc / (1 + (Lf_jnc / K_dL_jnc) + (Hf_jnc / K_dH_jnc))

    #v_cyc_NCX_blk = k_1 * q_blk_E_1Na * p_E_1_NCX_blk - k_2 * q_blk(E_2Na) * p_E_2_NCX_blk
    #v_cyc_NCX_iz = k_1 * q_iz_E_1Na * p_E_1_NCX_iz - k_2 * q_iz(E_2Na) * p_E_2_NCX_iz


    I_NCX = I_NCX_iz + I_NCX_blk

    
    # from huvec2

    p_O_KATP = 0.8 / (1.0 + np.power((ATP_cyt / 0.1),2))
    chi_KATP = 0.0236 * np.power((K_o),0.24)
    I_KATP_K = G_KATP * (Vm - E_K) * p_O_KATP * chi_KATP

    # from huvec2
    alpha_Mg = 12.0 * np.exp(-0.025 * (Vm - E_K))
    beta_Mg = 28 * Mg_2_cyt * np.exp(0.025 * (Vm - E_K))
    f_O = alpha_Mg / (alpha_Mg + beta_Mg)
    f_B = beta_Mg / (alpha_Mg + beta_Mg)
    po_Mg = f_O * f_O * f_O
    po_Mg1 = 3.0 * f_O * f_O * f_B
    po_Mg2 = 3.0 * f_O * f_B * f_B
    alpha_SPM = 0.17 * np.exp(-0.07 * (Vm - E_K + 8 * Mg_2_cyt)) / (1.0 + 0.01 * np.exp(0.12 * (Vm - E_K + 8 * Mg_2_cyt)))
    beta_SPM = 0.28 * SPM * np.exp(0.15 * (Vm - E_K + 8 * Mg_2_cyt)) / (1.0 + 0.01 * np.exp(0.13 * (Vm - E_K + 8 * Mg_2_cyt)))

    dPb_spm2 = beta_SPM * po_Mg * (1 - Pb_spm2) - alpha_SPM * Pb_spm2
    po_mode1 = f_mode1 * (1 - Pb_spm2) * (po_Mg + (2.0 / 3.0) * po_Mg1 + (1.0 / 3.0) * po_Mg2)
    po_mode2 = (1 - f_mode1) / (1.0 + SPM / (40.0 * np.exp(-(Vm - E_K) / 9.1)))
    p_O_K1 = po_mode1 + po_mode2
    chi_K1 = np.power((K_o / 4.5),0.4) / (1.0 + np.exp(-(K_o - 2.2) / 0.6))
    I_K1_K = G_K1 * chi_K1 * (Vm - E_K) * p_O_K1

    # sodium_potassium_pump_y_gate
    k2 = 0.04
    k3 = 0.01
    k4 = 0.165
    Km_Nai = 4.05
    Km_Nao = 69.8
    Km_Ki = 32.88
    Km_Ko = 0.258
    Km_ATP = 0.094
    P_NaK = 21
    p_E1Na = 1 / (1 + (np.power((Km_Nai / Na_i),1.06)) * (1 + np.power((K_i / Km_Ki),1.12)))
    Nao_Eff = Na_o * np.exp(-(0.82 * F * Vm) / (R * T))
    p_E2Na = 1 / (1 + (np.power((Km_Nao / Nao_Eff),1.06)) * np.power((1 + (K_o / Km_Ko)),1.12))
    p_E1K = 1 / (1 + (np.power((Km_Ki / K_i),1.12)) * (1 + np.power((Na_i / Km_Nai),1.06)))
    p_E2K = 1 / (1 + (np.power((Km_Ko / K_o),1.12)) * (1 + np.power((Nao_Eff / Km_Nao),1.06)))
    k1 = 0.37 / (1 + np.exp(Km_ATP / ATPi))

    alpha_y = k2 * p_E2Na + k4 * p_E2K
    beta_y = k1 * p_E1Na + k3 * p_E1K
    dy7 = alpha_y * (1 - y7) - (beta_y * y7)
    
    # sodium_potassium_pump

    i_NaK = P_NaK * Cm * (k1 * p_E1Na * y7 - (k2 * p_E2Na * (1 - y7)))

    # ATP_production
    ProducingRate_Max = 0.003
    Adenosine_Total = 5
    dATPi =  (ProducingRate_Max * (Adenosine_Total - ATPi) + dATPdt - (i_NaK / (F * Vi))) + (i_SR_U / (4 * F * Vi))

    # sodium_current
    #i_Na_Na = P_Na * CF_Na * SCVDG.p_AP_Na * SCUSG.y2
    #i_Na_K = 0.1 * P_Na * CF_K * SCVDG.p_AP_Na * SCUSG.y2
    f_L = 0.1575
    P_Na = 6.756
    #I_Na = I_NaT + I_NaL

    # Transient component
    #I_NaT = (1 - f_L) * P_Na * (GHK_Na + 0.1 * GHK_K) * p_O_NaT3
    k_I2O = 0.0001312
    dp_C_NaT = k_OC * p_O_NaT3 + k_I2C * p_I_2_NaT + k_Isb * p_I_s_NaT - (k_Isf + f_C_Na * (k_C2O + k_C2I2)) * p_C_NaT
    dp_O_NaT3 = k_I2O * p_I_2_NaT + f_C_Na * k_C2O * p_C_NaT - (k_OC + k_I2O) * p_O_NaT3
    dp_I_2_NaT = f_C_Na * k_C2I2 * p_C_NaT + k_OI2 * p_O_NaT3 + k_Isb * p_I_s_NaT - (k_I2C + k_I2O + k_Isf) * p_I_2_NaT
    dp_I_s_NaT = k_Isf * p_I_2_NaT + k_Isf * p_C_NaT - 2 * k_Isb * p_I_s_NaT
    f_C_Na = 1 / (1 + np.exp(-(Vm + 48) / 7))
    k_C2O = 1 / (0.0025 * np.exp(Vm / -8) + 0.15 * np.exp(Vm / -100))
    k_OC = 1 / (30 * np.exp(Vm / 12) + 0.53 * np.exp(Vm / 50))
    k_OI2 = 1 / (0.0433 * np.exp(Vm / -27) + 0.34 * np.exp(Vm / -2000))
    k_C2I2 = 0.5 / (1 + ((k_I2O * k_OC) / (k_OI2 * k_C2O)))
    k_I2C = 0.5 - k_C2I2
    k_Isb = 1 / (300000 * np.exp(Vm / 10) + 50000 * np.exp(Vm / 16))
    k_Isf = 1 / (0.016 * np.exp(Vm / -9.9) + 8 * np.exp(Vm / -45))

    # Late component
    #I_NaL = f_L * P_Na * (GHK_Na + 0.1 * GHK_K) * p_O_NaL3 
    k_I1I2 = 0.00534
    k_OI1 = k_OI2
    k_I1O = 0.01
    k_I1C = k_I2C
    k_C2I1 = k_C2I2
    dp_C_NaL = k_OC * p_O_NaL3 + k_I1C * p_I_1_NaL + k_I2C * p_I_2_NaL - k_Isb * p_I_s_NaL - (k_Isf + f_C_Na * (k_C2O + k_C2I2 + k_C2I1)) * p_C_NaL
    dp_O_NaL3 = k_I1O * p_I_1_NaL + f_C_Na * k_C2O * p_C_NaL - (k_OC + k_OI1) * p_O_NaL3
    dp_I_1_NaL = k_OI1 * p_O_NaL3 + f_C_Na * k_C2I1 * p_C_NaL - (k_I1O + k_I1C + k_I1I2) * p_I_1_NaL 
    dp_I_2_NaL = f_C_Na * k_C2I2 * p_C_NaL + k_I1I2 * p_I_1_NaL + k_Isb * p_I_s_NaL - (k_I2C + k_Isf) * p_I_2_NaL
    dp_I_s_NaL = k_Isf * p_I_2_NaL + k_Isf * p_C_NaL - 2 * k_Isb * p_I_s_NaL


    # sodium_current_voltage_dependent_gate
    k_AI_AP = 0.0000875
    p_RI_Na = 1 - p_RP_Na - p_AP_Na - p_AI_Na
    k_RP_AP = 1 / ((0.1027 * np.exp(-Vm / 8)) + (0.25 * np.exp(-Vm / 50)))
    k_AP_RP = 1 / ((26 * np.exp(Vm / 17)) + (0.02 * np.exp(Vm / 800)))
    k_AP_AI = 1 / (0.8 * np.exp(-Vm / 400))
    k_AI_RI = 1 / ((1300 * np.exp(Vm / 20)) + (0.04 * np.exp(Vm / 800)))
    k_RI_AI = 1 / ((0.0001027 * np.exp(-Vm / 8)) + (5 * np.exp(-Vm / 400)))
    k_RP_RI = 0.01 / (1 + ((k_AI_AP * k_AP_RP * k_RI_AI) / (k_AP_AI * k_RP_AP * k_AI_RI)))
    k_RI_RP = 0.01 - k_RP_RI
    dp_RP_Na = p_AP_Na * k_AP_RP + p_RI_Na * k_RI_RP - (p_RP_Na * (k_RP_RI + k_RP_AP))
    dp_AP_Na = p_RP_Na * k_RP_AP + p_AI_Na * k_AI_AP - (p_AP_Na * (k_AP_RP + k_AP_AI))
    dp_AI_Na = p_RI_Na * k_RI_AI + p_AP_Na * k_AP_AI - (p_AI_Na * (k_AI_RI + k_AI_AP))

    # sodium_current_ultra_slow_gate 
    alpha_y = 1 / ((9000000000 * np.exp(Vm / 5)) + (8000 * np.exp(Vm / 100)))
    beta_y = 1 / ((0.014 * np.exp(-Vm / 5)) + (4000 * np.exp(-Vm / 100)))
    dy2 = alpha_y * (1 - y2) - (beta_y * y2)

    #GHK_Ca = ((Z_Ca * F * Vm) / (R * T)) * (Ca_i - (Ca_o * np.exp((-Z_Ca * F * Vm) / R * T) / (1 - np.exp(-Z_Ca * F * Vm / R * T)))) 
    #GHK_Na = ((Z_Na * F * Vm) / (R * T)) * (Na_i - (Na_o * np.exp((-Z_Na * F * Vm) / R * T) / (1 - np.exp(-Z_Na * F * Vm / R * T)))) 
    #GHK_K = ((Z_K * F * Vm) / (R * T)) * (K_i - (K_o * np.exp((-Z_K * F * Vm) / R * T) / (1 - np.exp(-Z_K * F * Vm / R * T)))) 
    #GHK_Ca_blk = ((Z_Ca_blk * F * Vm) / (R * T)) * (Ca_blk_i - (Ca_blk_o * np.exp((-Z_Ca_blk * F * Vm) / R * T) / (1 - np.exp(-Z_Ca_blk * F * Vm / R * T)))) 
    #GHK_Ca_iz = ((Z_Ca_iz * F * Vm) / (R * T)) * (Ca_iz_i - (Ca_iz_o * np.exp((-Z_Ca_iz * F * Vm) / R * T) / (1 - np.exp(-Z_Ca_iz * F * Vm / R * T))))
    #GHK_Ca_jnc = ((Z_Ca_jnc * F * Vm) / (R * T)) * (Ca_jnc_i - (Ca_jnc_o * np.exp((-Z_Ca_jnc * F * Vm) / R * T) / (1 - np.exp(-Z_Ca_jnc * F * Vm / R * T))))
    #GHK_Na_blk = ((Z_Na_blk * F * Vm) / (R * T)) * (Na_blk_i - (Na_blk_o * np.exp((-Z_Na_blk * F * Vm) / R * T) / (1 - np.exp(-Z_Na_blk * F * Vm / R * T))))
    #GHK_Na_iz = ((Z_Na_iz * F * Vm) / (R * T)) * (Na_iz_i - (Na_iz_o * np.exp((-Z_Na_iz * F * Vm) / R * T) / (1 - np.exp(-Z_Na_iz * F * Vm / R * T))))
    #GHK_Na_jnc = ((Z_Na_jnc * F * Vm) / (R * T)) * (Na_jnc_i - (Na_jnc_o * np.exp((-Z_Na_jnc * F * Vm) / R * T) / (1 - np.exp(-Z_Na_jnc * F * Vm / R * T))))
    #GHK_K_blk = ((Z_K_blk * F * Vm) / (R * T)) * (K_blk_i - (K_blk_o * np.exp((-Z_K_blk * F * Vm) / R * T) / (1 - np.exp(-Z_K_blk * F * Vm / R * T))))
    #GHK_K_iz = ((Z_K_iz * F * Vm) / (R * T)) * (K_iz_i - (K_iz_o * np.exp((-Z_K_iz * F * Vm) / R * T) / (1 - np.exp(-Z_K_iz * F * Vm / R * T))))
    #GHK_K_jnc = ((Z_K_jnc * F * Vm) / (R * T)) * (K_jnc_i - (K_jnc_o * np.exp((-Z_K_jnc * F * Vm) / R * T) / (1 - np.exp(-Z_K_jnc * F * Vm / R * T))))

    #ATP = 6
    #I_CaL_Ca_blk = FOI.f_CaL_blk * CF.P_CaL_Ca * GHK_Ca_blk * pO_CaL_blk * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_Ca_iz = FOI.f_CaL_iz * CF.P_CaL_Ca * GHK_Ca_iz * pO_CaL_iz * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_Ca_jnc = FOI.f_CaL_jnc * CF.P_CaL_Ca * GHK_Ca_jnc * pO_CaL_jnc * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_Na_blk = FOI.f_CaL_blk * CF.P_CaL_Na * GHK_Na_blk * pO_CaL_blk * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_Na_iz = FOI.f_CaL_iz * CF.P_CaL_Na * GHK_Na_iz * pO_CaL_iz * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_Na_jnc = FOI.f_CaL_jnc * CF.P_CaL_Na * GHK_Na_jnc * pO_CaL_jnc * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_K_blk = FOI.f_CaL_blk * CF.P_CaL_K * GHK_K_blk * pO_CaL_blk * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_K_iz = FOI.f_CaL_iz * CF.P_CaL_K * GHK_K_iz * pO_CaL_iz * (1 / (1 + (1.4 / ATP),3))
    #I_CaL_K_jnc = FOI.f_CaL_jnc * CF.P_CaL_K * GHK_K_jnc * pO_CaL_jnc * (1 / (1 + (1.4 / ATP),3))

    # p_O_LCC = Y_ooo + Y_ooc

    # ATPfactor = 1 / (np.power(1 + (1.4 / ATP),3))

    # GHK_Ca_LR = (2 * Vm) / (R * T / F) * (Ca_2_nd_LR - Ca_o * np.exp((-2 * Vm) / (R * T / F))) / (1 - np.exp((-2 * Vm) / (R * T / F)))
    # GHK_Ca_L0 = (2 * Vm) / (R * T / F) * (Ca_2_nd_L0 - Ca_o * np.exp((-2 * Vm) / (R * T / F))) / (1 - np.exp((-2 * Vm) / (R * T / F)))
    # GHK_Ca_iz = (2 * Vm) / (R * T / F) * (Ca_2_iz - Ca_o * np.exp((-2 * Vm) / (R * T / F))) / (1 - np.exp((-2 * Vm) / (R * T / F)))
    # GHK_Ca_blk = (2 * Vm) / (R * T / F) * (Ca_2_blk - Ca_o * np.exp((-2 * Vm) / (R * T / F))) / (1 - np.exp((-2 * Vm) / (R * T / F)))
    # GHK_Na = ((1 * F * Vm) / (R * T)) * (Na_i - (Na_o * np.exp((-1 * F * Vm) / (R * T)) / (1 - np.exp((-1 * F * Vm) / (R * T))))) 
    # GHK_K = ((1 * F * Vm) / (R * T)) * (K_i - (K_o * np.exp((-1 * F * Vm) / (R * T)) / (1 - np.exp((-1 * F * Vm) / (R * T))))) 
    # GHK_Ca = ((1 * F * Vm) / (R * T)) * (Ca_i - (Ca_o * np.exp((-1 * F * Vm) / (R * T)) / (1 - np.exp((-1 * F * Vm) / (R * T)))))  

    # I_CaL_Ca_blk = f_CaL_blk * P_CaL_Ca * GHK_Ca_blk * Y_oo_blk * ATPfactor
    # I_CaL_Ca_iz = f_CaL_iz * P_CaL_Ca * GHK_Ca_iz * Y_oo_iz * ATPfactor
    # I_CaL_Ca_LR = f_CaL_jnc * P_CaL_Ca * GHK_Ca_LR * Y_ooo * ATPfactor
    # I_CaL_Ca_L0 = f_CaL_jnc * P_CaL_Ca * GHK_Ca_L0 * Y_ooc * ATPfactor
    # I_CaL_Ca_jnc = I_CaL_Ca_LR + I_CaL_Ca_L0
    # I_CaL_Na_blk = f_CaL_blk * P_CaL_Na * GHK_Na * Y_oo_blk * ATPfactor
    # I_CaL_Na_iz = f_CaL_iz * P_CaL_Na * GHK_Na * Y_oo_iz * ATPfactor
    # I_CaL_Na_jnc = f_CaL_jnc * P_CaL_Na * GHK_Na * p_O_LCC * ATPfactor
    # I_CaL_K_blk = f_CaL_blk * P_CaL_K * GHK_K * Y_oo_blk * ATPfactor
    # I_CaL_K_iz = f_CaL_iz * P_CaL_K * GHK_K * Y_oo_iz * ATPfactor
    # I_CaL_K_jnc = f_CaL_jnc * P_CaL_K * GHK_K * p_O_LCC * ATPfactor

    # I_CaL = (I_CaL_Ca_jnc + I_CaL_Na_jnc + I_CaL_K_jnc) + (I_CaL_Ca_iz + I_CaL_Na_iz + I_CaL_K_iz) + (I_CaL_Ca_blk + I_CaL_Na_blk + I_CaL_K_blk)

    # Fraction of ICaL
    f_CaL_blk = 0.1
    f_CaL_iz = 0.15
    f_CaL_jnc = 0.75

    # Converting factors
    P_CaL_Ca = 14.21
    P_CaL_Na = 0.0000185 * P_CaL_Ca
    P_CaL_K = 0.001468 * P_CaL_Ca

    # V-gate
    a0 = 1 / (5.55 * np.exp(-(Vm + 2) / 7) + 0.35 * np.exp(-(Vm + 2) / 3500))
    a1 = 1 / (4.35 * np.exp(-(Vm + 2) / 15) + 1.35 * np.exp((Vm + 2) / 100))

    # Ca2+ gate
    tL0 = 1885
    tL1 = 1495
    K_L = 0.00022
    ipsilon0 = (Ca_nd * a0) / (tL0 * K_L)
    ipsilon1 = (a0 / tL1) + (6 * (np.exp((Vm + 48) / -4.5) + 0.0225)) / (tL1 * (np.exp((Vm + 48) / -4.5) + 1)) 

    # L_type_Ca_channel_voltage_dependent_gate
    k_AP_AI = 0.004
    k_AI_AP = 0.001
    p_RI_CaL = 1 - p_AP_CaL - p_RP_CaL - p_AI_CaL
    k_RP_AP = 1 / ((0.27 * np.exp(-Vm / 5.9)) + (1.5 * np.exp(-Vm / 65)))
    k_AP_RP = 1 / ((480 * np.exp(Vm / 7)) + (2.2 * np.exp(Vm / 65)))
    k_RI_AI = 1 / ((0.0018 * np.exp(-Vm / 7.4)) + (2 * np.exp(-Vm / 100)))
    k_AI_RI = 1 / ((2200000 * np.exp(Vm / 7.4)) + (11 * np.exp(Vm / 100)))
    k_RP_RI = 0.04 / (1 + (k_AI_AP * k_AP_RP * k_RI_AI) / (k_AP_AI * k_RP_AP * k_AI_RI))
    k_RI_RP = 0.04 - k_RP_RI
    dp_RP_CaL = p_AP_CaL * k_AP_RP + p_RI_CaL * k_RI_RP - (p_RP_CaL * (k_RP_RI + k_RP_AP))
    dp_AP_CaL = p_RP_CaL * k_RP_AP + p_AI_CaL * k_AI_AP - (p_AP_CaL * (k_AP_RP + k_AP_AI))
    dp_AI_CaL = p_RI_CaL * k_RI_AI + p_AP_CaL * k_AP_AI - (p_AI_CaL * (k_AI_RI + k_AI_AP))



    # L_type_Ca_channel_ultra_slow_gate
    alpha_y = 1 / ((250000 * np.exp(Vm / 9)) + (58 * np.exp(Vm / 65)))
    beta_y = 1 / ((1800 * np.exp(-Vm / 14)) + (66 * np.exp(-Vm / 65)))
    dy3 = alpha_y * (1 - y3) - (beta_y * y3)

    # time_independent_potassium_current
    P_K1_0 = 1.353
    G_K1 = (P_K1_0 * np.power((K_o / 4.5),0.4)) / (1 + np.exp(-(K_o - 2.2) / 0.6))
    I_K1_K_cyt = G_K1 * (Vm - E_K) * p_O_K1
    p_O_K1 = po_mode1 + po_mode2
    fmode1 = 0.9
    po_mode1 = fmode1 * (1 - Pb_spm3) * (po_Mg + (2/3) * po_Mg1 + (1/3) * po_Mg2)
    SPM = 5
    po_mode2 = (1 - fmode1) / (1 + (SPM / (40 * np.exp(-(Vm - E_K) / 9.1))))
    mu = (0.75 * np.exp(0.035 * (Vm - E_K - 10))) / (1 + np.exp(0.015 * (Vm - E_K - 140)))
    lambda2 = (3 * np.exp(-0.048 * (Vm - E_K - 10)) * (1 + np.exp(0.064 * (Vm - E_K - 38)))) / (1 + np.exp(0.03 * (Vm - E_K - 70)))
    fB = mu / (mu + lambda2)
    fO = lambda2 / (mu + lambda2)
    fO2 = 2 * np.power(fO,2) * np.power(fB,2)
    fO3 = (8/3) * np.power(fO,3)* fB
    fO4 = np.power(fO,4)


    #I_K1 = G_K1 * (Vm - E_K) * (fO4 + fO3 + fO2) * TIPCY.y

    # The Mg2+-block in the mode1
    Mg_cyt = 0.8
    po_Mg = fO * fO * fO
    po_Mg1 = 3 * fO * fO * fB
    po_Mg2 = 3 * fO * fB * fB
    a_Mg = 12 * np.exp(-0.025 * (Vm - E_K))
    B_Mg = 28 * Mg_cyt * np.exp(0.025 * (Vm - E_K))
    fO = a_Mg / (a_Mg + B_Mg)
    fB = B_Mg / (a_Mg + B_Mg)


    # The SPM-block in the mode1
    a_SPM = (0.17 * np.exp(-0.07 * ((Vm - E_K) + 8 * Mg_cyt))) / (1 + 0.01 * np.exp(0.12 * ((Vm - E_K) + 8 * Mg_cyt)))
    B_SPM = (0.28 * SPM * np.exp(0.15 * ((Vm - E_K) + 8 * Mg_cyt))) / (1 + 0.01 * np.exp(0.13 * ((Vm - E_K) + 8 * Mg_cyt)))
    dPb_spm3 = B_SPM * po_Mg * (1 - Pb_spm3) - a_SPM * Pb_spm3

    # Pb_spm = (B_SPM * po_Mg) / (a_SPM + (B_SPM * po_Mg))

    # time_independent_potassium_current_y_gate
    #[TIPCY]

    #alpha_y = 1 / ((8000 * np.exp((Vm - E_K - 97) / 8.5)) + (7 * np.exp((Vm - E_K - 97) / 300)))
    #beta_y =  fO^4 / ((0.00014 * np.exp(-(Vm - E_K - 97) / 9.1)) + (0.2 * np.exp(-(Vm - E_K - 97) / 500)))
    #dy = alpha_y * (1 - y) - (beta_y * y)

    # rapid_time_dependent_potassium_current
    g_Kr= 0.01
    #I_Kr = g_Kr * (Vm - E_K) * TGP.p_O_Kr * (K_o / 4.5)^0.2 
    #I_Kr = g_Kr * (Vm - E_K) * p_O_Kr * (K_o / 4.5)^0.2 

    # three gating parameters, y1, y2 and y3
    a_y1 = 1 / (35 * np.exp(-Vm / 10.5) + 75 * np.exp(-Vm / 100))
    B_y1 = 1 / (470 * np.exp(Vm / 8.3) + 220 * np.exp(Vm / 29))
    a_y2 = 1 / (350 * np.exp(-Vm / 10.5) + 300 * np.exp(-Vm / 100))
    B_y2 = 1 / (1850 * np.exp(Vm / 8.3) + 2200 * np.exp(Vm / 29))
    a_y3 = 1 / (0.015 * np.exp(Vm / 6) + 7 * np.exp(Vm / 60))
    B_y3 = 1 / (0.114 * np.exp(-Vm / 9.2) + 2.3 * np.exp(-Vm / 1000))
    dy10 = a_y2 * (1 - y10) - B_y2 * y10
    dy11 = a_y3 * (1 - y11) - B_y3 * y11
    p_O_Kr = (0.6 * y9 + 0.4 + y10) * y11
    dy9 = a_y1 * (1 - y9) - B_y1 * y9


    # rapid_time_dependent_potassium_current_y1_gate
    #[RTDPC1]

    #alpha_y1 = 1 / ((20 * np.exp(-Vm / 11.5)) + (5 * np.exp(-Vm / 300)))
    #beta_y1 = (1 / ((160 * np.exp(Vm / 28)) + (200 * np.exp(Vm / 1000)))) + (1 / (2500 * np.exp(Vm / 20)))
    #dy1 = alpha_y1 * (1 - y1) - (beta_y1 * y1)

    # rapid_time_dependent_potassium_current_y2_gate
    #[RTDPC2]

    #alpha_y2 = 1 / ((200 * np.exp(-Vm / 13)) + (20 * np.exp(-Vm / 300)))
    #beta_y2 = (1 / ((1600 * np.exp(Vm / 28)) + (2000 * np.exp(Vm / 1000)))) + (1 / (10000 * np.exp(Vm / 20)))
    #dy2 = alpha_y2 * (1 - y2) - (beta_y2 * y2)

    # rapid_time_dependent_potassium_current_y3_gate
    #[RTDPC3]

    #alpha_y3 = 1 / ((10 * np.exp(Vm / 17)) + (2.5 * np.exp(Vm / 300)))
    #beta_y3 = 1 / ((0.35 * np.exp(-Vm / 17)) + (2 * np.exp(-Vm / 150)))
    #dy3 = alpha_y3 * (1 - y3) - (beta_y3 * y3)

 

    # from huvec2
    #[STDPC]

    #para_Xs1_infinity = 1 / (1 + np.exp(-(Vm + 11.60) / 8.932))
    #tau_Xs1 = 817.3 + 1 / (0.0002326 * np.exp((Vm + 48.28) / 17.80) + 0.001292 * np.exp(-(Vm + 210.0) / 230.0))
    #dpara_Xs1 = (para_Xs1_infinity - para_Xs1) / tau_Xs1
    #para_Xs2_infinity = para_Xs1_infinity
    #tau_Xs2 = 1 / (0.01 * np.exp((Vm - 50) / 20) + 0.0193 * np.exp(-(Vm + 66.54) / 31))
    #dpara_Xs2 = (para_Xs2_infinity - para_Xs2) / tau_Xs2
    #para_RKs_blk = 1 + 0.6 / (1 + (0.000038 / Ca_2_blk)^1.4)
    #para_RKs_iz = 1 + 0.6 / (1 + (0.000038 / Ca_2_iz)^1.4)
    #p_O_Ks_blk = para_Xs1 * para_Xs2 * para_RKs_blk
    #p_O_Ks_iz = para_Xs1 * para_Xs2 * para_RKs_iz
    #I_Ks_K_blk = f_Ks_blk * P_Ks_K * GHK_K * p_O_Ks_blk
    #I_Ks_K_iz = f_Ks_iz * P_Ks_K * GHK_K * p_O_Ks_iz
    #I_Ks_Na_blk = f_Ks_blk * P_Ks_Na * GHK_Na * p_O_Ks_blk
    #I_Ks_Na_iz = f_Ks_iz * P_Ks_Na * GHK_Na * p_O_Ks_iz
    #I_Ks = I_Ks_K_blk + I_Ks_K_iz + I_Ks_Na_blk + I_Ks_Na_iz

    # Voltage gate
    a_1 = 1 / (150 * np.exp(-Vm / 25) + 900 * np.exp(-Vm / 200))
    B_1 = 1 / (1000 * np.exp(Vm / 13) + 220 * np.exp(Vm / 50))
    dp_O_v = a_1 * (1 - p_O_v) - B_1 * p_O_v


    # Calcium gate
    B_2 = 0.00296
    a_3 = 0.0003
    B_3 = 0.03
    a_2blk = 2.24 * Ca_blk
    p_blk_C_1 = 1 - p_blk_C_2 - p_blk_O_c
    dp_blk_C_2 = a_3 * p_blk_C_1 - B_3 * p_blk_C_2 - a_2blk * p_blk_C_2 - B_2 * p_blk_O_c
    dp_blk_O_c = a_2blk * p_blk_C_2 - B_2 * p_blk_O_c

    a_2iz = 2.24 * Ca_iz
    dp_iz_O_c = a_2iz * p_iz_C_2 - B_2 * p_iz_O_c
    p_iz_C_1 = 1 - p_iz_C_2 - p_iz_O_c
    dp_iz_C_2 = a_3 * p_iz_C_1 - B_3 * p_iz_C_2 - a_2iz * p_iz_C_2 - B_2 * p_iz_O_c


    # Endo
    P_Kto_K = 0.01825
    P_Kto_Na = 0.09 * P_Kto_K
    p_O_Kto = y_1Kto * y_2Kto
    a_y1Kto = 1 / (9 * np.exp(-Vm / 20))
    B_y1Kto = 1 / (2.1 * np.exp(Vm / 60))
    a_y2Kto = 1 / (950 * np.exp(Vm / 500))
    B_y2Kto = (1 / (40 * np.exp(-Vm / 9))) + 13 * np.exp(-Vm / 1000)
    dy_1Kto = a_y1Kto * (1 - y_1Kto) - B_y1Kto * y_1Kto
    dy_2Kto = a_y2Kto * (1 - y_2Kto) - B_y2Kto * y_2Kto


    # Transient outward K+ current
    I_Kto_K = P_Kto_K * GHK_K * p_O_Kto
    I_Kto_Na = P_Kto_Na * GHK_Na * p_O_Kto
    I_Kto = I_Kto_Na + I_Kto_K


    # from huvec2
    #[TOK]

    #a_infinity = 1 / (1 + np.exp(-(Vm - 14.34) / 14.82))
    #tau_a = 1.0515 / (1 / (1.2089 * (1 + np.exp(-(Vm - 18.41) / 29.38))) + 3.5 / (1 + np.exp((Vm + 100) / 29.38)))
    #da = (a_infinity - a) / tau_a
    #i_infinity = 1 / (1 + np.exp((Vm + 43.94) / 5.711))
    #tau_i_fast = 4.562 + 1 / (0.3933 * np.exp(-(Vm + 100) / 100) + 0.08004 * np.exp((Vm + 50) / 16.59))
    #tau_i_slow = 23.62 + 1 / (0.001416 * np.exp(-(Vm + 96.52) / 59.05) + 0.0000000000017808 * np.exp((Vm + 114.1) / 8.079))
    #A_i_fast = 1 / (1 + np.exp((Vm - 213.6) / 151.2))
    #A_i_slow = 1 - A_i_fast
    #di_fast = (i_infinity - i_fast) / tau_i_fast
    #di_slow = (i_infinity - i_slow) / tau_i_slow
    #i = A_i_fast * i_fast + A_i_slow * i_slow
    #p_O_Kto = a * i
    #I_Kto = G_Kto * p_O_Kto * (Vm - E_K)

    # Voltage-dependent potassium current
    I_Kpl = P_Kpl * (Vm / (1 - np.exp(-Vm / 13))) / GHK_K
    P_Kpl = 0.000043 * np.power((K_o / 5.4),0.16)

    # Calcium-activated background cation current
    f_I_Ca_iz = 0.1
    f_I_Ca_blk = 0.9
    P_I_Ca_Na = 0.00273
    p_O_blk = 1 / (1 + np.power((0.0012 / Ca_blk),3))
    p_O_iz = 1 / (1 + np.power((0.0012 / Ca_iz),3)) 
    P_I_Ca_K = P_I_Ca_Na
    I_I_Ca_Na_blk = P_I_Ca_Na * f_I_Ca_blk * GHK_Na * p_O_blk
    I_I_Ca_K_blk = P_I_Ca_K * f_I_Ca_blk * GHK_K * p_O_blk
    I_I_Ca_Na_iz = P_I_Ca_Na * f_I_Ca_iz * GHK_Na * p_O_iz
    I_I_Ca_K_iz = P_I_Ca_K * f_I_Ca_iz * GHK_K * p_O_iz
    I_I_Ca = I_I_Ca_Na_iz + I_I_Ca_K_iz + I_I_Ca_Na_blk + I_I_Ca_K_blk

    

    # ATP-sensitive potassium current
    P_O_KATP = 0.8 / (1 + np.power((ATPi / 0.1),2))
    X_IKATP = 0.0236 * np.power((K_o),0.24)
    G_KATP = 17.674
    I_KATP = G_KATP * (Vm - E_K) * P_O_KATP * X_IKATP

    # Sarcoplasmic reticulum Ca2+ pump (SERCA) current 
    a1 = 25900 * MgATP
    a2 = 2540 / (1 + np.power((K_dCai / Ca_blk),1.7))
    a3 = 5.35 / (1 + np.power((Ca_2_SRup / K_dCasr),1.7))
    B1 = 0.1972 / (1 + np.power((Ca_blk / K_dCai),1.7))
    B2 = (25435 * MgADP) / (1 + np.power((K_dCasr / Ca_2_SRup),1.7))
    B3 = 149* Pi
    Pi = 0.50872066859173026
    v_cyc = (6.86 * (a1 * a2 * a3 - B1 * B2 * B3)) / (a2 * a3 + B1 * a3 + B1 * B2 + a1 * a3 + B2 * a1 + B2 * B3 + a1 * a2 + B3 * B1 + B3 * a2) 

    # slow_time_dependent_potassium_current_y1_gate
    #[STDPC1]

    #alpha_y1 = 1 / ((85 * np.exp(-Vm / 10.5)) + (370 * np.exp(-Vm / 62)))
    #beta_y1 = 1 / ((1450 * np.exp(Vm / 20)) + (260 * np.exp(Vm / 100)))
    #dy1 = alpha_y1 * (1 - y1) - (beta_y1 * y1)

    # slow_time_dependent_potassium_current_y2_gate
    #[STDPC2]

    #alpha_y2 = 3.7 * Ca_i
    #beta_y2 = 0.004444
    #dy2 = alpha_y2 * (1 - y2) - (beta_y2 * y2)

    # transient_outward_current
    #[TOC]

    #P_to_K = 0.033
    #P_to_Na = 0.00297
    #i_to_K = P_to_K * CF_K * TOC1.y1,3 * TOC2.y2
    #i_to_Na = P_to_Na * CF_Na * TOC1.y1,3 * TOC2.y2
    #i_to = i_to_Na + i_to_K

    # transient_outward_current_y1_gate
    #[TOC1]

    #alpha_y1 = 1 / ((11 * np.exp(-Vm / 28)) + (0.2 * np.exp(-Vm / 400)))
    #beta_y1 = 1 / ((4.4 * np.exp(Vm / 16)) + (0.2 * np.exp(Vm / 500)))
    #dy1 = alpha_y1 * (1 - y1) - (beta_y1 * y1)

    # transient_outward_current_y2_gate
    #[TOC2]

    #alpha_y2 = (0.0038 * np.exp(-(Vm + 13.5) / 11.3)) / (1 + (0.051335 * np.exp(-(Vm + 13.5) / 11.3)))
    #beta_y2 = (0.0038 * np.exp((Vm + 13.5) / 11.3)) / (1 + (0.067083 * np.exp((Vm + 13.5) / 11.3)))
    #dy2 = alpha_y2 * (1 - y2) - (beta_y2 * y2)

    # background_NSC_current
    P_bNSC = 0.385
    i_bNSC_K = 0.4 * P_bNSC * CF_K
    i_bNSC_Na = P_bNSC * CF_Na
    i_bNSC = i_bNSC_K + i_bNSC_Na

    # background_Kpl_current
    P_Kpl = 0.00011 * np.power((K_o / 5.4),0.16)
    if Vm == -3:
        i_Kpl = P_Kpl * CF_K * 13.0077
    else:
        i_Kpl = (P_Kpl * CF_K * (Vm + 3)) / (1 - np.exp(-(Vm + 3) / 13))
    # background_lCa_current
    P_lCa = 0.11
    p_open = 1 / (1 + np.power((0.0012 / Ca_i),3))
    i_lCa_K = P_lCa * CF_K * p_open
    i_lCa_Na = P_lCa * CF_Na * p_open

    i_lCa = i_lCa_K + i_lCa_Na

    # background_KATP_current
    P_KATP = 0.0236
    N = 2333
    gamma = P_KATP * N * np.power((K_o / 1),0.24)
    p_open = 0.8 / (1 + np.power((ATPi / 0.1),2))
    i_KATP = gamma * (Vm - E_K) * p_open 

    # sodium_calcium_exchanger_y_gate
    alpha_y = k2 * p_E2Na + k4 * p_E2Ca 
    beta_y = k1 * p_E1Na + k3 * p_E1Ca
    dy6 = alpha_y * (1 - y6) - (beta_y * y6)

    I_tot_Ca = I_tot_Ca_jnc + I_tot_Ca_iz + I_tot_Ca_blk
    I_tot_Na = (I_CaL_Na_jnc + I_CaL_Na_iz + I_CaL_Na_blk) + (I_NCX_Na_iz + I_NCX_Na_blk) + (I_Ks_Na_iz + I_Ks_Na_blk) + I_NaT_Na + I_NaL_Na + I_NaK_Na + I_Kto_Na + I_bNSC_Na + (I_LCCa_Na_iz + I_LCCa_Na_blk)
    I_tot_K = (I_CaL_K_jnc + I_CaL_K_iz + I_CaL_K_blk) + I_NaT_K + I_NaL_K + I_K1_K + I_Kr_K + (I_Ks_K_iz + I_Ks_K_blk) + I_Kto_K + I_Kpl + I_NaK_K + I_KATP_K + I_bNSC_K + (I_LCCa_K_iz + I_LCCa_K_blk)
    I_tot_cell = I_tot_Na + I_tot_Ca + I_tot_K
 
    if current_time >= 50 and current_time <= 53:
        i_stim = -12
        dVm = -(I_tot_cell + i_stim)
    else:
        i_stim = 0
        dVm = -(I_tot_cell + i_stim)
    # print(dVm)
    dK_i = ((I_tot_K + i_stim) * Cm) / (V_cyt * F)
    # update time
    current_time += dt

    # integrate
    Vm_Next = Vm + dVm*dt
    Na_i2_Next = Na_i + dNa_i*dt
    K_i2_Next = K_i + dK_i*dt
    Ca_2_tot_jnc2_Next = Ca_2_tot_jnc + dCa_2_tot_jnc*dt
    Ca_2_tot_iz2_Next = Ca_2_tot_iz + dCa_2_tot_iz*dt
    Ca_2_tot_blk2_Next = Ca_2_tot_blk + dCa_2_tot_blk*dt
    Ca_jnc2_Next = Ca_jnc + dCa_jnc*dt
    Ca_iz2_Next = Ca_iz + dCa_iz*dt
    Ca_blk2_Next = Ca_blk + dCa_blk*dt
    Ca_2_SRup2_Next = Ca_2_SRup + dCa_2_SRup*dt
    Ca_2_tot_SRrl2_Next = Ca_2_tot_SRrl + dCa_2_tot_SRrl*dt
    TnChCa2_Next = TnChCa2 + dTnChCa2*dt
    CaMCa2_Next = CaMCa2 + dCaMCa2*dt
    SRCa2_Next = SRCa2 + dSRCa2*dt
    TnChCa3_Next = TnChCa3 + dTnChCa3*dt
    CaMCa3_Next = CaMCa3 + dCaMCa3*dt
    SRCa3_Next = SRCa3 + dSRCa3*dt
    # Lb_jnc = 0.032868783807565416
    # Lb_iz = 0.011563775702743251
    # Hb_jnc = 0.226665087396588
    # Hb_iz = 0.099458321251370288
    # internal_ion_concentrations 
    Ca_Total2_Next = Ca_Total2 + dCa_Total2*dt
    ATPi2_Next = ATPi + dATPi*dt
    p_RP_Na2_Next = p_RP_Na + dp_RP_Na*dt
    p_AP_Na2_Next = p_AP_Na + dp_AP_Na*dt
    p_AI_Na2_Next = p_AI_Na + dp_AI_Na*dt
    y2_Next = y2 + dy2*dt
    p_RP_CaL2_Next = p_RP_CaL + dp_RP_CaL*dt
    p_AP_CaL2_Next = p_AP_CaL + dp_AP_CaL*dt
    p_AI_CaL2_Next = p_AI_CaL + dp_AI_CaL*dt
    p_U2_Next = p_U + dp_U*dt
    p_UCa2_Next = p_UCa + dp_UCa*dt
    p_C2_Next = p_C + dp_C*dt
    y3_Next = y3 + dy3*dt
    y4_Next = y4 + dy4*dt
    y5_Next = y5 + dy5*dt
    y6_Next = y6 + dy6*dt
    y7_Next = y7 + dy7*dt
    y8_Next = y8 + dy8*dt
    p_open_RyR2_Next = p_open_RyR * dp_open_RyR*dt
    p_close_RyR2_Next = p_close_RyR * dp_close_RyR*dt
    Ca_Total3_Next = Ca_Total3 * dCa_Total3*dt
    Caup2_Next = Caup * dCaup*dt
    X2_Next = X * dX*dt
    pCa2_Next = pCa * dpCa*dt
    pCaCB2_Next = pCaCB * dpCaCB*dt
    pCB2_Next = pCB * dpCB*dt
    p_O_NaT2_Next = p_O_NaT2 * dp_O_NaT2*dt
    p_O_NaT3_Next = p_O_NaT3 * dp_O_NaT3*dt
    p_I_2_NaT2_Next = p_I_2_NaT * dp_I_2_NaT*dt
    p_I_s_NaT2_Next = p_I_s_NaT * dp_I_s_NaT*dt
    p_O_NaL2_Next = p_O_NaL2 * dp_O_NaL2*dt
    p_O_NaL3_Next = p_O_NaL3 * dp_O_NaL3*dt
    p_I_1_NaL2_Next = p_I_1_NaL * dp_I_1_NaL*dt
    p_I_2_NaL2_Next = p_I_2_NaL + dp_I_2_NaL*dt
    p_I_s_NaL2_Next = p_I_s_NaL + dp_I_s_NaL*dt
    chi_r_fast2_Next = chi_r_fast + dchi_r_fast*dt
    chi_r_slow2_Next = chi_r_slow + dchi_r_slow*dt
    P_72_Next = P_7 + dP_7*dt
    P_8_132_Next = P_8_13 + dP_8_13*dt
    P_1_62_Next = P_1_6 + dP_1_6*dt
    p_E_1_NCX_blk2_Next = p_E_1_NCX_blk + dp_E_1_NCX_blk*dt
    p_I_1_NCX_blk2_Next = p_I_1_NCX_blk + dp_I_1_NCX_blk*dt
    p_I_2_NCX_blk2_Next = p_I_2_NCX_blk + dp_I_2_NCX_blk*dt
    p_E_1_NCX_iz2_Next = p_E_1_NCX_iz + dp_E_1_NCX_iz*dt
    p_I_1_NCX_iz2_Next = p_I_1_NCX_iz + dp_I_1_NCX_iz*dt
    p_I_2_NCX_iz2_Next = p_I_2_NCX_iz + dp_I_2_NCX_iz*dt
    Y_ooo2_Next = Y_ooo + dY_ooo*dt
    Y_ooc2_Next = Y_ooc + dY_ooc*dt
    Y_coo2_Next = Y_coo + dY_coo*dt
    Y_coc2_Next = Y_coc + dY_coc*dt
    Y_cco2_Next = Y_cco + dY_cco*dt
    Y_oco2_Next = Y_oco + dY_oco*dt
    Y_occ2_Next = Y_occ + dY_occ*dt
    Y_co_iz2_Next = Y_co_iz + dY_co_iz*dt
    Y_oo_iz2_Next = Y_oo_iz + dY_oo_iz*dt
    Y_oc_iz2_Next = Y_oc_iz + dY_oc_iz*dt
    Y_co_blk2_Next = Y_co_blk + dY_co_blk*dt
    Y_oo_blk2_Next = Y_oo_blk + dY_oo_blk*dt
    Y_oc_blk2_Next = Y_oc_blk + dY_oc_blk*dt
    TSCa_32_Next = TSCa_3 + dTSCa_3*dt
    TSCa_3W2_Next = TSCa_3W + dTSCa_3W*dt
    TSCa_3S2_Next = TSCa_3S + dTSCa_3S*dt
    TS_S2_Next = TS_S + dTS_S*dt
    TS_W2_Next = TS_W + dTS_W*dt
    hw2_Next = hw + dhw*dt
    hp2_Next = hp + dhp*dt
    Pb_spm2_Next = Pb_spm2 + dPb_spm2*dt
    Pb_spm3_Next = Pb_spm3 + dPb_spm3*dt
    # random 
    p_C_NaT2_Next = p_C_NaT + dp_C_NaT*dt
    p_C_NaL2_Next = p_C_NaL + dp_C_NaL*dt
    y9_Next = y9 + dy9*dt
    y10_Next = y10 + dy10*dt
    y11_Next = y11 + dy11*dt
    p_O_v2_Next = p_O_v + dp_O_v*dt
    p_blk_O_c2_Next = p_blk_O_c + dp_blk_O_c*dt
    p_blk_C_22_Next = p_blk_C_2 + dp_blk_C_2*dt
    p_iz_O_c2_Next = p_iz_O_c + dp_iz_O_c*dt
    p_iz_C_22_Next = p_iz_C_2 + dp_iz_C_2*dt
    y_1Kto2_Next = y_1Kto + dy_1Kto*dt
    y_2Kto2_Next = y_2Kto + dy_2Kto*dt

    # update values
    Vm = Vm_Next
    Na_i = Na_i2_Next 
    K_i = K_i2_Next
    Ca_2_tot_jnc = Ca_2_tot_jnc2_Next
    Ca_2_tot_iz = Ca_2_tot_iz2_Next
    Ca_2_tot_blk = Ca_2_tot_blk2_Next
    Ca_jnc = Ca_jnc2_Next
    Ca_iz = Ca_iz2_Next
    Ca_blk = Ca_blk2_Next
    Ca_2_SRup = Ca_2_SRup2_Next
    Ca_2_tot_SRrl = Ca_2_tot_SRrl2_Next
    TnChCa2 = TnChCa2_Next
    CaMCa2 = CaMCa2_Next
    SRCa2 = SRCa2_Next
    TnChCa3 = TnChCa3_Next
    CaMCa3 = CaMCa3_Next
    SRCa3 = SRCa3_Next
    # Lb_jnc = 0.032868783807565416
    # Lb_iz = 0.011563775702743251
    # Hb_jnc = 0.226665087396588
    # Hb_iz = 0.099458321251370288
    # internal_ion_concentrations 
    Ca_Total2 = Ca_Total2_Next
    ATPi = ATPi2_Next
    p_RP_Na = p_RP_Na2_Next
    p_AP_Na = p_AP_Na2_Next
    p_AI_Na = p_AI_Na2_Next
    y2 = y2_Next
    p_RP_CaL = p_RP_CaL2_Next
    p_AP_CaL = p_AP_CaL2_Next
    p_AI_CaL = p_AI_CaL2_Next
    p_U = p_U2_Next
    p_UCa = p_UCa2_Next
    p_C = p_C2_Next
    y3 = y3_Next
    y4 = y4_Next
    y5 = y5_Next
    y6 = y6_Next
    y7 = y7_Next
    y8 = y8_Next
    p_open_RyR = p_open_RyR2_Next
    p_close_RyR = p_close_RyR2_Next
    Ca_Total3 = Ca_Total3_Next
    Caup = Caup2_Next
    X = X2_Next
    pCa = pCa2_Next
    pCaCB = pCaCB2_Next
    pCB = pCB2_Next
    p_O_NaT2 = p_O_NaT2_Next
    p_O_NaT3 = p_O_NaT3_Next
    p_I_2_NaT = p_I_2_NaT2_Next
    p_I_s_NaT = p_I_s_NaT2_Next
    p_O_NaL2 = p_O_NaL2_Next
    p_O_NaL3 = p_O_NaL3_Next
    p_I_1_NaL = p_I_1_NaL2_Next
    p_I_2_NaL = p_I_2_NaL2_Next
    p_I_s_NaL = p_I_s_NaL2_Next
    chi_r_fast = chi_r_fast2_Next
    chi_r_slow = chi_r_slow2_Next
    P_7 = P_72_Next
    P_8_13 = P_8_132_Next
    P_1_6 = P_1_62_Next
    p_E_1_NCX_blk = p_E_1_NCX_blk2_Next
    p_I_1_NCX_blk = p_I_1_NCX_blk2_Next
    p_I_2_NCX_blk = p_I_2_NCX_blk2_Next
    p_E_1_NCX_iz = p_E_1_NCX_iz2_Next
    p_I_1_NCX_iz = p_I_1_NCX_iz2_Next
    p_I_2_NCX_iz = p_I_2_NCX_iz2_Next
    Y_ooo = Y_ooo2_Next
    Y_ooc = Y_ooc2_Next
    Y_coo = Y_coo2_Next
    Y_coc = Y_coc2_Next
    Y_cco = Y_cco2_Next
    Y_oco = Y_oco2_Next
    Y_occ = Y_occ2_Next
    Y_co_iz = Y_co_iz2_Next
    Y_oo_iz = Y_oo_iz2_Next
    Y_oc_iz = Y_oc_iz2_Next
    Y_co_blk = Y_co_blk2_Next
    Y_oo_blk = Y_oo_blk2_Next
    Y_oc_blk = Y_oc_blk2_Next
    TSCa_3 = TSCa_32_Next
    TSCa_3W = TSCa_3W2_Next
    TSCa_3S = TSCa_3S2_Next
    TS_S = TS_S2_Next
    TS_W = TS_W2_Next
    hw = hw2_Next
    hp = hp2_Next
    Pb_spm2 = Pb_spm2_Next
    Pb_spm3 = Pb_spm3_Next
    # random 
    p_C_NaT = p_C_NaT2_Next
    p_C_NaL = p_C_NaL2_Next
    y9 = y9_Next
    y10 = y10_Next
    y11 = y11_Next
    p_O_v = p_O_v2_Next
    p_blk_O_c = p_blk_O_c2_Next
    p_blk_C_2 = p_blk_C_22_Next
    p_iz_O_c = p_iz_O_c2_Next
    p_iz_C_2 = p_iz_C_22_Next
    y_1Kto = y_1Kto2_Next
    y_2Kto = y_2Kto2_Next

    times.append(current_time)
    Vm2.append(Vm_Next)

# Euler method call
plt.figure()
plt.title('Action potential')
plt.plot(times, Vm2)
# plt.xlim([0,1001])
# plt.ylim([-100,60])
plt.show()
