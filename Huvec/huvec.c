#include <stdio.h>
#include <math.h>

double R = 8.3143;
double T = 310;
double F = 96.4867;
double Vm = -91.4466885079348;
double dCaMCadt;
double k_on_CaM = 34;
double Ca_2_blk;
double B_tot_CaM = 0.024;
double CaMCa = 0.000228581865602447;
double k_off_CaM = 0.238;
double dTnChCadt;
double k_on_TnCh = 2.37;
double B_tot_TnCh = 0.12;
double TnChCa = 0.110742559707052;
double k_off_TnCh = 0.000032;
double dSRCadt;
double k_on_SR = 100;
double B_tot_SR = 0.0171;
double SRCa = 0.00172960014640511;
double k_off_SR = 0.06;
double L_free_iz;
double L_bound_iz = 0.0075621764602356;
double B_tot_L_iz = 0.6078;
double Ca_2_iz;
double k_off_L_iz = 1.3;
double k_on_L_iz = 100;
double K_dL_iz;
double H_free_iz;
double H_bound_iz = 0.0769149150028914;
double B_tot_H_iz = 0.2178;
double k_off_H_iz = 0.03;
double k_on_H_iz = 100;
double K_dH_iz;
double Ca_2_tot_iz = 0.084640522722006;
double L_free_jnc;
double B_tot_L_jnc = 1.1095;
double Ca_2_jnc = 0;
double k_off_L_jnc = 1.3;
double k_on_L_jnc = 100;
double K_dL_jnc;
double H_free_jnc;
double B_tot_H_jnc = 0.398;
double k_off_H_jnc = 0.03;
double k_on_H_jnc = 100;
double K_dH_jnc;
double Ca_2_tot_jnc = 0.207176351449979;
double B_tot_CSQN = 3;
double Ca_2_tot_SRrl = 2.21876221622152;
double k_off_CSQN = 65;
double k_on_CSQN = 100;
double K_d_CSQN_Ca;
double Ca_2_SRrl;
double J_Ca_jnciz;
double G_dCa_jnciz = 3395.88;
double Sc_Cell = 1;
double J_Ca_izblk;
double G_dCa_izblk = 3507.78;
double J_trans_SR;
double P_trans = 4.8037;
double Ca_2_SRup = 0.761077662687456;
double alpha_plus;
double alpha_minus;
double epsilon_plus_iz;
double epsilon_plus_blk;
double epsilon_plus_iz_loc;
double epsilon_plus_blk_loc;
double epsilon_minus;
double Ca_2_iz_loc;
double Ca_2_blk_loc;
double T_L = 147.51;
double K_L = 0.0044;
double ATPfactor;
double ATP = 6;
double p_O_LCC;
double P_CaL_Ca = 14.21;
double P_CaL_Na;
double P_CaL_K;
double GHK_Ca_LR;
double GHK_Ca_L0;
double GHK_Ca_iz;
double GHK_Ca_blk;
double GHK_Na;
double GHK_K;
double Ca_2_nd_LR;
double Ca_2_nd_L0;
double Cao = 1.8;
double Nao = 140;
double Ko = 4.5;
double Nai = 6.66894310282034;
double Ki = 139.238265011042;
double I_CaL_Ca_blk;
double I_CaL_Ca_iz;
double I_CaL_Ca_LR;
double I_CaL_Ca_L0;
double I_CaL_Na_blk;
double I_CaL_Na_iz;
double I_CaL_Na_jnc;
double I_CaL_K_blk;
double I_CaL_K_iz;
double I_CaL_K_jnc;
double f_CaL_blk = 0.1;
double f_CaL_iz = 0.15;
double f_CaL_jnc = 0.75;
double I_CaL;
double I_Na;
double I_NaT;
double I_NaL;
double k_OC;
double p_O_NaT = 0.000000706725155695262;
double k_I2C;
double p_I_2_NaT = 0.0117704053067285;
double k_Isb;
double p_I_s_NaT = 0.304002781414015;
double k_Isf;
double f_C_Na;
double k_C2O;
double p_C_NaT;
double k_OI2;
double k_I2O = 0.0001312;
double k_C2I2;
double k_I1I2 = 0.00534;
double k_OI1;
double k_I1O = 0.01;
double k_I1C;
double k_C2I1;
double f_LSM = 0.13125;
double P_Na = 8.1072;
double dp_O_NaTdt;
double dp_I_2_NaTdt;
double dp_I_s_NaTdt;
double p_O_NaL = 0.00000295214591324261;
double p_I_1_NaL = 0.00254273877063925;
double p_I_2_NaL = 0.0118261382165599;
double p_I_s_NaL = 0.303220346353844;
double dp_O_NaLdt;
double p_C_NaL;
double dp_I_1_NaLdt;
double dp_I_2_NaLdt;
double dp_I_s_NaLdt;
double alpha_Mg;
double E_K;
double E_Ca_jnc;
double E_Ca_blk;
double E_Ca_iz;
double E_Na;
double beta_Mg;
double Mg_2_cyt = 0.8;
double f_O;
double f_B;
double po_Mg;
double po_Mg1;
double po_Mg2;
double alpha_SPM;
double beta_SPM;
double SPM = 5.0;
double G_K1 = 1.353;
double f_mode1 = 0.9;
double dPb_spmdt;
double Pb_spm = 0.594875991179992;
double po_mode1;
double po_mode2;
double p_O_K1;
double chi_K1;
double I_K1;
double chi_r_infinity;
double tau_chi_r_fast;
double tau_chi_r_slow;
double A_chi_r_fast;
double A_chi_r_slow;
double dchi_r_fastdt;
double chi_r_fast = 0.00000486210633393005;
double dchi_r_slowdt;
double chi_r_slow = 0.437041249050081;
double chi_r;
double R_Kr;
double p_O_Kr;
double G_Kr = 0.0166;
double chi_Kr;
double I_Kr;
double f_Ks_iz = 0.1;
double f_Ks_blk = 0.9;
double P_Ks_K = 0.002782;
double P_Ks_Na;
double para_Xs1_infinity;
double tau_Xs1;
double dpara_Xs1dt;
double para_Xs1 = 0.277482694590328;
double para_Xs2_infinity;
double para_Xs2 = 0.000131110342877451;
double tau_Xs2;
double dpara_Xs2dt;
double para_RKs_blk;
double para_RKs_iz;
double p_O_Ks_blk;
double p_O_Ks_iz;
double I_Ks_K_blk;
double I_Ks_K_iz;
double I_Ks_Na_blk;
double I_Ks_Na_iz;
double I_Ks;
double a_infinity;
double tau_a;
double dadt;
double a = 0.000793627635934239;
double i_infinity;
double tau_i_fast;
double tau_i_slow;
double A_i_fast;
double A_i_slow;
double di_fastdt;
double i_fast = 0.999756080468878;
double di_slowdt;
double i_slow = 0.575995954010486;
double i;
double p_O_Kto;
double G_Kto = 0.0312;
double I_Kto;
double p_O_Kpl;
double chi_Kpl;
double P_Kpl = 0.0000172;
double I_Kpl;
double f_Cab_blk = 0.9;
double P_Cab = 0.00006822;
double I_Cab_blk;
double f_Cab_iz = 0.1;
double I_Cab_iz;
double I_Cab;
double P_bNSC_K = 0.00014;
double I_bNSC_K;
double P_bNSC_Na = 0.00035;
double I_bNSC_Na;
double I_bNSC;
double p_O_blk;
double p_O_iz;
double f_l_Ca_blk = 0.9;
double P_l_Ca_Na = 0.00273;
double I_l_Ca_Na_blk;
double f_l_Ca_iz = 0.1;
double I_l_Ca_Na_iz;
double P_l_Ca_K;
double I_l_Ca_K_blk;
double I_l_Ca_K_iz;
double I_l_Ca;
double ATP_cyt = 6.67701543987464;
double p_O_KATP;
double chi_KATP;
double G_KATP = 17.674;
double I_KATP;
double delta_Nai = -0.14;
double K_d_Nai_0 = 5;
double K_d_Nai;
double Nai_bar;
double delta_Nao = 0.44;
double K_d_Nao_0 = 26.8;
double K_d_Nao;
double Nao_bar;
double delta_Ki = -0.14;
double K_d_Ki_0 = 18.8;
double K_d_Ki;
double Ki_bar;
double delta_Ko = 0.23;
double K_d_Ko_0 = 0.8;
double K_d_Ko;
double Ko_bar;
double K_d_MgATP = 0.6;
double MgATP_cyt = 6.631643709767415;
double MgATP_bar;
double k_1_plus = 0.72;
double k_1_minus = 0.08;
double k_2_plus = 0.08;
double k_2_minus = 0.008;
double k_3_plus = 4;
double k_3_minus = 8000;
double k_4_plus = 0.3;
double k_4_minus = 0.2;
double alpha_1_plus;
double alpha_2_plus;
double alpha_3_plus;
double alpha_4_plus;
double alpha_1_minus;
double alpha_2_minus;
double Pi = 0.50872066859173026;
double H = 0.0001;
double alpha_3_minus;
double alpha_4_minus;
double dP_7dt;
double P_7 = 0.0831770174499825;
double P_8_13 = 0.281082409575779;
double P_1_6 = 0.435289193632868;
double P_14_15;
double dP_8_13dt;
double dP_1_6dt;
double V_step1;
double V_step2;
double V_step3;
double V_step4;
double v_cyc_NaK;
double Amp_NaK = 25.178;
double I_NaK;
double Stoi_NaK_Na = 3;
double Stoi_NaK_K = -2;
double I_NaK_Na;
double I_NaK_K;
double f_Caina_blk;
double K_m_act = 0.004;
double alpha_1_on = 0.002;
double alpha_1_off = 0.0015;
double alpha_2_on = 0.00006;
double alpha_2_off = 0.02;
double beta_1_on = 0.0012;
double beta_1_off = 0.0000005;
double beta_2_on = 0.18;
double beta_2_off = 0.0002;
double K_m_Nai = 20.74854;
double K_m_Nao = 87.5;
double K_m_Cai = 0.0184;
double K_m_Cao = 1.38;
double q_blk_E_1_Na;
double alpha_1_blk;
double f_Caina_iz;
double q_iz_E_1_Na;
double alpha_1_iz;
double beta_1_blk;
double beta_1_iz;
double alpha_2_blk;
double alpha_2_iz;
double beta_2_blk;
double beta_2_iz;
double k_1;
double k_2;
double k_3 = 1.0;
double k_4 = 1.0;
double q_blk_E_1_Ca;
double q_iz_E_1_Ca;
double q_E_2_Na;
double q_E_2_Ca;
double alpha_E;
double beta_E_blk;
double beta_E_iz;
double p_E_1_NCX_blk = 0.111872123711613;
double p_I_2_NCX_blk = 0.684869019924837;
double p_I_1_NCX_blk = 0.203023555446362;
double p_E_2_NCX_blk;
double dp_E_1_NCX_blkdt;
double dp_I_1_NCX_blkdt;
double dp_I_2_NCX_blkdt;
double p_E_2_NCX_iz;
double p_E_1_NCX_iz = 0.238718640001014;
double p_I_1_NCX_iz = 0.13771129457898;
double p_I_2_NCX_iz = 0.622892868847556;
double dp_E_1_NCX_izdt;
double dp_I_1_NCX_izdt;
double dp_I_2_NCX_izdt;
double v_cyc_NCX_blk;
double v_cyc_NCX_iz;
double Amp_NCX = 30.53;
double f_NCX_blk = 0.9;
double I_NCX_blk;
double f_NCX_iz = 0.1;
double I_NCX_iz;
double I_NCX;
double I_NCX_Na_blk;
double I_NCX_Na_iz;
double I_NCX_Ca_blk;
double I_NCX_Ca_iz;
double K_m = 0.0005;
double Amp_PMCA = 0.19;
double f_PMCA_blk = 0.9;
double f_PMCA_iz = 0.1;
double I_PMCA_blk;
double I_PMCA_iz;
double I_PMCA;
double p_O_RyR;
double Ca_2_nd_00;
double delta_RTF;
double J_L = 0.000913;
double J_R = 0.02;
double g_D = 0.065;
double f_L;
double f_R;
double Ca_2_nd_0R;
double Q_10 = 3;
double k_co_L0;
double k_co_00;
double k_oc;
double f_t_00;
double f_t_L0;
double sloc0 = 0.1;
double f_n = 7;
double k_rco_0;
double k_rco_L;
double N_RyR = 10;
double p_C_0;
double p_C_L;
double k_roc_0;
double k_roc_L;
double k_co_0R;
double k_co_LR;
double epsilon_plus_00;
double epsilon_plus_L0;
double epsilon_plus_0R;
double epsilon_plus_LR;
double Y_ooo = 0.00000172489315884865;
double Y_ooc = 0.00000142034754677507;
double Y_occ = 0.0000000249594301562175;
double Y_coc = 0.992110534408681;
double Y_ccc;
double Y_coo = 0.0000138422676498755;
double Y_cco = 0.0000000953816272498217;
double Y_oco = 0.00000000000156949238162028;
double dY_ooodt;
double dY_oocdt;
double dY_coodt;
double dY_cocdt;
double dY_ccodt;
double dY_ocodt;
double dY_occdt;
double p_O_RyR_base = 0.000075;
double p_O_RyR_t;
double P_RyR = 5191;
double J_Ca_rel;
double K_dCai = 0.0027;
double alpha_1;
double alpha_2;
double K_dCasr = 1.378;
double alpha_3;
double beta_1;
double MgADP_cyt = 0.025978226605534577;
double beta_2;
double beta_3;
double v_cyc;
double Amp_SERCA = 106.4448;
double J_SERCA;
double I_tot_Ca_jnc;
double I_tot_Ca_iz;
double I_tot_Ca_blk;
double I_tot_Ca;
double I_NaT_Na;
double I_NaT_K;
double I_NaL_Na;
double I_NaL_K;
double I_Kto_Na = 0;
double I_tot_Na;
double I_tot_K;
double I_tot_cell;
double I_app = 0;
double dVmdt;
double Y_cc_iz;
double Y_co_iz = 0.992251726297519;
double Y_oo_iz = 0.00000314564543512061;
double Y_oc_iz = 0.000000024556270151713;
double dY_co_izdt;
double dY_oo_izdt;
double dY_oc_izdt;
double Y_cc_blk;
double Y_co_blk = 0.992424981547859;
double Y_oo_blk = 0.00000314619469048683;
double Y_oc_blk = 0.0000000240070147854924;
double dY_co_blkdt;
double dY_oo_blkdt;
double dY_oc_blkdt;
double Cm;
double V_cell;
double V_jnc;
double V_iz;
double V_blk;
double V_SRt;
double V_SRrl;
double V_SRup;
double V_cyt;
double dCa_2_tot_jncdt;
double dCa_2_tot_izdt;
double dCa_2_tot_blkdt;
double dCa_2_SRupdt;
double dCa_2_tot_SRrldt;
double dNaidt;
double dKidt;
double Ca_2_tot_blk = 0.11279654524634;
double halfSL = 0.91;
double hp = 0.00600014761511324;
double TS_tot = 23;
double propFh = 28000;
double Za = 0.0023;
double Yv = 1.5;
double rate_g;
double Yd = 0.0333;
double Yc = 1;
double Lc = 1.2;
double Yvd;
double rate_gd;
double Zb = 0.1397;
double Yb = 0.1816;
double rate_f = 0.0023;
double convertF = 15;
double eqvhalfSL = 1.15;
double TS;
double dTSCa_3dt;
double TSCa_3 = 0.00899891910620064;
double Zp = 0.2095;
double Yp = 0.1397;
double dTSCa_3Wdt;
double TSCa_3W = 0.000369547640656701;
double Zr = 7.2626;
double Yr = 0.1397;
double dTSCa_3Sdt;
double TSCa_3S = 0.000153834503967436;
double Zq = 0.3724;
double Yq = 0.2328;
double dTS_Sdt;
double TS_S = 0.000876347322180234;
double hwr = 0.0001;
double dTS_Wdt;
double TS_W = 0.000492054058977473;
double dhwdt;
double rate_B = 0.5;
double hw = 0.000100147615113241;
double hpr = 0.006;
double dhpdt;
double dt = 0.01;

void euler() {
	CaMCa += dCaMCadt * dt;
	TnChCa += dTnChCadt * dt;
	SRCa += dSRCadt * dt;
	p_O_NaT += dp_O_NaTdt * dt;
	p_I_2_NaT += dp_I_2_NaTdt * dt;
	p_I_s_NaT += dp_I_s_NaTdt * dt;
	p_O_NaL += dp_O_NaLdt * dt;
	p_I_1_NaL += dp_I_1_NaLdt * dt;
	p_I_2_NaL += dp_I_2_NaLdt * dt;
	p_I_s_NaL += dp_I_s_NaLdt * dt;
	chi_r_fast += dchi_r_fastdt * dt;
	chi_r_slow += dchi_r_slowdt * dt;
	para_Xs1 += dpara_Xs1dt * dt;
	para_Xs2 += dpara_Xs2dt * dt;
	i_fast += di_fastdt * dt;
	i_slow += di_slowdt * dt;
	P_7 += dP_7dt * dt;
	P_8_13 += dP_8_13dt * dt;
	P_1_6 += dP_1_6dt * dt;
	p_E_1_NCX_blk += dp_E_1_NCX_blkdt * dt;
	p_I_1_NCX_blk += dp_I_1_NCX_blkdt * dt;
	p_I_2_NCX_blk += dp_I_2_NCX_blkdt * dt;
	p_E_1_NCX_iz += dp_E_1_NCX_izdt * dt;
	p_I_1_NCX_iz += dp_I_1_NCX_izdt * dt;
	p_I_2_NCX_iz += dp_I_2_NCX_izdt * dt;
	Y_ooo += dY_ooodt * dt;
	Y_ooc += dY_oocdt * dt;
	Y_coo += dY_coodt * dt;
	Y_coc += dY_cocdt * dt;
	Y_cco += dY_ccodt * dt;
	Y_oco += dY_ocodt * dt;
	Y_occ += dY_occdt * dt;
	Vm += dVmdt * dt;
	Y_co_iz += dY_co_izdt * dt;
	Y_oo_iz += dY_oo_izdt * dt;
	Y_oc_iz += dY_oc_izdt * dt;
	Y_co_blk += dY_co_blkdt * dt;
	Y_oo_blk += dY_oo_blkdt * dt;
	Y_oc_blk += dY_oc_blkdt * dt;
	Ca_2_tot_jnc += dCa_2_tot_jncdt * dt;
	Ca_2_tot_iz += dCa_2_tot_izdt * dt;
	Ca_2_tot_blk += dCa_2_tot_blkdt * dt;
	Ca_2_SRup += dCa_2_SRupdt * dt;
	Ca_2_tot_SRrl += dCa_2_tot_SRrldt * dt;
	Nai += dNaidt * dt;
	Ki += dKidt * dt;
	TSCa_3 += dTSCa_3dt * dt;
	TSCa_3W += dTSCa_3Wdt * dt;
	TSCa_3S += dTSCa_3Sdt * dt;
	TS_S += dTS_Sdt * dt;
	TS_W += dTS_Wdt * dt;
	hw += dhwdt * dt;
	hp += dhpdt * dt;
	Pb_spm += dPb_spmdt * dt;
	a += dadt * dt;
}

void contraction() {
	TS = TS_tot - TSCa_3 - TSCa_3W - TSCa_3S - TS_S - TS_W;
	rate_g = Za + Yv * (1 - exp(-propFh * pow(hw - hwr, 2)));
	rate_gd = Yd + Yc * pow(halfSL - Lc, 2) + Yvd * (1 - exp(-propFh * pow(hw - hwr, 2)));
	dTSCa_3dt = Yb * TS * pow(Ca_2_blk * 1000, 3) - Zb * TSCa_3 + rate_g * TSCa_3W - rate_f * exp(-convertF * pow(halfSL - eqvhalfSL, 2)) * TSCa_3;
	dTSCa_3Wdt = rate_f * exp(-convertF * pow(halfSL - eqvhalfSL, 2)) * TSCa_3 - rate_g * TSCa_3W + Zp * TSCa_3S - Yp * TSCa_3W;
	dTSCa_3Sdt = Yp * TSCa_3W - Zp * TSCa_3S + Zr * TS_S * pow(Ca_2_blk * 1000, 3) - Yr * TSCa_3S;
	dTS_Sdt = Yr * TSCa_3S - Zr * TS_S * pow(Ca_2_blk * 1000, 3) + Zq * TS_W - Yq * TS_S;
	dTS_Wdt = Yq * TS_S - Zq * TS_W - rate_gd * TS_W;
	dhwdt = -rate_B * (hw - hwr);
	dhpdt = -rate_B * (hp - hpr);
}

void ionConcentration() {
	dCa_2_tot_jncdt = -I_tot_Ca_jnc * Cm / (V_jnc * 2 * F) + J_Ca_rel / V_jnc - J_Ca_jnciz / V_jnc;
	dCa_2_tot_izdt = -I_tot_Ca_iz * Cm / (V_iz * 2 * F) + J_Ca_jnciz / V_iz - J_Ca_izblk / V_iz;
	dCa_2_tot_blkdt = -I_tot_Ca_blk * Cm / (V_blk * 2 * F) - J_SERCA / V_blk + J_Ca_izblk / V_blk;
	dCa_2_SRupdt = J_SERCA / V_SRup - J_trans_SR / V_SRup;
	dCa_2_tot_SRrldt = J_trans_SR / V_SRrl - J_Ca_rel / V_SRrl;
	dNaidt = -I_tot_Na * Cm / (V_cyt * F);
	dKidt = -(I_tot_K + I_app) * Cm / (V_cyt * F);
}

void membranePotential() {
	I_tot_K = (I_CaL_K_jnc + I_CaL_K_iz + I_CaL_K_blk) + I_NaT_K + I_NaL_K + I_K1 + I_Kr + (I_Ks_K_iz + I_Ks_K_blk) + I_Kto + I_Kpl + I_NaK_K + I_KATP + I_bNSC_K + (I_l_Ca_K_iz + I_l_Ca_K_blk);
	I_tot_Na = (I_CaL_Na_jnc + I_CaL_Na_iz + I_CaL_Na_blk) + (I_NCX_Na_iz + I_NCX_Na_blk) + (I_Ks_Na_iz + I_Ks_Na_blk) + I_NaT_Na + I_NaL_Na + I_NaK_Na + I_Kto_Na + I_bNSC_Na + (I_l_Ca_Na_iz + I_l_Ca_Na_blk);
	I_tot_Ca_blk = I_CaL_Ca_blk + I_PMCA_blk + I_NCX_Ca_blk + I_Cab_blk;
	I_tot_Ca_iz = I_CaL_Ca_iz + I_PMCA_iz + I_NCX_Ca_iz + I_Cab_iz;
	I_tot_Ca_jnc = I_CaL_Ca_LR + I_CaL_Ca_L0;
	I_tot_Ca = I_tot_Ca_jnc + I_tot_Ca_iz + I_tot_Ca_blk;
	I_tot_cell = I_tot_Na + I_tot_Ca + I_tot_K;
	dVmdt = -(I_tot_cell + I_app);
}

void SERCA() {
	alpha_1 = 25900 * MgATP_cyt;
	alpha_2 = 2540 / (1 + pow(K_dCai / Ca_2_blk, 1.7));
	alpha_3 = 5.35 / (1 + pow(Ca_2_SRup / K_dCasr, 1.7));
	beta_1 = 0.1972 / (1 + pow(Ca_2_blk / K_dCai, 1.7));
	beta_2 = 25435 * MgADP_cyt / (1 + pow(K_dCasr / Ca_2_SRup, 1.7));
	beta_3 = 149 * Pi;
	v_cyc = 6.86 * (alpha_1 * alpha_2 * alpha_3 - beta_1 * beta_2 * beta_3) / (alpha_2 * alpha_3 + beta_1 * alpha_3 + beta_1 * beta_2 + alpha_1 * alpha_3 + beta_2 * alpha_1 + beta_2 * beta_3 + alpha_1 * alpha_2 + beta_3 * beta_1 + beta_3 * alpha_2);
	J_SERCA = Amp_SERCA * v_cyc / (2 * F) * Sc_Cell;
}

void CaRU() {
	p_O_RyR = Y_ooo + Y_coo + Y_cco + Y_oco;
	p_O_RyR_t = p_O_RyR + p_O_RyR_base;

	J_Ca_rel = P_RyR * p_O_RyR_t * (Ca_2_SRrl - Ca_2_jnc) * Sc_Cell;

	Ca_2_nd_00 = Ca_2_jnc;
	Ca_2_nd_L0 = (Ca_2_nd_00 + f_L * delta_RTF * Vm * exp(-delta_RTF * Vm) / (1 - exp(-delta_RTF * Vm)) * Cao) / (1 + f_L * delta_RTF * Vm / (1 - exp(-delta_RTF * Vm)));
	Ca_2_nd_0R = (Ca_2_nd_00 + f_R * Ca_2_SRrl) / (1 + f_R);
	Ca_2_nd_LR = (Ca_2_nd_00 + f_R * Ca_2_SRrl + f_L * delta_RTF * Vm * exp(-delta_RTF * Vm) / (1 - exp(-delta_RTF * Vm)) * Cao) / (1 + f_R + f_L * delta_RTF * Vm / (1 - exp(-delta_RTF * Vm)));

	epsilon_plus_00 = (Ca_2_nd_00 * alpha_plus) / (T_L * K_L);
	epsilon_plus_L0 = (Ca_2_nd_L0 * alpha_plus) / (T_L * K_L);
	epsilon_plus_0R = (Ca_2_nd_0R * alpha_plus) / (T_L * K_L);
	epsilon_plus_LR = (Ca_2_nd_LR * alpha_plus) / (T_L * K_L);

	k_co_00 = Q_10 * 0.4 / (1 + pow(0.025 / Ca_2_nd_00, 2.7));
	k_co_L0 = Q_10 * 0.4 / (1 + pow(0.025 / Ca_2_nd_L0, 2.7));
	k_co_0R = Q_10 * 0.4 / (1 + pow(0.025 / Ca_2_nd_0R, 2.7));
	k_co_LR = Q_10 * 0.4 / (1 + pow(0.025 / Ca_2_nd_LR, 2.7));
	f_t_00 = k_co_00 / (k_co_00 + k_oc);
	f_t_L0 = k_co_L0 / (k_co_L0 + k_oc);
	k_rco_0 = f_n * f_t_00 * k_co_0R * (sloc0 + Ca_2_SRrl);
	k_rco_L = f_n * f_t_L0 * k_co_LR * (sloc0 + Ca_2_SRrl);
	p_C_0 = k_oc / (k_oc + f_t_00 * (k_rco_0 / (f_n * f_t_00)));
	p_C_L = k_oc / (k_oc + f_t_00 * (k_rco_L / (f_n * f_t_L0)));
	k_roc_0 = k_oc * pow(p_C_0, (N_RyR - 1) * 0.74);
	k_roc_L = k_oc * pow(p_C_L, (N_RyR - 1) * 0.74);

	Y_ccc = 1 - (Y_ooo + Y_ooc + Y_coo + Y_coc + Y_cco + Y_oco + Y_occ);
	dY_ooodt = k_rco_L * Y_ooc + alpha_plus * Y_coo + epsilon_minus * Y_oco - (k_roc_L + alpha_minus + epsilon_plus_LR) * Y_ooo;
	dY_oocdt = alpha_plus * Y_coc + k_roc_L * Y_ooo + epsilon_minus * Y_occ - (alpha_minus + k_rco_L + epsilon_plus_L0) * Y_ooc;
	dY_coodt = k_rco_0 * Y_coc + alpha_minus * Y_ooo + epsilon_minus * Y_cco - (k_roc_0 + alpha_plus + epsilon_plus_0R) * Y_coo;
	dY_cocdt = k_roc_0 * Y_coo + alpha_minus * Y_ooc + epsilon_minus * Y_ccc - (k_rco_0 + alpha_plus + epsilon_plus_00) * Y_coc;
	dY_ccodt = k_rco_0 * Y_ccc + alpha_minus * Y_oco + epsilon_plus_0R * Y_coo - (k_roc_0 + alpha_plus + epsilon_minus) * Y_cco;
	dY_ocodt = k_rco_0 * Y_occ + alpha_plus * Y_cco + epsilon_plus_LR * Y_ooo - (k_roc_0 + alpha_minus + epsilon_minus) * Y_oco;
	dY_occdt = alpha_plus * Y_ccc + k_roc_0 * Y_oco + epsilon_plus_L0 * Y_ooc - (alpha_minus + k_rco_0 + epsilon_minus) * Y_occ;
}

void currentPMCA() {
	I_PMCA_blk = f_PMCA_blk * Amp_PMCA * pow(Ca_2_blk, 1.6) / (pow(K_m, 1.6) + pow(Ca_2_blk, 1.6));
	I_PMCA_iz = f_PMCA_iz * Amp_PMCA * pow(Ca_2_iz, 1.6) / (pow(K_m, 1.6) + pow(Ca_2_iz, 1.6));
	I_PMCA = I_PMCA_iz + I_PMCA_blk;
}

void currentNCX() {
	f_Caina_blk = Ca_2_blk / (Ca_2_blk + K_m_act);
	f_Caina_iz = Ca_2_iz / (Ca_2_iz + K_m_act);
	q_blk_E_1_Na = 1.0 / (1.0 + pow(K_m_Nai / Nai, 3) * (1.0 + Ca_2_blk / K_m_Cai));
	q_iz_E_1_Na = 1.0 / (1.0 + pow(K_m_Nai / Nai, 3) * (1.0 + Ca_2_iz / K_m_Cai));
	q_blk_E_1_Ca = 1.0 / (1.0 + K_m_Cai / Ca_2_blk * (1.0 + pow(Nai / K_m_Nai, 3)));
	q_iz_E_1_Ca = 1.0 / (1.0 + K_m_Cai / Ca_2_iz * (1.0 + pow(Nai / K_m_Nai, 3)));
	q_E_2_Na = 1.0 / (1.0 + pow(K_m_Nao / Nao, 3) * (1.0 + Cao / K_m_Cao));
	q_E_2_Ca = 1.0 / (1.0 + K_m_Cao / Cao * (1.0 + pow(Nao / K_m_Nao, 3)));
	alpha_1_blk = q_blk_E_1_Na * (f_Caina_blk * alpha_1_on + (1 - f_Caina_blk) * alpha_1_off);
	alpha_1_iz = q_iz_E_1_Na * (f_Caina_iz * alpha_1_on + (1 - f_Caina_iz) * alpha_1_off);
	beta_1_blk = f_Caina_blk * beta_1_on + (1 - f_Caina_blk) * beta_1_off;
	beta_1_iz = f_Caina_iz * beta_1_on + (1 - f_Caina_iz) * beta_1_off;
	alpha_2_blk = f_Caina_blk * alpha_2_on + (1 - f_Caina_blk) * alpha_2_off;
	alpha_2_iz = f_Caina_iz * alpha_2_on + (1 - f_Caina_iz) * alpha_2_off;
	beta_2_blk = f_Caina_blk * beta_2_on + (1 - f_Caina_blk) * beta_2_off;
	beta_2_iz = f_Caina_iz * beta_2_on + (1 - f_Caina_iz) * beta_2_off;
	k_1 = exp((0.32 * F * Vm) / (R * T));
	k_2 = exp(((0.32 - 1) * F * Vm) / (R * T));
	alpha_E = k_2 * q_E_2_Na + k_4 * q_E_2_Ca;
	beta_E_blk = k_1 * q_blk_E_1_Na + k_3 * q_blk_E_1_Ca;
	beta_E_iz = k_1 * q_iz_E_1_Na + k_3 * q_iz_E_1_Ca;
	p_E_2_NCX_blk = 1 - p_E_1_NCX_blk - p_I_1_NCX_blk - p_I_2_NCX_blk;
	p_E_2_NCX_iz = 1 - p_E_1_NCX_iz - p_I_1_NCX_iz - p_I_2_NCX_iz;
	v_cyc_NCX_blk = k_1 * q_blk_E_1_Na * p_E_1_NCX_blk - k_2 * q_E_2_Na * p_E_2_NCX_blk;
	v_cyc_NCX_iz = k_1 * q_iz_E_1_Na * p_E_1_NCX_iz - k_2 * q_E_2_Na * p_E_2_NCX_iz;
	dp_E_1_NCX_blkdt = p_E_2_NCX_blk * alpha_E + p_I_1_NCX_blk * beta_1_blk + p_I_2_NCX_blk * beta_2_blk - p_E_1_NCX_blk * (beta_E_blk + alpha_1_blk + alpha_2_blk);
	dp_I_1_NCX_blkdt = p_E_1_NCX_blk * alpha_1_blk - p_I_1_NCX_blk * beta_1_blk;
	dp_I_2_NCX_blkdt = p_E_1_NCX_blk * alpha_2_blk - p_I_2_NCX_blk * beta_2_blk;
	dp_E_1_NCX_izdt = p_E_2_NCX_iz * alpha_E + p_I_1_NCX_iz * beta_1_iz + p_I_2_NCX_iz * beta_2_iz - p_E_1_NCX_iz * (beta_E_iz + alpha_1_iz + alpha_2_iz);
	dp_I_1_NCX_izdt = p_E_1_NCX_iz * alpha_1_iz - p_I_1_NCX_iz * beta_1_iz;
	dp_I_2_NCX_izdt = p_E_1_NCX_iz * alpha_2_iz - p_I_2_NCX_iz * beta_2_iz;
	I_NCX_blk = f_NCX_blk * Amp_NCX * v_cyc_NCX_blk;
	I_NCX_iz = f_NCX_iz * Amp_NCX * v_cyc_NCX_iz;
	I_NCX = I_NCX_iz + I_NCX_blk;
	I_NCX_Na_blk = 3 * I_NCX_blk;
	I_NCX_Na_iz = 3 * I_NCX_iz;
	I_NCX_Ca_blk = -2 * I_NCX_blk;
	I_NCX_Ca_iz = -2 * I_NCX_iz;
}

void currentNaK() {
	K_d_Nai = K_d_Nai_0 * exp((delta_Nai * F * Vm) / (R * T));
	K_d_Nao = K_d_Nao_0 * exp((delta_Nao * F * Vm) / (R * T));
	K_d_Ki = K_d_Ki_0 * exp((delta_Ki * F * Vm) / (R * T));
	K_d_Ko = K_d_Ko_0 * exp((delta_Ko * F * Vm) / (R * T));
	Nai_bar = Nai / K_d_Nai;
	Nao_bar = Nao / K_d_Nao;
	Ki_bar = Ki / K_d_Ki;
	Ko_bar = Ko / K_d_Ko;
	MgATP_bar = MgATP_cyt / K_d_MgATP;
	alpha_1_plus = k_1_plus * pow(Nai_bar, 3) / (pow(1 + Nai_bar, 3) + pow(1 + Ki_bar, 2) - 1);
	alpha_3_plus = k_3_plus * pow(Ko_bar, 2) / (pow(1 + Nao_bar, 3) + pow(1 + Ko_bar, 2) - 1);
	alpha_4_plus = k_4_plus * MgATP_bar / (1 + MgATP_bar);
	alpha_1_minus = k_1_minus * MgADP_cyt;
	alpha_2_minus = k_2_minus * pow(Nao_bar, 3) / (pow(1 + Nao_bar, 3) + pow(1 + Ko_bar, 2) - 1);
	alpha_3_minus = k_3_minus * Pi * H / (1 + MgATP_bar);
	alpha_4_minus = k_4_minus * pow(Ki_bar, 2) / (pow(1 + Nai_bar, 3) + pow(1 + Ki_bar, 2) - 1);
	P_14_15 = 1 - P_1_6 - P_7 - P_8_13;
	V_step1 = alpha_1_plus * P_1_6 - alpha_1_minus * P_7;
	V_step2 = alpha_2_plus * P_7 - alpha_2_minus * P_8_13;
	V_step3 = alpha_3_plus * P_8_13 - alpha_3_minus * P_14_15;
	V_step4 = alpha_4_plus * P_14_15 - alpha_4_minus * P_1_6;
	v_cyc_NaK = V_step2;
	dP_1_6dt = -alpha_1_plus * P_1_6 + alpha_1_minus * P_7 + alpha_4_plus * P_14_15 - alpha_4_minus * P_1_6;
	dP_7dt = -alpha_2_plus * P_7 + alpha_2_minus * P_8_13 + alpha_1_plus * P_1_6 - alpha_1_minus * P_7;
	dP_8_13dt = -alpha_3_plus * P_8_13 + alpha_3_minus * P_14_15 + alpha_2_plus * P_7 - alpha_2_minus * P_8_13;
	I_NaK = Amp_NaK * v_cyc_NaK;
	I_NaK_Na = Stoi_NaK_Na * I_NaK;
	I_NaK_K = Stoi_NaK_K * I_NaK;
}

void currentKATP() {
	p_O_KATP = 0.8 / (1.0 + pow(ATP_cyt / 0.1, 2));
	chi_KATP = 0.0236 * pow(Ko, 0.24);
	I_KATP = G_KATP * (Vm - E_K) * p_O_KATP * chi_KATP;
}

void currentlCa() {
	p_O_blk = 1.0 / (1.0 + pow(0.0012 / Ca_2_blk, 3));
	p_O_iz = 1.0 / (1.0 + pow(0.0012 / Ca_2_iz, 3));
	I_l_Ca_Na_blk = P_l_Ca_Na * f_l_Ca_blk * GHK_Na * p_O_blk;
	I_l_Ca_Na_iz = P_l_Ca_Na * f_l_Ca_iz * GHK_Na * p_O_iz;
	I_l_Ca_K_blk = P_l_Ca_K * f_l_Ca_blk * GHK_K * p_O_blk;
	I_l_Ca_K_iz = P_l_Ca_K * f_l_Ca_iz * GHK_K * p_O_iz;
	I_l_Ca = I_l_Ca_Na_iz + I_l_Ca_K_iz + I_l_Ca_Na_blk + I_l_Ca_K_blk;
}

void currentbNSC() {
	I_bNSC_K = P_bNSC_K * GHK_K;
	I_bNSC_Na = P_bNSC_Na * GHK_Na;
	I_bNSC = I_bNSC_K + I_bNSC_Na;
}

void currentCab() {
	I_Cab_blk = P_Cab * f_Cab_blk * GHK_Ca_blk;
	I_Cab_iz = P_Cab * f_Cab_iz * GHK_Ca_iz;
	I_Cab = I_Cab_iz + I_Cab_blk;
}

void currentKpl() {
	p_O_Kpl = Vm / (1 - exp(-Vm / 13.0));
	chi_Kpl = pow(Ko / 5.4, 0.16);
	I_Kpl = P_Kpl * chi_Kpl * p_O_Kpl * GHK_K;
}

void currentKto() {
	a_infinity = 1 / (1 + exp(-(Vm - 14.34) / 14.82));
	tau_a = 1.0515 / (1 / (1.2089 * (1 + exp(-(Vm - 18.41) / 29.38))) + 3.5 / (1 + exp((Vm + 100) / 29.38)));
	dadt = (a_infinity - a) / tau_a;
	i_infinity = 1 / (1 + exp((Vm + 43.94) / 5.711));
	tau_i_fast = 4.562 + 1 / (0.3933 * exp(-(Vm + 100) / 100) + 0.08004 * exp((Vm + 50) / 16.59));
	tau_i_slow = 23.62 + 1 / (0.001416 * exp(-(Vm + 96.52) / 59.05) + 0.000000017808 * exp((Vm + 114.1) / 8.079));
	A_i_fast = 1 / (1 + exp((Vm - 213.6) / 151.2));
	A_i_slow = 1 - A_i_fast;
	di_fastdt = (i_infinity - i_fast) / tau_i_fast;
	di_slowdt = (i_infinity - i_slow) / tau_i_slow;
	i = A_i_fast * i_fast + A_i_slow * i_slow;
	p_O_Kto = a * i;
	I_Kto = G_Kto * p_O_Kto * (Vm - E_K);
}

void currentKs() {
	para_Xs1_infinity = 1 / (1 + exp(-(Vm + 11.60) / 8.932));
	tau_Xs1 = 817.3 + 1 / (0.0002326 * exp((Vm + 48.28) / 17.80) + 0.001292 * exp(-(Vm + 210.0) / 230.0));
	dpara_Xs1dt = (para_Xs1_infinity - para_Xs1) / tau_Xs1;
	para_Xs2_infinity = para_Xs1_infinity;
	tau_Xs2 = 1 / (0.01 * exp((Vm - 50) / 20) + 0.0193 * exp(-(Vm + 66.54) / 31));
	dpara_Xs2dt = (para_Xs2_infinity - para_Xs2) / tau_Xs2;
	para_RKs_blk = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_blk, 1.4));
	para_RKs_iz = 1 + 0.6 / (1 + pow(0.000038 / Ca_2_iz, 1.4));
	p_O_Ks_blk = para_Xs1 * para_Xs2 * para_RKs_blk;
	p_O_Ks_iz = para_Xs1 * para_Xs2 * para_RKs_iz;
	I_Ks_K_blk = f_Ks_blk * P_Ks_K * GHK_K * p_O_Ks_blk;
	I_Ks_K_iz = f_Ks_iz * P_Ks_K * GHK_K * p_O_Ks_iz;
	I_Ks_Na_blk = f_Ks_blk * P_Ks_Na * GHK_Na * p_O_Ks_blk;
	I_Ks_Na_iz = f_Ks_iz * P_Ks_Na * GHK_Na * p_O_Ks_iz;
	I_Ks = I_Ks_K_blk + I_Ks_K_iz + I_Ks_Na_blk + I_Ks_Na_iz;
}

void currentKr() {
	chi_r_infinity = 1 / (1 + exp(-(Vm + 8.337) / 6.789));
	tau_chi_r_fast = 12.98 + 1 / (0.3652 * exp((Vm - 31.66) / 3.869) + 0.00004123 * exp(-(Vm - 47.78) / 20.38));
	tau_chi_r_slow = 1.865 + 1 / (0.06629 * exp((Vm - 34.70) / 7.355) + 0.00001128 * exp(-(Vm - 19.74) / 25.94));
	A_chi_r_fast = 1 / (1 + exp((Vm + 4.81) / 38.21));
	A_chi_r_slow = 1 - A_chi_r_fast;
	dchi_r_fastdt = (chi_r_infinity - chi_r_fast) / tau_chi_r_fast;
	dchi_r_slowdt = (chi_r_infinity - chi_r_slow) / tau_chi_r_slow;
	chi_r = A_chi_r_fast * chi_r_fast + A_chi_r_slow * chi_r_slow;
	R_Kr = 1 / ((1 + exp((Vm + 55) / 75)) * (1 + exp((Vm - 10) / 30)));
	p_O_Kr = chi_r * R_Kr;
	chi_Kr = sqrt(Ko / 4.5);
	I_Kr = G_Kr * chi_Kr * (Vm - E_K) * p_O_Kr;
}

void currentK1() {
	alpha_Mg = 12.0 * exp(-0.025 * (Vm - E_K));
	beta_Mg = 28 * Mg_2_cyt * exp(0.025 * (Vm - E_K));
	f_O = alpha_Mg / (alpha_Mg + beta_Mg);
	f_B = beta_Mg / (alpha_Mg + beta_Mg);
	po_Mg = f_O * f_O * f_O;
	po_Mg1 = 3.0 * f_O * f_O * f_B;
	po_Mg2 = 3.0 * f_O * f_B * f_B;
	alpha_SPM = 0.17 * exp(-0.07 * (Vm - E_K + 8 * Mg_2_cyt)) / (1.0 + 0.01 * exp(0.12 * (Vm - E_K + 8 * Mg_2_cyt)));
	beta_SPM = 0.28 * SPM * exp(0.15 * (Vm - E_K + 8 * Mg_2_cyt)) / (1.0 + 0.01 * exp(0.13 * (Vm - E_K + 8 * Mg_2_cyt)));

	dPb_spmdt = beta_SPM * po_Mg * (1 - Pb_spm) - alpha_SPM * Pb_spm;
	po_mode1 = f_mode1 * (1 - Pb_spm) * (po_Mg + (2.0 / 3.0) * po_Mg1 + (1.0 / 3.0) * po_Mg2);
	po_mode2 = (1 - f_mode1) / (1.0 + SPM / (40.0 * exp(-(Vm - E_K) / 9.1)));
	p_O_K1 = po_mode1 + po_mode2;
	chi_K1 = pow(Ko / 4.5, 0.4) / (1.0 + exp(-(Ko - 2.2) / 0.6));
	I_K1 = G_K1 * chi_K1 * (Vm - E_K) * p_O_K1;
}

void currentNa() {
	f_C_Na = 1 / (1 + exp(-(Vm + 48) / 7));
	k_C2O = 1 / (0.0025 * exp(Vm / -8.0) + 0.15 * exp(Vm / -100.0));
	k_OC = 1 / (30 * exp(Vm / 12.0) + 0.53 * exp(Vm / 50.0));
	k_OI2 = 1 / (0.0433 * exp(Vm / -27.0) + 0.34 * exp(Vm / -2000.0));
	k_C2I2 = 0.5 / (1.0 + (k_I2O * k_OC) / (k_OI2 * k_C2O));
	k_I2C = 0.5 - k_C2I2;
	k_Isb = 1 / (300000.0 * exp(Vm / 10.0) + 50000 * exp(Vm / 16.0));
	k_Isf = 1 / (0.016 * exp(Vm / -9.9) + 8.0 * exp(Vm / -45.0));
	k_OI1 = k_OI2;
	k_I1C = k_I2C;
	k_C2I1 = k_C2I2;
	p_C_NaT = 1.0 - p_I_2_NaT - p_I_s_NaT - p_O_NaT;
	dp_O_NaTdt = k_I2O * p_I_2_NaT + f_C_Na * k_C2O * p_C_NaT - (k_OC + k_OI2) * p_O_NaT;
	dp_I_2_NaTdt = f_C_Na * k_C2I2 * p_C_NaT + k_OI2 * p_O_NaT + k_Isb * p_I_s_NaT - (k_I2C + k_I2O + k_Isf) * p_I_2_NaT;
	dp_I_s_NaTdt = k_Isf * p_I_2_NaT + k_Isf * p_C_NaT - 2 * k_Isb * p_I_s_NaT;
	p_C_NaL = 1.0 - p_I_2_NaL - p_I_s_NaL - p_I_1_NaL - p_O_NaL;
	dp_O_NaLdt = k_I1O * p_I_1_NaL + f_C_Na * k_C2O * p_C_NaL - (k_OC + k_OI1) * p_O_NaL;
	dp_I_1_NaLdt = k_OI1 * p_O_NaL + f_C_Na * k_C2I1 * p_C_NaL - (k_I1O + k_I1C + k_I1I2) * p_I_1_NaL;
	dp_I_2_NaLdt = f_C_Na * k_C2I2 * p_C_NaL + k_I1I2 * p_I_1_NaL + k_Isb * p_I_s_NaL - (k_I2C + k_Isf) * p_I_2_NaL;
	dp_I_s_NaLdt = k_Isf * p_I_2_NaL + k_Isf * p_C_NaL - 2 * k_Isb * p_I_s_NaL;

	I_NaT_Na = (1 - f_LSM) * P_Na * GHK_Na * p_O_NaT;
	I_NaT_K = (1 - f_LSM) * P_Na * 0.1 * GHK_K * p_O_NaT;
	I_NaT = I_NaT_Na + I_NaT_K;
	I_NaL_Na = f_LSM * P_Na * GHK_Na * p_O_NaL;
	I_NaL_K = f_LSM * P_Na * 0.1 * GHK_K * p_O_NaL;
	I_NaL = I_NaL_Na + I_NaL_K;
	I_Na = I_NaT + I_NaL;
}

void currentCaL() {
	alpha_plus = 1 / (3.734 * exp(-Vm / 8.5) + 0.35 * exp(-Vm / 3500));
	alpha_minus = 1 / (4.65 * exp(Vm / 15) + 1.363 * exp(Vm / 100));

	epsilon_plus_iz = (Ca_2_iz * alpha_plus) / (T_L * K_L);
	epsilon_plus_blk = (Ca_2_blk * alpha_plus) / (T_L * K_L);
	Ca_2_iz_loc = (Ca_2_iz + f_L * delta_RTF * Vm * exp(-delta_RTF * Vm) / (1 - exp(-delta_RTF * Vm)) * Cao) / (1 + f_L * delta_RTF * Vm / (1 - exp(-delta_RTF * Vm)));
	Ca_2_blk_loc = (Ca_2_blk + f_L * (delta_RTF * Vm * exp(-delta_RTF * Vm)) / (1 - exp(-delta_RTF * Vm)) * Cao) / (1 + f_L * delta_RTF * Vm / (1 - exp(-delta_RTF * Vm)));
	epsilon_plus_iz_loc = (Ca_2_iz_loc * alpha_plus) / (T_L * K_L);
	epsilon_plus_blk_loc = (Ca_2_blk_loc * alpha_plus) / (T_L * K_L);
	epsilon_minus = 1 / (8084 * exp(Vm / 10) + 158 * exp(Vm / 1000)) + 1 / (134736 * exp(-Vm / 5) + 337 * exp(-Vm / 2000));

	p_O_LCC = Y_ooo + Y_ooc;

	ATPfactor = 1 / (1 + pow(1.4 / ATP, 3));

	E_Ca_jnc = (R * T / F) / 2 * log(Cao / Ca_2_nd_L0);
	E_Ca_blk = (R * T / F) / 2 * log(Cao / Ca_2_blk);
	E_Ca_iz = (R * T / F) / 2 * log(Cao / Ca_2_iz);
	E_K = (R * T / F) / 1 * log(Ko / Ki);
	E_Na = (R * T / F) / 1 * log(Nao / Nai);

	GHK_Ca_LR = 2 * Vm / (R * T / F) * (Ca_2_nd_LR - Cao * exp(-2 * Vm / (R * T / F))) / (1 - exp(-2 * Vm / (R * T / F)));
	GHK_Ca_L0 = 2 * Vm / (R * T / F) * (Ca_2_nd_L0 - Cao * exp(-2 * Vm / (R * T / F))) / (1 - exp(-2 * Vm / (R * T / F)));
	GHK_Ca_iz = 2 * Vm / (R * T / F) * (Ca_2_iz - Cao * exp(-2 * Vm / (R * T / F))) / (1 - exp(-2 * Vm / (R * T / F)));
	GHK_Ca_blk = 2 * Vm / (R * T / F) * (Ca_2_blk - Cao * exp(-2 * Vm / (R * T / F))) / (1 - exp(-2 * Vm / (R * T / F)));
	GHK_Na = 1 * Vm / (R * T / F) * (Nai - Nao * exp(-1 * Vm / (R * T / F))) / (1 - exp(-1 * Vm / (R * T / F)));
	GHK_K = 1 * Vm / (R * T / F) * (Ki - Ko * exp(-1 * Vm / (R * T / F))) / (1 - exp(-1 * Vm / (R * T / F)));

	I_CaL_Ca_blk = f_CaL_blk * P_CaL_Ca * GHK_Ca_blk * Y_oo_blk * ATPfactor;
	I_CaL_Ca_iz = f_CaL_iz * P_CaL_Ca * GHK_Ca_iz * Y_oo_iz * ATPfactor;
	I_CaL_Ca_LR = f_CaL_jnc * P_CaL_Ca * GHK_Ca_LR * Y_ooo * ATPfactor;
	I_CaL_Ca_L0 = f_CaL_jnc * P_CaL_Ca * GHK_Ca_L0 * Y_ooc * ATPfactor;
	I_CaL_Na_blk = f_CaL_blk * P_CaL_Na * GHK_Na * Y_oo_blk * ATPfactor;
	I_CaL_Na_iz = f_CaL_iz * P_CaL_Na * GHK_Na * Y_oo_iz * ATPfactor;
	I_CaL_Na_jnc = f_CaL_jnc * P_CaL_Na * GHK_Na * p_O_LCC * ATPfactor;
	I_CaL_K_blk = f_CaL_blk * P_CaL_K * GHK_K * Y_oo_blk * ATPfactor;
	I_CaL_K_iz = f_CaL_iz * P_CaL_K * GHK_K * Y_oo_iz * ATPfactor;
	I_CaL_K_jnc = f_CaL_jnc * P_CaL_K * GHK_K * p_O_LCC * ATPfactor;

	I_CaL = (I_CaL_Ca_LR + I_CaL_Ca_L0 + I_CaL_Na_jnc + I_CaL_K_jnc) + (I_CaL_Ca_iz + I_CaL_Na_iz + I_CaL_K_iz) + (I_CaL_Ca_blk + I_CaL_Na_blk + I_CaL_K_blk);

	Y_cc_iz = 1 - (Y_co_iz + Y_oo_iz + Y_oc_iz);
	dY_co_izdt = epsilon_minus * Y_cc_iz + alpha_minus * Y_oo_iz - (epsilon_plus_iz + alpha_plus) * Y_co_iz;
	dY_oo_izdt = alpha_plus * Y_co_iz + epsilon_minus * Y_oc_iz - (alpha_minus + epsilon_plus_iz_loc) * Y_oo_iz;
	dY_oc_izdt = epsilon_plus_iz_loc * Y_oo_iz + alpha_plus * Y_cc_iz - (epsilon_minus + alpha_minus) * Y_oc_iz;

	Y_cc_blk = 1 - (Y_co_blk + Y_oo_blk + Y_oc_blk);
	dY_co_blkdt = epsilon_minus * Y_cc_blk + alpha_minus * Y_oo_blk - (epsilon_plus_blk + alpha_plus) * Y_co_blk;
	dY_oo_blkdt = alpha_plus * Y_co_blk + epsilon_minus * Y_oc_blk - (alpha_minus + epsilon_plus_blk_loc) * Y_oo_blk;
	dY_oc_blkdt = epsilon_plus_blk_loc * Y_oo_blk + alpha_plus * Y_cc_blk - (epsilon_minus + alpha_minus) * Y_oc_blk;
}

void boundaryDiffusion() {
	J_Ca_jnciz = G_dCa_jnciz * (Ca_2_jnc - Ca_2_iz) * Sc_Cell;
	J_Ca_izblk = G_dCa_izblk * (Ca_2_iz - Ca_2_blk) * Sc_Cell;
	J_trans_SR = P_trans * (Ca_2_SRup - Ca_2_SRrl) * Sc_Cell;
}

void bulkSpace() {
	dCaMCadt = k_on_CaM * Ca_2_blk * (B_tot_CaM - CaMCa) - k_off_CaM * CaMCa;
	dTnChCadt = k_on_TnCh * Ca_2_blk * (B_tot_TnCh - TnChCa) - k_off_TnCh * TnChCa;
	dSRCadt = k_on_SR * Ca_2_blk * (B_tot_SR - SRCa) - k_off_SR * SRCa;
	Ca_2_blk = Ca_2_tot_blk - (CaMCa + TnChCa + SRCa + 3 * (TSCa_3 + TSCa_3W + TSCa_3S) / 1000);
}

void intermediateZone() {
	L_free_iz = B_tot_L_iz - L_bound_iz;
	H_free_iz = B_tot_H_iz - H_bound_iz;
	for(int n = 0; n < 10; n++) {
		Ca_2_iz = Ca_2_tot_iz / (1 + L_free_iz / K_dL_iz + H_free_iz / K_dH_iz);
		L_free_iz = B_tot_L_iz / (1 + Ca_2_iz / K_dL_iz);
		H_free_iz = B_tot_H_iz / (1 + Ca_2_iz / K_dH_iz);
	}
	L_bound_iz = B_tot_L_iz - L_free_iz;
	H_bound_iz = B_tot_H_iz - H_free_iz;
}

void junctionalSpace() {
	for(int n = 0; n < 10; n++) {
		Ca_2_jnc = Ca_2_tot_jnc / (1 + L_free_jnc / K_dL_jnc + H_free_jnc / K_dH_jnc);
		L_free_jnc = B_tot_L_jnc / (1 + Ca_2_jnc / K_dL_jnc);
		H_free_jnc = B_tot_H_jnc / (1 + Ca_2_jnc / K_dH_jnc);
	}
}

void releaseSiteOfSR() {
	double a = 1;
	double b = B_tot_CSQN - Ca_2_tot_SRrl + K_d_CSQN_Ca;
	double c = -K_d_CSQN_Ca * Ca_2_tot_SRrl;
	Ca_2_SRrl = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
}

void current() {
	currentCaL();
	currentNa();
	currentK1();
	currentKr();
	currentKs();
	currentKto();
	currentKpl();
	currentCab();
	currentbNSC();
	currentlCa();
	currentKATP();
	currentNaK();
	currentNCX();
	currentPMCA();
}

void ca2Buffer() {
	bulkSpace();
	intermediateZone();
	junctionalSpace();
	releaseSiteOfSR();
}

void calculate() {
	ca2Buffer();
	boundaryDiffusion();
	current();
	CaRU();
	SERCA();
	membranePotential();
	contraction();
	ionConcentration();
	euler();
}

int main(void) {
	K_dL_iz = k_off_L_iz / k_on_L_iz;
	K_dH_iz = k_off_H_iz / k_on_H_iz;
	K_dL_jnc = k_off_L_jnc / k_on_L_jnc;
	K_dH_jnc = k_off_H_jnc / k_on_H_jnc;
	K_d_CSQN_Ca = k_off_CSQN / k_on_CSQN;
	P_CaL_Na = 0.0000185 * P_CaL_Ca;
	P_CaL_K = 0.000367 * P_CaL_Ca;
	P_Ks_Na = 0.04 * P_Ks_K;
	P_l_Ca_K = P_l_Ca_Na;
	alpha_2_plus = k_2_plus;
	delta_RTF = 2 * F / (R * T);
	f_L = J_L / g_D;
	f_R = J_R / g_D;
	k_oc = Q_10 * 0.5664;
	Cm = 192.46 * Sc_Cell;
	V_cell = 120 * (37.62 * Sc_Cell) * 8.4;
	V_jnc = 0.008 * V_cell;
	V_iz = 0.035 * V_cell;
	V_blk = 0.68 * V_cell;
	V_SRt = 0.06 * V_cell;
	V_SRrl = 0.2 * V_SRt;
	V_SRup = 0.8 * V_SRt;
	V_cyt = V_jnc + V_iz + V_blk;
	Yvd = Yv;

	FILE *output;
	output = fopen("data.csv", "w");
	if(output == NULL) {
		printf("I/O error.\n");
		return -1;
	}

	for(double t = 0; t < 1000; t += dt) {
		fprintf(output, "%lf, %lf\n", t, Vm);
		if(50 <= t && t < 53) {
			I_app = -12;
		} else {
			I_app = 0;
		}
		calculate();
	}

	fclose(output);

	return 0;
}
