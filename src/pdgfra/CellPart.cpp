// CellPart.cpp : implementation file
// ICC model (headless port of 07_GiSim/PDGFa/CellPart.cpp, MFC removed)
// Changes vs. the original are tagged [headless].

/*
2015. 05
ICC project

*/
#include "CellPart.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>

/////////////////////////////////////////////////////////////////////////////
// CCellPart

CCellPart::CCellPart() : CCellPart(1u) {}

CCellPart::CCellPart(unsigned s)
{
	Reseed(s);
	// ------ initialize calculation 
	timeX = 0;
	timeX0 = 0;
	m_dDt = 5.000000e-003;
	ElapsedTime = 0;
	iRepeat = 1;	
	InitialCellVariables();
}

/////////////////////////////////////////////////////////////////////////////
// CCellPart commands


void CCellPart::Na_Current()
{	
	double dActCurve = 1.0 / (1.0 + exp(-(Vm + 57.99642) / 7.15595));		// my own conversion of above equation			
	double dActTau = (1.0 / (0.00196*exp(-Vm / 4.35885) + 8.05823*exp(Vm / 312.03936)) + 1.0 / (0.01875*exp(-Vm / 16.75279) + 100.0*exp(Vm / 11.34255)));

	double dInactCurve = 1.0 / (1.0 + exp((Vm + 84.74079) / 5.5584));
	double dInactTau = 1.0*1.0 / (0.00868*exp(-Vm / 1000.0) + 1.05102*exp(Vm / 195.94687)) + 1.0 / (6.22342E-7*exp(-Vm / 5.539) + 100.0*exp(Vm / 9.12919));		// remrem

	double dInactCurve2 = dInactCurve;
	double dInactTau2 = 200.0*1.0 / (0.012*exp(-Vm / 18.33336) + 0.11219*exp(Vm / 91.62613)) + 1.0 / (5.28434E-4*exp(-Vm / 14.96535) + 47.46058*exp(Vm / 7.33849));	// remrem

	double dNaAct_dt = (dActCurve - NaAct) / dActTau;
	double dNaInact_dt = (dInactCurve - NaInact) / dInactTau;
	double dNaInact2_dt = (dInactCurve2 - NaInact2) / dInactTau2;

	//----------- Integrate ------------------------------------------------//
	NaAct = TwoState(NaAct, dNaAct_dt, dActCurve / dActTau, (1.0 - dActCurve) / dActTau, m_dDt);
	NaInact = TwoState(NaInact, dNaInact_dt, dInactCurve / dInactTau, (1.0 - dInactCurve) / dInactTau, m_dDt);
	NaInact2 = TwoState(NaInact2, dNaInact2_dt, dInactCurve2 / dInactTau2, (1.0 - dInactCurve2) / dInactTau2, m_dDt);

	//---------- Current Amplitude ------------------------------------------//	
	double PNa = SF_Na*2.0*Cap*0.091;					// maximal INa conductance, 4.55 nS
	double ENa;
	ENa = 49.2; //mRTF * log(Nae/Nai);						//remrem REVERSAL POTENTIAL FOR Na	, 					// should be restored
	INa_Na = PNa * (0.9*NaAct*NaAct*NaAct*NaInact + 0.1*NaAct*NaAct*NaAct*NaInact2) * (Vm - ENa);
	INa_K = 0;
	INa = INa_Na + INa_K;
}
void CCellPart::CaL_Current()
{		
	double dAlfaAct = 0.0054;
	if (Vm != -30.0)	dAlfaAct = 0.002175*(Vm + 30.0) / (1.0 - exp(-(Vm + 30.0) / 2.5));         // alpha, activation parameter for L-type Ca channel

	double dBetaAct = 0.0016;
	if (Vm != 0.0)	dBetaAct = 0.0006315*Vm / (exp(Vm / 2.5) - 1.0);

	double dAlfaInact = 9.9977e-006;
	if (Vm != -34.0) dAlfaInact = 1.775e-6*(Vm + 34.0) / (exp((Vm + 34.0) / 5.633) - 1.0);
	double dBetaInact = 0.427*Cai*(Vm + 64.0) / (1.0 + exp(-(Vm + 44.0) / 4.16));       //beta

	double dCaAct_dt = dAlfaAct*(1 - CaAct) - dBetaAct*CaAct;
	double dCaInact_dt = dAlfaInact*(1 - CaInact) - dBetaInact*CaInact;

	//----------- Integrate ------------------------------------------------//
	CaAct = TwoState(CaAct, dCaAct_dt, dAlfaAct, dBetaAct, m_dDt);
	CaInact = TwoState(CaInact, dCaInact_dt, dAlfaInact, dBetaInact, m_dDt);

	//---------- Current Amplitude ------------------------------------------//	
	double PCaL = SF_CaL*57.6;
	ICaL_Ca = PCaL*CaAct*CaInact*CF_Ca;
	ICaL_Na = PCaL*5.0e-5*CaAct*CaInact*CF_Na;
	ICaL_K = PCaL*1.0e-3*CaAct*CaInact*CF_K;
	ICaL = (ICaL_Ca + ICaL_Na + ICaL_K);
	PoCaL = CaAct*CaInact;
}
void CCellPart::VDDR_Current()		
{
	double dY1_dt, dY2_dt;
	double alpha_Y1, beta_Y1, alpha_Y2, beta_Y2;
	double PVDDR;

	double DTinfinite = 1.0 / (1.0 + exp(-(Vm + 26.0) / 6.0));
	double tauDT = 0.6 + 5.4 / (1.0 + exp(0.03*(Vm + 100.0)));
	alpha_Y1 = DTinfinite / tauDT;
	beta_Y1 = (1.0 - DTinfinite) / tauDT;

	double FTinfinite = 1.0 / (1 + exp((Vm + 66.0) / 6.0));
	double tauFT = 10.0 + 400.0 / (1.0 + exp(0.02*(Vm + 65.0)));
	alpha_Y2 = FTinfinite / tauFT;
	beta_Y2 = (1.0 - FTinfinite) / tauFT;

	dY1_dt = alpha_Y1*(1 - YVDDR1) - beta_Y1*YVDDR1;
	dY2_dt = alpha_Y2*(1 - YVDDR2) - beta_Y2*YVDDR2;

	//----------- Integrate ------------------------------------------------//
	YVDDR1 = TwoState(YVDDR1, dY1_dt, alpha_Y1, beta_Y1, m_dDt);
	YVDDR2 = TwoState(YVDDR2, dY2_dt, alpha_Y2, beta_Y2, m_dDt);

	PVDDR = 0.5*15.84;
	IVDDR = SF_VDDR*PVDDR*CF_Ca*YVDDR1*YVDDR2;
}
void CCellPart::ClCa_Current()		// 
{			
	double ECatZero = 1.39e-3;		// 1.39 uM	
	double kcclca = 0.44*0.01248;		// 0.01248 mV-1
	double Vhclca = -100.0;			// -100 mV
	double kvclca = 0.0156;			// 0.0156 mV-1
	double EC50 = ECatZero*exp(-kcclca*Vm);

	double Dcclca = 40.0*250.0e-3;		// um^2 s-1
	double Dmclca = 40.0*75.0e-3;		// um^2 s-1
	double Kmclca = 10.0*0.1*1.0e-3;			// 1 uM
	double Bmclca = 50.0e-3;		// 50 uM
	double rclca = 30.0*0.05;			// 0.05 um
	double m_phi = 3.141592653;

	double dCai = CaPMU;
	double C2clca = Dcclca*dCai - Dmclca*Bmclca*Kmclca / (Kmclca + dCai);

	double nsoc = 20.0;			// 250
	double sigma = JSOC*Vl_cyt / nsoc;

	//sigma = -IVDDR / mF / 250.0;	
	double term1 = -Dcclca*Kmclca + sigma / 2.0 / m_phi / rclca + C2clca;
	//double term11 = -Dcclca*Kmclca + C2clca;
	double term2 = 4.0*Dcclca*Dmclca*Bmclca*Kmclca;
	double term3 = Dcclca*Kmclca + sigma / 2.0 / m_phi / rclca + C2clca;
	//double term33 = Dcclca*Kmclca + C2clca;
	CaClCa = (term1 + sqrt(term3*term3 + term2)) / 2.0 / Dcclca;
	//monitorVar1 = (term11 + sqrt(term33*term33 + term2)) / 2.0 / Dcclca;
	//---  Model I -------------------------------------------------------------------------------------//
	double OClCainf = 1.0 / (1.0 + exp((Vhclca - Vm)*kvclca)) / (1.0 + pow(EC50 / CaClCa, 2.0));
	double t1clca = 0.6*0.08163e3*exp(-0.57*CaClCa);			// s unit??
	double t2clca = 0.6*0.07617e3*exp(-0.05374*CaClCa);
	double t3clca = 1.9*70.3e3*exp(0.153*CaClCa);
	double tauClCa = t1clca + t2clca*exp(Vm / t3clca*1000.0);
	OClCa = TwoState(OClCa, (OClCainf - OClCa) / tauClCa, OClCainf / tauClCa, (1.0 - OClCainf) / tauClCa, m_dDt);
	double gAclca = 20.0;			// 20 nS

	// Electromotive force based
	//double Cli = 78.0;				// 78 mM
	//double Cle = 166.0;				// 166 mM
	//double ECl = mRTF*log(ClPMU / Cle);
	//double ECl = -20.2;				// mRTF*log(Cli/Cle);
	//IClCa = SF_ClCa*gAclca*OClCa*(Vm - ECl);			

	// Nernst equation based
	IClCa = 0.263*SF_ClCa*gAclca*OClCa*CF_K;		

}
void CCellPart::SOC_Current()
{	
	//double gSOC = 24.0*1.4*0.2*0.1;				// 0.1 nS, remremrem
	//double KSOC = 200.0e-3;			// 200 uM, remremrem
	//OSOC = pow(KSOC, 8.0) / (pow(KSOC, 8.0) + pow(CaSR*0.021, 8.0));		// remremrem	
	//ISOC = SF_SOC*gSOC*OSOC*CF_Ca;
	//JSOC = -ISOC / mF / Vl_cyt;

	double gSOC = 1.4*0.2*0.1;				// 0.1 nS, remremrem
	double KSOC = 200.0e-3;			// 200 uM, remremrem
	OSOC = pow(KSOC, 8.0) / (pow(KSOC, 8.0) + pow(CaSR*0.021, 8.0));		// remremrem	

	// Electromotive force based
	//double ECa = mRTF / 2.0*log(Cae / CaPMU);
	//ISOC = SF_SOC*gSOC*OSOC*(Vm - ECa);	

	// Nernst equation based
	ISOC = SF_SOC*gSOC*OSOC*20.0*CF_CaPMU;
	
	JSOC = -ISOC / mF / Vl_cyt;	
}
void CCellPart::NaLCN_Current()
{		
	// adjusted for repetitive firing, remrem
	PNaLCN = SF_NaLCN*0.06*1.5*1.5;

	INaLCN_Na = PNaLCN * CF_Na * (1.1);
	INaLCN_K = 0.9 * PNaLCN * CF_K * (1.1);
	INaLCN = INaLCN_Na + INaLCN_K;
}
void CCellPart::Cab_Current()
{		
	double PCab = 15 * 0.04 * SF_Cab;			// modified to increase the amplitude of low frequency stimulation-induced Ca2+ transient
	ICab = PCab * CF_Ca;
}
void CCellPart::Clb_Current()
{
	double PClb =  0.1*SF_Clb;			// 
	IClb = PClb*CF_Cl;
}
void CCellPart::AutonomousInward_Current()
{	
	/*------------ Autonomous inward Channel ------------------------------------*/
	double KmCai = 0.003;	// 0.003 (original value)		
	double BasalFraction = 0.0017;	// 0.0017 (original value)
	double PAIO = BasalFraction + (1.0 - BasalFraction) / (1.0 + KmCai / Cai);
	double PAI = SF_AutonomousInward*15.61;
	IAutonomousInward_Na = 0.4*PAI*PAIO*CF_Na;
	IAutonomousInward_K = 0.36*PAI*PAIO*CF_K;
	IAutonomousInward_Ca = PAI*PAIO*CF_Ca;
	IAutonomousInward = IAutonomousInward_Na + IAutonomousInward_K + IAutonomousInward_Ca;
}
double CCellPart::K1_Current(double Ioni)
{	
	/*---------- Source ---------------------------------------------------------------------------------------------------------//
	Biophysical J 81(6) 3046. PMID: 11720973
	Sandeep V. Pandit, Robert B Clark, Wayne R Giles, and Semahat S. Demir. A mathematical model of action
	potential heterogeneity in adult rat left ventricular myocytes.
	--------------------------------------------------------------------------------------------------------------------------------*/
	//-------------- Current Amplitude ------------------------------------------------//
	double GK1;
	double EK = mRTF * log(Ke / Ioni);			//REVERSAL POTENTIAL FOR K

	// Sandeep et al.
	GK1 = 1.0*0.3*SF_K1 * 60;				// 0.024 microsiemen, modified on 20110807
	double dOut = GK1 * (Vm - EK - 1.73) / (1.0 + exp(1.613*mFRT*(Vm - EK - 1.73)))*pow(Ke / 5.4, 0.4)
		+ (48.0 / (exp((Vm + 37.0) / 25.0) + exp(-(Vm + 37.0) / 25.0)) + 10.0)*(0.0001 / (1.0 + exp(-(Vm - EK - 76.77) / 17.0)));	
	return dOut;
}

void CCellPart::Kr_Current()		// MCKay Model
{		
	double ssAct = 0.1215 + 0.85599 / (1.0 + exp(-(Vm + 24.91989) / 7.32374));
	double tauAct = 1.0 / (0.01356*exp(-Vm / 59.46462) + 93.16276*exp(Vm / 1.60905)) + 1.0 / (0.05504*exp(-Vm / 15.9844) + 0.05177*exp(Vm / 39.0501));	// from my fitting to original current

	double ssInact = 0.23363 + 0.74654 / (1.0 + exp((Vm + 19.28558) / 11.33559));
	double tauInact = 1.0 / (0.01451*exp(-Vm / 0.63999) + 17.22878*exp(Vm / 138.03218)) + 1.0 / (0.26807*exp(-Vm / 75.95043) + 1.39618E-4*exp(Vm / 8.61787));	// from my fitting to original current

	//----------- Integrate ------------------------------------------------//
	YKr1 = TwoState(YKr1, (ssAct - YKr1) / tauAct, ssAct / tauAct, (1.0 - ssAct) / tauAct, m_dDt);
	YKr2 = TwoState(YKr2, (ssInact - YKr2) / tauInact, ssInact / tauInact, (1.0 - ssInact) / tauInact, m_dDt);

	double GKr;
	double EK;

	EK = mRTF * log(Ke / Ki);			//REVERSAL POTENTIAL FOR K	remrem

	GKr = SF_Kr*12.579*pow(Ke / 5.4, 0.2);	// 10.4775 for first model, 12.579 for second model	
	IKr = GKr*(Vm - EK)*YKr1*YKr2;
}


void CCellPart::Kto_Current()
{	
	/*---------- Source ---------------------------------------------------------------------------------------------------------//
	Biophysical J 81(6) 3046.
	Sandeep V. Pandit, Robert B Clark, Wayne R Giles, and Semahat S. Demir. A mathematical model of action
	potential heterogeneity in adult rat left ventricular myocytes.
	--------------------------------------------------------------------------------------------------------------------------------*/
	double KtoActCurve = 1.0 / (1.0 + exp(-(Vm + 1.5) / 25.81));		//----- Nie and Meng, 2005 (16168948). Arch Biochem & Biophys 442, 191
	double KtoInactCurve = 1.0 / (1.0 + exp((Vm + 59.43) / 11.90));	//----- Nie and Meng, 2005 (16168948). Arch Biochem & Biophys 442, 191
	double KtoInact2Curve = KtoInactCurve;						//----- Nie and Meng, 2005 (16168948). Arch Biochem & Biophys 442, 191
	double KtoActTau = 12.3*exp(-(Vm + 70)*(Vm + 70) / 53.4 / 53.4) + 1.15;		// Apkon & Nerbonne, 1991

	double KtoInactTau = 195.8*exp(-(Vm + 36)*(Vm + 36) / 18.4 / 18.4) + 8.7*exp(-(Vm - 57.3)*(Vm - 57.3) / 23 / 23) + 26.7;		// Wettwer et al., 1993
	double KtoInact2Tau = 2380 * exp(-(Vm + 70)*(Vm + 70) / 27.7 / 27.7) + 36.3; // Shimoni et al., 1995

	double dKtoAct_dt = (KtoActCurve - KtoAct) / KtoActTau;
	double dKtoInact_dt = (KtoInactCurve - KtoInact) / KtoInactTau;
	double dKtoInact2_dt = (KtoInact2Curve - KtoInact2) / KtoInact2Tau;

	//----------- Integrate ------------------------------------------------//
	KtoAct = TwoState(KtoAct, dKtoAct_dt, KtoActCurve / KtoActTau, (1.0 - KtoActCurve) / KtoActTau, m_dDt);
	KtoInact = TwoState(KtoInact, dKtoInact_dt, KtoInactCurve / KtoInactTau, (1.0 - KtoInactCurve) / KtoInactTau, m_dDt);
	KtoInact2 = TwoState(KtoInact2, dKtoInact2_dt, KtoInact2Curve / KtoInact2Tau, (1.0 - KtoInact2Curve) / KtoInact2Tau, m_dDt);

	//---------- Current Amplitude ------------------------------------------//	
	double EK;
	EK = mRTF * log(Ke / Ki);					//REVERSAL POTENTIAL FOR K
	double GK = SF_Kto*17.08;					// 0.035 microsiemen (maximum conductance for IKto), modified on 20110807, endocardiac modification
	IKto_K = GK*KtoAct*(0.886*KtoInact + 0.114*KtoInact2)*(Vm - EK);
	IKto_Na = 0;
	IKto = IKto_K;
}


sixelement CCellPart::SK_Current(sixelement dPoSK, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dCai)
{
	/*------------------ IKCa -----------------------------------------------------
	<large conductance Ca-dependent K channel>
	Bolton TB, Imaizumi Y (1996) Spontaneous transient outward currents in smooth
	muscle cells. Cell Calcium 20(2):141-152
	1. BK channel density: 1-5 channels per square micrometer
	2. unitary BK channel current: 1.6 pA
	3. half activation Ca2+: 500 nM (Po of single channel being 0.5)
	4. Duration: 100 - 200 msec

	Herrera GM, Heppner TJ, Nelson MT (2001) Voltage dependence of the coupling of
	Ca2+ sparks to BKCa channels in urinary bladder smooth muscle. AJP 280, C481-C490
	1. coupled to calcium release events from ryanodine receptor
	2. In rat cerebral arterial smooth muscles, Ca++ spark (10-100 micromolar nearby
	BKCa channels) increases the open probability of BKCa channels
	10,000 - 1,000,000 fold
	3. average frequency of BKCa transients: 1.23 Hz, average amplitude: 43.2 pA
	4. simultaneous recording of Ca2+ spark and BKCa current (p 486)
	Duration: Ca2+ spark (96-112 msec), BKCa (16-30 msec)
	5. channel activity: 10^6 fold increase from -50 mV to -20 mV
	6. free Ca2+ increase in the vicinity of BKCa: about 30 fold
	7. Ca2+ sparks increases [Ca2+] from 250 nM at rest to 7.5 micromole at -20 mV

	Frey BW, Carl A, Publicover NG (1998) Charybdoxin block of Ca2+-activated K+
	channels in colonic muscle depends on membrane potential dynamics. AJP C673-680

	Koh SD, Dick GM, Sanders KM (1997) Small-conductance Ca2+-dependent K+ channels
	activated by ATP in murine colonic smooth muscle. AJP 273:C2010-2021
	1. Single channel conductance: 250 pS (maximum Po: 0.12)

	<small conductance Ca-dependent K channel>
	Fujita et al (2001) Expression of Ca-activated K channels, SK3, in the
	interstitial cells of Cajal in the gastrointestinal tract, AJP 281, C1727-C1733
	1. Localized on processes of ICC
	2. SK4: human stomach, small intestine, colon, rat colon, SK3: rat intestine

	Koh SD, Dick GM, Sanders KM (1997) Small-conductance Ca2+-dependent K+ channels
	activated by ATP in murine colonic smooth muscle. AJP 273:C2010-2021
	1. unitary conductance = 0.35 pA at -60 mV (5.3 pS) with 140K/140K, 2.9 pS with 5K(pipette)/140K(bath)
	2. maximum Po = 0.26 + 0.07, at that time mean open time = 29 msec
	3. Ca2+ dependence I: 2 msec open time at 1 nM Ca2+, 53 msec at 1 micromole Ca2+
	4. Ca2+ dependence II: EC50 of NPo = 0.49 micromole
	5. Small amounts of free Ca2+ may preferentially cause activation of SK channels.
	Higher concentrations of Ca2+, perhaps supplemented by entry of Ca+ or more
	potent stimulation of Ca2+ release, may facilitate nonselective cation chanels
	and activate Ca2+-dependent Cl- conductances and BK channels
	------------------------------------------------------------------------------*/
	
	// Calculate Ca near SK channel
	double ECatZero = 1.39e-3;		// 1.39 uM	
	double kcSK = 0.44*0.01248;		// 0.01248 mV-1
	double VhSK = 20.0;				// -100 mV
	double kvSK = 0.0156;			// 0.0156 mV-1
	double EC50 = ECatZero*exp(-kcSK*Vm);

	double DcSK = 40.0*250.0e-3;		// um^2 s-1
	double DmSK = 40.0*75.0e-3;		// um^2 s-1
	double KmSK = 10.0*0.1*1.0e-3;			// 1 uM
	double BmSK = 50.0e-3;		// 50 uM
	double rSK = 30.0*0.05;			// 0.05 um
	double m_phi = 3.141592653;	
	double C2SK = DcSK*dCai - DmSK*BmSK*KmSK / (KmSK + dCai);

	double nsoc = 20.0;			// 250
	double sigma = JSOC*Vl_cyt / nsoc;
		
	double term1 = -DcSK*KmSK + sigma / 2.0 / m_phi / rSK + C2SK;	
	double term2 = 4.0*DcSK*DmSK*BmSK*KmSK;
	double term3 = DcSK*KmSK + sigma / 2.0 / m_phi / rSK + C2SK;	
	double CaSK = (term1 + sqrt(term3*term3 + term2)) / 2.0 / DcSK;
	
	//----------- State transition ---------------------------------------------------------//	
	double slopeK = -17.0;								// -17 mV, slope factor
	double hillCoeff = 2.0;								// 2, Hill coefficient
	double setCa = 2.5e-2;			// 최근값 = 1.5e-2, 0.5e-2					// 2.5 uM, Po is half-maximal Po in Vm=0 with this Ca2+ concentration
	double setCaTau = 1.0;								// 1 uM, Cao in tau equation, remremrem
	double hillCoeffTau = 1.25;							// 1.25, Hill Coefficient in tau equation															
	if (skCaSource == 1) { CaSK = CaPMU; setCa = skKd; hillCoeff = skHill; }	// [modified PDGFa] SK senses submembrane (release) Ca2+
	double VmSK = skVoltIndep ? 0.0 : Vm;		// [modified PDGFa]
	double ssActSK = 1.0 / (1.0 + exp(VmSK / slopeK - hillCoeff*log(CaSK / setCa)));			    //steady-state open of activation		
	double tauActSK = 1000.0 / (1380.1*exp(-VmSK / 34.9 - hillCoeffTau*log(CaSK / setCaTau)) + 35.2*exp(VmSK / 79.5 + hillCoeffTau*log(CaSK / setCaTau)));
	/*kc2c1SK = 0.95*(1000.0*CaSK)*(1000.0*CaSK)*0.001*exp((Vm - 20) / 30.0);
	kc1c2SK = 1270.0*0.001;
	kc1oSK = 2905.0*0.001;
	koc1SK = 180.0*0.001;*/
	kc2c1SK = 0;
	kc1c2SK = 0;
	kc1oSK = ssActSK / tauActSK;
	koc1SK = (1.0 - ssActSK) / tauActSK;

	// [modified PDGFa] deterministic (mean-field) SK gate. The stochastic scheme below lets a channel change state at
	// most once per step; when the rates exceed 1/dt (tauActSK can be ~1e-4 ms) every channel flips each step and the
	// open fraction locks near 0.5 regardless of Ca2+. Here Po relaxes to ssActSK with an SK3-like time constant.
	if (skDeterministic)
	{
		poSKdet += (ssActSK - poSKdet) * (1.0 - exp(-m_dDt / tauSKdet));
		double ErevDet = skUseEK ? mRTF * log(Ke / Ki) : -70.0;
		dPoSK.a = 0.0;
		dPoSK.b = 1.0 - poSKdet;
		dPoSK.c = poSKdet;
		dPoSK.e = SF_SK * SF_NoSK * 2.9e-3 * poSKdet * (Vm - ErevDet);
		return dPoSK;
	}

	//----------- Stochastic 3-state Model 2016. 10. 27 ----------------------------------------//					
	int NN = 0;
	int tState = 1;
	int unit0 = 0;	int unit1 = 0;	int unit2 = 0;
	for (NN = 0; NN<iChannelNumber; NN++)
	{
		tState = iState[NN];
		if (timeX >= dStayingTimePredict[NN] + dTimePredict[NN])
		{
			// [headless] random numbers are drawn only when a transition occurs
			// (the original drew two rand() per channel per step and discarded
			// them otherwise; distributionally identical, ~10x faster).
			// dSwitch is not used by the SK scheme and is not drawn.
			double dStay = fabs(log(double(Rand() + 0.9) / double(RAND_MAX + 0.9)));  // it is same as abs(log(1.0/double(rand())*double(RAND_MAX)))
			switch (tState)
			{
			case 0:				// c2 (close 2)
				iState[NN] = 1;						// in case, c1
				dStayingTimePredict[NN] = dStay / (kc1c2SK + kc1oSK);
				dTimePredict[NN] = timeX;
				break;
			case 1:				// c1 (close 1)
				//if (dSwitch > kc1oSK / (kc1c2SK + kc1oSK))		// in case, it goes to c2 state
				//{
				//	iState[NN] = 0;
				//	dStayingTimePredict[NN] = dStay / (kc2c1SK);
				//	dTimePredict[NN] = timeX;
				//}
				//else											// in case, it goes to open
				//{
					iState[NN] = 2;
					dStayingTimePredict[NN] = dStay / (koc1SK);
					dTimePredict[NN] = timeX;
				//}
				break;			
			case 2:				// open
				iState[NN] = 1;						// in case, c1
				dStayingTimePredict[NN] = dStay / (kc1c2SK + kc1oSK);
				dTimePredict[NN] = timeX;
				break;
			}
		}

		tState = iState[NN];
		switch (tState)
		{
		case 0:				// close 2
			unit0++;
			break;
		case 1:				// close 1
			unit1++;
			break;
		case 2:				// open
			unit2++;
			break;		
		}
	}
	// [headless] removed unused 'unitAll' (triplet-open count) loop for SK: result was never used

	dPoSK.a = double(unit0) / double(iChannelNumber);
	dPoSK.b = double(unit1) / double(iChannelNumber);
	dPoSK.c = double(unit2) / double(iChannelNumber);	
		
	// for 5.4 mM K+ concentration	
	//double GISK = 0.05*SF_SK*1.0;		// 최근값 = 20.0, 0.5
	double GISK = SF_SK;		// 최근값 = 20.0, 0.5
	//double dSingleChannelConductance = 234.0*0.001;			// 234 pS for equimolar 150mM K+, 100~150 pS for physiological solution
	double dSingleChannelConductance = 2.9*0.001;				// unitary conductance = 0.35 pA at - 60 mV(5.3 pS) with 140K / 140K, 2.9 pS with 5K(pipette) / 140K(bath)
	double Erev = skUseEK ? mRTF * log(Ke / Ki) : -70.0;	// [modified PDGFa]
	dPoSK.e = GISK*SF_NoSK*dSingleChannelConductance*dPoSK.c*(Vm - Erev);		// remrem

	return dPoSK;
}
sixelement CCellPart::initializeSTOCHASTICSK(sixelement dPoSK, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dCai)
{
	//------ SK channel (stochastic model) --------------//				
	// [headless] srand((unsigned)time(NULL)) removed: per-object RNG (Reseed)
	int NN = 0;

	// calculate initial stay time of each state
	int unitA = 0;
	int unitB = 0;
	int unitC = 0;	

	//----------- State transition ---------------------------------------------------------//	
	kc2c1SK = 1.2; // 0.95*(1000.0*dCai)*(1000.0*dCai)*0.001*exp((Vm - 20) / 30.0);
	kc1c2SK = 1270.0*0.001;
	kc1oSK = 2905.0*0.001;
	koc1SK = 180.0*0.001;

	for (NN = 0; NN<iChannelNumber; NN++)
	{
		if (NN<int(double(iChannelNumber)*dPoSK.a))		// close 2 state
		{
			iState[NN] = 0;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / (kc2c1SK);
			unitA++;
		}
		else if (NN<int(double(iChannelNumber)*(dPoSK.a + dPoSK.b)))		// close  1 state
		{
			iState[NN] = 1;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / (kc1c2SK + kc1oSK);
			unitB++;
		}
		else
		{
			iState[NN] = 2;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / (koc1SK);
			unitC++;
		}
	}
	dPoSK.a = double(unitA) / double(iChannelNumber);
	dPoSK.b = double(unitB) / double(iChannelNumber);	
	dPoSK.c = 1.0 - dPoSK.a - dPoSK.b;
		
	return dPoSK;
}
void CCellPart::initializeCurrentTimeSTOCHASTICSK(int iChannelNumber, double* dTimePredict)
{
	int NN;
	for (NN = 0; NN<iChannelNumber; NN++)
	{
		dTimePredict[NN] = 0;
	}
}
void CCellPart::NaCa_Current()
{	
	// rate constants array
	double KNae = 87.5;				// KmValue of external Na (mM)
	double KCae =  1.38;				// KmValue of external Ca (mM)
	double KNai = 0.5*8.75;
	double KCai = 0.00138;

	double kCa_f = 1002.0e-3;     // Calcium side E1-E2 inside rate
	double kCa_b = 1002.0e-3;    // Calcium side E1-E2 outside rate
	double kNa_f = exp(0.32*Vm*mFRT);				// Na side E1-E2 outside rate
	double kNa_b = exp((0.32 - 1.0)*Vm*mFRT);		// Na side E1-E2 inside rate
	
	double kion_on = 1.0e2;		// forward binding rate per mM ion ( 10^8 liter/mol/s = 10^2 mM-1 ms-1)

	double kNae_off = kion_on*KNae;					// external Na-off rate, ms-1 
	double kCae_off = kion_on*KCae;					// external Ca-off rate
	double kNai_off = kion_on*KNai;			// cytosol Na-off rate
	double kCai_off = kion_on*KCai;			// cytosol Ca-off rate

	double dCai = CaPMU;
	double dNai = NaPMU;

	using Eigen::MatrixXf;
	MatrixXf mE(10, 10);
	mE << -kNa_b - kNae_off, kion_on*Nae, 0, 0, 0, 0, 0, 0, 0, kNa_f,
		kNae_off, -kion_on*Nae - kNae_off, kion_on*Nae, 0, 0, 0, 0, 0, 0, 0,
		0, kNae_off, -kNae_off - kion_on*Nae, kion_on*Nae, 0, 0, 0, 0, 0, 0,
		0, 0, kNae_off, -kion_on*Nae - kion_on*Cae, kCae_off, 0, 0, 0, 0, 0,
		0, 0, 0, kion_on*Cae, -kCae_off - kCa_f, kCa_b, 0, 0, 0, 0,
		0, 0, 0, 0, kCa_f, -kCa_b - kCai_off, kion_on*dCai, 0, 0, 0,
		0, 0, 0, 0, 0, kCai_off, -kion_on*dNai - kion_on*dCai, kNai_off, 0, 0,
		0, 0, 0, 0, 0, 0, kion_on*dNai, -kion_on*dNai - kNai_off, kNai_off, 0,
		0, 0, 0, 0, 0, 0, 0, kion_on*dNai, -kion_on*dNai - kNai_off, kNai_off,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
	
	MatrixXf rE(10, 1);
	rE << 0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		Xtotal;

	MatrixXf aE(10, 1);
	aE = mE.inverse()*rE;

	X2Na3 = aE(0, 0);
	X2Na2 = aE(1, 0);
	X2Na1 = aE(2, 0);
	X2 = aE(3, 0);
	X2Ca = aE(4, 0);
	X1Ca = aE(5, 0);
	X1 = aE(6, 0);
	X1Na1 = aE(7, 0);
	X1Na2 = aE(8, 0);
	X1Na3 = Xtotal - (X2Na3 + X2Na2 + X2Na1 + X2 + X2Ca + X1Ca + X1 + X1Na1 + X1Na2);
	
	double v1 = kNae_off*X2Na3;
	double vn1 = kion_on*Nae*X2Na2;
	double v2 = kNae_off*X2Na2;
	double vn2 = kion_on*Nae*X2Na1;
	double v3 = kNae_off*X2Na1;
	double vn3 = kion_on*Nae*X2;
	double v4 = kion_on*Cae*X2;
	double vn4 = kCae_off*X2Ca;
	double v5 = kCa_f*X2Ca;
	double vn5 = kCa_b*X1Ca;
	double v6 = kCai_off*X1Ca;
	double vn6 = kion_on*dCai*X1;
	double v7 = kion_on*dNai*X1;
	double vn7 = kNai_off*X1Na1;
	double v8 = kion_on*dNai*X1Na1;
	double vn8 = kNai_off*X1Na2;
	double v9 = kion_on*dNai*X1Na2;
	double vn9 = kNai_off*X1Na3;

	double JX2_X1Ca = v4*v5 / (v5 + vn4);
	double JX1Ca_X2 = vn4*vn5 / (v5 + vn4);
	double JX2_X1 = JX2_X1Ca*v6 / (v6 + JX1Ca_X2);				// Ca influx
	double JX1_X2 = JX1Ca_X2*vn6 / (v6 + JX1Ca_X2);				// Ca efflux
	
	double PNaCa = SF_NaCa*8.458;

	//if (timeX > 2400)
	//	PNaCa *= 5.0;

	YNaCa = YNaCa;  											// Ca influx is negative of Na influx
	INaCa = PNaCa*(JX2_X1 - JX1_X2);

	/*if (timeX > 1000 && INaCa > 0)
		INaCa = 0;*/

	/*testVar1 = PNaCa*(JX2_X1 - JX1_X2);
	testVar2 = PNaCa*(k1*pE1Na*YNaCa - k2*pE2Na*(1.0 - YNaCa)); 
	INaCa = testVar2;*/
}
double CCellPart::NaK_Current(double dNai, double dKi)
{	
	double ampINaK = 25.63 * 0.8;         //32.033232000000005 * 0.8
	double k1p_NaK = 0.72;
	double k1m_NaK = 0.08;
	double k2p_NaK = 0.08;
	double k2m_NaK = 0.008;
	double k3p_NaK = 4;
	double k3m_NaK = 8000;
	double k4p_NaK = 0.3;
	double k4m_NaK = 0.2;

	double KdKe0_NaK = 0.8;             //dissociation const mM for [K]o at 0 mV
	double KdKi0_NaK = 18.8;            //dissociation const mM for [K]i at 0 mV
	double KdNae0_NaK = 26.8;           //dissociation const mM for [Na]o at 0 mV
	double KdNai0_NaK = 5;              //dissociation const mM for [Na]i at 0 mV
	double KdMgATP_NaK = 0.6;

	double stoiNa_NaK = 3;
	double stoiK_NaK = -2;
	double dV_Nae_NaK = 0.44;           //functional depth of the binding sites
	double dV_Nai_NaK = -0.14;          //functional depth of the binding sites
	double dV_Ke_NaK = 0.23;            //functional depth of the binding sites
	double dV_Ki_NaK = -0.14;           //functional depth of the bin	
	
	double ATPMg = 6.0;
	double kd_MgATP = 0.6;				// mM
	double Pifree_cyt = 1.9;			// mM concentration
	double H_cyt = 0.0001;				// mM/ms,  rate of cytosol H concentration
		
	
	double MgADP = 0.01;
	double fVm_NaK = mFRT * Vm;													 //VF/R/T
	double KdNae_NaK = KdNae0_NaK * exp(dV_Nae_NaK*fVm_NaK);					 //relative Kd 
	double KdNai_NaK = KdNai0_NaK * exp(dV_Nai_NaK*fVm_NaK);					 //relative Kd 
	double KdKe_NaK = KdKe0_NaK * exp(dV_Ke_NaK*fVm_NaK);						 //relative Kd 
	double KdKi_NaK = KdKi0_NaK * exp(dV_Ki_NaK*fVm_NaK);						 //relative Kd 
	double Nai_ = dNai / KdNai_NaK;												 //relative concentration referring to the dissociation constant 
	double Nae_ = Nae / KdNae_NaK;												 //relative concentration referring to the dissociation constant
	double Ki_ = dKi / KdKi_NaK;							                         //relative concentration referring to the dissociation constant
	double Ke_ = Ke / KdKe_NaK;													 //relative concentration referring to the dissociation constant
	double MgATP_ = ATPMg / kd_MgATP;							                 //relative concentration of ATP referring to the dissociation constant
																					 // transition rates in the forward direction
	double a1_plus_ = (k1p_NaK * pow(dNai, 3.0)) / (pow(1 + Nai_, 3) + pow(1 + Ki_, 2) - 1);
	double a2_plus_ = k2p_NaK;
	double a3_plus_ = k3p_NaK * pow(Ke_, 2) / (pow(1 + Nae_, 3) + pow(1 + Ke_, 2) - 1);
	double a4_plus_ = k4p_NaK * MgATP_ / (1 + MgATP_);

	// transition rates in the reverse direction
	double a1_minus_ = k1m_NaK * MgADP;
	double a2_minus_ = k2m_NaK * pow(Nae_, 3) / (pow(1 + Nae_, 3) + pow(1 + Ke_, 2) - 1);
	double a3_minus_ = k3m_NaK * Pifree_cyt * H_cyt / (1 + MgATP_);
	double a4_minus_ = k4m_NaK * pow(dKi, 2) / (pow(1 + Nai_, 3) + pow(1 + Ki_, 2) - 1);

	double vcyc_NaK;        //any step when calculating the steady state distribution

	double denomi;
	double numer;
	denomi = (a1_minus_ + a1_plus_) * a2_minus_ * a3_minus_ + a1_plus_ * a2_plus_ * (a3_plus_ + a3_minus_) + a2_plus_ * a3_plus_ * (a4_plus_ + a4_minus_) + (a2_plus_ + a2_minus_) * a3_minus_ * a4_minus_ + (a1_minus_ + a1_plus_) * a3_plus_ * a4_plus_ + a1_minus_ * (a3_plus_ + a3_minus_) * a4_minus_ + a1_plus_ * (a2_plus_ + a2_minus_) * a4_plus_ + a1_minus_ * a2_minus_ * (a4_plus_ + a4_minus_);
	numer = a1_plus_ * a2_plus_ * a3_plus_ * a4_plus_ - a1_minus_ * a2_minus_ * a3_minus_ * a4_minus_;

	vcyc_NaK = (numer / denomi);

	double dOut = ampINaK * vcyc_NaK * 165.0;					// current  density (uA/cm2)

	return dOut;	

	//INaK_Na = stoiNa_NaK * INaK;						// calculate ion flux   '3 * INaK
	//INaK_K = stoiK_NaK * INaK;
	
	// Original NaK model
	//double PNaK = 0.7*SF_NaK*39.5;
	//double KmNai = 3.5;
	//double Km2Nai = 18.0;
	//double KmKe = 0.27;
	//INaK = 0.7*PNaK / (1.0 + pow(KmNai / Nai, 1.36))*(1.0 - ((Vm + 50.0) / 250.0)*((Vm + 50.0) / 250.0)) / (1.0 + KmKe / Ke)
	//	+ 7.0*PNaK / (1.0 + pow(Km2Nai / Nai, 3.0))*(1.0 - ((Vm + 50.0) / 250.0)*((Vm + 50.0) / 250.0)) / (1.0 + KmKe / Ke);
}

// Biophysical Journal   Volume 81   December 2001   3029–3051
// A Mathematical Model of Action Potential Heterogeneity in Adult Rat Left Ventricular Myocytes
// Sandeep V. Pandit, Robert B. Clark, Wayne R. Giles, and Semahat S. Demir
void CCellPart::PMCA()
{
	// original ICC model
	double PPMCA = 927.4;
	double KmCai = 0.008;
	IPMCA = SF_PMCA * PPMCA / (1 + (KmCai / Cai)*(KmCai / Cai));
}
void CCellPart::NKCC()
{	
	//Nae = 140.0;
	//NaPMU = 5.0;
	//Ke = 5.0;
	//KPMU = 130.0;
	//Cle = 140.0;
	//ClPMU = 80.0;

	/*A quantitative description of the Na-K-2Cl cotransporter  and its conformity to experimental data
	Benjamin BA and Johnson  EA, AJP 1997:273:F473-F482*/
	double KNa = 15.8e-3;   //HeLa cell
	double KCl = 216.9e-3;
	double KK = 9.22e-3;
	double kfull_f = 1002.0e-3;
	double kfull_b = 2255.0e-3;
	double kempty_f = 37767.0e-3;	
	double kion_on = 1.0e2;
	double kNa_off = kion_on / KNa;
	double kK_off = kion_on / KK;
	double kCl_off = kion_on / KCl;
	double kempty_b = KNa*KCl*KK*KCl*kfull_f*kempty_f / KNa / KCl / KK / KCl / kfull_b;				// 16783.0e-3

	double dNai = NaPMU;
	double dKi = KPMU;
	double dCli = ClPMU;

	using Eigen::MatrixXf;
	MatrixXf mE(10, 10);
	mE << -kempty_b - kion_on*Nae, kNa_off, 0, 0, 0, 0, 0, 0, 0, kempty_f,
		kion_on*Nae, -kNa_off - kion_on*Cle, kCl_off, 0, 0, 0, 0, 0, 0, 0,
		0, kion_on*Cle, -kCl_off - kion_on*Ke, kK_off, 0, 0, 0, 0, 0, 0,
		0, 0, kion_on*Ke, -kK_off - kion_on*Cle, kCl_off, 0, 0, 0, 0, 0,
		0, 0, 0, kion_on*Cle, -kCl_off - kfull_f, kfull_b, 0, 0, 0, 0,
		0, 0, 0, 0, kfull_f, -kfull_b - kNa_off, kion_on*dNai, 0, 0, 0,
		0, 0, 0, 0, 0, kNa_off, -kion_on*dNai - kCl_off, kion_on*dCli, 0, 0,
		0, 0, 0, 0, 0, 0, kCl_off, -kion_on*dCli - kK_off, kion_on*dKi, 0,
		0, 0, 0, 0, 0, 0, 0, kK_off, -kion_on*dKi - kCl_off, kion_on*dCli,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1;		

	MatrixXf rE(10, 1);
	rE << 0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		Etotal;

	MatrixXf aE(10, 1);
	aE = mE.inverse()*rE;

	E1 = aE(0, 0);
	E1Na = aE(1, 0);
	E1NaCl = aE(2, 0);
	E1NaClK = aE(3, 0);
	E1NaClKCl = aE(4, 0);
	E2NaClKCl = aE(5, 0);
	E2ClKCl = aE(6, 0);
	E2KCl = aE(7, 0);
	E2Cl = aE(8, 0);
	E2 = Etotal - (E1 + E1Na + E1NaCl + E1NaClK + E1NaClKCl + E2NaClKCl + E2ClKCl + E2KCl + E2Cl);

	double v1 = kion_on*Nae*E1;
	double vn1 = kion_on / KNa*E1Na;
	double v2 = kion_on*Cle*E1Na;
	double vn2 = kion_on / KCl*E1NaCl;
	double v3 = kion_on*Ke*E1NaCl;
	double vn3 = kion_on / KK*E1NaClK;
	double v4 = kion_on*Cle*E1NaClK;
	double vn4 = kion_on / KCl*E1NaClKCl;
	double v5 = kfull_f*E1NaClKCl;
	double vn5 = kfull_b*E2NaClKCl;
	double v6 = kion_on / KNa*E2NaClKCl;
	double vn6 = kion_on*dNai*E2ClKCl;
	double v7 = kion_on / KCl*E2ClKCl;
	double vn7 = kion_on*dCli*E2KCl;
	double v8 = kion_on / KK*E2KCl;
	double vn8 = kion_on*dKi*E2Cl;
	double v9 = kion_on / KCl*E2Cl;
	double vn9 = kion_on*dCli*E2;

	double PNKCC = SF_NKCC*1.6*1500.0*0.07e-3;
	/*if (timeX > 20000)
		PNKCC = 0;*/
	double JE1_E1NaCl = v1*v2 / (v2 + vn1);
	double JE1NaCl_E1 = vn2*vn1 / (v2 + vn1);
	double JE1_E1NaClK = JE1_E1NaCl*v3 / (v3 + JE1NaCl_E1);
	double JE1NaClK_E1 = JE1NaCl_E1*vn3 / (v3 + JE1NaCl_E1);
	double JE1_E1NaClKCl = JE1_E1NaClK*v4 / (v4 + JE1NaClK_E1);
	double JE1NaClKCl_E1 = JE1NaClK_E1*vn4 / (v4 + JE1NaClK_E1);
	double JE1_E2NaClKCl = JE1_E1NaClKCl*v5 / (v5 + JE1NaClKCl_E1);
	double JE2NaClKCl_E1 = JE1NaClKCl_E1*vn5 / (v5 + JE1NaClKCl_E1);
	double JNKCC_Na_influx = JE1_E2NaClKCl*v6 / (v6 + JE2NaClKCl_E1);
	double JNKCC_Na_efflux = JE2NaClKCl_E1*vn6 / (v6 + JE2NaClKCl_E1);
	JNKCC_Na = PNKCC*(JNKCC_Na_influx - JNKCC_Na_efflux);
	/*if (timeX <5000)
		JNKCC_Na = 0.001*PNKCC*(JNKCC_Na_influx - JNKCC_Na_efflux);*/
	
	// Following codes are not neccessary. JNKCC_Na = JNKCC_K, JNKCC_Cl = JNKCC_Na + JNKCC_K
	/*double JE1NaCl_E1NaClKCl = v3*v4 / (v4 + vn3);
	double JE1NaClKCl_E1NaCl = vn3*vn4 / (v4 + vn3);
	double JE1NaCl_E2NaClKCl = JE1NaCl_E1NaClKCl*v5 / (v5 + JE1NaClKCl_E1NaCl);
	double JE2NaClKCl_E1NaCl = JE1NaClKCl_E1NaCl*vn5 / (v5 + JE1NaClKCl_E1NaCl);
	double JE1NaCl_E2ClKCl = JE1NaCl_E2NaClKCl*v6 / (v6 + JE2NaClKCl_E1NaCl);
	double JE2ClKCl_E1NaCl = JE2NaClKCl_E1NaCl*vn6 / (v6 + JE2NaClKCl_E1NaCl);
	double JE1NaCl_E2KCl = JE1NaCl_E2ClKCl*v7 / (v7 + JE2ClKCl_E1NaCl);
	double JE2KCl_E1NaCl = JE2ClKCl_E1NaCl*vn7 / (v7 + JE2ClKCl_E1NaCl);
	double JNKCC_K_influx = JE1NaCl_E2KCl*v8 / (v8 + JE2KCl_E1NaCl);
	double JNKCC_K_efflux = JE2KCl_E1NaCl*vn8 / (v8 + JE2KCl_E1NaCl);
	JNKCC_K = PNKCC*(JNKCC_K_influx - JNKCC_K_efflux);*/
	JNKCC_K = JNKCC_Na;
	/*if (timeX <5000)
		JNKCC_K = 0.001*PNKCC*(JNKCC_K_influx - JNKCC_K_efflux);*/

	// Following codes are not neccessary. JNKCC_Na = JNKCC_K, JNKCC_Cl = JNKCC_Na + JNKCC_K
	/*double JE1Na_E1NaClK = v2*v3 / (v3 + vn2);
	double JE1NaClK_E1Na = vn2*vn3 / (v3 + vn2);
	double JE1Na_E1NaClKCl = JE1Na_E1NaClK*v4 / (v4 + JE1NaClK_E1Na);
	double JE1NaClKCl_E1Na = JE1NaClK_E1Na*vn4 / (v4 + JE1NaClK_E1Na);
	double JE1Na_E2NaClKCl = JE1Na_E1NaClKCl*v5 / (v5 + JE1NaClKCl_E1Na);
	double JE2NaClKCl_E1Na = JE1NaClKCl_E1Na*vn5 / (v5 + JE1NaClKCl_E1Na);
	double JE1Na_E2ClKCl = JE1Na_E2NaClKCl*v6 / (v6 + JE2NaClKCl_E1Na);
	double JE2ClKCl_E1Na = JE2NaClKCl_E1Na*vn6 / (v6 + JE2NaClKCl_E1Na);
	double JNKCC_Cl_influx1 = JE1Na_E2ClKCl*v7 / (v7 + JE2ClKCl_E1Na);
	double JNKCC_Cl_efflux1 = JE2ClKCl_E1Na*vn7 / (v7 + JE2ClKCl_E1Na);

	double JE1NaClK_E2NaClKCl = v4*v5 / (v5 + vn4);
	double JE2NaClKCl_E1NaClK = vn4*vn5 / (v5 + vn4);
	double JE1NaClK_E2ClKCl = JE1NaClK_E2NaClKCl*v6 / (v6 + JE2NaClKCl_E1NaClK);
	double JE2ClKCl_E1NaClK = JE2NaClKCl_E1NaClK*vn6 / (v6 + JE2NaClKCl_E1NaClK);
	double JE1NaClK_E2KCl = JE1NaClK_E2ClKCl*v7 / (v7 + JE2ClKCl_E1NaClK);
	double JE2KCl_E1NaClK = JE2ClKCl_E1NaClK*vn7 / (v7 + JE2ClKCl_E1NaClK);
	double JE1NaClK_E2Cl = JE1NaClK_E2KCl*v8 / (v8 + JE2KCl_E1NaClK);
	double JE2Cl_E1NaClK = JE2KCl_E1NaClK*vn8 / (v8 + JE2KCl_E1NaClK);
	double JNKCC_Cl_influx2 = JE1NaClK_E2Cl*v9 / (v9 + JE2Cl_E1NaClK);
	double JNKCC_Cl_efflux2 = JE2Cl_E1NaClK*vn9 / (v9 + JE2Cl_E1NaClK);
	JNKCC_Cl = PNKCC*(JNKCC_Cl_influx1 - JNKCC_Cl_efflux1 + JNKCC_Cl_influx2 - JNKCC_Cl_efflux2);*/
	//JNKCC_Cl = 2.0*PNKCC*(JNKCC_Na_influx - JNKCC_Na_efflux);
	JNKCC_Cl = 2.0*JNKCC_Na;
}

sixelement CCellPart::IP3R_Current(sixelement dPoIP3R, int iChannelNumber, int *iState, double* dStayingTimePredict, double* dTimePredict, double dIP3, double dCai)
{	
	// IP3  rates
	kucIP3R = 400.0*dIP3;
	kcuIP3R = 400.0*0.9434e-3;			// d1: 0.13, d3: 0.9434 uM	
	kcoIP3R = 20.0*dCai;
	kocIP3R = 20.0*0.08234e-3;			// d5: 0.08234 uM
	kobIP3R = 0.2*dCai;
	kboIP3R = 0.2*1.049e-3;				// d2: 1.049, d4: 0.1445 uM		//*/

	//double dCIP3R_dt, dOIP3R_dt, dUIP3R_dt, dBIP3R_dt;
	double PIP3R;		
		
	//----------- (1) Stochastic 4-state Model  ----------------------------------------//				
	int NN = 0;	
	int tState = 1;	
	int unit0 = 0;	int unit1 = 0;	int unit2 = 0;	int unit3 = 0;
	for (NN=0; NN<iChannelNumber; NN++)
	{	
		tState = iState[NN];
		if (timeX >= dStayingTimePredict[NN]+dTimePredict[NN])
		{
			// [headless] draw random numbers only on a transition (see SK_Current)
			double dStay = fabs(log(double(Rand()+0.9)/double(RAND_MAX+0.9)));  // it is same as abs(log(1.0/double(rand())*double(RAND_MAX)))
			double dSwitch = double(Rand()+1)/double(RAND_MAX+1);
			switch (tState)
			{
				case 0:				// close
					if (dSwitch > kcoIP3R/(kcoIP3R+kcuIP3R))		// in case, it goes to unavailable				
					{
						iState[NN] = 2;					
						dStayingTimePredict[NN] = dStay/(kucIP3R);
						dTimePredict[NN] = timeX;							
					}
					else											// in case, it goes to open
					{
						iState[NN] = 1;
						dStayingTimePredict[NN] = dStay/(kocIP3R+kobIP3R);	
						dTimePredict[NN] = timeX;							
					}								
					break;
				case 1:				// open
					if (dSwitch > kocIP3R/(kobIP3R+kocIP3R)) // in case, blocked				
					{
						iState[NN] = 3;					
						dStayingTimePredict[NN] = dStay/(kboIP3R);					
						dTimePredict[NN] = timeX;								
					}
					else									// in case, close
					{
						iState[NN] = 0;
						dStayingTimePredict[NN] = dStay/(kcuIP3R+kcoIP3R);	
						dTimePredict[NN] = timeX;						
					}				
					break;			
				case 2:				// unavailable
					iState[NN] = 0;						// in case, close
					dStayingTimePredict[NN] = dStay/(kcuIP3R+kcoIP3R);								
					dTimePredict[NN] = timeX;
					break;
				case 3:				// blocked
					iState[NN] = 1;						// in case, open
					dStayingTimePredict[NN] = dStay/(kocIP3R+kobIP3R);								
					dTimePredict[NN] = timeX;
					break;
			}		
		}	

		tState = iState[NN];
		switch (tState)
		{
			case 0:				// close
				 unit0++;			
				break;
			case 1:				// open
				 unit1++;		
				break;			
			case 2:				// unavailable
				unit2++;
				break;
			case 3:				// blocked
				 unit3++;
				break;
		}	
	}
	int unitAll = 0;
	for (NN = 0; NN < iChannelNumber - 2; NN = NN + 3)
	{
		if (iState[NN] == 1 && iState[NN + 1] == 1 && iState[NN + 2] == 1)
		{
			unitAll++;
		}
	}

	dPoIP3R.a = double(unit0)/double(iChannelNumber);
	dPoIP3R.b = double(unit3)/double(iChannelNumber);
	dPoIP3R.c = double(unit2)/double(iChannelNumber);
	dPoIP3R.d = double(unit1)/double(iChannelNumber);
	
	powerOIP3R = dPoIP3R.d*dPoIP3R.d*dPoIP3R.d;

	//poweraOIP3R = pow(aOIP3R, 3.0);								// not appropriate method
	poweraOIP3R = 3.0*double(unitAll) / double(iChannelNumber);	// revised method

	PIP3R = 0.2*217.2;	
	dPoIP3R.e = PIP3R*(CaiSR - dCai)*poweraOIP3R;	
		
	return dPoIP3R;
	
	//----------- (2) Colligative 4-state Model  ----------------------------------------//
	//dCIP3R_dt = kucIP3R*aUIP3R + kocIP3R*aOIP3R - (kcuIP3R + kcoIP3R)*aCIP3R;
	//dOIP3R_dt = kcoIP3R*aCIP3R + kboIP3R*aBIP3R - (kocIP3R + kobIP3R)*aOIP3R;
	//dUIP3R_dt = kcuIP3R*aCIP3R - kucIP3R*aUIP3R;	
	//dBIP3R_dt = kobIP3R*aOIP3R - kboIP3R*aBIP3R;		
	//					
	//aCIP3R = Euler(aCIP3R, dCIP3R_dt, m_dDt);
	//aBIP3R = Euler(aBIP3R, dBIP3R_dt, m_dDt);
	//aUIP3R = Euler(aUIP3R, dUIP3R_dt, m_dDt);
	//aOIP3R = 1.0 - aCIP3R - aBIP3R - aUIP3R;					

	//PIP3R = 217.2;
	//IIP3R = PIP3R*(CaiSR - Cai)*pow(aOIP3R, 3.0);				

	//----------- (3) De Young  Model  ----------------------------------------//	
	//double capIP3 = 15.0e-6;
	//double Kd1 = 145.0e-6;
	//double Kd2 = 542.0e-6;
	//double d1 = Kd1 - capIP3;
	//double d2 = 1.049e-3;
	//double d3 = 943.4e-6; // (Kd2 - capIP3)*(1.0e-3 + d2) - d1*d2;
	//double d4 = 144.5e-6;
	//double d5 = 82.34e-6;
	////poweraOIP3R0 = pow(Cai*IP3*d2 / (Cai*IP3 + IP3*d2 + d1*d2 + Cai*d3) / (Cai + d5), 3.0);

	//double a1 = 400.0;		//uM-1 S-1 = mM-1 ms-1
	//double a2 = 0.2;
	//double a3 = 300.0;
	//double a4 = 0.2;
	//double a5 = 20.0;
	//double b1 = 0.13e-3*a1;
	//double b2 = 1.049e-3*a2;
	//double b3 = 943.4e-6*a3;
	//double b4 = 144.5e-6*a4;
	//double b5 = 82.34e-6*a5;

	//double dS000dt = b1*S100 + b4*S001 + b5*S010 - (a4*Cai + a1*IP3 + a5*Cai)*S000;
	//double dS100dt = b2*S101 + a1*IP3*S000 + b5*S110 - (a2*Cai + b1 + a5*Cai)*S100;
	//double dS101dt = a2*Cai*S100 + a3*IP3*S001 + b5*S111 - (b2 + b3 + a5*Cai)*S101;
	//double dS001dt = a4*Cai*S000 + b3*S101 + b5*S011 - (b4 + a3*IP3 + a5*Cai)*S001;
	//
	//double dS010dt = b1*S110 + b4*S011 + a5*Cai*S000 - (a4*Cai + a1*IP3 + b5)*S010;
	//double dS110dt = b2*S111 + a1*IP3*S010 + a5*Cai*S100 - (a2*Cai + b1 + b5)*S110;
	//double dS111dt = a2*Cai*S110 + a3*IP3*S011 + a5*Cai*S101 - (b2 + b3 + b5)*S111;
	//double dS011dt = a4*Cai*S010 + b3*S111 + a5*Cai*S001 - (b4 + a3*IP3 + b5)*S011;
	//
	//S000 = Euler(S000, dS000dt, m_dDt);
	//S100 = Euler(S100, dS100dt, m_dDt);
	//S101 = Euler(S101, dS101dt, m_dDt);
	//S001 = Euler(S001, dS001dt, m_dDt);
	//
	//S010 = Euler(S010, dS010dt, m_dDt);
	//S110 = Euler(S110, dS110dt, m_dDt);
	//S111 = Euler(S111, dS111dt, m_dDt);
	//S011 = Euler(S011, dS011dt, m_dDt);

	//poweraOIP3R = S110*S110*S110;

	//PIP3R = 217.2;		
	//IIP3R = PIP3R*(CaiSR-Cai)*poweraOIP3R;							// stochastic	

	//----------- (4) Original 3-state model  ----------------------------------------//	
	//double kic = 0.0025*pow(IP3/0.002,3.0)*CaSR;
	//double kci = 0.000849;													// 0.000849
	//double KmCaPMU = 0.000057*2.0;					// 0.000057
	//double kco = 0.0147*pow(CaPMU/KmCaPMU,3.0);								//0.147, 12.4*pow(Caf[L],1.0);							// 280000*pow(Caf[L]/0.00006,1.0)*pow(10, -10)
	//double koi = 0.02; 
	//double koc = 0;
	//double kio = 0;

	//// Colligative property
	//double dtC_dt = kic*tI + koc*tO - (kci + kco)*tC;
	//double dtO_dt = kco*tC + kio*tI - (koc + koi)*tO;
	//double dtI_dt = koi*tO + kci*tC - (kio + kic)*tI;
	//
	//twoelement ThreeStateSet;
	//twoelement dThreeStateSet;
	//sixelement RateConstSet;
	//twoelement pThreeStateSet;
	//
	//ThreeStateSet.a = tO;
	//ThreeStateSet.b = tC;
	//dThreeStateSet.a = dtO_dt;
	//dThreeStateSet.b = dtC_dt;

	//RateConstSet.a = koi;
	//RateConstSet.b = kio;
	//RateConstSet.c = kic;
	//RateConstSet.d = kci;
	//RateConstSet.e = koc;
	//RateConstSet.f = kco;

	//pThreeStateSet = ThreeState(ThreeStateSet, dThreeStateSet, RateConstSet, m_dDt);
	//
	//tO = pThreeStateSet.a;
	//tC = pThreeStateSet.b;
	//tI = 1.0 - tO - tC;		
	//PIP3R = 217.2;
	//IIP3R = PIP3R*(CaiSR-CaPMU)*tO;		

}
// initialize stochastic process (IP3 receptor channel)
sixelement CCellPart::initializeSTOCHASTICIP3R(sixelement dPoIP3R, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dIP3, double dCai)
{
	//------ IP3 Receptor (stochastic model) --------------//				
	// [headless] srand((unsigned)time(NULL)) removed: per-object RNG (Reseed)
	int NN;

	// calculate initial stay time of each state
	int unitA = 0;
	int unitB = 0;
	int unitC = 0;
	int unitD = 0;

	// initialize channels be sampled
	// IP3  rates
	kucIP3R = 400.0*dIP3;
	kcuIP3R = 400.0*0.9434e-3;			// d1: 0.13, d3: 0.9434 uM	CIP3R
	kcoIP3R = 20.0*dCai;
	kocIP3R = 20.0*0.08234e-3;			// d5: 0.08234 uM
	kobIP3R = 0.2*dCai;
	kboIP3R = 0.2*1.049e-3;				// d2: 1.049, d4: 0.1445 uM		//*/

	for (NN = 0; NN<iChannelNumber; NN++)
	{
		if (NN<int(double(iChannelNumber)*dPoIP3R.a))
		{
			iState[NN] = 0;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / (kcuIP3R + kcoIP3R);
			unitA++;
		}
		else if (NN<int(double(iChannelNumber)*(dPoIP3R.a + dPoIP3R.d)))
		{
			iState[NN] = 1;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / (kocIP3R + kobIP3R);
			unitB++;
		}
		else if (NN<int(double(iChannelNumber)*(dPoIP3R.a + dPoIP3R.d + dPoIP3R.b)))
		{
			iState[NN] = 3;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / kboIP3R;
			unitD++;
		}
		else
		{
			iState[NN] = 2;
			dTimePredict[NN] = 0;
			dStayingTimePredict[NN] = fabs(log(double(Rand() % RAND_MAX) / double(RAND_MAX))) / kucIP3R;
			unitC++;
		}
	}
	dPoIP3R.a = double(unitA) / double(iChannelNumber);
	dPoIP3R.d = double(unitB) / double(iChannelNumber);
	dPoIP3R.c = double(unitC) / double(iChannelNumber);
	dPoIP3R.b = 1.0 - dPoIP3R.a - dPoIP3R.d - dPoIP3R.c;
		
	return dPoIP3R;
}
void CCellPart::initializeCurrentTimeSTOCHASTICIP3R(int iChannelNumber, double* dTimePredict)
{
	int NN;
	for (NN = 0; NN<iChannelNumber; NN++)
	{
		dTimePredict[NN] = 0;
	}
}
void CCellPart::SERCA()
{		
	// SERCA 1
	double KmCao_SERCA, KmCai_SERCA;
	double alpha_YSERCA, beta_YSERCA;
	double k1SERCA, k2SERCA, k3SERCA, k4SERCA;
	double pE1CaSERCA, pE2CaSERCA, pE1SERCA, pE2SERCA;
	double dYSERCA_dt;


	k1SERCA = 0.01;
	k2SERCA = 1.0;
	k3SERCA = 1.0;
	k4SERCA = 0.01;

	KmCao_SERCA = 0.86*1.5*0.002; //1.5*0.002;				// mM	
	KmCai_SERCA = 0.08; //0.08;					// mM

	//KmCao_SERCA = 0.004;				// mM			remremrem
	//KmCai_SERCA = 0.1;					// mM			remremrem

	pE1CaSERCA = 1.0 / (1.0 + KmCai_SERCA / CaUP);
	pE2CaSERCA = 1.0 / (1.0 + KmCao_SERCA / Cai);
	pE1SERCA = 1.0 / (1.0 + CaUP / KmCai_SERCA);
	pE2SERCA = 1.0 / (1.0 + Cai / KmCao_SERCA);

	alpha_YSERCA = k2SERCA * pE2CaSERCA + k4SERCA * pE2SERCA;
	beta_YSERCA = k1SERCA * pE1CaSERCA + k3SERCA * pE1SERCA;

	dYSERCA_dt = alpha_YSERCA*(1 - YSERCA) - beta_YSERCA*YSERCA;

	//----------- Integrate ------------------------------------------------//
	YSERCA = TwoState(YSERCA, dYSERCA_dt, alpha_YSERCA, beta_YSERCA, m_dDt);

	double PSERCA = 0.1*7728.0;			// remremrem
	ISERCA = SF_SERCA*1.0*PSERCA * (k2SERCA*pE2CaSERCA*(1 - YSERCA) - k1SERCA*pE1CaSERCA*YSERCA);
}

void CCellPart::InitialCellVariables()
{ 
	// from file
	SF_Na = 0.000000e+000;
	INa = 1.000000e-012;
	INa_Na = 1.000000e-012;
	INa_K = 1.000000e-012;
	NaAct = 0.000000e+000;
	NaInact = 9.900000e-001;
	NaInact2 = 1.000000e-002;
	SF_NaLCN = 2.000000e+000;
	PNaLCN = 1.350000e-001;
	INaLCN = -5.888044e+001;
	INaLCN_Na = -6.026994e+001;
	INaLCN_K = 1.389503e+000;
	SF_CaL = 0*1.000000e+000;
	CaAct = 8.980201e-008;
	CaInact = 6.804725e-002;
	PoCaL = 6.110780e-009;
	ICaL = -3.459157e-006;
	ICaL_Na = -7.142717e-009;
	ICaL_K = 3.659398e-009;
	ICaL_Ca = -3.455674e-006;
	SF_VDDR = 0*1.000000e+000;
	YVDDR1 = 4.428786e-004;
	YVDDR2 = 7.261682e-001;
	IVDDR = -5.001390e-002;
	SF_Kto = 0*0.000000e+000;
	KtoAct = 2.191519e-003;
	KtoInact = 9.842542e-001;
	KtoInact2 = 6.421196e-001;
	IKto = 0.000000e+000;
	IKto_Na = 0.000000e+000;
	IKto_K = 0.000000e+000;
	SF_Kr = 0*0.000000e+000;
	YKr1 = 1.928000e-003;
	YKr2 = 2.086000e-001;
	YKr3 = 9.677000e-001;
	IKr = 0.000000e+000;
	SF_Ks = 0*0.000000e+000;
	SF_K1 = 0*1.000000e+000;
	nK1 = 8.632682e-001;
	IK1 = 6.981454e+001;
	SF_SK = 1.000000e+00;
	kc2c1SK = 0.000000e+00;
	kc1c2SK = 0.000000e+00;
	kc1oSK = 2.033294e+00;
	koc1SK = 5.472752e+05;
	ISK = 0.000000e+00;
	SF_NoSK = 20000;
	PoSK = { 0.000000e+00,9.480000e-01,5.200000e-02,0.000000e+00, 0, 0};
	SF_ClCa = 0.000000e+000;
	CaClCa = 1.505251e-004;
	OClCa = 3.249719e-003;
	IClCa = -3.883260e+000;
	SF_SOC = 1.000000e+000;
	OSOC = 7.908945e-001;
	ISOC = -4.525868e+000;
	JSOC = 6.555787e-005;
	SF_NaCa = 1.000000e+000;
	Xtotal = 1.000000e+000;
	X2Na3 = 3.566759e-002;
	X2Na2 = 2.229218e-002;
	X2Na1 = 1.393255e-002;
	X2 = 8.707782e-003;
	X2Ca = 1.135154e-002;
	X1Ca = 1.046436e-002;
	X1 = 3.873480e-002;
	X1Na1 = 9.327766e-002;
	X1Na2 = 2.246276e-001;
	X1Na3 = 5.409440e-001;
	YNaCa = 9.798973e-001;
	INaCa = 7.518478e-003;
	SF_NaK = 1.000000e+000;
	YNaK = 5.983000e-001;
	INaK = 2.022200e+001;
	SF_PMCA = 1.000000e+000;
	IPMCA = 1.151047e-001;
	SF_NKCC = 1.000000e+000;
	JNKCC_Na = 4.055627e-003;
	JNKCC_K = 4.055627e-003;
	JNKCC_Cl = 8.111254e-003;
	Etotal = 1.000000e+000;
	E1 = 1.579359e-003;
	E1Na = 3.489729e-003;
	E1NaCl = 1.059168e-001;
	E1NaClK = 5.271158e-003;
	E1NaClKCl = 1.600116e-001;
	E2NaClKCl = 6.039509e-002;
	E2ClKCl = 6.046385e-001;
	E2KCl = 3.186887e-002;
	E2Cl = 2.548814e-002;
	E2 = 1.340750e-003;
	SF_Cab = 0.000000e+000;
	ICab = 0.000000e+000;
	SF_Clb = 0*1.000000e+000;
	IClb = -2.268677e+001;
	SF_AutonomousInward = 0.000000e+000;
	IAutonomousInward = 0.000000e+000;
	IAutonomousInward_K = 0.000000e+000;
	IAutonomousInward_Na = -0.000000e+000;
	IAutonomousInward_Ca = -0.000000e+000;
	SF_IP3R = 1.000000e+000;
	PoIP3R = {0, 1.16e-1, 8.84e-1, 0, 0, 0};
	PoIP3pmuR = {0, 1.16e-1, 8.84e-1, 0, 0, 0};
	powerOIP3R = 5.717870e-004;
	powerOIP3R0 = 7.290000e-004;
	kucIP3R = 1.358491e-001;
	kcuIP3R = 3.773600e-001;
	kcoIP3R = 1.782624e-003;
	kocIP3R = 1.646800e-003;
	kobIP3R = 1.782624e-005;
	kboIP3R = 2.098000e-004;
	IIP3R = 3.268562e+000;
	aUIP3R = 4.609766e-001;
	aCIP3R = 1.659768e-001;
	aOIP3R = 2.148620e-001;
	aBIP3R = 1.581846e-001;
	S000 = 2.484670e-001;
	S100 = 5.210876e-001;
	S101 = 1.901195e-002;
	S001 = 6.579107e-002;
	S010 = 4.235619e-002;
	S110 = 8.882987e-002;
	S111 = 3.240970e-003;
	S011 = 1.121541e-002;
	poweraOIP3R = 7.009339e-004;
	poweraOIP3R0 = 8.448703e-004;
	SF_NoIP3R = 750;
	tO = 1.223854e-004;
	tC = 5.057212e-002;
	tI = 9.493055e-001;
	SF_SRT = 1.000000e+000;
	Itr = 5.456317e+000;
	SF_SRL = 1.000000e+000;
	Ileak = 2.006184e+000;
	SF_SERCA = 1.000000e+000;
	YSERCA = 4.226115e-001;
	YSERCA2 = 4.563000e-001;
	ISERCA = 9.748728e+000;
	Nai = 6.320874e+000;
	Nae = 1.400000e+002;
	NaPMU = 6.321518e+000;
	Ki = 1.356016e+002;
	Ke = 5.400000e+000;
	KPMU = 1.356023e+002;
	Cai = 8.913122e-005;
	Cae = 1.800000e+000;
	CaPMUt = 1.717251e-004;
	CaPMU = 1.433250e-004;
	Cat = 1.071804e-004;
	CaSR = 8.064770e+000;
	CaiSR = 1.517202e+000;
	CaUP = 1.783366e+000;
	Cli = 8.729039e+001;
	Cle = 1.400000e+002;
	ClPMU = 8.746454e+001;
	CMDNf = 4.819509e-004;
	CMDNPMU = 4.716000e-004;
	CSQNf = 3.452440e+000;
	TRPNf = 0.000000e+000;
	IP3 = 3.396232e-004;
	PIP2 = 1.364663e-003;
	IP3total = 3.300000e-003;
	CF_Na = -4.058582e+002;
	CF_K = 1.039658e+001;
	CF_Ca = -9.817787e+000;
	CF_Cl = -2.268677e+002;
	CF_NaPMU = -4.058581e+002;
	CF_KPMU = 1.039672e+001;
	CF_CaPMU = -9.817786e+000;
	CF_ClPMU = -2.273740e+002;
	testVar1 = 1.000000e-001;
	testVar2 = 1.000000e-001;	
	RatioCaLtoSR = 0.000000e+000;
	RatioNaCatoSR = 0.000000e+000;
	monitorVar1 = 0;
	monitorVar2 = 0;
	Ilin = 0.000000e+000;
	Iext = 0.000000e+000;
	Itot = 1.769826e-002;
	Vtest = 0.000000e+000;
	SF_AccessResistance = 0.000000e+000;
	SF_PipetteVolumeRatio = 1.000000e+000;
	Vm = -7.233222e+001;
	dVm_dt = -7.079303e-004;
	dtCai = 1.000000e-003;
	dtCai0 = 0.000000e+000;
	dtIP3 = 1.000000e-002;
	dtIP30 = 0.000000e+000;
	NaFromMito = 0.000000e+000;
	CaFromMito = 0.000000e+000;
	KFromMito = 0.000000e+000;
	ATPFromMito = 0.000000e+000;
	ADPFromMito = 0.000000e+000;

	InitializeStochasticChannels();	
}

void CCellPart::InitializeStochasticChannels()
{
	//------ initialize stochastic array 
	int SF_NoIP3nonpmuR = int(fractionPMU_IIP3R*double(SF_NoIP3R));
	int SF_NoIP3pmuR = SF_NoIP3R - SF_NoIP3nonpmuR;	

	// [headless] the original set all pointers to NULL *before* the delete[]
	// block, so every call leaked the previous arrays. DeepArray::allocate
	// frees the old buffer before allocating the new one.
	iStateIP3R.allocate(SF_NoIP3nonpmuR);
	dTimePredictIP3R.allocate(SF_NoIP3nonpmuR);
	dStayingTimePredictIP3R.allocate(SF_NoIP3nonpmuR);
	iStateIP3pmuR.allocate(SF_NoIP3pmuR);
	dTimePredictIP3pmuR.allocate(SF_NoIP3pmuR);
	dStayingTimePredictIP3pmuR.allocate(SF_NoIP3pmuR);

	iStateSK.allocate(SF_NoSK);
	dTimePredictSK.allocate(SF_NoSK);
	dStayingTimePredictSK.allocate(SF_NoSK);
	
	// initialize all dynamic arrays
	PoIP3R = initializeSTOCHASTICIP3R(PoIP3R, SF_NoIP3nonpmuR, iStateIP3R, dStayingTimePredictIP3R, dTimePredictIP3R, IP3, Cai);
	PoIP3pmuR = initializeSTOCHASTICIP3R(PoIP3pmuR, SF_NoIP3pmuR, iStateIP3pmuR, dStayingTimePredictIP3pmuR, dTimePredictIP3pmuR, IP3, CaPMU);

	PoSK = initializeSTOCHASTICSK(PoSK, SF_NoSK, iStateSK, dStayingTimePredictSK, dTimePredictSK, Cai);	
}
void CCellPart::StoreVariables()
{
	dAllParameters[0] = NaAct;
	dAllParameters[1] = NaInact;
	dAllParameters[2] = NaInact2;
	dAllParameters[3] = CaAct;
	dAllParameters[4] = CaInact;
	dAllParameters[5] = YVDDR1;
	dAllParameters[6] = YVDDR2;
	dAllParameters[7] = KtoAct;
	dAllParameters[8] = KtoInact;
	dAllParameters[9] = KtoInact2;
	dAllParameters[10] = YKr1;
	dAllParameters[11] = YKr2;
	dAllParameters[12] = YKr3;
	dAllParameters[13] = PoSK.c;
	dAllParameters[14] = OClCa;
	dAllParameters[15] = YNaCa;
	dAllParameters[16] = YNaK;
	dAllParameters[17] =  Etotal;
	dAllParameters[18] =  E1;
	dAllParameters[19] =  E1Na;
	dAllParameters[20] =  E1NaCl;
	dAllParameters[21] =  E1NaClK;
	dAllParameters[22] =  E1NaClKCl;
	dAllParameters[23] =  E2NaClKCl;
	dAllParameters[24] =  E2ClKCl;
	dAllParameters[25] =  E2KCl;
	dAllParameters[26] =  E2Cl;
	dAllParameters[27] =  E2;
	dAllParameters[28] = aUIP3R;
	dAllParameters[29] = aCIP3R;
	dAllParameters[30] = aOIP3R;
	dAllParameters[31] = aBIP3R;
	dAllParameters[32] = S000;
	dAllParameters[33] = S100;
	dAllParameters[34] = S101;
	dAllParameters[35] = S001;
	dAllParameters[36] = S010;
	dAllParameters[37] = S110;
	dAllParameters[38] = S111;
	dAllParameters[39] = S011;
	dAllParameters[40] = tO;
	dAllParameters[41] = tC;
	dAllParameters[42] = tI;
	dAllParameters[43] = YSERCA;
	dAllParameters[44] = YSERCA2;
	dAllParameters[45] = Nai;
	dAllParameters[46] = NaPMU;
	dAllParameters[47] = Ki;	
	dAllParameters[48] = KPMU;
	dAllParameters[49] = Cai;
	dAllParameters[50] = CaPMU;
	dAllParameters[51] = Cat;
	dAllParameters[52] = CaSR;
	dAllParameters[53] = CaiSR;
	dAllParameters[54] = CaUP;
	dAllParameters[55] = Cli;
	dAllParameters[56] = ClPMU;
	dAllParameters[57] = CMDNf;
	dAllParameters[58] = CSQNf;
	dAllParameters[59] = TRPNf;
	dAllParameters[60] = IP3;
	dAllParameters[61] = PIP2;
	dAllParameters[62] = Vm;	
	// [headless] CaPMUt and CMDNPMU were not saved/restored in the original, so a
	// failed step (negative concentration) left CaPMUt at its failed value and
	// the retry could loop forever. Added to the store/retrieve set.
	dAllParameters[63] = CaPMUt;
	dAllParameters[64] = CMDNPMU;
}
void CCellPart::RetrieveVariables()
{
NaAct = dAllParameters[0];
NaInact = dAllParameters[1];
NaInact2 = dAllParameters[2];
CaAct = dAllParameters[3];
CaInact = dAllParameters[4];
YVDDR1 = dAllParameters[5];
YVDDR2 = dAllParameters[6];
KtoAct = dAllParameters[7];
KtoInact = dAllParameters[8];
KtoInact2 = dAllParameters[9];
YKr1 = dAllParameters[10];
YKr2 = dAllParameters[11];
YKr3 = dAllParameters[12];
PoSK.c = dAllParameters[13];
OClCa = dAllParameters[14];
YNaCa = dAllParameters[15];
YNaK = dAllParameters[16];
Etotal = dAllParameters[17];
E1 = dAllParameters[18];
E1Na = dAllParameters[19];
E1NaCl = dAllParameters[20];
E1NaClK = dAllParameters[21];
E1NaClKCl = dAllParameters[22];
E2NaClKCl = dAllParameters[23];
E2ClKCl = dAllParameters[24];
E2KCl = dAllParameters[25];
E2Cl = dAllParameters[26];
E2 = dAllParameters[27];
aUIP3R = dAllParameters[28];
aCIP3R = dAllParameters[29];
aOIP3R = dAllParameters[30];
aBIP3R = dAllParameters[31];
S000 = dAllParameters[32];
S100 = dAllParameters[33];
S101 = dAllParameters[34];
S001 = dAllParameters[35];
S010 = dAllParameters[36];
S110 = dAllParameters[37];
S111 = dAllParameters[38];
S011 = dAllParameters[39];
tO = dAllParameters[40];
tC = dAllParameters[41];
tI = dAllParameters[42];
YSERCA = dAllParameters[43];
YSERCA2 = dAllParameters[44];
Nai = dAllParameters[45];
NaPMU = dAllParameters[46];
Ki = dAllParameters[47];
KPMU = dAllParameters[48];
Cai = dAllParameters[49];
CaPMU = dAllParameters[50];
Cat = dAllParameters[51];
CaSR = dAllParameters[52];
CaiSR = dAllParameters[53];
CaUP = dAllParameters[54];
Cli = dAllParameters[55];
ClPMU = dAllParameters[56];
CMDNf = dAllParameters[57];
CSQNf = dAllParameters[58];
TRPNf = dAllParameters[59];
IP3 = dAllParameters[60];
PIP2 = dAllParameters[61];
Vm = dAllParameters[62];
CaPMUt = dAllParameters[63];
CMDNPMU = dAllParameters[64];
}
bool CCellPart::RunCellModel(bool AP, double InjectAmplitude)
{
	StoreVariables();	

	if (Ca_Buffer() == false)			// 농도가 마이너스가 되면 time step 줄여서 다시 시작
	{
		RetrieveVariables();
		return false;
	}		

	Constant_Field();	

	int SF_NoIP3nonpmuR = int(fractionPMU_IIP3R*double(SF_NoIP3R));
	ip3Ag += (ip3AgTarget - ip3Ag) * m_dDt / tauIP3Ag;		// [modified PDGFa] agonist-evoked IP3
	double IP3eff = IP3 + ip3Ag;
	PoIP3R = IP3R_Current(PoIP3R, SF_NoIP3nonpmuR, iStateIP3R, dStayingTimePredictIP3R, dTimePredictIP3R, IP3eff, Cai);
	nonPMU_IP3R = PoIP3R.e;
	PoIP3pmuR = IP3R_Current(PoIP3pmuR, SF_NoIP3R - SF_NoIP3nonpmuR, iStateIP3pmuR, dStayingTimePredictIP3pmuR, dTimePredictIP3pmuR, IP3eff, CaPMU);
	PMU_IP3R = PoIP3pmuR.e;
	IIP3R = nonPMU_IP3R + PMU_IP3R;
		
	PoSK = SK_Current(PoSK, SF_NoSK, iStateSK, dStayingTimePredictSK, dTimePredictSK, Cai);	
	ISK = PoSK.e;
	
	CaL_Current();	
	VDDR_Current();
	
	PMU_IK1 = fractionPMU_IK1*K1_Current(KPMU);		
	nonPMU_IK1 = (1.0-fractionPMU_IK1)*K1_Current(Ki);
	IK1 = PMU_IK1 + nonPMU_IK1;

	NaCa_Current();	

	PMU_INaK = fractionPMU_INaK*NaK_Current(NaPMU, KPMU);
	nonPMU_INaK = (1.0-fractionPMU_INaK)*NaK_Current(Nai, Ki);
	INaK = PMU_INaK + nonPMU_INaK;

	NaLCN_Current();
	//Na_Current();
	AutonomousInward_Current();
	ClCa_Current();
	SOC_Current();
	//Kr_Current();			
	
	// does not consider
	//Cab_Current();
	Clb_Current();
	//Kto_Current();		
	
	SERCA();		
	PMCA();	
	NKCC();
	if (SR_Concentration() == false)		// 농도가 마이너스가 되면 time step 줄여서 다시 시작
	{
		RetrieveVariables();
		return false;
	}	
	
	if (Ion_Concentration() == false)		// 농도가 마이너스가 되면 time step 줄여서 다시 시작
	{
		RetrieveVariables();
		return false;
	}	

	if (IP3metabolism() == false)
	{
		RetrieveVariables();
		return false;
	}	

	Membrane_Potential(AP, InjectAmplitude);
	if (AP == false) Vm = InjectAmplitude;	// should be restored remrem

	return true;
}
// Save Cell Parameters
void CCellPart::SaveCellParamters()
{
		// [headless] Win32 CreateFile("temp.cel")/AfxMessageBox removed; plain C stdio only.
		FILE* AsciiFile = fopen("allparameters.dat", "w+t");
		if (AsciiFile == NULL) return;
				
			fprintf(AsciiFile, "SF_Na = %e;\n", SF_Na);
			fprintf(AsciiFile, "INa = %e;\n", INa);
			fprintf(AsciiFile, "INa_Na = %e;\n", INa_Na);
			fprintf(AsciiFile, "INa_K = %e;\n", INa_K);
			fprintf(AsciiFile, "NaAct = %e;\n", NaAct);
			fprintf(AsciiFile, "NaInact = %e;\n", NaInact);
			fprintf(AsciiFile, "NaInact2 = %e;\n", NaInact2);

			fprintf(AsciiFile, "SF_NaLCN = %e;\n", SF_NaLCN);
			fprintf(AsciiFile, "PNaLCN = %e;\n", PNaLCN);
			fprintf(AsciiFile, "INaLCN = %e;\n", INaLCN);
			fprintf(AsciiFile, "INaLCN_Na = %e;\n", INaLCN_Na);
			fprintf(AsciiFile, "INaLCN_K = %e;\n", INaLCN_K);
			
			fprintf(AsciiFile, "SF_CaL = %e;\n", SF_CaL);
			fprintf(AsciiFile, "CaAct = %e;\n", CaAct);
			fprintf(AsciiFile, "CaInact = %e;\n", CaInact);
			fprintf(AsciiFile, "PoCaL = %e;\n", PoCaL);
			fprintf(AsciiFile, "ICaL = %e;\n", ICaL);
			fprintf(AsciiFile, "ICaL_Na = %e;\n", ICaL_Na);
			fprintf(AsciiFile, "ICaL_K = %e;\n", ICaL_K);
			fprintf(AsciiFile, "ICaL_Ca = %e;\n", ICaL_Ca);
			
			fprintf(AsciiFile, "SF_VDDR = %e;\n", SF_VDDR);
			fprintf(AsciiFile, "YVDDR1 = %e;\n", YVDDR1);
			fprintf(AsciiFile, "YVDDR2 = %e;\n", YVDDR2);
			fprintf(AsciiFile, "IVDDR = %e;\n", IVDDR);
			
			fprintf(AsciiFile, "SF_Kto = %e;\n", SF_Kto);
			fprintf(AsciiFile, "KtoAct = %e;\n", KtoAct);
			fprintf(AsciiFile, "KtoInact = %e;\n", KtoInact);
			fprintf(AsciiFile, "KtoInact2 = %e;\n", KtoInact2);
			fprintf(AsciiFile, "IKto = %e;\n", IKto);
			fprintf(AsciiFile, "IKto_Na = %e;\n", IKto_Na);
			fprintf(AsciiFile, "IKto_K = %e;\n", IKto_K);
			
			fprintf(AsciiFile, "SF_Kr = %e;\n", SF_Kr);
			fprintf(AsciiFile, "YKr1 = %e;\n", YKr1);
			fprintf(AsciiFile, "YKr2 = %e;\n", YKr2);
			fprintf(AsciiFile, "YKr3 = %e;\n", YKr3);
			fprintf(AsciiFile, "IKr = %e;\n", IKr);
			
			fprintf(AsciiFile, "SF_Ks = %e;\n", SF_Ks);
		
			fprintf(AsciiFile, "SF_K1 = %e;\n", SF_K1);
			fprintf(AsciiFile, "nK1 = %e;\n", nK1);
			fprintf(AsciiFile, "IK1 = %e;\n", IK1);
			
			fprintf(AsciiFile, "SF_SK = %e;\n", SF_SK);
			fprintf(AsciiFile, "kc2c1SK = %e;\n", kc2c1SK);
			fprintf(AsciiFile, "kc1c2SK = %e;\n", kc1c2SK);
			fprintf(AsciiFile, "kc1oSK = %e;\n", kc1oSK);
			fprintf(AsciiFile, "koc1SK = %e;\n", koc1SK);
			fprintf(AsciiFile, "ISK = %e;\n", ISK);
			fprintf(AsciiFile, "SF_NoSK = %d;\n", SF_NoSK);
			fprintf(AsciiFile, "PoSK = {%e,%e,%e,%e,%e,%e};\n", PoSK.a, PoSK.b, PoSK.c, PoSK.d, PoSK.e, PoSK.f);
					
			fprintf(AsciiFile, "SF_ClCa = %e;\n", SF_ClCa);
			fprintf(AsciiFile, "CaClCa = %e;\n", CaClCa);
			fprintf(AsciiFile, "OClCa = %e;\n", OClCa);
			fprintf(AsciiFile, "IClCa = %e;\n", IClCa);

			fprintf(AsciiFile, "SF_SOC = %e;\n", SF_SOC);
			fprintf(AsciiFile, "OSOC = %e;\n", OSOC);
			fprintf(AsciiFile, "ISOC = %e;\n", ISOC);
			fprintf(AsciiFile, "JSOC = %e;\n", JSOC);
			
			fprintf(AsciiFile, "SF_NaCa = %e;\n", SF_NaCa);
			fprintf(AsciiFile, "Xtotal = %e;\n", Xtotal);
			fprintf(AsciiFile, "X2Na3 = %e;\n", X2Na3);
			fprintf(AsciiFile, "X2Na2 = %e;\n", X2Na2);
			fprintf(AsciiFile, "X2Na1 = %e;\n", X2Na1);
			fprintf(AsciiFile, "X2 = %e;\n", X2);
			fprintf(AsciiFile, "X2Ca = %e;\n", X2Ca);
			fprintf(AsciiFile, "X1Ca = %e;\n", X1Ca);
			fprintf(AsciiFile, "X1 = %e;\n", X1);
			fprintf(AsciiFile, "X1Na1 = %e;\n", X1Na1);
			fprintf(AsciiFile, "X1Na2 = %e;\n", X1Na2);
			fprintf(AsciiFile, "X1Na3 = %e;\n", X1Na3);			
			fprintf(AsciiFile, "YNaCa = %e;\n", YNaCa);
			fprintf(AsciiFile, "INaCa = %e;\n", INaCa);
			
			fprintf(AsciiFile, "SF_NaK = %e;\n", SF_NaK);
			fprintf(AsciiFile, "YNaK = %e;\n", YNaK);
			fprintf(AsciiFile, "INaK = %e;\n", INaK);
			
			fprintf(AsciiFile, "SF_PMCA = %e;\n", SF_PMCA);
			fprintf(AsciiFile, "IPMCA = %e;\n", IPMCA);
			
			fprintf(AsciiFile, "SF_NKCC = %e;\n", SF_NKCC);
			fprintf(AsciiFile, "JNKCC_Na = %e;\n", JNKCC_Na);
			fprintf(AsciiFile, "JNKCC_K = %e;\n", JNKCC_K);
			fprintf(AsciiFile, "JNKCC_Cl = %e;\n", JNKCC_Cl);
			fprintf(AsciiFile, "Etotal = %e;\n", Etotal);
			fprintf(AsciiFile, "E1 = %e;\n", E1);
			fprintf(AsciiFile, "E1Na = %e;\n", E1Na);
			fprintf(AsciiFile, "E1NaCl = %e;\n", E1NaCl);
			fprintf(AsciiFile, "E1NaClK = %e;\n", E1NaClK);
			fprintf(AsciiFile, "E1NaClKCl = %e;\n", E1NaClKCl);
			fprintf(AsciiFile, "E2NaClKCl = %e;\n", E2NaClKCl);
			fprintf(AsciiFile, "E2ClKCl = %e;\n", E2ClKCl);
			fprintf(AsciiFile, "E2KCl = %e;\n", E2KCl);
			fprintf(AsciiFile, "E2Cl = %e;\n", E2Cl);
			fprintf(AsciiFile, "E2 = %e;\n", E2);

			fprintf(AsciiFile, "SF_Cab = %e;\n", SF_Cab);
			fprintf(AsciiFile, "ICab = %e;\n", ICab);			
			fprintf(AsciiFile, "SF_Clb = %e;\n", SF_Clb);
			fprintf(AsciiFile, "IClb = %e;\n", IClb);

			fprintf(AsciiFile, "SF_AutonomousInward = %e;\n", SF_AutonomousInward);
			fprintf(AsciiFile, "IAutonomousInward = %e;\n", IAutonomousInward);
			fprintf(AsciiFile, "IAutonomousInward_K = %e;\n", IAutonomousInward_K);
			fprintf(AsciiFile, "IAutonomousInward_Na = %e;\n", IAutonomousInward_Na);
			fprintf(AsciiFile, "IAutonomousInward_Ca = %e;\n", IAutonomousInward_Ca);

			fprintf(AsciiFile, "SF_IP3R = %e;\n", SF_IP3R);
			fprintf(AsciiFile, "PoIP3R = {%e,%e,%e,%e,%e,%e};\n", PoIP3R.a, PoIP3R.b, PoIP3R.c, PoIP3R.d, PoIP3R.e, PoIP3R.f);
			fprintf(AsciiFile, "PoIP3pmuR = {%e,%e,%e,%e,%e,%e};\n", PoIP3pmuR.a, PoIP3pmuR.b, PoIP3pmuR.c, PoIP3pmuR.d, PoIP3pmuR.e, PoIP3pmuR.f);
			fprintf(AsciiFile, "powerOIP3R = %e;\n", powerOIP3R);
			fprintf(AsciiFile, "powerOIP3R0 = %e;\n", powerOIP3R0);
			fprintf(AsciiFile, "kucIP3R = %e;\n", kucIP3R);
			fprintf(AsciiFile, "kcuIP3R = %e;\n", kcuIP3R);
			fprintf(AsciiFile, "kcoIP3R = %e;\n", kcoIP3R);
			fprintf(AsciiFile, "kocIP3R = %e;\n", kocIP3R);
			fprintf(AsciiFile, "kobIP3R = %e;\n", kobIP3R);
			fprintf(AsciiFile, "kboIP3R = %e;\n", kboIP3R);
			fprintf(AsciiFile, "IIP3R = %e;\n", IIP3R);

			fprintf(AsciiFile, "aUIP3R = %e;\n", aUIP3R);
			fprintf(AsciiFile, "aCIP3R = %e;\n", aCIP3R);
			fprintf(AsciiFile, "aOIP3R = %e;\n", aOIP3R);
			fprintf(AsciiFile, "aBIP3R = %e;\n", aBIP3R);

			fprintf(AsciiFile, "S000 = %e;\n", S000);
			fprintf(AsciiFile, "S100 = %e;\n", S100);
			fprintf(AsciiFile, "S101 = %e;\n", S101);
			fprintf(AsciiFile, "S001 = %e;\n", S001);
			fprintf(AsciiFile, "S010 = %e;\n", S010);
			fprintf(AsciiFile, "S110 = %e;\n", S110);
			fprintf(AsciiFile, "S111 = %e;\n", S111);
			fprintf(AsciiFile, "S011 = %e;\n", S011);
			
			fprintf(AsciiFile, "poweraOIP3R = %e;\n", poweraOIP3R);
			fprintf(AsciiFile, "poweraOIP3R0 = %e;\n", poweraOIP3R0);
			fprintf(AsciiFile, "SF_NoIP3R = %d;\n", SF_NoIP3R);
			fprintf(AsciiFile, "tO = %e;\n", tO);
			fprintf(AsciiFile, "tC = %e;\n", tC);
			fprintf(AsciiFile, "tI = %e;\n", tI);

			fprintf(AsciiFile, "SF_SRT = %e;\n", SF_SRT);
			fprintf(AsciiFile, "Itr = %e;\n", Itr);
			fprintf(AsciiFile, "SF_SRL = %e;\n", SF_SRL);
			fprintf(AsciiFile, "Ileak = %e;\n", Ileak);
			
			fprintf(AsciiFile, "SF_SERCA = %e;\n", SF_SERCA);
			fprintf(AsciiFile, "YSERCA = %e;\n", YSERCA);
			fprintf(AsciiFile, "YSERCA2 = %e;\n", YSERCA2);
			fprintf(AsciiFile, "ISERCA = %e;\n", ISERCA);
			
			fprintf(AsciiFile, "Nai = %e;\n", Nai);
			fprintf(AsciiFile, "Nae = %e;\n", Nae);
			fprintf(AsciiFile, "NaPMU = %e;\n", NaPMU);
			fprintf(AsciiFile, "Ki = %e;\n", Ki);
			fprintf(AsciiFile, "Ke = %e;\n", Ke);
			fprintf(AsciiFile, "KPMU = %e;\n", KPMU);
			fprintf(AsciiFile, "Cai = %e;\n", Cai);
			fprintf(AsciiFile, "Cae = %e;\n", Cae);
			fprintf(AsciiFile, "CaPMUt = %e;\n", CaPMUt);
			fprintf(AsciiFile, "CaPMU = %e;\n", CaPMU);
			fprintf(AsciiFile, "Cat = %e;\n", Cat);
			fprintf(AsciiFile, "CaSR = %e;\n", CaSR);
			fprintf(AsciiFile, "CaiSR = %e;\n", CaiSR);
			fprintf(AsciiFile, "CaUP = %e;\n", CaUP);
			fprintf(AsciiFile, "Cli = %e;\n", Cli);
			fprintf(AsciiFile, "Cle = %e;\n", Cle);
			fprintf(AsciiFile, "ClPMU = %e;\n", ClPMU);
			
			fprintf(AsciiFile, "CMDNf = %e;\n", CMDNf);
			fprintf(AsciiFile, "CMDNPMU = %e;\n", CMDNPMU);
			fprintf(AsciiFile, "CSQNf = %e;\n", CSQNf);
			fprintf(AsciiFile, "TRPNf = %e;\n", TRPNf);

			fprintf(AsciiFile, "IP3 = %e;\n", IP3);
			fprintf(AsciiFile, "PIP2 = %e;\n", PIP2);
			fprintf(AsciiFile, "IP3total = %e;\n", IP3total);

			fprintf(AsciiFile, "CF_Na = %e;\n", CF_Na);
			fprintf(AsciiFile, "CF_K = %e;\n", CF_K);
			fprintf(AsciiFile, "CF_Ca = %e;\n", CF_Ca);
			fprintf(AsciiFile, "CF_Cl = %e;\n", CF_Cl);
			fprintf(AsciiFile, "CF_NaPMU = %e;\n", CF_NaPMU);
			fprintf(AsciiFile, "CF_KPMU = %e;\n", CF_KPMU);			
			fprintf(AsciiFile, "CF_CaPMU = %e;\n", CF_CaPMU);
			fprintf(AsciiFile, "CF_ClPMU = %e;\n", CF_ClPMU);

			fprintf(AsciiFile, "testVar1 = %e;\n", testVar1);
			fprintf(AsciiFile, "testVar2 = %e;\n", testVar2);
			fprintf(AsciiFile, "RatioCaLtoSR = %e;\n", RatioCaLtoSR);
			fprintf(AsciiFile, "RatioNaCatoSR = %e;\n", RatioNaCatoSR);
			
			fprintf(AsciiFile, "Ilin = %e;\n", Ilin);
			fprintf(AsciiFile, "Iext = %e;\n", Iext);
			fprintf(AsciiFile, "Itot = %e;\n", Itot);
			fprintf(AsciiFile, "Vtest = %e;\n", Vtest);
			fprintf(AsciiFile, "SF_AccessResistance = %e;\n", SF_AccessResistance);
			fprintf(AsciiFile, "SF_PipetteVolumeRatio = %e;\n", SF_PipetteVolumeRatio);
			fprintf(AsciiFile, "Vm = %e;\n", Vm);
			fprintf(AsciiFile, "dVm_dt = %e;\n", dVm_dt);
			
			fprintf(AsciiFile, "dtCai = %e;\n", dtCai);
			fprintf(AsciiFile, "dtCai0 = %e;\n", dtCai0);
			fprintf(AsciiFile, "dtIP3 = %e;\n", dtIP3);
			fprintf(AsciiFile, "dtIP30 = %e;\n", dtIP30);

			fprintf(AsciiFile, "NaFromMito = %e;\n", NaFromMito);
			fprintf(AsciiFile, "CaFromMito = %e;\n", CaFromMito);
			fprintf(AsciiFile, "KFromMito = %e;\n", KFromMito);
			fprintf(AsciiFile, "ATPFromMito = %e;\n", ATPFromMito);
			fprintf(AsciiFile, "ADPFromMito = %e;\n", ADPFromMito);

		fclose(AsciiFile);		
}
void CCellPart::Membrane_Potential(bool AP, double Iinj)
{	
	Itot = ISOC + ICaL + IVDDR + IK1 + ISK + IAutonomousInward + IClCa + IClb + IKr + INaK + INaCa + INaLCN + INa + 0*IPMCA;		// total current	
	
	Iext = Iinj; 	
	
	if (AP) dVm_dt = -(Itot + Iext)/Cap;	// insert temporary code, injection current, 2007.09.18
	else dVm_dt = 1.0 / Cap * ((Iinj - Vm) / (0.001*SF_AccessResistance) - Itot); 

	Vm = Euler(Vm, dVm_dt, m_dDt);		
}
void CCellPart::Constant_Field()
{	
	// experimental cell condition for patch clamp
	//double Cli = 135.0;			// 26.0 mM, http://www.ncbi.nlm.nih.gov/pubmed/?term=Zhu+Y+Parsons+SP PMID: 20403100, remrem
	//double Cle = 135.0;			// 130.0 mM	remrem

	// real cell condition
	//double Cli = 26.0;			// 26.0 mM, http://www.ncbi.nlm.nih.gov/pubmed/?term=Zhu+Y+Parsons+SP PMID: 20403100, remrem
	//double Cle = 130.0;			// 130.0 mM	remrem

	if (Vm == 0)			// modified code
	{
		CF_Na = Nai - Nae;
		CF_K = Ki - Ke;
		CF_Ca = Cai - Cae;
		CF_Cl = Cle - Cli;
		CF_NaPMU = NaPMU - Nae;
		CF_KPMU = KPMU - Ke;
		CF_CaPMU = CaPMU - Cae;
		CF_ClPMU = Cle - ClPMU;
	}
	else			
	{		
		CF_Na = mFRT * Vm * (Nai - Nae * exp(-Vm * mFRT)) / (1 - exp(-Vm*mFRT));
		CF_K = mFRT * Vm * (Ki - Ke * exp(-Vm * mFRT)) / (1 - exp(-Vm*mFRT));
		CF_Ca = 2.0 * mFRT * Vm * (Cai - Cae * exp(-2.0 * Vm * mFRT)) / (1 - exp(-2.0 * Vm *mFRT));
		CF_Cl = mFRT * Vm * (Cli - Cle * exp(Vm * mFRT)) / (1 - exp(Vm*mFRT));
		CF_NaPMU = mFRT * Vm * (NaPMU - Nae * exp(-Vm * mFRT)) / (1 - exp(-Vm*mFRT));
		CF_KPMU = mFRT * Vm * (KPMU - Ke * exp(-Vm * mFRT)) / (1 - exp(-Vm*mFRT));
		CF_CaPMU = 2.0 * mFRT * Vm * (CaPMU - Cae * exp(-2.0 * Vm * mFRT)) / (1 - exp(-2.0 * Vm *mFRT));
		CF_ClPMU = mFRT * Vm * (ClPMU - Cle * exp(Vm * mFRT)) / (1 - exp(Vm*mFRT));
	}
}
bool CCellPart::Ca_Buffer()
{	
	int N = 1;
    do
	{
            Cai = Cat / (1.0 + (TRPNf/Km_TRPN + CMDNf/Km_CMDN));            
            CMDNf = CMDNt / (1.0 + (Cai/Km_CMDN));  
			TRPNf = TRPNt / (1.0 + (Cai/Km_TRPN));
            N = N + 1;
	} while(N<=10);
            N = 1;
    do
	{
			CaiSR = CaSR / (1.0 + (CSQNf/Km_CSQN));
            CSQNf = CSQNt / (1.0 + (CaiSR/Km_CSQN));
            N = N + 1;            
    } while(N<10);		

	do
	{
		CaPMU = CaPMUt / (1.0 + (CMDNPMU / Km_CMDN));
		CMDNPMU = CMDNt / (1.0 + (CaPMU / Km_CMDN));
		N = N + 1;
	} while (N<10);
	
	if (Cai >= 0 && CaiSR >= 0 && CaPMU >= 0) return true;
	else return false;

	/*double bmyo = CMDNt - Cat + Km_CMDN;
	double cmyo = Km_CMDN * Cat;

	Cai = (pow(bmyo*bmyo+4*cmyo, 0.5) - bmyo) * 0.5;

	bmyo = CSQNt - CaSR + Km_CSQN;
	cmyo = Km_CSQN * CaSR;
	CaiSR = (pow(bmyo*bmyo+4*cmyo, 0.5) - bmyo) * 0.5;	*/		
}
bool CCellPart::SR_Concentration()
{
	Itr = SF_SRT*10.0*20.5*(CaUP - CaiSR);
	Ileak = SF_SRL*1.7*2.5*0.45*(CaUP - Cai);
	double dCaSR_dt;		
	dCaSR_dt = (Itr-IIP3R)/2.0/mF/Vl_SR;		
	CaSR = Euler(CaSR, dCaSR_dt, m_dDt);		

	double dCaUP_dt;		
	dCaUP_dt = (ISERCA-Itr-Ileak)/2.0/mF/Vl_UP;		
	CaUP = Euler(CaUP, dCaUP_dt, m_dDt);		
	if (CaUP>5.0)		CaUP=5.0;

	if (CaSR >= 0 && CaUP >= 0) return true;
	else return false;
}
bool CCellPart::IP3metabolism()
{
	double PIP2_IP3, IP3_PIP2, IP3_IP4, IP4_PIP2, IP4_IP3, PIP2_IP4;
	
	/* original
	PIP2_IP3 = 0.2*exp((Vm+48.5)/18.1);		// kab	
	IP3_PIP2 = 0.5*exp((Vm+100.0)/(-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035*exp((Vm+100)/(-25.5));					//kca		
	PIP2_IP4 = 0.000;							
	*/

	PIP2_IP3 = 0.2*exp((Vm + 48.5) / 18.1)*(Cai) / (Cai + 0.00001);		// kab	
	PIP2_IP3 *= (1.0 + agonistScale);		// [headless] neural/agonist drive (default 0 -> unchanged)
	IP3_PIP2 = 0.5*exp((Vm + 100.0) / (-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035*exp((Vm + 100) / (-25.5));					//kca		
	PIP2_IP4 = 0.000;

	/*/ remrem for IClCa
	PIP2_IP3 = 0.2*exp((Vm+48.5)/18.1)*(Cai)/(Cai+0.00001);		// kab	
	IP3_PIP2 = 0.5*exp((Vm+100.0)/(-48.5));						// kba	
	IP3_IP4 = 0.004;											// kbc		
	IP4_PIP2 = 0.0005;											// kca	
	PIP2_IP4 = 0.000;							
	*/

	/* original remrem for new adjustment
	PIP2_IP3 = 0.02*exp((Vm+48.5)/18.1)*(Cai)/(Cai+0.00001);		// kab	
	IP3_PIP2 = 0.05*exp((Vm+100.0)/(-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.00035;											//kca		
	PIP2_IP4 = 0.000;							
	*/

	/* remrem for no voltage-dependency, V1, 2003.12.6
	PIP2_IP3 = 0.2*exp((Vm+48.5)/18.1)*(Cai)/(Cai+0.00001);		// kab	
	IP3_PIP2 = 0.5*exp((Vm+100.0)/(-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035*exp((-35.0+100)/(-25.5));					//kca		
	PIP2_IP4 = 0.000;													
	*/

	/* remrem for no voltage-dependency, V2
	PIP2_IP3 = 0.02*(Cai)/(Cai+0.00001); //0.2*exp((Vm+48.5)/18.1)*(Cai)/(Cai+0.00001);		// kab	
	IP3_PIP2 = 0.5*exp((Vm+100.0)/(-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035*exp((-35.0+100)/(-25.5));					//kca		
	PIP2_IP4 = 0.000;													
	*/

	/* remrem for no voltage-dependency, V3
	PIP2_IP3 = 0.01*(Cai)/(Cai+0.0007); //0.2*exp((Vm+48.5)/18.1)*(Cai)/(Cai+0.00001);		// kab	
	IP3_PIP2 = 0.05; //0.5*exp((Vm+100.0)/(-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035*exp((-35.0+100)/(-25.5));					//kca		
	PIP2_IP4 = 0.000;													
	*/

	/* remrem for no voltage-dependency, V4
	PIP2_IP3 = 0.01*(Cai)/(Cai+0.0007);														// kab	
	IP3_PIP2 = 0.05;																		//kba	
	IP3_IP4 = 5.0*0.004;																		//kbc			
	IP4_PIP2 = 0.0035*exp((-35.0+100)/(-25.5));												//kca		
	PIP2_IP4 = 0.000;													
	*/

	/* remrem for ca-dependent IP3 production
	double alphaIP3 = 0;
	double v4IP3 = 0.0002; // 0.2 /uM/s or 0.05 /uM/s
	double k4IP3 = 0.0011;	// 1.1 uM, dissociation constant
	double IrIP3 = 1.0;   // 1 /s, loss of IP3
	double dIP3_dt = (Cai + (1.0-alphaIP3)*k4IP3)/(Cai+k4IP3)*v4IP3 - IrIP3*IP3;
	*/

	double dIP3_dt = PIP2*PIP2_IP3 - IP3_IP4*IP3 - IP3_PIP2*IP3;
	double dPIP2_dt = IP4_PIP2*(IP3total-IP3-PIP2) + IP3_PIP2*IP3- (PIP2_IP4+PIP2_IP3)*PIP2;
		

	//----------- Integrate ------------------------------------------------//
	twoelement ThreeStateSet;
	twoelement dThreeStateSet;
	sixelement RateConstSet;
	twoelement pThreeStateSet;

	double fIP3, fPIP2;		// fraction of IP3 and PIP2 to total IPs(IP4, IP3, PIP2)
	fIP3 = IP3/IP3total;
	fPIP2 = PIP2/IP3total;
	ThreeStateSet.a = fIP3;
	ThreeStateSet.b = fPIP2;
	dThreeStateSet.a = dIP3_dt/IP3total;
	dThreeStateSet.b = dPIP2_dt/IP3total;
	
	RateConstSet.a = IP3_IP4;
	RateConstSet.b = 0;
	RateConstSet.c = IP4_PIP2;
	RateConstSet.d = 0;
	RateConstSet.e = IP3_PIP2;
	RateConstSet.f = PIP2_IP3;

	pThreeStateSet = ThreeState(ThreeStateSet, dThreeStateSet, RateConstSet, m_dDt);
	
	fIP3 = pThreeStateSet.a;
	fPIP2 = pThreeStateSet.b;

	IP3 = fIP3*IP3total;
	PIP2 = fPIP2*IP3total;
		
	//IP3 = Euler(IP3, dIP3_dt, m_dDt); //remrem for ca-dependent IP3 production
	if (fPIP2 >= 0 && fIP3 >= 0 && 1.0-fPIP2-fIP3 >= 0) return true;
	else return false;
}
bool CCellPart::Ion_Concentration()
{	
	double dNa_dt, dK_dt, dCa_dt, dCl_dt;
	double dNaPMU_dt, dKPMU_dt, dCaPMU_dt, dClPMU_dt;
	double diffNa = 3.0e-1;
	double diffK = 3.0e-1;
	double diffCa = 2.0e-1;
	double diffCl = 2.0e-3;
	
	double Vl_PipettelVol = Vl_cyt*SF_PipetteVolumeRatio;
	dNa_dt = -(IAutonomousInward_Na + INaLCN_Na + INa_Na + ICaL_Na + 3.0*nonPMU_INaK) / mF / Vl_PipettelVol + diffNa*(NaPMU - Nai);
	dK_dt = -(nonPMU_IK1 + ICaL_K + IAutonomousInward_K + IKr + ISK + INaLCN_K + INa_K - 2.0*nonPMU_INaK) / mF / Vl_PipettelVol + diffK*(KPMU - Ki);
	if (gapCarriedByK) dK_dt += -Iext / mF / Vl_PipettelVol;	// [modified PDGFa] gap-junction current carried by K+ (as in the mCSMC model)
	dCa_dt = -(ICaL_Ca + IVDDR + IAutonomousInward_Ca + ISERCA + IPMCA - Ileak - PMU_IP3R) / 2.0 / mF / Vl_PipettelVol + diffCa*(CaPMU - Cai);
	/*if (timeX > 3500)
		dCa_dt = -(0*ISOC + ICaL_Ca + IVDDR + IAutonomousInward_Ca + ISERCA + IPMCA - Ileak - IIP3R) / 2.0 / mF / Vl_PipettelVol + diffCa*(CaPMU - Cai);
	else
		dCa_dt = -(ISOC + ICaL_Ca + IVDDR + IAutonomousInward_Ca + ISERCA + IPMCA - Ileak - IIP3R) / 2.0 / mF / Vl_PipettelVol + diffCa*(CaPMU - Cai);*/
	dCl_dt = IClb / mF / Vl_PipettelVol + diffCl*(ClPMU - Cli);

	dNaPMU_dt = -(3.0*INaCa + 3.0*PMU_INaK) / mF / Vl_PMU + JNKCC_Na + diffNa*(Nai - NaPMU)*Vl_PipettelVol/Vl_PMU;
	dKPMU_dt = -(PMU_IK1 - 2.0*PMU_INaK) / mF / Vl_PMU + JNKCC_K + diffK*(Ki - KPMU)*Vl_PipettelVol / Vl_PMU;
	dCaPMU_dt = -(ISOC - 2.0*INaCa - nonPMU_IP3R) / 2.0 / mF / Vl_PMU + diffCa*(Cai - CaPMU)*Vl_PipettelVol / Vl_PMU;
	dClPMU_dt = IClCa / mF / Vl_PMU + JNKCC_Cl + diffCl*(Cli - ClPMU)*Vl_PipettelVol / Vl_PMU;

	Nai = Euler(Nai, dNa_dt, m_dDt);
	Ki = Euler(Ki, dK_dt, m_dDt);
	Cat = Euler(Cat, dCa_dt, m_dDt);
	Cli = Euler(Cli, dCl_dt, m_dDt);
	NaPMU = Euler(NaPMU, dNaPMU_dt, m_dDt);
	KPMU = Euler(KPMU, dKPMU_dt, m_dDt);
	CaPMUt = Euler(CaPMUt, dCaPMU_dt, m_dDt);
	ClPMU = Euler(ClPMU, dClPMU_dt, m_dDt);

	if (Nai >= 0 && Ki >= 0 && Cat >= 0 && Cli >=0 && NaPMU>=0 && KPMU>=0 && CaPMUt>=0 && ClPMU>=0) return true;
	else return false;
}


double CCellPart::Euler(double dInput, double dSlope, double dTime)
{	
	return dInput+dSlope*dTime;	
}

double CCellPart::TwoState(double dInput, double dSlope, double Alpha, double Beta, double dTime)
{
	return dInput+dSlope/(Alpha+Beta)*(1-exp(-dTime*(Alpha+Beta)));
}

fourelement CCellPart::Get_State_Coeff(sixelement kVal)
{
	double State_bf, State_sbf;
	fourelement State_Co;
	double State_c1, State_c2, State_c3, State_d1, State_d2, State_d3;
	State_c1 = kVal.e - kVal.c;
	State_c2 = -(kVal.c + kVal.d + kVal.f);
	State_c3 = kVal.c;
	State_d2 = kVal.f - kVal.b;
	State_d1 = -(kVal.b + kVal.a + kVal.e);
	State_d3 = kVal.b;
	State_Co.c = (State_c2*State_d3 - State_c3*State_d2) / (State_c1*State_d2 - State_c2*State_d1);
	State_Co.d = (State_c1*State_d3 - State_c3*State_d1) / (State_c2*State_d1 - State_c1*State_d2);
	State_bf = kVal.a + kVal.b + kVal.c + kVal.d + kVal.e + kVal.f;
	State_sbf = pow(State_bf, 2.0) - 4 * ((kVal.c + kVal.d + kVal.f)*(kVal.b + kVal.a + kVal.e) - (kVal.e - kVal.c)*(kVal.f - kVal.b));
	State_Co.a = (State_bf + pow(fabs(State_sbf), 0.5)) / 2.0;
	State_Co.b = (State_bf - pow(fabs(State_sbf), 0.5)) / 2.0;
	return State_Co;
}

twoelement CCellPart::ThreeState(twoelement dInput, twoelement dSlope, sixelement RateConst, double dTime)
{
	fourelement	State_Coeff;
	double	stO, ndO, stA, ndA;
	twoelement dOutput;

	State_Coeff = Get_State_Coeff(RateConst);

	stO = (dSlope.a + State_Coeff.b * dInput.a - State_Coeff.c * State_Coeff.b) / (State_Coeff.b - State_Coeff.a);
	ndO = (dSlope.a + State_Coeff.a * dInput.a - State_Coeff.c * State_Coeff.a) / (State_Coeff.a - State_Coeff.b);

	stA = (dSlope.b + State_Coeff.b * dInput.b - State_Coeff.d * State_Coeff.b) / (State_Coeff.b - State_Coeff.a);
	ndA = (dSlope.b + State_Coeff.a * dInput.b - State_Coeff.d * State_Coeff.a) / (State_Coeff.a - State_Coeff.b);

	dOutput.a = stO*exp(-State_Coeff.a * dTime) + ndO*exp(-State_Coeff.b * dTime) + State_Coeff.c;
	dOutput.b = stA*exp(-State_Coeff.a * dTime) + ndA*exp(-State_Coeff.b * dTime) + State_Coeff.d;

	return dOutput;
}