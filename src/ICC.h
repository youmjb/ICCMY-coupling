#pragma once
typedef struct {
	double a, b;
}twoelement;
typedef struct {
	double a, b, c, d;
}fourelement;
typedef struct {
	double a, b, c, d, e, f;
}sixelement;

class ICC
{
public:
	const double mR = 8.3143;
	const double mT = 309.15;
	const double mF = 96.4867;	// Faraday constant

	const double mRTF = mR * mT / mF;        //  RT/F (mV)
	const double mFRT = mF / mR / mT;

	const double Vl_cyt = 715.5;			    // cytosolic volume of ICC
	double Vl_UP = 0.03 * 715.5;		// single SR compartment
	double Vl_SR = 0.01 * 715.5;		// single SR compartment
	const double Vl_PMU = 0.05 * 715.5;		// pacemaker unit compartment

	const double Mgi = 1.0;

	// cell capacitance, Katsube, Y., H. Yokoshiki, L. Nguyen, M. Yamamoto, and N. Sperelakis. 1998. 
	// L-type Ca2+ currents in ventricular myocytes from neonatal and adult rats. Can. J. Physiol. Pharmacol
	const double Cap = 25.0;

	const double Km_TRPN = 0.0005;		// 0.0005 mM
	const double Km_CMDN = 0.00238;		// 2.38e-3;	// Km value of calmodulin (MM), 
	const double Km_CSQN = 8.00e-1;		// Km value of calsequestrin (MM)

	const double Kf_TRPN = 39;		// Kf value of troponin (mM/ms)
	const double Kf_CMDN = 100;		// Kf value of calmodulin (mM/ms)
	const double Kf_CSQN = 0;		// Kf value of calsequestrin (mM/ms)

	const double Kb_TRPN = 0.03;	// Kb value of troponin (ms-1)
	const double Kb_CMDN = 0.238;	// Kb value of calmodulin (ms-1)
	const double Kb_CSQN = 0;		// Kb value of calsequestrin (ms-1)

	const double TRPNt = 0.0;		// 0.049
	const double CMDNt = 0.0005;	// 0.0005, calmodulin concentration (mM)		
	const double CSQNt = 10000e-3;	// calsequestrin concentration (mM)

	// ------------------ Integration tools -------------------------------------------------//
	double Euler(double dInput, double dSlope, double dTime);
	double TwoState(double dInput, double dSlope, double Alpha, double Beta, double dTime);
	twoelement ThreeState(twoelement dInput, twoelement dSlope, sixelement RateConst, double dTime);
	fourelement Get_State_Coeff(sixelement kVal);

	double Igap = 0;

	// [Coupling] multiplies the PIP2->IP3 production rate: <1 removes pacemaking (ICC-IM), >1 mimics cholinergic Gq drive
	double ip3ProdScale = 1.0;
	double vddrH = -66.0;	// [Coupling] VDDR inactivation midpoint (mV); -66 = original
	double ip3VK = 1.0, ip3Vref = -72.0;	// [Coupling] voltage dependence of PIP2->IP3: exponent ((ip3Vref+48.5) + ip3VK*(Vm-ip3Vref))/18.1; 1 = original, 0 = fixed at the rest value
	double lumExp = 0.0, lumRef = 2.0;		// luminal (store) Ca2+ sensitivity of IP3R release: x (CaiSR/lumRef)^lumExp (Coupling copy; 0 = original)

	// Na leak channel
	void NaLCN_Current();
	double SF_NaLCN;
	double PNaLCN;
	double INaLCN;
	double INaLCN_Na;
	double INaLCN_K;

	// CaL-channel
	void CaL_Current();
	double SF_CaL;
	double CaAct;
	double CaInact;
	double PoCaL;
	double ICaL;
	double ICaL_Na;
	double ICaL_K;
	double ICaL_Ca;

	// VDDR channel
	void VDDR_Current();
	double SF_VDDR;
	double YVDDR1;
	double YVDDR2;
	double IVDDR;

	// PMCA	
	void PMCA();
	double SF_PMCA;
	double IPMCA;

	// K1_Current
	void K1_Current();
	double SF_K1;
	double nK1;
	double IK1;

	// Ca-activated Chloride channel
	void ClCa_Current();
	double SF_ClCa;
	double CaClCa;
	double OClCa;
	double IClCa;

	// SOC 
	void SOC_Current();
	double SF_SOC;
	double OSOC;
	double ISOC;
	double JSOC;

	// Na-Ca exchanger	
	void NaCa_Current();
	double SF_NaCa;
	double Xtotal;
	double X2Na3;
	double X2Na2;
	double X2Na1;
	double X2;
	double X2Ca;
	double X1Ca;
	double X1;
	double X1Na1;
	double X1Na2;
	double X1Na3;
	double YNaCa;
	double INaCa;

	// Na-K pump
	void NaK_Current();
	double SF_NaK;
	double YNaK;
	double INaK;

	// NKCC
	void NKCC();
	double SF_NKCC;
	double JNKCC_Na;
	double JNKCC_K;
	double JNKCC_Cl;
	double Etotal;
	double E1;
	double E1Na;
	double E1NaCl;
	double E1NaClK;
	double E1NaClKCl;
	double E2NaClKCl;
	double E2ClKCl;
	double E2KCl;
	double E2Cl;
	double E2;

	// background Cl-channel
	void Clb_Current();
	double SF_Clb;
	double IClb;

	// ----------------- Intracellular ion flux ------------------------//
	// IP3 receptor channel (4-state model)
	double SF_IP3R;
	double UIP3R;
	double CIP3R;
	double OIP3R;
	double BIP3R;
	double UIP3R0;
	double CIP3R0;
	double OIP3R0;
	double BIP3R0;		// store value
	double powerOIP3R;
	double powerOIP3R0;				// all subnits
	double kucIP3R;
	double kcuIP3R;
	double kcoIP3R;
	double kocIP3R;
	double kobIP3R;
	double kboIP3R;	// rate constants
	double IIP3R;

	double aUIP3R;
	double aCIP3R;
	double aOIP3R;
	double aBIP3R;
	double S000;
	double S100;
	double S101;
	double S001;
	double S010;
	double S110;
	double S111;
	double S011;

	double poweraOIP3R;
	double poweraOIP3R0;			// all subnits
	int	   SF_NoIP3R;							// number of IP3receptor channel (channel density is fixed)	
	void IP3R_Current();
	void IP3rates();
	double tO;
	double tC;
	double tI;

	// SR transfer, leak
	double SF_SRT;
	double Itr;
	double SF_SRL;
	double Ileak;

	// SERCA 
	void SERCA();
	double SF_SERCA;
	double YSERCA;
	double YSERCA2;
	double ISERCA;

	// SR
	void SR_Concentration();

	// ion concentration
	void Ion_Concentration();
	double Nai;
	double Nae;
	double NaPMU;
	double Ki;
	double Ke;
	double KPMU;
	double Cai;
	double Cae;
	double CaPMUt;
	double CaPMU;
	double Cat;
	double CaSR;
	double CaiSR;
	double CaUP;
	double Cli;
	double Cle;
	double ClPMU;

	// Ca buffers
	void Ca_Buffer();
	double CMDNf;
	double CMDNPMU;
	double CSQNf;
	double TRPNf;

	// IP3 metabolism
	void IP3metabolism();
	double IP3;
	double PIP2;
	double IP3total;

	// Constant field
	void Constant_Field();
	double CF_Na;
	double CF_K;
	double CF_Ca;
	double CF_Cl;
	double CF_NaPMU;
	double CF_KPMU;
	double CF_CaPMU;
	double CF_ClPMU;

	// Run cell model
	void InitialCellVariables();
	void RunCellModel(BOOL AP, double InjectAmplitude);
	void Membrane_Potential(BOOL AP, double Iinj);
	
	// ---------------- Clamp mode variables ----------------------------//
	double Ilin;
	double Iext;
	double Itot;
	double Vtest;
	double SF_AccessResistance;
	double SF_PipetteVolumeRatio;
	double Vm;
	double dVm_dt;
	double timeX;
	double timeX0;
	double m_dDt;
	double ElapsedTime;
	int	   iRepeat;

	// fixed ion concentration
	double dtCai;
	double dtCai0;
	double dtIP3;
	double dtIP30;
};
