// CellPart.h : headless port of the MFC CCellPart (PDGFa project)
// Numerics kept identical to 07_GiSim/PDGFa/CellPart.cpp except where noted
// in comments tagged [headless].
#pragma once
#include <random>
#include <cstdint>
#include <cstdlib>

// [headless] owning heap array with deep-copy semantics. Converts implicitly
// to T* so the original channel code (iState[NN], passing to T* params) is
// unchanged, while the enclosing class gets correct copy/assign for free.
template <class T>
struct DeepArray
{
	T* p = nullptr;
	int n = 0;
	DeepArray() = default;
	DeepArray(const DeepArray& o) { copyFrom(o); }
	DeepArray& operator=(const DeepArray& o) { if (this != &o) { reset(); copyFrom(o); } return *this; }
	~DeepArray() { reset(); }
	void allocate(int k) { reset(); if (k > 0) { p = new T[k]; n = k; } }
	void reset() { delete[] p; p = nullptr; n = 0; }
	operator T*() const { return p; }
private:
	void copyFrom(const DeepArray& o) { if (o.p && o.n > 0) { p = new T[o.n]; n = o.n; for (int i = 0; i < n; ++i) p[i] = o.p[i]; } }
};

//
// 2015. 05. 20. Mouse intestine ICC
/////////////////////////////////////////////////////////////////////////////
// CCellPart document

const double mR = 8.3143;
const double mT = 309.15;
const double mF	=  96.4867;	// Faraday constant

const double mRTF = mR * mT / mF;        //  RT/F (mV)
const double mFRT = mF / mR / mT; 

const double Vl_cyt = 715.5;			    // cytosolic volume of ICC
const double Vl_UP = 0.03*715.5;		// single SR compartment
const double Vl_SR = 0.01*715.5;		// single SR compartment
const double Vl_PMU = 0.05*715.5;		// pacemaker unit compartment

const double Mgi = 1.0;		

// cell capacitance, Katsube, Y., H. Yokoshiki, L. Nguyen, M. Yamamoto, and N. Sperelakis. 1998. 
// L-type Ca2+ currents in ventricular myocytes from neonatal and adult rats. Can. J. Physiol. Pharmacol
const double Cap	=	25.0;		

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
const double CMDNt	= 0.0005;	// 0.0005, calmodulin concentration (mM)		
const double CSQNt	= 10000e-3;	// calsequestrin concentration (mM)

typedef struct {
	double a, b;
}twoelement;
typedef struct {
	double a, b, c, d;
}fourelement;
typedef struct {
	double a, b, c, d, e, f;
}sixelement;

class CCellPart
{

// Attributes
public:	
	// ------------------ Integration tools -------------------------------------------------//
	double Euler(double dInput, double dSlope, double dTime);
	double TwoState(double dInput, double dSlope, double Alpha, double Beta, double dTime);
	twoelement ThreeState(twoelement dInput, twoelement dSlope, sixelement RateConst, double dTime);
	fourelement Get_State_Coeff(sixelement kVal);

	// Na-channel
	void Na_Current();
	double SF_Na;
	double INa;
	double INa_Na;
	double INa_K;
	double NaAct;
	double NaInact;
	double NaInact2;

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
		
	// to-channel
	void Kto_Current();
	double SF_Kto;
	double KtoAct;
	double KtoInact;
	double KtoInact2;
	double IKto;
	double IKto_Na;
	double IKto_K;
	
	// Kr-channel
	void Kr_Current();
	double SF_Kr;
	double YKr1;
	double YKr2;
	double YKr3;
	double IKr;	

	// Ks-channel
	double SF_Ks;

	// K1_Current
	double K1_Current(double Ioni);		
	double SF_K1;
	double nK1;
	double IK1;
	double fractionPMU_IK1 = 0.1;
	double PMU_IK1 = 0;
	double nonPMU_IK1 = 0;
	
	// SK channel (3-state model)
	double SF_SK;
	double kc2c1SK;
	double kc1c2SK;
	double kc1oSK;
	double koc1SK;
	double ISK = 0;
		
	int	   SF_NoSK;						// number of SK channel (channel density is fixed)	
	DeepArray<int> iStateSK;					// 0: close2, 1:close1, 2:open
	DeepArray<double> dTimePredictSK;				// current time
	DeepArray<double> dStayingTimePredictSK;		// staying time
	sixelement PoSK;					// C2, C1, O
	sixelement SK_Current(sixelement dPoSK, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dCai);
	sixelement initializeSTOCHASTICSK(sixelement dPoSK, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dCai);
	void initializeCurrentTimeSTOCHASTICSK(int iChannelNumber, double* dTimePredict);		// initialize only current time
	
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
	double NaK_Current(double dNai, double dKi);
	double SF_NaK;
	double YNaK;
	double INaK;
	double fractionPMU_INaK = 0.1;
	double PMU_INaK = 0;
	double nonPMU_INaK = 0;

	// PMCA	
	void PMCA();
	double SF_PMCA;
	double IPMCA;	

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

	// background Ca-channel
	void Cab_Current();
	double SF_Cab;
	double ICab;		

	// background Cl-channel
	void Clb_Current();
	double SF_Clb;
	double IClb;

	// AutonomousInward_Current
	void AutonomousInward_Current();
	double SF_AutonomousInward;
	double IAutonomousInward;
	double IAutonomousInward_K;
	double IAutonomousInward_Na;
	double IAutonomousInward_Ca;
	
	// ----------------- Intracellular ion flux ------------------------//
	// IP3 receptor channel (4-state model)
	double SF_IP3R;	
	double powerOIP3R;
	double powerOIP3R0;				// all subnits
	double kucIP3R;
	double kcuIP3R;
	double kcoIP3R;
	double kocIP3R;
	double kobIP3R;
	double kboIP3R;	// rate constants
	double IIP3R;
	double fractionPMU_IIP3R = 0.3;
	double PMU_IP3R = 0;
	double nonPMU_IP3R = 0;
		
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
	double poweraOIP3R0;						// all subnits
	int	   SF_NoIP3R;							// number of IP3receptor channel (channel density is fixed)	
	DeepArray<int> iStateIP3R;							// 0: close, 1:open, 2:unavailable, 3: blocked (Ca-inhibited)
	DeepArray<double> dTimePredictIP3R;					// current time
	DeepArray<double> dStayingTimePredictIP3R;			// staying time
	int	   SF_NoIP3pmuR;						// number of IP3receptor channel (channel density is fixed)	
	DeepArray<int> iStateIP3pmuR;						// 0: close, 1:open, 2:unavailable, 3: blocked (Ca-inhibited)
	DeepArray<double> dTimePredictIP3pmuR;				// current time
	DeepArray<double> dStayingTimePredictIP3pmuR;			// staying time
	sixelement PoIP3R, PoIP3pmuR;				// C, B, U, O
	sixelement IP3R_Current(sixelement dPoIP3R, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dIP3, double dCai);
	sixelement initializeSTOCHASTICIP3R(sixelement dPoIP3R, int iChannelNumber, int* iState, double* dStayingTimePredict, double* dTimePredict, double dIP3, double dCai);	// initialize stochastic process (IP3 receptor channel)
	void initializeCurrentTimeSTOCHASTICIP3R(int iChannelNumber, double* dTimePredict);		// initialize only current time

	void InitializeStochasticChannels();


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
	bool SR_Concentration();

	// ion concentration
	bool Ion_Concentration();
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
	bool Ca_Buffer();
	double CMDNf;
	double CMDNPMU;
	double CSQNf;
	double TRPNf;

	// IP3 metabolism
	bool IP3metabolism();
	double IP3;
	double PIP2;
	double IP3total;
	// [headless] agonist (neural) drive: PIP2->IP3 rate is multiplied by (1+agonistScale)
	double agonistScale = 0.0;
	// [modified PDGFa] switches toward the literature SK3/PDGFRa+ phenotype (all off = original model)
	bool skUseEK = false;		// SK reversal = EK (Nernst) instead of fixed -70 mV
	bool skVoltIndep = false;	// SK3 activation depends on Ca only (Vm term evaluated at 0 mV)
	// [modified PDGFa] purinergic (P2Y1 -> PLC -> IP3) drive: an agonist-evoked IP3 signal added to the IP3 seen by IP3R.
	// Kept separate from the basal PIP2/IP3/IP4 metabolism, whose PIP2 pool is depleted at steady state.
	double ip3AgTarget = 0.0;	// mM, set by the agonist
	double ip3Ag = 0.0;			// mM, first-order follower of ip3AgTarget
	double tauIP3Ag = 300.0;	// ms
	// [modified PDGFa] Ca2+ sensed by SK: 0 = original (SOC-entry nanodomain + bulk, half-activation 25 uM, Hill 2)
	//                                    1 = submembrane (PMU) Ca2+, which receives IP3R release; SK3-like Kd and Hill
	int skCaSource = 0;
	double skKd = 4.0e-4;		// mM (0.4 uM)
	double skHill = 4.0;
	bool skDeterministic = false;	// [modified PDGFa] mean-field SK gate instead of the 20000-channel stochastic scheme
	double tauSKdet = 10.0;			// ms, SK3-like activation time constant
	double poSKdet = 0.0;
	bool gapCarriedByK = false;	// [modified PDGFa] charge the external (gap-junction) current to the K+ balance
	
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
	bool RunCellModel(bool AP, double InjectAmplitude);
	void Membrane_Potential(bool AP, double Iinj);
	void StoreVariables();
	void RetrieveVariables();
	void SaveCellParamters();
	
	// for monitoring
	double testVar1;
	double testVar2;
	double RatioCaLtoSR;
	double RatioNaCatoSR;
	double monitorVar1;
	double monitorVar2;
	double dAllParameters[100];

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

	// from mitochondria
	double NaFromMito;
	double CaFromMito;
	double KFromMito;
	double ATPFromMito;
	double ADPFromMito;	

// Operations
public:
	CCellPart();
	explicit CCellPart(unsigned seed);
	// [headless] all members are value types (DeepArray deep-copies the
	// stochastic channel arrays, std::mt19937 copies its state), so the
	// defaulted copy yields a fully independent cell. Call Reseed() on the
	// copy to make its stochastic trajectory diverge.
	CCellPart(const CCellPart& o) = default;
	CCellPart& operator=(const CCellPart& o) = default;
	~CCellPart() = default;

	// [headless] per-object RNG replacing global rand()/srand(time(NULL)).
	// Rand() returns ints uniformly in [0, RAND_MAX], like the CRT rand().
	unsigned seed = 1;
	std::mt19937 rng;
	void Reseed(unsigned s) { seed = s; rng.seed(s); }
	inline int Rand() { return int(rng() % (unsigned(RAND_MAX) + 1u)); }
};
