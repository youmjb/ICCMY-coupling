#pragma once
#define _USE_MATH_DEFINES // for C
#include <math.h>
typedef struct {
	double a1, a2, i1, i2, INa, IK, ICa, Itot;
} ode1;
typedef struct {
	double a, o, supf, deep;
} odesr;
typedef struct {
	double i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
} matin;
typedef struct {
	double c0, c1, c2, c3;
} matc;
class mCSMC
{
public:
	const double RTF = 26.7257;							 // RT/F (mV)
	const double RTF2 = 26.7257*0.5;							 // RT/F (mV), for Ca
	const double FRT = 0.037417;							// F/RT
	const double RT = 8.314 * 310.0;						// RT
	const double FFRT = 96.4867 * 96.4867 / 8.3143 / 310.0;
	const double F = 96.4867;   // C/mmol
	const double Cm = 26.4; //80 pF, Cell capacitance
	const double Vcell = 5e3; // 5 pL = 5e3 fL, cell volume
	const double betaCa = 1.0 / 41.0;   // calcium buffering power

	const double KmTRPN = 0.0005;		//DISSOCIATION CONSTANT OF TROPONIN Luo & Rudy, Circ Res,74:1071-1096,1994	
	const double KmCMDN = 0.00238;		//DISSOCIATION CONSTANT OF CALMODULIN Luo & Rudy, Circ Res,74:1071-1096,1994
	const double KmCSQN = 0.8;		//DISSOCIATION CALSEQUESTRIN OF TROPONIN Luo & Rudy, Circ Res,74:1071-1096,1994
	const double TRPNt = 0.049;		//TOTAL CONCENTRATION OF TROPONIN, mM/L
	const double CMDNt = 0.0005;		//TOTAL CONCENTRATIONT OF CALMODULIN, mM/L
	const double CSQNt = 10.0;		//TOTAL CONCENTRATIONT OF CALSEQUESTRIN, mM/L
	double TRPNf = 0;		//FREE CONCENTRATION OF TROPONIN, mM/L  necessary to start calculation "buffer"
	double CMDNf = 0;		//FREE CONCENTRATION OF CALMODULIN, mM/L	
	double CSQNf = 0;		//FREE CONCENTRATION OF CALMODULIN, mM/L	

	// variables
	double Cae = 2; //2mM, calcium concentration
	double Nae = 140; //140mM, sodium concentration
	double Ke = 5.4; // 6mM, K concentration
	double Caitotal = 0.0035; // 0.0035
	double Cai = 3.7920519569088907e-05; //100nM, calcium concentration
	double Casupf = 0.0035;
	double Caisupf = 3.7920519569088907e-05;
	double Cadeep = 0.0035;
	double Caideep = 3.7920519569088907e-05;
	double Carel = 13.1;
	double Cairel = 4.4662920295621760;
	double Caup = 4.61;                  // 4.61, initial Ca concentration in uptake pool (C[4][6][L])
	const double Caupmax = 5.0;			// maximum Ca concentration in uptake site (Caupmax)

	double Nai = 10; //10mM, sodium concentration
	double Ki = 140; // 140mM, K concentration

	// temperature coefficient
	double Q10Ca = 1.25; // 1.25, Q10 for calcium currents
	double Q10Na = 1.7; // 1.7, Q10 for sodium currents
	double Q10K = 1.44; // 1.44, Q10 for potassium currents

	// Gap junction current
	double Igap = 0;

	// stimulus current
	double Istim = 0;

	// Integration time step
	double delT = 0.01;

	// Voltage
	double Vm = -50.0;

	double time = 0;

	ode1 IRK(ode1 input, double voltage, double Ex, double Xe, double Xi, double dt);
	ode1 IRK1 = { 0,0,1,1,0,0,0,0 };
	ode1 IKfi(ode1 input, double voltage, double Ex, double dt);
	ode1 IKfi1 = { 0,0,1,1,0,0,0,0 };
	ode1 ICaL(ode1 input, double voltage, double CFCa, double CFK, double CFNa, double Xi, double dt);
	ode1 ICaL1 = { 0,0,1,1,0,0,0,0 };
	ode1 ICaT(ode1 input, double voltage, double Ex, double dt);
	ode1 ICaT1 = { 0,0,1,1,0,0,0,0 };
	ode1 IKCa(ode1 input, double voltage, double CFNa, double CFK, double Xi, double dt);
	ode1 IKCa1 = { 0,0,1,1,0,0,0,0 };
	ode1 INCX(ode1 input, double voltage, double X1e, double X1i, double X2e, double X2supf, double dt);
	ode1 INCX1 = { 0,0,1,1,0,0,0,0 };
	ode1 INaK(ode1 input, double voltage, double X1);
	ode1 INaK1 = { 0,0,1,1,0,0,0,0 };
	ode1 INSLC(ode1 input, double voltage, double CFNa, double CFK);
	ode1 INSLC1 = { 0,0,1,1,0,0,0,0 };
	odesr IRelRyr(odesr input, double xrel, double xfree, double xsupf, double xdeep, double ical, double dt);
	odesr IRelRyr1 = { 0,1,0,0 };
	odesr Iup(odesr input, double xup, double xsupf, double xdeep, double dt);
	odesr Iup1 = { 0,1,0,0 };


	matc getfourmatrix(matin i1, matc i2);

	void CaBuffer();

	double KCF, NaCF, CaCF;
	
	double sf1= 1.0;

	bool baICC = false;

	void InitializeVar();
	// Integrate model
	void IntegrateCurrents(int celltype);
	// Calculate Na, K, Ca ions
	void CalculateIon();
	// integrate all components
	void Integrate(int celltype, double ct, bool bICC);
};

