#pragma once
#define _USE_MATH_DEFINES // for C
#include <math.h>
typedef struct {
	double a1, a2, i1, i2, INa, IK, ICa, Itot;
} odes;

class SmoothM
{
public:
	
	const double RTF = 26.7257;							 // RT/F (mV)
	const double FRT = 0.037417;							// F/RT
	const double RT = 8.314 * 310.0;						// RT
	const double FFRT = 96.4867 * 96.4867 / 8.3143 / 310.0;

	// cell geometry
	double subV = M_PI * (6e-4 * 6e-4 - 5e-4 * 5e-4) * 200e-4;            // 0.2~1.2 um thickness, 
	double CV = 200e-4 * M_PI * 6e-4 * 6e-4;							            // 200um x 12um = 200 x phi x (12/2)^2, volume of core compartment
	double APVOL = 1e4 * CV;								                            // 1000~100000, apparent volume of core compartment which is 1000-100000 x greater than its physical volume
	double CompLoss = 7e-6;                                                          // 1e-6 cm2 s-1, slower than the rate of Ca-diffusion through water (7e-6 cm2 s-1)
	double SACC = M_PI*200e-4*12e-4;                                            // um2, surface area of core compartment, um = 1e-4 cm, 7.54e-5 cm2 = 75.4 pF
	double Cai = 0.1;                                // 0.1 uM, micromole unit
	double Cacore = 0.1;                           // 0.1 uM, micromole unit
	const double F = 96.4867;                     // 96.4867 C/mmol
	const double Cm = 50;							// cell capacitance, pF
	// Initialize all variables
	void InitializeVar();

	// Current (function, current amplitude, conductance)
	odes SmoothM::ICaL(odes input, double sf, double voltage, double Xe, double Xi, double dt);
	odes SmoothM::IKs(odes input, int type, double sf, double voltage, double Xe, double Xi, double X2e, double X2i, double dt);
	odes SmoothM::IKCa(odes input, double sf, double voltage, double Xi);
	odes SmoothM::Ileak(odes input, double sf, double voltage);

	odes ICaL1 = { 0,0,1,1,0,0,0,0 };
	odes IKs1 = { 0,0,1,1,0,0,0,0 };
	odes IKCa1 = { 0,0,1,1,0,0,0,0 };
	odes Ileak1 = { 0,0,1,1,0,0,0,0 };

	double sfical1 = 1.0;
	double sfiks1 = 1.0;
	double sfikca1 = 1.0;
	double sfileak1 = 1.0;
	double tauCaLScale = 1.0;	// ICaL activation time-constant scale (Coupling copy only)

	// Gap junction current
	double Igap = 0;

	// stimulus current
	double Istim = 0;

	// Voltage
	double Vm = -80.0;

	// Ion
	double Ke = 5;
	double Ki = 140;
	double Cae = 1.8;
	double Nae = 140.0;
	double Nai = 7.68;

	// equilibrium potential
	double ENa = 75;
	double EK = -75;
	double ECa = 100;

	// Integration time step
	double delT = 0.01;

	// Integrate model
	void IntegrateCurrents(int celltype);
	// Calculate Na, K, Ca ions
	void CalculateIon();
	// integrate all components
	void Integrate(int celltype);
};

