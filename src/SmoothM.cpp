#include "pch.h"
#include "SmoothM.h"
odes SmoothM::ICaL(odes input, double sf, double voltage, double Xe, double Xi, double dt)
{
	odes output = input;

	// activation
	double alfa = 13 * (voltage + 15) / (1 - exp(-(voltage + 15) / 3));
	double beta = 31 * exp(-0.04 * (voltage + 15));
	double a1inf = alfa / (alfa + beta);
	double a1tau = tauCaLScale * 1e3 / (alfa + beta);	// tauCaLScale (Coupling copy): < 1 speeds activation (original 17 ms at -30 mV)
	output.a1 = a1inf - (a1inf - input.a1) * exp(-dt / a1tau);

	// fast inactivation
	double i1inf = 1.0 / (1.0 + exp((voltage + 40.5) / 13.1));
	double i1tau = 28 / (1 + exp(-(voltage + 20) / 1));
	output.i1 = i1inf - (i1inf - input.i1) * exp(-dt / i1tau);

	// slow inactivation
	double i2inf = i1inf;
	double i2tau =300 / (1 + exp(-(voltage + 60) / 2));
	output.i2 = i2inf - (i2inf - input.i2) * exp(-dt / i2tau);

	// Ca-dependent inactivation
	double KCa = 9;           // 0.1, 1, 10 uM
	double hCa = 1 / (1 + Xi / KCa);

	// current amplitude
	double GCaL = 65;									// 1800 nS, 
	double ECaK = 87;									   // (60, 40) mV, reversal potential
	output.ICa = sf * GCaL * output.a1 * output.i1 * output.i2 * hCa * (voltage - ECaK);
	output.Itot = output.ICa;

	return output;
}
odes SmoothM::IKs(odes input, int type, double sf, double voltage, double Xe, double Xi, double X2e, double X2i, double dt)
{
	odes output = input;

	// activation
	double alfa = -1.126 * (voltage + 8.821) / (exp(-0.386 * (voltage + 8.821)) - 1);
	double beta = 17.02 * exp(-0.022 * (voltage + 60));
	double a1inf = alfa / (alfa + beta);
	double a1tau = 1e3 / (alfa + beta);
	output.a1 = a1inf - (a1inf - input.a1) * exp(-dt / a1tau);

	// current amplitude
	double GKs = 3.1;                 // 540 nS
	output.IK = sf*GKs * pow(output.a1, 2) * (voltage - EK);
	output.Itot = output.IK;

	return output;
}
odes SmoothM::IKCa(odes input, double sf, double voltage, double Xi)
{
	odes output = input;
	double kvca = -65;                          // (-140,-50) mV
	double kn = kvca * Xi + 119.44;
	double skca = 14.372 * Xi + 11.566;
	double hkca = 1 / (1 + exp((kn - voltage) / skca));
	double gkca = 0.7;                   // 480 nS
	output.IK = sf*gkca * hkca * (voltage - EK);
	output.Itot = output.IK;

	return output;
}
odes SmoothM::Ileak(odes input, double sf, double voltage)
{
	odes output = input;
	double gleak = 5e-3;                     // 4~40 uS/cm2,  20 uS/uF, 20 pS/pF, 20e-3 nS/pF
	double Eleak = -46.3;                       // (-60, -20) mV, automaticity factor
	output.Itot = sf*gleak *Cm* (voltage-Eleak);     // pA unit
	return output;
}
// Initialize all variables
void SmoothM::InitializeVar()
{
	ICaL1 = { 0,0,1,1,0,0,0,0 };
	IKs1 = { 0,0,1,1,0,0,0,0 };
	IKCa1 = { 0,0,1,1,0,0,0,0 };
	Ileak1 = { 0,0,1,1,0,0,0,0 };

	sfiks1 = 1.0;
	sfical1 = 1.0;
	sfikca1 = 1.0;
	sfileak1 = 1.0;
	
	// Gap junction current
	Igap = 0;

	// stimulus current
	Istim = 0;

	// Voltage
	Vm = -50;

	Ke = 5.40e+00;
	Ki = 1.09e+02;
	Cae = 2.50e+00;
	Cai = 0.1;
	Cacore = 0.1;
	Nai = 7.92e+00;
	
	// equilibrium potential
	ENa = 75;
	EK = -85;
	ECa = 100;

	// Integration time step
	delT = 0.01;
}
// Integrate model, celltype: 0: epicardial type, 1: M caridiac, 1: endocardial type
void SmoothM::IntegrateCurrents(int celltype)
{
	/*if (celltype == 0) Vc = 50 * 100 * 20 * 8 / 2;
	else if(celltype == 0) Vc = 20 * 100 * 20 * 8 / 2;
	else Vc = 70 * 100 * 20 * 8 / 2;*/

	// equilibrium potential
	ENa = RTF * log(Nae / Nai);
	EK = RTF * log(Ke / Ki);
	ECa = 0.5 * RTF * log(Cae / (1e3*Cai));
	ICaL1 = ICaL(ICaL1, sfical1, Vm, Cae, Cai, delT);
	IKs1 = IKs(IKs1, celltype, sfiks1, Vm, Ke, Ki, Nae, Nai, delT);
	IKCa1 = IKCa(IKCa1, sfikca1, Vm, Cai);
	Ileak1 = Ileak(Ileak1, sfileak1, Vm);
}

// Calculate Na, K, Ca ions
void SmoothM::CalculateIon()
{
	// calculate cytoplasmic total calcium change
	// pA/cm3/(C/mmole) = 1e3*pA/L/C*mmole = 1e3*1e-12*mM/s = 1e-9 mM/s = 1e-9 uM/ms
	//double dCasubdt = -1e-9*(ICaL1.ICa) / 2 / subV / F;			
	//double dCacoredt = 1e-3 * Cai * CompLoss / SACC * CV / APVOL;          // uM x cm2/s /cm2 = uM/s = 1e-3 uM/ms
	//Cacore += dCacoredt * delT;	
	//Cai += (dCasubdt - 1000*dCacoredt * APVOL / CV) * delT;

	// A model for Ca-dependent Ca inactivation by Standen and Stanfield
	double djCadt = -1e-9 * (ICaL1.ICa) / 2 / subV / F;
	double djCadifdt = 1e-1 * Cai * CompLoss / SACC;
	double djCaleakdt = 1e-1 * Cacore * CompLoss / SACC;
	Cai += (djCadt - djCadifdt + djCaleakdt) * delT;

	// voltage
	double dvmdt =  -(
		ICaL1.Itot 
		+ IKs1.Itot 
		+ IKCa1.Itot 
		+ Ileak1.Itot
		+ Istim + Igap)/Cm;       // mV/ms = pA/pF, Q = VC , mV/ms = mA/mF = mQ/ms/mF = mV*mF/ms/mF = mV/ms 
	Vm += dvmdt * delT;
}
void SmoothM::Integrate(int celltype)
{
	IntegrateCurrents(celltype);
	CalculateIon();
}


