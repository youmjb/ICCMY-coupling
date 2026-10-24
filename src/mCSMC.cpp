#include "pch.h"
#include "mCSMC.h"

ode1 mCSMC::IRK(ode1 input, double voltage, double Ex, double Xe, double Xi, double dt)
{
	ode1 output = input;

	double minf = 1.0 / (1.0 + 0.1 * exp((voltage - Ex) / 15.0) + 0.048 * exp((voltage - Ex) / 7.0));     //7APPROXIMATED FROM STEADY STATE i-V
	double rtau = 2.5 / (4.7 * exp((voltage - Ex) / 28.7)) + 2.5 / (6.0 * exp((Ex - voltage) / 25.8));      //1/TAU
	double mtau = 1.0 / rtau;
	output.a1 = minf - (minf - input.a1) * exp(-dt / mtau);

	double lambda = 26.0 * exp((Ex - voltage - 6.0) / 52.0);     //Lamda, 120                               'adjust half block = -50 mV
	const double Mgi = 1.0;
	double myu = Mgi * 1.85 * exp((voltage - Ex - 6.0) / 15.0);  //MYU, M^-1,msec^-1  single channel data maximum slope
	double hinf = lambda / (lambda + myu);

	double GIRK1 = 0.3024;
	output.IK = GIRK1 * (voltage - Ex) * pow((Xe / 5.4), 0.62) * output.a1 * hinf;
	output.Itot = output.IK;

	return output;
}
ode1 mCSMC::IKfi(ode1 input, double voltage, double Ex, double dt)
{
	ode1 output = input;

	double ACh = 0.2;
	double mf = 0;
	if (ACh == 0.0) mf = 0.165 * 0.001 * 1.5;      //spontaneous opening   /msec"
	else mf = 12.32 / (1.0 + 4.2 / ACh) * 0.001;
	double mb = 10.0 * exp(0.0133 * (voltage + 40.0)) * 0.001; //beta of IKACh
	double minf = mf / (mf + mb);
	double mtau = 1.0 / (mf + mb);
	output.a1 = minf - (minf - input.a1) * exp(-dt / mtau);
	
	output.IK = 0.2700 / 2.0 * 26.4 * output.a1 * (voltage - Ex) / (1.0 + exp((voltage + 20.0) / 20.0));
	output.Itot = output.IK;

	return output;
}
ode1 mCSMC::ICaL(ode1 input, double voltage, double CFCa, double CFK, double CFNa, double Xi, double dt)
{
	ode1 output;

	// forward rate constant of activation gate
	double mSpeed = 0.15; //0.04;
	double hSpeed = 0.005; //0.005;
	double df = 0;
	if (voltage == -30.0) df = mSpeed * 0.02792845;
	else	df = mSpeed * 0.01045 * (voltage + 30.0) / (1.0 - exp(-(voltage + 30.0) / 2.5));         //+ .03125 * E / (1 - EXP(-E / 4.8)) 'alpha, activation parameter for L-type Ca channel
	// backward rate constant of activation gate
	double db = 0;
	if (voltage == 0.0)		db = mSpeed * 0.01052639;
	else		db = mSpeed * 0.00421 * voltage / (exp(voltage / 2.5) - 1.0);                 //beta
	double dinf = df / (df + db);
	double dtau = 1 / (df + db);
	output.a1 = dinf - (dinf - input.a1) * exp(-dt / dtau);

	// forward rate constant of inactivation gate
	double ff = 0;
	if (voltage == -34.0)		ff = hSpeed * 0.00200021;
	else		ff = hSpeed * 0.000355 * (voltage + 34.0) / (exp((voltage + 34.0) / 5.633) - 1.0);
	// backward rate constant of inactivation gate
	double fb = hSpeed * Xi * pow(10.0, 4.5) * 0.0007 * (voltage + 64.0) / (1.0 + exp(-(voltage + 44.0) / 4.16));       //beta
	
	double finf = ff / (ff + fb);
	double ftau = 1.0 / (ff + fb);
	output.i1 = finf - (finf - input.i1) * exp(-dt / ftau);
	double PCa = 96.0;
	output.ICa = 2.0 * PCa * output.a1 * output.i1 * CFCa;
	output.IK = 0.002 * PCa * output.a1 * output.i1 * CFK;
	output.INa = 0.0001 * PCa * output.a1 * output.i1 * CFNa;
	output.Itot = output.ICa + output.INa + output.IK;
	return output;
}
ode1 mCSMC::ICaT(ode1 input, double voltage, double Ex, double dt)
{
	ode1 output = input;

	double gcat = 0*1.72; //1.72nS, maximum conductance
	double vdcat = -34; //-48mV, half-activation potential
	double kdcat = 5.0; //6.5mV, slope potential for activation
	double vfcat = -29; //-73mV, half-inactivation potential
	double kfcat = 4; //6mV, slope potential of inactivation

	double dcatinf = 1 / (1 + exp(-(voltage - vdcat) / kdcat));
	double taudcat = 0.04439 - 0.02824 / (1 + exp(-(voltage + 26.156) / 5.536));
	taudcat *= 1e3;
	output.a1 = dcatinf - (dcatinf - input.a1) * exp(-dt / taudcat);

	double fcatinf = 1 / (1 + exp((voltage - vfcat) / kfcat));
	double taufcat = 0.034756 - 0.02762 / (1 + exp(-(voltage + 16.786) / 10.58));
	taufcat *= 1e3;
	output.i1 = fcatinf - (fcatinf - input.i1) * exp(-dt / taufcat);

	output.ICa = gcat * output.a1 * output.i1 * (voltage - Ex);
	output.Itot = output.ICa;

	return output;
}
ode1 mCSMC::IKCa(ode1 input, double voltage, double CFNa, double CFK, double Xi, double dt)
{
	ode1 output;

	double minf = 1.0 / (1.0 + exp((voltage + 25.0) / (-6.0))) / (1 + pow(0.0002 / Xi, 2.0));        //steady-state open of activation
	double mtau = 9.0 / (1.0 + exp((voltage + 5.0) / 12.0)) + 0.5;      //time constant of activation gate
	output.a1 = minf - (minf - input.a1) * exp(-dt / mtau);

	double finf = 0.4 / (1.0 + exp((voltage + 0.0) / 10.0)) + 0.6;   //steady-state open of inactivation gate, Mathematical model
	double ftau = 47.0 / (1.0 + exp((voltage + 60.0) / 10.0)) + 300.0;    //time constant of inactivation gate
	output.i1 = finf - (finf - input.i1) * exp(-dt / ftau);

	double mCaDep;
	mCaDep = 1500 / (1 + 0.005 / Xi) + 500 / (1 + 0.0005 / Xi);

	output.ICa = 0;
	double GIKCaK = 0.006;
	double RGIKCaNa = 0.09;
	output.IK = mCaDep * output.a1 * output.i1 * GIKCaK * CFK;
	output.INa = mCaDep * output.a1 * output.i1 * GIKCaK * RGIKCaNa * CFNa;
	output.Itot = output.INa + output.IK + output.ICa;

	return output;
}
ode1 mCSMC::INCX(ode1 input, double voltage, double X1e, double X1i, double X2e, double X2supf, double dt)   // Nae, Nai, Cae, Casupf
{
	ode1 output = input;
	double KNAI = 670;				// dissociation constant for Na at the inner side of the membrane (NaCa exchanger)
	double KNAO = 670000;			// dissociation constant for Na at the outer side of the membrane (NaCa exchanger)
	double KCaI = 0.00138;			// dissociation constant for Ca at the inner side of the membrane (NaCa exchanger) 
	double KCaO = 1.38;			// dissociation constant for Ca at the outer side of the membrane (NaCa exchanger)
	double BasalRate = 250.0;		// basal rate constant for the Na-Ca exchanger (BasalRate)
	double PNAO = 1.0 / (1.0 + ((1.0 + X2e / KCaO) / (pow(X1e, 3.0) / KNAO))); //'PROBABILITY OF Na BINDING
	double PNAI = 1.0 / (1.0 + ((1.0 + X2supf / KCaI) / (pow(X1i, 3.0) / KNAI))); //'PROBABILITY OF Na BINDING
	double PCAO = 1.0 / (1.0 + ((1.0 + pow(X1e, 3.0) / KNAO) / (X2e / KCaO))); //'PROBABILITY OF Na BINDING
	double PCAI = 1.0 / (1.0 + ((1.0 + pow(X1i, 3.0) / KNAI) / (X2supf / KCaI))); //'PROBABILITY OF Na BINDING

	double Partition_v = 0.32;		// Partition_v of Na-K pump
	double K30 = exp(Partition_v * voltage / RTF) * PNAI;
	double K03 = exp((Partition_v - 1.0) * voltage / RTF) * PNAO;
	double mb = (K30 + PCAI);    //RATE FROM PXI
	double mf = (K03 + PCAO);   //RATE TO   PXI
	double minf = mf / (mf + mb);
	double mtau = 1.0 / (mf + mb);
	output.a1 = minf - (minf - input.a1) * exp(-dt / mtau);

	double INaCaMax = 8.45856;
	output.Itot = INaCaMax * (output.a1 * K30 - (1.0 - output.a1) * K03);
	output.INa = 3.0 * output.Itot;
	output.ICa = -2.0 * output.Itot;

	return output;
}
ode1 mCSMC::INaK(ode1 input, double voltage, double X1)
{
	ode1 output = input;

	double RectifyingFactor = 0.2*50;
	double INaKmax = 2.965248;
	double KmNa = 11.0;
	output.Itot = (1.0 / (1.0 + pow(KmNa / X1, 1.36))) * (1.0 - pow((voltage - RectifyingFactor) / 125.0 / 2.0, 2.0)) * INaKmax;
	output.INa = 3.0 * output.Itot;
	output.IK = -2.0 * output.Itot;

	return output;
}
ode1 mCSMC::INSLC(ode1 input, double voltage, double CFNa, double CFK)
{
	ode1 output = input;
	double GILNSCNa = 0.7*0.030294;

	output.INa = GILNSCNa * CFNa;
	output.IK = 1.73*GILNSCNa * CFK;
	output.Itot = output.INa + output.IK;

	return output;
}
// xrel = Cafrel, xfree = Caf, xsupf = Cafsupf, ip3 = IP3
odesr mCSMC::IRelRyr(odesr input, double xrel, double xfree, double xsupf, double xdeep, double ical, double dt)
{
	odesr output = input;

	// inactivation to close
	double a = 2.0 * 0.000377 * pow(xrel / 1.0, 2.0);		// 0.000377*pow([Ca]rel, 2.0)
	// close to inactivation
	double b = 0.000849;									// 0.000849
	// close to open
	double c = 280000 * pow(xfree, 2.0);					// 280000*pow([Ca]cyt,2.0)
	// open to inactivation
	double d = 2 * 0.08 / (1 + 0.36 / xrel);					// 0.08/(1+0.36/[Ca]rel)

	double term1 = a + b + c + d;
	double term2 = a * c + b * d + c * d + a * d;

	double y0 = input.a;  double x0 = 1.0 - input.a - input.o;

	double alfa = 0; double beta = 0; double r1 = 0; double r2 = 0;
	matc mout = { 0,0,0,0 };
	if (term1 * term1 - 4 * term2 > 0)
	{
		alfa = 0.5 * (-term1 + sqrt(term1 * term1 - 4 * term2));
		beta = 0.5 * (-term1 - sqrt(term1 * term1 - 4 * term2));
		r2 = (b * d + c * d) / term2;
		r1 = a * d / term2;
		
		double k1 = alfa + b + c;
		double k2 = -a;
		double k3 = beta + b + c;
		double k4 = -a;
		double A1 = y0 - r1;
		double B1 = x0 - r2;
		double C2 = (k1 * A1 + k2 * B1) / (k1 - k3 * k2 / k4);
		double C1 = A1 - C2;
		double C3 = -k1 * C1 / k2;
		double C4 = -k3 * C2 / k4;
		double PU = C3 * exp(alfa * dt) + C4 * exp(beta * dt) + r2;
		output.a = C1 * exp(alfa * dt) + C2 * exp(beta * dt) + r1;
		output.o = 1.0 - PU - output.a;
	}
	else
	{
		alfa = 0.5 * (-term1);
		beta = 0.5 * sqrt(fabs(term1 * term1 - 4 * term2));
		r2 = (b * d + c * d) / term2;
		r1 = a * d / term2;
		
		double k1 = alfa + b + c;
		double k2 = beta;
		double k3 = -a;
		double k4 = -beta;
		double k5 = alfa + b + c;
		double k6 = -a;
		double A1 = y0 - r1;
		double B1 = x0 - r2;
		double C1 = A1;
		double C3 = B1;
		double C2 = -(k1 * C1 + k3 * C3) / k2;
		double C4 = -(k4 * C1 + k5 * C2) / k6;

		output.a = exp(alfa * dt) * (C1 * cos(beta * dt) + C2 * sin(beta * dt)) + r1;
		double PU = exp(alfa * dt) * (C3 * cos(beta * dt) + C4 * sin(beta * dt)) + r2;
		output.o = 1.0 - output.a - PU;
	}
	
	// current amplitude
	double FARyR = 2071.2079153847521;
	double supfSR = 0.6;

	output.supf = supfSR * FARyR * (xrel - xsupf) * output.o;				//global or local Ca gradient?
	output.deep = (1 - supfSR) * FARyR * (xrel - xdeep) * output.o;		//global or local Ca gradient?

	return output;
}
matc mCSMC::getfourmatrix(matin i1, matc i2)
{
	matc output;
	//x + 3y - 2z + 5w = 7
	//3x + 5y + 6z + 7w = 9
	//2x + 4y + 3z + 8w = 5
	//4x + 11y + 2z + 8w = 3

	/*double mx[4][5] = { 
	 { i1.i0,i1.i1,i1.i2,i1.i3,i2.c0 },
	 { i1.i4,i1.i5,i1.i6,i1.i7,i2.c1 },
	 { i1.i8,i1.i9,i1.i10,i1.i11,i2.c2 },
	 { i1.i12,i1.i13,i1.i14,i1.i15,i2.c3 } 
	};*/

	double mx[4][5] = {
	 {1, 1, 0, 0,0.074515042074067006 },
	 { 0,0,1,1,-0.072004604133259736 },
	 { -0.016136804395722015,0,-0.016087375908733433,0,0 },
	 { 0,-0.14713057467190940,0,-0.016087375908733433,0 }
	};



	double tmp = 0;
	int i, j, dl, count;
	int N = 4;
	int P = N + 1;
	for (i = 0; i < N; i++)
	{
		// set all element before the diagonal line as 0
		dl = i;
		count = 0;
		while (count < dl)      // proceed until the diagonal line appears
		{
			// get multiplying factor 
			tmp = mx[i][count];
			// make all the elements before diagonal line as 0
			for (j = 0; j < P; j++)    mx[i][j] -= tmp * mx[count][j];
			count += 1;
		}
		if (count == dl)
		{
			// set diagonal element as 1
			tmp = mx[i][count];
			for (j = 0; j < P; j++)    mx[i][j] /= tmp;
		}
	}
	for (i = N - 1; i >= 0; i--)
	{
		// 
		dl = i;
		count = N - 1;
		while (count > dl)
		{
			// get multiplying factor
			tmp = mx[i][count];
			// make all the elements after diagonal line as 0
			for (j = 0; j < P; j++)    mx[i][j] -= tmp * mx[count][j];
			count -= 1;
		}
	}

	// solution of the equation
	output.c0 = mx[0][4];
	output.c1 = mx[1][4];
	output.c2 = mx[2][4];
	output.c3 = mx[3][4];

	return output;
}
odesr mCSMC::Iup(odesr input, double xup, double xsupf, double xdeep, double dt)
{
	odesr output = input;

	double Km_Cao = 0.002;
	double Km_Cai = 0.08;
	double supfK1 = 0.01 / (1 + Km_Cai / xup);
	double supfK2 = 1.0 / (1 + Km_Cao / xsupf);
	double supfK3 = 1 / (1 + xup / Km_Cai);
	double supfK4 = 0.01 / (1 + xsupf / Km_Cao);
	double deepK1 = 0.01 / (1 + Km_Cai / xup);
	double deepK2 = 1.0 / (1 + Km_Cao / xdeep);
	double deepK3 = 1 / (1 + xup / Km_Cai);
	double deepK4 = 0.01 / (1 + xdeep / Km_Cao);
	double Aup = 2761610.5538463364;
	double supfSR = 0.6;
	
	double mf1 = supfK2 + supfK4;
	double mb1 = supfK1 + supfK3;
	double minf1 = mf1 / (mf1 + mb1);
	double mtau1 = 1.0 / (mf1 + mb1);
	output.a = minf1 - (minf1 - input.a) * exp(-dt / mtau1);

	double ff1 = deepK2 + deepK4;
	double fb1 = deepK1 + deepK3;
	double finf1 = ff1 / (ff1 + fb1);
	double ftau1 = 1.0 / (ff1 + fb1);
	output.o = finf1 - (finf1 - input.o) * exp(-dt / ftau1);

	output.supf = 0.0015 * Aup * supfSR * (supfK2 * (1 - output.a) - supfK1 * output.a);	//0.001
	output.deep = 0.003 * Aup * (1 - supfSR) * (deepK2 * (1 - output.o) - deepK1 * output.o);  //0.002	

	return output;
}
void mCSMC::CaBuffer()
{
	int N = 1;
	do
	{
		Cai = Caitotal / (1.0 + (TRPNf / KmTRPN) + (CMDNf / KmCMDN));
		TRPNf = TRPNt / (1.0 + (Cai / KmTRPN));
		CMDNf = CMDNt / (1.0 + (Cai / KmCMDN));

		N = N + 1;
	} while (N <= 10);
	N = 1;
	do
	{
		Cairel = Carel / (1.0 + (CSQNf / KmCSQN));
		CSQNf = CSQNt / (1.0 + (Cairel / KmCSQN));
		N = N + 1;
	} while (N < 10);

	Caisupf = Casupf / (1.0 + (TRPNf / KmTRPN) + (CMDNf / KmCMDN));
	Caideep = Cadeep / (1.0 + (TRPNf / KmTRPN) + (CMDNf / KmCMDN));
}
void mCSMC::InitializeVar()
{
	// variables
	Cae = 2; //2mM, calcium concentration
	Nae = 140; //140mM, sodium concentration
	Ke = 5.4; // 6mM, K concentration
	Caitotal = 0.0035; // 0.0035
	Cai = 3.7920519569088907e-05; //100nM, calcium concentration
	Casupf = 0.0035;
	Caisupf = 3.7920519569088907e-05;
	Cadeep = 0.0035;
	Caideep = 3.7920519569088907e-05;
	Carel = 13.1;
	Cairel = 4.4662920295621760;
	Caup = 4.61;                  // 4.61, initial Ca concentration in uptake pool 

	Nai = 10; //10mM, sodium concentration
	Ki = 140; // 140mM, K concentration

	// Calcium buffer
	TRPNf = 0;		//FREE CONCENTRATION OF TROPONIN, mM/L  necessary to start calculation "buffer"
	CMDNf = 0;		//FREE CONCENTRATION OF CALMODULIN, mM/L	
	CSQNf = 0;		//FREE CONCENTRATION OF CALMODULIN, mM/L	

	// temperature coefficient
	Q10Ca = 1.25; // 1.25, Q10 for calcium currents
	Q10Na = 1.7; // 1.7, Q10 for sodium currents
	Q10K = 1.44; // 1.44, Q10 for potassium currents

	// Gap junction current
	Igap = 0;

	// stimulus current
	Istim = 0;

	// Voltage
	Vm = -50.0;

	// Integration time step
	delT = 0.01;

	IRK1 = { 0.1143,0.9985,0.1549,1,0,0,0,0 };
	IKfi1 = { 0.1143,0.9985,0.1549,1,0,0,0,0 };
	ICaL1 = { 0.00438,0,0.7483,0.9935,0,0,0,0 };
	ICaT1 = { 0.4242,0,0.02114,1,0,0,0,0 };
	IKCa1 = { 0.02392,0,0.07126,1,0,0,0,0 };
	INCX1 = { 0,0,1,1,0,0,0,0 };
	INaK1 = { 0,0,1,1,0,0,0,0 };
	INSLC1 = { 0,0,1,1,0,0,0,0 };
	IRelRyr1 = { 1,0,0,0 };
}
void mCSMC::IntegrateCurrents(int celltype)
{
	// equilibrium potential
	double ENa = RTF * log(Nae / Nai);
	double EK = RTF * log(Ke / Ki);
	double ECa = 0.5 * RTF * log(Cae / Cai);

	// constant field for ions
	double RectifyingFactor = 50.0;
	// Ca2+
	if (Vm == RectifyingFactor) CaCF = 0.99906;
	else CaCF = ((Vm - RectifyingFactor) / RTF2) / (1 - exp(-(Vm - RectifyingFactor) / RTF2));
	CaCF = CaCF * (Cai * exp(RectifyingFactor / RTF2) - Cae * exp(-(Vm - RectifyingFactor) / RTF2));

	// K+
	if (Vm == RectifyingFactor) KCF = 0.99906;
	else KCF = ((Vm - RectifyingFactor) / RTF) / (1.0 - exp(-(Vm - RectifyingFactor) / RTF));
	KCF = KCF * (Ki * exp(RectifyingFactor / RTF) - Ke * exp(-(Vm - RectifyingFactor) / RTF));

	// Na+
	if (Vm == RectifyingFactor) NaCF = 0.99906;
	else	NaCF = ((Vm - RectifyingFactor) / RTF) / (1.0 - exp(-(Vm - RectifyingFactor) / RTF));
	NaCF = NaCF * (Nai * exp(RectifyingFactor / RTF) - Nae * exp(-(Vm - RectifyingFactor) / RTF));
	   
	IRK1 = IRK(IRK1, Vm, EK, Ke, Ki, delT);   // 
	IKfi1 = IKfi(IKfi1, Vm, EK, delT);   // 
	ICaL1 = ICaL(ICaL1, Vm, CaCF, KCF, NaCF, Cai, delT); // Cai
	ICaT1 = ICaT(ICaT1, Vm, ECa, delT);
	IKCa1 = IKCa(IKCa1, Vm, NaCF, KCF, Caisupf, delT);		// Caisupf
	INCX1 = INCX(INCX1, Vm, Nae, Nai, Cae, Caisupf, delT);   // Caisupf
	INaK1 = INaK(INaK1, Vm, Nai);
	INSLC1 = INSLC(INSLC1, Vm, NaCF, KCF);

	IRelRyr1 = IRelRyr(IRelRyr1, Cairel, Cai, Caisupf, Caideep, ICaL1.ICa, delT);
	Iup1 = Iup(Iup1, Caup, Caisupf, Caideep, delT);

}

// Calculate Na, K, Ca ions
void mCSMC::CalculateIon()
{
	double PMCA = 1.8 * 0.0012 * 2761610.5538463364 / (1 + pow(0.08 / (Caisupf), 2.0));		// [modified]
	double Ileaksupf = 0.6 * 0.45 * (Caup - Caisupf);
	double Ileakdeep = (1 - 0.6) * 0.45 * (Caup - Caideep);
	double Itr = 20.507009063215367 * (Caup - Cairel);
	double Vup = 28.621670111997314;
	double Vrel = 7.1554175279993286;

	// calculate uptake pool
	double dCaupdt = (Iup1.deep + Iup1.supf - Itr - Ileaksupf - Ileakdeep) / 2.0 / Vup / F;			//change at uptake site
	Caup = Caup + dCaupdt * delT;
	if (Caup > Caupmax) Caup = Caupmax;

	// calculate release pool
	double dCareldt = (Itr - IRelRyr1.deep - IRelRyr1.supf) / 2.0 / Vrel / F;	//change at releasing site
	Carel = Carel + dCareldt * delT;


	//-------- Determine Ca concentration at superficial or deep cytosol ------------///		
	double DiffConst = 0.5;		//diffusion constants
	//net current at superficial SR
	double IBySRsupf = -IRelRyr1.supf  + Iup1.supf - Ileaksupf;	
	//net current at deep SR
	double IBySRdeep = -IRelRyr1.deep  + Iup1.deep - Ileakdeep;

	// calculate Cai total
	double ICanet = (ICaL1.ICa
		+ INCX1.ICa
		+ PMCA);

	//superficial SR Ca concentration change
	double Vsupf = 0.1; double Vdeep = 0.9;
	double Vtotal = 715.54175279993285;
	double dCasupfdt = (-ICanet / 2.0 / F - DiffConst * (Casupf - Cadeep)) / (Vsupf * Vtotal) - IBySRsupf / 2.0 / F / (Vsupf * Vtotal);	
	Casupf = Casupf + dCasupfdt * delT;

	//deep SR Ca concentration change
	double dCadeepdt = DiffConst * (Casupf - Cadeep) / (Vdeep * Vtotal) - IBySRdeep / 2.0 / F / (Vdeep * Vtotal);	
	Cadeep = Cadeep + dCadeepdt * delT;
	
    // whole cytosolic Ca concentration change
	double dCaitotaldt = -(ICanet + IBySRsupf + IBySRdeep) / 2.0 / Vtotal / F;	//calculate ion concentration changes, d[I]/dt
	Caitotal = Caitotal + dCaitotaldt * delT;

	// calculate K
	double IKnet = (IRK1.IK
		+ ICaL1.IK
		+ IKCa1.IK
		+ IKfi1.IK
		+ INaK1.IK
		+ INSLC1.IK
		+ Istim + Igap);
	double dKdt = -IKnet / Vtotal / F;
	Ki = Ki + dKdt * delT;

	// calculate Na
	double INanet = (IKCa1.INa 
		+ INSLC1.INa
		+ INCX1.INa
		+ INaK1.INa
		+ ICaL1.INa);
	double dNadt = -INanet / Vtotal / F;
	Nai = Nai + dNadt * delT;
	   
	// calculate voltage
	double Iadd = Istim + Igap;
	
	double dvmdt = -(
		IRK1.Itot
		+ ICaL1.Itot
		+ IKCa1.Itot
		+ IKfi1.Itot
		+ INSLC1.Itot
		+ INCX1.Itot
		+ INaK1.Itot
		+ Iadd) / Cm;       // mV/ms = pA/pF, Q = VC , mV/ms = mA/mF = mQ/ms/mF = mV*mF/ms/mF = mV/ms 
	Vm += dvmdt * delT;
}
void mCSMC::Integrate(int celltype, double ct, bool bICC)
{
	time = ct;
	baICC = bICC;
	CaBuffer();
	IntegrateCurrents(celltype);
	CalculateIon();
}