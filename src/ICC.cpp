#include "pch.h"
#include "ICC.h"

void ICC::CaL_Current()
{
	double dAlfaAct = 0.0054;
	if (Vm != -30.0)	dAlfaAct = 0.002175 * (Vm + 30.0) / (1.0 - exp(-(Vm + 30.0) / 2.5));         // alpha, activation parameter for L-type Ca channel

	double dBetaAct = 0.0016;
	if (Vm != 0.0)	dBetaAct = 0.0006315 * Vm / (exp(Vm / 2.5) - 1.0);

	double dAlfaInact = 9.9977e-006;
	if (Vm != -34.0) dAlfaInact = 1.775e-6 * (Vm + 34.0) / (exp((Vm + 34.0) / 5.633) - 1.0);
	double dBetaInact = 0.427 * Cai * (Vm + 64.0) / (1.0 + exp(-(Vm + 44.0) / 4.16));       //beta

	double dCaAct_dt = dAlfaAct * (1 - CaAct) - dBetaAct * CaAct;
	double dCaInact_dt = dAlfaInact * (1 - CaInact) - dBetaInact * CaInact;

	//----------- Integrate ------------------------------------------------//
	CaAct = TwoState(CaAct, dCaAct_dt, dAlfaAct, dBetaAct, m_dDt);
	CaInact = TwoState(CaInact, dCaInact_dt, dAlfaInact, dBetaInact, m_dDt);

	//---------- Current Amplitude ------------------------------------------//	
	double PCaL = SF_CaL * 57.6;
	ICaL_Ca = PCaL * CaAct * CaInact * CF_Ca;
	ICaL_Na = PCaL * 5.0e-5 * CaAct * CaInact * CF_Na;
	ICaL_K = PCaL * 1.0e-3 * CaAct * CaInact * CF_K;
	ICaL = (ICaL_Ca + ICaL_Na + ICaL_K);
	PoCaL = CaAct * CaInact;
}
void ICC::VDDR_Current()
{
	double dY1_dt, dY2_dt;
	double alpha_Y1, beta_Y1, alpha_Y2, beta_Y2;
	double PVDDR;

	double DTinfinite = 1.0 / (1.0 + exp(-(Vm + 26.0) / 6.0));
	double tauDT = 0.6 + 5.4 / (1.0 + exp(0.03 * (Vm + 100.0)));
	alpha_Y1 = DTinfinite / tauDT;
	beta_Y1 = (1.0 - DTinfinite) / tauDT;

	double FTinfinite = 1.0 / (1 + exp((Vm - vddrH) / 6.0));
	double tauFT = 10.0 + 400.0 / (1.0 + exp(0.02 * (Vm + 65.0)));
	alpha_Y2 = FTinfinite / tauFT;
	beta_Y2 = (1.0 - FTinfinite) / tauFT;

	dY1_dt = alpha_Y1 * (1 - YVDDR1) - beta_Y1 * YVDDR1;
	dY2_dt = alpha_Y2 * (1 - YVDDR2) - beta_Y2 * YVDDR2;

	//----------- Integrate ------------------------------------------------//
	YVDDR1 = TwoState(YVDDR1, dY1_dt, alpha_Y1, beta_Y1, m_dDt);
	YVDDR2 = TwoState(YVDDR2, dY2_dt, alpha_Y2, beta_Y2, m_dDt);

	PVDDR = 0.5 * 15.84;
	IVDDR = SF_VDDR * PVDDR * CF_Ca * YVDDR1 * YVDDR2;
}
// Biophysical Journal   Volume 81   December 2001   3029–3051
// A Mathematical Model of Action Potential Heterogeneity in Adult Rat Left Ventricular Myocytes
// Sandeep V. Pandit, Robert B. Clark, Wayne R. Giles, and Semahat S. Demir
void ICC::PMCA()
{
	// original ICC model, electroneutral
	double PPMCA = 927.4;
	double KmCai = 0.008;
	IPMCA = SF_PMCA * PPMCA / (1 + (KmCai / Cai) * (KmCai / Cai));
}
void ICC::ClCa_Current()		// 
{
	double ECatZero = 1.39e-3;		// 1.39 uM	
	double kcclca = 0.44 * 0.01248;		// 0.01248 mV-1
	double Vhclca = -100.0;			// -100 mV
	double kvclca = 0.0156;			// 0.0156 mV-1
	double EC50 = ECatZero * exp(-kcclca * Vm);

	double Dcclca = 250.0e-3;		// um^2 s-1
	double Dmclca = 75.0e-3;		// um^2 s-1
	double Kmclca = 0.1 * 1.0e-3;			// 1 uM
	double Bmclca = 50.0e-3;		// 50 uM
	double rclca = 0.05;			// 0.05 um
	double m_phi = 3.141592653;

	double dCai = CaPMU;
	double C2clca = Dcclca * dCai - Dmclca * Bmclca * Kmclca / (Kmclca + dCai);

	double nsoc = 250.0;
	double sigma = JSOC * Vl_cyt / nsoc;

	double term1 = -Dcclca * Kmclca + sigma / 2.0 / m_phi / rclca + C2clca;
	double term2 = 4.0 * Dcclca * Dmclca * Bmclca * Kmclca;
	double term3 = Dcclca * Kmclca + sigma / 2.0 / m_phi / rclca + C2clca;
	CaClCa = (term1 + sqrt(term3 * term3 + term2)) / 2.0 / Dcclca;

	//---  Model I -------------------------------------------------------------------------------------//
	double OClCainf = 1.0 / (1.0 + exp((Vhclca - Vm) * kvclca)) / (1.0 + pow(EC50 / CaClCa, 2.0));
	double t1clca = 0.6 * 0.08163e3 * exp(-0.57 * CaClCa);			// s unit??
	double t2clca = 0.6 * 0.07617e3 * exp(-0.05374 * CaClCa);
	double t3clca = 1.9 * 70.3e3 * exp(0.153 * CaClCa);
	double tauClCa = t1clca + t2clca * exp(Vm / t3clca * 1000.0);
	OClCa = TwoState(OClCa, (OClCainf - OClCa) / tauClCa, OClCainf / tauClCa, (1.0 - OClCainf) / tauClCa, m_dDt);
	double gAclca = 20.0;			// 20 nS

	double ECl = mRTF * log(ClPMU / Cle);
	IClCa = SF_ClCa * gAclca * OClCa * (Vm - ECl);
}
void ICC::SOC_Current()
{
	double gSOC = 1.4 * 0.2 * 0.1;				// 0.1 nS, remremrem
	double KSOC = 200.0e-3;			// 200 uM, remremrem
	OSOC = pow(KSOC, 8.0) / (pow(KSOC, 8.0) + pow(CaSR * 0.021, 8.0));		// remremrem	
	double ECa = mRTF / 2.0 * log(Cae / CaPMU);
	ISOC = SF_SOC * gSOC * OSOC * (Vm - ECa);
	JSOC = -ISOC / mF / Vl_cyt;
}
void ICC::NaLCN_Current()
{
	// adjusted for repetitive firing, remrem
	PNaLCN = SF_NaLCN * 0.06 * 1.5 * 1.5;

	INaLCN_Na = PNaLCN * CF_Na * (1.1);
	INaLCN_K = 0.9 * PNaLCN * CF_K * (1.1);
	INaLCN = INaLCN_Na + INaLCN_K;
}
void ICC::Clb_Current()
{
	double PClb = 0.1 * SF_Clb;			// 
	IClb = PClb * CF_Cl;
}
void ICC::K1_Current()
{
	/*---------- Source ---------------------------------------------------------------------------------------------------------//
	Biophysical J 81(6) 3046. PMID: 11720973
	Sandeep V. Pandit, Robert B Clark, Wayne R Giles, and Semahat S. Demir. A mathematical model of action
	potential heterogeneity in adult rat left ventricular myocytes.
	--------------------------------------------------------------------------------------------------------------------------------*/
	//-------------- Current Amplitude ------------------------------------------------//
	double GK1;
	double EK = mRTF * log(Ke / Ki);			//REVERSAL POTENTIAL FOR K

	// Sandeep et al.
	GK1 = 1.0 * 0.3 * SF_K1 * 60;				// 0.024 microsiemen, modified on 20110807
	IK1 = GK1 * (Vm - EK - 1.73) / (1.0 + exp(1.613 * mFRT * (Vm - EK - 1.73))) * pow(Ke / 5.4, 0.4)
		+ (48.0 / (exp((Vm + 37.0) / 25.0) + exp(-(Vm + 37.0) / 25.0)) + 10.0) * (0.0001 / (1.0 + exp(-(Vm - EK - 76.77) / 17.0)));
}
void ICC::NaCa_Current()
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

	double mx[10][11] = {
	{ -kNa_b - kNae_off, kion_on * Nae, 0, 0, 0, 0, 0, 0, 0, kNa_f,0},
	{kNae_off, -kion_on * Nae - kNae_off, kion_on * Nae, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, kNae_off, -kNae_off - kion_on * Nae, kion_on * Nae, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, kNae_off, -kion_on * Nae - kion_on * Cae, kCae_off, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, kion_on * Cae, -kCae_off - kCa_f, kCa_b, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, kCa_f, -kCa_b - kCai_off, kion_on * dCai, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, kCai_off, -kion_on * dNai - kion_on * dCai, kNai_off, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, kion_on * dNai, -kion_on * dNai - kNai_off, kNai_off, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, kion_on * dNai, -kion_on * dNai - kNai_off, kNai_off, 0},
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Xtotal}
	};
	double tmp = 0;
	int i, j, dl, count;
	int N = 10;
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

	X2Na3 = mx[0][10];
	X2Na2 = mx[1][10];
	X2Na1 = mx[2][10];
	X2 = mx[3][10];
	X2Ca = mx[4][10];
	X1Ca = mx[5][10];
	X1 = mx[6][10];
	X1Na1 = mx[7][10];
	X1Na2 = mx[8][10];
	X1Na3 = mx[9][10];

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

	YNaCa = YNaCa;  											// Ca influx is negative of Na influx
	INaCa = PNaCa*(JX2_X1 - JX1_X2);
}
void ICC::NaK_Current()
{
	double PNaK = 0.7 * SF_NaK * 39.5;
	double KmNai = 3.5;
	double Km2Nai = 18.0;
	double KmKe = 0.27;
	INaK = 0.7 * PNaK / (1.0 + pow(KmNai / Nai, 1.36)) * (1.0 - ((Vm + 50.0) / 250.0) * ((Vm + 50.0) / 250.0)) / (1.0 + KmKe / Ke)
		+ 7.0 * PNaK / (1.0 + pow(Km2Nai / Nai, 3.0)) * (1.0 - ((Vm + 50.0) / 250.0) * ((Vm + 50.0) / 250.0)) / (1.0 + KmKe / Ke);
}
void ICC::NKCC()
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

	double mx[10][11] = {
	{ -kempty_b - kion_on * Nae, kNa_off, 0, 0, 0, 0, 0, 0, 0, kempty_f, 0},
	{kion_on * Nae, -kNa_off - kion_on * Cle, kCl_off, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, kion_on * Cle, -kCl_off - kion_on * Ke, kK_off, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, kion_on * Ke, -kK_off - kion_on * Cle, kCl_off, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, kion_on * Cle, -kCl_off - kfull_f, kfull_b, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, kfull_f, -kfull_b - kNa_off, kion_on * dNai, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, kNa_off, -kion_on * dNai - kCl_off, kion_on * dCli, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, kCl_off, -kion_on * dCli - kK_off, kion_on * dKi, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, kK_off, -kion_on * dKi - kCl_off, kion_on * dCli, 0},
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, Etotal}
	};
	double tmp = 0;
	int i, j, dl, count;
	int N = 10;
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
	
	E1 = mx[0][10];
	E1Na = mx[1][10];
	E1NaCl = mx[2][10];
	E1NaClK = mx[3][10];
	E1NaClKCl = mx[4][10];
	E2NaClKCl = mx[5][10];
	E2ClKCl = mx[6][10];
	E2KCl = mx[7][10];
	E2Cl = mx[8][10];
	E2 = mx[9][10];

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
	JNKCC_K = JNKCC_Na;
	JNKCC_Cl = 2.0*JNKCC_Na;
}
void ICC::IP3rates()
{
	kucIP3R = 400.0 * IP3;
	kcuIP3R = 400.0 * 0.9434e-3;			// d1: 0.13, d3: 0.9434 uM	
	kcoIP3R = 20.0 * Cai;
	kocIP3R = 20.0 * 0.08234e-3;			// d5: 0.08234 uM
	kobIP3R = 0.2 * Cai;
	kboIP3R = 0.2 * 1.049e-3;				// d2: 1.049, d4: 0.1445 uM		//*/
}
void ICC::IP3R_Current()
{
	IP3rates();
	double dCIP3R_dt, dOIP3R_dt, dUIP3R_dt, dBIP3R_dt;
	double PIP3R;

	//----------- (2) Colligative 4-state Model  ----------------------------------------//
	dCIP3R_dt = kucIP3R*aUIP3R + kocIP3R*aOIP3R - (kcuIP3R + kcoIP3R)*aCIP3R;
	dOIP3R_dt = kcoIP3R*aCIP3R + kboIP3R*aBIP3R - (kocIP3R + kobIP3R)*aOIP3R;
	dUIP3R_dt = kcuIP3R*aCIP3R - kucIP3R*aUIP3R;	
	dBIP3R_dt = kobIP3R*aOIP3R - kboIP3R*aBIP3R;		
						
	aCIP3R = Euler(aCIP3R, dCIP3R_dt, m_dDt);
	aBIP3R = Euler(aBIP3R, dBIP3R_dt, m_dDt);
	aUIP3R = Euler(aUIP3R, dUIP3R_dt, m_dDt);
	aOIP3R = 1.0 - aCIP3R - aBIP3R - aUIP3R;					

	PIP3R = 217.2;
	IIP3R = SF_IP3R*PIP3R*(CaiSR - Cai)*pow(aOIP3R, 3.0) * (lumExp > 0 ? pow((CaiSR > 0.0 ? CaiSR : 0.0) / lumRef, lumExp) : 1.0);		// SF_IP3R (default 1) added in the Coupling copy (used by gisim psd)				
}
void ICC::SERCA()
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

	KmCao_SERCA = 0.86 * 1.5 * 0.002; //1.5*0.002;				// mM	
	KmCai_SERCA = 0.08; //0.08;					// mM

	pE1CaSERCA = 1.0 / (1.0 + KmCai_SERCA / CaUP);
	pE2CaSERCA = 1.0 / (1.0 + KmCao_SERCA / Cai);
	pE1SERCA = 1.0 / (1.0 + CaUP / KmCai_SERCA);
	pE2SERCA = 1.0 / (1.0 + Cai / KmCao_SERCA);

	alpha_YSERCA = k2SERCA * pE2CaSERCA + k4SERCA * pE2SERCA;
	beta_YSERCA = k1SERCA * pE1CaSERCA + k3SERCA * pE1SERCA;

	dYSERCA_dt = alpha_YSERCA * (1 - YSERCA) - beta_YSERCA * YSERCA;

	//----------- Integrate ------------------------------------------------//
	YSERCA = TwoState(YSERCA, dYSERCA_dt, alpha_YSERCA, beta_YSERCA, m_dDt);

	double PSERCA = 0.1* 7728.0;			// remremrem
	ISERCA = SF_SERCA * 1.0 * PSERCA * (k2SERCA * pE2CaSERCA * (1 - YSERCA) - k1SERCA * pE1CaSERCA * YSERCA);
}
void ICC::InitialCellVariables()
{
	// from file
	SF_NaLCN = 1.000000e+000;
	PNaLCN = 1.350000e-001;
	INaLCN = -5.888044e+001;
	INaLCN_Na = -6.026994e+001;
	INaLCN_K = 1.389503e+000;
	SF_CaL = 1.000000e+000;
	CaAct = 8.980201e-008;
	CaInact = 6.804725e-002;
	PoCaL = 6.110780e-009;
	ICaL = -3.459157e-006;
	ICaL_Na = -7.142717e-009;
	ICaL_K = 3.659398e-009;
	ICaL_Ca = -3.455674e-006;
	SF_VDDR = 1.000000e+000;
	YVDDR1 = 4.428786e-004;
	YVDDR2 = 7.261682e-001;
	IVDDR = -5.001390e-002;
	SF_K1 = 1.000000e+000;
	nK1 = 8.632682e-001;
	IK1 = 6.981454e+001;
	SF_ClCa = 1.000000e+000;
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
	SF_Clb = 1.000000e+000;
	IClb = -2.268677e+001;
	SF_IP3R = 1.000000e+000;
	UIP3R = 8.840000e-001;
	CIP3R = 0.000000e+000;
	OIP3R = 0.000000e+000;
	BIP3R = 1.160000e-001;
	UIP3R0 = 6.083333e-001;
	CIP3R0 = 1.900000e-001;
	OIP3R0 = 8.300000e-002;
	BIP3R0 = 1.186667e-001;
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

	Vl_UP = 0.03 * 715.5;		// single SR compartment
	Vl_SR = 0.01 * 715.5;		// single SR compartment
}
void ICC::RunCellModel(BOOL AP, double InjectAmplitude)
{
	Ca_Buffer();
	Constant_Field();

	IP3R_Current();

	CaL_Current();
	VDDR_Current();
	K1_Current();
	NaCa_Current();
	NaK_Current();
	NaLCN_Current();
	ClCa_Current();
	SOC_Current();
	Clb_Current();

	SERCA();
	PMCA();
	NKCC();
	SR_Concentration();

	Ion_Concentration();
	IP3metabolism();

	Membrane_Potential(AP, InjectAmplitude);
	if (AP == FALSE) Vm = InjectAmplitude;	// should be restored remrem
}
void ICC::Membrane_Potential(BOOL AP, double Iinj)
{	
	Itot = ISOC + ICaL + IVDDR + IK1 + IClCa + IClb + INaK + INaCa + INaLCN + 0*IPMCA;		// total current	

	Iext = Iinj;

	if (AP) dVm_dt = -(Itot + Iext) / Cap;	// insert temporary code, injection current, 2007.09.18
	else dVm_dt = 1.0 / Cap * ((Iinj - Vm) / (0.001 * SF_AccessResistance) - Itot);

	Vm = Euler(Vm, dVm_dt, m_dDt);
}
void ICC::Constant_Field()
{	
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
		CF_Na = mFRT * Vm * (Nai - Nae * exp(-Vm * mFRT)) / (1 - exp(-Vm * mFRT));
		CF_K = mFRT * Vm * (Ki - Ke * exp(-Vm * mFRT)) / (1 - exp(-Vm * mFRT));
		CF_Ca = 2.0 * mFRT * Vm * (Cai - Cae * exp(-2.0 * Vm * mFRT)) / (1 - exp(-2.0 * Vm * mFRT));
		CF_Cl = mFRT * Vm * (Cli - Cle * exp(Vm * mFRT)) / (1 - exp(Vm * mFRT));
		CF_NaPMU = mFRT * Vm * (NaPMU - Nae * exp(-Vm * mFRT)) / (1 - exp(-Vm * mFRT));
		CF_KPMU = mFRT * Vm * (KPMU - Ke * exp(-Vm * mFRT)) / (1 - exp(-Vm * mFRT));
		CF_CaPMU = 2.0 * mFRT * Vm * (CaPMU - Cae * exp(-2.0 * Vm * mFRT)) / (1 - exp(-2.0 * Vm * mFRT));
		CF_ClPMU = mFRT * Vm * (ClPMU - Cle * exp(Vm * mFRT)) / (1 - exp(Vm * mFRT));
	}
}
void ICC::Ca_Buffer()
{	
	/*int N = 1;
	do
	{
		Cai = Cat / (1.0 + (TRPNf / Km_TRPN + CMDNf / Km_CMDN));
		CMDNf = CMDNt / (1.0 + (Cai / Km_CMDN));
		TRPNf = TRPNt / (1.0 + (Cai / Km_TRPN));
		N = N + 1;
	} while (N <= 10);
	N = 1;
	do
	{
		CaiSR = CaSR / (1.0 + (CSQNf / Km_CSQN));
		CSQNf = CSQNt / (1.0 + (CaiSR / Km_CSQN));
		N = N + 1;
	} while (N < 10);

	do
	{
		CaPMU = CaPMUt / (1.0 + (CMDNPMU / Km_CMDN));
		CMDNPMU = CMDNt / (1.0 + (CaPMU / Km_CMDN));
		N = N + 1;
	} while (N < 10);

	if (Cai >= 0 && CaiSR >= 0 && CaPMU >= 0) return true;
	else return false;*/

	double bmyo = CMDNt - Cat + Km_CMDN;
	double cmyo = Km_CMDN * Cat;
	Cai = (sqrt(bmyo * bmyo + 4 * cmyo) - bmyo) * 0.5;

	bmyo = CSQNt - CaSR + Km_CSQN;
	cmyo = Km_CSQN * CaSR;
	CaiSR = (sqrt(bmyo * bmyo + 4 * cmyo) - bmyo) * 0.5;

	bmyo = CMDNt - CaPMUt + Km_CMDN;
	cmyo = Km_CMDN * CaPMUt;
	CaPMU = (sqrt(bmyo * bmyo + 4 * cmyo) - bmyo) * 0.5;
}
void ICC::SR_Concentration()
{
	Itr = SF_SRT * 10.0 * 20.5 * (CaUP - CaiSR);
	Ileak = SF_SRL * 1.7 * 2.5 * 0.45 * (CaUP - Cai);
	double dCaSR_dt;
	dCaSR_dt = (Itr - IIP3R) / 2.0 / mF / Vl_SR;
	CaSR = Euler(CaSR, dCaSR_dt, m_dDt);

	double dCaUP_dt;
	dCaUP_dt = (ISERCA - Itr - Ileak) / 2.0 / mF / Vl_UP;
	CaUP = Euler(CaUP, dCaUP_dt, m_dDt);
	if (CaUP > 5.0)		CaUP = 5.0;
}
void ICC::IP3metabolism()
{
	double PIP2_IP3, IP3_PIP2, IP3_IP4, IP4_PIP2, IP4_IP3, PIP2_IP4;

	PIP2_IP3 = ip3ProdScale * 0.2 * exp(((ip3Vref + 48.5) + ip3VK * (Vm - ip3Vref)) / 18.1) * (Cai) / (Cai + 0.00001);		// kab	
	IP3_PIP2 = 0.5 * exp((Vm + 100.0) / (-28.5));						//kba	
	IP3_IP4 = 0.004;											//kbc			
	IP4_PIP2 = 0.0035 * exp((Vm + 100) / (-25.5));					//kca		
	PIP2_IP4 = 0.000;

	double dIP3_dt = PIP2 * PIP2_IP3 - IP3_IP4 * IP3 - IP3_PIP2 * IP3;
	double dPIP2_dt = IP4_PIP2 * (IP3total - IP3 - PIP2) + IP3_PIP2 * IP3 - (PIP2_IP4 + PIP2_IP3) * PIP2;

	//----------- Integrate ------------------------------------------------//
	twoelement ThreeStateSet;
	twoelement dThreeStateSet;
	sixelement RateConstSet;
	twoelement pThreeStateSet;

	double fIP3, fPIP2;		// fraction of IP3 and PIP2 to total IPs(IP4, IP3, PIP2)
	fIP3 = IP3 / IP3total;
	fPIP2 = PIP2 / IP3total;
	ThreeStateSet.a = fIP3;
	ThreeStateSet.b = fPIP2;
	dThreeStateSet.a = dIP3_dt / IP3total;
	dThreeStateSet.b = dPIP2_dt / IP3total;

	RateConstSet.a = IP3_IP4;
	RateConstSet.b = 0;
	RateConstSet.c = IP4_PIP2;
	RateConstSet.d = 0;
	RateConstSet.e = IP3_PIP2;
	RateConstSet.f = PIP2_IP3;

	pThreeStateSet = ThreeState(ThreeStateSet, dThreeStateSet, RateConstSet, m_dDt);

	fIP3 = pThreeStateSet.a;
	fPIP2 = pThreeStateSet.b;

	IP3 = fIP3 * IP3total;
	PIP2 = fPIP2 * IP3total;
}
void ICC::Ion_Concentration()
{	
	double dNa_dt, dK_dt, dCa_dt, dCl_dt;
	double dNaPMU_dt, dKPMU_dt, dCaPMU_dt, dClPMU_dt;
	double diffNa = 3.0e-1;
	double diffK = 3.0e-1;
	double diffCa = 2.0e-1;
	double diffCl = 2.0e-3;

	double Vl_PipettelVol = Vl_cyt * SF_PipetteVolumeRatio;
	dNa_dt = -(INaLCN_Na + ICaL_Na + 3.0 * INaK) / mF / Vl_PipettelVol + diffNa * (NaPMU - Nai);
	dK_dt = -(IK1 + ICaL_K + INaLCN_K - 2.0 * INaK) / mF / Vl_PipettelVol + diffK * (KPMU - Ki);
	dCa_dt = -(ISOC + ICaL_Ca + IVDDR + ISERCA + IPMCA - Ileak - IIP3R) / 2.0 / mF / Vl_PipettelVol + diffCa * (CaPMU - Cai);
	
	dCl_dt = IClb / mF / Vl_PipettelVol + diffCl * (ClPMU - Cli);

	dNaPMU_dt = -(3.0 * INaCa) / mF / Vl_PMU + JNKCC_Na + diffNa * (Nai - NaPMU) * Vl_PipettelVol / Vl_PMU;
	dKPMU_dt = JNKCC_K + diffK * (Ki - KPMU) * Vl_PipettelVol / Vl_PMU;
	dCaPMU_dt = -(-2.0 * INaCa) / 2.0 / mF / Vl_PMU + diffCa * (Cai - CaPMU) * Vl_PipettelVol / Vl_PMU;
	dClPMU_dt = IClCa / mF / Vl_PMU + JNKCC_Cl + diffCl * (Cli - ClPMU) * Vl_PipettelVol / Vl_PMU;

	Nai = Euler(Nai, dNa_dt, m_dDt);
	Ki = Euler(Ki, dK_dt, m_dDt);
	Cat = Euler(Cat, dCa_dt, m_dDt);
	Cli = Euler(Cli, dCl_dt, m_dDt);
	NaPMU = Euler(NaPMU, dNaPMU_dt, m_dDt);
	KPMU = Euler(KPMU, dKPMU_dt, m_dDt);
	CaPMUt = Euler(CaPMUt, dCaPMU_dt, m_dDt);
	ClPMU = Euler(ClPMU, dClPMU_dt, m_dDt);
}


double ICC::Euler(double dInput, double dSlope, double dTime)
{	
	return dInput+dSlope*dTime;	
}

double ICC::TwoState(double dInput, double dSlope, double Alpha, double Beta, double dTime)
{
	return dInput+dSlope/(Alpha+Beta)*(1-exp(-dTime*(Alpha+Beta)));
}

fourelement ICC::Get_State_Coeff(sixelement kVal)
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

twoelement ICC::ThreeState(twoelement dInput, twoelement dSlope, sixelement RateConst, double dTime)
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