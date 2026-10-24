#include <memory>
// Headless driver: pacemaker ICC-MY bidirectionally coupled to a downstream SIP compartment.
//
// Topology (2-compartment, step 1):
//   ICC-MY (1 cell, ICC model)  <-- g_c -->  SIP compartment = M isopotential mCSMC cells
//
// Gap junction current (outward-positive, pA, g in nS):
//   ICC side : Igap_I = g_c * (V_I - V_S)
//   SIP side : total current g_c*(V_S - V_I) is shared by M cells -> per-cell Igap_S = g_c*(V_S - V_I)/M
// Tonic PDGFRa-like SK conductance on the SIP compartment (per SMC, pA); keeps SMCs quiescent near -65 mV:
//   I_K0 = gK0 * (V_S + 80)      (global gK0, default 0.2 nS)
// Neural input on the SIP compartment (per SMC, pA):
//   I_n = g_n * (V_S - E_n)     E_n = -80 mV (inhibitory, K-like) or 0 mV (excitatory, NSCC-like)
//
// Usage:
//   gisim single_icc  T dt out.csv
//   gisim single_smc  T dt out.csv
//   gisim pair  gc M gn En T dt out.csv [gK0]    (time series of one run)
//   gisim sweep M T dt summary.csv [gK0]         (grid over gc x signed neural input, OpenMP parallel)
//        gnSigned < 0 : inhibitory (E=-80), > 0 : excitatory (E=0)

#include "pch.h"
#include "ICC.h"
#include "mCSMC.h"
#include "SmoothM.h"		// minimal SMC model (small intestine, per the model author)
#include "PohSMC.h"		// Poh et al. 2012 human jejunal SMC (generated from CellML; smcModel=2)
#include "pdgfra/PDGFRaCell.h"	// modified PDGFa (original model, SK = EK and voltage-independent)
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

struct Trace {
	std::vector<double> t, vI, vS, caI;
};

struct Summary {
	double freqI, ampI, maxI, minI;		// ICC-MY: cycles/min, peak-to-peak (mV), max, min
	double freqS, ampS, maxS, minS;		// SIP compartment
	double meanIgapI;					// mean load current on ICC (pA)
};

// Cycles/min from upward crossings of the midpoint with hysteresis; oscillation requires amp > 1 mV.
static void Analyze(const std::vector<double>& t, const std::vector<double>& v, double& freq, double& amp, double& vmax, double& vmin)
{
	vmax = *std::max_element(v.begin(), v.end());
	vmin = *std::min_element(v.begin(), v.end());
	amp = vmax - vmin;
	freq = 0;
	if (amp < 1.0) return;
	double hi = vmin + 0.6 * amp, lo = vmin + 0.4 * amp;
	bool armed = v[0] < lo;
	std::vector<double> up;
	for (size_t i = 1; i < v.size(); i++)
	{
		if (!armed && v[i] < lo) armed = true;
		if (armed && v[i] >= hi) { up.push_back(t[i]); armed = false; }
	}
	if (up.size() >= 2) freq = 60000.0 * (up.size() - 1) / (up.back() - up.front());
}

// Pre-equilibrated uncoupled cells, shared by all runs
static double gK0 = 0.2;
static ICC   gIccRest;
static double gIccSOC = 1.0, gIccANO1 = 1.0;	// keys iccSOC, iccANO1: scale the store-operated Ca2+ entry (SF_SOC) and ANO1 (SF_ClCa) of ICC-MY; 1 = Youm 2019
static double gIccRestSF[5] = {1, 1, 1, 1, 1};	// keys iccNaLCN, iccK1, iccClb, iccNaK, iccNKCC (iccRest sets all): scale the resting conductances and the steady transporters of ICC-MY together (low own conductance with the same balance); 1 = Youm 2019
static double gIccVddr = 1.0, gIccIp3VK = 1.0, gIccVH = -66.0;	// + iccVH: VDDR inactivation midpoint	// ICC-MY variants (keys iccVDDR, ip3VK): VDDR scale and voltage dependence of IP3 production; 1 = Youm 2019
static mCSMC gSmcRest;
static SmoothM gSmsRest;	// small-intestine SMC prototype (equilibrated like gSmcRest)
static PohSMC gPohRest;	// Poh et al. 2012 prototype, equilibrated uncoupled with no added current (its own rest)
// SmoothM conductance scales from environment variables SMS_CAL, SMS_KS, SMS_KCA, SMS_LEAK and the ICaL activation
// time-constant scale SMS_TAUCAL (all default 1) so that every mode
// (single cell, cable, sheet) uses the same spike calibration. Always initialise SmoothM through SmsInit().
static double SmsEnv(const char* k) { const char* v = getenv(k); return v ? atof(v) : 1.0; }
static void SmsInit(SmoothM& c)
{
	c.InitializeVar();
	c.sfical1 = SmsEnv("SMS_CAL"); c.sfiks1 = SmsEnv("SMS_KS"); c.sfikca1 = SmsEnv("SMS_KCA"); c.sfileak1 = SmsEnv("SMS_LEAK"); c.tauCaLScale = SmsEnv("SMS_TAUCAL");
}

static void Equilibrate(double Tpre, double dt)
{
	gIccRest.InitialCellVariables();
	gIccRest.m_dDt = dt;
	gIccRest.SF_VDDR = gIccVddr; gIccRest.ip3VK = gIccIp3VK; gIccRest.vddrH = gIccVH;
	gIccRest.SF_NaLCN = gIccRestSF[0]; gIccRest.SF_K1 = gIccRestSF[1]; gIccRest.SF_Clb = gIccRestSF[2]; gIccRest.SF_NaK = gIccRestSF[3]; gIccRest.SF_NKCC = gIccRestSF[4]; gIccRest.SF_SOC = gIccSOC; gIccRest.SF_ClCa = gIccANO1;
	gSmcRest.InitializeVar();
	gSmcRest.delT = dt;
	for (double t = 0; t < Tpre; t += dt)
	{
		gIccRest.RunCellModel(TRUE, 0.0);
		gSmcRest.Igap = gK0 * (gSmcRest.Vm + 80.0);
		gSmcRest.Integrate(0, t, false);
	}
	SmsInit(gSmsRest);
	gSmsRest.delT = dt;
	for (double t = 0; t < Tpre; t += dt) { gSmsRest.Igap = gK0 * (gSmsRest.Vm + 80.0); gSmsRest.Integrate(0); }
	gPohRest.InitializeVar(); gPohRest.delT = dt;
	for (double t = 0; t < Tpre; t += dt) { gPohRest.Igap = 0; gPohRest.Integrate(); }
}

// Runs the coupled pair; analysis uses the second half of T. Trace is filled only if tr != NULL.
static Summary RunPair(double gc, double M, double gn, double En, double T, double dt, Trace* tr, double sample = 1.0)
{
	ICC icc(gIccRest);
	mCSMC smc(gSmcRest);
	// explicit Euler on the gap term is stable only for gc*dt/Cap_ICC < 1; keep a 5x margin
	if (gc > 0) dt = std::min(dt, 0.2 * icc.Cap / gc);
	icc.m_dDt = dt;
	smc.delT = dt;

	std::vector<double> ta, vIa, vSa;
	double sumIgap = 0; long nIgap = 0;
	long step = 0, every = std::max(1L, (long)(sample / dt + 0.5));
	for (double t = 0; t < T; t += dt, step++)
	{
		double vI = icc.Vm, vS = smc.Vm;
		double IgapI = gc * (vI - vS);
		smc.Igap = gc * (vS - vI) / M + gK0 * (vS + 80.0) + gn * (vS - En);	// all extra currents use the external-current slot
		icc.RunCellModel(TRUE, IgapI);
		smc.Integrate(0, t, false);

		if (step % every == 0)
		{
			if (t >= 0.5 * T)
			{
				ta.push_back(t); vIa.push_back(icc.Vm); vSa.push_back(smc.Vm);
				sumIgap += IgapI; nIgap++;
			}
			if (tr) { tr->t.push_back(t); tr->vI.push_back(icc.Vm); tr->vS.push_back(smc.Vm); tr->caI.push_back(icc.Cai); }
		}
	}
	Summary s;
	Analyze(ta, vIa, s.freqI, s.ampI, s.maxI, s.minI);
	Analyze(ta, vSa, s.freqS, s.ampS, s.maxS, s.minS);
	s.meanIgapI = nIgap ? sumIgap / nIgap : 0;
	return s;
}

// ---------------------------------------------------------------------------------------------
// Step 2: ICC-MY node coupled to a 1D SIP cable through the muscle thickness.
//   Layer i (0..L-1) holds m isopotential SMCs per ICC-MY; layers are coupled by gSS per SMC pair.
//   Per-SMC external current in layer i:
//     gSS*(2V_i - V_{i-1} - V_{i+1})  (sealed end at L-1)  + [i==0] gc*(V_0 - V_I)/m  + gK0*(V_i+80) + gn*(V_i-E_n)
// ---------------------------------------------------------------------------------------------

// Slow-wave frequency with SMC spikes filtered out: 250 ms moving average (samples at 1 ms), then Analyze.
static double SlowWaveFreq(const std::vector<double>& t, const std::vector<double>& v)
{
	const int w = 250;
	if ((int)v.size() <= w) return 0;
	std::vector<double> ts, vs;
	double s = 0;
	for (int i = 0; i < (int)v.size(); i++)
	{
		s += v[i];
		if (i >= w) s -= v[i - w];
		if (i >= w - 1) { vs.push_back(s / w); ts.push_back(t[i]); }
	}
	double f, a, mx, mn;
	Analyze(ts, vs, f, a, mx, mn);
	return f;
}

struct CableSummary {
	double freqI, ampI, freqS0, meanIgapI;
	std::vector<double> ampL;
	double lambdaE;		// layers until amplitude falls to amp(0)/e (linear interpolation; L if never)
	int depth20;		// number of consecutive layers from 0 with amplitude >= 20 mV
};

static CableSummary RunCable(double gc, double m, double gSS, int L, double gn, double En, double T, double dt, const char* tracePath)
{
	ICC icc(gIccRest);
	std::vector<mCSMC> smc(L, gSmcRest);
	if (gc > 0) dt = std::min(dt, 0.2 * icc.Cap / gc);
	if (gSS > 0) dt = std::min(dt, 0.5 * smc[0].Cm / (2 * gSS));
	icc.m_dDt = dt;
	for (auto& c : smc) c.delT = dt;

	FILE* ft = tracePath ? fopen(tracePath, "w") : NULL;
	const int show[] = { 0, 1, 2, 5, 10, 20, L - 1 };
	if (ft) { fprintf(ft, "t_ms,V_ICC"); for (int k : show) fprintf(ft, ",V_L%d", k); fprintf(ft, "\n"); }

	std::vector<double> ta, vIa, V(L);
	std::vector<std::vector<double>> vL(L);
	double sumIgap = 0; long nIgap = 0;
	long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
	for (double t = 0; t < T; t += dt, step++)
	{
		double vI = icc.Vm;
		for (int i = 0; i < L; i++) V[i] = smc[i].Vm;
		double IgapI = gc * (vI - V[0]);
		for (int i = 0; i < L; i++)
		{
			double Ic = 0;
			if (i > 0)     Ic += gSS * (V[i] - V[i - 1]);
			if (i < L - 1) Ic += gSS * (V[i] - V[i + 1]);
			if (i == 0)    Ic += gc * (V[0] - vI) / m;
			smc[i].Igap = Ic + gK0 * (V[i] + 80.0) + gn * (V[i] - En);
		}
		icc.RunCellModel(TRUE, IgapI);
		for (int i = 0; i < L; i++) smc[i].Integrate(0, t, false);

		if (step % every == 0)
		{
			if (t >= 0.5 * T)
			{
				ta.push_back(t); vIa.push_back(icc.Vm);
				for (int i = 0; i < L; i++) vL[i].push_back(smc[i].Vm);
				sumIgap += IgapI; nIgap++;
			}
			if (ft) { fprintf(ft, "%.1f,%.3f", t, icc.Vm); for (int k : show) fprintf(ft, ",%.3f", smc[k].Vm); fprintf(ft, "\n"); }
		}
	}
	if (ft) fclose(ft);

	CableSummary s;
	double f, mx, mn;
	Analyze(ta, vIa, f, s.ampI, mx, mn);
	s.freqI = s.ampI >= 5.0 ? SlowWaveFreq(ta, vIa) : 0;
	s.freqS0 = SlowWaveFreq(ta, vL[0]);
	s.meanIgapI = nIgap ? sumIgap / nIgap : 0;
	s.ampL.resize(L);
	for (int i = 0; i < L; i++)
		s.ampL[i] = *std::max_element(vL[i].begin(), vL[i].end()) - *std::min_element(vL[i].begin(), vL[i].end());
	double target = s.ampL[0] / exp(1.0);
	s.lambdaE = L;
	for (int i = 1; i < L; i++)
		if (s.ampL[i] <= target) { s.lambdaE = (i - 1) + (s.ampL[i - 1] - target) / (s.ampL[i - 1] - s.ampL[i]); break; }
	s.depth20 = 0;
	while (s.depth20 < L && s.ampL[s.depth20] >= 20.0) s.depth20++;
	return s;
}

// ---------------------------------------------------------------------------------------------
// Step 3: active SIP. Each cable layer (per ICC-MY) holds three nodes:
//   m SMC (mCSMC + gK0), fP*m PDGFRa cells (mean-field SK), fI*m ICC-IM (ICC model with ip3ProdScale = sIM)
//   SMC-PDGFRa and SMC-ICC-IM gap conductance per cell pair: gSP, gSI (nS). Layer-to-layer coupling via SMC only.
// Neural inputs:
//   inhibitory (purinergic) : drive u in [0,1] raises PDGFRa Ca -> SK3 opens
//   excitatory (cholinergic): ICC-IM ip3ProdScale = sIM*(1+a)  (Gq -> IP3 -> Ca release -> ANO1)
// ---------------------------------------------------------------------------------------------

// Mean-field PDGFRa+ cell. SK3 is voltage independent, Hill-type Ca activation.
struct PDGFRa {
	static constexpr double Cm = 25.0;			// pF
	static constexpr double gSKmax = 58.0;		// nS (20000 channels x 2.9 pS, as in the PDGFa model)
	static constexpr double Kd = 4.0e-4;		// mM (0.4 uM)
	static constexpr double nH = 4.0;
	static constexpr double tauSK = 10.0;		// ms
	static constexpr double tauCa = 150.0;		// ms
	static constexpr double CaRest = 1.0e-4;	// mM
	static constexpr double dCaMax = 1.5e-3;	// mM, Ca rise at full purinergic drive
	static constexpr double gbg = 0.3, Ebg = -20.0, EK = -80.0;
	double V = -65.0, Ca = CaRest, Po = 0.0;
	double ISK = 0;
	void Step(double dt, double Iext, double u)
	{
		Ca += (CaRest + u * dCaMax - Ca) * dt / tauCa;
		double cn = pow(Ca, nH);
		Po += (cn / (cn + pow(Kd, nH)) - Po) * dt / tauSK;
		ISK = gSKmax * Po * (V - EK);
		V += -(ISK + gbg * (V - Ebg) + Iext) / Cm * dt;
	}
};

// Hypothetical fast-recovering ICC-IM (phenomenological): voltage-triggered Ca entry -> ANO1, lumped into one gate.
//   I = gAno * a * r * (V - ECl) + gL * (V - EL)
//   a: fast activation, a_inf = 1/(1+exp(-(V-Vth)/3)), tau 20 ms
//   r: inactivation,   r_inf = 1/(1+exp((V+62)/3)),    tau = tauR (recovery time constant); little a*r window current at rest
// Silent at rest (Vth above rest); fires when depolarized past Vth; refractory ~ tauR.
struct FastIM {
	static constexpr double Cm = 25.0, gL = 1.0, EL = -70.0, ECl = -20.0, tauA = 20.0;
	double V = -68.0, a = 0.0, r = 1.0;
	void Init(double Vth) { a = 1.0 / (1.0 + exp(-(V - Vth) / 3.0)); r = 1.0 / (1.0 + exp((V + 62.0) / 3.0)); }
	void Step(double dt, double Iext, double Vth, double tauR, double gAno)
	{
		a += (1.0 / (1.0 + exp(-(V - Vth) / 3.0)) - a) * dt / tauA;
		r += (1.0 / (1.0 + exp((V + 62.0) / 3.0)) - r) * dt / tauR;
		V += -(gAno * a * r * (V - ECl) + gL * (V - EL) + Iext) / Cm * dt;
	}
};

// Small-intestine ICC-DMP (imModel 2): no voltage-dependent Ca2+ channels. Each cell fires stochastic Ca2+ transients
// (Poisson, rate independent of Vm) that open ANO1 and produce STICs; neural input only changes the event rate.
//   Cm dV/dt = -(gA(t) (V - ECl) + gL (V - EL) + Iext),  gA = gU * sum over events of a difference of exponentials (peak gU)
// Parameters are placeholders until calibrated to ICC-DMP Ca2+ imaging / STIC data.
struct StochDMP {
	double V = -70.0, s1 = 0.0, s2 = 0.0;	// s1 decays with tauD, s2 with tauR; gA = gU * norm * (s1 - s2)
	std::mt19937 rng;
	void Step(double dt, double Iext, double rate_per_ms, double Cm, double gL, double EL, double gU, double tauR, double tauD, double norm)
	{
		if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < rate_per_ms * dt) { s1 += 1.0; s2 += 1.0; }
		s1 -= s1 * dt / tauD; s2 -= s2 * dt / tauR;
		double gA = gU * norm * (s1 - s2);
		V += -(gA * (V + 20.0) + gL * (V - EL) + Iext) / Cm * dt;
	}
};

struct ActiveParams {
	double m = 3, gSS = 30, fP = 0.1, fI = 0.1, gSP = 3, gSI = 3, sIM = 0.03, vIM = 10;	// vIM: ICC-IM SF_VDDR (excitable, non-pacemaking)
	int L = 30;
	int imModel = 0;						// 0: ICC model as ICC-IM, 1: FastIM, 2: StochDMP (ICC-DMP only in the deepest dmpN layers; fI = DMP cells per SMC there)
	int dmpN = 1, dmpK = 10;				// StochDMP: number of deep layers holding ICC-DMP; independent cells simulated per layer
	double dmpRate = 0.5, dmpG = 1.0, dmpTr = 50, dmpTd = 300, dmpCm = 20, dmpGL = 1.0, dmpEL = -70;	// Hz per cell, nS peak per event, ms, pF, nS, mV; rate x (1 + aI)
	double vth = -52, tauR = 1000, gAno = 10;	// FastIM parameters; cholinergic gain a lowers Vth by 3*a mV
	int pdModel = 0;						// 0: mean-field PDGFRa, 1: modified PDGFa v2b (one per layer, own seed, dt 0.02 ms sub-steps)
	double ip3AgMax = 2e-3;					// mM agonist IP3 at purinergic drive uP = 1 (pdModel 1)
	double skKd = 2.0e-3, skHill = 4;		// SK3 sensitivity to submembrane Ca2+ (pdModel 1); Kd 2 uM fitted to Kito et al. 2014 FLC data
	double pdNaLCN = 2.5;					// NaLCN scale in the PDGFa cell (original 2.0); 2.5 with Kd 2 uM fits FLC in situ (Kito et al. 2014)
	double uOn = -1, uOff = -1;				// ms; if uOn < uOff the purinergic drive uP is applied only in [uOn, uOff) (pdModel 1)
	double wavePer = 0, waveHW = 670, waveRise = 60, waveFall = 400, waveRest = -70, waveAmp = 50;	// wavePer > 0 (ms): ICC-MY Vm prescribed as a trapezoid (half-width waveHW, linear rise/fall ms, first upstroke at 1 s); ICC model not stepped
	double iccGL = 0, iccEL = -72;			// extra passive conductance on ICC-MY (nS, real units for the scaled cell; may be negative) and its reversal (mV); sensitivity of the divider to the ICC-MY own conductance
	double iccInj = 0, injOn = -1, injOff = -1;	// pA (outward positive) injected into ICC-MY in [injOn, injOff) ms, for input-resistance tests
	double smcInj = 0;						// pA (outward positive) injected into SMC layer 1 in the same window (syncytial input conductance)
	int injN = 1;							// smcInj goes into each of SMC layers 0..injN-1 (propagation tests)
	int dumpPeak = 0;						// 1: print peak Vm of every SMC layer over the second half of the run
	int traceAll = 0;						// 1: also write <trace>.layers.csv with V_SMC of every layer (1-ms samples)
	double iccScale = 1.0;					// ICC-MY size: all currents and Cap scaled by s (implemented by dividing external current by s); 0.72 -> Cm 18 pF
	double gK0run = -1;						// SMC tonic K conductance during the run (-1: use global gK0, which is also used for equilibration)
	int seed = 1;
	int smcModel = 0;						// 0: mouse colonic SMC (mCSMC), 1: small-intestine SMC (SmoothM), 2: Poh et al. 2012 human jejunal SMC
	// extra dihydropyridine-resistant, voltage-dependent inward current in every SMC (T-type-like; test of active amplification of the
	// slow wave): I = gT m h (V - ET), dm/dt = (minf - m)/tauM, dh/dt = (hinf - h)/tauH, minf = 1/(1+exp(-(V-vmT)/kmT)), hinf = 1/(1+exp((V-vhT)/khT))
	double gT = 0, ET = 40, vmT = -50, kmT = 5, vhT = -70, khT = 6, tauM = 5, tauH = 100;
	double gKATP = 0, EKATP = -90;
	// longitudinal-muscle (LM) column on the other side of ICC-MY: LL layers of mL SMC per ICC-MY, coupled to ICC-MY by gcL and to each other
	// by gSSL; same SMC model, tonic K and K_ATP as the circular column, no FLC or T-type current (default LL 0 = no LM)
	double gIF = 0, fP0 = -1;	// direct ICC-MY to FLC coupling: gIF (nS, per ICC-MY) to the FLC of layer 1, whose density fP0 (FLC per SMC; -1: fP) may differ from deeper layers (pdModel 1)
	int LL = 0; double gcL = 0, mL = 0.3, gSSL = 30;			// extra K+ conductance in every SMC (nS) and its reversal (mV): K_ATP opener (pinacidil) test
};

static ICC gIccImRest;
static PDGFRa gPdgfRest;
static PDGFRaCell gPdgfaRest;

static void EquilibrateActive(const ActiveParams& P, double Tpre, double dt)
{
	gIccImRest.InitialCellVariables();
	gIccImRest.m_dDt = dt;
	gIccImRest.ip3ProdScale = P.sIM;
	gIccImRest.SF_VDDR = P.vIM;
	for (double t = 0; t < Tpre; t += dt) gIccImRest.RunCellModel(TRUE, 0.0);
	// PDGFRa alone relaxes to its own rest; it is re-equilibrated inside the syncytium during the run
	for (double t = 0; t < Tpre; t += dt) gPdgfRest.Step(dt, 0.0, 0.0);
	if (P.pdModel == 1)
	{
		// modified PDGFa equilibrated 10 s while held near the syncytium rest (-65 mV through 50 nS)
		gPdgfaRest = PDGFRaCell(1);
		gPdgfaRest.Init();
		gPdgfaRest.SetModified(true, true, P.pdNaLCN, 1.0);	// SK = EK, voltage-independent, deterministic gate; gap current carried by K+
		gPdgfaRest.SetSKCaSource(1, P.skKd, P.skHill);	// SK senses submembrane (IP3R release) Ca2+
		for (double t = 0; t < 10000; t += 0.02) gPdgfaRest.Step(0.02, 50.0 * (gPdgfaRest.Vm() + 65.0));
	}
}

static inline void StepSmc(mCSMC& c, double t) { c.Integrate(0, t, false); }
static inline void StepSmc(SmoothM& c, double) { c.Integrate(0); }
static inline void StepSmc(PohSMC& c, double) { c.Integrate(); }

// Prescribed ICC-MY waveform (wavePer > 0): trapezoid from waveRest to waveRest + waveAmp, period wavePer, first upstroke at 1000 ms;
// the flat top lasts waveHW - (waveRise + waveFall) / 2 so that the half-amplitude width equals waveHW.
static double IccWave(const ActiveParams& P, double t)
{
	double x = fmod(t - 1000.0 + 100.0 * P.wavePer, P.wavePer);
	double flat = std::max(0.0, P.waveHW - 0.5 * (P.waveRise + P.waveFall));
	double f = x < P.waveRise ? x / P.waveRise : x < P.waveRise + flat ? 1.0 : x < P.waveRise + flat + P.waveFall ? 1.0 - (x - P.waveRise - flat) / P.waveFall : 0.0;
	return P.waveRest + P.waveAmp * f;
}

template <class SMC>
static CableSummary RunActiveT(const ActiveParams& P, const SMC& smcRest, double gc, double uP, double aI, double T, double dt, const char* tracePath)
{
	const int L = P.L;
	ICC icc(gIccRest);
	std::vector<SMC> smc(L, smcRest);
	std::vector<ICC> im(P.imModel == 0 ? L : 0, gIccImRest);
	std::vector<FastIM> fim(P.imModel == 1 ? L : 0);
	// StochDMP: dmpK independent ICC-DMP cells in each of the deepest dmpN layers (dmp[i] empty elsewhere)
	std::vector<std::vector<StochDMP>> dmp(L);
	const double dmpNorm = [&] { double tp = P.dmpTr * P.dmpTd / (P.dmpTd - P.dmpTr) * log(P.dmpTd / P.dmpTr); return 1.0 / (exp(-tp / P.dmpTd) - exp(-tp / P.dmpTr)); }();
	const double dmpRate = P.dmpRate * std::max(0.0, 1.0 + aI) / 1000.0;	// per ms
	if (P.imModel == 2)
		for (int i = std::max(0, L - P.dmpN); i < L; i++)
		{
			dmp[i].resize(P.dmpK);
			for (int k = 0; k < P.dmpK; k++) { dmp[i][k].V = P.dmpEL; dmp[i][k].rng.seed(7919u * P.seed + 101u * i + k); }
		}
	const double vthEff = P.vth - 3.0 * aI;
	for (auto& c : fim) c.Init(vthEff);
	std::vector<PDGFRa> pd(L, gPdgfRest);
	std::vector<PDGFRaCell> pdf;
	if (P.pdModel == 1) for (int i = 0; i < L; i++) { pdf.push_back(gPdgfaRest); pdf.back().Reseed(1000u * P.seed + i); pdf.back().SetIP3Agonist(P.uOn < P.uOff ? 0.0 : uP * P.ip3AgMax); }
	const double gK0eff = P.gK0run >= 0 ? P.gK0run : gK0;
	if (gc + P.gIF > 0) dt = std::min(dt, 0.2 * icc.Cap * P.iccScale / (gc + P.gIF));
	if (P.gSS > 0) dt = std::min(dt, 0.5 * smc[0].Cm / (2 * P.gSS));
	icc.m_dDt = dt;
	for (auto& c : smc) c.delT = dt;
	for (auto& c : im) { c.m_dDt = dt; c.ip3ProdScale = P.sIM * (1.0 + aI); }

	const int LL0 = P.gcL > 0 ? P.LL : 0;
	FILE* ft = tracePath ? fopen(tracePath, "w") : NULL;
	const int show[] = { 0, std::min(5, L - 1), std::min(10, L - 1), std::min(20, L - 1) };	// clamped for short cables (column names give the actual layer)
	if (ft) { fprintf(ft, "t_ms,V_ICC"); for (int k : show) fprintf(ft, ",V_SMC%d,V_IM%d,V_P%d", k, k, k); if (P.imModel == 2) fprintf(ft, ",V_SMClast,V_DMPmean,V_DMP0"); if (LL0 > 0) fprintf(ft, ",V_LM0"); fprintf(ft, "\n"); }
	FILE* fa = (tracePath && P.traceAll) ? fopen((std::string(tracePath) + ".layers.csv").c_str(), "w") : NULL;
	if (fa) { fprintf(fa, "t_ms,V_ICC"); for (int i = 0; i < L; i++) fprintf(fa, ",V_SMC%d", i); fprintf(fa, "\n"); }

	std::vector<double> ta, vIa, V(L), VI(L), VP(L);
	const int LL = P.gcL > 0 ? P.LL : 0;
	std::vector<SMC> lm(LL, smcRest);
	std::vector<double> VLm(LL);
	for (auto& c : lm) c.delT = dt;
	auto mInfT = [&](double v) { return 1.0 / (1.0 + exp(-(v - P.vmT) / P.kmT)); };
	auto hInfT = [&](double v) { return 1.0 / (1.0 + exp((v - P.vhT) / P.khT)); };
	std::vector<double> mT(L), hT(L);
	for (int i = 0; i < L; i++) { mT[i] = mInfT(smc[i].Vm); hT[i] = hInfT(smc[i].Vm); }
	std::vector<std::vector<double>> vL(L);
	double sumIgap = 0; long nIgap = 0;
	long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
	for (double t = 0; t < T; t += dt, step++)
	{
		double vI = icc.Vm;
		for (int i = 0; i < L; i++) { V[i] = smc[i].Vm; VI[i] = P.imModel == 1 ? fim[i].V : P.imModel == 0 ? im[i].Vm : 0.0;
			if (P.imModel == 2 && !dmp[i].empty()) { VI[i] = 0; for (auto& c : dmp[i]) VI[i] += c.V / dmp[i].size(); } VP[i] = P.pdModel ? pdf[i].Vm() : pd[i].V; }
		for (int i = 0; i < LL; i++) VLm[i] = lm[i].Vm;
		double IgapI = gc * (vI - V[0]) + (LL > 0 ? P.gcL * (vI - VLm[0]) : 0.0) + P.gIF * (vI - VP[0]);
		for (int i = 0; i < LL; i++)
		{
			double Ic = 0;
			if (i > 0)      Ic += P.gSSL * (VLm[i] - VLm[i - 1]);
			if (i < LL - 1) Ic += P.gSSL * (VLm[i] - VLm[i + 1]);
			if (i == 0)     Ic += P.gcL * (VLm[0] - vI) / P.mL;
			lm[i].Igap = Ic + gK0eff * (VLm[i] + 80.0) + P.gKATP * (VLm[i] - P.EKATP);
		}
		for (int i = 0; i < L; i++)
		{
			double Ic = 0;
			if (i > 0)     Ic += P.gSS * (V[i] - V[i - 1]);
			if (i < L - 1) Ic += P.gSS * (V[i] - V[i + 1]);
			if (i == 0)    Ic += gc * (V[0] - vI) / P.m;
			Ic += (i == 0 && P.fP0 >= 0 ? P.fP0 : P.fP) * P.gSP * (V[i] - VP[i]);
			if (P.imModel != 2 || !dmp[i].empty()) Ic += P.fI * P.gSI * (V[i] - VI[i]);	// mean DMP Vm: exact, since the gap current is linear
			smc[i].Igap = Ic + gK0eff * (V[i] + 80.0) + P.gKATP * (V[i] - P.EKATP);
			if (P.gT > 0)
			{
				smc[i].Igap += P.gT * mT[i] * hT[i] * (V[i] - P.ET);
				mT[i] += dt * (mInfT(V[i]) - mT[i]) / P.tauM;
				hT[i] += dt * (hInfT(V[i]) - hT[i]) / P.tauH;
			}
		}
		if (P.smcInj != 0 && t >= P.injOn && t < P.injOff) for (int i = 0; i < std::min(P.injN, L); i++) smc[i].Igap += P.smcInj / P.m;
		if (P.wavePer > 0) icc.Vm = IccWave(P, t + dt);
		else icc.RunCellModel(TRUE, (IgapI + P.iccGL * (vI - P.iccEL) + ((t >= P.injOn && t < P.injOff) ? P.iccInj : 0.0)) / P.iccScale);
		if (P.pdModel == 1 && P.fP > 0)
		{
			if (P.uOn < P.uOff) { double a = (t >= P.uOn && t < P.uOff) ? uP * P.ip3AgMax : 0.0; for (auto& c : pdf) c.SetIP3Agonist(a); }
			const int nsub = std::max(1, (int)(dt / 0.02 + 0.5));
			const double dts = dt / nsub;
#pragma omp parallel for schedule(static)
			for (int i = 0; i < L; i++)
				for (int k = 0; k < nsub; k++) pdf[i].Step(dts, P.gSP * (pdf[i].Vm() - V[i])
					+ (i == 0 && P.gIF > 0 ? P.gIF * (pdf[i].Vm() - vI) / (P.m * (P.fP0 >= 0 ? P.fP0 : P.fP)) : 0.0));
		}
		for (int i = 0; i < LL; i++) StepSmc(lm[i], t);
		for (int i = 0; i < L; i++)
		{
			StepSmc(smc[i], t);
			if (P.fI > 0)
			{
				if (P.imModel == 1) fim[i].Step(dt, P.gSI * (VI[i] - V[i]), vthEff, P.tauR, P.gAno);
				else if (P.imModel == 0) im[i].RunCellModel(TRUE, P.gSI * (VI[i] - V[i]));
				else for (auto& c : dmp[i]) c.Step(dt, P.gSI * (c.V - V[i]), dmpRate, P.dmpCm, P.dmpGL, P.dmpEL, P.dmpG, P.dmpTr, P.dmpTd, dmpNorm);
			}
			if (P.fP > 0 && P.pdModel == 0) pd[i].Step(dt, P.gSP * (VP[i] - V[i]), (P.uOn < P.uOff) ? ((t >= P.uOn && t < P.uOff) ? uP : 0.0) : uP);	// uOn/uOff pulse also for the mean-field FLC
		}

		if (step % every == 0)
		{
			if (t >= 0.5 * T)
			{
				ta.push_back(t); vIa.push_back(icc.Vm);
				for (int i = 0; i < L; i++) vL[i].push_back(smc[i].Vm);
				sumIgap += IgapI; nIgap++;
			}
			if (ft)
			{
				fprintf(ft, "%.1f,%.3f", t, icc.Vm);
				for (int k : show) fprintf(ft, ",%.3f,%.3f,%.3f", smc[k].Vm, P.imModel == 1 ? fim[k].V : (P.imModel == 0 ? (P.fI > 0 ? im[k].Vm : 0.0) : VI[k]), P.pdModel ? pdf[k].Vm() : pd[k].V);
				if (P.imModel == 2) fprintf(ft, ",%.3f,%.3f,%.3f", smc[L - 1].Vm, VI[L - 1], dmp[L - 1].empty() ? 0.0 : dmp[L - 1][0].V);
				if (LL > 0) fprintf(ft, ",%.3f", lm[0].Vm);
				fprintf(ft, "\n");
			}
			if (fa) { fprintf(fa, "%.1f,%.3f", t, icc.Vm); for (int i = 0; i < L; i++) fprintf(fa, ",%.3f", smc[i].Vm); fprintf(fa, "\n"); }
		}
	}
	if (ft) fclose(ft);
	if (fa) fclose(fa);
	if (P.dumpPeak)
	{
		printf("peak by layer:");
		for (int i = 0; i < L; i++) printf(" %.1f", *std::max_element(vL[i].begin(), vL[i].end()));
		printf("\n");
	}

	CableSummary s;
	double f, mx, mn;
	Analyze(ta, vIa, f, s.ampI, mx, mn);
	s.freqI = s.ampI >= 5.0 ? SlowWaveFreq(ta, vIa) : 0;
	s.freqS0 = SlowWaveFreq(ta, vL[0]);
	s.meanIgapI = nIgap ? sumIgap / nIgap : 0;
	s.ampL.resize(L);
	for (int i = 0; i < L; i++)
		s.ampL[i] = *std::max_element(vL[i].begin(), vL[i].end()) - *std::min_element(vL[i].begin(), vL[i].end());
	double target = s.ampL[0] / exp(1.0);
	s.lambdaE = L;
	for (int i = 1; i < L; i++)
		if (s.ampL[i] <= target) { s.lambdaE = (i - 1) + (s.ampL[i - 1] - target) / (s.ampL[i - 1] - s.ampL[i]); break; }
	s.depth20 = 0;
	while (s.depth20 < L && s.ampL[s.depth20] >= 20.0) s.depth20++;
	return s;
}

static CableSummary RunActive(const ActiveParams& P, double gc, double uP, double aI, double T, double dt, const char* tracePath)
{
	if (P.smcModel == 2) return RunActiveT(P, gPohRest, gc, uP, aI, T, dt, tracePath);
	return P.smcModel == 1 ? RunActiveT(P, gSmsRest, gc, uP, aI, T, dt, tracePath) : RunActiveT(P, gSmcRest, gc, uP, aI, T, dt, tracePath);
}

static void PrintCable(const CableSummary& s)
{
	printf("ICC-MY %.2f cpm amp %.1f | L0 %.2f cpm | amp by layer:", s.freqI, s.ampI, s.freqS0);
	for (size_t i = 0; i < s.ampL.size(); i++) if (i < 3 || i % 5 == 0 || i == s.ampL.size() - 1) printf(" [%zu]%.1f", i, s.ampL[i]);
	printf(" | lambdaE %.1f depth20 %d | Igap %.1f pA\n", s.lambdaE, s.depth20, s.meanIgapI);
}

static void WriteTrace(const char* path, const Trace& tr)
{
	FILE* f = fopen(path, "w");
	if (!f) { fprintf(stderr, "cannot open %s\n", path); return; }
	fprintf(f, "t_ms,V_ICC,V_SIP,Cai_ICC\n");
	for (size_t i = 0; i < tr.t.size(); i++) fprintf(f, "%.2f,%.4f,%.4f,%.6e\n", tr.t[i], tr.vI[i], tr.vS[i], tr.caI[i]);
	fclose(f);
}

static void PrintSummary(const Summary& s)
{
	printf("ICC-MY: %.2f cpm, amp %.2f mV (%.2f..%.2f) | SIP: %.2f cpm, amp %.2f mV (%.2f..%.2f) | mean Igap on ICC %.2f pA\n",
		s.freqI, s.ampI, s.minI, s.maxI, s.freqS, s.ampS, s.minS, s.maxS, s.meanIgapI);
}

int main(int argc, char** argv)
{
	if (argc < 2) { fprintf(stderr, "see header of main.cpp for usage\n"); return 1; }
	std::string mode = argv[1];
	const double Tpre = 20000;	// ms, uncoupled equilibration

	if (mode == "single_icc" || mode == "single_smc")
	{
		double T = atof(argv[2]), dt = atof(argv[3]);
		ICC icc; icc.InitialCellVariables(); icc.m_dDt = dt;
		mCSMC smc; smc.InitializeVar(); smc.delT = dt;
		Trace tr; long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			if (mode == "single_icc") icc.RunCellModel(TRUE, 0.0); else smc.Integrate(0, t, false);
			if (step % every == 0)
			{
				tr.t.push_back(t);
				tr.vI.push_back(icc.Vm); tr.vS.push_back(smc.Vm); tr.caI.push_back(icc.Cai);
			}
		}
		WriteTrace(argv[4], tr);
		double f, a, mx, mn;
		std::vector<double> th(tr.t.begin() + tr.t.size() / 2, tr.t.end());
		std::vector<double> vh = (mode == "single_icc") ? std::vector<double>(tr.vI.begin() + tr.vI.size() / 2, tr.vI.end())
		                                                : std::vector<double>(tr.vS.begin() + tr.vS.size() / 2, tr.vS.end());
		Analyze(th, vh, f, a, mx, mn);
		printf("%s: %.2f cpm, amp %.2f mV, range %.2f..%.2f mV\n", mode.c_str(), f, a, mn, mx);
		return 0;
	}

	if (mode == "pair")
	{
		double gc = atof(argv[2]), M = atof(argv[3]), gn = atof(argv[4]), En = atof(argv[5]);
		double T = atof(argv[6]), dt = atof(argv[7]);
		if (argc > 9) gK0 = atof(argv[9]);
		Equilibrate(Tpre, dt);
		Trace tr;
		Summary s = RunPair(gc, M, gn, En, T, dt, &tr);
		WriteTrace(argv[8], tr);
		PrintSummary(s);
		return 0;
	}

	if (mode == "sweep")
	{
		double M = atof(argv[2]), T = atof(argv[3]), dt = atof(argv[4]);
		const char* out = argv[5];
		if (argc > 6) gK0 = atof(argv[6]);
		std::vector<double> gcs = { 0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500 };
		std::vector<double> gns = { -2, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 2 };	// signed, nS per SMC
		Equilibrate(Tpre, dt);
		int n = (int)(gcs.size() * gns.size());
		std::vector<Summary> res(n);
#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < n; k++)
		{
			double g = gns[k % gns.size()];
			res[k] = RunPair(gcs[k / gns.size()], M, fabs(g), g < 0 ? -80.0 : 0.0, T, dt, NULL);
		}
		FILE* f = fopen(out, "w");
		fprintf(f, "gc_nS,M,gK0,gnSigned,freqICC_cpm,ampICC_mV,maxICC,minICC,freqSIP_cpm,ampSIP_mV,maxSIP,minSIP,meanIgapICC_pA\n");
		for (int k = 0; k < n; k++)
		{
			const Summary& s = res[k];
			fprintf(f, "%g,%g,%g,%g,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", gcs[k / gns.size()], M, gK0, gns[k % gns.size()],
				s.freqI, s.ampI, s.maxI, s.minI, s.freqS, s.ampS, s.maxS, s.minS, s.meanIgapI);
		}
		fclose(f);
		printf("wrote %s (%d runs)\n", out, n);
		return 0;
	}
	if (mode == "sms")
	{
		// gisim sms gc gK T dt out.csv : SmoothM (small-intestine SMC) alone (gc = 0) or coupled to one ICC-MY through gc (nS);
		// gK adds a tonic K+ conductance toward -80 mV (nS). SmoothM is equilibrated 20 s alone first.
		double gc = atof(argv[2]), gK = atof(argv[3]), T = atof(argv[4]), dt = atof(argv[5]);
		SmoothM c; SmsInit(c); c.delT = dt;
		for (double t = 0; t < 20000; t += dt) { c.Igap = gK * (c.Vm + 80.0); c.Integrate(0); }
		Equilibrate(Tpre, dt);
		ICC icc(gIccRest); icc.m_dDt = dt;
		if (gc > 0) { dt = std::min(dt, 0.2 * icc.Cap / gc); icc.m_dDt = dt; c.delT = dt; }
		Trace tr; long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			double vI = icc.Vm, vS = c.Vm;
			c.Igap = gc * (vS - vI) + gK * (vS + 80.0);
			if (gc > 0) icc.RunCellModel(TRUE, gc * (vI - vS));
			c.Integrate(0);
			if (step % every == 0) { tr.t.push_back(t); tr.vI.push_back(gc > 0 ? icc.Vm : 0.0); tr.vS.push_back(c.Vm); tr.caI.push_back(c.Cai); }
		}
		WriteTrace(argv[6], tr);
		std::vector<double> th(tr.t.begin() + tr.t.size() / 2, tr.t.end()), vh(tr.vS.begin() + tr.vS.size() / 2, tr.vS.end());
		double f, a, mx, mn;
		Analyze(th, vh, f, a, mx, mn);
		printf("SMS gc %.1f gK %.2f: %.2f cpm, amp %.1f mV, range %.1f..%.1f mV\n", gc, gK, f, a, mn, mx);
		return 0;
	}

	if (mode == "sheet")
	{
		// gisim sheet nx ny gSS gK side amp T dt out_prefix : SMC-only sheet (SmoothM, nearest-neighbour gap gSS nS, tonic K+ gK nS
		// toward -80 mV, sealed edges). ny = 1 gives a 1D chain with identical cells. At t = 20 ms, amp pA (inward, entered positive)
		// for 10 ms into a central side x side square (side cells in 1D). Writes <prefix>_peak.csv and <prefix>_tact.csv (ny rows x nx
		// columns: peak Vm, first time > -30 mV or -1) and prints peak Vm along +x from the centre; "PROPAGATES" if the cell 20 away spikes.
		// Optional [gy] [sidey]: anisotropy. gSS couples along x (longitudinal), gy along y (transverse, default gSS); the stimulus is
		// side (x) by sidey (y) cells (default square). In 2D "PROPAGATES" requires the cells 20 away along +x and along +y to spike.
		int nx = atoi(argv[2]), ny = atoi(argv[3]); double g = atof(argv[4]), gK = atof(argv[5]);
		int side = atoi(argv[6]); double amp = atof(argv[7]), T = atof(argv[8]), dt = atof(argv[9]);
		std::string pre = argv[10];
		const double gy = argc > 11 ? atof(argv[11]) : g;
		const int sidey = argc > 12 ? atoi(argv[12]) : side;
		SmoothM rest; SmsInit(rest); rest.delT = dt;
		for (double t = 0; t < 20000; t += dt) { rest.Igap = gK * (rest.Vm + 80.0); rest.Integrate(0); }
		const int n = nx * ny, cx = nx / 2, cy = ny / 2;
		std::vector<SmoothM> c(n, rest);
		// Optional PDGFRa+ cells (env FLC_FP > 0): one modified PDGFa cell per SMC, weighted fP, gap FLC_GSP (default 3 nS), as in the
		// SIP cable. The stimulus then starts at SHEET_TSTIM ms (default 2000) so that SMC and FLC settle together; peaks and activation
		// times are only recorded from the stimulus onset. dt must be <= 0.02 ms (PDGFa step).
		const double fP = getenv("FLC_FP") ? atof(getenv("FLC_FP")) : 0.0, gSP = getenv("FLC_GSP") ? atof(getenv("FLC_GSP")) : 3.0;
		const double tStim = getenv("SHEET_TSTIM") ? atof(getenv("SHEET_TSTIM")) : (fP > 0 ? 2000.0 : 20.0);
		std::vector<PDGFRaCell> fl;
		std::vector<double> VP(n, 0.0);
		if (fP > 0)
		{
			PDGFRaCell proto(1); proto.Init(); proto.SetModified(true, true, 2.5, 1.0); proto.SetSKCaSource(1, 2.0e-3, 4);
			for (double t = 0; t < 10000; t += 0.02) proto.Step(0.02, 50.0 * (proto.Vm() + 65.0));
			// the sheet is uniform, so equilibrate one SMC + fP-weighted FLC pair (FLC_EQ ms, default 20000) and copy it everywhere
			const double teq = getenv("FLC_EQ") ? atof(getenv("FLC_EQ")) : 20000.0;
			SmoothM sp = rest; double vs = 0, vp = 0, nlast = 0;
			for (double t = 0; t < teq; t += dt)
			{
				const double v = sp.Vm, vq = proto.Vm();
				sp.Igap = fP * gSP * (v - vq) + gK * (v + 80.0); sp.Integrate(0);
				proto.Step(dt, gSP * (vq - v));
				if (t >= teq - 2000) { vs += v; vp += vq; nlast++; }
			}
			fprintf(stderr, "SMC+FLC pair after %.0f ms: SMC %.1f mV, FLC %.1f mV (mean of last 2 s)\n", teq, vs / nlast, vp / nlast);
			c = std::vector<SmoothM>(n, sp);	// SmoothM has const members (no copy assignment), so rebuild the vector
			fl.assign(n, proto);
			for (int i = 0; i < n; i++) fl[i].Reseed(1000u + i);
			// identical copies fire STHs nearly in step (the ensemble mean then sits about 3 mV depolarised for several seconds), so
			// advance every SMC+FLC pair on its own for a random 0..FLC_DESYNC ms (default 3000, STH period about 1.3 s)
			const double tds = getenv("FLC_DESYNC") ? atof(getenv("FLC_DESYNC")) : 3000.0;
#pragma omp parallel for schedule(dynamic, 16)
			for (int i = 0; i < n; i++)
			{
				std::mt19937 rg(424242u + i);
				const double tend = std::uniform_real_distribution<double>(0.0, tds)(rg);
				for (double t = 0; t < tend; t += dt)
				{
					const double v = c[i].Vm, vq = fl[i].Vm();
					c[i].Igap = fP * gSP * (v - vq) + gK * (v + 80.0); c[i].Integrate(0);
					fl[i].Step(dt, gSP * (vq - v));
				}
			}
		}
		std::vector<double> V(n), pk(n, -200), ta(n, -1);
		std::vector<char> stim(n, 0);
		for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++)
			if (abs(x - cx) * 2 < side + (side % 2 == 0) && (ny == 1 || abs(y - cy) * 2 < sidey + (sidey % 2 == 0))) stim[y * nx + x] = 1;
		int nstim = 0; for (char s : stim) nstim += s;
		const int i10 = cy * nx + std::min(cx + side / 2 + 10, nx - 1), i15 = cy * nx + std::min(cx + side / 2 + 15, nx - 1);	// 10 and 15 cells beyond the stimulus edge
		std::vector<double> tr10, tr15; long stp = 0; const long every = std::max(1L, (long)(0.1 / dt + 0.5));
		for (double t = 0; t < T; t += dt, stp++)
		{
			if (stp % every == 0) { tr10.push_back(c[i10].Vm); tr15.push_back(c[i15].Vm); }
			for (int i = 0; i < n; i++) V[i] = c[i].Vm;
			if (fP > 0) for (int i = 0; i < n; i++) VP[i] = fl[i].Vm();
			const bool on = t >= tStim && t < tStim + 10, rec = t >= tStim;
#pragma omp parallel for schedule(static)
			for (int i = 0; i < n; i++)
			{
				int x = i % nx, y = i / nx; double Ig = 0;
				if (x > 0) Ig += g * (V[i] - V[i - 1]);
				if (x < nx - 1) Ig += g * (V[i] - V[i + 1]);
				if (y > 0) Ig += gy * (V[i] - V[i - nx]);
				if (y < ny - 1) Ig += gy * (V[i] - V[i + nx]);
				if (fP > 0) { Ig += fP * gSP * (V[i] - VP[i]); fl[i].Step(dt, gSP * (VP[i] - V[i])); }
				c[i].Igap = Ig + gK * (V[i] + 80.0) - ((on && stim[i]) ? amp : 0.0);
				c[i].Integrate(0);
				if (rec && c[i].Vm > pk[i]) pk[i] = c[i].Vm;
				if (rec && ta[i] < 0 && c[i].Vm > -30) ta[i] = t - tStim;
			}
		}
		for (int k = 0; k < 2; k++)
		{
			FILE* fo = fopen((pre + (k ? "_tact.csv" : "_peak.csv")).c_str(), "w");
			for (int y = 0; y < ny; y++) { for (int x = 0; x < nx; x++) fprintf(fo, x ? ",%.2f" : "%.2f", k ? ta[y * nx + x] : pk[y * nx + x]); fprintf(fo, "\n"); }
			fclose(fo);
		}
		{
			FILE* fo = fopen((pre + "_trace.csv").c_str(), "w");
			fprintf(fo, "t_ms,V_edge10,V_edge15\n");
			for (size_t i = 0; i < tr10.size(); i++) fprintf(fo, "%.1f,%.3f,%.3f\n", i * every * dt, tr10[i], tr15[i]);
			fclose(fo);
		}
		const bool px = cx + 20 < nx && pk[cy * nx + cx + 20] > -20, py = ny == 1 || (cy + 20 < ny && pk[(cy + 20) * nx + cx] > -20);
		printf("sheet %dx%d gSS %.0f gK %.2f side %d (%d cells) amp %.0f: peak at d =", nx, ny, g, gK, side, nstim, amp);
		for (int d : {0, 1, 2, 3, 5, 10, 15, 20}) if (cx + d < nx) printf(" %d:%.1f", d, pk[cy * nx + cx + d]);
		if (ny > 1) printf(" | gy %.1f sidey %d, +y d=20: %.1f", gy, sidey, cy + 20 < ny ? pk[(cy + 20) * nx + cx] : NAN);
		printf("%s\n", (px && py) ? "  PROPAGATES" : (px ? "  x-only" : ""));
		return 0;
	}

	if (mode == "flcsingle")
	{
		// gisim flcsingle T_ms out.csv [skKd=mM] [ip3=mM] [ip3On=ms] [inj=pA] : isolated modified PDGFa cell as used in the tissue runs
		// (SK = EK, voltage independent, deterministic gate, NaLCN 2.5, SK senses submembrane Ca with Kd 2 uM, Hill 4), equilibrated 10 s
		// near -65 mV and 10 s free, then T_ms free with optional SK Kd (apamin: 1000), agonist IP3 from ip3On, constant current inj
		// (pA, outward positive). Writes t, Vm, ISK (pA), CaPMU (uM), Cai (uM) every 1 ms and prints mean/min/max Vm.
		double T = atof(argv[2]), dt = 0.02, kd = 2.0e-3, ip3 = 0, ip3On = 0, inj = 0;
		for (int a = 4; a < argc; a++) { std::string kv = argv[a]; size_t e = kv.find('='); if (e == std::string::npos) continue;
			std::string k = kv.substr(0, e); double v = atof(kv.c_str() + e + 1);
			if (k == "skKd") kd = v; else if (k == "ip3") ip3 = v; else if (k == "ip3On") ip3On = v; else if (k == "inj") inj = v; }
		PDGFRaCell p(1); p.Init();
		p.SetModified(true, true, 2.5, 1.0); p.SetSKCaSource(1, 2.0e-3, 4);
		for (double t = 0; t < 10000; t += dt) p.Step(dt, 50.0 * (p.Vm() + 65.0));
		for (double t = 0; t < 10000; t += dt) p.Step(dt, 0.0);
		p.SetSKCaSource(1, kd, 4);
		FILE* fo = fopen(argv[3], "w"); fprintf(fo, "t_ms,V,ISK,CaPMU_uM,Cai_uM\n");
		double sum = 0, mn = 1e9, mx = -1e9; long n = 0, k = 0, every = (long)(1.0 / dt + 0.5);
		for (double t = 0; t < T; t += dt, k++)
		{
			p.SetIP3Agonist(t >= ip3On ? ip3 : 0.0);
			p.Step(dt, inj);
			if (k % every == 0) fprintf(fo, "%.1f,%.3f,%.3f,%.5f,%.5f\n", t, p.Vm(), p.ISK(), 1e3 * p.CaPMU(), 1e3 * p.Cai());
			if (t >= ip3On) { sum += p.Vm(); n++; mn = std::min(mn, p.Vm()); mx = std::max(mx, p.Vm()); }
		}
		fclose(fo);
		printf("FLC single: skKd %g mM, ip3 %g mM from %g ms, inj %g pA: Vm mean %.2f, min %.2f, max %.2f mV (from ip3On)\n", kd, ip3, ip3On, inj, sum / n, mn, mx);
		return 0;
	}

	if (mode == "flcclamp")
	{
		// gisim flcclamp T_ms v1 v2 ... : modified PDGFa cell (as in the SIP runs: SK = EK, voltage independent, deterministic gate,
		// NaLCN 2.5, SK senses submembrane Ca, Kd 2 uM, Hill 4) equilibrated 10 s near -65 mV, then clamped (1000 nS) at each voltage
		// for T_ms; mean clamp current (includes spontaneous STH/SK events) and its SD. Gives the time-averaged FLC I-V (leak seen by SMC).
		double T = atof(argv[2]), dt = 0.02, gcl = 1000;
		PDGFRaCell p(1); p.Init();
		p.SetModified(true, true, 2.5, 1.0); p.SetSKCaSource(1, 2.0e-3, 4);
		for (double t = 0; t < 10000; t += dt) p.Step(dt, 50.0 * (p.Vm() + 65.0));
		PDGFRaCell p0 = p;
		for (double t = 0; t < 20000; t += dt) p0.Step(dt, 0.0);
		double vfree = 0; for (double t = 0; t < 20000; t += dt) { p0.Step(dt, 0.0); vfree += p0.Vm() * dt / 20000; }
		printf("FLC Cap %.1f pF, free-running mean Vm %.1f mV (20 s)\n", p.Cap(), vfree);
		for (int j = 3; j < argc; j++)
		{
			double vh = atof(argv[j]); PDGFRaCell c = p; double s = 0, s2 = 0; long n = 0;
			for (double t = 0; t < 3000; t += dt) c.Step(dt, gcl * (c.Vm() - vh));
			for (double t = 0; t < T; t += dt) { c.Step(dt, gcl * (c.Vm() - vh)); double I = -gcl * (c.Vm() - vh); s += I; s2 += I * I; n++; }	// membrane current = -clamp current
			double m = s / n;
			printf("Vh %6.1f mV: mean membrane current %8.2f pA (SD %.2f), mean ISK-free check Vm %.2f\n", vh, m, sqrt(std::max(0.0, s2 / n - m * m)), c.Vm());
		}
		return 0;
	}

	if (mode == "psd")
	{
		// gisim psd T dt out.csv [key=value ...]
		// Poststimulus depolarization (Koh et al. 2022 PNAS): N ICC (store-based ICC model, size s) each coupled by gc to one passive
		// SIP node (Csip, Gsip, Esip; all totals for the population). Nitrergic stimulation of frequency f for dur ms from tOn gives a drive
		// u = f/(f + f50). Two NO effects with first-order kinetics:
		//   n (cGMP, tauNon/tauNoff): inhibits IP3R release in ICC, SF_IP3R = 1 - aRel*n, and opens a hyperpolarizing conductance gNO*n
		//      (EK) in the SIP node (the nitrergic IJP is preserved in ICC-specific Ano1 knockdown, so it is not carried by ANO1);
		//   p (cGMP -> PDE3A -> cAMP/PKA -> PLN, tauPon/tauPoff): stimulates SERCA, SF_SERCA = 1 + bSerca*p.
		// ICCs are desynchronized by running each for a random extra time (0..desync ms) before coupling.
		// lum, lumRef: optional luminal Ca2+ sensitivity of IP3R release, x (CaSR/lumRef)^lum (0 = original ICC model).
		// Writes t_ms, V_SIP, V_ICCmean, I_ANO1 (total, pA), CaSR_mean, n, p every ms.
		std::map<std::string, double> P = { {"N", 30}, {"s", 0.72}, {"gc", 5}, {"ip3", 0.3}, {"Csip", 6000}, {"Gsip", 60}, {"Esip", -45},
			{"gNO", 45}, {"EK", -85}, {"tOn", 30000}, {"dur", 5000}, {"f", 5}, {"f50", 2}, {"aRel", 0.9}, {"bSerca", 1.0},
			{"tauNon", 300}, {"tauNoff", 1000}, {"tauPon", 2000}, {"tauPoff", 8000}, {"seed", 1}, {"desync", 20000}, {"pre", 20000}, {"lum", 0}, {"lumRef", 2.0} };
		double T = atof(argv[2]), dt = atof(argv[3]);
		for (int j = 5; j < argc; j++)
		{
			std::string a = argv[j]; size_t e = a.find('='); if (e == std::string::npos) continue;
			P[a.substr(0, e)] = atof(a.c_str() + e + 1);
		}
		const int N = (int)P["N"]; const double s = P["s"], gc = P["gc"];
		ICC proto; proto.InitialCellVariables(); proto.m_dDt = dt; proto.ip3ProdScale = P["ip3"]; proto.lumExp = P["lum"]; proto.lumRef = P["lumRef"];
		for (double t = 0; t < P["pre"]; t += dt) proto.RunCellModel(TRUE, 0.2 * (proto.Vm - P["Esip"]) / s);	// settle near the SIP potential
		std::mt19937 rng((unsigned)P["seed"]); std::uniform_real_distribution<double> U(0.0, 1.0);
		std::vector<ICC> cells(N, proto);
		for (auto& c : cells) { double tx = U(rng) * P["desync"]; for (double t = 0; t < tx; t += dt) c.RunCellModel(TRUE, 0.2 * (c.Vm - P["Esip"]) / s); }
		double V = P["Esip"], n = 0, p = 0;
		FILE* f = fopen(argv[4], "w");
		fprintf(f, "t_ms,V_SIP,V_ICCmean,I_ANO1,CaSR_mean,n,p\n");
		long step = 0, ev = std::max(1L, (long)(1.0 / dt + 0.5));
		const double u1 = P["f"] / (P["f"] + P["f50"]);
		for (double t = 0; t < T; t += dt, step++)
		{
			double u = (t >= P["tOn"] && t < P["tOn"] + P["dur"]) ? u1 : 0.0;
			n += dt * (u - n) / (u > n ? P["tauNon"] : P["tauNoff"]);
			p += dt * (u - p) / (u > p ? P["tauPon"] : P["tauPoff"]);
			double Igap = 0, Iano = 0, vm = 0, sr = 0;
			for (auto& c : cells)
			{
				c.SF_IP3R = 1.0 - P["aRel"] * n; c.SF_SERCA = 1.0 + P["bSerca"] * p;
				double ig = gc * (c.Vm - V);						// outward from the ICC into the SIP node
				c.RunCellModel(TRUE, ig / s);
				Igap += ig; Iano += s * c.IClCa; vm += c.Vm / N; sr += c.CaiSR / N;
			}
			V += -(P["Gsip"] * (V - P["Esip"]) + P["gNO"] * n * (V - P["EK"]) - Igap) / P["Csip"] * dt;
			if (step % ev == 0) fprintf(f, "%.0f,%.3f,%.3f,%.2f,%.5f,%.4f,%.4f\n", t, V, vm, Iano, sr, n, p);
		}
		fclose(f);
		printf("psd done: N %d, gc %.2f, ip3 %.3f, f %.1f Hz, dur %.0f ms, aRel %.2f, bSerca %.2f\n", N, gc, P["ip3"], P["f"], P["dur"], P["aRel"], P["bSerca"]);
		return 0;
	}

	if (mode == "iccclamp")
	{
		// gisim iccclamp gc Vh T dt out.csv [s=0.72 dV=-5 on=200 every=1000]
		// Sucrose-gap voltage clamp of a SIP node (Connor, Prosser & Weems 1974): the clamped muscle is an ideal voltage source Vc(t);
		// one ICC-MY (size scale s) hangs on it through g_c and is NOT clamped. Every 'every' ms, Vc steps by dV for 'on' ms
		// (slope-conductance probe). Writes t_ms, Vc, V_ICC, I_node (membrane-current convention of the node: -g_c (V_ICC - Vc),
		// inward negative, pA), gANO1 (ICC ANO1 conductance, nS, scaled by s). 20 s coupled settling before recording.
		double gc = atof(argv[2]), Vh = atof(argv[3]), T = atof(argv[4]), dt = atof(argv[5]);
		double s = 0.72, dV = -5, on = 200, every = 1000;
		for (int j = 7; j < argc; j++)
		{
			std::string a = argv[j]; size_t e = a.find('='); if (e == std::string::npos) continue;
			std::string k = a.substr(0, e); double v = atof(a.c_str() + e + 1);
			if (k == "s") s = v; else if (k == "dV") dV = v; else if (k == "on") on = v; else if (k == "every") every = v;
		}
		Equilibrate(Tpre, dt);
		ICC icc(gIccRest);
		if (gc > 0) dt = std::min(dt, 0.2 * icc.Cap * s / gc);
		icc.m_dDt = dt;
		FILE* f = fopen(argv[6], "w");
		fprintf(f, "t_ms,Vc,V_ICC,I_node,gANO1\n");
		long step = 0, ev = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = -20000; t < T; t += dt, step++)
		{
			double ph = fmod(t + 1e6 * every, every);
			double vc = Vh + ((t >= 0 && ph >= every - on) ? dV : 0.0);
			double Ig = gc * (icc.Vm - vc);
			icc.RunCellModel(TRUE, Ig / s);
			if (t >= 0 && step % ev == 0)
			{
				double den = icc.Vm - icc.mRTF * log(icc.ClPMU / icc.Cle);	// ECl as in ICC::IClCa
				fprintf(f, "%.1f,%.2f,%.4f,%.4f,%.5f\n", t, vc, icc.Vm, -Ig, fabs(den) > 1e-6 ? s * icc.IClCa / den : 0.0);
			}
		}
		fclose(f);
		printf("iccclamp gc %.3g nS, Vh %.1f mV, dt %.4g ms done\n", gc, Vh, dt);
		return 0;
	}

	if (mode == "iccprc")
	{
		// gisim iccprc amp dur nphase dt out.csv [s=0.72 vth=-45 ncyc=6 gL=0 EL=-65.5 trace=-1]
		// Phase response curve of the isolated ICC-MY (size scale s). Reference = an upstroke (Vm crossing vth upward) after 20 s
		// equilibration + 10 cycles; T0 = mean free period. For each phase k/nphase a rectangular current pulse (amp pA, entered positive =
		// inward/depolarizing, scaled like a coupling current: Iinj = -amp/s) of dur ms starts at k/nphase*T0 after the reference.
		// Runs (ncyc + 0.5) T0 per phase; writes every upstroke after the reference: phase (-1 = control), event no., time, peak Vm.
		double amp = atof(argv[2]), dur = atof(argv[3]); int np = atoi(argv[4]); double dt = atof(argv[5]);
		double s = 0.72, vth = -45, gL = 0, EL = -65.5, trph = -1; int ncyc = 6;
		for (int j = 7; j < argc; j++)
		{
			std::string a = argv[j]; size_t e = a.find('='); if (e == std::string::npos) continue;
			std::string k = a.substr(0, e); double v = atof(a.c_str() + e + 1);
			if (k == "s") s = v; else if (k == "vth") vth = v; else if (k == "ncyc") ncyc = (int)v;
			else if (k == "gL") gL = v; else if (k == "EL") EL = v; else if (k == "trace") trph = v;
		}
		// gL (nS) toward EL: a passive stand-in for the muscle column, g_c*G_load/(g_c+G_load); its current is scaled by 1/s like g_c.
		// trace = phase: also writes the Vm trace (1-ms samples) of that phase's run to out.csv.trace.csv
		Equilibrate(Tpre, dt);
		ICC icc(gIccRest); icc.m_dDt = dt;
		for (double t = 0; t < 20000; t += dt) icc.RunCellModel(TRUE, gL * (icc.Vm - EL) / s);
		// free run: collect upstrokes, keep the state at the 10th as the reference
		std::vector<double> up; std::unique_ptr<ICC> refp; double vPrev = icc.Vm;
		for (double t = 0; up.size() < 20; t += dt)
		{
			icc.RunCellModel(TRUE, gL * (icc.Vm - EL) / s);
			if (vPrev < vth && icc.Vm >= vth) { up.push_back(t); if (up.size() == 10) refp.reset(new ICC(icc)); }
			vPrev = icc.Vm;
			if (t > 600000) { fprintf(stderr, "no rhythm\n"); return 1; }
		}
		double T0 = (up[19] - up[9]) / 10.0;
		auto run = [&](double ton, double I, std::vector<double>& tn, std::vector<double>& pk, FILE* tf = nullptr)
		{
			ICC c(*refp); double vp = c.Vm; tn.clear(); pk.clear(); bool in = false; long st = 0, ev = std::max(1L, (long)(1.0 / dt + 0.5));
			for (double t = dt; t < (ncyc + 0.5) * T0; t += dt, st++)
			{
				c.RunCellModel(TRUE, ((t >= ton && t < ton + dur) ? -I / s : 0.0) + gL * (c.Vm - EL) / s);
				if (tf && st % ev == 0) fprintf(tf, "%.1f,%.3f\n", t, c.Vm);
				if (vp < vth && c.Vm >= vth && t > 0.2 * T0) { tn.push_back(t); pk.push_back(c.Vm); in = true; }	// skip the reference upstroke
				if (in) { if (c.Vm > pk.back()) pk.back() = c.Vm; if (c.Vm < vth) in = false; }
				vp = c.Vm;
			}
		};
		std::vector<double> t0n, p0n; run(-1e9, 0.0, t0n, p0n);
		FILE* f = fopen(argv[6], "w");
		fprintf(f, "phase,event,t_ms,peak_mV\n");
		for (size_t n = 0; n < t0n.size(); n++) fprintf(f, "-1,%zu,%.1f,%.2f\n", n + 1, t0n[n], p0n[n]);
		printf("iccprc amp %.1f pA dur %.0f ms: T0 %.1f ms (%.2f cpm); control events %zu\n", amp, dur, T0, 60000.0 / T0, t0n.size());
		for (int k = 0; k < np; k++)
		{
			double ph = (double)k / np; std::vector<double> tn, pk;
			FILE* tf = nullptr;
			if (trph >= 0 && fabs(ph - trph) < 0.5 / np) { tf = fopen((std::string(argv[6]) + ".trace.csv").c_str(), "w"); fprintf(tf, "t_ms,V\n"); }
			run(ph * T0, amp, tn, pk, tf);
			if (tf) fclose(tf);
			for (size_t n = 0; n < tn.size(); n++) fprintf(f, "%.4f,%zu,%.1f,%.2f\n", ph, n + 1, tn[n], pk[n]);
		}
		fclose(f);
		return 0;
	}

	if (mode == "smsfit")
	{
		// gisim smsfit gK [out.csv] : spike metrics of an isolated SmoothM (conductance scales from SMS_* env) with tonic K+ gK (nS).
		// Rin from a -5 pA step; 10-ms rheobase by bisection (all-or-none check: peak > -10 mV); then a 1.5x rheobase pulse:
		// peak, amplitude from rest, APD50 (width at half amplitude), time above -30 mV, max dV/dt, Vm 200 ms after onset.
		double gK = atof(argv[2]), dt = 0.01;
		SmoothM rest; SmsInit(rest); rest.delT = dt;
		for (double t = 0; t < 20000; t += dt) { rest.Igap = gK * (rest.Vm + 80.0); rest.Integrate(0); }
		auto run = [&](double amp, double dur, std::vector<double>* tr) {
			SmoothM c = rest; double pk = -200;
			for (double t = 0; t < 1500; t += dt)
			{
				c.Igap = gK * (c.Vm + 80.0) - ((t >= 100 && t < 100 + dur) ? amp : 0.0);
				c.Integrate(0); pk = std::max(pk, c.Vm);
				if (tr && fmod(t, 0.1) < dt) tr->push_back(c.Vm);
			}
			return pk; };
		std::vector<double> h; run(-5, 1000, &h);
		double Rin = (h[(size_t)(1000 / 0.1)] - rest.Vm) / -5.0 * 1000;	// MOhm, at 900 ms into the step
		double lo = 0, hi = 2000;
		if (run(hi, 10, NULL) < -10) { printf("SMS CaL %.2f Ks %.2f KCa %.2f leak %.2f gK %.2f: rest %.1f mV, Rin %.0f MOhm, not excitable up to 2 nA\n", rest.sfical1, rest.sfiks1, rest.sfikca1, rest.sfileak1, gK, rest.Vm, Rin); return 0; }
		for (int k = 0; k < 30; k++) { double mid = 0.5 * (lo + hi); (run(mid, 10, NULL) > -10 ? hi : lo) = mid; }
		std::vector<double> v; double pk = run(1.5 * hi, 10, &v);
		double amp = pk - rest.Vm, half = rest.Vm + 0.5 * amp, a50 = 0, a30 = 0, dvmax = 0;
		for (size_t i = 1100; i < v.size(); i++) { if (v[i] > half) a50 += 0.1; if (v[i] > -30) a30 += 0.1; dvmax = std::max(dvmax, (v[i] - v[i - 1]) / 0.1); }
		printf("SMS CaL %.2f Ks %.2f KCa %.2f leak %.2f gK %.2f: rest %.1f mV, Rin %.0f MOhm, rheobase(10 ms) %.0f pA; at 1.5x: peak %.1f mV, amp %.1f mV, APD50 %.0f ms, >-30 mV %.0f ms, dV/dt %.1f V/s, V(+200 ms) %.1f\n",
			rest.sfical1, rest.sfiks1, rest.sfikca1, rest.sfileak1, gK, rest.Vm, Rin, hi, pk, amp, a50, a30, dvmax, v[3000]);
		if (argc > 3) { FILE* fo = fopen(argv[3], "w"); fprintf(fo, "t_ms,V\n"); for (size_t i = 0; i < v.size(); i++) fprintf(fo, "%.1f,%.3f\n", 0.1 * i - 100, v[i]); fclose(fo); }
		return 0;
	}

	if (mode == "pohref")
	{
		// gisim pohref dtSub out.csv : the Poh et al. 2012 model with its ORIGINAL ICC drive (Gcouple 2.6 nS to a prescribed ICC
		// waveform, 10-s period, 40 s) from the CellML initial state, for comparison with paper_R1/check_poh_ref.py (LSODA).
		double h = atof(argv[2]);
		PohSMC c; c.delT = h; c.dtSub = h;
		FILE* fo = fopen(argv[3], "w"); fprintf(fo, "t_ms,V\n");
		long every = std::max(1L, (long)(1.0 / h + 0.5)), k = 0;
		for (double t = 0; t < 40000.0 - 1e-9; t += h, k++)
		{
			if (k % every == 0) fprintf(fo, "%.1f,%.6f\n", t, c.Vm);
			double loc = t - 10000.0 * floor(t / 10000.0);
			double vicc = loc < 300.0 ? -57.0 + 33.5 * loc / 300.0
				: -57.0 + 33.5 * (1.0 + exp(-12000.0 / 1200.0)) / (1.0 + exp((loc - 300.0 - 6000.0) / 600.0));
			c.Igap = 2.6 * (c.Vm - vicc);
			c.Integrate();
		}
		fclose(fo);
		return 0;
	}

	if (mode == "smsstep")
	{
		// gisim smsstep gK dur_ms out.csv amp1 amp2 ... : isolated SmoothM with tonic K+ gK (nS, toward -80 mV), equilibrated 20 s,
		// then one inward current pulse (amp pA, entered positive) of dur_ms per amplitude, from the same rest state.
		// Regenerative (spike) test: peak Vm, Vm at pulse end, and whether Vm keeps rising after the pulse ends.
		double gK = atof(argv[2]), dur = atof(argv[3]), dt = 0.01;
		SmoothM rest; SmsInit(rest); rest.delT = dt;
		for (double t = 0; t < 20000; t += dt) { rest.Igap = gK * (rest.Vm + 80.0); rest.Integrate(0); }
		int na = argc - 5;
		std::vector<std::vector<double>> vv(na);
		for (int j = 0; j < na; j++)
		{
			double amp = atof(argv[5 + j]);
			SmoothM c = rest; double vEnd = 0, vPk = -200, tPk = 0;
			for (double t = 0; t < 1000; t += dt)
			{
				c.Igap = gK * (c.Vm + 80.0) - ((t >= 100 && t < 100 + dur) ? amp : 0.0);
				c.Integrate(0);
				if (t < 100 + dur && t + dt >= 100 + dur) vEnd = c.Vm;
				if (c.Vm > vPk) { vPk = c.Vm; tPk = t; }
				if (fmod(t, 0.5) < dt) vv[j].push_back(c.Vm);
			}
			printf("gK %.2f dur %.0f ms amp %6.1f pA: rest %.1f, end of pulse %.1f, peak %.1f mV at %.1f ms after onset%s\n",
				gK, dur, amp, rest.Vm, vEnd, vPk, tPk - 100, tPk > 100 + dur + 1 ? "  (rises after pulse: regenerative)" : "");
		}
		FILE* fo = fopen(argv[4], "w");
		fprintf(fo, "t_ms"); for (int j = 0; j < na; j++) fprintf(fo, ",A%s", argv[5 + j]); fprintf(fo, "\n");
		for (size_t i = 0; i < vv[0].size(); i++) { fprintf(fo, "%.1f", 0.5 * i); for (int j = 0; j < na; j++) fprintf(fo, ",%.3f", vv[j][i]); fprintf(fo, "\n"); }
		fclose(fo);
		return 0;
	}

	if (mode == "single_aux")
	{
		// gisim single_aux T dt out.csv : isolated FastIM and mean-field PDGFRa (at rest, and with purinergic drive u=1 from T/3 to 2T/3)
		double T = atof(argv[2]), dt = atof(argv[3]);
		ActiveParams P;
		FastIM f; f.Init(P.vth);
		PDGFRa p0, p1;
		FILE* fo = fopen(argv[4], "w");
		fprintf(fo, "t_ms,V_fastIM,V_pdgf_rest,V_pdgf_drive,ISK_drive_pA\n");
		long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			f.Step(dt, 0.0, P.vth, P.tauR, P.gAno);
			p0.Step(dt, 0.0, 0.0);
			p1.Step(dt, 0.0, (t > T / 3 && t < 2 * T / 3) ? 1.0 : 0.0);
			if (step % every == 0) fprintf(fo, "%.1f,%.3f,%.3f,%.3f,%.3f\n", t, f.V, p0.V, p1.V, p1.ISK);
		}
		fclose(fo);
		return 0;
	}

	if (mode == "iccim_drive")
	{
		// gisim iccim_drive ip3Scale gload period_ms width_ms dV T dt [SF_SOC SF_VDDR SF_CaL SF_ClCa]
		// load voltage = -65 mV, stepped to -65+dV for width_ms every period_ms (first pulse at 10 s); prints events per pulse
		double s = atof(argv[2]), gl = atof(argv[3]), per = atof(argv[4]), wid = atof(argv[5]), dV = atof(argv[6]);
		double T = atof(argv[7]), dt = atof(argv[8]);
		ICC c; c.InitialCellVariables(); c.m_dDt = dt; c.ip3ProdScale = s;
		if (argc > 9)  c.SF_SOC = atof(argv[9]);
		if (argc > 10) c.SF_VDDR = atof(argv[10]);
		if (argc > 11) c.SF_CaL = atof(argv[11]);
		if (argc > 12) c.SF_ClCa = atof(argv[12]);
		std::vector<double> ev; bool above = false; int pulses = 0;
		for (double t = 0; t < T; t += dt)
		{
			bool on = t >= 10000 && fmod(t - 10000, per) < wid;
			if (on && fmod(t - 10000, per) < dt) pulses++;
			c.RunCellModel(TRUE, gl * (c.Vm - (on ? -65 + dV : -65)));
			if (!above && c.Vm > -35) { above = true; ev.push_back(t); }
			if (above && c.Vm < -55) above = false;
		}
		int before = 0, after = 0;
		for (double e : ev) (e < 10000 ? before : after)++;
		printf("s %.3f SOC %.2f VDDR %.2f CaL %.2f ClCa %.2f dV %.0f: spontaneous events (0-10 s) %d, events during drive %d / %d pulses; first events:",
			s, c.SF_SOC, c.SF_VDDR, c.SF_CaL, c.SF_ClCa, dV, before, after, pulses);
		for (size_t i = 0; i < ev.size() && i < 8; i++) printf(" %.1f", ev[i] / 1000);
		printf("\n");
		return 0;
	}

	if (mode == "dmp")
	{
		// gisim dmp ip3Scale SF_VDDR SF_CaL gload T dt out.csv : one ICC with voltage-dependent Ca2+ entry scaled (0 = ICC-DMP candidate),
		// loaded by gload (nS) toward -65 mV (a crude stand-in for the syncytium). Ca2+ events = Cai upward crossings of
		// (min + max)/2 over the second half; prints event count, mean interval and CV of intervals (CV ~0: clock-like, ~1: Poisson-like).
		double s = atof(argv[2]), T = atof(argv[6]), dt = atof(argv[7]), gl = atof(argv[5]);
		ICC c; c.InitialCellVariables(); c.m_dDt = dt; c.ip3ProdScale = s;
		c.SF_VDDR = atof(argv[3]); c.SF_CaL = atof(argv[4]);
		Trace tr; long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			c.RunCellModel(TRUE, gl * (c.Vm + 65.0));
			if (step % every == 0) { tr.t.push_back(t); tr.vI.push_back(c.Vm); tr.vS.push_back(c.IClCa); tr.caI.push_back(c.Cai); }
		}
		WriteTrace(argv[8], tr);
		size_t h = tr.t.size() / 2;
		double cmin = 1e9, cmax = -1e9, vmin = 1e9, vmax = -1e9;
		for (size_t i = h; i < tr.t.size(); i++) { cmin = std::min(cmin, tr.caI[i]); cmax = std::max(cmax, tr.caI[i]); vmin = std::min(vmin, tr.vI[i]); vmax = std::max(vmax, tr.vI[i]); }
		std::vector<double> ev;
		double th = 0.5 * (cmin + cmax);
		if (cmax - cmin > 0.05 * cmax) for (size_t i = h + 1; i < tr.t.size(); i++) if (tr.caI[i - 1] < th && tr.caI[i] >= th) ev.push_back(tr.t[i]);
		double mi = 0, sd = 0; int n = (int)ev.size() - 1;
		for (int i = 0; i < n; i++) mi += (ev[i + 1] - ev[i]) / n;
		for (int i = 0; i < n; i++) sd += (ev[i + 1] - ev[i] - mi) * (ev[i + 1] - ev[i] - mi) / std::max(1, n - 1);
		printf("s %.3f VDDR %.2f CaL %.2f gload %.1f: Vm %.1f..%.1f mV, Cai %.2e..%.2e mM, events %d (2nd half), interval %.2f s, CV %.3f\n",
			s, c.SF_VDDR, c.SF_CaL, gl, vmin, vmax, cmin, cmax, (int)ev.size(), mi / 1000, n > 1 ? sqrt(sd) / mi : NAN);
		return 0;
	}

	if (mode == "iccim")
	{
		// gisim iccim ip3Scale gload Vload T dt out.csv : one ICC loaded by gload (nS) toward Vload (mV)
		double s = atof(argv[2]), gl = atof(argv[3]), Vl = atof(argv[4]), T = atof(argv[5]), dt = atof(argv[6]);
		ICC c; c.InitialCellVariables(); c.m_dDt = dt; c.ip3ProdScale = s;
		if (argc > 8) c.SF_VDDR = atof(argv[8]);	// optional: ICC-IM uses SF_VDDR = 10
		Trace tr; long step = 0, every = std::max(1L, (long)(1.0 / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			c.RunCellModel(TRUE, gl * (c.Vm - Vl));
			if (step % every == 0) { tr.t.push_back(t); tr.vI.push_back(c.Vm); tr.vS.push_back(c.IClCa); tr.caI.push_back(c.Cai); }
		}
		WriteTrace(argv[7], tr);
		std::vector<double> th(tr.t.begin() + tr.t.size() / 2, tr.t.end()), vh(tr.vI.begin() + tr.vI.size() / 2, tr.vI.end());
		double f, a, mx, mn;
		Analyze(th, vh, f, a, mx, mn);
		printf("ip3Scale %.2f gload %.1f: %.2f cpm, amp %.1f mV, range %.1f..%.1f\n", s, gl, f, a, mn, mx);
		return 0;
	}

	if (mode == "cluster")
	{
		// gisim cluster T dt out.csv [key=value ...]: N ICC (no muscle) coupled as a chain (topo 0) or all-to-all (topo 1) by gII nS;
		// electrode in cell 0. Current clamp: inj pA (outward positive) in [injOn, injOff). Voltage clamp (vc 1): electrode conductance
		// 1000/Rs nS (Rs in MOhm) to Vh, stepped by dV mV in [tOn, tOn + tDur). Cell size: currents and Cap x s (s 0.22 = 5.5 pF).
		// ICC keys as in 'active': iccVDDR, iccVH, ip3VK, iccRest (iccNaLCN, iccK1, iccClb, iccNaK, iccNKCC). smp: sampling interval (ms).
		double T = atof(argv[2]), dt = atof(argv[3]); const char* out = argv[4];
		int N = 1, topo = 0, vc = 0;
		double gII = 1, s = 0.22, inj = 0, injOn = -1, injOff = -1, Rs = 10, Vh = -70, dV = -5, tOn = 100, tDur = 100, smp = 1;
		for (int a = 5; a < argc; a++)
		{
			std::string kv = argv[a]; size_t e = kv.find('='); if (e == std::string::npos) continue;
			std::string k = kv.substr(0, e); double v = atof(kv.c_str() + e + 1);
			if (k == "N") N = (int)v; else if (k == "topo") topo = (int)v; else if (k == "vc") vc = (int)v; else if (k == "gII") gII = v;
			else if (k == "s") s = v; else if (k == "inj") inj = v; else if (k == "injOn") injOn = v; else if (k == "injOff") injOff = v;
			else if (k == "Rs") Rs = v; else if (k == "Vh") Vh = v; else if (k == "dV") dV = v; else if (k == "tOn") tOn = v; else if (k == "tDur") tDur = v;
			else if (k == "smp") smp = v; else if (k == "iccVDDR") gIccVddr = v; else if (k == "ip3VK") gIccIp3VK = v; else if (k == "iccVH") gIccVH = v;
			else if (k == "iccRest") { for (double& x : gIccRestSF) x = v; } else if (k == "iccNaLCN") gIccRestSF[0] = v; else if (k == "iccK1") gIccRestSF[1] = v;
			else if (k == "iccClb") gIccRestSF[2] = v; else if (k == "iccNaK") gIccRestSF[3] = v; else if (k == "iccNKCC") gIccRestSF[4] = v;
			else if (k == "iccSOC") gIccSOC = v; else if (k == "iccANO1") gIccANO1 = v;
		}
		Equilibrate(Tpre, 0.1);
		std::vector<ICC> c(N, gIccRest);
		for (auto& x : c) x.m_dDt = dt;
		FILE* f = fopen(out, "w"); fprintf(f, "t_ms,V0,Vlast,Iel_pA\n");
		std::vector<double> V(N), I(N);
		long step = 0, every = std::max(1L, (long)(smp / dt + 0.5));
		for (double t = 0; t < T; t += dt, step++)
		{
			for (int i = 0; i < N; i++) V[i] = c[i].Vm;
			for (int i = 0; i < N; i++)
			{
				I[i] = 0;
				if (topo == 1) { for (int j = 0; j < N; j++) if (j != i) I[i] += gII * (V[i] - V[j]); }
				else { if (i > 0) I[i] += gII * (V[i] - V[i - 1]); if (i < N - 1) I[i] += gII * (V[i] - V[i + 1]); }
			}
			double Iel = 0;	// outward from cell 0 through the electrode (pA)
			if (vc) Iel = (V[0] - (Vh + ((t >= tOn && t < tOn + tDur) ? dV : 0.0))) * 1000.0 / Rs;
			else if (t >= injOn && t < injOff) Iel = inj;
			I[0] += Iel;
			for (int i = 0; i < N; i++) c[i].RunCellModel(TRUE, I[i] / s);
			if (step % every == 0) fprintf(f, "%.4f,%.4f,%.4f,%.5f\n", t, V[0], V[N - 1], -Iel);
		}
		fclose(f);
		printf("cluster N %d topo %d gII %.3g s %.3g vc %d done\n", N, topo, gII, s, vc);
		return 0;
	}

	if (mode == "active" || mode == "asweep")
	{
		// gisim active gc uP aI T dt out.csv  [key=value ...]
		// gisim asweep T dt summary.csv       [key=value ...]
		// keys: m gSS fP fI gSP gSI sIM vIM L gK0 imModel vth tauR gAno
		ActiveParams P;
		int first = (mode == "active") ? 8 : 5;
		for (int a = first; a < argc; a++)
		{
			std::string kv = argv[a];
			size_t e = kv.find('=');
			if (e == std::string::npos) continue;
			std::string k = kv.substr(0, e); double v = atof(kv.c_str() + e + 1);
			if (k == "m") P.m = v; else if (k == "gSS") P.gSS = v; else if (k == "fP") P.fP = v; else if (k == "fI") P.fI = v;
			else if (k == "gSP") P.gSP = v; else if (k == "gSI") P.gSI = v; else if (k == "sIM") P.sIM = v; else if (k == "vIM") P.vIM = v; else if (k == "imModel") P.imModel = (int)v; else if (k == "vth") P.vth = v; else if (k == "tauR") P.tauR = v; else if (k == "gAno") P.gAno = v; else if (k == "pdModel") P.pdModel = (int)v; else if (k == "gK0run") P.gK0run = v; else if (k == "seed") P.seed = (int)v; else if (k == "ip3AgMax") P.ip3AgMax = v; else if (k == "skKd") P.skKd = v; else if (k == "skHill") P.skHill = v; else if (k == "pdNaLCN") P.pdNaLCN = v; else if (k == "smcModel") P.smcModel = (int)v; else if (k == "uOn") P.uOn = v; else if (k == "uOff") P.uOff = v; else if (k == "iccInj") P.iccInj = v; else if (k == "iccGL") P.iccGL = v; else if (k == "wavePer") P.wavePer = v; else if (k == "waveHW") P.waveHW = v; else if (k == "waveRise") P.waveRise = v; else if (k == "waveFall") P.waveFall = v; else if (k == "waveRest") P.waveRest = v; else if (k == "waveAmp") P.waveAmp = v; else if (k == "iccEL") P.iccEL = v; else if (k == "injOn") P.injOn = v; else if (k == "injOff") P.injOff = v; else if (k == "smcInj") P.smcInj = v; else if (k == "iccScale") P.iccScale = v; else if (k == "injN") P.injN = (int)v; else if (k == "dumpPeak") P.dumpPeak = (int)v; else if (k == "traceAll") P.traceAll = (int)v; else if (k == "gT") P.gT = v; else if (k == "ET") P.ET = v; else if (k == "vmT") P.vmT = v; else if (k == "kmT") P.kmT = v; else if (k == "vhT") P.vhT = v; else if (k == "khT") P.khT = v; else if (k == "tauM") P.tauM = v; else if (k == "tauH") P.tauH = v; else if (k == "gKATP") P.gKATP = v; else if (k == "iccVDDR") gIccVddr = v; else if (k == "ip3VK") gIccIp3VK = v; else if (k == "iccVH") gIccVH = v; else if (k == "iccRest") { for (double& x : gIccRestSF) x = v; } else if (k == "iccNaLCN") gIccRestSF[0] = v; else if (k == "iccK1") gIccRestSF[1] = v; else if (k == "iccClb") gIccRestSF[2] = v; else if (k == "iccNaK") gIccRestSF[3] = v; else if (k == "iccNKCC") gIccRestSF[4] = v; else if (k == "iccSOC") gIccSOC = v; else if (k == "iccANO1") gIccANO1 = v; else if (k == "EKATP") P.EKATP = v; else if (k == "LL") P.LL = (int)v; else if (k == "gcL") P.gcL = v; else if (k == "mL") P.mL = v; else if (k == "gSSL") P.gSSL = v; else if (k == "gIF") P.gIF = v; else if (k == "fP0") P.fP0 = v;
			else if (k == "dmpN") P.dmpN = (int)v; else if (k == "dmpK") P.dmpK = (int)v; else if (k == "dmpRate") P.dmpRate = v; else if (k == "dmpG") P.dmpG = v; else if (k == "dmpTr") P.dmpTr = v; else if (k == "dmpTd") P.dmpTd = v; else if (k == "dmpCm") P.dmpCm = v; else if (k == "dmpGL") P.dmpGL = v; else if (k == "dmpEL") P.dmpEL = v;
			else if (k == "L") P.L = (int)v; else if (k == "gK0") gK0 = v;
		}
		if (mode == "active")
		{
			double gc = atof(argv[2]), uP = atof(argv[3]), aI = atof(argv[4]), T = atof(argv[5]), dt = atof(argv[6]);
			Equilibrate(Tpre, dt);
			EquilibrateActive(P, Tpre, dt);
			PrintCable(RunActive(P, gc, uP, aI, T, dt, argv[7]));
			return 0;
		}
		double T = atof(argv[2]), dt = atof(argv[3]);
		const char* out = argv[4];
		std::vector<double> gcs = { 0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500 };
		// neural axis: negative = purinergic drive uP, positive = cholinergic ICC-IM IP3 gain aI
		std::vector<double> col = { -1, -0.5, -0.3, -0.2, -0.1, 0, 0.2, 0.5, 1, 2, 5 };
		Equilibrate(Tpre, dt);
		EquilibrateActive(P, Tpre, dt);
		int n = (int)(gcs.size() * col.size());
		std::vector<CableSummary> res(n);
#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < n; k++)
		{
			double c = col[k % col.size()];
			res[k] = RunActive(P, gcs[k / col.size()], c < 0 ? -c : 0.0, c > 0 ? c : 0.0, T, dt, NULL);
		}
		FILE* f = fopen(out, "w");
		fprintf(f, "gc_nS,m,gSS,fP,fI,gSP,gSI,sIM,vIM,imModel,vth,tauR,gAno,gK0,neural,freqICC_cpm,ampICC_mV,freqL0_cpm,lambdaE_layers,depth20_layers,meanIgapICC_pA");
		for (int i = 0; i < P.L; i++) fprintf(f, ",amp_L%d", i);
		fprintf(f, "\n");
		for (int k = 0; k < n; k++)
		{
			const CableSummary& s = res[k];
			fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%.3f,%.3f,%.3f,%.3f,%d,%.3f", gcs[k / col.size()], P.m, P.gSS, P.fP, P.fI, P.gSP, P.gSI, P.sIM, P.vIM, P.imModel, P.vth, P.tauR, P.gAno, gK0,
				col[k % col.size()], s.freqI, s.ampI, s.freqS0, s.lambdaE, s.depth20, s.meanIgapI);
			for (int i = 0; i < P.L; i++) fprintf(f, ",%.3f", s.ampL[i]);
			fprintf(f, "\n");
		}
		fclose(f);
		printf("wrote %s (%d runs)\n", out, n);
		return 0;
	}

	if (mode == "cable")
	{
		// gisim cable gc m gSS L gn En T dt out.csv [gK0]
		double gc = atof(argv[2]), m = atof(argv[3]), gSS = atof(argv[4]);
		int L = atoi(argv[5]);
		double gn = atof(argv[6]), En = atof(argv[7]), T = atof(argv[8]), dt = atof(argv[9]);
		if (argc > 11) gK0 = atof(argv[11]);
		Equilibrate(Tpre, dt);
		PrintCable(RunCable(gc, m, gSS, L, gn, En, T, dt, argv[10]));
		return 0;
	}

	if (mode == "csweep")
	{
		// gisim csweep m gSS L T dt summary.csv [gK0]
		double m = atof(argv[2]), gSS = atof(argv[3]);
		int L = atoi(argv[4]);
		double T = atof(argv[5]), dt = atof(argv[6]);
		const char* out = argv[7];
		if (argc > 8) gK0 = atof(argv[8]);
		std::vector<double> gcs = { 0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500 };
		std::vector<double> gns = { -2, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 2 };
		Equilibrate(Tpre, dt);
		int n = (int)(gcs.size() * gns.size());
		std::vector<CableSummary> res(n);
#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < n; k++)
		{
			double g = gns[k % gns.size()];
			res[k] = RunCable(gcs[k / gns.size()], m, gSS, L, fabs(g), g < 0 ? -80.0 : 0.0, T, dt, NULL);
		}
		FILE* f = fopen(out, "w");
		fprintf(f, "gc_nS,m,gSS,L,gK0,gnSigned,freqICC_cpm,ampICC_mV,freqL0_cpm,lambdaE_layers,depth20_layers,meanIgapICC_pA");
		for (int i = 0; i < L; i++) fprintf(f, ",amp_L%d", i);
		fprintf(f, "\n");
		for (int k = 0; k < n; k++)
		{
			const CableSummary& s = res[k];
			fprintf(f, "%g,%g,%g,%d,%g,%g,%.3f,%.3f,%.3f,%.3f,%d,%.3f", gcs[k / gns.size()], m, gSS, L, gK0, gns[k % gns.size()],
				s.freqI, s.ampI, s.freqS0, s.lambdaE, s.depth20, s.meanIgapI);
			for (int i = 0; i < L; i++) fprintf(f, ",%.3f", s.ampL[i]);
			fprintf(f, "\n");
		}
		fclose(f);
		printf("wrote %s (%d runs)\n", out, n);
		return 0;
	}

	fprintf(stderr, "unknown mode\n");
	return 1;
}
