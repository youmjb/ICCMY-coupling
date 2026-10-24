// pdgfra_test.cpp : tests for the headless PDGFRa cell (PDGFRaCell wrapper)
// usage: pdgfra_test [part] [dt_ms] [T_s]
//   part = a | b | c | d | all (default all); dt default 0.02 ms (part a runs 0.02/0.05/0.1
//   unless dt is given); T_s = part-(a) simulated duration (default 30 s)
#include "PDGFRaCell.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <omp.h>

struct Stats
{
	double n = 0, sVm = 0, sVm2 = 0, vmin = 1e9, vmax = -1e9, sISK = 0, cmin = 1e9, cmax = -1e9, sCai = 0, sIP3 = 0;
	void add(const PDGFRaCell& c)
	{
		double v = c.Vm();
		n++; sVm += v; sVm2 += v * v; sISK += c.ISK(); sCai += c.Cai(); sIP3 += c.IP3();
		if (v < vmin) vmin = v; if (v > vmax) vmax = v;
		if (c.Cai() < cmin) cmin = c.Cai(); if (c.Cai() > cmax) cmax = c.Cai();
	}
	double mVm() const { return sVm / n; }
	double sdVm() const { double m = mVm(); return std::sqrt(std::fmax(0.0, sVm2 / n - m * m)); }
	double mISK() const { return sISK / n; }
	double mCai() const { return sCai / n; }
	double mIP3() const { return sIP3 / n; }
};

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// run for T ms, collecting stats (sampled every step) when t >= tStat
static void run(PDGFRaCell& c, double dt, double T, double Iext, Stats* st, double tStatFrom, bool* ok, bool perSecond = false)
{
	long long n = (long long)std::llround(T / dt);
	Stats sec;
	double t0 = c.Time();
	for (long long i = 0; i < n; ++i)
	{
		if (!c.Step(dt, Iext)) { if (ok) *ok = false; }
		double t = c.Time() - t0;
		if (st && t >= tStatFrom) st->add(c);
		if (perSecond)
		{
			sec.add(c);
			if (((i + 1) % (long long)std::llround(1000.0 / dt)) == 0)
			{
				std::printf("    t=%5.1f s  Vm mean %7.2f [%7.2f,%7.2f]  ISK %7.3f pA  Cai %.3e..%.3e mM  IP3 %.3e\n",
					c.Time() / 1000.0, sec.mVm(), sec.vmin, sec.vmax, sec.mISK(), sec.cmin, sec.cmax, sec.mIP3());
				sec = Stats();
			}
		}
	}
}

static void partA(double dtOnly, double Ts)
{
	std::printf("=== (a) single cell, Iext=0, %.0f s; stats over last %.0f s ===\n", Ts, Ts / 2);
	std::vector<double> dts = { 0.02, 0.05, 0.1 };
	if (dtOnly > 0) dts = { dtOnly };
	for (double dt : dts)
	{
		PDGFRaCell c(12345);
		c.Init();
		Stats st; bool ok = true;
		double w0 = now();
		run(c, dt, Ts * 1000.0, 0.0, &st, Ts * 500.0, &ok, true);
		double w = now() - w0;
		std::printf("  dt=%.3f ms: Vm min %.3f max %.3f mean %.3f (sd %.3f) mV | mean ISK %.4f pA | Cai %.4e..%.4e mM | IP3 mean %.4e | retries %lld | ok=%d | wall %.2f s per sim s\n",
			dt, st.vmin, st.vmax, st.mVm(), st.sdVm(), st.mISK(), st.cmin, st.cmax, st.mIP3(), c.Retries(), (int)ok, w / Ts);
		std::fflush(stdout);
	}
}

static PDGFRaCell equilibrated(double dt, double Teq_ms, unsigned seed)
{
	PDGFRaCell c(seed);
	c.Init();
	bool ok = true;
	run(c, dt, Teq_ms, 0.0, nullptr, 0, &ok);
	return c;
}

static void partB(double dt)
{
	std::printf("=== (b) input resistance: equilibrate 15 s, then Iext = 0/-10/+10 pA for 10 s (stats last 5 s), dt=%.3f ===\n", dt);
	PDGFRaCell base = equilibrated(dt, 15000.0, 777);
	double I[3] = { 0.0, -10.0, 10.0 };
	Stats st[3];
	#pragma omp parallel for schedule(static)
	for (int k = 0; k < 3; ++k)
	{
		PDGFRaCell c(base);   // same seed/RNG state -> only Iext differs
		bool ok = true;
		run(c, dt, 10000.0, I[k], &st[k], 5000.0, &ok);
	}
	for (int k = 0; k < 3; ++k)
		std::printf("  Iext=%+6.1f pA: Vm mean %.3f [%.3f,%.3f] mV, ISK %.4f pA, Cai %.4e mM\n", I[k], st[k].mVm(), st[k].vmin, st[k].vmax, st[k].mISK(), st[k].mCai());
	double Rin_dep = (st[1].mVm() - st[0].mVm()) / (-10.0) * -1.0;  // dV/dI_inward
	double Rin_hyp = (st[2].mVm() - st[0].mVm()) / (10.0) * -1.0;
	double Rin = (st[1].mVm() - st[2].mVm()) / 20.0;          // mV/pA = GOhm
	std::printf("  Rin (from -10 pA) %.3f GOhm, (from +10 pA) %.3f GOhm, symmetric %.3f GOhm (= %.1f MOhm)\n", Rin_dep, Rin_hyp, Rin, Rin * 1000.0);
}

static void partC(double dt)
{
	std::printf("=== (c) agonist: equilibrate 15 s, then SetAgonist(s) for 15 s (stats last 10 s), dt=%.3f ===\n", dt);
	PDGFRaCell base = equilibrated(dt, 15000.0, 4242);
	double S[4] = { 0.0, 1.0, 3.0, 0.0 };  // last = agonist 3 with SF_SK=0 (SK block)
	Stats st[4];
	#pragma omp parallel for schedule(static)
	for (int k = 0; k < 4; ++k)
	{
		PDGFRaCell c(base);
		c.SetAgonist(k == 3 ? 3.0 : S[k]);
		if (k == 3) c.SetSKScale(0.0);
		bool ok = true;
		run(c, dt, 15000.0, 0.0, &st[k], 5000.0, &ok);
	}
	const char* lab[4] = { "s=0 (control)", "s=1", "s=3", "s=3, SK blocked" };
	for (int k = 0; k < 4; ++k)
		std::printf("  %-16s Vm mean %.3f [%.3f,%.3f] mV  dVm %+.3f | ISK %.4f pA (x%.2f) | Cai %.4e mM | IP3 %.4e mM\n", lab[k],
			st[k].mVm(), st[k].vmin, st[k].vmax, st[k].mVm() - st[0].mVm(), st[k].mISK(), st[k].mISK() / st[0].mISK(), st[k].mCai(), st[k].mIP3());
}

static void partD(double dt)
{
	std::printf("=== (d) copy semantics: equilibrate 10 s, copy, step 3 s, dt=%.3f ===\n", dt);
	PDGFRaCell base = equilibrated(dt, 10000.0, 99);
	PDGFRaCell A(base), B(base);     // same RNG state
	PDGFRaCell C(base), D(base);
	C.Reseed(1001); D.Reseed(2002);
	PDGFRaCell* cells[4] = { &A, &B, &C, &D };
	long long n = (long long)std::llround(3000.0 / dt);
	std::vector<double> vm[4];
	for (auto& v : vm) v.resize((size_t)n);
	#pragma omp parallel for schedule(static)
	for (int k = 0; k < 4; ++k)
		for (long long i = 0; i < n; ++i) { cells[k]->Step(dt, 0.0); vm[k][(size_t)i] = cells[k]->Vm(); }
	auto maxdiff = [&](int a, int b) { double m = 0; for (long long i = 0; i < n; ++i) m = std::fmax(m, std::fabs(vm[a][i] - vm[b][i])); return m; };
	auto ms = [&](int a, double& mean, double& sd) { double s = 0, s2 = 0; for (double v : vm[a]) { s += v; s2 += v * v; } mean = s / n; sd = std::sqrt(std::fmax(0.0, s2 / n - mean * mean)); };
	std::printf("  same seed (A vs B): max|dVm| = %.3e mV  -> %s\n", maxdiff(0, 1), maxdiff(0, 1) == 0.0 ? "IDENTICAL" : "DIFFERENT");
	std::printf("  reseeded  (C vs D): max|dVm| = %.3e mV;  A vs C: %.3e mV\n", maxdiff(2, 3), maxdiff(0, 2));
	for (int k = 0; k < 4; ++k) { double m, s; ms(k, m, s); std::printf("  cell %c: Vm mean %.3f sd %.3f mV, final Vm %.4f\n", "ABCD"[k], m, s, vm[k][n - 1]); }
	// the base cell must be unaffected by stepping the copies
	std::printf("  base cell time %.1f ms (should be 10000.0), Vm %.4f\n", base.Time(), base.Vm());
}

int main(int argc, char** argv)
{
	std::string part = argc > 1 ? argv[1] : "all";
	double dt = argc > 2 ? std::atof(argv[2]) : 0.0;
	double Ts = argc > 3 ? std::atof(argv[3]) : 30.0;
	std::printf("pdgfra_test: part=%s dt=%g, OpenMP threads=%d\n", part.c_str(), dt, omp_get_max_threads());
	double dtB = dt > 0 ? dt : 0.02;
	if (part == "a" || part == "all") partA(dt, Ts);
	if (part == "b" || part == "all") partB(dtB);
	if (part == "c" || part == "all") partC(dtB);
	if (part == "d" || part == "all") partD(dtB);
	return 0;
}
