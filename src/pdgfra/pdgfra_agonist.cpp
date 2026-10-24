// Modified PDGFa cell held near the syncytium potential; purinergic agonist (IP3 signal) applied from tOn to tOff.
// Reports mean current into the syncytium (I < 0 = inhibitory), ISK and Cai before, during and after the agonist.
// Usage: pdgfra_agonist Vhold_mV g_nS ip3Ag_mM tOn_s tOff_s T_s seed out.csv [skKd_uM skHill]
//   skKd_uM > 0 switches SK to sense submembrane (PMU) Ca2+ with that Kd and Hill coefficient
#include "PDGFRaCell.h"
#include "CellPart.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
	double Vh = atof(argv[1]), g = atof(argv[2]), ag = atof(argv[3]);
	double tOn = atof(argv[4]) * 1000, tOff = atof(argv[5]) * 1000, T = atof(argv[6]) * 1000;
	unsigned seed = (unsigned)atoi(argv[7]);
	const double dt = 0.02;
	PDGFRaCell c(seed);
	c.Init();
	c.SetModified(true, true, 2.0, 1.0);
	double kd = argc > 9 ? atof(argv[9]) : 0.0, hill = argc > 10 ? atof(argv[10]) : 4.0;
	if (kd > 0) c.SetSKCaSource(1, kd * 1e-3, hill);
	FILE* f = fopen(argv[8], "w");
	fprintf(f, "t_ms,Vm,Icoupling_pA,ISK_pA,Cai_mM,CaPMU_mM,Ki_mM,Nai_mM\n");
	double sI[3] = { 0 }, sSK[3] = { 0 }, sCa[3] = { 0 }; long n[3] = { 0 };
	long step = 0, every = (long)(1.0 / dt + 0.5);
	for (double t = 0; t < T; t += dt, step++)
	{
		c.SetIP3Agonist(t >= tOn && t < tOff ? ag : 0.0);
		double I = g * (c.Vm() - Vh);
		c.Step(dt, I);
		if (step % every == 0)
		{
			fprintf(f, "%.1f,%.3f,%.3f,%.3f,%.4e,%.4e,%.3f,%.3f\n", t, c.Vm(), I, c.ISK(), c.Cai(), c.CaPMU(), c.Raw().Ki, c.Raw().Nai);
			// windows: baseline 5 s..tOn, agonist tOn+2 s..tOff, washout tOff+5 s..T
			int w = (t > 5000 && t < tOn) ? 0 : (t > tOn + 2000 && t < tOff) ? 1 : (t > tOff + 5000) ? 2 : -1;
			if (w >= 0) { sI[w] += I; sSK[w] += c.ISK(); sCa[w] += c.Cai(); n[w]++; }
		}
	}
	fclose(f);
	const char* lab[3] = { "baseline", "agonist ", "washout " };
	for (int w = 0; w < 3; w++)
		printf("Kd %.2f uM | ip3Ag %.4f mM | %s: I->syncytium %+7.1f pA  ISK %6.1f pA  mean Cai %.3f uM\n", kd, ag, lab[w], sI[w] / n[w], sSK[w] / n[w], sCa[w] / n[w] * 1e3);
	return 0;
}
