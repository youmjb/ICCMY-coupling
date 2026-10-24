// PDGFa cell (original or modified) either isolated (g = 0) or held near a syncytium potential through
// a large coupling conductance g. The coupling current I = g*(Vm - Vhold) is what this cell would inject
// into the syncytium: I > 0 depolarizes it (excitatory load), I < 0 hyperpolarizes it (inhibitory).
// Usage: pdgfra_hold Vhold_mV g_nS T_s seed skUseEK skVoltIndep sfNaLCN sfSOC out.csv
// Prints time averages over t > 5 s: Vm, I into syncytium, ISK, SK conductance, and SK event rate.
#include "PDGFRaCell.h"
#include "CellPart.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
	double Vh = atof(argv[1]), g = atof(argv[2]), T = atof(argv[3]) * 1000.0;
	unsigned seed = (unsigned)atoi(argv[4]);
	bool ek = atoi(argv[5]) != 0, vi = atoi(argv[6]) != 0;
	double sfNa = atof(argv[7]), sfSOC = atof(argv[8]);
	const double dt = 0.02;
	PDGFRaCell c(seed);
	c.Init();
	c.SetModified(ek, vi, sfNa, sfSOC);
	CCellPart& raw = c.Raw();
	FILE* f = fopen(argv[9], "w");
	fprintf(f, "t_ms,Vm,Icoupling_pA,ISK_pA,gSK_nS,Cai_mM\n");
	double sV = 0, sI = 0, sSK = 0, sG = 0; long n = 0; int events = 0; bool inEv = false;
	long step = 0, every = (long)(1.0 / dt + 0.5);
	for (double t = 0; t < T; t += dt, step++)
	{
		double I = g * (c.Vm() - Vh);
		c.Step(dt, I);
		if (step % every == 0)
		{
			double Erev = ek ? mRTF * log(raw.Ke / raw.Ki) : -70.0;
			double drive = c.Vm() - Erev;
			double gSK = fabs(drive) > 1.0 ? c.ISK() / drive : 0.0;
			fprintf(f, "%.1f,%.3f,%.3f,%.3f,%.4f,%.4e\n", t, c.Vm(), I, c.ISK(), gSK, c.Cai());
			if (t > 5000)
			{
				sV += c.Vm(); sI += I; sSK += c.ISK(); sG += gSK; n++;
				if (!inEv && c.Cai() > 5e-4) { inEv = true; events++; }
				if (inEv && c.Cai() < 2e-4) inEv = false;
			}
		}
	}
	fclose(f);
	double dur = (T - 5000) / 60000.0;
	printf("EK %d Vind %d NaLCN %.2f SOC %.2f | %s %4.0f mV: Vm %6.1f  I->syncytium %+6.1f pA  ISK %5.1f pA  gSK %.3f nS  Ca events %.1f /min\n",
		ek, vi, sfNa, sfSOC, g > 0 ? "held" : "free", g > 0 ? Vh : 0.0, sV / n, sI / n, sSK / n, sG / n, events / dur);
	return 0;
}
