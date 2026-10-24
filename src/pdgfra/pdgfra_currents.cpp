// Modified PDGFa v2b cell held near the syncytium potential: records the individual membrane currents to find
// what carries the depolarizing current between spontaneous Ca2+ events.
// Usage: pdgfra_currents Vhold_mV g_nS T_s seed sfNaLCN sfSOC out.csv   (SK: EK, voltage-independent, PMU Ca, Kd 1.5 uM, Hill 4)
// Currents are membrane currents of the cell (outward positive, pA). Icoupling = g*(Vm-Vhold) is the current sent into
// the syncytium (> 0 depolarizes it).
#include "PDGFRaCell.h"
#include "CellPart.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
	double Vh = atof(argv[1]), g = atof(argv[2]), T = atof(argv[3]) * 1000.0;
	unsigned seed = (unsigned)atoi(argv[4]);
	double sfNa = atof(argv[5]), sfSOC = atof(argv[6]);
	const double dt = 0.02;
	PDGFRaCell c(seed);
	c.Init();
	c.SetModified(true, true, sfNa, sfSOC);
	c.SetSKCaSource(1, 1.5e-3, 4.0);
	CCellPart& r = c.Raw();
	FILE* f = fopen(argv[7], "w");
	fprintf(f, "t_ms,Vm,Icoupling,ISK,ISOC,INaLCN,INaK,INaCa,IClb,IK1,IAuto,Itot,CaSR_mM,Cai_mM\n");
	// averages over t > 5 s, split into "between events" (ISK < 50 pA) and "during events"
	const int NC = 9;
	double sq[2][NC] = { { 0 } }; long nq[2] = { 0, 0 };
	long step = 0, every = (long)(1.0 / dt + 0.5);
	for (double t = 0; t < T; t += dt, step++)
	{
		double I = g * (c.Vm() - Vh);
		c.Step(dt, I);
		if (step % every == 0)
		{
			double v[NC] = { I, r.ISK, r.ISOC, r.INaLCN, r.INaK, r.INaCa, r.IClb, r.IK1, r.IAutonomousInward };
			fprintf(f, "%.1f,%.3f", t, c.Vm());
			for (int k = 0; k < NC; k++) fprintf(f, ",%.3f", v[k]);
			fprintf(f, ",%.3f,%.4e,%.4e\n", r.Itot, r.CaSR, r.Cai);
			if (t > 5000)
			{
				int q = r.ISK < 50.0 ? 0 : 1;
				for (int k = 0; k < NC; k++) sq[q][k] += v[k];
				nq[q]++;
			}
		}
	}
	fclose(f);
	const char* name[NC] = { "Icoupling", "ISK", "ISOC", "INaLCN", "INaK", "INaCa", "IClb", "IK1", "IAuto" };
	printf("NaLCN x%.2f SOC x%.2f | time between events %.0f%%\n", sfNa, sfSOC, 100.0 * nq[0] / (nq[0] + nq[1]));
	for (int q = 0; q < 2; q++)
	{
		printf("  %s:", q == 0 ? "between events" : "during events ");
		for (int k = 0; k < NC; k++) printf(" %s %+.1f", name[k], nq[q] ? sq[q][k] / nq[q] : 0.0);
		printf("\n");
	}
	return 0;
}
