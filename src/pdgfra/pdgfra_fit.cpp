// Calibration runs for the modified PDGFa (FLC / PDGFRa+) cell, in an "in situ-like" setting:
// the cell is coupled through g_nS to a circular SMC held at Vhold (Kito et al. 2014: CSMC -63.5 mV, FLC -49.5 mV).
// Usage: pdgfra_fit Vhold_mV g_nS T_s seed sfNaLCN skKd_uM skScale out.csv
//   skScale multiplies the SK conductance (0 = apamin). Writes a 1 ms trace; the analysis is done in Python.
#include "PDGFRaCell.h"
#include "CellPart.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
	double Vh = atof(argv[1]), g = atof(argv[2]), T = atof(argv[3]) * 1000.0;
	unsigned seed = (unsigned)atoi(argv[4]);
	double sfNa = atof(argv[5]), kd = atof(argv[6]), sk = atof(argv[7]);
	const double dt = 0.02;
	PDGFRaCell c(seed);
	c.Init();
	c.SetModified(true, true, sfNa, 1.0);	// SK = EK, voltage-independent, deterministic gate; gap current carried by K+
	c.SetSKCaSource(1, kd * 1e-3, 4.0);		// SK senses submembrane (IP3R release) Ca2+
	c.SetSKScale(sk);
	FILE* f = fopen(argv[8], "w");
	fprintf(f, "t_ms,Vm,ISK,Cai_mM\n");
	long step = 0, every = (long)(1.0 / dt + 0.5);
	for (double t = 0; t < T; t += dt, step++)
	{
		c.Step(dt, g * (c.Vm() - Vh));
		if (step % every == 0) fprintf(f, "%.1f,%.3f,%.3f,%.4e\n", t, c.Vm(), c.ISK(), c.Cai());
	}
	fclose(f);
	return 0;
}
