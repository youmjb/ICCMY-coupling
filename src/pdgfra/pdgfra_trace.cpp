// Isolated original PDGFa-project cell (stochastic SK/IP3R), Iext = 0: writes a 1 ms-sampled trace.
// Usage: pdgfra_trace T_s dt_ms seed out.csv
#include "PDGFRaCell.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
	double T = argc > 1 ? atof(argv[1]) * 1000.0 : 60000.0;
	double dt = argc > 2 ? atof(argv[2]) : 0.02;	// the GUI default; SK current is dt-dependent
	unsigned seed = argc > 3 ? (unsigned)atoi(argv[3]) : 1;
	const char* out = argc > 4 ? argv[4] : "pdgfra_trace.csv";
	PDGFRaCell c(seed);
	c.Init();
	FILE* f = fopen(out, "w");
	fprintf(f, "t_ms,Vm,Cai_mM,ISK_pA\n");
	long step = 0, every = (long)(1.0 / dt + 0.5);
	for (double t = 0; t < T; t += dt, step++)
	{
		c.Step(dt, 0.0);
		if (step % every == 0) fprintf(f, "%.1f,%.3f,%.4e,%.3f\n", t, c.Vm(), c.Cai(), c.ISK());
	}
	fclose(f);
	return 0;
}
