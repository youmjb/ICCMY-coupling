"""Model vs Kito & Suzuki 2007 (mouse ileum, nifedipine): apamin (SK block) depolarized CSMC by 3.1 mV and ICC-MY by 3.3 mV
(Table 1; unpaired tissues). Existing runs: control out/v4win/m0.3_gSS30_gc{gc}_u0, SK block out/paper_cal/apa_gc{gc} (skKd 1000).
Resting potential = 5th percentile of V over the second half (ICC-MY) and median between slow waves (CSMC), as in analyze_paper_cal.py;
amplitudes event-aligned (amp_event.py)."""
import os, sys, csv
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import amp_event as AE


def load(f):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0")]


def rest(f):
    t, vi, vs = load(f)
    h = t >= 20000; t, vi, vs = t[h], vi[h], vs[h]
    near = np.convolve((vi >= -65).astype(float), np.ones(601), mode="same") > 0
    return np.percentile(vi, 5), np.median(vs[~near]), np.median(vi[~near])


print("g_c | ICC-MY rest ctrl -> SK block (d) | CSMC rest ctrl -> SK block (d) | transfer dICC/dCSMC | amp ICC, CSMC ctrl -> SK block")
for gc in ("1", "1.5", "2"):
    c, a = f"out/v4win/m0.3_gSS30_gc{gc}_u0", f"out/paper_cal/apa_gc{gc}"
    ic, sc, ic2 = rest(c + ".csv"); ia, sa, ia2 = rest(a + ".csv")
    mc, ma = AE.measure(c + ".csv"), AE.measure(a + ".csv")
    print(f"{gc:>3} | {ic:6.1f} -> {ia:6.1f} ({ia - ic:+.1f}); between waves {ic2:6.1f} -> {ia2:6.1f} ({ia2 - ic2:+.1f}) | "
          f"{sc:6.1f} -> {sa:6.1f} ({sa - sc:+.1f}) | {(ia - ic) / (sa - sc):.2f} | {mc['icc']:.1f}, {mc['smc']:.1f} -> {ma['icc']:.1f}, {ma['smc']:.1f}")
print("Kito & Suzuki 2007 Table 1: CSMC -70.8 -> -67.7 (+3.1, n 15), primary 36.2 -> 36.0; ICC-MY -71.7 -> -68.4 (+3.3, n 7), amplitude 61.3 -> 58.8")
