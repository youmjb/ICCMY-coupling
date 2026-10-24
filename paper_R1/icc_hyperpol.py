"""Isolated ICC-MY model (Youm et al. 2019 equations; gisim 'active' with g_c 0, so the muscle column is not connected):
does steady hyperpolarization enlarge the pacemaker potential, as pinacidil does in situ (Kito, Ward & Sanders 2005: ICC-MY -9.6 mV,
amplitude ~56 -> ~65 mV, dV/dt up, frequency unchanged)?
Two ways: outward current iccInj (pA, whole run) and an added K+ conductance iccGL (nS) with iccEL -90 mV (K_ATP-like).
Cell size: iccScale 1 (2019 cell) and 0.72 (18 pF, as in the paper).
Usage: python paper_R1/icc_hyperpol.py run | analyze       Runs: out/paper_icc_hyp/<scale>_<mode><value>.csv (40 s)"""
import os, sys, csv, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_icc_hyp"
SC = ["1", "0.72"]
INJ = ["0", "1", "2", "3", "5", "8", "15", "25", "35", "50"]
GK = ["0.5", "1", "2", "4"]
JOBS = [(s, "i", v) for s in SC for v in INJ] + [(s, "g", v) for s in SC for v in GK]


def name(s, mode, v):
    return f"{OUT}/s{s}_{mode}{v}"


def run(j):
    s, mode, v = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    k = [f"iccScale={s}", "L=2", "m=0.3", "fI=0", "fP=0", "smcModel=1", "seed=1"]
    k += ([f"iccInj={v}", "injOn=0", "injOff=1e9"] if mode == "i" else [f"iccGL={v}", "iccEL=-90"])
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", "0", "0", "0", "40000", "0.1", f + ".csv"] + k, check=True, stdout=lg)


def stats(f):
    r = list(csv.DictReader(open(f)))
    t, v = (np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC"))
    h = t >= 20000; t, v = t[h], v[h]
    rest, peak = np.percentile(v, 2), np.percentile(v, 99.5)
    mid = rest + 0.5 * (peak - rest)
    up = np.where((v[:-1] < mid) & (v[1:] >= mid))[0]
    amp, dv = [], []
    for a, b in zip(up[:-1], up[1:]):
        seg = v[a:b]; base = v[max(0, a - 200):max(1, a - 50)].mean()
        amp.append(seg.max() - base); dv.append(np.diff(v[max(0, a - 100):a + 100]).max() / (t[1] - t[0]))
    cpm = 60000.0 / np.mean(np.diff(t[up])) if len(up) > 1 else 0.0
    return rest, (np.mean(amp) if amp else 0.0), peak, (np.mean(dv) if dv else 0.0), cpm, len(up)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(20) as ex:
            list(ex.map(run, JOBS))
    else:
        for s in SC:
            r0 = stats(name(s, "i", "0") + ".csv")
            print(f"iccScale {s}: control rest {r0[0]:.1f} mV, amplitude {r0[1]:.1f} mV, peak {r0[2]:.1f}, dV/dt {r0[3]:.2f} V/s, {r0[4]:.1f} min-1")
            for mode, vals, lab in (("i", INJ[1:], "pA"), ("g", GK, "nS K")):
                for v in vals:
                    r = stats(name(s, mode, v) + ".csv")
                    print(f"   {v:>3} {lab:5s} rest {r[0]:6.1f} ({r[0] - r0[0]:+5.1f}) | amplitude {r[1]:5.1f} (x{r[1] / r0[1]:.2f}) | "
                          f"peak {r[2]:6.1f} | dV/dt {r[3]:5.2f} (x{r[3] / r0[3]:.2f}) | {r[4]:5.1f} min-1 (n {r[5]})")
        print("in situ (Kito 2005, pinacidil): ICC-MY -9.6 mV, amplitude x~1.16, dV/dt up, frequency unchanged")
