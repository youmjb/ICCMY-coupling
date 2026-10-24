"""Calibration of the active-SMC variant (T-type-like current, gisim keys gT vmT kmT vhT) to the mouse control state only:
CSMC resting potential about -70 mV (Kito & Suzuki 2007 Table 1: -69.1 to -72.7), muscle/ICC-MY amplitude ratio about 0.60
(Kito & Suzuki 2003: 34.1/56.5; 2007: 37/61) and CSMC dV/dt max below 1 V/s (0.85-0.92). The drug responses are NOT used here;
they are predicted afterwards (active_drugs.py).
Usage: python paper_R1/active_cal.py run | analyze       Runs: out/paper_active/cal_*.csv (20 s, calibrated column m 0.3, seed 1)"""
import os, sys, csv, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
OUT = "out/paper_active"
B = "m=0.3 fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1 kmT=4".split()
GRID = dict(gc=["0.5", "1", "1.5"], gT=["2", "4", "8"], vmT=["-50", "-45"], vhT=["-65", "-60"], gK=["0.4", "1"])
JOBS = list(itertools.product(*GRID.values()))


def name(gc, gT, vmT, vhT, gK):
    return f"{OUT}/cal_gc{gc}_T{gT}_m{vmT}_h{vhT}_k{gK}"


def run(j):
    gc, gT, vmT, vhT, gK = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + B
                       + [f"gT={gT}", f"vmT={vmT}", f"vhT={vhT}", f"gK0run={gK}"], check=True, stdout=lg)


def summary(f):
    """CSMC rest, ICC-MY rest, event-aligned amplitudes, CSMC dV/dt max (V/s) over the second half"""
    import amp_event as AE
    r = list(csv.DictReader(open(f)))
    t, vi, vs = (np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0"))
    h = t >= t[-1] / 2; t, vi, vs = t[h], vi[h], vs[h]
    near = np.convolve((vi >= -65).astype(float), np.ones(601), mode="same") > 0
    sr = np.median(vs[~near]) if (~near).any() else np.percentile(vs, 5)
    a = AE.measure(f)
    dv = np.diff(vs).max() / (t[1] - t[0])
    return dict(sr=sr, ir=np.percentile(vi, 5), a=a, dv=dv)


def score(s):
    if s["a"] is None: return 1e9
    q = s["a"]["smc"] / s["a"]["icc"]
    return abs(s["sr"] + 70) / 2 + abs(q - 0.60) / 0.05 + max(0, s["dv"] - 1.0) / 0.2


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(24) as ex:
            list(ex.map(run, JOBS))
    else:
        rows = []
        for j in JOBS:
            f = name(*j) + ".csv"
            if not os.path.exists(f): continue
            s = summary(f); rows.append((score(s), j, s))
        rows.sort(key=lambda r: r[0])
        print("score | g_c gT vmT vhT gK0 | CSMC rest | ICC rest | ICC amp | CSMC amp | ratio | CSMC dV/dt (V/s)")
        for sc, j, s in rows[:25]:
            a = s["a"] or {"icc": float("nan"), "smc": float("nan")}
            print(f"{sc:6.2f} | {' '.join(j)} | {s['sr']:6.1f} | {s['ir']:6.1f} | {a['icc']:5.1f} | {a['smc']:5.1f} | "
                  f"{a['smc'] / a['icc']:.2f} | {s['dv']:.2f}")
        print("targets: CSMC rest about -70 mV, ratio 0.60, dV/dt < 1 V/s")
