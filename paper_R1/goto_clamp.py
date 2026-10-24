"""Input resistance of the isolated model ICC measured as Goto, Matsuoka & Noma (2004) did: voltage clamp at a holding potential of -80 mV,
test pulses of -20, +10 and +20 mV; resting input resistance from the steady current change (no time-dependent current), slight inward
rectification expected. Recorded: 0.76 +/- 0.41 GOhm (n 3), 25.2 pF.
Model: single ICC ('cluster' mode, N 1), electrode series resistance 2 MOhm, held at -80 mV from t = 0 after 20 s of free running.
A baseline run (no pulse) gives the holding current; the first quiet window after 0.5 s (no spontaneous inward current, holding current within
+/-2 pA of its median for 600 ms) receives a 300-ms pulse in each test run. dI = I_test - I_baseline averaged over 200-300 ms of the pulse.
Cells: published model (25 pF, unmodified) and the cell of this study (18 pF, upstroke/clock variant), resting currents x f.
Usage: python paper_R1/goto_clamp.py run | analyze      Runs: out/paper_goto/"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_goto"
CELL = {"pub": (["s=1"], 25.0), "var": (["s=0.72", "iccVDDR=16", "iccVH=-80", "ip3VK=0"], 18.0)}
F = ["1", "0.5", "0.3", "0.2", "0.1"]
DV = ["-20", "10", "20"]
T, DT = "8000", "0.005"
COMMON = ["N=1", "vc=1", "Rs=2", "Vh=-80", "smp=0.1"]


def name(c, f, tag): return f"{OUT}/{c}_f{f}_{tag}.csv"


def load(p):
    r = np.loadtxt(p, delimiter=",", skiprows=1); return r[:, 0], r[:, 1], r[:, 3]


def quiet_window(c, f):
    t, v, i = load(name(c, f, "base"))
    for t0 in np.arange(800, 7300, 50):                           # the variant clock runs at -80 mV: first quiet window after ~0.5 s
        w = (t >= t0 - 300) & (t < t0 + 600)
        if np.ptp(i[w]) < 4.0: return t0
    return None


def run_base(job):
    c, f = job; out = name(c, f, "base")
    if not os.path.exists(out + ".log"):
        with open(out + ".log", "w") as lg:
            subprocess.run(["./bin/gisim.exe", "cluster", T, DT, out] + CELL[c][0] + COMMON + [f"iccRest={f}", "dV=0"], stdout=lg)


def run_test(job):
    c, f, dv = job; out = name(c, f, f"dv{dv}")
    if os.path.exists(out + ".log"): return
    t0 = quiet_window(c, f)
    if t0 is None: return
    with open(out + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "cluster", T, DT, out] + CELL[c][0] + COMMON +
                       [f"iccRest={f}", f"dV={dv}", f"tOn={t0}", "tDur=300"], stdout=lg)
    open(out + ".t0", "w").write(f"{t0}\n")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        base = [(c, f) for c in CELL for f in F]
        with ThreadPoolExecutor(20) as ex: list(ex.map(run_base, base))
        with ThreadPoolExecutor(30) as ex: list(ex.map(run_test, [(c, f, d) for c, f in base for d in DV]))
    else:
        print("cell f | holding current at -80 mV (pA) | R_in from -20 / +10 / +20 mV (MOhm) | conductance per pF (-20 mV) | tau = C R (ms)")
        for c in CELL:
            C = CELL[c][1]
            for f in F:
                try:
                    tb, vb, ib = load(name(c, f, "base"))
                    t0 = float(open(name(c, f, "dv-20") + ".t0").read())
                except Exception as e:
                    print(f"  {c} f {f}: no quiet window or failed ({e})"); continue
                hold = np.median(ib[(tb > t0 - 300) & (tb < t0)])
                Rs = []
                for dv in DV:
                    t, v, i = load(name(c, f, f"dv{dv}"))
                    w = (t >= t0 + 200) & (t < t0 + 300)
                    dI = np.mean(i[w] - np.interp(t[w], tb, ib))
                    dVm = np.mean(v[w] - np.interp(t[w], tb, vb))           # membrane step after series-resistance drop
                    Rs.append(dVm / dI * 1000.0)                           # mV / pA -> GOhm *1000 -> MOhm
                g = 1000.0 / Rs[0]
                print(f"  {c} f {f:>3}: hold {hold:7.1f} pA | {Rs[0]:6.0f} / {Rs[1]:6.0f} / {Rs[2]:6.0f} MOhm | {g / C:.3f} nS/pF | {C * Rs[0] / 1000:.1f} ms")
        print("Goto 2004: 0.76 +/- 0.41 GOhm (n 3; 0.7 GOhm in the example cell, slightly larger current for -20 than for +20 mV), 25.2 pF -> 0.052 nS/pF, ~19 ms")
