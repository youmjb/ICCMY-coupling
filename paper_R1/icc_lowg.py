"""Isolated ICC-MY variant (18 pF; VDDR x16, inactivation midpoint -80 mV, IP3 production voltage-independent) with its resting
conductances and steady transporters scaled together by f (key iccRest: NaLCN, K1, Clb, Na/K pump, NKCC). Own conductance from the
rest shift with a small outward current (chord, 5 pA). Question: can the own conductance fall to ~0.1-0.3 nS and the cell still pace?
Usage: python paper_R1/icc_lowg.py run | analyze      Runs: out/paper_lowg/iso_f<f>_i<pA>.csv (60 s; last 30 s analysed)"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from icc_hyperpol import stats
OUT = "out/paper_lowg"
VAR = ["iccScale=0.72", "iccVDDR=16", "iccVH=-80", "ip3VK=0"]
F = ["1", "0.5", "0.3", "0.2", "0.1", "0.05", "0.03"]
INJ = ["0", "5", "35"]
JOBS = [(f, i) for f in F for i in INJ]


def name(f, i): return f"{OUT}/iso_f{f}_i{i}"


def run(j):
    f, i = j; n = name(f, i)
    if os.path.exists(n + ".log"): return
    k = VAR + [f"iccRest={f}", "L=2", "m=0.3", "fI=0", "fP=0", "smcModel=1", "seed=1", f"iccInj={i}", "injOn=0", "injOff=1e9"]
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", "0", "0", "0", "60000", "0.1", n + ".csv"] + k, stdout=lg)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(21) as ex: list(ex.map(run, JOBS))
    else:
        for f in F:
            try:
                r = {i: stats(name(f, i) + ".csv") for i in INJ}
            except Exception as e:
                print(f"f {f}: failed ({e})"); continue
            r0, r5, r35 = r["0"], r["5"], r["35"]
            g = 5.0 / (r0[0] - r5[0]) if r0[0] != r5[0] else float("nan")
            print(f"f {f:>4}: rest {r0[0]:6.1f} mV, amp {r0[1]:5.1f}, dV/dt {r0[3]:5.2f}, rate {r0[4]:4.1f} min-1 (n {r0[5]}) | "
                  f"g_own(5 pA) {g:5.2f} nS | 35 pA: dV {r35[0] - r0[0]:+5.1f}, amp x{r35[1] / max(r0[1], 1e-9):.2f}, "
                  f"dV/dt x{r35[3] / max(r0[3], 1e-9):.2f}, rate {r35[4]:4.1f}")
