"""Contribution of the resting-conductance factor f alone: the original 2019 ICC model (VDDR x1, midpoint -66 mV, voltage-dependent IP3
production; 18 pF) in the passive calibrated column (tonic K 0.4 nS, nifedipine), g_c 1.75 nS, resting currents x f; pinacidil (K_ATP,
E -84 mV) interpolated to the recorded muscle change of -10.3 mV. Compare with the same column and the upstroke/clock variant
(lowg_coupled.py: transfer 0.40, 0.77, 0.97 for f 1, 0.3, 0.1). Seed 1, 40 s.
Usage: python paper_R1/f_only.py run | analyze      Runs: out/paper_fonly/"""
import os, sys, subprocess, importlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
os.environ["FM_MODEL"] = "A"
import final_model as FM
OUT = "out/paper_fonly"
BASE = FM.BASE + ["iccVDDR=1", "iccVH=-66", "ip3VK=1"]          # undo the upstroke/clock variant of the A base
F = ["1", "0.3", "0.1"]
PIN = ["1", "2", "4", "8"]
JOBS = [(f, c) for f in F for c in ["ctrl"] + [f"pin{p}" for p in PIN]]


def name(f, c): return f"{OUT}/f{f}_{c}"


def run(j):
    f, c = j; n = name(f, c)
    if os.path.exists(n + ".log"): return
    a = ["./bin/gisim.exe", "active", "1.75", "0", "0", "40000", "0.1", n + ".csv"] + BASE + [f"iccRest={f}", "seed=1"] + \
        ([f"gKATP={c[3:]}"] if c.startswith("pin") else [])
    with open(n + ".log", "w") as lg: subprocess.run(a, stdout=lg, env=dict(os.environ, SMS_CAL="0"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(15) as ex: list(ex.map(run, JOBS))
    else:
        for f in F:
            c0 = summary(name(f, "ctrl") + ".csv"); rows = []
            for p in PIN:
                c1 = summary(name(f, f"pin{p}") + ".csv")
                if c1 and c1["a"]:
                    dM = c1["sr"] - c0["sr"]; rows.append((dM, (c1["ir"] - c0["ir"]) / dM, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
            rows.sort(); x = [r[0] for r in rows]
            s = (f"original ICC, f {f}: ICC rest {c0['ir']:.1f} amp {c0['a']['icc']:.1f} | muscle rest {c0['sr']:.1f} amp {c0['a']['smc']:.1f} | "
                 f"ratio {c0['a']['smc'] / c0['a']['icc']:.2f}")
            if rows and x[0] <= -10.3 <= x[-1]:
                t = [np.interp(-10.3, x, [r[i] for r in rows]) for i in (1, 2, 3)]
                s += f" | at -10.3 mV: transfer {t[0]:.2f}, ICC x{t[1]:.2f}, muscle x{t[2]:.2f}"
            else: s += f" | dM range {x[0]:.1f} to {x[-1]:.1f}"
            print(s)
        print("variant upstroke/clock, same column (lowg_coupled): transfer 0.40 / 0.77 / 0.97, ICC x1.01 / x1.13 / x1.20 for f 1 / 0.3 / 0.1")
