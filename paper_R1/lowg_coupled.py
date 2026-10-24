"""Reference column (final_model.BASE, nifedipine) with the ICC-MY resting conductances and steady transporters scaled by f (iccRest).
Control amplitude ratio against g_c, and pinacidil (SMC K_ATP, E -84 mV) at each g_c: does a low own conductance of ICC-MY give the
recorded transfer (0.93) with a larger ICC-MY amplitude (x1.17) at the control ratio of the pinacidil set (0.52)? Seed 1, 40 s.
Usage: python paper_R1/lowg_coupled.py run | analyze      Runs: out/paper_lowg/col_f<f>_<cond>_gc<gc>.csv"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
from final_model import BASE
OUT = "out/paper_lowg"
F = ["1", "0.3", "0.2", "0.1", "0.05"]
GC = ["0.5", "0.75", "1", "1.25", "1.5", "1.75", "2", "2.5", "3"]
PIN = ["2", "4", "8"]
JOBS = [(f, c, gc) for f in F for gc in GC for c in ["ctrl"] + [f"pin{p}" for p in PIN]]
PIN_LOW = ["0.5", "1", "1.5"]                                       # sizing to the recorded -10.3 mV at g_c 1.75
JOBS += [(f, f"pin{p}", "1.75") for f in F for p in PIN_LOW]


def name(f, c, gc): return f"{OUT}/col_f{f}_{c}_gc{gc}"


def run(j):
    f, c, gc = j; n = name(*j)
    if os.path.exists(n + ".log"): return
    k = BASE + [f"iccRest={f}", "seed=1"] + ([f"gKATP={c[3:]}"] if c.startswith("pin") else [])
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", n + ".csv"] + k, stdout=lg,
                       env=dict(os.environ, SMS_CAL="0"))


def S(*j):
    try:
        x = summary(name(*j) + ".csv"); return x if x and x["a"] else None
    except Exception:
        return None


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, JOBS))
    else:
        for f in F:
            print(f"\n=== iccRest f {f}")
            for gc in GC:
                c0 = S(f, "ctrl", gc)
                if not c0: print(f"  g_c {gc:>4}: control failed"); continue
                line = (f"  g_c {gc:>4}: ICC rest {c0['ir']:6.1f} amp {c0['a']['icc']:5.1f} | muscle rest {c0['sr']:6.1f} "
                        f"amp {c0['a']['smc']:5.1f} | ratio {c0['a']['smc'] / c0['a']['icc']:.2f} || pin:")
                for p in PIN:
                    c1 = S(f, f"pin{p}", gc)
                    if not c1: line += f"  [{p}: fail]"; continue
                    dS, dI = c1["sr"] - c0["sr"], c1["ir"] - c0["ir"]
                    line += (f"  [{p}: dM {dS:+5.1f} T {dI / dS:.2f} ICCx{c1['a']['icc'] / c0['a']['icc']:.2f} "
                             f"Mx{c1['a']['smc'] / c0['a']['smc']:.2f}]")
                print(line)
