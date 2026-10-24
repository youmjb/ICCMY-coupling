"""Low own conductance of ICC-MY (iccRest f 0.2, 0.1) with the active muscle C1 (regenerative inward current g_T 2 nS, vmT -50,
vhT -65, kmT 4; SMC tonic K 1 nS; figA_runs 'act'), nifedipine. Can the muscle's voltage-dependent component, removed by
hyperpolarization, give the recorded fall of the muscle slow wave with pinacidil (x0.22) while the transfer stays high?
Seed 1, 40 s. Usage: python paper_R1/lowg_active.py run | analyze      Runs: out/paper_lowg/act_f<f>_<cond>_gc<gc>.csv"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
from final_model import BASE
OUT = "out/paper_lowg"
F = ["0.2", "0.1"]
ACT = ["gT=2", "vmT=-50", "vhT=-65", "kmT=4", "gK0run=1"]
GC = ["0.75", "1", "1.5", "2", "3"]
PIN = ["1", "2", "4", "8", "16"]
JOBS = [(f, c, gc) for f in F for gc in GC for c in ["ctrl"] + [f"pin{p}" for p in PIN]]


def name(f, c, gc): return f"{OUT}/act_f{f}_{c}_gc{gc}"


def run(j):
    f, c, gc = j; n = name(*j)
    if os.path.exists(n + ".log"): return
    a = BASE + ACT + [f"iccRest={f}", "seed=1"] + ([f"gKATP={c[3:]}"] if c.startswith("pin") else [])
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", n + ".csv"] + a, stdout=lg,
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
            print(f"\n=== f {f}, active muscle C1")
            for gc in GC:
                c0 = S(f, "ctrl", gc)
                if not c0: print(f"  g_c {gc}: control failed"); continue
                s = (f"  g_c {gc:>4}: ICC rest {c0['ir']:6.1f} amp {c0['a']['icc']:5.1f} | muscle rest {c0['sr']:6.1f} amp {c0['a']['smc']:5.1f} "
                     f"| ratio {c0['a']['smc'] / c0['a']['icc']:.2f} || pin:")
                for p in PIN:
                    c1 = S(f, f"pin{p}", gc)
                    if not c1: s += f" [{p}: fail]"; continue
                    dM = c1["sr"] - c0["sr"]
                    s += f" [{p}: dM {dM:+5.1f} T {(c1['ir'] - c0['ir']) / dM:.2f} ICCx{c1['a']['icc'] / c0['a']['icc']:.2f} Mx{c1['a']['smc'] / c0['a']['smc']:.2f}]"
                print(s)
