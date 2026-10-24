"""Low own conductance of ICC-MY (iccRest f 0.2, 0.1) with the muscle tonic K+ conductance raised (gK0run) so that both cells rest
near the recorded -70 mV. If the muscle rests closer to E_K, the recorded 10-mV hyperpolarization needs more K_ATP, which may shunt
the muscle slow wave (recorded x0.22). Control ratio against g_c and pinacidil at each g_c; seed 1, 40 s.
Usage: python paper_R1/lowg_rest.py run | analyze      Runs: out/paper_lowg/rest_f<f>_k<gK0>_<cond>_gc<gc>.csv"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
from final_model import BASE
OUT = "out/paper_lowg"
F = ["0.2", "0.1"]
K = ["0.8", "1.2", "1.6"]
GC = ["2", "3", "4"]
PIN = ["1", "2", "4", "8"]
JOBS = [(f, k, c, gc) for f in F for k in K for gc in GC for c in ["ctrl"] + [f"pin{p}" for p in PIN]]


def name(f, k, c, gc): return f"{OUT}/rest_f{f}_k{k}_{c}_gc{gc}"


def run(j):
    f, k, c, gc = j; n = name(*j)
    if os.path.exists(n + ".log"): return
    a = BASE + [f"iccRest={f}", f"gK0run={k}", "seed=1"] + ([f"gKATP={c[3:]}"] if c.startswith("pin") else [])
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", n + ".csv"] + a, stdout=lg,
                       env=dict(os.environ, SMS_CAL="0"))


def S(*j):
    try:
        x = summary(name(*j) + ".csv"); return x if x and x["a"] else None
    except Exception:
        return None


def sized(f, k, gc, target=-10.3):
    c0 = S(f, k, "ctrl", gc); rows = []
    for p in PIN:
        c1 = S(f, k, f"pin{p}", gc)
        if c0 and c1:
            dM = c1["sr"] - c0["sr"]
            rows.append((dM, (c1["ir"] - c0["ir"]) / dM, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
    rows.sort()
    if not rows: return None
    x = [r[0] for r in rows]
    return [float(np.interp(target, x, [r[i] for r in rows])) for i in (1, 2, 3)] + [min(x), max(x)]


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, JOBS))
    else:
        for f in F:
            for k in K:
                print(f"\n=== f {f}, muscle tonic K {k} nS")
                for gc in GC:
                    c0 = S(f, k, "ctrl", gc)
                    if not c0: print(f"  g_c {gc}: control failed"); continue
                    z = sized(f, k, gc)
                    s = (f"  g_c {gc}: ICC rest {c0['ir']:6.1f} amp {c0['a']['icc']:5.1f} | muscle rest {c0['sr']:6.1f} amp {c0['a']['smc']:5.1f} "
                         f"| ratio {c0['a']['smc'] / c0['a']['icc']:.2f}")
                    if z: s += f" || pin at -10.3: T {z[0]:.2f} ICCx{z[1]:.2f} Mx{z[2]:.2f} (dM {z[3]:.1f} to {z[4]:.1f})"
                    print(s)
