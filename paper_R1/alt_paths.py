"""Option C: can the drug transfer be reached WITHOUT lowering the own conductance of ICC-MY?
(1) Direct coupling of ICC-MY to FLC: original-conductance model (final_model 'A' BASE: ICC-MY upstroke/clock variant, f = 1, passive
    muscle, nifedipine) with ICC-MY coupled by gIF to the FLC of layer 1, whose density fP0 (FLC per SMC) is raised there.
(2) Scaling check: model C with ICC-MY and column scaled 10-fold (iccScale 7.2, m 3, g_c x10). With the resting currents x0.1 this ICC has
    the ABSOLUTE own conductance of the original 18-pF model but 10-fold active currents and capacitance; if the result equals model C,
    the data constrain the own conductance relative to the cell's active currents, not its absolute value.
Seed 1, 40 s. Usage: python paper_R1/alt_paths.py run | analyze      Runs: out/paper_alt/"""
import os, sys, subprocess, importlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
OUT = "out/paper_alt"


def base(model):
    os.environ["FM_MODEL"] = model
    import final_model; importlib.reload(final_model); return list(final_model.BASE)


BA, BC = base("A"), base("C")
DR = {"ctrl": [], "pin2": ["gKATP=2"], "pin4": ["gKATP=4"], "pin8": ["gKATP=8"], "apa": ["skKd=1000"]}
FLC = [(fp, g, gc) for fp in ("1", "3") for g in ("1", "3", "10") for gc in ("1", "1.75")]
JOBS = [("flc", f"fp{fp}_g{g}", gc, BA + [f"fP0={fp}", f"gIF={g}"], d) for fp, g, gc in FLC for d in DR]
JOBS += [("scale", "x10", "12.5", BC + ["iccScale=7.2", "m=3"], d) for d in ("ctrl", "pin2", "pin3", "pin4", "apa")]
DR["pin3"] = ["gKATP=3"]


def name(kind, tag, gc, d): return f"{OUT}/{kind}_{tag}_gc{gc}_{d}"


def run(j):
    kind, tag, gc, keys, d = j; n = name(kind, tag, gc, d)
    if os.path.exists(n + ".log"): return
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", n + ".csv"] + keys + ["seed=1"] + DR[d], stdout=lg,
                       env=dict(os.environ, SMS_CAL="0"))


def S(*a):
    try:
        x = summary(name(*a) + ".csv"); return x if x and x["a"] else None
    except Exception:
        return None


def line(kind, tag, gc, pins):
    c0 = S(kind, tag, gc, "ctrl")
    if not c0: return f"{kind} {tag} g_c {gc}: control failed"
    s = (f"{kind} {tag} g_c {gc}: ICC rest {c0['ir']:.1f} amp {c0['a']['icc']:.1f} | muscle rest {c0['sr']:.1f} amp {c0['a']['smc']:.1f} | "
         f"ratio {c0['a']['smc'] / c0['a']['icc']:.2f} ||")
    rows = []
    for p in pins:
        c1 = S(kind, tag, gc, p)
        if not c1: s += f" [{p}: fail]"; continue
        dM = c1["sr"] - c0["sr"]; T = (c1["ir"] - c0["ir"]) / dM
        if p == "apa": s += f" [apa: dM {dM:+.1f} T {T:.2f} ICCx{c1['a']['icc'] / c0['a']['icc']:.2f}]"
        else:
            rows.append((dM, T, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
            s += f" [{p}: dM {dM:+.1f} T {T:.2f} ICCx{rows[-1][2]:.2f} Mx{rows[-1][3]:.2f}]"
    rows.sort()
    if rows and rows[0][0] <= -10.3 <= rows[-1][0]:
        v = [np.interp(-10.3, [r[0] for r in rows], [r[i] for r in rows]) for i in (1, 2, 3)]
        s += f" => at -10.3 mV: T {v[0]:.2f} ICCx{v[1]:.2f} Mx{v[2]:.2f}"
    return s


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, JOBS))
    else:
        for fp, g, gc in FLC: print(line("flc", f"fp{fp}_g{g}", gc, ["pin2", "pin4", "pin8", "apa"]))
        print(line("scale", "x10", "12.5", ["pin2", "pin3", "pin4", "apa"]))
        print("reference: original-conductance model g_c 1.75 (Table 2): T 0.40, ICC x1.01; model C g_c 1.23: T 0.79, ICC x1.14, muscle x0.51")
