"""Which resting current makes the published ICC model leakier than freshly prepared ICC? Scale only the inward-rectifier K+ current
(iccK1) or all resting currents and their transporters (iccRest), and compare input resistance, resting potential and pacemaker output.
(1) Isolated cell: published model (25 pF, no changes) and the cell of this study (18 pF, upstroke/clock variant); chord conductance from a
    5-pA outward current (60 s, last 30 s). Goto et al. 2004: 25.2 pF, 0.76 +/- 0.41 GOhm, -72.5 mV.
(2) Coupled: model D base with iccRest=1 and iccK1 = x (active muscle, nifedipine), g_c 1.25 nS, pinacidil 2/4/6 nS, seed 1.
Usage: python paper_R1/k1_only.py run | analyze      Runs: out/paper_k1/"""
import os, sys, subprocess, importlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from icc_hyperpol import stats
from active_cal import summary
OUT = "out/paper_k1"
os.environ["FM_MODEL"] = "D"
import final_model as FM
CELL = {"pub": ["iccScale=1"], "var": ["iccScale=0.72", "iccVDDR=16", "iccVH=-80", "ip3VK=0"]}
SCALE = [("rest", "1"), ("K1", "0.5"), ("K1", "0.3"), ("K1", "0.2"), ("K1", "0.1"), ("rest", "0.3")]
ISO = [(c, k, v, i) for c in CELL for k, v in SCALE for i in ("0", "5")]
COL = [(k, v, d) for k, v in (("K1", "0.5"), ("K1", "0.3"), ("K1", "0.2")) for d in ("ctrl", "pin2", "pin4", "pin6")]


def iso_name(c, k, v, i): return f"{OUT}/iso_{c}_{k}{v}_i{i}"
def col_name(k, v, d): return f"{OUT}/col_{k}{v}_{d}"


def key(k, v): return [f"iccRest={v}"] if k == "rest" else ["iccRest=1", f"iccK1={v}"]


def run(job):
    if job[0] == "iso":
        _, c, k, v, i = job; n = iso_name(c, k, v, i)
        if os.path.exists(n + ".log"): return
        a = ["./bin/gisim.exe", "active", "0", "0", "0", "60000", "0.1", n + ".csv"] + CELL[c] + key(k, v) + \
            ["L=2", "m=0.3", "fI=0", "fP=0", "smcModel=1", "seed=1", f"iccInj={i}", "injOn=0", "injOff=1e9"]
        env = os.environ
    else:
        _, k, v, d = job; n = col_name(k, v, d)
        if os.path.exists(n + ".log"): return
        a = ["./bin/gisim.exe", "active", "1.25", "0", "0", "40000", "0.1", n + ".csv"] + FM.BASE + key(k, v) + ["seed=1"] + \
            ([f"gKATP={d[3:]}"] if d.startswith("pin") else [])
        env = dict(os.environ, SMS_CAL="0")
    with open(n + ".log", "w") as lg: subprocess.run(a, stdout=lg, env=env)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, [("iso",) + j for j in ISO] + [("col",) + j for j in COL]))
    else:
        print("isolated cell: scale | rest (mV) amp (mV) rate (min-1) | own conductance (nS), R_in (MOhm), tau (ms), per pF (nS/pF)")
        for c, C in (("pub", 25.0), ("var", 18.0)):
            for k, v in SCALE:
                a, b = stats(iso_name(c, k, v, "0") + ".csv"), stats(iso_name(c, k, v, "5") + ".csv")
                paced = a[0] < -50
                if not paced: print(f"  {c} {k} x{v}: rest {a[0]:.1f} (no stable rest)"); continue
                g = 5.0 / (a[0] - b[0])
                print(f"  {c} {k} x{v}: {a[0]:6.1f} {a[1]:5.1f} {a[4]:5.1f} | {g:.2f} nS, {1000 / g:.0f} MOhm, {C / g:.1f} ms, {g / C:.3f} nS/pF")
        print("Goto 2004: 25.2 pF, 760 MOhm (tau 19 ms, 0.052 nS/pF), rest -72.5 mV")
        print("\ncoupled (model D base, g_c 1.25, seed 1): pinacidil interpolated to -10.3 mV")
        for k, v in (("K1", "0.5"), ("K1", "0.3"), ("K1", "0.2")):
            c0 = summary(col_name(k, v, "ctrl") + ".csv"); rows = []
            for d in ("pin2", "pin4", "pin6"):
                c1 = summary(col_name(k, v, d) + ".csv")
                if c1 and c1["a"]:
                    dM = c1["sr"] - c0["sr"]; rows.append((dM, (c1["ir"] - c0["ir"]) / dM, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
            rows.sort(); x = [r[0] for r in rows]
            s = (f"  {k} x{v}: ICC rest {c0['ir']:.1f} amp {c0['a']['icc']:.1f} | muscle rest {c0['sr']:.1f} amp {c0['a']['smc']:.1f} | "
                 f"ratio {c0['a']['smc'] / c0['a']['icc']:.2f}")
            if rows and x[0] <= -10.3 <= x[-1]:
                t = [np.interp(-10.3, x, [r[i] for r in rows]) for i in (1, 2, 3)]
                s += f" | at -10.3 mV: T {t[0]:.2f} ICCx{t[1]:.2f} Mx{t[2]:.2f}"
            else: s += f" | dM range {x[0]:.1f} to {x[-1]:.1f}"
            print(s)
        print("model D (all resting x0.3) at g_c 1.19: T 0.55 ICCx1.06 Mx0.40")
