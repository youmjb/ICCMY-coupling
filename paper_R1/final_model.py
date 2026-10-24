"""Final reference model for plan A (RESTRUCTURE_PLAN_A.md): the paper's passive calibrated column (m 0.3, g_SS 30, L 30, SMC tonic K
0.4 nS, FLC as in the paper) in nifedipine (SMC L-type blocked, env SMS_CAL=0), with the ICC-MY variant (VDDR x16, inactivation
midpoint -80 mV, IP3 production voltage-independent; icc_variant.py). Calibrated to the mouse controls only; drugs sized by the recorded
CSMC resting change only; everything else is a prediction. 40 s runs, seeds 1-3. Data: kito_data.md.
Data sets: 'K05' Kito, Ward & Sanders 2005 (control ratio 29.3/56.1 = 0.52; pinacidil), 'K07a' Kito & Suzuki 2007 apamin set
(36.2/61.3 = 0.59; apamin, phenylephrine), 'K07b' 4-AP set (37.0/60.0 = 0.62; 4-AP).
Drugs: pinacidil = SMC K_ATP (E -84 mV); apamin = FLC SK block (skKd 1000, not sized); 4-AP = SMC tonic K reduced;
phenylephrine = small SMC K conductance (E -84 mV), as a stand-in for its SK-mediated hyperpolarization.
Usage: python paper_R1/final_model.py run | analyze      Runs: out/paper_final/<cond>_gc<gc>_s<seed>.csv"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
OUT = "out/paper_final"
BASE = ("m=0.3 gSS=30 L=30 fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 gK0run=0.4 "
        "iccVDDR=16 iccVH=-80 ip3VK=0 EKATP=-84").split()
SEEDS = ["1", "2", "3"]
GC_CTRL = ["1.5", "1.75", "2", "2.25", "2.5", "2.75", "3"]
GC_DRUG = ["1.75", "2.25", "2.75"]
DRUGS = {f"pin{g}": [f"gKATP={g}"] for g in ("1", "1.5", "2", "3", "4")}
DRUGS["apa"] = ["skKd=1000"]
DRUGS.update({f"4ap{x}": [f"gK0run={0.4 * float(x):.3f}"] for x in ("0.75", "0.5", "0.25", "0")})
DRUGS.update({f"pe{g}": [f"gKATP={g}"] for g in ("0.1", "0.2", "0.4")})
MODEL = os.environ.get("FM_MODEL", "A")      # A: reference model of draft A1; C: ICC-MY own conductance x0.1 + active muscle C1
if MODEL in ("C", "D"):                       # lowg_active.py: iccRest 0.1 (C) or 0.3 (D: own conductance per pF of freshly prepared ICC,
    OUT = "out/paper_final" + MODEL            # Goto et al. 2004), SMC regenerative current g_T 2 nS, tonic K 1 nS
    BASE = BASE + ["gT=2", "vmT=-50", "vhT=-65", "kmT=4", "gK0run=1", "iccRest=" + ("0.1" if MODEL == "C" else "0.3")]
    GC_CTRL = ["1", "1.25", "1.5", "1.75", "2"]
    GC_DRUG = ["1.25", "1.5", "1.75"]
    DRUGS = {f"pin{g}": [f"gKATP={g}"] for g in ("2", "3", "4", "6")}
    DRUGS["apa"] = ["skKd=1000"]
    DRUGS.update({f"4ap{x}": [f"gK0run={1.0 * float(x):.3f}"] for x in ("0.75", "0.5", "0.25", "0")})
    DRUGS.update({f"pe{g}": [f"gKATP={g}"] for g in ("0.2", "0.4", "0.8")})
JOBS = [("ctrl", gc, s) for gc in GC_CTRL for s in SEEDS] + [(d, gc, s) for d in DRUGS for gc in GC_DRUG for s in SEEDS]
DATA = {  # control ratio, drug, recorded CSMC change, recorded: dICC, ICC x, CSMC x, ICC dV/dt x, CSMC dV/dt x
    "K05 pinacidil": (0.52, "pin", -10.3, -9.6, 1.17, 0.22, 1.52, 0.25),
    "K07a apamin": (0.59, "apa", +3.1, +3.3, 0.96, 1.00, 0.89, 1.05),
    "K07a phenylephrine": (0.59, "pe", -2.4, -1.6, 1.04, 1.07, 1.01, 0.99),
    "K07b 4-AP": (0.62, "4ap", +3.1, +3.9, 0.94, 1.22, 0.83, 1.31),
}


def name(c, gc, s):
    return f"{OUT}/{c}_gc{gc}_s{s}"


def run(j):
    c, gc, s = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", f + ".csv"] + BASE + [f"seed={s}"]
                       + DRUGS.get(c, []), check=True, stdout=lg, env=dict(os.environ, SMS_CAL="0"))


def S(c, gc, s):
    try:
        return summary(name(c, gc, s) + ".csv")
    except Exception:
        return None


def ctrl_ratio(gc):
    r = [x["a"]["smc"] / x["a"]["icc"] for s in SEEDS if (x := S("ctrl", gc, s)) and x["a"]]
    return np.mean(r), np.min(r), np.max(r)


def effect(d, gc, s):
    """drug effect relative to the same seed's control: dCSMC, dICC, ICC x, CSMC x, ICC dV/dt x (not saved), CSMC dV/dt x"""
    c0, c1 = S("ctrl", gc, s), S(d, gc, s)
    if not (c0 and c1 and c0["a"] and c1["a"]): return None
    return (c1["sr"] - c0["sr"], c1["ir"] - c0["ir"], c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"], c1["dv"] / c0["dv"])


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex:
            list(ex.map(run, JOBS))
    else:
        print("control ratio by g_c (mean, range over seeds 1-3):")
        cr = {gc: ctrl_ratio(gc) for gc in GC_CTRL}
        for gc in GC_CTRL: print(f"   g_c {gc}: {cr[gc][0]:.3f} ({cr[gc][1]:.3f}-{cr[gc][2]:.3f})")
        g = np.array([float(x) for x in GC_CTRL]); q = np.array([cr[x][0] for x in GC_CTRL])
        print("\ndrug levels at the drug g_c values (mean over seeds; dCSMC dICC transfer | ICC x | CSMC x | CSMC dV/dt x):")
        for d in DRUGS:
            for gc in GC_DRUG:
                e = [x for s in SEEDS if (x := effect(d, gc, s))]
                if not e: print(f"   {d:7s} g_c {gc}: no result"); continue
                m = np.mean(e, axis=0); lo, hi = np.min(e, axis=0), np.max(e, axis=0)
                print(f"   {d:7s} g_c {gc}: {m[0]:+5.1f} {m[1]:+5.1f} {m[1] / m[0]:.2f} | x{m[2]:.2f} ({lo[2]:.2f}-{hi[2]:.2f}) | "
                      f"x{m[3]:.2f} ({lo[3]:.2f}-{hi[3]:.2f}) | x{m[4]:.2f}")
        print("\nprediction for each data set (g_c from the control ratio; drug size interpolated to the recorded CSMC change):")
        for lab, (ratio, dk, dS, rI, rIx, rSx, rIdv, rSdv) in DATA.items():
            gfit = float(np.interp(ratio, q, g))
            gd = min(GC_DRUG, key=lambda x: abs(float(x) - gfit))
            lv = [k for k in DRUGS if k.startswith(dk)]
            rows = []
            for k in lv:
                e = [x for s in SEEDS if (x := effect(k, gd, s))]
                if e: rows.append(np.mean(e, axis=0))
            rows = sorted(rows, key=lambda r: r[0])
            dSs = np.array([r[0] for r in rows])
            if dk == "apa": p = rows[0]
            else: p = np.array([np.interp(dS, dSs, [r[i] for r in rows]) for i in range(5)])
            print(f"   {lab:20s} g_c fit {gfit:.2f} (runs at {gd}): CSMC {p[0]:+.1f} (rec {dS:+.1f}) | ICC {p[1]:+.1f} (rec {rI:+.1f}), "
                  f"transfer {p[1] / p[0]:.2f} (rec {rI / dS:.2f}) | ICC x{p[2]:.2f} (rec {rIx:.2f}) | CSMC x{p[3]:.2f} (rec {rSx:.2f}) | "
                  f"CSMC dV/dt x{p[4]:.2f} (rec {rSdv:.2f})")
