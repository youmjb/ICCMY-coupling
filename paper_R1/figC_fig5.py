"""Model-C Figure 5: a low own conductance of ICC-MY, not stronger coupling, reconciles transfer, pacemaker amplitude and input resistance.
A  transfer of the muscle resting change to ICC-MY against the ICC-MY amplitude change at the pinacidil level (CSMC -10.3 mV):
   coupling and muscle variants with the original own conductance (grey; figA_runs.py, out/paper_figA), the own-conductance series
   (iccRest f, passive column, g_c 1.75 nS; lowg_coupled.py) and model C (final_model.py with FM_MODEL=C, seeds 1-3).
C  circular-muscle amplitude change at the pinacidil level for the same models (control ratio closest to 0.52).
B  input resistance of ICC-MY in a passive ICC-MY lattice (41 x 41, each node loaded by its column) against ICC-MY-ICC-MY coupling,
   with Cai et al. (260 +/- 38 MOhm, mean +/- SEM).
Output: fig/<M>_fig5 (.pdf/.tif/.png), fig<M>_fig5.out.txt with M = FM_MODEL (C, default, or D). With D, model D (resting currents x0.3)
is the reference and model C is shown as what full agreement requires."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
MODEL = os.environ.setdefault("FM_MODEL", "C")
import figstyle; figstyle.apply()
from active_cal import summary
import final_model as FM
from fig_lowg import sized, gin, GOWN
OUT = []
FA = "out/paper_figA"


def log(s): print(s); OUT.append(s)


def at_pin(k):                                                       # as figA.py
    try:
        c0 = summary(f"{FA}/{k}_ctrl.csv")
    except Exception:
        return None
    rows = []
    for p in ("p1", "p2", "p3", "p4", "p6", "p8"):
        if not os.path.exists(f"{FA}/{k}_{p}.csv"): continue
        try:
            s = summary(f"{FA}/{k}_{p}.csv")
        except Exception:
            continue
        if s["a"]: rows.append((s["sr"] - c0["sr"], s["ir"] - c0["ir"], s["a"]["icc"] / c0["a"]["icc"], s["a"]["smc"] / c0["a"]["smc"]))
    rows.sort()
    if len(rows) < 2 or not (rows[0][0] <= -10.3 <= rows[-1][0]): return None
    v = [np.interp(-10.3, [r[0] for r in rows], [r[i] for r in rows]) for i in range(4)]
    return v[1] / v[0], v[2], v[3], c0["a"]["smc"] / c0["a"]["icc"]


def model_c_pin():                                                   # as figA.fm_predict, pinacidil set (ratio 0.52)
    q = np.array([FM.ctrl_ratio(x)[0] for x in FM.GC_CTRL]); g = np.array([float(x) for x in FM.GC_CTRL])
    gfit = float(np.interp(0.52, q, g)); gd = min(FM.GC_DRUG, key=lambda x: abs(float(x) - gfit))
    rows = sorted([np.mean(e, axis=0) for k in FM.DRUGS if k.startswith("pin")
                   if (e := [x for s in FM.SEEDS if (x := FM.effect(k, gd, s))])], key=lambda r: r[0])
    p = [np.interp(-10.3, [r[0] for r in rows], [r[i] for r in rows]) for i in range(4)]
    return p[1] / p[0], p[2], p[3], gfit


FAMS = [("$g_c$ only", "gc", ["1.75", "2.75", "4", "6", "8"], "o"), ("+ longitudinal muscle", "lm2", ["1.75", "2.75", "4"], "D"),
        ("long transmural length constant", "lam", ["1.75", "2.75", "4"], "^"), ("active muscle", "act", ["1", "1.5", "2"], "v")]
fig = plt.figure(figsize=(figstyle.WIDTH, 5.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], height_ratios=[1.35, 0.8], wspace=0.35, hspace=0.45, left=0.2, right=0.97)
a, c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]); b = fig.add_subplot(gs[1, :])
bars = []
for lab, key, gcs, mk in FAMS:
    pts = [v for gc in gcs if (v := at_pin(f"{key}_{gc}"))]
    if not pts: continue
    a.plot([v[0] for v in pts], [v[1] for v in pts], "-", marker=mk, color="0.6", mfc="white", ms=4, lw=0.7,
           label=lab if key == "gc" else None)
    best = min(pts, key=lambda v: abs(v[3] - 0.52)); bars.append((lab.replace("long transmural length constant", "long λ"), best[2], "0.85"))
    log(f"{lab}: transfer {', '.join(f'{v[0]:.2f}' for v in pts)}; best-ratio CSMC x{best[2]:.2f}")
a.plot([], [], "-", color="0.6", lw=0.7, label="variants of the\noriginal-conductance model")
ser = [(f, sized(f)[0]) for f in ("1", "0.3", "0.2", "0.1", "0.05")]
a.plot([v[0] for _, v in ser], [v[1] for _, v in ser], "-", marker="s", color="k", ms=4, lw=1.0, label="resting conductances lowered")
for f, v in ser:
    a.annotate(f"×{f}", (v[0], v[1]), xytext=(6, -9), textcoords="offset points", fontsize=figstyle.SMALL, ha="left")
    log(f"own conductance x{f} (g_c 1.75, passive muscle): transfer {v[0]:.2f}, ICC x{v[1]:.2f}, CSMC x{v[2]:.2f}")
bars.append(("resting conductances ×0.1", dict(ser)["0.1"][2], "0.45"))
import re
PAT = r"K05 pinacidil\s+g_c fit ([\d.]+).*?transfer ([\d.]+) \(rec.*?ICC x([\d.]+) \(rec.*?CSMC x([\d.]+) \(rec"
def pin_of(m):
    g = re.search(PAT, open(f"paper_R1/final_model_{m}.out.txt", encoding="utf-8").read()).groups()
    return float(g[1]), float(g[2]), float(g[3]), float(g[0])
for m, mk, fc, lab in (("D", "P", "black", "model D (×0.3, active muscle)"), ("C", "P", "white", "model C (×0.1, active muscle)")):
    if m == "D" and MODEL != "D": continue
    mc = pin_of(m)
    a.plot(mc[0], mc[1], mk, mfc=fc, mec="k", mew=1.0, ms=10, zorder=5, label=lab)
    bars.append((f"model {m}", mc[2], "black" if m == MODEL else "0.3"))
    log(f"model {m} (g_c {mc[3]:.2f}, seeds 1-3): transfer {mc[0]:.2f}, ICC x{mc[1]:.2f}, CSMC x{mc[2]:.2f}")
a.plot(0.93, 1.17, "*", color="k", ms=11, mfc="white", mew=1.0, label="recorded (pinacidil)")
a.axhline(1, color="0.85", lw=0.5, zorder=0)
a.set_xlim(0.15, 1.12); a.set_ylim(0.8, 1.3)
a.set_xlabel("transfer of the muscle resting change"); a.set_ylabel("ICC-MY amplitude change (fold)")
a.legend(frameon=False, fontsize=figstyle.SMALL, loc="upper left")
a.text(-0.17, 1.02, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
b.barh(range(len(bars)), [x[1] for x in bars], color=[x[2] for x in bars], edgecolor="black", lw=0.5)
b.axvline(0.22, color="black", lw=1.0, ls="--"); b.text(0.235, -0.62, "recorded 0.22", fontsize=figstyle.SMALL, va="bottom")
b.set_yticks(range(len(bars))); b.set_yticklabels([x[0] for x in bars], fontsize=figstyle.SMALL); b.invert_yaxis()
b.set_xlabel("muscle amplitude change (fold)"); b.set_xlim(0, 1.05)
b.text(-0.25, 1.04, "C", transform=b.transAxes, fontsize=10, fontweight="bold")
gii = np.logspace(-1, 1.3, 30)
for gown, ls in ((1.95, "-"), (0.97, "--"), (0.19, ":")):          # own conductance measured in the column: original, model D, model C
    r = [1000.0 / gin(g, gown) for g in gii]
    c.plot(gii, r, ls, color="k", lw=1, label={1.95: "original-conductance", 0.97: "model D", 0.19: "model C"}[gown] + f", own {gown} nS")
    log(f"R_in, own conductance {gown} nS: alone {1000 / gin(0, gown):.0f} MOhm; ICC-MY coupling for 260 MOhm "
        f"{np.interp(-260, [-x for x in r], gii):.2f} nS")
c.axhspan(222, 298, color="0.88", lw=0); c.axhline(260, color="0.5", lw=0.6)
c.text(0.11, 300, "recorded in situ", fontsize=figstyle.SMALL, va="bottom", color="0.3")
c.set_xscale("log"); c.set_xlabel("coupling between neighbouring ICC-MY, g$_{II}$ (nS)"); c.set_ylabel("input resistance of ICC-MY (MΩ)")
c.legend(frameon=False, fontsize=figstyle.SMALL, loc="upper right")
c.text(-0.25, 1.02, "B", transform=c.transAxes, fontsize=10, fontweight="bold")
figstyle.save(fig, f"paper_R1/fig/{MODEL}_fig5"); plt.close(fig)
open(f"paper_R1/fig{MODEL}_fig5.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
