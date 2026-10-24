"""Plan-A figures (grayscale, journal form via figstyle.py). Output: fig/A_fig3, A_fig4, A_fig5 (.pdf/.tif/.png) and figA.out.txt.
Fig. 3  Drug experiments (Kito 2005, 2007; kito_data.md) against the final reference model (final_model.py):
        A resting change of ICC-MY against that of the circular muscle; B predicted against recorded amplitude changes.
Fig. 4  Isolated ICC-MY hyperpolarized by current: Youm 2019 against the variant (VDDR x16, midpoint -80 mV, IP3 V-independent);
        A amplitude, B dV/dt max, C rate, relative to rest; pinacidil in situ (-9.6 mV) for comparison.
Fig. 5  Coupling trade-off at the pinacidil level (CSMC -10.3 mV): A muscle-to-ICC-MY transfer against the ICC-MY amplitude change for
        four model families; B CSMC amplitude change. Runs: figA_runs.py (out/paper_figA)."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import figstyle
from active_cal import summary
import final_model as FM
from icc_hyperpol import stats
figstyle.apply()
OUTTXT = []
FA = "out/paper_figA"
PRE = FM.MODEL + "_" if FM.MODEL in ("C", "D") else "A_"         # FM_MODEL=C: model C (final_model.py), outputs C_fig3, C_fig4


def log(s):
    print(s); OUTTXT.append(s)


# ---------- recorded data (kito_data.md) ----------
REC = {  # label: dCSMC, dICC, CSMC x, ICC x
    "Pinacidil": (-10.3, -9.6, 0.22, 1.17),
    "Phenylephrine": (-2.4, -1.6, 37.0 / 34.7, 62.7 / 60.4),
    "Apamin": (3.1, 3.3, 36.0 / 36.2, 58.8 / 61.3),
    "4-AP": (3.1, 3.9, 45.3 / 37.0, 56.5 / 60.0),
}
MK = {"Pinacidil": "o", "Phenylephrine": "^", "Apamin": "s", "4-AP": "D"}
SDC = lambda s0, s1, rho=0.8: (s0 * s0 + s1 * s1 - 2 * rho * s0 * s1) ** 0.5   # SD of a change from the SDs before and in the drug
CI95 = {  # label: half-width of the 95 % interval of the mean change, CSMC and ICC-MY (SD/sqrt(n) x 1.96; kito_data.md)
    "Pinacidil": (1.96 * 1.0 / 11 ** 0.5, 1.96 * 2.0 / 8 ** 0.5),
    "Phenylephrine": (1.96 * 2.2 / 9 ** 0.5, 1.96 * 1.3 / 9 ** 0.5),
    "Apamin": (1.96 * SDC(2.8, 2.5) / 15 ** 0.5, 1.96 * SDC(2.9, 2.8) / 7 ** 0.5),      # within-tissue correlation 0.8 assumed
    "4-AP": (1.96 * SDC(1.2, 1.5) / 12 ** 0.5, 1.96 * SDC(3.4, 4.0) / 8 ** 0.5),
}
JIT = {"Pinacidil": 0, "Phenylephrine": 0, "Apamin": -0.18, "4-AP": 0.18}   # x offset (mV) so that apamin and 4-AP (both +3.1) do not overlap


def fm_predict(ratio, dk, dS):
    q = np.array([FM.ctrl_ratio(x)[0] for x in FM.GC_CTRL]); g = np.array([float(x) for x in FM.GC_CTRL])
    gfit = float(np.interp(ratio, q, g)); gd = min(FM.GC_DRUG, key=lambda x: abs(float(x) - gfit))
    rows = []
    for k in [k for k in FM.DRUGS if k.startswith(dk)]:
        e = [x for s in FM.SEEDS if (x := FM.effect(k, gd, s))]
        if e: rows.append(np.mean(e, axis=0))
    rows.sort(key=lambda r: r[0])
    if dk == "apa": return gfit, rows[0]
    return gfit, np.array([np.interp(dS, [r[0] for r in rows], [r[i] for r in rows]) for i in range(5)])


PRED = {}
for lab, (ratio, dk) in {"Pinacidil": (0.52, "pin"), "Phenylephrine": (0.59, "pe"), "Apamin": (0.59, "apa"), "4-AP": (0.62, "4ap")}.items():
    gfit, p = fm_predict(ratio, dk, REC[lab][0]); PRED[lab] = p
    log(f"Fig3 {lab}: g_c {gfit:.2f}; model dCSMC {p[0]:+.1f} dICC {p[1]:+.1f} CSMC x{p[3]:.2f} ICC x{p[2]:.2f} | "
        f"recorded {REC[lab][0]:+.1f} {REC[lab][1]:+.1f} x{REC[lab][2]:.2f} x{REC[lab][3]:.2f}")

fig, (a, b) = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 2.9))
lim = 12
a.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.7, ls="--", zorder=0)
a.axhline(0, color="0.85", lw=0.5, zorder=0); a.axvline(0, color="0.85", lw=0.5, zorder=0)
for lab in REC:
    dS, dI = REC[lab][:2]; p = PRED[lab]
    a.errorbar(dS + JIT[lab], dI, xerr=CI95[lab][0], yerr=CI95[lab][1], fmt="none", ecolor="black", elinewidth=0.6, capsize=1.5, zorder=1)
    a.plot(dS + JIT[lab], dI, MK[lab], color="black", ms=6, label=f"{lab}, recorded", zorder=3)
    log(f"Fig3 {lab}: 95 % interval of the recorded mean change, CSMC +-{CI95[lab][0]:.2f}, ICC +-{CI95[lab][1]:.2f} mV")
    a.plot(p[0] + JIT[lab], p[1], MK[lab], mfc="white", mec="black", ms=6)
    a.annotate("", xy=(p[0] + JIT[lab], p[1]), xytext=(dS + JIT[lab], dI), arrowprops=dict(arrowstyle="-", color="0.6", lw=0.6, ls=":"))
a.plot([], [], "o", mfc="white", mec="black", label=f"model {FM.MODEL} (open)" if FM.MODEL in ("C", "D") else "model (open)")
if FM.MODEL in ("C", "D"):                                # the same test with the original own conductance (draft A1 model), grey
    import re
    PAT = (r"(K0\S+ (\S+))\s+g_c fit .*?CSMC ([+-][\d.]+) \(rec .*?\| ICC ([+-][\d.]+) \(rec .*?ICC x([\d.]+) \(rec .*?CSMC x([\d.]+) \(rec")
    OLD = {m.group(2): [float(m.group(k)) for k in (3, 4, 5, 6)] for m in re.finditer(PAT, open("paper_R1/final_model.out.txt", encoding="utf-8").read())}
    for lab in REC:
        o = OLD[{"Pinacidil": "pinacidil", "Phenylephrine": "phenylephrine", "Apamin": "apamin", "4-AP": "4-AP"}[lab]]
        a.plot(o[0] + JIT[lab], o[1], MK[lab], mfc="none", mec="0.6", ms=4.5, mew=0.6)
        log(f"Fig3 {lab}: original own conductance dCSMC {o[0]:+.1f} dICC {o[1]:+.1f} ICC x{o[2]:.2f} CSMC x{o[3]:.2f}")
    a.plot([], [], "o", mfc="none", mec="0.6", ms=4.5, mew=0.6, label="original-conductance model")
    if FM.MODEL == "D":                                   # model C (resting currents x0.1): what full agreement requires
        MC = {m.group(2): [float(m.group(k)) for k in (3, 4, 5, 6)] for m in re.finditer(PAT, open("paper_R1/final_model_C.out.txt", encoding="utf-8").read())}
        for lab in REC:
            o = MC[{"Pinacidil": "pinacidil", "Phenylephrine": "phenylephrine", "Apamin": "apamin", "4-AP": "4-AP"}[lab]]
            a.plot(o[0] + JIT[lab], o[1], MK[lab], color="0.45", ms=3.5)
            log(f"Fig3 {lab}: model C dCSMC {o[0]:+.1f} dICC {o[1]:+.1f} ICC x{o[2]:.2f} CSMC x{o[3]:.2f}")
        a.plot([], [], "o", color="0.45", ms=3.5, label="resting conductances ×0.1 (model C)")
a.set_xlim(-lim, 6); a.set_ylim(-lim, 6)
a.set_xlabel("change of circular-muscle resting potential (mV)"); a.set_ylabel("change of ICC-MY resting potential (mV)")
a.text(-6.8, -7.1, "identity (complete transfer)", fontsize=figstyle.SMALL, color="0.4", rotation=45, rotation_mode="anchor", ha="center", va="top")
a.legend(frameon=False, loc="upper left", fontsize=figstyle.SMALL)
a.text(-0.2, 1.02, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
b.plot([0, 1.4], [0, 1.4], color="0.6", lw=0.7, ls="--", zorder=0)
for lab in REC:
    p = PRED[lab]
    b.plot(REC[lab][2], p[3], MK[lab], color="black", ms=6)
    b.plot(REC[lab][3], p[2], MK[lab], mfc="white", mec="black", ms=6)
b.plot([], [], "o", color="black", label="circular muscle"); b.plot([], [], "o", mfc="white", mec="black", label="ICC-MY")
b.set_xlim(0, 1.4); b.set_ylim(0, 1.4)
b.set_xlabel("recorded amplitude change (fold)"); b.set_ylabel("predicted amplitude change (fold)")
b.annotate("pinacidil,\nmuscle", xy=(0.22, PRED["Pinacidil"][3]), xytext=(0.3, 1.05), fontsize=figstyle.SMALL,
           arrowprops=dict(arrowstyle="->", lw=0.5))
b.legend(frameon=False, loc="lower right", fontsize=figstyle.SMALL)
b.text(-0.2, 1.02, "B", transform=b.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(); figstyle.save(fig, f"paper_R1/fig/{PRE}fig3"); plt.close(fig)

# ---------- Fig. 4: isolated ICC-MY ----------
fig, ax = plt.subplots(1, 3, figsize=(figstyle.WIDTH, 2.3))
for m, ls, lab in (("orig", "--", "published model"), ("var", "-", "this study")):
    r = [stats(f"{FA}/iso_{m}_i{i}.csv") for i in ("0", "5", "10", "15", "25", "35", "50")]
    r0 = r[0]
    dV = [x[0] - r0[0] for x in r]
    for k, (col, ylab) in enumerate(((1, "amplitude (fold)"), (3, "dV/dt max (fold)"), (4, "rate (fold)"))):
        ax[k].plot(dV, [x[col] / r0[col] for x in r], ls, color="black", marker="o", ms=3, label=lab)
    log(f"Fig4 {lab}: control rest {r0[0]:.1f} amp {r0[1]:.1f} dV/dt {r0[3]:.2f} rate {r0[4]:.1f}; at {dV[-2]:.1f} mV: amp x{r[-2][1] / r0[1]:.2f} "
        f"dV/dt x{r[-2][3] / r0[3]:.2f} rate x{r[-2][4] / r0[4]:.2f}")
for k, (val, ylab) in enumerate(((1.17, "amplitude (fold)"), (1.52, "dV/dt max (fold)"), (1.0, "rate (fold)"))):
    ax[k].plot(-9.6, val, "*", color="black", ms=9, label="in situ, pinacidil")
    ax[k].set_xlabel("ΔV (mV)"); ax[k].set_ylabel(ylab); ax[k].axhline(1, color="0.85", lw=0.5, zorder=0)
    ax[k].text(-0.3, 1.02, "ABC"[k], transform=ax[k].transAxes, fontsize=10, fontweight="bold")
ax[2].legend(*ax[0].get_legend_handles_labels(), frameon=False, fontsize=figstyle.SMALL, loc="lower right")   # empty corner of C
fig.tight_layout(); figstyle.save(fig, f"paper_R1/fig/{PRE}fig4"); plt.close(fig)

# ---------- Fig. 5: trade-off ----------
def at_pin(k):
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


FAMS = [("$g_c$ only", "gc", ["1.75", "2.75", "4", "6", "8"], "o", "-"), ("+ longitudinal muscle, $g_{cL}$ = $g_c$", "lm1", ["1.75", "2.75", "4"], "s", "--"),
        ("+ longitudinal muscle, $g_{cL}$ = 2 $g_c$", "lm2", ["1.75", "2.75", "4"], "D", "--"), ("long transmural length constant", "lam", ["1.75", "2.75", "4"], "^", ":"),
        ("active muscle", "act", ["1", "1.5", "2"], "v", "-.")]
if FM.MODEL not in ("C", "D"):                           # models C, D: Fig. 5 is made by figC_fig5.py
    fig, (a, b) = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 2.9), gridspec_kw={"width_ratios": [1.3, 1]})
    labs, cs = [], []
    for lab, key, gcs, mk, ls in FAMS:
        pts = [(gc, at_pin(f"{key}_{gc}")) for gc in gcs]
        pts = [(gc, v) for gc, v in pts if v]
        if not pts: continue
        a.plot([v[0] for _, v in pts], [v[1] for _, v in pts], ls, marker=mk, color="black", mfc="white", ms=5, lw=0.8, label=lab)
        for gc, v in pts:
            log(f"Fig5 {lab} g_c {gc}: control ratio {v[3]:.2f}, transfer {v[0]:.2f}, ICC x{v[1]:.2f}, CSMC x{v[2]:.2f}")
        best = min(pts, key=lambda t: abs(t[1][3] - 0.52))
        labs.append(lab.replace("+ longitudinal muscle", "+ LM").replace("long transmural length constant", "long lambda")); cs.append(best[1][2])
    a.plot(0.93, 1.17, "*", color="black", ms=11, label="in situ (pinacidil)")
    a.set_xlabel("muscle-to-ICC-MY transfer of the resting change"); a.set_ylabel("ICC-MY amplitude change (fold)")
    a.set_xlim(0.2, 1.0); a.axhline(1, color="0.85", lw=0.5, zorder=0)
    a.legend(frameon=False, fontsize=figstyle.SMALL, loc="upper left", bbox_to_anchor=(0.0, 1.02))
    a.set_ylim(0.8, 1.22)
    a.text(-0.17, 1.02, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
    b.barh(range(len(cs)), cs, color="0.7", edgecolor="black", lw=0.5)
    b.axvline(0.22, color="black", lw=1.0, ls="--"); b.text(0.23, -0.62, "in situ 0.22", fontsize=figstyle.SMALL, va="bottom")
    b.set_yticks(range(len(cs))); b.set_yticklabels(labs, fontsize=figstyle.SMALL); b.invert_yaxis()
    b.set_xlabel("circular-muscle amplitude change (fold)\nat control ratio closest to 0.52"); b.set_xlim(0, 1)
    b.text(-0.75, 1.02, "B", transform=b.transAxes, fontsize=10, fontweight="bold")
    fig.tight_layout(); figstyle.save(fig, "paper_R1/fig/A_fig5"); plt.close(fig)
open("paper_R1/figA.out.txt" if PRE == "A_" else f"paper_R1/fig{FM.MODEL}.out.txt", "w", encoding="utf-8").write("\n".join(OUTTXT) + "\n")
