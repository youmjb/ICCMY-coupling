"""Figures 1 and 3-6 of the R1 paper, assembled in the style of fig2_assemble.py (7.2 in wide, Arial 7 pt, 300 dpi + PDF).
No new simulations and no second copy of the analysis: the values come from the existing analysis scripts, imported as modules
(importing them also re-writes their draft panels in paper_R1/fig/ and prints their tables).
  Fig. 1  fig1_schematic.draw (A, B) + analyze_paper_cal (C traces, D model vs Kito z-scores)
  Fig. 3  analytic R_in vs Cai 2017 (A, B; same formulas as fig_rin_cai.py) + analyze_paper_a C (direct R_in, own conductance)
  Fig. 4  analyze_loadrule (A collapse on G_load, B criteria)
  Fig. 5  analyze_paper_a B (IJP fraction vs g_c at three ICC-MY rests) + analyze_paper_sens (fraction vs own conductance)
  Fig. 6  analyze_paper_sens (one-at-a-time g_c* and IJP fraction)
  Fig. S1 analyze_paper_dmp (ICC-DMP input); S2 fig_supp_smc (SMC model); S3 analyze_paper_flc (FLC model)
Usage: python paper_R1/figs_assemble.py [1 3 4 5 6 S1 S2 S3]   (default: all). Output: paper_R1/fig/fig<N>.png/.pdf, paper_R1/figs.out.txt
"""
import os, sys, io, contextlib
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(os.path.dirname(HERE))
import figstyle
STYLE = figstyle.STYLE
ORANGE, BLUE, GREEN, GREY, PURPLE = "#D9822B", "#2B5DA8", "#3A8D4F", "0.88", "#7B4FA0"
G1, GOWN, GOWN_SIM = 5.37, 3.6, (2.3, 3.1)          # column conductance per m (nS); assumed and simulated ICC-MY own conductance
CAI, CAI_SEM, CAI_N = 260.0, 38.0, 13
OUT = []


def quiet_import(name):
    """import an analysis script without its printed table (the module still writes its draft panel)"""
    with contextlib.redirect_stdout(io.StringIO()):
        mod = __import__(name)
    plt.rcParams.update(plt.rcParamsDefault); plt.rcParams.update(STYLE)
    return mod


def title(a, s):
    a.set_title(s, loc="left", fontweight="bold")


def save(fig, n):
    figstyle.save(fig, f"paper_R1/fig/fig{n}"); plt.close(fig)
    OUT.append(f"wrote paper_R1/fig/fig{n}.png, .pdf")


def fig1():
    import fig1_schematic as sch
    sch.ICC, sch.SMC, sch.FLC, sch.LONG, sch.RED, sch.BOX, sch.INK, sch.MUTED = "0.42", "0.88", "0.70", "0.93", "k", "0.94", "0.1", "0.42"     # grayscale (main figure)
    cal = quiet_import("analyze_paper_cal")
    fig = plt.figure(figsize=(figstyle.WIDTH, 6.6))
    sch.draw(fig.add_axes([0.0, 0.50, 0.40, 0.49]), fig.add_axes([0.40, 0.50, 0.60, 0.49]), k=0.66,
             title_a="A", title_b="B")
    t, vi, vs, vp = cal.load("out/v4win/m0.3_gSS30_gc1.5_u0.csv")
    ta, _, _, vpa = cal.load("out/paper_cal/apa_gc1.5.csv")
    w = (t >= 30000) & (t < 38000); wa = (ta >= 30000) & (ta < 38000); x = (t[w] - 30000) / 1000
    for r, (v, lab, c) in enumerate(((vi, "ICC-MY", "k"), (vs, "muscle layer 1", "k"), (vp, "FLC, layer 1", "k"))):
        a = fig.add_axes([0.075, 0.345 - r * 0.135, 0.40, 0.105])
        a.plot(x, v[w], color=c, lw=0.9, label=lab)
        if r == 2:
            a.plot(x, vpa[wa], color="0.55", lw=0.8, label="FLC, SK blocked")
            a.set_xlabel("time (s)")
        else:
            a.set_xticklabels([])
        a.set_xlim(0, 8); a.set_ylabel("V$_m$ (mV)")
        a.legend(frameon=False, loc="lower right" if r == 2 else "upper center", ncol=1 if r == 2 else 2, handlelength=1.2, borderaxespad=0)
        if r == 0: title(a, "C")
    a = fig.add_axes([0.70, 0.065, 0.29, 0.385])
    keys = list(cal.KITO.keys()); y = np.arange(len(keys))[::-1]
    a.axvspan(-1, 1, color=GREY, lw=0); a.axvline(0, color="0.5", lw=0.6)
    for g, mk, c, fc in zip(cal.GCS, ("o", "s", "^"), ("0.45", "k", "k"), ("0.45", "k", "white")):
        z = np.array([(cal.M[g][k] - cal.KITO[k][0]) / cal.KITO[k][1] for k in keys])
        a.plot(np.clip(z, -8, 8), y, mk, color=c, mfc=fc, mec=c, ms=3.2, label=f"rabbit [3], g$_c$ {g} nS", mew=0.7)
        OUT.append(f"Fig1D rabbit g_c {g}: " + ", ".join(f"{k} {zz:+.1f}" for k, zz in zip(keys, z)))
    # mouse small intestine, Kito & Suzuki 2003 (35 C; ICC-MY n 39, circular muscle n 34; primary component), model at g_c 2 nS
    MOUSE = {"ICC-MY RMP (mV)": (-69.6, 4.5), "ICC-MY amplitude (mV)": (56.5, 5.3), "ICC-MY rate (min$^{-1}$)": (26.1, 2.4),
             "ICC-MY half-width (s)": (1.00, 0.16), "ICC-MY dV/dt max (V s$^{-1}$)": (1.51, 0.30), "CSMC RMP (mV)": (-69.9, 3.9),
             "CSMC amplitude (mV)": (34.1, 4.6), "CSMC half-width (s)": (1.07, 0.20), "CSMC dV/dt max (V s$^{-1}$)": (0.54, 0.10)}
    zm = [((cal.M[2][k] - MOUSE[k][0]) / MOUSE[k][1], yy) for k, yy in zip(keys, y) if k in MOUSE]
    a.plot([np.clip(z, -8, 8) for z, _ in zm], [yy for _, yy in zm], "D", mfc="white", mec="k", mew=0.9, ms=3.4,
           label="mouse [4], g$_c$ 2 nS")
    OUT.append("Fig1D mouse g_c 2: " + ", ".join(f"{k} {(cal.M[2][k] - MOUSE[k][0]) / MOUSE[k][1]:+.1f}" for k in keys if k in MOUSE))
    a.set_yticks(y); a.set_yticklabels([k.replace("CSMC", "muscle") for k in keys], fontsize=6.5)
    a.set_xlim(-8.5, 8.5); a.set_xticks([-8, -4, -1, 0, 1, 4, 8])
    a.set_xlabel("(model − recorded mean) / recorded SD")
    a.legend(frameon=True, framealpha=1, edgecolor="none", loc="lower left", handletextpad=0.2)
    fig.text(0.535, 0.462, "D", fontweight="bold", fontsize=8)
    save(fig, 1)


def fig3():
    pa = quiet_import("analyze_paper_a")
    sd, t95 = CAI_SEM * np.sqrt(CAI_N), stats.t.ppf(0.975, CAI_N - 1)
    fig, ax = plt.subplots(2, 2, figsize=(figstyle.WIDTH, 5.0), gridspec_kw=dict(hspace=0.5, wspace=0.32, left=0.08, right=0.985, top=0.95, bottom=0.09))
    ax = ax.ravel()
    gc = np.logspace(-1, np.log10(30), 200)
    MC = [(0.1, ":"), (0.3, "-"), (1, "--"), (3, "-.")]                 # grayscale: line style per m

    def cai_bands(a, lab=True):
        a.axhspan(CAI - sd, CAI + sd, color="0.94", lw=0, label="Cai mean ± SD" if lab else None)
        a.axhspan(CAI - t95 * CAI_SEM, CAI + t95 * CAI_SEM, color="0.86", lw=0, label="95 % CI of mean" if lab else None)
        a.axhspan(CAI - CAI_SEM, CAI + CAI_SEM, color="0.74", lw=0, label="mean ± SEM" if lab else None)
        a.axhline(CAI, color="k", lw=0.6)

    def logx(a):
        a.set_xscale("log"); a.set_xticks([0.1, 0.3, 1, 3, 10, 30]); a.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"])
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    a = ax[0]; cai_bands(a)
    for m, c in MC:
        G = m * G1; a.plot(gc, 1000 / (GOWN + gc * G / (gc + G)), color="k", ls=c, label=f"m {m}")
    a.axhline(1000 / GOWN, color="k", ls=":", lw=0.6); a.text(0.105, 1000 / GOWN + 5, "isolated ICC-MY, 278 MΩ", fontsize=6.5)
    logx(a); a.set_ylim(80, 420); a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("ICC-MY R$_{in}$, analytic (MΩ)")
    a.legend(frameon=False, loc="lower left", ncol=2, handlelength=1.2, columnspacing=0.8); title(a, "A")

    a = ax[1]
    labs = ["mean", "mean − SEM", "95 % CI\nlower limit", "mean − SD"]
    rins = [CAI, CAI - CAI_SEM, CAI - t95 * CAI_SEM, CAI - sd]
    for i, r in enumerate(rins):
        cap = 1000 / r - GOWN; lo, hi = 1000 / r - GOWN_SIM[1], 1000 / r - GOWN_SIM[0]
        a.bar(i - 0.12, cap, color="0.6", width=0.45)
        a.text(i - 0.12, cap + 0.08, f"{cap:.2f}", fontsize=6.5, ha="center", va="bottom")
        a.plot([i + 0.2, i + 0.2], [lo, hi], color="k", lw=2.2, solid_capstyle="butt")
        a.text(i + 0.27, (lo + hi) / 2, f"{lo:.1f}–\n{hi:.1f}", fontsize=6.5, va="center")
        OUT.append(f"Fig3B {labs[i].replace(chr(10), ' ')} R_in {r:.0f} MOhm: cap {cap:.2f} nS (g_own 3.6); {lo:.2f}-{hi:.2f} (g_own 3.1-2.3)")
    ge = 1 * 0.3 * G1 / (1 + 0.3 * G1)
    a.axhline(ge, color="k", ls="--", lw=0.8, label=f"calibrated column at g$_c$ 1 nS ({ge:.2f} nS)")
    a.set_xticks(range(4)); a.set_xticklabels([f"{l}\n{r:.0f} MΩ" for l, r in zip(labs, rins)], fontsize=6.5)
    a.set_xlim(-0.5, 3.7); a.set_ylim(0, 6.2)
    a.set_ylabel("largest effective load g$_c$G/(g$_c$+G) (nS)")
    a.plot([], [], "s", color="0.6", label="g$_{own}$ 3.6 nS (assumed)"); a.plot([], [], color="k", lw=2.2, label="g$_{own}$ 2.3–3.1 nS (simulated)")
    a.legend(frameon=False, loc="upper left"); title(a, "B")

    a = ax[2]; cai_bands(a, lab=False)
    for m, c, ls_, fc in ((0.3, "k", "-", "k"), (1, "k", "--", "white")):
        a.plot(gc, 1000 / (GOWN + gc * m * G1 / (gc + m * G1)), color=c, ls=ls_, lw=0.9, label=f"analytic, m {m}")
        for (mm, g), (rr, an) in pa.C.items():
            if mm == m and rr:
                v = [x[0] for x in rr]; x0 = max(g, 0.1) * (1.06 if m == 1 else 1)
                a.plot([x0] * len(v), v, "o", color="0.55", alpha=0.5, ms=2.2, mew=0)
                a.plot([x0], [np.mean(v)], "o", color=fc, mfc=fc, mec="k", mew=0.7, ms=4)
                OUT.append(f"Fig3C m {m} g_c {g}: R_in {np.mean(v):.0f} MOhm (n {len(v)}), analytic {an:.0f}")
    logx(a); a.set_ylim(80, 420); a.set_xlabel("g$_c$ (nS; 0 plotted at 0.1)"); a.set_ylabel("ICC-MY R$_{in}$, +10 pA pulse (MΩ)")
    a.legend(frameon=False, loc="lower left"); title(a, "C")

    a = ax[3]
    cols = {0: "white", 1: "0.85", 1.5: "0.65", 2: "0.45", 5: "0.25", 20: "k"}          # grayscale ramp with g_c
    for (m, g), (rr, an) in pa.C.items():
        if rr:
            a.plot([x[4] for x in rr], [x[2] for x in rr], "o" if m == 0.3 else "^", color=cols[g], mec="k", ms=3.4, mew=0.4,
                   label=f"m {m}, g$_c$ {g}")
    a.axhline(GOWN, color="0.5", ls=":", lw=0.8); a.text(-66.6, GOWN + 0.03, "assumed 3.6 nS", fontsize=6.5, color="0.4", va="bottom", ha="right")
    a.set_xlabel("ICC-MY V$_m$ before the pulse (mV)"); a.set_ylabel("ICC-MY own slope conductance (nS)")
    a.set_ylim(1.6, None); a.legend(frameon=False, loc="lower left", ncol=3, fontsize=6.5, handletextpad=0.1, columnspacing=0.6); title(a, "D")
    gall = [x[2] for (m, g), (rr, an) in pa.C.items() if g > 0 for x in rr]
    OUT.append(f"Fig3D coupled g_own {min(gall):.2f}-{max(gall):.2f} nS; isolated "
               f"{np.mean([x[2] for x in pa.C[(0.3, 0)][0]]):.2f} nS")
    save(fig, 3)


def fig4():
    lr = quiet_import("analyze_loadrule")
    fig, ax = plt.subplots(1, 3, figsize=(figstyle.WIDTH, 2.8), gridspec_kw=dict(width_ratios=[1, 1.15, 1], wspace=0.55, left=0.075, right=0.985, top=0.91, bottom=0.17))
    a = ax[0]
    for g, c, mk in zip(lr.GSS, ("0.55", "k", "k"), ("o", "s", "^")):
        p = [x for x in lr.pts if x[1] == g]
        a.plot([x[2] for x in p], [x[3] for x in p], mk + "-", color=c, mfc="white" if mk == "^" else c, ms=3.2, label=f"g$_{{SS}}$ {g} nS")
    a.axhline(24.4, color="0.5", ls=":", lw=0.8); a.text(4.2, 25.8, "rabbit CSMC, 24.4 mV", fontsize=6.5, color="0.4")
    a.set_xscale("log"); a.set_xticks([0.3, 1, 3, 10]); a.set_xticklabels(["0.3", "1", "3", "10"])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xlabel("G$_{load}$ = m·√(g$_{SS}$·g$_m$) (nS)"); a.set_ylabel("best muscle layer-1 amplitude\nwith ICC-MY protected (mV)")
    a.legend(frameon=False, loc="upper right"); title(a, "A")
    a = ax[1]                                        # one row per ICC-MY amplitude criterion; 9 rate/muscle combinations each
    AMP = [0.9, 0.8, 0.7]
    for y, amp in enumerate(AMP):
        sub = [x for x in lr.crit if abs(x[0] - amp) < 1e-9]
        for k, x in enumerate(sub):
            dy = (k - 4) * 0.045
            main = abs(x[1] - 0.15) < 1e-9 and x[2] == 20
            a.plot([x[3], x[4]], [y + dy] * 2, color="0.3" if main else "0.82", lw=1.2 if main else 0.6, zorder=1 if not main else 3)
            a.plot([x[3]], [y + dy], "o", color="k", ms=4 if main else 2.2, zorder=4)
            a.plot([x[4]], [y + dy], "x", color="0.4", ms=4.5 if main else 2.6, mew=1.0 if main else 0.6, zorder=4)
        OUT.append(f"Fig4B ICC-MY >= {amp:.0%}: largest with window {min(x[3] for x in sub):.2f}-{max(x[3] for x in sub):.2f} nS, "
                   f"smallest without {min(x[4] for x in sub):.2f}-{max(x[4] for x in sub):.2f} nS (9 rate/muscle criteria)")
    a.set_yticks(range(3)); a.set_yticklabels([f"ICC-MY ≥ {amp:.0%}" for amp in AMP]); a.set_ylim(-0.5, 2.5); a.invert_yaxis()
    a.set_xscale("log"); a.set_xlim(0.8, 10); a.set_xticks([1, 1.5, 2, 3, 5, 10]); a.set_xticklabels(["1", "1.5", "2", "3", "5", "10"])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xlabel("G$_{load}$ (nS):  ● largest with a window,  × smallest without")
    a.text(0.98, 0.02, "bold: rate ±15 %, muscle ≥ 20 mV\nthin: rate ±10–20 %, muscle 15–25 mV", transform=a.transAxes,
           fontsize=6.5, ha="right", va="bottom", color="0.35")
    title(a, "B")
    a = ax[2]                                        # ICC loss: each remaining ICC-MY serves m = 0.3/f, g_c per ICC-MY fixed
    il = quiet_import("icc_loss"); rows = il.loss()
    for sp, col, mk in (("rabbit", "k", "o"), ("mouse", "0.5", "s")):
        r = [d for d in rows if d["sp"] == sp]
        x = [100 * (1 - d["f"]) for d in r]
        a.plot(x, [d["smc"] for d in r], mk + "-", color=col, ms=3, label="muscle layer 1" if sp == "rabbit" else None)
        a.plot(x, [d["icc"] for d in r], mk + "--", color=col, mfc="white", ms=3, lw=0.8, label="ICC-MY" if sp == "rabbit" else None)
        a.text(1, r[0]["smc"] + 2.2, f"g$_c$ {r[0]['gc']} nS", fontsize=6.5, color=col)
        # passive divider with constant ICC-MY amplitude: A(f)/A(1) = f (g_c + G0) / (f g_c + G0), G0 = G_load of the intact column
        gcv, G0 = r[0]["gc"], r[0]["G"]; ff = np.linspace(0.08, 1, 60)
        a.plot(100 * (1 - ff), r[0]["smc"] * ff * (gcv + G0) / (ff * gcv + G0), ":", color=col, lw=1.0,
               label="passive divider" if sp == "rabbit" else None)
        fh = G0 / (gcv + 2 * G0)
        OUT.append(f"Fig4C {sp} passive divider (G0 {G0:.2f}, g_c {gcv}): half amplitude at {100 * (1 - fh):.0f} % lost")
        for d in r:
            OUT.append(f"Fig4C {sp} g_c {d['gc']}: lost {100 * (1 - d['f']):.0f} %, G_load {d['G']:.2f} nS, ICC-MY {d['icc']:.1f} mV, muscle {d['smc']:.1f} mV, rate {d['cpm_icc']:.1f}/{d['cpm_l0']:.1f}")
    a.set_xlabel("ICC-MY lost (%)"); a.set_ylabel("slow-wave amplitude (mV)"); a.set_xlim(-3, 93); a.set_ylim(0, 58)
    a.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.6), fontsize=6.5, handlelength=1.8); title(a, "C")
    save(fig, 4)


def figS4():
    """transmural decay of the slow wave (revision 3), saved as Fig. S4"""
    il = quiet_import("icc_loss"); tm = il.transmural()
    fig, ax = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 2.7), gridspec_kw=dict(wspace=0.32, left=0.08, right=0.985, top=0.9, bottom=0.17))
    a = ax[0]
    for name, col, ls in (("L30_rabbit", BLUE, "-"), ("L30_mouse", PURPLE, "-"), ("L7_rabbit", BLUE, ":"), ("L7_mouse", PURPLE, ":")):
        if name not in tm: continue
        d = tm[name]; k = np.arange(1, d["L"] + 1)
        a.plot(k, d["amp"], ls, marker="o", ms=2.2, color=col, label=f"{d['L']} layers, m {d['m']}, g$_c$ {d['gc']} nS")
        a.plot(k, d["passive"], "--", color=col, lw=0.6, alpha=0.6)
        OUT.append(f"FigS4A {name}: " + ", ".join(f"L{i + 1} {v:.1f}" for i, v in enumerate(d["amp"]) if i in (0, 4, 9, 19, 29) or i == d["L"] - 1))
    a.set_xlabel("muscle layer (1 = next to ICC-MY)"); a.set_ylabel("slow-wave amplitude (mV)"); a.set_xlim(0, 31); a.set_ylim(0, None)
    a.legend(frameon=False, loc="upper right"); title(a, "A")
    a = ax[1]
    import csv as _csv
    p = "out/paper_rev3/tm_L30_rabbit.csv.layers.csv"
    with open(p) as fh:
        r = _csv.reader(fh); h = next(r); arr = np.array([[float(x) for x in row] for row in r])
    t, vi = arr[:, 0], arr[:, 1]
    up = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; u = up[t[up] >= 30000][0]
    sl = slice(u - 500, u + 2800); tt = (t[sl] - t[u]) / 1000
    a.plot(tt, vi[sl], color=ORANGE, lw=0.9, label="ICC-MY")
    for k, c in zip((1, 5, 10, 20, 30), plt.cm.Blues(np.linspace(0.95, 0.4, 5))):
        a.plot(tt, arr[sl, 1 + k], color=c, lw=0.9, label=f"layer {k}")
    a.set_xlabel("time from ICC-MY upstroke (s)"); a.set_ylabel("V$_m$ (mV)"); a.legend(frameon=False, loc="upper right", fontsize=6.5)
    title(a, "B")
    save(fig, "S4")


def fig5():
    pa = quiet_import("analyze_paper_a"); se = quiet_import("analyze_paper_sens")
    fig, ax = plt.subplots(1, 3, figsize=(figstyle.WIDTH, 2.6), gridspec_kw=dict(wspace=0.45, left=0.095, right=0.985, top=0.91, bottom=0.18))
    a = ax[0]
    gg = np.linspace(0.8, 2.2, 50); a.plot(gg, gg / (gg + GOWN), color="0.5", lw=0.9, label="g$_c$/(g$_c$ + 3.6 nS)")
    for lab, c, dx, mk in (("0", "k", -0.04, "o"), ("-10", "0.5", 0, "s"), ("-20", "white", 0.04, "^")):
        xs, ms, rests = [], [], []
        for g in (1, 1.5, 2):
            rr, r1 = pa.B.get((g, lab), ([], np.nan))
            if not rr: continue
            q = [x[3] for x in rr]
            a.plot([g + dx] * len(q), q, mk, color="0.6", alpha=0.5, ms=2.2, mew=0)
            xs.append(g + dx); ms.append(np.mean(q)); rests.append(r1)
            OUT.append(f"Fig5A iccInj {lab} g_c {g}: rest {r1:.1f} mV, fraction {np.mean(q):.3f} (n {len(q)})")
        a.plot(xs, ms, mk + "-", color="k", mfc=c, ms=4, mec="k", mew=0.6, lw=0.8, label=f"ICC-MY rest {np.mean(rests):.0f} mV".replace("-", "−"))
    a.set_xticks([1, 1.5, 2]); a.set_ylim(0, 0.75)
    a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("ICC-MY hyperpolarization /\nmuscle layer-1 IJP")
    a.legend(frameon=False, loc="lower right"); title(a, "A")
    a = ax[1]
    xx = np.linspace(2, 8, 50); a.plot(xx, 1.5 / (1.5 + xx), color="0.5", lw=0.9, label="g$_c$/(g$_c$ + g$_{own}$), g$_c$ 1.5 nS")
    for c, g in {"base": 0.0, "iccGL-0.6": -0.6, "iccGL1.8": 1.8, "iccGL3.6": 3.6}.items():
        q = [x[2] for x in se.res[c]["ijp"]]
        a.plot([se.GOWN0 + g] * len(q), q, "o", color="0.6", alpha=0.5, ms=2.2, mew=0)
        a.plot([se.GOWN0 + g], [np.mean(q)], "o", color="0.55" if c != "base" else "k", mec="k", mew=0.5, ms=4.5)
        OUT.append(f"Fig5B g_own {se.GOWN0 + g:.2f}: fraction {np.mean(q):.3f} (n {len(q)})")
    a.set_ylim(0, 0.75)
    a.set_xlabel("ICC-MY own conductance at rest (nS)"); a.set_ylabel("IJP fraction at g$_c$ 1.5 nS")
    a.legend(frameon=False, loc="upper right"); title(a, "B")
    import network_divider as nd                     # passive ICC-MY network: local vs widespread muscle IJP
    d = nd.data(); a = ax[2]
    for (lab, _), c, ls_ in zip(nd.PATCH, ("k", "k", "k", "0.55"), ("-", "--", "-.", "-")):
        a.plot(nd.GII, d[lab], color=c, ls=ls_, lw=1.0 if c == "k" else 1.8, label=lab)
        OUT.append(f"Fig5C {lab}: fraction " + ", ".join(f"{g:.1f} nS {v:.3f}" for g, v in zip(nd.GII[::6], d[lab][::6])))
    a.set_xscale("log"); a.set_xticks([0.1, 1, 10, 100]); a.set_xticklabels(["0.1", "1", "10", "100"]); a.set_ylim(0, 0.35); a.set_xlabel("ICC-MY–ICC-MY coupling g$_{II}$ (nS)")
    a.set_ylabel("IJP fraction, passive network\n(g$_c$ 1.5 nS, g$_{own}$ 3.6 nS)")
    a.legend(frameon=False, loc="lower left", title="muscle IJP in", title_fontsize=6.5); title(a, "C")
    save(fig, 5)


LABELS = {"base": "calibrated model", "gSS15": "g$_{SS}$ 15 nS", "gSS60": "g$_{SS}$ 60 nS", "L15": "L 15 layers", "L45": "L 45 layers",
          "m0.15": "m 0.15", "m0.45": "m 0.45", "gK0run0.2": "tonic K$^+$ 0.2 nS", "gK0run0.6": "tonic K$^+$ 0.6 nS",
          "fP0.05": "f$_P$ 0.05", "fP0.15": "f$_P$ 0.15", "gSP1.5": "g$_{SP}$ 1.5 nS", "gSP4.5": "g$_{SP}$ 4.5 nS",
          "iccGL-0.6": "ICC-MY g$_{own}$ −0.6 nS", "iccGL1.8": "ICC-MY g$_{own}$ +1.8 nS", "iccGL3.6": "ICC-MY g$_{own}$ +3.6 nS",
          "seed2": "seed 2", "seed3": "seed 3"}


def fig6():
    se = quiet_import("analyze_paper_sens")
    names = [c for c, _ in se.CONDS]; y = np.arange(len(names))[::-1]
    fig, ax = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 3.6), sharey=True, gridspec_kw=dict(wspace=0.08, left=0.17, right=0.985, top=0.93, bottom=0.14))
    a = ax[0]
    for yy, c in zip(y, names):
        (glo, flo), (gm, fm), (ghi, fhi) = se.res[c]["gstar"]
        col = "k" if c == "base" else "0.45"
        a.plot([glo, ghi], [yy, yy], color=col, lw=3, alpha=0.3, solid_capstyle="butt")
        a.plot([gm], [yy], "o", color=col, ms=3.5)
        for g, f in ((glo, flo), (ghi, fhi)):
            if f: a.text(g, yy + 0.3, f, fontsize=6.5, ha="center")
        OUT.append(f"Fig6A {c}: g_c* {fm}{gm:.2f} ({flo}{glo:.2f}-{fhi}{ghi:.2f})")
    (blo, _), _, (bhi, _) = se.res["base"]["gstar"]
    a.axvspan(blo, bhi, color=GREY, lw=0, zorder=0)
    a.set_xscale("log"); a.set_xlim(0.45, 5.6); a.set_xticks([0.5, 1, 1.5, 2, 3, 5]); a.set_xticklabels(["0.5", "1", "1.5", "2", "3", "5"])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_yticks(y); a.set_yticklabels([LABELS[c] for c in names])
    a.set_xlabel("g$_c$ giving muscle layer 1 = 24.4 mV (nS)\n(bar, 21.7–27.1 mV)"); title(a, "A")
    a = ax[1]
    base = np.mean([x[2] for x in se.res["base"]["ijp"]])
    a.axvline(base, color="0.5", lw=0.6, ls=":")
    for yy, c in zip(y, names):
        q = [x[2] for x in se.res[c]["ijp"]]
        col = "k" if c == "base" else "0.45"
        a.plot(q, [yy] * len(q), "o", color=col, alpha=0.3, ms=2.5, mew=0)
        a.plot([np.mean(q)], [yy], "o", color=col, ms=4)
        a.text(0.455, yy, f"n {len(q)}", fontsize=6.5, va="center", color="0.45")
        OUT.append(f"Fig6B {c}: fraction {np.mean(q):.3f} (n {len(q)})")
    a.set_xlim(0.15, 0.47); a.set_xlabel("ICC-MY hyperpolarization / muscle IJP\nat g$_c$ 1.5 nS"); title(a, "B")
    save(fig, 6)


def figS1():
    dm = quiet_import("analyze_paper_dmp")
    fig, ax = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 2.8), gridspec_kw=dict(width_ratios=[1.3, 1], wspace=0.3, left=0.08, right=0.985, top=0.91, bottom=0.16))
    SITES = (("V_SMClast", PURPLE, "muscle next to ICC-DMP"), ("V_SMC0", BLUE, "muscle layer 1"), ("V_ICC", ORANGE, "ICC-MY"))
    a = ax[0]
    d0, d9 = dm.load(dm.path(10, 0, 1)), dm.load(dm.path(10, 9, 1))
    w = (d0["t_ms"] >= 40000) & (d0["t_ms"] < 48000); tt = d0["t_ms"][w] / 1000 - 40
    for k, c, lab in SITES:
        a.plot(tt, d0[k][w], color=c, lw=0.7, ls="--"); a.plot(tt, d9[k][w], color=c, lw=1.0, label=lab)
    a.plot([], [], color="0.4", lw=0.7, ls="--", label="DMP baseline"); a.plot([], [], color="0.4", lw=1.0, label="ACh ×10")
    a.set_xlim(0, 8); a.set_ylim(-75, 3); a.set_xlabel("time (s)"); a.set_ylabel("V$_m$ (mV)")
    a.legend(frameon=False, loc="upper center", ncol=3, handlelength=1.6, columnspacing=0.8, borderaxespad=0); title(a, "A")
    a = ax[1]
    x = np.arange(len(dm.LS)); bw = 0.26
    for j, (_, c, lab) in enumerate(SITES):
        vals = np.array([[dm.D[(L, s)][j] for s in dm.SEEDS] for L in dm.LS])
        a.bar(x + (j - 1) * bw, vals.mean(axis=1), bw, color=c, alpha=0.8, label=lab)
        for k in range(len(dm.LS)):
            a.plot([x[k] + (j - 1) * bw] * len(dm.SEEDS), vals[k], "o", color="k", ms=1.8, mew=0)
            OUT.append(f"FigS1B L {dm.LS[k]} {lab}: {vals[k].mean():.2f} mV (seeds {', '.join(f'{v:.2f}' for v in vals[k])})")
    a.axhline(0, color="k", lw=0.6)
    a.set_xticks(x); a.set_xticklabels([f"{L} layers" for L in dm.LS])
    a.set_ylim(0, 15.5); a.set_ylabel("depolarization by ACh ×10 (mV)")
    a.legend(frameon=False, loc="upper center", ncol=3, fontsize=6.5, handlelength=1, columnspacing=0.6, borderaxespad=0); title(a, "B")
    save(fig, "S1_DMP_removed")


def figS2():
    sm = quiet_import("fig_supp_smc")
    fig, ax = plt.subplots(2, 2, figsize=(figstyle.WIDTH, 5.4), gridspec_kw=dict(hspace=0.45, wspace=0.3, left=0.08, right=0.985, top=0.95, bottom=0.08))
    ax = ax.ravel()
    vv, g = sm.vv, sm.g
    p = ax[0].get_position(); ax[0].remove()                     # panel A = steady states (left 60 %) + time constants (right)
    a = fig.add_axes([p.x0, p.y0, p.width * 0.58, p.height]); a2 = fig.add_axes([p.x0 + p.width * 0.74, p.y0, p.width * 0.26, p.height])
    a.plot(vv, g[:, 0], color="tab:red", label="d$_\\infty$")
    a.plot(vv, g[:, 2], color="tab:red", ls="--", label="f$_\\infty$")
    a.plot(vv, g[:, 5] ** 2, color=BLUE, label="n$_\\infty^2$")
    for cai, ls in ((0.1, ":"), (1.0, "-.")):
        a.plot(vv, [sm.sm.kca_h(v, cai) for v in vv], color=GREEN, ls=ls, label=f"h$_\\infty$, {cai:g} µM Ca$^{{2+}}$")
    a.set_xticks([-80, -40, 0, 40]); a.set_xlabel("V (mV)"); a.set_ylabel("steady-state gate")
    for x_, y_, lab, c, ha in ((-32, 0.95, "d$_\\infty$", "tab:red", "right"), (-74, 0.82, "f$_\\infty$", "tab:red", "left"),
                               (14, 0.52, "n$_\\infty^2$", BLUE, "left"), (2, 0.27, "h$_\\infty$, 1 µM", GREEN, "left"),
                               (2, 0.045, "h$_\\infty$, 0.1 µM", GREEN, "left")):
        a.text(x_, y_, lab, color=c, fontsize=6.5, ha=ha, va="center")         # direct labels instead of a legend
    for col, c, ls, lab in ((1, "tab:red", "-", "d"), (3, "tab:red", "--", "f$_1$"), (4, "tab:red", ":", "f$_2$"), (6, BLUE, "-", "n")):
        a2.semilogy(vv, g[:, col], color=c, ls=ls, lw=0.8, label=lab)
    a2.set_ylim(1e-2, 1e3); a2.set_xticks([-80, -40, 0, 40]); a2.set_xlabel("V (mV)"); a2.set_ylabel("time constant (ms)", labelpad=1)
    a2.legend(frameon=False, loc="lower right", fontsize=6.5, handlelength=1.2, borderaxespad=0.1)
    title(a, "A")
    a = ax[1]
    a.plot(sm.steps, sm.pk, "o-", color="tab:red", ms=2.5, label="peak I$_{CaL}$")
    a.plot(sm.steps, sm.kend, "s-", color=BLUE, ms=2.5, label="I$_K$ + I$_{KCa}$ at 500 ms")
    a.axhline(0, color="k", lw=0.6); a.set_xlabel("test potential (mV); holding −70 mV"); a.set_ylabel("current density (pA/pF)")
    a.legend(frameon=False, loc="upper left"); title(a, "B")
    i0 = int(np.argmin(sm.pk)); OUT.append(f"FigS2B peak I_CaL {sm.pk[i0]:.1f} pA/pF at {sm.steps[i0]} mV; K+ at +40 mV {sm.kend[-1]:.1f} pA/pF")
    a = ax[2]
    t = sm.t
    for k, col in (("A5", "0.7"), ("A10", "0.45"), ("A20", ORANGE), ("A40", "tab:red")):
        a.plot(t - 100, [float(x[k]) for x in sm.r], color=col, lw=0.9, label=f"+{k[1:]} pA")
    a.axvspan(0, 900, color="0.94", lw=0, zorder=0)
    a.set_xlabel("time from current onset (ms)"); a.set_ylabel("V$_m$ (mV)"); a.legend(frameon=False, loc="upper right")
    title(a, "C")
    a = ax[3]
    for gk, col in ((0.0, "0.55"), (0.4, "k")):
        a.plot(sm.vq, [sm.iss(v, gk) for v in sm.vq], color=col, label=f"tonic K$^+$ {gk} nS")
    a.axhline(0, color="k", lw=0.6); a.set_ylim(-15, 25)
    a.set_xlabel("V (mV)"); a.set_ylabel("steady-state current (pA)"); a.legend(frameon=False, loc="upper left")
    title(a, "D")
    save(fig, "S1")


def figS3():
    fl = quiet_import("analyze_paper_flc")
    fig, ax = plt.subplots(2, 2, figsize=(figstyle.WIDTH, 5.0), gridspec_kw=dict(hspace=0.45, wspace=0.3, left=0.08, right=0.985, top=0.95, bottom=0.09))
    ax = ax.ravel()
    a = ax[0]
    for f, c, lab in (("ctrl", "k", "control"), ("apamin", BLUE, "SK blocked"), ("agonist", "tab:red", "IP$_3$ agonist from 10 sec")):
        t, v = fl.load(f"out/paper_flc/{f}.csv", ("t_ms", "V")); a.plot(t / 1000, v, color=c, lw=0.7, label=lab)
    a.set_xlabel("time (s)"); a.set_ylabel("V$_m$ (mV)")
    a.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.9), fontsize=6.5, handlelength=1.2); title(a, "A")
    a = ax[1]
    t, v = fl.load("out/v4win/m0.3_gSS30_gc1.5_u0.csv", ("t_ms", "V_P0")); w = (t >= 30000) & (t < 38000)
    ta, va = fl.load("out/paper_cal/apa_gc1.5.csv", ("t_ms", "V_P0")); wa = (ta >= 30000) & (ta < 38000)
    a.axhspan(-49.5 - 5.3, -49.5 + 5.3, color=GREY, lw=0, zorder=0); a.text(0.1, -54.4, "rabbit FLC RMP ± SD", fontsize=6.5, color="0.35", va="bottom")
    a.plot((t[w] - 30000) / 1000, v[w], color=GREEN, lw=0.9, label="FLC, layer 1")
    a.plot((ta[wa] - 30000) / 1000, va[wa], color="0.5", lw=0.9, label="FLC, SK blocked")
    a.set_xlim(0, 8); a.set_xlabel("time (s)"); a.set_ylabel("V$_m$ (mV)"); a.legend(frameon=False, loc="lower right"); title(a, "B")
    a = ax[2]
    cmap = {"k": "k", "tab:purple": PURPLE}
    a.axhspan(fl.KITO - fl.KSD, fl.KITO + fl.KSD, color=GREY, lw=0, zorder=0)
    for lab, pre, ip, c in fl.SETS:
        a.plot(fl.GCS, fl.R[lab]["smc"], "o-", color=cmap[c], ms=3, label=lab.replace("modified PDGFa FLC", "full FLC model").replace("minimal mean-field FLC", "minimal FLC model"))
        OUT.append(f"FigS3C {lab}: g_c* {fl.R[lab]['gs'][1]:.2f} ({fl.R[lab]['gs'][0]:.2f}-{fl.R[lab]['gs'][2]:.2f}) nS")
    a.set_xscale("log"); a.set_xticks(fl.GCS); a.set_xticklabels([f"{g:g}" for g in fl.GCS]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("muscle layer-1 slow wave (mV)"); a.legend(frameon=False, loc="upper left"); title(a, "C")
    a = ax[3]
    for k, (lab, pre, ip, c) in enumerate(fl.SETS):
        q = [x[2] for x in fl.R[lab]["rr"]]; sI = np.mean([x[0] for x in fl.R[lab]["rr"]])
        a.plot([k] * len(q), q, "o", color=cmap[c], alpha=0.35, ms=2.5, mew=0); a.plot([k], [np.mean(q)], "o", color=cmap[c], ms=5)
        a.text(k + 0.1, np.mean(q), f"{np.mean(q):.2f} (n {len(q)})\nmuscle IJP {sI:.1f} mV", fontsize=6.5, va="center")
        OUT.append(f"FigS3D {lab}: fraction {np.mean(q):.3f} (n {len(q)}), muscle IJP {sI:.2f} mV")
    a.set_xticks([0, 1]); a.set_xticklabels(["full FLC model", "minimal FLC model"]); a.set_xlim(-0.4, 1.7); a.set_ylim(0, 0.6)
    a.set_ylabel("ICC-MY hyperpolarization /\nmuscle IJP at g$_c$ 1.5 nS"); title(a, "D")
    save(fig, "S2")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    plt.rcParams.update(STYLE)
    want = sys.argv[1:] or ["1", "3", "4", "5", "6", "S2", "S3", "S4"]   # old S1 (ICC-DMP) dropped in revision 1; S2/S3 now saved as S1/S2
    for n in want: {"1": fig1, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "S1": figS1, "S2": figS2, "S3": figS3, "S4": figS4}[n]()
    open("paper_R1/figs.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
    print("\n".join(OUT))
