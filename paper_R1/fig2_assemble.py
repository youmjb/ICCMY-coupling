"""Figure 2 of the R1 paper, assembled from saved runs (no new simulations). Every number is recomputed from the logs/CSVs:
A traces (out/v4win, calibrated column), B amplitudes vs g_c, C g_c vs G_load across 8 geometries (ratio criterion SMC/ICC = 0.43;
out/v4win, out/paper_sens, out/paper_anat), D ICC-MY vs muscle amplitude along g_c for growing loads, E robustness of g_c* (absolute
criterion, SMC layer 1 = 24.4 +/- 2.7 mV) to the ICC-MY waveform (out/paper_wave), L-type block (out/paper_nif), the SMC model
(out/paper_poh) and the FLC model (out/paper_flcmf). Output: paper_R1/fig/fig2.png (+ .pdf) and fig2.out.txt (printed values).
"""
import os, re, csv, math
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle; figstyle.apply()
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
KITO, KSD, KICC, RATIO, GM = 24.4, 2.7, 56.9, 24.4 / 56.9, 0.96
ORANGE, BLUE, GREY = "#D9822B", "#2B5DA8", "0.88"


def rd(f):
    x = PAT.search(AE.log(f + ".log")); return float(x.group(2)), float(x.group(4))   # ICC amp, SMC L1 amp


def curve(pre, gcs):
    return [(g,) + rd(pre(g)) for g in gcs if os.path.exists(pre(g) + ".log")]


def cross(xs, ys, target):
    for (x0, y0), (x1, y1) in zip(zip(xs[:-1], ys[:-1]), zip(xs[1:], ys[1:])):
        if (y0 - target) * (y1 - target) <= 0 and y1 != y0:
            return math.exp(math.log(x0) + (target - y0) / (y1 - y0) * (math.log(x1) - math.log(x0)))
    return float("nan")


def gload(m, gss, L):
    lam = math.sqrt(gss / GM); return m * math.sqrt(gss * GM) * math.tanh(L / lam)


def gstar(c):
    g = [x[0] for x in c]; s = [x[2] for x in c]
    return [cross(g, s, KITO + d) for d in (-KSD, 0, KSD)]


out = []
fig = plt.figure(figsize=(figstyle.WIDTH, 5.4))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42, left=0.07, right=0.985, top=0.95, bottom=0.08)

# A traces
a = fig.add_subplot(gs[0, 0])
for k, (g, col, ls_) in enumerate(((0.5, "0.55", "-"), (1, "k", "--"), (3, "k", "-"))):
    r = list(csv.DictReader(open(f"out/v4win/m0.3_gSS30_gc{g}_u0.csv")))
    t = np.array([float(x["t_ms"]) for x in r]); vi = np.array([float(x["V_ICC"]) for x in r]); vs = np.array([float(x["V_SMC0"]) for x in r])
    up = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; up = up[t[up] >= 30000][0]
    sl = slice(up - 500, up + 2500); x = (t[sl] - t[up]) / 1000
    if k == 0: a.plot(x, vi[sl], color="k", lw=0.8, ls=":", label="ICC-MY")
    a.plot(x, vs[sl], color=col, ls=ls_, lw=1.0, label=f"muscle, g$_c$ {g} nS")
a.set_xlabel("time (sec)"); a.set_ylabel("V$_m$ (mV)"); a.legend(frameon=False, loc="upper right", handlelength=1.2)
a.set_title("A", loc="left", fontweight="bold")

# B amplitude vs g_c (calibrated column)
a = fig.add_subplot(gs[0, 1])
c = curve(lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [0.5, 1, 1.5, 2, 3, 5, 10, 20])
g = [x[0] for x in c]
a.axhspan(KITO - KSD, KITO + KSD, color="0.84", lw=0); a.axhspan(KICC - 6.4, KICC + 6.4, color="0.93", lw=0)
a.axhspan(34.1 - 4.6, 34.1 + 4.6, color="0.90", lw=0)
a.plot(g, [x[1] for x in c], "o-", color="k", mfc="white", ms=2.8, label="ICC-MY")
a.plot(g, [x[2] for x in c], "o-", color="k", ms=2.8, label="muscle layer 1")
gb = gstar(c[:6])
a.axvspan(gb[0], gb[2], color="k", alpha=0.10, lw=0)
a.set_xscale("log"); a.set_xticks([0.5, 1, 2, 5, 10, 20]); a.set_xticklabels(["0.5", "1", "2", "5", "10", "20"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("slow-wave amplitude (mV)"); a.legend(frameon=False, loc="lower right")
a.text(0.52, KITO, "rabbit CSMC", fontsize=6.5, color="0.4", va="center"); a.text(0.52, 34.1 + 4.6 + 0.5, "mouse CSMC", fontsize=6.5, color="0.4")
a.text(0.52, KICC + 6.4 + 1, "ICC-MY (rabbit; mouse 56.5 ± 5.3)", fontsize=6.5, color="0.4")
a.set_title("B", loc="left", fontweight="bold")
out.append(f"B: calibrated column g_c* {gb[1]:.2f} ({gb[0]:.2f}-{gb[2]:.2f}) nS")

# C g_c vs G_load (ratio criterion)
GEOM = [(0.3, 30, 30, lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [0.5, 1, 1.5, 2, 3, 5], "k"),
        (0.3, 15, 30, lambda g: f"out/paper_sens/gSS15/gc{g}_u0", [0.5, 1, 1.5, 2, 3], BLUE),
        (0.3, 60, 30, lambda g: f"out/paper_sens/gSS60/gc{g}_u0", [0.5, 1, 1.5, 2, 3, 4, 5], BLUE),
        (0.15, 30, 30, lambda g: f"out/paper_sens/m0.15/gc{g}_u0", [0.5, 1, 1.5, 2, 3], BLUE),
        (0.45, 30, 30, lambda g: f"out/paper_sens/m0.45/gc{g}_u0", [0.5, 1, 1.5, 2, 3, 4], BLUE),
        (0.3, 30, 15, lambda g: f"out/paper_sens/L15/gc{g}_u0", [0.5, 1, 1.5, 2, 3], BLUE),
        (0.3, 30, 45, lambda g: f"out/paper_sens/L45/gc{g}_u0", [0.5, 1, 1.5, 2, 3], BLUE),
        (1.3, 30, 7, lambda g: f"out/paper_anat/L7m1.3_gc{g}", [1, 2, 3, 4, 6, 8, 10, 12], "tab:red"),
        (3, 30, 8, lambda g: f"out/paper_anat/gSS30_gc{g}", [3, 6, 10, 15, 25, 40], "tab:red"),
        (3, 300, 8, lambda g: f"out/paper_anat/gSS300_gc{g}", [5, 10, 15, 20, 30, 50], "tab:red")]
a = fig.add_subplot(gs[0, 2])
GG = np.logspace(-0.3, 1.45, 50)
import sys; sys.path.insert(0, "paper_R1"); import null_model          # passive null model (same geometries)
RMOUSE, PURP = 34.1 / 56.5, "#7B4FA0"
for target, rt, col, mk in (("rabbit", RATIO, "k", "o"), ("mouse", RMOUSE, "0.5", "s")):
    NULL = {(null_model.GEO[d["name"]]): d["dc"] for d in null_model.results(target)[0]}
    ratios = []
    for m, gss, L, pre, gcs, _ in GEOM:
        cc = curve(pre, gcs + ([5] if (m, gss, L) == (0.45, 30, 30) else [])); G = gload(m, gss, L)
        gr = cross([x[0] for x in cc], [x[2] / x[1] for x in cc], rt)
        gn = NULL[(m, gss, L)]
        a.plot([G, G], [gr, gn], color="0.8", lw=0.6, zorder=1)
        a.plot([G], [gn], mk, mfc="white", mec=col, mew=0.8, ms=3.2, zorder=2)
        a.plot([G], [gr], mk, color=col, mec="k" if (m, gss, L) == (0.3, 30, 30) else col, mew=0.6, ms=3.5, zorder=3)
        if gr == gr: ratios.append(gr / G)
        out.append(f"C {target}: m {m} gSS {gss} L {L}: G_load {G:.2f} nS, g_c(ratio {rt:.3f}) {gr:.2f} nS, ratio {gr / G:.2f}; passive null {gn:.2f} nS ({gr / gn:.2f})")
    a.fill_between(GG, min(ratios) * GG, max(ratios) * GG, color=col, alpha=0.12, lw=0, zorder=0)
    out.append(f"C {target}: g_c/G_load range {min(ratios):.2f}-{max(ratios):.2f} (n {len(ratios)})")
a.set_xscale("log"); a.set_yscale("log")
a.set_xticks([1, 3, 10]); a.set_xticklabels(["1", "3", "10"]); a.set_yticks([1, 3, 10]); a.set_yticklabels(["1", "3", "10"])
a.set_xlabel("load per ICC-MY, G$_{load}$ (nS)"); a.set_ylabel("g$_c$ giving the recorded amplitude ratio (nS)")
a.plot([], [], "o", color="k", ms=3, label="rabbit, r 0.43"); a.plot([], [], "s", color="0.5", ms=3, label="mouse, r 0.60")
a.plot([], [], "o", mfc="white", mec="0.3", ms=3, label="passive cable"); a.legend(frameon=False, loc="lower right")
a.set_title("C", loc="left", fontweight="bold")

# D ICC-MY vs muscle amplitude along g_c
a = fig.add_axes([0.07, 0.08, 0.25, 0.36])
a.axhspan(KITO - KSD, KITO + KSD, color="0.84", lw=0); a.axvspan(KICC - 6.4, KICC + 6.4, color="0.93", lw=0)
for (m, gss, L, pre, gcs, _), col, mk in zip([next(g for g in GEOM if g[:3] == k) for k in ((0.3, 30, 30), (1.3, 30, 7), (3, 30, 8))],["k", "0.45", "k"], ["o", "s", "^"]):
    cc = curve(pre, gcs); G = gload(m, gss, L)
    a.plot([x[1] for x in cc], [x[2] for x in cc], mk + "-", color=col, mfc="white" if mk == "^" else col, ms=2.8, label=f"G$_{{load}}$ {G:.1f} nS")
a.set_xlabel("ICC-MY amplitude (mV)"); a.set_ylabel("muscle layer-1 amplitude (mV)")
a.legend(frameon=False, loc="upper left", title="g$_c$ increasing →", title_fontsize=6.5)
a.set_title("D", loc="left", fontweight="bold")

# E robustness forest (absolute criterion)
a = fig.add_axes([0.575, 0.08, 0.41, 0.36])
GC6 = [0.5, 1, 1.5, 2, 3, 5]
rows = [("calibrated model", lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", GC6, "k"),
        ("ICC-MY half-width 1.5 sec", lambda g: f"out/paper_wave/W4_gc{g}", GC6, "0.35"),
        ("ICC-MY half-width 2.4 sec, 10 min$^{-1}$", lambda g: f"out/paper_wave/W5_gc{g}", GC6, "0.35"),
        ("ICC-MY upstroke 5 msec", lambda g: f"out/paper_wave/W6_gc{g}", GC6, "0.35"),
        ("rabbit-like ICC-MY wave (−65 mV, 57 mV)", lambda g: f"out/paper_wave/W7_gc{g}", GC6, "0.35"),
        ("L-type block 90 %", lambda g: f"out/paper_nif/cal0.1/gc{g}_u0", GC6, "0.35"),
        ("L-type block 100 %", lambda g: f"out/paper_nif/cal0/gc{g}_u0", GC6, "0.35"),
        ("SMC model of Poh et al. 2012", lambda g: f"out/paper_poh/gK0_gc{g}", GC6, "0.35"),
        ("minimal FLC model", lambda g: f"out/paper_flcmf/gc{g}_u0", [0.5, 1, 1.5, 2, 3], "0.35")]
y = np.arange(len(rows))[::-1]
blo, _, bhi = gstar(curve(rows[0][1], rows[0][2]))          # calibrated-model range (+/- 1 SD of the recorded amplitude)
for yy, (lab, pre, gcs, col) in zip(y, rows):
    lo, mid, hi = gstar(curve(pre, gcs))
    a.plot([lo if not math.isnan(lo) else 0.5, hi], [yy, yy], color=col, lw=3, alpha=0.35, solid_capstyle="butt")
    a.plot([mid], [yy], "o", color=col, ms=3.5)
    out.append(f"E: {lab}: g_c* {mid:.2f} ({lo:.2f}-{hi:.2f}) nS")
a.axvspan(blo, bhi, color=GREY, lw=0, zorder=0)
a.set_yticks(y); a.set_yticklabels([r[0] for r in rows])
a.set_xscale("log"); a.set_xticks([0.5, 0.7, 1, 1.5, 2, 3]); a.set_xticklabels(["0.5", "0.7", "1", "1.5", "2", "3"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.set_xlim(0.45, 3)
a.set_xlabel("g$_c$ giving muscle layer 1 = 24.4 ± 2.7 mV (nS), calibrated column")
a.set_title("E", loc="left", fontweight="bold")

figstyle.save(fig, "paper_R1/fig/fig2")
open("paper_R1/fig2.out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out)); print("wrote paper_R1/fig/fig2.png, .pdf")
