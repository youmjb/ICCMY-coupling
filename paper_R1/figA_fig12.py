"""Plan-A Figures 1 and 2 (grayscale, figstyle). Output: fig/A_fig1, fig/A_fig2 (.pdf/.tif/.png), figA_fig12.out.txt.
Fig. 1  A anatomy (fig1_schematic.draw panel A, grays); B model of plan A: ICC-MY (variant properties), circular column with FLC,
        optional longitudinal column, drug targets and the two recorded quantities.
Fig. 2  A-C as revision-4 Fig. 2 A-C (original ICC-MY model, out/v4win, out/paper_sens, out/paper_anat, null_model);
        D reference model of plan A (ICC-MY variant, nifedipine; out/paper_final, seeds 1-3) amplitude ratio against g_c with the control
        ratios of the three mouse data sets, and the original model (out/v4win) for comparison."""
import os, sys, re, csv, math
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import figstyle; figstyle.apply()
import fig1_schematic as sch
import amp_event as AE
OUT = []
MODEL = os.environ.get("FM_MODEL", "A")          # C: model C (ICC-MY own conductance x0.1, active muscle); outputs C_fig1, C_fig2
PRE = MODEL + "_" if MODEL in ("C", "D") else "A_"          # D: model D (resting currents x0.3, freshly prepared ICC)

# ======================= Figure 1 =======================
sch.ICC, sch.SMC, sch.FLC, sch.INK, sch.MUTED = "0.35", "0.88", "0.70", "black", "0.4"
sch.LONG, sch.RED, sch.BOX = "0.93", "black", "0.95"
fig = plt.figure(figsize=(figstyle.WIDTH, 3.6))
a = fig.add_axes([0.0, 0.02, 0.36, 0.93]); b = fig.add_axes([0.38, 0.02, 0.62, 0.93])
sch.draw(a, plt.figure().add_subplot(), k=0.75, title_a="", title_b="")   # panel A only (B drawn into a throwaway axes)
plt.close(2)
a.text(0.0, 9.9, "A", fontsize=10, fontweight="bold", va="top")
# --- B ---
b.set_xlim(0, 16); b.set_ylim(0, 10); b.axis("off")
b.text(0, 9.9, "B", fontsize=10, fontweight="bold", va="top")
fs, fsm = 7, 6.5
# longitudinal column (variant)
for i, y in enumerate((9.3, 8.8)):
    b.add_patch(Rectangle((1.4, y - 0.18), 2.2, 0.36, fc="white", ec="0.45", lw=0.6, ls="--"))
b.text(2.5, 9.3, "", ha="center", va="center", fontsize=fsm, color="0.4")
b.text(3.75, 9.05, "longitudinal column\n(variant, g$_{cL}$)", fontsize=fsm, color="0.4", va="center")
b.plot([2.5, 2.5], [8.62, 8.15], color="0.45", lw=1.0, ls=(0, (2, 1.5)))
# ICC-MY
b.add_patch(FancyBboxPatch((1.2, 7.25), 2.6, 0.85, boxstyle="round,pad=0.05", fc="0.35", ec="black", lw=0.8))
b.text(2.5, 7.68, "ICC-MY (18 pF)", ha="center", va="center", fontsize=fs, color="white", fontweight="bold")
b.text(4.05, 7.95, "store clock, insensitive to rest potential", fontsize=fsm, va="center")
b.text(4.05, 7.5, "upstroke: VDDR, inactivation midpoint −80 mV", fontsize=fsm, va="center")
if MODEL in ("C", "D"):
    b.text(4.05, 7.05, "resting conductances ×0.1" if MODEL == "C" else
           "resting conductances ×0.3", fontsize=fsm, va="center", fontweight="bold")
    b.text(6.15, 6.15, "SMC: voltage-dependent\ninward current", fontsize=fsm, va="center", color="0.25")
b.plot([2.5, 2.5], [7.2, 6.55], color="black", lw=1.4, ls=(0, (2, 1.5)))
b.text(2.65, 6.88, "g$_c$", fontsize=fs + 0.5, va="center", fontweight="bold")
rows = [(6.15, "layer 1"), (5.2, "layer 2"), (4.25, "layer 3"), (1.6, "layer 30")]
for y, lab in rows:
    b.add_patch(Rectangle((1.2, y - 0.35), 2.6, 0.7, fc="0.88", ec="black", lw=0.7))
    b.text(2.5, y, "m SMC", ha="center", va="center", fontsize=fsm)
    b.add_patch(Rectangle((4.6, y - 0.27), 1.2, 0.54, fc="0.70", ec="black", lw=0.7))
    b.text(5.2, y, "FLC", ha="center", va="center", fontsize=fsm)
    b.plot([3.8, 4.6], [y, y], color="black", lw=0.7)
    b.text(1.1, y, lab, ha="right", va="center", fontsize=fsm, color="0.4")
for (y1, _), (y2, _) in zip(rows[:2], rows[1:3]):
    b.plot([2.5, 2.5], [y1 - 0.35, y2 + 0.35], color="black", lw=1.0)
b.plot([2.5] * 3, [2.7, 3.05, 3.4], "o", color="black", ms=2)
b.text(2.5, 0.95, "circular column: 30 layers, sealed end", ha="center", fontsize=fsm, color="0.4")
# drug targets
def drug(txt, xy, xytext):
    b.annotate(txt, xy=xy, xytext=xytext, fontsize=fsm, va="center",
               arrowprops=dict(arrowstyle="-|>", color="black", lw=0.7, shrinkA=0, shrinkB=1))
b.annotate("pinacidil: K$_{ATP}$ ↑\n4-AP: K$^+$ ↓\nphenylephrine: K$^+$ ↑", xy=(3.8, 4.15), xytext=(6.6, 3.0), fontsize=fsm,
           va="center", arrowprops=dict(arrowstyle="-|>", color="black", lw=0.7, shrinkA=2, shrinkB=1))
drug("apamin: SK ↓", (5.8, 5.2), (6.6, 5.2))
b.text(6.6, 1.9, "drugs act on every layer\n(and on the longitudinal column in the variant)", fontsize=fsm, color="0.4", va="center")
# recorded quantities
def tag(x, y):
    b.add_patch(plt.Circle((x, y), 0.2, fc="black", ec="none")); b.text(x, y, "R", ha="center", va="center", fontsize=6, color="white", fontweight="bold")
tag(0.75, 7.68); tag(0.75, 6.65)
b.add_patch(FancyBboxPatch((10.3, 3.9), 5.5, 2.9, boxstyle="round,pad=0.1", fc="0.95", ec="none"))
b.text(10.45, 6.65, "Recorded and compared", fontsize=fs, fontweight="bold", va="top")
b.text(10.45, 6.05, "R  resting potential and slow-wave\n    amplitude of ICC-MY and of\n    circular muscle (layer 1)", fontsize=fsm, va="top")
b.text(10.45, 4.95, "control: ratio → g$_c$/G$_{load}$\ndrugs: resting transfer,\namplitude changes", fontsize=fsm, va="top", color="0.25")
figstyle.save(fig, f"paper_R1/fig/{PRE}fig1"); plt.close(fig)

# ======================= Figure 2 =======================
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
KITO, KSD, KICC, RATIO, GM = 24.4, 2.7, 56.9, 24.4 / 56.9, 0.96


def rd(f):
    x = PAT.search(AE.log(f + ".log")); return float(x.group(2)), float(x.group(4))


def curve(pre, gcs):
    return [(g,) + rd(pre(g)) for g in gcs if os.path.exists(pre(g) + ".log")]


def cross(xs, ys, target):
    for (x0, y0), (x1, y1) in zip(zip(xs[:-1], ys[:-1]), zip(xs[1:], ys[1:])):
        if (y0 - target) * (y1 - target) <= 0 and y1 != y0:
            return math.exp(math.log(x0) + (target - y0) / (y1 - y0) * (math.log(x1) - math.log(x0)))
    return float("nan")


def gload(m, gss, L):
    lam = math.sqrt(gss / GM); return m * math.sqrt(gss * GM) * math.tanh(L / lam)


fig = plt.figure(figsize=(figstyle.WIDTH, 5.2))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.32, left=0.08, right=0.985, top=0.95, bottom=0.09)
# A traces
a = fig.add_subplot(gs[0, 0])
for k, (g, col, ls_) in enumerate(((0.5, "0.55", "-"), (1, "k", "--"), (3, "k", "-"))):
    r = list(csv.DictReader(open(f"out/v4win/m0.3_gSS30_gc{g}_u0.csv")))
    t = np.array([float(x["t_ms"]) for x in r]); vi = np.array([float(x["V_ICC"]) for x in r]); vs = np.array([float(x["V_SMC0"]) for x in r])
    up = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; up = up[t[up] >= 30000][0]
    sl = slice(up - 500, up + 2500); x = (t[sl] - t[up]) / 1000
    if k == 0: a.plot(x, vi[sl], color="k", lw=0.8, ls=":", label="ICC-MY")
    a.plot(x, vs[sl], color=col, ls=ls_, lw=1.0, label=f"muscle, g$_c$ {g} nS")
a.set_xlabel("time (s)"); a.set_ylabel("V$_m$ (mV)"); a.legend(frameon=False, loc="upper right", handlelength=1.2)
a.set_title("A", loc="left", fontweight="bold")
# B amplitude vs g_c
a = fig.add_subplot(gs[0, 1])
c = curve(lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [0.5, 1, 1.5, 2, 3, 5, 10, 20])
g = [x[0] for x in c]
a.axhspan(KITO - KSD, KITO + KSD, color="0.84", lw=0); a.axhspan(KICC - 6.4, KICC + 6.4, color="0.93", lw=0)
a.axhspan(34.1 - 4.6, 34.1 + 4.6, color="0.90", lw=0)
a.plot(g, [x[1] for x in c], "o-", color="k", mfc="white", ms=2.8, label="ICC-MY")
a.plot(g, [x[2] for x in c], "o-", color="k", ms=2.8, label="muscle layer 1")
a.set_xscale("log"); a.set_xticks([0.5, 1, 2, 5, 10, 20]); a.set_xticklabels(["0.5", "1", "2", "5", "10", "20"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("slow-wave amplitude (mV)"); a.legend(frameon=False, loc="lower right")
a.text(0.52, KITO, "rabbit CSMC", fontsize=6.5, color="0.4", va="center"); a.text(0.52, 34.1 + 4.6 + 0.5, "mouse CSMC", fontsize=6.5, color="0.4")
a.text(0.52, KICC + 6.4 + 1, "ICC-MY (rabbit; mouse 56.5 ± 5.3)", fontsize=6.5, color="0.4")
a.set_title("B", loc="left", fontweight="bold")
# C g_c vs G_load
GEOM = [(0.3, 30, 30, lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [0.5, 1, 1.5, 2, 3, 5]),
        (0.3, 15, 30, lambda g: f"out/paper_sens/gSS15/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
        (0.3, 60, 30, lambda g: f"out/paper_sens/gSS60/gc{g}_u0", [0.5, 1, 1.5, 2, 3, 4, 5]),
        (0.15, 30, 30, lambda g: f"out/paper_sens/m0.15/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
        (0.45, 30, 30, lambda g: f"out/paper_sens/m0.45/gc{g}_u0", [0.5, 1, 1.5, 2, 3, 4, 5]),
        (0.3, 30, 15, lambda g: f"out/paper_sens/L15/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
        (0.3, 30, 45, lambda g: f"out/paper_sens/L45/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
        (1.3, 30, 7, lambda g: f"out/paper_anat/L7m1.3_gc{g}", [1, 2, 3, 4, 6, 8, 10, 12]),
        (3, 30, 8, lambda g: f"out/paper_anat/gSS30_gc{g}", [3, 6, 10, 15, 25, 40]),
        (3, 300, 8, lambda g: f"out/paper_anat/gSS300_gc{g}", [5, 10, 15, 20, 30, 50])]
a = fig.add_subplot(gs[1, 0])
GG = np.logspace(-0.3, 1.45, 50)
import null_model
for target, rt, col, mk in (("rabbit", RATIO, "k", "o"), ("mouse", 34.1 / 56.5, "0.5", "s")):
    NULL = {(null_model.GEO[d["name"]]): d["dc"] for d in null_model.results(target)[0]}
    ratios = []
    for m, gss, L, pre, gcs in GEOM:
        cc = curve(pre, gcs); G = gload(m, gss, L)
        gr = cross([x[0] for x in cc], [x[2] / x[1] for x in cc], rt); gn = NULL[(m, gss, L)]
        a.plot([G, G], [gr, gn], color="0.8", lw=0.6, zorder=1)
        a.plot([G], [gn], mk, mfc="white", mec=col, mew=0.8, ms=3.2, zorder=2)
        a.plot([G], [gr], mk, color=col, ms=3.5, zorder=3)
        if gr == gr: ratios.append(gr / G)
    a.fill_between(GG, min(ratios) * GG, max(ratios) * GG, color=col, alpha=0.12, lw=0, zorder=0)
    OUT.append(f"Fig2C {target}: g_c/G_load {min(ratios):.2f}-{max(ratios):.2f} (n {len(ratios)})")
a.set_xscale("log"); a.set_yscale("log")
a.set_xticks([1, 3, 10]); a.set_xticklabels(["1", "3", "10"]); a.set_yticks([1, 3, 10]); a.set_yticklabels(["1", "3", "10"])
a.set_xlabel("load per ICC-MY, G$_{load}$ (nS)"); a.set_ylabel("g$_c$ giving the recorded ratio (nS)")
a.plot([], [], "o", color="k", ms=3, label="rabbit, r 0.43"); a.plot([], [], "s", color="0.5", ms=3, label="mouse, r 0.60")
a.plot([], [], "o", mfc="white", mec="0.3", ms=3, label="passive cable"); a.legend(frameon=False, loc="lower right")
a.set_title("C", loc="left", fontweight="bold")
# D reference model ratio vs g_c
import final_model as FM
a = fig.add_subplot(gs[1, 1])
gg = [float(x) for x in FM.GC_CTRL]; q = [FM.ctrl_ratio(x) for x in FM.GC_CTRL]
X0, Y0 = (0.9, 0.38) if MODEL in ("C", "D") else (1.4, 0.44)
a.plot(gg, [x[0] for x in q], "o-", color="k", ms=3, label={"C": "model C (resting conductances ×0.1, active muscle)", "D": "model D (×0.3, active muscle)"}.get(MODEL, "reference model (ICC-MY variant, nifedipine)"))
co = curve(lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [1, 1.5, 2, 3])
a.plot([x[0] for x in co], [x[2] / x[1] for x in co], "o--", color="0.5", mfc="white", ms=3, label="published model (18 pF)")
for lab, r in (("pinacidil set", 0.52), ("apamin and\nphenylephrine sets", 0.59), ("4-AP set", 0.62)):
    gfit = float(np.interp(r, [x[0] for x in q], gg))
    a.plot([X0, gfit], [r, r], color="0.6", lw=0.6, ls=":"); a.plot([gfit, gfit], [Y0, r], color="0.6", lw=0.6, ls=":")
    a.text(gfit + 0.03, r - 0.004, f"{lab} ({r:.2f})", fontsize=6.5, va="top", color="0.3")
    OUT.append(f"Fig2D {lab}: control ratio {r:.2f} -> g_c {gfit:.2f} nS")
a.set_xlim(X0, 3.1); a.set_ylim(Y0, 0.68)
a.set_xlabel("g$_c$ (nS), calibrated column"); a.set_ylabel("muscle / ICC-MY amplitude ratio")
a.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.04), fontsize=6.5)
a.set_title("D", loc="left", fontweight="bold")
for gc_, r_ in zip(gg, q): OUT.append(f"Fig2D reference model g_c {gc_}: ratio {r_[0]:.3f} ({r_[1]:.3f}-{r_[2]:.3f})")
for x in co: OUT.append(f"Fig2D original model g_c {x[0]}: ratio {x[2] / x[1]:.3f}")
figstyle.save(fig, f"paper_R1/fig/{PRE}fig2"); plt.close(fig)
open(f"paper_R1/fig{PRE}fig12.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
print("\n".join(OUT))
