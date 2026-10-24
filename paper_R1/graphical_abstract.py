"""Graphical abstract for KJPP (4:3, 2000 x 1500 px at 300 dpi, sans-serif 10-14 pt, RGB).
Left: model schematic with the model's own slow waves (out/v4win/m0.3_gSS30_gc1_u0: ICC-MY 52.9 mV, muscle layer 1 23.2 mV at g_c 1 nS).
Right: g_c at the recorded amplitude ratio (0.43) against G_load for the ten column geometries (paper_R1/criteria.out.txt, group "geometry");
bands = event-aligned g_c/G_load ranges (rabbit 0.66-0.75, mouse 1.23-1.49); shading = loads at which one 18-pF ICC-MY no longer reproduced both recorded rabbit amplitudes (reproduced up to 2.4 nS, not at 5.9 nS).
Output: paper_R1/fig/graphical_abstract.png / .tif (RGB) / .pdf
"""
import os, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse
from PIL import Image
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "mathtext.default": "regular", "pdf.fonttype": 42,
                     "axes.linewidth": 0.8})
ICC, MUS, FLC, INK, MUTED, RED = "#D9822B", "#2B5DA8", "#3A8D4F", "#1B2440", "#5F6673", "#C0392B"
W, H = 2000 / 300, 1500 / 300

# --- data ---
r = list(csv.DictReader(open("out/v4win/m0.3_gSS30_gc1_u0.csv")))
t = np.array([float(x["t_ms"]) for x in r]); vi = np.array([float(x["V_ICC"]) for x in r]); vs = np.array([float(x["V_SMC0"]) for x in r])
up = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; up = up[t[up] >= 30000][0]
sl = slice(up - 600, up + 2700); tt = (t[sl] - t[up]) / 1000
log = AE.log("out/v4win/m0.3_gSS30_gc1_u0.log")
amp_i, amp_s = float(re.search(r"amp ([\d.]+)", log).group(1)), float(re.search(r"\[0\]([\d.]+)", log).group(1))
pts = []
for line in open("paper_R1/criteria.out.txt", encoding="utf-8"):
    f = [x.strip() for x in line.split("|")]
    if len(f) > 5 and f[1] == "geometry":
        pts.append((f[0], float(f[2]), float(f[5]), float(f[9])))   # name, G_load, g_c at rabbit 0.429, at mouse 0.604

fig = plt.figure(figsize=(W, H))
fig.text(0.03, 0.955, "How strongly is the intestinal pacemaker coupled to the muscle?", fontsize=13, fontweight="bold", va="top", color=INK)

# --- schematic ---
AX = [0.0, 0.03, 0.52, 0.85]
a = fig.add_axes(AX); a.set_xlim(0, 10); a.set_ylim(0.8, 9.3); a.axis("off")
sx, sy = AX[2] * W / 10, AX[3] * H / 8.5                     # inches per data unit, for round FLC symbols
a.add_patch(FancyBboxPatch((1.1, 8.0), 3.0, 1.1, boxstyle="round,pad=0.08", fc=ICC, ec=INK, lw=0.8))
a.text(2.6, 8.55, "pacemaker\n(ICC-MY)", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
a.plot([2.6, 2.6], [8.0, 6.95], color=INK, lw=1.6, ls=(0, (2, 1.5)))
a.text(2.85, 7.45, "g$_c$", fontsize=12, fontweight="bold", color=INK, va="center")
ys = [6.4, 5.4, 4.4, 1.6]
for y in ys:
    a.add_patch(Rectangle((1.1, y - 0.4), 3.0, 0.8, fc="#D5E3F0", ec=INK, lw=0.7))
    a.add_patch(Ellipse((4.55, y), 0.16 / sx, 0.16 / sy, fc=FLC, ec=INK, lw=0.6))
    a.plot([4.1, 4.55 - 0.08 / sx], [y, y], color=INK, lw=0.7)
for y1, y2 in zip(ys[:2], ys[1:3]):
    a.plot([2.6, 2.6], [y1 - 0.4, y2 + 0.4], color=INK, lw=1.0)
a.plot([2.6] * 3, [2.6, 3.0, 3.4], "o", color=INK, ms=2.5)
a.text(2.6, 6.4, "muscle layer 1", ha="center", va="center", fontsize=10, color=INK)
a.text(4.55, 6.4 + 0.12 / sy, "FLC", ha="center", va="bottom", fontsize=10, color=FLC)
a.annotate("", xy=(0.75, 1.2), xytext=(0.75, 6.8), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.9))
a.text(0.6, 4.0, "muscle load G$_{load}$", rotation=90, ha="right", va="center", fontsize=10, color=MUTED)
a.annotate("inhibitory\nnerve input", xy=(4.75, 5.3), xytext=(5.3, 3.2), fontsize=10, color=RED, ha="left",
           arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.0))
a.text(5.3, 1.75, "25–41 % reaches the\npacemaker if widespread,\nless if local", fontsize=10, color=RED, va="center")
# slow waves of the model at g_c 1 nS
for v, col, y0, amp in ((vi, ICC, 8.0, amp_i), (vs, MUS, 5.75, amp_s)):
    ax = a.inset_axes([0.52, (y0 - 0.8) / 8.5, 0.46, 0.14], transform=a.transAxes)
    ax.plot(tt, v[sl], color=col, lw=1.4); ax.set_ylim(-75, -15); ax.axis("off")
    ax.text(1.0, 0.95, f"{amp:.0f} mV", transform=ax.transAxes, ha="right", va="top", fontsize=11, fontweight="bold", color=col)
a.text(7.5, 7.3, "slow wave shrinks\nacross g$_c$", ha="center", va="center", fontsize=10, color=MUTED)

# --- result: g_c vs G_load ---
b = fig.add_axes([0.64, 0.50, 0.33, 0.36])
GG = np.logspace(np.log10(0.6), np.log10(30), 60)
b.axvspan(2.4, 5.9, color="#FBEFED", lw=0, zorder=0); b.axvspan(5.9, 30, color="#F3D9D5", lw=0, zorder=0)
PURP = "#7B4FA0"
b.fill_between(GG, 0.66 * GG, 0.75 * GG, color=MUS, alpha=0.15, lw=0, zorder=1)
b.fill_between(GG, 1.23 * GG, 1.49 * GG, color=PURP, alpha=0.15, lw=0, zorder=1)
for name, G, gr, gm in pts:
    b.plot(G, gr, "o", color=MUS, ms=3.5, zorder=3); b.plot(G, gm, "s", color=PURP, ms=3.5, zorder=3)
b.set_xscale("log"); b.set_yscale("log"); b.set_xlim(0.6, 30); b.set_ylim(0.3, 30)
b.set_xticks([1, 3, 10, 30]); b.set_xticklabels(["1", "3", "10", "30"]); b.set_yticks([0.3, 1, 3, 10]); b.set_yticklabels(["0.3", "1", "3", "10"])
b.minorticks_off(); b.tick_params(labelsize=10)
b.set_xlabel("muscle load G$_{load}$ (nS)", fontsize=10); b.set_ylabel("coupling g$_c$ (nS)", fontsize=10)
b.text(0.65, 21, "mouse: 1.23–1.49 G$_{load}$", fontsize=10, color=PURP)
b.text(0.65, 11.5, "rabbit: 0.66–0.75 G$_{load}$", fontsize=10, color=MUS)
b.text(6.3, 0.36, "too large for\none model\nICC-MY", fontsize=10, color=RED, va="bottom")

# --- ICC-MY loss (icc_loss.py; uniform loss, g_c per remaining ICC-MY fixed at the rabbit value) ---
import sys; sys.path.insert(0, "paper_R1"); import icc_loss
lr = [d for d in icc_loss.loss() if d["sp"] == "rabbit"]
c = fig.add_axes([0.64, 0.10, 0.33, 0.25])
c.plot([100 * (1 - d["f"]) for d in lr], [100 * d["smc"] / lr[0]["smc"] for d in lr], "o-", color=MUS, ms=3.5, lw=1.4)
gcv, G0 = lr[0]["gc"], lr[0]["G"]; ff = np.linspace(0.08, 1, 60)  # passive divider, as in Fig. 4C
c.plot(100 * (1 - ff), 100 * ff * (gcv + G0) / (ff * gcv + G0), "--", color=INK, lw=1.0)
c.axhline(50, color=MUTED, lw=0.7, ls=":")
c.set_xlim(-3, 93); c.set_ylim(0, 110); c.set_xticks([0, 30, 60, 90]); c.set_yticks([0, 50, 100]); c.tick_params(labelsize=10)
c.set_xlabel("ICC-MY lost (%)", fontsize=10); c.set_ylabel("muscle wave (%)", fontsize=10)
c.text(1, 4, "model (dots) follows\npassive divider (dashed)", fontsize=10, color=INK, ha="left", va="bottom")

for ext in ("png", "tif"):
    fig.savefig(f"paper_R1/fig/graphical_abstract.{ext}", dpi=300, facecolor="white")
    Image.open(f"paper_R1/fig/graphical_abstract.{ext}").convert("RGB").save(f"paper_R1/fig/graphical_abstract.{ext}", dpi=(300, 300))
fig.savefig("paper_R1/fig/graphical_abstract.pdf", facecolor="white")
print("points:", pts)
print("amplitudes:", amp_i, amp_s)
