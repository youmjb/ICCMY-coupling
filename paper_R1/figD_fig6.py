"""Model-D Figure 6: how tight the ICC-MY membrane must be, and what isolated-cell recordings measure.
A  Transfer of the muscle resting change to ICC-MY at the pinacidil level against the input resistance of the isolated published 25-pF
   cell measured as Goto et al. 2004 did (holding -80 mV, -20-mV step; bar to the +10-mV step; goto_clamp.py). Passive muscle, g_c 1.75 nS (open,
   lowg_coupled/fig_lowg.sized, resting currents x1, 0.3, 0.2, 0.1); active muscle (filled): model D (x0.3) and model C (x0.1) with g_c set
   by the control ratio (final_model_D/C.out.txt). Bands: tau of freshly prepared ICC (Goto et al. 2004: 25.2 pF x 0.76 +/- 0.41 GOhm) and
   of cultured ICC (Koh et al. 1998: 29 pF x 250 MOhm); dashed line, recorded transfer (0.93).
B  Time constant measured in clusters of N = 1, 3, 6 ICC coupled all-to-all by 3 nS (cluster_pulse.py) for resting currents x1, 0.3,
   0.2 and 0.1: it does not depend on the number of coupled cells.
Output: fig/D_fig6 (.pdf/.tif/.png), figD_fig6.out.txt."""
import os, sys, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import figstyle; figstyle.apply()
import cluster_pulse as CP
from fig_lowg import sized
OUT = []


def log(s): print(s); OUT.append(s)


def tau(f, n="1", topo="0", g="1"): return CP.measure((f, n, topo, g))[1]


PAT = r"K05 pinacidil\s+g_c fit ([\d.]+).*?transfer ([\d.]+) \(rec.*?ICC x([\d.]+) \(rec"
pin = {m: [float(x) for x in re.search(PAT, open(f"paper_R1/final_model_{m}.out.txt", encoding="utf-8").read()).groups()] for m in ("C", "D")}
TG = 25.2 * 0.76; TGlo, TGhi = 25.2 * (0.76 - 0.41), 25.2 * (0.76 + 0.41)     # Goto et al. 2004 (ms)
TK = 29 * 0.25                                                                 # Koh et al. 1998 (ms)

GT = {}                                                                         # Goto-protocol R_in of the published 25-pF cell (goto_clamp.py)
for line in open("paper_R1/goto_clamp.out.txt", encoding="utf-8"):
    m = re.match(r"\s+pub f\s+([\d.]+): hold .*?\|\s+(-?\d+) /\s+(-?\d+) /\s+(-?\d+) MOhm", line)
    if m: GT[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
fig, (a, b) = plt.subplots(1, 2, figsize=(figstyle.WIDTH, 2.9), gridspec_kw={"width_ratios": [1.25, 1]})
a.axvspan(760 - 410, 760 + 410, color="0.88", lw=0); a.axvline(760, color="0.5", lw=0.8)
a.text(800, 1.06, "freshly prepared ICC", fontsize=figstyle.SMALL, color="0.3", va="top")
a.axhline(0.93, color="k", lw=0.8, ls="--"); a.text(9000, 0.895, "recorded (pinacidil)", fontsize=figstyle.SMALL, ha="right")
fs = ["1", "0.3", "0.2", "0.1"]
xs = [GT[f][0] for f in fs]; xe = [GT[f][1] for f in fs]; ys = [sized(f)[0][0] for f in fs]
for x, x2, y in zip(xs, xe, ys):
    a.plot([x, x2], [y, y], color="0.55", lw=2.5, solid_capstyle="butt", zorder=1)
a.plot(xs, ys, "o-", color="k", mfc="white", ms=4.5, lw=0.9, label="passive muscle", zorder=3)
for f, x, y in zip(fs, xs, ys):
    a.annotate(f"×{f}", (x, y), xytext=(-4, 6), textcoords="offset points", fontsize=figstyle.SMALL, ha="right")
    log(f"A resting currents x{f}: R_in (-80 mV hold) -20 mV step {GT[f][0]:.0f}, +10 mV step {GT[f][1]:.0f} MOhm; transfer (passive muscle, g_c 1.75) {y:.2f}")
for m, f in (("D", "0.3"), ("C", "0.1")):
    a.plot(GT[f][0], pin[m][1], "P", color="k", ms=8, mfc="k" if m == "D" else "white", label=f"model {m} (active muscle)", zorder=4)
    log(f"A model {m}: R_in -20 mV step {GT[f][0]:.0f} MOhm, transfer {pin[m][1]:.2f}, ICC x{pin[m][2]:.2f} (g_c {pin[m][0]:.2f})")
a.plot([], [], color="0.55", lw=2.5, label="range to the +10-mV step")
log("A Goto 2004 R_in 760 +/- 410 MOhm (25.2 pF)")
a.set_xscale("log"); a.set_xlim(50, 10000); a.set_ylim(0.2, 1.1)
a.set_xlabel("input resistance of the isolated 25-pF cell,\nheld at −80 mV, −20-mV step (MΩ)"); a.set_ylabel("transfer of the muscle resting change")
a.legend(frameon=False, fontsize=figstyle.SMALL, loc="lower right", bbox_to_anchor=(1.0, 0.08))
a.text(-0.17, 1.02, "A", transform=a.transAxes, fontsize=10, fontweight="bold")
for f, mk in zip(fs, ("o", "s", "^", "D")):
    t = [tau(f, n, "1" if n != "1" else "0", "3" if n != "1" else "1") for n in ("1", "3", "6")]
    b.plot([1, 3, 6], t, "-", marker=mk, color="k", mfc="white", ms=4, lw=0.8, label=f"×{f}")
    log(f"B x{f}: tau N1/N3/N6 (all-to-all 3 nS) {t[0]:.0f}/{t[1]:.0f}/{t[2]:.0f} ms")
b.axhspan(TGlo, TGhi, color="0.88", lw=0); b.axhline(TG, color="0.5", lw=0.8); b.axhline(TK, color="0.5", lw=0.8, ls=":")
b.set_yscale("log"); b.set_ylim(2, 3000); b.set_xticks([1, 3, 6]); b.set_xlim(0.5, 6.5)
b.set_xlabel("number of ICC in the cluster"); b.set_ylabel("measured time constant (ms)")
b.legend(frameon=False, fontsize=figstyle.SMALL, loc="upper center", ncol=4, columnspacing=0.8, handlelength=1.5, title="resting conductances", title_fontsize=figstyle.SMALL)
b.text(-0.22, 1.02, "B", transform=b.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout(); figstyle.save(fig, "paper_R1/fig/D_fig6"); plt.close(fig)
open("paper_R1/figD_fig6.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
