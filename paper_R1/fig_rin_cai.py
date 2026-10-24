"""Passive in-situ R_in of ICC-MY vs g_c, compared with Cai et al. 2017 (mean +/- SEM, n = 13).

R_in = 1000 / (g_own + g_c*G/(g_c+G)), G = m*G1 (same analytic estimate as R1 Fig. 2B).
G1 = 5.37 nS per column (out/iccfit/A_*.csv, m 1); g_own = 5 nS x 0.72 = 3.6 nS (assumed).
Cai 2017 (J Neurogastroenterol Motil 23:117, checked against the PDF): R_in 260 +/- 38 MOhm,
means +/- SEM, n = 13 cells, room temperature (~23 C), mice 9-15 days old, LMMP with circular
muscle locally removed. SD = SEM*sqrt(13) = 137 MOhm.
Output: paper_R1/fig/rin_cai.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams["font.family"] = "Arial"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "fig")
os.makedirs(OUT, exist_ok=True)

G1, GOWN = 5.37, 3.6
MEAN, SEM, N = 260.0, 38.0, 13
SD = SEM * np.sqrt(N)
T95 = stats.t.ppf(0.975, N - 1)

gc = np.logspace(-1, np.log10(50), 300)
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4))

a.axhspan(MEAN - SD, MEAN + SD, color="0.92", label="Cai mean ± SD")
a.axhspan(MEAN - T95 * SEM, MEAN + T95 * SEM, color="0.80", label="Cai 95% CI of mean")
a.axhspan(MEAN - SEM, MEAN + SEM, color="0.65", label="Cai mean ± SEM")
a.axhline(MEAN, color="k", lw=0.8)
a.axhline(1000 / GOWN, color="k", ls=":", lw=0.8)
a.text(0.11, 1000 / GOWN + 4, "isolated ICC-MY (278 MΩ)", fontsize=8)
for m, c in [(0.1, "tab:green"), (0.3, "tab:blue"), (1, "tab:orange"), (3, "tab:red")]:
    G = m * G1
    a.plot(gc, 1000 / (GOWN + gc * G / (gc + G)), color=c, label=f"m = {m}")
a.axvspan(1, 2, color="tab:blue", alpha=0.08)
a.set_xscale("log")
a.set_ylim(80, 420)
a.set_xlabel("g$_c$ (nS)")
a.set_ylabel("ICC-MY input resistance (MΩ)")
a.set_title("A  Passive in-situ R$_{in}$ vs Cai 2017", loc="left", fontsize=10)
a.legend(fontsize=7, loc="lower left", ncol=2)

# B: what the measurement allows, expressed as the effective load g_c*G/(g_c+G)
labels = ["mean", "mean − SEM", "95% CI\nlower", "mean − SD"]
rins = [MEAN, MEAN - SEM, MEAN - T95 * SEM, MEAN - SD]
loads = [1000 / r - GOWN for r in rins]
b.bar(range(4), loads, color=["0.3", "0.5", "0.65", "0.8"])
for i, (r, L) in enumerate(zip(rins, loads)):
    b.text(i, L + 0.08, f"{L:.2f} nS\n(R$_{{in}}$ {r:.0f} MΩ)", ha="center", fontsize=8)
for m, c in [(0.3, "tab:blue"), (1, "tab:orange")]:
    G = m * G1
    for g, ls in [(1, "-"), (2, "--")]:
        v = g * G / (g + G)
        b.axhline(v, color=c, ls=ls, lw=0.9, label=f"m {m}, g$_c$ {g}: {v:.2f} nS")
b.set_xticks(range(4))
b.set_xticklabels(labels, fontsize=8)
b.set_ylabel("max. effective load g$_c$G/(g$_c$+G) (nS)")
b.set_ylim(0, 5.3)
b.set_title("B  Upper bound on the load seen by ICC-MY", loc="left", fontsize=10)
b.legend(fontsize=7, loc="upper left")

fig.tight_layout()
fig.savefig(os.path.join(OUT, "rin_cai.png"), dpi=150)
print("SD", round(SD, 1), "t95", round(T95, 3), "CI", round(MEAN - T95 * SEM, 1), round(MEAN + T95 * SEM, 1))
for l, r, L in zip(labels, rins, loads):
    print(l.replace("\n", " "), round(r, 1), "MOhm -> load", round(L, 2), "nS")
