"""Plan-A Fig. S3: the revision-4 Fig. 1C-D (model traces and calibration z-scores) as a standalone figure. Same data and code as
figs_assemble.fig1 (original ICC-MY model, out/v4win and out/paper_cal). Output: fig/A_figS3 (.pdf/.tif/.png)."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import figs_assemble as FA
import figstyle
cal = FA.quiet_import("analyze_paper_cal")
fig = plt.figure(figsize=(figstyle.WIDTH, 3.6))
t, vi, vs, vp = cal.load("out/v4win/m0.3_gSS30_gc1.5_u0.csv")
ta, _, _, vpa = cal.load("out/paper_cal/apa_gc1.5.csv")
w = (t >= 30000) & (t < 38000); wa = (ta >= 30000) & (ta < 38000); x = (t[w] - 30000) / 1000
for r, (v, lab) in enumerate(((vi, "ICC-MY"), (vs, "muscle layer 1"), (vp, "FLC, layer 1"))):
    a = fig.add_axes([0.075, 0.70 - r * 0.29, 0.40, 0.23])
    a.plot(x, v[w], color="k", lw=0.9, label=lab)
    if r == 2:
        a.plot(x, vpa[wa], color="0.55", lw=0.8, label="FLC, SK blocked"); a.set_xlabel("time (s)")
    else:
        a.set_xticklabels([])
    a.set_xlim(0, 8); a.set_ylabel("V$_m$ (mV)")
    a.legend(frameon=False, loc="lower right" if r == 2 else "upper center", ncol=1 if r == 2 else 2, handlelength=1.2, borderaxespad=0)
    if r == 0: a.set_title("A", loc="left", fontweight="bold")
a = fig.add_axes([0.70, 0.12, 0.29, 0.83])
keys = list(cal.KITO.keys()); y = np.arange(len(keys))[::-1]
a.axvspan(-1, 1, color=FA.GREY, lw=0); a.axvline(0, color="0.5", lw=0.6)
for g, mk, c, fc in zip(cal.GCS, ("o", "s", "^"), ("0.45", "k", "k"), ("0.45", "k", "white")):
    z = np.array([(cal.M[g][k] - cal.KITO[k][0]) / cal.KITO[k][1] for k in keys])
    a.plot(np.clip(z, -8, 8), y, mk, color=c, mfc=fc, mec=c, ms=3.2, label=f"rabbit [3], g$_c$ {g} nS", mew=0.7)
MOUSE = {"ICC-MY RMP (mV)": (-69.6, 4.5), "ICC-MY amplitude (mV)": (56.5, 5.3), "ICC-MY rate (min$^{-1}$)": (26.1, 2.4),
         "ICC-MY half-width (s)": (1.00, 0.16), "ICC-MY dV/dt max (V s$^{-1}$)": (1.51, 0.30), "CSMC RMP (mV)": (-69.9, 3.9),
         "CSMC amplitude (mV)": (34.1, 4.6), "CSMC half-width (s)": (1.07, 0.20), "CSMC dV/dt max (V s$^{-1}$)": (0.54, 0.10)}
zm = [((cal.M[2][k] - MOUSE[k][0]) / MOUSE[k][1], yy) for k, yy in zip(keys, y) if k in MOUSE]
a.plot([np.clip(z, -8, 8) for z, _ in zm], [yy for _, yy in zm], "D", mfc="white", mec="k", mew=0.9, ms=3.4, label="mouse [4], g$_c$ 2 nS")
a.set_yticks(y); a.set_yticklabels([k.replace("CSMC", "muscle") for k in keys], fontsize=6.5)
a.set_xlim(-8.5, 8.5); a.set_xticks([-8, -4, -1, 0, 1, 4, 8]); a.set_xlabel("(model − recorded mean) / recorded SD")
a.legend(frameon=True, framealpha=1, edgecolor="none", loc="lower left", handletextpad=0.2)
a.set_title("B", loc="left", fontweight="bold")
figstyle.save(fig, "paper_R1/fig/A_figS3"); plt.close(fig)
print("wrote paper_R1/fig/A_figS3")
