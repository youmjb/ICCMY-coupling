"""Anatomy check: g_c giving the Kito CSMC amplitude (24.4 +/- 2.7 mV) relative to the load G_load for the calibrated column
(L 30, m 0.3, g_SS 30; out/v4win), the one-at-a-time m/g_SS changes (out/paper_sens) and short wide columns (L 8, m 3, g_SS 30 and
300; out/paper_anat). G_load = m*sqrt(g_SS*g_m)*tanh(L/lambda), lambda = sqrt(g_SS/g_m) layers, g_m 0.96 nS (sealed finite cable).
Output: paper_R1/fig/anat.png and printed table.
"""
import os, re, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GM, KITO, KSD = 0.96, 24.4, 2.7


def gload(m, gss, L):
    lam = math.sqrt(gss / GM)
    return m * math.sqrt(gss * GM) * math.tanh(L / lam)


def curve(pre, gcs):
    out = []
    for g in gcs:
        f = pre(g) + ".log"
        x = PAT.search(AE.log(f)) if os.path.exists(f) else None
        if x: out.append((g, float(x.group(2)), float(x.group(4))))
    return out


def cross(c, target):
    for (g0, _, y0), (g1, _, y1) in zip(c[:-1], c[1:]):
        if y0 <= target <= y1:
            return math.exp(math.log(g0) + (target - y0) / (y1 - y0) * (math.log(g1) - math.log(g0)))
    return float("nan")


CASES = [("L30 m0.3 gSS30 (model)", 0.3, 30, 30, lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", [0.5, 1, 1.5, 2, 3, 5]),
         ("L30 m0.3 gSS15", 0.3, 15, 30, lambda g: f"out/paper_sens/gSS15/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
         ("L30 m0.3 gSS60", 0.3, 60, 30, lambda g: f"out/paper_sens/gSS60/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
         ("L30 m0.15 gSS30", 0.15, 30, 30, lambda g: f"out/paper_sens/m0.15/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
         ("L30 m0.45 gSS30", 0.45, 30, 30, lambda g: f"out/paper_sens/m0.45/gc{g}_u0", [0.5, 1, 1.5, 2, 3]),
         ("L7 m1.3 gSS30 (mouse-like, same membrane)", 1.3, 30, 7, lambda g: f"out/paper_anat/L7m1.3_gc{g}", [1, 2, 3, 4, 6, 8]),
         ("L8 m3 gSS30", 3, 30, 8, lambda g: f"out/paper_anat/gSS30_gc{g}", [3, 6, 10, 15, 25, 40]),
         ("L8 m3 gSS300", 3, 300, 8, lambda g: f"out/paper_anat/gSS300_gc{g}", [5, 10, 15, 20, 30, 50])]
rows = []
print("case | G_load (nS) | g_c* 21.7/24.4/27.1 | g_c*/G_load | ICC amp at g_c* (nearest) | g_c at SMC/ICC = 0.43 (Kito ratio) | /G_load")
for lab, m, gss, L, pre, gcs in CASES:
    c = curve(pre, gcs); G = gload(m, gss, L)
    gs = [cross(c, KITO + s) for s in (-KSD, 0, KSD)]
    near = min(c, key=lambda x: abs(x[0] - gs[1])) if c and not math.isnan(gs[1]) else (0, float("nan"), 0)
    rat = [(g, a, s1 / a) for g, a, s1 in c]
    gr = cross([(g, 0, r) for g, a, r in rat], 24.4 / 56.9)
    rows.append((lab, G, gs, gr, c))
    print(f"{lab} | {G:.2f} | {' / '.join(f'{x:.2f}' for x in gs)} | {gs[1] / G:.2f} ({gs[0] / G:.2f}-{gs[2] / G:.2f}) | {near[1]:.1f} | {gr:.2f} | {gr / G:.2f}")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
a = ax[0]
GG = np.logspace(-0.3, 1.5, 50)
a.fill_between(GG, 0.43 * GG, 0.66 * GG, color="0.9", label="g$_c$ = 0.43–0.66 G$_{load}$")
for lab, G, gs, gr, c in rows:
    col = "k" if "model" in lab else ("tab:red" if lab.startswith(("L8", "L7")) else "tab:blue")
    a.plot([G], [gr], "o", color=col, ms=6)
    if not math.isnan(gs[1]): a.plot([G], [gs[1]], "x", color=col, ms=6)
    a.annotate(lab, (G, gr), xytext=(5, -9), textcoords="offset points", fontsize=6)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("G$_{load}$ = m·√(g$_{SS}$g$_m$)·tanh(L/λ) (nS)"); a.set_ylabel("g$_c$ (nS)")
a.set_title("A  g$_c$ matching the Kito ratio SMC/ICC 0.43 (●) or SMC 24.4 mV (×)", loc="left", fontsize=9)
a.legend(fontsize=7, loc="upper left")
a = ax[1]
for lab, G, gs, gr, c in rows:
    if not c: continue
    col = "k" if "model" in lab else ("tab:red" if lab.startswith(("L8", "L7")) else None)
    if col is None: continue
    a.plot([x[1] for x in c], [x[2] for x in c], "o-", color=col, ms=3, label=lab)
a.axhspan(KITO - KSD, KITO + KSD, color="0.9", zorder=0); a.axvspan(56.9 - 6.4, 56.9 + 6.4, color="0.93", zorder=0)
a.set_xlabel("ICC-MY amplitude (mV)"); a.set_ylabel("SMC layer-1 amplitude (mV)")
a.set_title("B  Larger columns: muscle never reaches 24 mV before ICC-MY collapses", loc="left", fontsize=9)
a.legend(fontsize=6)
fig.tight_layout(); fig.savefig("paper_R1/fig/anat.png", dpi=150)
print("wrote paper_R1/fig/anat.png")
