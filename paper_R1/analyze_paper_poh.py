"""Step 3: g_c estimate with the Poh et al. 2012 human jejunal SMC model in place of SmoothM (run_paper_poh.ps1).
Compares SMC layer-1 amplitude vs g_c, the g_c giving the Kito CSMC amplitude (24.4 +/- 2.7 mV) and the SMC/ICC-MY ratio at g_c 1.5,
for SmoothM (out/v4win) and Poh (tonic K+ 0 and 0.4 nS). Output: paper_R1/fig/poh.png and printed table.
"""
import os, re, csv
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GCS = [0.5, 1, 1.5, 2, 3, 5]
KITO, KSD = 24.4, 2.7
SETS = [("SmoothM (tonic K$^+$ 0.4 nS)", lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", "k"),
        ("Poh 2012 (no tonic K$^+$)", lambda g: f"out/paper_poh/gK0_gc{g}", "tab:red"),
        ("Poh 2012 (tonic K$^+$ 0.4 nS)", lambda g: f"out/paper_poh/gK0.4_gc{g}", "tab:orange")]


def logf(pre):
    return tuple(float(x) for x in PAT.search(AE.log(pre + ".log")).groups())


def load(f):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0")]


def cross(y, target):
    for k in range(len(y) - 1):
        if y[k] <= target <= y[k + 1]:
            a, b = np.log(GCS[k]), np.log(GCS[k + 1])
            return float(np.exp(a + (target - y[k]) / (y[k + 1] - y[k]) * (b - a)))
    return np.nan


R = {}
print("model | uncoupled SMC rest | SMC L1 at g_c 0.5/1/1.5/2/3/5 | ICC amp at 1.5 | g_c* 21.7/24.4/27.1 | SMC/ICC at 1.5 | SMC rest at 1.5")
for lab, pre, c in SETS:
    d = [logf(pre(g)) for g in GCS]
    smc = [x[3] for x in d]
    gs = [cross(smc, KITO + s) for s in (-KSD, 0, KSD)]
    t0, i0, s0 = load(pre(0) + ".csv"); r0 = np.median(s0[t0 >= 20000])
    t, i, s = load(pre(1.5) + ".csv"); h = t >= 20000
    bw = ~(np.convolve((i[h] >= -65).astype(float), np.ones(601), mode="same") > 0)
    r15 = np.median(s[h][bw])
    R[lab] = dict(smc=smc, gs=gs, c=c)
    print(f"{lab} | {r0:.1f} | {' '.join(f'{x:.1f}' for x in smc)} | {d[2][1]:.1f} | {' / '.join(f'{x:.2f}' for x in gs)} | {d[2][3] / d[2][1]:.3f} | {r15:.1f}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(width_ratios=[1.3, 1, 0.9]))
a = ax[0]
for lab, pre, c in SETS:
    t, i, s = load(pre(1.5) + ".csv")
    k = np.where((i[:-1] < -45) & (i[1:] >= -45))[0]; k = k[t[k] >= 25000][0]
    sl = slice(k - 1000, k + 4000); x = (t[sl] - t[k]) / 1000
    a.plot(x, i[sl], color=c, lw=0.8, ls="--"); a.plot(x, s[sl], color=c, lw=1.7, label=lab)
a.set_xlabel("time from ICC-MY upstroke (s), g$_c$ 1.5 nS"); a.set_ylabel("Vm (mV)")
a.set_title("A  ICC-MY (dashed) and SMC layer 1 (solid)", loc="left", fontsize=10); a.legend(fontsize=7)
a = ax[1]
for lab, pre, c in SETS:
    a.plot(GCS, R[lab]["smc"], "o-", color=c, ms=4, label=lab)
a.axhspan(KITO - KSD, KITO + KSD, color="0.9", zorder=0)
a.set_xscale("log"); a.set_xticks([0.5, 1, 2, 3, 5]); a.set_xticklabels(["0.5", "1", "2", "3", "5"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("SMC layer-1 slow wave (mV)")
a.set_title("B  Amplitude vs g$_c$ (band: Kito ± SD)", loc="left", fontsize=10); a.legend(fontsize=7)
a = ax[2]
y = np.arange(len(SETS))[::-1]
for yy, (lab, pre, c) in zip(y, SETS):
    lo, mid, hi = R[lab]["gs"]
    a.plot([lo, hi], [yy, yy], color=c, lw=5, alpha=0.35, solid_capstyle="butt"); a.plot([mid], [yy], "o", color=c)
a.axvspan(0.7, 1.6, color="0.9", zorder=0)
a.set_yticks(y); a.set_yticklabels([s[0] for s in SETS], fontsize=7)
a.set_xscale("log"); a.set_xticks([0.5, 0.7, 1, 1.6, 2, 3]); a.set_xticklabels(["0.5", "0.7", "1", "1.6", "2", "3"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.set_xlim(0.45, 3.5)
a.set_xlabel("g$_c$ giving SMC layer 1 = 24.4 mV (shade: paper's range)")
a.set_title("C  Estimate of g$_c$", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/poh.png", dpi=150)
print("wrote paper_R1/fig/poh.png")
