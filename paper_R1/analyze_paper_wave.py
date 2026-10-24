"""Plateau-duration test (paper_R1/run_paper_wave.ps1): prescribed ICC-MY trapezoids W1-W7 vs the model ICC-MY (out/v4win).
Per waveform: SMC layer-1 amplitude vs g_c, g_c giving 24.4 mV (Kito CSMC, band +/- SD 2.7), SMC/ICC amplitude ratio,
SMC layer-1 half-width at g_c 1.5. Output: paper_R1/fig/wave.png and printed table.
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
WAVES = [("model", "model ICC-MY (3.3 s, HW 0.66 s)", "k", "-"),
         ("W1", "W1 3.33 s, HW 0.67, rise 60 ms", "0.5", "--"), ("W2", "W2 3.33 s, HW 1.5", "tab:purple", "--"),
         ("W3", "W3 6 s, HW 0.67", "tab:blue", "-"), ("W4", "W4 6 s, HW 1.5", "tab:green", "-"),
         ("W5", "W5 6 s, HW 2.4", "tab:orange", "-"), ("W6", "W6 = W5, rise 5 ms", "tab:red", "-"),
         ("W7", "W7 = W6, -65 mV, 57 mV (Kito-like)", "tab:brown", "-")]


def pre(w, gc):
    return f"out/v4win/m0.3_gSS30_gc{gc}_u0" if w == "model" else f"out/paper_wave/{w}_gc{gc}"


def logf(w, gc):
    return tuple(float(x) for x in PAT.search(AE.log(pre(w, gc) + ".log")).groups())


def load(f):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0")]


def halfwidth(t, vi, v):
    k = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]
    k = k[t[k] >= 20000]
    out = []
    for a, b in zip(k[:-1], k[1:]):
        if a < 500: continue
        seg = v[a - 500:b - 500]; half = seg.min() + (seg.max() - seg.min()) / 2
        out.append((seg > half).sum() / 1000.0)
    return np.mean(out)


def cross(y, target):
    for k in range(len(y) - 1):
        if y[k] <= target <= y[k + 1]:
            a, b = np.log(GCS[k]), np.log(GCS[k + 1])
            return float(np.exp(a + (target - y[k]) / (y[k + 1] - y[k]) * (b - a)))
    return np.nan


reg = open("out/paper_wave/regress_gc1.5.csv").read() == open("out/v4win/m0.3_gSS30_gc1.5_u0.csv").read()
print("regression (wavePer 0 vs v4win g_c 1.5) identical:", reg)
R = {}
print("wave | ICC amp | SMC L1 at g_c 0.5/1/1.5/2/3/5 | g_c* 21.7/24.4/27.1 | SMC/ICC at 1.5 | SMC HW at 1.5 (s) | cpm at 1.5")
for w, lab, _, _ in WAVES:
    d = [logf(w, g) for g in GCS]
    smc = [x[3] for x in d]
    gs = [cross(smc, KITO + s) for s in (-KSD, 0, KSD)]
    t, vi, vs = load(pre(w, 1.5) + ".csv")
    hw = halfwidth(t, vi, vs)
    R[w] = dict(smc=smc, gs=gs, ampI=d[2][1], hw=hw, ratio=[x[3] / x[1] for x in d])
    print(f"{w:5s} | {d[2][1]:5.1f} | {' '.join(f'{x:5.1f}' for x in smc)} | {' / '.join(f'{x:.2f}' for x in gs)} | "
          f"{d[2][3] / d[2][1]:.3f} | {hw:.2f} | {d[2][0]:.2f}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
a = ax[0]
for w, c, off in (("model", "k", 0), ("W1", "0.5", 0), ("W5", "tab:orange", 0), ("W7", "tab:brown", 0)):
    t, vi, vs = load(pre(w, 1.5) + ".csv")
    k = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; k = k[t[k] >= 25000][0]
    s = slice(k - 1000, k + 5000); x = (t[s] - t[k]) / 1000
    a.plot(x, vi[s], color=c, lw=0.9, ls="--"); a.plot(x, vs[s], color=c, lw=1.8, label=dict((q[0], q[1]) for q in WAVES)[w])
a.set_xlabel("time from ICC-MY upstroke (s), g$_c$ 1.5 nS"); a.set_ylabel("Vm (mV)")
a.set_title("A  ICC-MY (dashed) and SMC layer 1 (solid)", loc="left", fontsize=10); a.legend(fontsize=6, loc="upper right")
a = ax[1]
for w, lab, c, ls in WAVES:
    a.plot(GCS, R[w]["smc"], marker="o", ms=3, color=c, ls=ls, label=w if w != "model" else "model")
a.axhspan(KITO - KSD, KITO + KSD, color="0.9", zorder=0)
a.set_xscale("log"); a.set_xticks([0.5, 1, 2, 3, 5]); a.set_xticklabels(["0.5", "1", "2", "3", "5"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("SMC layer-1 slow wave (mV)")
a.set_title("B  Amplitude vs g$_c$ (band: Kito ± SD)", loc="left", fontsize=10); a.legend(fontsize=6, ncol=2)
a = ax[2]
names = [q[0] for q in WAVES]; y = np.arange(len(names))[::-1]
for yy, w in zip(y, names):
    lo, mid, hi = R[w]["gs"]
    a.plot([lo, hi], [yy, yy], color="tab:blue", lw=4, alpha=0.35, solid_capstyle="butt"); a.plot([mid], [yy], "o", color="tab:blue" if w != "model" else "k")
a.axvspan(1, 2, color="0.9", zorder=0)
a.set_yticks(y); a.set_yticklabels([q[1] for q in WAVES], fontsize=7)
a.set_xscale("log"); a.set_xticks([0.5, 1, 1.5, 2, 3]); a.set_xticklabels(["0.5", "1", "1.5", "2", "3"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.set_xlim(0.45, 3.5)
a.set_xlabel("g$_c$ giving SMC layer 1 = 24.4 mV")
a.set_title("C  Coupling fixed by the amplitude", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/wave.png", dpi=150)
print("wrote paper_R1/fig/wave.png")
