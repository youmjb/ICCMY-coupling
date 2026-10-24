"""Nifedipine condition (paper_R1/run_paper_nif.ps1): SmoothM L-type Ca2+ scaled by SMS_CAL 1 (control = out/v4win, out/v4ijp),
0.1 and 0 (out/paper_nif/cal<level>/). Same metrics as analyze_paper_sens.py (g_c giving SMC layer-1 = Kito 24.4 +/- 2.7 mV;
IJP fraction at g_c 1.5 with the 200-ms guarded clean-window test). Output: paper_R1/fig/nif.png, printed table.
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
GCS = [0, 0.5, 1, 1.5, 2, 3, 5]
T0 = [20000, 20450, 20900, 21350, 21800, 22250, 22700, 23150]
KITO, KSD = 24.4, 2.7
LEV = [("1", "control (SMS_CAL 1)", "k"), ("0.1", "90 % L-type block", "tab:orange"), ("0", "full L-type block", "tab:red")]


def path(lev, name):
    if lev == "1":
        return f"out/v4win/m0.3_gSS30_{name}" if name.startswith("gc") else f"out/v4ijp/gc1.5_{name[4:]}"
    return f"out/paper_nif/cal{lev}/{name}"


def logf(lev, gc):
    m = PAT.search(AE.log(path(lev, f"gc{gc}_u0") + ".log"))
    return tuple(float(x) for x in m.groups())


def load(f, cols=("t_ms", "V_ICC", "V_SMC0")):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in cols]


def rest(t, v):
    return np.percentile(v[t >= t[-1] / 2], 5)


def cross(gcs, y, target):
    for k in range(len(y) - 1):
        if y[k] <= target <= y[k + 1]:
            a, b = np.log(gcs[k]), np.log(gcs[k + 1])
            return float(np.exp(a + (target - y[k]) / (y[k + 1] - y[k]) * (b - a)))
    return np.nan


t, ib, _ = load(path("1", "ijp_ctrl") + ".csv"); RB = rest(t, ib)
R = {}
print("level | uncoupled ICC amp/cpm, SMC amp | SMC L1 at g_c 0.5/1/1.5/2/3/5 | ICC amp at 1.5 | g_c* 21.7/24.4/27.1 | SMC rest | IJP n, SMC, ICC, ratio")
for lev, lab, _ in LEV:
    d = [logf(lev, g) for g in GCS]
    smc = [x[3] for x in d[1:]]
    gs = [cross(GCS[1:], smc, KITO + s) for s in (-KSD, 0, KSD)]
    tc, ic, sc = load(path(lev, "ijp_ctrl") + ".csv")
    rc = rest(tc, ic); thr = -65.0 + (rc - RB)
    srest = np.percentile(sc[tc >= tc[-1] / 2], 5)
    rr = []
    for t0 in T0:
        tt, i, s = load(path(lev, f"ijp_p{t0}") + ".csv")
        n = min(len(tt), len(tc)); w = (tt[:n] >= t0) & (tt[:n] < t0 + 1000); g = (tt[:n] >= t0) & (tt[:n] < t0 + 1200)
        if ic[:n][g].max() < thr and i[:n][g].max() < thr:
            dS, dI = -(s[:n] - sc[:n])[w].min(), -(i[:n] - ic[:n])[w].min()
            rr.append((dS, dI, dI / dS))
    R[lev] = dict(d=d, smc=smc, gs=gs, rr=rr)
    ij = f"{len(rr)}, {np.mean([x[0] for x in rr]):.2f}, {np.mean([x[1] for x in rr]):.2f}, {np.mean([x[2] for x in rr]):.3f}" if rr else "0"
    print(f"{lab:20s} | {d[0][1]:.1f}/{d[0][0]:.2f}, {d[0][3]:.1f} | {' '.join(f'{x:.1f}' for x in smc)} | {d[3][1]:.1f} | "
          f"{' / '.join(f'{x:.2f}' for x in gs)} | {srest:.1f} | {ij}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(width_ratios=[1.4, 1, 1]))
a = ax[0]
for k, (lev, lab, c) in enumerate(LEV):
    tt, i, s = load(path(lev, "gc1.5_u0") + ".csv")
    w = (tt >= 30000) & (tt < 37000)
    a.plot((tt[w] - 30000) / 1000, i[w], color=c, lw=0.9, ls="--")
    a.plot((tt[w] - 30000) / 1000, s[w], color=c, lw=1.6, label=lab)
a.set_xlabel("time (s), g$_c$ 1.5 nS"); a.set_ylabel("Vm (mV)")
a.set_title("A  SMC layer 1 (solid) and ICC-MY (dashed)", loc="left", fontsize=10); a.legend(fontsize=7)
a = ax[1]
for lev, lab, c in LEV:
    a.plot(GCS[1:], R[lev]["smc"], "o-", color=c, label=lab)
a.axhspan(KITO - KSD, KITO + KSD, color="0.9", zorder=0); a.text(0.52, 25, "Kito CSMC (in nifedipine)", fontsize=7)
a.set_xscale("log"); a.set_xticks([0.5, 1, 2, 3, 5]); a.set_xticklabels(["0.5", "1", "2", "3", "5"])
a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("SMC layer-1 slow wave (mV)")
a.set_title("B  Amplitude vs g$_c$", loc="left", fontsize=10); a.legend(fontsize=7)
a = ax[2]
for k, (lev, lab, c) in enumerate(LEV):
    q = [x[2] for x in R[lev]["rr"]]
    if q:
        a.plot([k] * len(q), q, "o", color=c, alpha=0.4); a.plot([k], [np.mean(q)], "o", color=c, ms=9)
        a.text(k + 0.1, np.mean(q), f"{np.mean(q):.2f} (n {len(q)})", fontsize=7, va="center")
a.set_xticks(range(3)); a.set_xticklabels(["control", "90 %", "full"]); a.set_xlim(-0.5, 2.8); a.set_ylim(0, 0.5)
a.set_xlabel("L-type block in SMC"); a.set_ylabel("ICC-MY hyperpol. / SMC layer-1 IJP (g$_c$ 1.5)")
a.set_title("C  IJP fraction", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/nif.png", dpi=150)
print("wrote paper_R1/fig/nif.png")
