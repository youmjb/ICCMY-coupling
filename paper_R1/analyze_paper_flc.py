"""FLC supplement (S3): single-cell behaviour (out/paper_flc, gisim flcsingle) and robustness of the main results to the FLC model
(out/paper_flcmf: minimal mean-field FLC, pdModel=0) vs the modified PDGFa FLC (out/v4win, out/v4ijp).
IJP windows use the 200-ms guarded clean test (as analyze_paper_sens.py). Output: paper_R1/fig/supp_flc.png and printed table.
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
GCS = [0.5, 1, 1.5, 2, 3]
T0 = [20000, 20450, 20900, 21350, 21800, 22250, 22700, 23150]
KITO, KSD = 24.4, 2.7


def load(f, cols):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in cols]


def cross(y, target):
    for k in range(len(y) - 1):
        if y[k] <= target <= y[k + 1]:
            a, b = np.log(GCS[k]), np.log(GCS[k + 1])
            return float(np.exp(a + (target - y[k]) / (y[k + 1] - y[k]) * (b - a)))
    return np.nan


def ijp(pre):
    tc, ic, sc = load(pre + "_ctrl.csv", ("t_ms", "V_ICC", "V_SMC0"))
    rr = []
    for t0 in T0:
        t, i, s = load(f"{pre}_p{t0}.csv", ("t_ms", "V_ICC", "V_SMC0"))
        n = min(len(t), len(tc)); w = (t[:n] >= t0) & (t[:n] < t0 + 1000); g = (t[:n] >= t0) & (t[:n] < t0 + 1200)
        if ic[:n][g].max() < -65 and i[:n][g].max() < -65:
            dS, dI = -(s[:n] - sc[:n])[w].min(), -(i[:n] - ic[:n])[w].min()
            rr.append((dS, dI, dI / dS))
    return rr


SETS = [("modified PDGFa FLC", lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", "out/v4ijp/gc1.5", "k"),
        ("minimal mean-field FLC", lambda g: f"out/paper_flcmf/gc{g}_u0", "out/paper_flcmf/ijp", "tab:purple")]
R = {}
print("FLC model | SMC L1 at g_c 0.5/1/1.5/2/3 | g_c* 21.7/24.4/27.1 | FLC rest in tissue (g_c 1.5) | IJP n, SMC, ICC, ratio")
for lab, pre, ip, c in SETS:
    smc = [float(PAT.search(AE.log(pre(g) + ".log")).group(4)) for g in GCS]
    gs = [cross(smc, KITO + s) for s in (-KSD, 0, KSD)]
    t, vp = load(pre(1.5) + ".csv", ("t_ms", "V_P0")); vrest = np.median(vp[t >= 20000])
    rr = ijp(ip)
    R[lab] = dict(smc=smc, gs=gs, rr=rr, c=c)
    print(f"{lab} | {' '.join(f'{x:.1f}' for x in smc)} | {' / '.join(f'{x:.2f}' for x in gs)} | {vrest:.1f} | "
          f"{len(rr)}, {np.mean([x[0] for x in rr]):.2f}, {np.mean([x[1] for x in rr]):.2f}, {np.mean([x[2] for x in rr]):.3f}")

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5)); ax = ax.ravel()
a = ax[0]
for f, c, lab in (("ctrl", "k", "control"), ("apamin", "tab:blue", "SK blocked (apamin)"), ("agonist", "tab:red", "IP$_3$ agonist from 10 s")):
    t, v = load(f"out/paper_flc/{f}.csv", ("t_ms", "V")); a.plot(t / 1000, v, color=c, lw=0.8, label=lab)
a.set_xlabel("time (s)"); a.set_ylabel("Vm (mV)"); a.legend(fontsize=7)
a.set_title("A  Isolated FLC (no coupling)", loc="left", fontsize=10)
a = ax[1]
t, v = load("out/v4win/m0.3_gSS30_gc1.5_u0.csv", ("t_ms", "V_P0")); w = (t >= 30000) & (t < 38000)
ta, va = load("out/paper_cal/apa_gc1.5.csv", ("t_ms", "V_P0")); wa = (ta >= 30000) & (ta < 38000)
a.plot((t[w] - 30000) / 1000, v[w], color="k", lw=0.9, label="FLC in tissue (layer 1), g$_c$ 1.5")
a.plot((ta[wa] - 30000) / 1000, va[wa], color="tab:blue", lw=0.9, label="same, apamin")
a.axhspan(-49.5 - 5.3, -49.5 + 5.3, color="0.9", zorder=0); a.text(0.1, -43.5, "Kito FLC RMP", fontsize=7)
a.set_xlabel("time (s)"); a.set_ylabel("Vm (mV)"); a.legend(fontsize=7, loc="lower right")
a.set_title("B  FLC coupled to the SMC syncytium", loc="left", fontsize=10)
a = ax[2]
for lab, pre, ip, c in SETS:
    a.plot(GCS, R[lab]["smc"], "o-", color=c, ms=4, label=lab)
a.axhspan(KITO - KSD, KITO + KSD, color="0.9", zorder=0)
a.set_xscale("log"); a.set_xticks(GCS); a.set_xticklabels([f"{g:g}" for g in GCS]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("g$_c$ (nS)"); a.set_ylabel("SMC layer-1 slow wave (mV)"); a.legend(fontsize=7)
a.set_title("C  Amplitude vs g$_c$ with two FLC models", loc="left", fontsize=10)
a = ax[3]
for k, (lab, pre, ip, c) in enumerate(SETS):
    q = [x[2] for x in R[lab]["rr"]]
    a.plot([k] * len(q), q, "o", color=c, alpha=0.4); a.plot([k], [np.mean(q)], "o", color=c, ms=9)
    a.text(k + 0.08, np.mean(q), f"{np.mean(q):.2f} (n {len(q)}); SMC IJP {np.mean([x[0] for x in R[lab]['rr']]):.1f} mV", fontsize=7, va="center")
a.set_xticks([0, 1]); a.set_xticklabels([s[0] for s in SETS], fontsize=8); a.set_xlim(-0.4, 1.9); a.set_ylim(0, 0.6)
a.set_ylabel("ICC-MY hyperpol. / SMC layer-1 IJP (g$_c$ 1.5)")
a.set_title("D  IJP fraction with two FLC models", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/supp_flc.png", dpi=150)
print("wrote paper_R1/fig/supp_flc.png")
