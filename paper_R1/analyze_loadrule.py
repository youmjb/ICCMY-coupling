"""Load rule checks from existing out/v4win runs (no new runs) + direct column conductance at gSS 10/100 (out/paper_g1).
1. sqrt(gSS) scaling: G1(gSS) = 50 pA / dV of V_SMC0 (20.8-21.0 s, inj - ctrl), m 1, g_c 0, as R1 Fig. 2B at gSS 30
   (out/iccfit/A_*.csv); compared with sqrt(gSS * 0.96 nS).
2. Collapse: best protected SMC layer-1 amplitude vs G_load = m*sqrt(gSS*g_m) for gSS 10/30/100 (same criteria as plot_v4win.py).
3. Criteria sensitivity: largest G_load with a window, for ICC amplitude >= 70/80/90 % of uncoupled, rate within 10/15/20 %,
   SMC layer 1 >= 15/20/25 mV (layer-0 locking within 5 % kept).
Output: paper_R1/fig/loadrule.png
"""
import os, re, csv, itertools
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GM = 0.96
MS, GSS, GCS = [0.1, 0.3, 1, 3, 10], [10, 30, 100], [0.5, 1, 2, 3, 5, 10, 20, 50]


def v4(m, g, gc):
    f = f"out/v4win/m{m}_gSS{g}_gc{gc}_u0.log"
    if not os.path.exists(f): return None
    x = PAT.search(AE.log(f))
    return tuple(float(v) for v in x.groups()) if x else None


def col(f, k="V_SMC0"):
    r = list(csv.DictReader(open(f)))
    return np.array([float(x["t_ms"]) for x in r]), np.array([float(x[k]) for x in r])


print("1. column conductance G1 (m 1):")
G1 = {}
for g, pre in ((10, "out/paper_g1/gSS10"), (30, "out/iccfit/A"), (100, "out/paper_g1/gSS100")):
    if not os.path.exists(f"{pre}_inj.csv"): print("   missing", pre); continue
    t, c = col(f"{pre}_ctrl.csv"); t2, i = col(f"{pre}_inj.csv"); n = min(len(t), len(t2))
    w = (t[:n] >= 20800) & (t[:n] < 21000)
    G1[g] = 50 / abs((i[:n] - c[:n])[w].mean())
    print(f"   gSS {g}: G1 {G1[g]:.2f} nS, sqrt(gSS*{GM}) {np.sqrt(g * GM):.2f} nS")

ref = v4(0.3, 30, 0)   # uncoupled ICC-MY: cpm, amp, L0 cpm, L1


def window(a, amp=0.8, rate=0.15, l1=20):
    return a[1] >= amp * ref[1] and abs(a[0] / ref[0] - 1) <= rate and abs(a[2] / a[0] - 1) <= 0.05 and a[3] >= l1


print("\n2. best protected SMC L1 amplitude (L1 criterion dropped) vs G_load")
pts = []
for m, g in itertools.product(MS, GSS):
    rr = [(gc, v4(m, g, gc)) for gc in GCS]; rr = [x for x in rr if x[1]]
    ok = [a[3] for gc, a in rr if window(a, l1=-1)]
    pts.append((m, g, m * np.sqrt(g * GM), max(ok) if ok else np.nan))
    print(f"   m {m:4} gSS {g:3}: G_load {pts[-1][2]:6.2f} nS, best L1 {pts[-1][3]:.1f}")

print("\n3. largest G_load with a window (and smallest without)")
crit = []
for amp, rate, l1 in itertools.product([0.7, 0.8, 0.9], [0.10, 0.15, 0.20], [15, 20, 25]):
    has = []
    for m, g in itertools.product(MS, GSS):
        rr = [v4(m, g, gc) for gc in GCS]
        has.append((m * np.sqrt(g * GM), any(a and window(a, amp, rate, l1) for a in rr)))
    gin = max([G for G, h in has if h], default=np.nan); gout = min([G for G, h in has if not h], default=np.nan)
    crit.append((amp, rate, l1, gin, gout))
    print(f"   ICC >= {amp:.0%}, rate +/-{rate:.0%}, L1 >= {l1} mV: largest with window {gin:.2f}, smallest without {gout:.2f} nS")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
a = ax[0]
for g, c in zip(GSS, ("tab:green", "tab:blue", "tab:red")):
    p = [x for x in pts if x[1] == g]
    a.plot([x[2] for x in p], [x[3] for x in p], "o-", color=c, label=f"gSS {g} nS")
a.axhline(24.4, color="0.5", ls=":"); a.text(0.32, 25.5, "Kito CSMC 24.4 mV", fontsize=7, color="0.4")
a.set_xscale("log"); a.set_xlabel("G$_{load}$ = m·√(gSS·g$_m$) (nS)"); a.set_ylabel("best SMC layer-1 amplitude with ICC-MY protected (mV)")
a.set_title("A  Different gSS collapse on G$_{load}$", loc="left", fontsize=10); a.legend(fontsize=8)
a = ax[1]
lab = [f"{x[0]:.0%} / ±{x[1]:.0%} / {x[2]}" for x in crit]
yy = np.arange(len(crit))[::-1]
for y, x in zip(yy, crit):
    a.plot([x[3], x[4]], [y, y], color="0.6", lw=1)
    a.plot([x[3]], [y], "o", color="tab:blue", ms=5); a.plot([x[4]], [y], "x", color="tab:red", ms=5)
a.set_yticks(yy); a.set_yticklabels(lab, fontsize=6); a.set_xscale("log")
a.set_xticks([1.5, 2, 3, 5, 10]); a.set_xticklabels(["1.5", "2", "3", "5", "10"]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_xlabel("G$_{load}$ (nS): ● largest with a window, × smallest without")
a.set_title("B  Threshold vs criteria (ICC amp / rate / SMC L1)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/loadrule.png", dpi=150)
print("wrote paper_R1/fig/loadrule.png")
