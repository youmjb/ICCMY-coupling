"""Supplement: ICC-DMP (excitatory, cholinergic) input. ACh x10 (aI 9) minus the DMP baseline (aI 0, same L and seed) -- one
baseline for every site (R1 Fig. 5 mixed DMP and no-DMP baselines). Seeds 1 (out/dmp1) and 2-3 (out/paper_dmp); L 10/15/30 layers.
Sites (second half of the 60-s run, t >= 30 s):
  SMC next to ICC-DMP = mean of V_SMClast (deepest layer);
  SMC layer 1 and ICC-MY = median between slow waves (ICC-MY below -65 mV and > 300 ms from any sample above it).
ICC-DMP parameters are placeholders (HANDOFF gap 13); ICC-MY rate is not reported (voltage-insensitive pacemaker, gap 3).
Output: paper_R1/fig/supp_dmp.png and printed table.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LS, SEEDS = [10, 15, 30], [1, 2, 3]


def path(L, a, s):
    if s == 1:
        return f"out/dmp1/{'' if L == 30 else f'L{L}_'}f3_a{a}.csv"
    return f"out/paper_dmp/L{L}_a{a}_s{s}.csv"


def load(f):
    r = list(csv.DictReader(open(f)))
    return {k: np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0", "V_SMClast")}


def sites(f):
    d = load(f); h = d["t_ms"] >= 30000
    vi = d["V_ICC"][h]
    near = np.convolve((vi >= -65).astype(float), np.ones(601), mode="same") > 0
    bw = ~near
    return np.array([d["V_SMClast"][h].mean(), np.median(d["V_SMC0"][h][bw]), np.median(vi[bw])])


D = {}
print("L | seed | dV deep SMC | dV SMC layer 1 | dV ICC-MY (mV, ACh - DMP baseline) | ICC-MY/SMC L1")
for L in LS:
    for s in SEEDS:
        D[(L, s)] = sites(path(L, 9, s)) - sites(path(L, 0, s))
        v = D[(L, s)]
        print(f"{L} | {s} | {v[0]:6.2f} | {v[1]:6.2f} | {v[2]:6.2f} | {v[2] / v[1] if abs(v[1]) > 0.2 else float('nan'):.2f}")
    m = np.mean([D[(L, s)] for s in SEEDS], axis=0)
    print(f"{L} | mean | {m[0]:6.2f} | {m[1]:6.2f} | {m[2]:6.2f} | {m[2] / m[1] if abs(m[1]) > 0.2 else float('nan'):.2f}")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw=dict(width_ratios=[1.3, 1]))
a = ax[0]
d0, d9 = load(path(10, 0, 1)), load(path(10, 9, 1))
w = (d0["t_ms"] >= 40000) & (d0["t_ms"] < 48000); tt = d0["t_ms"][w] / 1000 - 40
for k, c, lab in (("V_SMClast", "tab:purple", "SMC next to ICC-DMP"), ("V_SMC0", "tab:blue", "SMC layer 1"), ("V_ICC", "tab:orange", "ICC-MY")):
    a.plot(tt, d0[k][w], color=c, lw=1.0, ls="--"); a.plot(tt, d9[k][w], color=c, lw=1.6, label=lab)
a.set_xlabel("time (s)"); a.set_ylabel("Vm (mV)")
a.set_title("A  10 layers, DMP baseline (dashed) vs ACh x10 (solid)", loc="left", fontsize=10); a.legend(fontsize=7, loc="upper right")
a = ax[1]
x = np.arange(3); bw = 0.26
for j, (c, lab) in enumerate((("tab:purple", "SMC next to ICC-DMP"), ("tab:blue", "SMC layer 1"), ("tab:orange", "ICC-MY"))):
    vals = np.array([[D[(L, s)][j] for s in SEEDS] for L in LS])
    a.bar(x + (j - 1) * bw, vals.mean(axis=1), bw, color=c, alpha=0.75, label=lab)
    for k in range(3):
        a.plot([x[k] + (j - 1) * bw] * 3, vals[k], "k.", ms=4)
a.axhline(0, color="k", lw=0.6)
a.set_xticks(x); a.set_xticklabels([f"{L} layers" for L in LS])
a.set_ylabel("depolarization by ACh (mV)"); a.legend(fontsize=7)
a.set_title("B  Three sites, 3 seeds (dots)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/supp_dmp.png", dpi=150)
print("wrote paper_R1/fig/supp_dmp.png")
