"""Priority 2 analysis (run by paper_R1/run_paper_a.ps1).
A  out/paper_rmp : slow-wave amplitudes vs ICC-MY rest, shifted by a sustained inward current (iccInj).
B  out/paper_rijp_i10, _i20, paper_rijp (iccInj -10/-20/-30 pA): purinergic IJP ratio ICC-MY / SMC layer 1, vs out/v4ijp (iccInj 0).
C  out/paper_rin : ICC-MY input resistance from +10 pA, 300 ms pulses vs the analytic estimate (R1 Fig. 2B) and Cai 2017.
"SMC layer 1" = column V_SMC0 / log index [0], as in make_report_figs.py.
Output: paper_R1/fig/rmp_rin.png and a printed table.
"""
import os, re, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
FIG = os.path.join("paper_R1", "fig")
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
T0 = [20000, 20450, 20900, 21350, 21800, 22250, 22700, 23150]
G1, GOWN = 5.37, 3.6


def load(f, cols=("t_ms", "V_ICC", "V_SMC0")):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[c]) for x in r]) for c in cols]


def rest(t, v):
    """ICC-MY rest = 5th percentile of Vm over the second half of the run."""
    return np.percentile(v[t >= t[-1] / 2], 5)


# ---------- A: RMP shift ----------
GCS, INJ = [0, 1, 1.5, 2, 5], [0, -10, -20, -30, -40]
A = {}
print("A  g_c | iccInj | ICC rest | ICC peak | ICC amp | cpm | SMC L1 amp | SMC/ICC")
for gc in GCS:
    for i in INJ:
        pre = f"out/paper_rmp/gc{gc}_inj{i}"
        m = PAT.search(AE.log(pre + ".log"))
        cpm, ampI, _, ampS = (float(x) for x in m.groups())
        t, vi, _ = load(pre + ".csv")
        h = t >= t[-1] / 2
        A[(gc, i)] = (rest(t, vi), vi[h].max(), ampI, cpm, ampS)
        r = A[(gc, i)]
        print(f"   {gc:4} | {i:4} | {r[0]:6.1f} | {r[1]:6.1f} | {r[2]:5.1f} | {r[3]:5.2f} | {r[4]:5.1f} | {r[4] / r[2]:.3f}")


# ---------- B: IJP with depolarized ICC-MY ----------
def ijp(pre, thr):
    tc, ic, sc = load(f"{pre}_ctrl.csv")
    out = []
    for t0 in T0:
        t, i, s = load(f"{pre}_p{t0}.csv")
        n = min(len(t), len(tc)); w = (t[:n] >= t0) & (t[:n] < t0 + 1000)
        if ic[:n][w].max() < thr and i[:n][w].max() < thr:
            dS, dI = -(s[:n] - sc[:n])[w].min(), -(i[:n] - ic[:n])[w].min()
            out.append((t0, dS, dI, dI / dS))
    return out, rest(tc, ic)


B = {}
RIJP = [("0", "out/v4ijp"), ("-10", "out/paper_rijp_i10"), ("-20", "out/paper_rijp_i20"), ("-30", "out/paper_rijp")]
print("\nB  g_c | iccInj | ICC rest | clean n | SMC L1 IJP | ICC-MY | ratio (mean, range)")
for gc in [1, 1.5, 2]:
    tc, ic, _ = load(f"out/v4ijp/gc{gc}_ctrl.csv"); r0 = rest(tc, ic)
    for lab, d in RIJP:
        if not os.path.exists(f"{d}/gc{gc}_ctrl.csv"):
            continue
        tc, ic, _ = load(f"{d}/gc{gc}_ctrl.csv"); r1 = rest(tc, ic)
        rr, _ = ijp(f"{d}/gc{gc}", -65.0 + (r1 - r0))   # same margin above rest as the -65 mV criterion at iccInj 0
        B[(gc, lab)] = (rr, r1)
        if rr:
            q = [x[3] for x in rr]
            print(f"   {gc:4} | {lab:>4} | {r1:6.1f} | {len(rr)} | {np.mean([x[1] for x in rr]):5.2f} | {np.mean([x[2] for x in rr]):5.2f} | {np.mean(q):.3f} ({min(q):.3f}-{max(q):.3f})")
        else:
            print(f"   {gc:4} | {lab:>4} | {r1:6.1f} | 0 clean windows")

# ---------- C: direct R_in ----------
CONDS = [(0.3, 0), (0.3, 1), (0.3, 1.5), (0.3, 2), (0.3, 5), (0.3, 20), (1, 1), (1, 2), (1, 5)]
C = {}
print("\nC  m | g_c | clean n | R_in direct (MOhm, mean, range) | analytic | dV ratio 150/300 ms | g_own nS | column G nS (m*G1) | ICC Vm before")
for m, gc in CONDS:
    pre = f"out/paper_rin/m{m}_gc{gc}"
    tc, ic, sc = load(f"{pre}_ctrl.csv")
    rr = []
    for t0 in T0:
        t, i, s = load(f"{pre}_i{t0}.csv")
        n = min(len(t), len(tc)); w = (t[:n] >= t0 - 50) & (t[:n] < t0 + 400)
        if ic[:n][w].max() >= -65 or i[:n][w].max() >= -65:
            continue
        d = i[:n] - ic[:n]
        late = d[(t[:n] >= t0 + 250) & (t[:n] < t0 + 300)].mean()
        mid = d[(t[:n] >= t0 + 125) & (t[:n] < t0 + 175)].mean()
        wl = (t[:n] >= t0 + 250) & (t[:n] < t0 + 300)
        dSm = (s[:n] - sc[:n])[wl].mean()
        ig = gc * (late - dSm)                       # extra gap current ICC-MY -> layer 0 (pA, negative)
        gown = (10 + ig) / -late                     # ICC-MY own slope conductance (nS): injected 10 pA = own + gap
        gcol = ig / dSm if gc > 0 else np.nan        # conductance of the m-SMC column seen at layer 0 (nS)
        vpre = ic[:n][(t[:n] >= t0 - 20) & (t[:n] < t0)].mean()
        rr.append((-late / 10 * 1000, mid / late, gown, gcol, vpre))
    G = m * G1
    an = 1000 / (GOWN + gc * G / (gc + G))
    C[(m, gc)] = (rr, an)
    if rr:
        v = [x[0] for x in rr]
        print(f"   {m} | {gc:4} | {len(rr)} | {np.mean(v):6.1f} ({min(v):.1f}-{max(v):.1f}) | {an:6.1f} | {np.mean([x[1] for x in rr]):.3f} | {np.mean([x[2] for x in rr]):.2f} | {np.nanmean([x[3] for x in rr]):.2f} ({min(x[3] for x in rr):.2f}-{max(x[3] for x in rr):.2f}) ({G:.2f}) | {np.mean([x[4] for x in rr]):.1f}")
    else:
        print(f"   {m} | {gc:4} | 0 clean windows | analytic {an:.1f}")

# ---------- figure ----------
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5)); ax = ax.ravel()
cols = {0: "0.5", 1: "tab:blue", 1.5: "tab:green", 2: "tab:orange", 5: "tab:red"}
a = ax[0]
for gc in GCS:
    ii = [i for i in INJ if A[(gc, i)][3] > 0 and A[(gc, i)][0] < -40]   # oscillating, rest below -40 mV
    x = [A[(gc, i)][0] for i in ii]
    a.plot(x, [A[(gc, i)][2] for i in ii], "-o", color=cols[gc], ms=4, label=f"g$_c$ {gc}: ICC-MY")
    if gc > 0:
        a.plot(x, [A[(gc, i)][4] for i in ii], "--s", color=cols[gc], ms=4, label=f"g$_c$ {gc}: SMC layer 1")
a.axhspan(24.4 - 2.7, 24.4 + 2.7, color="0.9", zorder=0)
a.text(-74.5, 25.5, "Kito CSMC 24.4 ± 2.7 (SD)", fontsize=7)
a.axvspan(-65.1 - 2.8, -62.2 + 2.8, color="tab:purple", alpha=0.07)
a.text(-63.6, 57, "recorded ICC-MY rest\n(Kito, Cai)", fontsize=7, ha="center", color="tab:purple")
a.set_xlim(-75, -40); a.set_ylim(0, 62)
a.set_xlabel("ICC-MY rest (mV), shifted by iccInj 0 to -40 pA"); a.set_ylabel("slow-wave amplitude (mV)")
a.set_title("A  Amplitudes vs ICC-MY rest (m 0.3)", loc="left", fontsize=10)
a.legend(fontsize=6, ncol=2, loc="lower left")

a = ax[1]
gg = np.linspace(0.8, 2.2, 50); a.plot(gg, gg / (gg + GOWN), color="0.5", label="g$_c$/(g$_c$+3.6 nS)")
ccol = {"0": "tab:blue", "-10": "tab:green", "-20": "tab:orange", "-30": "tab:red"}
dx = {"0": -0.06, "-10": -0.02, "-20": 0.02, "-30": 0.06}
seen = set()
for (gc, lab), (rr, r1) in B.items():
    if rr:
        q = [x[3] for x in rr]
        a.scatter([gc + dx[lab]] * len(q), q, color=ccol[lab], alpha=0.4, s=12)
        a.scatter([gc + dx[lab]], [np.mean(q)], color=ccol[lab], edgecolor="k", s=45,
                  label=None if lab in seen else f"iccInj {lab} pA")
        a.annotate(f"{r1:.0f}", (gc + dx[lab], np.mean(q)), xytext=(4, 3), textcoords="offset points", fontsize=6)
        seen.add(lab)
a.set_yscale("log")
a.set_xlabel("g$_c$ (nS)   (numbers: ICC-MY rest, mV)"); a.set_ylabel("ICC-MY hyperpol. / SMC layer-1 IJP")
a.set_title("B  IJP fraction reaching ICC-MY", loc="left", fontsize=10); a.legend(fontsize=7)

a = ax[2]
a.axhspan(260 - 2.179 * 38, 260 + 2.179 * 38, color="0.88", label="Cai 95% CI of mean")
gg = np.logspace(-1, np.log10(30), 100)
for m, c in ((0.3, "tab:blue"), (1, "tab:orange")):
    a.plot(gg, 1000 / (GOWN + gg * m * G1 / (gg + m * G1)), color=c, label=f"analytic, m {m}")
    for (mm, gc), (rr, an) in C.items():
        if mm == m and rr:
            v = [x[0] for x in rr]; x0 = max(gc, 0.1)
            a.scatter([x0] * len(v), v, color=c, alpha=0.4, s=15)
            a.scatter([x0], [np.mean(v)], color=c, edgecolor="k", s=50)
a.set_xscale("log"); a.set_xlabel("g$_c$ (nS)  (g$_c$ 0 plotted at 0.1)"); a.set_ylabel("ICC-MY R$_{in}$ (MΩ)")
a.set_title("C  Direct R$_{in}$ (+10 pA, 300 ms; dots) vs analytic", loc="left", fontsize=10); a.legend(fontsize=7)

a = ax[3]
for (m, gc), (rr, an) in C.items():
    if rr:
        a.scatter([x[4] for x in rr], [x[2] for x in rr], s=18, color=cols.get(gc, "k") if gc != 20 else "k",
                  marker="o" if m == 0.3 else "^", label=f"m {m}, g$_c$ {gc}")
a.axhline(GOWN, color="0.5", ls=":"); a.text(-71.5, GOWN + 0.03, "assumed 3.6 nS", fontsize=7, color="0.4")
a.set_xlabel("ICC-MY Vm just before the pulse (mV)"); a.set_ylabel("ICC-MY own slope conductance (nS)")
a.set_title("D  ICC-MY own conductance falls as it depolarizes", loc="left", fontsize=10); a.legend(fontsize=6, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "rmp_rin.png"), dpi=150)
print("wrote", os.path.join(FIG, "rmp_rin.png"))
