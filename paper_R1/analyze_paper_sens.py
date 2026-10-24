"""Priority 3: one-at-a-time sensitivity (paper_R1/run_paper_sens.ps1 -> out/paper_sens/<cond>/; base = out/v4win + out/v4ijp).
Per condition:
  g_c* = g_c at which the SMC layer-1 slow wave equals the Kito CSMC amplitude 24.4 mV (band 21.7-27.1 = +/- SD),
         by linear interpolation in log g_c over 0.5-3 nS ("<0.5" / ">3" if outside);
  SMC/ICC amplitude ratio at g_c 1.5 (Kito 24.4/56.9 = 0.43);
  IJP fraction ICC-MY / SMC layer 1 at g_c 1.5 (clean windows: ICC-MY below rest + 6 mV for the whole 1-s window,
         i.e. the -65 mV criterion of R1 shifted with the condition's rest; rest = 5th percentile of the control run;
         the clean test covers [t0, t0+1200) ms, 200 ms beyond the measurement window, to exclude delayed-upstroke artefacts).
Note iccGL -0.6: the ICC-MY does not oscillate uncoupled or at g_c 0.5 (0 cpm); it paces from g_c 1 on.
"SMC layer 1" = log index [0] / column V_SMC0.
Output: paper_R1/fig/sens.png, printed table (saved by the caller to analyze_paper_sens.out.txt).
"""
import os, re, csv
import numpy as np
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GCS = [0.5, 1, 1.5, 2, 3]
T0 = [20000, 20450, 20900, 21350, 21800, 22250, 22700, 23150]
KITO, KSD = 24.4, 2.7
CONDS = [("base", "v4 base")] + [(f"{k}{v}", f"{lab} {v}") for k, lab, vals in [
    ("gSS", "gSS", ["15", "60"]), ("L", "L", ["15", "45"]), ("m", "m", ["0.15", "0.45"]), ("gK0run", "gK0", ["0.2", "0.6"]),
    ("fP", "fP", ["0.05", "0.15"]), ("gSP", "gSP", ["1.5", "4.5"]), ("iccGL", "ICC +gL", ["-0.6", "1.8", "3.6"]),
    ("seed", "seed", ["2", "3"])] for v in vals]


def logf(c, gc):
    f = f"out/v4win/m0.3_gSS30_gc{gc}_u0.log" if c == "base" else f"out/paper_sens/{c}/gc{gc}_u0.log"
    m = PAT.search(AE.log(f))
    return tuple(float(x) for x in m.groups()) if m else None   # cpm, ICC amp, L0 cpm, SMC L1 amp


def ijp_pre(c):
    return "out/v4ijp/gc1.5" if c == "base" else f"out/paper_sens/{c}/ijp"


def load(f):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0")]


def rest(t, v):
    return np.percentile(v[t >= t[-1] / 2], 5)


def cross(gcs, y, target):
    """g_c where y(g_c) crosses target (log-linear), assuming y increases with g_c."""
    if y[0] > target: return 0.5, "<"
    if y[-1] < target: return float(gcs[-1]), ">"
    for k in range(len(y) - 1):
        if y[k] <= target <= y[k + 1]:
            a, b = np.log(gcs[k]), np.log(gcs[k + 1])
            return float(np.exp(a + (target - y[k]) / (y[k + 1] - y[k]) * (b - a))), ""


tb, tib, _ = load(ijp_pre("base") + "_ctrl.csv"); RB = rest(tb, tib)
res = {}
print("cond | uncoupled ICC amp, cpm | SMC L1 at g_c 0.5/1/1.5/2/3 | g_c* (21.7 / 24.4 / 27.1) | SMC/ICC at 1.5 | ICC rest | IJP n | SMC IJP | ICC IJP | ratio")
for c, lab in CONDS:
    d0 = logf(c, 0)
    gcs = [g for g in GCS + [4, 5] if os.path.exists((f"out/v4win/m0.3_gSS30_gc{g}_u0" if c == "base" else f"out/paper_sens/{c}/gc{g}_u0") + ".log")]
    d = [logf(c, g) for g in gcs]
    smc = [x[3] for x in d]
    gstar = [cross(gcs, smc, KITO + s) for s in (-KSD, 0, KSD)]
    ratio = d[2][3] / d[2][1]
    pre = ijp_pre(c)
    tc, ic, sc = load(pre + "_ctrl.csv"); rc = rest(tc, ic)
    thr = -65.0 + (rc - RB)
    rr = []
    for t0 in T0:
        t, i, s = load(f"{pre}_p{t0}.csv")
        n = min(len(t), len(tc)); w = (t[:n] >= t0) & (t[:n] < t0 + 1000)
        g = (t[:n] >= t0) & (t[:n] < t0 + 1200)   # 200-ms guard: a delayed next upstroke (timing shift) must not enter the window
        if ic[:n][g].max() < thr and i[:n][g].max() < thr:
            dS, dI = -(s[:n] - sc[:n])[w].min(), -(i[:n] - ic[:n])[w].min()
            rr.append((dS, dI, dI / dS))
    res[c] = dict(lab=lab, d0=d0, smc=smc, gstar=gstar, ratio=ratio, rest=rc, ijp=rr)
    gs = " / ".join(f"{f}{g:.2f}" for g, f in gstar)
    ij = (f"{len(rr)} | {np.mean([x[0] for x in rr]):.2f} | {np.mean([x[1] for x in rr]):.2f} | {np.mean([x[2] for x in rr]):.3f} "
          f"({min(x[2] for x in rr):.3f}-{max(x[2] for x in rr):.3f})") if rr else "0 clean windows"
    print(f"{lab:12s} | {d0[1]:.1f}, {d0[0]:.2f} | {' '.join(f'{x:.1f}' for x in smc)} | {gs} | {ratio:.3f} | {rc:.1f} | {ij}")

# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
names = [c for c, _ in CONDS]
y = np.arange(len(names))[::-1]
a = ax[0]
for yy, c in zip(y, names):
    (glo, flo), (gm, fm), (ghi, fhi) = res[c]["gstar"]
    a.plot([glo, ghi], [yy, yy], color="tab:blue", lw=4, alpha=0.35, solid_capstyle="butt")
    a.plot([gm], [yy], "o", color="tab:blue" if c != "base" else "k")
    for g, f in ((glo, flo), (ghi, fhi)):
        if f: a.text(g, yy + 0.25, f, fontsize=7, ha="center")
a.axvspan(1, 2, color="0.9", zorder=0)
a.set_xscale("log"); a.set_xlim(0.45, 3.3); a.set_xticks([0.5, 1, 1.5, 2, 3]); a.set_xticklabels(["0.5", "1", "1.5", "2", "3"]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a.set_yticks(y); a.set_yticklabels([res[c]["lab"] for c in names], fontsize=8)
a.set_xlabel("g$_c$ giving SMC layer-1 slow wave = Kito 24.4 mV\n(bar: 21.7–27.1 mV, ± SD)")
a.set_title("A  Coupling fixed by the amplitude", loc="left", fontsize=10)

a = ax[1]
for yy, c in zip(y, names):
    rr = res[c]["ijp"]
    if rr:
        q = [x[2] for x in rr]
        a.plot(q, [yy] * len(q), "o", color="tab:red", alpha=0.35, ms=4)
        a.plot([np.mean(q)], [yy], "o", color="tab:red" if c != "base" else "k", ms=7)
        a.text(max(q) + 0.01, yy - 0.15, f"n {len(q)}", fontsize=6)
a.axvline(np.mean([x[2] for x in res["base"]["ijp"]]), color="k", lw=0.7, ls=":")
a.set_yticks(y); a.set_yticklabels([])
a.set_xlabel("ICC-MY hyperpol. / SMC layer-1 IJP at g$_c$ 1.5 nS")
a.set_title("B  IJP fraction reaching ICC-MY", loc="left", fontsize=10)

a = ax[2]
GOWN0 = 3.74   # isolated ICC-MY slope conductance at -72.5 mV (out/paper_rin, m 0.3 g_c 0)
gl = {"base": 0.0, "iccGL-0.6": -0.6, "iccGL1.8": 1.8, "iccGL3.6": 3.6}
xx = np.linspace(2, 8, 50); a.plot(xx, 1.5 / (1.5 + xx), color="0.5", label="g$_c$/(g$_c$+G$_{own}$), g$_c$ 1.5")
for c, g in gl.items():
    rr = res[c]["ijp"]
    if rr:
        q = [x[2] for x in rr]
        a.plot([GOWN0 + g] * len(q), q, "o", color="tab:purple", alpha=0.35, ms=4)
        a.plot([GOWN0 + g], [np.mean(q)], "o", color="tab:purple", ms=7)
a.set_xlabel("ICC-MY own conductance at rest (nS)\n(3.74 isolated + iccGL)"); a.set_ylabel("IJP fraction at g$_c$ 1.5 nS")
a.set_title("C  Own-conductance axis", loc="left", fontsize=10); a.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join("paper_R1", "fig", "sens.png"), dpi=150)
print("wrote paper_R1/fig/sens.png")
