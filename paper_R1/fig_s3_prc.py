"""Fig. S3 (revision 1): phase response of the model ICC-MY and the current needed to trigger a slow wave, against its load.

Data: column runs out/paper_prc_col/events_s<seed>.csv (prc_column.py; gisim active, seeds 1-3, 200-msec pulses at 9 phases) and the
passive-load runs out/paper_prc/load_gL<gL>_a<amp>.csv + prc_a<amp>_d200.csv (run_paper_prc.sh; gisim iccprc, 20 or 40 phases).
Early wave = peak > -30 mV within the pulse + 300 msec and > 0.1 T0 before the control wave. Threshold = smallest tested amplitude
giving early waves at >= 1/3 of the phases (reported as the interval between the largest failing and the smallest passing amplitude).
Effective load of a column = g_c G_load / (g_c + G_load).
A  ICC-MY traces, column g_c 1.17 nS, seed 1: control and 200-msec pulses of 30 and 40 pA at phase 0.7 (prc_traces.py), with the
   upstroke times; the shift is read from the waves after the pulse, matched to the control waves by order after the reference
B  lasting phase shift against pulse phase, column g_c 1.17 nS (mean of seeds; band = range of the 1-pA pulses over seeds)
C  threshold current against effective load: passive stand-in (line) and real columns (symbols, three seeds)
Output: paper_R1/fig/figS3.png/.pdf, paper_R1/fig_s3_prc.out.txt
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle; figstyle.apply()
BLUE, PURPLE, RED, GREY = "#2B5DA8", "#7B4FA0", "#C0392B", "0.45"
DUR, FULL, SEEDS = 200.0, -30.0, (1, 2, 3)
COND = {"rabbit": (1.17, 1.61, BLUE, "o"), "mouse": (2.33, 1.61, PURPLE, "s"), "large": (3.0, 5.9, RED, "^")}
OUT = []


def analyse(rows, ctrl_rows, ref, T0, phases_on):
    """rows: {phase: [(t, peak)]}; returns {phase: (early, shift)}"""
    ct = [t for t, p in ctrl_rows if t > ref + 0.05 * T0]
    res = {}
    for ph, ev in rows.items():
        on = phases_on(ph)
        ev = [e for e in ev if e[0] > ref + 0.05 * T0]
        early = any(p > FULL and t < ct[0] - 0.1 * T0 for t, p in ev if on <= t <= on + DUR + 300)
        full = [t for t, p in ev if p > FULL]
        m = int(np.argmin([abs(full[-1] - x) for x in ct])); s = ((ct[m] - full[-1]) / T0 + 0.5) % 1.0 - 0.5
        res[ph] = (early, s)
    return res


def threshold(byamp):
    """byamp: {amp: {phase: (early, shift)}} -> (largest failing, smallest passing) amplitude"""
    amps = sorted(a for a in byamp if a > 1)
    ok = [a for a in amps if np.mean([e for e, s in byamp[a].values()]) >= 1 / 3]
    lo = max([a for a in amps if a not in ok and (not ok or a < min(ok))], default=0)
    return lo, (min(ok) if ok else np.nan)


# ---------- real columns ----------
col = {}     # (cond, seed) -> {amp: {phase: (early, shift)}}
for seed in SEEDS:
    d = {}
    for r in csv.DictReader(open(f"out/paper_prc_col/events_s{seed}.csv")):
        d.setdefault((r["cond"], float(r["amp_pA"])), {}).setdefault(float(r["phase"]), []).append((float(r["t_ms"]), float(r["peak_mV"])))
    for c in COND:
        ctrl = [e for e in d[(c, -1.0)][-1.0] if e[0] >= 20000]; ref = ctrl[0][0]; T0 = (ctrl[-1][0] - ref) / (len(ctrl) - 1)
        col[(c, seed)] = {a: analyse(rows, ctrl, ref, T0, lambda ph: ref + ph * T0) for (cc, a), rows in d.items() if cc == c and a > 0}
        col[(c, seed)]["T0"] = T0

# ---------- passive stand-in ----------
def load_iccprc(path):
    d = {}
    for r in csv.DictReader(open(path)):
        d.setdefault(float(r["phase"]), []).append((float(r["t_ms"]), float(r["peak_mV"])))
    ctrl = d.pop(-1.0); T0 = (ctrl[-1][0] - ctrl[0][0]) / (len(ctrl) - 1)
    return analyse(d, ctrl, 0.0, T0, lambda ph: ph * T0)

pas = {}
for gl in ("0", "0.6", "0.91", "2"):
    pas[float(gl)] = {a: load_iccprc(f"out/paper_prc/load_gL{gl}_a{a}.csv") for a in (25, 30, 35, 40, 45, 60, 80)}
for a in (10, 20):
    pas[0.0][a] = load_iccprc(f"out/paper_prc/prc_a{a}_d200.csv")

OUT.append("Threshold current (pA; largest failing .. smallest passing tested amplitude; early waves at >= 1/3 of phases)")
OUT.append("passive stand-in gL (nS): " + "; ".join(f"{g}: {threshold(v)[0]}..{threshold(v)[1]}" for g, v in pas.items()))
for c, (gc, G, col_, mk) in COND.items():
    ge = gc * G / (gc + G)
    th = [threshold({a: v for a, v in col[(c, s)].items() if a != "T0"}) for s in SEEDS]
    OUT.append(f"column {c} (g_c {gc}, G_load {G}, effective {ge:.2f} nS): " + "; ".join(f"seed {s}: {lo}..{hi}" for s, (lo, hi) in zip(SEEDS, th))
               + f"; T0 {np.mean([col[(c, s)]['T0'] for s in SEEDS]):.0f} msec")
for c in COND:
    noise = [s for seed in SEEDS for e, s in col[(c, seed)][1.0].values()]
    OUT.append(f"1-pA pulses ({c}): lasting shift {min(noise):+.3f}..{max(noise):+.3f} T0 over seeds")
    sub = [s for seed in SEEDS for a in (20, 30) for e, s in col[(c, seed)][float(a)].values()]
    OUT.append(f"  20-30 pA ({c}): lasting shift {min(sub):+.3f}..{max(sub):+.3f} T0; max advance above threshold "
               f"{max(s for seed in SEEDS for a in (50, 60) for e, s in col[(c, seed)][float(a)].values()):+.2f} T0")

# ---------- figure ----------
fig = plt.figure(figsize=(figstyle.WIDTH, 4.6))
t_ = fig.add_axes([0.08, 0.66, 0.88, 0.28]); a = fig.add_axes([0.08, 0.09, 0.36, 0.43]); b = fig.add_axes([0.56, 0.09, 0.33, 0.43])
ON = float(open("out/paper_prc_col/trace_onset.txt").read())         # pulse onset written by prc_traces.py
for f, cl, ls, lab in (("ctrl", "0.6", "-", "control"), ("a30", BLUE, "--", "30 pA"), ("a40", "k", "-", "40 pA")):
    r = list(csv.DictReader(open(f"out/paper_prc_col/trace_{f}.csv")))
    tt = np.array([float(x["t_ms"]) for x in r]); vv = np.array([float(x["V_ICC"]) for x in r])
    k = (tt >= 20000) & (tt <= 34000)
    t_.plot((tt[k] - ON) / 1000, vv[k], ls, color=cl, lw=0.8, label=lab)
    up = tt[np.where((vv[:-1] < -45) & (vv[1:] >= -45))[0]]; up = up[(up >= 20000) & (up <= 34000)]
    t_.plot((up - ON) / 1000, np.full(len(up), {"ctrl": -10, "a30": -6, "a40": -2}[f]), "|", color=cl, ms=5, mew=1)
    OUT.append(f"trace {f}: upstrokes (sec after pulse onset) " + ", ".join(f"{(u - ON) / 1000:.2f}" for u in up))
t_.axvspan(0, 0.2, color="#F3D9D5", lw=0)
t_.set_xlabel("time from pulse onset (sec)"); t_.set_ylabel("ICC-MY V$_m$ (mV)"); t_.set_xlim(-4.5, 9.5)
t_.legend(frameon=False, loc="center left"); t_.set_title("A", loc="left", fontweight="bold")
ph = sorted(col[("rabbit", 1)][1.0])
nz = np.array([[col[("rabbit", s)][1.0][p][1] for p in ph] for s in SEEDS])
a.fill_between(ph, nz.min(0), nz.max(0), color="0.85", lw=0, label="1 pA (noise)")
for amp, ls in ((20.0, ":"), (30.0, "--"), (40.0, "-"), (60.0, "-")):
    y = np.array([[col[("rabbit", s)][amp][p][1] for p in ph] for s in SEEDS])
    a.plot(ph, y.mean(0), ls, color=BLUE if amp < 60 else "k", marker="o", ms=2.5, label=f"{amp:.0f} pA")
a.axhline(0, color=GREY, lw=0.5)
a.set_xlabel("pulse phase (fraction of the period)"); a.set_ylabel("lasting phase shift (period; advance > 0)")
a.set_xlim(0, 1); a.set_ylim(-0.15, 0.55); a.legend(frameon=False, loc="upper left")
a.set_title("B", loc="left", fontweight="bold")

gs = sorted(pas)
lo = [threshold(pas[g])[0] for g in gs]; hi = [threshold(pas[g])[1] for g in gs]
b.fill_between(gs, lo, hi, color="0.85", lw=0, label="passive stand-in")
b.plot(gs, [(l + h) / 2 for l, h in zip(lo, hi)], color=GREY, lw=0.8)
for c, (gc, G, colr, mk) in COND.items():
    ge = gc * G / (gc + G)
    for i, s in enumerate(SEEDS):
        l, h = threshold({k: v for k, v in col[(c, s)].items() if k != "T0"})
        x = ge + (i - 1) * 0.03
        b.errorbar(x, (l + h) / 2, yerr=[[(h - l) / 2], [(h - l) / 2]], fmt=mk, color=colr, ms=3.5, capsize=1.5, lw=0.7,
                   label={"rabbit": "g$_c$ 1.17 nS", "mouse": "g$_c$ 2.33 nS", "large": "L 7, m 1.3, g$_c$ 3 nS"}[c] if i == 0 else None)
b.set_xlabel("effective load g$_c$G$_{load}$/(g$_c$+G$_{load}$) (nS)"); b.set_ylabel("current to trigger a slow wave (pA)")
b.set_xlim(-0.1, 2.2); b.set_ylim(0, 90); b.legend(frameon=False, loc="upper left")
r = b.secondary_yaxis("right", functions=(lambda i: i / 50.0, lambda g: g * 50.0)); r.set_ylabel("g$_{II}$ for a 50-mV neighbour (nS)")
b.set_title("C", loc="left", fontweight="bold")
figstyle.save(fig, "paper_R1/fig/figS3")
open("paper_R1/fig_s3_prc.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
print("\n".join(OUT))
