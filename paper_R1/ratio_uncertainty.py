"""Uncertainty of the target amplitude ratio (reviewer point 1-2) and its propagation to g_c / G_load, for both targets.
Rabbit ileum (Kito et al. 2014, Table 1): ICC-MY upstroke 56.9 +/- 6.4 mV (n 10), circular SMC 24.4 +/- 2.7 mV (n 38).
Mouse small intestine (Kito & Suzuki 2003, 35 C): ICC-MY pacemaker potential 56.5 +/- 5.3 mV (n 39), circular SMC slow wave
(primary component) 34.1 +/- 4.6 mV (n 34). Means +/- SD from different cells (unpaired) in both.
(a) ratio of population means: Monte Carlo over the two sample means (SEM = SD/sqrt(n)) -> 95 % interval;
(b) cell-to-cell spread: ratio of single cells drawn from the two distributions -> central 68 % and 95 %.
Propagation: passive divider g_c/Y = r/(1-r), and the full model (calibrated column, out/v4win amplitudes interpolated in log g_c),
as g_c/G_load with G_load 1.61 nS.
"""
import os, re
import numpy as np
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rng = np.random.default_rng(1)
N = 200000
TARGETS = [("rabbit ileum", 56.9, 6.4, 10, 24.4, 2.7, 38), ("mouse small intestine", 56.5, 5.3, 39, 34.1, 4.6, 34)]
PAT = re.compile(r"amp ([\d.]+) \| L0 [\d.]+ cpm \| amp by layer: \[0\]([\d.]+)")
gcs = [0.5, 1, 1.5, 2, 3, 5]
rat = []
for g in gcs:
    a, s = (float(x) for x in PAT.search(AE.log(f"out/v4win/m0.3_gSS30_gc{g}_u0.log")).groups())
    rat.append(s / a)
G = 1.61


def model_gc(r):
    return np.exp(np.interp(r, rat, np.log(gcs), left=np.nan, right=np.nan))


def q(x, p):
    return np.nanpercentile(x, p)


print("model ratio at g_c", dict(zip(gcs, [round(x, 3) for x in rat])))
for name, mi, si, ni, ms, ss, ns in TARGETS:
    print(f"== {name}: point ratio {ms / mi:.3f}")
    r_mean = rng.normal(ms, ss / np.sqrt(ns), N) / rng.normal(mi, si / np.sqrt(ni), N)
    r_cell = rng.normal(ms, ss, N) / rng.normal(mi, si, N)
    for lab, r in (("(a) ratio of means", r_mean), ("(b) single-cell ratio", r_cell)):
        pv = r / (1 - r); gm = model_gc(r) / G
        print(f"  {lab}: r 95 % {q(r, 2.5):.3f}-{q(r, 97.5):.3f} (68 % {q(r, 16):.3f}-{q(r, 84):.3f}); "
              f"passive g_c/Y 95 % {q(pv, 2.5):.2f}-{q(pv, 97.5):.2f}; "
              f"full model g_c/G_load 95 % {q(gm, 2.5):.2f}-{q(gm, 97.5):.2f} (68 % {q(gm, 16):.2f}-{q(gm, 84):.2f}); "
              f"outside simulated g_c 0.5-5 nS: {np.isnan(gm).mean():.1%}")
