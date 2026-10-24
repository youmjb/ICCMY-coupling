"""Revision 1 runs (run_paper_rev1.ps1, run_paper_rev1b.ps1).
(1) Load-limit generality, mouse-like column L 7, m 1.3 (G_load 5.9 nS): base (out/paper_anat/L7m1.3_gc*), ICC-MY own conductance
    -0.6 / +1.8 nS (iccGL) and a 25-pF ICC-MY (iccScale 1.0). For each: muscle layer-1 amplitude and ICC-MY amplitude vs g_c, the
    largest muscle amplitude, and the ICC-MY amplitude where the muscle first reaches 24.4 mV (linear interpolation).
(2) Time-step convergence at the calibrated column (dt 0.1 = out/v4win; 0.05, 0.025, 0.0125 = out/paper_rev1): amplitudes and
    g_c at the recorded ratio 0.429 where it is bracketed.
Output: printed (saved to analyze_rev1.out.txt).
"""
import os, re
import numpy as np
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY [\d.]+ cpm amp ([\d.]+) \| L0 [\d.]+ cpm \| amp by layer: \[0\]([\d.]+)")
R = 24.4 / 56.9


def rd(pre):
    if not os.path.exists(pre + ".log"): return None
    m = PAT.search(AE.log(pre + ".log"))
    return (float(m.group(1)), float(m.group(2))) if m else None


print("(1) load limit, L 7 m 1.3 (G_load 5.9 nS): g_c: ICC-MY / muscle layer-1 amplitude (mV)")
VAR = [("18 pF, base", lambda g: f"out/paper_anat/L7m1.3_gc{g}", [1, 2, 3, 4, 6, 8, 10, 12]),
       ("18 pF, own conductance -0.6 nS", lambda g: f"out/paper_rev1/L7_gLm0.6_gc{g}", [1, 2, 3, 4, 6, 8, 12]),
       ("18 pF, own conductance +1.8 nS", lambda g: f"out/paper_rev1/L7_gL1.8_gc{g}", [1, 2, 3, 4, 6, 8, 12]),
       ("25 pF (iccScale 1.0)", lambda g: f"out/paper_rev1/L7_s1.0_gc{g}", [1, 2, 3, 4, 6, 8, 12])]
for lab, pre, gcs in VAR:
    v = [(g,) + rd(pre(g)) for g in gcs if rd(pre(g))]
    s = [x[2] for x in v]; i = [x[1] for x in v]
    reach = next((k for k in range(1, len(v)) if s[k - 1] < 24.4 <= s[k]), None)
    if reach is not None:
        f = (24.4 - s[reach - 1]) / (s[reach] - s[reach - 1]); icc = i[reach - 1] + f * (i[reach] - i[reach - 1])
        txt = f"muscle reaches 24.4 mV with ICC-MY at {icc:.1f} mV ({icc / 56.9:.0%} of the recorded 56.9 mV)"
    else:
        k = int(np.argmax(s)); txt = f"muscle never reaches 24.4 mV (max {s[k]:.1f} mV at g_c {v[k][0]}, ICC-MY {i[k]:.1f} mV)"
    print(f"  {lab}: " + "; ".join(f"{g}: {a:.1f}/{b:.1f}" for g, a, b in v) + f" -> {txt}")

print("(2) time step, calibrated column: dt (ms) | g_c (nS): ICC-MY / muscle / ratio | g_c at ratio 0.429")
for dt, pre in (("0.1", "out/v4win/m0.3_gSS30_gc{}_u0"), ("0.05", "out/paper_rev1/dt0.05_gc{}"),
                ("0.025", "out/paper_rev1/dt0.025_gc{}"), ("0.0125", "out/paper_rev1/dt0.0125_gc{}")):
    v = [(g,) + rd(pre.format(g)) for g in (0.5, 0.75, 1, 1.5, 2) if rd(pre.format(g))]
    rat = [(g, b / a) for g, a, b in v]
    cross = next((np.exp(np.log(g0) + (R - r0) / (r1 - r0) * np.log(g1 / g0)) for (g0, r0), (g1, r1) in zip(rat[:-1], rat[1:]) if r0 <= R <= r1), None)
    print(f"  {dt} | " + "; ".join(f"{g}: {a:.1f}/{b:.1f}/{b / a:.3f}" for g, a, b in v) + (f" | {cross:.2f}" if cross else " | not bracketed"))
