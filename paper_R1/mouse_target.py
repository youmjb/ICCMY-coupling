"""Species of the amplitude-ratio target (user point, 2026-10-15). The ICC-MY model is a mouse small-intestine model, but the
target ratio came from rabbit ileum (Kito et al. 2014: 24.4/56.9 = 0.429). Kito & Suzuki 2003 (J Physiol 553:803, mouse small
intestine, 35 C; checked in the PMC full text) give ICC-MY pacemaker potentials 56.5 +/- 5.3 mV (n 39) and circular-muscle slow
waves 34.1 +/- 4.6 mV (n 34; primary component) -> ratio 0.604. Here: g_c at the mouse ratio for every saved condition (same
curves and load definition as analyze_criteria.py), the passive divider value, and the uncertainty of the ratio of means.
Output: printed (saved to mouse_target.out.txt).
"""
import io, contextlib, math, sys
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import analyze_criteria as ac

R_MOUSE, R_RABBIT = 34.1 / 56.5, 24.4 / 56.9
print(f"ratios: rabbit {R_RABBIT:.3f}, mouse {R_MOUSE:.3f}; passive g_c/Y: rabbit {R_RABBIT / (1 - R_RABBIT):.3f}, mouse {R_MOUSE / (1 - R_MOUSE):.3f}")
rng = np.random.default_rng(1); N = 200000
rm = rng.normal(34.1, 4.6 / math.sqrt(34), N) / rng.normal(56.5, 5.3 / math.sqrt(39), N)
print("mouse ratio of means 95 %:", np.round(np.percentile(rm, [2.5, 97.5]), 3))
print("case | group | G_load | g_c rabbit (0.429) | /G | g_c mouse (0.604) | /G | max simulated ratio")
summ = {}
for lab, m, gss, L, pre, gcs, grp in ac.CASES:
    c = ac.curve(pre, gcs)
    if not c: continue
    G = ac.gload(m, gss, L); g = [x[0] for x in c]; r = [x[2] / x[1] for x in c]
    gr, gm = ac.cross(g, r, R_RABBIT), ac.cross(g, r, R_MOUSE)
    summ.setdefault(grp, []).append(gm / G)
    print(f"{lab} | {grp} | {G:.2f} | {gr:.2f} | {gr / G:.2f} | {gm:.2f} | {gm / G:.2f} | {max(r):.3f}")
for grp, v in summ.items():
    x = [t for t in v if not math.isnan(t)]
    print(f"SUMMARY {grp} mouse g_c/G_load: {min(x):.2f}-{max(x):.2f} (n {len(x)} of {len(v)})")
