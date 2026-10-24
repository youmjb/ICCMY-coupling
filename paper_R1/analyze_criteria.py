"""QC follow-up (qc_report_2026-10-14 A1-A3): g_c under three criteria for every saved condition, with ONE definition of the load,
G_load = m*sqrt(g_SS*g_m)*tanh(L/lambda), lambda = sqrt(g_SS/g_m) layers, g_m = 0.96 nS (sealed cable of L layers).
Criteria (Kito et al. 2014, rabbit ileum, Table 1):
  abs   : muscle layer-1 amplitude = 24.4 mV (CSMC), band 21.7-27.1 (+/- SD)
  r_up  : muscle / ICC-MY amplitude = 24.4 / 56.9 = 0.429 (ICC-MY upstroke amplitude)
  r_pl  : muscle / ICC-MY amplitude = 24.4 / 49.7 = 0.491 (ICC-MY plateau amplitude)
The model's ICC-MY amplitude (log value) is peak minus minimum, i.e. the upstroke amplitude, so r_up is the like-for-like ratio.
Output: paper_R1/criteria.out.txt (table) -- no new simulations.
"""
import os, re, math
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GM = 0.96
R_UP, R_PL = 24.4 / 56.9, 24.4 / 49.7
R_MOUSE = 34.1 / 56.5          # Kito & Suzuki 2003, mouse small intestine (CSMC primary component / ICC-MY pacemaker potential)
R_MPL = 31.1 / 56.5            # same, CSMC plateau component (3rd review, 4)


def rd(f):
    x = PAT.search(AE.log(f + ".log")); return float(x.group(2)), float(x.group(4))


def curve(pre, gcs):
    return [(g,) + rd(pre(g)) for g in gcs if os.path.exists(pre(g) + ".log")]


def cross(xs, ys, t):
    for (x0, y0), (x1, y1) in zip(zip(xs[:-1], ys[:-1]), zip(xs[1:], ys[1:])):
        if (y0 - t) * (y1 - t) <= 0 and y1 != y0:
            return math.exp(math.log(x0) + (t - y0) / (y1 - y0) * (math.log(x1) - math.log(x0)))
    return float("nan")


def gload(m, gss, L):
    lam = math.sqrt(gss / GM); return m * math.sqrt(gss * GM) * math.tanh(L / lam)


GC6 = [0.5, 1, 1.5, 2, 3, 5]
GC5 = [0.5, 1, 1.5, 2, 3]
CASES = [  # label, m, g_SS, L, run prefix, g_c list, group
    ("calibrated column", 0.3, 30, 30, lambda g: f"out/v4win/m0.3_gSS30_gc{g}_u0", GC6, "geometry"),
    ("g_SS 15", 0.3, 15, 30, lambda g: f"out/paper_sens/gSS15/gc{g}_u0", GC5, "geometry"),
    ("g_SS 60", 0.3, 60, 30, lambda g: f"out/paper_sens/gSS60/gc{g}_u0", GC5 + [4, 5], "geometry"),
    ("m 0.15", 0.15, 30, 30, lambda g: f"out/paper_sens/m0.15/gc{g}_u0", GC5, "geometry"),
    ("m 0.45", 0.45, 30, 30, lambda g: f"out/paper_sens/m0.45/gc{g}_u0", GC5 + [4, 5], "geometry"),
    ("L 15", 0.3, 30, 15, lambda g: f"out/paper_sens/L15/gc{g}_u0", GC5, "geometry"),
    ("L 45", 0.3, 30, 45, lambda g: f"out/paper_sens/L45/gc{g}_u0", GC5, "geometry"),
    ("L 7, m 1.3", 1.3, 30, 7, lambda g: f"out/paper_anat/L7m1.3_gc{g}", [1, 2, 3, 4, 6, 8, 10, 12], "geometry"),
    ("L 8, m 3, g_SS 30", 3, 30, 8, lambda g: f"out/paper_anat/gSS30_gc{g}", [3, 6, 10, 15, 25, 40], "geometry"),
    ("L 8, m 3, g_SS 300", 3, 300, 8, lambda g: f"out/paper_anat/gSS300_gc{g}", [5, 10, 15, 20, 30, 50], "geometry"),
    ("ICC-MY trapezoid HW 0.67 s, 3.3 s", 0.3, 30, 30, lambda g: f"out/paper_wave/W1_gc{g}", GC6, "robustness"),
    ("ICC-MY HW 1.5 s, 6 s", 0.3, 30, 30, lambda g: f"out/paper_wave/W4_gc{g}", GC6, "robustness"),
    ("ICC-MY HW 2.4 s, 6 s", 0.3, 30, 30, lambda g: f"out/paper_wave/W5_gc{g}", GC6, "robustness"),
    ("ICC-MY upstroke 5 ms", 0.3, 30, 30, lambda g: f"out/paper_wave/W6_gc{g}", GC6, "robustness"),
    ("Kito-like ICC-MY wave", 0.3, 30, 30, lambda g: f"out/paper_wave/W7_gc{g}", GC6, "robustness"),
    ("L-type block 90 %", 0.3, 30, 30, lambda g: f"out/paper_nif/cal0.1/gc{g}_u0", GC6, "robustness"),
    ("L-type block 100 %", 0.3, 30, 30, lambda g: f"out/paper_nif/cal0/gc{g}_u0", GC6, "robustness"),
    ("SMC model Poh 2012", 0.3, 30, 30, lambda g: f"out/paper_poh/gK0_gc{g}", GC6, "robustness"),
    ("minimal FLC model", 0.3, 30, 30, lambda g: f"out/paper_flcmf/gc{g}_u0", GC5, "robustness"),
    ("tonic K 0.2", 0.3, 30, 30, lambda g: f"out/paper_sens/gK0run0.2/gc{g}_u0", GC5, "robustness"),
    ("tonic K 0.6", 0.3, 30, 30, lambda g: f"out/paper_sens/gK0run0.6/gc{g}_u0", GC5, "robustness"),
    ("ICC-MY +gL 1.8 nS", 0.3, 30, 30, lambda g: f"out/paper_sens/iccGL1.8/gc{g}_u0", GC5 + [4, 5], "robustness"),
]
lines = ["case | group | G_load (nS) | g_c abs 24.4 (21.7-27.1) | /G | g_c ratio 0.429 | /G | g_c ratio 0.491 | /G | g_c mouse ratio 0.604 | /G | g_c mouse plateau 0.550 | /G"]
summ = {}
for lab, m, gss, L, pre, gcs, grp in CASES:
    c = curve(pre, gcs)
    if not c: lines.append(f"{lab} | missing"); continue
    G = gload(m, gss, L); g = [x[0] for x in c]
    ab = [cross(g, [x[2] for x in c], 24.4 + d) for d in (-2.7, 0, 2.7)]
    ru = cross(g, [x[2] / x[1] for x in c], R_UP); rp = cross(g, [x[2] / x[1] for x in c], R_PL)
    rm = cross(g, [x[2] / x[1] for x in c], R_MOUSE); rq = cross(g, [x[2] / x[1] for x in c], R_MPL)
    summ.setdefault(grp, []).append((ab[1] / G, ru / G, rp / G, rm / G, rq / G))
    lines.append(f"{lab} | {grp} | {G:.2f} | {ab[1]:.2f} ({ab[0]:.2f}-{ab[2]:.2f}) | {ab[1] / G:.2f} | {ru:.2f} | {ru / G:.2f} | {rp:.2f} | {rp / G:.2f} | {rm:.2f} | {rm / G:.2f} | {rq:.2f} | {rq / G:.2f}")
for grp, v in summ.items():
    for k, name in enumerate(("abs/G", "r_up/G", "r_pl/G", "mouse/G", "mouse plateau/G")):
        x = [t[k] for t in v if not math.isnan(t[k])]
        lines.append(f"SUMMARY {grp} {name}: {min(x):.2f}-{max(x):.2f} (n {len(x)} of {len(v)})")
open("paper_R1/criteria.out.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
