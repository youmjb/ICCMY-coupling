"""Event-aligned amplitudes (3rd review, 3-2): does the g_c estimate depend on measuring amplitude as max - min over the window?

For every saved run (same cases and runs as analyze_criteria.py, plus the SK-blocked calibrated runs) and every ICC-MY upstroke
(mid-level crossing, as amp_event.py) in the second half of the run:
  baseline = mean Vm over 50-200 msec before the ICC-MY upstroke (ICC-MY and muscle layer 1)
  peak     = maximum within 2 sec after the upstroke
  trough   = minimum between this upstroke and the next (muscle layer 1)
  event amplitude = peak - baseline (mean over events); window amplitude = max - min over the half (as in the run logs)
g_c at each target ratio is read by log-interpolation as in analyze_criteria.py, with the event-aligned ratio
(mean layer-1 event amplitude / mean ICC-MY event amplitude). Targets: rabbit 24.4/56.9, mouse primary 34.1/56.5, mouse plateau
31.1/56.5 (Kito & Suzuki 2003). SK contribution at the calibrated column recomputed with both measures.
Output: paper_R1/event_amplitude.out.txt
"""
import os, sys, io, csv, math, contextlib
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import analyze_criteria as AC
import amp_event as AE
os.chdir(os.path.dirname(HERE))
TARGETS = (("rabbit", 24.4 / 56.9), ("mouse", 34.1 / 56.5), ("mouse plateau", 31.1 / 56.5))
APA = [(1, "out/paper_cal/apa_gc1"), (1.5, "out/paper_cal/apa_gc1.5"), (2, "out/paper_cal/apa_gc2"),
       (2.5, "out/paper_rev1/apa_gc2.5"), (3, "out/paper_rev1/apa_gc3")]
CACHE = {}


def measure(pre):
    """window and event-aligned amplitudes of ICC-MY and layer 1 for one run (event detection as in amp_event.py)"""
    m = AE.measure(pre + ".csv")
    return dict(win_i=m["win_icc"], win_s=m["win_smc"], ev_i=m["icc"], ev_s=m["smc"], n=m["n"], base=m["base"], peak=m["peak"],
                trough=m["trough"])


def gc_at(gcs, rat, target):
    return AC.cross(gcs, rat, target)


out = ["Event-aligned vs window (max - min) amplitudes; g_c (nS) at each target ratio",
       "case | G_load | " + " | ".join(f"{n} window / event" for n, _ in TARGETS) + " | layer-1 amp at g_c 1 nS window / event (mV)"]
summ = {}
for lab, m, gss, L, pre, gcs, grp in AC.CASES:
    runs = [g for g in gcs if os.path.exists(pre(g) + ".csv")]
    if len(runs) < 2: out.append(f"{lab} | traces missing"); continue
    ms = [measure(pre(g)) for g in runs]; G = AC.gload(m, gss, L)
    rw = [x["win_s"] / x["win_i"] for x in ms]; re_ = [x["ev_s"] / x["ev_i"] for x in ms]
    cells = []
    for n, tr in TARGETS:
        a, b = gc_at(runs, rw, tr), gc_at(runs, re_, tr)
        cells.append(f"{a:.2f} / {b:.2f}"); summ.setdefault((grp, n), []).append((a / G, b / G, b / a))
    one = next((x for g, x in zip(runs, ms) if g == 1), None)
    out.append(f"{lab} | {G:.2f} | " + " | ".join(cells) + (f" | {one['win_s']:.1f} / {one['ev_s']:.1f}" if one else " | -"))
for (grp, n), v in summ.items():
    w = [x[0] for x in v if x[0] == x[0]]; e = [x[1] for x in v if x[1] == x[1]]; q = [x[2] for x in v if x[2] == x[2]]
    out.append(f"SUMMARY {grp} {n}: g_c/G_load window {min(w):.2f}-{max(w):.2f}, event {min(e):.2f}-{max(e):.2f}; event/window g_c {min(q):.2f}-{max(q):.2f}")

# ---- calibrated column with and without SK ----
out.append("\nCalibrated column (m 0.3, g_SS 30, L 30): layer 1 per g_c, full vs SK blocked")
out.append("g_c | condition | window amp | event amp | baseline | peak | most negative trough (mV) | ICC-MY window / event")
full = {g: measure(f"out/v4win/m0.3_gSS30_gc{g}_u0") for g in (0.5, 1, 1.5, 2, 3)}
nosk = {g: measure(f) for g, f in APA if os.path.exists(f + ".csv")}
for g in sorted(set(full) | set(nosk)):
    for nm, d in (("full", full), ("SK blocked", nosk)):
        if g in d:
            x = d[g]
            out.append(f"{g} | {nm} | {x['win_s']:.1f} | {x['ev_s']:.1f} | {x['base']:.1f} | {x['peak']:.1f} | {x['trough']:.1f} | {x['win_i']:.1f} / {x['ev_i']:.1f}")
import null_model as NM
for n, tr in TARGETS:
    Yd = NM.column_admittance(0.3, 30, 30, 0).real; dc = tr / (1 - tr) * Yd
    row = []
    for meas in (("win_s", "win_i"), ("ev_s", "ev_i")):
        gf = sorted(full); gn = sorted(nosk)
        f_ = gc_at(gf, [full[g][meas[0]] / full[g][meas[1]] for g in gf], tr)
        n_ = gc_at(gn, [nosk[g][meas[0]] / nosk[g][meas[1]] for g in gn], tr)
        row.append(f"{'window' if meas[0] == 'win_s' else 'event'}: passive {dc:.2f} -> no SK {n_:.2f} -> full {f_:.2f} nS, "
                   f"full/passive {f_ / dc:.2f}, SK share {(n_ - f_) / (dc - f_):.0%}")
    out.append(f"{n} (ratio {tr:.3f}): " + " || ".join(row))
open("paper_R1/event_amplitude.out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out))
