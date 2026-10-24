"""PRC pilot of the isolated model ICC-MY (out/paper_prc/prc_a*_d200.csv from run_paper_prc.sh; gisim iccprc).

Per pulse amplitude and phase: whether the pulse triggered an event (upstroke within the pulse or 300 msec after it), its peak, the
number of extra slow waves, and the lasting phase shift of the clock = (control time - perturbed time)/T0 of the last event, matched to
the nearest control event (positive = advance). A full slow wave = peak above -30 mV (control peak -14 mV).
Load series: an early full slow wave = peak > -30 mV within the pulse + 300 msec and > 0.1 T0 before the control event.
Output: paper_R1/prc_analysis.out.txt
"""
import csv, os, sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DUR, FULL = 200.0, -30.0
out = []


def load(a, path=None):
    d = {}
    for r in csv.DictReader(open(path or f"out/paper_prc/prc_a{a}_d200.csv")):
        d.setdefault(float(r["phase"]), []).append((float(r["t_ms"]), float(r["peak_mV"])))
    return d


for a in (10, 20, 50, -20, -50):
    d = load(a); ctrl = [t for t, p in d.pop(-1.0)]; T0 = (ctrl[-1] - ctrl[0]) / (len(ctrl) - 1)
    out.append(f"\n== pulse {a:+d} pA, {DUR:.0f} msec (T0 {T0:.0f} msec) ==")
    out.append("phase | triggered (peak mV) | extra waves | lasting shift (T0)")
    shifts, trig = [], []
    for ph, ev in sorted(d.items()):
        ton = ph * T0
        tr = [p for t, p in ev if ton <= t <= ton + DUR + 300]
        full = [t for t, p in ev if p > FULL]
        last = full[-1]; m = int(np.argmin([abs(last - c) for c in ctrl]))
        sh = (ctrl[m] - last) / T0; extra = len(full) - (m + 1)
        shifts.append(sh); trig.append(tr[0] if tr else None)
        out.append(f"{ph:.3f} | {'%.1f' % tr[0] if tr else '-':>6} | {extra:+d} | {sh:+.3f}")
    s = np.array(shifts); tp = [p for p in trig if p is not None]
    out.append(f"summary {a:+d} pA: lasting shift {s.min():+.3f} to {s.max():+.3f} T0; "
               f"pulse triggered an event at {len(tp)}/{len(trig)} phases"
               + (f" (peak {min(tp):.1f} to {max(tp):.1f} mV)" if tp else ""))

# ---- load series: passive conductance gL standing in for the muscle column ----
out.append("\n== load series (20 phases, 200-msec pulses): phases where the pulse triggered an early full slow wave / max lasting advance (T0) ==")
AMPS = (25, 30, 35, 40, 45, 60, 80)
out.append("gL nS | T0 msec | " + " | ".join(f"{a} pA" for a in AMPS))
for gl in ("0", "0.6", "0.91", "2"):
    row = []
    for a in AMPS:
        d = load(a, f"out/paper_prc/load_gL{gl}_a{a}.csv"); ctrl = [t for t, p in d.pop(-1.0)]; T0 = (ctrl[-1] - ctrl[0]) / (len(ctrl) - 1)
        n, adv = 0, []
        for ph, ev in d.items():
            ton = ph * T0
            if any(p > FULL and t < ctrl[0] - 0.1 * T0 for t, p in ev if ton <= t <= ton + DUR + 300): n += 1
            full = [t for t, p in ev if p > FULL]; m = int(np.argmin([abs(full[-1] - c) for c in ctrl])); adv.append((ctrl[m] - full[-1]) / T0)
        row.append(f"{n}/{len(d)}, {max(adv):+.2f}")
    out.append(f"{gl} | {T0:.0f} | " + " | ".join(row))

open("paper_R1/prc_analysis.out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(l for l in out if l.startswith("summary")))
print("\n".join(out[next(i for i, l in enumerate(out) if "load series" in l):]))
