"""PRC of the model ICC-MY coupled to its real muscle column (check of the passive-load PRC pilot, prc_analysis.py).

gisim active (final configuration, seed 1, dt 0.1 msec, 40 sec) with a 200-msec current pulse into ICC-MY (iccInj = -amp pA, i.e.
inward/depolarizing, scaled by 1/iccScale like the coupling current). Reference = first ICC-MY upstroke (-45 mV) after 20 sec in the
control run; T0 = mean control period after it; pulse onset = reference + phase*T0. The same seed is used for every run.
Conditions: the calibrated column (m 0.3, L 30; G_load 1.61) at the rabbit and mouse g_c (1.17, 2.33 nS), and a large
load (L 7, m 1.3; G_load 5.9 nS) at g_c 3 nS. Noise check: 1-pA pulses.
Per run only the ICC-MY upstroke times and peaks are kept (out/paper_prc_col/events.csv); traces are deleted.
Usage: python paper_R1/prc_column.py run|analyze [seed]   Output: out/paper_prc_col/events_s<seed>.csv, paper_R1/prc_column_s<seed>.out.txt
"""
import csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_prc_col"
BASE = "fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 gK0run=0.4 iccScale=0.72".split()
COND = {"rabbit": (1.17, ["m=0.3"]), "mouse": (2.33, ["m=0.3"]), "large": (3.0, ["m=1.3", "L=7"])}   # revision 4: current calibrated g_c (event-aligned); earlier runs at 0.96/2.08 in out/paper_prc_col/old_gc0.96_2.08
AMPS = (1, 20, 30, 35, 40, 45, 50, 60, 70, 80)
SEEDS = (1, 2, 3)
PHASES = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
T, DT, VTH, DUR, FULL = 40000, 0.1, -45.0, 200.0, -30.0


def events(path):
    t, v = [], []
    for r in csv.DictReader(open(path)):
        t.append(float(r["t_ms"])); v.append(float(r["V_ICC"]))
    t, v = np.array(t), np.array(v)
    ev = []
    for i in np.where((v[:-1] < VTH) & (v[1:] >= VTH))[0]:
        j = i + 1
        while j < len(v) and v[j] >= VTH: j += 1
        ev.append((t[i + 1], v[i + 1:j].max() if j > i + 1 else v[i + 1]))
    return ev


def sim(name, gc, keys, tag):
    path = f"{OUT}/{name}_{tag}.csv"
    keys = keys + [f"seed={SEED}"]
    subprocess.run(["./bin/gisim.exe", "active", str(gc), "0", "0", str(T), str(DT), path] + BASE + keys,
                   check=True, stdout=subprocess.DEVNULL)
    ev = events(path); os.remove(path)
    return ev


def run():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    with ThreadPoolExecutor(3) as ex:
        ctrl = dict(zip(COND, ex.map(lambda n: sim(n, COND[n][0], COND[n][1], "ctrl"), COND)))
    jobs = []
    for n, (gc, keys) in COND.items():
        ev = [e for e in ctrl[n] if e[0] >= 20000]; ref = ev[0][0]; T0 = (ev[-1][0] - ref) / (len(ev) - 1)
        rows += [(n, -1, -1, t, p) for t, p in ctrl[n]]
        for a in AMPS:
            for ph in PHASES:
                on = ref + ph * T0
                jobs.append((n, a, ph, gc, keys + [f"iccInj={-a}", f"injOn={on:.1f}", f"injOff={on + DUR:.1f}"]))
    with ThreadPoolExecutor(28) as ex:
        res = list(ex.map(lambda j: sim(j[0], j[3], j[4], f"a{j[1]}_p{j[2]}"), jobs))
    for j, ev in zip(jobs, res):
        rows += [(j[0], j[1], j[2], t, p) for t, p in ev]
    with open(f"{OUT}/events_s{SEED}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cond", "amp_pA", "phase", "t_ms", "peak_mV"]); w.writerows(rows)


def analyze():
    d = {}
    for r in csv.DictReader(open(f"{OUT}/events_s{SEED}.csv")):
        d.setdefault((r["cond"], float(r["amp_pA"]), float(r["phase"])), []).append((float(r["t_ms"]), float(r["peak_mV"])))
    out = ["Column PRC: per condition and amplitude, phases with an early full slow wave triggered by the pulse / lasting shift range (T0)",
           "(early = peak > -30 mV within the pulse + 300 msec and > 0.1 T0 before the control event; shift > 0 = advance)"]
    for n in COND:
        c = [e for e in d[(n, -1, -1)] if e[0] >= 20000]; ref = c[0][0]; ct = [t for t, p in c[1:]]
        T0 = (c[-1][0] - ref) / (len(c) - 1)
        out.append(f"\n== {n}: g_c {COND[n][0]} nS, {' '.join(COND[n][1])}; T0 {T0:.0f} msec, control peak {np.mean([p for t, p in c]):.1f} mV ==")
        out.append("amp pA | early waves | lasting shift min..max | per phase (shift)")
        for a in AMPS:
            k, sh, per = 0, [], []
            for ph in PHASES:
                on = ref + ph * T0
                ev = [e for e in d[(n, a, ph)] if e[0] > ref + 0.05 * T0]
                if any(p > FULL and t < ct[0] - 0.1 * T0 for t, p in ev if on <= t <= on + DUR + 300): k += 1
                full = [t for t, p in ev if p > FULL]
                m = int(np.argmin([abs(full[-1] - x) for x in ct])); s = (ct[m] - full[-1]) / T0
                s = (s + 0.5) % 1.0 - 0.5                              # wrap: an extra wave near the run end can pair with the wrong control wave
                sh.append(s); per.append(f"{ph:.1f}:{s:+.2f}")
            out.append(f"{a:>3} | {k}/{len(PHASES)} | {min(sh):+.3f}..{max(sh):+.3f} | " + " ".join(per))
    open(f"paper_R1/prc_column_s{SEED}.out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    run() if sys.argv[1] == "run" else analyze()
