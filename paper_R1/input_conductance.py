"""Three input conductances of the muscle column per ICC-MY (3rd review, 3-1), for the ten column geometries:
  G_load      continuum formula m*sqrt(g_SS*g_m)*tanh(L/lambda), g_m 0.96 nS (used to normalise g_c)
  Y_discrete  exact input conductance of the discrete passive ladder (layer 1 of L layers, g_m per cell, g_SS between layers)
  G_model     input conductance of the full model column at rest, measured: ICC-MY uncoupled (g_c 0), +/-10 pA, or +/-4 mV x Y_discrete for large loads (outward/inward),
              injected into muscle layer 1 for 5 sec (20-25 sec of a 30-sec run, seeds 1-3, mean); G = 2I / (V(-I) - V(+I)), V =
              mean layer-1 Vm over 22-25 sec; also with SK blocked (apamin = skKd 1000).
Usage: python paper_R1/input_conductance.py [run]    Output: out/paper_gin/*.txt, paper_R1/input_conductance.out.txt
"""
import os, sys, csv, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.chdir(os.path.dirname(HERE))
import null_model as NM
OUT = "out/paper_gin"
BASE = "fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 gK0run=0.4 iccScale=0.72".split()
GM = 0.96
SEEDS = (1, 2, 3)
# geometries (no extra keys) and the tonic K+ variants of the calibrated column (3rd review, 3-3: they change the SMC membrane)
CASES = [(n, g, []) for n, g in NM.GEO.items()] + [("tonic K 0.2", (0.3, 30, 30), ["gK0run=0.2"]), ("tonic K 0.6", (0.3, 30, 30), ["gK0run=0.6"])]


def gload(m, gss, L):
    lam = np.sqrt(gss / GM); return m * np.sqrt(gss * GM) * np.tanh(L / lam)


def sim(tag, keys):
    txt = f"{OUT}/{tag}.txt"
    if os.path.exists(txt): return float(open(txt).read())
    path = f"{OUT}/{tag}.csv"
    subprocess.run(["./bin/gisim.exe", "active", "0", "0", "0", "30000", "0.1", path] + BASE + keys, check=True, stdout=subprocess.DEVNULL)
    with open(path) as fh:
        v = [float(x["V_SMC0"]) for x in csv.DictReader(fh) if 22000 <= float(x["t_ms"]) < 25000]
    os.remove(path); val = float(np.mean(v)); open(txt, "w").write(f"{val}")
    return val


def amp(m, gss, L):
    """injected current (pA): 10 pA, or about 4 mV x the passive input conductance for large loads (FLC noise is ~0.5 mV)"""
    return max(10, int(round(4 * NM.column_admittance(m, gss, L, 0).real)))


def run_all():
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for name, (m, gss, L), var in CASES:
        geo = [f"m={m}", f"gSS={gss}", f"L={L}"] + var
        for sk, extra in (("full", []), ("nosk", ["skKd=1000"])):
            for seed in SEEDS:
                for i in (amp(m, gss, L), -amp(m, gss, L)):
                    tag = f"{name.replace(' ', '').replace(',', '_')}_{sk}_i{i}" + ("" if seed == 1 else f"_s{seed}")
                    jobs.append((tag, geo + extra + [f"smcInj={i}", "injOn=20000", "injOff=25000", f"seed={seed}"]))
    with ThreadPoolExecutor(28) as ex:
        list(ex.map(lambda j: sim(*j), jobs))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run": run_all()
    out = ["geometry | G_load formula | Y discrete passive | G_model full (rest) | G_model SK blocked | discrete/formula | full/discrete | injected pA | SD of full over seeds"]
    for name, (m, gss, L), var in CASES:
        t = name.replace(' ', '').replace(',', '_')
        G = gload(m, gss, L); Y = NM.column_admittance(m, gss, L, 0).real
        a = amp(m, gss, L)
        sfx = lambda sd: "" if sd == 1 else f"_s{sd}"
        g = {sk: np.mean([2.0 * a / (sim(f"{t}_{sk}_i{-a}{sfx(sd)}", []) - sim(f"{t}_{sk}_i{a}{sfx(sd)}", [])) for sd in SEEDS]) for sk in ("full", "nosk")}
        g["sd"] = np.std([2.0 * a / (sim(f"{t}_full_i{-a}{sfx(sd)}", []) - sim(f"{t}_full_i{a}{sfx(sd)}", [])) for sd in SEEDS])
        out.append(f"{name} | {G:.2f} | {Y:.2f} | {g['full']:.2f} | {g['nosk']:.2f} | {Y / G:.3f} | {g['full'] / Y:.2f} | {a} | {g['sd']:.2f}")
    # passive divider over the measured resting input conductance vs the full model (g_c at the target ratio, criteria.out.txt)
    full = {}
    for line in open("paper_R1/criteria.out.txt", encoding="utf-8"):
        f = [x.strip() for x in line.split("|")]
        if len(f) == 13 and f[0] in [c[0] for c in CASES]: full[f[0]] = (float(f[5]), float(f[9]))
    out.append("\npassive divider over Y (discrete, g_m 0.96) and over G_model (measured) vs full model: g_c (nS) at rabbit 0.429 | mouse 0.604")
    out.append("geometry | rabbit: Y-null / G_model-null / full | full/G_model-null | mouse: Y-null / G_model-null / full | full/G_model-null")
    q = []
    for line in out[1:len(CASES) + 1]:
        f = [x.strip() for x in line.split("|")]; n = f[0]; Y, Gm = float(f[2]), float(f[3])
        if n not in full: continue
        cells = []
        for k, r in enumerate((24.4 / 56.9, 34.1 / 56.5)):
            a_, b_, c_ = r / (1 - r) * Y, r / (1 - r) * Gm, full[n][k]
            cells.append(f"{a_:.2f} / {b_:.2f} / {c_:.2f} | {c_ / b_:.2f}"); q.append((k, c_ / a_, c_ / b_))
        out.append(f"{n} | " + " | ".join(cells))
    for k, sp in ((0, "rabbit"), (1, "mouse")):
        a_ = [x[1] for x in q if x[0] == k and x[1] == x[1]]; b_ = [x[2] for x in q if x[0] == k and x[2] == x[2]]
        out.append(f"SUMMARY {sp}: full / Y-null {min(a_):.2f}-{max(a_):.2f}; full / G_model-null {min(b_):.2f}-{max(b_):.2f}")
    open("paper_R1/input_conductance.out.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("\n".join(out))
