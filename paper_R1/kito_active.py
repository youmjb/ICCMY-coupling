"""Small test: can an active, dihydropyridine-resistant inward current in the SMC (T-type-like, gisim keys gT, vmT, ...) reproduce the
pinacidil result of Kito, Ward & Sanders 2005 (CSMC slow wave x0.2 for -10 mV, ICC-MY ~x1.16) that the passive column misses?
Grid: g_c x gT, control (gK0run 0.4) and pinacidil (gK0run 4, K only in SMC = scenario A of kito_pinacidil.py).
Usage: [VHT=-60] python paper_R1/kito_active.py run | analyze      Runs: out/paper_kito/act_gc<gc>_T<gT>_k<gK>.csv"""
import os, sys, csv, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
OUT = "out/paper_kito"
B = "m=0.3 fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1".split()
VH = os.environ.get("VHT", "-70")                # T inactivation midpoint (mV); run 2 uses -60
GC = ["0.5", "1", "1.5", "2.33"] if VH == "-70" else ["0.25", "0.5", "1"]
GT = ["0", "0.5", "1", "2"] if VH == "-70" else ["4", "8"]
GK = ["0.4", "4"]
JOBS = [(gc, gt, gk) for gc in GC for gt in GT for gk in GK]


def name(gc, gt, gk):
    return f"{OUT}/act_gc{gc}_T{gt}_k{gk}" + ("" if VH == "-70" else f"_vh{VH}")


def run(j):
    gc, gt, gk = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", f + ".csv"] + B
                       + [f"gK0run={gk}", f"gT={gt}", f"vhT={VH}"], check=True, stdout=lg)


def rest(f):
    r = list(csv.DictReader(open(f)))
    t, vi, vs = (np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0"))
    h = t >= 20000; vi, vs = vi[h], vs[h]
    near = np.convolve((vi >= -65).astype(float), np.ones(601), mode="same") > 0
    return np.percentile(vi, 5), (np.median(vs[~near]) if (~near).any() else np.percentile(vs, 5))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(16) as ex:
            list(ex.map(run, JOBS))
    else:
        import amp_event as AE
        print("g_c | gT | control: ICC rest, CSMC rest, ICC amp, CSMC amp, ratio | pinacidil: dCSMC rest, dICC rest, ICC amp x, CSMC amp x, ratio | ratio change")
        for gc in GC:
            for gt in GT:
                c, p = name(gc, gt, "0.4"), name(gc, gt, "4")
                if not (os.path.exists(c + ".csv") and os.path.exists(p + ".csv")): continue
                rc, rp = rest(c + ".csv"), rest(p + ".csv"); ac, ap = AE.measure(c + ".csv"), AE.measure(p + ".csv")
                if ac is None or ap is None: print(gc, gt, "no slow waves"); continue
                q0, q1 = ac["smc"] / ac["icc"], ap["smc"] / ap["icc"]
                print(f"{gc:>4} | {gt:>3} | {rc[0]:6.1f} {rc[1]:6.1f} {ac['icc']:5.1f} {ac['smc']:5.1f} {q0:.2f} | "
                      f"{rp[1] - rc[1]:+5.1f} {rp[0] - rc[0]:+5.1f} {ap['icc'] / ac['icc']:.2f} {ap['smc'] / ac['smc']:.2f} {q1:.2f} | x{q1 / q0:.2f}")
        print("Kito 2005: CSMC -10.3 mV, ICC -9.6 mV; ICC x~1.16, CSMC x~0.2; ratio ~0.50 -> ~0.09 (x~0.18). Mouse control ratio target 0.60.")
