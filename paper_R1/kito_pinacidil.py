"""Model vs Kito, Ward & Sanders 2005 (mouse small intestine, nifedipine): pinacidil 10 uM (K_ATP opener) hyperpolarized CSMC by
10.3 +/- 1.0 mV (n 11) and ICC-MY by 9.6 +/- 2.0 mV (n 8); CSMC slow wave about 28 -> 6 mV, ICC-MY pacemaker potential about
56 -> 65 mV (Fig. 9A, read from the bar chart; unpaired tissues). The drug is modelled as a K+ conductance added to the SMC
(tonic gK0run, E -80 mV) and, in scenario B, also to ICC-MY (iccGL, iccEL -90 mV); its size is chosen from the resting potential
only, and the amplitudes are the prediction.
Usage: python paper_R1/kito_pinacidil.py run | analyze
Runs: out/paper_kito/pin_gc<gc>_k<gK>_i<gI>.csv (calibrated column m 0.3, seed 1, 40 s, dt 0.1 ms)."""
import os, sys, csv, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
OUT = "out/paper_kito"
B = "m=0.3 fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1".split()
GC = ["2.33", "1.17"]
GK = ["0.4", "1", "2", "4", "8"]                 # SMC tonic K (nS); 0.4 = control
GI = ["0", "1", "3"]                             # extra ICC-MY K conductance (nS, E -90 mV); 0 = scenario A
JOBS = [(gc, gk, gi) for gc in GC for gk in GK for gi in GI if not (gk == "0.4" and gi != "0")]


def name(gc, gk, gi):
    return f"{OUT}/pin_gc{gc}_k{gk}_i{gi}"


def run(j):
    gc, gk, gi = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", f + ".csv"] + B
                       + [f"gK0run={gk}", f"iccGL={gi}", "iccEL=-90"], check=True, stdout=lg)


def rest(f):
    r = list(csv.DictReader(open(f)))
    t, vi, vs = (np.array([float(x[k]) for x in r]) for k in ("t_ms", "V_ICC", "V_SMC0"))
    h = t >= 20000; vi, vs = vi[h], vs[h]
    near = np.convolve((vi >= -65).astype(float), np.ones(601), mode="same") > 0
    return np.percentile(vi, 5), np.median(vs[~near]) if (~near).any() else np.percentile(vs, 5)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(16) as ex:
            list(ex.map(run, JOBS))
    else:
        import amp_event as AE
        print("g_c | gK SMC | gK ICC | ICC rest (d) | CSMC rest (d) | ICC amp | CSMC amp | ratio")
        for gc in GC:
            r0 = rest(name(gc, "0.4", "0") + ".csv"); a0 = AE.measure(name(gc, "0.4", "0") + ".csv")
            for gk in GK:
                for gi in GI:
                    f = name(gc, gk, gi)
                    if not os.path.exists(f + ".csv"): continue
                    r = rest(f + ".csv"); a = AE.measure(f + ".csv")
                    if a is None: print(gc, gk, gi, "no slow waves"); continue
                    print(f"{gc} | {gk:>3} | {gi} | {r[0]:6.1f} ({r[0] - r0[0]:+5.1f}) | {r[1]:6.1f} ({r[1] - r0[1]:+5.1f}) | "
                          f"{a['icc']:5.1f} ({a['icc'] / a0['icc']:.2f}x) | {a['smc']:5.1f} ({a['smc'] / a0['smc']:.2f}x) | {a['smc'] / a['icc']:.2f}")
        print("Kito 2005: CSMC -10.3 mV, slow wave ~28 -> ~6 mV (0.2x); ICC-MY -9.6 mV, ~56 -> ~65 mV (1.16x); ratio ~0.50 -> ~0.09")
