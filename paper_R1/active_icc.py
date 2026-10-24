"""Step 1-2 after active_drugs.py: active-SMC candidate C1 (gT 2, vmT -50, vhT -65, SMC tonic K 1 nS) with
(1) pinacidil as a separate SMC K_ATP conductance with E -90 mV (gKATP 2/4/8 nS; the recorded CSMC change -10.3 mV picks the size),
(2) a reduced ICC-MY resting conductance (iccGL 0/-1/-2/-3 nS at iccEL -72 mV, near ICC-MY rest), testing whether the near-complete
transfer of bath-drug resting changes from muscle to ICC-MY (apamin: CSMC +3.1, ICC-MY +3.3 mV; pinacidil: -10.3 / -9.6 mV) needs a
smaller ICC-MY own conductance relative to g_c. g_c 1/1.5/2 nS so that the control ratio can be re-matched.
Usage: python paper_R1/active_icc.py run | analyze      Runs: out/paper_active/icc_gc<gc>_L<iccGL>_<cond>.csv (20 s)"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import B, OUT, summary
C1 = ["gT=2", "vmT=-50", "vhT=-65", "gK0run=1", "iccEL=-72"]
GC = ["1", "1.5", "2"]
GL = ["0", "-1", "-2", "-3"]
COND = {"ctrl": [], "apa": ["skKd=1000"], "pin2": ["gKATP=2"], "pin4": ["gKATP=4"], "pin8": ["gKATP=8"]}
JOBS = list(itertools.product(GC, GL, COND))


def name(gc, gl, c):
    return f"{OUT}/icc_gc{gc}_L{gl}_{c}"


def run(j):
    gc, gl, c = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + B + C1 + [f"iccGL={gl}"] + COND[c],
                       check=True, stdout=lg)


def safe(f):
    try:
        return summary(f)
    except ValueError:                                     # NaN in the trace: the run diverged
        return None


if __name__ == "__main__":
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(24) as ex:
            list(ex.map(run, JOBS))
    else:
        for gc in GC:
            for gl in GL:
                s0 = safe(name(gc, gl, "ctrl") + ".csv"); a0 = s0 and s0["a"]
                if not a0: print(gc, gl, "control: no slow waves or diverged"); continue
                print(f"g_c {gc} iccGL {gl}: ctrl CSMC rest {s0['sr']:.1f} ICC rest {s0['ir']:.1f} | ICC {a0['icc']:.1f} CSMC {a0['smc']:.1f} "
                      f"ratio {a0['smc'] / a0['icc']:.2f} dV/dt {s0['dv']:.2f}")
                for c in list(COND)[1:]:
                    f = name(gc, gl, c) + ".csv"
                    if not os.path.exists(f): continue
                    s = safe(f); a = s and s["a"]
                    if not a: print(f"     {c:5s} no slow waves or diverged"); continue
                    dS, dI = s["sr"] - s0["sr"], s["ir"] - s0["ir"]
                    print(f"     {c:5s} dCSMC {dS:+5.1f} dICC {dI:+5.1f} (transfer {dI / dS:.2f}) | ICC x{a['icc'] / a0['icc']:.2f} "
                          f"CSMC x{a['smc'] / a0['smc']:.2f} dV/dt x{s['dv'] / s0['dv']:.2f}")
        print("recorded: apamin dCSMC +3.1 dICC +3.3 (transfer ~1.06), ICC x0.96, CSMC x1.00 | "
              "pinacidil dCSMC -10.3 dICC -9.6 (0.93), ICC x~1.16, CSMC x~0.2 | control ratio 0.60")
