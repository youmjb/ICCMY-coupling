"""ICC-MY variant (icc_variant.py: VDDR x16, inactivation midpoint -80 mV, IP3 production voltage-independent) coupled to the
calibrated column (m 0.3, L 30), with a passive muscle (gT 0, SMC tonic K 0.4 nS) or the active muscle C1 (gT 2, vmT -50, vhT -65,
SMC tonic K 1 nS). g_c is swept so that the mouse control ratio 0.60 can be re-matched; drugs are sized by the CSMC resting change
only: apamin (skKd 1000, no sizing), pinacidil (SMC K_ATP, E -90 mV, gKATP 1/2/4 nS; recorded CSMC -10.3 mV).
Recorded (mouse): apamin CSMC +3.1, ICC-MY +3.3 mV, CSMC x1.00, ICC-MY x0.96; pinacidil CSMC -10.3, ICC-MY -9.6 mV,
CSMC x~0.2, ICC-MY x~1.16.
Usage: python paper_R1/coupled_var.py run | analyze      Runs: out/paper_active/cv_<muscle>_gc<gc>_<cond>.csv (20 s)"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import B, OUT, summary
ICCV = ["iccVDDR=16", "iccVH=-80", "ip3VK=0"]
MUS = {"pas": ["gK0run=0.4"], "act": ["gT=2", "vmT=-50", "vhT=-65", "gK0run=1"]}
GC = ["1", "1.5", "2.33", "3", "4"]
COND = {"ctrl": [], "apa": ["skKd=1000"], "pin1": ["gKATP=1"], "pin2": ["gKATP=2"], "pin4": ["gKATP=4"]}
JOBS = list(itertools.product(MUS, GC, COND))


def name(mu, gc, c):
    return f"{OUT}/cv_{mu}_gc{gc}_{c}"


def run(j):
    mu, gc, c = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + B + ICCV + MUS[mu] + COND[c],
                       check=True, stdout=lg)


def safe(f):
    try:
        return summary(f)
    except ValueError:
        return None


if __name__ == "__main__":
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(25) as ex:
            list(ex.map(run, JOBS))
    else:
        for mu in MUS:
            for gc in GC:
                s0 = safe(name(mu, gc, "ctrl") + ".csv"); a0 = s0 and s0["a"]
                if not a0: print(mu, gc, "control: no slow waves or diverged"); continue
                print(f"{mu} g_c {gc}: ctrl CSMC rest {s0['sr']:.1f} ICC rest {s0['ir']:.1f} | ICC {a0['icc']:.1f} CSMC {a0['smc']:.1f} "
                      f"ratio {a0['smc'] / a0['icc']:.2f} dV/dt {s0['dv']:.2f}")
                for c in list(COND)[1:]:
                    s = safe(name(mu, gc, c) + ".csv"); a = s and s["a"]
                    if not a: print(f"     {c:5s} no slow waves or diverged"); continue
                    dS, dI = s["sr"] - s0["sr"], s["ir"] - s0["ir"]
                    print(f"     {c:5s} dCSMC {dS:+5.1f} dICC {dI:+5.1f} (transfer {dI / dS:.2f}) | ICC x{a['icc'] / a0['icc']:.2f} "
                          f"CSMC x{a['smc'] / a0['smc']:.2f} | ratio {a['smc'] / a['icc']:.2f}")
        print("recorded: apamin +3.1/+3.3 (transfer ~1.06), ICC x0.96, CSMC x1.00 | pinacidil -10.3/-9.6 (0.93), ICC x~1.16, CSMC x~0.2 | "
              "control ratio 0.60")
