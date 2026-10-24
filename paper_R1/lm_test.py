"""Longitudinal muscle on the other side of ICC-MY (gisim keys LL, gcL, mL, gSSL): does coupling ICC-MY to both muscle layers,
which a bath drug affects together, give the near-complete transfer of resting changes from muscle to ICC-MY seen in situ
(apamin +3.1 -> +3.3 mV; pinacidil -10.3 -> -9.6 mV) and the larger fall of the CSMC slow wave (x~0.2)?
Cousins et al. (2003): coupling of the ICC-MY network to the longitudinal layer about equal to that to the circular layer.
Set-up as the best case of pin_longlambda.py: passive SMC, L-type blocked (SMS_CAL=0, nifedipine), ICC-MY variant, K_ATP E -84 mV,
circular column m 0.3, L 7, g_SS 300; LM column LL 7, mL 0.3, g_SSL 300, gcL = 0, 1 or 2 x g_c. The LM has no FLC, so apamin acts on
the circular layer only.
Usage: python paper_R1/lm_test.py run | analyze      Runs: out/paper_active/lm_gc<gc>_r<gcL/gc>_<cond>.csv (20 s)"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import OUT, summary
BASE = ("fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1 gK0run=0.4 iccVDDR=16 iccVH=-80 ip3VK=0 EKATP=-84 "
        "m=0.3 gSS=300 L=7 LL=7 mL=0.3 gSSL=300").split()
GC = ["1.5", "2", "3", "4"]
R = ["0", "1", "2"]
COND = {"ctrl": [], "apa": ["skKd=1000"], "pin1": ["gKATP=1"], "pin2": ["gKATP=2"], "pin4": ["gKATP=4"]}
JOBS = list(itertools.product(GC, R, COND))


def name(gc, r, c):
    return f"{OUT}/lm_gc{gc}_r{r}_{c}"


def run(j):
    gc, r, c = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + BASE + [f"gcL={float(gc) * float(r)}"]
                       + COND[c], check=True, stdout=lg, env=dict(os.environ, SMS_CAL="0"))


def safe(f):
    try:
        return summary(f)
    except ValueError:
        return None


if __name__ == "__main__":
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(30) as ex:
            list(ex.map(run, JOBS))
    else:
        for gc in GC:
            for r in R:
                s0 = safe(name(gc, r, "ctrl") + ".csv"); a0 = s0 and s0["a"]
                if not a0: print(gc, r, "control: no slow waves or diverged"); continue
                print(f"g_c {gc} gcL {r}x: ctrl CSMC rest {s0['sr']:.1f} ICC rest {s0['ir']:.1f} | ICC {a0['icc']:.1f} CSMC {a0['smc']:.1f} "
                      f"ratio {a0['smc'] / a0['icc']:.2f}")
                for c in list(COND)[1:]:
                    s = safe(name(gc, r, c) + ".csv"); a = s and s["a"]
                    if not a: print(f"     {c:5s} no slow waves or diverged"); continue
                    dS, dI = s["sr"] - s0["sr"], s["ir"] - s0["ir"]
                    print(f"     {c:5s} dCSMC {dS:+5.1f} dICC {dI:+5.1f} (transfer {dI / dS:.2f}) | ICC x{a['icc'] / a0['icc']:.2f} "
                          f"CSMC x{a['smc'] / a0['smc']:.2f} | ratio {a['smc'] / a['icc']:.2f}")
        print("recorded: apamin +3.1/+3.3 (transfer ~1.06), ICC x0.96, CSMC x1.00 | pinacidil -10.3/-9.6 (0.93), ICC x~1.16, CSMC x~0.2 | ratio 0.60")
