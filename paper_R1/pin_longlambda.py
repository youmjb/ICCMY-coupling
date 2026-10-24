"""Hypothesis (c): the pinacidil fall of the muscle slow wave (x~0.2 for -10.3 mV, Kito, Ward & Sanders 2005) is passive if
(1) E_K is about -84 mV (Krebs with 5.9 mM K+), so that -10 mV needs a K_ATP conductance several times the resting one, and
(2) the transmural length constant is long (as the depth profile of Suzuki, Prosser & Dahms 1986 suggests), so that the column load
grows with the membrane conductance itself rather than with its square root.
Passive muscle (gT 0, SMC tonic K 0.4 nS), ICC-MY variant of icc_variant.py, pinacidil = SMC K_ATP with E -84 mV (gKATP 1-16 nS).
Columns: 'cal' m 0.3, g_SS 30, L 30 (paper); 'thin300' and 'thin3000' m 0.3, L 7, g_SS 300 / 3000 nS (long lambda, mouse-like).
g_c is swept so that each column can be matched to the mouse control ratio 0.60.
Usage: [NIF=1] python paper_R1/pin_longlambda.py run | analyze     Runs: out/paper_active/ll_<col>_gc<gc>_<cond>.csv (20 s)"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import OUT, summary
BASE = "fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1 kmT=4 gK0run=0.4 iccVDDR=16 iccVH=-80 ip3VK=0 EKATP=-84".split()
COL = {"cal": ["m=0.3", "gSS=30", "L=30"], "thin300": ["m=0.3", "gSS=300", "L=7"], "thin3000": ["m=0.3", "gSS=3000", "L=7"]}
NIF = os.environ.get("NIF", "0") == "1"          # 1: SMC L-type blocked (env SMS_CAL=0), as in the nifedipine recordings
GC = ["1", "2", "3", "4", "6"] if not NIF else ["1", "1.5", "2", "3", "4", "6"]
COND = {"ctrl": []}
COND.update({f"pin{g}": [f"gKATP={g}"] for g in ("1", "2", "4", "8", "16")})
if NIF:
    COL = {k: v for k, v in COL.items() if k != "cal"}
    COND = {k: v for k, v in COND.items() if k in ("ctrl", "pin1", "pin2", "pin4")}
JOBS = list(itertools.product(COL, GC, COND))


def name(col, gc, c):
    return f"{OUT}/ll_{col}_gc{gc}_{c}" + ("_nif" if NIF else "")


def run(j):
    col, gc, c = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    env = dict(os.environ, SMS_CAL="0") if NIF else None
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + BASE + COL[col] + COND[c],
                       check=True, stdout=lg, env=env)


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
        for col in COL:
            for gc in GC:
                s0 = safe(name(col, gc, "ctrl") + ".csv"); a0 = s0 and s0["a"]
                if not a0: print(col, gc, "control: no slow waves or diverged"); continue
                print(f"{col} g_c {gc}: ctrl CSMC rest {s0['sr']:.1f} ICC rest {s0['ir']:.1f} | ICC {a0['icc']:.1f} CSMC {a0['smc']:.1f} "
                      f"ratio {a0['smc'] / a0['icc']:.2f}")
                for c in list(COND)[1:]:
                    s = safe(name(col, gc, c) + ".csv"); a = s and s["a"]
                    if not a: print(f"     {c:6s} no slow waves or diverged"); continue
                    dS, dI = s["sr"] - s0["sr"], s["ir"] - s0["ir"]
                    print(f"     {c:6s} dCSMC {dS:+5.1f} dICC {dI:+5.1f} (transfer {dI / dS:.2f}) | ICC x{a['icc'] / a0['icc']:.2f} "
                          f"CSMC x{a['smc'] / a0['smc']:.2f} | ratio {a['smc'] / a['icc']:.2f}")
        print("recorded pinacidil: CSMC -10.3, ICC-MY -9.6 mV (transfer 0.93); ICC-MY x~1.16, CSMC x~0.2; control ratio 0.60")
