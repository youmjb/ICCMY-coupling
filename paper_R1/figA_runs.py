"""Runs for the plan-A figures under one set of conditions (ICC-MY variant, nifedipine SMS_CAL=0, K_ATP E -84 mV, 40 s, seed 1).
(1) Trade-off families for pinacidil (control + gKATP 1/2/3/4 nS, the -10.3 mV level found by interpolation):
    gc   : the final model column (m 0.3, g_SS 30, L 30) at g_c 1.75-8 nS
    lm   : the same plus a longitudinal-muscle column (LL 7, mL 0.3, g_SSL 30) coupled by gcL = g_c or 2 g_c
    lam  : long transmural lambda (L 7, g_SS 3000)
    act  : active muscle C1 (gT 2, vmT -50, vhT -65, kmT 4, SMC tonic K 1 nS)
(2) Isolated ICC-MY (g_c 0) hyperpolarized by outward current, Youm 2019 vs the variant, for Fig. 5.
Usage: python paper_R1/figA_runs.py run      Runs: out/paper_figA/"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_figA"
ICCV = ["iccVDDR=16", "iccVH=-80", "ip3VK=0"]
COM = "fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 iccScale=0.72 seed=1 EKATP=-84".split()
COL = ["m=0.3", "gSS=30", "L=30"]
FAM = {}
for gc in ("1.75", "2.75", "4", "6", "8"): FAM[f"gc_{gc}"] = (gc, COL + ["gK0run=0.4"])
for gc in ("1.75", "2.75", "4"):
    for r in ("1", "2"): FAM[f"lm{r}_{gc}"] = (gc, COL + ["gK0run=0.4", "LL=7", "mL=0.3", "gSSL=30", f"gcL={float(gc) * float(r)}"])
for gc in ("1.75", "2.75", "4"): FAM[f"lam_{gc}"] = (gc, ["m=0.3", "gSS=3000", "L=7", "gK0run=0.4"])
for gc in ("1", "1.5", "2"): FAM[f"act_{gc}"] = (gc, COL + ["gT=2", "vmT=-50", "vhT=-65", "kmT=4", "gK0run=1"])
PIN = {"ctrl": [], "p1": ["gKATP=1"], "p2": ["gKATP=2"], "p3": ["gKATP=3"], "p4": ["gKATP=4"], "p6": ["gKATP=6"], "p8": ["gKATP=8"]}
ISO = {"orig": [], "var": ICCV}
INJ = ["0", "5", "10", "15", "25", "35", "50"]


def job_fam(k, c):
    gc, keys = FAM[k]; f = f"{OUT}/{k}_{c}"
    return f, ["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", f + ".csv"] + COM + ICCV + keys + PIN[c]


def job_iso(m, i):
    f = f"{OUT}/iso_{m}_i{i}"
    return f, ["./bin/gisim.exe", "active", "0", "0", "0", "40000", "0.1", f + ".csv", "iccScale=0.72", "L=2", "m=0.3", "fI=0", "fP=0",
               "smcModel=1", "seed=1", f"iccInj={i}", "injOn=0", "injOff=1e9"] + ISO[m]


def run(fc):
    f, cmd = fc
    if os.path.exists(f + ".log"): return
    with open(f + ".log", "w") as lg:
        subprocess.run(cmd, check=True, stdout=lg, env=dict(os.environ, SMS_CAL="0", OMP_NUM_THREADS="1"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    jobs = [job_fam(k, c) for k in FAM for c in PIN if c not in ("p6", "p8") or k.startswith("act")] + [job_iso(m, i) for m in ISO for i in INJ]
    print(len(jobs), "runs")
    with ThreadPoolExecutor(30) as ex:
        list(ex.map(run, jobs))
