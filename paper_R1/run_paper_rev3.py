"""Revision 3 runs (2026-10-16): ICC loss and transmural decay. Final configuration (fI 0, fP 0.1, gSP 3, pdModel 1, smcModel 1,
gK0run 0.4, iccScale 0.72, seed 1), 40 sec, dt 0.1 msec, no neural input.
(1) ICC loss: a fraction f of the ICC-MY remains, uniformly, so each remaining ICC-MY serves m = 0.3/f muscle cells per layer
    (L 30, g_SS 30), with g_c per ICC-MY kept at the calibrated value (rabbit 1.17 nS, mouse 2.33 nS)
    -> out/paper_rev3/loss_gc<gc>_f<f>.csv/.log
(2) Transmural decay: every layer written (traceAll=1) for the calibrated column (L 30, m 0.3) and the mouse-like column
    (L 7, m 1.3) at their rabbit and mouse g_c -> out/paper_rev3/tm_<name>.csv/.log (+ .layers.csv)
Usage: python paper_R1/run_paper_rev3.py
"""
import os, subprocess
from concurrent.futures import ThreadPoolExecutor

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_rev3"
BASE = "fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 gK0run=0.4 iccScale=0.72 seed=1".split()
FRAC = (1.0, 0.8, 2 / 3, 0.5, 0.4, 0.3, 0.2, 0.1)
GC = {"rabbit": 1.17, "mouse": 2.33}
TM = [("L30_rabbit", 1.17, ["m=0.3", "L=30", "gSS=30"]), ("L30_mouse", 2.33, ["m=0.3", "L=30", "gSS=30"]),
      ("L7_rabbit", 4.22, ["m=1.3", "L=7", "gSS=30"]), ("L7_mouse", 8.26, ["m=1.3", "L=7", "gSS=30"])]


def run(tag, gc, keys):
    pre = f"{OUT}/{tag}"
    if os.path.exists(pre + ".log") and os.path.getsize(pre + ".log"): return
    with open(pre + ".log", "w") as log:
        subprocess.run(["./bin/gisim.exe", "active", str(gc), "0", "0", "40000", "0.1", pre + ".csv"] + BASE + keys, check=True, stdout=log)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    jobs = [(f"loss_gc{gc}_f{f:.3f}", gc, [f"m={0.3 / f:.4f}", "L=30", "gSS=30"]) for gc in GC.values() for f in FRAC]
    jobs += [(f"tm_{n}", gc, k + ["traceAll=1"]) for n, gc, k in TM]
    with ThreadPoolExecutor(24) as ex:
        list(ex.map(lambda j: run(*j), jobs))
    print("done", len(jobs))
