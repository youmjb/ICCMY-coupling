"""Example ICC-MY traces for Fig. S3A (replaces run_prc_traces.sh; revision 4 uses the current rabbit g_c).
Column g_c 1.17 nS (m 0.3, L 30), seed 1. The reference upstroke and period T0 come from the control run of prc_column.py
(out/paper_prc_col/events_s1.csv, condition "rabbit"); 200-msec pulses at phase 0.7: 30 pA (below threshold) and 40 pA (above).
Writes out/paper_prc_col/trace_{ctrl,a30,a40}.csv and trace_onset.txt (pulse onset, msec) read by fig_s3_prc.py.
"""
import os, csv, subprocess
from concurrent.futures import ThreadPoolExecutor

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GC, PHASE = 1.17, 0.7
B = "m=0.3 fI=0 fP=0.1 gSP=3 pdModel=1 smcModel=1 gK0run=0.4 iccScale=0.72 seed=1".split()
ev = [float(r["t_ms"]) for r in csv.DictReader(open("out/paper_prc_col/events_s1.csv"))
      if r["cond"] == "rabbit" and r["amp_pA"] == "-1" and float(r["t_ms"]) >= 20000]
ref, T0 = ev[0], (ev[-1] - ev[0]) / (len(ev) - 1)
on = round(ref + PHASE * T0)
open("out/paper_prc_col/trace_onset.txt", "w").write(f"{on}\n")
jobs = [("ctrl", []), ("a30", ["iccInj=-30", f"injOn={on}", f"injOff={on + 200}"]), ("a40", ["iccInj=-40", f"injOn={on}", f"injOff={on + 200}"])]
with ThreadPoolExecutor(3) as ex:
    list(ex.map(lambda j: subprocess.run(["./bin/gisim.exe", "active", str(GC), "0", "0", "40000", "0.1",
                                          f"out/paper_prc_col/trace_{j[0]}.csv"] + B + j[1], check=True, stdout=subprocess.DEVNULL), jobs))
print(f"reference {ref:.0f} msec, T0 {T0:.0f} msec, pulse onset {on} msec")
