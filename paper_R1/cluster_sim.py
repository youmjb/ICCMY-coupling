"""Could the 'single' cultured ICC of Koh, Sanders & Ward 1998 (29 pF, 250 MOhm, -53 mV, pacing) have been small clusters of coupled
ICC with a low own conductance? N ICC of 5.5 pF each (Zhu et al. 2009 single-cell capacitance; s 0.22), coupled as a chain or all-to-all
by gII, with the upstroke/clock variant (VDDR x16, midpoint -80 mV, IP3 voltage-independent) and resting currents x f; no muscle.
Measurements as an experimenter would make them with an electrode in cell 0:
  cc : current clamp, I = 0 (30 s): resting potential (5th percentile), slow-wave amplitude and rate (last 15 s)
  inj: +5 pA outward from 0 s: input resistance from the shift of the rest
  vc : voltage clamp (Rs 10 MOhm, Vh -70 mV, -5 mV step at 100 ms for 100 ms): capacitance from the charge in the first 1 ms
       (what a fast transient fit sees) and in the first 50 ms (all cells), and input resistance from the steady current.
Usage: python paper_R1/cluster_sim.py run | analyze      Runs: out/paper_cluster/"""
import os, sys, csv, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_cluster"
VAR = ["s=0.22", "iccVDDR=16", "iccVH=-80", "ip3VK=0"]
F = ["1", "0.3", "0.2", "0.1"]
CFG = [("1", "0", "1")] + [(n, topo, g) for n in ("3", "6") for topo in ("0", "1") for g in ("0.3", "1", "3", "10")]
PROT = {"cc": ("30000", "0.02", []), "inj": ("30000", "0.02", ["inj=5", "injOn=0", "injOff=1e9"]),
        "vc": ("250", "0.002", ["vc=1", "smp=0.01", "Rs=10", "Vh=-70", "dV=-5", "tOn=100", "tDur=100"])}
JOBS = [(f, n, topo, g, p) for f in F for n, topo, g in CFG for p in PROT]


def name(f, n, topo, g, p): return f"{OUT}/f{f}_N{n}_t{topo}_g{g}_{p}.csv"


def run(j):
    f, n, topo, g, p = j; out = name(*j); T, dt, extra = PROT[p]
    if os.path.exists(out + ".log"): return
    with open(out + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "cluster", T, dt, out] + VAR + [f"iccRest={f}", f"N={n}", f"topo={topo}", f"gII={g}"] + extra, stdout=lg)


def load(path):
    r = np.loadtxt(path, delimiter=",", skiprows=1)
    return r[:, 0], r[:, 1], r[:, 3]


def cc_stats(path):
    t, v, _ = load(path); h = t >= 15000; t, v = t[h], v[h]
    rest, peak = np.percentile(v, 5), np.percentile(v, 99.5)
    mid = rest + 0.5 * (peak - rest); up = np.where((v[:-1] < mid) & (v[1:] >= mid))[0]
    paced = peak - rest > 10 and len(up) > 1
    return rest, (peak - rest if paced else 0.0), (60000.0 / np.mean(np.diff(t[up])) if paced else 0.0)


def vc_stats(path, dV=-5.0):
    t, v, i = load(path)
    base = i[(t > 80) & (t < 100)].mean(); ss = i[(t > 180) & (t < 200)].mean()
    def charge(tw):
        w = (t >= 100) & (t < 100 + tw); return np.trapezoid(i[w] - ss, t[w])
    return charge(1.0) / dV, charge(50.0) / dV, 1000.0 * dV / (ss - base)   # pF, pF, MOhm


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, JOBS))
    else:
        print("f | N topo gII | rest (mV) amp (mV) rate (min-1) | R_in from +5 pA (MOhm) | C 1 ms / C 50 ms (pF) | R_in VC (MOhm)")
        for f in F:
            for n, topo, g in CFG:
                try:
                    r0 = cc_stats(name(f, n, topo, g, "cc")); r5 = cc_stats(name(f, n, topo, g, "inj")); c1, c50, rv = vc_stats(name(f, n, topo, g, "vc"))
                except Exception as e:
                    print(f"{f} | {n} {topo} {g} | failed ({e})"); continue
                rin = 1000.0 * 5.0 / (r0[0] - r5[0]) if r0[0] != r5[0] else float("nan")
                print(f"{f:>4} | N {n} {'chain' if topo == '0' else 'all  '} g {g:>4} | {r0[0]:6.1f} {r0[1]:5.1f} {r0[2]:5.1f} | {rin:7.0f} | "
                      f"{c1:5.1f} / {c50:5.1f} | {rv:6.0f}")
        print("recorded: Koh 1998 cultured 'single' ICC -53 mV, 29 pF, 250 MOhm, pacing; Zhu 2009 single ICC 5.5 pF, -51 mV")
