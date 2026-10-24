"""Input resistance and apparent capacitance of ICC clusters measured between slow waves, as with a current pulse in current clamp.
For each cluster_sim.py configuration the 30-s current-clamp run (cc) gives a quiet diastolic interval after 15 s; the same deterministic
run is repeated with a 2-pA outward (hyperpolarizing) pulse into cell 0 for 400 ms from 300 ms after a repolarization. R_in = dV at
400 ms / 2 pA; tau from the time to 63 % of dV; apparent C = tau / R_in. Cells that do not pace are pulsed at 20 s.
Usage: python paper_R1/cluster_pulse.py run | analyze      Runs: out/paper_cluster/pulse_*.csv"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import cluster_sim as CS
I0, DUR = 2.0, 400.0


def onset(f, n, topo, g):
    t, v, _ = CS.load(CS.name(f, n, topo, g, "cc"))
    h = t >= 15000; t, v = t[h], v[h]
    rest, peak = np.percentile(v, 5), np.percentile(v, 99.5)
    if peak - rest < 10: return 20000.0
    mid = rest + 0.5 * (peak - rest)
    down = np.where((v[:-1] >= mid) & (v[1:] < mid))[0]; up = np.where((v[:-1] < mid) & (v[1:] >= mid))[0]
    for d in down:
        t0 = t[d] + 600.0
        nxt = up[up > d]
        if len(nxt) and t[nxt[0]] > t0 + DUR + 200: return t0
    return None


def pname(f, n, topo, g): return f"{CS.OUT}/pulse_f{f}_N{n}_t{topo}_g{g}.csv"


def run(cfg):
    f, n, topo, g = cfg; out = pname(*cfg)
    if os.path.exists(out + ".log"): return
    t0 = onset(f, n, topo, g)
    if t0 is None: return
    with open(out + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "cluster", f"{t0 + DUR + 100:.0f}", "0.02", out] + CS.VAR
                       + [f"iccRest={f}", f"N={n}", f"topo={topo}", f"gII={g}", f"inj={I0}", f"injOn={t0:.1f}", f"injOff={t0 + DUR:.1f}"], stdout=lg)
    open(out + ".t0", "w").write(f"{t0}\n")


def measure(cfg):
    f, n, topo, g = cfg
    t0 = float(open(pname(*cfg) + ".t0").read())
    tc, vc, _ = CS.load(CS.name(f, n, topo, g, "cc")); tp, vp, _ = CS.load(pname(*cfg))
    w = (tp >= t0 - 5) & (tp <= t0 + DUR)
    d = vp[w] - np.interp(tp[w], tc, vc); tt = tp[w] - t0
    dV = d[tt <= DUR][-1]
    k = np.where(d <= 0.632 * dV)[0] if dV < 0 else np.where(d >= 0.632 * dV)[0]
    tau = tt[k[0]] if len(k) else float("nan")
    R = abs(dV) / I0 * 1000.0
    return R, tau, tau / R * 1000.0, dV


CFGS = [(f, n, topo, g) for f in CS.F for n, topo, g in CS.CFG]
if __name__ == "__main__":
    if sys.argv[1] == "run":
        with ThreadPoolExecutor(30) as ex: list(ex.map(run, CFGS))
    else:
        print("f | N topo gII | true C (pF) | R_in (MOhm) | tau (ms) | apparent C (pF) | rest (mV) paced")
        for cfg in CFGS:
            f, n, topo, g = cfg
            try:
                R, tau, C, dV = measure(cfg); r0 = CS.cc_stats(CS.name(f, n, topo, g, "cc"))
            except Exception as e:
                print(f"{f} | {n} {topo} {g} | failed ({e})"); continue
            print(f"{f:>4} | N {n} {'chain' if topo == '0' else 'all  '} g {g:>4} | {5.5 * int(n):5.1f} | {R:6.0f} | {tau:6.1f} | {C:6.1f} | "
                  f"{r0[0]:6.1f} {'yes' if r0[1] > 0 else 'no'}")
        print("recorded: Koh 1998 'single' cultured ICC 29 pF, 250 MOhm (tau ~7 ms), -53 mV, pacing")
