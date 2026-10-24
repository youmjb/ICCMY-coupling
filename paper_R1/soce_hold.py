"""Published ICC model (2019, 25 pF, unmodified): does it still generate spontaneous inward currents at a holding potential of -80 mV when
store-operated Ca2+ entry (SOC/Orai, iccSOC) is removed while ANO1 is kept? Also free running (I = 0) and the ANO1-removed control.
SOC and ANO1 scales apply from the start of the 20-s uncoupled equilibration. Voltage clamp: 'cluster' mode, N 1, Rs 2 MOhm, 60 s.
Usage: python paper_R1/soce_hold.py run | analyze      Runs: out/paper_soce/"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = "out/paper_soce"
CASES = {"ctrl": [], "soc0.5": ["iccSOC=0.5"], "soc0.25": ["iccSOC=0.25"], "soc0": ["iccSOC=0"], "ano0": ["iccANO1=0"]}
T = "60000"


def name(c, m): return f"{OUT}/{c}_{m}.csv"


def run(job):
    c, m = job; out = name(c, m)
    if os.path.exists(out + ".log"): return
    extra = ["vc=1", "Rs=2", "Vh=-80", "dV=0", "smp=0.5"] if m == "vc" else ["smp=1"]
    with open(out + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "cluster", T, "0.005" if m == "vc" else "0.02", out, "N=1", "s=1"] + CASES[c] + extra, stdout=lg)


def events(t, y, inward=True, thr=50.0):
    """events: excursions beyond the quiet level (10th/90th percentile) by more than thr; returns (count, mean peak, mean half-width)"""
    q = np.percentile(y, 90 if inward else 10)
    s = (q - y) if inward else (y - q)
    on = np.where((s[:-1] < thr) & (s[1:] >= thr))[0]
    pk, hw = [], []
    for a in on:
        e = a + np.argmax(s[a:] < thr) if np.any(s[a:] < thr) else len(s) - 1
        seg = s[max(0, a - 2000):e + 1]; tt = t[max(0, a - 2000):e + 1]
        p = seg.max(); idx = np.where(seg >= p / 2)[0]
        pk.append(p); hw.append(tt[idx[-1]] - tt[idx[0]])
    return q, len(on), (np.mean(pk) if pk else 0.0), (np.mean(hw) if hw else 0.0)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        with ThreadPoolExecutor(10) as ex: list(ex.map(run, [(c, m) for c in CASES for m in ("vc", "cc")]))
    else:
        print("case | voltage clamp at -80 mV (last 50 s): quiet holding current, events > 50 pA inward, mean peak, half-width | "
              "free running (last 50 s): rest, peak, events > 20 mV")
        for c in CASES:
            r = np.loadtxt(name(c, "vc"), delimiter=",", skiprows=1); t, i = r[:, 0], r[:, 3]; h = t >= 10000
            q, n, pk, hw = events(t[h], i[h], inward=True, thr=50.0)
            r2 = np.loadtxt(name(c, "cc"), delimiter=",", skiprows=1); t2, v = r2[:, 0], r2[:, 1]; h2 = t2 >= 10000
            qv, nv, pv, hv = events(t2[h2], v[h2], inward=False, thr=20.0)
            print(f"  {c:8s} | hold {q:7.1f} pA, {n:2d} events, peak {pk:6.0f} pA, half-width {hw:5.0f} ms | "
                  f"rest {np.percentile(v[h2], 5):6.1f} mV, max {v[h2].max():6.1f} mV, {nv} events ({60 * nv / 50:.1f} min-1)")
