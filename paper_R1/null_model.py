"""Passive null model for the amplitude-ratio estimate of g_c (reviewer point 1-1), and what the full model adds.

The muscle column is a passive sealed cable of L layers, m cells per layer, each cell = SMC + f_P FLC share:
membrane conductance g_m 0.96 nS, capacitance c_m = 50 pF (SMC) + 0.1 x 25 pF (FLC) = 52.5 pF; layers coupled by g_SS per cell pair.
The ICC-MY voltage is imposed (the model ICC-MY waveform at the calibrated point), so layer 1 follows
V1(w) = g_c Vicc(w) / (g_c + Y(w)), Y(w) = admittance of the discrete column at layer 1.
  DC null : Y = Y(0)             -> g_c = r/(1-r) * Y(0)
  AC null : Y(w) with capacitance, event-aligned amplitude ratio (peak - pre-upstroke baseline) of the filtered model ICC-MY waveform
Targets: rabbit ileum r = 24.4/56.9 (Kito et al. 2014) and mouse small intestine r = 34.1/56.5 (Kito & Suzuki 2003).
Full model: g_c at the target ratio from criteria.out.txt. Decomposition at the calibrated column: SK blocked (apamin = skKd 1000;
out/paper_cal/apa_gc1-2, out/paper_rev1/apa_gc2.5-3) isolates the FLC contribution; the rest of the gap to the null model comes
from the SMC membrane.
Output: null_model.out.txt (run as a script); results(target) is used by fig2_assemble.py.
"""
import os, csv, re
import numpy as np
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GM, CM = 0.96, 52.5                        # nS, pF per cell (SMC + FLC share)
R = 24.4 / 56.9                            # rabbit ileum, Kito et al. 2014
R_MOUSE = 34.1 / 56.5                      # mouse small intestine, Kito & Suzuki 2003
TARGETS = (("rabbit", R, 5), ("mouse", R_MOUSE, 9))   # criteria.out.txt column with g_c at each target
GEO = {"calibrated column": (0.3, 30, 30), "g_SS 15": (0.3, 15, 30), "g_SS 60": (0.3, 60, 30), "m 0.15": (0.15, 30, 30),
       "m 0.45": (0.45, 30, 30), "L 15": (0.3, 30, 15), "L 45": (0.3, 30, 45), "L 7, m 1.3": (1.3, 30, 7),
       "L 8, m 3, g_SS 30": (3, 30, 8), "L 8, m 3, g_SS 300": (3, 300, 8)}
PAT = re.compile(r"amp ([\d.]+) \| L0 [\d.]+ cpm \| amp by layer: \[0\]([\d.]+)")
APA = [(1, "out/paper_cal/apa_gc1"), (1.5, "out/paper_cal/apa_gc1.5"), (2, "out/paper_cal/apa_gc2"),
       (2.5, "out/paper_rev1/apa_gc2.5"), (3, "out/paper_rev1/apa_gc3")]


def column_admittance(m, gss, L, w):
    """input admittance (nS) at layer 1 of a sealed passive cable at angular frequency w (rad/ms); pF x rad/ms = nS"""
    y = m * (GM + 1j * w * CM)
    Y = y
    for _ in range(L - 1):
        Y = y + (m * gss * Y) / (m * gss + Y)
    return Y


def ratio_ac(gc, m, gss, L, vicc):
    V = np.fft.rfft(vicc - vicc.mean()); w = 2 * np.pi * np.fft.rfftfreq(len(vicc), 1.0)
    H = np.array([gc / (gc + column_admittance(m, gss, L, x)) for x in w])
    v1 = np.fft.irfft(V * H, len(vicc))
    # event-aligned, as amp_event.py: peak within 2 sec after the upstroke minus mean over 50-200 msec before it (the 5th upstroke)
    u = 4 * len(vicc) // 8
    ev = lambda v: v[u:u + 2000].max() - v[u - 200:u - 50].mean()
    return ev(v1) / ev(vicc)


def solve(f, target, lo=0.01, hi=200):
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) < target else (lo, mid)
    return np.sqrt(lo * hi)


def num(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


def results(target="rabbit"):
    """per-geometry passive and full-model g_c for one target, and the SK-blocked estimate at the calibrated column"""
    _, Rt, col = next(t for t in TARGETS if t[0] == target)
    r = list(csv.DictReader(open(os.path.join(ROOT, "out/v4win/m0.3_gSS30_gc1_u0.csv"))))
    t = np.array([float(x["t_ms"]) for x in r]); vi = np.array([float(x["V_ICC"]) for x in r])
    up = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]; up = up[t[up] >= 30000]
    vicc = np.tile(vi[up[0]:up[1]], 8)
    rows = []
    for line in open(os.path.join(ROOT, "paper_R1/criteria.out.txt"), encoding="utf-8"):
        f = [x.strip() for x in line.split("|")]
        if len(f) > col and f[1] == "geometry":
            m, gss, L = GEO[f[0]]; Yd = column_admittance(m, gss, L, 0).real
            ac = solve(lambda g: ratio_ac(g, m, gss, L, vicc), Rt)
            rows.append(dict(name=f[0], G=float(f[2]), Ydc=Yd, dc=Rt / (1 - Rt) * Yd, ac=ac, full=num(f[col])))
    runs = [(g, f) for g, f in APA if os.path.exists(os.path.join(ROOT, f + ".log")) and os.path.getsize(os.path.join(ROOT, f + ".log"))]
    gcs = [g for g, _ in runs]
    rat = []
    for _, f in runs:
        a, s = (float(x) for x in PAT.search(AE.log(os.path.join(ROOT, f + ".log"))).groups())
        rat.append(s / a)
    nosk = float(np.exp(np.interp(Rt, rat, np.log(gcs)))) if min(rat) <= Rt <= max(rat) else float("nan")
    return rows, nosk, dict(zip(gcs, rat))


if __name__ == "__main__":
    out = []
    for target, Rt, _ in TARGETS:
        rows, nosk, rat = results(target)
        out.append(f"== target {target}, ratio {Rt:.3f} (passive g_c/Y = {Rt / (1 - Rt):.3f})")
        out.append("geometry | G_load formula (nS) | Y_dc discrete (nS) | g_c DC null | g_c AC null | g_c full model | full/DC")
        for d in rows:
            out.append(f"{d['name']:20s} | {d['G']:5.2f} | {d['Ydc']:5.2f} | {d['dc']:5.2f} | {d['ac']:5.2f} | {d['full']:5.2f} | {d['full'] / d['dc']:.2f}")
        cal = rows[0]
        out.append(f"calibrated column, SK blocked: layer-1/ICC-MY ratio {', '.join(f'{k}: {v:.3f}' for k, v in rat.items())} -> g_c at {Rt:.3f} = {nosk:.2f} nS")
        gap = cal["dc"] - cal["full"]
        out.append(f"decomposition: passive {cal['dc']:.2f} -> SMC membrane active, no SK {nosk:.2f} -> full {cal['full']:.2f} nS; "
                   f"FLC/SK share {(nosk - cal['full']) / gap:.0%}, SMC membrane share {(cal['dc'] - nosk) / gap:.0%}")
        vals = [d["full"] / d["dc"] for d in rows if d["full"] == d["full"]]
        out.append(f"full/passive across geometries: {min(vals):.2f}-{max(vals):.2f} (n {len(vals)})")
    open(os.path.join(ROOT, "paper_R1/null_model.out.txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("\n".join(out))
