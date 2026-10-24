"""Fig. 1 calibration re-check in the final v4 configuration (m 0.3, gSS 30, L 30, gK0run 0.4, fP 0.1, gSP 3, iccScale 0.72),
g_c 1 / 1.5 / 2 nS, against Kito et al. 2014 (rabbit ileum, 35 C, nifedipine; Table 1 means +/- SD; text values checked 2026-10-14).
Runs: control = out/v4win/m0.3_gSS30_gc<gc>_u0 (40 s, no input); apamin = out/paper_cal/apa_gc<gc> (skKd = 1000 mM, SK never
activated); fIJP = out/v4ijp (0.5 s purinergic pulse, 200-ms guarded clean windows as in analyze_paper_sens.py).
All metrics over the second half of the run (20-40 s). "CSMC" = SMC layer 1 (V_SMC0); "FLC" = the PDGFa cell of layer 1 (V_P0).
Between-wave samples: ICC-MY below -65 mV and > 300 ms from any sample above it.
Output: paper_R1/fig/cal.png and the printed table.
"""
import os, re, csv
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amp_event as AE                     # event-aligned amplitudes (revision 2); replaces the window values of the logs

plt.rcParams["font.family"] = "Arial"
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAT = re.compile(r"ICC-MY ([\d.]+) cpm amp ([\d.]+) \| L0 ([\d.]+) cpm \| amp by layer: \[0\]([\d.]+)")
GCS = [1, 1.5, 2]
T0 = [20000, 20450, 20900, 21350, 21800, 22250, 22700, 23150]
# (label, Kito mean, SD, note)
KITO = {
    "ICC-MY RMP (mV)": (-65.1, 2.8), "ICC-MY amplitude (mV)": (56.9, 6.4), "ICC-MY rate (min$^{-1}$)": (10.2, 1.1),
    "ICC-MY half-width (sec)": (2.4, 0.3), "ICC-MY dV/dt max (V sec$^{-1}$)": (10.6, 3.4),
    "CSMC RMP (mV)": (-63.5, 2.4), "CSMC amplitude (mV)": (24.4, 2.7), "CSMC half-width (sec)": (2.3, 0.4),
    "CSMC dV/dt max (V sec$^{-1}$)": (0.4, 0.1), "CSMC fIJP (mV)": (9.0, 1.5),
    "FLC RMP (mV)": (-49.5, 5.3), "FLC slow wave (mV)": (14.4, 4.7),
    "FLC apamin depol. (mV)": (9.5, 1.7),
}
RANGE_ONLY = {"CSMC STH max (mV)": (1.1, 8.8), "FLC STH amplitude (mV)": (4.2, 34.2), "FLC STH rate (min$^{-1}$)": (24, 60)}   # Kito ranges
# STH metrics are method dependent in the model: FLC SK hyperpolarizations come in multi-second clusters after each slow wave
# (down to about -75 mV) rather than as discrete ~17.5-mV events at 43/min, so they are reported, not scored.


def load(f, cols=("t_ms", "V_ICC", "V_SMC0", "V_P0")):
    r = list(csv.DictReader(open(f)))
    return [np.array([float(x[k]) for x in r]) for k in cols]


def between(t, vi):
    up = vi >= -65
    near = np.convolve(up.astype(float), np.ones(601), mode="same") > 0      # within 300 ms of an above-threshold sample (1-ms samples)
    return ~near


def waves(t, v, vi):
    """per-cycle half-width (s) and max dV/dt (V/s) of v, cycles defined by ICC-MY upstrokes through -45 mV"""
    k = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]
    hw, dv = [], []
    for a, b in zip(k[:-1], k[1:]):
        seg = v[a - 500:b - 500] if a >= 500 else None
        if seg is None: continue
        lo, hi = seg.min(), seg.max(); half = lo + (hi - lo) / 2
        hw.append((seg > half).sum() / 1000.0)
        dv.append(np.diff(seg).max())                                     # mV/ms = V/s
    return np.mean(hw), np.mean(dv)


def cycle_avg(t, v, vi):
    k = np.where((vi[:-1] < -45) & (vi[1:] >= -45))[0]
    segs = [v[a - 1000:a + 2500] for a in k if a >= 1000 and a + 2500 <= len(v)]
    m = np.mean(segs, axis=0)
    return m.max() - m.min()


def sths(v, mask, thr):
    """hyperpolarizing transients on the detrended trace (v minus a 1-s running median, which removes the slow wave):
    minima deeper than thr mV, >= 250 ms apart, inside mask; returns depths and events per minute of mask time"""
    r = v - median_filter(v, size=1001, mode="nearest")
    p, _ = find_peaks(-r, height=thr, distance=250)
    keep = mask[p]
    return -r[p][keep], keep.sum() / (mask.sum() / 60000.0)


def ijp(gc):
    pre = f"out/v4ijp/gc{gc}"
    tc, ic, sc = load(pre + "_ctrl.csv", ("t_ms", "V_ICC", "V_SMC0"))
    out = []
    for t0 in T0:
        t, i, s = load(f"{pre}_p{t0}.csv", ("t_ms", "V_ICC", "V_SMC0"))
        n = min(len(t), len(tc)); w = (t[:n] >= t0) & (t[:n] < t0 + 1000); g = (t[:n] >= t0) & (t[:n] < t0 + 1200)
        if ic[:n][g].max() < -65 and i[:n][g].max() < -65:
            out.append(-(s[:n] - sc[:n])[w].min())
    return np.mean(out), len(out)


M = {}
for gc in GCS:
    f = f"out/v4win/m0.3_gSS30_gc{gc}_u0"
    cpm, ampI, _, ampS = (float(x) for x in PAT.search(AE.log(f + ".log")).groups())
    t, vi, vs, vp = load(f + ".csv")
    h = t >= 20000; t, vi, vs, vp = t[h], vi[h], vs[h], vp[h]
    bw = between(t, vi)
    hwI, dvI = waves(t, vi, vi); hwS, dvS = waves(t, vs, vi)
    sa, _ = sths(vs, bw, 1.0)
    fa, fr = sths(vp, np.ones_like(vp, bool), 4.0)
    ta, va_i, va_s, va_p = load(f"out/paper_cal/apa_gc{gc}.csv")
    ha = ta >= 20000
    fa0, fr0 = sths(va_p[ha], np.ones(ha.sum(), bool), 4.0)
    ij, nij = ijp(gc)
    M[gc] = {
        "ICC-MY RMP (mV)": np.percentile(vi, 5), "ICC-MY amplitude (mV)": ampI, "ICC-MY rate (min$^{-1}$)": cpm,
        "ICC-MY half-width (sec)": hwI, "ICC-MY dV/dt max (V sec$^{-1}$)": dvI,
        "CSMC RMP (mV)": np.median(vs[bw]), "CSMC amplitude (mV)": ampS, "CSMC half-width (sec)": hwS, "CSMC dV/dt max (V sec$^{-1}$)": dvS,
        "CSMC fIJP (mV)": ij,
        "FLC RMP (mV)": np.median(vp[bw]), "FLC slow wave (mV)": cycle_avg(t, vp, vi), "FLC STH rate (min$^{-1}$)": fr,
        "FLC apamin depol. (mV)": va_p[ha].mean() - vp.mean(),
        "CSMC STH max (mV)": sa.max() if len(sa) else 0.0, "FLC STH amplitude (mV)": np.median(fa) if len(fa) else 0.0,
        "_csmc_sth_n": len(sa), "_csmc_sth_mean": sa.mean() if len(sa) else 0.0, "_flc_sth_apamin_rate": fr0, "_ijp_n": nij,
        "_flc_mean": vp.mean(), "_flc_apamin_mean": va_p[ha].mean(),
    }

print("metric | Kito | g_c 1 | 1.5 | 2   (z = (model - mean)/SD)")
for k, (mu, sd) in KITO.items():
    print(f"{k:26s} | {mu:6.1f} ± {sd:4.1f} | " + " | ".join(f"{M[g][k]:7.2f} (z {(M[g][k] - mu) / sd:+5.1f})" for g in GCS))
for k, (lo, hi) in RANGE_ONLY.items():
    print(f"{k:26s} | range {lo}-{hi} | " + " | ".join(f"{M[g][k]:7.2f}" for g in GCS))
for g in GCS:
    print(f"g_c {g}: CSMC STH n {M[g]['_csmc_sth_n']} mean {M[g]['_csmc_sth_mean']:.1f} mV; FLC STH rate with apamin "
          f"{M[g]['_flc_sth_apamin_rate']:.1f}/min; FLC mean {M[g]['_flc_mean']:.1f} -> apamin {M[g]['_flc_apamin_mean']:.1f}; fIJP clean n {M[g]['_ijp_n']}")

# ---------- figure ----------
fig = plt.figure(figsize=(14, 7.5))
gsp = fig.add_gridspec(3, 2, width_ratios=[1.35, 1])
t, vi, vs, vp = load("out/v4win/m0.3_gSS30_gc1.5_u0.csv")
ta, _, _, vpa = load("out/paper_cal/apa_gc1.5.csv")
w = (t >= 30000) & (t < 38000); x = (t[w] - 30000) / 1000
for r, (v, lab, c) in enumerate(((vi, "ICC-MY", "tab:orange"), (vs, "CSMC (SMC layer 1)", "tab:blue"), (vp, "FLC (layer 1)", "tab:green"))):
    a = fig.add_subplot(gsp[r, 0])
    a.plot(x, v[w], color=c, lw=1.1, label=lab)
    if r == 2:
        a.plot(x, vpa[(ta >= 30000) & (ta < 38000)], color="0.5", lw=0.9, label="FLC, apamin (SK blocked)")
        a.set_xlabel("time (sec), g$_c$ 1.5 nS")
    a.set_ylabel("mV"); a.legend(fontsize=7, loc="upper right")
    if r == 0: a.set_title("A  Model traces in the final configuration", loc="left", fontsize=10)
a = fig.add_subplot(gsp[:, 1])
keys = list(KITO.keys()); y = np.arange(len(keys))[::-1]
a.axvspan(-1, 1, color="0.9"); a.axvline(0, color="0.5", lw=0.8)
for g, mk, c in zip(GCS, ("o", "s", "^"), ("tab:blue", "k", "tab:red")):
    z = [(M[g][k] - KITO[k][0]) / KITO[k][1] for k in keys]
    a.plot(np.clip(z, -6, 6), y, mk, color=c, ms=5, label=f"g$_c$ {g} nS", alpha=0.8)
a.set_yticks(y); a.set_yticklabels(keys, fontsize=8); a.set_xlim(-6.5, 6.5)
a.set_xlabel("(model − Kito mean) / Kito SD   (clipped at ±6)")
a.set_title("B  Model vs Kito et al. 2014 (shaded ±1 SD)", loc="left", fontsize=10); a.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig("paper_R1/fig/cal.png", dpi=150)
print("wrote paper_R1/fig/cal.png")
