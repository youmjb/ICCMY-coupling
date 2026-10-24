"""ICC loss and transmural decay (revision 3; runs from run_paper_rev3.py, event-aligned amplitudes as amp_event.py).

loss(): for each g_c (rabbit 1.17, mouse 2.33 nS, kept per remaining ICC-MY) and remaining ICC-MY fraction f (m = 0.3/f):
        G_load, ICC-MY and layer-1 amplitudes, ICC-MY and layer-1 rates, and the passive-ladder prediction of the layer-1/ICC-MY
        ratio at the same g_c (DC divider, g_m 0.96 nS).
transmural(): event-aligned amplitude of every muscle layer (peak within 2 sec after the ICC-MY upstroke minus the mean over
        50-200 msec before it), for the calibrated column (L 30, m 0.3) and the mouse-like column (L 7, m 1.3) at their rabbit and
        mouse g_c, with the passive-ladder decay (DC, g_m 0.96 nS) for comparison.
Output (run as a script): paper_R1/icc_loss.out.txt
"""
import os, re, sys, csv
import numpy as np

if hasattr(sys.stdout, "reconfigure"): sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
import amp_event as AE

OUT = os.path.join(ROOT, "out", "paper_rev3")
FRAC = (1.0, 0.8, 2 / 3, 0.5, 0.4, 0.3, 0.2, 0.1)
GC = {"rabbit": 1.17, "mouse": 2.33}
TM = {"L30_rabbit": (1.17, 0.3, 30), "L30_mouse": (2.33, 0.3, 30), "L7_rabbit": (4.22, 1.3, 7), "L7_mouse": (8.26, 1.3, 7)}
GM, GSS = 0.96, 30.0
ICC_UNCOUPLED = 56.3          # event-aligned amplitude of the uncoupled model ICC-MY (out/v4win/m0.3_gSS30_gc0_u0)


def gload(m, L=30):
    lam = np.sqrt(GSS / GM); return m * np.sqrt(GSS * GM) * np.tanh(L / lam)


def ladder(m, L):
    """input admittance from each layer to the sealed end (per column, DC): Y[k] for layers k = 0..L-1"""
    Y = [0.0] * L; Y[L - 1] = m * GM
    for k in range(L - 2, -1, -1):
        Y[k] = m * GM + (m * GSS * Y[k + 1]) / (m * GSS + Y[k + 1])
    return Y


def loss():
    rows = []
    for sp, gc in GC.items():
        for f in FRAC:
            pre = f"{OUT}/loss_gc{gc}_f{f:.3f}"
            if not os.path.exists(pre + ".csv"): continue
            ev = AE.measure(pre + ".csv"); log = open(pre + ".log").read()
            ci, cl = (float(x) for x in re.search(r"ICC-MY ([\d.]+) cpm .*\| L0 ([\d.]+) cpm", log).groups())
            m = 0.3 / f; Y = ladder(m, 30)[0]
            rows.append(dict(sp=sp, gc=gc, f=f, m=m, G=gload(m), icc=ev["icc"], smc=ev["smc"], cpm_icc=ci, cpm_l0=cl,
                             passive=gc / (gc + Y)))
    return rows


def transmural():
    res = {}
    for name, (gc, m, L) in TM.items():
        p = f"{OUT}/tm_{name}.csv.layers.csv"
        if not os.path.exists(p): continue
        with open(p) as fh:
            r = csv.reader(fh); h = next(r); a = np.array([[float(x) for x in row] for row in r])
        a = a[len(a) // 2:]; t, vi = a[:, 0], a[:, 1]; dt = t[1] - t[0]
        mid = vi.min() + 0.5 * (vi.max() - vi.min())
        up = np.where((vi[:-1] < mid) & (vi[1:] >= mid))[0] + 1
        up = up[(t[up] - t[0] > 300) & (t[-1] - t[up] > 2100)]
        n200, n50, n2s = int(200 / dt), int(50 / dt), int(2000 / dt)
        amp = [float(np.mean([a[u:u + n2s, 2 + k].max() - a[u - n200:u - n50, 2 + k].mean() for u in up])) for k in range(L)]
        Y = ladder(m, L); dec = [1.0]
        for k in range(1, L):                       # DC voltage transfer layer k-1 -> k in the passive ladder
            dec.append(dec[-1] * (m * GSS) / (m * GSS + Y[k]))
        res[name] = dict(gc=gc, m=m, L=L, amp=amp, passive=[amp[0] * d for d in dec])
    return res


if __name__ == "__main__":
    out = ["(1) ICC loss: remaining ICC-MY fraction f, m = 0.3/f, g_c per ICC-MY fixed",
           "species g_c | f | m | G_load | ICC-MY amp (% uncoupled) | layer-1 amp | layer-1/ICC-MY (passive DC) | ICC-MY / layer-1 rate (cpm)"]
    rows = loss()
    for d in rows:
        out.append(f"{d['sp']} {d['gc']} | {d['f']:.2f} | {d['m']:.3f} | {d['G']:.2f} | {d['icc']:.1f} ({d['icc'] / ICC_UNCOUPLED:.0%}) | "
                   f"{d['smc']:.1f} | {d['smc'] / d['icc']:.3f} ({d['passive']:.3f}) | {d['cpm_icc']:.1f} / {d['cpm_l0']:.1f}")
    for sp in GC:
        r = [d for d in rows if d["sp"] == sp]
        if not r: continue
        base = r[0]["smc"]
        for thr, lab in ((0.75 * base, "75 % of intact"), (0.5 * base, "half of intact")):
            fs = [d["f"] for d in r]; ys = [d["smc"] for d in r]
            cross = next((f0 + (thr - y0) / (y1 - y0) * (f1 - f0) for (f0, y0), (f1, y1) in zip(zip(fs, ys), zip(fs[1:], ys[1:]))
                          if (y0 - thr) * (y1 - thr) <= 0 and y1 != y0), float("nan"))
            out.append(f"  {sp}: layer-1 amplitude reaches {lab} at f = {cross:.2f} ({1 - cross:.0%} of ICC-MY lost)")
    out.append("\n(2) transmural decay: event-aligned amplitude by layer (mV); passive DC ladder from layer 1 in parentheses")
    for name, d in transmural().items():
        a, p = d["amp"], d["passive"]
        idx = [k for k in (0, 1, 2, 4, 6, 9, 14, 19, 24, 29) if k < d["L"]]
        out.append(f"{name} (g_c {d['gc']}, m {d['m']}, L {d['L']}): " + ", ".join(f"L{k + 1} {a[k]:.1f} ({p[k]:.1f})" for k in idx))
        half = next((k + 1 for k in range(d["L"]) if a[k] < 0.5 * a[0]), None)
        out.append(f"   halves by layer {half}; deepest layer {a[-1]:.1f} mV = {a[-1] / a[0]:.0%} of layer 1; mean over layers {np.mean(a):.1f} mV")
    open(os.path.join(HERE, "icc_loss.out.txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("\n".join(out))
