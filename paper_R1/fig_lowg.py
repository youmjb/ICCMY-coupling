"""Preview figure: a low own conductance of ICC-MY (iccRest f) reconciles the drug data (lowg_coupled.py runs, seed 1, g_c 1.75 nS,
control ratio 0.52; pinacidil sized by interpolation to the recorded muscle change of -10.3 mV) and the in situ input resistance
(passive ICC-MY lattice, 41 x 41, each node loaded by its column). Output: fig/lowg_preview.png, fig_lowg.out.txt."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
import figstyle; figstyle.apply()
from lowg_coupled import S, F, GC
from icc_hyperpol import stats
OUT = []
PINS = ["0.5", "1", "1.5", "2", "4", "8"]
GOWN = {"1": 2.91, "0.5": 1.51, "0.3": 0.90, "0.2": 0.72, "0.1": None, "0.05": None}   # isolated cell, 5 pA chord (icc_lowg.py)


def sized(f, gc="1.75", target=-10.3):
    c0 = S(f, "ctrl", gc); rows = []
    for p in PINS:
        c1 = S(f, f"pin{p}", gc)
        if c1 and c0:
            dM, dI = c1["sr"] - c0["sr"], c1["ir"] - c0["ir"]
            rows.append((dM, dI / dM, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
    rows.sort()
    x = np.array([r[0] for r in rows])
    return [float(np.interp(target, x, [r[k] for r in rows])) for k in (1, 2, 3)], (x.min(), x.max())


def gin(gii, gown, gc=1.74, gl=1.61, n=41):
    c = n // 2; A = lil_matrix((n * n, n * n)); g = gown + gc * gl / (gc + gl)
    for i in range(n):
        for j in range(n):
            k = i * n + j; A[k, k] += g
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n: A[k, k] += gii; A[k, ii * n + jj] -= gii
    b = np.zeros(n * n); b[c * n + c] = 1.0
    return 1.0 / spsolve(A.tocsr(), b)[c * n + c]


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 3, figsize=(figstyle.WIDTH, 2.5))
    # A: transfer against ICC-MY amplitude change at the pinacidil level
    a = ax[0]
    for f, mk in zip(F, ("o", "s", "^", "D", "v")):
        (T, xi, xm), rng = sized(f)
        a.plot(xi, T, mk, color="k", mfc="white" if f == "1" else "k", ms=4, label=f"f {f}")
        OUT.append(f"A f {f}: transfer {T:.2f}, ICC-MY x{xi:.2f}, muscle x{xm:.2f} (pinacidil dM range {rng[0]:.1f} to {rng[1]:.1f})")
    # the reference model's trade-off (f 1, g_c 0.5-3, pinacidil 2 nS) for comparison
    tr = [(S("1", "pin2", g)["ir"] - S("1", "ctrl", g)["ir"]) / (S("1", "pin2", g)["sr"] - S("1", "ctrl", g)["sr"]) for g in GC]
    ia = [S("1", "pin2", g)["a"]["icc"] / S("1", "ctrl", g)["a"]["icc"] for g in GC]
    a.plot(ia, tr, "-", color="0.6", lw=0.8, label="f 1, g$_c$ 0.5–3")
    a.plot(1.17, 0.93, "*", color="k", ms=9, label="recorded")
    a.set_xlabel("ICC-MY amplitude (fold)"); a.set_ylabel("resting transfer"); a.set_ylim(0, 1.05)
    a.legend(frameon=False, fontsize=6, loc="lower right"); a.set_title("A", loc="left", fontweight="bold")
    # B: control amplitudes and rest in the column at g_c 1.75 against f
    a = ax[1]; ff = [float(f) for f in F]
    a.plot(ff, [S(f, "ctrl", "1.75")["a"]["icc"] for f in F], "o-", color="k", mfc="white", ms=3.5, label="ICC-MY amplitude")
    a.plot(ff, [S(f, "ctrl", "1.75")["a"]["smc"] for f in F], "o-", color="k", ms=3.5, label="muscle amplitude")
    a.plot(ff, [-S(f, "ctrl", "1.75")["ir"] for f in F], "s--", color="0.5", ms=3, label="−ICC-MY rest")
    a.axhspan(56.1 - 5.3, 56.1 + 5.3, color="0.9", lw=0); a.axhspan(29.3 - 4.6, 29.3 + 4.6, color="0.9", lw=0)
    a.set_xscale("log"); a.set_xlabel("resting-current scale f"); a.set_ylabel("mV"); a.legend(frameon=False, fontsize=6)
    a.set_title("B", loc="left", fontweight="bold")
    for f in F:
        c = S(f, "ctrl", "1.75"); r = stats(f"out/paper_lowg/col_f{f}_ctrl_gc1.75.csv")
        OUT.append(f"B f {f}: ICC rest {c['ir']:.1f}, amp {c['a']['icc']:.1f}, muscle rest {c['sr']:.1f}, amp {c['a']['smc']:.1f}, rate {r[4]:.1f}")
    # C: in situ input resistance of ICC-MY in a lattice
    a = ax[2]; gii = np.logspace(-1, 1.3, 30)
    for gown, ls in ((2.91, "-"), (0.72, "--"), (0.13, ":")):
        r = [1000.0 / gin(g, gown) for g in gii]
        a.plot(gii, r, ls, color="k", lw=1, label=f"g$_{{own}}$ {gown} nS")
        OUT.append(f"C g_own {gown}: R_in alone {1000 / gin(0, gown):.0f} MOhm; g_II for 260 MOhm {np.interp(-260, [-x for x in r], gii):.2f} nS")
    a.axhspan(222, 298, color="0.9", lw=0)                              # Cai et al. mean ± SEM (Fig. S7); a.axhline(260, color="0.5", lw=0.6)
    a.set_xscale("log"); a.set_xlabel("ICC-MY–ICC-MY coupling (nS)"); a.set_ylabel("R$_{in}$ of ICC-MY (MΩ)")
    a.legend(frameon=False, fontsize=6); a.set_title("C", loc="left", fontweight="bold")
    fig.tight_layout()
    figstyle.save(fig, "paper_R1/fig/lowg_preview")
    open("paper_R1/fig_lowg.out.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n"); print("\n".join(OUT))
