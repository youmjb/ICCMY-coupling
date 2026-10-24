"""Supplementary Fig. S2: SmoothM single-cell behaviour (paper_R1/smoothm_py.py, verified against the C++ code).
A steady-state gating and time constants; B voltage clamp from -70 mV (500-ms steps): peak I_CaL and end-of-step K+ current,
as densities (Cm 50 pF); C current clamp with tonic K+ 0.4 nS (C++ smsstep output, out/paper_smc/step_gK0.4.csv);
D quasi-steady-state I-V with and without the tonic K+ conductance (gates and [Ca]i at steady state for each V).
Printed numbers go to the supplement text. Output: paper_R1/fig/supp_smc.png
"""
import os, sys, csv, math, copy
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); os.chdir(ROOT)
import smoothm_py as sm

plt.rcParams["font.family"] = "Arial"
fig, ax = plt.subplots(2, 2, figsize=(11, 8)); ax = ax.ravel()

# A gating
vv = np.linspace(-90, 40, 261)
g = np.array([sm.gates_inf_tau(v) for v in vv])
a = ax[0]
a.plot(vv, g[:, 0], color="tab:red", label="I$_{CaL}$ activation d$_\\infty$")
a.plot(vv, g[:, 2], color="tab:red", ls="--", label="I$_{CaL}$ inactivation f$_\\infty$")
a.plot(vv, g[:, 5] ** 2, color="tab:blue", label="I$_{K}$ n$_\\infty^2$")
for cai, ls in ((0.1, ":"), (1.0, "-.")):
    a.plot(vv, [sm.kca_h(v, cai) for v in vv], color="tab:green", ls=ls, label=f"I$_{{KCa}}$ h$_\\infty$, Ca$_i$ {cai} µM")
a.set_xlabel("V (mV)"); a.set_ylabel("steady state"); a.legend(fontsize=7, loc="center left")
a2 = a.inset_axes([0.62, 0.55, 0.35, 0.4])
a2.semilogy(vv, g[:, 1], color="tab:red"); a2.semilogy(vv, g[:, 3], color="tab:red", ls="--"); a2.semilogy(vv, g[:, 4], color="tab:red", ls=":")
a2.semilogy(vv, g[:, 6], color="tab:blue"); a2.set_ylim(1e-2, 1e3)
a2.set_title("τ (ms); f1 τ → 0 below −40 mV", fontsize=7); a2.tick_params(labelsize=6)
a.set_title("A  Gating", loc="left", fontsize=10)

# B voltage clamp
dt = 0.01
hold = sm.SmoothM()
for _ in range(int(3000 / dt)): hold.step(dt, vclamp=-70.0)
steps = np.arange(-60, 45, 5); pk, kend = [], []
for vs in steps:
    c = copy.deepcopy(hold); p = 0.0
    for i in range(int(500 / dt)):
        ical, iks, ikca, _ = c.step(dt, vclamp=float(vs))
        p = min(p, ical)
    pk.append(p / sm.CM); kend.append((iks + ikca) / sm.CM)
pk, kend = np.array(pk), np.array(kend)
a = ax[1]
a.plot(steps, pk, "o-", color="tab:red", ms=3, label="peak I$_{CaL}$")
a.plot(steps, kend, "s-", color="tab:blue", ms=3, label="I$_K$ + I$_{KCa}$ at 500 ms")
a.axhline(0, color="k", lw=0.6); a.set_xlabel("test potential (mV), holding −70 mV"); a.set_ylabel("current density (pA/pF)")
a.set_title("B  Voltage clamp", loc="left", fontsize=10); a.legend(fontsize=8)
i0 = int(np.argmin(pk))
print(f"B: peak I_CaL {pk[i0]:.2f} pA/pF at {steps[i0]} mV; I_K at +40 mV {kend[-1]:.2f} pA/pF")

# C current clamp (C++)
r = list(csv.DictReader(open("out/paper_smc/step_gK0.4.csv")))
t = np.array([float(x["t_ms"]) for x in r])
a = ax[2]
for k, col in (("A5", "0.7"), ("A10", "0.5"), ("A20", "tab:orange"), ("A40", "tab:red")):
    a.plot(t - 100, [float(x[k]) for x in r], color=col, label=f"+{k[1:]} pA")
a.axvspan(0, 900, color="0.93", zorder=0)
a.set_xlabel("time from current onset (ms)"); a.set_ylabel("V (mV)"); a.legend(fontsize=8)
a.set_title("C  Current clamp, tonic K$^+$ 0.4 nS (rest −67 mV, R$_{in}$ 1.54 GΩ)", loc="left", fontsize=10)

# D quasi-steady-state I-V
def iss(v, gk):
    d_inf, _, f_inf, _, _, n_inf, _ = sm.gates_inf_tau(v)
    c = sm.SmoothM(); c.d, c.f1, c.f2, c.n = d_inf, f_inf, f_inf, n_inf
    ical = sm.GCAL * d_inf * f_inf * f_inf * (v - sm.ECAK)                   # Ca-dependent factor ~1 at 0.1 uM; [Ca]i set below
    k = 1e-1 * sm.COMPLOSS / sm.SACC
    c.cai = 0.1 + (-1e-9 * ical / 2 / sm.SUBV / sm.F) / k
    ical, iks, ikca, ileak = c.currents(v)
    return ical + iks + ikca + ileak + gk * (v + 80.0)
vq = np.linspace(-90, 0, 181)
a = ax[3]
for gk, col in ((0.0, "0.5"), (0.4, "k")):
    ii = np.array([iss(v, gk) for v in vq])
    a.plot(vq, ii, color=col, label=f"tonic K$^+$ {gk} nS")
    z = np.where(np.diff(np.sign(ii)) != 0)[0]
    print(f"D: gK {gk}: zero crossings at " + ", ".join(f"{vq[j]:.1f}" for j in z) + " mV; slope at rest " +
          (f"{(ii[z[0] + 1] - ii[z[0]]) / (vq[1] - vq[0]):.2f} nS" if len(z) else "-"))
a.axhline(0, color="k", lw=0.6); a.set_ylim(-15, 25)
a.text(-88, 20, "no tonic K$^+$: only an unstable crossing (−28 mV);\nthe isolated cell oscillates (period ~10 s)\n0.1 nS or more: stable rest (−56 mV at 0.1, −67 mV at 0.4)",
       fontsize=7, va="top")
a.set_xlabel("V (mV)"); a.set_ylabel("steady-state current (pA)"); a.legend(fontsize=8)
a.set_title("D  Quasi-steady-state I–V", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("paper_R1/fig/supp_smc.png", dpi=150)
print("wrote paper_R1/fig/supp_smc.png")
