"""Preview: recorded against predicted drug responses for the reference model of draft A1 (final_model.out.txt) and model C
(ICC-MY own conductance x0.1 + active muscle C1; final_model_C.out.txt). Values parsed from the 'prediction' lines.
Output: fig/modelC_preview.png"""
import os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import figstyle; figstyle.apply()
PAT = (r"(K0\S+ \S+)\s+g_c fit .*?CSMC ([+-][\d.]+) \(rec ([+-][\d.]+)\) \| ICC ([+-][\d.]+) \(rec ([+-][\d.]+)\), transfer ([\d.]+) "
       r"\(rec ([\d.]+)\) \| ICC x([\d.]+) \(rec ([\d.]+)\) \| CSMC x([\d.]+) \(rec ([\d.]+)\)")


def load(f):
    return {m.group(1): [float(x) for x in m.groups()[1:]] for m in re.finditer(PAT, open(f, encoding="utf-8").read())}


A, C = load("final_model.out.txt"), load("final_model_C.out.txt")
MK = {"K05 pinacidil": "o", "K07a phenylephrine": "^", "K07a apamin": "s", "K07b 4-AP": "D"}
fig, ax = plt.subplots(1, 3, figsize=(figstyle.WIDTH, 2.5))
for a, (i_m, i_r), lab, lim in zip(ax, ((4, 5), (6, 7), (8, 9)), ("resting transfer", "ICC-MY amplitude (fold)", "muscle amplitude (fold)"),
                                    ((0.2, 1.4), (0.85, 1.25), (0.1, 1.3))):
    for k, mk in MK.items():
        a.plot(A[k][i_r], A[k][i_m], mk, mfc="white", mec="0.5", ms=4.5)
        a.plot(C[k][i_r], C[k][i_m], mk, color="k", ms=4.5)
        a.plot([A[k][i_r]] * 2, [A[k][i_m], C[k][i_m]], ":", color="0.6", lw=0.7)
    a.plot(lim, lim, "--", color="0.5", lw=0.7); a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel(f"recorded {lab}"); a.set_ylabel(f"model {lab}")
for k, mk in MK.items(): ax[0].plot([], [], mk, color="k", ms=4, label=k.split()[1])
ax[0].plot([], [], "o", mfc="white", mec="0.5", ms=4, label="draft A1 model"); ax[0].plot([], [], "o", color="k", ms=4, label="model C")
ax[0].legend(frameon=False, fontsize=5.5, loc="upper left")
for a, t in zip(ax, "ABC"): a.set_title(t, loc="left", fontweight="bold")
fig.tight_layout(); figstyle.save(fig, "fig/modelC_preview")
for k in MK: print(k, "A:", A[k][4], A[k][6], A[k][8], "| C:", C[k][4], C[k][6], C[k][8], "| rec:", C[k][5], C[k][7], C[k][9])
