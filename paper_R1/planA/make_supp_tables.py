"""New supplementary tables for plan A, computed from the run outputs (no hand transcription). Writes planA/supp_tables_A.md, or with
FM_MODEL=C planA/supp_tables_C.md (drug levels of model C, own-conductance tests as Table S9, no note on the resting conductance)."""
import os, sys
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, "paper_R1")
from icc_hyperpol import stats
import final_model as FM
from active_cal import summary
def own_conductance_table():
    """Table S9 (model C only): own conductance of ICC-MY lowered (iccRest f), isolated and in the column (icc_lowg.py, lowg_coupled.py),
    and with the muscle variants that lead to model C (lowg_rest.py, lowg_active.py); seed 1, pinacidil interpolated to -10.3 mV."""
    import icc_lowg as IL, lowg_coupled as LC, lowg_rest as LR, lowg_active as LA
    from fig_lowg import sized
    T = ["**Table S9.** Own conductance of ICC-MY lowered by scaling its resting conductances and steady transporters together by f (seed 1). "
         "Isolated cell: own conductance from a 5-pA outward current, resting potential, amplitude and rate. Column: passive muscle (tonic K⁺ 0.4 nS), "
         "g\\_c 1.75 nS; at the pinacidil level (muscle −10.3 mV, K\\_ATP interpolated): transfer of the resting change and amplitude changes. "
         "Lower part: the same with the muscle variants that lead to model C (tonic K⁺ raised; active muscle), at the g\\_c closest to a control ratio of 0.52; "
         "—, −10.3 mV not reached. Main Fig. 5A.\n",
         "| f | Isolated: own conductance (nS) | Isolated: rest (mV) | Isolated: amplitude (mV) | Isolated: rate (min⁻¹) | Column: control ratio | Column: ICC-MY rest (mV) | Column: ICC-MY amplitude (mV) | Transfer | ICC-MY amplitude (fold) | Muscle amplitude (fold) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for f in LC.F:
        r0, r5 = stats(IL.name(f, "0") + ".csv"), stats(IL.name(f, "5") + ".csv")
        paced = r0[0] < -50
        g = f"{5.0 / (r0[0] - r5[0]):.2f}" if paced else "—"
        iso = f"{r0[0]:.1f} | {r0[1]:.1f} | {r0[4]:.1f}" if paced else f"{r0[0]:.1f} (no pacing) | — | —"
        c0 = LC.S(f, "ctrl", "1.75"); (tr, xi, xm), _ = sized(f)
        T.append(f"| {f} | {g} | {iso} | {c0['a']['smc'] / c0['a']['icc']:.2f} | {c0['ir']:.1f} | {c0['a']['icc']:.1f} | {tr:.2f} | {xi:.2f} | {xm:.2f} |".replace("-", "−"))
    T.append("\n| Muscle variant | f | g\\_c (nS) | Control ratio | ICC-MY rest (mV) | Muscle rest (mV) | ICC-MY amplitude (mV) | Transfer | ICC-MY amplitude (fold) | Muscle amplitude (fold) |")
    T.append("|---|---|---|---|---|---|---|---|---|---|")

    def interp(rows):
        rows.sort()
        x = [r[0] for r in rows]
        if not rows or not (x[0] <= -10.3 <= x[-1]): return None
        return [float(np.interp(-10.3, x, [r[i] for r in rows])) for i in (1, 2, 3)]

    for lab, f, k, gcs, Sf, pins in [("tonic K⁺ 1.2 nS", fv, "1.2", LR.GC, None, LR.PIN) for fv in LR.F] + \
                                     [("active muscle (model C muscle)", fv, None, LA.GC, None, LA.PIN) for fv in LA.F]:
        get = (lambda c, gc: LR.S(f, k, c, gc)) if k else (lambda c, gc: LA.S(f, c, gc))
        best = min((gc for gc in gcs if get("ctrl", gc)), key=lambda gc: abs(get("ctrl", gc)["a"]["smc"] / get("ctrl", gc)["a"]["icc"] - 0.52))
        c0 = get("ctrl", best); rows = []
        for p in pins:
            c1 = get(f"pin{p}", best)
            if c1:
                dM = c1["sr"] - c0["sr"]
                rows.append((dM, (c1["ir"] - c0["ir"]) / dM, c1["a"]["icc"] / c0["a"]["icc"], c1["a"]["smc"] / c0["a"]["smc"]))
        v = interp(rows)
        tail = f"{v[0]:.2f} | {v[1]:.2f} | {v[2]:.2f}" if v else "— | — | —"
        T.append(f"| {lab} | {f} | {best} | {c0['a']['smc'] / c0['a']['icc']:.2f} | {c0['ir']:.1f} | {c0['sr']:.1f} | {c0['a']['icc']:.1f} | {tail} |".replace("-", "−").replace("— | — | —", "— | — | —"))
    T.append("\nRecorded: transfer 0.93, ICC-MY ×1.17, muscle ×0.22. Seeds 1–3: model D (f 0.3) 0.55, ×1.06, ×0.40; model C (f 0.1) 0.79, ×1.14, ×0.51.")
    return T


def cluster_table():
    """Table S10 (model D): simulated current-clamp recordings from isolated ICC and clusters (cluster_pulse.py)."""
    import cluster_pulse as CP, cluster_sim as CS
    T = ["**Table S10.** Simulated recordings from isolated ICC (5.5 pF each, upstroke and clock as in models C and D, no muscle): N cells "
         "coupled as a chain or all-to-all by g\\_II; 2-pA hyperpolarizing pulse for 400 ms between slow waves in cell 0. Input resistance from the "
         "voltage change at 400 ms, time constant τ from the time to 63 % of it, apparent capacitance = τ / input resistance. Recorded: cultured "
         "ICC [21] 29 pF, 250 MΩ (τ ≈ 7 ms), −53 mV; freshly prepared ICC [19] 25.2 pF, 0.76 ± 0.41 GΩ (τ ≈ 19 ms), −72.5 mV. With f = 0.1 the cells settle near −21 mV, so that τ there is not that of a resting cell. "
         "Main Fig. 6.\n",
         "| Resting currents (f) | N | Coupling | g\\_II (nS) | True C (pF) | Input resistance (MΩ) | τ (ms) | Apparent C (pF) | Rest (mV) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for f, n, topo, g in CP.CFGS:
        if n != "1" and g not in ("1", "3", "10"): continue
        R, tau, C, dV = CP.measure((f, n, topo, g)); r0 = CS.cc_stats(CS.name(f, n, topo, g, "cc"))
        T.append(f"| {f} | {n} | {'—' if n == '1' else ('chain' if topo == '0' else 'all-to-all')} | {'—' if n == '1' else g} | "
                 f"{5.5 * int(n):.1f} | {R:.0f} | {tau:.0f} | {C:.1f} | {r0[0]:.1f} |".replace("-", "−").replace("—", "—"))
    return T


L = []
H = "out/paper_icc_hyp"

# Table S5: isolated ICC-MY variants (current 35 pA)
L.append("**Table S5.** Isolated ICC-MY (18 pF) in control and hyperpolarized by 35 pA outward current: VDDR scale, VDDR inactivation midpoint "
         "and voltage dependence of IP₃ production (1, original; 0, fixed at its −72 mV value). Fold changes relative to the control of the same variant.\n")
L.append("| VDDR scale | Midpoint (mV) | IP₃ V-dependence | Rest (mV) | Amplitude (mV) | dV/dt (V/s) | Rate (min⁻¹) | ΔV at 35 pA (mV) | Amplitude (fold) | dV/dt (fold) | Rate at 35 pA (min⁻¹) |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
rows = [(vd, "-66", vk, f"{H}/var_V{vd}_K{vk}_i{{i}}") for vd in ("1", "2", "4", "8") for vk in ("1", "0.5", "0")]
rows += [(vd, vh, "0", f"{H}/var_V{vd}_K0_i{{i}}_H{vh}") for vh in ("-70", "-75", "-80") for vd in ("2", "4", "8", "16")]
for vd, vh, vk, pat in rows:
    r0 = stats(pat.format(i="0") + ".csv"); r1 = stats(pat.format(i="35") + ".csv")
    L.append(f"| {vd} | {vh.replace('-', '−')} | {vk} | {r0[0]:.1f} | {r0[1]:.1f} | {r0[3]:.2f} | {r0[4]:.1f} | {r1[0] - r0[0]:+.1f} | "
             f"{r1[1] / r0[1]:.2f} | {r1[3] / r0[3]:.2f} | {r1[4]:.1f} |".replace("-", "−").replace("+", "+"))
L.append("\nIn situ (pinacidil, [5]): −9.6 mV, amplitude ×1.17, dV/dt ×1.52, rate unchanged. Reference variant: VDDR ×16, midpoint −80 mV, IP₃ V-dependence 0.\n")

# Table S6: isolated ICC-MY with an added K+ conductance (E -90 mV)
L.append("**Table S6.** Isolated ICC-MY hyperpolarized by an added K⁺ conductance (reversal −90 mV) instead of current: a direct action of a K⁺-channel "
         "opener on ICC-MY shunts the upstroke.\n")
L.append("| Model | Added conductance (nS) | ΔV (mV) | Amplitude (fold) | dV/dt (fold) | Rate at 35 pA (min⁻¹) |")
L.append("|---|---|---|---|---|---|")
for sc in ("1", "0.72"):
    r0 = stats(f"{H}/s{sc}_i0.csv")
    for g in ("0.5", "1", "2", "4"):
        r = stats(f"{H}/s{sc}_g{g}.csv")
        L.append(f"| original, {'25 pF' if sc == '1' else '18 pF'} | {g} | {r[0] - r0[0]:+.1f} | {r[1] / r0[1]:.2f} | {r[3] / r0[3]:.2f} | {r[4]:.1f} |".replace("-", "−"))
L.append("")

# Table S7: drug sizing in the reference model
MOD = {"C": "model C", "D": "model D"}.get(FM.MODEL, "the reference model")
L.append(f"**Table S7.** Drug levels in {MOD} (mean of seeds 1–3, each relative to the control of the same seed): change of the muscle "
         "and ICC-MY resting potentials, their ratio (transfer) and the slow-wave amplitude changes. The level matching the recorded muscle change was "
         "found by linear interpolation (main Table 2).\n")
L.append("| Drug level | g\\_c (nS) | ΔV muscle (mV) | ΔV ICC-MY (mV) | Transfer | ICC-MY amplitude (fold) | Muscle amplitude (fold) |")
L.append("|---|---|---|---|---|---|---|")
lab = lambda d: (f"K\_ATP {d[3:]} nS" if d.startswith("pin") else "SK block" if d == "apa" else
                 f"tonic K⁺ ×{d[3:]}" if d.startswith("4ap") else f"K⁺ {d[2:]} nS (phenylephrine)")
for d in FM.DRUGS:
    for gc in FM.GC_DRUG:
        e = [x for s in FM.SEEDS if (x := FM.effect(d, gc, s))]
        if not e: continue
        m = np.mean(e, axis=0)
        L.append(f"| {lab(d)} | {gc} | {m[0]:+.1f} | {m[1]:+.1f} | {m[1] / m[0]:.2f} | {m[2]:.2f} | {m[3]:.2f} |".replace("-", "−"))
L.append("")

# Table S8: coupling and muscle variants at the pinacidil level
L.append("**Table S8.** Coupling and muscle variants at the pinacidil level (muscle −10.3 mV, K\\_ATP reversal −84 mV, ICC-MY variant, nifedipine, seed 1; "
         "values interpolated between K\\_ATP 1–8 nS). Main Fig. 5.\n")
L.append("| Variant | g\\_c (nS) | Control ratio | Transfer | ICC-MY amplitude (fold) | Muscle amplitude (fold) |")
L.append("|---|---|---|---|---|---|")
for line in open("paper_R1/figA.out.txt", encoding="utf-8"):
    if line.startswith("Fig5"):
        name, rest = line[5:].split(" g_c ")
        gc, vals = rest.split(":", 1)
        v = [x.strip() for x in vals.split(",")]
        nm = name.replace("$g_c$", "g\\_c").replace("$g_{cL}$", "g\\_cL").replace("$", "")
        L.append(f"| {nm} | {gc} | {v[0].split()[-1]} | {v[1].split()[-1]} | "
                 f"{v[2].split('x')[-1]} | {v[3].split('x')[-1].strip()} |")
L.append("\nRecorded: transfer 0.93, ICC-MY ×1.17, muscle ×0.22.\n")

# ICC-MY resting conductance reduced (text)
if FM.MODEL in ("C", "D"): L.extend(own_conductance_table()); L.append("")
if FM.MODEL == "D": L.extend(cluster_table()); L.append("")
L.append("**Reducing the resting conductance of ICC-MY.** A negative passive conductance at −72 mV (−1 to −3 nS) was added to ICC-MY of the original-conductance model in "
         "the active-muscle column. At −1 nS the cell either latched near +20 mV (g\\_c 1 nS) or, at g\\_c 1.5–2 nS, paced "
         "but exploded with apamin; at −2 and −3 nS it did not pace or the integration diverged. A lower own conductance therefore needs a different "
         "balance of the resting currents, not a simple leak change"
         + (": scaling the resting conductances and the transporters that balance them together (Table S9) lowers the own conductance while the "
            "isolated cell keeps pacing down to f = 0.2.\n" if FM.MODEL in ("C", "D") else ".\n"))
open(f"paper_R1/planA/supp_tables_{FM.MODEL}.md" if FM.MODEL in ("C", "D") else "paper_R1/planA/supp_tables_A.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L[:12]))
