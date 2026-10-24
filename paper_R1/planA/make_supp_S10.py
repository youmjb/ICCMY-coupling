"""Supplement S10 for draft D (planA/supp_S10_D.md): the published ICC model against freshly prepared ICC (Goto et al. 2004).
Table S14 waveform and passive properties (model2019_vs_goto.md values, recomputed here where the run exists); Table S15 input resistance
with the Goto protocol (goto_clamp.out.txt); Table S16 spontaneous currents at -80 mV with SOCE reduced (out/paper_soce).
Citation numbers follow draft D1: Goto et al. 2004 [19], Youm et al. 2019 [15]."""
import os, re, sys
import numpy as np
os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
L = ["## Supplement S10. The published ICC model against freshly prepared ICC",
     "",
     "The ICC model of [15] (25 pF, without the changes made in this study) was compared with the freshly prepared ICC on which it was built "
     "[19], using the recording conditions of [19] where possible. The model reproduces the resting potential and the capacitance; it differs "
     "in the input resistance and its rectification, the rate of rise and peak of the pacemaker potential, and the spontaneous currents under "
     "voltage clamp at a constant potential (Tables S14–S16). The comparison of the input resistance (Table S15) sets the range of resting-"
     "conductance scale factors used for models D and C.",
     ""]
L.append("**Table S14.** The published ICC model (isolated, 25 pF, unmodified) against freshly prepared ICC [19] (35–36 °C, perforated patch).\n")
L.append("| Quantity | Freshly prepared ICC [19] | Model |")
L.append("|---|---|---|")
for q, rec, mod in [
        ("Capacitance", "25.2 ± 4.7 pF (n = 5)", "25 pF"),
        ("Resting potential", "−72.5 ± 6.1 mV (n = 7)", "−72.8 mV"),
        ("Peak of the pacemaker potential", "−1.6 ± 3.8 mV", "−14.4 mV"),
        ("Frequency", "26.8 ± 6.2 min⁻¹", "17.2 min⁻¹"),
        ("Duration at 50 % repolarization", "553 ± 174 ms", "775 ms"),
        ("Maximum rate of rise", "6.7 ± 2.5 V/s at −26.3 mV", "1.24 V/s at −36 mV"),
        ("Input resistance, holding −80 mV (Table S15)", "0.76 ± 0.41 GΩ (n = 3)", "77 MΩ (−20-mV step) to 271 MΩ (+20-mV step)"),
        ("Rectification", "slight inward rectification", "conductance 3.5 times larger for −20 than for +20 mV"),
        ("Spontaneous currents at a constant holding potential (Table S16)", "mostly no full-size inward currents; small currents < 100 pA, tens of ms",
         "full-size inward currents of about 440 pA every 3.6 s at −80 mV")]:
    L.append(f"| {q} | {rec} | {mod} |")
L.append("")
L.append("**Table S15.** Input resistance of the isolated model ICC measured with the protocol of [19]: voltage clamp at −80 mV (series "
         "resistance 2 MΩ), 300-ms steps of −20, +10 and +20 mV in a quiet interval, steady current change against a run without a step. "
         "Resting conductances and steady transporters scaled by f. Published cell: 25 pF, unmodified; cell of this study: 18 pF with the "
         "upstroke and clock of variant (i) (no quiet interval at f = 1 and 0.1, because its clock fires continuously at −80 mV). "
         "Recorded [19]: 0.76 ± 0.41 GΩ (0.7 GΩ in the example cell), 25.2 pF.\n")
L.append("| Cell | f | Holding current at −80 mV (pA) | Input resistance, −20 mV (MΩ) | +10 mV (MΩ) | +20 mV (MΩ) |")
L.append("|---|---|---|---|---|---|")
for line in open("paper_R1/goto_clamp.out.txt", encoding="utf-8"):
    m = re.match(r"\s+(pub|var) f\s+([\d.]+): hold\s+(-?[\d.]+) pA \|\s+(-?\d+) /\s+(-?\d+) /\s+(-?\d+) MOhm", line)
    if not m: continue
    cell = "published (25 pF)" if m.group(1) == "pub" else "this study (18 pF)"
    r20 = m.group(6) if float(m.group(6)) > 0 else "— (active current)"
    L.append(f"| {cell} | {m.group(2)} | {m.group(3)} | {m.group(4)} | {m.group(5)} | {r20} |".replace("-", "−").replace("— (active", "— (active"))
L.append("")
L.append("**Table S16.** Spontaneous inward currents of the published model (25 pF) voltage-clamped at −80 mV for 60 s (last 50 s analysed), "
         "with store-operated Ca²⁺ entry (SOC) or ANO1 scaled from the start of the run; free running for comparison. Events: inward "
         "deviations of more than 30 pA from the median holding current.\n")
L.append("| Condition | Median holding current (pA) | Events per min | Peak (pA) | Half-width (ms) | Free running |")
L.append("|---|---|---|---|---|---|")
for c, lab in (("ctrl", "control"), ("soc0.5", "SOC ×0.5"), ("soc0.25", "SOC ×0.25"), ("soc0", "SOC ×0"), ("ano0", "ANO1 ×0")):
    r = np.loadtxt(f"out/paper_soce/{c}_vc.csv", delimiter=",", skiprows=1); t, i = r[:, 0], r[:, 3]; h = t >= 10000; t, i = t[h], i[h]
    base = np.median(i); thr = base - 30
    on = np.where((i[:-1] >= thr) & (i[1:] < thr))[0]; pk, hw = [], []
    for a in on:
        w = (t >= t[a] - 500) & (t <= t[a] + 3000); ii, tt = i[w], t[w]
        loc = np.median(i[(t >= t[a] - 1500) & (t < t[a] - 200)]) if t[a] - 1500 > t[0] else base
        p = loc - ii.min(); pk.append(p); x = np.where(ii <= loc - p / 2)[0]; hw.append(tt[x[-1]] - tt[x[0]])
    r2 = np.loadtxt(f"out/paper_soce/{c}_cc.csv", delimiter=",", skiprows=1); t2, v = r2[:, 0], r2[:, 1]; v = v[t2 >= 10000]
    paced = v.max() - v.min() > 20
    fr = f"pacing, {v.min():.1f} to {v.max():.1f} mV" if paced else f"no pacing, {np.median(v):.1f} mV"
    L.append(f"| {lab} | {base:.1f} | {60 * len(on) / 50:.1f} | {np.mean(pk) if pk else 0:.0f} | {np.mean(hw) if hw else 0:.0f} | {fr} |"
             .replace("-", "−").replace("| 0 | 0 |", "| — | — |"))
L.append("")
open("paper_R1/planA/supp_S10_D.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L[-8:]))
