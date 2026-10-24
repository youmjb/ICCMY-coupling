"""95 % intervals of the recorded resting-potential transfer (ICC-MY change / circular-muscle change), Methods 'Data' and Table 2.

The two changes were recorded in different tissues, so each mean change is drawn independently from a normal distribution with its
recorded mean and standard error s/sqrt(n) (n = animals); the interval is the 2.5th-97.5th percentile of the ratio of the draws.
Where only the SDs before (s0) and in the drug (s1) are reported (apamin, 4-AP; Kito & Suzuki 2007 Table 1), the SD of the change is
s = sqrt(s0^2 + s1^2 - 2 rho s0 s1) with an assumed within-tissue correlation rho.
Data: kito_data.md (Supplement S9). Output: transfer_intervals.out.txt.   Usage: python paper_R1/transfer_intervals.py"""
import os
import numpy as np

N, SEED = 1_000_000, 1
SDC = lambda s0, s1, rho: np.sqrt(s0 * s0 + s1 * s1 - 2 * rho * s0 * s1)
# drug: (ICC mean change, ICC SD of change, n ICC, muscle mean change, muscle SD of change, n muscle); SDs of change may depend on rho
DATA = {
    "pinacidil": lambda rho: (-9.6, 2.0, 8, -10.3, 1.0, 11),                              # Kito et al. 2005 (text)
    "phenylephrine": lambda rho: (-1.6, 1.3, 9, -2.4, 2.2, 9),                            # Kito & Suzuki 2007 (text)
    "apamin": lambda rho: (+3.3, SDC(2.9, 2.8, rho), 7, +3.1, SDC(2.8, 2.5, rho), 15),     # Kito & Suzuki 2007 Table 1
    "4-AP": lambda rho: (+3.9, SDC(3.4, 4.0, rho), 8, +3.1, SDC(1.2, 1.5, rho), 12),       # Kito & Suzuki 2007 Table 1
}


def interval(mi, si, ni, mm, sm, nm, rng):
    r = rng.normal(mi, si / np.sqrt(ni), N) / rng.normal(mm, sm / np.sqrt(nm), N)
    return np.percentile(r, [2.5, 97.5])


lines = [f"95 % intervals of the transfer (Monte Carlo, {N} draws, seed {SEED}, percentile interval)"]
for drug, f in DATA.items():
    rhos = (None,) if drug in ("pinacidil", "phenylephrine") else (0.8, 0.95, 0.0)
    for rho in rhos:
        rng = np.random.default_rng(SEED)
        mi, si, ni, mm, sm, nm = f(rho)
        lo, hi = interval(mi, si, ni, mm, sm, nm, rng)
        tag = "" if rho is None else f" (rho {rho})"
        lines.append(f"{drug}{tag}: transfer {mi / mm:.2f}, 95 % interval {lo:.2f}-{hi:.2f}")
out = "\n".join(lines)
print(out)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "transfer_intervals.out.txt"), "w", encoding="utf-8").write(out + "\n")
