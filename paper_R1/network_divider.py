"""Passive ICC-MY network divider for inhibitory input (reviewer point 1-6).

Square lattice (N x N) of ICC-MY nodes coupled to their four neighbours by g_II (unknown; swept). Each node has its own
conductance g_own and is coupled by g_c to layer 1 of its own muscle column. In columns that receive the muscle IJP, layer 1 is
treated as a voltage source dV_m = 1 (the muscle IJP the full model reports); in the other columns the passive column loads the
node through g_c*G_load/(g_c + G_load). The fraction is the hyperpolarization of the central ICC-MY relative to dV_m.
Widespread input (every column driven, e.g. field stimulation) gives g_c/(g_c + g_own) for any g_II, i.e. the single-column result;
local input (one column or a patch) is diluted by current flowing into neighbouring ICC-MY.
Values: g_c 1.5 nS, g_own 3.6 nS, G_load 1.61 nS (calibrated column). Output: printed table; data() used by figs_assemble.py.
"""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

GC, GOWN, GLOAD, N = 1.5, 3.6, 1.61, 41
GII = np.logspace(-1, 2, 25)
PATCH = [("one column", 0), ("3 × 3 columns", 1), ("7 × 7 columns", 3), ("widespread", None)]


def fraction(gii, half, gc=GC, gown=GOWN, gload=GLOAD, n=N):
    c = n // 2
    A = lil_matrix((n * n, n * n)); b = np.zeros(n * n)
    gpass = gc * gload / (gc + gload)
    for i in range(n):
        for j in range(n):
            k = i * n + j
            driven = half is None or (abs(i - c) <= half and abs(j - c) <= half)
            A[k, k] += gown + (gc if driven else gpass)
            if driven:
                b[k] += gc
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:          # sealed edge
                    A[k, k] += gii; A[k, ii * n + jj] -= gii
    v = spsolve(A.tocsr(), b)
    return v[c * n + c]


def data():
    return {lab: np.array([fraction(g, h) for g in GII]) for lab, h in PATCH}


if __name__ == "__main__":
    d = data()
    print("single-node divider g_c/(g_c+g_own):", round(GC / (GC + GOWN), 3))
    print("g_II (nS) | " + " | ".join(lab for lab, _ in PATCH))
    for k, g in enumerate(GII):
        if k % 4 == 0:
            print(f"{g:7.2f} | " + " | ".join(f"{d[lab][k]:.3f}" for lab, _ in PATCH))
