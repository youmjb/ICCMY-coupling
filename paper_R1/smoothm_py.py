"""Python transcription of SmoothM (src/SmoothM.cpp/.h, Coupling copy) for the supplementary model description and voltage-clamp
figures. Same equations, constants, initial values and explicit-Euler / exponential-gate scheme as the C++ code (all SMS_* scales 1).
Units: V mV, t ms, I pA, g nS, Cai uM. Check: run as a script -> compares a current-clamp step with out/paper_smc/step_gK0.4.csv.
"""
import math
import numpy as np

RTF, CM = 26.7257, 50.0
KE, KI = 5.40, 109.0
EK = RTF * math.log(KE / KI)                       # -80.3 mV (recomputed each step in C++ from the same constants)
SUBV = math.pi * (6e-4 ** 2 - 5e-4 ** 2) * 200e-4  # cm^3, 1-um submembrane shell of a 200 x 12 um cylinder
SACC = math.pi * 200e-4 * 12e-4                    # cm^2
COMPLOSS, F = 7e-6, 96.4867
GCAL, ECAK, KCA = 65.0, 87.0, 9.0
GKS, GKCA, GLEAK, ELEAK = 3.1, 0.7, 5e-3 * CM, -46.3


def cal_rates(v):
    x = v + 15.0
    a = 13 * x / (1 - math.exp(-x / 3)) if abs(x) > 1e-9 else 13 * 3.0
    b = 31 * math.exp(-0.04 * x)
    return a, b


def ks_rates(v):
    x = v + 8.821
    a = -1.126 * x / (math.exp(-0.386 * x) - 1) if abs(x) > 1e-9 else 1.126 / 0.386
    b = 17.02 * math.exp(-0.022 * (v + 60))
    return a, b


def gates_inf_tau(v):
    """steady states and time constants (ms) of every gate at voltage v"""
    a, b = cal_rates(v); d_inf, d_tau = a / (a + b), 1e3 / (a + b)
    f_inf = 1 / (1 + math.exp((v + 40.5) / 13.1))
    f1_tau = 28 / (1 + math.exp(-(v + 20) / 1))
    f2_tau = 300 / (1 + math.exp(-(v + 60) / 2))
    a, b = ks_rates(v); n_inf, n_tau = a / (a + b), 1e3 / (a + b)
    return d_inf, d_tau, f_inf, f1_tau, f2_tau, n_inf, n_tau


def kca_h(v, cai):
    kn = -65 * cai + 119.44
    s = 14.372 * cai + 11.566
    return 1 / (1 + math.exp((kn - v) / s))


class SmoothM:
    def __init__(self):
        self.v, self.cai = -50.0, 0.1
        self.d, self.f1, self.f2, self.n = 0.0, 1.0, 1.0, 0.0

    def currents(self, v=None):
        v = self.v if v is None else v
        ical = GCAL * self.d * self.f1 * self.f2 / (1 + self.cai / KCA) * (v - ECAK)
        iks = GKS * self.n ** 2 * (v - EK)
        ikca = GKCA * kca_h(v, self.cai) * (v - EK)
        ileak = GLEAK * (v - ELEAK)
        return ical, iks, ikca, ileak

    def step(self, dt, igap=0.0, vclamp=None):
        v = self.v if vclamp is None else vclamp
        if vclamp is not None: self.v = v
        d_inf, d_tau, f_inf, f1_tau, f2_tau, n_inf, n_tau = gates_inf_tau(v)
        self.d = d_inf - (d_inf - self.d) * math.exp(-dt / d_tau)
        self.f1 = f_inf - (f_inf - self.f1) * math.exp(-dt / f1_tau)
        self.f2 = f_inf - (f_inf - self.f2) * math.exp(-dt / f2_tau)
        self.n = n_inf - (n_inf - self.n) * math.exp(-dt / n_tau)
        ical, iks, ikca, ileak = self.currents(v)
        k = 1e-1 * COMPLOSS / SACC                                   # per ms; core [Ca] fixed at 0.1 uM
        self.cai += (-1e-9 * ical / 2 / SUBV / F - k * self.cai + k * 0.1) * dt
        if vclamp is None:
            self.v += -(ical + iks + ikca + ileak + igap) / CM * dt
        return ical, iks, ikca, ileak


def equilibrate(gk, T=20000.0, dt=0.01):
    c = SmoothM()
    for _ in range(int(round(T / dt))):
        c.step(dt, igap=gk * (c.v + 80.0))
    return c


if __name__ == "__main__":
    import csv, copy
    gk, dt = 0.4, 0.01
    rest = equilibrate(gk)
    print(f"rest {rest.v:.2f} mV (C++ smsstep: -67.0)")
    ref = list(csv.DictReader(open("out/paper_smc/step_gK0.4.csv")))
    for amp in (10.0, 20.0):
        c = copy.deepcopy(rest); v = []
        for i in range(int(1000 / dt)):
            t = i * dt
            c.step(dt, igap=gk * (c.v + 80.0) - (amp if 100 <= t < 1100 else 0.0))
            if i % 50 == 0: v.append(c.v)
        r = np.array([float(x[f"A{amp:g}"]) for x in ref])
        n = min(len(r), len(v))
        print(f"{amp:g} pA step: max |python - C++| = {np.max(np.abs(np.array(v[:n]) - r[:n])):.3f} mV over {n} samples")
