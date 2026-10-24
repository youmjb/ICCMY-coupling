"""Own (non-junctional) resting conductance of ICC-MY in the coupled column: a small steady outward current I into ICC-MY shifts its rest
by dVi and layer 1 by dVm; g_own = (I - g_c (dVi - dVm)) / dVi (chord between slow waves). Model C and model D (final_model, FM_MODEL=C/D, g_c 1.25)
and the draft-A1 reference model (g_c 1.75). Seed 1, 40 s. Usage: python paper_R1/gown_coupled.py run | analyze"""
import os, sys, subprocess, importlib
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import summary
OUT = "out/paper_lowg"
INJ = ["0", "5", "-5"]                                          # outward (+) and inward (-) 5 pA


def base(model):
    os.environ["FM_MODEL"] = model
    import final_model; importlib.reload(final_model); return final_model.BASE


CFG = {"C": ("1.25", base("C")), "A": ("1.75", base("A")), "D": ("1.25", base("D"))}
JOBS = [(m, i) for m in CFG for i in INJ]


def name(m, i): return f"{OUT}/gown_{m}_i{i}"


def run(j):
    m, i = j; n = name(m, i); gc, b = CFG[m]
    if os.path.exists(n + ".log"): return
    with open(n + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "40000", "0.1", n + ".csv"] + b
                       + ["seed=1", f"iccInj={i}", "injOn=0", "injOff=1e9"], stdout=lg, env=dict(os.environ, SMS_CAL="0"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "run":
        with ThreadPoolExecutor(4) as ex: list(ex.map(run, JOBS))
    else:
        for m, (gc, _) in CFG.items():
            a, b = summary(name(m, "0") + ".csv"), summary(name(m, "5") + ".csv")
            dVi, dVm = b["ir"] - a["ir"], b["sr"] - a["sr"]
            g = (-5.0 - float(gc) * (dVi - dVm)) / dVi
            print(f"model {m} (g_c {gc}): dV ICC {dVi:+.2f} mV, layer 1 {dVm:+.2f} mV -> own conductance {g:.2f} nS; "
                  f"passive transfer g_c/(g_c+g_own) {float(gc) / (float(gc) + g):.2f}")
            if os.path.exists(name(m, "-5") + ".csv"):
                c = summary(name(m, "-5") + ".csv"); dVi2, dVm2 = c["ir"] - a["ir"], c["sr"] - a["sr"]
                g2 = (5.0 - float(gc) * (dVi2 - dVm2)) / dVi2
                print(f"    inward 5 pA: dV ICC {dVi2:+.2f} mV, layer 1 {dVm2:+.2f} mV -> own conductance {g2:.2f} nS (depolarizing direction)")
