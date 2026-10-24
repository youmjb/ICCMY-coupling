"""Isolated ICC-MY variants (iccScale 0.72): VDDR scaled up (key iccVDDR) and the voltage dependence of IP3 production reduced
(key ip3VK; 1 = Youm 2019, 0 = production fixed at its -72 mV value). Does steady hyperpolarization by outward current then
enlarge the pacemaker potential and speed its upstroke without slowing it, as in situ (Kito, Ward & Sanders 2005, pinacidil:
-9.6 mV, amplitude x~1.16, dV/dt up, frequency unchanged)? Mouse ICC-MY in situ control (Kito & Suzuki 2003/2007): rest about -70 mV,
amplitude 56-61 mV, dV/dt 1.5-2.5 V/s, 25-26 min-1.
Usage: [VH=-75] python paper_R1/icc_variant.py run | analyze      Runs: out/paper_icc_hyp/var_V<vddr>_K<vk>_i<inj>.csv (40 s)"""
import os, sys, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from icc_hyperpol import stats, OUT
VH = os.environ.get("VH", "-66")               # VDDR inactivation midpoint (key iccVH); -66 = Youm 2019
VD = ["1", "2", "4", "8"] if VH == "-66" else ["2", "4", "8", "16"]
VK = ["1", "0.5", "0"] if VH == "-66" else ["0"]
INJ = ["0", "25", "35"]
JOBS = list(itertools.product(VD, VK, INJ))


def name(vd, vk, inj):
    return f"{OUT}/var_V{vd}_K{vk}_i{inj}" + ("" if VH == "-66" else f"_H{VH}")


def run(j):
    vd, vk, inj = j; f = name(*j)
    if os.path.exists(f + ".log"): return
    k = ["iccScale=0.72", "L=2", "m=0.3", "fI=0", "fP=0", "smcModel=1", "seed=1", f"iccVDDR={vd}", f"ip3VK={vk}",
         f"iccInj={inj}", "injOn=0", "injOff=1e9", f"iccVH={VH}"]
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", "0", "0", "0", "40000", "0.1", f + ".csv"] + k, check=True, stdout=lg)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        with ThreadPoolExecutor(24) as ex:
            list(ex.map(run, JOBS))
    else:
        print(f"VDDR inactivation midpoint {VH} mV"); print("VDDR x | IP3 V-dep | control: rest, amplitude, dV/dt, rate | hyperpolarized by 25 / 35 pA: d rest, amplitude x, dV/dt x, rate")
        for vd in VD:
            for vk in VK:
                r0 = stats(name(vd, vk, "0") + ".csv")
                if r0[5] < 2: print(f"{vd:>2} | {vk:>3} | no pacemaking (rest {r0[0]:.1f})"); continue
                line = f"{vd:>2} | {vk:>3} | {r0[0]:6.1f} {r0[1]:5.1f} {r0[3]:5.2f} {r0[4]:5.1f} |"
                for inj in INJ[1:]:
                    r = stats(name(vd, vk, inj) + ".csv")
                    line += (f" {r[0] - r0[0]:+5.1f} x{r[1] / r0[1]:.2f} x{r[3] / r0[3]:.2f} {r[4]:5.1f} |" if r[5] >= 2
                             else f" {r[0] - r0[0]:+5.1f} no pacemaking |")
                print(line)
        print("in situ: control -70 mV, 56-61 mV, 1.5-2.5 V/s, 25-26 min-1; pinacidil -9.6 mV, amplitude x~1.16, dV/dt up (~x1.5), rate unchanged")
