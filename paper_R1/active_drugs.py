"""Out-of-sample drug predictions for the active-SMC candidates of active_cal.py (control state fitted only).
Each drug is sized by the recorded CSMC resting-potential change only; amplitudes and the ICC-MY response are predictions.
  pinacidil (Kito, Ward & Sanders 2005): SMC K+ added (gK0run base + add, add 1/2/4/8 nS); pick the add closest to -10.3 mV.
      Recorded: CSMC -10.3, ICC-MY -9.6 mV; CSMC amplitude x~0.2, ICC-MY x~1.16 (Fig. 9A, read from the chart).
  apamin (Kito & Suzuki 2007): SK block in FLC (skKd 1000), no fitting.
      Recorded: CSMC +3.1, ICC-MY +3.3 mV; CSMC primary 36.2 -> 36.0 (x1.00), ICC-MY 61.3 -> 58.8 (x0.96).
  4-AP (Kito & Suzuki 2007): SMC K+ reduced (gK0run 0.85/0.75/0.5/0.25/0 of base); pick the one closest to +3.1 mV.
      Recorded: CSMC +3.1, ICC-MY +3.9 mV; CSMC primary 37.0 -> 45.3 (x1.22), dV/dt x1.31; ICC-MY 60.0 -> 56.5 (x0.94).
Usage: python paper_R1/active_drugs.py run | analyze     Runs: out/paper_active/drug_<cand>_<drug>.csv (20 s)"""
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, "paper_R1")
from active_cal import B, OUT, name as cal_name, summary
CANDS = {"C1": ("1.5", "2", "-50", "-65", "1"), "C2": ("1.5", "4", "-45", "-65", "1"), "C3": ("1", "2", "-50", "-60", "1"),
         "P": ("2.33", "0", "-50", "-65", "0.4")}               # P: the passive calibrated column (gT 0) for comparison
DRUGS = {f"pin{a}": {"gK0run": None, "add": a} for a in ("1", "2", "4", "8")}
DRUGS.update({"apa": {"skKd": "1000"}})
DRUGS.update({f"4ap{x}": {"scale": x} for x in ("0.85", "0.75", "0.5", "0.25", "0")})


def keys(c, d):
    gc, gT, vmT, vhT, gK = CANDS[c]; spec = DRUGS[d]; k = {"gT": gT, "vmT": vmT, "vhT": vhT, "gK0run": gK}
    if "add" in spec: k["gK0run"] = str(float(gK) + float(spec["add"]))
    if "scale" in spec: k["gK0run"] = str(float(gK) * float(spec["scale"]))
    if "skKd" in spec: k["skKd"] = spec["skKd"]
    return gc, [f"{a}={b}" for a, b in k.items()]


def run(j):
    c, d = j; f = f"{OUT}/drug_{c}_{d}"
    if os.path.exists(f + ".log"): return
    gc, k = keys(c, d)
    with open(f + ".log", "w") as lg:
        subprocess.run(["./bin/gisim.exe", "active", gc, "0", "0", "20000", "0.1", f + ".csv"] + B + k, check=True, stdout=lg)


def control_file(c):
    if c == "P":
        f = f"{OUT}/drug_P_ctrl"
        if not os.path.exists(f + ".log"):
            with open(f + ".log", "w") as lg:
                subprocess.run(["./bin/gisim.exe", "active", "2.33", "0", "0", "20000", "0.1", f + ".csv"] + B + ["gK0run=0.4"], check=True, stdout=lg)
        return f + ".csv"
    return cal_name(*CANDS[c]) + ".csv"


if __name__ == "__main__":
    if sys.argv[1] == "run":
        os.environ["OMP_NUM_THREADS"] = "1"
        control_file("P")
        with ThreadPoolExecutor(24) as ex:
            list(ex.map(run, [(c, d) for c in CANDS for d in DRUGS]))
    else:
        for c in CANDS:
            s0 = summary(control_file(c)); a0 = s0["a"]
            print(f"{c} {CANDS[c]}: control CSMC rest {s0['sr']:.1f}, ICC rest {s0['ir']:.1f}, ICC {a0['icc']:.1f}, CSMC {a0['smc']:.1f} mV, "
                  f"ratio {a0['smc'] / a0['icc']:.2f}, dV/dt {s0['dv']:.2f}")
            for d in DRUGS:
                f = f"{OUT}/drug_{c}_{d}.csv"
                if not os.path.exists(f): continue
                s = summary(f); a = s["a"]
                if a is None: print(f"   {d:7s} no slow waves"); continue
                print(f"   {d:7s} dCSMC {s['sr'] - s0['sr']:+5.1f}  dICC {s['ir'] - s0['ir']:+5.1f} | ICC x{a['icc'] / a0['icc']:.2f}  "
                      f"CSMC x{a['smc'] / a0['smc']:.2f}  dV/dt x{s['dv'] / s0['dv']:.2f} | ratio {a['smc'] / a['icc']:.2f}")
        print("recorded: pinacidil dCSMC -10.3 dICC -9.6, ICC x~1.16, CSMC x~0.2 | apamin +3.1/+3.3, ICC x0.96, CSMC x1.00 | "
              "4-AP +3.1/+3.9, ICC x0.94, CSMC x1.22, dV/dt x1.31")
