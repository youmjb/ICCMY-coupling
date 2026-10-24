"""Reference solution of the ORIGINAL Poh et al. 2012 CellML model (own I_stim drive, 4 slow waves at 10-s period), solved with
scipy LSODA from Myokit-generated Python expressions. Used to verify src/PohSMC.h (see check_poh_cpp in main.cpp, mode 'pohref').
Usage: <python with myokit + scipy> paper_R1/check_poh_ref.py poh_2012.cellml out.csv
"""
import sys, math
import numpy as np
import myokit, myokit.formats
from scipy.integrate import solve_ivp

m = myokit.formats.importer("cellml").model(sys.argv[1])
states = list(m.states()); idx = {s.qname(): i for i, s in enumerate(states)}
w = myokit.formats.ewriter("python")
def lhs(x):
    v = x.var()
    if isinstance(x, myokit.Derivative): return f"dy[{idx[v.qname()]}]"
    if v.is_state(): return f"y[{idx[v.qname()]}]"
    if v.binding() == "time": return "t"
    return v.qname().replace(".", "__")
w.set_lhs_function(lhs)
body = []
for comp, eqs in m.solvable_order().items():
    for eq in eqs:
        v = eq.lhs.var()
        if v.binding() == "time": continue
        body.append(("    " + (f"dy[{idx[v.qname()]}]" if isinstance(eq.lhs, myokit.Derivative) else v.qname().replace('.', '__'))
                     + " = " + w.ex(eq.rhs)))
src = "def rhs(t, y):\n    dy = [0.0] * %d\n%s\n    return dy\n" % (len(states), "\n".join(body))
ns = {"math": math, "exp": math.exp, "log": math.log}
try:
    import numpy
    ns["numpy"] = numpy
except ImportError:
    pass
exec(src, ns)
y0 = [s.initial_value(as_float=True) for s in states]
T = 40000.0
sol = solve_ivp(ns["rhs"], (0, T), y0, method="LSODA", rtol=1e-8, atol=1e-10, max_step=1.0, t_eval=np.arange(0, T, 1.0))
np.savetxt(sys.argv[2], np.c_[sol.t, sol.y[0]], delimiter=",", header="t_ms,V", comments="")
print("ok", sol.status, "V range %.2f .. %.2f" % (sol.y[0].min(), sol.y[0].max()))
