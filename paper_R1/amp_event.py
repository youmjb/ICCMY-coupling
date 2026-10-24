"""Event-aligned slow-wave amplitudes, the primary amplitude measure of the paper from revision 2 (3rd review, 3-2).

The run logs of gisim report window amplitudes (max - min over the second half of the run), which include FLC/SK-driven troughs
that occur between slow waves. Here, for every ICC-MY upstroke in the second half of the run (crossing of the ICC-MY mid-level,
min + 0.5 (max - min)), the amplitude of ICC-MY and of muscle layer 1 is peak (within 2 sec) minus baseline (mean over 50-200 msec
before the upstroke), averaged over events.

log(path) returns the text of a run log with the ICC-MY amplitude ("ICC-MY .. cpm amp X") and the layer-1 amplitude ("[0]Y")
replaced by the event-aligned values computed from the run's trace (<prefix>.csv), so the analysis scripts keep their parsers.
Runs without a trace or without detectable upstrokes keep the window values (counted in FALLBACK); a negative event amplitude
(no wave in an uncoupled layer) is written as 0. Results are cached in
out/amp_event_cache.json (keyed by trace path, size and mtime).
"""
import os, re, csv, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_F = os.path.join(ROOT, "out", "amp_event_cache.json")
try:
    _cache = json.load(open(CACHE_F))
except (OSError, ValueError):
    _cache = {}
FALLBACK = []


def measure(trace):
    """event-aligned amplitudes {icc, smc, n, base, peak, trough} of one trace, or None"""
    p = trace if os.path.isabs(trace) else os.path.join(ROOT, trace)
    st = os.stat(p); key = f"{os.path.relpath(p, ROOT)}|{st.st_size}|{int(st.st_mtime)}"
    if key in _cache: return _cache[key]
    r = csv.reader(open(p)); h = next(r)
    if "V_SMC0" not in h: _cache[key] = None; return None
    iI, iS = h.index("V_ICC"), h.index("V_SMC0")
    a = np.array([[float(x[0]), float(x[iI]), float(x[iS])] for x in r]); a = a[len(a) // 2:]
    t, vi, vs = a[:, 0], a[:, 1], a[:, 2]; dt = t[1] - t[0]
    mid = vi.min() + 0.5 * (vi.max() - vi.min())
    up = np.where((vi[:-1] < mid) & (vi[1:] >= mid))[0] + 1
    up = up[(t[up] - t[0] > 300) & (t[-1] - t[up] > 2100)]
    res = None
    if vi.max() - vi.min() > 1 and len(up):
        n200, n50, n2s = int(200 / dt), int(50 / dt), int(2000 / dt)
        ei = [vi[u:u + n2s].max() - vi[u - n200:u - n50].mean() for u in up]
        es = [vs[u:u + n2s].max() - vs[u - n200:u - n50].mean() for u in up]
        res = dict(icc=float(np.mean(ei)), smc=float(np.mean(es)), n=int(len(up)),
                   base=float(np.mean([vs[u - n200:u - n50].mean() for u in up])), peak=float(np.mean([vs[u:u + n2s].max() for u in up])),
                   trough=float(vs.min()), win_icc=float(vi.max() - vi.min()), win_smc=float(vs.max() - vs.min()))
    _cache[key] = res
    return res


def save():
    os.makedirs(os.path.dirname(CACHE_F), exist_ok=True)
    json.dump(_cache, open(CACHE_F, "w"))


def log(path):
    """log text with event-aligned ICC-MY and layer-1 amplitudes (see module doc)"""
    p = path if os.path.isabs(path) else os.path.join(ROOT, path)
    txt = open(p).read()
    tr = p[:-4] + ".csv"
    m = measure(tr) if os.path.exists(tr) else None
    if m is None:
        FALLBACK.append(os.path.relpath(p, ROOT)); return txt
    txt = re.sub(r"(ICC-MY [\d.]+ cpm amp )([\d.]+)", lambda x: x.group(1) + f"{max(m['icc'], 0.0):.2f}", txt, count=1)
    txt = re.sub(r"(amp by layer: \[0\])([\d.]+)", lambda x: x.group(1) + f"{max(m['smc'], 0.0):.2f}", txt, count=1)
    return txt


import atexit
atexit.register(save)
