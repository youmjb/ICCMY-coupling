import os
"""Digitize horizontal bar charts (Kito, Ward & Sanders 2005 Fig. 9A, 9D). For every row below the axis, the contiguous dark run that
starts at the left axis gives the bar end (white bars: their top and bottom borders; black and hatched bars: every row). Rows are
grouped into bars; the bar value is the median run end over the border rows; the whisker end (mean + SD) is the rightmost dark pixel of
the horizontal whisker line, contiguous with the bar end, in the bar's middle rows."""
import sys
import numpy as np
from PIL import Image
sys.stdout.reconfigure(encoding="utf-8")
S = os.environ.get("KITO2005_RENDER_DIR", "kito2005_render")   # folder with the 400-dpi rendering of Fig. 9 of Kito et al. 2005 (not distributed)


def run_end(row, x0):
    x = x0
    while x + 1 < len(row) and (row[x + 1] or (x + 2 < len(row) and row[x + 2])): x += 1
    return x


def analyze(fn, ticks, labels):
    a = np.asarray(Image.open(f"{S}\\{fn}").convert("L")) < 140
    H, W = a.shape
    ax = int(np.argmax(a[: H // 3].sum(1)))
    xs = np.where(a[ax])[0]; xL = xs.min()
    tk = [x for x in range(xL, xs.max() + 1) if a[ax - 6:ax - 1, x].sum() >= 4]
    tk = [x for i, x in enumerate(tk) if i == 0 or x - tk[i - 1] > 3]
    p = np.polyfit(tk, ticks[: len(tk)], 1)
    left = xL + 3                                   # just right of the left axis line
    ends = np.array([run_end(a[y], left) if a[y, left] else left for y in range(H)])
    long_rows = [y for y in range(ax + 8, H) if ends[y] - left > 15]
    groups, g = [], [long_rows[0]]
    for y in long_rows[1:]:
        if y - g[-1] > 8: groups.append(g); g = [y]
        else: g.append(y)
    groups.append(g)
    res = []
    for g, lab in zip(groups, labels):
        e = np.array([ends[y] for y in g])
        bar = np.median(e[e >= np.percentile(e, 25)])
        # whisker: rightmost dark pixel contiguous from the bar end on the middle rows
        ym = range(g[0] + len(g) // 3, g[-1] - len(g) // 3 + 1)
        wk = max(run_end(a[y], int(bar)) if a[y, int(bar)] else int(bar) for y in ym)
        v, w = np.polyval(p, bar), np.polyval(p, wk)
        res.append((lab, v, w - v))
        print(f"   {lab:22s} rows {g[0]}-{g[-1]}: {v:6.2f}  (SD {w - v:5.2f})")
    return res


print("Fig. 9A amplitude (mV); ticks", )
analyze("f9A.png", [0, 10, 20, 30, 40, 50, 60, 70], ["ICC-MY control", "ICC-MY pinacidil", "ICC-MY pin+glib", "CM control", "CM pinacidil", "CM pin+glib"])
print("Fig. 9D dV/dt max (V/s)")
analyze("f9D.png", [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], ["ICC-MY control", "ICC-MY pinacidil", "ICC-MY pin+glib", "CM control", "CM pinacidil", "CM pin+glib"])
