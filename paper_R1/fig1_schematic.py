"""Fig. 1 schematic for the R1 paper (pure drawing, no data).
A  Anatomy: longitudinal muscle, ICC-MY network, circular muscle; one ICC-MY serves one column of circular muscle
   (Publicover 1995 argument); the column width is set by the ICC-MY density (m SMC per ICC-MY per layer).
B  Model: ICC-MY node (own conductance g_own) coupled by g_c to layer 1 of a 30-layer cable; each layer = m SMC + fP*m FLC,
   layers coupled by g_SS; purinergic input acts on FLC (SK); recording/injection sites used as constraints.
Output: paper_R1/fig/fig1_schematic.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse, FancyArrowPatch

plt.rcParams["font.family"] = "Arial"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig", "fig1_schematic.png")
ICC, SMC, FLC, INK, MUTED = "#D9822B", "#C9DCE6", "#BFD9C4", "#1B2440", "#6B7280"
LONG, RED, BOX = "#E8E2D6", "tab:red", "#F4F1EA"          # longitudinal muscle fill, purinergic arrow (figs_assemble sets grays for the main figure)

def draw(a, b, k=1.0, title_a="A  One ICC-MY drives one column of circular muscle", title_b="B  Model"):
    """Draw panel A (anatomy) into axes a and panel B (model) into axes b; k scales all font sizes."""
    # ---------------- A: anatomy ----------------
    a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
    a.text(0, 9.7, title_a, fontsize=max(11*k, 8.5), fontweight="bold", va="top")
    a.add_patch(Rectangle((0.5, 7.6), 9, 1.0, fc=LONG, ec="none"))
    a.text(0.6, 8.1, "longitudinal muscle", fontsize=max(8*k, 6.5), va="center", color=MUTED)
    a.add_patch(Rectangle((0.5, 1.0), 9, 5.9, fc=SMC, ec="none", alpha=0.6))
    a.text(0.6, 1.3, "circular muscle", fontsize=max(8*k, 6.5), color=MUTED)
    xs = [1.4, 3.2, 5.0, 6.8, 8.6]
    for i, x in enumerate(xs):
        a.add_patch(Ellipse((x, 7.25), 1.0, 0.45, fc=ICC, ec=INK, lw=0.8))
        if i < len(xs) - 1:
            a.plot([x + 0.5, xs[i + 1] - 0.5], [7.25, 7.25], color=ICC, lw=2.5)
    a.text(0.6, 6.6, "ICC-MY\nnetwork", fontsize=max(8*k, 6.5), va="top", ha="left")
    x0 = xs[2]
    a.add_patch(Rectangle((x0 - 0.9, 1.0), 1.8, 5.9, fc="none", ec=INK, lw=1.2, ls="--"))
    for y in (6.6, 6.0, 5.4, 4.8):
        a.add_patch(Rectangle((x0 - 0.75, y - 0.2), 1.5, 0.4, fc="white", ec=INK, lw=0.6))
    a.plot([x0] * 3, [3.6, 3.95, 4.3], "o", color=INK, ms=2.5*k)
    a.add_patch(Rectangle((x0 - 0.75, 1.3), 1.5, 0.4, fc="white", ec=INK, lw=0.6))
    a.annotate("", xy=(x0 - 0.9, 0.55), xytext=(x0 + 0.9, 0.55), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8))
    a.text(x0, 0.15, "column width set by ICC-MY density", ha="center", fontsize=max(8*k, 6.5), va="center")
    a.text(x0 + 1.05, 4.0, "column\nmodelled\nin B", fontsize=max(8*k, 6.5), color=INK, va="center")

    # ---------------- B: model ----------------
    b.set_xlim(0, 14); b.set_ylim(0, 10); b.axis("off")
    b.text(0, 9.7, title_b, fontsize=max(11*k, 8.5), fontweight="bold", va="top")
    b.add_patch(FancyBboxPatch((1.0, 8.0), 3.0, 1.0, boxstyle="round,pad=0.05", fc=ICC, ec=INK, lw=0.8))
    b.text(2.5, 8.5, "ICC-MY (18 pF)", ha="center", va="center", fontsize=max(10*k, 7), color="white", fontweight="bold")
    b.text(4.3, 8.75, "g$_{own}$ ≈ 3.7 nS", fontsize=max(8*k, 6.5), va="center")
    b.text(4.3, 8.25, "18 min$^{-1}$", fontsize=max(8*k, 6.5), va="center", color=MUTED)
    b.plot([2.5, 2.5], [8.0, 7.05], color=INK, lw=1.5, ls=(0, (2, 1.5)))
    b.text(2.7, 7.5, "g$_c$", fontsize=max(9*k, 7), va="center", fontweight="bold")
    rows = [(6.6, "layer 1"), (5.5, "layer 2"), (4.4, "layer 3"), (1.6, "layer 30")]
    for y, lab in rows:
        b.add_patch(Rectangle((1.0, y - 0.4), 3.0, 0.8, fc=SMC, ec=INK, lw=0.8))
        b.text(2.5, y, "m SMC", ha="center", va="center", fontsize=max(9*k, 6.5))
        b.add_patch(Rectangle((5.0, y - 0.3), 1.4, 0.6, fc=FLC, ec=INK, lw=0.8))
        b.text(5.7, y, "FLC", ha="center", va="center", fontsize=max(8*k, 6.5))
        b.plot([4.0, 5.0], [y, y], color=INK, lw=0.8)
        b.text(0.9, y, lab, ha="right", va="center", fontsize=max(8*k, 6.5), color=MUTED)
    for (y1, _), (y2, _) in zip(rows[:2], rows[1:3]):
        b.plot([2.5, 2.5], [y1 - 0.4, y2 + 0.4], color=INK, lw=1.2)
    b.plot([2.5] * 3, [2.8, 3.2, 3.6], "o", color=INK, ms=2.5*k)
    b.text(0.9, 1.0, "sealed end", fontsize=max(7.5*k, 6.5), color=MUTED, ha="right")
    b.text(6.6, 6.6, "SK channels", fontsize=max(8*k, 6.5), va="center")
    b.annotate("purinergic\ninput", xy=(6.4, 5.5), xytext=(7.0, 4.6), fontsize=max(8*k, 6.5), va="center",
               arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.0), color=RED)
    b.text(2.5, 0.75, "SMC + tonic K$^+$", fontsize=max(8*k, 6.5), va="center", ha="center", color=MUTED)
    # constraints
    b.add_patch(FancyBboxPatch((9.6, 5.9), 4.2, 3.5, boxstyle="round,pad=0.1", fc=BOX, ec="none"))
    b.text(9.75, 9.2, "Constraints", fontsize=max(9*k, 7), fontweight="bold", va="top")
    b.text(10.35, 8.45, "layer 1 / ICC-MY\namplitude ratio\n→ g$_c$/G$_{load}$", fontsize=max(8*k, 6.5), va="top")
    b.text(10.35, 6.85, "ICC-MY R$_{in}$\n(mouse): check", fontsize=max(8*k, 6.5), va="top")
    b.text(9.75, 3.0, "G$_{load}$: load per\nICC-MY (1.61 nS)", fontsize=max(8*k, 6.5), va="top", color=MUTED)
    def tag(x, y, n):
        b.add_patch(plt.Circle((x, y), 0.22, fc=INK, ec="none"))
        b.text(x, y, n, ha="center", va="center", fontsize=max(7.5*k, 6.5), color="white", fontweight="bold")
    # where each constraint is read: 1 = ICC-MY and layer-1 amplitudes, 2 = current injected into ICC-MY
    tag(0.7, 8.75, "1"); tag(0.7, 8.25, "2"); tag(1.0, 7.2, "1")
    tag(9.95, 8.3, "1"); tag(9.95, 6.7, "2")



if __name__ == "__main__":
    fig = plt.figure(figsize=(12, 5.6))
    draw(fig.add_axes([0.02, 0.05, 0.40, 0.88]), fig.add_axes([0.45, 0.05, 0.54, 0.88]))
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)
