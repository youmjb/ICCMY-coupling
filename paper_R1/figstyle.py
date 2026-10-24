"""Common style and output for the submission figures (revision 3 figure check).

KJPP: figures >= 300 dpi as PDF/JPEG/DOCX/PPTX (TIFF accepted for the graphical abstract), RGB, title in the legend only.
Common journal practice: sans-serif (Arial) throughout, also in math text; >= 6.5 pt at final size; double-column width 17.4 cm.
save(fig, path_without_ext) writes .pdf (TrueType fonts embedded), .tif (RGB, LZW, 300 dpi) and .png (RGB, 300 dpi, for Word).
"""
import matplotlib
from PIL import Image

WIDTH = 6.85                       # inches = 17.4 cm, double column
STYLE = {"font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
         "ytick.labelsize": 7, "legend.fontsize": 6.5, "legend.title_fontsize": 6.5, "axes.spines.top": False,
         "axes.spines.right": False, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
         "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
         "mathtext.sf": "Arial", "mathtext.default": "regular", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
SMALL = 6.5                        # smallest text size used in annotations


def apply():
    matplotlib.rcParams.update(STYLE)


def save(fig, stem, dpi=300):
    fig.savefig(stem + ".pdf")
    import io
    buf = io.BytesIO()                                   # render in memory and write each file once (Dropbox locks fresh files)
    fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
    buf.seek(0)
    im = Image.open(buf).convert("RGB")
    im.save(stem + ".png", dpi=(dpi, dpi))
    im.save(stem + ".tif", dpi=(dpi, dpi), compression="tiff_lzw")
