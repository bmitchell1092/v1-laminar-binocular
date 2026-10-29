"""House style for the paper figures (docs/09). Every figure script calls `apply()` and uses these colours.

Colours (validated with the dataviz palette validator, light surface):
* contrast is ordinal -> one-hue blue ramp, light -> dark, lightest step clears 2:1 on white;
* conditions at one contrast: dominant eye blue, other eye orange, both eyes ink, linear-summation reference
  muted gray (dashed);
* models: data ink, fixed-delay (feedforward) orange, recurrent aqua (aqua is < 3:1 on white -> always
  direct-labelled);
* monkeys are encoded by marker shape, not colour;
* CSD: diverging blue (sink) <-> red (source) with a gray midpoint.
Sizes follow journal columns (single 85 mm, 1.5 column 114 mm, double 174 mm); text 7 pt; vector PDF with
editable text plus a 300 dpi PNG.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .config import out_root

MM = 1 / 25.4
WIDTH = {"single": 85 * MM, "onehalf": 114 * MM, "double": 174 * MM}

INK, INK2, MUTED, RULE = "#0b0b0b", "#52514e", "#8a8985", "#d9d8d3"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
CONTRAST = {0.055: "#86b6ef", 0.11: "#3987e5", 0.225: "#256abf", 0.45: "#184f95", 0.9: "#0d366b"}
COND = {"dominant": BLUE, "other": ORANGE, "both": INK, "linear": MUTED}
MODEL = {"data": INK, "fixed-delay": ORANGE, "recurrent": AQUA}
MONKEY_MARKER = {"B": "o", "J": "s"}
CSD_CMAP = LinearSegmentedColormap.from_list("csd", ["#0d366b", "#3987e5", "#f0efec", "#e34948", "#8f1f1e"])


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.titlepad": 4,
        "axes.linewidth": 0.6, "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.spines.top": False,
        "axes.spines.right": False, "axes.facecolor": "white", "figure.facecolor": "white",
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "lines.linewidth": 1.2, "lines.markersize": 3.5, "patch.linewidth": 0.6,
        "legend.frameon": False, "legend.handlelength": 1.6, "legend.borderaxespad": 0.2,
        "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def figure(width: str = "double", height_mm: float = 100.0):
    return plt.figure(figsize=(WIDTH[width], height_mm * MM), layout="constrained")


def panel_label(ax, letter: str, x: float = -0.02, y: float = 1.0):
    """Bold panel letter at the top-left corner, outside the axes (figure-relative offset in axes units)."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(x * 100 - 14, 6 + y - 1),
                textcoords="offset points", fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")


def direct_label(ax, x, y, text, color, dx: float = 3, **kw):
    ax.annotate(text, xy=(x, y), xytext=(dx, 0), textcoords="offset points", color=color, fontsize=6.5,
                va="center", ha="left", **kw)


def contrast_label(c: float) -> str:
    return f"{100 * c:g}%"


def log_contrast_axis(ax, contrasts=(0.055, 0.11, 0.225, 0.45, 0.9), label="contrast"):
    ax.set_xscale("log")
    ax.set_xticks(contrasts)
    ax.set_xticklabels([contrast_label(c) for c in contrasts])
    ax.minorticks_off()
    ax.set_xlabel(label)


def save(fig, name: str) -> Path:
    out = out_root() / "figures"
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        (out / f"{name}.{ext}").unlink(missing_ok=True)      # an image viewer's handle can block overwriting
    fig.savefig(out / f"{name}.pdf")
    fig.savefig(out / f"{name}.png", dpi=300)
    return out / f"{name}.png"


def boot_ci(m: np.ndarray, n: int = 2000, seed: int = 0, axis: int = 0, ci: float = 95, stat=np.nanmean):
    """Statistic (mean by default) over rows (penetrations) and its bootstrap CI."""
    rng = np.random.default_rng(seed)
    m = np.asarray(m, float)
    b = np.array([stat(m[rng.integers(0, len(m), len(m))], axis=axis) for _ in range(n)])
    lo, hi = np.nanpercentile(b, [(100 - ci) / 2, 100 - (100 - ci) / 2], axis=0)
    return stat(m, axis=axis), lo, hi


def ci_line(ax, x, Y, color, ls="-", lw=1.1, alpha=0.18, stat=np.nanmean, **kw):
    """Penetration statistic with its 95% bootstrap band; Y [penetration, x]."""
    m, lo, hi = boot_ci(Y, stat=stat)
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, lw=0)
    ax.plot(x, m, color=color, ls=ls, lw=lw, **kw)
    return m


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "p n/a"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def zero_line(ax, y=0.0, **kw):
    ax.axhline(y, color=kw.pop("color", RULE), lw=kw.pop("lw", 0.6), zorder=kw.pop("zorder", 0), **kw)


def line_key(ax, items, x=0.6, y=0.95, dy=0.085, length=0.08, fontsize=5.5):
    """Compact legend in axes fraction: items = [(label, color, linestyle[, marker[, marker face]])]."""
    for k, it in enumerate(items):
        lab, col, ls = it[:3]
        mk = it[3] if len(it) > 3 else None
        mfc = it[4] if len(it) > 4 else col
        yy = y - k * dy
        ax.plot([x, x + length], [yy, yy], color=col, ls=ls, lw=1.0, transform=ax.transAxes, clip_on=False)
        if mk:
            ax.plot([x + length / 2], [yy], color=col, marker=mk, ms=2.8, mfc=mfc, transform=ax.transAxes,
                    clip_on=False)
        ax.text(x + length + 0.03, yy, lab, transform=ax.transAxes, va="center", fontsize=fontsize, color=INK)


def stim_bar(ax, x0, x1, y, lw=1.5):
    """Thick ink bar marking the stimulus (data coordinates)."""
    ax.plot([x0, x1], [y, y], color=INK, lw=lw, solid_capstyle="butt", clip_on=False)
