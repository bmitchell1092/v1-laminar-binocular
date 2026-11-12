"""Figure 1: the experiment, an example, and response onset follows drive (docs/09; rebuilt 2026-11-02 for readers
outside the field, docs/11; 2026-11-11: the onset-vs-drive panel of the simplified Fig 2 is Fig 1's final panel, and the
contrast responses and ocular dominance moved to the supplement, figS_monocular.py).

Usage: .venv\\Scripts\\python scripts\\figures\\fig1_setup.py [example penetration, default 221212_J_eA]
(a) the setup (schematics): a mirror stereoscope (top view; two mirrors per eye shift each eye's line of sight onto
    its own half of one monitor, divided by a black septum), so each eye sees its own grating at its own contrast; a
    32-contact laminar probe spans the depth of V1; below a, one trial (three conditions of ~317 ms with ~317 ms
    between them; restored under a without moving it, Blake's V6 comment);
(b) the 6 x 6 dichoptic contrast design (each eye 0, 5.5, 11, 22.5, 45 or 90%), with a grating drawn at each
    contrast (Michelson contrast on mid grey) so the contrasts can be seen; each channel's axes are later reoriented
    to its dominant and other eye;
(c) example penetration, channel mean with each channel aligned to its own response onset (half-max, dominant eye at
    90%), oriented to its dominant eye and divided by its dominant-eye peak at 90%: both eyes at 45% each, each eye
    alone at 45%, and what linear summation of contrast predicts, the dominant eye's response at 45% x (1 + a) with
    the penetration's median eye weight a (v1lb.binocular.response_at); the figure names the monkey, not the session;
(d) response onset against response strength for all 35 conditions of the design (check H, v4_drive_latency.py):
    half-max onset of each condition's own peak and that peak, per channel; medians over channels per penetration,
    then over penetrations; onsets relative to each channel's dominant eye at 90% plus that onset's population median;
    grey band, the least-squares line of onset on log peak through the monocular conditions, over their range; the
    statistics (checks H, G and U) are in the text and caption only (Blake's V6 comment: legend only in d).
crfs() (monocular contrast responses) and ocular_dominance() are drawn in Supplementary figS_monocular.py.
The example flash CSD with the laminar alignment is in Supplementary Fig S1 (figS1_methods.py).
Inputs: figure_data.py (population.npz, population_meta.parquet), run_crf.py, run_dichoptic.py.
Writes <out_root>/figures/fig1_setup.pdf and .png; prints the statistics shown.
"""

import glob
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
H_MM = 116
POS = {"stereo": (2, 66, 40, 38), "probe": (44, 66, 22, 38),          # (x, y, w, h) in mm; a top-aligned with b
       "matrix": (97, 69, 35, 35), "key": (137, 69, 37, 35),
       "example": (14, 10, 64, 36), "drive": (108, 10, 62, 36),
       "trial": (2, 54, 64, 10)}                  # trial timeline under a, a not moved (Blake, V6 comment)
QC = REPO / "docs" / "qc"
SWATCH_CMAP = "gray"


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def grating(c: float, n: int = 64, cycles: float = 2.5) -> np.ndarray:
    """Vertical sinusoidal grating of Michelson contrast c on mid grey (values 0-1)."""
    x = np.linspace(0, 2 * np.pi * cycles, n)
    return np.tile(0.5 + 0.5 * c * np.sin(x), (n, 1))


def swatch(ax, c, x0, x1, y0, y1):
    ax.imshow(grating(c), cmap=SWATCH_CMAP, vmin=0, vmax=1, extent=(x0, x1, y0, y1), interpolation="bilinear",
              clip_on=False, zorder=3)
    ax.add_patch(mpatches.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=fs.RULE, lw=0.4, clip_on=False,
                                    zorder=4))


def stereoscope(ax):
    """Top view: two mirrors per eye shift each line of sight onto its own half of one monitor; each half shows
    that eye's grating."""
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.25, 2.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(mpatches.Rectangle((-0.95, 1.2), 1.9, 0.07, color=fs.INK2, lw=0))
    ax.text(0, 1.31, "monitor", ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    ax.plot([0, 0], [1.2, 0.5], color=fs.INK, lw=2.2, solid_capstyle="butt")                 # septum
    for s, c, lab in ((-1, 0.45, "left eye: 45%"), (1, 0.11, "right eye: 11%")):
        ex, mx = s * 0.17, s * 0.6
        sight = dict(color=fs.MUTED, lw=0.6, ls=(0, (2, 1.5)))
        ax.plot([ex, ex], [0.1, 0.35], **sight)
        ax.plot([ex, mx], [0.35, 0.35], **sight)
        ax.plot([mx, mx], [0.35, 1.2], **sight)
        for x in (ex, mx):                                                                     # 45-degree mirrors
            ax.plot([x - s * 0.09, x + s * 0.09], [0.26, 0.44], color=fs.INK, lw=1.6, solid_capstyle="butt")
        ax.add_patch(mpatches.Circle((ex, 0.02), 0.08, fc="white", ec=fs.INK, lw=0.8))
        swatch(ax, c, mx - 0.24, mx + 0.24, 1.42, 1.9)
        ax.text(mx, 1.94, lab, ha="center", va="bottom", fontsize=5.2, color=fs.INK)
    ax.text(0, -0.1, "eyes", ha="center", va="top", fontsize=5.5, color=fs.INK2)
    ax.text(0.62, 0.2, "mirrors", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def probe(ax):
    """Side view of V1 with the laminar probe: 32 contacts, 0.1 mm apart, spanning the cortical depth."""
    ax.set_xlim(0, 1)
    ax.set_ylim(3.6, -0.75)
    ax.axis("off")
    ax.add_patch(mpatches.Rectangle((0.05, 0), 0.9, 2.4, fc="#f0efec", ec="none"))           # grey matter
    ax.plot([0.05, 0.95], [0, 0], color=fs.INK2, lw=0.6)
    ax.text(0.95, -0.05, "top of cortex", ha="right", va="bottom", fontsize=5.2, color=fs.INK2)
    ax.text(0.95, 2.47, "white matter", ha="right", va="top", fontsize=5.2, color=fs.INK2)
    ax.text(0.92, 1.2, "V1", ha="right", va="center", fontsize=6.5, color=fs.INK)
    x0 = 0.32
    ax.add_patch(mpatches.Rectangle((x0 - 0.035, -0.6), 0.07, 0.6 + 3.25, fc=fs.INK2, ec="none"))
    z = -0.2 + 0.1 * np.arange(32)
    ax.scatter(np.full(32, x0), z, s=1.2, color="white", zorder=3, lw=0)
    ax.annotate("32 contacts,\n0.1 mm apart", xy=(x0 + 0.04, 1.6), xytext=(0.42, 1.95), fontsize=5.2, color=fs.INK,
                va="top", linespacing=1.0, arrowprops=dict(arrowstyle="-", lw=0.4, color=fs.INK2))


def trial(ax):
    """One trial: fixation, then three conditions of ~317 ms each with ~317 ms between them."""
    ax.set_xlim(-450, 1950)
    ax.set_ylim(-0.6, 1.25)
    ax.axis("off")
    ax.plot([-400, 1900], [0.3, 0.3], color=fs.MUTED, lw=0.6)
    ax.text(-400, 0.45, "fixate", fontsize=5.5, color=fs.INK2, va="bottom")
    for k in range(3):
        x0 = k * 634
        ax.add_patch(mpatches.Rectangle((x0, 0.3), 317, 0.42, fc=fs.INK, ec="none"))
        ax.text(x0 + 158, 0.8, f"condition {k + 1}", ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    ax.plot([1583, 1900], [0.05, 0.05], color=fs.INK, lw=0.8)
    ax.text(1741, -0.08, "317 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def design(ax, kax):
    cls = np.zeros((6, 6), int)
    for i in range(6):
        for j in range(6):
            cls[i, j] = 0 if i == j == 0 else 1 if (i == 0 or j == 0) else 2 if i == j else 3
    cols = {0: "white", 1: fs.INK2, 2: fs.INK, 3: fs.RULE}
    for i in range(6):
        for j in range(6):
            ax.add_patch(mpatches.Rectangle((j + 0.07, i + 0.07), 0.86, 0.86, fc=cols[cls[i, j]],
                                            ec=fs.MUTED if cls[i, j] == 0 else "none", lw=0.5))
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    for k, c in enumerate(CONTRASTS):                    # each contrast as it looked
        swatch(ax, c, k + 0.14, k + 0.86, -0.92, -0.2)
        swatch(ax, c, -0.92, -0.2, k + 0.14, k + 0.86)
    for k, c in enumerate(CONTRASTS):                    # numbers beyond the swatches
        lab = "0" if c == 0 else f"{100 * c:g}"
        ax.text(k + 0.5, -1.02, lab, ha="center", va="top", fontsize=5.5, color=fs.INK, clip_on=False)
        ax.text(-1.02, k + 0.5, lab, ha="right", va="center", fontsize=5.5, color=fs.INK, clip_on=False)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("contralateral eye contrast (%)", labelpad=25)
    ax.set_ylabel("ipsilateral eye contrast (%)", labelpad=31)

    kax.set_xlim(0, 1)
    kax.set_ylim(0, 1)
    kax.axis("off")
    for k, (lab, y) in enumerate((("one eye", 0.95), ("both eyes, equal", 0.85), ("both eyes, unequal", 0.75))):
        kax.add_patch(mpatches.Rectangle((0.0, y - 0.035), 0.1, 0.07, fc=cols[k + 1], ec="none"))
        kax.text(0.14, y, lab, va="center", fontsize=5.5, color=fs.INK)
    kax.text(0.0, 0.62, "~30 presentations per\ncondition, random order", va="top", fontsize=5.5, color=fs.INK2,
             linespacing=1.0)
    kax.text(0.0, 0.38, "analysed per channel as\ndominant eye (the eye that\ndrives it more) × other eye",
             va="top", fontsize=5.5, color=fs.INK2, linespacing=1.0)


def example_raw(ax, pop, k, monkey):
    """Example penetration on the stimulus clock: both eyes and each eye alone at 45% (simplified set)."""
    t = pop["t"]
    Y = pop["flash"][k]
    i45 = list(CONTRASTS).index(0.45)
    sel = (t >= -20) & (t <= 340)
    items = [("both eyes, 45% each", fs.COND["both"], Y[i45, i45]),
             ("dominant eye alone, 45%", fs.COND["dominant"], Y[i45, 0]),
             ("other eye alone, 45%", fs.COND["other"], Y[0, i45])]
    for lab, col, y in items[::-1]:
        ax.plot(t[sel], y[sel], color=col, lw=1.2 if col == fs.INK else 1.0)
    fs.zero_line(ax)
    top = 1.1 * np.nanmax(Y[i45, i45][sel])
    ax.set_xlim(-20, 400)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_ylim(-0.1, top * 1.02)
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("MUA (dominant-eye\npeak = 1)")
    fs.line_key(ax, [(lab, col, "-") for lab, col, _ in items], x=0.5, y=0.8, dy=0.12, length=0.07)
    ax.set_title(f"example penetration, monkey {monkey}", fontsize=6.5)
    return {lab: float(np.nanmax(y[sel])) for lab, _, y in items}


def example(ax, pop, k, monkey, a):
    rel = pop["rel"]
    Y = pop["onset"][k]
    i45 = list(CONTRASTS).index(0.45)
    sel = (rel >= -20) & (rel <= 200)
    items = [("both eyes, 45% each", fs.COND["both"], "-", Y[i45, i45]),
             ("dominant eye alone, 45%", fs.COND["dominant"], "-", Y[i45, 0]),
             ("other eye alone, 45%", fs.COND["other"], "-", Y[0, i45]),
             (f"{fs.SUM}\n(45% + weight × 45%)", fs.COND["linear"], fs.SUM_LS, response_at(Y[:, 0], 0.45 * (1 + a)))]
    for lab, col, ls, y in items[::-1]:
        ax.plot(rel[sel], y[sel], color=col, ls=ls, lw=1.2 if col == fs.INK else 1.0)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye\npeak = 1)")
    top = 1.08 * np.nanmax(Y[5, 0][sel])
    ax.set_ylim(-0.1, top)
    fs.line_key(ax, [(lab, col, ls) for lab, col, ls, _ in items], x=0.42, y=0.95, dy=0.125, length=0.07)
    ax.set_title(f"example penetration, monkey {monkey}", fontsize=6.5)
    ax.text(23, top * 0.97, "first 20 ms", ha="left", va="top", fontsize=5.5, color=fs.INK2)
    e = (rel >= 0) & (rel < 20)
    late = (rel >= 50) & (rel < 200)
    return {w: {lab.split(" (")[0]: float(y[m].mean()) for lab, _, _, y in items}
            for w, m in (("0-20 ms", e), ("50-200 ms", late))}


def crfs(ax, pop):
    t = pop["t"]
    w = (t >= 40) & (t < 250)
    F = pop["flash"][:, :, :, w].mean(axis=3)
    c = CONTRASTS[1:]
    out = {}
    for y, col, lab in ((F[:, 1:, 0], fs.COND["dominant"], "dominant\neye"), (F[:, 0, 1:], fs.COND["other"], "other\neye")):
        m, lo, hi = fs.boot_ci(y)
        ax.fill_between(c, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(c, m, color=col, marker="o", ms=2.8, lw=1.1)
        fs.direct_label(ax, c[-1], m[-1], lab, col, dx=4, linespacing=0.9)
        out[lab.replace("\n", " ")] = (m, lo, hi)
    fs.log_contrast_axis(ax, c)
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"])
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel("MUA (dominant-eye\npeak = 1)")
    ax.set_ylim(0, None)
    ax.set_xlim(0.045, 1.0)
    ax.set_title("each eye alone, 40–250 ms", fontsize=6.5)
    ax.text(0.05, 0.98 * ax.get_ylim()[1], f"mean, 95% CI,\n{F.shape[0]} penetrations", fontsize=5.5,
            color=fs.INK2, va="top", linespacing=1.0)
    ratio = F[:, 0, 5] / F[:, 5, 0]
    return out, ratio


def ocular_dominance(ax):
    D = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "dichoptic" / "*.parquet"))], ignore_index=True)
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[R.condition.str.startswith("mono") & R.driven].groupby(["penetration", "site", "site_id"]).size()
    D = D.join(drv.rename("nd"), on=["penetration", "site", "site_id"])
    D = D[D.nd.notna() & (D.window == "full")]
    isolated = (D.isi_contamination < 0.5) & (D.snr >= 3)                     # lenient isolation (05 C14)
    D = D[(D.site == "mua") | ((D.site == "unit") & isolated)]
    odi = {"mua": [], "unit": []}
    for _, g in D.groupby(["penetration", "site", "site_id"]):
        mm = g.set_index(["contrast_contra", "contrast_ipsi"])["mean"]
        try:
            rc, ri = mm[(0.9, 0.0)], mm[(0.0, 0.9)]
        except KeyError:
            continue
        if max(rc, ri) <= 0:
            continue
        odi[g.site.iloc[0]].append((rc - ri) / (abs(rc) + abs(ri)))
    bins = np.linspace(-1, 1, 21)
    hm = np.histogram(odi["mua"], bins)[0] / len(odi["mua"])
    hu = np.histogram(odi["unit"], bins)[0] / len(odi["unit"])
    ax.bar(bins[:-1], hm, width=np.diff(bins), align="edge", color=fs.RULE, lw=0)
    ax.stairs(hu, bins, color=fs.INK, lw=0.9)
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, 0, 1], ["−1", "0", "1"])
    ax.set_xlabel("ocular dominance index", labelpad=1)
    ax.set_ylabel("fraction")
    ax.set_title("most sites respond to both eyes", fontsize=6.5)
    top = max(hm.max(), hu.max()) * 1.3
    ax.set_ylim(0, top)
    ax.text(-0.97, top * 0.62, "ipsilateral\neye only", fontsize=5.5, color=fs.INK2, ha="left", va="top",
            linespacing=0.95)
    ax.text(0.97, top * 0.62, "contralateral\neye only", fontsize=5.5, color=fs.INK2, ha="right", va="top",
            linespacing=0.95)
    ax.add_patch(mpatches.Rectangle((-0.95, top * 0.9), 0.14, top * 0.06, fc=fs.RULE, ec="none"))
    ax.text(-0.76, top * 0.93, f"MUA, n = {len(odi['mua'])}", fontsize=5.5, va="center", color=fs.INK)
    ax.plot([-0.95, -0.81], [top * 0.8] * 2, color=fs.INK, lw=0.9)
    ax.text(-0.76, top * 0.8, f"units, n = {len(odi['unit'])}", fontsize=5.5, va="center", color=fs.INK)
    return {k: np.asarray(v) for k, v in odi.items()}


def drive_latency(ax, show_added=True, notes=True):
    """Response onset against response strength (peak) for all 35 conditions (check H, v4_drive_latency.py)."""
    X = pd.read_csv(REPO / "docs" / "qc" / "v4_H_cells.csv")
    C = pd.read_parquet(out_root() / "binocular" / "v4_H_cells.parquet")
    A = C[(C.split == "all") & (C.c_dom == 0.9) & (C.c_oth == 0)]
    base = float(A.groupby("penetration").onset.median().median())       # dominant eye at 90%, absolute onset
    X["y"], X["y_add"] = X.onset_rel_median + base, X.onset_added_rel_median + base
    dom = X[X.c_oth == 0].sort_values("c_dom")
    oth = X[X.c_dom == 0].sort_values("c_oth")
    both = X[(X.c_dom > 0) & (X.c_oth > 0)]
    dio = both[np.isclose(both.c_dom, both.c_oth)]
    mono = pd.concat([dom, oth])
    slope, icpt = np.polyfit(np.log(mono.peak_median), mono.y, 1)
    xs = np.linspace(mono.peak_median.min(), mono.peak_median.max(), 50)     # only over the monocular range
    ax.plot(xs, icpt + slope * np.log(xs), color=fs.RULE, lw=2.4, zorder=1, solid_capstyle="round")
    if show_added:
        for _, r in dio.iterrows():
            ax.plot([r.peak_median, r.peak_added_median], [r.y, r.y_add], color=fs.RULE, lw=0.5,
                    ls=(0, (1.5, 1.5)), zorder=2)
        ax.scatter(both.peak_added_median, both.y_add, s=6, marker="D", facecolor="white", edgecolor="#b5b4af",
                   lw=0.5, zorder=3)
    ax.scatter(both.peak_median, both.y, s=9, marker="D", color=fs.INK, lw=0, zorder=4)
    ax.scatter(dio.peak_median, dio.y, s=22, marker="D", color=fs.INK, lw=0, zorder=5)
    ax.plot(dom.peak_median, dom.y, color=fs.COND["dominant"], lw=0.8, marker="o", ms=3.2, zorder=6)
    ax.plot(oth.peak_median, oth.y, color=fs.COND["other"], lw=0.8, marker="o", ms=3.2, zorder=6)
    for _, r in dom.iterrows():
        up = r.c_dom >= 0.45                                       # the both-eyes points crowd below these two
        ax.text(r.peak_median - 0.025, r.y + (0.35 if up else -0.25), f"{100 * r.c_dom:g}", fontsize=5, ha="right",
                va="bottom" if up else "top", color=fs.COND["dominant"])
    ax.set_xlim(0.2, 1.9 if show_added else 1.3)
    ax.set_yticks(np.arange(30, 46, 5))
    ax.set_xlabel("response strength (peak; dominant eye at 90% = 1)")
    ax.set_ylabel("response onset\n(ms after stimulus onset)")
    ax.text(0.26, icpt + slope * np.log(0.26) + 0.5, "stronger drive, earlier start", fontsize=5.5, color=fs.INK2,
            va="bottom")
    keys = [("o", fs.COND["dominant"], fs.COND["dominant"], "dominant eye alone (contrast, %)"),
            ("o", fs.COND["other"], fs.COND["other"], "other eye alone"),
            ("D", fs.INK, fs.INK, "both eyes (large: equal contrasts)")]
    if show_added:
        keys.append(("D", "white", "#b5b4af", "if the two eyes' responses simply added"))
    kx = 1.05 if notes else 0.42
    for k, (mk, fc, ec, lab) in enumerate(keys):
        y = 0.95 - k * 0.09
        ax.plot([kx], [y], marker=mk, ms=3.0, mfc=fc, mec=ec, mew=0.6, transform=ax.transAxes, clip_on=False)
        ax.text(kx + 0.03, y, lab, transform=ax.transAxes, va="center", fontsize=5.5, color=fs.INK)
    H = pd.read_csv(REPO / "docs" / "qc" / "v4_H.csv").set_index("test")
    rho = stats.spearmanr(X.peak_median, X.y).statistic
    h1 = H.loc["vs_added"]
    h3 = H.loc["rho"]
    if not notes:
        return dict(base=base, rho=rho, slope=slope, H=H)
    ax.text(1.05, 0.38, f"35 conditions, 26 penetrations\nrank correlation over conditions {rho:.2f}\n"
            f"per penetration: median {h3['median']:.2f}, {int(h3.n_in_direction)}/{int(h3.n)} negative\n"
            f"both eyes start {abs(h1['median']):.1f} ms earlier than the sum\nof the two monocular responses would "
            f"({int(h1.n_in_direction)}/{int(h1.n)}, {fs.fmt_p(h1.p)})",
            transform=ax.transAxes, va="top", fontsize=5.5, color=fs.INK2, linespacing=1.15)
    return dict(base=base, rho=rho, slope=slope, H=H)


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    k = list(pop["pens"]).index(pen)

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")                     # axes placed in mm: the in-panel keys stay inside their panels
    sax, pax, dax, kax = (axes_mm(fig, *POS[n]) for n in ("stereo", "probe", "matrix", "key"))
    stereoscope(sax)
    probe(pax)
    trial(axes_mm(fig, *POS["trial"]))
    design(dax, kax)
    bax, cax = (axes_mm(fig, *POS[n]) for n in ("example", "drive"))
    ex = example_raw(bax, pop, k, meta.loc[pen, "monkey"])
    dl = drive_latency(cax, show_added=False, notes=False)
    U = pd.read_csv(QC / "v4_U.csv").set_index("test")
    G = pd.read_csv(QC / "v4_G.csv").set_index("test")
    u4, g45 = U.loc["U4 rank correlation of onset and peak < 0"], G.loc["G1 both eyes earlier, 0.45"]
    fs.panel_label(sax, "a", x=0.04)
    for ax, letter in zip((dax, bax, cax), "bcd"):
        fs.panel_label(ax, letter, x=-0.18 if ax is dax else -0.02 if ax is bax else -0.25)
    print(fs.save(fig, "fig1_setup"))
    print(f"(d) rank correlation over the 35 conditions {dl['rho']:.3f}; hierarchical per penetration "
          f"{u4.estimate:.3f} (p {u4.p_one_sided:.2g}); both eyes - dominant at 45% {g45['median']:.2f} ms")
    print(f"(c) example: monkey {meta.loc[pen, 'monkey']}, {int(meta.loc[pen, 'n_channels'])} channels, "
          f"median eye weight {meta.loc[pen, 'a']:.2f}, onset {meta.loc[pen, 'onset']:.1f} ms after stimulus onset; "
          "peaks " + "; ".join(f"{lab} {v:.2f}" for lab, v in ex.items()))
    print(f"    median eye weight a over penetrations {meta.a.median():.2f} (IQR {meta.a.quantile(.25):.2f}-"
          f"{meta.a.quantile(.75):.2f}); {len(meta)} penetrations, {int(meta.n_channels.sum())} channels "
          f"(B {int((meta.monkey == 'B').sum())}, J {int((meta.monkey == 'J').sum())})")


if __name__ == "__main__":
    main()
