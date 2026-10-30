"""Figure 1: binocular combination measured with millisecond and laminar precision (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig1_setup.py [example penetration, default 221212_J_eA]
(a) mirror stereoscope (schematic, top view): two mirrors per eye shift each eye's line of sight onto its own half of
    one monitor, divided by a black septum; (b) the 6 x 6 dichoptic contrast design and the trial sequence;
(c) example penetration: flash-evoked CSD (bad contacts interpolated, Hamming-smoothed), registered to the grating
    recording by the alignment's depth shift, with the top of cortex and the L4/5 boundary below the flash sink
    (laminar alignment v2);
(d) example penetration, channel mean (each channel oriented to its dominant eye, divided by its dominant-eye peak
    at .9): each eye alone at .45, both eyes at .45, and the dominant eye at .9, about what linear summation of
    contrast predicts ((1 + a) c with a ~ 0.9);
(e) monocular contrast response per eye, 40-250 ms, penetration means (95% bootstrap CI over penetrations);
(f) ocular dominance index at .9 for driven MUA channels and isolated single units.
Photodiode timing (onset jitter by whole frames) is no longer a panel (standard practice; docs/01).
Inputs: figure_data.py, datasets, flash CSD, laminar alignment v2, run_crf.py, run_dichoptic.py.
Writes <out_root>/figures/fig1_setup.pdf and .png.
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

from v1lb import figstyle as fs  # noqa: E402
from v1lb import laminar as lam  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import metadata_dir, out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]


def stereoscope(ax):
    """Top-down schematic: two mirrors per eye shift each line of sight outward onto its own half of one monitor,
    and a black septum divides the monitor between the eyes."""
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.25, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(mpatches.Rectangle((-0.95, 1.2), 1.9, 0.07, color=fs.INK2, lw=0))
    ax.text(0, 1.31, "monitor", ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    ax.plot([0, 0], [1.2, 0.5], color=fs.INK, lw=2.2, solid_capstyle="butt")                 # septum
    ax.text(0.04, 0.82, "septum", ha="left", va="center", fontsize=5.5, color=fs.INK2, rotation=90)
    for s in (-1, 1):
        ex, mx = s * 0.17, s * 0.6
        sight = dict(color=fs.MUTED, lw=0.6, ls=(0, (2, 1.5)))
        ax.plot([ex, ex], [0.1, 0.35], **sight)
        ax.plot([ex, mx], [0.35, 0.35], **sight)
        ax.plot([mx, mx], [0.35, 1.2], **sight)
        for x in (ex, mx):                                                                     # 45-degree mirrors
            ax.plot([x - s * 0.09, x + s * 0.09], [0.26, 0.44], color=fs.INK, lw=1.6, solid_capstyle="butt")
        ax.add_patch(mpatches.Circle((ex, 0.02), 0.08, fc="white", ec=fs.INK, lw=0.8))
        ax.text(mx, 1.1, "left half" if s < 0 else "right half", ha="center", va="top", fontsize=5,
                color=fs.INK2)
    ax.text(0, -0.1, "eyes", ha="center", va="top", fontsize=5.5, color=fs.INK2)
    ax.text(0.6, 0.18, "mirrors", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def design(ax, tax):
    cls = np.zeros((6, 6), int)
    for i in range(6):
        for j in range(6):
            cls[i, j] = 0 if i == j == 0 else 1 if (i == 0 or j == 0) else 2 if i == j else 3
    cols = {0: "white", 1: fs.INK2, 2: fs.INK, 3: fs.RULE}
    for i in range(6):
        for j in range(6):
            ax.add_patch(mpatches.Rectangle((j + 0.06, i + 0.06), 0.88, 0.88, fc=cols[cls[i, j]],
                                            ec=fs.MUTED if cls[i, j] == 0 else "none", lw=0.5))
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ticks = np.arange(6) + 0.5
    labs = ["0"] + [f"{100 * c:g}" for c in CONTRASTS[1:]]
    ax.set_xticks(ticks, labs)
    ax.set_yticks(ticks, labs)
    ax.tick_params(length=0, pad=1.5, labelsize=5.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("contralateral eye contrast (%)", labelpad=1.5)
    ax.set_ylabel("ipsilateral eye contrast (%)", labelpad=1.5)
    for k, (lab, y) in enumerate((("one eye", 5.3), ("both eyes, equal", 4.2), ("both eyes, unequal", 3.1))):
        ax.add_patch(mpatches.Rectangle((6.4, y - 0.3), 0.6, 0.6, fc=cols[k + 1], ec="none", clip_on=False))
        ax.text(7.2, y, lab, va="center", fontsize=6, clip_on=False)
    ax.text(6.4, 1.6, "~30 presentations\nper cell, random\norder", va="center", fontsize=6, color=fs.INK2,
            clip_on=False, linespacing=1.0)

    # trial: fixation, then three random cells (~317 ms on, ~317 ms off)
    tax.set_xlim(-450, 1950)
    tax.set_ylim(-0.4, 1.2)
    tax.axis("off")
    tax.plot([-400, 1900], [0.3, 0.3], color=fs.MUTED, lw=0.6)
    tax.text(-400, 0.45, "fixate", fontsize=5.5, color=fs.INK2, va="bottom")
    for k in range(3):
        x0 = k * 634
        tax.add_patch(mpatches.Rectangle((x0, 0.3), 317, 0.42, fc=fs.INK, ec="none"))
        tax.text(x0 + 158, 0.8, f"cell {k + 1}", ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    tax.plot([1583, 1900], [0.05, 0.05], color=fs.INK, lw=0.8)
    tax.text(1741, -0.05, "317 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def flash_csd(ax, cax, pen, al):
    # the flash recording is registered to the grating recording (laminar alignment v2): contact + shift
    reg = pd.read_csv(REPO / "docs" / "qc" / "flash_landmarks_v2.csv")
    reg = reg[(reg.penetration == pen) & reg.usable].iloc[0]
    z = np.load(out_root() / "flash_csd" / f"{pen}_{reg.evp}.npz")
    lfp, t, con = z["lfp"], z["t_ms"], z["contacts"]
    bad = lam.bad_channels(lfp[(t >= 20) & (t <= 150)])
    csd = lam.csd(lfp, smooth=True, bad=bad)
    sgn = 1 if al.orientation == "ascending" else -1
    order = np.argsort(sgn * con)
    C = csd[:, order].T
    C -= np.nanmean(C[:, (t >= -50) & (t < 0)], axis=1, keepdims=True)
    depth = sgn * (con[order] + reg["shift"] - al.top_of_cortex) * 0.1     # mm below the top of cortex
    ok = np.isfinite(C).all(axis=1)
    lim = np.nanpercentile(np.abs(C[ok][:, (t > 20) & (t < 150)]), 98)
    im = ax.pcolormesh(t, depth[ok], C[ok], cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, shading="nearest",
                       rasterized=True)
    ax.set_xlim(-10, 120)
    ax.set_ylim(min(depth[ok].max() + 0.05, 2.8), min(depth[ok].min() - 0.05, -0.15))
    l45 = (al.l45 - al.top_of_cortex) * 0.1
    for y, lab in ((-0.05, "top of cortex"), (l45, "L4/5")):
        ax.axhline(y, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
        ax.text(117, y - 0.05, lab, ha="right", va="bottom", fontsize=6, color=fs.INK,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
    ax.set_xlabel("time from flash (ms)")
    ax.set_ylabel("depth below top of cortex (mm)")
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)


def photodiode(ax):
    d = out_root() / "datasets"
    lat = pd.concat([pd.read_parquet(f).pd_on_lat_ms for f in glob.glob(str(d / "*_presentations.parquet"))])
    lat = lat[np.isfinite(lat)]
    ax.hist(lat, bins=np.arange(5, 40, 0.5), color=fs.INK2, lw=0)
    ax.set_yscale("log")
    ax.set_ylim(0.8, 2e5)
    ax.set_xlim(5, 40)
    ax.set_xticks([10, 20, 30, 40])
    ax.set_xlabel("photodiode onset − event code (ms)")
    ax.set_ylabel("presentations")
    ax.annotate("", xy=(21, 4e4), xytext=(29.5, 4e4), arrowprops=dict(arrowstyle="<->", lw=0.6, color=fs.INK,
                                                                       shrinkA=0, shrinkB=0))
    ax.text(30.5, 4e4, "1 frame\n(8.3 ms)", ha="left", va="center", fontsize=5.5, color=fs.INK, linespacing=0.9)
    ax.text(6, 1.2e5, f"n = {len(lat):,}", ha="left", va="top", fontsize=6, color=fs.INK2)


def example(ax, pop, k):
    t = pop["t"]
    Y = pop["flash"][k]
    i45, i9 = list(CONTRASTS).index(0.45), 5
    traces = ((Y[i9, 0], fs.COND["linear"], (0, (3, 1.5)), "dominant eye, 90%\n(≈ linear sum)"),
              (Y[i45, 0], fs.COND["dominant"], "-", "dominant eye"),
              (Y[0, i45], fs.COND["other"], "-", "other eye"),
              (Y[i45, i45], fs.COND["both"], "-", "both eyes"))
    for y, col, ls, lab in traces:
        ax.plot(t, y, color=col, ls=ls, lw=1.0 if col != fs.INK else 1.2)
    ax.axhline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_xlim(-20, 250)
    ax.set_xticks([0, 100, 200])
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    top = 1.05 * np.nanmax(Y[i9, 0])
    ax.set_ylim(-0.1, top)
    ax.text(245, top, "45% contrast\none penetration", ha="right", va="top", fontsize=5.5, color=fs.INK2,
            linespacing=1.0)
    for (y, col, ls, lab), yy in zip(traces, (0.76, 0.64, 0.55, 0.46)):
        ax.plot([108, 126], [yy * top] * 2, color=col, ls=ls, lw=1.0)
        ax.text(130, yy * top, lab, va="center", fontsize=5.5, color=fs.INK, linespacing=0.9)


def crfs(ax, pop):
    t = pop["t"]
    w = (t >= 40) & (t < 250)
    F = pop["flash"][:, :, :, w].mean(axis=3)
    c = CONTRASTS[1:]
    for y, col, lab in ((F[:, 1:, 0], fs.COND["dominant"], "dominant\neye"), (F[:, 0, 1:], fs.COND["other"], "other\neye")):
        m, lo, hi = fs.boot_ci(y)
        ax.fill_between(c, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(c, m, color=col, marker="o", ms=2.8, lw=1.1)
        fs.direct_label(ax, c[-1], m[-1], lab, col, dx=4, linespacing=0.9)
    fs.log_contrast_axis(ax, c)
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"])
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel("MUA, 40–250 ms")
    ax.set_ylim(0, None)
    ax.set_xlim(0.045, 1.0)
    ax.text(0.05, 0.98 * ax.get_ylim()[1], f"{F.shape[0]} penetrations\nmean, 95% CI", fontsize=5.5,
            color=fs.INK2, va="top", linespacing=1.0)


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
    ax.set_xticks([-1, 0, 1], ["−1\nipsi", "0", "1\ncontra"])
    ax.set_xlabel("ocular dominance index", labelpad=1)
    ax.set_ylabel("fraction of sites")
    top = max(hm.max(), hu.max()) * 1.3
    ax.set_ylim(0, top)
    ax.add_patch(mpatches.Rectangle((-0.95, top * 0.9), 0.14, top * 0.06, fc=fs.RULE, ec="none"))
    ax.text(-0.76, top * 0.93, f"MUA, n = {len(odi['mua'])}", fontsize=5.5, va="center", color=fs.INK)
    ax.plot([-0.95, -0.81], [top * 0.8] * 2, color=fs.INK, lw=0.9)
    ax.text(-0.76, top * 0.8, f"units, n = {len(odi['unit'])}", fontsize=5.5, va="center", color=fs.INK)


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration").loc[pen]
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    k = list(pop["pens"]).index(pen)

    fig = fs.figure("double", 118)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.1, 1], hspace=0.03)
    t1, t2, t3 = top.subfigures(1, 3, width_ratios=[0.7, 1.45, 1.15], wspace=0.02)
    sax = t1.subplots()
    stereoscope(sax)
    dax, tax = t2.subplots(2, 1, gridspec_kw=dict(height_ratios=[4.2, 1]))
    design(dax, tax)
    cax_ = t3.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 0.05]))
    flash_csd(cax_[0], cax_[1], pen, al)
    ba = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    example(ba[0], pop, k)
    crfs(ba[1], pop)
    ocular_dominance(ba[2])
    for ax, letter in zip((sax, dax, cax_[0], ba[0], ba[1], ba[2]), "abcdef"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig1_setup"))


if __name__ == "__main__":
    main()
