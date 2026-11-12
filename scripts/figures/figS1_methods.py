"""Supplementary Figure S1: how the linearity index is read, and the laminar alignment (docs/09, figure set for V3;
the stereoscope schematic moved to Fig 1a and the equivalent-contrast readout came here from Fig 2b on 2026-11-02).

Usage: .venv\\Scripts\\python scripts\\figures\\figS1_methods.py [example penetration, default 221212_J_eA]
(a) how the equivalent contrast is read (population means, illustration), each window on its own response axis: the
    binocular response at 45% is placed on the dominant eye's contrast response measured in the same window and read
    off as the dominant-eye contrast that gives the same response; the triangle marks what linear summation of
    contrast predicts, 45% x (1 + median eye weight). The linearity index is (equivalent contrast / 45% - 1) / eye
    weight: 1 when the contrasts add, 0 for the dominant eye alone;
(b) the laminar alignment (v4; 2026-11-12): monocular grating CSD averaged over monkey J's tier-1 penetrations,
    each aligned to its layer 4C/5 anchor and scaled by its squish factor (tier1_csd); the single-session flash CSD
    with the v2 landmarks (flash_csd) is kept for reference but no longer drawn.
Inputs: figure_data.py (population.npz, population_meta.parquet); <out_root>/flash_csd/*.npz,
docs/qc/flash_landmarks_v2.csv, metadata/laminar_alignment_v2.csv.
Writes <out_root>/figures/figS1_methods.pdf and .png; prints the readout and the landmarks shown.
"""

import sys
import warnings
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb import laminar as lam  # noqa: E402
from v1lb.binocular import CONTRASTS, invert_crf  # noqa: E402
from v1lb.config import metadata_dir, out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]


def equivalent_contrast_sketch(ax, pop, w, title, a_med):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= w[0]) & (rel < w[1])
    crf = O[:, :, 0, :][:, :, sel].mean(axis=2).mean(axis=0)              # dominant eye, population mean
    i = list(CONTRASTS).index(0.45)
    rb = O[:, i, i, :][:, sel].mean(axis=1).mean()
    kc, cens = invert_crf(crf, np.array([rb]))
    kc = float(kc[0]) if cens[0] == 0 else CONTRASTS[-1]
    c = CONTRASTS[1:]
    ymax = 1.12 * max(crf.max(), rb)
    ax.plot(c, crf[1:], color=fs.COND["dominant"], marker="o", ms=2.6, lw=1.0)
    ax.plot([0.45], [rb], color=fs.INK, marker="D", ms=3.4, zorder=5)
    ax.plot([0.45, kc], [rb, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([kc, kc], [0, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([(1 + a_med) * 0.45], [0], marker="^", ms=3.8, color=fs.MUTED, clip_on=False, zorder=6)
    fs.log_contrast_axis(ax, (0.055, 0.11, 0.225, 0.45, 0.9))
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"], fontsize=5.5)
    ax.set_xlim(0.045, 1.05)
    ax.set_ylim(0, ymax)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xlabel("dominant-eye contrast (%)")
    ax.set_title(title, fontsize=6.5, color=fs.INK)
    ax.text(0.05, 0.97, f"both eyes at 45%\n= dominant eye at {100 * kc:.0f}%", transform=ax.transAxes,
            ha="left", va="top", fontsize=6, color=fs.INK, linespacing=1.0)
    return kc / 0.45


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
    for y, lab, dy, va in ((-0.05, "top of cortex", 0.04, "top"), (l45, "L4/5", -0.05, "bottom")):
        ax.axhline(y, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
        ax.text(117, y + dy, lab, ha="right", va=va, fontsize=6, color=fs.INK,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
    ax.set_xlabel("time from flash (ms)")
    ax.set_ylabel("depth below top of cortex (mm)")
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return dict(l45_mm=l45, shift_contacts=float(reg["shift"]), bad_contacts=int(np.sum(bad)),
                depth_range=(float(depth[ok].min()), float(depth[ok].max())))


def tier1_csd(ax, cax, monkey="J"):
    """Average monocular grating CSD (one eye at >= 45%, eyes pooled; the laminar audit's cache) over the tier-1
    penetrations of one monkey (laminar alignment v4), each aligned to its 4C/5 anchor, scaled by its squish factor s
    and divided by its own peak |CSD| (20-120 ms) before averaging (Blake, 2026-11-12: examples from a clean J session
    or the tier-1 average of one monkey)."""
    V = pd.read_csv(metadata_dir() / "laminar_alignment_v5.csv")          # v5: anchor at the initial sink's bottom
    V = V[(V.tier == 1) & (V.monkey == monkey)]
    grid = np.round(np.arange(-0.6, 1.45, 0.05), 3)
    stack = []
    for r in V.itertuples():
        M = np.load(out_root() / "laminar_audit" / f"{r.penetration}_mono.npz")
        t = M["t"]
        c = (M["csd_contra"] + M["csd_ipsi"]) / 2
        c = c - np.nanmean(c[t < 0], axis=0)
        sign = 1 if r.orientation == "ascending" else -1
        h = -sign * (M["contact"] - r.anchor_contact) * r.spacing_mm / r.s
        ok = np.isfinite(c).all(axis=0)
        o = np.argsort(h[ok])
        C = np.stack([np.interp(grid, h[ok][o], c[k, ok][o], left=np.nan, right=np.nan) for k in range(t.size)])
        stack.append(C / np.nanmax(np.abs(C[(t >= 20) & (t < 120)])))
    A = np.nanmean(np.stack(stack), axis=0).T                       # [height, time]
    lim = np.nanpercentile(np.abs(A[:, (t > 20) & (t < 120)]), 98)
    im = ax.pcolormesh(t, grid, A, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, shading="nearest", rasterized=True)
    ax.set_xlim(-10, 120)
    ax.set_ylim(grid[0], grid[-1])
    for y, lab in ((1.24, "top of cortex"), (0.59, "top of layer 4"), (0.0, "layer 4C/5")):     # scaled heights
        ax.axhline(y, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
        ax.text(117, y + 0.03, lab, ha="right", va="bottom", fontsize=6, color=fs.INK,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("height above the layer 4C/5\nboundary (mm, scaled)")
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return dict(n=len(V), penetrations=list(V.penetration))


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet")
    a_med = float(meta.a.median())
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration").loc[pen]
    fig = fs.figure("double", 64)
    sa, sb = fig.subfigures(1, 2, width_ratios=[1.25, 1], wspace=0.03)
    ax_a1, ax_a2 = sa.subplots(1, 2)
    k_early = equivalent_contrast_sketch(ax_a1, pop, (0, 20), "first 20 ms after onset", a_med)
    k_late = equivalent_contrast_sketch(ax_a2, pop, (50, 200), "50–200 ms after onset", a_med)
    ax_a1.set_ylabel("MUA (dominant-eye peak = 1)")
    for k, (mk, col, lab) in enumerate((("D", fs.INK, "both eyes, 45% each"), ("o", fs.COND["dominant"], "dominant eye alone"),
                                        ("^", fs.MUTED, fs.SUM))):
        y = 0.75 - k * 0.075
        ax_a1.plot([0.07], [y], marker=mk, color=col, ms=2.8, transform=ax_a1.transAxes, clip_on=False)
        ax_a1.text(0.12, y, lab, transform=ax_a1.transAxes, va="center", fontsize=5.5, color=fs.INK)
    cax_ = sb.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 0.06]))
    info = tier1_csd(cax_[0], cax_[1], "J")
    cax_[0].set_title(f"monkey J, mean of {info['n']} tier-1 penetrations; one eye at 45 or 90%", loc="right",
                      fontsize=6, color=fs.INK2)
    for ax, letter in zip((ax_a1, cax_[0]), "ab"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS1_methods"))
    print(f"(a) illustration: equivalent contrast / 45%: early {k_early:.2f}, late {k_late:.2f}; median eye weight "
          f"{a_med:.2f}; {len(meta)} penetrations")
    print(f"(b) tier-1 average grating CSD, monkey J: {info['n']} penetrations {info['penetrations']}")


if __name__ == "__main__":
    main()
