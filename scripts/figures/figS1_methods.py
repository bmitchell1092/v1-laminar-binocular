"""Supplementary Figure S1: apparatus and laminar alignment (docs/09, figure set for V3; formerly Fig 1a and 1c).

Usage: .venv\\Scripts\\python scripts\\figures\\figS1_methods.py [example penetration, default 221212_J_eA]
(a) mirror stereoscope (schematic, top view): two mirrors per eye shift each eye's line of sight onto its own half of
    one monitor, divided by a black septum;
(b) example penetration (the one in Fig 1b): flash-evoked CSD (bad contacts interpolated, Hamming-smoothed),
    registered to the grating recording by the alignment's depth shift, with the top of cortex and the L4/5
    boundary below the flash sink (laminar alignment v2). The figure names the monkey, not the session.
Inputs: <out_root>/flash_csd/*.npz, docs/qc/flash_landmarks_v2.csv, metadata/laminar_alignment_v2.csv.
Writes <out_root>/figures/figS1_methods.pdf and .png; prints the landmarks shown.
"""

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


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration").loc[pen]
    fig = fs.figure("onehalf", 62)
    sa, sb = fig.subfigures(1, 2, width_ratios=[0.8, 1.2], wspace=0.03)
    sax = sa.subplots()
    stereoscope(sax)
    sax.set_anchor("N")
    cax_ = sb.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 0.06]))
    info = flash_csd(cax_[0], cax_[1], pen, al)
    cax_[0].set_title(f"example penetration, monkey {pen.split('_')[1]}", loc="right", fontsize=6, color=fs.INK2)
    for ax, letter in zip((sax, cax_[0]), "ab"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS1_methods"))
    print(f"(b) L4/5 boundary {info['l45_mm']:.2f} mm below the top of cortex; flash recording shifted by "
          f"{info['shift_contacts']:+.0f} contacts; {info['bad_contacts']} bad contacts interpolated; depths "
          f"{info['depth_range'][0]:.2f} to {info['depth_range'][1]:.2f} mm")


if __name__ == "__main__":
    main()
