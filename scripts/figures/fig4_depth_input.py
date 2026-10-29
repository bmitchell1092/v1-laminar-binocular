"""Figure 4: the same sequence runs at every depth and starts in the earliest synaptic input (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_depth_input.py [example penetration, default 221212_J_eA]
(a) example penetration: grating-evoked CSD (presentations with either eye >= .45), the earliest current sink
    (dichoptic_csd.py) and the L4/5 landmark (laminar alignment v2);
(b) linearity L (window-specific eye weights; dioptic .11 and .225 pooled) by depth below the top of cortex, in
    the early (0-20 ms after each channel's onset) and late (50-200) epochs; median of penetration
    medians, 95% bootstrap CI (as fig3_depth.py);
(c) equivalent contrast k of the earliest sink (CSD) and of MUA at the same channels, windows timed from the sink
    onset, at .225 and .45; median of penetrations with 95% bootstrap CI (dichoptic_csd.py);
(d-f) channels above, at (within 2 contacts) and below the earliest sink (laminar_timing.py): response onset on the
    flash clock; time from each channel's own onset to normalization (first window with L < 0.5); and that time on
    the flash clock (15 ms windows, so times are quantized). Lines = penetrations (values beyond the axis are
    drawn at its edge as triangles); Wilcoxon tests against the sink group, uncorrected.
Inputs: datasets, laminar alignment v2, k_windows.parquet, csd_windows.parquet, laminar_timing.parquet,
docs/qc/dichoptic_csd_sinks.csv.
Writes <out_root>/figures/fig4_depth_input.pdf and .png; prints the statistics shown.
"""

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import metadata_dir, out_root  # noqa: E402
from v1lb.dataset import clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
EPOCHS = {"early": (0, 20), "transition": (20, 50), "late": (50, 200)}
EPOCH_COLOR = {"early": fs.MUTED, "late": fs.INK}
LIMITS = {"onset": (20, 65), "t_rel": (0, 100), "t_abs": (45, 125)}       # values beyond are drawn at the edge
DEPTH_BINS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 32)]     # contacts below the top
CSD_WINDOWS = [(0, 10), (10, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150)]
GROUPS = ["above sink", "at sink", "below sink"]


def example_csd(ax, cax, pen, sink):
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration").loc[pen]
    ds = open_penetration(pen)
    d = out_root() / "datasets"
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E) & ((P.contrast_contra >= 0.45) | (P.contrast_ipsi >= 0.45)).to_numpy()
    t = ds.time.values
    inc = ds.in_cortex.values
    C = ds.csd.values[keep].astype(float).mean(axis=0)[:, inc]
    C -= C[(t >= -50) & (t < 0)].mean(axis=0)
    depth = ds.depth_below_toc.values[inc] * 0.1
    ok = np.isfinite(C).all(axis=0)
    lim = np.nanpercentile(np.abs(C[(t > 20) & (t < 150)][:, ok]), 98)
    im = ax.pcolormesh(t, depth[ok], C[:, ok].T, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, shading="nearest",
                       rasterized=True)
    ax.set_xlim(-10, 150)
    ax.set_ylim(depth[ok].max() + 0.05, depth[ok].min() - 0.05)
    ys = sink.sink_depth_below_toc * 0.1
    ax.plot([sink.sink_onset_ms], [ys], marker=">", ms=4, color=fs.INK, clip_on=False)
    ax.text(sink.sink_onset_ms - 3, ys, "earliest\nsink", ha="right", va="center", fontsize=5.5, color=fs.INK,
            linespacing=0.9, bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    l45 = (al.l45 - al.top_of_cortex) * 0.1
    ax.axhline(l45, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
    ax.text(147, l45 - 0.04, "L4/5", ha="right", va="bottom", fontsize=5.5, color=fs.INK,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth below top of cortex (mm)")
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)


def depth_profile(ax, W):
    W = W[W.c.isin([0.11, 0.225])].copy()
    W["epoch"] = np.select([W.w0 < 20, W.w0 < 50, W.w0 >= 50], list(EPOCHS), default="")
    W["bin"] = pd.cut(W.depth, [b[0] for b in DEPTH_BINS] + [DEPTH_BINS[-1][1]], right=False,
                      labels=range(len(DEPTH_BINS)))
    per = W.groupby(["penetration", "bin", "epoch"], observed=True).L_w.median().unstack(["epoch", "bin"])
    y = np.array([(a + b - 1) / 2 for a, b in DEPTH_BINS]) * 0.1
    for k, ep in enumerate(EPOCH_COLOR):
        M = per[ep].reindex(columns=range(len(DEPTH_BINS))).to_numpy()
        m, lo, hi = fs.boot_ci(M, stat=np.nanmedian)
        n = np.isfinite(M).sum(axis=0)
        good = n >= 5
        yy = y[good] + (k - 0.5) * 0.04
        ax.errorbar(m[good], yy, xerr=[m[good] - lo[good], hi[good] - m[good]], color=EPOCH_COLOR[ep], lw=1.1,
                    elinewidth=0.6, marker="o", ms=2.6, capsize=0)
    ax.axvline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_ylim(3.0, -0.1)
    ax.set_xlim(-0.8, 2.0)
    ax.set_xlabel("linearity index L")
    ax.set_ylabel("depth below top of cortex (mm)")
    fs.line_key(ax, [(f"{ep}, {EPOCHS[ep][0]}–{EPOCHS[ep][1]} ms", EPOCH_COLOR[ep], "-", "o") for ep in EPOCH_COLOR],
                x=0.6, y=0.97, dy=0.07, length=0.07)


def sink_vs_mua(ax, R, c, a_sink):
    centres = np.array([(a + b) / 2 for a, b in CSD_WINDOWS])
    for sig, col, ls, mfc, dx in (("csd_sink", fs.INK, (0, (3, 1.5)), "white", -1.5),
                                  ("mua_sink_clock", fs.INK, "-", fs.INK, 1.5)):
        per = R[(R.signal == sig) & (R.c == c)].groupby(["penetration", "w0"]).k.median().unstack()
        per = per.reindex(columns=[w for w, _ in CSD_WINDOWS])
        m, lo, hi = fs.boot_ci(per.to_numpy(), stat=np.nanmedian)
        ax.errorbar(centres + dx, m, [m - lo, hi - m], color=col, ls=ls, lw=1.0, elinewidth=0.6, marker="o", ms=2.8,
                    mfc=mfc, capsize=0)
    ax.axhline(1 + a_sink, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax, 1.0)
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 3.2)
    ax.set_xlabel("time from sink onset (ms)")
    ax.set_title(fs.contrast_label(c), fontsize=6.5)


def grouped(ax, T, var, ylabel):
    piv = T.pivot_table(index="penetration", columns="grp", values=var).reindex(columns=GROUPS)
    lo, hi = LIMITS[var]
    for _, r in piv.iterrows():
        v = r.values.clip(lo, hi)
        ax.plot(range(3), v, color=fs.RULE, lw=0.5, zorder=1)
        for k, (x, raw) in enumerate(zip(v, r.values)):
            out_of = np.isfinite(raw) and (raw > hi or raw < lo)
            ax.plot([k], [x], marker=("^" if raw > hi else "v") if out_of else "o", ms=2.4 if out_of else 1.5,
                    color=fs.MUTED if out_of else fs.RULE, zorder=2, clip_on=False)
    med = piv.median()
    ax.plot(range(3), med.values, color=fs.INK, lw=1.4, marker="o", ms=3.2, zorder=3)
    ax.set_xticks(range(3), ["above", "at", "below"])
    ax.set_xlabel("position relative to earliest sink")
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.4, 2.4)
    ax.set_ylim(lo, hi)
    out = {}
    for k, g in ((0, "above sink"), (2, "below sink")):
        ok = piv[["at sink", g]].dropna()
        d = ok[g] - ok["at sink"]
        p = stats.wilcoxon(ok[g], ok["at sink"]).pvalue if (d != 0).any() and len(ok) > 5 else np.nan
        ax.text(k, hi + 0.05 * (hi - lo), f"{d.median():+.1f} ms vs at\n{fs.fmt_p(p)}, n = {len(ok)}", ha="center",
                va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
        out[g] = dict(n=len(ok), median_diff=d.median(), p=p)
    ax.set_ylim(lo, hi)
    return med, out


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    R = pd.read_parquet(out_root() / "binocular" / "csd_windows.parquet")
    T = pd.read_parquet(out_root() / "binocular" / "laminar_timing.parquet")
    a_sink = float(R[(R.signal == "csd_sink") & (R.w0 < 20)].groupby("penetration").a.median().median())

    fig = fs.figure("double", 118)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.05, 1], hspace=0.04)
    ta = top.subplots(1, 5, gridspec_kw=dict(width_ratios=[1.05, 0.05, 0.95, 0.75, 0.75]))
    example_csd(ta[0], ta[1], pen, S.loc[pen])
    depth_profile(ta[2], W)
    sink_vs_mua(ta[3], R, 0.225, a_sink)
    sink_vs_mua(ta[4], R, 0.45, a_sink)
    ta[3].set_ylabel("equivalent contrast k")
    ta[4].tick_params(labelleft=False)
    ta[4].text(153, 1 + a_sink, "linear\nsum", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    ta[4].text(153, 1.0, "dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2,
               linespacing=0.95)
    fs.line_key(ta[3], [("current sink (CSD)", fs.INK, (0, (3, 1.5)), "o", "white"),
                        ("MUA, same channels", fs.INK, "-", "o")], x=0.3, y=0.96, dy=0.07, length=0.1)
    ba = bot.subplots(1, 3)
    m_on, s_on = grouped(ba[0], T, "onset", "response onset (ms after stimulus)")
    m_rel, s_rel = grouped(ba[1], T, "t_rel", "normalization (ms after own onset)")
    m_abs, s_abs = grouped(ba[2], T, "t_abs", "normalization (ms after stimulus)")
    for ax, letter in zip((ta[0], ta[2], ta[3], ba[0], ba[1], ba[2]), "abcdef"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig4_depth_input"))
    print(f"sink eye weight a (0-20 ms) median {a_sink:.2f}")
    for lab, m, s in (("onset", m_on, s_on), ("t_rel", m_rel, s_rel), ("t_abs", m_abs, s_abs)):
        print(lab, m.round(1).to_dict(), {g: {k: round(v, 3) for k, v in d.items()} for g, d in s.items()})


if __name__ == "__main__":
    main()
