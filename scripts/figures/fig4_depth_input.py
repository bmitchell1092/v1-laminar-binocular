"""Figure 4: the same sequence runs at every depth and starts in the earliest synaptic input (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_depth_input.py
(a) population CSD aligned to each penetration's earliest current sink (presentations with either eye >= .45; each
    penetration scaled by its sink's peak, then averaged; sink_data.py), with the median position of the L4/5
    boundary relative to the sink;
(b) linearity index by depth below the top of cortex (window-specific eye weights; dioptic .11 and .225 pooled),
    early (0-20 ms after each channel's onset) and late (50-200 ms); median of penetration medians, 95% bootstrap
    CI (as fig3_depth.py);
(c) binocular gain (both eyes minus the dominant eye alone, at the same contrast) of the sink current and of MUA on
    the same channels, on the sink-onset clock, each scaled by its dominant-eye peak at .9; median and 95% bootstrap
    CI over penetrations (medians: a few penetrations with small sinks dominate the mean; sink_data.py);
(d) the laminar sequence (laminar_timing.py): for channels above, at (within 2 contacts) and below the earliest sink,
    the response onset and the time normalization is reached (first 15 ms window with linearity index < 0.5), on the
    stimulus clock; dots = penetrations, bars = medians with interquartile ranges.
Inputs: sink_data.py (sink.npz), k_windows.parquet, laminar_timing.parquet, docs/qc/dichoptic_csd_sinks.csv.
Writes <out_root>/figures/fig4_depth_input.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
EPOCHS = {"early": (0, 20), "late": (50, 200)}
EPOCH_COLOR = {"early": fs.MUTED, "late": fs.INK}
DEPTH_BINS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 32)]     # contacts below the top
GROUPS = ["above sink", "at sink", "below sink"]
GROUP_LABEL = {"above sink": "above", "at sink": "at sink", "below sink": "below"}


def sink_csd(ax, cax, Z, S):
    M = np.nanmean(Z["maps"], axis=0)
    off, t = Z["offsets"] * 0.1, Z["t_map"]
    lim = np.nanpercentile(np.abs(M[:, (t > 20) & (t < 150)]), 99)
    im = ax.pcolormesh(t, off, M, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, shading="gouraud", rasterized=True)
    ax.set_ylim(off.max(), off.min())
    ax.set_xlim(-10, 150)
    l45 = np.nanmedian(S.sink_depth_l45) * 0.1                # L4/5 lies this far below the sink (median)
    ax.axhline(l45, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
    ax.text(147, l45 - 0.03, "L4/5 (median)", ha="right", va="bottom", fontsize=5.5, color=fs.INK,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.plot([-10], [0], marker=">", ms=4, color=fs.INK, clip_on=False, zorder=5)
    ax.text(-6, -0.05, "earliest sink", ha="left", va="bottom", fontsize=5.5, color=fs.INK,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to earliest sink (mm)")
    ax.text(0.97, 0.03, f"{Z['maps'].shape[0]} penetrations", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.5, color=fs.INK2, bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return l45


def depth_profile(ax, W):
    W = W[W.c.isin([0.11, 0.225])].copy()
    W["epoch"] = np.select([W.w0 < 20, W.w0 < 50, W.w0 >= 50], ["early", "transition", "late"], default="")
    W["bin"] = pd.cut(W.depth, [b[0] for b in DEPTH_BINS] + [DEPTH_BINS[-1][1]], right=False,
                      labels=range(len(DEPTH_BINS)))
    per = W.groupby(["penetration", "bin", "epoch"], observed=True).L_w.median().unstack(["epoch", "bin"])
    y = np.array([(a + b - 1) / 2 for a, b in DEPTH_BINS]) * 0.1
    for k, ep in enumerate(EPOCH_COLOR):
        M = per[ep].reindex(columns=range(len(DEPTH_BINS))).to_numpy()
        m, lo, hi = fs.boot_ci(M, stat=np.nanmedian)
        good = np.isfinite(M).sum(axis=0) >= 5
        ax.errorbar(m[good], y[good] + (k - 0.5) * 0.04, xerr=[m[good] - lo[good], hi[good] - m[good]],
                    color=EPOCH_COLOR[ep], lw=1.1, elinewidth=0.6, marker="o", ms=2.6, capsize=0)
    ax.axvline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_ylim(3.0, -0.1)
    ax.set_xlim(-0.8, 2.0)
    ax.set_xlabel("linearity index")
    ax.set_ylabel("depth below top of cortex (mm)")
    fs.line_key(ax, [(f"{ep}, {EPOCHS[ep][0]}–{EPOCHS[ep][1]} ms", EPOCH_COLOR[ep], "-", "o") for ep in EPOCH_COLOR],
                x=0.6, y=0.97, dy=0.07, length=0.07)


def sink_gain(ax, Z, c):
    rel = Z["rel"]
    i = list(CONTRASTS).index(c)
    sel = (rel >= -10) & (rel <= 150)
    for sig, ls in (("csd", (0, (3, 1.5))), ("mua", "-")):
        X = Z[sig]
        fs.ci_line(ax, rel[sel], (X[:, i, i] - X[:, i, 0])[:, sel], fs.INK, ls=ls, alpha=0.12, stat=np.nanmedian)
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xlim(-10, 150)
    ax.set_xlabel("time from sink onset (ms)")
    ax.set_title(fs.contrast_label(c), fontsize=6.5)


def sequence(ax, T):
    rng = np.random.default_rng(3)
    out = {}
    for row, g in enumerate(GROUPS):
        G = T[T.grp == g]
        y = row
        for var, mfc in (("onset", "white"), ("t_abs", fs.INK)):
            v = G[var].dropna()
            ax.scatter(v.clip(15, 130), y + rng.uniform(-0.18, 0.18, len(v)), s=4, color=fs.RULE, zorder=1, lw=0)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [y, y], color=fs.INK, lw=0.7, zorder=2)
            ax.plot([md], [y], marker="o", ms=4.5, mfc=mfc, mec=fs.INK, mew=0.8, zorder=4)
            out[(g, var)] = md
        ax.plot([out[(g, "onset")], out[(g, "t_abs")]], [y, y], color=fs.INK, lw=1.6, zorder=3)
    ax.set_yticks(range(3), [GROUP_LABEL[g] for g in GROUPS])
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(15, 130)
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.tick_params(axis="y", length=0)
    ax.set_ylabel("channels relative to earliest sink")
    for k, (mfc, lab) in enumerate((("white", "response onset"), (fs.INK, "normalized (linearity index < 0.5)"))):
        ax.plot([0.5], [0.97 - k * 0.07], marker="o", ms=4, mfc=mfc, mec=fs.INK, mew=0.8, transform=ax.transAxes,
                clip_on=False)
        ax.text(0.53, 0.97 - k * 0.07, lab, transform=ax.transAxes, va="center", fontsize=5.5, color=fs.INK)
    return out


def tests(T):
    rows = []
    for var in ("onset", "t_abs"):
        piv = T.pivot_table(index="penetration", columns="grp", values=var)
        for g in ("above sink", "below sink"):
            ok = piv[["at sink", g]].dropna()
            d = ok[g] - ok["at sink"]
            rows.append(dict(variable=var, group=g, n=len(ok), median_diff=d.median(),
                             p=stats.wilcoxon(ok[g], ok["at sink"]).pvalue))
    return pd.DataFrame(rows)


def main():
    fs.apply()
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    T = pd.read_parquet(out_root() / "binocular" / "laminar_timing.parquet")
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))

    fig = fs.figure("double", 112)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.05, 1], hspace=0.05)
    ta = top.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.25, 0.05, 1]))
    l45 = sink_csd(ta[0], ta[1], Z, S)
    depth_profile(ta[2], W)
    ba = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[0.75, 0.75, 1.3]))
    sink_gain(ba[0], Z, 0.225)
    sink_gain(ba[1], Z, 0.45)
    ba[1].sharey(ba[0])
    ba[0].set_ylabel("both eyes − dominant eye\n(dominant-eye peak = 1)")
    ba[1].tick_params(labelleft=False)
    fs.line_key(ba[0], [("current sink (CSD)", fs.INK, (0, (3, 1.5))), ("MUA, same channels", fs.INK, "-")],
                x=0.3, y=0.95, dy=0.08, length=0.13)
    med = sequence(ba[2], T)
    for ax, letter in zip((ta[0], ta[2], ba[0], ba[2]), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig4_depth_input"))
    print(f"L4/5 median {l45:.2f} mm below the sink; sinks {len(S)}")
    print({k: round(v, 1) for k, v in med.items()})
    print(tests(T).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
