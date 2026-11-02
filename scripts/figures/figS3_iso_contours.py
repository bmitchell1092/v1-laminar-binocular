"""Supplementary Figure S3: iso-response contours at every target contrast, early, during the transition and late
(docs/09, V3 analysis C; the main-text version, 45% early and late without fitted contours, is Fig 2c).

Usage: .venv\\Scripts\\python scripts\\figures\\figS3_iso_contours.py
Contours from v3_iso_contours.py: for a target equal to the dominant eye's response at 11, 22.5 or 45% in a window
(0-20, 20-50 or 50-200 ms after each penetration's response onset), each other-eye contrast gives the dominant-eye
contrast at which the binocular response reaches the target. Axes are in percent contrast, as in Fig 2c; the other
eye's contrast is multiplied by the eye weight fitted in the same window, so linear summation is the diagonal and
"the dominant eye alone" the horizontal line at the target contrast.
(a) population contours (population-average window responses), 95% penetration-bootstrap CIs of the dominant-eye
    contrast; filled = measured; open circle on the horizontal axis = the other eye alone already reaches the target
    (the contour meets the axis at or before that point); open triangle = CI beyond the axis; lines = fitted Minkowski
    contours (exponent 1 = linear summation, large = the dominant eye alone);
(b) per-penetration Minkowski exponents (log axis; the fit is bounded to 0.5-20, values at a bound are drawn on
    it); lines = penetrations, black = median; Wilcoxon signed-rank test, late vs early.
Inputs: <out_root>/binocular/v3_C_population.parquet, v3_C_population_boot.parquet, v3_C_penetration.parquet.
Writes <out_root>/figures/figS3_iso_contours.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
TARGETS = (0.11, 0.225, 0.45)
WINDOWS = {"early": "0–20 ms", "mid": "20–50 ms", "late": "50–200 ms"}
SHADE = {"early": fs.MUTED, "mid": fs.INK2, "late": fs.INK}
XMAX, YMAX = 1.55, 1.33                               # axis limits as multiples of the target contrast (Fig 2c)


def minkowski(y, m):
    return (1 - np.clip(y, 0, 1) ** m) ** (1 / m)


def contour_panels(axs, Pop, Boot):
    out = []
    yy = np.linspace(0, 1, 101)
    for ax, cs in zip(axs, TARGETS):
        T = 100 * cs
        ax.plot([0, T], [T, 0], color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=1)
        ax.plot([0, XMAX * T], [T, T], color=fs.RULE, lw=0.9, zorder=0)
        for k, w in enumerate(WINDOWS):
            col = SHADE[w]
            p = Pop[(Pop.window == w) & (Pop.cstar == cs)].sort_values("c_N")
            b = Boot[(Boot.window == w) & (Boot.cstar == cs)].groupby("c_N").x
            lo = T * b.quantile(0.025).reindex(p.c_N).to_numpy()
            hi = T * b.quantile(0.975).reindex(p.c_N).to_numpy()
            r = p.iloc[0]
            ax.plot(T * yy, T * minkowski(yy, r.m), color=col, lw=0.8, zorder=2)
            show = p.y.to_numpy() <= XMAX
            xo, yd = T * (p.y.to_numpy()[show] + (k - 1) * 0.025), T * p.x.to_numpy()[show]
            lo, hi = np.minimum(lo[show], YMAX * T), np.minimum(hi[show], YMAX * T)
            st = p.status.to_numpy()[show]
            ax.vlines(xo, lo, hi, color=col, lw=0.6, zorder=3)
            ok = st == "ok"
            ax.scatter(xo[ok], yd[ok], s=9, color=col, zorder=4, lw=0)
            ax.scatter(xo[~ok], np.minimum(yd[~ok], YMAX * T), s=9, facecolor="white", edgecolor=col, lw=0.7,
                       zorder=4, clip_on=False)
            top = hi >= YMAX * T
            ax.scatter(xo[top], np.full(top.sum(), YMAX * T), marker="^", s=8, facecolor="white", edgecolor=col,
                       lw=0.6, zorder=4, clip_on=False)
            out.append(dict(target=cs, window=w, eye_weight=r.a, m=r.m, m_lo=r.m_lo, m_hi=r.m_hi,
                            suppression=r.suppression, below_axis=int((st == "below_axis").sum()),
                            censored=int((st == "censored").sum()), ci_beyond_axis=int(top.sum())))
        ax.set_xlim(0, XMAX * T)
        ax.set_ylim(0, YMAX * T)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"target: dominant eye alone at {fs.contrast_label(cs)}", fontsize=6.5)
        ax.text(T * 0.36, T * 0.70, "linear sum", rotation=-45, rotation_mode="anchor", fontsize=5.5,
                color=fs.INK2, ha="center", va="bottom")
        ax.text(XMAX * T, T * 1.02, "dominant eye alone", fontsize=5.5, color=fs.INK2, ha="right", va="bottom")
    axs[0].set_ylabel("dominant-eye contrast needed (%)")
    fs.line_key(axs[0], [(lab, SHADE[w], "-", "o") for w, lab in WINDOWS.items()], x=0.66, y=0.52, dy=0.08,
                length=0.1, fontsize=5.5)
    return pd.DataFrame(out)


def exponent_panel(ax, Pen):
    out = []
    gap = 3.6
    for g, cs in enumerate(TARGETS):
        M = Pen[Pen.cstar == cs].pivot(index="penetration", columns="window", values="m")[list(WINDOWS)].dropna()
        xs = g * gap + np.arange(3)
        for _, row in M.iterrows():
            ax.plot(xs, row.to_numpy(), color=fs.RULE, lw=0.5, zorder=1)
        ax.plot(xs, M.median().to_numpy(), color=fs.INK, lw=1.4, marker="o", ms=3, zorder=3)
        p = stats.wilcoxon(M.late, M.early).pvalue
        ax.text(xs[1], 26, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        ax.text(xs[1], 0.3, f"target {fs.contrast_label(cs)}", ha="center", va="top", fontsize=6.5, color=fs.INK)
        mk = M.index.str[7]
        out.append(dict(target=cs, n=len(M), early=M.early.median(), mid=M["mid"].median(), late=M.late.median(),
                        k_late_gt_early=int((M.late > M.early).sum()), p=p,
                        change_B=(M.late - M.early)[mk == "B"].median(), change_J=(M.late - M.early)[mk == "J"].median(),
                        at_upper_bound_late=int((M.late >= 20).sum()), at_upper_bound_early=int((M.early >= 20).sum())))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(0.45, 24)
    ax.set_yticks([0.5, 1, 2, 5, 10, 20], ["0.5", "1", "2", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_xticks([g * gap + i for g in range(len(TARGETS)) for i in range(3)],
                  ["0–20", "20–50", "50–200"] * len(TARGETS), fontsize=5.5)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, (len(TARGETS) - 1) * gap + 2.6)
    ax.set_ylabel("Minkowski exponent")
    ax.set_xlabel("window after response onset (ms)", labelpad=12)
    ax.text(ax.get_xlim()[1], 1, " linear sum", ha="left", va="center", fontsize=5.5, color=fs.INK2)
    return pd.DataFrame(out)


def main():
    fs.apply()
    b = out_root() / "binocular"
    Pop = pd.read_parquet(b / "v3_C_population.parquet")
    Boot = pd.read_parquet(b / "v3_C_population_boot.parquet")
    Pen = pd.read_parquet(b / "v3_C_penetration.parquet")
    fig = fs.figure("double", 118)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.05, 1], hspace=0.04)
    axs = top.subplots(1, 3)
    A = contour_panels(axs, Pop, Boot)
    top.supxlabel("other-eye contrast × eye weight (%)", fontsize=7)
    left, _ = bot.subfigures(1, 2, width_ratios=[1.6, 0.4])
    ax_b = left.subplots()
    B = exponent_panel(ax_b, Pen)
    fs.panel_label(axs[0], "a")
    fs.panel_label(ax_b, "b")
    print(fs.save(fig, "figS3_iso_contours"))
    pd.set_option("display.width", 220)
    print("(a) population contours:")
    print(A.round(3).to_string(index=False))
    print("(b) per-penetration exponents (medians; Wilcoxon late vs early, two-sided; monkey medians of late - early):")
    print(B.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
