"""Supplementary Figure S2: the linearity index is a valid measure of the binocular combination rule (docs/09, V3
analyses A1-A4).

Usage: .venv\\Scripts\\python scripts\\figures\\figS2_validity.py
(a) fixed-rule null through the identical readout (v3_null_model_index.py): each penetration's fitted fixed-delay
    model with the pool weight set to the drive weight, so that the binocular response equals the dominant eye's
    response at (1 + a)c at every moment, with realistic transient shapes; read out noise-free and with the data's
    noise sample, (odd - even) / 2, added and subtracted; the data in the same channel-average readout. Linearity
    index in onset-aligned windows (window centres), median over penetrations; 95% bootstrap CI of the data median;
(b) cross-validated dominance (v3_dominance_cv.py): dominant eye, eye weight and onset chosen on one half of the
    presentations, the index read from the other half (both directions averaged); per penetration, early (0-20 ms
    after onset) vs late (50-200 ms), per-channel readout of Fig 2d;
(c) ocular-dominance terciles at 45% (|ODI| from the selection half; index from the other half);
(d) isolated single units (lenient set; v3_units_index.py), per unit, same readout as Fig 2d.
(b-d): lines = penetrations (units in d), coloured line = median; Wilcoxon signed-rank tests, early vs late; values
beyond the axis are drawn at its edge as open triangles and counted.
Inputs: <out_root>/binocular/v3_A1_L.parquet, v3_A2_penetrations.parquet, v3_A3_penetrations.parquet,
v3_A4_units.parquet; docs/qc/v3_A2.csv, v3_A3.csv.
Writes <out_root>/figures/figS2_validity.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
DIOPTIC = (0.11, 0.225, 0.45)
WIN = {0: 10, 10: 20, 20: 35, 35: 50, 50: 75, 75: 100, 100: 150}       # channel-average readout windows (start: end)
EARLY_W, LATE_W = (0, 10), (50, 75, 100)                              # A1 early / late window starts
EARLY, LATE = [0, 10], [50, 75, 100, 150]                             # per-channel readout (Fig 2d)
LO, HI = -1.5, 3.0                                                    # axis limits of the paired panels
NULL_COL = fs.INK2


def null_panel(axs, L):
    out = []
    for ax, c in zip(axs, DIOPTIC):
        g = L[L.c == c].pivot_table(index=["penetration", "w0"], columns="source", values="L").reset_index()
        x = np.array([(w0 + WIN[w0]) / 2 for w0 in WIN])
        for src, col, ls, lw in (("null", NULL_COL, (0, (3, 1.5)), 1.0), ("null+noise", fs.MUTED, "-", 0.7),
                                 ("null-noise", fs.MUTED, "-", 0.7)):
            m = g.groupby("w0")[src].median().reindex(list(WIN)).to_numpy()
            ax.plot(x, m, color=col, ls=ls, lw=lw, marker="o" if "noise" in src else None, ms=1.8, zorder=2)
        D = g.pivot(index="penetration", columns="w0", values="data")[list(WIN)].to_numpy()
        m, lo, hi = fs.boot_ci(D, stat=np.nanmedian)
        ax.fill_between(x, lo, hi, color=fs.INK, alpha=0.12, lw=0)
        ax.plot(x, m, color=fs.INK, lw=1.2, marker="o", ms=2.4, zorder=3)
        ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
        fs.zero_line(ax)
        ax.set_xlim(0, 130)
        ax.set_xticks([0, 50, 100])
        ax.set_ylim(-1.2, 2.0)
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
        for src in ("data", "null", "null+noise", "null-noise"):
            med = g.groupby("w0")[src].median()
            out.append(dict(c=c, source=src, **{f"w{k}": round(float(v), 2) for k, v in med.items()}))
    for ax in axs[1:]:
        ax.tick_params(labelleft=False)
    axs[0].set_ylabel("linearity index")
    fs.line_key(axs[0], [("data", fs.INK, "-", "o"), ("fixed rule, noise-free", NULL_COL, (0, (3, 1.5))),
                         ("fixed rule + data noise (both signs)", fs.MUTED, "-", "o")], x=0.08, y=0.2, dy=0.09)
    axs[-1].text(132, 1, "linear\nsum", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    axs[-1].text(132, 0, "dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2,
                 linespacing=0.95)
    return pd.DataFrame(out)


def null_drop(L):
    """Early-minus-late drop at 45% per penetration: data vs each noisy null (paired one-sided Wilcoxon)."""
    g = L[L.c == 0.45].pivot_table(index=["penetration", "w0"], columns="source", values="L")
    e = g[g.index.get_level_values("w0").isin(EARLY_W)].groupby("penetration").mean()
    l_ = g[g.index.get_level_values("w0").isin(LATE_W)].groupby("penetration").mean()
    drop = (e - l_)
    rows = []
    for src in ("null", "null+noise", "null-noise"):
        d = drop[["data", src]].dropna()
        p = stats.wilcoxon(d.data, d[src], alternative="greater").pvalue
        rows.append(dict(null=src, n=len(d), data_drop=d.data.median(), null_drop=d[src].median(), p_one_sided=p))
    return pd.DataFrame(rows)


def paired(ax, groups, gap=2.6):
    """groups = [(label, colour, DataFrame with early/late columns)]; paired lines per item, median in colour."""
    out = []
    for k, (lab, col, d) in enumerate(groups):
        d = d[["early", "late"]].dropna()
        e, l_ = d.early.to_numpy(), d.late.to_numpy()
        x0, x1 = k * gap, k * gap + 1
        ec, lc = np.clip(e, LO, HI), np.clip(l_, LO, HI)
        for a, b in zip(ec, lc):
            ax.plot([x0, x1], [a, b], color=fs.RULE, lw=0.5, zorder=1)
        for xx, v, vc in ((x0, e, ec), (x1, l_, lc)):
            hi, lo = v > HI, v < LO
            ax.scatter(np.full(hi.sum(), xx), vc[hi], marker="^", s=8, facecolor="white", edgecolor=fs.INK2,
                       lw=0.5, zorder=2, clip_on=False)
            ax.scatter(np.full(lo.sum(), xx), vc[lo], marker="v", s=8, facecolor="white", edgecolor=fs.INK2,
                       lw=0.5, zorder=2, clip_on=False)
        p = stats.wilcoxon(e, l_).pvalue
        ax.plot([x0, x1], [np.median(e), np.median(l_)], color=col, lw=1.6, marker="o", ms=3.2, zorder=3)
        n_clip = int(((e > HI) | (e < LO)).sum() + ((l_ > HI) | (l_ < LO)).sum())
        ax.text((x0 + x1) / 2, HI + 0.2, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        ax.text((x0 + x1) / 2, LO - 0.6, lab, ha="center", va="top", fontsize=6, color=fs.INK, linespacing=0.95)
        ax.text((x0 + x1) / 2, -0.2, f"{n_clip} off axis" if n_clip else "", ha="center", va="top", fontsize=5,
                color=fs.INK2, transform=blended_transform_factory(ax.transData, ax.transAxes))
        out.append(dict(group=lab.replace("\n", " "), n=len(d), early=np.median(e), late=np.median(l_),
                        frac_drop=float((e > l_).mean()), p=p, clipped=n_clip))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for k in range(len(groups)) for x in (k * gap, k * gap + 1)], ["early", "late"] * len(groups),
                  fontsize=5.5)
    ax.set_xlim(-0.6, (len(groups) - 1) * gap + 1.6)
    ax.set_ylim(LO - 0.1, HI + 0.75)
    ax.set_yticks([-1, 0, 1, 2, 3])
    ax.spines["left"].set_bounds(LO - 0.1, HI + 0.1)
    ax.tick_params(axis="x", length=0)
    return pd.DataFrame(out)


def unit_early_late(U):
    g = U[U.variant == "all"]
    per = g.groupby(["penetration", "cluster_id", "direction", "c", "w0"]).L.first().unstack("w0")
    el = pd.DataFrame({"early": per[EARLY].mean(axis=1), "late": per[LATE].mean(axis=1)})
    return el.groupby(["penetration", "cluster_id", "c"]).mean()


def main():
    fs.apply()
    b = out_root() / "binocular"
    L = pd.read_parquet(b / "v3_A1_L.parquet")
    P2 = pd.read_parquet(b / "v3_A2_penetrations.parquet")
    P3 = pd.read_parquet(b / "v3_A3_penetrations.parquet")
    U = pd.read_parquet(b / "v3_A4_units.parquet")
    A2 = pd.read_csv(REPO / "docs" / "qc" / "v3_A2.csv")
    A3 = pd.read_csv(REPO / "docs" / "qc" / "v3_A3.csv")

    fig = fs.figure("double", 112)
    top, bot = fig.subfigures(2, 1, height_ratios=[0.9, 1], hspace=0.05)
    axs_a = top.subplots(1, 3, sharey=True)
    Sa = null_panel(axs_a, L)
    top.supxlabel("time from response onset (ms, window centre)", fontsize=7)
    drops = null_drop(L)
    ax_b, ax_c, ax_d = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[3, 3, 2.1]))

    cv = P2.xs("cv", level="variant")
    Sb = paired(ax_b, [(fs.contrast_label(c), fs.CONTRAST[c], cv.xs(c, level="c")) for c in DIOPTIC])
    sh = A2[(A2.variant == "cv vs same") & (A2.c == 0.45)].iloc[0]
    ax_b.set_ylabel("linearity index")
    ax_b.set_title("dominant eye chosen on the other half", fontsize=6.5)
    ax_b.text(0.02, 0.03, f"late index at 45%, cross-validated − same half:\n{sh.shift_late:+.3f} "
              f"[{sh.lo:.3f}, {sh.hi:.3f}]", transform=ax_b.transAxes, fontsize=5.2, color=fs.INK2, va="bottom",
              linespacing=1.0)

    t45 = P3.xs(0.45, level="c")
    names = {0: "most\nbinocular", 1: "middle", 2: "most\nmonocular"}
    rng_lab = A3[A3.c == 0.45].set_index("tercile").abs_odi_range
    Sc = paired(ax_c, [(names[t], fs.CONTRAST[0.45], t45.xs(t, level="tercile")) for t in (0, 1, 2)])
    print("(c) |ODI| ranges of the terciles:", rng_lab.to_dict())
    ax_c.set_title("ocular-dominance terciles, both eyes 45%", fontsize=6.5)
    ax_c.tick_params(labelleft=False)

    el = unit_early_late(U)
    Sd = paired(ax_d, [(fs.contrast_label(c), fs.CONTRAST[c], el.xs(c, level="c")) for c in (0.225, 0.45)])
    shown = el[el.index.get_level_values("c").isin([0.225, 0.45])].dropna()
    n_pen = shown.index.get_level_values("penetration").nunique()
    n_units = shown.reset_index()[["penetration", "cluster_id"]].drop_duplicates().shape[0]
    ax_d.set_title(f"{n_units} isolated units, {n_pen} penetrations", fontsize=6.5)
    print(f"(d) units shown: {n_units} units in {n_pen} penetrations")
    ax_d.tick_params(labelleft=False)

    for ax, letter in zip((axs_a[0], ax_b, ax_c, ax_d), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS2_validity"))
    pd.set_option("display.width", 220)
    print("(a) medians over penetrations per window start:")
    print(Sa.to_string(index=False))
    print("(a) early-minus-late drop at 45% (data vs null; paired one-sided Wilcoxon):")
    print(drops.round(4).to_string(index=False))
    for lab, S in (("(b) cross-validated", Sb), ("(c) OD terciles at 45%", Sc), ("(d) lenient units", Sd)):
        print(lab)
        print(S.round(4).to_string(index=False))
    print(f"(b) selection-bias shift of the late index at 45% (cv - same): {sh.shift_late:+.3f} [{sh.lo:.3f}, "
          f"{sh.hi:.3f}]")


if __name__ == "__main__":
    main()
