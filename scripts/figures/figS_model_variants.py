"""Supplementary figure: model variants that add interocular gain control before binocular normalization
(docs/09, model tests 2-3; formerly Fig 5d-e, moved when Figs 5 and 6 were merged into one model figure).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_model_variants.py
(a) two-stage vs single-stage: per penetration, median % reduction of held-out error over channels, for the
    registered raw-error test (fig5_channel_test.py) and the interaction-loss test (fig5_interaction_test.py);
    one-sided Wilcoxon over penetrations; values above the axis are drawn at its edge (open triangles);
(b) the registered capture test (fig5_interaction_test.py): late linearity index (50-150 ms after onset) of each
    model's fit minus the data, median over penetrations with its 90% bootstrap CI; a model captures the collapse only
    if every CI lies inside the +-0.25 band; positive = the model is too linear.
Inputs: <out_root>/binocular/fig5_channel_test.parquet, fig5_interaction_test.parquet,
fig5_interaction_test_L.parquet.
Writes <out_root>/figures/figS_model_variants.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
DIOPTIC = [0.11, 0.225, 0.45]
LATE = (50, 75, 100)                   # window starts, as registered for the interaction test
MARGIN = 0.25                          # registered equivalence margin
GAIN_MAX = 8.0                         # % axis limit in (a)
MODELS = (("full", "single-stage", "o", fs.MODEL["fixed-delay"]), ("full2", "two-stage", "s", "white"),
          ("ds", "Ding–Sperling", "^", "white"))



def two_stage(ax, C, I):
    rows = []
    for x, (lab, D) in enumerate((("raw\nerror", C), ("interaction\nloss", I))):
        r = np.log(D.cv_full / D.cv_full2)
        per = r.groupby([D.penetration, D.monkey]).median()
        pct = 100 * (1 - np.exp(-per))
        rng = np.random.default_rng(x)
        for mk, v in pct.groupby(level="monkey"):
            xs = x + rng.uniform(-0.15, 0.15, len(v))
            inside = v <= GAIN_MAX
            ax.scatter(xs[inside], v[inside], s=8, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                       lw=0.4, zorder=3)
            ax.scatter(xs[~inside], np.full((~inside).sum(), GAIN_MAX), s=10, marker="^", facecolor="white",
                       edgecolor=fs.INK, lw=0.5, zorder=3, clip_on=False)
        ax.plot([x - 0.25, x + 0.25], [pct.median()] * 2, color=fs.INK, lw=1.4, zorder=4)
        p = stats.wilcoxon(per.values, alternative="greater").pvalue
        ax.text(x, GAIN_MAX * 1.06, f"{int((per > 0).sum())}/{len(per)}, {fs.fmt_p(p)}\nmedian {pct.median():.1f}%",
                ha="center", va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
        rows.append(dict(test=lab.replace("\n", " "), n=len(per), positive=int((per > 0).sum()), median_pct=pct.median(),
                         p=p, clipped=int((pct > GAIN_MAX).sum())))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks([0, 1], ["raw\nerror", "interaction\nloss"], fontsize=5.5)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-2, GAIN_MAX)
    ax.set_ylabel("two-stage gain\n(% less held-out error)")
    return pd.DataFrame(rows)


def capture(ax, L):
    late = L[L.w0.isin(LATE)].groupby(["penetration", "channel", "source", "c"]).L.mean().unstack("source")
    rng = np.random.default_rng(0)
    ax.axhspan(-MARGIN, MARGIN, color=fs.RULE, alpha=0.6, lw=0, zorder=0)
    out = []
    for k, c in enumerate(DIOPTIC):
        g = late.xs(c, level="c")
        for j, (m, lab, mk, mfc) in enumerate(MODELS):
            per = (g[m] - g["data"]).dropna().groupby(level="penetration").median()
            boot = [np.median(rng.choice(per.values, len(per))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            x = k + (j - 1) * 0.22
            ax.errorbar([x], [per.median()], [[per.median() - lo], [hi - per.median()]], color=fs.MODEL["fixed-delay"],
                        marker=mk, mfc=mfc, ms=3.2, elinewidth=0.8, capsize=0, lw=0)
            out.append(dict(c=c, model=m, n=len(per), median=per.median(), lo=lo, hi=hi,
                            inside=(lo > -MARGIN) and (hi < MARGIN)))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks(range(3), ["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_ylabel("late linearity index,\nmodel − data")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.6, 0.8)
    ax.text(-0.45, MARGIN - 0.02, "equivalence margin", ha="left", va="top", fontsize=5.5, color=fs.INK2)
    ax.text(-0.45, 0.78, "model too linear", fontsize=5.5, color=fs.INK2, va="top")
    ax.text(-0.45, -0.58, "model too suppressive", fontsize=5.5, color=fs.INK2, va="bottom")
    fs.line_key(ax, [(lab, fs.MODEL["fixed-delay"], "none", mk, mfc) for _, lab, mk, mfc in MODELS], x=0.68, y=0.96,
                dy=0.075, length=0.05)
    return pd.DataFrame(out)


def main():
    fs.apply()
    b = out_root() / "binocular"
    C = pd.read_parquet(b / "fig5_channel_test.parquet")
    I = pd.read_parquet(b / "fig5_interaction_test.parquet")
    L = pd.read_parquet(b / "fig5_interaction_test_L.parquet")

    fig = fs.figure("double", 52)
    ax_a, ax_b = fig.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 1.5]))
    T = two_stage(ax_a, C, I)
    K = capture(ax_b, L)
    for ax, letter in zip((ax_a, ax_b), "ab"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS_model_variants"))
    print(T.round(4).to_string(index=False))
    print(K.round(3).to_string(index=False))

if __name__ == "__main__":
    main()
