"""Figure 5: a delayed binocular normalization pool explains the direction, not the strength, of the collapse
(docs/09, model tests 1-3).

Usage: .venv\\Scripts\\python scripts\\figures\\fig5_models.py
(a) the single-stage delayed normalization model (v1lb.dynorm 'full'); the two-stage variant adds interocular gain
    control on each eye first;
(b) held-out (odd/even) error of static vs delayed normalization per penetration (fig5_model.py);
(c) fitted pool weight w vs excitatory weight a of the non-dominant eye per penetration (fig5_model.py);
(d) two-stage vs single-stage: per penetration, median % reduction of held-out error over channels, for the
    registered raw-error test (fig5_channel_test.py) and the interaction-loss test (fig5_interaction_test.py);
    one-sided Wilcoxon over penetrations; values above the axis are drawn at its edge (open triangles);
(e) the registered capture test (fig5_interaction_test.py): late L (50-150 ms after onset) of each model's fit minus
    the data, median over penetrations with its 90% bootstrap CI; a model captures the collapse only if every CI
    lies inside the +-0.25 band. None does; positive = the model is too linear.
Inputs: <out_root>/binocular/fig5_model.parquet, fig5_channel_test.parquet, fig5_interaction_test.parquet,
fig5_interaction_test_L.parquet.
Writes <out_root>/figures/fig5_models.pdf and .png; prints the statistics shown.
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
GAIN_MAX = 8.0                         # % axis limit in (d)
MODELS = (("full", "single-stage", "o", fs.MODEL["fixed-delay"]), ("full2", "two-stage", "s", "white"),
          ("ds", "Ding–Sperling", "^", "white"))


def schematic(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec=fs.INK2, lw=0.6)
    arrow = dict(arrowstyle="-|>", lw=0.6, color=fs.INK2, mutation_scale=6, shrinkA=0, shrinkB=0)

    def node(x, y, s):
        ax.text(x, y, s, ha="center", va="center", fontsize=6, color=fs.INK, bbox=box)

    def arr(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=arrow)

    ax.text(0.0, 0.74, "dominant\neye $c_D$", ha="left", va="center", fontsize=6, color=fs.COND["dominant"],
            linespacing=0.95)
    ax.text(0.0, 0.3, "other\neye $c_N$", ha="left", va="center", fontsize=6, color=fs.COND["other"],
            linespacing=0.95)
    node(0.43, 0.74, "drive  $E = c_D + a\\,c_N$")
    node(0.35, 0.3, "pool  $c_D + w\\,c_N$")
    node(0.675, 0.3, "low-pass $\\tau_2$")
    node(0.89, 0.52, "$R = \\dfrac{E^n}{\\sigma^n + P^n}$")
    arr(0.15, 0.74, 0.255, 0.74)
    arr(0.15, 0.71, 0.225, 0.36)
    arr(0.15, 0.3, 0.215, 0.3)
    arr(0.15, 0.33, 0.255, 0.69)
    arr(0.485, 0.3, 0.565, 0.3)
    arr(0.61, 0.74, 0.79, 0.58)
    arr(0.775, 0.33, 0.8, 0.44)
    ax.text(0.71, 0.7, "fast", fontsize=5.5, color=fs.INK2, ha="center")
    ax.text(0.675, 0.17, "delayed pool $P(t)$", fontsize=5.5, color=fs.INK2, ha="center", va="top")
    ax.text(0.0, 0.04, "two-stage: each eye is first divided by $S + c_{self} + \\omega\\,c_{other}$ (own delay)",
            fontsize=5.5, color=fs.INK2, ha="left", va="center")


def static_vs_delayed(ax, M):
    for mk, g in M.groupby("monkey"):
        ax.scatter(g.cv_static, g.cv_full, s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                   lw=0.4, zorder=3)
    lo = min(M.cv_full.min(), M.cv_static.min()) / 1.5
    hi = max(M.cv_full.max(), M.cv_static.max()) * 1.5
    ax.plot([lo, hi], [lo, hi], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("held-out error, static")
    ax.set_ylabel("held-out error, delayed")
    worst = M.loc[(M.cv_full / M.cv_static).idxmax()]
    ax.annotate(f"{worst.penetration[:6]} (few trials)", xy=(worst.cv_static, worst.cv_full), xytext=(-4, 0),
                textcoords="offset points", ha="right", va="center", fontsize=5.5, color=fs.INK2)
    better = int((M.cv_full < M.cv_static).sum())
    p = stats.wilcoxon(M.cv_full, M.cv_static).pvalue
    ax.text(0.97, 0.06, f"delayed better\n{better}/{len(M)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
    return better, p


def weights(ax, M):
    for mk, g in M.groupby("monkey"):
        ax.scatter(g.p_a, g.p_w, s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4,
                   zorder=3)
    top = max(1.6, float(M[["p_a", "p_w"]].max().max()) * 1.06)
    ax.plot([0, top], [0, top], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal")
    ax.set_xlabel("excitatory weight a")
    ax.set_ylabel("pool weight w")
    n = int((M.p_w > M.p_a).sum())
    p = stats.wilcoxon(M.p_w, M.p_a).pvalue
    ax.text(0.04, 0.96, f"w > a\n{n}/{len(M)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.5, color=fs.INK, linespacing=1.0)
    for k, (mk, lab) in enumerate((("o", "monkey B"), ("s", "monkey J"))):
        y = top * (0.2 - k * 0.08)
        ax.scatter([top * 0.64], [y], s=9, marker=mk, facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4)
        ax.text(top * 0.69, y, lab, va="center", fontsize=5.5, color=fs.INK)
    return n, p


def two_stage(ax, C, I):
    rows = []
    for x, (lab, D) in enumerate((("raw error\n(registered)", C), ("interaction\nloss", I))):
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
    ax.set_xticks([0, 1], ["raw error\n(registered)", "interaction\nloss"], fontsize=5.5)
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
    ax.set_ylabel("late L, model − data")
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
    M = pd.read_parquet(b / "fig5_model.parquet")
    C = pd.read_parquet(b / "fig5_channel_test.parquet")
    I = pd.read_parquet(b / "fig5_interaction_test.parquet")
    L = pd.read_parquet(b / "fig5_interaction_test_L.parquet")

    fig = fs.figure("double", 100)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1], hspace=0.05)
    ta = top.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.7, 1, 1]))
    schematic(ta[0])
    nb, pb = static_vs_delayed(ta[1], M)
    nw, pw = weights(ta[2], M)
    ba = bot.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 1.5]))
    T = two_stage(ba[0], C, I)
    K = capture(ba[1], L)
    for ax, letter in zip((ta[0], ta[1], ta[2], ba[0], ba[1]), "abcde"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig5_models"))
    print(f"delayed beats static {nb}/{len(M)} (p {pb:.2g}); w > a {nw}/{len(M)} (p {pw:.2g}); "
          f"median w {M.p_w.median():.2f}, a {M.p_a.median():.2f}, tau2 {M.p_tau2.median():.0f} ms")
    print(T.round(4).to_string(index=False))
    print(K.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
