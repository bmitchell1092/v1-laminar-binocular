"""Supplementary Figure S5: model checks (docs/09, figure set for V3; absorbs the former figS_model_variants.py and
panel g of the former fig5_models.py).

Usage: .venv\\Scripts\\python scripts\\figures\\figS5_model_checks.py
(a) two-stage vs single-stage (fixed-delay) normalization: per penetration, median % reduction of held-out error over
    channels, with the raw error (fig5_channel_test.py) and with the interaction loss (fig5_interaction_test.py);
    one-sided Wilcoxon over penetrations; values above the axis are drawn at its edge (open triangles);
(b) does a fitted model reproduce the late collapse? Late linearity index (windows starting 50-100 ms after onset) of
    each channel-level fit minus the data (fig5_interaction_test.py), median over penetrations with its 90% bootstrap
    CI, against a +-0.25 band; positive = the model is too linear;
(c) the same for the penetration-level fits of the three dynamic models (fig5_cond_model.py, fig5_org_model.py); at
    11% the data's own noise moves the index by +-0.4-0.7 (docs/09), so that contrast is drawn open;
(d) parameter recovery: pool weight minus drive weight of the other eye, refitted ('full') from synthetic data with
    equal weights plus the data's noise (100 sets x 25 penetrations, v3_param_recovery.py) and from the data
    (fig5_cond.parquet, full fits); bars = medians;
(e) separate pool delays for the two eyes' terms (model 'sep', P = c_D v(tau2) + w c_N v(tau2N); v3_interocular_delay.py):
    all-presentation fits per penetration; held-out comparison with the shared-delay model over 24 penetrations
    (one-sided Wilcoxon on log error ratios; one penetration, 211008, has no usable even half).
Inputs: <out_root>/binocular/fig5_channel_test.parquet, fig5_interaction_test.parquet,
fig5_interaction_test_L.parquet, fig5_cond_L.parquet, fig5_org_L.parquet, fig5_cond.parquet, v3_B3.parquet,
v3_B4_fits.parquet.
Writes <out_root>/figures/figS5_model_checks.pdf and .png; prints the statistics shown.
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
DIOPTIC = [0.11, 0.225, 0.45]
LATE = (50, 75, 100)                   # window starts
MARGIN = 0.25                          # equivalence band
GAIN_MAX = 8.0                         # % axis limit in (a)
VARIANTS = (("full", "single-stage (fixed-delay)", "o", fs.MODEL["fixed-delay"]), ("full2", "two-stage", "s", "white"),
            ("ds", "Ding–Sperling", "^", "white"))
DYNAMIC = {"full": ("fixed-delay", fs.MODEL["fixed-delay"]), "cond": ("conductance", fs.MODEL["conductance"]),
           "org": ("recurrent", fs.MODEL["recurrent"])}


def band_axes(ax, ylim):
    ax.axhspan(-MARGIN, MARGIN, color=fs.RULE, alpha=0.6, lw=0, zorder=0)
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks(range(3), ["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(*ylim)
    ax.text(-0.45, ylim[1] - 0.02, "model too linear", fontsize=5.5, color=fs.INK2, va="top")
    ax.text(-0.45, ylim[0] + 0.02, "model too suppressive", fontsize=5.5, color=fs.INK2, va="bottom")


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
        ax.text(x, GAIN_MAX * 1.04, f"{int((per > 0).sum())}/{len(per)}, {fs.fmt_p(p)}\nmedian {pct.median():.1f}%",
                ha="center", va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
        rows.append(dict(test=lab.replace("\n", " "), n=len(per), positive=int((per > 0).sum()),
                         median_pct=pct.median(), p=p, clipped=int((pct > GAIN_MAX).sum())))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks([0, 1], ["raw\nerror", "interaction\nloss"], fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-2, GAIN_MAX)
    ax.set_ylabel("two-stage vs single-stage\n(% less held-out error)")
    return pd.DataFrame(rows)


def capture(ax, L):
    late = L[L.w0.isin(LATE)].groupby(["penetration", "channel", "source", "c"]).L.mean().unstack("source")
    rng = np.random.default_rng(0)
    out = []
    for k, c in enumerate(DIOPTIC):
        g = late.xs(c, level="c")
        for j, (m, lab, mk, mfc) in enumerate(VARIANTS):
            per = (g[m] - g["data"]).dropna().groupby(level="penetration").median()
            boot = [np.median(rng.choice(per.values, len(per))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            x = k + (j - 1) * 0.22
            ax.errorbar([x], [per.median()], [[per.median() - lo], [hi - per.median()]], color=fs.MODEL["fixed-delay"],
                        marker=mk, mfc=mfc, ms=3.2, elinewidth=0.8, capsize=0, lw=0)
            out.append(dict(c=c, model=m, n=len(per), median=per.median(), lo=lo, hi=hi,
                            inside=(lo > -MARGIN) and (hi < MARGIN)))
    band_axes(ax, (-0.8, 0.8))
    ax.set_ylabel("late linearity index,\nmodel − data")
    fs.line_key(ax, [(lab, fs.MODEL["fixed-delay"], "none", mk, mfc) for _, lab, mk, mfc in VARIANTS], x=0.55,
                y=0.95, dy=0.075, length=0.05)
    ax.text(0.97, 0.06, "channel-level fits", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5,
            color=fs.INK2)
    return pd.DataFrame(out)


def dynamic_late(ax, L):
    late = L[L.w0.isin(LATE)].groupby(["penetration", "source", "c"]).L.mean().unstack("source")
    rng = np.random.default_rng(0)
    out = []
    for k, c in enumerate(DIOPTIC):
        g = late.xs(c, level="c")
        for j, m in enumerate(DYNAMIC):
            d = (g[m] - g["data"]).dropna()
            boot = [np.median(rng.choice(d.values, len(d))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            x = k + (j - 1) * 0.24
            col = DYNAMIC[m][1]
            ax.errorbar([x], [d.median()], [[d.median() - lo], [hi - d.median()]], color=col, marker="o",
                        mfc=col if c != 0.11 else "white", ms=3.2, elinewidth=0.8, capsize=0, lw=0)
            out.append(dict(c=c, model=m, n=len(d), median=d.median(), lo=lo, hi=hi,
                            inside=(lo > -MARGIN) and (hi < MARGIN)))
    band_axes(ax, (-0.8, 0.8))
    ax.text(0, -0.55, "noise-\nlimited", ha="center", va="bottom", fontsize=5.5, color=fs.INK2, linespacing=1.0)
    fs.line_key(ax, [(lab, col, "none", "o") for lab, col in DYNAMIC.values()], x=0.66, y=0.95, dy=0.075,
                length=0.05)
    ax.text(0.97, 0.06, "penetration-level fits", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5,
            color=fs.INK2)
    return pd.DataFrame(out)


def recovery(ax, B4, F):
    rng = np.random.default_rng(3)
    syn = B4.w_minus_a.to_numpy()
    ax.scatter(rng.uniform(-0.2, 0.2, len(syn)), syn, s=2, color=fs.MUTED, alpha=0.35, lw=0, zorder=2)
    ax.plot([-0.27, 0.27], [np.median(syn)] * 2, color=fs.INK, lw=1.4, zorder=4)
    d = F.full_p_w - F.full_p_a
    for mk, idx in F.groupby("monkey").groups.items():
        ax.scatter(1 + rng.uniform(-0.15, 0.15, len(idx)), d[idx], s=9, marker=fs.MONKEY_MARKER[mk],
                   facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4, zorder=3)
    ax.plot([0.73, 1.27], [d.median()] * 2, color=fs.INK, lw=1.4, zorder=4)
    fs.zero_line(ax, color=fs.INK2)
    set_med = B4.groupby("set").w_minus_a.median()
    lo, hi = np.percentile(set_med, [5, 95])
    ax.set_xticks([0, 1], ["synthetic,\nequal weights", "data"], fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.25, 0.82)
    ax.set_ylabel("pool weight − drive weight\nof the other eye")
    ax.text(0, ax.get_ylim()[1], f"median of 100 sets\n{np.median(set_med):.3f} [{lo:.3f}, {hi:.3f}]", ha="center",
            va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)
    ax.text(1, ax.get_ylim()[1], f"median {d.median():.2f}\n{int((d > 0).sum())}/{len(d)} > 0", ha="center",
            va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)
    return dict(data_median=float(d.median()), data_positive=int((d > 0).sum()), n=len(d),
                synthetic_set_median=float(np.median(set_med)), synthetic_5=float(lo), synthetic_95=float(hi),
                synthetic_fit_median_abs=float(np.median(np.abs(syn))))


def separate_delays(ax, B3):
    for mk, g in B3.groupby("monkey"):
        ax.scatter(g.sep_p_tau2, g.sep_p_tau2N, s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2,
                   edgecolor=fs.INK, lw=0.4, zorder=3)
    lo, hi = 8, 250
    ax.plot([lo, hi], [lo, hi], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ticks = [10, 20, 50, 100, 200]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter([str(t) for t in ticks]))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_aspect("equal")
    ax.set_xlabel("pool delay, dominant-eye term (ms)")
    ax.set_ylabel("pool delay, other-eye term (ms)")
    dt = B3.sep_p_tau2N - B3.sep_p_tau2
    p_dt = stats.wilcoxon(dt, alternative="greater").pvalue
    held = B3[~B3.empty_half]
    r = np.log(held.cv_full / held.cv_sep)
    p_r = stats.wilcoxon(r, alternative="greater").pvalue
    ax.text(0.04, 0.80, f"other-eye term slower by\n{dt.median():.1f} ms (median), {int((dt > 0).sum())}/{len(dt)}",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)
    ax.text(0.96, 0.04, f"held out, separate delays\nbetter in {int((r > 0).sum())}/{len(r)}, {fs.fmt_p(p_r)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
    return dict(n=len(B3), median_dtau=float(dt.median()), dtau_positive=int((dt > 0).sum()), p_dtau=p_dt,
                median_tau2=float(B3.sep_p_tau2.median()), median_tau2N=float(B3.sep_p_tau2N.median()),
                held_n=len(r), held_better=int((r > 0).sum()), held_p=p_r,
                held_median_pct=float(100 * (1 - np.exp(-r.median()))))


def main():
    fs.apply()
    b = out_root() / "binocular"
    C = pd.read_parquet(b / "fig5_channel_test.parquet")
    I = pd.read_parquet(b / "fig5_interaction_test.parquet")
    LI = pd.read_parquet(b / "fig5_interaction_test_L.parquet")
    LD = pd.concat([pd.read_parquet(b / "fig5_cond_L.parquet"),
                    pd.read_parquet(b / "fig5_org_L.parquet").query("source == 'org'")])
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    B3 = pd.read_parquet(b / "v3_B3.parquet")
    B4 = pd.read_parquet(b / "v3_B4_fits.parquet")

    fig = fs.figure("double", 108)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1.05], hspace=0.05)
    ax_a, ax_b, ax_c = top.subplots(1, 3, gridspec_kw=dict(width_ratios=[0.75, 1.2, 1.2]))
    ax_d, ax_e, ax_x = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
    ax_x.axis("off")
    T = two_stage(ax_a, C, I)
    K = capture(ax_b, LI)
    D = dynamic_late(ax_c, LD)
    R = recovery(ax_d, B4, F)
    S = separate_delays(ax_e, B3)
    for k, (mk, lab) in enumerate((("o", "monkey B"), ("s", "monkey J"))):
        ax_x.scatter([0.05], [0.9 - k * 0.08], s=9, marker=mk, facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4,
                     transform=ax_x.transAxes)
        ax_x.text(0.1, 0.9 - k * 0.08, lab, va="center", fontsize=5.5, color=fs.INK, transform=ax_x.transAxes)
    for ax, letter in zip((ax_a, ax_b, ax_c, ax_d, ax_e), "abcde"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS5_model_checks"))
    print("(a) two-stage gain:")
    print(T.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("(b) channel-level late index, model - data (90% CI):")
    print(K.round(3).to_string(index=False))
    print("(c) penetration-level late index, model - data (90% CI):")
    print(D.round(3).to_string(index=False))
    print(f"(d) parameter recovery: {R}")
    print(f"(e) separate delays: {S}")


if __name__ == "__main__":
    main()
