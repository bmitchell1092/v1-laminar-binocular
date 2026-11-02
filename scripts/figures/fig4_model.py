"""Figure 4: binocular normalization is delayed and weights the other eye more than it drives the response (docs/09,
figure set for V3; rebuilt 2026-11-02 for readers outside the field, docs/11).

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_model.py
Normalization models fitted per penetration to the channel-averaged 6 x 6 time courses (0-219 ms after stimulus
onset; v1lb.dynorm; dominant eye on axis 0, each channel scaled by its dominant-eye peak at 90%):
* static: R = g E^n / (sigma^n + P^n) with drive E = (c_D + a c_N) u(t) and pool P = (c_D + w c_N) u(t), where a is the
  other eye's drive weight and w its pool weight;
* fixed-delay ('full'): the same, with the pool low-passed by tau2 (the "delayed" model of panels c and d);
* conductance ('cond'): R relaxes toward the fixed-delay value with time constant tau_m / (1 + (P / sigma)^n), so the
  pool both divides the response and speeds it up; tau_m -> 0 gives the fixed-delay model.
(a) illustration of delayed normalization (dynorm.step_responses; median fitted drive weight a, pool weight w and
    pool time constant tau2 of the fixed-delay fits; n = 2, sigma = 0.25 and tau1 = 5 ms hand-set, delta = 0):
    the drive (both eyes add at once), the pool (both eyes, lagging) and the response drive^n / (sigma^n + pool^n),
    for both eyes at 45% each and the dominant eye alone at 45%;
(b) data and conductance fit at 45% in the format of the former Fig 2a, on the stimulus clock: dominant eye at 90%,
    each eye alone at 45% and both eyes at 45% (the text's comparison with 90% rests on this panel; every data trace
    is solid and every fit dashed, and the key says so); data = mean over penetrations of the channel-averaged
    responses with 95% bootstrap CI (v3_A1_population.npz, all presentations), fit = mean over penetrations of each
    penetration's all-presentation conductance fit (fig5_cond.parquet, predicted here with dynorm.predict);
(c) what the other eye adds (both eyes minus the dominant eye alone) at 11 and 45%, data (mean, 95% bootstrap CI,
    5 ms running mean for display) against the static fit (v3_B1_params.parquet, static_p_*) and the fixed-delay
    fit (fig5_cond.parquet, full_p_*), each predicted per penetration and averaged; window means printed at 11, 22.5
    and 45% (replaces the linearity index over time, whose model readouts are in v3_B1.csv);
(d) held-out (odd/even) error of the binocular interaction R(c_D, c_N) - R(c_D, 0) with the monocular time courses
    (interaction loss), static vs delayed, each model fitted with that loss (v3_B2.parquet from
    v3_interocular_delay.py); one-sided Wilcoxon on log error ratios;
(e) pool weight vs drive weight of the other eye in the conductance fit (fig5_cond.parquet); two-sided Wilcoxon.
The held-out comparison of the three dynamic models (former panel f) is Supplementary Fig S5f.
Held-out comparisons use 24 penetrations: one penetration (211008, early B) has no condition with two presentations
in its even half, so its even-half fits never see data (docs/09).
Inputs: <out_root>/binocular/v3_A1_population.npz, v3_B1_params.parquet, v3_B2.parquet, v3_B3.parquet (empty-half
flag), fig5_cond.parquet, fig5_model.parquet.
Writes <out_root>/figures/fig4_model.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
T = np.arange(0, 220)
SHOWN = [0.11, 0.45]                    # (c) contrasts drawn; 22.5% printed
PRINTED = [0.11, 0.225, 0.45]
C_WINDOWS = {"30-60 ms": (30, 60), "70-140 ms": (70, 140), "150-219 ms": (150, 220)}     # stimulus clock
COND_KEYS = ("g", "sigma", "n", "delta", "tau1", "tau2", "a", "w", "tau_m")
FULL_KEYS = ("g", "sigma", "n", "delta", "tau1", "tau2", "a", "w", "beta")
MODELS = {"static": ("no delay (static)", fs.MODEL["static"], (0, (3, 1.5))),
          "full": ("delayed pool (fixed delay)", fs.MODEL["fixed-delay"], "-")}
H_MM = 122
POS = {"a1": (12, 104, 62, 9), "a2": (12, 89, 62, 9), "a3": (12, 68, 62, 15), "b": (98, 70, 50, 43),
       "c1": (14, 12, 30, 40), "c2": (48, 12, 30, 40), "d": (96, 12, 36, 40), "e": (142, 12, 30, 40)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def timing_sketch(axs, a, w, tau2):
    """Drive, pool and response of delayed normalization for both eyes at 45% and the dominant eye alone."""
    t = np.arange(-20, 201, 1.0)
    u, v = dynorm.step_responses(t, 0.0, 5.0, tau2)
    c, n, sig = 0.45, 2.0, 0.25
    drive = {"both": (c + a * c) * u, "dom": c * u}
    pool = {"both": (c + w * c) * v, "dom": c * v}
    resp = {k: drive[k] ** n / (sig ** n + pool[k] ** n) for k in drive}
    rows = (("drive: both eyes add, at once", drive),
            (f"pool: both eyes, lagging (time constant {tau2:.0f} ms)", pool),
            ("response = drive² / (σ² + pool²)", resp))
    for ax, (title, d) in zip(axs, rows):
        ax.plot(t, d["dom"], color=fs.COND["dominant"], lw=1.0)
        ax.plot(t, d["both"], color=fs.COND["both"], lw=1.2)
        ax.set_xlim(-20, 200)
        ax.set_ylim(0, max(d["both"].max(), d["dom"].max()) * 1.15)
        ax.set_yticks([])
        ax.set_title(title, fontsize=6, loc="left", pad=2)
        ax.tick_params(labelbottom=False)
    axs[2].tick_params(labelbottom=True)
    axs[2].set_xticks([0, 50, 100, 150, 200])
    axs[2].set_xlabel("time from response onset (ms)")
    fs.line_key(axs[2], [("both eyes, 45% each", fs.COND["both"], "-"), ("dominant eye alone", fs.COND["dominant"], "-")],
                x=0.62, y=0.9, dy=0.2, length=0.07)
    r = resp["both"]
    k_pk = int(np.argmax(r))
    axs[2].annotate("eyes add", xy=(t[k_pk], r[k_pk]), xytext=(t[k_pk] + 22, r[k_pk] * 0.98), fontsize=5.5,
                    color=fs.INK, va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=fs.INK2))
    k = int(np.searchsorted(t, 140))
    axs[2].annotate("pool arrives: both eyes ≈ dominant eye", xy=(140, r[k]), xytext=(200, r.max() * 0.33),
                    fontsize=5.5, color=fs.INK, ha="right", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=fs.INK2))
    axs[0].text(1.0, 1.32, "illustration: median fitted weights and delay", transform=axs[0].transAxes, ha="right",
                va="bottom", fontsize=5.5, color=fs.INK2)
    late = t >= 120
    return dict(peak_both_over_dom=float(r.max() / resp["dom"].max()),
                late_both_over_dom=float(r[late].mean() / resp["dom"][late].mean()))


def fit_vs_data(ax, Y, pens, F):
    """Data (solid, CI) and the conductance fit (dashed) for four conditions at 45%."""
    i, j = list(CONTRASTS).index(0.45), list(CONTRASTS).index(0.9)
    conds = [("dominant eye, 90%", fs.COND["linear"], (j, 0)), ("dominant eye, 45%", fs.COND["dominant"], (i, 0)),
             ("other eye, 45%", fs.COND["other"], (0, i)), ("both eyes, 45% each", fs.COND["both"], (i, i))]
    cd = np.array([CONTRASTS[c[2][0]] for c in conds])
    cn = np.array([CONTRASTS[c[2][1]] for c in conds])
    P = np.array([dynorm.predict({k: F.loc[p, f"cond_p_{k}"] for k in COND_KEYS}, T.astype(float), cd, cn)
                  for p in pens])                                               # [pen, condition, time]
    for k, (lab, col, (a, b)) in enumerate(conds):
        fs.ci_line(ax, T, Y[:, a, b], col, lw=1.2 if col == fs.INK else 1.0, alpha=0.15)
        ax.plot(T, P[:, k].mean(axis=0), color=col, ls=(0, (2.2, 1.2)), lw=0.9, zorder=4)
    fs.zero_line(ax)
    ax.set_xlim(0, 219)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    fs.line_key(ax, [(lab, col, "-") for lab, col, _ in conds], x=0.40, y=0.97, dy=0.075, length=0.08)
    fs.line_key(ax, [("solid: data (95% CI)", fs.INK, "-"), ("dashed: model fit", fs.INK, (0, (2.2, 1.2)))],
                x=0.40, y=0.62, dy=0.075, length=0.08)
    ax.set_title("delayed normalization (conductance form) fits the data", fontsize=6.5)
    # numbers for the text: both eyes relative to the dominant eye at 90% (early) and 45% (late), data and fit
    D = {lab: Y[:, a, b].mean(axis=0) for lab, _, (a, b) in conds}
    M = {lab: P[:, k].mean(axis=0) for k, (lab, _, _) in enumerate(conds)}
    out = {}
    for src, X in (("data", D), ("fit", M)):
        pk_both = float(X["both eyes, 45% each"].max())
        pk90 = float(X["dominant eye, 90%"].max())
        late = (T >= 120) & (T < 200)
        out[src] = dict(peak_both_over_peak_dom90=pk_both / pk90,
                        late_both_over_dom45=X["both eyes, 45% each"][late].mean() / X["dominant eye, 45%"][late].mean(),
                        late_both_over_dom90=X["both eyes, 45% each"][late].mean() / X["dominant eye, 90%"][late].mean())
    return pd.DataFrame(out).T


def model_added(Y, pens, S, F):
    """What the other eye adds, data and the static and fixed-delay fits: dict c -> {source: [pen, time]}."""
    out = {}
    for c in PRINTED:
        i = list(CONTRASTS).index(c)
        cd, cn = np.array([c, c]), np.array([c, 0.0])
        d = {"data": Y[:, i, i] - Y[:, i, 0]}
        for src, tab, pre, keys, static in (("static", S, "static_p_", FULL_KEYS, True),
                                            ("full", F, "full_p_", FULL_KEYS, False)):
            P = np.array([dynorm.predict({k: tab.loc[p, f"{pre}{k}"] for k in keys}, T.astype(float), cd, cn,
                                         static=static) for p in pens])
            d[src] = P[:, 0] - P[:, 1]
        out[c] = d
    return out


def added_vs_models(axs, Ad):
    for k, (ax, c) in enumerate(zip(axs, SHOWN)):
        d = Ad[c]
        fs.ci_line(ax, T, fs.smooth(d["data"]), fs.MODEL["data"], lw=1.2)
        for src in ("static", "full"):
            lab, col, ls = MODELS[src]
            ax.plot(T, d[src].mean(axis=0), color=col, ls=ls, lw=1.0, zorder=3)
        fs.zero_line(ax)
        ax.set_xlim(0, 219)
        ax.set_ylim(-0.06, 0.32)
        ax.set_xticks([0, 100, 200])
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
        if k:
            ax.tick_params(labelleft=False)
    axs[0].set_ylabel("added by the other eye (MUA)\n= both eyes − dominant eye alone")
    fig = axs[0].figure
    x0 = axs[0].get_position().x0
    x1 = axs[1].get_position().x1
    fig.text((x0 + x1) / 2, axs[0].get_position().y0 - 0.042, "time from stimulus onset (ms)", ha="center",
             va="top", fontsize=7, color=fs.INK)
    fs.line_key(axs[1], [("data (95% CI)", fs.MODEL["data"], "-"), MODELS["static"], MODELS["full"]],
                x=0.25, y=0.95, dy=0.085, length=0.12)
    rows = []
    for c in PRINTED:
        for name, (w0, w1) in C_WINDOWS.items():
            s = (T >= w0) & (T < w1)
            rows.append(dict(c=c, window=name, **{src: float(v[:, s].mean()) for src, v in Ad[c].items()}))
    return pd.DataFrame(rows)


def interaction(ax, B):
    for mk, g in B.groupby("monkey"):
        ax.scatter(g.cv_static, g.cv_full, s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                   lw=0.4, zorder=3)
    lo = min(B.cv_full.min(), B.cv_static.min()) / 1.4
    hi = max(B.cv_full.max(), B.cv_static.max()) * 1.4
    ax.plot([lo, hi], [lo, hi], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ticks = [t for t in (0.5, 1, 2, 5, 10) if lo <= t <= hi]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter([f"{t:g}" for t in ticks]))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_aspect("equal")
    ax.set_xlabel("prediction error, no delay")
    ax.set_ylabel("prediction error, delayed pool")
    ax.set_title("held-out binocular interaction", fontsize=6.5)
    r = np.log(B.cv_static / B.cv_full)
    better = int((r > 0).sum())
    p = stats.wilcoxon(r, alternative="greater").pvalue
    pct = 100 * (1 - np.exp(-r.median()))
    ax.text(0.97, 0.05, f"delay better\n{better}/{len(B)}, {fs.fmt_p(p)}\nmedian {pct:.0f}% less error",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
    for k, (mk, lab) in enumerate((("o", "monkey B"), ("s", "monkey J"))):
        y = 0.93 - k * 0.08
        ax.scatter([0.07], [y], s=9, marker=mk, facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4, transform=ax.transAxes)
        ax.text(0.12, y, lab, va="center", fontsize=5.5, color=fs.INK, transform=ax.transAxes)
    by = {mk: (int((r[B.monkey == mk] > 0).sum()), int((B.monkey == mk).sum()), float(r[B.monkey == mk].median()))
          for mk in ("B", "J")}
    return dict(n=len(B), better=better, p=p, median_log_ratio=float(r.median()), median_pct_less=pct, by_monkey=by)


def weights(ax, F):
    a, w = F.cond_p_a, F.cond_p_w
    for mk, idx in F.groupby("monkey").groups.items():
        ax.scatter(a[idx], w[idx], s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.MODEL["conductance"],
                   edgecolor="white", lw=0.3, zorder=3)
    top = max(1.6, float(max(a.max(), w.max())) * 1.06)
    ax.plot([0, top], [0, top], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal")
    ax.set_xlabel("other eye's weight in the drive")
    ax.set_ylabel("other eye's weight in the pool")
    ax.set_title("the other eye normalizes\nmore than it drives", fontsize=6.5, linespacing=1.0)
    n = int((w > a).sum())
    p = stats.wilcoxon(w, a).pvalue
    ax.text(0.04, 0.96, f"pool > drive\n{n}/{len(F)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.5, color=fs.INK, linespacing=1.0)
    ax.text(0.96, 0.04, "conductance fits", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5,
            color=fs.INK2)
    return dict(n=len(F), pool_gt_drive=n, p=p, median_w=float(w.median()), median_a=float(a.median()),
                median_w_minus_a=float((w - a).median()))


def main():
    fs.apply()
    pd.set_option("display.width", 200)
    b = out_root() / "binocular"
    z = np.load(b / "v3_A1_population.npz")
    pens, Y = [str(p) for p in z["pens"]], z["Y"]
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    S = pd.read_parquet(b / "v3_B1_params.parquet").set_index("penetration")
    B2 = pd.read_parquet(b / "v3_B2.parquet")
    empty = pd.read_parquet(b / "v3_B3.parquet").empty_half
    M = pd.read_parquet(b / "fig5_model.parquet").set_index("penetration")

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")                     # axes placed in mm: in-panel keys stay inside their panels
    ax_a = [axes_mm(fig, *POS[n]) for n in ("a1", "a2", "a3")]
    ax_b, ax_c1, ax_c2, ax_d, ax_e = (axes_mm(fig, *POS[n]) for n in ("b", "c1", "c2", "d", "e"))
    sk = timing_sketch(ax_a, float(F.full_p_a.median()), float(F.full_p_w.median()), float(F.full_p_tau2.median()))
    fit = fit_vs_data(ax_b, Y, pens, F)
    Ad = model_added(Y, pens, S, F)
    Wc = added_vs_models([ax_c1, ax_c2], Ad)
    inter = interaction(ax_d, B2)
    wts = weights(ax_e, F)
    for ax, letter, y in ((ax_a[0], "a", 8), (ax_b, "b", 1.0), (ax_c1, "c", 1.0), (ax_d, "d", 1.0), (ax_e, "e", 1.0)):
        fs.panel_label(ax, letter, y=y)
    print(fs.save(fig, "fig4_model"))

    print(f"(a) sketch with median fixed-delay a {F.full_p_a.median():.2f}, w {F.full_p_w.median():.2f}, tau2 "
          f"{F.full_p_tau2.median():.1f} ms: {sk}")
    print(f"    median fitted pool delay tau2 (conductance) {F.cond_p_tau2.median():.1f} ms, fixed-delay "
          f"{F.full_p_tau2.median():.1f} ms; tau_m median {F.cond_p_tau_m.median():.2f} ms; n = {len(F)}")
    print("(b) data vs conductance fit (means over penetrations; late = 120-200 ms after stimulus onset):")
    print(fit.round(3).to_string())
    print("(c) what the other eye adds, window means (stimulus clock): data, static, fixed-delay:")
    print(Wc.round(4).to_string(index=False))
    print(f"(d) interaction loss, delayed vs static: {inter}")
    r_raw = np.log(M.drop(index=empty[empty].index).cv_static / M.drop(index=empty[empty].index).cv_full)
    print(f"    (not plotted) total held-out error, delayed vs static: {int((r_raw > 0).sum())}/{len(r_raw)}, "
          f"one-sided p {stats.wilcoxon(r_raw, alternative='greater').pvalue:.2g}, median "
          f"{100 * (1 - np.exp(-r_raw.median())):.1f}% less")
    print(f"(e) conductance weights: {wts}")
    print(f"    by monkey: " + "; ".join(f"{mk} {int((g.cond_p_w > g.cond_p_a).sum())}/{len(g)}"
                                       for mk, g in F.groupby("monkey")))


if __name__ == "__main__":
    main()
