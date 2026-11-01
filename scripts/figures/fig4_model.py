"""Figure 4: binocular normalization is delayed, weights the other eye more than it drives the response, and acts
like a conductance (docs/09, figure set for V3; replaces the former fig5_models.py, whose rise-time and late-index
panels moved to figS6_dynamics.py and figS5_model_checks.py).

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_model.py
Normalization models fitted per penetration to the channel-averaged 6 x 6 time courses (0-219 ms after stimulus
onset; v1lb.dynorm; dominant eye on axis 0, each channel scaled by its dominant-eye peak at 90%):
* static: R = g E^n / (sigma^n + P^n) with drive E = (c_D + a c_N) u(t) and pool P = (c_D + w c_N) u(t), where a is the
  other eye's drive weight and w its pool weight;
* fixed-delay ('full'): the same, with the pool low-passed by tau2;
* conductance ('cond'): R relaxes toward the fixed-delay value with time constant tau_m / (1 + (P / sigma)^n), so the
  pool both divides the response and speeds it up; tau_m -> 0 gives the fixed-delay model;
* recurrent ('org'): binocular ORGaNICs, whose normalization signal is built recurrently from the responses.
(a) the model in words, and the colour key;
(b) data and conductance fit at 45% in the format of Fig 2a, on the stimulus clock: dominant eye at 90%, each eye
    alone at 45% and both eyes at 45%; data = mean over penetrations of the channel-averaged responses with 95%
    bootstrap CI (v3_A1_population.npz, all presentations), fit = mean over penetrations of each penetration's
    all-presentation conductance fit (fig5_cond.parquet, predicted here with dynorm.predict);
(c) linearity index over time (channel-average readout fig5_model.windowed_L, windows 0-10 ... 100-150 ms after the
    response onset, plotted at their centres) of the data and of the static and delayed (fixed-delay) fits; median
    over penetrations, data with its 95% bootstrap CI (v3_A1_L.parquet from v3_null_model_index.py); CIs beyond the
    axis are cut at its edge;
(d) held-out (odd/even) error of the binocular interaction R(c_D, c_N) - R(c_D, 0) with the monocular time courses
    (interaction loss), static vs delayed, each model fitted with that loss (v3_B2.parquet from
    v3_interocular_delay.py); one-sided Wilcoxon on log error ratios;
(e) pool weight vs drive weight of the other eye in the conductance fit (fig5_cond.parquet); two-sided Wilcoxon;
(f) % less held-out error per penetration for three pairings of the dynamic models (least-squares fits, odd/even
    cross-validation; fig5_cond.parquet, fig5_org.parquet); one-sided Wilcoxon on log error ratios; values beyond the
    axis are drawn at its edge (open triangles).
Held-out comparisons use 24 penetrations: one penetration (211008, early B) has no condition with two presentations
in its even half, so its even-half fits never see data (docs/09).
Inputs: <out_root>/binocular/v3_A1_population.npz, v3_A1_L.parquet, v3_B2.parquet, v3_B3.parquet (empty-half flag),
fig5_cond.parquet, fig5_org.parquet, fig5_model.parquet.
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
DIOPTIC = [0.11, 0.225, 0.45]
WINDOWS = {0: 10, 10: 20, 20: 35, 35: 50, 50: 75, 75: 100, 100: 150}
GAIN_LIM = (-15.0, 25.0)               # % axis limits in (f)
L_LIM = (-1.3, 2.0)                    # linearity-index axis limits in (c)
COND_KEYS = ("g", "sigma", "n", "delta", "tau1", "tau2", "a", "w", "tau_m")
MODELS = {"static": ("static", fs.MODEL["static"], (0, (3, 1.5))),
          "full": ("fixed-delay", fs.MODEL["fixed-delay"], "-"),
          "cond": ("conductance", fs.MODEL["conductance"], "-"),
          "org": ("recurrent", fs.MODEL["recurrent"], "-")}


def schematic(ax, tau2):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec=fs.INK2, lw=0.6)
    arrow = dict(arrowstyle="-|>", lw=0.6, color=fs.INK2, mutation_scale=6, shrinkA=1, shrinkB=1)

    def node(x, y, s, **kw):
        return ax.text(x, y, s, ha="center", va="center", fontsize=6, color=fs.INK, bbox=kw.pop("bbox", box),
                       linespacing=1.05, **kw)

    def arr(src, dst):
        """Arrow from a node centre (or a point) to a node centre, clipped at the boxes."""
        is_node = hasattr(src, "get_bbox_patch")
        ax.annotate("", xy=dst.get_position(), xytext=src.get_position() if is_node else src,
                    arrowprops=dict(arrow, patchA=src.get_bbox_patch() if is_node else None,
                                    patchB=dst.get_bbox_patch()))

    yd, yp, bus = 0.86, 0.58, 0.34
    for y, lab, col in ((0.77, "dominant eye", fs.COND["dominant"]), (0.67, "other eye", fs.COND["other"])):
        ax.text(0.0, y + 0.012, lab, ha="left", va="bottom", fontsize=6, color=fs.INK)
        ax.plot([0.0, bus], [y, y], color=col, lw=1.4, solid_capstyle="butt")
    ax.plot([bus, bus], [yp, yd], color=fs.INK2, lw=0.6)
    drive = node(0.47, yd, "drive")
    pool = node(0.47, yp, "pool")
    delay = node(0.73, yp, f"delay\n{tau2[0]:.0f}–{tau2[1]:.0f} ms")
    resp = node(0.87, 0.80, "response\n= drive$^n$ /\n($\\sigma^n$ + pool$^n$)")
    arr((bus, yd), drive)
    arr((bus, yp), pool)
    arr(drive, resp)
    arr(pool, delay)
    arr(delay, resp)
    ax.text(0.37, yd + 0.075, "dominant + drive weight × other", fontsize=5.3, color=fs.INK2, ha="left", va="bottom")
    ax.text(0.37, yp - 0.095, "dominant + pool weight × other", fontsize=5.3, color=fs.INK2, ha="left", va="top")
    ax.text(0.0, 0.335, "conductance: the delayed pool also shortens the\nresponse's time constant",
            fontsize=5.5, color=fs.INK, ha="left", va="center", linespacing=1.0)
    items = [("data", fs.MODEL["data"], "-"),
             ("static: pool without delay", *MODELS["static"][1:]),
             ("fixed-delay", *MODELS["full"][1:]),
             ("conductance", *MODELS["cond"][1:]),
             ("recurrent: pool built from the responses", *MODELS["org"][1:])]
    fs.line_key(ax, items, x=0.0, y=0.215, dy=0.055, length=0.07)


def fit_vs_data(ax, Y, pens, F):
    """Fig 2a format: data (solid, CI) and the conductance fit (dashed) for four conditions at 45%."""
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
    fs.line_key(ax, [("data (95% CI)", fs.INK2, "-"), ("conductance fit", fs.INK2, (0, (2.2, 1.2)))], x=0.40,
                y=0.62, dy=0.075, length=0.08)
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


def index_time_course(axs, L):
    xs = {w0: (w0 + w1) / 2 for w0, w1 in WINDOWS.items()}
    rows, clipped = [], 0
    for ax, c in zip(axs, DIOPTIC):
        g = L[np.isclose(L.c, c)]
        ax.axhline(1, color=fs.COND["linear"], lw=0.6, ls=(0, (3, 1.5)), zorder=0)
        fs.zero_line(ax)
        for src in ("static", "full"):
            m = g[g.source == src].groupby("w0").L.median()
            lab, col, ls = MODELS[src]
            ax.plot([xs[w] for w in m.index], m.values, color=col, ls=ls, lw=1.0, zorder=3)
            rows += [dict(c=c, source=src, w0=w, median=v) for w, v in m.items()]
        d = g[g.source == "data"].pivot(index="penetration", columns="w0", values="L")
        med, lo, hi = fs.boot_ci(d.to_numpy(), stat=np.nanmedian)
        x = np.array([xs[w] for w in d.columns])
        clipped += int((lo < L_LIM[0]).sum() + (hi > L_LIM[1]).sum())
        ax.errorbar(x, med, [med - np.maximum(lo, L_LIM[0]), np.minimum(hi, L_LIM[1]) - med], color=fs.MODEL["data"],
                    marker="o", ms=2.6, lw=1.1, elinewidth=0.6, capsize=0, zorder=4)
        rows += [dict(c=c, source="data", w0=w, median=v, lo=a, hi=b, n=int(d[w].notna().sum()))
                 for w, v, a, b in zip(d.columns, med, lo, hi)]
        ax.set_xlim(0, 150)
        ax.set_ylim(*L_LIM)
        ax.set_xticks([0, 50, 100])
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=5.8, loc="center", pad=2)
    for ax in axs[1:]:
        ax.tick_params(labelleft=False)
    axs[0].set_ylabel("linearity index")
    axs[1].set_xlabel("time from response onset (ms, window centre)")
    # direct labels at the right end of the last panel (45%)
    last = {s: float(g) for s, g in pd.DataFrame(rows).query("c == 0.45 and w0 == 100").set_index("source")
            ["median"].items()}
    for y, lab, col in ((1.0, "linear sum", fs.INK2), (last["static"], "static", MODELS["static"][1]),
                        (0.17, "dominant\neye alone", fs.INK2), (last["full"] - 0.1, "fixed-delay", fs.INK)):
        axs[2].text(153, y, lab, ha="left", va="center", fontsize=5.3, color=col, linespacing=0.95, clip_on=False)
    return pd.DataFrame(rows), clipped


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
    ax.set_xlabel("held-out interaction error, static")
    ax.set_ylabel("held-out interaction error, delayed")
    r = np.log(B.cv_static / B.cv_full)
    better = int((r > 0).sum())
    p = stats.wilcoxon(r, alternative="greater").pvalue
    pct = 100 * (1 - np.exp(-r.median()))
    ax.text(0.97, 0.05, f"delayed better\n{better}/{len(B)}, {fs.fmt_p(p)}\nmedian {pct:.0f}% less error",
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
    ax.set_xlabel("drive weight of the other eye")
    ax.set_ylabel("pool weight of the other eye")
    n = int((w > a).sum())
    p = stats.wilcoxon(w, a).pvalue
    ax.text(0.04, 0.96, f"pool > drive\n{n}/{len(F)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.5, color=fs.INK, linespacing=1.0)
    ax.text(0.96, 0.04, "conductance fits", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5,
            color=fs.INK2)
    return dict(n=len(F), pool_gt_drive=n, p=p, median_w=float(w.median()), median_a=float(a.median()),
                median_w_minus_a=float((w - a).median()))


def held_out(ax, F):
    pairs = (("cond", "full"), ("cond", "org"), ("org", "full"))
    rows = []
    for x, (ch, ref) in enumerate(pairs):
        r = np.log(F[f"cv_{ref}"] / F[f"cv_{ch}"])
        pct = 100 * (1 - np.exp(-r))
        rng = np.random.default_rng(x)
        for mk, idx in F.groupby("monkey").groups.items():
            v = pct[idx].to_numpy()
            xs = x + rng.uniform(-0.15, 0.15, len(v))
            inside = (v >= GAIN_LIM[0]) & (v <= GAIN_LIM[1])
            ax.scatter(xs[inside], v[inside], s=8, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                       lw=0.4, zorder=3)
            for edge, mkr in ((GAIN_LIM[1], "^"), (GAIN_LIM[0], "v")):
                out = (v > edge) if edge > 0 else (v < edge)
                ax.scatter(xs[out], np.full(out.sum(), edge), s=10, marker=mkr, facecolor="white", edgecolor=fs.INK,
                           lw=0.5, zorder=3, clip_on=False)
        ax.plot([x - 0.25, x + 0.25], [pct.median()] * 2, color=MODELS[ch][1], lw=1.6, zorder=4)
        p = stats.wilcoxon(r, alternative="greater").pvalue
        ax.text(x, GAIN_LIM[1] + 1.5, f"{int((r > 0).sum())}/{len(r)}\n{fs.fmt_p(p)}", ha="center", va="bottom",
                fontsize=5.5, color=fs.INK, linespacing=1.0)
        by = {mk: f"{int((r[F.monkey == mk] > 0).sum())}/{int((F.monkey == mk).sum())}, median "
                  f"{100 * (1 - np.exp(-r[F.monkey == mk].median())):.1f}%" for mk in ("B", "J")}
        rows.append(dict(pair=f"{MODELS[ch][0]} vs {MODELS[ref][0]}", n=len(r), better=int((r > 0).sum()),
                         median_pct=pct.median(), median_log_ratio=r.median(), p=p, B=by["B"], J=by["J"],
                         clipped=int(((pct < GAIN_LIM[0]) | (pct > GAIN_LIM[1])).sum())))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks(range(3), ["conductance\nvs\nfixed-delay", "conductance\nvs\nrecurrent", "recurrent\nvs\nfixed-delay"],
                  fontsize=5.5)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(*GAIN_LIM)
    ax.set_ylabel("% less held-out error\n(first model vs second)")
    return pd.DataFrame(rows)


def main():
    fs.apply()
    b = out_root() / "binocular"
    z = np.load(b / "v3_A1_population.npz")
    pens, Y = [str(p) for p in z["pens"]], z["Y"]
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    O = pd.read_parquet(b / "fig5_org.parquet").set_index("penetration")
    F = F.join(O[[c for c in O.columns if c.startswith(("cv_org", "org_"))]])
    L = pd.read_parquet(b / "v3_A1_L.parquet")
    B2 = pd.read_parquet(b / "v3_B2.parquet")
    empty = pd.read_parquet(b / "v3_B3.parquet").empty_half
    M = pd.read_parquet(b / "fig5_model.parquet").set_index("penetration")
    held = F.drop(index=empty[empty].index)

    fig = fs.figure("double", 118)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1.02], hspace=0.05)
    gt = top.add_gridspec(1, 5, width_ratios=[1.55, 1.15, 0.62, 0.62, 0.62], wspace=0.06)
    ax_a, ax_b = top.add_subplot(gt[0]), top.add_subplot(gt[1])
    ax_c = [top.add_subplot(gt[2])]
    ax_c += [top.add_subplot(gt[k], sharey=ax_c[0]) for k in (3, 4)]
    ax_d, ax_e, ax_f = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[1, 1, 1.1]))

    schematic(ax_a, (float(F.full_p_tau2.median()), float(F.cond_p_tau2.median())))
    fit = fit_vs_data(ax_b, Y, pens, F)
    tc, n_clip = index_time_course(ax_c, L)
    inter = interaction(ax_d, B2)
    wts = weights(ax_e, F)
    H = held_out(ax_f, held)
    for ax, letter in zip((ax_a, ax_b, ax_c[0], ax_d, ax_e, ax_f), "abcdef"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig4_model"))

    print(f"(a) median fitted pool delay tau2 (conductance) {F.cond_p_tau2.median():.1f} ms, fixed-delay "
          f"{F.full_p_tau2.median():.1f} ms; tau_m median {F.cond_p_tau_m.median():.2f} ms; n = {len(F)}")
    print("(b) data vs conductance fit (means over penetrations; late = 120-200 ms after stimulus onset):")
    print(fit.round(3).to_string())
    print(f"(c) linearity index medians (clipped CI ends: {n_clip}):")
    print(tc.pivot_table(index=["c", "w0"], columns="source", values="median").round(2).to_string())
    print(tc[tc.source == "data"][["c", "w0", "median", "lo", "hi", "n"]].round(2).to_string(index=False))
    print(f"(d) interaction loss, delayed vs static: {inter}")
    r_raw = np.log(M.drop(index=empty[empty].index).cv_static / M.drop(index=empty[empty].index).cv_full)
    print(f"    (not plotted) total held-out error, delayed vs static: {int((r_raw > 0).sum())}/{len(r_raw)}, "
          f"one-sided p {stats.wilcoxon(r_raw, alternative='greater').pvalue:.2g}, median "
          f"{100 * (1 - np.exp(-r_raw.median())):.1f}% less")
    print(f"(e) conductance weights: {wts}")
    print(f"    by monkey: " + "; ".join(f"{mk} {int((g.cond_p_w > g.cond_p_a).sum())}/{len(g)}"
                                       for mk, g in F.groupby("monkey")))
    print(f"(f) held-out comparisons (n = {len(held)}; excluded: {list(empty[empty].index)}):")
    print(H.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
