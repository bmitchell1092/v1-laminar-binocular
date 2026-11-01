"""Figure 2: the two eyes add at response onset, then normalize within tens of milliseconds (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_add_then_normalize.py
(a) population responses at 45% aligned to each channel's onset: each eye alone at 45%, both eyes at 45%, and the
    dominant eye at 90% (about what linear summation of contrast predicts); penetration means, 95% bootstrap CI;
(b) how the equivalent contrast is read (population means, illustration), on one shared response axis: the
    binocular response at 45% is placed on the dominant eye's contrast response measured in the same window and
    read off as the dominant-eye contrast that gives the same response;
(c) iso-response contours across the dichoptic conditions (v3_iso_contours.py), for targets equal to the dominant
    eye's response at 22.5% and 45%, early (0-20 ms after each penetration's onset) and late (50-200 ms): for each
    other-eye contrast, the dominant-eye contrast at which the binocular response reaches the target, from the
    population-average window responses, with 95% penetration-bootstrap CIs. Axes are fractions of the target
    contrast, the other eye's contrast multiplied by the eye weight of the same window, so linear summation is the
    diagonal and "the dominant eye alone" the L-shaped rule; lines = fitted Minkowski contours; open circle on the
    horizontal axis = the other eye alone already reaches the target; open triangle = CI beyond the axis. The note
    gives the exponent early -> late for the population and the per-penetration medians (two-sided Wilcoxon);
(d) per penetration, early (0-20 ms after onset) vs late (50-200 ms) linearity index (mean of the windows' channel
    medians, as in equivalent_contrast.py), including 5.5%; Wilcoxon signed-rank tests; values beyond the axis are
    drawn at its edge as open triangles and counted. At 90% the equivalent contrast cannot be measured (the
    dominant eye has no higher contrast to match).
The linearity index over time (formerly panel c) is Fig 4c, with the model predictions.
Inputs: figure_data.py (population.npz), <out_root>/binocular/k_windows.parquet, v3_C_population.parquet,
v3_C_population_boot.parquet, v3_C_penetration.parquet.
Writes <out_root>/figures/fig2_add_then_normalize.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, invert_crf  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
PAIRED = [0.055, 0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts: 0-20 ms and 50-200 ms after onset
CONTOUR_TARGETS = (0.225, 0.45)
EPOCH = {"early": fs.MUTED, "late": fs.INK}           # docs/11: epoch grey ramp
XMAX, YMAX = 1.5, 1.6                                 # contour axes (fractions of the target contrast)


def timecourses(ax, pop):
    rel, O = pop["rel"], pop["onset"]
    i = list(CONTRASTS).index(0.45)
    sel = (rel >= -20) & (rel <= 200)
    items = [("dominant eye, 90%  (≈ linear sum)", fs.COND["linear"], (0, (3, 1.5)), O[:, 5, 0]),
             ("dominant eye, 45%", fs.COND["dominant"], "-", O[:, i, 0]),
             ("other eye, 45%", fs.COND["other"], "-", O[:, 0, i]),
             ("both eyes, 45% each", fs.COND["both"], "-", O[:, i, i])]
    for lab, col, ls, Y in items:
        fs.ci_line(ax, rel[sel], Y[:, sel], col, ls=ls, lw=1.2 if col == fs.INK else 1.0)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    fs.line_key(ax, [(lab, col, ls) for lab, col, ls, _ in items], x=0.42, y=0.93)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.1)
    ax.text(10, ax.get_ylim()[1] * 0.99, "first 20 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def equivalent_contrast_sketch(ax, pop, w, title, a_med, ymax):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= w[0]) & (rel < w[1])
    crf = O[:, :, 0, :][:, :, sel].mean(axis=2).mean(axis=0)              # dominant eye, population mean
    i = list(CONTRASTS).index(0.45)
    rb = O[:, i, i, :][:, sel].mean(axis=1).mean()
    kc, cens = invert_crf(crf, np.array([rb]))
    kc = float(kc[0]) if cens[0] == 0 else CONTRASTS[-1]
    c = CONTRASTS[1:]
    ax.plot(c, crf[1:], color=fs.COND["dominant"], marker="o", ms=2.6, lw=1.0)
    ax.plot([0.45], [rb], color=fs.INK, marker="D", ms=3.4, zorder=5)
    ax.plot([0.45, kc], [rb, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([kc, kc], [0, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([(1 + a_med) * 0.45], [0], marker="^", ms=3.8, color=fs.MUTED, clip_on=False, zorder=6)
    fs.log_contrast_axis(ax, (0.055, 0.11, 0.225, 0.45, 0.9))
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"], fontsize=5.5)
    ax.set_xlim(0.045, 1.05)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("dominant-eye contrast (%)")
    ax.set_title(title, fontsize=6.5, color=fs.INK)
    ax.text(0.05, 0.97, f"both eyes at 45%\n= dominant eye at {100 * kc:.0f}%", transform=ax.transAxes,
            ha="left", va="top", fontsize=6, color=fs.INK, linespacing=1.0)
    return kc / 0.45


def minkowski(y, m):
    """Dominant-eye contrast (fraction of target) on the Minkowski contour of exponent m, for weighted other-eye
    contrast y (fraction of target) in [0, 1]."""
    return (1 - np.clip(y, 0, 1) ** m) ** (1 / m)


def iso_contours(axs, Pop, Boot, Pen):
    """Population iso-response contours (v3_iso_contours.py) at 22.5% and 45%, early vs late, with 95%
    penetration-bootstrap CIs of the dominant-eye contrast needed; per-penetration exponents in the note."""
    out = []
    for ax, cs in zip(axs, CONTOUR_TARGETS):
        yy = np.linspace(0, 1, 101)
        ax.plot([0, 1], [1, 0], color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=1)
        ax.plot([0, 1, 1], [1, 1, 0], color=fs.RULE, lw=0.9, zorder=0)
        for k, w in enumerate(("early", "late")):
            col = EPOCH[w]
            p = Pop[(Pop.window == w) & (Pop.cstar == cs)].sort_values("c_N")
            b = Boot[(Boot.window == w) & (Boot.cstar == cs)].groupby("c_N").x
            lo = b.quantile(0.025).reindex(p.c_N).to_numpy()
            hi = b.quantile(0.975).reindex(p.c_N).to_numpy()
            m = float(p.m.iloc[0])
            ax.plot(yy, minkowski(yy, m), color=col, lw=0.8, zorder=2)
            show = p.y.to_numpy() <= XMAX
            xo, yd = p.y.to_numpy()[show] + (k - 0.5) * 0.025, p.x.to_numpy()[show]
            lo, hi = np.minimum(lo[show], YMAX), np.minimum(hi[show], YMAX)
            bound = p.status.to_numpy()[show] != "ok"
            ax.vlines(xo, lo, hi, color=col, lw=0.6, zorder=3)
            ax.scatter(xo[~bound], yd[~bound], s=10, color=col, zorder=4, lw=0)
            ax.scatter(xo[bound], np.minimum(yd[bound], YMAX), s=10, facecolor="white", edgecolor=col, lw=0.7,
                       zorder=4, clip_on=False)
            ax.scatter(xo[hi >= YMAX], np.full((hi >= YMAX).sum(), YMAX), marker="^", s=9, facecolor="white",
                       edgecolor=col, lw=0.6, zorder=4, clip_on=False)
            r = p.iloc[0]
            out.append(dict(target=cs, window=w, eye_weight=r.a, m=r.m, m_lo=r.m_lo, m_hi=r.m_hi,
                            suppression=r.suppression, n_bound=int(bound.sum()), n_ci_at_top=int((hi >= YMAX).sum())))
        g = Pen[Pen.cstar == cs].pivot(index="penetration", columns="window", values="m")[["early", "late"]].dropna()
        pw = stats.wilcoxon(g.late, g.early).pvalue
        out[-1].update(pen_n=len(g), pen_early=g.early.median(), pen_late=g.late.median(),
                       pen_k=int((g.late > g.early).sum()), pen_p=pw)
        pe, pl = out[-2], out[-1]
        ax.text(0.0, 1.03, f"exponent, early → late:\npopulation {pe['m']:.1f} → {pl['m']:.1f}; penetrations "
                f"{g.early.median():.1f} → {g.late.median():.1f}\n({int((g.late > g.early).sum())}/{len(g)} larger, "
                f"{fs.fmt_p(pw)})", transform=ax.transAxes, fontsize=5, color=fs.INK2, va="bottom", ha="left",
                linespacing=1.05, in_layout=False)
        ax.set_xlim(0, XMAX + 0.05)
        ax.set_ylim(0, YMAX)
        ax.set_xticks([0, 0.5, 1.0, 1.5])
        ax.set_yticks([0, 0.5, 1.0, 1.5])
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"target: dominant eye at {fs.contrast_label(cs)}", fontsize=6.5, pad=22)
        ax.text(0.27, 0.46, "linear sum", rotation=-45, fontsize=5, color=fs.INK2, ha="center", va="center",
                rotation_mode="anchor", in_layout=False)
        ax.text(1.045, 0.45, "dominant eye alone", rotation=90, fontsize=5, color=fs.INK2, ha="left",
                va="center", in_layout=False)
    axs[0].set_ylabel("dominant-eye contrast\n(fraction of target)")
    axs[1].tick_params(labelleft=False)
    fs.line_key(axs[1], [("0–20 ms", EPOCH["early"], "-", "o"), ("50–200 ms", EPOCH["late"], "-", "o")],
                x=0.55, y=0.95, dy=0.085, length=0.1, fontsize=5)
    return pd.DataFrame(out)


def early_late(ax, W):
    out = []
    lo_clip, hi_clip = -1.5, 3.0
    gap = 2.6
    for g, c in enumerate(PAIRED):
        per = W[W.c == c].groupby(["penetration", "w0"]).L.median().unstack()
        e, l_ = per[EARLY].mean(axis=1), per[LATE].mean(axis=1)
        ok = e.notna() & l_.notna()
        e, l_ = e[ok], l_[ok]
        p = stats.wilcoxon(e, l_).pvalue
        x0, x1 = g * gap, g * gap + 1
        for a, b in zip(e.clip(lo_clip, hi_clip), l_.clip(lo_clip, hi_clip)):
            ax.plot([x0, x1], [a, b], color=fs.RULE, lw=0.5, zorder=1)
        for xx, v in ((x0, e.to_numpy()), (x1, l_.to_numpy())):
            for sel, mk, yv in ((v > hi_clip, "^", hi_clip), (v < lo_clip, "v", lo_clip)):
                ax.scatter(np.full(sel.sum(), xx), np.full(sel.sum(), yv), marker=mk, s=8, facecolor="white",
                           edgecolor=fs.INK2, lw=0.5, zorder=2, clip_on=False)
        n_clip = int(((e > hi_clip) | (e < lo_clip)).sum() + ((l_ > hi_clip) | (l_ < lo_clip)).sum())
        ax.plot([x0, x1], [e.median(), l_.median()], color=fs.CONTRAST[c], lw=1.6, marker="o", ms=3.2, zorder=3)
        ax.text((x0 + x1) / 2, hi_clip + 0.2, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        ax.text((x0 + x1) / 2, lo_clip - 0.55, fs.contrast_label(c), ha="center", va="top", fontsize=6.5,
                color=fs.INK)
        ax.text((x0 + x1) / 2, -0.21, f"{n_clip} off axis" if n_clip else "", ha="center", va="top", fontsize=5,
                color=fs.INK2, transform=blended_transform_factory(ax.transData, ax.transAxes))
        out.append(dict(c=c, n=int(ok.sum()), early=e.median(), late=l_.median(), drop=(e - l_).median(),
                        frac_drop=float((e > l_).mean()), p=p, off_axis=n_clip))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for g in range(len(PAIRED)) for x in (g * gap, g * gap + 1)], ["early", "late"] * len(PAIRED),
                  fontsize=5.5)
    ax.set_xlim(-0.6, (len(PAIRED) - 1) * gap + 1.6)
    ax.set_ylim(lo_clip - 0.1, hi_clip + 0.75)
    ax.set_yticks([-1, 0, 1, 2, 3])
    ax.spines["left"].set_bounds(lo_clip - 0.1, hi_clip + 0.1)
    ax.set_ylabel("linearity index")
    ax.tick_params(axis="x", length=0)
    ax.text(x1 + 0.7, 1, "linear\nsum", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95,
            clip_on=False)
    ax.text(x1 + 0.7, 0, "dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2,
            linespacing=0.95, clip_on=False)
    return pd.DataFrame(out)


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    Pop = pd.read_parquet(out_root() / "binocular" / "v3_C_population.parquet")
    Boot = pd.read_parquet(out_root() / "binocular" / "v3_C_population_boot.parquet")
    Pen = pd.read_parquet(out_root() / "binocular" / "v3_C_penetration.parquet")
    a_med = float(meta.a.median())

    fig = fs.figure("double", 122)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1], hspace=0.04)
    ax_a, ax_b1, ax_b2 = top.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    timecourses(ax_a, pop)
    rel, O = pop["rel"], pop["onset"]
    ymax = 1.08 * O[:, 5, 0, :][:, (rel >= 0) & (rel < 20)].mean(axis=1).mean()
    k_early = equivalent_contrast_sketch(ax_b1, pop, (0, 20), "first 20 ms after onset", a_med, ymax)
    k_late = equivalent_contrast_sketch(ax_b2, pop, (50, 200), "50–200 ms after onset", a_med, ymax)
    ax_b1.set_ylabel("MUA (dominant-eye peak = 1)")
    ax_b2.tick_params(labelleft=False)
    for k, (mk, col, lab) in enumerate((("D", fs.INK, "both eyes, 45% each"), ("o", fs.COND["dominant"], "dominant eye alone"),
                                        ("^", fs.MUTED, "linear sum, (1 + eye weight) × 45%"))):
        y = 0.62 - k * 0.075
        ax_b2.plot([0.06], [y], marker=mk, color=col, ms=2.8, transform=ax_b2.transAxes, clip_on=False)
        ax_b2.text(0.11, y, lab, transform=ax_b2.transAxes, va="center", fontsize=5.5, color=fs.INK)
    left, right = bot.subfigures(1, 2, width_ratios=[1.25, 1.0], wspace=0.03)
    axs_c = left.subplots(1, 2)
    C = iso_contours(axs_c, Pop, Boot, Pen)
    left.supxlabel("other-eye contrast × eye weight (fraction of target)", fontsize=7)
    ax_d = right.subplots()
    S = early_late(ax_d, W)
    for ax, letter in zip((ax_a, ax_b1, ax_d), "abd"):
        fs.panel_label(ax, letter)
    fs.panel_label(axs_c[0], "c", y=33)                  # above the title and the exponent note
    print(fs.save(fig, "fig2_add_then_normalize"))
    pd.set_option("display.width", 220)
    print(f"illustration equivalent contrast / 45%: early {k_early:.2f}, late {k_late:.2f}; median eye weight "
          f"{a_med:.2f}; {len(meta)} penetrations")
    print("(c) iso-response contours (population exponent [95% CI], suppression fraction, bounds; per-penetration "
          "medians, k/n larger late, two-sided Wilcoxon):")
    print(C.round(3).to_string(index=False))
    print("(d) early vs late per penetration:")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
