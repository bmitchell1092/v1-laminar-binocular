"""Figure 2: the two eyes add at response onset, then normalize within tens of milliseconds (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_add_then_normalize.py
(a) population responses at 45% aligned to each channel's onset: each eye alone at 45%, both eyes at 45%, and the
    dominant eye at 90% (about what linear summation of contrast predicts); penetration means, 95% bootstrap CI;
(b) how the equivalent contrast is read (population means, illustration), each window on its own response axis: the
    binocular response at 45% is placed on the dominant eye's contrast response measured in the same window and
    read off as the dominant-eye contrast that gives the same response;
(c) iso-response contour across the dichoptic conditions (v3_iso_contours.py) for a target equal to the dominant
    eye's response at 45%, early (0-20 ms after each penetration's onset) and late (50-200 ms): for each other-eye
    contrast, the dominant-eye contrast at which the binocular response reaches the target, from the
    population-average window responses, with 95% penetration-bootstrap CIs, in percent contrast. The other eye's
    contrast is multiplied by the eye weight of the same window, so linear summation is the diagonal and "the
    dominant eye alone" the horizontal line at the target; open circle = the other eye alone already reaches the
    target; open triangle = CI beyond the axis. The note gives the per-penetration test of the fitted Minkowski
    exponent, early vs late (two-sided Wilcoxon); the fitted contours and the 11 and 22.5% targets are in
    figS3_iso_contours.py (statistics for 22.5% are still printed here);
(d) per penetration, early (0-20 ms after onset) vs late (50-200 ms) linearity index (mean of the windows' channel
    medians, as in equivalent_contrast.py), including 5.5%; Wilcoxon signed-rank tests; values beyond the axis are
    drawn at its edge as open triangles (counted in the printout, for the legend). At 90% the equivalent contrast cannot be measured (the
    dominant eye has no higher contrast to match).
The linearity index over time (formerly panel c) is Fig 4c, with the model predictions.
Inputs: figure_data.py (population.npz), <out_root>/binocular/k_windows.parquet, v3_C_population.parquet,
v3_C_population_boot.parquet, v3_C_penetration.parquet.
Writes <out_root>/figures/fig2_add_then_normalize.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, invert_crf  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
PAIRED = [0.055, 0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts: 0-20 ms and 50-200 ms after onset
CONTOUR_TARGET = 0.45                                 # shown; 22.5% is printed and in Fig S3
EPOCH = {"early": fs.MUTED, "late": fs.INK}           # docs/11: epoch grey ramp


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


def equivalent_contrast_sketch(ax, pop, w, title, a_med):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= w[0]) & (rel < w[1])
    crf = O[:, :, 0, :][:, :, sel].mean(axis=2).mean(axis=0)              # dominant eye, population mean
    i = list(CONTRASTS).index(0.45)
    rb = O[:, i, i, :][:, sel].mean(axis=1).mean()
    kc, cens = invert_crf(crf, np.array([rb]))
    kc = float(kc[0]) if cens[0] == 0 else CONTRASTS[-1]
    c = CONTRASTS[1:]
    ymax = 1.12 * max(crf.max(), rb)
    ax.plot(c, crf[1:], color=fs.COND["dominant"], marker="o", ms=2.6, lw=1.0)
    ax.plot([0.45], [rb], color=fs.INK, marker="D", ms=3.4, zorder=5)
    ax.plot([0.45, kc], [rb, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([kc, kc], [0, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([(1 + a_med) * 0.45], [0], marker="^", ms=3.8, color=fs.MUTED, clip_on=False, zorder=6)
    fs.log_contrast_axis(ax, (0.055, 0.11, 0.225, 0.45, 0.9))
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"], fontsize=5.5)
    ax.set_xlim(0.045, 1.05)
    ax.set_ylim(0, ymax)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xlabel("dominant-eye contrast (%)")
    ax.set_title(title, fontsize=6.5, color=fs.INK)
    ax.text(0.05, 0.97, f"both eyes at 45%\n= dominant eye at {100 * kc:.0f}%", transform=ax.transAxes,
            ha="left", va="top", fontsize=6, color=fs.INK, linespacing=1.0)
    return kc / 0.45


def contour_stats(Pop, Pen, cs) -> list[dict]:
    """Population exponent per window and the per-penetration early vs late test of the exponent."""
    out = []
    for w in ("early", "late"):
        r = Pop[(Pop.window == w) & (Pop.cstar == cs)].iloc[0]
        out.append(dict(target=cs, window=w, eye_weight=r.a, m=r.m, m_lo=r.m_lo, m_hi=r.m_hi,
                        suppression=r.suppression))
    g = Pen[Pen.cstar == cs].pivot(index="penetration", columns="window", values="m")[["early", "late"]].dropna()
    out[-1].update(pen_n=len(g), pen_early=g.early.median(), pen_late=g.late.median(),
                   pen_k=int((g.late > g.early).sum()), pen_p=stats.wilcoxon(g.late, g.early).pvalue)
    return out


def iso_contour(ax, Pop, Boot, Pen):
    """Population iso-response contour for the 45% target, early vs late, in percent contrast: the dominant-eye
    contrast needed to reach the target at each weighted other-eye contrast, with 95% penetration-bootstrap CIs."""
    T = 100 * CONTOUR_TARGET
    xmax, ymax = 1.55 * T, 1.33 * T
    ax.plot([0, T], [T, 0], color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=1)
    ax.plot([0, xmax], [T, T], color=fs.RULE, lw=0.9, zorder=0)
    for k, w in enumerate(("early", "late")):                 # late last: its final point carries the label
        col = EPOCH[w]
        p = Pop[(Pop.window == w) & (Pop.cstar == CONTOUR_TARGET)].sort_values("c_N")
        b = Boot[(Boot.window == w) & (Boot.cstar == CONTOUR_TARGET)].groupby("c_N").x
        lo = T * b.quantile(0.025).reindex(p.c_N).to_numpy()
        hi = T * b.quantile(0.975).reindex(p.c_N).to_numpy()
        xv, yv = T * p.y.to_numpy() + (k - 0.5) * 0.8, T * p.x.to_numpy()
        bound = p.status.to_numpy() != "ok"
        ax.vlines(xv, lo, np.minimum(hi, ymax), color=col, lw=0.6, zorder=3)
        ax.plot(xv, yv, color=col, lw=1.1 if w == "late" else 1.0, zorder=3)
        ax.scatter(xv[~bound], yv[~bound], s=11, color=col, zorder=4, lw=0)
        ax.scatter(xv[bound], yv[bound], s=11, facecolor="white", edgecolor=col, lw=0.7, zorder=4, clip_on=False)
        top = hi > ymax
        ax.scatter(xv[top], np.full(top.sum(), ymax), marker="^", s=9, facecolor="white", edgecolor=col, lw=0.6,
                   zorder=4, clip_on=False)
        last = (xv[~bound][-1], yv[~bound][-1])
    ax.text(1.11 * T, 0.18 * T, "first 20 ms", fontsize=6.5, color=fs.INK, ha="center", va="bottom")
    ax.annotate("50–200 ms", last, xytext=(0, 7), textcoords="offset points", fontsize=6.5, color=fs.INK,
                ha="center", va="bottom")
    ax.text(T * 0.60, T * 0.42, "linear sum", rotation=-45, rotation_mode="anchor", fontsize=5.5, color=fs.INK2,
            ha="center", va="bottom")
    ax.text(T * 0.85, T + 1, "dominant eye alone", fontsize=5.5, color=fs.INK2, ha="center", va="bottom")
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_yticks([0, 20, 40, 60])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("other-eye contrast × eye weight (%)")
    ax.set_ylabel("dominant-eye contrast needed (%)")
    ax.set_title(f"same response as the dominant eye alone at {fs.contrast_label(CONTOUR_TARGET)}", fontsize=6.5)
    st = contour_stats(Pop, Pen, CONTOUR_TARGET)[-1]
    ax.text(0.03, 0.03, f"late contour further\nfrom the linear sum:\n{st['pen_k']}/{st['pen_n']} penetrations, "
            f"{fs.fmt_p(st['pen_p'])}", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5, color=fs.INK,
            linespacing=1.05)


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
        ax.text((x0 + x1) / 2, hi_clip + 0.95, f"both eyes {fs.contrast_label(c)}", ha="center", va="bottom",
                fontsize=6.5, color=fs.INK)
        ax.text((x0 + x1) / 2, hi_clip + 0.3, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        out.append(dict(c=c, n=int(ok.sum()), early=e.median(), late=l_.median(), drop=(e - l_).median(),
                        frac_drop=float((e > l_).mean()), p=p, off_axis=n_clip))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for g in range(len(PAIRED)) for x in (g * gap, g * gap + 1)], ["0–20", "50–200"] * len(PAIRED),
                  fontsize=5.5)
    ax.set_xlabel("window after onset (ms)")
    ax.set_xlim(-0.6, (len(PAIRED) - 1) * gap + 1.6)
    ax.set_ylim(lo_clip - 0.1, hi_clip + 1.55)
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
    k_early = equivalent_contrast_sketch(ax_b1, pop, (0, 20), "first 20 ms after onset", a_med)
    k_late = equivalent_contrast_sketch(ax_b2, pop, (50, 200), "50–200 ms after onset", a_med)
    ax_b1.set_ylabel("MUA (dominant-eye peak = 1)")
    for k, (mk, col, lab) in enumerate((("D", fs.INK, "both eyes, 45% each"), ("o", fs.COND["dominant"], "dominant eye alone"),
                                        ("^", fs.MUTED, "linear-sum prediction"))):
        y = 0.75 - k * 0.075
        ax_b1.plot([0.07], [y], marker=mk, color=col, ms=2.8, transform=ax_b1.transAxes, clip_on=False)
        ax_b1.text(0.12, y, lab, transform=ax_b1.transAxes, va="center", fontsize=5.5, color=fs.INK)
    ax_c, ax_d = bot.subplots(1, 2, gridspec_kw=dict(width_ratios=[0.85, 1.3]))
    iso_contour(ax_c, Pop, Boot, Pen)
    C = pd.DataFrame([r for cs in (0.225, CONTOUR_TARGET) for r in contour_stats(Pop, Pen, cs)])
    S = early_late(ax_d, W)
    for ax, letter in zip((ax_a, ax_b1, ax_c, ax_d), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig2_add_then_normalize"))
    pd.set_option("display.width", 220)
    print(f"illustration equivalent contrast / 45%: early {k_early:.2f}, late {k_late:.2f}; median eye weight "
          f"{a_med:.2f}; {len(meta)} penetrations")
    print("(c) iso-response contours (population exponent [95% CI], suppression fraction; per-penetration "
          "medians, k/n larger late, two-sided Wilcoxon; 45% shown, 22.5% in Fig S3):")
    print(C.round(3).to_string(index=False))
    print("(d) early vs late per penetration:")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
