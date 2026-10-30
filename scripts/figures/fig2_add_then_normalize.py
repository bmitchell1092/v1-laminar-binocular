"""Figure 2: the two eyes add at response onset, then normalize within tens of milliseconds (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_add_then_normalize.py
(a) population responses at 45% aligned to each channel's onset: each eye alone at 45%, both eyes at 45%, and the
    dominant eye at 90% (about what linear summation of contrast predicts); penetration means, 95% bootstrap CI;
(b) how the equivalent contrast is read (population means, illustration), on one shared response axis: the
    binocular response at 45% is placed on the dominant eye's contrast response measured in the same window and
    read off as the dominant-eye contrast that gives the same response;
(c) linearity index over time (1 = linear summation, 0 = binocular equals the dominant eye alone) of the
    population-average response, read as in (b): sliding 15 ms windows every 5 ms from each channel's onset, one eye
    weight from the whole response; open triangles = binocular above the dominant eye at 90% (a lower bound); 95% CI
    from 500 bootstrap resamples of penetrations. Averaging responses before the readout avoids the heavy tails of
    per-channel ratios; the per-penetration statistics are in (d);
(d) per penetration, early (0-20 ms after onset) vs late (50-200 ms) linearity index (mean of the windows' channel
    medians, as in equivalent_contrast.py), including 5.5%; Wilcoxon signed-rank tests; values beyond the axis are
    drawn at its edge. At 90% the equivalent contrast cannot be measured (the dominant eye has no higher contrast
    to match); Fig 3a shows the raw binocular gain at 90%.
Inputs: figure_data.py (population.npz), <out_root>/binocular/k_windows.parquet.
Writes <out_root>/figures/fig2_add_then_normalize.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight, invert_crf  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
DIOPTIC = [0.11, 0.225, 0.45]
PAIRED = [0.055, 0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts: 0-20 ms and 50-200 ms after onset
SLIDE = 15                                            # sliding-window width (ms), as in normalization_timing.py


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
    ax.text(10, ax.get_ylim()[1] * 0.98, "first\n20 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2,
            linespacing=0.95)


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


def population_L(O, rel, c, starts):
    """Linearity index of the population-average response (as in panel b) in sliding windows; O [pen, 6, 6, rel].
    One eye weight from the whole response (0-200 ms after onset), as for the per-penetration statistic in (d): a
    weight refitted in the first windows, before the other eye has fully responded, would inflate the index.
    Returns the index and whether k was top-censored (binocular above the dominant eye at 90%: a lower bound)."""
    i = list(CONTRASTS).index(c)
    M = np.nanmean(O, axis=0)
    whole = (rel >= 0) & (rel < 200)
    a = eye_weight(M[:, 0, whole].mean(axis=1), M[0, :, whole].mean(axis=1))
    out, censored = np.full(len(starts), np.nan), np.zeros(len(starts), bool)
    for k, s0 in enumerate(starts):
        sel = (rel >= s0) & (rel < s0 + SLIDE)
        d, b = M[:, 0, sel].mean(axis=1), M[i, i, sel].mean()
        kk, cens = equivalent_contrast(d[:, None], np.array([b]), c)
        censored[k] = cens[0] == 1
        kk = CONTRASTS[-1] / c if censored[k] else float(kk[0])
        out[k] = (kk - 1) / a
    return out, censored


def linearity_time(axs, pop, n_boot=500):
    """One small panel per contrast (shared axes), 0-120 ms after onset."""
    rel, O = pop["rel"], pop["onset"]
    starts = np.arange(0, 106, 5)
    x = starts + SLIDE / 2
    rng = np.random.default_rng(0)
    out = {}
    for ax, c in zip(axs, PAIRED):
        m, cens = population_L(O, rel, c, starts)
        B = np.array([population_L(O[rng.integers(0, len(O), len(O))], rel, c, starts)[0] for _ in range(n_boot)])
        lo, hi = np.nanpercentile(B, [2.5, 97.5], axis=0)
        ax.fill_between(x, lo, hi, color=fs.CONTRAST[c], alpha=0.2, lw=0)
        ax.plot(x, m, color=fs.CONTRAST[c], lw=1.4)
        ax.plot(x[cens], m[cens], ls="none", marker="^", ms=3, mfc="white", mec=fs.CONTRAST[c], mew=0.7)
        ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
        fs.zero_line(ax)
        ax.set_xlim(0, 120)
        ax.set_xticks([0, 50, 100])
        ax.set_ylim(-1.2, 2.0)
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
        out[c] = (x, m)
    for ax in axs[1:]:
        ax.tick_params(labelleft=False)
    axs[0].set_ylabel("linearity index")
    axs[-1].text(123, 1, "linear\nsum", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    axs[-1].text(123, 0, "dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2,
                 linespacing=0.95)
    return out


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
        ax.plot([x0, x1], [e.median(), l_.median()], color=fs.CONTRAST[c], lw=1.6, marker="o", ms=3.2, zorder=3)
        ax.text((x0 + x1) / 2, hi_clip + 0.15, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        ax.text((x0 + x1) / 2, lo_clip - 0.55, fs.contrast_label(c), ha="center", va="top", fontsize=6.5,
                color=fs.INK)
        out.append(dict(c=c, n=int(ok.sum()), early=e.median(), late=l_.median(), drop=(e - l_).median(),
                        frac_drop=float((e > l_).mean()), p=p))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for g in range(len(PAIRED)) for x in (g * gap, g * gap + 1)], ["early", "late"] * len(PAIRED),
                  fontsize=5.5)
    ax.set_xlim(-0.6, (len(PAIRED) - 1) * gap + 1.6)
    ax.set_ylim(lo_clip - 0.1, hi_clip + 0.1)
    ax.set_ylabel("linearity index")
    ax.tick_params(axis="x", length=0)
    return pd.DataFrame(out)


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    a_med = float(meta.a.median())

    fig = fs.figure("double", 112)
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
                                        ("^", fs.MUTED, "linear sum, (1 + a)·45%"))):
        y = 0.62 - k * 0.075
        ax_b2.plot([0.06], [y], marker=mk, color=col, ms=2.8, transform=ax_b2.transAxes, clip_on=False)
        ax_b2.text(0.11, y, lab, transform=ax_b2.transAxes, va="center", fontsize=5.5, color=fs.INK)
    left, right = bot.subfigures(1, 2, width_ratios=[1.45, 1.0], wspace=0.02)
    axs_c = left.subplots(1, len(PAIRED), sharey=True)
    curves = linearity_time(axs_c, pop)
    left.supxlabel("time from response onset (ms, window centre)", fontsize=7)
    ax_c = axs_c[0]
    ax_d = right.subplots()
    S = early_late(ax_d, W)
    for ax, letter in zip((ax_a, ax_b1, ax_c, ax_d), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig2_add_then_normalize"))
    print(f"illustration k: early {k_early:.2f}, late {k_late:.2f}; median eye weight a {a_med:.2f}; "
          f"{len(meta)} penetrations")
    for c, (x, m) in curves.items():
        print(f"sliding L at {c}: first window {m[0]:.2f}; min {np.nanmin(m):.2f} at {x[np.nanargmin(m)]:.1f} ms")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
