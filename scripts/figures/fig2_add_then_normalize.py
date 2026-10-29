"""Figure 2: the two eyes add at response onset, then normalize within tens of milliseconds (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_add_then_normalize.py
(a) population responses at 45% aligned to each channel's onset: each eye alone, both eyes, and the dominant eye
    at 90% (about what linear summation of contrast predicts); penetration means, 95% bootstrap CI;
(b) how the equivalent contrast k is read (population means, illustration): the binocular response at 45% is
    placed on the dominant eye's contrast response in the same window; early it sits near the linear-summation
    contrast (1 + a) c, late near c;
(c) linearity index L = (k - 1) / a per onset-aligned window (1 = linear summation, 0 = binocular equals the
    dominant eye alone), median of penetration medians with 95% bootstrap CI (equivalent_contrast.py);
(d) per penetration, early (0-20 ms after onset) vs late (50-200 ms) L (mean of the windows' channel medians, as
    in equivalent_contrast.py), Wilcoxon signed-rank test; lines beyond the axis limits are clipped.
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
from v1lb.binocular import CONTRASTS, invert_crf  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
DIOPTIC = [0.11, 0.225, 0.45]
WINDOWS = [(0, 10), (10, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150), (150, 200)]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts: 0-20 ms and 50-200 ms after onset


def timecourses(ax, pop):
    rel, O = pop["rel"], pop["onset"]
    i = list(CONTRASTS).index(0.45)
    sel = (rel >= -20) & (rel <= 200)
    items = [("dominant eye, 90%", fs.COND["linear"], (0, (3, 1.5)), O[:, 5, 0]),
             ("dominant eye", fs.COND["dominant"], "-", O[:, i, 0]),
             ("other eye", fs.COND["other"], "-", O[:, 0, i]),
             ("both eyes", fs.COND["both"], "-", O[:, i, i])]
    for lab, col, ls, Y in items:
        fs.ci_line(ax, rel[sel], Y[:, sel], col, ls=ls, lw=1.2 if col == fs.INK else 1.0)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    fs.line_key(ax, [(lab, col, ls) for lab, col, ls, _ in items], x=0.5, y=0.93)
    ax.text(0.5, 0.56, "45% contrast; dashed ≈ linear sum", transform=ax.transAxes, fontsize=5.5, color=fs.INK2)
    ax.text(10, ax.get_ylim()[1] * 0.98, "early", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def equivalent_contrast_sketch(ax, pop, w, title, a_med):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= w[0]) & (rel < w[1])
    crf = O[:, :, 0, :][:, :, sel].mean(axis=2).mean(axis=0)              # dominant eye, population mean
    i = list(CONTRASTS).index(0.45)
    rb = O[:, i, i, :][:, sel].mean(axis=1).mean()
    kc, cens = invert_crf(crf, np.array([rb]))
    kc = float(kc[0]) if cens[0] == 0 else CONTRASTS[-1]
    c = CONTRASTS[1:]
    ax.plot(c, crf[1:], color=fs.COND["dominant"], marker="o", ms=2.6, lw=1.0)
    ax.plot([0.45], [crf[i]], color=fs.COND["dominant"], marker="o", ms=2.6)
    ax.plot([0.45], [rb], color=fs.INK, marker="D", ms=3.2, zorder=5)
    ax.plot([0.45, kc], [rb, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    ax.plot([kc, kc], [0, rb], color=fs.INK, lw=0.6, ls=(0, (2, 1.5)))
    lin = (1 + a_med) * 0.45
    ax.plot([lin], [0], marker="^", ms=3.5, color=fs.MUTED, clip_on=False, zorder=6)
    fs.log_contrast_axis(ax, (0.055, 0.11, 0.225, 0.45, 0.9))
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"], fontsize=5.5)
    ax.set_xlim(0.045, 1.05)
    ax.set_ylim(0, None)
    ax.set_xlabel("dominant-eye contrast (%)")
    ax.set_title(title, fontsize=6.5, color=fs.INK)
    ax.text(kc, rb * 0.5, f" k = {kc / 0.45:.1f}", fontsize=5.5, va="center", ha="left" if kc < 0.6 else "right",
            color=fs.INK)
    return kc / 0.45


def linearity_time(ax, W):
    centres = np.array([(a + b) / 2 for a, b in WINDOWS])
    for k, c in enumerate(DIOPTIC):
        per = W[W.c == c].groupby(["penetration", "w0"]).L.median().unstack().reindex(columns=[w for w, _ in WINDOWS])
        m, lo, hi = fs.boot_ci(per.to_numpy(), stat=np.nanmedian)
        x = centres + (k - 1) * 2.2
        ax.errorbar(x, m, [m - lo, hi - m], color=fs.CONTRAST[c], lw=1.1, elinewidth=0.6, capsize=0, marker="o",
                    ms=2.8)
    fs.line_key(ax, [(fs.contrast_label(c), fs.CONTRAST[c], "-", "o") for c in DIOPTIC], x=0.8, y=0.97, dy=0.075,
                length=0.06)
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.text(203, 1, "linear\nsummation", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    ax.text(203, 0, "dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    ax.set_xlim(0, 200)
    ax.set_ylim(-1.0, 1.75)
    ax.set_xlabel("time from response onset (ms, window centre)")
    ax.set_ylabel("linearity index L")


def early_late(ax, W):
    out = []
    lo_clip, hi_clip = -1.5, 3.0
    for g, c in enumerate(DIOPTIC):
        per = W[W.c == c].groupby(["penetration", "w0"]).L.median().unstack()
        e, l_ = per[EARLY].mean(axis=1), per[LATE].mean(axis=1)
        ok = e.notna() & l_.notna()
        e, l_ = e[ok], l_[ok]
        p = stats.wilcoxon(e, l_).pvalue
        x0, x1 = g * 2.6, g * 2.6 + 1
        for a, b in zip(e.clip(lo_clip, hi_clip), l_.clip(lo_clip, hi_clip)):
            ax.plot([x0, x1], [a, b], color=fs.RULE, lw=0.5, zorder=1)
        ax.plot([x0, x1], [e.median(), l_.median()], color=fs.CONTRAST[c], lw=1.6, marker="o", ms=3.2, zorder=3)
        ax.text((x0 + x1) / 2, hi_clip + 0.15, fs.fmt_p(p), ha="center", va="bottom", fontsize=5.5, color=fs.INK)
        ax.text((x0 + x1) / 2, lo_clip - 0.55, fs.contrast_label(c), ha="center", va="top", fontsize=6.5,
                color=fs.INK)
        n_clip = int(((e < lo_clip) | (e > hi_clip) | (l_ < lo_clip) | (l_ > hi_clip)).sum())
        out.append(dict(c=c, n=int(ok.sum()), early=e.median(), late=l_.median(), drop=(e - l_).median(),
                        frac_drop=float((e > l_).mean()), p=p, clipped=n_clip))
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for g in range(3) for x in (g * 2.6, g * 2.6 + 1)], ["early", "late"] * 3, fontsize=5.5)
    ax.set_xlim(-0.6, 2 * 2.6 + 1.6)
    ax.set_ylim(lo_clip - 0.1, hi_clip + 0.1)
    ax.set_ylabel("linearity index L")
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
    k_early = equivalent_contrast_sketch(ax_b1, pop, (0, 20), "0–20 ms after onset", a_med)
    k_late = equivalent_contrast_sketch(ax_b2, pop, (50, 200), "50–200 ms after onset", a_med)
    ax_b1.set_ylabel("MUA (dominant-eye peak = 1)")
    for k, (mk, col, lab) in enumerate((("D", fs.INK, "both eyes at 45%"), ("o", fs.COND["dominant"], "dominant eye"),
                                        ("^", fs.MUTED, "linear sum, (1 + a)·45%"))):
        y = 0.95 - k * 0.075
        ax_b1.plot([0.06], [y], marker=mk, color=col, ms=2.8, transform=ax_b1.transAxes, clip_on=False)
        ax_b1.text(0.11, y, lab, transform=ax_b1.transAxes, va="center", fontsize=5.5, color=fs.INK)
    ax_c, ax_d = bot.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.45, 1]))
    linearity_time(ax_c, W)
    S = early_late(ax_d, W)
    for ax, letter in zip((ax_a, ax_b1, ax_c, ax_d), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig2_add_then_normalize"))
    print(f"illustration k: early {k_early:.2f}, late {k_late:.2f}; median eye weight a {a_med:.2f}; "
          f"{len(meta)} penetrations")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
