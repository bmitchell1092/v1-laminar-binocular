"""Figure 2: the two eyes add at response onset, then normalize within tens of milliseconds (docs/09).

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_add_then_normalize.py
Built for readers outside the field (2026-11-02, docs/11 "shared reading"): what each hypothesis predicts, then the
data in the same format, then what the other eye adds.
(a) illustration (hand-set curves): if the eyes' contrasts add, the response to both eyes follows the dominant eye's
    response at the summed contrast (grey dashed); if the eyes normalize each other, it follows the dominant eye
    alone (blue);
(b) population responses at 45% aligned to each channel's onset: both eyes at 45% each, the dominant eye alone at
    45%, and the summation prediction, the dominant eye's response at 45% x (1 + a) (each penetration's median eye
    weight a; v1lb.binocular.response_at, log-contrast interpolation between 45 and 90% at each time point);
    penetration means, 95% bootstrap CI, 5 ms running mean for display;
(c) the response the other eye adds, both eyes minus the dominant eye alone at 45%, against what summation predicts
    it adds (prediction minus the dominant eye alone); same means and CIs. The windows of the dip-and-recovery check
    (docs/09 F; scripts/v4_dip_recovery.py, docs/qc/v4_F.csv) are printed: 0-20, 30-100, 130-200 ms;
(d) per penetration, early (0-20 ms after onset) vs late (50-200 ms) linearity index (mean of the windows' channel
    medians, as in equivalent_contrast.py), including 5.5%; Wilcoxon signed-rank tests; values beyond the axis are
    drawn at its edge as open triangles (counted in the printout, for the legend). At 90% the equivalent contrast
    cannot be measured (the dominant eye has no higher contrast to match).
Moved out on 2026-11-02: the equivalent-contrast readout (former b) to Supplementary Fig S1, the iso-response contour
(former c) to Supplementary Fig S3 (which has all three targets and the per-penetration test).
Inputs: figure_data.py (population.npz, population_meta.parquet), <out_root>/binocular/k_windows.parquet,
docs/qc/v4_F.csv.
Writes <out_root>/figures/fig2_add_then_normalize.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, other_eye_added  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
PAIRED = [0.055, 0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts: 0-20 ms and 50-200 ms after onset
F_WINDOWS = {"0-20 ms": (0, 20), "30-100 ms": (30, 100), "130-200 ms": (130, 200)}
H_MM = 122
POS = {"a1": (14, 95, 36, 14), "a2": (14, 67, 36, 14), "b": (68, 64, 44, 48), "c": (130, 64, 41, 48),   # (x, y, w, h)
       "d": (16, 9, 144, 38)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def predictions(axs):
    """Two hand-set sketches in the format of panel b: what each hypothesis predicts."""
    t = np.arange(-20, 201)
    s = np.where(t < 0, 0, (t / 9) ** 2 * np.exp(-(t / 9 - 1) * 2) * 0.75 + 0.3 * (1 - np.exp(-np.maximum(t, 0) / 12)))
    s = np.where(t < 0, 0, s)
    titles = ("if the eyes' contrasts add:\nblack follows the dashed line",
              "if the eyes normalize each other:\nblack follows blue")
    for ax, title, gain in zip(axs, titles, (1.3, 1.0)):
        ax.plot(t, 1.3 * s, color=fs.COND["linear"], ls=fs.SUM_LS, lw=1.0, zorder=3)
        ax.plot(t, s, color=fs.COND["dominant"], lw=1.0)
        ax.plot(t, gain * s + 0.05, color=fs.COND["both"], lw=1.2)
        ax.set_title(title, fontsize=6, loc="left", pad=2, linespacing=1.05)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-20, 200)
        ax.set_ylim(-0.05, 1.45)
        ax.set_ylabel("response", fontsize=6, labelpad=2)
    axs[1].set_xlabel("time", fontsize=6, labelpad=2)
    axs[0].text(0.0, 1.75, "what to look for (illustration)", transform=axs[0].transAxes, ha="left", va="bottom",
                fontsize=6, color=fs.INK2)


def timecourses(ax, rel, dom, both, pred_resp):
    sel = (rel >= -20) & (rel <= 200)
    items = [("both eyes, 45% each", fs.COND["both"], "-", both, 1.2),
             ("dominant eye alone, 45%", fs.COND["dominant"], "-", dom, 1.0),
             (f"{fs.SUM}\n(45% + weight × 45%)", fs.COND["linear"], fs.SUM_LS, pred_resp, 1.0)]
    for lab, col, ls, Y, lw in items[::-1]:
        fs.ci_line(ax, rel[sel], fs.smooth(Y)[:, sel], col, ls=ls, lw=lw)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_ylim(-0.05, 1.0)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    fs.line_key(ax, [(lab, col, ls) for lab, col, ls, _, _ in items], x=0.38, y=0.95, dy=0.09)
    ax.text(10, 0.99, "first 20 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def added(ax, rel, A, P):
    sel = (rel >= -20) & (rel <= 200)
    fs.added_panel(ax, rel[sel], fs.smooth(A)[:, sel], fs.smooth(P)[:, sel], label_added="measured",
                   key_xy=(0.50, 0.95))
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_ylim(-0.1, 0.34)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("added by the other eye (MUA)\n= both eyes − dominant eye alone")
    out = {}
    for name, (w0, w1) in F_WINDOWS.items():
        s = (rel >= w0) & (rel < w1)
        out[name] = dict(added=float(A[:, s].mean()), predicted=float(P[:, s].mean()))
    return pd.DataFrame(out).T


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
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=fs.SUM_LS, zorder=0)
    fs.zero_line(ax)
    ax.set_xticks([x for g in range(len(PAIRED)) for x in (g * gap, g * gap + 1)], ["0–20", "50–200"] * len(PAIRED),
                  fontsize=5.5)
    ax.set_xlabel("window after onset (ms)")
    ax.set_xlim(-0.6, (len(PAIRED) - 1) * gap + 1.6)
    ax.set_ylim(lo_clip - 0.1, hi_clip + 1.55)
    ax.set_yticks([-1, 0, 1, 2, 3])
    ax.spines["left"].set_bounds(lo_clip - 0.1, hi_clip + 0.1)
    ax.set_ylabel(fs.LINEARITY)
    ax.tick_params(axis="x", length=0)
    ax.text(x1 + 0.7, 1, "1: contrasts\nadd", ha="left", va="center", fontsize=5.5, color=fs.INK2, linespacing=0.95,
            clip_on=False)
    ax.text(x1 + 0.7, 0, "0: dominant\neye alone", ha="left", va="center", fontsize=5.5, color=fs.INK2,
            linespacing=0.95, clip_on=False)
    return pd.DataFrame(out)


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    F = pd.read_csv(REPO / "docs" / "qc" / "v4_F.csv")
    O, rel = pop["onset"], pop["rel"]
    a = meta.set_index("penetration").loc[[str(p) for p in pop["pens"]], "a"].to_numpy()
    i45 = list(CONTRASTS).index(0.45)
    A, P = other_eye_added(O, a, 0.45)
    dom = O[:, i45, 0]

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")                     # axes placed in mm: in-panel keys stay inside their panels
    ax_a = [axes_mm(fig, *POS[n]) for n in ("a1", "a2")]
    ax_b, ax_c, ax_d = (axes_mm(fig, *POS[n]) for n in "bcd")
    predictions(ax_a)
    timecourses(ax_b, rel, dom, O[:, i45, i45], dom + P)
    Fw = added(ax_c, rel, A, P)
    S = early_late(ax_d, W)
    for ax, letter, y in ((ax_a[0], "a", 23), (ax_b, "b", 1.0), (ax_c, "c", 1.0), (ax_d, "d", 1.0)):
        fs.panel_label(ax, letter, y=y)
    print(fs.save(fig, "fig2_add_then_normalize"))
    pd.set_option("display.width", 220)
    print(f"{len(meta)} penetrations; median eye weight {np.median(a):.2f} (summation predicts 45% x "
          f"{1 + np.median(a):.2f} = {45 * (1 + np.median(a)):.0f}% for the median penetration)")
    print("(c) population means of what the other eye adds and of the summation prediction (dominant-eye peak = 1):")
    print(Fw.round(4).to_string())
    print("    per-penetration check (docs/09 F, 45%):")
    print(F[np.isclose(F.c, 0.45)][["test", "n", "median", "n_in_direction", "p", "p_odd", "p_even", "holds"]]
          .round(4).to_string(index=False))
    print("(d) early vs late per penetration:")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
