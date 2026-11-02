"""Figure 1: design and an example (docs/09, figure set for V3).

Usage: .venv\\Scripts\\python scripts\\figures\\fig1_setup.py [example penetration, default 221212_J_eA]
(a) the 6 x 6 dichoptic contrast design (each eye 0, 5.5, 11, 22.5, 45 or 90%) and the trial sequence (three
    conditions per trial, ~317 ms on and ~317 ms off); each channel's axes are later reoriented to its dominant
    and other eye;
(b) example penetration in the Fig 2a format, channel mean with each channel aligned to its own response onset
    (half-max, dominant eye at 90%), oriented to its dominant eye and divided by its dominant-eye peak at 90%: each
    eye alone at 45%, both eyes at 45%, and the dominant eye at 90%, about what linear summation of contrast
    predicts ((1 + a) c with eye weight a ~ 0.6-0.9); the figure names the monkey, not the session;
(c) monocular contrast response per eye, mean MUA 40-250 ms after stimulus onset in units of each channel's
    dominant-eye peak, penetration means (95% bootstrap CI over penetrations);
(d) ocular dominance index at 90% for driven MUA channels and isolated single units.
The stereoscope schematic and the example flash CSD moved to Supplementary Fig S1 (figS1_methods.py).
Inputs: figure_data.py (population.npz, population_meta.parquet), run_crf.py, run_dichoptic.py.
Writes <out_root>/figures/fig1_setup.pdf and .png; prints the statistics shown.
"""

import glob
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
H_MM = 76
POS = {"matrix": (11, 37, 31, 31), "key": (44, 37, 21, 31), "trial": (2, 12, 62, 14),    # (x, y, w, h) in mm
       "example": (78, 12, 44, 58), "crf": (137, 49, 25, 22), "odi": (137, 11, 33, 22)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def design(ax, kax, tax):
    cls = np.zeros((6, 6), int)
    for i in range(6):
        for j in range(6):
            cls[i, j] = 0 if i == j == 0 else 1 if (i == 0 or j == 0) else 2 if i == j else 3
    cols = {0: "white", 1: fs.INK2, 2: fs.INK, 3: fs.RULE}
    for i in range(6):
        for j in range(6):
            ax.add_patch(mpatches.Rectangle((j + 0.07, i + 0.07), 0.86, 0.86, fc=cols[cls[i, j]],
                                            ec=fs.MUTED if cls[i, j] == 0 else "none", lw=0.5))
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ticks = np.arange(6) + 0.5
    labs = ["0"] + [f"{100 * c:g}" for c in CONTRASTS[1:]]
    ax.set_xticks(ticks, labs)
    ax.set_yticks(ticks, labs)
    ax.tick_params(length=0, pad=1.5, labelsize=5.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("contralateral eye contrast (%)", labelpad=1.5)
    ax.set_ylabel("ipsilateral eye contrast (%)", labelpad=1.5)

    # key and notes beside the matrix
    kax.set_xlim(0, 1)
    kax.set_ylim(0, 1)
    kax.axis("off")
    for k, (lab, y) in enumerate((("one eye", 0.95), ("both eyes,\nequal", 0.82), ("both eyes,\nunequal", 0.66))):
        kax.add_patch(mpatches.Rectangle((0.0, y - 0.035), 0.14, 0.07, fc=cols[k + 1], ec="none"))
        kax.text(0.2, y, lab, va="center", fontsize=5.5, color=fs.INK, linespacing=0.95)
    kax.text(0.0, 0.53, "~30 presentations\nper condition,\nrandom order", va="top", fontsize=5.5,
             color=fs.INK2, linespacing=1.0)
    kax.text(0.0, 0.25, "analysed per\nchannel as\ndominant eye ×\nother eye", va="top", fontsize=5.5,
             color=fs.INK2, linespacing=1.0)

    # trial: fixation, then three conditions in random order (~317 ms on, ~317 ms off)
    tax.set_xlim(-450, 1950)
    tax.set_ylim(-0.6, 1.25)
    tax.axis("off")
    tax.plot([-400, 1900], [0.3, 0.3], color=fs.MUTED, lw=0.6)
    tax.text(-400, 0.45, "fixate", fontsize=5.5, color=fs.INK2, va="bottom")
    for k in range(3):
        x0 = k * 634
        tax.add_patch(mpatches.Rectangle((x0, 0.3), 317, 0.42, fc=fs.INK, ec="none"))
        tax.text(x0 + 158, 0.8, f"condition {k + 1}", ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    tax.plot([1583, 1900], [0.05, 0.05], color=fs.INK, lw=0.8)
    tax.text(1741, -0.08, "317 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2)


def example(ax, pop, k, monkey):
    rel = pop["rel"]
    Y = pop["onset"][k]
    i45 = list(CONTRASTS).index(0.45)
    sel = (rel >= -20) & (rel <= 200)
    items = [("dominant eye, 90%\n(≈ linear sum)", fs.COND["linear"], (0, (3, 1.5)), Y[5, 0]),
             ("dominant eye, 45%", fs.COND["dominant"], "-", Y[i45, 0]),
             ("other eye, 45%", fs.COND["other"], "-", Y[0, i45]),
             ("both eyes, 45% each", fs.COND["both"], "-", Y[i45, i45])]
    for lab, col, ls, y in items:
        ax.plot(rel[sel], y[sel], color=col, ls=ls, lw=1.2 if col == fs.INK else 1.0)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-20, 200)
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    top = 1.08 * np.nanmax(Y[5, 0][sel])
    ax.set_ylim(-0.1, top)
    for (lab, col, ls, _), yy in zip(items, (0.92, 0.83, 0.765, 0.70)):
        ax.plot([0.40, 0.48], [yy, yy], color=col, ls=ls, lw=1.0, transform=ax.transAxes)
        ax.text(0.51, yy, lab, transform=ax.transAxes, va="center", fontsize=5.5, color=fs.INK, linespacing=0.95)
    ax.text(0.40, 0.62, f"example penetration,\nmonkey {monkey}", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.5, color=fs.INK2, linespacing=1.0)
    ax.text(23, top, "first 20 ms", ha="left", va="top", fontsize=5.5, color=fs.INK2)
    e = (rel >= 0) & (rel < 20)
    late = (rel >= 50) & (rel < 200)
    return {w: {lab.split("\n")[0]: float(y[m].mean()) for lab, _, _, y in items}
            for w, m in (("0-20 ms", e), ("50-200 ms", late))}


def crfs(ax, pop):
    t = pop["t"]
    w = (t >= 40) & (t < 250)
    F = pop["flash"][:, :, :, w].mean(axis=3)
    c = CONTRASTS[1:]
    out = {}
    for y, col, lab in ((F[:, 1:, 0], fs.COND["dominant"], "dominant\neye"), (F[:, 0, 1:], fs.COND["other"], "other\neye")):
        m, lo, hi = fs.boot_ci(y)
        ax.fill_between(c, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(c, m, color=col, marker="o", ms=2.8, lw=1.1)
        fs.direct_label(ax, c[-1], m[-1], lab, col, dx=4, linespacing=0.9)
        out[lab.replace("\n", " ")] = (m, lo, hi)
    fs.log_contrast_axis(ax, c)
    ax.set_xticklabels(["5.5", "11", "22.5", "45", "90"])
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    ax.set_ylim(0, None)
    ax.set_xlim(0.045, 1.0)
    ax.set_title("40–250 ms after stimulus onset", fontsize=6, color=fs.INK)
    ax.text(0.05, 0.98 * ax.get_ylim()[1], f"mean, 95% CI,\n{F.shape[0]} penetrations", fontsize=5.5,
            color=fs.INK2, va="top", linespacing=1.0)
    ratio = F[:, 0, 5] / F[:, 5, 0]
    return out, ratio


def ocular_dominance(ax):
    D = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "dichoptic" / "*.parquet"))], ignore_index=True)
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[R.condition.str.startswith("mono") & R.driven].groupby(["penetration", "site", "site_id"]).size()
    D = D.join(drv.rename("nd"), on=["penetration", "site", "site_id"])
    D = D[D.nd.notna() & (D.window == "full")]
    isolated = (D.isi_contamination < 0.5) & (D.snr >= 3)                     # lenient isolation (05 C14)
    D = D[(D.site == "mua") | ((D.site == "unit") & isolated)]
    odi = {"mua": [], "unit": []}
    for _, g in D.groupby(["penetration", "site", "site_id"]):
        mm = g.set_index(["contrast_contra", "contrast_ipsi"])["mean"]
        try:
            rc, ri = mm[(0.9, 0.0)], mm[(0.0, 0.9)]
        except KeyError:
            continue
        if max(rc, ri) <= 0:
            continue
        odi[g.site.iloc[0]].append((rc - ri) / (abs(rc) + abs(ri)))
    bins = np.linspace(-1, 1, 21)
    hm = np.histogram(odi["mua"], bins)[0] / len(odi["mua"])
    hu = np.histogram(odi["unit"], bins)[0] / len(odi["unit"])
    ax.bar(bins[:-1], hm, width=np.diff(bins), align="edge", color=fs.RULE, lw=0)
    ax.stairs(hu, bins, color=fs.INK, lw=0.9)
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, 0, 1], ["−1", "0", "1"])
    ax.set_xlabel("ocular dominance index", labelpad=1)
    ax.set_ylabel("fraction")
    top = max(hm.max(), hu.max()) * 1.3
    ax.set_ylim(0, top)
    ax.text(-0.97, top * 0.62, "ipsilateral\neye only", fontsize=5.5, color=fs.INK2, ha="left", va="top",
            linespacing=0.95)
    ax.text(0.97, top * 0.62, "contralateral\neye only", fontsize=5.5, color=fs.INK2, ha="right", va="top",
            linespacing=0.95)
    ax.add_patch(mpatches.Rectangle((-0.95, top * 0.9), 0.14, top * 0.06, fc=fs.RULE, ec="none"))
    ax.text(-0.76, top * 0.93, f"MUA, n = {len(odi['mua'])}", fontsize=5.5, va="center", color=fs.INK)
    ax.plot([-0.95, -0.81], [top * 0.8] * 2, color=fs.INK, lw=0.9)
    ax.text(-0.76, top * 0.8, f"units, n = {len(odi['unit'])}", fontsize=5.5, va="center", color=fs.INK)
    return {k: np.asarray(v) for k, v in odi.items()}


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    k = list(pop["pens"]).index(pen)

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")                     # axes placed in mm: the in-panel keys stay inside their panels
    dax, kax, tax = (axes_mm(fig, *POS[n]) for n in ("matrix", "key", "trial"))
    design(dax, kax, tax)
    bax, cax, oax = (axes_mm(fig, *POS[n]) for n in ("example", "crf", "odi"))
    ex = example(bax, pop, k, meta.loc[pen, "monkey"])
    crf, ratio = crfs(cax, pop)
    odi = ocular_dominance(oax)
    for ax, letter in zip((dax, bax, cax, oax), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig1_setup"))

    print(f"(b) example: monkey {meta.loc[pen, 'monkey']}, {int(meta.loc[pen, 'n_channels'])} channels, "
          f"median eye weight {meta.loc[pen, 'a']:.2f}, onset {meta.loc[pen, 'onset']:.1f} ms after stimulus onset")
    for w, d in ex.items():
        print(f"    mean {w} after onset: " + "; ".join(f"{lab} {v:.2f}" for lab, v in d.items()))
    for lab, (m, lo, hi) in crf.items():
        print(f"(c) {lab}: " + "; ".join(f"{100 * c:g}% {mm:.3f} [{l_:.3f}, {h:.3f}]"
                                         for c, mm, l_, h in zip(CONTRASTS[1:], m, lo, hi)))
    print(f"    other / dominant eye at 90% per penetration: median {np.median(ratio):.2f} "
          f"(IQR {np.percentile(ratio, 25):.2f}-{np.percentile(ratio, 75):.2f}); median eye weight a over "
          f"penetrations {meta.a.median():.2f} (IQR {meta.a.quantile(.25):.2f}-{meta.a.quantile(.75):.2f}); "
          f"{len(meta)} penetrations, {int(meta.n_channels.sum())} channels (B {int((meta.monkey == 'B').sum())}, "
          f"J {int((meta.monkey == 'J').sum())})")
    for s, v in odi.items():
        print(f"(d) {s}: n {len(v)}, median ODI {np.median(v):+.2f}, median |ODI| {np.median(np.abs(v)):.2f}, "
              f"|ODI| < 0.5 {np.mean(np.abs(v) < 0.5):.0%}, contralateral-dominant {np.mean(v > 0):.0%}, "
              f"|ODI| > 0.9 {np.mean(np.abs(v) > 0.9):.0%}")


if __name__ == "__main__":
    main()
