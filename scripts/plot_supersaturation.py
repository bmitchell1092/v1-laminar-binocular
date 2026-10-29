"""Figure for the supersaturation assessment: <out_root>/crf/supersaturation_assessment.png

Usage: .venv\\Scripts\\python scripts\\plot_supersaturation.py
Row 1: cross-validated SI vs depth below the top of cortex (sustained window), per condition; MUA sites and single
units (strict / lenient); MUA running median. Row 2: fraction of driven sites with a held-out downturn
(p < 0.05) against the 5% expected without supersaturation; normalised CRFs of single units with a downturn.
"""

import glob

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v1lb import crf  # noqa: E402
from v1lb.config import out_root  # noqa: E402

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
COL = {"mua": "#2a78d6", "unit_lenient": "#eb6834", "unit_strict": "#1baf7a"}     # categorical slots 1-3
NAME = {"mua": "MUA", "unit_lenient": "units, lenient", "unit_strict": "units, strict"}
CONDS = ["mono_contra", "mono_ipsi", "dioptic"]


def load():
    R = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(out_root() / "crf" / "*_crf.parquet")))],
                  ignore_index=True)
    strict = R.well_isolated.fillna(False).astype(bool)
    lenient = (R.site == "unit") & (R.isi_contamination < 0.5) & (R.snr >= 3)
    R["kind"] = np.where(R.site == "mua", "mua", np.where(strict, "unit_strict", np.where(lenient, "unit_lenient", "")))
    D = R[R.driven & (R.kind != "")]
    return pd.concat([D, D[D.kind == "unit_strict"].assign(kind="unit_lenient")], ignore_index=True)


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def main():
    D = load()
    fig, axs = plt.subplots(2, 3, figsize=(14, 9), facecolor=SURFACE)
    for ax, cond in zip(axs[0], CONDS):
        style(ax)
        g = D[(D.window == "sustained") & (D.condition == cond)]
        for kind in ("mua", "unit_lenient", "unit_strict"):
            k = g[(g.kind == kind) & g.depth_below_toc.notna()]
            x = k.si_cv.clip(-1, 1)
            ax.scatter(x, k.depth_below_toc + (0.15 if kind != "mua" else 0), s=10 if kind == "mua" else 26,
                       color=COL[kind], alpha=0.35 if kind == "mua" else 0.9, lw=0,
                       label=f"{NAME[kind]} (n={len(k)})", zorder=2 if kind == "mua" else 3)
        m = g[g.kind == "mua"].groupby(g.depth_below_toc // 2 * 2).si_cv.median()
        ax.plot(m.values, m.index + 1, color=COL["mua"], lw=2, zorder=4)
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_ylim(32, -3)
        ax.set_xlim(-1.05, 1.05)
        ax.set_title(f"{cond.replace('_', ' ')}, sustained (150-250 ms)", fontsize=10, color=INK)
        ax.set_xlabel("cross-validated supersaturation index (clipped at ±1)", fontsize=8, color=INK2)
        ax.set_ylabel("contacts below top of cortex (100 µm)", fontsize=8, color=INK2)
        ax.legend(fontsize=7, frameon=False, loc="lower left")

    ax = axs[1, 0]
    style(ax)
    rows = []
    for (kind, window, cond), g in D.groupby(["kind", "window", "condition"]):
        rows.append(dict(kind=kind, window=window, cond=cond, frac=(g.p_downturn < 0.05).mean(), n=len(g)))
    F = pd.DataFrame(rows)
    labels = [f"{w[:5]}\n{c.replace('mono_', 'm-').replace('dioptic', 'diop')}" for w in ("transient", "sustained") for c in CONDS]
    xs = np.arange(len(labels))
    for i, kind in enumerate(("mua", "unit_lenient", "unit_strict")):
        f = F[F.kind == kind].set_index(["window", "cond"])
        vals = [f.frac.get((w, c), np.nan) for w in ("transient", "sustained") for c in CONDS]
        ax.bar(xs + (i - 1) * 0.27, vals, width=0.25, color=COL[kind], label=NAME[kind])
    ax.axhline(0.05, color=INK, lw=1, ls="--")
    ax.text(len(labels) - 0.5, 0.055, "5% expected if no site supersaturates", fontsize=7, color=INK2, ha="right", va="bottom")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("fraction of driven sites with held-out downturn (p < .05)", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("how often responses fall at high contrast", fontsize=10, color=INK)

    for ax, kind in zip(axs[1, 1:], ("unit_strict", "mua")):
        style(ax)
        g = D[(D.kind == kind) & (D.p_downturn < 0.05)].sort_values("p_downturn").drop_duplicates(["penetration", "site_id"])
        g = g.head(12)
        cols = [f"r{c:g}" for c in crf.CONTRASTS]
        for _, row in g.iterrows():
            v = row[cols].to_numpy(float)
            ax.plot(np.arange(6), v / np.nanmax(v), color=COL[kind], lw=1.2, alpha=0.8, marker="o", ms=3)
        ax.set_xticks(range(6))
        ax.set_xticklabels([f"{c:g}" for c in crf.CONTRASTS], fontsize=7)
        ax.set_xlabel("contrast", fontsize=8, color=INK2)
        ax.set_ylabel("response / peak (blank-subtracted)", fontsize=8, color=INK2)
        ax.set_title(f"{NAME[kind]} with a held-out downturn (up to 12 shown, n={len(g)})", fontsize=10, color=INK)
    fig.suptitle("Supersaturation in V1: cross-validated statistics, 26 penetrations (B: 8, J: 18)", fontsize=12, color=INK)
    fig.tight_layout()
    fn = out_root() / "crf" / "supersaturation_assessment.png"
    fig.savefig(fn, dpi=90, facecolor=SURFACE)
    print(fn)


if __name__ == "__main__":
    main()
