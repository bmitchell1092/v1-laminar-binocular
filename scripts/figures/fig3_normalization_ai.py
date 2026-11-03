"""Figure 3 (simplified set, docs/09 check V): binocular summation gives way to normalization within ~30 ms, more at
and above the input layer, and more when the other eye's contrast is higher.

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_normalization_ai.py
Additivity index AI = both eyes / (dominant eye + other eye); the dashed reference is the index that summation of the
two eyes' contrasts predicts: the dominant eye's response at c_D + a c_O over the same denominator (each channel's or
penetration's eye weight a). Time from each channel's response onset.
(a) AI over time at 45% in each eye, one trace per 0.2 mm depth bin relative to the earliest sink (population: mean
    over penetrations of the channel means; ratio of the means; >= 10 penetrations; 5 ms running mean); grey where
    the index falls below the summation reference; the thin rule under each trace is AI = 0.5;
(b) AI 30-100 ms after onset per penetration and depth group (above / at, within 0.2 mm / below the sink), data and
    summation reference; medians and interquartile ranges; check V2 (below - at, one-sided Wilcoxon);
(c, d) AI against the other eye's contrast with the dominant eye at 11%, 0-20 ms (c) and 130-200 ms (d) after
    onset; penetration means and 95% bootstrap CI; check V3 (all cells with the other eye's contrast higher).
Inputs: spacetime_data.py (spacetime.npz, spacetime_meta.parquet); docs/qc/v4_V.csv, v4_V_cells.csv, v4_V_groups.csv.
Writes <out_root>/figures/fig3_normalization_ai.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 104
BIN = 2
EDGES = np.arange(-9, 14, BIN)
MIN_PEN = 10
CS4 = (0.11, 0.225, 0.45, 0.9)
POS = {"a": (20, 10, 62, 86), "b": (104, 58, 40, 36), "c": (104, 10, 28, 32), "d": (138, 10, 28, 32)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def depth_ai(X, meta, c=0.45):
    i, j = list(CONTRASTS).index(c), CS4.index(c)
    D, O, B = X[:, i], X[:, 6 + j], X[:, 10 + j]
    P = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    b = np.digitize(meta.rel_sink.to_numpy(), EDGES) - 1
    centres = (EDGES[:-1] + BIN / 2) * 0.1
    pen = meta.penetration.to_numpy()
    out, z = {k: [] for k in "DOBP"}, []
    for k in range(len(centres)):
        sel = b == k
        pens = np.unique(pen[sel])
        if len(pens) < MIN_PEN:
            continue
        z.append(centres[k])
        for name, x in zip("DOBP", (D, O, B, P)):
            out[name].append(np.nanmean([np.nanmean(x[sel & (pen == p)], axis=0) for p in pens], axis=0))
    M = {k: fs.smooth(np.array(v)) for k, v in out.items()}
    den = M["D"] + M["O"]
    ok = den > 0.12
    return np.array(z), np.where(ok, M["B"] / den, np.nan), np.where(ok, M["P"] / den, np.nan)


def ridgeline(ax, rel, z, AI, AS, gap=0.5):
    sel = (rel >= 0) & (rel <= 150)
    t = rel[sel]
    for k in range(len(z))[::-1]:
        y0 = -k * gap
        a, p = AI[k, sel] + y0, AS[k, sel] + y0
        ax.fill_between(t, a, p, where=a < p, color=fs.INK2, alpha=0.35, lw=0, zorder=2 + 3 * k)
        ax.plot(t, p, color=fs.MUTED, lw=0.7, ls=fs.SUM_LS, zorder=3 + 3 * k)
        ax.plot(t, a, color=fs.INK, lw=0.9, zorder=3 + 3 * k)
        ax.plot(t, np.full(t.size, y0 + 0.5), color=fs.RULE, lw=0.4, zorder=0)
        ax.text(-4, y0 + 0.65, f"{z[k]:+.1f}", ha="right", va="center", fontsize=5.5, color=fs.INK2)
    ax.plot([152, 152], [0.5, 1.0], color=fs.INK, lw=0.8, clip_on=False)               # scale bar: 0.5 of the index
    ax.text(154, 0.75, "0.5", fontsize=5.5, va="center", color=fs.INK)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 150)
    ax.set_ylim(-gap * (len(z) - 1) + 0.35, 1.35)
    ax.set_xlabel("time from response onset (ms)")
    ax.text(-26, -gap * (len(z) - 1) / 2 + 0.6, "depth relative to earliest sink (mm)", rotation=90, ha="center",
            va="center", fontsize=7, color=fs.INK)
    ax.set_title("additivity index, both eyes 45%", fontsize=6.5, loc="left")
    fs.line_key(ax, [("measured", fs.INK, "-"), (fs.SUM, fs.MUTED, fs.SUM_LS)], x=0.5, y=1.0, dy=0.035, length=0.08)


def by_group(ax, G, V):
    rng = np.random.default_rng(7)
    D = G[G.window == "dip"]
    rows = []
    for k, g in enumerate(("above", "at", "below")):
        x = D[D.grp == g]
        for col, colr, dy, mfc in (("AI", fs.INK, -0.12, fs.INK), ("AI_sum", fs.MUTED, 0.12, "white")):
            v = x[col].dropna()
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            rows.append(dict(grp=g, series=col, n=len(v), median=md))
    ax.set_yticks(range(3), ["above", "at sink", "below"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(0.3, 0.9)
    ax.set_xlabel("additivity index, 30–100 ms")
    v2 = V.set_index("test").loc["V2 (AI - AI_sum) below - at the sink > 0, 45%, 30-100 ms"]
    ax.text(1.03, 0.98, f"measured\n{fs.SUM} (open)\n\nbelow − at sink:\n{int(v2.n_in_direction)}/{int(v2.n)} "
            f"penetrations,\n{fs.fmt_p(v2.p)}", transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=fs.INK2, linespacing=1.1)
    return pd.DataFrame(rows)


def by_other(axs, C, V):
    rows = []
    for ax, w, title in zip(axs, ("early", "late"), ("0–20 ms after onset", "130–200 ms")):
        X = C[np.isclose(C.c_dom, 0.11) & (C.window == w)]
        for col, colr, ls, mfc in (("AI", fs.INK, "-", fs.INK), ("AI_sum", fs.MUTED, fs.SUM_LS, "white")):
            P = X.pivot(index="penetration", columns="c_oth", values=col).reindex(columns=list(CS4))
            m, lo, hi = fs.boot_ci(P.to_numpy())
            ax.errorbar(CS4, m, [m - lo, hi - m], color=colr, ls=ls, lw=1.0, marker="o", ms=3.0, mfc=mfc, mec=colr,
                        elinewidth=0.6, capsize=0)
            rows += [dict(window=w, series=col, c_oth=c, mean=mm) for c, mm in zip(CS4, m)]
        fs.log_contrast_axis(ax, CS4, label="")
        ax.set_xticklabels(["11", "22.5", "45", "90"])
        ax.set_xlim(0.095, 1.05)
        ax.set_ylim(0.35, 1.0)
        ax.set_title(title, fontsize=6.5)
    axs[0].set_ylabel("additivity index")
    axs[1].tick_params(labelleft=False)
    axs[0].text(1.12, -0.3, "other eye's contrast (%), dominant eye 11%", transform=axs[0].transAxes, ha="center",
                fontsize=7)
    v3 = V.set_index("test").loc["V3 other eye stronger, 130-200 ms: AI - AI_sum < 0"]
    axs[1].text(0.97, 0.97, f"other eye stronger,\nall cells: {int(v3.n_in_direction)}/{int(v3.n)},\n{fs.fmt_p(v3.p)}",
                transform=axs[1].transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK2, linespacing=1.1)
    return pd.DataFrame(rows)


def main():
    fs.apply()
    Z = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    V = pd.read_csv(QC / "v4_V.csv")
    C, G = pd.read_csv(QC / "v4_V_cells.csv"), pd.read_csv(QC / "v4_V_groups.csv")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    z, AI, AS = depth_ai(Z["onset"], meta)
    ridgeline(ax["a"], Z["rel"], z, AI, AS)
    B = by_group(ax["b"], G, V)
    O = by_other([ax["c"], ax["d"]], C, V)
    for k, x in (("a", -0.2), ("b", -0.32), ("c", -0.45)):
        fs.panel_label(ax[k], k, x=x)
    print(fs.save(fig, "fig3_normalization_ai"))
    pd.set_option("display.width", 200)
    rel = Z["rel"]
    for k in range(len(z)):
        e, d = (rel >= 0) & (rel < 20), (rel >= 30) & (rel < 100)
        print(f"(a) {z[k]:+.1f} mm: 0-20 ms {np.nanmean(AI[k, e]):.2f} (sum {np.nanmean(AS[k, e]):.2f}); "
              f"30-100 ms {np.nanmean(AI[k, d]):.2f} (sum {np.nanmean(AS[k, d]):.2f})")
    print(B.round(3).to_string(index=False))
    print(O.round(3).pivot_table(index=["window", "series"], columns="c_oth", values="mean").to_string())
    print(V.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
