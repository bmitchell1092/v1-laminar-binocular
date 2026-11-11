"""Figure 2 (docs/09 checks Q and R): the eyes' contrasts add at onset; afterwards the binocular response falls short of
summation in proportion to the other eye's contrast, and stays short when the other eye is the stronger one.

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_dichoptic_time.py
Shortfall from summation = both eyes - "if the contrasts add" (the dominant eye's response at c_D + a c_O, the
penetration's median eye weight a, interpolated in log contrast; v1lb.binocular.response_at), with each channel
aligned to its own response onset (half-max, dominant eye at 90%) and scaled to its dominant-eye peak at 90%;
penetration means of channel means. Cells with both contrasts >= 11% and c_D + a c_O <= 90% (check Q).
(a) one cell of the design, dominant eye 22.5% and other eye 90%: both eyes, "if the contrasts add" and the dominant
    eye alone at 22.5%, means over the penetrations in which summation can be predicted for this cell; orange where
    both eyes exceed the prediction, grey where they fall short; 5 ms running mean for display; the three analysis
    windows marked;
(b) the shortfall in each window against the other eye's contrast, mean over the dominant-eye contrasts per
    penetration, then mean over penetrations with 95% bootstrap CI; the slope on the other eye's contrast from the
    hierarchical model (check Q1, Q2) per window, and the held-out confirmation (check R2);
(c) the design as a matrix in each window: rows the dominant eye's contrast, columns the other eye's; mean over
    penetrations (v4_Q_cells.csv), in % of the dominant-eye peak; outlined, matched contrasts; hatched, cells where
    summation can be predicted in fewer than 10 penetrations.
Statistics: checks Q (all 26 penetrations, hierarchical) and R (eight held-out cross-recordings).
Inputs: <out_root>/figures/data/population.npz, population_meta.parquet; docs/qc/v4_Q_cells.csv, v4_Q.csv, v4_R.csv.
Writes <out_root>/figures/fig2_dichoptic_time.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 100
CD = (0.11, 0.225, 0.45)
CO = (0.11, 0.225, 0.45, 0.9)
EXAMPLE = (0.225, 0.9)                                                   # (dominant, other) for panel a
WIN = {"early": ("0–20 ms", (0, 20)), "dip": ("30–100 ms", (30, 100)), "late": ("130–200 ms", (130, 200))}
WIN_STYLE = {"early": (fs.MUTED, "o", fs.MUTED), "dip": (fs.INK2, "s", fs.INK2), "late": (fs.INK, "D", "white")}
MIN_PEN = 10
# grey (below summation) - off-white - orange (above), as figstyle.ADDED_CMAP but with an off-white midpoint
DIV = LinearSegmentedColormap.from_list("shortfall", ["#3d3c39", "#8a8985", "#f4f3f0", "#f2a674", "#c94a14"])
LIM = 6.0                                                                # colour limit, % of the dominant-eye peak
POS = {"a": (16, 58, 92, 34), "b": (128, 58, 40, 34),
       "c0": (16, 10, 30, 27), "c1": (52, 10, 30, 27), "c2": (88, 10, 30, 27), "cb": (124, 10, 2, 27)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def example_cell(ax, pop, meta, cd, co):
    """Both eyes, the summation prediction and the dominant eye alone for one cell, over the penetrations in which
    summation can be predicted (c_D + a c_O <= 90%)."""
    rel, O = pop["rel"], pop["onset"]
    ci = {c: list(CONTRASTS).index(c) for c in CONTRASTS}
    B, P, D = [], [], []
    for k, p in enumerate(pop["pens"]):
        a = float(meta.loc[str(p), "a"])
        if cd + a * co > CONTRASTS[-1]:
            continue
        B.append(O[k, ci[cd], ci[co]])
        P.append(response_at(O[k, :, 0], cd + a * co))
        D.append(O[k, ci[cd], 0])
    sel = (rel >= -10) & (rel <= 200)
    x = rel[sel]
    b, p, d = (fs.smooth(np.nanmean(np.array(Y), axis=0))[sel] for Y in (B, P, D))
    ax.fill_between(x, b, p, where=b >= p, color=fs.ORANGE, alpha=0.55, lw=0, interpolate=True, zorder=1)
    ax.fill_between(x, b, p, where=b < p, color=fs.MUTED, alpha=0.55, lw=0, interpolate=True, zorder=1)
    ax.plot(x, d, color=fs.COND["dominant"], lw=0.9, zorder=2)
    ax.plot(x, p, color=fs.MUTED, ls=fs.SUM_LS, lw=1.0, zorder=3)
    ax.plot(x, b, color=fs.INK, lw=1.1, zorder=4)
    top = 1.02 * max(b.max(), p.max())
    for key, (lab, (w0, w1)) in WIN.items():
        ax.axvspan(w0, w1, color=fs.RULE, alpha=0.3, lw=0, zorder=0)
        ax.text((w0 + w1) / 2, top * 1.06, lab, ha="center", va="bottom", fontsize=5.5, color=fs.INK2)
    fs.zero_line(ax, color=fs.RULE)
    ax.set_xlim(-10, 200)
    ax.set_ylim(-0.03, top * 1.2)
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("MUA (dominant-eye peak = 1)")
    ax.set_title(f"one cell: dominant eye {fs.contrast_label(cd)}, other eye {fs.contrast_label(co)}",
                 fontsize=6.5, loc="left", pad=11)
    fs.line_key(ax, [("both eyes, measured", fs.INK, "-"), (fs.SUM, fs.MUTED, fs.SUM_LS),
                     (f"dominant eye alone, {fs.contrast_label(cd)}", fs.COND["dominant"], "-")],
                x=0.5, y=0.76, dy=0.1, length=0.06)
    i0 = int(np.argmax(b))
    ax.annotate("more than\nsummation", (x[i0] + 4, (b[i0] + p[i0]) / 2), xytext=(x[i0] + 16, b[i0] * 0.98),
                fontsize=5.5, color=fs.INK2, va="center", linespacing=1.05,
                arrowprops=dict(arrowstyle="-", color=fs.INK2, lw=0.5, shrinkA=1, shrinkB=0))
    ax.text(112, 0.13, "less than summation", fontsize=5.5, color=fs.INK2, va="top")
    ax.annotate("", (118, 0.27), xytext=(118, 0.14),
                arrowprops=dict(arrowstyle="-", color=fs.INK2, lw=0.5, shrinkA=0, shrinkB=0))
    win = {k: (float(np.mean((b - p)[(x >= w0) & (x < w1)]))) for k, (_, (w0, w1)) in WIN.items()}
    return dict(n=len(B), peak_both=float(b.max()), peak_dom=float(d.max()), **{f"S_{k}": v for k, v in win.items()})


def by_other(ax, Q, T, Tr):
    """Window means against the other eye's contrast (mean over the dominant-eye contrasts per penetration)."""
    rows = []
    for k, (key, (lab, _)) in enumerate(WIN.items()):
        V = (Q[Q.window == key].groupby(["unit", "c_oth"]).S.mean().unstack().reindex(columns=list(CO))) * 100
        m, lo, hi = fs.boot_ci(V.to_numpy())
        col, mk, mfc = WIN_STYLE[key]
        x = np.arange(len(CO)) + (k - 1) * 0.09
        ax.errorbar(x, m, [m - lo, hi - m], color=col, marker=mk, ms=2.8, mfc=mfc, mec=col, mew=0.7, lw=1.0,
                    elinewidth=0.6, capsize=0, zorder=3)
        slope = T.loc[f"Q x_O slope, {key}"]
        ax.text(3.3, m[-1] + LABEL_DY[key], f"{lab}\n{100 * slope.estimate:+.1f} per doubling".replace("-", "−"),
                fontsize=5.2, color=fs.INK2, va="center", linespacing=1.05)
        rows += [dict(window=key, c_oth=c, n=int(V[c].notna().sum()), mean=mm, lo=l_, hi=h)
                 for c, mm, l_, h in zip(CO, m, lo, hi)]
    fs.zero_line(ax, color=fs.MUTED)
    ax.set_ylim(-12, 11)
    ax.set_xticks(range(len(CO)), [f"{100 * c:g}" for c in CO])
    ax.set_xlim(-0.35, 3.25)
    ax.set_xlabel("other eye's contrast (%)")
    ax.set_ylabel("both eyes − if the contrasts add\n(% of dominant-eye peak)")
    ax.set_title("all cells, by the other eye's contrast", fontsize=6.5, loc="left", pad=11)
    q, r = T.loc["Q late shortfall, other eye stronger (c_O > c_D)"], \
        Tr.loc["R late shortfall, other eye stronger (c_O > c_D)"]
    ax.text(0.03, 0.03, (f"130–200 ms, other eye stronger: {100 * q.estimate:+.1f}, {fs.fmt_p(q.p_one_sided)}\n"
                         f"held-out recordings: {100 * r.estimate:+.1f}, {fs.fmt_p(r.p_one_sided)}").replace("-", "−"),
            transform=ax.transAxes, fontsize=5.2, color=fs.INK2, va="bottom", linespacing=1.1)
    return pd.DataFrame(rows)


LABEL_DY = {"early": 0.0, "dip": -1.6, "late": 1.6}


def matrices(axs, cax, Q):
    norm = TwoSlopeNorm(0, -LIM, LIM)
    rows = []
    for ax, (key, (lab, _)) in zip(axs, WIN.items()):
        W = Q[Q.window == key]
        M = np.full((len(CD), len(CO)), np.nan)
        for i, cd in enumerate(CD):
            for j, co in enumerate(CO):
                v = W[np.isclose(W.c_dom, cd) & np.isclose(W.c_oth, co)].S * 100
                rows.append(dict(window=key, c_dom=cd, c_oth=co, n=len(v), mean=v.mean(), n_below=int((v < 0).sum())))
                if len(v) >= MIN_PEN:
                    M[i, j] = v.mean()
        im = ax.imshow(M, cmap=DIV, norm=norm, origin="lower", aspect="auto")
        for i, cd in enumerate(CD):
            for j, co in enumerate(CO):
                if np.isnan(M[i, j]):
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="//////", ec=fs.RULE,
                                               lw=0))
                    continue
                ax.text(j, i, f"{M[i, j]:+.1f}".replace("-", "−"), ha="center", va="center", fontsize=5,
                        color="white" if abs(M[i, j]) > 0.62 * LIM else fs.INK)
                if cd == co:
                    ax.add_patch(plt.Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fill=False, ec=fs.INK, lw=0.7))
        ax.set_xticks(range(len(CO)), [f"{100 * c:g}" for c in CO])
        ax.set_yticks(range(len(CD)), [f"{100 * c:g}" for c in CD])
        ax.set_title(f"{lab} after onset", fontsize=6.5)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    axs[0].set_ylabel("dominant eye's\ncontrast (%)")
    axs[1].set_xlabel("other eye's contrast (%)")
    for ax in axs[1:]:
        ax.tick_params(labelleft=False)
    cb = plt.colorbar(im, cax=cax, ticks=[-LIM, -LIM / 2, 0, LIM / 2, LIM])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.ax.set_yticklabels([f"{t:+g}".replace("-", "−") if t else "0" for t in cb.get_ticks()])
    cb.set_label("both eyes − if the contrasts add\n(% of dominant-eye peak)", fontsize=6)
    cax.text(14, 1.0, "outlined: matched\ncontrasts\n\nhatched: summation\npredictable in fewer\nthan 10 penetrations",
             fontsize=5.2, color=fs.INK2, transform=cax.transAxes, va="top", ha="left", linespacing=1.1)
    return pd.DataFrame(rows)


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    Q = pd.read_csv(QC / "v4_Q_cells.csv")
    Q = Q[Q.split == "all"]
    T = pd.read_csv(QC / "v4_Q.csv").set_index("test")
    Tr = pd.read_csv(QC / "v4_R.csv").set_index("test")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    ex = example_cell(ax["a"], pop, meta, *EXAMPLE)
    B = by_other(ax["b"], Q, T, Tr)
    C = matrices([ax["c0"], ax["c1"], ax["c2"]], ax["cb"], Q)
    fs.panel_label(ax["a"], "a", x=-0.13, y=1.12)
    fs.panel_label(ax["b"], "b", x=-0.42, y=1.12)
    fs.panel_label(ax["c0"], "c", x=-0.55, y=1.15)
    print(fs.save(fig, "fig2_dichoptic_time"))
    pd.set_option("display.width", 220)
    print(f"(a) dominant {EXAMPLE[0]}, other {EXAMPLE[1]}: {ex}")
    print("(b) window means by the other eye's contrast (% of dominant-eye peak):")
    print(B.round(2).to_string(index=False))
    print("(c) cell means (% of dominant-eye peak), n penetrations, n below summation:")
    print(C.round(2).pivot_table(index=["window", "c_dom"], columns="c_oth", values="mean").to_string())
    print(C.pivot_table(index=["window", "c_dom"], columns="c_oth", values=["n", "n_below"]).to_string())
    print(T.round(4).to_string())
    print(Tr.round(4).to_string())


if __name__ == "__main__":
    main()
