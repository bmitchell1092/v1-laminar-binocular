"""Figure 4 (simplified set, docs/09 checks K, N, P, U): a delayed normalization pool, weighted toward the other eye.

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_delayed_pool.py
(a) what the model does (illustration): drive, pool and response for both eyes at 45% and the dominant eye alone, with
    the median fitted weights, pool delay, exponent and semi-saturation of the conductance fits per depth group;
(b) held-out error, % less with a delayed pool than without, per penetration and depth group (above / at, within
    0.2 mm / below the earliest sink; v4_model_groups.py, both cross-validation directions summed); medians and
    interquartile ranges; check N1;
(c) the other eye's weight in the drive (a) and in the pool (w), conductance fits, per penetration and depth group;
    medians and interquartile ranges; check U6 (w - a > 0, hierarchical).
Inputs: <out_root>/binocular/v4_N_groups.parquet; docs/qc/v4_N.csv, v4_U.csv.
Writes <out_root>/figures/fig4_delayed_pool.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig4_model import timing_sketch  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 66
GROUPS = ("above", "at", "below")
LABEL = {"above": "above", "at": "at sink", "below": "below"}
POS = {"a1": (12, 46, 52, 8), "a2": (12, 32, 52, 8), "a3": (12, 11, 52, 14), "b": (90, 11, 30, 42),
       "c": (128, 11, 30, 42)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def rows_panel(ax, series, xlabel, xlim, ref=None):
    rng = np.random.default_rng(8)
    out = []
    for k, g in enumerate(GROUPS):
        for lab, (v_all, colr, dy, mfc) in series.items():
            v = v_all[g].dropna()
            ax.scatter(v.clip(*xlim), k + dy + rng.uniform(-0.06, 0.06, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9, zorder=2)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            out.append(dict(series=lab, grp=g, n=len(v), median=md))
    if ref is not None:
        ax.axvline(ref, color=fs.MUTED, lw=0.6)
    ax.set_yticks(range(3), [LABEL[g] for g in GROUPS])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    return pd.DataFrame(out)


def main():
    fs.apply()
    N = pd.read_parquet(out_root() / "binocular" / "v4_N_groups.parquet")
    Tn = pd.read_csv(QC / "v4_N.csv").set_index("test")
    U = pd.read_csv(QC / "v4_U.csv").set_index("test")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    sk = timing_sketch([ax["a1"], ax["a2"], ax["a3"]], float(N.cond_p_a.median()), float(N.cond_p_w.median()),
                       float(N.cond_p_tau2.median()), n=float(N.cond_p_n.median()), sig=float(N.cond_p_sigma.median()))
    by = lambda col: N.pivot(index="penetration", columns="grp", values=col)    # noqa: E731
    N["pct_delay"] = 100 * (1 - N.cv_full / N.cv_static)
    B = rows_panel(ax["b"], {"delay": (by("pct_delay"), fs.INK, 0.0, fs.INK)},
                   "held-out error, % less\nwith a delayed pool", (-10, 100), ref=0)
    ax["b"].set_title("needed at every depth", fontsize=6.5)
    C = rows_panel(ax["c"], {"drive (a)": (by("cond_p_a"), fs.MUTED, 0.12, "white"),
                             "pool (w)": (by("cond_p_w"), fs.INK, -0.12, fs.INK)},
                   "the other eye's weight\n(dominant eye = 1)", (0, 2.0), ref=1.0)
    ax["c"].set_title("other eye weighted more in the pool", fontsize=6.5)
    ax["c"].tick_params(labelleft=False)
    ax["c"].text(0.97, 0.03, "filled: pool\nopen: drive", transform=ax["c"].transAxes, ha="right", va="bottom",
                 fontsize=5.5, color=fs.INK2, linespacing=1.0)
    for k, letter, x in (("a1", "a", -0.02), ("b", "b", -0.4), ("c", "c", -0.08)):
        fs.panel_label(ax[k], letter, x=x)
    print(fs.save(fig, "fig4_delayed_pool"))
    pd.set_option("display.width", 200)
    print("(a)", sk)
    print(B.round(3).to_string(index=False))
    print(C.round(3).to_string(index=False))
    print(Tn.loc[[t for t in Tn.index if t.startswith("N1")], ["n", "median", "n_in_direction", "p", "holds"]].round(4))
    print(U.loc[[t for t in U.index if t.startswith("U6")], ["estimate", "p_one_sided", "holds"]].round(4))


if __name__ == "__main__":
    main()
