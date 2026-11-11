"""Figure 4 (simplified set, docs/09 checks N, U and their v3 rerun Y): a delayed normalization pool, weighted toward
the other eye.

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_delayed_pool.py [sink]
(a) what the model does (illustration): drive, pool and response for both eyes at 45% and the dominant eye alone, with
    the median fitted weights, pool delay, exponent and semi-saturation of the conductance fits per depth group;
(b) held-out error, % less with a delayed pool than without, per penetration and depth group (superficial / input,
    0.5 mm above the input/deep boundary of laminar alignment v3 / deep, to 0.7 mm below; y_model_v3.py, both
    cross-validation directions summed); medians and interquartile ranges; check Y1 (N1 on v3);
(c) the other eye's weight in the drive (a) and in the pool (w), conductance fits, per penetration and depth group;
    medians and interquartile ranges; check Y4 (U6', w - a > 0, hierarchical).
`sink` as argument: the earlier version, groups above / at, within 0.2 mm / below the earliest sink (checks N, U6).
Inputs: <out_root>/binocular/y_N_groups.parquet; docs/qc/y_N.csv, y_U.csv (v4_N_groups.parquet, v4_N.csv, v4_U.csv
for `sink`).
Writes <out_root>/figures/fig4_delayed_pool.pdf and .png (fig4_delayed_pool_sink for `sink`); prints the statistics
shown.
"""

import sys
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
DEPTH = {
    "v3": dict(groups=("superficial", "input", "deep"), label={"superficial": "superficial", "input": "input",
                                                               "deep": "deep"},
               files=("y_N_groups.parquet", "y_N.csv", "y_U.csv"), u6="U6'"),
    "sink": dict(groups=("above", "at", "below"), label={"above": "above", "at": "at sink", "below": "below"},
                 files=("v4_N_groups.parquet", "v4_N.csv", "v4_U.csv"), u6="U6"),
}
POS = {"a1": (12, 46, 52, 8), "a2": (12, 32, 52, 8), "a3": (12, 11, 52, 14), "b": (90, 11, 30, 42),
       "c": (128, 11, 30, 42)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def rows_panel(ax, series, xlabel, xlim, groups, label, ref=None):
    rng = np.random.default_rng(8)
    out = []
    for k, g in enumerate(groups):
        for lab, (v_all, colr, dy, mfc) in series.items():
            v = v_all[g].dropna()
            ax.scatter(v.clip(*xlim), k + dy + rng.uniform(-0.06, 0.06, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9, zorder=2)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            out.append(dict(series=lab, grp=g, n=len(v), median=md))
    if ref is not None:
        ax.axvline(ref, color=fs.MUTED, lw=0.6)
    ax.set_yticks(range(3), [label[g] for g in groups])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    return pd.DataFrame(out)


def main():
    mode = "sink" if len(sys.argv) > 1 and sys.argv[1] == "sink" else "v3"
    cfg = DEPTH[mode]
    gl = (cfg["groups"], cfg["label"])
    fs.apply()
    f_groups, f_n, f_u = cfg["files"]
    N = pd.read_parquet(out_root() / "binocular" / f_groups)
    Tn = pd.read_csv(QC / f_n).set_index("test")
    U = pd.read_csv(QC / f_u).set_index("test")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    sk = timing_sketch([ax["a1"], ax["a2"], ax["a3"]], float(N.cond_p_a.median()), float(N.cond_p_w.median()),
                       float(N.cond_p_tau2.median()), n=float(N.cond_p_n.median()), sig=float(N.cond_p_sigma.median()))
    by = lambda col: N.pivot(index="penetration", columns="grp", values=col)    # noqa: E731
    N["pct_delay"] = 100 * (1 - N.cv_full / N.cv_static)
    B = rows_panel(ax["b"], {"delay": (by("pct_delay"), fs.INK, 0.0, fs.INK)},
                   "held-out error, % less\nwith a delayed pool", (-10, 100), *gl, ref=0)
    ax["b"].set_title("needed at every depth", fontsize=6.5)
    C = rows_panel(ax["c"], {"drive (a)": (by("cond_p_a"), fs.MUTED, 0.12, "white"),
                             "pool (w)": (by("cond_p_w"), fs.INK, -0.12, fs.INK)},
                   "the other eye's weight\n(dominant eye = 1)", (0, 2.0), *gl, ref=1.0)
    ax["c"].set_title("other eye weighted more in the pool", fontsize=6.5)
    ax["c"].tick_params(labelleft=False)
    ax["c"].text(0.97, 0.03, "filled: pool\nopen: drive", transform=ax["c"].transAxes, ha="right", va="bottom",
                 fontsize=5.5, color=fs.INK2, linespacing=1.0)
    for k, letter, x in (("a1", "a", -0.02), ("b", "b", -0.4), ("c", "c", -0.08)):
        fs.panel_label(ax[k], letter, x=x)
    print(fs.save(fig, "fig4_delayed_pool" + ("_sink" if mode == "sink" else "")))
    pd.set_option("display.width", 200)
    print("(a)", sk)
    print(B.round(3).to_string(index=False))
    print(C.round(3).to_string(index=False))
    print(Tn.loc[[t for t in Tn.index if t.startswith("N1")], ["n", "median", "n_in_direction", "p", "holds"]].round(4))
    print(U.loc[[t for t in U.index if t.startswith(cfg["u6"] + " ")], ["estimate", "p_one_sided", "holds"]].round(4))


if __name__ == "__main__":
    main()
