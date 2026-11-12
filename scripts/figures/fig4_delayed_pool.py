"""Figure 4: what the normalization requires: a delayed pool, weighted toward the other eye, and what a single-stage
pool misses across depth (docs/09 checks Y and AA; redesigned 2026-11-12 after Blake's V5 comment that the held-out
error panel did not show what the delay does).

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_delayed_pool.py
(a) what the model does (illustration): drive, pool and response for both eyes at 45% and the dominant eye alone, with
    the median fitted weights, pool delay, exponent and semi-saturation of the conductance fits per depth group;
(b) why a delay: both eyes - contrast summation at 45% in each eye over time from response onset, for the data, the
    conductance fit and the fit without a delay (static pool); each array on its own onset and with its own eye weight
    and contrast-summation prediction (aa_model_misses_v3.deficit); per penetration the mean over its depth groups,
    then mean over penetrations; 5 ms running mean; orange where the data exceed contrast summation, grey where
    they fall short (as Figs 2 and 3);
(c) why the pool weights the other eye more: the same difference 30-100 ms after onset against the other eye's
    contrast (mean over dominant-eye contrasts 11, 22.5 and 45% where contrast summation can be predicted, as Fig 2b),
    for the data, the conductance fit and the fit with equal weights in drive and pool (w = a);
(d) the other eye's weight in the drive (a) and in the pool (w), conductance fits, per penetration and depth group;
    check Y4 (U6', w - a > 0, hierarchical);
(e) the shortfall at 45%, 30-100 ms, per depth group: data, the conductance fit and the two-stage fit (check AA,
    docs/qc/aa_groups.csv), with the deep - input tests (aa_model_misses.csv).
The static and w = a fits are the fixed-delay model without its delay or with w = a; at this scale the fixed-delay and
conductance fits are indistinguishable (printed).
Inputs: <out_root>/binocular/y_N_pred.npz, y_N_groups.parquet; docs/qc/aa_groups.csv, aa_model_misses.csv.
Writes <out_root>/figures/fig4_delayed_pool.pdf and .png; prints the statistics shown.
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
from v1lb.binocular import CONTRASTS, eye_weight, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 118
GROUPS = ("superficial", "input", "deep")
CD = (0.11, 0.225, 0.45)
CO = (0.11, 0.225, 0.45, 0.9)
REL = np.arange(-20, 181)                                   # ms from each array's own response onset
FITS = {"data": ("data", fs.INK, "-"), "conductance": ("cond", fs.VIOLET, "-"),
        "no delay (static pool)": ("static", fs.INK2, (0, (2.2, 1.2))),
        "equal weights (w = a)": ("w_eq_a", fs.AQUA, (0, (2.2, 1.2))), "fixed delay": ("full", fs.AQUA, "-")}
TWO_STAGE = fs.MAGENTA
POS = {"a1": (12, 102, 52, 7), "a2": (12, 89, 52, 7), "a3": (12, 68, 52, 13), "b": (88, 68, 80, 40),
       "c": (14, 11, 42, 38), "d": (82, 11, 34, 38), "e": (122, 11, 46, 38)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def measures(Y, T, cd, cn):
    """One group's array [condition, time] -> its onset, eye weight, both eyes - contrast summation at 45% in each eye
    on the onset clock (REL), and the 30-100 ms difference per (c_D, c_O) cell where summation can be predicted."""
    C = list(CONTRASTS)
    f = np.zeros((6, 6, T.size))
    for k, (x, y) in enumerate(zip(cd, cn)):
        f[C.index(x), C.index(y)] = Y[k]
    dom, oth = f[:, 0], f[0, :]
    j = np.flatnonzero((T >= 10) & (dom[-1] >= 0.5 * dom[-1].max()))
    on = T[j[0]]
    wt = T >= 40
    a = eye_weight(dom[:, wt].mean(axis=1), oth[:, wt].mean(axis=1))
    i = C.index(0.45)
    S45 = np.interp(REL + on, T, f[i, i] - response_at(dom, 0.45 * (1 + a)), left=np.nan, right=np.nan)
    sel = (T >= on + 30) & (T < on + 100)
    cells = {(x, y): float(np.mean(f[C.index(x), C.index(y), sel] - response_at(dom, x + a * y)[sel]))
             for x in CD for y in CO if x + a * y <= 0.9}
    return on, a, S45, cells


def collect():
    Z = np.load(out_root() / "binocular" / "y_N_pred.npz", allow_pickle=True)
    cd, cn, T = Z["cd"], Z["cn"], Z["t"]
    curves, cells = {}, []
    for name, (key, _, _) in FITS.items():
        per = {}
        for q, pen in enumerate(Z["penetration"]):
            on, a, S45, cc = measures(Z[key][q], T, cd, cn)
            per.setdefault(pen, []).append(S45)
            cells += [dict(fit=name, penetration=pen, c_dom=x, c_oth=y, S=v) for (x, y), v in cc.items()]
        curves[name] = np.array([np.nanmean(v, axis=0) for v in per.values()]) * 100     # [penetration, REL], %
    Cells = pd.DataFrame(cells)
    Cells["S"] *= 100
    return curves, Cells


def over_time(ax, curves):
    sm = lambda y: fs.smooth(y[None])[0]                                             # noqa: E731
    m = sm(np.nanmean(curves["data"], axis=0))
    ax.fill_between(REL, m, 0, where=m >= 0, color=fs.ORANGE, alpha=0.55, lw=0, interpolate=True, zorder=1)
    ax.fill_between(REL, m, 0, where=m < 0, color=fs.MUTED, alpha=0.55, lw=0, interpolate=True, zorder=1)
    ax.plot(REL, m, color=fs.INK, lw=1.2, zorder=4)
    for name in ("conductance", "no delay (static pool)"):
        _, col, ls = FITS[name]
        ax.plot(REL, sm(np.nanmean(curves[name], axis=0)), color=col, ls=ls, lw=1.1, zorder=3)
    fs.zero_line(ax, color=fs.MUTED)
    ax.axvspan(30, 100, color=fs.RULE, alpha=0.3, lw=0, zorder=0)
    ax.set_xlim(0, 150)
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("both eyes − contrast summation\n(% of dominant-eye peak)")
    ax.set_title("why a delay: both eyes at 45%", fontsize=6.5, loc="left")
    fs.line_key(ax, [("data", fs.INK, "-"), ("conductance fit (delayed pool)", fs.VIOLET, "-"),
                     ("fit without a delay (static pool)", fs.INK2, FITS["no delay (static pool)"][2])],
                x=0.42, y=0.95, dy=0.09, length=0.07)
    ax.text(148, 0.4, "contrast summation", fontsize=5.5, color=fs.INK2, va="bottom", ha="right")
    w = lambda Y, lo_, hi_: np.nanmean(Y[:, (REL >= lo_) & (REL < hi_)], axis=1)       # noqa: E731
    return pd.DataFrame([dict(fit=n, early_0_20=np.median(w(curves[n], 0, 20)), dip_30_100=np.median(w(curves[n], 30, 100)),
                              late_130_180=np.median(w(curves[n], 130, 180)), n=len(curves[n])) for n in curves])


def by_other(ax, Cells):
    rows = []
    for k, name in enumerate(("data", "conductance", "equal weights (w = a)")):
        _, col, ls = FITS[name]
        V = Cells[Cells.fit == name].groupby(["penetration", "c_oth"]).S.mean().unstack().reindex(columns=list(CO))
        m, lo, hi = fs.boot_ci(V.to_numpy())
        x = np.arange(len(CO)) + (k - 1) * 0.08
        if name == "data":
            ax.errorbar(x, m, [m - lo, hi - m], color=col, marker="o", ms=3, lw=1.1, elinewidth=0.6, capsize=0, zorder=4)
        else:
            ax.plot(x, m, color=col, ls=ls, marker="o", ms=2.6, mfc="white", mec=col, mew=0.8, lw=1.0, zorder=3)
        rows += [dict(fit=name, c_oth=c, n=int(V[c].notna().sum()), mean=v) for c, v in zip(CO, m)]
    fs.zero_line(ax, color=fs.MUTED)
    ax.set_xticks(range(len(CO)), [f"{100 * c:g}" for c in CO])
    ax.set_xlim(-0.35, 3.35)
    ax.set_xlabel("other eye's contrast (%)")
    ax.set_ylabel("both eyes − contrast summation,\n30–100 ms (% of dominant-eye peak)")
    ax.set_title("why the pool weights the other eye more", fontsize=6.5, loc="left")
    fs.line_key(ax, [("data", fs.INK, "-", "o"), ("conductance fit", fs.VIOLET, "-", "o", "white"),
                     ("equal weights in drive\nand pool (w = a)", fs.AQUA, FITS["equal weights (w = a)"][2], "o", "white")],
                x=0.05, y=0.30, dy=0.11, length=0.08)
    return pd.DataFrame(rows)


def rows_panel(ax, series, groups, xlim, ref=None, ylim=(3.5, -1.4)):
    rng = np.random.default_rng(8)
    out = []
    for k, g in enumerate(groups):
        for lab, (v_all, colr, dy, mfc, mk) in series.items():
            v = v_all[g].dropna()
            ax.scatter(v.clip(*xlim), k + dy + rng.uniform(-0.05, 0.05, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9, zorder=2)
            ax.plot([md], [k + dy], marker=mk, ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            out.append(dict(series=lab, grp=g, n=len(v), median=md))
    if ref is not None:
        ax.axvline(ref, color=fs.MUTED, lw=0.6)
    ax.set_yticks(range(len(groups)), groups)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    return pd.DataFrame(out)


def main():
    fs.apply()
    N = pd.read_parquet(out_root() / "binocular" / "y_N_groups.parquet")
    G = pd.read_csv(QC / "aa_groups.csv")
    A = pd.read_csv(QC / "aa_model_misses.csv").set_index("test")
    curves, Cells = collect()
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    sk = timing_sketch([ax["a1"], ax["a2"], ax["a3"]], float(N.cond_p_a.median()), float(N.cond_p_w.median()),
                       float(N.cond_p_tau2.median()), n=float(N.cond_p_n.median()), sig=float(N.cond_p_sigma.median()))
    Bt = over_time(ax["b"], curves)
    Co = by_other(ax["c"], Cells)

    by = lambda D, col: D.pivot(index="penetration", columns="grp", values=col)         # noqa: E731
    Dw = rows_panel(ax["d"], {"drive (a)": (by(N, "cond_p_a"), fs.MUTED, 0.13, "white", "o"),
                              "pool (w)": (by(N, "cond_p_w"), fs.INK, -0.13, fs.INK, "o")}, GROUPS, (0, 2.0), ref=1.0)
    ax["d"].set_xlabel("the other eye's weight\n(dominant eye = 1)")
    ax["d"].set_title("the other eye's weights", fontsize=6.5, loc="left")
    ax["d"].text(0.97, 0.03, "filled: pool\nopen: drive", transform=ax["d"].transAxes, ha="right", va="bottom",
                 fontsize=5.5, color=fs.INK2, linespacing=1.0)

    G[["data", "conductance", "two-stage"]] *= 100
    De = rows_panel(ax["e"], {"data": (by(G, "data"), fs.INK, -0.2, fs.INK, "o"),
                              "conductance": (by(G, "conductance"), fs.VIOLET, 0.0, "white", "o"),
                              "two-stage": (by(G, "two-stage"), TWO_STAGE, 0.2, "white", "D")}, GROUPS, (-16, 8), ref=0)
    ax["e"].tick_params(labelleft=False)
    ax["e"].set_xlabel("both eyes − contrast summation,\n45%, 30–100 ms (% of peak)")
    ax["e"].set_title("what a single-stage pool misses", fontsize=6.5, loc="left")
    d0 = A.loc["(reported) data deep - input"]
    d1 = A.loc["AA2 conductance: model deep - input > 0"]
    d2 = A.loc["AA2 two-stage: model deep - input > 0"]
    ax["e"].text(0.98, 0.03, (f"deep − input\ndata {100 * d0['median']:+.1f}, {fs.fmt_p(d0.p)}\n"
                              f"conductance {100 * d1['median']:+.1f}, {fs.fmt_p(d1.p)}\n"
                              f"two-stage {100 * d2['median']:+.1f}, {fs.fmt_p(d2.p)}").replace("-", "−"),
                 transform=ax["e"].transAxes, ha="right", va="bottom", fontsize=5.2, color=fs.INK2, linespacing=1.1)
    fs.line_key(ax["e"], [("data", fs.INK, "none", "o"), ("conductance fit", fs.VIOLET, "none", "o", "white"),
                          ("two-stage fit", TWO_STAGE, "none", "D", "white")], x=0.03, y=1.0, dy=0.075, length=0.06)

    for k, letter, x in (("a1", "a", -0.02), ("b", "b", -0.2), ("c", "c", -0.42), ("d", "d", -0.5), ("e", "e", -0.06)):
        fs.panel_label(ax[k], letter, x=x)
    print(fs.save(fig, "fig4_delayed_pool"))
    pd.set_option("display.width", 200)
    print("(a)", sk)
    print("(b) medians over penetrations, % of dominant-eye peak\n", Bt.round(2).to_string(index=False))
    print("(c)\n", Co.pivot(index="fit", columns="c_oth", values="mean").round(2).to_string())
    print(Co.pivot(index="fit", columns="c_oth", values="n").to_string())
    print("(d)\n", Dw.round(3).to_string(index=False))
    print("(e)\n", De.round(2).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
