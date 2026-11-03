"""Figure 4 (space-time figure set, docs/09 checks K and N): delayed normalization across time and depth.

Usage: .venv\\Scripts\\python scripts\\figures\\fig4_model_spacetime.py
(a) what the model does (illustration, as the former Fig 4a): drive, pool and response for both eyes at 45% and the
    dominant eye alone, with the median fitted weights and pool delay (conductance fits per depth group);
(b) what the other eye adds at 45%, depth x time on the response clock: data, no-delay fits and conductance fits of
    each channel (v4_model_depth.py; all presentations); population per 0.2 mm depth bin (>= 10 penetrations);
(c) data minus each model;
(d) held-out error, % less than without a delay, per penetration and depth group (above / at, within 0.2 mm / below
    the earliest sink; group fits, v4_model_groups.py, both cross-validation directions summed); medians and
    interquartile ranges; one-sided Wilcoxon (check N1);
(e) the same for weighting the other eye in the pool separately from the drive (fixed delay vs w = a, N2) and for
    the conductance form (conductance vs fixed delay, N3);
(f) the other eye's extra weight in the pool (pool weight - drive weight, conductance fits, K2/N4);
(g) the pool delay (conductance fits);
(h) both eyes - "if the contrasts add" at 45%, 30-100 ms after the group's onset, in the data and in the conductance
    fits of the same groups (descriptive).
Inputs: <out_root>/binocular/v4_K_channels.parquet, v4_K_pred.npz, v4_K_cond_channels.parquet, v4_K_cond_pred.npz,
v4_N_groups.parquet, v4_N_pred.npz.
Writes <out_root>/figures/fig4_model_spacetime.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from fig4_model import timing_sketch  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, eye_weight, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
H_MM = 160
REL = np.arange(-10, 151)
BIN = 2
EDGES = np.arange(-9, 14, BIN)
MIN_PEN = 10
GROUPS = ("above", "at", "below")
GROUP_LABEL = {"above": "above", "at": "at sink", "below": "below"}
RES_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("residual", [fs.INK2, "#ffffff", fs.VIOLET])
POS = {"a1": (12, 141, 46, 9), "a2": (12, 126, 46, 9), "a3": (12, 106, 46, 15),
       "b1": (84, 106, 24, 40), "b2": (111, 106, 24, 40), "b3": (138, 106, 24, 40), "bcb": (165, 106, 2, 40),
       "c1": (16, 58, 26, 30), "c2": (45, 58, 26, 30), "ccb": (74, 58, 2, 30),
       "d": (104, 60, 26, 28), "e": (144, 60, 26, 28),
       "f": (22, 12, 30, 27), "g": (76, 12, 30, 27), "h": (132, 12, 36, 27)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def added_maps(K, X, T, cd, cn, c=0.45):
    """Per depth bin, the penetration-averaged response the other eye adds, on each channel's response clock."""
    iD = np.flatnonzero(np.isclose(cd, c) & (cn == 0))[0]
    iB = np.flatnonzero(np.isclose(cd, c) & np.isclose(cn, c))[0]
    added = {k: np.stack([np.interp(on + REL, T, v[q, iB] - v[q, iD]) for q, on in enumerate(K.onset)])
             for k, v in X.items()}
    b = np.digitize(K.rel_sink.to_numpy(), EDGES) - 1
    centres = (EDGES[:-1] + BIN / 2) * 0.1
    pen = K.penetration.to_numpy()
    out = {k: np.full((len(centres), REL.size), np.nan) for k in added}
    for j in range(len(centres)):
        sel = b == j
        pens = np.unique(pen[sel])
        if len(pens) < MIN_PEN:
            continue
        for k, a in added.items():
            out[k][j] = np.nanmean([a[sel & (pen == p)].mean(axis=0) for p in pens], axis=0)
    ok = np.isfinite(out["data"]).any(axis=1)
    return centres[ok], {k: fs.smooth(v[ok]) for k, v in out.items()}


def show_map(ax, M, z, cmap, lim, title):
    im = ax.imshow(M, aspect="auto", origin="lower", cmap=cmap, vmin=-lim, vmax=lim, interpolation="bilinear",
                   extent=(REL[0] - 0.5, REL[-1] + 0.5, z[0] - 0.1, z[-1] + 0.1), rasterized=True)
    ax.set_ylim(z[-1] + 0.1, z[0] - 0.1)
    ax.set_xticks([0, 50, 100, 150])
    ax.axvline(0, color=fs.INK, lw=0.6)
    ax.axhline(0, color=fs.INK, lw=0.5, ls=(0, (2, 2)))
    ax.set_title(title, fontsize=6.5)
    return im


def colorbar(fig, cax, im, label, ticks):
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label(label, fontsize=6)


def group_dots(ax, series: dict, xlabel, ref=0.0, xlim=None, clip=None):
    """series: {label: (pd.Series indexed by (penetration, grp) or DataFrame with grp column, colour, dy)}."""
    rng = np.random.default_rng(3)
    rows = []
    for lab, (S, col, dy) in series.items():
        for k, g in enumerate(GROUPS):
            v = S[S.grp == g].value.dropna()
            y = k + dy
            shown = v.clip(*clip) if clip else v
            ax.scatter(shown, y + rng.uniform(-0.07, 0.07, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [y, y], color=col, lw=0.9, zorder=2)
            ax.plot([md], [y], marker="o", ms=3.2, color=col, zorder=3)
            rows.append(dict(series=lab, grp=g, n=len(v), median=md, n_pos=int((v > 0).sum()),
                             p_greater=stats.wilcoxon(v, alternative="greater").pvalue if len(v) > 5 else np.nan))
    if ref is not None:
        ax.axvline(ref, color=fs.MUTED, lw=0.6)
    ax.set_yticks(range(3), [GROUP_LABEL[g] for g in GROUPS])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    return pd.DataFrame(rows)


def long(N, col):
    return pd.DataFrame(dict(penetration=N.penetration, grp=N.grp, value=N[col]))


def group_deficit(D, Pc, cd, cn, T, c=0.45, w=(30, 100)):
    """Both eyes - "if the contrasts add" in w ms after the group's onset, for data and model [group, cond, t]."""
    out = []
    for q in range(D.shape[0]):
        r = {}
        for name, Y in (("data", D[q]), ("conductance", Pc[q])):
            full = np.zeros((6, 6, T.size))
            for k, (a_, b_) in enumerate(zip(cd, cn)):
                full[list(CONTRASTS).index(a_), list(CONTRASTS).index(b_)] = Y[k]
            dom, oth = full[:, 0], full[0, :]
            j = np.flatnonzero((T >= 10) & (dom[-1] >= 0.5 * dom[-1].max()))
            on = T[j[0]]
            wt = (T >= 40)
            a = eye_weight(dom[:, wt].mean(axis=1), oth[:, wt].mean(axis=1))
            i = list(CONTRASTS).index(c)
            pred = response_at(dom, c * (1 + a))
            sel = (T >= on + w[0]) & (T < on + w[1])
            r[name] = float(np.mean(full[i, i, sel] - pred[sel]))
        out.append(r)
    return pd.DataFrame(out)


def main():
    fs.apply()
    b = out_root() / "binocular"
    K = pd.read_parquet(b / "v4_K_channels.parquet")
    Zk = np.load(b / "v4_K_pred.npz")
    Kc = pd.read_parquet(b / "v4_K_cond_channels.parquet")
    Zc = np.load(b / "v4_K_cond_pred.npz")
    key = lambda d: list(zip(d.penetration, d.site_id))           # noqa: E731
    pos = {k: i for i, k in enumerate(key(Kc))}
    order = [pos[k] for k in key(K)]
    X = {"data": Zk["data"], "static": Zk["static"], "cond": Zc["cond"][order]}
    z, M = added_maps(K, X, Zk["t"], Zk["cd"], Zk["cn"])
    N = pd.read_parquet(b / "v4_N_groups.parquet")
    Zn = np.load(b / "v4_N_pred.npz", allow_pickle=True)

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    sk = timing_sketch([ax["a1"], ax["a2"], ax["a3"]], float(N.cond_p_a.median()), float(N.cond_p_w.median()),
                       float(N.cond_p_tau2.median()))
    titles = {"data": "data", "static": "no delay (static)", "cond": "delayed pool (conductance)"}
    for k, name in zip(("b1", "b2", "b3"), ("data", "static", "cond")):
        im = show_map(ax[k], M[name], z, fs.ADDED_CMAP, 0.25, titles[name])
        if k != "b1":
            ax[k].tick_params(labelleft=False)
    ax["b1"].set_ylabel("depth relative to\nearliest sink (mm)")
    ax["b2"].set_xlabel("time from response onset (ms)")
    colorbar(fig, ax["bcb"], im, "added by the other eye\n(both eyes 45%)", [-0.2, -0.1, 0, 0.1, 0.2])
    for k, name in zip(("c1", "c2"), ("static", "cond")):
        im2 = show_map(ax[k], M["data"] - M[name], z, RES_CMAP, 0.12, f"data − {titles[name].split(' (')[0]}")
    fig.text(44.5 / 174, 48 / H_MM, "time from response onset (ms)", ha="center", va="center", fontsize=7)
    ax["c2"].tick_params(labelleft=False)
    ax["c1"].set_ylabel("depth relative to\nearliest sink (mm)")
    colorbar(fig, ax["ccb"], im2, "data − model", [-0.1, 0, 0.1])

    pct = lambda worse, better: 100 * (1 - N[f"cv_{better}"] / N[f"cv_{worse}"])   # noqa: E731
    N = N.assign(pct_delay=pct("static", "full"), pct_weight=pct("w=a", "full"), pct_cond=pct("full", "cond"),
                 w_minus_a=N.cond_p_w - N.cond_p_a)
    Rd = group_dots(ax["d"], {"delay": (long(N, "pct_delay"), fs.INK, 0.0)},
                    "held-out error, % less\nwith a delayed pool", xlim=(-10, 100))
    Re = group_dots(ax["e"], {"pool weighting": (long(N, "pct_weight"), fs.MODEL["fixed-delay"], -0.15),
                              "conductance": (long(N, "pct_cond"), fs.MODEL["conductance"], 0.15)},
                    "held-out error, % less", xlim=(-2, 12), clip=(-2, 12))
    ax["e"].tick_params(labelleft=False)
    ax["e"].text(1.0, 1.13, "separate pool weight (vs w = a)", transform=ax["e"].transAxes,
                 ha="right", va="bottom", fontsize=5.5, color=fs.MODEL["fixed-delay"])
    ax["e"].text(1.0, 1.03, "conductance (vs fixed delay)", transform=ax["e"].transAxes, ha="right",
                 va="bottom", fontsize=5.5, color=fs.MODEL["conductance"])
    Rf = group_dots(ax["f"], {"w - a": (long(N, "w_minus_a"), fs.INK, 0.0)},
                    "other eye: pool weight −\ndrive weight", xlim=(-0.4, 1.0), clip=(-0.4, 1.0))
    Rg = group_dots(ax["g"], {"tau2": (long(N, "cond_p_tau2"), fs.INK, 0.0)}, "pool delay (ms)", ref=None,
                    xlim=(0, 120), clip=(0, 120))
    ax["g"].tick_params(labelleft=False)
    Gd = group_deficit(Zn["data"], Zn["cond"], Zn["cd"], Zn["cn"], Zn["t"])
    Gd["grp"] = Zn["grp"]
    Gd["penetration"] = Zn["penetration"]
    rows = []
    for k, g in enumerate(GROUPS):
        for name, col, dy in (("data", fs.INK, -0.15), ("conductance", fs.MODEL["conductance"], 0.15)):
            v = Gd[Gd.grp == g][name]
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax["h"].plot([q1, q3], [k + dy] * 2, color=col, lw=0.9)
            ax["h"].plot([md], [k + dy], marker="o", ms=3.2, color=col)
            rows.append(dict(source=name, grp=g, n=len(v), median=md))
    ax["h"].axvline(0, color=fs.MUTED, lw=0.8, ls=fs.SUM_LS)
    ax["h"].set_yticks(range(3), [GROUP_LABEL[g] for g in GROUPS])
    ax["h"].tick_params(axis="y", length=0)
    ax["h"].set_ylim(2.6, -0.6)
    ax["h"].set_xlabel("both eyes − if the contrasts add,\n45%, 30–100 ms after onset")
    ax["h"].text(1.0, 1.13, "data", transform=ax["h"].transAxes, fontsize=5.5, color=fs.INK, ha="right",
                 va="bottom")
    ax["h"].text(1.0, 1.03, "conductance fits", transform=ax["h"].transAxes, fontsize=5.5,
                 color=fs.MODEL["conductance"], ha="right", va="bottom")
    for k, letter in (("a1", "a"), ("b1", "b"), ("c1", "c"), ("d", "d"), ("e", "e"), ("f", "f"), ("g", "g"),
                      ("h", "h")):
        fs.panel_label(ax[k], letter, x={"d": -0.18, "f": -0.18, "h": -0.18, "e": -0.3}.get(k, -0.02))
    print(fs.save(fig, "fig4_model_spacetime"))
    print(f"(b, c) {K.penetration.nunique()} penetrations, {len(K)} channels; depth bins {np.round(z, 1).tolist()}")
    for name in ("static", "cond"):
        for w0, w1 in ((0, 20), (30, 100), (130, 150)):
            sel = (REL >= w0) & (REL < w1)
            print(f"   data - {name}, {w0}-{w1} ms: mean over depth {np.nanmean((M['data'] - M[name])[:, sel]):+.3f}")
    print(f"(a) sketch with a {N.cond_p_a.median():.2f}, w {N.cond_p_w.median():.2f}, tau2 {N.cond_p_tau2.median():.0f} ms;"
          f" {sk if sk is not None else ''}")
    for lab, R in (("(d) % less error with a delay", Rd), ("(e) pool weighting / conductance", Re),
                   ("(f) w - a", Rf), ("(g) pool delay", Rg)):
        print(lab)
        print(R.round(4).to_string(index=False))
    print("(h) both eyes - if the contrasts add, 45%, 30-100 ms, medians over groups:")
    print(pd.DataFrame(rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
