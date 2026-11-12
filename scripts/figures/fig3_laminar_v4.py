"""Figure 3 on laminar alignment v4 (docs/09 AD; Blake 2026-11-12: tiers 1 + 2 as the main laminar set, coarse
compartments in the main figure): normalization is strongest in and above layer 4 and weaker in layers 5-6.

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_laminar_v4.py [tier1|all]
Penetrations: tiers 1 + 2 of metadata/laminar_alignment_v4.csv (12; `tier1` or `all` for the sensitivity sets).
Compartments (ad_compartments_v4.label): height h above the 4C/5 anchor, scaled by the penetration's squish factor s:
supragranular (layers 2/3) h > 0.59 s, granular (layer 4) 0 < h <= 0.59 s, infragranular (layers 5-6) -0.49 s <= h < 0.
Additivity index AI = both eyes / (dominant eye + other eye); dashed reference, the index predicted by contrast
summation (the dominant eye's response at c (1 + a), each channel's eye weight a). Channels: check X's per-channel
arrays (the 50 um probes recomputed; ad_compartments_v4.load_spacetime).
(a) AI over time at 45% in each eye per compartment: population = mean over penetrations (>= 2 channels in the
    compartment) of channel means, ratio of the means, 5 ms running mean; orange where the index exceeds contrast
    summation, grey where it falls below.
(b) AI 30-100 ms after onset per penetration and compartment (docs/qc/ad_groups.csv), measured and summation;
    medians and interquartile ranges; AD1 (infragranular - granular, one-sided Wilcoxon).
(c) MUA for both eyes at 45% each against scaled height in 0.1 mm bins (>= MIN_PEN penetrations), each bin scaled to
    its own peak, with the median onset (median over penetrations of each penetration's median channel onset) for both
    eyes, the dominant eye alone and the other eye alone at 45%; stimulus clock.
Writes <out_root>/figures/fig3_laminar_v4[_tier1|_all].pdf and .png; prints the values shown.
"""

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ad_compartments_v4 import AD_PFX, ALIGN, alignment, load_channels, load_spacetime, relabel  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import metadata_dir  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 104
COARSE = (("supragranular", "supragranular (layers 2/3)"), ("granular", "granular (layer 4)"),
          ("infragranular", "infragranular (layers 5–6)"))
SETS = {"main": (1, 2), "tier1": (1,), "all": (0, 1, 2)}
SET_NAME = {"main": "tiers 1+2", "tier1": "tier 1", "all": "all anchored"}
CS4 = (0.11, 0.225, 0.45, 0.9)
MIN_CH, MIN_PEN, BIN_MM = 2, 8, 0.1
ROWS = {"both": 12, "dominant": list(CONTRASTS).index(0.45), "other": 8}      # space-time rows at 45%
POS = {"a": (22, 10, 58, 84),
       "b": (104, 60, 40, 34), "c": (104, 10, 50, 36), "ccb": (157, 10, 2, 36)}
FINE_BOUNDS = ((0.59, "-"), (0.28, ":"), (0.0, "-"), (-0.22, ":"))      # layer boundaries (scaled mm), check AD
FINE_LABELS = ((0.9, "L2/3"), (0.435, "L4A/B"), (0.14, "L4C"), (-0.11, "L5"), (-0.355, "L6"))
LAYERS4 = (("L2/3", ("L2/3",)), ("L4", ("L4A/B", "L4C")), ("L5", ("L5",)), ("L6", ("L6",)))


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def comp_ai(X, meta, comp, c=0.45):
    i, j = list(CONTRASTS).index(c), CS4.index(c)
    sel = (meta.coarse == comp).to_numpy()
    pen = meta.penetration.to_numpy()
    pens = [p for p in np.unique(pen[sel]) if np.sum(sel & (pen == p)) >= MIN_CH]
    P = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    M = {}
    for name, x in zip("DOBP", (X[:, i], X[:, 6 + j], X[:, 10 + j], P)):
        M[name] = fs.smooth(np.nanmean([np.nanmean(x[sel & (pen == p)], axis=0) for p in pens], axis=0)[None])[0]
    den = M["D"] + M["O"]
    ok = den > 0.12
    return np.where(ok, M["B"] / den, np.nan), np.where(ok, M["P"] / den, np.nan), len(pens)


def panel_a(axs, X, meta, rel):
    t_sel = (rel >= 0) & (rel <= 150)
    t = rel[t_sel]
    d = (rel >= 30) & (rel < 100)
    rows = []
    for ax, (comp, title) in zip(axs, COARSE):
        AI, AS, n = comp_ai(X, meta, comp)
        a, p = AI[t_sel], AS[t_sel]
        ax.fill_between(t, a, p, where=a >= p, color=fs.ORANGE, alpha=0.55, lw=0, interpolate=True)
        ax.fill_between(t, a, p, where=a < p, color=fs.MUTED, alpha=0.55, lw=0, interpolate=True)
        ax.plot(t, p, color=fs.MUTED, lw=0.9, ls=fs.SUM_LS)
        ax.plot(t, a, color=fs.INK, lw=1.1)
        ax.text(0.99, 0.95, f"30–100 ms: {np.nanmean(AI[d]):.2f}; {fs.SUM}, {np.nanmean(AS[d]):.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK2)
        ax.set_title(f"{title}; {n} penetrations", fontsize=6.5, loc="left")
        ax.axvspan(30, 100, color=fs.RULE, alpha=0.3, lw=0, zorder=0)
        ax.set_xlim(0, 150)
        ax.set_ylim(0.35, 1.05)
        ax.set_yticks([0.4, 0.6, 0.8, 1.0])
        ax.set_ylabel("additivity index")
        if ax is not axs[-1]:
            ax.tick_params(labelbottom=False)
        rows.append(dict(compartment=comp, n_pen=n, AI_dip=np.nanmean(AI[d]), AI_sum_dip=np.nanmean(AS[d]),
                         AI_early=np.nanmean(AI[(rel >= 0) & (rel < 20)]),
                         AI_sum_early=np.nanmean(AS[(rel >= 0) & (rel < 20)])))
    axs[-1].set_xlabel("time from response onset (ms)")
    fs.line_key(axs[0], [("measured, both eyes 45%", fs.INK, "-")], x=0.0, y=1.36, dy=0.13, length=0.08)
    fs.line_key(axs[0], [(fs.SUM, fs.MUTED, fs.SUM_LS)], x=0.52, y=1.36, dy=0.13, length=0.08)
    return pd.DataFrame(rows)


PROFILE_BIN = 0.2                                   # mm (scaled) for the depth profile of panel a
CI = (16, 84)                                       # bootstrap 68% interval (about +-1 s.e.), stated in the caption


def gap_profile(C, window):
    """AI - AI_sum at 45% in a window, per penetration and 0.1 mm bin of scaled height (ratio of channel means, as
    x_depth_v3.ai_groups); median over penetrations and bootstrap 95% CI, bins with >= MIN_PEN penetrations."""
    X = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == window) & (C.coarse != "")].copy()
    X["bin"] = np.round(np.floor(X.h_scaled / PROFILE_BIN + 1e-9) * PROFILE_BIN + PROFILE_BIN / 2, 2)
    g = X.assign(DO=X.D + X.O).groupby(["penetration", "bin"])[["B", "P_sum", "DO"]].mean()
    g = g[g.DO > 0.12]
    per = ((g.B - g.P_sum) / g.DO).unstack()
    per = per.loc[:, per.notna().sum() >= MIN_PEN]
    rng = np.random.default_rng(3)
    boot = np.array([np.nanmedian(per.to_numpy()[rng.integers(0, len(per), len(per))], axis=0) for _ in range(2000)])
    lo, hi = np.nanpercentile(boot, CI, axis=0)
    return pd.DataFrame(dict(height=per.columns.to_numpy(), median=per.median().to_numpy(), lo=lo, hi=hi,
                             n=per.notna().sum().to_numpy()))


def panel_a(ax, C):
    """Depth profile of the gap to contrast summation: 0-20 ms (orange) and 30-100 ms (black)."""
    out = []
    for window, col, lab in (("early", fs.ORANGE, "0–20 ms after onset"), ("dip", fs.INK, "30–100 ms")):
        P = gap_profile(C, window)
        ax.fill_betweenx(P.height, P.lo, P.hi, color=col, alpha=0.18, lw=0)
        ax.plot(P["median"], P.height, color=col, lw=1.3, marker="o", ms=2.2)
        out.append(P.assign(window=window))
    ax.axvline(0, color=fs.MUTED, lw=0.8, ls=fs.SUM_LS)
    ax.text(0.005, 1.32, fs.SUM, fontsize=5.5, color=fs.INK2, ha="left", va="center")
    z = pd.concat(out).height
    for yb, ls in FINE_BOUNDS:
        ax.axhline(yb, color=fs.INK2 if ls == "-" else fs.RULE, lw=0.6, ls=ls, zorder=0)
    for y, lab in FINE_LABELS:
        if z.min() - 0.05 <= y <= z.max() + 0.05:
            ax.text(0.97, y, lab, transform=ax.get_yaxis_transform(), fontsize=5.5, ha="right", va="center",
                    color=fs.INK2)
    ax.set_ylim(z.min() - 0.08, z.max() + 0.12)
    ax.set_xlim(-0.25, 0.18)
    ax.set_xlabel("additivity index − contrast summation\n(both eyes at 45%)")
    ax.set_ylabel("height above the layer 4C/5 boundary (mm, scaled)")
    ax.set_title("where the eyes add and where they normalize", fontsize=6.5, loc="left")
    fs.line_key(ax, [("0–20 ms after onset", fs.ORANGE, "-"), ("30–100 ms", fs.INK, "-")], x=0.03, y=0.10, dy=0.04,
                length=0.08)
    return pd.concat(out)


def panel_b(ax, C, T, set_name):
    """The index 30-100 ms in layers 2/3, 4, 5 and 6 (fine labels, layer 4 pooled)."""
    from x_depth_v3 import ai_groups
    X = C.copy()
    X["layer"] = X.fine.map({f: lab for lab, fs_ in LAYERS4 for f in fs_}).fillna("")
    G = ai_groups(X, grp="layer").reset_index()
    rows = []
    for k, (lab, _) in enumerate(LAYERS4):
        x = G[G.layer == lab]
        for col, colr, dy, mfc in (("AI", fs.INK, -0.12, fs.INK), ("AI_sum", fs.MUTED, 0.12, "white")):
            v = x[col].dropna()
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            rows.append(dict(layer=lab, series=col, n=len(v), median=md))
    ax.set_yticks(range(len(LAYERS4)), ["layers 2/3", "layer 4", "layer 5", "layer 6"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(LAYERS4) - 0.4, -0.6)
    ax.set_xlim(0.3, 0.9)
    ax.set_xlabel("additivity index, 30–100 ms")
    r = T.set_index("test").loc[f"AD1 (AI - AI_sum) infragranular - granular > 0 [{set_name}]"]
    ax.text(1.03, 0.98, f"measured\n{fs.SUM} (open)\n\nlayers 5–6 − layer 4:\n{int(r.n_in_direction)}/{int(r.n)} "
            f"penetrations,\n{fs.fmt_p(r.p)}", transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=fs.INK2, linespacing=1.1)
    return pd.DataFrame(rows)


def panel_b_coarse(ax, G, T, set_name):
    rows = []
    for k, (comp, _) in enumerate(COARSE):
        x = G[G.compartment == comp]
        for col, colr, dy, mfc in (("AI", fs.INK, -0.12, fs.INK), ("AI_sum", fs.MUTED, 0.12, "white")):
            v = x[col].dropna()
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            rows.append(dict(compartment=comp, series=col, n=len(v), median=md))
    ax.set_yticks(range(3), ["layers 2/3", "layer 4", "layers 5–6"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(0.3, 0.9)
    ax.set_xlabel("additivity index, 30–100 ms")
    r = T.set_index("test").loc[f"AD1 (AI - AI_sum) infragranular - granular > 0 [{set_name}]"]
    ax.text(1.03, 0.98, f"measured\n{fs.SUM} (open)\n\nlayers 5–6 − layer 4:\n{int(r.n_in_direction)}/{int(r.n)} "
            f"penetrations,\n{fs.fmt_p(r.p)}", transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=fs.INK2, linespacing=1.1)
    return pd.DataFrame(rows)


def panel_c(ax, cax, S, t, meta, L):
    hb = np.round(np.floor(meta.h_scaled.to_numpy() / BIN_MM + 1e-9) * BIN_MM + BIN_MM / 2, 2)
    pen = meta.penetration.to_numpy()
    use = (meta.coarse != "").to_numpy()
    z, curves = [], {k: [] for k in ROWS}
    for b in np.unique(hb[use]):
        sel = use & (hb == b)
        pens = np.unique(pen[sel])
        if len(pens) < MIN_PEN:
            continue
        z.append(b)
        for name, r in ROWS.items():
            curves[name].append(np.nanmean([np.nanmean(S[sel & (pen == p), r], axis=0) for p in pens], axis=0))
    z = np.array(z)
    C = {k: np.array(v) for k, v in curves.items()}
    w = (t >= 0) & (t < 100)
    img = C["both"] / np.nanmax(C["both"][:, w], axis=1, keepdims=True)
    im = ax.imshow(img, aspect="auto", origin="lower", cmap="Greys", vmin=0, vmax=1, interpolation="bilinear",
                   extent=(t[0] - 0.5, t[-1] + 0.5, z[0] - BIN_MM / 2, z[-1] + BIN_MM / 2), rasterized=True)
    A = L[(L.split == "all") & (L.coarse != "")].copy()
    A["bin"] = np.round(np.floor(A.h_scaled / BIN_MM + 1e-9) * BIN_MM + BIN_MM / 2, 2)
    per = A.groupby(["penetration", "bin"])[["dom_45", "oth_45", "bin_45"]].median().reset_index()
    med = per.groupby("bin")[["dom_45", "oth_45", "bin_45"]].median().reindex(z)
    on = {"both": med.bin_45.to_numpy(), "dominant": med.dom_45.to_numpy(), "other": med.oth_45.to_numpy()}
    for name, col, lab in (("other", fs.COND["other"], "other eye alone"),
                           ("dominant", fs.COND["dominant"], "dominant eye alone"), ("both", fs.INK, "both eyes")):
        ax.plot(on[name], z, color="white", lw=2.4, zorder=3)
        ax.plot(on[name], z, color=col, lw=1.1, marker="o", ms=1.8, zorder=4, label=lab)
    ax.set_ylim(z[0] - BIN_MM / 2, z[-1] + BIN_MM / 2)
    ax.set_xlim(0, 50)
    for yb, ls in FINE_BOUNDS:
        ax.axhline(yb, color=fs.INK2 if ls == "-" else fs.MUTED, lw=0.5, ls=(0, (2, 2)) if ls == "-" else ":")
    for y, lab in FINE_LABELS:
        if z[0] <= y <= z[-1] + 0.05:
            ax.text(1.0, y, lab, fontsize=5.2, va="center", ha="left", color=fs.INK2)
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("height above the layer 4C/5\nboundary (mm, scaled)")
    ax.set_title("onset at each depth, 45% in each eye", fontsize=6.5)
    fs.line_key(ax, [("both eyes", fs.INK, "-"), ("dominant eye alone", fs.COND["dominant"], "-"),
                     ("other eye alone", fs.COND["other"], "-")], x=0.22, y=0.97, dy=0.075, length=0.07)
    cb = ax.figure.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("MUA, both eyes\n(each depth's peak = 1)", fontsize=6)
    return pd.DataFrame(dict(height_mm=z, **{f"onset_{k}": v for k, v in on.items()}))


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    tiers = SETS[which]
    fs.apply()
    V = alignment()
    V = V[V.s.notna()]
    on, stim, M, rel, t = load_spacetime(V)
    keep_idx = np.flatnonzero(M.penetration.isin(V[V.tier.isin(tiers)].penetration).to_numpy())
    on, stim, M = on[keep_idx], stim[keep_idx], M.iloc[keep_idx].reset_index(drop=True)
    meta = relabel(M.assign(_i=np.arange(len(M))), V)
    idx = meta._i.to_numpy()
    on, stim, meta = on[idx], stim[idx], meta.reset_index(drop=True)
    _, L = load_channels(V)
    L = relabel(L[L.penetration.isin(meta.penetration.unique())], V)
    G = pd.read_csv(QC / f"{AD_PFX[ALIGN]}groups.csv")
    G = G[(G.scheme == "coarse") & G.tier.isin(tiers)]
    T = pd.read_csv(QC / f"{AD_PFX[ALIGN]}tests.csv")
    C0, _ = load_channels(V)
    C = relabel(C0[C0.penetration.isin(V[V.tier.isin(tiers)].penetration)], V)
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    A = panel_a(ax["a"], C)
    B = panel_b(ax["b"], C, T, SET_NAME[which])
    O = panel_c(ax["c"], ax["ccb"], stim, t, meta, L)
    for k, lab, x in (("a", "a", -0.2), ("b", "b", -0.32), ("c", "c", -0.32)):
        fs.panel_label(ax[k], lab, x=x)
    print(fs.save(fig, "fig3_laminar_v4" + ("" if which == "main" else f"_{which}")))
    pd.set_option("display.width", 200)
    print(f"{meta.penetration.nunique()} penetrations; channels per compartment:",
          meta[meta.coarse != ""].coarse.value_counts().to_dict())
    print(A.round(3).to_string(index=False))
    print(B.round(3).to_string(index=False))
    print(O.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
