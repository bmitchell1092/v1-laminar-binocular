"""Figure 3 (docs/09 checks V, X, G): the laminar profile. Summation gives way to normalization within ~30 ms, more
in and above the input layer than below it; responses start first just below the input/deep boundary, and the second
eye advances the onset at every depth. (2026-11-11: panels c-d, the index against the other eye's contrast, are
superseded by Fig 2, fig2_dichoptic_time.py; the laminar onset, formerly Fig 2b of the simplified set, is panel c.)

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_normalization_ai.py
Additivity index AI = both eyes / (dominant eye + other eye); the dashed reference is the index that summation of the
two eyes' contrasts predicts: the dominant eye's response at c_D + a c_O over the same denominator (each channel's or
penetration's eye weight a). Time from each channel's response onset.
(a) AI over time at 45% in each eye in three bands of five contacts (0.5 mm) from the input/deep boundary of laminar
    alignment v3: 0.5-1.0 mm above, the 0.5 mm above (which holds the earliest sink in 17/22 penetrations, check Z1)
    and the 0.5 mm below (Blake, comment on V4, 2026-11-11: three traces matched in contacts instead of the
    ridgeline, which moved to figS_ridgeline.py); population: mean over penetrations (>= 2 channels in the band) of
    the channel means, ratio of the means, 5 ms running mean; grey where the index falls below contrast summation;
(b) AI 30-100 ms after onset per penetration and depth group (superficial / input, 0.5 mm above the boundary /
    deep, to 0.7 mm below), data and summation reference; medians and interquartile ranges; check X2 (deep - input,
    one-sided Wilcoxon);
(c) the first 50 ms at each depth: MUA for both eyes at 45% each, each depth bin scaled to its own peak (0-100 ms),
    with the median onset at each depth for both eyes, the dominant eye alone and the other eye alone at 45%
    (fig2_latency_drive.laminar_onset; checks G, X3, X4), stimulus clock; tick, the median depth of the earliest
    sink (check Z1, docs/qc/z_onset_sink.csv).
`sink` as argument: the earlier version, depth relative to the earliest sink (groups above / at / below; check V2).
Inputs: x_depth_v3.py (spacetime_v3.npz, docs/qc/x_depth_v3.csv, x_groups.csv) or spacetime_data.py (spacetime.npz;
docs/qc/v4_V.csv, v4_V_groups.csv); for (c) <out_root>/binocular/x_latency.parquet (v4_G_channels.parquet for
`sink`).
Writes <out_root>/figures/fig3_normalization_ai.pdf and .png (fig3_normalization_ai_sink for `sink`); prints the
statistics shown.
"""

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig2_latency_drive import DEPTH as ONSET_DEPTH, laminar_onset  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 104
BIN = 2                                          # contacts (0.2 mm)
MIN_PEN = 10
DEPTH = {
    "v3": dict(data="spacetime_v3", rel="rel_anchor_mm", per_contact=10, edges=np.arange(-16, 9, BIN),
               label="depth from input/deep boundary (mm; negative: above)",
               groups=("superficial", "input", "deep"),
               names=("superficial", "input", "deep"), files=("x_groups.csv", "x_depth_v3.csv"),
               test="X2 (AI - AI_sum) deep - input > 0, 45%, 30-100 ms", contrast="deep − input"),
    "sink": dict(data="spacetime", rel="rel_sink", per_contact=1, edges=np.arange(-9, 14, BIN),
                 label="depth relative to earliest sink (mm)", groups=("above", "at", "below"),
                 names=("above", "at sink", "below"), files=("v4_V_groups.csv", "v4_V.csv"),
                 test="V2 (AI - AI_sum) below - at the sink > 0, 45%, 30-100 ms", contrast="below − at sink"),
}
CS4 = (0.11, 0.225, 0.45, 0.9)
# (a) three bands of five contacts each (Blake, comment on V4: one trace above, at and below the sink, matched in
# contacts); the input band holds the earliest sink in 17 of 22 penetrations (check Z1). mm from the boundary, (lo, hi]
BANDS = {"superficial": (-1.0, -0.5, "0.5–1.0 mm above the boundary"),
         "input": (-0.5, 0.0, "0.5 mm above the boundary (earliest sink)"),
         "deep": (0.0, 0.5, "0.5 mm below the boundary")}
MIN_CH = 2                                       # channels per penetration in a band
POS = {"a0": (20, 68, 62, 24), "a1": (20, 39, 62, 24), "a2": (20, 10, 62, 24),
       "b": (104, 58, 40, 36), "c": (104, 10, 50, 34), "ccb": (157, 10, 2, 34)}


POS_A_SINK = (20, 10, 62, 86)


def axes_box(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return x / W, y / H, w / W, h / H


def axes_mm(fig, x, y, w, h):
    return fig.add_axes(axes_box(fig, x, y, w, h))


def depth_ai(X, meta, cfg, c=0.45):
    i, j = list(CONTRASTS).index(c), CS4.index(c)
    D, O, B = X[:, i], X[:, 6 + j], X[:, 10 + j]
    P = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    edges = cfg["edges"]
    b = np.digitize(meta[cfg["rel"]].to_numpy() * cfg["per_contact"], edges) - 1
    centres = (edges[:-1] + BIN / 2) * 0.1
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


def band_ai(X, meta, lo, hi, c=0.45):
    """AI and its summation reference over time for one depth band: ratio of the population means (over penetrations
    with >= MIN_CH channels in the band) of channel means, 5 ms running mean."""
    i, j = list(CONTRASTS).index(c), CS4.index(c)
    rel = np.round(meta.rel_anchor_mm.to_numpy(), 3)
    sel = (rel > lo) & (rel <= hi)
    pen = meta.penetration.to_numpy()
    pens = [p for p in np.unique(pen[sel]) if np.sum(sel & (pen == p)) >= MIN_CH]
    P = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    M = {}
    for name, x in zip("DOBP", (X[:, i], X[:, 6 + j], X[:, 10 + j], P)):
        M[name] = fs.smooth(np.nanmean([np.nanmean(x[sel & (pen == p)], axis=0) for p in pens], axis=0)[None])[0]
    den = M["D"] + M["O"]
    ok = den > 0.12
    n_ch = np.median([np.sum(sel & (pen == p)) for p in pens])
    return np.where(ok, M["B"] / den, np.nan), np.where(ok, M["P"] / den, np.nan), len(pens), n_ch


def bands(axs, X, meta, rel):
    """(a) one trace per band: measured index (solid), contrast-summation reference (dashed), grey where below."""
    t_sel = (rel >= 0) & (rel <= 150)
    t = rel[t_sel]
    d = (rel >= 30) & (rel < 100)
    rows = []
    for ax, (name, (lo, hi, where)) in zip(axs, BANDS.items()):
        AI, AS, n, n_ch = band_ai(X, meta, lo, hi)
        a, p = AI[t_sel], AS[t_sel]
        ax.fill_between(t, a, p, where=a < p, color=fs.INK2, alpha=0.35, lw=0)
        ax.plot(t, p, color=fs.MUTED, lw=0.9, ls=fs.SUM_LS)
        ax.plot(t, a, color=fs.INK, lw=1.1)
        gap = np.nanmean(AI[d]) - np.nanmean(AS[d])
        ax.text(0.99, 0.95, f"30–100 ms: {np.nanmean(AI[d]):.2f}; {fs.SUM}, {np.nanmean(AS[d]):.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK2)
        ax.set_title(f"{name}: {where}; {n} penetrations", fontsize=6.5, loc="left")
        ax.axvspan(30, 100, color=fs.RULE, alpha=0.3, lw=0, zorder=0)
        ax.set_xlim(0, 150)
        ax.set_ylim(0.35, 1.05)
        ax.set_yticks([0.4, 0.6, 0.8, 1.0])
        ax.set_ylabel("additivity index")
        if ax is not axs[-1]:
            ax.tick_params(labelbottom=False)
        rows.append(dict(band=name, n_pen=n, channels_median=n_ch, AI_dip=np.nanmean(AI[d]), AI_sum_dip=np.nanmean(AS[d]),
                         gap=gap, AI_early=np.nanmean(AI[(rel >= 0) & (rel < 20)]),
                         AI_sum_early=np.nanmean(AS[(rel >= 0) & (rel < 20)])))
    axs[-1].set_xlabel("time from response onset (ms)")
    fs.line_key(axs[0], [("measured, both eyes 45%", fs.INK, "-")], x=0.0, y=1.36, dy=0.13, length=0.08)
    fs.line_key(axs[0], [(fs.SUM, fs.MUTED, fs.SUM_LS)], x=0.52, y=1.36, dy=0.13, length=0.08)
    return pd.DataFrame(rows)


def ridgeline(ax, rel, z, AI, AS, label, gap=0.5):
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
    ax.text(-26, -gap * (len(z) - 1) / 2 + 0.6, label, rotation=90, ha="center",
            va="center", fontsize=7, color=fs.INK)
    ax.set_title("additivity index, both eyes 45%", fontsize=6.5, loc="left")
    fs.line_key(ax, [("measured", fs.INK, "-"), (fs.SUM, fs.MUTED, fs.SUM_LS)], x=0.5, y=1.0, dy=0.035, length=0.08)


def by_group(ax, G, V, cfg):
    D = G[G.window == "dip"] if "window" in G else G
    rows = []
    for k, g in enumerate(cfg["groups"]):
        x = D[D.grp == g]
        for col, colr, dy, mfc in (("AI", fs.INK, -0.12, fs.INK), ("AI_sum", fs.MUTED, 0.12, "white")):
            v = x[col].dropna()
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9)
            ax.plot([md], [k + dy], marker="o", ms=3.6, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            rows.append(dict(grp=g, series=col, n=len(v), median=md))
    ax.set_yticks(range(3), list(cfg["names"]))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(0.3, 0.9)
    ax.set_xlabel("additivity index, 30–100 ms")
    v2 = V.set_index("test").loc[cfg["test"]]
    ax.text(1.03, 0.98, f"measured\n{fs.SUM} (open)\n\n{cfg['contrast']}:\n{int(v2.n_in_direction)}/{int(v2.n)} "
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
    version = sys.argv[1] if len(sys.argv) > 1 else "v3"
    cfg = DEPTH[version]
    fs.apply()
    Z = dict(np.load(out_root() / "figures" / "data" / f"{cfg['data']}.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / f"{cfg['data']}_meta.parquet")
    if version == "v3":                                     # first sheet only (the file also holds sink-anchor channels)
        keep = (meta.grp != "").to_numpy()
        Z["onset"], meta = Z["onset"][keep], meta[keep].reset_index(drop=True)
    G, V = pd.read_csv(QC / cfg["files"][0]), pd.read_csv(QC / cfg["files"][1])
    C, Vc = pd.read_csv(QC / "v4_V_cells.csv"), pd.read_csv(QC / "v4_V.csv")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    if version == "v3":
        Ab = bands([ax["a0"], ax["a1"], ax["a2"]], Z["onset"], meta, Z["rel"])
    else:                                                   # the sink version keeps the ridgeline
        for k in ("a1", "a2"):
            ax[k].remove()
        ax["a"] = ax.pop("a0")
        ax["a"].set_position(ax["a"].get_position().from_bounds(*np.array(axes_box(fig, *POS_A_SINK))))
        z, AI, AS = depth_ai(Z["onset"], meta, cfg)
        ridgeline(ax["a"], Z["rel"], z, AI, AS, cfg["label"])
        Ab = None
    B = by_group(ax["b"], G, V, cfg)
    oc = ONSET_DEPTH[version]
    Zs = np.load(out_root() / "figures" / "data" / f"{oc['data']}.npz")
    ms = pd.read_parquet(out_root() / "figures" / "data" / f"{oc['data']}_meta.parquet")
    O = laminar_onset(ax["c"], ax["ccb"], Zs, ms, oc)
    ax["c"].set_ylabel(cfg["label"].replace(" (mm;", "\n(mm;"))
    if version == "v3":                                     # median depth of the earliest sink (check Z1)
        zs = pd.read_csv(QC / "z_onset_sink.csv").set_index("test").iloc[0]["median"]
        ax["c"].plot([0, 3], [zs, zs], color=fs.INK, lw=1.2, clip_on=False, zorder=5)
        ax["c"].text(3.8, zs, "earliest sink (median)", fontsize=5.5, va="center", color=fs.INK2)
    pa = "a0" if version == "v3" else "a"
    for k, lab, x in ((pa, "a", -0.2), ("b", "b", -0.32), ("c", "c", -0.32)):
        fs.panel_label(ax[k], lab, x=x)
    print(fs.save(fig, "fig3_normalization_ai" + ("" if version == "v3" else "_sink")))
    pd.set_option("display.width", 200)
    rel = Z["rel"]
    if Ab is not None:
        print("(a) bands:")
        print(Ab.round(3).to_string(index=False))
        return
    for k in range(len(z)):
        e, d = (rel >= 0) & (rel < 20), (rel >= 30) & (rel < 100)
        print(f"(a) {z[k]:+.1f} mm: 0-20 ms {np.nanmean(AI[k, e]):.2f} (sum {np.nanmean(AS[k, e]):.2f}); "
              f"30-100 ms {np.nanmean(AI[k, d]):.2f} (sum {np.nanmean(AS[k, d]):.2f})")
    print(B.round(3).to_string(index=False))
    print("(c) onset per depth bin (ms after stimulus onset):")
    print(O.round(1).to_string(index=False))
    print(V.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
