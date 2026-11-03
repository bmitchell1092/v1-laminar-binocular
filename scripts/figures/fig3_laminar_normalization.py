"""Figure 3 (revised set, docs/09 checks L and S): space and time, and the laminar organization of normalization.

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_laminar_normalization.py
Population: per depth bin (0.2 mm, relative to each penetration's earliest sink), the mean over penetrations (bins
with >= 10) of the mean over its channels; 5 ms running mean for display (spacetime_data.py).
(a) what the other eye adds at 45% in each eye (both eyes - dominant eye alone), depth x time, on the stimulus clock
    with each depth's median response onset drawn, and on each channel's own response-onset clock; values below 0
    are drawn faint (a response below the dominant eye alone was not found, check F);
(b) the same at 45% as one trace per depth on the response clock: what the other eye adds (ink) and what summation
    predicts it adds (grey dashed); orange = what it adds, grey = how far it falls short of summation;
(c) both eyes - "if the contrasts add" in 30-100 ms after onset, per penetration and depth group (above / at,
    within 0.2 mm / below the earliest sink), at 22.5 and 45%; medians and interquartile ranges; check L (below - at,
    one-sided Wilcoxon) and its depth-reference variants (check S).
Supplementary Fig S8 (figS8_spacetime_maps.py) uses this script's maps of both eyes - "if the contrasts add" at
three contrasts, its descriptive normalization time and the CSD depth landmark.
Inputs: spacetime_data.py (spacetime.npz, spacetime_meta.parquet); docs/qc/v4_L_groups.csv, v4_S.csv.
Writes <out_root>/figures/fig3_laminar_normalization.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
H_MM = 124
CS = (0.11, 0.225, 0.45)
ROWS = {c: (list(CONTRASTS).index(c), 6 + i, 10 + i) for i, c in enumerate(CS)}   # dominant, other, both
BIN = 2                                                     # contacts (0.2 mm)
EDGES = np.arange(-9, 14, BIN)
MIN_PEN = 10
LIM = 0.12
POS = {"a1": (16, 76, 66, 38), "a2": (92, 76, 66, 38), "acb": (161, 76, 2, 38),
       "b": (22, 10, 64, 52), "c": (106, 14, 32, 44)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def by_depth(X, meta, c):
    """Per depth bin: penetration-averaged dominant-eye, both-eyes and summation-prediction responses."""
    iD, _, iB = ROWS[c]
    pred = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    sig = {"D": X[:, iD], "B": X[:, iB], "P": pred}
    b = np.digitize(meta.rel_sink.to_numpy(), EDGES) - 1
    centres = (EDGES[:-1] + BIN / 2) * 0.1
    out = {k: np.full((len(centres), X.shape[-1]), np.nan) for k in sig}
    npen = np.zeros(len(centres), int)
    pen = meta.penetration.to_numpy()
    for k in range(len(centres)):
        sel = b == k
        pens = np.unique(pen[sel])
        npen[k] = len(pens)
        if len(pens) < MIN_PEN:
            continue
        for name, x in sig.items():
            out[name][k] = np.nanmean([np.nanmean(x[sel & (pen == p)], axis=0) for p in pens], axis=0)
    ok = npen >= MIN_PEN
    D, B, P = (fs.smooth(out[k][ok]) for k in ("D", "B", "P"))
    return centres[ok], npen[ok], dict(added=B - D, pred=P - D, deficit=B - P)


def norm_time(rel, added, pred, t_max=100):
    out = np.full(added.shape[0], np.nan)
    for k in range(added.shape[0]):
        ok = (rel >= 5) & (rel <= t_max) & (pred[k] > 0.01)
        hit = np.flatnonzero(ok & (added[k] < 0.5 * pred[k]))
        if hit.size:
            out[k] = rel[hit[0]]
    return out


def maps(axs, cax, rel, res, times):
    for k, (ax, c) in enumerate(zip(axs, CS)):
        z, _, r = res[c]
        im = ax.imshow(r["deficit"], aspect="auto", origin="lower", cmap=fs.ADDED_CMAP, vmin=-LIM, vmax=LIM,
                       interpolation="bilinear", extent=(rel[0] - 0.5, rel[-1] + 0.5, z[0] - 0.1, z[-1] + 0.1),
                       rasterized=True)
        ax.plot(times[c], z, color=fs.INK, lw=0.9, marker="o", ms=1.8)
        ax.set_ylim(z[-1] + 0.1, z[0] - 0.1)
        ax.set_xlim(-10, 150)
        ax.set_xticks([0, 50, 100, 150])
        ax.axvline(0, color=fs.INK, lw=0.6)
        ax.axhline(0, color=fs.INK, lw=0.5, ls=(0, (2, 2)))
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
        if k:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("depth relative to\nearliest sink (mm)")
            ax.text(-7, z[0] + 0.02, "normalized", fontsize=5.5, va="top", color=fs.INK,
                    bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    axs[1].set_xlabel("time from response onset (ms)")
    cb = axs[0].figure.colorbar(im, cax=cax, ticks=[-0.1, -0.05, 0, 0.05, 0.1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("both eyes − if the contrasts add\n(grey: normalized)", fontsize=6)


def ridgeline(ax, rel, z, r, gap=0.21):
    sel = (rel >= -10) & (rel <= 150)
    t = rel[sel]
    for k in range(len(z))[::-1]:
        y0 = -k * gap
        a, p = r["added"][k, sel] + y0, r["pred"][k, sel] + y0
        ax.fill_between(t, y0, a, where=a > y0, color=fs.COND["other"], alpha=0.22, lw=0, zorder=1 + 3 * k)
        ax.fill_between(t, a, p, where=a < p, color=fs.INK2, alpha=0.38, lw=0, zorder=2 + 3 * k)
        ax.plot(t, p, color=fs.MUTED, lw=0.7, ls=fs.SUM_LS, zorder=3 + 3 * k)
        ax.plot(t, a, color=fs.INK, lw=0.8, zorder=3 + 3 * k)
        ax.plot(t, np.full(t.size, y0), color=fs.RULE, lw=0.4, zorder=0)
        ax.text(-13, y0, f"{z[k]:+.1f}", ha="right", va="center", fontsize=5.5, color=fs.INK2)
    ax.axvline(0, color=fs.INK, lw=0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-10, 150)
    ax.set_ylim(-gap * (len(z) - 0.4), 0.62)
    ax.set_xlabel("time from response onset (ms)")
    ax.text(-34, -gap * (len(z) - 1) / 2, "depth relative to earliest sink (mm)", rotation=90, ha="center",
            va="center", fontsize=7, color=fs.INK)
    fs.line_key(ax, [("added by the other eye", fs.INK, "-"), (fs.SUM, fs.MUTED, fs.SUM_LS)], x=0.42, y=1.0,
                dy=0.042, length=0.08)
    ax.text(0.42, 0.9, "grey: short of summation", transform=ax.transAxes, fontsize=5.5, color=fs.INK2,
            va="center")
    ax.set_title("both eyes 45%", fontsize=6.5, loc="left")


def strength(ax, G):
    rng = np.random.default_rng(2)
    rows = []
    for c, dy in ((0.225, 0.17), (0.45, -0.17)):
        piv = G[(G.split == "all") & np.isclose(G.c, c) & (G.window == "dip")].pivot(
            index="penetration", columns="grp", values="deficit")
        for k, g in enumerate(("above", "at", "below")):
            v = piv[g].dropna()
            y = k + dy
            ax.scatter(v.clip(-0.25, 0.15), y + rng.uniform(-0.06, 0.06, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [y, y], color=fs.CONTRAST[c], lw=0.9, zorder=2)
            ax.plot([md], [y], marker="o", ms=3.4, color=fs.CONTRAST[c], zorder=3)
            if g == "above":
                ax.text(-0.245, y, fs.contrast_label(c), ha="left", va="center", fontsize=5.5, color=fs.CONTRAST[c])
            rows.append(dict(c=c, grp=g, n=len(v), median=md, n_below=int((v < 0).sum()),
                             p_less=stats.wilcoxon(v, alternative="less").pvalue))
        d = (piv["below"] - piv["at"]).dropna()
        rows.append(dict(c=c, grp="below - at", n=len(d), median=d.median(), n_below=int((d > 0).sum()),
                         p_less=stats.wilcoxon(d, alternative="greater").pvalue))
    ax.axvline(0, color=fs.MUTED, lw=0.8, ls=fs.SUM_LS)
    ax.set_yticks(range(3), ["above", "at sink", "below"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(-0.25, 0.15)
    ax.set_xlabel("both eyes − if the contrasts add,\n30–100 ms after onset")
    R = pd.DataFrame(rows)
    r45 = R[(R.c == 0.45) & (R.grp == "below - at")].iloc[0]
    ax.text(0.03, 0.99, f"45%, below − at sink:\n{r45.n_below}/{r45.n} penetrations, {fs.fmt_p(r45.p_less)}",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)
    S = pd.read_csv(REPO / "docs" / "qc" / "v4_S.csv")
    s = {k: S[S.test.str.startswith(k)].iloc[0] for k in ("S1", "S2", "S3", "S4")}
    ax.text(1.1, 0.5, "other depth references\n(45%, 30–100 ms):\n"
            f"depth below the sink:\n   {s['S1'].estimate:+.3f} per mm, {fs.fmt_p(s['S1'].p_one_sided)}\n"
            f"depth below L4/5:\n   {s['S2'].estimate:+.3f} per mm, {fs.fmt_p(s['S2'].p_one_sided)}\n"
            f"without deep sinks:\n   {int(s['S3'].n_in_direction)}/{int(s['S3'].n_penetrations)}, "
            f"{fs.fmt_p(s['S3'].p_one_sided)}\n"
            f"sink chosen on the other half:\n   {int(s['S4'].n_in_direction)}/{int(s['S4'].n_penetrations)}, "
            f"{fs.fmt_p(s['S4'].p_one_sided)}", transform=ax.transAxes, ha="left", va="center", fontsize=5.5,
            color=fs.INK2, linespacing=1.2)
    return R


def added_by_depth(X, meta, c=0.45):
    """Per depth bin: what the other eye adds (both eyes - dominant eye alone) and the median response onset."""
    i = list(CONTRASTS).index(c)
    added = X[:, ROWS[c][2]] - X[:, i]
    b = np.digitize(meta.rel_sink.to_numpy(), EDGES) - 1
    centres = (EDGES[:-1] + BIN / 2) * 0.1
    pen = meta.penetration.to_numpy()
    out, onset = np.full((len(centres), X.shape[-1]), np.nan), np.full(len(centres), np.nan)
    for k in range(len(centres)):
        sel = b == k
        pens = np.unique(pen[sel])
        if len(pens) < MIN_PEN:
            continue
        out[k] = np.nanmean([np.nanmean(added[sel & (pen == p)], axis=0) for p in pens], axis=0)
        onset[k] = np.nanmedian([np.median(meta.onset.to_numpy()[sel & (pen == p)]) for p in pens])
    ok = np.isfinite(out).any(axis=1)
    return centres[ok], fs.smooth(out[ok]), onset[ok]


def added_map(ax, t, z, A, onset, clock):
    im = ax.imshow(A, aspect="auto", origin="lower", cmap=fs.ADDED_SOFT, vmin=-0.25, vmax=0.25,
                   interpolation="bilinear", extent=(t[0] - 0.5, t[-1] + 0.5, z[0] - 0.1, z[-1] + 0.1),
                   rasterized=True)
    ax.set_ylim(z[-1] + 0.1, z[0] - 0.1)
    ax.axhline(0, color=fs.INK, lw=0.5, ls=(0, (2, 2)))
    if clock == "stim":
        ax.plot(onset, z, color=fs.INK, lw=0.9, marker="o", ms=1.8)
        ax.set_xlim(0, 200)
        ax.set_xlabel("time from stimulus onset (ms)")
        ax.set_title("aligned to stimulus onset", fontsize=6.5)
        ax.text(onset[0] + 4, z[0] + 0.02, "response onset\nat each depth", fontsize=5.5, va="top", color=fs.INK,
                linespacing=1.0)
    else:
        ax.axvline(0, color=fs.INK, lw=0.8)
        ax.set_xlim(-10, 200)
        ax.set_xlabel("time from response onset (ms)")
        ax.set_title("aligned to each site's own response onset", fontsize=6.5)
    return im


def csd_landmark(ax, cax, Z, meta):
    """Population CSD aligned on the earliest sink (sink_data.py), with the median MUA half-max onset at each depth."""
    M = np.nanmean(Z["maps"], axis=0)
    off, t = Z["offsets"] * 0.1, Z["t_map"]
    lim = np.nanpercentile(np.abs(M[:, (t > 20) & (t < 150)]), 99)
    dt, dz = t[1] - t[0], off[1] - off[0]
    im = ax.imshow(M, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, interpolation="bilinear", aspect="auto", origin="lower",
                   extent=(t[0] - dt / 2, t[-1] + dt / 2, off[0] - dz / 2, off[-1] + dz / 2), rasterized=True)
    o = np.arange(-8, 13)
    per = meta.assign(o=np.round(meta.rel_sink).astype(int)).groupby(["penetration", "o"]).onset.median()
    on = per.groupby(level="o").agg(["median", "count"]).reindex(o)
    ok = on["count"].to_numpy() >= MIN_PEN
    ax.plot(on["median"].to_numpy()[ok], (o * 0.1)[ok], color="white", lw=2.2, zorder=3)
    ax.plot(on["median"].to_numpy()[ok], (o * 0.1)[ok], color=fs.INK, lw=0.9, marker="o", ms=1.6, zorder=4)
    ax.set_ylim(off.max(), off.min())
    ax.set_xlim(-10, 150)
    ax.plot([-10], [0], marker=">", ms=4, color=fs.INK, clip_on=False, zorder=5)
    ax.text(44, (o * 0.1)[ok][0] + 0.02, "MUA half-max onset at each depth\n(not a timing comparison with the CSD)",
            ha="left", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to\nearliest sink (mm)")
    ax.set_title(f"current-source density, {Z['maps'].shape[0]} penetrations", fontsize=6.5)
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return pd.DataFrame(dict(offset_mm=o * 0.1, onset=on["median"].to_numpy(), n=on["count"].to_numpy()))[ok]


def timing(ax, res, times):
    for c in CS:
        z = res[c][0]
        ax.plot(times[c], z, color=fs.CONTRAST[c], lw=1.0, marker="o", ms=2.4)
        ha, dx = {0.11: ("left", 1.0), 0.225: ("right", -1.0), 0.45: ("center", 0.0)}[c]
        ax.text(times[c][0] + dx, z[0] - 0.08, fs.contrast_label(c), ha=ha, va="bottom", fontsize=5.5,
                color=fs.CONTRAST[c])
    ax.set_ylim(1.3, -1.1)
    ax.set_xlim(0, 100)
    ax.set_xlabel("normalized at\n(ms after each depth's onset)")
    ax.set_ylabel("depth relative to\nearliest sink (mm)")


def main():
    fs.apply()
    Z = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    G = pd.read_csv(REPO / "docs" / "qc" / "v4_L_groups.csv")
    rel = Z["rel"]
    z45, npen, r45 = by_depth(Z["onset"], meta, 0.45)
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    peaks = {}
    for k, clock, tk in (("a1", "stim", "t"), ("a2", "onset", "rel")):
        z, A, onset = added_by_depth(Z[clock], meta)
        im = added_map(ax[k], Z[tk], z, A, onset, clock)
        peaks[clock] = dict(zip(np.round(z, 1), Z[tk][np.nanargmax(np.where(np.isfinite(A), A, -9), axis=1)]))
    ax["a1"].set_ylabel("depth relative to\nearliest sink (mm)")
    ax["a2"].tick_params(labelleft=False)
    cb = fig.colorbar(im, cax=ax["acb"], ticks=[-0.2, -0.1, 0, 0.1, 0.2])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("added by the other eye\n(both eyes − dominant eye alone, 45%)", fontsize=6)
    ridgeline(ax["b"], rel, z45, r45)
    R = strength(ax["c"], G)
    for k, letter in (("a1", "a"), ("b", "b"), ("c", "c")):
        fs.panel_label(ax[k], letter, x=-0.14 if k == "b" else -0.3 if k == "c" else -0.02)
    print(fs.save(fig, "fig3_laminar_normalization"))
    print(f"{meta.penetration.nunique()} penetrations, {len(meta)} channels; depth bins and penetrations per bin:",
          dict(zip(np.round(z45, 1), npen)))
    print("(a) time of the peak of what the other eye adds, per depth bin:", peaks)
    print("(c) both eyes - if the contrasts add, 30-100 ms (check L):")
    print(R.round(4).to_string(index=False))
    print(pd.read_csv(REPO / "docs" / "qc" / "v4_S.csv").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
