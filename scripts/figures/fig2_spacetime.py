"""Figure 2 (space-time figure set, docs/09): the space-time resolution of binocular summation.

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_spacetime.py
(a) the linearity index (share of the other eye's contrast added; 1 = the contrasts add, 0 = dominant eye alone) in
    windows after each channel's response onset, for both eyes at 11, 22.5 and 45%: channel-average readout (the
    data of Supplementary Fig S2a, v3_null_model_index.py), median over penetrations and 95% bootstrap CI;
(b) the population current-source density (presentations with either eye >= 45%, aligned on each penetration's
    earliest sink and scaled by its peak; sink_data.py), with the median MUA half-max onset at each depth (dominant
    eye at 90%) drawn on it (penetrations with that depth, >= 10);
(c, d) what the other eye adds at 45% (both eyes - dominant eye alone), depth x time, on the stimulus clock with
    each depth's median onset drawn (c) and on each channel's own response-onset clock (d); population per 0.2 mm
    depth bin (>= 10 penetrations), 5 ms running mean for display (spacetime_data.py).
Inputs: <out_root>/binocular/v3_A1_L.parquet; figures/data/sink.npz, spacetime.npz, spacetime_meta.parquet.
Writes <out_root>/figures/fig2_spacetime.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
H_MM = 104
CS = (0.11, 0.225, 0.45)
WIN = {0: 10, 10: 20, 20: 35, 35: 50, 50: 75, 75: 100, 100: 150}
EARLY_W, LATE_W = (0, 10), (50, 75, 100)
BIN = 2
EDGES = np.arange(-9, 14, BIN)
MIN_PEN = 10
POS = {"a": (16, 60, 52, 36), "b": (88, 60, 46, 36), "bcb": (136, 60, 2, 36),
       "c": (16, 11, 66, 36), "d": (92, 11, 66, 36), "dcb": (161, 11, 2, 36)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def index_over_time(ax, L):
    x = np.array([(w0 + w1) / 2 for w0, w1 in WIN.items()])
    rows = []
    for c in CS:
        g = L[(L.source == "data") & np.isclose(L.c, c)].pivot(index="penetration", columns="w0", values="L")[list(WIN)]
        m, lo, hi = fs.boot_ci(g.to_numpy(), stat=np.nanmedian)
        ax.fill_between(x, lo, hi, color=fs.CONTRAST[c], alpha=0.09, lw=0)
        ax.plot(x, m, color=fs.CONTRAST[c], lw=1.1, marker="o", ms=2.4)
        early, late = g[list(EARLY_W)].mean(axis=1), g[list(LATE_W)].mean(axis=1)
        d = (early - late).dropna()
        rows.append(dict(c=c, n=len(d), early=float(early.median()), late=float(late.median()),
                         n_drop=int((d > 0).sum()), p=stats.wilcoxon(d, alternative="greater").pvalue,
                         **{f"w{k}": round(float(v), 2) for k, v in zip(WIN, m)}))
    ax.axhline(1, color=fs.MUTED, lw=0.8, ls=fs.SUM_LS, zorder=0)
    fs.zero_line(ax)
    ax.set_xlim(0, 130)
    ax.set_ylim(-1.0, 1.6)
    ax.set_xlabel("time from response onset (ms, window centre)")
    ax.set_ylabel(fs.LINEARITY)
    ax.text(129, 1.04, "1: contrasts add", ha="right", va="bottom", fontsize=5.5, color=fs.INK2)
    ax.text(129, 0.04, "0: dominant eye alone", ha="right", va="bottom", fontsize=5.5, color=fs.INK2)
    for k, c in enumerate(CS):
        ax.text(0.98, 0.32 - k * 0.075, f"both eyes {fs.contrast_label(c)}", transform=ax.transAxes, ha="right",
                fontsize=5.5, color=fs.CONTRAST[c])
    return pd.DataFrame(rows)


def onset_by_depth(meta):
    off = np.arange(-8, 13)
    per = meta.assign(o=np.round(meta.rel_sink).astype(int)).groupby(["penetration", "o"]).onset.median()
    on = per.groupby(level="o").agg(["median", "count"]).reindex(off)
    return off, on


def csd_map(ax, cax, Z, meta):
    M = np.nanmean(Z["maps"], axis=0)
    off, t = Z["offsets"] * 0.1, Z["t_map"]
    lim = np.nanpercentile(np.abs(M[:, (t > 20) & (t < 150)]), 99)
    dt, dz = t[1] - t[0], off[1] - off[0]
    im = ax.imshow(M, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, interpolation="bilinear", aspect="auto", origin="lower",
                   extent=(t[0] - dt / 2, t[-1] + dt / 2, off[0] - dz / 2, off[-1] + dz / 2), rasterized=True)
    o, on = onset_by_depth(meta)
    ok = on["count"].to_numpy() >= MIN_PEN
    ax.plot(on["median"].to_numpy()[ok], (o * 0.1)[ok], color="white", lw=2.2, zorder=3)
    ax.plot(on["median"].to_numpy()[ok], (o * 0.1)[ok], color=fs.INK, lw=0.9, marker="o", ms=1.6, zorder=4)
    ax.set_ylim(off.max(), off.min())
    ax.set_xlim(-10, 150)
    ax.set_xticks([0, 50, 100, 150])
    ax.plot([-10], [0], marker=">", ms=4, color=fs.INK, clip_on=False, zorder=5)
    ax.text(44, (o * 0.1)[ok][0] + 0.02, "spiking onset\nat each depth", ha="left", va="top", fontsize=5.5,
            color=fs.INK, linespacing=1.0, bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to\nearliest sink (mm)")
    ax.set_title(f"current-source density, {Z['maps'].shape[0]} penetrations", fontsize=6.5)
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return pd.DataFrame(dict(offset_mm=o * 0.1, onset=on["median"].to_numpy(), n=on["count"].to_numpy()))[ok]


def added_by_depth(X, meta, c=0.45):
    i = list(CONTRASTS).index(c)
    added = X[:, 10 + 2] - X[:, i]                               # both eyes 45% (row 12) - dominant eye 45%
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
    im = ax.imshow(A, aspect="auto", origin="lower", cmap=fs.ADDED_CMAP, vmin=-0.25, vmax=0.25,
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


def main():
    fs.apply()
    L = pd.read_parquet(out_root() / "binocular" / "v3_A1_L.parquet")
    Zs = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    Zt = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    A = index_over_time(ax["a"], L)
    O = csd_map(ax["b"], ax["bcb"], Zs, meta)
    res = {}
    for k, clock, tk in (("c", "stim", "t"), ("d", "onset", "rel")):
        z, M, onset = added_by_depth(Zt[clock], meta)
        im = added_map(ax[k], Zt[tk], z, M, onset, clock)
        res[clock] = (Zt[tk], z, M, onset)
    ax["c"].set_ylabel("depth relative to\nearliest sink (mm)")
    ax["d"].tick_params(labelleft=False)
    cb = fig.colorbar(im, cax=ax["dcb"], ticks=[-0.2, -0.1, 0, 0.1, 0.2])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("added by the other eye\n(both eyes − dominant eye alone, 45%)", fontsize=6)
    for k in ("a", "b", "c", "d"):
        fs.panel_label(ax[k], k, x=-0.3 if k == "a" else -0.02)
    print(fs.save(fig, "fig2_spacetime"))
    pd.set_option("display.width", 250)
    print("(a) linearity index, medians per window (window starts), and early (0-20) vs late (50-150 ms) per "
          "penetration, one-sided Wilcoxon:")
    print(A.round(4).to_string(index=False))
    print("(b) median MUA half-max onset (dominant eye 90%) by depth relative to the sink (ms):")
    print(O.round(1).to_string(index=False))
    for clock, (t, z, M, onset) in res.items():
        print(f"({'c' if clock == 'stim' else 'd'}) {clock} clock: time of the peak of what the other eye adds, "
              f"per depth bin:", dict(zip(np.round(z, 1), t[np.nanargmax(np.where(np.isfinite(M), M, -9), axis=1)])))


if __name__ == "__main__":
    main()
