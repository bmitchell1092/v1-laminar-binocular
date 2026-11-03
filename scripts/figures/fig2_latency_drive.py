"""Figure 2 (simplified set, docs/09): response onset follows drive, and the second eye is drive, at every depth.

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_latency_drive.py
(a) response onset (half-max of each condition's own peak) against response strength (that peak) for all 35 conditions
    of the design; medians over channels per penetration, then over penetrations (check H; fig1_setup.drive_latency);
(b) the first 50 ms at each depth: MUA for both eyes at 45% each, each depth bin scaled to its own peak (0-100 ms),
    with the median onset at each depth (check G: half-max of each channel's own peak, interpolated; median over the
    channels of a penetration in the bin, then over penetrations) for both eyes, the dominant eye alone and the other
    eye alone at 45%; population per 0.2 mm depth bin relative to the earliest sink (>= 10
    penetrations), stimulus clock (spacetime_data.py).
Statistics: checks G (binocular advance by depth group, laminar onset order), H and U (hierarchical).
Inputs: docs/qc/v4_H_cells.csv, v4_H.csv, v4_G.csv, v4_U.csv; <out_root>/binocular/v4_H_cells.parquet;
figures/data/spacetime.npz, spacetime_meta.parquet.
Writes <out_root>/figures/fig2_latency_drive.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig1_setup import drive_latency  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 66
BIN = 2
EDGES = np.arange(-9, 14, BIN)
MIN_PEN = 10
POS = {"a": (16, 11, 62, 48), "b": (100, 11, 56, 48), "bcb": (159, 11, 2, 48)}
ROWS = {"both": 12, "dominant": list(CONTRASTS).index(0.45), "other": 8}      # spacetime.npz rows at 45%


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def channel_onsets(z):
    """Per depth bin: median over penetrations of each penetration's median channel onset (check G estimator)."""
    C = pd.read_parquet(out_root() / "binocular" / "v4_G_channels.parquet")
    A = C[(C.split == "all") & C.rel_sink.notna()].copy()
    A["bin"] = np.digitize(A.rel_sink.to_numpy(), EDGES) - 1
    centres = np.round((EDGES[:-1] + BIN / 2) * 0.1, 2)
    per = A.groupby(["penetration", "bin"])[["dom_0.45", "oth_0.45", "bin_0.45"]].median().reset_index()
    per["depth"] = centres[per.bin.clip(0, len(centres) - 1)]
    med = per.groupby("depth")[["dom_0.45", "oth_0.45", "bin_0.45"]].median().reindex(np.round(z, 2))
    return {"both": med["bin_0.45"].to_numpy(), "dominant": med["dom_0.45"].to_numpy(),
            "other": med["oth_0.45"].to_numpy()}


def laminar_onset(ax, cax, Z, meta):
    t, X = Z["t"], Z["stim"]
    b = np.digitize(meta.rel_sink.to_numpy(), EDGES) - 1
    centres = (EDGES[:-1] + BIN / 2) * 0.1
    pen = meta.penetration.to_numpy()
    curves = {k: [] for k in ROWS}
    z = []
    for k in range(len(centres)):
        sel = b == k
        pens = np.unique(pen[sel])
        if len(pens) < MIN_PEN:
            continue
        z.append(centres[k])
        for name, r in ROWS.items():
            curves[name].append(np.nanmean([np.nanmean(X[sel & (pen == p), r], axis=0) for p in pens], axis=0))
    z = np.array(z)
    C = {k: np.array(v) for k, v in curves.items()}
    w = (t >= 0) & (t < 100)
    img = C["both"] / np.nanmax(C["both"][:, w], axis=1, keepdims=True)
    im = ax.imshow(img, aspect="auto", origin="lower", cmap="Greys", vmin=0, vmax=1, interpolation="bilinear",
                   extent=(t[0] - 0.5, t[-1] + 0.5, z[0] - 0.1, z[-1] + 0.1), rasterized=True)
    on = channel_onsets(z)
    for name, col, ls, lab in (("other", fs.COND["other"], "-", "other eye alone"),
                               ("dominant", fs.COND["dominant"], "-", "dominant eye alone"),
                               ("both", fs.INK, "-", "both eyes")):
        ax.plot(on[name], z, color="white", lw=2.4, zorder=3)
        ax.plot(on[name], z, color=col, lw=1.1, ls=ls, marker="o", ms=1.8, zorder=4, label=lab)
    ax.set_ylim(z[-1] + 0.1, z[0] - 0.1)
    ax.set_xlim(0, 50)
    ax.axhline(0, color=fs.INK2, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to\nearliest sink (mm)")
    ax.set_title("onset at each depth, 45% in each eye", fontsize=6.5)
    fs.line_key(ax, [("both eyes", fs.INK, "-"), ("dominant eye alone", fs.COND["dominant"], "-"),
                     ("other eye alone", fs.COND["other"], "-")], x=0.03, y=0.97, dy=0.075, length=0.07)
    cb = ax.figure.colorbar(im, cax=cax, ticks=[0, 0.5, 1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("MUA, both eyes\n(each depth's peak = 1)", fontsize=6)
    return pd.DataFrame(dict(depth_mm=z, **{f"onset_{k}": v for k, v in on.items()}))


def main():
    fs.apply()
    Z = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    G = pd.read_csv(QC / "v4_G.csv").set_index("test")
    U = pd.read_csv(QC / "v4_U.csv").set_index("test")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    dl = drive_latency(ax["a"], show_added=False, notes=False)
    u4 = U.loc["U4 rank correlation of onset and peak < 0"]
    g45 = G.loc["G1 both eyes earlier, 0.45"]
    ax["a"].text(0.97, 0.64, f"35 conditions; rank correlation {dl['rho']:.2f}\n"
                 f"per penetration {u4.estimate:.2f}, {fs.fmt_p(u4.p_one_sided)}\n"
                 f"both eyes {abs(g45['median']):.1f} ms earlier than\nthe dominant eye alone at 45%",
                 transform=ax["a"].transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK2, linespacing=1.15)
    O = laminar_onset(ax["b"], ax["bcb"], Z, meta)
    for k in ("a", "b"):
        fs.panel_label(ax[k], k, x=-0.22)
    print(fs.save(fig, "fig2_latency_drive"))
    pd.set_option("display.width", 200)
    print(f"(a) rank correlation over the 35 conditions {dl['rho']:.3f}; hierarchical per-penetration "
          f"{u4.estimate:.3f} (p {u4.p_one_sided:.2g}); both eyes - dominant at 45% {g45['median']:.2f} ms")
    print("(b) onset per depth bin (ms after stimulus onset):")
    print(O.round(1).to_string(index=False))
    print(G.loc[[t for t in G.index if t.startswith(("G3 both eyes earlier", "G4"))],
                ["n", "median", "n_in_direction", "p", "holds"]].round(4).to_string())


if __name__ == "__main__":
    main()
