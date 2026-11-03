"""Supplementary Figure S8: the depth axis and the shortfall from summation in space and time (moved out of the main
figures after the strength review, docs/09 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\figures\\figS8_spacetime_maps.py
(a) the depth landmark: population current-source density aligned on each penetration's earliest sink (presentations
    with either eye >= 45%, each penetration scaled by its sink's peak; sink_data.py), with the median MUA half-max
    onset at each depth (dominant eye at 90%; penetrations with that depth, >= 10). The two signals' onsets are not
    compared here (Supplementary Fig S7d);
(b) both eyes - "if the contrasts add" (the dominant eye at c (1 + a), each channel's own eye weight), depth x time on
    each channel's response clock, at 11, 22.5 and 45% in each eye; the line marks when each depth is normalized:
    the first time >= 5 ms after onset at which the other eye adds less than half of what summation predicts (not
    drawn where this is not reached by 100 ms);
(c) that time for the three contrasts (descriptive; not tested, and the earlier per-penetration timing measures were
    not reliable between halves of the presentations, docs/09 V3 analysis D).
Inputs: spacetime_data.py (spacetime.npz, spacetime_meta.parquet), sink_data.py (sink.npz).
Writes <out_root>/figures/figS8_spacetime_maps.pdf and .png; prints the values shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig3_laminar_normalization import CS, axes_mm, by_depth, csd_landmark, maps, norm_time, timing  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
H_MM = 112
POS = {"a": (18, 64, 50, 38), "acb": (70, 64, 2, 38), "c": (112, 64, 34, 38),
       "b1": (16, 10, 40, 40), "b2": (60, 10, 40, 40), "b3": (104, 10, 40, 40), "bcb": (148, 10, 2.2, 40)}


def main():
    fs.apply()
    Z = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    Zs = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    rel = Z["rel"]
    res = {c: by_depth(Z["onset"], meta, c) for c in CS}
    times = {c: norm_time(rel, res[c][2]["added"], res[c][2]["pred"]) for c in CS}
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    O = csd_landmark(ax["a"], ax["acb"], Zs, meta)
    maps([ax["b1"], ax["b2"], ax["b3"]], ax["bcb"], rel, res, times)
    timing(ax["c"], res, times)
    ax["c"].set_title("descriptive", fontsize=6.5, color=fs.INK2)
    for k in ("a", "b1", "c"):
        fs.panel_label(ax[k], k[0], x=-0.3 if k == "c" else -0.02)
    print(fs.save(fig, "figS8_spacetime_maps"))
    pd.set_option("display.width", 200)
    print("(a) median MUA half-max onset by depth relative to the sink (ms):")
    print(O.round(1).to_string(index=False))
    print("(b, c) normalized at (ms after onset) per depth bin:")
    for c in CS:
        print(f"   {fs.contrast_label(c)}:", dict(zip(np.round(res[c][0], 1), times[c])))


if __name__ == "__main__":
    main()
