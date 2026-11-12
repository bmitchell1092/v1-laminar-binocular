"""Supplementary Figure S2: the additivity index over time at every depth (formerly Fig 3a; moved to the supplement
after Blake's comment on draft V4, 2026-11-11, with three matched depth bands in its place in Fig 3a).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_ridgeline.py
AI at 45% in each eye over time from response onset, one trace per 0.2 mm bin of height above the layer 4C/5 anchor
of laminar alignment v4, scaled by each penetration's squish factor (2026-11-12; main laminar set, tiers 1 + 2;
>= 10 penetrations per bin; ratio of the population means of channel means; 5 ms running mean), with the
contrast-summation reference (fig3_normalization_ai.depth_ai, ridgeline). `v3` as argument: the earlier version on
v3 depth, all penetrations.
Inputs: ad_compartments_v4.load_spacetime (X's space-time arrays, 50 um probes recomputed).
Writes <out_root>/figures/figS_ridgeline.pdf and .png.
"""

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig3_normalization_ai import DEPTH, axes_mm, depth_ai, ridgeline  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import metadata_dir, out_root  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ad_compartments_v4 import alignment, load_spacetime, relabel  # noqa: E402

warnings.filterwarnings("ignore")


def main():
    version = sys.argv[1] if len(sys.argv) > 1 else "v4"
    fs.apply()
    if version == "v3":
        cfg = DEPTH["v3"]
        Z = dict(np.load(out_root() / "figures" / "data" / f"{cfg['data']}.npz"))
        meta = pd.read_parquet(out_root() / "figures" / "data" / f"{cfg['data']}_meta.parquet")
        keep = (meta.grp != "").to_numpy()
        X, meta, rel, sign, label = Z["onset"][keep], meta[keep].reset_index(drop=True), Z["rel"], 1, cfg["label"]
    else:
        V = alignment()
        V = V[V.s.notna()]
        on, _, M, rel, _ = load_spacetime(V)
        M = M.assign(_i=np.arange(len(M)))
        meta = relabel(M[M.penetration.isin(V[V.tier.isin((1, 2))].penetration)], V)
        meta = meta[meta.coarse != ""].reset_index(drop=True)
        X = on[meta._i.to_numpy()]
        meta["depth_v4"] = -meta.h_scaled                      # binned as depth (+ deeper), shown as height
        cfg = dict(rel="depth_v4", per_contact=10, edges=np.arange(-14, 9, 2))
        import fig3_normalization_ai as F3
        F3.MIN_PEN = 8                                            # 10 penetrations in the main set (v5): >= 8 per bin
        sign, label = -1, "height above the layer 4C/5 boundary (mm, scaled)"
    fig = fs.figure("onehalf", 100)
    fig.set_layout_engine("none")
    ax = axes_mm(fig, 22, 10, 70, 82)
    z, AI, AS = depth_ai(X, meta, cfg)
    z = sign * z
    ridgeline(ax, rel, z, AI, AS, label)
    print(fs.save(fig, "figS_ridgeline" + ("_v3" if version == "v3" else "")))
    e, d = (rel >= 0) & (rel < 20), (rel >= 30) & (rel < 100)
    print(f"{meta.penetration.nunique()} penetrations")
    for k in range(len(z)):
        print(f"{z[k]:+.1f} mm: 0-20 ms {np.nanmean(AI[k, e]):.2f} (sum {np.nanmean(AS[k, e]):.2f}); "
              f"30-100 ms {np.nanmean(AI[k, d]):.2f} (sum {np.nanmean(AS[k, d]):.2f})")

if __name__ == "__main__":
    main()
