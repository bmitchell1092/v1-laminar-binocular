"""Supplementary Figure S2: the additivity index over time at every depth (formerly Fig 3a; moved to the supplement
after Blake's comment on draft V4, 2026-11-11, with three matched depth bands in its place in Fig 3a).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_ridgeline.py
AI at 45% in each eye over time from response onset, one trace per 0.2 mm depth bin from the input/deep boundary of
laminar alignment v3 (>= 10 penetrations; ratio of the population means of channel means; 5 ms running mean), with the
contrast-summation reference (fig3_normalization_ai.depth_ai, ridgeline).
Inputs: <out_root>/figures/data/spacetime_v3.npz, spacetime_v3_meta.parquet (x_depth_v3.py).
Writes <out_root>/figures/figS_ridgeline.pdf and .png.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig3_normalization_ai import DEPTH, axes_mm, depth_ai, ridgeline  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")


def main():
    cfg = DEPTH["v3"]
    fs.apply()
    Z = dict(np.load(out_root() / "figures" / "data" / f"{cfg['data']}.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / f"{cfg['data']}_meta.parquet")
    keep = (meta.grp != "").to_numpy()
    X, meta = Z["onset"][keep], meta[keep].reset_index(drop=True)
    fig = fs.figure("onehalf", 100)
    fig.set_layout_engine("none")
    ax = axes_mm(fig, 22, 10, 70, 82)
    z, AI, AS = depth_ai(X, meta, cfg)
    ridgeline(ax, Z["rel"], z, AI, AS, cfg["label"])
    print(fs.save(fig, "figS_ridgeline"))
    rel = Z["rel"]
    e, d = (rel >= 0) & (rel < 20), (rel >= 30) & (rel < 100)
    for k in range(len(z)):
        print(f"{z[k]:+.1f} mm: 0-20 ms {np.nanmean(AI[k, e]):.2f} (sum {np.nanmean(AS[k, e]):.2f}); "
              f"30-100 ms {np.nanmean(AI[k, d]):.2f} (sum {np.nanmean(AS[k, d]):.2f})")


if __name__ == "__main__":
    main()
