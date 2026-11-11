"""Supplementary figure: each eye alone (moved from Fig 1 on 2026-11-11).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_monocular.py
(a) monocular contrast response per eye, mean MUA 40-250 ms after stimulus onset in units of each channel's
    dominant-eye peak, penetration means (95% bootstrap CI over penetrations);
(b) ocular dominance index at 90% for driven MUA channels and isolated single units.
Inputs: figure_data.py (population.npz, population_meta.parquet), run_crf.py, run_dichoptic.py.
Writes <out_root>/figures/figS_monocular.pdf and .png; prints the statistics shown.
"""

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fig1_setup import axes_mm, crfs, ocular_dominance  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet")
    fig = fs.figure("onehalf", 46)
    fig.set_layout_engine("none")
    a, b = axes_mm(fig, 14, 9, 34, 30), axes_mm(fig, 72, 9, 34, 30)
    crf, ratio = crfs(a, pop)
    odi = ocular_dominance(b)
    fs.panel_label(a, "a", x=-0.35)
    fs.panel_label(b, "b", x=-0.3)
    print(fs.save(fig, "figS_monocular"))
    for lab, (m, lo, hi) in crf.items():
        print(f"(a) {lab}: " + "; ".join(f"{100 * c:g}% {mm:.3f} [{l_:.3f}, {h:.3f}]"
                                         for c, mm, l_, h in zip(CONTRASTS[1:], m, lo, hi)))
    print(f"    other / dominant eye at 90% per penetration: median {np.median(ratio):.2f} "
          f"(IQR {np.percentile(ratio, 25):.2f}-{np.percentile(ratio, 75):.2f}); median eye weight a over "
          f"penetrations {meta.a.median():.2f} (IQR {meta.a.quantile(.25):.2f}-{meta.a.quantile(.75):.2f})")
    for s, v in odi.items():
        print(f"(b) {s}: n {len(v)}, median ODI {np.median(v):+.2f}, median |ODI| {np.median(np.abs(v)):.2f}, "
              f"|ODI| < 0.5 {np.mean(np.abs(v) < 0.5):.0%}, contralateral-dominant {np.mean(v > 0):.0%}, "
              f"|ODI| > 0.9 {np.mean(np.abs(v) > 0.9):.0%}")


if __name__ == "__main__":
    main()
