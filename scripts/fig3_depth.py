"""Figure 3: where is interocular normalization strongest? Depth below the top of cortex (primary) and
landmark-based compartments (secondary), per monkey (docs/09).

Usage: .venv\\Scripts\\python scripts\\fig3_depth.py
Input: <out_root>/binocular/k_windows.parquet (equivalent_contrast.py): per channel and onset-aligned window, the
linearity index with window-specific eye weights (L_w), dioptic .11 and .225 pooled (the primary contrasts).
Epochs: early 0-20, transition 20-50, late 50-200 ms after each channel's onset.
Compartments (penetrations with an L4/5 landmark): supragranular depth_l45 > 5, granular 0 < depth_l45 <= 5,
infragranular depth_l45 <= 0. Within-penetration comparisons (Wilcoxon), and a mixed model of late L on depth below
the top of cortex (linear + quadratic) with a random intercept per penetration, overall and per monkey.
Writes docs/qc/fig3_depth.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]


def main():
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    W = W[W.c.isin([0.11, 0.225])].copy()
    W["monkey"] = W.penetration.str[7]
    W["epoch"] = np.select([W.w0 < 20, W.w0 < 50, W.w0 >= 50], ["early", "transition", "late"], default="")
    W["comp"] = np.select([W.depth_l45 > 5, W.depth_l45 > 0, W.depth_l45 <= 0], ["supra", "granular", "infra"], default="")
    rows = []
    # (1) compartments, within penetration
    C = W[W.comp != ""].groupby(["penetration", "monkey", "comp", "epoch"]).L_w.median().unstack("epoch")
    print("compartment medians over penetrations (L, 1 = linear, 0 = invariance):")
    print(C.groupby(level="comp").median().round(2).to_string())
    for ep in ("early", "transition", "late"):
        P = C[ep].unstack("comp")
        for other in ("supra", "infra"):
            ok = P[["granular", other]].dropna()
            if len(ok) >= 6:
                d = ok["granular"] - ok[other]
                mk = ok.index.get_level_values("monkey")
                p = stats.wilcoxon(ok["granular"], ok[other]).pvalue
                rows.append(dict(test=f"{ep}: granular - {other}", n_pen=len(ok), median_diff=d.median(),
                                 frac_pos=(d > 0).mean(), frac_pos_B=(d[mk == "B"] > 0).mean(),
                                 frac_pos_J=(d[mk == "J"] > 0).mean(), p=p))
    # (2) depth below top of cortex, late L, mixed model overall and per monkey
    late = W[W.epoch == "late"].groupby(["penetration", "monkey", "site_id", "depth"]).L_w.median().reset_index()
    late = late[np.isfinite(late.L_w)].assign(d=lambda x: x.depth / 10, L=lambda x: x.L_w.clip(-3, 4))
    for label, g in [("all", late), ("B", late[late.monkey == "B"]), ("J", late[late.monkey == "J"])]:
        m = smf.mixedlm("L ~ d + I(d**2)", g, groups="penetration").fit()
        peak = -m.params["d"] / (2 * m.params["I(d ** 2)"]) if m.params["I(d ** 2)"] < 0 else np.nan
        rows.append(dict(test=f"late L ~ depth + depth^2 [{label}]", n_pen=g.penetration.nunique(),
                         median_diff=m.params["I(d ** 2)"], frac_pos=peak, frac_pos_B=m.params["d"],
                         frac_pos_J=m.pvalues["d"], p=m.pvalues["I(d ** 2)"]))
    T = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print("\n(rows 'granular - X': median_diff and fraction of penetrations with granular more linear; "
          "rows 'late L ~ depth': median_diff = quadratic coef, frac_pos = peak depth in mm, frac_pos_B = linear coef, "
          "frac_pos_J = its p, p = quadratic p)")
    print(T.round(3).to_string(index=False))
    T.to_csv(REPO / "docs" / "qc" / "fig3_depth.csv", index=False)


if __name__ == "__main__":
    main()
