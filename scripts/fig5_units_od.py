"""Figure 5 first look: binocular modulation vs ocular dominance, early vs late, in single units and MUA.

Usage: .venv\Scripts\python scripts\fig5_units_od.py
From the 6x6 matrices (run_dichoptic.py; transient 40-140, sustained 150-250 ms), for sites driven in a monocular
condition (run_crf.py): ODI = (R_contra - R_ipsi) / (|R_contra| + |R_ipsi|) at .9, full window;
binocular modulation BMI(c) = (R(c,c) - R_D(c)) / (|R(c,c)| + |R_D(c)|), R_D = dominant eye alone at c, for c = .45, .9.
Zhang et al. 2024 predict late BMI more negative with |ODI|; Dougherty et al. 2019 early facilitation.
Mixed model BMI ~ |ODI| with a random intercept per penetration, per epoch, kind and contrast.
Writes docs/qc/fig5_units_od.csv.
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]


def main():
    D = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "dichoptic" / "*.parquet"))], ignore_index=True)
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[R.condition.str.startswith("mono") & R.driven].groupby(["penetration", "site", "site_id"]).size()
    D = D.join(drv.rename("nd"), on=["penetration", "site", "site_id"])
    D = D[D.nd.notna()]
    D["kind"] = np.where(D.site == "mua", "mua",
                         np.where(D.well_isolated.fillna(False).astype(bool), "unit_strict",
                                  np.where((D.isi_contamination < 0.5) & (D.snr >= 3), "unit_lenient", "")))
    D = D[D.kind != ""]
    rows = []
    for (pen, site, sid), g in D.groupby(["penetration", "site", "site_id"]):
        m = g.set_index(["window", "contrast_contra", "contrast_ipsi"])["mean"]
        try:
            rc, ri = m[("full", 0.9, 0.0)], m[("full", 0.0, 0.9)]
        except KeyError:
            continue
        if abs(rc) + abs(ri) <= 0:
            continue
        odi = (rc - ri) / (abs(rc) + abs(ri))
        for w in ("transient", "sustained"):
            for c in (0.45, 0.9):
                try:
                    b = m[(w, c, c)]
                    dom = m[(w, c, 0.0)] if odi >= 0 else m[(w, 0.0, c)]
                except KeyError:
                    continue
                den = abs(b) + abs(dom)
                if den > 0:
                    rows.append(dict(penetration=pen, kind=g.kind.iloc[0], site_id=sid, odi=odi, abs_odi=abs(odi),
                                     window=w, c=c, bmi=(b - dom) / den, depth=g.depth_below_toc.iloc[0]))
    T = pd.DataFrame(rows)
    T = pd.concat([T, T[T.kind == "unit_strict"].assign(kind="unit_lenient")], ignore_index=True)
    out = []
    for (kind, w, c), g in T.groupby(["kind", "window", "c"]):
        g = g[np.isfinite(g.bmi)]
        rho, p_rho = stats.spearmanr(g.abs_odi, g.bmi)
        try:
            mm = smf.mixedlm("bmi ~ abs_odi", g, groups="penetration").fit()
            slope, p_mm = mm.params["abs_odi"], mm.pvalues["abs_odi"]
        except Exception:
            slope, p_mm = np.nan, np.nan
        out.append(dict(kind=kind, window=w, c=c, n=len(g), n_pen=g.penetration.nunique(), median_bmi=g.bmi.median(),
                        frac_facilitated=(g.bmi > 0).mean(), median_abs_odi=g.abs_odi.median(), spearman=rho,
                        p_spearman=p_rho, mixed_slope=slope, p_mixed=p_mm))
    O = pd.DataFrame(out)
    pd.set_option("display.width", 220)
    print(O.round(3).to_string(index=False))
    O.to_csv(REPO / "docs" / "qc" / "fig5_units_od.csv", index=False)


if __name__ == "__main__":
    main()
