"""Is supersaturation vs depth / laminar compartment answerable with this dataset? (honest assessment)

Usage: .venv\\Scripts\\python scripts\\assess_supersaturation.py
Needs <out_root>/crf/<pen>_crf.parquet (run_crf.py). Writes docs/qc/supersaturation_assessment_*.csv and prints.

Statistic: the cross-validated supersaturation index (si_cv, `v1lb.crf.cv_supersaturation_index`), which is <= 0
in expectation for any monotonic CRF, so a positive mean is evidence of supersaturation (the plain SI is not).
Sites: MUA per in-cortex channel; single units 'strict' (well_isolated: contamination < 0.1, amplitude cutoff
< 0.1, presence >= 0.9, SNR >= 3) and 'lenient' (contamination < 0.5, same otherwise). Only driven sites.
Compartments (penetrations with an L4/5): supragranular depth_l45 > 5 contacts, granular 0 < depth_l45 <= 5,
infragranular depth_l45 <= 0 (approximate: L4 ~0.5 mm above the L4/5 boundary). Depth: contacts below the top of
cortex, bins of 5.
Tests: (1) site-level mixed model si_cv ~ compartment + (1 | penetration); (2) conservative: per penetration,
mean si_cv per compartment, supragranular - granular and infragranular - granular differences tested across
penetrations (neighbouring MUA channels share signal, so sites are not independent). Power: one-sample t-test on
the per-penetration differences, from their observed SD and number of penetrations.
"""

import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from v1lb.config import out_root

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]


def compartment(d):
    return np.select([d > 5, d > 0, d <= 0], ["supra", "granular", "infra"], default="")


def power_one_sample(delta, sd, n, alpha=0.05):
    if n < 2 or not np.isfinite(sd) or sd <= 0:
        return np.nan
    df = n - 1
    ncp = delta / (sd / np.sqrt(n))
    tc = stats.t.ppf(1 - alpha / 2, df)
    p = stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)
    if not np.isfinite(p):                                        # nct overflows for very large ncp
        p = stats.norm.sf(tc - ncp) + stats.norm.cdf(-tc - ncp)
    return float(min(p, 1.0))


def main():
    R = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(out_root() / "crf" / "*_crf.parquet")))],
                  ignore_index=True)
    R["lenient"] = (R.site == "unit") & (R.isi_contamination < 0.5) & (R.snr >= 3)
    R["kind"] = np.where(R.site == "mua", "mua", np.where(R.well_isolated.fillna(False).astype(bool), "unit_strict",
                                                          np.where(R.lenient, "unit_lenient", "unit_other")))
    R["compartment"] = np.where(np.isfinite(R.depth_l45), compartment(R.depth_l45), "")
    R["depth_bin"] = (R.depth_below_toc // 5 * 5).astype("Int64")
    D = R[R.driven & (R.kind != "unit_other")].copy()
    D.loc[D.kind == "unit_strict", "kind"] = "unit_strict"
    # lenient set includes strict units
    D = pd.concat([D, D[D.kind == "unit_strict"].assign(kind="unit_lenient")], ignore_index=True)

    # ---- yield
    y = (R[R.kind != "unit_other"].assign(driven_any=R.driven)
         .groupby(["penetration", "kind", "site_id"]).driven_any.any().groupby(["penetration", "kind"]).sum()
         .unstack(fill_value=0))
    print("driven sites per penetration (any condition/window):\n", y.to_string(), "\ntotal:", y.sum().to_dict())

    rows, comp_rows = [], []
    for (kind, window, cond), g in D.groupby(["kind", "window", "condition"]):
        x = g.si_cv.dropna()
        w = stats.wilcoxon(x, alternative="greater").pvalue if len(x) >= 6 else np.nan
        rows.append(dict(kind=kind, window=window, condition=cond, n_sites=len(x), n_pen=g.penetration.nunique(),
                         si_cv_mean=x.mean(), si_cv_median=x.median(), p_greater_0=w,
                         frac_si_cv_gt_0_1=(x > 0.1).mean(),
                         frac_fit_supersat=((g.s > 1) & (g.s_p < 0.05)).mean(),
                         frac_peak_below_top=(g.peak_contrast_obs < 0.9).mean()))
        # compartment effects
        gc = g[g.compartment != ""]
        if gc.penetration.nunique() >= 3 and gc.compartment.nunique() >= 2:
            try:
                m = smf.mixedlm("si_cv ~ C(compartment, Treatment('granular'))", gc.dropna(subset=["si_cv"]),
                                groups="penetration").fit(reml=True)
                mm = {k.split("[T.")[-1].rstrip("]"): (m.params[k], m.pvalues[k]) for k in m.params.index if "T." in k}
            except Exception:
                mm = {}
            per = gc.groupby(["penetration", "compartment"]).si_cv.mean().unstack()
            for other in ("supra", "infra"):
                if other in per and "granular" in per:
                    d = (per[other] - per["granular"]).dropna()
                    sd = d.std(ddof=1)
                    comp_rows.append(dict(kind=kind, window=window, condition=cond, contrast=f"{other}-granular",
                                          n_pen=len(d), mean_diff=d.mean(), sd_diff=sd,
                                          p_pen=stats.ttest_1samp(d, 0).pvalue if len(d) >= 3 else np.nan,
                                          mixed_coef=mm.get(other, (np.nan, np.nan))[0],
                                          mixed_p=mm.get(other, (np.nan, np.nan))[1],
                                          power_0_05=power_one_sample(0.05, sd, len(d)),
                                          power_0_10=power_one_sample(0.10, sd, len(d)),
                                          power_0_20=power_one_sample(0.20, sd, len(d)),
                                          n_pen_for_80pct_at_0_10=next((n for n in range(3, 200)
                                                                        if power_one_sample(0.10, sd, n) >= 0.8), np.nan)))
    pop = pd.DataFrame(rows)
    comp = pd.DataFrame(comp_rows)
    pd.set_option("display.width", 250)
    print("\npopulation (driven sites):\n", pop.round(3).to_string(index=False))
    print("\ncompartment effects:\n", comp.round(3).to_string(index=False))
    dep = (D[D.depth_bin.notna()].groupby(["kind", "window", "condition", "depth_bin"]).si_cv
           .agg(["size", "mean"]).round(3))
    print("\nsi_cv by depth below top of cortex (contacts):\n", dep.unstack("depth_bin").to_string())
    qc = REPO / "docs" / "qc"
    pop.to_csv(qc / "supersaturation_assessment_population.csv", index=False)
    comp.to_csv(qc / "supersaturation_assessment_compartments.csv", index=False)
    dep.reset_index().to_csv(qc / "supersaturation_assessment_depth.csv", index=False)
    y.to_csv(qc / "supersaturation_assessment_yield.csv")


if __name__ == "__main__":
    main()
