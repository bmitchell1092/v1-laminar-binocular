"""First-pass feasibility of the top literature questions (reports/V1 binocular contrast normalization.md).

Usage: .venv\\Scripts\\python scripts\\binocular_feasibility.py
Needs <out_root>/dichoptic/*.parquet (run_dichoptic.py) and <out_root>/crf/*_crf.parquet (driven flags).

Per driven MUA channel (driven in any monocular condition, sustained or transient):
* ocular dominance ODI = (R_contra - R_ipsi) / (R_contra + R_ipsi) at .9 contrast, full window; balanced = |ODI| < 0.3
* octave ratio O(c) = R(c,c) / R_D(2c) with R_D the dominant eye's monocular response, c in {.055, .11, .225, .45}:
  early (linear) summation predicts 1 on balanced channels; normalization before summation < 1, constant;
  two-stage / Ding-Sperling < 1 and falling with contrast (all derived, see report Q2)
* additivity index AI(c) = R(c,c) / (R(c,0) + R(0,c)): 1 = linear sum, 0.5 = ocularity invariance for balanced
  channels (report Q7)
Summaries by depth below top of cortex (5-contact bins) and window, with bootstrap 95% CIs over penetrations
(penetration-level resampling, so neighbouring channels do not inflate precision).
Writes docs/qc/binocular_feasibility.csv.
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]
C = [0.055, 0.11, 0.225, 0.45]


def cell(m, cc, ci):
    v = m.get((cc, ci), np.nan)
    return v


def per_site(g: pd.DataFrame, odi_full: float) -> dict:
    m = g.set_index(["contrast_contra", "contrast_ipsi"])["mean"].to_dict()
    dom = "contra" if odi_full >= 0 else "ipsi"
    out = {}
    for c in C:
        rdd = cell(m, c, c)
        r2 = cell(m, 2 * c if c != 0.11 else 0.225, 0) if dom == "contra" else cell(m, 0, 2 * c if c != 0.11 else 0.225)
        if c == 0.45:
            r2 = cell(m, 0.9, 0) if dom == "contra" else cell(m, 0, 0.9)
        if c == 0.225:
            r2 = cell(m, 0.45, 0) if dom == "contra" else cell(m, 0, 0.45)
        if c == 0.055:
            r2 = cell(m, 0.11, 0) if dom == "contra" else cell(m, 0, 0.11)
        sum_mono = cell(m, c, 0) + cell(m, 0, c)
        out[f"oct_{c:g}"] = rdd / r2 if r2 > 0 else np.nan
        out[f"ai_{c:g}"] = rdd / sum_mono if sum_mono > 0 else np.nan
    return out


def boot_pen(df: pd.DataFrame, col: str, n: int = 2000, seed: int = 0):
    """Median of the per-penetration medians, with a bootstrap CI over penetrations."""
    per = df.groupby("penetration")[col].median().dropna()
    if len(per) < 3:
        return np.nan, np.nan, np.nan, len(per)
    rng = np.random.default_rng(seed)
    b = [np.median(rng.choice(per.values, len(per))) for _ in range(n)]
    return float(np.median(per)), *np.percentile(b, [2.5, 97.5]), len(per)


def main():
    D = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "dichoptic" / "*.parquet"))], ignore_index=True)
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono")].groupby(["penetration", "site_id"]).driven.any()
    D = D[D.site == "mua"].join(drv, on=["penetration", "site_id"])
    D = D[D.driven.fillna(False)]
    full = D[D.window == "full"].set_index(["penetration", "site_id", "contrast_contra", "contrast_ipsi"])["mean"]
    rows = []
    for (pen, sid), g in D.groupby(["penetration", "site_id"]):
        try:
            rc, ri = full.loc[(pen, sid, 0.9, 0.0)], full.loc[(pen, sid, 0.0, 0.9)]
        except KeyError:
            continue
        odi = (rc - ri) / (abs(rc) + abs(ri)) if (abs(rc) + abs(ri)) > 0 else np.nan
        for w, gw in g.groupby("window"):
            rows.append(dict(penetration=pen, site_id=sid, window=w, depth_below_toc=gw.depth_below_toc.iloc[0],
                             depth_l45=gw.depth_l45.iloc[0], odi=odi) | per_site(gw, odi))
    S = pd.DataFrame(rows)
    S["balanced"] = S.odi.abs() < 0.3
    S["depth_bin"] = (S.depth_below_toc // 5 * 5)
    print(f"driven MUA channels: {S.site_id.groupby(S.penetration).nunique().sum()} in {S.penetration.nunique()} "
          f"penetrations; balanced (|ODI|<0.3): {S[S.window == 'full'].balanced.mean():.0%}; "
          f"median |ODI| {S[S.window == 'full'].odi.abs().median():.2f}")
    out = []
    for w in ("transient", "sustained"):
        for subset, sel in (("balanced", S.balanced), ("all", S.balanced | ~S.balanced)):
            g = S[(S.window == w) & sel]
            for c in C:
                for stat in ("oct", "ai"):
                    med, lo, hi, n = boot_pen(g, f"{stat}_{c:g}")
                    out.append(dict(window=w, subset=subset, stat=stat, contrast=c, depth_bin="all", median=med, lo=lo,
                                    hi=hi, n_pen=n, n_sites=int(g[f"{stat}_{c:g}"].notna().sum())))
            if subset == "balanced":
                for db, gd in g.groupby("depth_bin"):
                    for c in (0.11, 0.225):
                        med, lo, hi, n = boot_pen(gd, f"oct_{c:g}")
                        out.append(dict(window=w, subset=subset, stat="oct", contrast=c, depth_bin=int(db), median=med,
                                        lo=lo, hi=hi, n_pen=n, n_sites=int(gd[f"oct_{c:g}"].notna().sum())))
    T = pd.DataFrame(out)
    pd.set_option("display.width", 200)
    print(T[T.depth_bin == "all"].round(3).to_string(index=False))
    print(T[T.depth_bin != "all"].round(3).to_string(index=False))
    T.to_csv(REPO / "docs" / "qc" / "binocular_feasibility.csv", index=False)


if __name__ == "__main__":
    main()
