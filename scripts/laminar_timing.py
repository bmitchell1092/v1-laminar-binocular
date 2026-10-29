"""Figure 4: when does each depth respond and normalize, relative to the earliest current sink? (docs/09)

Usage: .venv\\Scripts\\python scripts\\laminar_timing.py
Input: <out_root>/binocular/k_windows.parquet (equivalent_contrast.py) and docs/qc/dichoptic_csd_sinks.csv
(dichoptic_csd.py, usable sinks only). Dioptic .225 and .45 pooled. Channels are grouped by position relative to
their penetration's earliest grating sink: at the sink (within 2 contacts), above it, below it.
* Onset: each channel's half-max onset (dominant eye, .9); penetration median per group.
* Normalization time: the penetration-median L (window-specific eye weights) per onset-aligned window; the first
  window starting >= 20 ms after onset with L < 0.5, at its window centre (t_rel); on the flash clock, t_rel + the
  group's onset (t_abs).
Paired Wilcoxon tests (above vs sink, below vs sink) over penetrations.
Writes <out_root>/binocular/laminar_timing.parquet and docs/qc/laminar_timing.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]
GROUPS = ["above sink", "at sink", "below sink"]
CENTRES = {0: 5, 10: 15, 20: 27.5, 35: 42.5, 50: 62.5, 75: 87.5, 100: 125, 150: 175}


def main():
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    W = W[W.c.isin([0.225, 0.45]) & W.penetration.isin(S.index[S.usable])].copy()
    W["rel"] = W.depth - W.penetration.map(S.sink_depth_below_toc)          # contacts below the sink (+ = deeper)
    W["grp"] = np.select([W.rel.abs() <= 2, W.rel < -2, W.rel > 2], ["at sink", "above sink", "below sink"], default="")
    on = W.drop_duplicates(["penetration", "site_id"]).groupby(["penetration", "grp"]).onset.median().unstack()
    per = W.groupby(["penetration", "grp", "w0"]).L_w.median().unstack("w0")
    rows = []
    for (pen, grp), r in per.iterrows():
        below = [w for w in r.index if w >= 20 and np.isfinite(r[w]) and r[w] < 0.5]
        t_rel = CENTRES[below[0]] if below else np.nan
        rows.append(dict(penetration=pen, monkey=pen[7], grp=grp, onset=on.loc[pen, grp], t_rel=t_rel,
                         t_abs=t_rel + on.loc[pen, grp], **{f"L_{w}": r[w] for w in r.index}))
    T = pd.DataFrame(rows)
    T.to_parquet(out_root() / "binocular" / "laminar_timing.parquet")
    out = []
    for var in ("onset", "t_rel", "t_abs"):
        piv = T.pivot_table(index="penetration", columns="grp", values=var)
        for g in ("above sink", "below sink"):
            ok = piv[["at sink", g]].dropna()
            d = ok[g] - ok["at sink"]
            p = stats.wilcoxon(ok[g], ok["at sink"]).pvalue if len(ok) > 5 and (d != 0).any() else np.nan
            out.append(dict(variable=var, comparison=f"{g} - at sink", n_pen=len(ok), median_sink=ok["at sink"].median(),
                            median_other=ok[g].median(), median_diff=d.median(), frac_later=(d > 0).mean(), p=p))
    R = pd.DataFrame(out)
    pd.set_option("display.width", 200)
    print(R.round(3).to_string(index=False))
    R.to_csv(REPO / "docs" / "qc" / "laminar_timing.csv", index=False)


if __name__ == "__main__":
    main()
