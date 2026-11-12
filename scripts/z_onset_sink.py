"""Check Z (descriptive; Blake's comment on draft V4, 2026-11-11): is spiking onset earliest at the earliest sink?

Usage: .venv\\Scripts\\python scripts\\z_onset_sink.py
The V4 text said onsets were earliest 0.1-0.5 mm below the input/deep boundary. Blake asked whether the earliest sink
still marks the initial input and whether spiking begins earliest below it. Not registered: run after the question,
on the onsets of check X (laminar alignment v3, first-sheet channels, half-max onsets interpolated; x_depth_v3.py).
Z1 where the earliest monocular sink lies relative to the boundary, per penetration;
Z2 within penetrations, onset in a band below the earliest sink minus onset at the sink (+-0.15 mm), for 0.2-0.4 and
   0.5-0.8 mm below (bands of channel centres (0.15, 0.45] and (0.45, 0.85] mm), and the same above; dominant eye at 90% and 45%, both eyes at 45%; two-sided Wilcoxon over penetrations, with
   the medians of odd and even presentations;
Z3 the population onset profile in 0.2 mm bins from the sink and from the boundary (median over penetrations of each
   penetration's median channel onset, and relative to that penetration's median; bins with >= 10 penetrations).
Inputs: <out_root>/binocular/x_latency.parquet. Writes docs/qc/z_onset_sink.csv and z_onset_profile.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
COLS = {"dom_90": "dominant eye 90%", "dom_45": "dominant eye 45%", "bin_45": "both eyes 45%"}
AT = 0.15                                        # mm from the sink: "at" |d| <= AT
BANDS = {"0.2-0.4 mm": (0.15, 0.45), "0.5-0.8 mm": (0.45, 0.85)}
MIN_PEN = 10


def load():
    C = pd.read_parquet(out_root() / "binocular" / "x_latency.parquet")
    return C[C.rel_anchor_mm.notna() & C.rel_sink_mm.notna() & (C.grp != "")].copy()


def band_medians(d, col, sel):
    return d[sel].groupby("penetration")[col].median()


def paired(C, col, side, band):
    near, far = BANDS[band]
    rows = {}
    for split in ("all", "odd", "even"):
        d = C[C.split == split].dropna(subset=[col])
        z = d.rel_sink_mm if side == "below" else -d.rel_sink_mm
        at = band_medians(d, col, d.rel_sink_mm.abs() <= AT)
        other = band_medians(d, col, (z > near) & (z <= far))
        rows[split] = (other - at).dropna()
    x = rows["all"]
    return dict(test=f"Z2 onset {band} {side} - at the earliest sink, {COLS[col]}", n=len(x),
                median=x.median(), n_earlier=int((x < 0).sum()), p=stats.wilcoxon(x).pvalue,
                median_odd=rows["odd"].median(), median_even=rows["even"].median())


def profile(C, col, ref):
    d = C[C.split == "all"].dropna(subset=[col]).copy()
    d["bin"] = np.round(np.floor(d[ref] / 0.2 + 1e-9) * 0.2 + 0.1, 2)
    per = d.groupby(["penetration", "bin"])[col].median().reset_index()
    per["rel"] = per[col] - per.groupby("penetration")[col].transform("median")
    t = per.groupby("bin").agg(n=("rel", "size"), onset=(col, "median"), rel=("rel", "median")).reset_index()
    return t[t.n >= MIN_PEN].assign(series=col, reference=ref.replace("rel_", "").replace("_mm", ""))


def main():
    C = load()
    A = C[C.split == "all"]
    s = A.groupby("penetration").apply(lambda d: (d.rel_anchor_mm - d.rel_sink_mm).median())
    in_input = int(((s >= -0.5) & (s <= 0)).sum())
    rows = [dict(test="Z1 earliest sink relative to the boundary (mm; negative = above)", n=len(s), median=s.median(),
                 q25=s.quantile(0.25), q75=s.quantile(0.75), n_in_input_group=in_input)]
    for col in COLS:
        for side in ("below", "above"):
            for band in BANDS:
                rows.append(paired(C, col, side, band))
    T = pd.DataFrame(rows)
    T.to_csv(QC / "z_onset_sink.csv", index=False)
    P = pd.concat([profile(C, col, ref) for col in COLS for ref in ("rel_sink_mm", "rel_anchor_mm")])
    P.to_csv(QC / "z_onset_profile.csv", index=False)
    pd.set_option("display.width", 200)
    print(T.round(3).to_string(index=False))
    for (col, ref), g in P.groupby(["series", "reference"], sort=False):
        print(f"\n{col} by depth from the {ref}:")
        print(g.set_index("bin")[["n", "onset", "rel"]].round(2).T.to_string())


if __name__ == "__main__":
    main()
