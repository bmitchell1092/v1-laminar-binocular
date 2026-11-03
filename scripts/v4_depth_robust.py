"""V4 check S: depth robustness of the laminar difference L2 (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_depth_robust.py
Per channel (v4_L_channels.parquet): the shortfall from summation at 45%, 30-100 ms after onset (added - predicted).
S1 continuous depth relative to the earliest sink, shortfall ~ depth (mm) + penetration random intercept, slope > 0;
S2 the same with depth relative to the L4/5 landmark (laminar alignment v2: l45_minus_toc); S3 L2 (below - at the
sink) excluding penetrations whose earliest sink is > 0.5 mm below L4/5 or > 1.5 mm below the top of cortex; S4 L2
with the sink chosen on one half of the presentations (the dichoptic_csd.py rule, dataset CSD, either eye >= 45%) and
the shortfall measured on the other half, both directions averaged per penetration and group. One-sided tests.
Writes docs/qc/v4_S.csv.
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from v4_units_sink import earliest_sink

from v1lb.config import metadata_dir, out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]


def channel_table():
    C = pd.read_parquet(out_root() / "binocular" / "v4_L_channels.parquet")
    C = C[np.isclose(C.c, 0.45) & (C.window == "dip")].assign(S=lambda d: d.added - d.predicted)
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    A = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration")
    C["depth"] = C.rel_sink + C.penetration.map(S.sink_depth_below_toc)               # contacts below the top
    C["rel_l45"] = C.depth - C.penetration.map(A.l45_minus_toc)
    return C, S


def slope_test(D, x):
    D = D.dropna(subset=[x, "S"]).assign(xmm=lambda d: d[x] * 0.1)
    m = smf.mixedlm("S ~ xmm", D, groups=D.penetration).fit(reml=True)
    est, p = float(m.params["xmm"]), float(m.pvalues["xmm"])
    return dict(n_channels=len(D), n_penetrations=D.penetration.nunique(), estimate=est,
                p_one_sided=p / 2 if est > 0 else 1 - p / 2)


def l2(G):
    piv = G.pivot(index="penetration", columns="grp", values="S")
    d = (piv["below"] - piv["at"]).dropna()
    return dict(n_penetrations=len(d), estimate=float(d.median()), n_in_direction=int((d > 0).sum()),
                p_one_sided=stats.wilcoxon(d, alternative="greater").pvalue)


def cv_groups(C):
    """Per penetration and direction: the sink chosen on one half; groups and shortfall from the other."""
    rows = []
    for pen, g in C.groupby("penetration"):
        ds = open_penetration(pen)
        d = out_root() / "datasets"
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = clean_presentations(ds, P, E)
        t = ds.time.values
        cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
        idx = np.arange(len(P))
        csd = ds.csd.values
        incx = ds.in_cortex.values.astype(bool)
        base = (t >= -50) & (t < 0)
        for X, Y in (("odd", "even"), ("even", "odd")):
            sx = keep & (np.maximum(cc, ci) >= 0.45) & (idx % 2 == (1 if X == "odd" else 0))
            cm = csd[sx].astype(float).mean(axis=0)
            cm -= cm[base].mean(axis=0)
            ch = earliest_sink(cm, t, incx & np.isfinite(cm).all(axis=0))
            sink = float(ds.depth_below_toc.values[ch])
            h = g[g.split == Y].copy()
            rel = h.depth - sink
            h["grp"] = np.where(rel.abs() <= 2, "at", np.where(rel < -2, "above", "below"))
            for grp, v in h.groupby("grp").S:
                rows.append(dict(penetration=pen, direction=f"{X}>{Y}", grp=grp, S=float(v.mean())))
    R = pd.DataFrame(rows)
    return R.groupby(["penetration", "grp"]).S.mean().reset_index()


def main():
    C, S = channel_table()
    A = C[C.split == "all"]
    rows = [dict(test="S1 shortfall ~ depth relative to the sink (mm), slope > 0", **slope_test(A, "rel_sink")),
            dict(test="S2 shortfall ~ depth relative to L4/5 (mm), slope > 0", **slope_test(A, "rel_l45"))]
    deep = S.index[(S.sink_depth_l45 < -5) | (S.sink_depth_below_toc > 15)]
    G = A.groupby(["penetration", "grp"]).S.mean().reset_index()
    rows.append(dict(test=f"S3 L2 without deep-sink penetrations ({', '.join(deep)})",
                     **l2(G[~G.penetration.isin(deep)])))
    rows.append(dict(test="S4 L2 with a cross-validated sink", **l2(cv_groups(C))))
    rows.append(dict(test="(reference) L2 as registered", **l2(G)))
    T = pd.DataFrame(rows)
    T["holds"] = T.p_one_sided < .05
    T.to_csv(REPO / "docs" / "qc" / "v4_S.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
