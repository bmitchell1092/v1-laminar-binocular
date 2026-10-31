"""V3 analysis A4: the linearity index early and late in isolated single units (docs/09 "V3 analysis plan").

Usage: .venv\\Scripts\\python scripts\\v3_units_index.py [n_workers]
Units: driven in a monocular condition (run_crf.py, site 'unit'); lenient set = well isolated, or refractory
contamination < 0.5 and SNR >= 3 (as fig5_units_od.py); strict set = well isolated (v1lb.units). Spike counts in
1 ms bins around the photodiode onset of cleaned presentations, converted to spikes/s, then exactly the per-channel
readout of Fig 2d (v3_dominance_cv: condition means >= 3 presentations, blank-subtracted, Gaussian sigma 3 ms,
dominant eye at 90%, eye weight a, half-max onset, onset-aligned windows, early = 0-20, late = 50-200 ms).
Units are the sample (median over units, Wilcoxon early vs late). Also, for the same units, the binocular modulation
index BMI = (R(c,c) - R_D(c)) / (|R(c,c)| + |R_D(c)|) at 45% in the transient (40-140 ms) and sustained (150-250 ms)
windows after stimulus onset (the net-suppression observation of draft V2, fig5_units_od.py). A cross-validated
variant (selection on one half, index on the other) is reported as well (descriptive).
Writes <out_root>/binocular/v3_A4_units.parquet and docs/qc/v3_A4.csv.
"""

import glob
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v3_dominance_cv import EARLY, LATE, cell_means, selection, windowed  # noqa: E402

from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import bin_spikes, clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
CS = (0.225, 0.45)


def penetration(task):
    pen, driven = task
    t0 = time.time()
    upath = out_root() / "units" / f"{pen}_units.parquet"
    if not upath.exists():
        return pd.DataFrame(), pd.DataFrame()
    U = pd.read_parquet(upath)
    U["lenient"] = U.well_isolated.fillna(False).astype(bool) | ((U.isi_contamination < 0.5) & (U.snr >= 3))
    U = U[U.lenient.to_numpy() & np.array([(pen, str(c)) in driven for c in U.cluster_id])].reset_index(drop=True)
    if U.empty:
        return pd.DataFrame(), pd.DataFrame()
    S = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    win = (int(t[0]), int(t[-1]) + 1)
    ticks = [S.loc[S.cluster_id == cid, "tick"].to_numpy() for cid in U.cluster_id]
    rates = bin_spikes(ticks, P.pd_on_tick.to_numpy(), window_ms=win).astype(float) * 1000.0   # [pres, t, unit]
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    M = {"all": cell_means(rates, cc, ci, keep), "odd": cell_means(rates, cc, ci, keep & (idx % 2 == 1)),
         "even": cell_means(rates, cc, ci, keep & (idx % 2 == 0))}
    w_full = (t >= 40) & (t < 250)
    rows, bmi = [], []
    for q, u in U.iterrows():
        base = dict(penetration=pen, cluster_id=int(u.cluster_id), strict=bool(u.well_isolated),
                    depth=u.depth_below_toc, rate_hz=u.rate_hz)
        s = selection(M["all"], q, t, w_full)
        if s is None:
            continue
        for (c, w0), (L, cens) in windowed(M["all"], q, t, s["contra_dom"], s["a"], s["onset"]).items():
            rows.append(base | dict(variant="all", direction="all", c=c, w0=w0, L=L, cens=cens, a=s["a"], odi=s["odi"]))
        for X, Y in (("odd", "even"), ("even", "odd")):
            sX = selection(M[X], q, t, w_full)
            if sX is None:
                continue
            for (c, w0), (L, cens) in windowed(M[Y], q, t, sX["contra_dom"], sX["a"], sX["onset"]).items():
                rows.append(base | dict(variant="cv", direction=f"{X}->{Y}", c=c, w0=w0, L=L, cens=cens, a=sX["a"],
                                        odi=sX["odi"]))
        i45 = list(CONTRASTS).index(0.45)
        rd = M["all"][i45, 0, :, q] if s["contra_dom"] else M["all"][0, i45, :, q]
        rb = M["all"][i45, i45, :, q]
        for name, (w0, w1) in (("transient", (40, 140)), ("sustained", (150, 250))):
            sel = (t >= w0) & (t < w1)
            b, dm = rb[sel].mean(), rd[sel].mean()
            den = abs(b) + abs(dm)
            bmi.append(base | dict(window=name, bmi=(b - dm) / den if den > 0 else np.nan, odi=s["odi"]))
    print(f"{pen}: {len(U)} driven lenient units ({time.time() - t0:.0f}s)", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(bmi)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "unit") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        res = list(ex.map(penetration, [(p, driven) for p in pens]))
    C = pd.concat([r[0] for r in res], ignore_index=True)
    B = pd.concat([r[1] for r in res], ignore_index=True)
    C.to_parquet(out_root() / "binocular" / "v3_A4_units.parquet")

    rows = []
    for subset, sel in (("lenient", C.index == C.index), ("strict", C.strict)):
        for variant in ("all", "cv"):
            g = C[sel & (C.variant == variant)]
            per = g.groupby(["penetration", "cluster_id", "direction", "c", "w0"]).L.first().unstack("w0")
            el = pd.DataFrame({"early": per[EARLY].mean(axis=1), "late": per[LATE].mean(axis=1)})
            el = el.groupby(["penetration", "cluster_id", "c"]).mean()           # directions averaged per unit
            for c in CS:
                x = el.xs(c, level="c").dropna()
                if len(x) < 5:
                    rows.append(dict(subset=subset, variant=variant, c=c, n_units=len(x)))
                    continue
                rows.append(dict(subset=subset, variant=variant, c=c, n_units=len(x),
                                 n_pen=x.index.get_level_values("penetration").nunique(), early=x.early.median(),
                                 late=x.late.median(), drop=(x.early - x.late).median(),
                                 frac_drop=float((x.early > x.late).mean()), frac_late_below0=float((x.late < 0).mean()),
                                 p=stats.wilcoxon(x.early, x.late).pvalue))
    for subset, sel in (("lenient", B.index == B.index), ("strict", B.strict)):
        for w, g in B[sel].groupby("window"):
            g = g.dropna(subset=["bmi"])
            rows.append(dict(subset=subset, variant=f"BMI 45% {w}", c=0.45, n_units=len(g),
                             n_pen=g.penetration.nunique(), early=g.bmi.median(), late=np.nan, drop=np.nan,
                             frac_drop=float((g.bmi > 0).mean()), frac_late_below0=np.nan,
                             p=stats.wilcoxon(g.bmi).pvalue if len(g) >= 5 else np.nan))
    S = pd.DataFrame(rows)
    S.to_csv(REPO / "docs" / "qc" / "v3_A4.csv", index=False)
    pd.set_option("display.width", 240)
    print("(index rows: early/late medians over units, p = Wilcoxon early vs late; BMI rows: early = median BMI, "
          "frac_drop = fraction facilitated, p = Wilcoxon vs 0)")
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
