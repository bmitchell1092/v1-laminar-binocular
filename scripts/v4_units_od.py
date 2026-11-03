"""V4 check T: single units by ocular dominance, hierarchically (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_units_od.py [n_workers]
Lenient driven units (as v3_units_index.py): spike counts in 1 ms bins around the photodiode onset of clean
presentations, condition means (>= 3), blank subtracted, Gaussian sigma 3 ms (v3_dominance_cv.cell_means), the unit's
dominant eye, eye weight a, |ODI| and half-max onset from all presentations, scaled by the dominant-eye peak at 90%
(40-250 ms). At 45% in each eye: added A = R(c, c) - R_D(c) and predicted P = R_D(c (1 + a)) - R_D(c) in 0-20 and
30-100 ms after the unit's onset, for all / odd / even presentations. Mixed models with units as observations, a
random intercept per session and a penetration variance component, monkey as a fixed sum-coded effect.
T1 A(0-20) > 0; T2 A - P (30-100) < 0; T3 T1 in the more monocular half of units (|ODI| above the median).
Holds if one-sided p < .05 and the same sign on odd and even presentations.
Writes docs/qc/v4_T.csv and docs/qc/v4_T_units.csv.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from v3_dominance_cv import cell_means, selection

from v1lb.binocular import CONTRASTS, response_at
from v1lb.config import out_root
from v1lb.dataset import bin_spikes, clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
WINDOWS = {"early": (0, 20), "dip": (30, 100)}
C45 = 0.45


def penetration(task):
    pen, driven = task
    upath = out_root() / "units" / f"{pen}_units.parquet"
    if not upath.exists():
        return []
    U = pd.read_parquet(upath)
    U["lenient"] = U.well_isolated.fillna(False).astype(bool) | ((U.isi_contamination < 0.5) & (U.snr >= 3))
    U = U[U.lenient.to_numpy() & np.array([(pen, str(c)) in driven for c in U.cluster_id])].reset_index(drop=True)
    if U.empty:
        return []
    S = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    win = (int(t[0]), int(t[-1]) + 1)
    ticks = [S.loc[S.cluster_id == cid, "tick"].to_numpy() for cid in U.cluster_id]
    rates = bin_spikes(ticks, P.pd_on_tick.to_numpy(), window_ms=win).astype(float) * 1000.0
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    M = {"all": cell_means(rates, cc, ci, keep), "odd": cell_means(rates, cc, ci, keep & (idx % 2 == 1)),
         "even": cell_means(rates, cc, ci, keep & (idx % 2 == 0))}
    w_full = (t >= 40) & (t < 250)
    i = list(CONTRASTS).index(C45)
    rows = []
    for q, u in U.iterrows():
        s = selection(M["all"], q, t, w_full)
        if s is None:
            continue
        rd_all = M["all"][:, 0, :, q] if s["contra_dom"] else M["all"][0, :, :, q]
        scale = np.nanmax(rd_all[-1, w_full])
        if not np.isfinite(scale) or scale <= 0:
            continue
        for split, X in M.items():
            rd = (X[:, 0, :, q] if s["contra_dom"] else X[0, :, :, q]) / scale
            rb = X[i, i, :, q] / scale
            if np.isnan(rd).any() or np.isnan(rb).any():
                continue
            added = rb - rd[i]
            pred = response_at(rd, C45 * (1 + s["a"])) - rd[i]
            r = dict(penetration=pen, monkey=pen[7], session=pen[:8], cluster_id=int(u.cluster_id), split=split,
                     strict=bool(u.well_isolated) if pd.notna(u.well_isolated) else False, odi=abs(s["odi"]),
                     a=s["a"], depth=float(u.depth_below_toc))
            for wn, (w0, w1) in WINDOWS.items():
                sel = (t >= s["onset"] + w0) & (t < s["onset"] + w1)
                r[f"A_{wn}"] = float(np.nanmean(added[sel]))
                r[f"S_{wn}"] = float(np.nanmean((added - pred)[sel]))
            rows.append(r)
    return rows


def intercept(D, col):
    D = D.dropna(subset=[col])
    m = smf.mixedlm(f"{col} ~ 1 + C(monkey, Sum)", D, groups=D.session, re_formula="1",
                    vc_formula={"pen": "0 + C(penetration)"}).fit(reml=True)
    return float(m.params["Intercept"]), float(m.pvalues["Intercept"]), len(D), m


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "unit") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        Uu = pd.DataFrame([r for res in ex.map(penetration, [(p, driven) for p in pens]) for r in res])
    Uu.to_csv(REPO / "docs" / "qc" / "v4_T_units.csv", index=False)
    med = Uu[Uu.split == "all"].odi.median()
    tests = {"T1 the other eye adds at onset, A(0-20 ms) > 0": ("A_early", lambda d: d, 1),
             "T2 single-neuron normalization, A - P (30-100 ms) < 0": ("S_dip", lambda d: d, -1),
             f"T3 T1 in the more monocular half (|ODI| > {med:.2f})": ("A_early", lambda d: d[d.odi > med], 1),
             "(reported) A - P (0-20 ms)": ("S_early", lambda d: d, 1),
             "(reported) T2 in the more monocular half": ("S_dip", lambda d: d[d.odi > med], -1)}
    rows = []
    for name, (col, sub, want) in tests.items():
        res = {s: intercept(sub(Uu[Uu.split == s]), col) for s in ("all", "odd", "even")}
        est, p, n, m = res["all"]
        p1 = p / 2 if np.sign(est) == want else 1 - p / 2
        rows.append(dict(test=name, n_units=n, n_penetrations=sub(Uu[Uu.split == "all"]).penetration.nunique(),
                         estimate=est, p_one_sided=p1, monkey_effect=float(m.params.get("C(monkey, Sum)[S.B]",
                                                                                        np.nan)),
                         p_monkey=float(m.pvalues.get("C(monkey, Sum)[S.B]", np.nan)),
                         odd=res["odd"][0], even=res["even"][0],
                         holds=bool(p1 < .05 and np.sign(res["odd"][0]) == want and np.sign(res["even"][0]) == want)
                         if not name.startswith("(reported)") else np.nan))
    T = pd.DataFrame(rows)
    T.to_csv(REPO / "docs" / "qc" / "v4_T.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    A = Uu[Uu.split == "all"]
    print(f"{A.cluster_id.size} units in {A.penetration.nunique()} penetrations ({A.groupby('monkey').size().to_dict()});"
          f" median |ODI| {med:.2f}")
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
