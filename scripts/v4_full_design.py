"""V4 checks Q and R: the full 6 x 6 design over time, and its held-out confirmation (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_full_design.py [n_workers]
Q: per penetration, onset-aligned channel means as check F (v4_dip_recovery.penetration; all / odd / even
presentations). For every binocular cell with both contrasts >= 11% and c_D + a c_O <= 90%: shortfall
S = window mean of R(c_D, c_O) - R_D(c_D + a c_O) in 0-20, 30-100, 130-200 ms after onset. Linear mixed model
S ~ 0 + window + window:x_O + window:x_D + monkey (x = log2(c / 22.5%)), random intercept per session with a
penetration variance component (REML). Q1 the x_O slope < 0 in 30-100 and 130-200 ms; Q2 those slopes < the 0-20 ms
slope; Q3 the 130-200 ms shortfall over cells with c_O > c_D < 0. Same sign required in the odd and even refits and in
each monkey alone.
R: the dual-probe cross-recordings (output/datasets_cross; driven channels by one-way ANOVA across the dominant eye's
monocular contrasts, 40-250 ms, p < .01; >= 5 per recording): the same S and model with a recording random intercept.
R1 = Q1 in 130-200 ms; R2 = Q3.
Writes docs/qc/v4_Q.csv, v4_Q_cells.csv, v4_R.csv and v4_R_cells.csv.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from normalization_timing import eye_table
from v4_dip_recovery import REL, SPLITS, penetration

from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
CS = (0.11, 0.225, 0.45, 0.9)
CROSS = out_root().parent / "datasets_cross" if (out_root().parent / "datasets_cross").exists() else \
    out_root() / "datasets_cross"


def cross_recording(name):
    """Onset-aligned channel means [split, 6, 6, REL] and the median eye weight for one cross-recording."""
    d = CROSS
    ds = open_penetration(name, d)
    P = pd.read_parquet(d / f"{name}_presentations.parquet")
    E = eye_table(d, name, P)
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    incx = np.flatnonzero(ds.in_cortex.values)
    mua = ds.mua.values[:, :, incx].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    w = (t >= 40) & (t < 250)
    M = {}
    for split, sel, nmin in (("all", keep, 3), ("odd", keep & (idx % 2 == 1), 2), ("even", keep & (idx % 2 == 0), 2)):
        m = np.full((6, 6, t.size, len(incx)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= nmin:
                    m[i, j] = mua[s].mean(axis=0)
        M[split] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    win = mua[:, w].mean(axis=1)
    out, a_ = {s: [] for s in SPLITS}, []
    for q in range(len(incx)):
        m = M["all"][..., q]
        if np.isnan(m[:, 0]).any() or np.isnan(m[0, :]).any():
            continue
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        groups = [win[keep & (np.isclose(cc, c) & np.isclose(ci, 0) if dom_contra else np.isclose(ci, c) &
                              np.isclose(cc, 0)), q] for c in CONTRASTS]
        if min(len(g) for g in groups) < 3 or stats.f_oneway(*groups).pvalue >= 0.01:
            continue
        mo = m if dom_contra else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0:
            continue
        k = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] / pk >= 0.5))
        if not k.size:
            continue
        on = float(t[k[0]])
        a_.append(eye_weight((mo[:, 0] / pk)[:, w].mean(axis=1), (mo[0, :] / pk)[:, w].mean(axis=1)))
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            out[s].append(np.stack([np.interp(on + REL, t, y) for y in x.reshape(36, -1)]).reshape(6, 6, -1))
    if len(a_) < 5:
        return name, None, np.nan, len(a_)
    return name, np.stack([np.nanmean(out[s], axis=0) for s in SPLITS]), float(np.median(a_)), len(a_)


def shortfalls(unit, Y, a):
    rows = []
    for k, s in enumerate(SPLITS):
        for cd in CS:
            for co in CS:
                ce = cd + a * co
                if ce > CONTRASTS[-1]:
                    continue
                iD, iO = list(CONTRASTS).index(cd), list(CONTRASTS).index(co)
                d = Y[k, iD, iO] - response_at(Y[k, :, 0], ce)
                for wn, (w0, w1) in WINDOWS.items():
                    sel = (REL >= w0) & (REL < w1)
                    rows.append(dict(unit=unit, split=s, c_dom=cd, c_oth=co, window=wn, a=a,
                                     S=float(np.nanmean(d[sel]))))
    return rows


def fit(D: pd.DataFrame, groups: str, monkey: bool):
    D = D.dropna(subset=["S"]).assign(xO=np.log2(D.c_oth / 0.225), xD=np.log2(D.c_dom / 0.225))
    f = "S ~ 0 + C(window) + C(window):xO + C(window):xD" + (" + C(monkey, Sum)" if monkey else "")
    vc = {"pen": "0 + C(unit)"} if groups == "session" else None
    m = smf.mixedlm(f, D, groups=D[groups], re_formula="1", vc_formula=vc).fit(reml=True)
    return m


def slope(m, w):
    return m.params[f"C(window)[{w}]:xO"], m.pvalues[f"C(window)[{w}]:xO"]


def contrast(m, w):
    """x_O slope in window w minus that in the early window: estimate and two-sided p."""
    names = list(m.fe_params.index)
    L = np.zeros((1, len(names)))
    L[0, names.index(f"C(window)[{w}]:xO")] = 1
    L[0, names.index("C(window)[early]:xO")] = -1
    r = m.t_test(L)
    return float(np.ravel(r.effect)[0]), float(np.ravel(r.pvalue)[0])


def q3(D, groups, monkey):
    sub = D[(D.window == "late") & (D.c_oth > D.c_dom)].dropna(subset=["S"])
    f = "S ~ 1" + (" + C(monkey, Sum)" if monkey else "")
    vc = {"pen": "0 + C(unit)"} if groups == "session" else None
    m = smf.mixedlm(f, sub, groups=sub[groups], re_formula="1", vc_formula=vc).fit(reml=True)
    return float(m.params["Intercept"]), float(m.pvalues["Intercept"])


def one_sided(est, p2, want=-1):
    return p2 / 2 if np.sign(est) == want else 1 - p2 / 2


def summarise(C: pd.DataFrame, groups: str, monkey: bool, label: str):
    rows = []
    A = C[C.split == "all"]
    m = fit(A, groups, monkey)
    sub_fits = {s: fit(C[C.split == s], groups, monkey) for s in ("odd", "even")}
    per_monkey = {}
    if monkey:
        for mk in ("B", "J"):
            per_monkey[mk] = fit(A[A.monkey == mk], groups, False)
    for w in ("early", "dip", "late"):
        est, p = slope(m, w)
        row = dict(test=f"{label} x_O slope, {w}", n_obs=len(A[A.window == w]), estimate=est,
                   p_one_sided=one_sided(est, p), odd=slope(sub_fits["odd"], w)[0],
                   even=slope(sub_fits["even"], w)[0], **{f"monkey_{k}": slope(v, w)[0] for k, v in per_monkey.items()})
        if w != "early":
            signs = [row["odd"], row["even"]] + [row[f"monkey_{k}"] for k in per_monkey]
            row["holds"] = bool(row["p_one_sided"] < .05 and all(s < 0 for s in signs))
        rows.append(row)
    for w in ("dip", "late"):
        est, p = contrast(m, w)
        signs = [contrast(sub_fits[s], w)[0] for s in ("odd", "even")] + [contrast(v, w)[0] for v in per_monkey.values()]
        rows.append(dict(test=f"{label} x_O slope, {w} - early", estimate=est, p_one_sided=one_sided(est, p),
                         odd=signs[0], even=signs[1],
                         **{f"monkey_{k}": contrast(v, w)[0] for k, v in per_monkey.items()},
                         holds=bool(one_sided(est, p) < .05 and all(s < 0 for s in signs))))
    est, p = q3(A, groups, monkey)
    signs = [q3(C[C.split == s], groups, monkey)[0] for s in ("odd", "even")]
    pm = {k: q3(A[A.monkey == k], groups, False)[0] for k in per_monkey}
    rows.append(dict(test=f"{label} late shortfall, other eye stronger (c_O > c_D)",
                     n_obs=int(((A.window == "late") & (A.c_oth > A.c_dom)).sum()), estimate=est,
                     p_one_sided=one_sided(est, p), odd=signs[0], even=signs[1],
                     **{f"monkey_{k}": v for k, v in pm.items()},
                     holds=bool(one_sided(est, p) < .05 and all(s < 0 for s in signs + list(pm.values())))))
    # reported: the matched diagonal, late
    diag = A[(A.window == "late") & np.isclose(A.c_oth, A.c_dom)]
    rows.append(dict(test=f"(reported) {label} late shortfall, matched contrasts", n_obs=len(diag),
                     estimate=float(diag.S.mean())))
    return pd.DataFrame(rows), m


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    names = sorted(Path(f).name[:-3] for f in glob.glob(str(CROSS / "*__bminteroc*.nc")) if not f.endswith("_eye.nc"))
    with ProcessPoolExecutor(n_workers) as ex:
        main_done = list(ex.map(penetration, [(p, driven) for p in pens]))
        cross_done = list(ex.map(cross_recording, names))
    C = pd.DataFrame([r for d in main_done if d is not None for r in shortfalls(d[0], d[1], d[2])])
    C["monkey"] = C.unit.str[7]
    C["session"] = C.unit.str[:8]
    C.to_csv(REPO / "docs" / "qc" / "v4_Q_cells.csv", index=False)
    Q, mq = summarise(C, "session", True, "Q")
    Q.to_csv(REPO / "docs" / "qc" / "v4_Q.csv", index=False)
    X = pd.DataFrame([r for d in cross_done if d[1] is not None for r in shortfalls(d[0], d[1], d[2])])
    n_rec = {d[0]: d[3] for d in cross_done}
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(mq.summary())
    print(Q.round(4).to_string(index=False))
    print("cross-recordings, driven channels:", n_rec)
    if X.unit.nunique() >= 5:
        X["monkey"] = "J"
        X.to_csv(REPO / "docs" / "qc" / "v4_R_cells.csv", index=False)
        Rr, mr = summarise(X, "unit", False, "R")
        Rr.to_csv(REPO / "docs" / "qc" / "v4_R.csv", index=False)
        print(Rr.round(4).to_string(index=False))
    else:
        pd.DataFrame([dict(test="R inconclusive: fewer than 5 cross-recordings", n=X.unit.nunique())]).to_csv(
            REPO / "docs" / "qc" / "v4_R.csv", index=False)
        print("R inconclusive")


if __name__ == "__main__":
    main()
