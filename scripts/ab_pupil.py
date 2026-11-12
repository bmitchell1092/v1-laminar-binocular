"""Check AB: the pupil-size control (#23) on the current measures (docs/09, registered 2026-11-12).

Usage: .venv\\Scripts\\python scripts\\ab_pupil.py [n_workers]
Clean presentations split per penetration at the median of the pre-stimulus pupil signal (median over -200...0 ms;
the analog signal has a per-monkey offset, so only the within-penetration split is used): "lower" and "higher" halves,
and "IQR" (reported). Dominant eye, scale, half-max onset and eye weight from all clean presentations (as the odd/even
halves); conditions need >= 2 presentations in a subset (>= 3 for "all").
AB1 (F1) A - P < 0 at 45%, 30-100 ms: population channel mean per penetration (driven in-cortex MUA, as
    v4_dip_recovery.py), one-sided Wilcoxon.
AB2 (Q1, Q3) x_O slope of the shortfall < 0 in 30-100 and 130-200 ms; 130-200 ms shortfall where c_O > c_D < 0:
    v4_full_design.py's mixed model on the same population means.
AB3 (X1) (A - P) deep - input > 0 at 45%, 30-100 ms: X's channels and v3 groups (x_channels.parquet).
AB4 (X4) binocular onset at 45% - dominant eye's at 45% < 0: per penetration, the median over X's channels.
The control passes if AB1-AB4 hold (p < .05, predicted sign) in the lower and in the higher half. Reported: "all" (the
reference), "IQR", higher - lower (two-sided Wilcoxon; for AB2 the window:x_O:half interaction).
Writes docs/qc/ab_pupil.csv and docs/qc/ab_pupil_units.csv.
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

from v4_full_design import fit, one_sided, q3, slope
from v4_latency import onset as interp_onset

from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
SUBSETS = ("all", "lower", "higher", "IQR")
REL = np.arange(-30, 301)
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
CS = (0.11, 0.225, 0.45, 0.9)
I45 = list(CONTRASTS).index(0.45)


def subsets(ds, P, E):
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    pu = np.nanmedian(ds.pupil.values[:, (t >= -200) & (t < 0)].astype(float), axis=1)
    ok = keep & np.isfinite(pu)
    med = np.median(pu[ok])
    q1, q3_ = np.percentile(pu[ok], [25, 75])
    return {"all": keep, "lower": ok & (pu <= med), "higher": ok & (pu > med), "IQR": ok & (pu >= q1) & (pu <= q3_)}


def channel_arrays(ds, P, masks, chans):
    """Per channel: (onset-aligned [subset, 6, 6, REL] in dominant-eye orientation and peak units, eye weight,
    stimulus-clock onsets at 45% per subset), from condition means of the subset's presentations."""
    t = ds.time.values
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    M = {}
    for s, sel in masks.items():
        nmin = 3 if s == "all" else 2
        m = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                x = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if x.sum() >= nmin:
                    m[i, j] = mua[x].mean(axis=0)
        M[s] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    out = {}
    for q in range(len(chans)):
        m = M["all"][..., q]
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        mo = m if dom_contra else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0 or np.isnan(mo[:, 0]).any() or np.isnan(mo[0, :]).any():
            continue
        mo = mo / pk
        k = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] >= 0.5))
        if not k.size:
            continue
        on = float(t[k[0]])
        a = eye_weight(mo[:, 0][:, w].mean(axis=1), mo[0, :][:, w].mean(axis=1))
        Y, adv = [], {}
        for s in SUBSETS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            Y.append(np.stack([np.interp(on + REL, t, y) for y in x.reshape(36, -1)]).reshape(6, 6, -1))
            adv[s] = interp_onset(x[I45, I45], t) - interp_onset(x[I45, 0], t)
        out[q] = (np.stack(Y), a, adv)
    return out


def penetration(task):
    pen, driven, xgrp = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    masks = subsets(ds, P, E)
    labels = [str(ch) for ch in ds.channel.values]
    pop = [k for k, ch in enumerate(labels) if ds.in_cortex.values[k] and (pen, ch) in driven]
    res = {"pen": pen, "n": {s: int(m.sum()) for s, m in masks.items()}}
    if pop:
        A = channel_arrays(ds, P, masks, pop)
        if A:
            res["pop"] = (np.nanmean([v[0] for v in A.values()], axis=0), float(np.median([v[1] for v in A.values()])))
    xs = [k for k, ch in enumerate(labels) if ch in xgrp]
    if xs:
        A = channel_arrays(ds, P, masks, xs)
        rows = []
        for q, (Y, a, adv) in A.items():
            g = xgrp[labels[xs[q]]]
            for k, s in enumerate(SUBSETS):
                dom = Y[k, I45, 0]
                sel = (REL >= 30) & (REL < 100)
                S = float(np.nanmean((Y[k, I45, I45] - response_at(Y[k, :, 0], 0.45 * (1 + a)))[sel]))
                rows.append(dict(penetration=pen, site_id=labels[xs[q]], grp=g, subset=s, S=S, adv=adv[s],
                                 dom=float(np.nanmean(dom[sel]))))
        res["chan"] = rows
    return res


def pop_rows(pen, Y, a):
    """F1 at 45% and the full-design shortfalls (v4_full_design.shortfalls) per subset."""
    f1, cells = [], []
    sel = (REL >= 30) & (REL < 100)
    for k, s in enumerate(SUBSETS):
        dom = Y[k, I45, 0]
        f1.append(dict(penetration=pen, subset=s, S=float(np.nanmean((Y[k, I45, I45] - response_at(
            Y[k, :, 0], 0.45 * (1 + a)))[sel]))))
        for cd in CS:
            for co in CS:
                ce = cd + a * co
                if ce > CONTRASTS[-1]:
                    continue
                iD, iO = list(CONTRASTS).index(cd), list(CONTRASTS).index(co)
                dd = Y[k, iD, iO] - response_at(Y[k, :, 0], ce)
                for wn, (w0, w1) in WINDOWS.items():
                    ww = (REL >= w0) & (REL < w1)
                    cells.append(dict(unit=pen, split=s, c_dom=cd, c_oth=co, window=wn, a=a,
                                      S=float(np.nanmean(dd[ww]))))
    return f1, cells


def wil(v, alt, name, subset):
    v = pd.Series(v).dropna()
    want = {"less": -1, "greater": 1}.get(alt, np.sign(v.median()))
    p = stats.wilcoxon(v, alternative=alt).pvalue
    return dict(test=name, subset=subset, n=len(v), estimate=float(v.median()),
                n_in_direction=int((np.sign(v) == want).sum()), p=p,
                holds=bool(p < .05 and np.sign(v.median()) == want) if alt != "two-sided" else np.nan)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    X = pd.read_parquet(out_root() / "binocular" / "x_channels.parquet")
    X = X[(X.split == "all") & np.isclose(X.c, 0.45) & (X.window == "dip") & (X.grp != "")]
    xg = {p: dict(zip(g.site_id, g.grp)) for p, g in X.groupby("penetration")}
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        done = list(ex.map(penetration, [(p, driven, xg.get(p, {})) for p in pens]))

    F1, cells = [], []
    for r in done:
        if "pop" in r:
            f, c = pop_rows(r["pen"], *r["pop"])
            F1 += f
            cells += c
    F1 = pd.DataFrame(F1).pivot(index="penetration", columns="subset", values="S")
    Cells = pd.DataFrame(cells)
    Cells["monkey"] = Cells.unit.str[7]
    Cells["session"] = Cells.unit.str[:8]
    Ch = pd.DataFrame([x for r in done for x in r.get("chan", [])])
    Ch.to_csv(QC / "ab_pupil_units.csv", index=False)
    Sg = Ch.groupby(["penetration", "subset", "grp"]).S.mean().unstack("grp")
    X1 = (Sg["deep"] - Sg["input"]).unstack("subset")
    X4 = Ch.groupby(["penetration", "subset"]).adv.median().unstack("subset")

    rows = []
    for s in SUBSETS:
        rows.append(wil(F1[s], "less", "AB1 (F1) A - P < 0, 45%, 30-100 ms", s))
        m = fit(Cells[Cells.split == s], "session", True)
        for wn, lab in (("dip", "30-100"), ("late", "130-200")):
            est, p = slope(m, wn)
            rows.append(dict(test=f"AB2 (Q1) x_O slope < 0, {lab} ms", subset=s, n=Cells[Cells.split == s].unit.nunique(),
                             estimate=est, p=one_sided(est, p), holds=bool(one_sided(est, p) < .05 and est < 0)))
        est, p = q3(Cells[Cells.split == s], "session", True)
        rows.append(dict(test="AB2 (Q3) 130-200 ms shortfall, c_O > c_D, < 0", subset=s,
                         n=Cells[Cells.split == s].unit.nunique(), estimate=est, p=one_sided(est, p),
                         holds=bool(one_sided(est, p) < .05 and est < 0)))
        rows.append(wil(X1[s], "greater", "AB3 (X1) (A - P) deep - input > 0", s))
        rows.append(wil(X4[s], "less", "AB4 (X4) binocular advance at 45% < 0 (ms)", s))
    for name, T in (("AB1", F1), ("AB3", X1), ("AB4", X4)):
        rows.append(wil(T["higher"] - T["lower"], "two-sided", f"(reported) {name} higher - lower", "higher - lower"))
    H = Cells[Cells.split.isin(["lower", "higher"])].dropna(subset=["S"]).copy()
    H = H.assign(xO=np.log2(H.c_oth / 0.225), xD=np.log2(H.c_dom / 0.225), hi=(H.split == "higher").astype(float))
    m = smf.mixedlm("S ~ 0 + C(window) + C(window):xO + C(window):xD + C(window):hi + C(window):xO:hi + C(monkey, Sum)",
                    H, groups=H.session, re_formula="1", vc_formula={"pen": "0 + C(unit)"}).fit(reml=True)
    for wn, lab in (("dip", "30-100"), ("late", "130-200")):
        k = f"C(window)[{wn}]:xO:hi"
        rows.append(dict(test=f"(reported) AB2 x_O slope, {lab} ms, higher - lower", subset="higher - lower",
                         estimate=float(m.params[k]), p=float(m.pvalues[k])))
        k = f"C(window)[{wn}]:hi"
        rows.append(dict(test=f"(reported) AB2 shortfall at x = 0, {lab} ms, higher - lower", subset="higher - lower",
                         estimate=float(m.params[k]), p=float(m.pvalues[k])))
    T = pd.DataFrame(rows)
    passes = T[T.subset.isin(["lower", "higher"])].holds.astype(bool).all()
    T.to_csv(QC / "ab_pupil.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.max_rows", 200)
    n = pd.DataFrame({r["pen"]: r["n"] for r in done}).T
    print("presentations per subset (median over penetrations):", n.median().to_dict())
    print(T.round(4).to_string(index=False))
    print("control passes (AB1-AB4 hold in both halves):", passes)


if __name__ == "__main__":
    main()
