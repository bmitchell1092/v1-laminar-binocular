"""V4 check F: dip and partial recovery of the response the other eye adds (docs/09, registered 2026-11-02).

Usage: .venv\\Scripts\\python scripts\\v4_dip_recovery.py [n_workers]
Per penetration, onset-aligned channel-mean responses as figure_data.py (driven in-cortex MUA channels; dominant eye,
dominant-eye peak scale, half-max onset and eye weight a from all presentations), for all presentations and for odd
and even presentations (conditions with >= 3 presentations for all, >= 2 for a half). In each window after response
onset: added A = mean of R(c, c) - R(c, 0); predicted P = mean of R_D(c (1 + a)) - R(c, 0), R_D interpolated in log
contrast between the dominant eye's measured contrasts at each time point (a = the penetration's median eye weight).
Tests at 45% (11 and 22.5% reported), Wilcoxon signed-rank over penetrations, one-sided:
F1 dip, A - P < 0 in 30-100 ms; F2 recovery, A(130-200) - A(30-100) > 0; F3 partial, A - P < 0 in 130-200 ms;
F4 undershoot, A < 0 in 30-100 ms. A test holds if p < .05 on all presentations and the median difference has the
same sign in both halves.
Writes docs/qc/v4_F.csv (tests) and docs/qc/v4_F_windows.csv (per-penetration window means); prints both.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
REL = np.arange(-30, 301)
SPLITS = ("all", "odd", "even")
WINDOWS = {"early": (0, 20), "dip": (30, 100), "recovery": (130, 200)}
DIOPTIC = (0.11, 0.225, 0.45)


def penetration(task):
    """Onset-aligned channel means [split, 6, 6, REL] and the median eye weight, or None."""
    pen, driven = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return None
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    M = {}
    for split, sel, nmin in (("all", keep, 3), ("odd", keep & (idx % 2 == 1), 2), ("even", keep & (idx % 2 == 0), 2)):
        m = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= nmin:
                    m[i, j] = mua[s].mean(axis=0)
        M[split] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    out, a_ = {s: [] for s in SPLITS}, []
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
        a_.append(eye_weight(mo[:, 0][:, w].mean(axis=1), mo[0, :][:, w].mean(axis=1)))
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            out[s].append(np.stack([np.interp(on + REL, t, y) for y in x.reshape(36, -1)]).reshape(6, 6, -1))
    if not a_:
        return None
    return pen, np.stack([np.nanmean(out[s], axis=0) for s in SPLITS]), float(np.median(a_))


def window_rows(pen, Y, a):
    rows = []
    for k, s in enumerate(SPLITS):
        for c in DIOPTIC:
            i = list(CONTRASTS).index(c)
            dom = Y[k, i, 0]
            added = Y[k, i, i] - dom
            pred = response_at(Y[k, :, 0], c * (1 + a)) - dom
            for wname, (w0, w1) in WINDOWS.items():
                sel = (REL >= w0) & (REL < w1)
                rows.append(dict(penetration=pen, split=s, c=c, window=wname, a=a, added=float(np.nanmean(added[sel])),
                                 predicted=float(np.nanmean(pred[sel])), dominant=float(np.nanmean(dom[sel]))))
    return rows


def tests(Wd: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in DIOPTIC:
        per = {}
        for s in SPLITS:
            g = Wd[(Wd.split == s) & np.isclose(Wd.c, c)].pivot(index="penetration", columns="window",
                                                                values=["added", "predicted"])
            per[s] = {"F1 dip: added - predicted, 30-100 ms (< 0)": (g["added"]["dip"] - g["predicted"]["dip"], "less"),
                      "F2 recovery: added 130-200 - added 30-100 ms (> 0)":
                          (g["added"]["recovery"] - g["added"]["dip"], "greater"),
                      "F3 partial: added - predicted, 130-200 ms (< 0)":
                          (g["added"]["recovery"] - g["predicted"]["recovery"], "less"),
                      "F4 undershoot: added, 30-100 ms (< 0)": (g["added"]["dip"], "less")}
        for name, (d, alt) in per["all"].items():
            d = d.dropna()
            p = stats.wilcoxon(d, alternative=alt).pvalue
            want = -1 if alt == "less" else 1
            halves = {s: per[s][name][0].dropna() for s in ("odd", "even")}
            same = all(np.sign(h.median()) == want for h in halves.values())
            rows.append(dict(test=name, c=c, n=len(d), median=float(d.median()),
                             n_in_direction=int((np.sign(d) == want).sum()), p=p,
                             median_odd=float(halves["odd"].median()), median_even=float(halves["even"].median()),
                             p_odd=stats.wilcoxon(halves["odd"], alternative=alt).pvalue,
                             p_even=stats.wilcoxon(halves["even"], alternative=alt).pvalue,
                             holds=bool(p < 0.05 and same) if np.isclose(c, 0.45) else np.nan))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        done = [d for d in ex.map(penetration, [(p, driven) for p in pens]) if d is not None]
    Wd = pd.DataFrame([r for pen, Y, a in done for r in window_rows(pen, Y, a)])
    T = tests(Wd)
    Wd.to_csv(REPO / "docs" / "qc" / "v4_F_windows.csv", index=False)
    T.to_csv(REPO / "docs" / "qc" / "v4_F.csv", index=False)
    pd.set_option("display.width", 250)
    print(f"{len(done)} penetrations")
    print(Wd[Wd.split == "all"].groupby(["c", "window"])[["added", "predicted", "dominant"]].mean().round(4).to_string())
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
