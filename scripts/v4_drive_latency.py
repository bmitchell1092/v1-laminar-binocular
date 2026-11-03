"""V4 check H: does onset latency follow drive across the whole 6 x 6 design? (docs/09, registered 2026-11-03)

Usage: .venv\\Scripts\\python scripts\\v4_drive_latency.py [n_workers]
Channels, presentations, smoothing, scaling (dominant-eye peak at 90% = 1), splits and onset estimator as
v4_latency.py, for all 35 non-blank cells. Drive = a cell's peak response in 10-250 ms. Per channel, the monocular
drive-latency curve = least-squares line of onset on log peak through the eight monocular cells (each eye alone at
11-90%, >= 5 defined), read at a binocular peak up to 1.5 x the largest monocular peak. "If the responses added" =
the onset of D(c_D, t) + O(c_O, t). Binocular cells for H1-H2: both contrasts >= 11%.
H1 onset(both) - onset(D + O) < 0; H2 onset(both) - curve(peak(both)) within +-1 ms; H3 Spearman rho(peak, onset) < 0
per penetration. Reported: peak-based linearity index at 11, 22.5, 45% and binocular peak - dominant peak at c (1 + a).
Writes <out_root>/binocular/v4_H_cells.parquet (per channel and cell), docs/qc/v4_H.csv (tests) and
docs/qc/v4_H_cells.csv (population medians per cell, for the figure).
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

from v4_latency import MARGIN, SPLITS, boot_ci, onset

from v1lb.binocular import CONTRASTS, eye_weight, invert_crf
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
MONO = (0.11, 0.225, 0.45, 0.9)
DIOPTIC = (0.11, 0.225, 0.45)


def peak(x, t, lo=10, hi=250):
    w = (t >= lo) & (t < hi)
    return float(np.nanmax(x[w]))


def curve(mono: pd.DataFrame):
    """Least-squares line of onset on log peak through the monocular cells, or None."""
    ok = mono.onset.notna() & (mono.peak > 0)
    if ok.sum() < 5:
        return None
    b, a = np.polyfit(np.log(mono.peak[ok]), mono.onset[ok], 1)
    top = mono.peak[ok].max()
    return lambda p: np.where((p > 0) & (p <= 1.5 * top), a + b * np.log(np.maximum(p, 1e-9)), np.nan)


def penetration(task):
    pen, driven = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return []
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
    rows = []
    for q, k in enumerate(chans):
        m = M["all"][..., q]
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        mo = m if dom_contra else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0 or np.isnan(mo[:, 0]).any() or np.isnan(mo[0, :]).any():
            continue
        a = eye_weight((mo[:, 0] / pk)[:, w].mean(axis=1), (mo[0, :] / pk)[:, w].mean(axis=1))
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            for i, cd in enumerate(CONTRASTS):
                for j, co in enumerate(CONTRASTS):
                    if i == 0 and j == 0 or not np.isfinite(x[i, j]).all():
                        continue
                    summed = x[i, 0] + x[0, j] if i and j else np.full(t.size, np.nan)
                    rows.append(dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), split=s,
                                     a=a, c_dom=cd, c_oth=co, peak=peak(x[i, j], t), onset=onset(x[i, j], t),
                                     peak_added=peak(summed, t) if i and j else np.nan,
                                     onset_added=onset(summed, t) if i and j else np.nan))
    return rows


def per_channel(C: pd.DataFrame) -> pd.DataFrame:
    """Curve residuals, summed-response differences and the peak-based index, per channel and split."""
    out = []
    for (pen, site, s), g in C.groupby(["penetration", "site_id", "split"]):
        mono = g[((g.c_oth == 0) & g.c_dom.isin(MONO)) | ((g.c_dom == 0) & g.c_oth.isin(MONO))]
        f = curve(mono)
        b = g[g.c_dom.isin(MONO) & g.c_oth.isin(MONO)]
        r = dict(penetration=pen, site_id=site, split=s, monkey=g.monkey.iloc[0],
                 vs_added=float((b.onset - b.onset_added).median()),
                 vs_curve=float((b.onset - f(b.peak.to_numpy())).median()) if f is not None else np.nan,
                 rho=stats.spearmanr(g.peak, g.onset, nan_policy="omit").statistic if g.onset.notna().sum() >= 6
                 else np.nan)
        pd_ = g[g.c_oth == 0].set_index("c_dom").peak.reindex(CONTRASTS).to_numpy(copy=True)
        po_ = g[g.c_dom == 0].set_index("c_oth").peak.reindex(CONTRASTS).to_numpy(copy=True)
        pd_[0] = po_[0] = 0.0
        if np.isfinite(pd_).all() and np.isfinite(po_).all():
            a_pk = eye_weight(pd_, po_)
            for c in DIOPTIC:
                pb = g[np.isclose(g.c_dom, c) & np.isclose(g.c_oth, c)].peak
                if pb.empty:
                    continue
                kc, cens = invert_crf(pd_, np.array([pb.iloc[0]]))
                k = (CONTRASTS[-1] if cens[0] == 1 else kc[0]) / c
                pred = np.interp(np.log(c * (1 + a_pk)), np.log(CONTRASTS[1:]), np.maximum.accumulate(pd_[1:]))
                r.update({f"Lpeak_{c:g}": (k - 1) / a_pk if a_pk > 0 else np.nan, f"cens_{c:g}": int(cens[0] == 1),
                          f"peakdiff_{c:g}": float(pb.iloc[0] - pred)})
        out.append(r)
    return pd.DataFrame(out)


def test(H: pd.DataFrame, col: str, kind: str) -> dict:
    per = {s: H[H.split == s].groupby("penetration")[col].median().dropna() for s in SPLITS}
    v = per["all"]
    lo, hi = boot_ci(v)
    row = dict(test=col, n=len(v), median=float(v.median()), ci_lo=lo, ci_hi=hi,
               median_odd=float(per["odd"].median()), median_even=float(per["even"].median()))
    if kind == "less":
        row.update(p=stats.wilcoxon(v, alternative="less").pvalue, n_in_direction=int((v < 0).sum()),
                   holds=bool(stats.wilcoxon(v, alternative="less").pvalue < .05 and row["median_odd"] < 0
                              and row["median_even"] < 0))
    elif kind == "equiv":
        row.update(p=stats.wilcoxon(v).pvalue, n_in_direction=int((v.abs() <= MARGIN).sum()),
                   holds=bool(lo >= -MARGIN and hi <= MARGIN and abs(row["median_odd"]) <= MARGIN
                              and abs(row["median_even"]) <= MARGIN))
    elif kind == "rho":
        n_neg = int((v < 0).sum())
        p = stats.binomtest(n_neg, len(v), 0.5, alternative="greater").pvalue
        row.update(p=p, n_in_direction=n_neg, holds=bool(n_neg >= 0.8 * len(v) and p < .05 and row["median_odd"] < 0
                                                         and row["median_even"] < 0))
    else:
        row.update(p=stats.wilcoxon(v).pvalue, n_in_direction=int((v > 0).sum()), holds=np.nan)
    return row


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        C = pd.DataFrame([r for done in ex.map(penetration, [(p, driven) for p in pens]) for r in done])
    C.to_parquet(out_root() / "binocular" / "v4_H_cells.parquet")
    H = per_channel(C)
    rows = [dict(test_name="H1 both eyes - if the responses added", **test(H, "vs_added", "less")),
            dict(test_name="H2 both eyes - monocular drive curve", **test(H, "vs_curve", "equiv")),
            dict(test_name="H3 Spearman rho, peak vs onset", **test(H, "rho", "rho"))]
    for c in DIOPTIC:
        rows.append(dict(test_name=f"(reported) peak-based linearity index, {c:g}", **test(H, f"Lpeak_{c:g}", "")))
        rows.append(dict(test_name=f"(reported) binocular peak - dominant peak at c(1+a), {c:g}",
                         **test(H, f"peakdiff_{c:g}", "")))
        cens = H[H.split == "all"][f"cens_{c:g}"]
        rows[-2]["censored_fraction"] = float(cens.mean())
    T = pd.DataFrame(rows)
    T.to_csv(REPO / "docs" / "qc" / "v4_H.csv", index=False)
    # population medians per cell (median over channels per penetration, then over penetrations), onset relative to
    # the dominant eye at 90%
    A = C[C.split == "all"].copy()
    ref = A[(A.c_dom == 0.9) & (A.c_oth == 0)].set_index(["penetration", "site_id"]).onset
    A["onset_rel"] = A.onset - ref.reindex(pd.MultiIndex.from_frame(A[["penetration", "site_id"]])).to_numpy()
    A["onset_added_rel"] = A.onset_added - ref.reindex(pd.MultiIndex.from_frame(A[["penetration", "site_id"]])).to_numpy()
    per = A.groupby(["penetration", "c_dom", "c_oth"])[["peak", "onset_rel", "peak_added", "onset_added_rel"]].median()
    cells = per.groupby(["c_dom", "c_oth"]).agg(["median", "count"])
    cells.columns = [f"{a}_{b}" for a, b in cells.columns]
    cells = cells.reset_index()
    cells.to_csv(REPO / "docs" / "qc" / "v4_H_cells.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"{A.penetration.nunique()} penetrations, {A.site_id.groupby(A.penetration).nunique().sum()} channels")
    print(T.round(4).to_string(index=False))
    print(cells.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
