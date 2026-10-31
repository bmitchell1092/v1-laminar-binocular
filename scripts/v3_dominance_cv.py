"""V3 analyses A2 (cross-validated dominance) and A3 (ocular-dominance strata) (docs/09 "V3 analysis plan").

Usage: .venv\\Scripts\\python scripts\\v3_dominance_cv.py [n_workers]
Per-channel readout of Fig 2d (equivalent_contrast.py): driven in-cortex MUA channels, cleaned presentations,
condition means (>= 3 presentations), blank-subtracted, Gaussian sigma 3 ms; onset-aligned windows; k from the
dominant eye's monocular responses in the same window (top-censored k = 0.9 / c); L = (k - 1) / a with a from the
full 40-250 ms window; per penetration the median over channels; early = mean of windows 0-10, 10-20; late = mean
of 50-75, 75-100, 100-150, 150-200 ms after onset.
Variants (selection = dominant eye, eye weight a and onset; analysis = the responses the index is read from):
* all: selection and analysis on all presentations (reproduces Fig 2d);
* same: selection and analysis on the same half (odd or even);
* cv: selection on one half, analysis on the other (the registered A2 estimate);
* cv_dom: dominant eye from the other half, a and onset from the analysis half (dominance alone).
Half-based variants: per penetration, the two directions (odd -> even, even -> odd) are averaged.
A3: channels split into terciles of |ODI| measured on the selection half (pooled over channels and directions);
index from the analysis half (cv).
Writes <out_root>/binocular/v3_A2_channels.parquet, v3_A2_penetrations.parquet, v3_A3_penetrations.parquet and
docs/qc/v3_A2.csv, docs/qc/v3_A3.csv; prints the registered tests.
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
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent))
from equivalent_contrast import WINDOWS_REL, onset_halfmax  # noqa: E402

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
DIOPTIC = [0.055, 0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]
CRIT_C = 0.45
MARGIN = 0.25


def cell_means(mua, cc, ci, sel):
    M = np.full((6, 6, mua.shape[1], mua.shape[2]), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = sel & np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= 3:
                M[i, j] = mua[s].mean(axis=0)
    return gaussian_filter1d(M - M[0, 0], 3, axis=2)


def selection(M, q, t, w_full, contra_dom=None):
    """Dominant eye (unless given), eye weight a, ODI and onset of channel q from condition means M."""
    rc, ri = M[:, 0, :, q], M[0, :, :, q]
    if np.isnan(rc).any() or np.isnan(ri).any():
        return None
    mc, mi = rc[-1, w_full].mean(), ri[-1, w_full].mean()
    if contra_dom is None:
        contra_dom = mc >= mi
    rd, rn = (rc, ri) if contra_dom else (ri, rc)
    on = onset_halfmax(rd[-1], t)
    if not np.isfinite(on):
        return None
    return dict(contra_dom=bool(contra_dom), a=eye_weight(rd[:, w_full].mean(axis=1), rn[:, w_full].mean(axis=1)),
                odi=(mc - mi) / (abs(mc) + abs(mi)), onset=on)


def windowed(M, q, t, contra_dom, a, on):
    """L per (contrast, window start) read from condition means M of channel q, with the given selection."""
    rc, ri = M[:, 0, :, q], M[0, :, :, q]
    if np.isnan(rc).any() or np.isnan(ri).any():
        return {}
    rd = rc if contra_dom else ri
    out = {}
    for c in DIOPTIC:
        i = list(CONTRASTS).index(c)
        rb = M[i, i, :, q]
        for w0, w1 in WINDOWS_REL:
            sel = (t >= on + w0) & (t < on + w1)
            kw, cw = equivalent_contrast(rd[:, sel].mean(axis=1)[:, None], np.array([rb[sel].mean()]), c)
            kw = float(CONTRASTS[-1] / c) if cw[0] == 1 else float(kw[0])
            out[(c, w0)] = ((kw - 1) / a if np.isfinite(kw) else np.nan, int(cw[0]))
    return out


def penetration(task):
    pen, driven = task
    t0 = time.time()
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return pd.DataFrame()
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    M = {"all": cell_means(mua, cc, ci, keep), "odd": cell_means(mua, cc, ci, keep & (idx % 2 == 1)),
         "even": cell_means(mua, cc, ci, keep & (idx % 2 == 0))}
    w_full = (t >= 40) & (t < 250)
    rows = []

    def emit(q, variant, direction, s, Lw, s_sel):
        for (c, w0), (L, cens) in Lw.items():
            rows.append(dict(penetration=pen, site_id=str(ds.channel.values[chans[q]]), variant=variant,
                             direction=direction, c=c, w0=w0, L=L, cens=cens, a=s["a"], onset=s["onset"],
                             contra_dom=s["contra_dom"], odi_sel=s_sel["odi"], contra_dom_sel=s_sel["contra_dom"]))

    for q in range(len(chans)):
        s_all = selection(M["all"], q, t, w_full)
        if s_all is not None:
            emit(q, "all", "all", s_all, windowed(M["all"], q, t, s_all["contra_dom"], s_all["a"], s_all["onset"]),
                 s_all)
        for X, Y in (("odd", "even"), ("even", "odd")):
            sX, sY = selection(M[X], q, t, w_full), selection(M[Y], q, t, w_full)
            if sX is None or sY is None:
                continue
            direction = f"{X}->{Y}"
            emit(q, "same", direction, sY, windowed(M[Y], q, t, sY["contra_dom"], sY["a"], sY["onset"]), sY)
            emit(q, "cv", direction, sX, windowed(M[Y], q, t, sX["contra_dom"], sX["a"], sX["onset"]), sX)
            sD = selection(M[Y], q, t, w_full, contra_dom=sX["contra_dom"])
            if sD is not None:
                emit(q, "cv_dom", direction, sD, windowed(M[Y], q, t, sD["contra_dom"], sD["a"], sD["onset"]), sX)
    print(f"{pen}: {len(chans)} channels ({time.time() - t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def early_late(C: pd.DataFrame, by=("penetration",)) -> pd.DataFrame:
    """Per group x direction: median over channels per window, then early and late means; directions averaged."""
    per = C.groupby(list(by) + ["direction", "c", "w0"]).L.median().unstack("w0")
    out = pd.DataFrame({"early": per[EARLY].mean(axis=1), "late": per[LATE].mean(axis=1)})
    return out.groupby(list(by) + ["c"]).mean()


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        C = pd.concat(list(ex.map(penetration, [(p, driven) for p in pens])), ignore_index=True)
    out = out_root() / "binocular"
    C.to_parquet(out / "v3_A2_channels.parquet")

    # ---- reproduction check against the saved Fig 2d input (k_windows.parquet)
    W = pd.read_parquet(out / "k_windows.parquet")
    per_ref = W.groupby(["penetration", "c", "w0"]).L.median().unstack("w0")
    ref = pd.DataFrame({"early": per_ref[EARLY].mean(axis=1), "late": per_ref[LATE].mean(axis=1)})

    rows = []
    PEN = {}
    for v in ("all", "same", "cv", "cv_dom"):
        PEN[v] = early_late(C[C.variant == v])
    PEN["fig2d_saved"] = ref
    for v, D in PEN.items():
        for c in DIOPTIC:
            g = D.xs(c, level="c").dropna()
            p = stats.wilcoxon(g.early, g.late).pvalue
            rows.append(dict(analysis="A2", variant=v, c=c, n_pen=len(g), early=g.early.median(), late=g.late.median(),
                             drop=(g.early - g.late).median(), p_early_vs_late=p, shift_late=np.nan, lo=np.nan,
                             hi=np.nan, passed=np.nan))
    # registered criteria (45%): cv drop significant; |median late cv - median late same| <= 0.25
    rng = np.random.default_rng(0)
    for c in DIOPTIC:
        cv, same = PEN["cv"].xs(c, level="c"), PEN["same"].xs(c, level="c")
        j = cv.index.intersection(same.index)
        sh = (cv.loc[j, "late"] - same.loc[j, "late"]).dropna()
        boot = [np.median(rng.choice(sh.values, len(sh))) for _ in range(5000)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        gap = cv.loc[j, "late"].median() - same.loc[j, "late"].median()
        cvp = cv.dropna()
        p_cv = stats.wilcoxon(cvp.early, cvp.late).pvalue if c == CRIT_C else np.nan
        rows.append(dict(analysis="A2 shift (cv - same), late", variant="cv vs same", c=c, n_pen=len(sh),
                         early=np.nan, late=gap, drop=np.nan, p_early_vs_late=stats.wilcoxon(sh).pvalue,
                         shift_late=sh.median(), lo=lo, hi=hi,
                         passed=((p_cv < 0.05) and (abs(gap) <= MARGIN)) if c == CRIT_C else np.nan))
    # dominance flips between halves, by |ODI| tercile (descriptive)
    ch = C[(C.variant == "cv") & (C.c == 0.45) & (C.w0 == 0)].copy()
    same_dom = C[(C.variant == "same") & (C.c == 0.45) & (C.w0 == 0)].set_index(["penetration", "site_id", "direction"])
    ch = ch.join(same_dom.contra_dom.rename("contra_dom_analysis"), on=["penetration", "site_id", "direction"])
    ch["flip"] = ch.contra_dom != ch.contra_dom_analysis
    S = pd.DataFrame(rows)
    S.to_csv(REPO / "docs" / "qc" / "v3_A2.csv", index=False)
    pd.concat({k: v for k, v in PEN.items()}, names=["variant"]).to_parquet(out / "v3_A2_penetrations.parquet")

    # ---- A3: |ODI| terciles from the selection half, index from the analysis half (cv)
    cvC = C[C.variant == "cv"].copy()
    cvC["abs_odi"] = cvC.odi_sel.abs()
    unit = cvC.drop_duplicates(["penetration", "site_id", "direction"])
    edges = np.quantile(unit.abs_odi, [1 / 3, 2 / 3])
    cvC["tercile"] = np.digitize(cvC.abs_odi, edges)            # 0 = most binocular, 2 = most monocular
    T3 = early_late(cvC, by=("penetration", "tercile"))
    rows3 = []
    for (terc, c), g in T3.groupby(level=["tercile", "c"]):
        g = g.dropna()
        if c not in (0.225, 0.45) or len(g) < 6:
            continue
        nch = unit[np.digitize(unit.abs_odi, edges) == terc]
        rows3.append(dict(tercile=int(terc), abs_odi_range=f"{[0, *edges, 1][terc]:.2f}-{[0, *edges, 1][terc + 1]:.2f}",
                          c=c, n_pen=len(g), n_channel_halves=len(nch), early=g.early.median(), late=g.late.median(),
                          drop=(g.early - g.late).median(), frac_drop=float((g.early > g.late).mean()),
                          p=stats.wilcoxon(g.early, g.late).pvalue))
    S3 = pd.DataFrame(rows3)
    S3.to_csv(REPO / "docs" / "qc" / "v3_A3.csv", index=False)
    T3.to_parquet(out / "v3_A3_penetrations.parquet")

    pd.set_option("display.width", 240)
    print("A2 (early, late = medians over penetrations; p = Wilcoxon early vs late; shift rows: late = gap of medians, "
          "shift_late = median paired shift, lo/hi = 95% CI)")
    print(S.round(3).to_string(index=False))
    flips = ch.groupby(pd.cut(ch.odi_sel.abs(), [0, *edges, 1], include_lowest=True)).flip.agg(["mean", "size"])
    print("dominance flips between halves by |ODI| (selection half):")
    print(flips.round(3).to_string())
    print(f"overall flip fraction {ch.flip.mean():.3f} of {len(ch)} channel-halves")
    print("A3 (cv index by |ODI| tercile):")
    print(S3.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
