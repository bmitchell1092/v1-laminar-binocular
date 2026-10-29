"""Supplementary S1: do the core Fig 2/3 results survive the confound controls? (docs/09 register)

Usage: .venv\\Scripts\\python scripts\\fig2_robustness.py
Core statistics, recomputed per control subset of presentations (always: completed, no saccade >= 1 deg, no MUA artifact (05 C16)):
* claim 2: linearity L (window-specific eye weights) early (0-20 ms after each channel's onset) vs late (50-200);
  per penetration the channel median; Wilcoxon over penetrations, at dioptic .11, .225, .45;
* claim 3: slope of L on log contrast in the transition window (20-50 ms); fraction of penetrations negative.
Subsets: baseline; no saccade of any size in -200...500 ms (register #6); previous presentation absent or <= .11 in
both eyes (#8, adaptation); pre-stimulus pupil (-200...0 ms) inside its penetration's interquartile range (#23,
arousal); first and second half of each session (#23, drift).
Groupings on the baseline: per monkey (#18); without 211008 (#32); channels with eye weight a >= 0.3 (#26).
Writes docs/qc/fig2_robustness.csv.
"""

import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
DIOPTIC = [0.11, 0.225, 0.45]
EPOCHS = {"early": (0, 20), "transition": (20, 50), "late": (50, 200)}


def subsets(ds, P: pd.DataFrame, E: pd.DataFrame, pupil: np.ndarray) -> dict:
    base = clean_presentations(ds, P, E)
    prev = np.fmax(P.prev_contrast_L.to_numpy(float), P.prev_contrast_R.to_numpy(float))
    q1, q3 = np.nanpercentile(pupil[base], [25, 75])
    order = np.arange(len(P))
    return {"baseline": base,
            "no saccade at all": base & (E.n_saccades.to_numpy() == 0),
            "previous absent or <= .11": base & (~np.isfinite(prev) | (prev <= 0.11)),
            "pupil within IQR": base & (pupil >= q1) & (pupil <= q3),
            "first half of session": base & (order < len(P) / 2),
            "second half of session": base & (order >= len(P) / 2)}


def channel_rows(ds, mask, chans, pen, subset):
    t = ds.time.values
    P = pd.read_parquet(out_root() / "datasets" / f"{pen}_presentations.parquet")
    cc, ci = P.contrast_contra.to_numpy()[mask], P.contrast_ipsi.to_numpy()[mask]
    mua = ds.mua.values[mask][:, :, chans].astype(float)
    M = np.full((6, 6, t.size, len(chans)), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= 3:
                M[i, j] = mua[s].mean(axis=0)
    M = gaussian_filter1d(M - M[0, 0], 3, axis=2)
    w = (t >= 40) & (t < 250)
    rows = []
    for q in range(len(chans)):
        rc, ri = M[:, 0, :, q], M[0, :, :, q]
        if np.isnan(rc).any() or np.isnan(ri).any():
            continue
        rd, rn = (rc, ri) if rc[-1, w].mean() >= ri[-1, w].mean() else (ri, rc)
        pk = np.nanmax(rd[-1][(t >= 10) & (t < 250)])
        idx = np.flatnonzero((t >= 10) & (rd[-1] >= 0.5 * pk))
        if pk <= 0 or not idx.size:
            continue
        on = t[idx[0]]
        a_full = eye_weight(rd[:, w].mean(axis=1), rn[:, w].mean(axis=1))
        for c in DIOPTIC:
            i = list(CONTRASTS).index(c)
            for ep, (w0, w1) in EPOCHS.items():
                sel = (t >= on + w0) & (t < on + w1)
                d, n, b = rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1), M[i, i, sel, q].mean()
                a = eye_weight(d, n)
                k, cens = equivalent_contrast(d[:, None], np.array([b]), c)
                k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
                rows.append(dict(penetration=pen, subset=subset, channel=q, a_full=a_full, c=c, epoch=ep,
                                 L=(k - 1) / a if np.isfinite(k) else np.nan))
    return rows


def summarise(D: pd.DataFrame, label: str) -> list:
    out = []
    per = D.groupby(["penetration", "c", "epoch"]).L.median().unstack("epoch")
    for c in DIOPTIC:
        g = per.xs(c, level="c")[["early", "late"]].dropna()
        if len(g) >= 6:
            d = g.early - g.late
            out.append(dict(subset=label, stat=f"early-late drop, c={c}", n_pen=len(g), median=d.median(),
                            frac_pos=(d > 0).mean(), p=stats.wilcoxon(g.early, g.late).pvalue))
    tr = per["transition"].unstack("c").reindex(columns=DIOPTIC).dropna()
    if len(tr) >= 6:
        sl = tr.apply(lambda r: np.polyfit(np.log(DIOPTIC), r.values, 1)[0], axis=1)
        n_neg = int((sl < 0).sum())
        out.append(dict(subset=label, stat="transition slope on log contrast", n_pen=len(sl), median=sl.median(),
                        frac_pos=n_neg / len(sl), p=stats.binomtest(n_neg, len(sl), 0.5).pvalue))
    return out


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    rows = []
    for pen in sorted({p for p, _ in driven}):
        ds = open_penetration(pen)
        d = out_root() / "datasets"
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        t = ds.time.values
        pupil = ds.pupil.values[:, (t >= -200) & (t < 0)].astype(float).mean(axis=1)
        chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
        if not chans:
            continue
        for name, mask in subsets(ds, P, E, pupil).items():
            rows += channel_rows(ds, mask, chans, pen, name)
        print(pen, flush=True)
    D = pd.DataFrame(rows)
    D["monkey"] = D.penetration.str[7]
    out = []
    for name, g in D.groupby("subset", sort=False):
        out += summarise(g, name)
    b = D[D.subset == "baseline"]
    out += summarise(b[b.monkey == "B"], "baseline, monkey B")
    out += summarise(b[b.monkey == "J"], "baseline, monkey J")
    out += summarise(b[b.penetration != "211008_B_eA"], "baseline, without 211008")
    out += summarise(b[b.a_full >= 0.3], "baseline, channels with a >= 0.3")
    S = pd.DataFrame(out)
    pd.set_option("display.width", 200)
    print("(frac_pos: early-late rows = fraction of penetrations with a drop; slope rows = fraction negative)")
    print(S.round(4).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "fig2_robustness.csv", index=False)


if __name__ == "__main__":
    main()
