"""Shared inputs for the paper figures: per-penetration population time courses of the full 6 x 6 matrix.

Usage: .venv\\Scripts\\python scripts\\figures\\figure_data.py
Per driven in-cortex MUA channel (driven in a monocular condition, run_crf.py; the rule used by Fig 2 and the model
tests): condition means over completed presentations without a saccade >= 1 deg, blank (0/0) subtracted, Gaussian
sigma 2 ms; axes reoriented so that axis 0 is the dominant eye (larger monocular response at .9, 40-250 ms); divided
by the channel's dominant-eye peak at .9 (40-250 ms). Penetration = mean over its channels, either on the flash
clock or with each channel aligned to its own half-max onset (dominant eye, .9) before averaging.
Writes <out_root>/figures/data/population.npz (flash [pen, 6, 6, t], onset [pen, 6, 6, rel]) and
population_meta.parquet (per penetration: monkey, channels, median eye weight a, onset, stimulus duration).
"""

import glob
import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS, eye_weight
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
T = np.arange(-50, 451)
REL = np.arange(-30, 301)


def penetration(pen: str, driven: set):
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return None
    mua = ds.mua.values[keep][:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy()[keep], P.contrast_ipsi.to_numpy()[keep]
    M = np.full((6, 6, t.size, len(chans)), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= 3:
                M[i, j] = mua[s].mean(axis=0)
    M = gaussian_filter1d(M - M[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    flash, onset, a_, on_ = [], [], [], []
    for q in range(len(chans)):
        m = M[..., q]
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        mo = m if dom_contra else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0 or np.isnan(mo[:, 0]).any() or np.isnan(mo[0, :]).any():
            continue
        mo = mo / pk
        idx = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] >= 0.5))
        if not idx.size:
            continue
        on = float(t[idx[0]])
        flash.append(np.stack([np.interp(T, t, x) for x in mo.reshape(36, -1)]).reshape(6, 6, -1))
        onset.append(np.stack([np.interp(on + REL, t, x) for x in mo.reshape(36, -1)]).reshape(6, 6, -1))
        a_.append(eye_weight(mo[:, 0][:, w].mean(axis=1), mo[0, :][:, w].mean(axis=1)))
        on_.append(on)
    if not flash:
        return None
    dur = float(np.median((P.pd_off_tick - P.pd_on_tick)[keep & (P.pd_off_tick > 0).to_numpy()]) / 30.0)
    meta = dict(penetration=pen, monkey=pen[7], n_channels=len(flash), a=float(np.median(a_)),
                onset=float(np.median(on_)), dur=dur)
    return np.nanmean(flash, axis=0), np.nanmean(onset, axis=0), meta


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    F, O, meta = [], [], []
    for pen in sorted({p for p, _ in driven}):
        r = penetration(pen, driven)
        if r is None:
            continue
        F.append(r[0])
        O.append(r[1])
        meta.append(r[2])
        print(pen, r[2]["n_channels"], flush=True)
    out = out_root() / "figures" / "data"
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "population.npz", flash=np.stack(F), onset=np.stack(O), t=T, rel=REL,
                        pens=np.array([m["penetration"] for m in meta]))
    pd.DataFrame(meta).to_parquet(out / "population_meta.parquet")
    print(out / "population.npz")


if __name__ == "__main__":
    main()
