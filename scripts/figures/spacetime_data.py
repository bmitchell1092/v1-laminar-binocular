"""Shared inputs for the space-time figures: per-channel condition time courses with depth relative to the earliest
sink, on the stimulus clock and on each channel's response-onset clock.

Usage: .venv\\Scripts\\python scripts\\figures\\spacetime_data.py [n_workers]
Per driven in-cortex MUA channel (the figure_data.py rule) in penetrations with a usable earliest sink: condition
means over clean presentations, blank subtracted, Gaussian sigma 2 ms, axes reoriented so that the dominant eye comes
first, divided by the channel's dominant-eye peak at 90% (40-250 ms). Conditions kept: dominant eye alone at all six
contrasts (rows 0-5), other eye alone at 11, 22.5, 45, 90% (rows 6-9) and both eyes at the same four (rows 10-13).
Writes <out_root>/figures/data/spacetime.npz (stim [channel, 14, t], onset [channel, 14, rel], t, rel) and
spacetime_meta.parquet (penetration, monkey, channel, half-max onset, eye weight a, depth, depth relative to the sink
in contacts).
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS, eye_weight
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
T = np.arange(-20, 251)
REL = np.arange(-30, 251)
CS = (0.11, 0.225, 0.45, 0.9)


def penetration(task):
    pen, driven, sink = task
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
    ii = [list(CONTRASTS).index(c) for c in CS]
    out = dict(stim=[], onset=[], meta=[])
    for q, k in enumerate(chans):
        m = M[..., q]
        mo = m if m[-1, 0][w].mean() >= m[0, -1][w].mean() else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0 or np.isnan(mo[:, 0]).any() or np.isnan(mo[0, :]).any():
            continue
        mo = mo / pk
        idx = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] >= 0.5))
        if not idx.size:
            continue
        on = float(t[idx[0]])
        X = np.concatenate([mo[:, 0], mo[0, ii], mo[ii, ii]])
        out["stim"].append(np.stack([np.interp(T, t, x) for x in X]))
        out["onset"].append(np.stack([np.interp(on + REL, t, x) for x in X]))
        depth = float(ds.depth_below_toc.values[k])
        out["meta"].append(dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), onset=on,
                                a=eye_weight(mo[:, 0][:, w].mean(axis=1), mo[0, :][:, w].mean(axis=1)),
                                depth=depth, rel_sink=depth - sink))
    return out if out["meta"] else None


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    pens = sorted(p for p in {p for p, _ in driven} if p in S.index and S.usable[p])
    with ProcessPoolExecutor(n_workers) as ex:
        done = [r for r in ex.map(penetration, [(p, driven, float(S.sink_depth_below_toc[p])) for p in pens]) if r]
    out = out_root() / "figures" / "data"
    np.savez_compressed(out / "spacetime.npz", stim=np.concatenate([np.stack(r["stim"]) for r in done]).astype(np.float32),
                        onset=np.concatenate([np.stack(r["onset"]) for r in done]).astype(np.float32), t=T, rel=REL)
    meta = pd.DataFrame([m for r in done for m in r["meta"]])
    meta.to_parquet(out / "spacetime_meta.parquet")
    print(out / "spacetime.npz", len(meta), "channels,", meta.penetration.nunique(), "penetrations")


if __name__ == "__main__":
    main()
