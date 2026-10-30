"""Figure 4 inputs: population CSD aligned to each penetration's earliest current sink, and the sink current and
MUA at the sink channels for every contrast cell.

Usage: .venv\\Scripts\\python scripts\\figures\\sink_data.py
Sinks as chosen by dichoptic_csd.py (docs/qc/dichoptic_csd_sinks.csv, usable sinks only); presentations as everywhere
(completed, no saccade >= 1 deg, no MUA artifact). CSD from the dataset (bad contacts interpolated, 3-point Hamming
smoothing across depth, 05 C12), baseline -50...0 ms removed.
* map: mean CSD over presentations with either eye >= .45, channels -8...+12 contacts from the sink channel
  (+ = deeper), divided by the magnitude of the sink's own peak (20-70 ms), so each penetration counts equally;
* cells: sink current (-CSD over the sink channel and its two neighbours) and MUA over the same channels, condition
  means (blank-subtracted, Gaussian sigma 2 ms) on the sink-onset clock, axis 0 = dominant eye (larger response at .9
  in the first 50 ms after sink onset), each scaled by its dominant-eye peak at .9 in that interval.
Writes <out_root>/figures/data/sink.npz.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

REPO = Path(__file__).resolve().parents[2]
OFFSETS = np.arange(-8, 13)
T_MAP = np.arange(-10, 151)
REL = np.arange(-10, 201)


def cell_means(x, cc, ci):
    M = np.full((6, 6, x.shape[1]), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= 3:
                M[i, j] = x[s].mean(axis=0)
    return M - M[0, 0]


def main():
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    d = out_root() / "datasets"
    maps, csd_cells, mua_cells, pens = [], [], [], []
    for _, s in S.iterrows():
        pen = s.penetration
        ds = open_penetration(pen)
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = clean_presentations(ds, P, E)
        t = ds.time.values
        cc, ci = P.contrast_contra.to_numpy()[keep], P.contrast_ipsi.to_numpy()[keep]
        csd = ds.csd.values[keep].astype(float)
        csd -= np.nanmean(csd[:, (t >= -50) & (t < 0)], axis=1, keepdims=True)
        k = list(ds.channel.values.astype(str)).index(s.sink_channel)
        m = np.nanmean(csd[np.maximum(cc, ci) >= 0.45], axis=0)
        peak = abs(np.nanmin(m[(t >= 20) & (t < 70), k]))
        rows = []
        for o in OFFSETS:
            q = k + o
            ok = 0 <= q < m.shape[1] and ds.in_cortex.values[q]
            rows.append(np.interp(T_MAP, t, m[:, q]) / peak if ok else np.full(T_MAP.size, np.nan))
        maps.append(np.array(rows))
        chans = [q for q in (k - 1, k, k + 1) if 0 <= q < m.shape[1] and ds.in_cortex.values[q]]
        for sig, store in (("csd", csd_cells), ("mua", mua_cells)):
            x = -np.nanmean(csd[:, :, chans], axis=2) if sig == "csd" else ds.mua.values[keep][:, :, chans].astype(float).mean(axis=2)
            M = cell_means(gaussian_filter1d(x, 2, axis=1), cc, ci)
            Mi = np.stack([[np.interp(s.sink_onset_ms + REL, t, M[i, j]) for j in range(6)] for i in range(6)])
            w50 = (REL >= 0) & (REL < 50)
            Mi = Mi if np.nanmean(Mi[-1, 0, w50]) >= np.nanmean(Mi[0, -1, w50]) else Mi.transpose(1, 0, 2)
            store.append(Mi / np.nanmax(Mi[-1, 0, w50]))
        pens.append(pen)
        print(pen, flush=True)
    out = out_root() / "figures" / "data"
    np.savez_compressed(out / "sink.npz", maps=np.array(maps), csd=np.array(csd_cells), mua=np.array(mua_cells),
                        offsets=OFFSETS, t_map=T_MAP, rel=REL, pens=np.array(pens))
    print(out / "sink.npz", len(pens))


if __name__ == "__main__":
    main()
