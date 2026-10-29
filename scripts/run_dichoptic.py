"""Full 6 x 6 dichoptic response matrices per site (MUA channel or unit) and window.

Usage: .venv\\Scripts\\python scripts\\run_dichoptic.py [penetration ...]
Presentations: completed, no saccade >= 1 deg in -200...500 ms, no MUA artifact (05 C16). Responses: window means (MUA uV, units spikes/s)
minus the blank (0/0) mean, so cell (0, 0) is 0 by construction. Windows: transient 40-140, sustained 150-250,
full 40-250 ms after the photodiode flash.
Writes <out_root>/dichoptic/<pen>.parquet, long format: penetration, site, site_id, depth coordinates, unit metrics,
window, contrast_contra, contrast_ipsi, mean, sem, n.
"""

import sys
import time
import warnings

import numpy as np
import pandas as pd

from v1lb import crf
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
WINDOWS = dict(crf.WINDOWS, full=(40, 250))


def unit_rates(spikes, cluster_ids, onsets, window):
    out = np.zeros((len(onsets), len(cluster_ids)))
    for j, cid in enumerate(cluster_ids):
        tk = np.sort(spikes.loc[spikes.cluster_id == cid, "tick"].to_numpy())
        a = np.searchsorted(tk, onsets + window[0] * 30)
        b = np.searchsorted(tk, onsets + window[1] * 30)
        out[:, j] = (b - a) / ((window[1] - window[0]) / 1000)
    return out


def matrices(R: np.ndarray, Pk: pd.DataFrame, meta: pd.DataFrame, window: str) -> pd.DataFrame:
    blank = np.isclose(Pk.contrast_contra, 0) & np.isclose(Pk.contrast_ipsi, 0)
    R = R - np.nanmean(R[blank], axis=0)
    rows = []
    for (cc, ci), idx in Pk.groupby(["contrast_contra", "contrast_ipsi"]).groups.items():
        x = R[np.asarray(idx)]
        m, sd, n = np.nanmean(x, axis=0), np.nanstd(x, axis=0, ddof=1), np.sum(np.isfinite(x), axis=0)
        rows.append(meta.assign(window=window, contrast_contra=cc, contrast_ipsi=ci, mean=m, sem=sd / np.sqrt(n), n=n))
    return pd.concat(rows, ignore_index=True)


if __name__ == "__main__":
    d = out_root() / "datasets"
    out = out_root() / "dichoptic"
    out.mkdir(parents=True, exist_ok=True)
    pens = [a for a in sys.argv[1:] if not a.startswith("--")] or sorted(
        p.name.removesuffix("_presentation_eye.parquet") for p in d.glob("*_presentation_eye.parquet"))
    for pen in pens:
        t0 = time.time()
        ds = open_penetration(pen)
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = clean_presentations(ds, P, E)
        Pk = P[keep].reset_index(drop=True)
        t = ds.time.values
        incx = ds.in_cortex.values
        mua = ds.mua.values[keep][:, :, incx].astype(float)
        meta = pd.DataFrame(dict(penetration=pen, site="mua", site_id=ds.channel.values[incx].astype(str),
                                 contact=ds.contact.values[incx], depth_below_toc=ds.depth_below_toc.values[incx],
                                 depth_l45=ds.depth_l45.values[incx]))
        parts = [matrices(crf.window_mean(mua, t, w), Pk, meta, name) for name, w in WINDOWS.items()]
        upath = out_root() / "units" / f"{pen}_units.parquet"
        if upath.exists():
            U = pd.read_parquet(upath)
            S = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
            umeta = pd.DataFrame(dict(penetration=pen, site="unit", site_id=U.cluster_id.astype(str).to_numpy(),
                                      contact=U.peak_contact.fillna(U.contact).to_numpy(),
                                      depth_below_toc=U.depth_below_toc.to_numpy(), depth_l45=U.depth_l45.to_numpy(),
                                      well_isolated=U.well_isolated.to_numpy(),
                                      isi_contamination=U.isi_contamination.to_numpy(), snr=U.snr.to_numpy(),
                                      trough_to_peak_ms=U.trough_to_peak_ms.to_numpy(), rate_hz=U.rate_hz.to_numpy()))
            parts += [matrices(unit_rates(S, U.cluster_id.to_numpy(), Pk.pd_on_tick.to_numpy(), w), Pk, umeta, name)
                      for name, w in WINDOWS.items()]
        res = pd.concat(parts, ignore_index=True)
        res.to_parquet(out / f"{pen}.parquet")
        print(f"{pen}: {time.time() - t0:.0f}s, {res.site_id.nunique()} sites, cells n median "
              f"{int(res.n.median())}", flush=True)
