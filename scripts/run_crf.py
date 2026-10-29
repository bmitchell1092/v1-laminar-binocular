"""CRFs and supersaturation for every built dataset: MUA per in-cortex channel and every imported unit.

Usage: .venv\\Scripts\\python scripts\\run_crf.py [penetration ...]
Presentations: completed, no saccade >= 1 deg in -200...500 ms (run_saccades.py). Per site (channel or unit),
condition (mono_contra, mono_ipsi, dioptic) and window (transient 40-140, sustained 150-250 ms): blank-subtracted
means per contrast, driven (one-way ANOVA across the 6 contrasts p < 0.01 and peak > 0), SI with a 500-sample
bootstrap CI, the cross-validated SI (unbiased for monotonic CRFs; the
population statistic), the held-out downturn test (p_downturn; calibrated: <= 5% of
monotonic sites reach p < 0.05), and the Naka-Rushton fit with supersaturation exponent s (F-test vs s = 1).
Writes <out_root>/crf/<pen>_crf.parquet.
"""

import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from v1lb import crf
from v1lb.config import out_root
from v1lb.dataset import open_penetration

warnings.filterwarnings("ignore")


def unit_rates(spikes: pd.DataFrame, cluster_ids, onsets: np.ndarray, window) -> np.ndarray:
    """Spikes/s per presentation in [onset + window) for each unit -> [presentation, unit]."""
    out = np.zeros((len(onsets), len(cluster_ids)))
    for j, cid in enumerate(cluster_ids):
        tk = np.sort(spikes.loc[spikes.cluster_id == cid, "tick"].to_numpy())
        a = np.searchsorted(tk, onsets + window[0] * 30)
        b = np.searchsorted(tk, onsets + window[1] * 30)
        out[:, j] = (b - a) / ((window[1] - window[0]) / 1000)
    return out


def analyse(R: np.ndarray, Pk: pd.DataFrame, meta: pd.DataFrame, window: str) -> list[dict]:
    blank = np.isclose(Pk.contrast_contra, 0) & np.isclose(Pk.contrast_ipsi, 0)
    R = R - np.nanmean(R[blank], axis=0)
    rows = []
    for cond in crf.CONDITIONS:
        c, r = crf.crf_points(R, Pk, cond)
        m = crf.means_by_contrast(c, r)
        si = crf.supersaturation_index(m)
        boot = crf.bootstrap_si(c, r, 500)
        cv = crf.cv_supersaturation_index(c, r, 200)
        p_down = crf.downturn_p(c, r, 50)
        n_per = [int(np.isclose(c, k).sum()) for k in crf.CONTRASTS]
        for k in range(R.shape[1]):
            groups = [r[np.isclose(c, cc), k] for cc in crf.CONTRASTS]
            p_anova = stats.f_oneway(*groups).pvalue if min(len(g) for g in groups) > 1 else np.nan
            driven = bool(p_anova < 0.01 and np.nanmax(m[:, k]) > 0)
            row = dict(meta.iloc[k]) | dict(window=window, condition=cond, driven=driven, p_anova=p_anova,
                                            n_min=min(n_per), si=si[k], si_cv=cv[k], p_downturn=p_down[k], si_lo=np.nanpercentile(boot[:, k], 2.5),
                                            si_hi=np.nanpercentile(boot[:, k], 97.5),
                                            peak_contrast_obs=crf.CONTRASTS[int(np.nanargmax(m[:, k]))])
            row |= {f"r{cc:g}": m[i, k] for i, cc in enumerate(crf.CONTRASTS)}
            if driven:
                f = crf.fit_crf(c, r[:, k])
                row |= dict(rmax=f.rmax, c50=f.c50, n=f.n, s=f.s, s_p=f.p, peak_contrast_fit=f.peak_contrast)
            rows.append(row)
    return rows


if __name__ == "__main__":
    d = out_root() / "datasets"
    out = out_root() / "crf"
    out.mkdir(parents=True, exist_ok=True)
    pens = [a for a in sys.argv[1:] if not a.startswith("--")] or sorted(
        p.name.removesuffix("_presentation_eye.parquet") for p in d.glob("*_presentation_eye.parquet"))
    for pen in pens:
        t0 = time.time()
        ds = open_penetration(pen)
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
        Pk = P[keep].reset_index(drop=True)
        t = ds.time.values
        incx = ds.in_cortex.values
        mua = ds.mua.values[keep][:, :, incx].astype(float)
        meta = pd.DataFrame(dict(penetration=pen, site="mua", site_id=ds.channel.values[incx].astype(str),
                                 contact=ds.contact.values[incx], depth_below_toc=ds.depth_below_toc.values[incx],
                                 depth_l45=ds.depth_l45.values[incx]))
        rows = []
        for wname, w in crf.WINDOWS.items():
            rows += analyse(crf.window_mean(mua, t, w), Pk, meta, wname)
        upath = out_root() / "units" / f"{pen}_units.parquet"
        n_units = 0
        if upath.exists():
            U = pd.read_parquet(upath)
            S = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
            umeta = pd.DataFrame(dict(penetration=pen, site="unit", site_id=U.cluster_id.astype(str).to_numpy(),
                                      contact=U.peak_contact.fillna(U.contact).to_numpy(),
                                      depth_below_toc=U.depth_below_toc.to_numpy(), depth_l45=U.depth_l45.to_numpy(),
                                      group=U.group.to_numpy(), well_isolated=U.well_isolated.to_numpy(),
                                      isi_contamination=U.isi_contamination.to_numpy(),
                                      amplitude_cutoff=U.amplitude_cutoff.to_numpy(), snr=U.snr.to_numpy(),
                                      trough_to_peak_ms=U.trough_to_peak_ms.to_numpy()))
            for wname, w in crf.WINDOWS.items():
                rows += analyse(unit_rates(S, U.cluster_id.to_numpy(), Pk.pd_on_tick.to_numpy(), w), Pk, umeta, wname)
            n_units = len(U)
        res = pd.DataFrame(rows)
        res.to_parquet(out / f"{pen}_crf.parquet")
        dm = res[(res.site == "mua") & res.driven]
        print(f"{pen}: {time.time() - t0:.0f}s, kept {keep.sum()}/{len(P)}; MUA sites {incx.sum()}, units {n_units}; "
              f"driven MUA rows {len(dm)}, sig SI>0 {int((dm.si_lo > 0).sum())}", flush=True)
