"""V4 check J, sensitivity (descriptive, not registered; docs/09): MUA - sink current timing at the earliest sink
under different processing.

Usage: .venv\\Scripts\\python scripts\\v4_csd_sensitivity.py
CSD from condition-mean LFP (bad contacts interpolated, blank-subtracted) with or without 3-point Hamming smoothing
across contacts, at the sink contact alone or averaged with its two neighbours; MUA at the sink contact alone or the
same three; with or without the Gaussian sigma 2 ms; the current advanced by the measured 0.9 ms hardware-filter
delay (qc_nsx_clock.py). Conditions: both eyes 45%, both eyes 90%, and the single eye with the larger sink at 90%.
Half-max onset (10-150 ms after the flash, interpolated) and peak time (10-100 ms), per penetration, MUA - current.
Writes docs/qc/v4_J_sensitivity.csv (per penetration and variant); prints the medians.
"""

import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from v4_latency import onset

from v1lb import laminar as lam
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
HW_DELAY = 0.9


def peak_time(x, t, lo=10, hi=100):
    w = (t >= lo) & (t < hi)
    return float(t[w][np.nanargmax(x[w])])


def work(task):
    pen, ch_label = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values.astype(float)
    ch = list(ds.channel.values.astype(str)).index(ch_label)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    blank = keep & (cc == 0) & (ci == 0)
    conds = {"both 45%": keep & np.isclose(cc, .45) & np.isclose(ci, .45),
             "both 90%": keep & np.isclose(cc, .9) & np.isclose(ci, .9),
             "contra 90%": keep & np.isclose(cc, .9) & (ci == 0), "ipsi 90%": keep & (cc == 0) & np.isclose(ci, .9)}
    lfp = ds.lfp.values.astype(float)
    mua = ds.mua.values.astype(float)
    bad = ds.bad_lfp.values
    L0, M0 = lfp[blank].mean(axis=0), mua[blank].mean(axis=0)
    rows = []
    for name, sel in conds.items():
        L, M = lfp[sel].mean(axis=0) - L0, mua[sel].mean(axis=0) - M0
        for ham in (False, True):
            C = -lam.csd(L, smooth=ham, bad=bad)                          # sink positive
            C -= np.nanmean(C[(t >= -50) & (t < 0)], axis=0)
            nb = [k for k in (ch - 1, ch, ch + 1) if 0 <= k < C.shape[1] and np.isfinite(C[:, k]).all()]
            for span, cs in (("1 contact", [ch]), ("3 contacts", nb)):
                cur, mu = C[:, cs].mean(axis=1), M[:, cs].mean(axis=1)
                for sm in (0, 2):
                    c_, m_ = (gaussian_filter1d(cur, sm), gaussian_filter1d(mu, sm)) if sm else (cur, mu)
                    pc, pm = np.nanmax(c_[(t >= 10) & (t < 150)]), np.nanmax(m_[(t >= 10) & (t < 150)])
                    oc = onset(c_ / pc, t, lo=10, hi=150, floor=0) - HW_DELAY
                    om = onset(m_ / pm, t, lo=10, hi=150, floor=0)
                    rows.append(dict(penetration=pen, condition=name, hamming=ham, span=span, sigma=sm,
                                     d_onset=om - oc, d_peak=peak_time(m_, t) - (peak_time(c_, t) - HW_DELAY),
                                     sink_peak=pc))
    return rows


def main():
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    with ProcessPoolExecutor(6) as ex:
        R = pd.DataFrame([r for res in ex.map(work, list(zip(S.penetration, S.sink_channel))) for r in res])
    mono = R[R.condition.str.endswith("90%") & ~R.condition.str.startswith("both")]
    idx = mono.groupby(["penetration", "hamming", "span", "sigma"]).sink_peak.idxmax()
    R = pd.concat([R[R.condition.str.startswith("both")], mono.loc[idx].assign(condition="dominant eye 90%")])
    R.to_csv(REPO / "docs" / "qc" / "v4_J_sensitivity.csv", index=False)
    pd.set_option("display.width", 200)
    out = R.groupby(["condition", "hamming", "span", "sigma"]).agg(
        n=("d_onset", "count"), onset_median=("d_onset", "median"),
        mua_earlier=("d_onset", lambda v: int((v < 0).sum())), peak_median=("d_peak", "median"),
        mua_peak_earlier=("d_peak", lambda v: int((v < 0).sum())))
    print("MUA - current (ms; current corrected for the 0.9 ms hardware delay):")
    print(out.round(2).to_string())


if __name__ == "__main__":
    main()
