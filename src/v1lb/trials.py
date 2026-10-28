"""Match Blackrock (NEV) trials to MonkeyLogic stimulus-file trials by time.

The stimulus file stamps each trial's rows with `now` (local PC time) at the end
of the trial; NEV trials end at the last `18` (NSP clock, origin in UTC). Their
difference is constant within a session (time zone + PC clock offset, e.g.
-5 h - 2.2 s), so the correct pairing is the trial shift that makes it constant.
Index-based matching (the legacy method) is wrong whenever Blackrock missed a
trial, as in 211009_B_bminteroc001.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TrialMatch:
    shift: int            # stimulus-file trial row k + shift  <->  NEV trial k (row order)
    offset_s: float       # median (ML local time - NSP UTC time), seconds
    max_dev_s: float      # max |residual - median| over matched trials
    n_matched: int
    table: pd.DataFrame   # nev_trial, ml_trial, residual_s


def match_trials(nev_trials: pd.DataFrame, origin_utc: dt.datetime, stim: pd.DataFrame,
                 tick_rate: float = 30000.0, max_shift: int = 5, tolerance_s: float = 1.0) -> TrialMatch:
    """Pair NEV trials (from `events.split_trials`) with stimulus-file trials (`read_stim_text`).

    Raises if no shift gives residuals within `tolerance_s` of a constant.
    """
    ml = stim.groupby("trial", sort=True)["ml_time_local"].first()
    ml_trials = ml.index.to_numpy()
    epoch = np.datetime64("1970-01-01T00:00:00", "us")
    ml_s = (ml.to_numpy(dtype="datetime64[us]") - epoch) / np.timedelta64(1, "s")
    origin_s = (np.datetime64(origin_utc.replace(tzinfo=None), "us") - epoch) / np.timedelta64(1, "s")
    nev_s = origin_s + nev_trials["end_tick"].to_numpy() / tick_rate

    best = None
    for shift in range(-max_shift, max_shift + 1):
        k = np.arange(len(nev_s))
        j = k + shift
        ok = (j >= 0) & (j < len(ml_s))
        if ok.sum() < max(3, 0.5 * len(nev_s)):
            continue
        res = ml_s[j[ok]] - nev_s[k[ok]]
        med = float(np.median(res))
        dev = float(np.max(np.abs(res - med)))
        if best is None or dev < best[2]:
            best = (shift, med, dev, k[ok], j[ok], res)
    if best is None or best[2] > tolerance_s:
        raise ValueError(f"no trial shift within ±{max_shift} gives a constant ML-NEV time offset "
                         f"(best max deviation {None if best is None else round(best[2], 3)} s)")
    shift, med, dev, k, j, res = best
    table = pd.DataFrame({"nev_trial": nev_trials["nev_trial"].to_numpy()[k],
                          "ml_trial": ml_trials[j], "residual_s": res - med})
    return TrialMatch(shift, med, dev, len(k), table)
