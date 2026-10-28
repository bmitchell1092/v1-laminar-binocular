import datetime as dt

import numpy as np
import pandas as pd
import pytest

from v1lb.trials import match_trials

ORIGIN = dt.datetime(2021, 10, 9, 20, 0, 0, tzinfo=dt.timezone.utc)


def make(n_ml=40, nev_missing_first=0, offset_s=-5 * 3600 + 2.2, seed=0):
    rng = np.random.default_rng(seed)
    end_s = np.cumsum(rng.uniform(1.5, 4.0, n_ml)) + 10           # trial ends, seconds after origin
    ml_local = [ORIGIN.replace(tzinfo=None) + dt.timedelta(seconds=s + offset_s + rng.normal(0, 0.02))
                for s in end_s]
    stim = pd.DataFrame({"trial": np.repeat(np.arange(1, n_ml + 1), 3),
                         "ml_time_local": np.repeat(ml_local, 3)})
    nev = pd.DataFrame({"nev_trial": np.arange(n_ml - nev_missing_first),
                        "end_tick": (end_s[nev_missing_first:] * 30000).astype(np.int64)})
    return nev, stim


def test_index_aligned():
    nev, stim = make()
    m = match_trials(nev, ORIGIN, stim)
    assert m.shift == 0 and m.n_matched == 40
    assert m.offset_s == pytest.approx(-5 * 3600 + 2.2, abs=0.1)
    assert m.table["ml_trial"].tolist() == list(range(1, 41))


def test_blackrock_missed_first_trial():
    nev, stim = make(nev_missing_first=1)
    m = match_trials(nev, ORIGIN, stim)
    assert m.shift == 1
    assert m.table["ml_trial"].iloc[0] == 2           # NEV trial 0 is ML trial 2
    assert m.max_dev_s < 0.2


def test_no_consistent_offset_raises():
    nev, stim = make()
    nev["end_tick"] = np.random.default_rng(1).permutation(nev["end_tick"].to_numpy())
    with pytest.raises(ValueError):
        match_trials(nev, ORIGIN, stim)
