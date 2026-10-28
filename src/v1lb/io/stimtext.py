"""MonkeyLogic 2 stimulus text files (`.g<PARADIGM>Grating_di`, tab-delimited, header row).

T_binocular_gratings.m appends one row per *planned* presentation (3 per
trial) at the end of every trial, including aborted trials, with
TrialRecord.CurrentTrialNumber and `now` (MATLAB datenum, local PC time).
The column layout changed between sessions (e.g. 211008 has gaborfilter_on /
gabor_std, later files grating_phzdist / grating_posdist), so columns are taken
from the header rather than from a fixed list.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd


def datenum_to_datetime(x: float) -> dt.datetime:
    """MATLAB datenum -> naive datetime (same local clock it was written in)."""
    day = int(x)
    return dt.datetime.fromordinal(day - 366) + dt.timedelta(days=x - day)


def read_stim_text(path: str | Path) -> pd.DataFrame:
    """One row per planned presentation, with `trial` (ML trial number), `pres`
    (1-based position in the trial) and `ml_time_local` (datetime) added."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    if df.columns[0] != "trial":
        raise ValueError(f"{path.name}: expected a header row starting with 'trial'")
    df["trial"] = df["trial"].astype(int)
    if (np.diff(df["trial"].to_numpy()) < 0).any():
        raise ValueError(f"{path.name}: trial numbers are not monotonic")
    df["pres"] = df.groupby("trial").cumcount() + 1
    if "timestamp" in df.columns:
        df["ml_time_local"] = [datenum_to_datetime(v) for v in df["timestamp"].to_numpy(float)]
    return df
