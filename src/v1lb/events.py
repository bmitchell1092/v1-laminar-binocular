"""MonkeyLogic event codes recorded on the Blackrock digital port.

Codes arrive with bit 7 set (value = code + 128). Each ML trial is bracketed by
`9 9 9` ... `18 18 18`. bminteroc (T_binocular_gratings.m) sends:
116+block, 116+condition, 116+(trial mod 10) [not in 211008-211012]; 35, 11 fix spot
on / wait; 8 fixation acquired; 23/25/27 stimulus 1/2/3 on; 24/26 ISI scene
(= stimulus 1/2 off); 32, 36 end scene after stimulus 3 (there is no 28);
97 fixation broken; 96 reward.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TRIAL_START, TRIAL_END = 9, 18
STIM_ON = (23, 25, 27, 29, 31)
STIM_OFF = (24, 26, 28, 30, 32)
BREAK_FIXATION, REWARD = 97, 96
TRIAL_ID_RANGE = (116, 125)


def ml_codes(values: np.ndarray) -> np.ndarray:
    """Strip the strobe bit from raw NEV digital values."""
    values = np.asarray(values)
    if values.size and (values.min() < 128 or values.max() > 255):
        raise ValueError("digital values outside 128..255; strobe convention differs from this rig")
    return values.astype(np.int64) - 128


def _runs(codes: np.ndarray, value: int) -> list[tuple[int, int]]:
    """(first index, run length) of consecutive runs of `value`."""
    idx = np.flatnonzero(codes == value)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) != 1)
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    return list(zip(starts.tolist(), (ends - starts + 1).tolist()))


def split_trials(codes: np.ndarray, ticks: np.ndarray) -> pd.DataFrame:
    """One row per complete trial (a 9-triplet followed by an 18-triplet).

    Columns: nev_trial (0-based), i_start/i_end (code indices, inclusive),
    start_tick (first 9), end_tick (last 18), trial_digit (last 116-series code - 116,
    or -1 if absent), rewarded, broke_fixation.
    Incomplete brackets (recording started or stopped mid-trial, or a run of 9s/18s
    that is not exactly three long) are skipped; see `trial_bracket_problems`.
    """
    starts = [i for i, n in _runs(codes, TRIAL_START) if n == 3]
    ends = [i + n - 1 for i, n in _runs(codes, TRIAL_END) if n == 3]
    rows, e = [], 0
    for k, s in enumerate(starts):
        nxt = starts[k + 1] if k + 1 < len(starts) else len(codes)
        while e < len(ends) and ends[e] < s:
            e += 1
        if e == len(ends) or ends[e] >= nxt:
            continue
        seg = codes[s:ends[e] + 1]
        ids = seg[(seg >= TRIAL_ID_RANGE[0]) & (seg <= TRIAL_ID_RANGE[1])]
        rows.append(dict(i_start=s, i_end=ends[e], start_tick=int(ticks[s]), end_tick=int(ticks[ends[e]]),
                         trial_digit=int(ids[2] - 116) if ids.size >= 3 else -1,
                         rewarded=bool(np.any(seg == REWARD)),
                         broke_fixation=bool(np.any(seg == BREAK_FIXATION))))
    df = pd.DataFrame(rows, columns=["i_start", "i_end", "start_tick", "end_tick",
                                     "trial_digit", "rewarded", "broke_fixation"])
    df.insert(0, "nev_trial", np.arange(len(df)))
    return df


def trial_bracket_problems(codes: np.ndarray) -> list[str]:
    """Human-readable list of 9/18 runs that are not triplets."""
    out = []
    for value in (TRIAL_START, TRIAL_END):
        out += [f"run of {n} x {value} at code index {i}" for i, n in _runs(codes, value) if n != 3]
    return out


def presentations(codes: np.ndarray, ticks: np.ndarray, trials: pd.DataFrame) -> pd.DataFrame:
    """One row per stimulus onset code inside a trial.

    The offset is the first code in STIM_OFF after the onset and before the next
    onset. `completed` is False if there is no offset or fixation broke first.
    Times are event-code ticks; photodiode times are added elsewhere.
    """
    rows = []
    for t in trials.itertuples():
        seg = codes[t.i_start:t.i_end + 1]
        segt = ticks[t.i_start:t.i_end + 1]
        on_idx = np.flatnonzero(np.isin(seg, STIM_ON))
        for k, i in enumerate(on_idx):
            stop = on_idx[k + 1] if k + 1 < len(on_idx) else len(seg)
            off = next((j for j in range(i + 1, stop) if seg[j] in STIM_OFF), None)
            broke = next((j for j in range(i + 1, stop) if seg[j] == BREAK_FIXATION), None)
            completed = off is not None and (broke is None or broke > off)
            rows.append(dict(nev_trial=t.nev_trial, pres=(int(seg[i]) - 21) // 2, on_code=int(seg[i]),
                             on_tick=int(segt[i]),
                             off_code=int(seg[off]) if off is not None else -1,
                             off_tick=int(segt[off]) if off is not None else -1,
                             completed=completed))
    return pd.DataFrame(rows, columns=["nev_trial", "pres", "on_code", "on_tick",
                                       "off_code", "off_tick", "completed"])
