"""Flash-evoked LFP/CSD/MUA from the evp recordings (full-field flashes; one per ML trial).

evp event codes: 9 9 9, 23 (flash on), 24 (flash off), 18 18 18. Flash timing comes from the
photodiode exactly as for bminteroc (`photodiode.detect`). The landmark compared with the
laminar alignment is the initial (earliest) sink: its onset, and the contacts it spans in its
first 15 ms. Its bottom edge is the classic "sink bottom" estimate of the L4/L5 boundary.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import photodiode
from .config import recording_files
from .dataset import lfp_epochs
from .events import ml_events, presentations, split_trials
from .io import NevFile, NsxFile
from .laminar import csd as compute_csd
from .laminar import earliest_sink
from .signals import analog_mua_epochs, probe_labels


@dataclass
class FlashCSD:
    recording: str
    n_flashes: int
    t_ms: np.ndarray
    contacts: np.ndarray        # contact numbers, superficial -> deep
    lfp: np.ndarray             # [time, channel], baseline-corrected mean, uV
    csd: np.ndarray             # [time, channel], nA/mm^3, sinks negative
    mua: np.ndarray             # [time, channel], baseline-corrected analog MUA, uV
    sink_index: int             # channel index (superficial -> deep) of the earliest qualifying sink, -1 if none
    sink_onset_ms: float
    sink_top: int               # channel index range of the initial sink (inclusive), -1 if none
    sink_bottom: int


def initial_sink_extent(csd: np.ndarray, t_ms: np.ndarray, onset_ms: float, seed: int,
                        span_ms: float = 15.0, rel: float = 0.5) -> tuple[int, int]:
    """Contiguous channels around `seed` whose mean CSD over [onset, onset + span_ms] is below
    `rel` x the seed channel's value (i.e. part of the same sink). Returns (top, bottom) indices."""
    w = (t_ms >= onset_ms) & (t_ms <= onset_ms + span_ms)
    m = np.nanmean(csd[w], axis=0)
    core = m[seed]
    if not np.isfinite(core) or core >= 0:
        return seed, seed
    top = bottom = seed
    while top - 1 >= 0 and np.isfinite(m[top - 1]) and m[top - 1] <= rel * core:
        top -= 1
    while bottom + 1 < m.size and np.isfinite(m[bottom + 1]) and m[bottom + 1] <= rel * core:
        bottom += 1
    return top, bottom


def flash_csd(session: str, recording: str, bank: str, orientation: str, window_ms=(-100, 300),
              sink_window_ms=(20, 90)) -> FlashCSD:
    files = recording_files(session, recording)
    nev = NevFile(files[".nev"])
    ticks, values, _ = nev.digital_events()
    ticks, codes, _ = ml_events(ticks, values)
    pres = presentations(codes, ticks, split_trials(codes, ticks))
    ns2, ns6 = NsxFile(files[".ns2"]), NsxFile(files[".ns6"])
    pres, _ = photodiode.detect(ns6, pres)
    onsets = pres.loc[pres.pd_on_tick > 0, "pd_on_tick"].to_numpy()

    labels = probe_labels(ns2, bank)
    contacts = np.arange(1, len(labels) + 1) if orientation == "ascending" else np.arange(len(labels), 0, -1)
    ordered = [labels[c - 1] for c in contacts]
    t_ms = np.arange(*window_ms).astype(float)
    base = (t_ms >= -50) & (t_ms < 0)
    lfp = lfp_epochs(ns2, ordered, onsets, window_ms)
    lfp = (lfp - lfp[:, base].mean(axis=1, keepdims=True)).mean(axis=0)
    csd = compute_csd(lfp)
    mua = analog_mua_epochs(ns6, ordered, onsets, pre_ms=-window_ms[0], post_ms=window_ms[1])
    mua = (mua - mua[:, base].mean(axis=1, keepdims=True)).mean(axis=0)
    onset, first = earliest_sink(csd, t_ms, window_ms=sink_window_ms)
    if first >= 0:
        top, bottom = initial_sink_extent(csd, t_ms, float(onset[first]), first)
        on = float(onset[first])
    else:
        top = bottom = -1
        on = float("nan")
    return FlashCSD(recording, len(onsets), t_ms, contacts, lfp, csd, mua, first, on, top, bottom)
