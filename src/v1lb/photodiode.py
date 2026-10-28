"""Photodiode timing of stimulus presentations (analog input `ainp1` in the .ns6).

What the signal looks like (measured, docs/05 A2-A2c):
* The ViewPixx backlight is strobed: while the photodiode box is white, the
  diode sees one ~3.5 ms flash per 8.33 ms frame. Flashes drive `ainp1` from
  ~+980 mV (dark) to ~-334 mV. The gray background (fixation scene) dips only
  to ~830-910 mV, and the black ISI box stays near dark.
* The first flash follows the stimulus-on code by ~20.8 ms in most files, by
  ~12.5 ms (one frame less) in some, and switches between the two within a
  few files. Event codes alone are therefore not a reliable stimulus time.

Definitions used here (all NSP ticks):
* pd_on: first flash rising through the half-way level after the on code.
* pd_n_frames: number of consecutive flashes (one per frame) from pd_on.
* pd_off: one frame after the last flash, i.e. when the next flash would
  have come. pd_off - pd_on = pd_n_frames x frame period.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .io.nsx import NsxFile

TICKS_PER_MS = 30


@dataclass(frozen=True)
class PhotodiodeLevels:
    dark: float        # signal with the box black / no flash (physical units, mV)
    white: float       # flash extreme with the box white
    polarity: int      # +1 if flashes raise the signal, -1 if they lower it
    frame_ticks: float # flash-to-flash interval

    @property
    def threshold(self) -> float:
        return 0.5 * (self.dark + self.white)


def flash_onsets(x: np.ndarray, ticks: np.ndarray, levels: PhotodiodeLevels) -> np.ndarray:
    """NSP ticks where the signal crosses the half-way level in the flash direction."""
    s = levels.polarity * (x - levels.threshold) > 0
    idx = np.flatnonzero(s[1:] & ~s[:-1]) + 1
    return ticks[idx]


def estimate_levels(nsx: NsxFile, on_ticks: np.ndarray, channel: str = "ainp1",
                    n: int = 200) -> PhotodiodeLevels:
    """Session-level dark/white levels and frame period from a sample of presentations."""
    ms = TICKS_PER_MS
    sample = np.asarray(on_ticks)[np.linspace(0, len(on_ticks) - 1, min(n, len(on_ticks))).astype(int)]
    dark, white = [], []
    for on in sample:
        pre, _ = nsx.read(on - 60 * ms, on - 2 * ms, [channel])
        post, _ = nsx.read(on + 5 * ms, on + 60 * ms, [channel])
        d = np.nanmedian(pre)
        dark.append(d)
        dev = post[:, 0] - d
        white.append(post[np.nanargmax(np.abs(dev)), 0])
    d, w = float(np.median(dark)), float(np.median(white))
    levels = PhotodiodeLevels(d, w, 1 if w > d else -1, 250.0)
    intervals = []
    for on in sample[:50]:
        x, t = nsx.read(on, on + 200 * ms, [channel])
        f = flash_onsets(x[:, 0], t, levels)
        intervals.extend(np.diff(f)[np.diff(f) < 1.5 * 250])
    frame = float(np.median(intervals)) if intervals else 250.0
    return PhotodiodeLevels(d, w, levels.polarity, frame)


def detect(nsx: NsxFile, pres: pd.DataFrame, levels: PhotodiodeLevels | None = None,
           channel: str = "ainp1", first_flash_ms: tuple[float, float] = (5, 60),
           max_stim_ms: float = 500, tail_ms: float = 60) -> tuple[pd.DataFrame, PhotodiodeLevels]:
    """Add photodiode columns to a presentations table (from `events.presentations`).

    Adds pd_on_tick, pd_off_tick, pd_n_frames, pd_on_lat_ms (pd_on - on code),
    pd_off_lat_ms (pd_off - off code) and pd_flag ('' = ok, 'no_flash', 'gap').
    """
    ms = TICKS_PER_MS
    if levels is None:
        levels = estimate_levels(nsx, pres["on_tick"].to_numpy(), channel)
    on = pres["on_tick"].to_numpy()
    off_code = pres["off_tick"].to_numpy()
    order = np.argsort(on)
    next_on = np.full(len(on), np.iinfo(np.int64).max)
    next_on[order[:-1]] = on[order[1:]]
    frame = levels.frame_ticks
    out = {k: [] for k in ("pd_on_tick", "pd_off_tick", "pd_n_frames", "pd_flag")}
    for i in range(len(pres)):
        end = (off_code[i] if off_code[i] > 0 else on[i] + max_stim_ms * ms) + tail_ms * ms
        end = min(end, next_on[i] - 2 * ms)
        x, t = nsx.read(on[i] + first_flash_ms[0] * ms, end, [channel])
        f = flash_onsets(x[:, 0], t, levels)
        if f.size == 0 or f[0] > on[i] + first_flash_ms[1] * ms:
            out["pd_on_tick"].append(-1); out["pd_off_tick"].append(-1)
            out["pd_n_frames"].append(0); out["pd_flag"].append("no_flash")
            continue
        gaps = np.flatnonzero(np.diff(f) > 1.5 * frame)
        run = f[: gaps[0] + 1] if gaps.size else f
        out["pd_on_tick"].append(int(run[0]))
        out["pd_off_tick"].append(int(round(run[-1] + frame)))
        out["pd_n_frames"].append(int(run.size))
        out["pd_flag"].append("gap" if gaps.size else "")
    res = pres.copy()
    for k, v in out.items():
        res[k] = v
    valid = res["pd_on_tick"] > 0
    res["pd_on_lat_ms"] = np.where(valid, (res["pd_on_tick"] - res["on_tick"]) / ms, np.nan)
    res["pd_off_lat_ms"] = np.where(valid & (res["off_tick"] > 0),
                                    (res["pd_off_tick"] - res["off_tick"]) / ms, np.nan)
    return res, levels
