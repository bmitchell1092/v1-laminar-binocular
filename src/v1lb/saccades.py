"""Saccade and microsaccade detection on calibrated eye traces (Engbert & Kliegl 2003; Engbert & Mergenthaler 2006).

Velocity from a 5-sample moving window; per-component threshold = `lam` x a median-based SD of the
velocity; a saccade is a run of samples outside the elliptic threshold lasting at least `min_dur_ms`;
runs closer than `merge_ms` are merged (overshoots). Samples flagged invalid (blinks, signal loss) are never
part of a saccade and break runs. Events above `max_peak_velocity` (macaque saccades peak near 1000 deg/s
even at 30 deg) or `max_amp_deg` are tracker artefacts (signal jumps) and are dropped.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def velocity(x: np.ndarray, fs: float) -> np.ndarray:
    """(x[n+2] + x[n+1] - x[n-1] - x[n-2]) * fs / 6; NaN at the two samples at each end."""
    v = np.full(len(x), np.nan)
    v[2:-2] = (x[4:] + x[3:-1] - x[1:-3] - x[:-4]) * fs / 6.0
    return v


def _median_sd(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(np.sqrt(max(np.median(v ** 2) - np.median(v) ** 2, 1e-12)))


def detect(x: np.ndarray, y: np.ndarray, fs: float = 1000.0, lam: float = 6.0, min_dur_ms: float = 6.0,
           merge_ms: float = 20.0, valid: np.ndarray | None = None, max_amp_deg: float = 30.0,
           max_peak_velocity: float = 1500.0) -> pd.DataFrame:
    """Saccades in one trace (deg). Returns onset/offset sample indices, duration (ms), peak velocity (deg/s),
    amplitude (deg; largest excursion within the saccade), dx, dy (deg; offset minus onset position)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if valid is not None:
        ok &= np.asarray(valid, bool)
    vx, vy = velocity(x, fs), velocity(y, fs)
    good = ok & np.isfinite(vx) & np.isfinite(vy)
    for k in (-2, -1, 1, 2):                               # the velocity window must not touch invalid samples
        good &= np.roll(ok, k)
    cols = ["onset", "offset", "duration_ms", "peak_velocity", "amplitude", "dx", "dy"]
    if good.sum() < 10:
        return pd.DataFrame(columns=cols)
    ex, ey = lam * _median_sd(vx[good]), lam * _median_sd(vy[good])
    above = good & ((vx / ex) ** 2 + (vy / ey) ** 2 > 1)
    d = np.diff(np.concatenate([[0], above.astype(int), [0]]))
    on, off = np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1      # inclusive runs
    if len(on) == 0:
        return pd.DataFrame(columns=cols)
    runs = [[on[0], off[0]]]
    merge = int(round(merge_ms * fs / 1000))
    for a, b in zip(on[1:], off[1:]):
        if a - runs[-1][1] <= merge and good[runs[-1][1]:a].all():
            runs[-1][1] = b
        else:
            runs.append([a, b])
    rows = []
    min_n = int(np.ceil(min_dur_ms * fs / 1000))
    for a, b in runs:
        if b - a + 1 < min_n:
            continue
        seg_x, seg_y = x[a:b + 1], y[a:b + 1]
        amp = float(np.hypot(seg_x.max() - seg_x.min(), seg_y.max() - seg_y.min()))
        peak_v = float(np.nanmax(np.hypot(vx[a:b + 1], vy[a:b + 1])))
        if amp > max_amp_deg or peak_v > max_peak_velocity:
            continue
        rows.append(dict(onset=int(a), offset=int(b), duration_ms=(b - a + 1) * 1000 / fs,
                         peak_velocity=peak_v, amplitude=amp,
                         dx=float(x[b] - x[a]), dy=float(y[b] - y[a])))
    return pd.DataFrame(rows, columns=cols)


def valid_mask(x: np.ndarray, y: np.ndarray, pupil: np.ndarray | None = None, max_deg: float = 20.0,
               pupil_mad: float = 6.0, pad_ms: float = 50.0, fs: float = 1000.0) -> np.ndarray:
    """Samples usable for saccade detection: gaze within `max_deg` of fixation and pupil within `pupil_mad`
    robust SDs of its median (blinks), widened by `pad_ms` on each side."""
    ok = np.isfinite(x) & np.isfinite(y) & (np.abs(x) < max_deg) & (np.abs(y) < max_deg)
    if pupil is not None and np.isfinite(pupil).any():
        med = np.nanmedian(pupil)
        mad = np.nanmedian(np.abs(pupil - med)) / 0.6745
        ok &= np.abs(pupil - med) < pupil_mad * mad
    pad = int(round(pad_ms * fs / 1000))
    return np.convolve(~ok, np.ones(2 * pad + 1), "same") == 0


def detect_continuous(ticks: np.ndarray, x: np.ndarray, y: np.ndarray, pupil: np.ndarray | None = None,
                      segment_s: float = 60.0, fs: float = 1000.0, **kw) -> pd.DataFrame:
    """Saccades in a whole-recording trace on the NSP clock (thresholds re-estimated per `segment_s`, so slow
    changes in tracker noise do not bias detection). Adds onset_tick / offset_tick."""
    valid = valid_mask(x, y, pupil, fs=fs)
    seg = int(segment_s * fs)
    out = []
    for s0 in range(0, len(x), seg):
        d = detect(x[s0:s0 + seg], y[s0:s0 + seg], fs=fs, valid=valid[s0:s0 + seg], **kw)
        out.append(d.assign(onset=d.onset + s0, offset=d.offset + s0))
    s = pd.concat(out, ignore_index=True)
    s["onset_tick"] = np.asarray(ticks)[s.onset.to_numpy(int)]
    s["offset_tick"] = np.asarray(ticks)[s.offset.to_numpy(int)]
    return s


def around_onsets(sacc: pd.DataFrame, onsets: np.ndarray, window_ms=(-200, 500), micro_deg: float = 1.0,
                  ticks_per_ms: int = 30) -> pd.DataFrame:
    """Per stimulus onset (NSP ticks): saccades / microsaccades starting within `window_ms`, time from the
    end of the last saccade before onset (ms; the fixation age) and to the next saccade onset after it."""
    on = np.sort(sacc.onset_tick.to_numpy())
    order = np.argsort(sacc.onset_tick.to_numpy())
    off = sacc.offset_tick.to_numpy()[order]
    amp = sacc.amplitude.to_numpy()[order]
    rows = []
    for t in np.asarray(onsets):
        a, b = np.searchsorted(on, [t + window_ms[0] * ticks_per_ms, t + window_ms[1] * ticks_per_ms])
        prev = np.searchsorted(on, t) - 1
        nxt = np.searchsorted(on, t)
        rows.append(dict(n_saccades=int(b - a), n_micro=int(np.sum(amp[a:b] < micro_deg)),
                         n_large=int(np.sum(amp[a:b] >= micro_deg)),
                         since_saccade_ms=(t - off[prev]) / ticks_per_ms if prev >= 0 else np.nan,
                         to_saccade_ms=(on[nxt] - t) / ticks_per_ms if nxt < len(on) else np.nan))
    return pd.DataFrame(rows)
