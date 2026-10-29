"""Eye position and pupil from the Blackrock analog inputs, calibrated against MonkeyLogic.

The eye tracker's analog outputs were split to MonkeyLogic and to Blackrock `ainp*` channels;
Central gives them no labels. MonkeyLogic stores the calibrated eye position (degrees) per
trial in the .bhv2 (`AnalogData.Eye`, time 0 = trial start = first `9` code). Regressing the
Blackrock channels on it identifies which channels are X and Y and gives an affine
calibration from mV to degrees. Pupil is not in the .bhv2; it is the remaining channel that
carries signal (blink transients, slow drift) and does not track gaze.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import recording_files
from .events import ml_events, split_trials
from .io import NevFile, NsxFile, read_stim_text
from .io.bhv2 import read_bhv2
from .trials import match_trials


@dataclass
class EyeCalibration:
    x_channel: str
    y_channel: str
    affine: np.ndarray                 # 2 x 3: [deg_x, deg_y] = A @ [mv_x, mv_y, 1]
    r2: tuple[float, float]            # fit quality for deg_x, deg_y
    channel_r2: dict = field(default_factory=dict)   # per analog channel: (r2 vs ML x, r2 vs ML y)
    live_channels: list = field(default_factory=list)  # channels carrying signal (not flat noise)
    fixation: tuple = (0.0, 0.0)       # fixation cross in MonkeyLogic screen degrees (stereoscope: left-half centre)


def calibrate(session: str, recording: str, max_trials: int = 200, flat_sd_mv: float = 5.0) -> EyeCalibration:
    files = recording_files(session, recording)
    ns2 = NsxFile(files[".ns2"])
    ain = [lab for lab in ns2.labels if lab.startswith("ainp") and lab != "ainp1"]   # ainp1 = photodiode
    nev = NevFile(files[".nev"])
    ticks, values, _ = nev.digital_events()
    ticks, codes, _ = ml_events(ticks, values)
    trials = split_trials(codes, ticks)
    match = match_trials(trials, nev.origin_utc, read_stim_text(files[".gBMINTEROCGrating_di"]))
    bhv = read_bhv2(files[".bhv2"])
    mc = bhv["MLConfig"]
    ppd = abs(float(np.ravel(mc["PixelsPerDegree"])[0]))   # stored as [ppd_x, -ppd_y]
    width_px = float(str(mc["Resolution"]).split("x")[0])
    fixation = (-0.25 * width_px / ppd, 0.0)      # T_binocular_gratings: fix1.Target = [-0.25*scrsize(1) 0]
    X, Y = [], []
    for _, m in match.table.iloc[:max_trials].iterrows():
        t = bhv[f"Trial{int(m.ml_trial)}"]
        eye = np.asarray(t["AnalogData"]["Eye"], float)
        if eye.ndim != 2 or eye.shape[0] < 10:
            continue
        dt_ms = float(t["AnalogData"]["SampleInterval"])
        start = int(trials.loc[trials.nev_trial == m.nev_trial, "start_tick"].iloc[0])
        sample_ticks = start + np.round(np.arange(eye.shape[0]) * dt_ms * 30).astype(np.int64)
        x, tk = ns2.read(sample_ticks[0], sample_ticks[-1] + 30, ain)
        idx = np.clip(np.searchsorted(tk, sample_ticks), 0, len(tk) - 1)
        X.append(x[idx])
        Y.append(eye)
    X, Y = np.concatenate(X), np.concatenate(Y)
    ok = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)
    X, Y = X[ok], Y[ok]
    live = [lab for i, lab in enumerate(ain) if np.std(X[:, i]) > flat_sd_mv]
    # drop blinks / off-screen samples: they move every channel at once and inflate all correlations
    keep = (np.abs(Y) < 20).all(axis=1)
    for lab in live:
        v = X[:, ain.index(lab)]
        med = np.median(v)
        mad = np.median(np.abs(v - med)) / 0.6745 + 1e-9
        keep &= np.abs(v - med) < 5 * mad
    X, Y = X[keep], Y[keep]
    channel_r2 = {lab: tuple(float(np.corrcoef(X[:, i], Y[:, k])[0, 1] ** 2) if np.std(X[:, i]) > 0 else 0.0
                             for k in (0, 1)) for i, lab in enumerate(ain)}

    def fit(xc, yc):
        D = np.column_stack([X[:, ain.index(xc)], X[:, ain.index(yc)], np.ones(len(X))])
        A, *_ = np.linalg.lstsq(D, Y, rcond=None)
        pred = D @ A
        return A.T, tuple(float(1 - np.var(Y[:, k] - pred[:, k]) / np.var(Y[:, k])) for k in (0, 1))

    pairs = [(a, b) for a in live for b in live if a != b]
    xc, yc = max(pairs, key=lambda p: (lambda r: r[0] + r[1])(fit(*p)[1]))
    A, r2 = fit(xc, yc)
    return EyeCalibration(xc, yc, A, r2, channel_r2, live, fixation)


PUPIL_CHANNEL = "ainp5"   # verified by the pupillary light reflex to full-field flashes (B and J; 05 C11)


def to_degrees(cal: EyeCalibration, mv_x: np.ndarray, mv_y: np.ndarray, relative: bool = True
               ) -> tuple[np.ndarray, np.ndarray]:
    """mV -> degrees; relative=True gives gaze relative to the fixation cross (the stereoscope puts
    it at the centre of the left half of the screen, about -10 to -12 deg in screen coordinates)."""
    A = cal.affine
    x = A[0, 0] * mv_x + A[0, 1] * mv_y + A[0, 2]
    y = A[1, 0] * mv_x + A[1, 1] * mv_y + A[1, 2]
    if relative:
        x, y = x - cal.fixation[0], y - cal.fixation[1]
    return x, y


def continuous(ns2: NsxFile, cal: EyeCalibration, start: int | None = None, stop: int | None = None):
    """Whole-recording eye trace on the NS2 grid: (ticks, x_deg, y_deg, pupil_mv)."""
    chans = [cal.x_channel, cal.y_channel] + ([PUPIL_CHANNEL] if PUPIL_CHANNEL in ns2.labels else [])
    x, ticks = ns2.read(ns2.t_start if start is None else start, ns2.t_stop if stop is None else stop, chans)
    ex, ey = to_degrees(cal, x[:, 0], x[:, 1])
    pupil = x[:, 2] if x.shape[1] > 2 else np.full(len(ticks), np.nan, np.float32)
    return ticks, ex.astype(np.float32), ey.astype(np.float32), pupil.astype(np.float32)
