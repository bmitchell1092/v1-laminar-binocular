"""Continuous signals on the NSP clock: LFP (.ns2) and threshold-crossing MUA spikes (.ns6).

Hardware filtering (from the NSx headers): .ns2 = 0.3 Hz (1st order) to 500 Hz
(4th order), 1 kHz; .ns6 = 0.3 Hz to 7.5 kHz, 30 kHz. Electrodes are in uV.
"""

from __future__ import annotations

import re

import numpy as np
from scipy import signal

from .io.nsx import NsxFile


def probe_labels(nsx: NsxFile, bank: str) -> list[str]:
    """Electrode labels of one probe bank ('A' -> eA01, eA02, ...), sorted by contact number.

    Contact number = label number. Whether number increases with depth is
    decided by laminar analysis, not assumed here.
    """
    pat = re.compile(rf"^e{bank}(\d+)$")
    labs = [(int(m.group(1)), lab) for lab in nsx.labels if (m := pat.match(lab))]
    return [lab for _, lab in sorted(labs)]


def read_lfp(ns2: NsxFile, labels: list[str], start: int | None = None, stop: int | None = None
             ) -> tuple[np.ndarray, np.ndarray]:
    """LFP [n_samples, n_channels] in uV and sample ticks (NSP clock)."""
    start = ns2.t_start if start is None else start
    stop = ns2.t_stop if stop is None else stop
    return ns2.read(start, stop, labels)


def threshold_spikes(ns6: NsxFile, labels: list[str], start: int | None = None, stop: int | None = None,
                     highpass_hz: float = 500.0, k_sd: float = 4.0, window_s: float = 1.0,
                     chunk_s: float = 30.0, pad_s: float = 0.05, dead_ms: float = 1.0) -> list[np.ndarray]:
    """Multi-unit threshold crossings per channel (NSP ticks), following Davis et al. 2023:
    high-pass at `highpass_hz` (4th-order Butterworth, zero-phase), negative crossings of
    -k_sd x SD, with SD computed in consecutive `window_s` windows; `dead_ms` dead time.
    """
    start = ns6.t_start if start is None else start
    stop = ns6.t_stop if stop is None else stop
    fs = ns6.fs
    per = ns6.period
    sos = signal.butter(4, highpass_hz, "highpass", fs=fs, output="sos")
    win = int(round(window_s * fs))
    chunk = int(round(chunk_s * fs)) // win * win          # whole SD windows per chunk
    pad = int(round(pad_s * fs))
    dead = int(round(dead_ms * fs / 1000))
    out: list[list[np.ndarray]] = [[] for _ in labels]
    t = start
    while t < stop:
        t1 = min(t + chunk * per, stop)
        x, ticks = ns6.read(t - pad * per, t1 + pad * per, labels)
        valid = ~np.isnan(x)
        x = np.where(valid, x, 0.0)
        y = signal.sosfiltfilt(sos, x, axis=0)[pad:pad + (t1 - t) // per]
        ok = valid[pad:pad + y.shape[0]]
        tk = ticks[pad:pad + y.shape[0]]
        nwin = int(np.ceil(y.shape[0] / win))
        for w in range(nwin):
            seg, segok = y[w * win:(w + 1) * win], ok[w * win:(w + 1) * win]
            sd = np.array([seg[segok[:, c], c].std() if segok[:, c].any() else np.inf for c in range(seg.shape[1])])
            below = seg < -k_sd * sd
            cross = below[1:] & ~below[:-1]
            for c in range(seg.shape[1]):
                idx = np.flatnonzero(cross[:, c]) + 1
                if idx.size > 1:
                    keep = np.r_[True, np.diff(idx) > dead]
                    idx = idx[keep]
                out[c].append(tk[w * win + idx])
        t = t1
    return [np.concatenate(o) if o else np.empty(0, np.int64) for o in out]


def analog_mua_epochs(ns6: NsxFile, labels: list[str], event_ticks: np.ndarray,
                      pre_ms: float = 100, post_ms: float = 200, pad_ms: float = 50,
                      band: tuple[float, float] = (500.0, 5000.0), lowpass_hz: float = 200.0) -> np.ndarray:
    """Analog MUA envelope around events: [events, time (1 kHz, -pre..post-1 ms), channels], uV.

    One documented chain (Super & Roelfsema 2005 style): 4th-order Butterworth band-pass
    `band`, full-wave rectification, 4th-order low-pass `lowpass_hz`, all zero-phase,
    then every 30th sample (1 kHz). Edges are padded by `pad_ms` before filtering.
    """
    fs = ns6.fs
    per = ns6.period
    bp = signal.butter(4, band, "bandpass", fs=fs, output="sos")
    lp = signal.butter(4, lowpass_hz, "lowpass", fs=fs, output="sos")
    step = int(round(fs / 1000))
    pre, post, pad = (int(round(v * fs / 1000)) for v in (pre_ms, post_ms, pad_ms))
    out = []
    for t in np.asarray(event_ticks):
        x, _ = ns6.read(t - (pre + pad) * per, t + (post + pad) * per, labels)
        x = np.nan_to_num(x)
        env = signal.sosfiltfilt(lp, np.abs(signal.sosfiltfilt(bp, x, axis=0)), axis=0)
        out.append(env[pad:pad + pre + post:step])
    return np.stack(out) if out else np.empty((0, (pre + post) // step, len(labels)))
