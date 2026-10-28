"""Laminar markers: spike-phase coupling (primary), spectrolaminar crossing and earliest CSD sink.

Marker hierarchy (docs/06): SPI vs depth decides probe orientation and the
input/deep boundary (Davis et al. 2023, eLife 12:e84512); the L4/L5 boundary is
cross-checked with the spectrolaminar crossing (Mendoza-Halliday et al. 2024)
and the earliest stimulus-evoked CSD sink.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import interpolate, signal


# ----------------------------------------------------------------------------- phase
def bandpass(x: np.ndarray, fs: float, band: tuple[float, float], order: int = 4) -> np.ndarray:
    sos = signal.butter(order, band, "bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x, axis=0)


def analytic(x: np.ndarray) -> np.ndarray:
    return signal.hilbert(x, axis=0)


def generalized_phase(x: np.ndarray, fs: float, lp: float, nwin: int = 3) -> np.ndarray:
    """Generalized phase (Davis et al. 2020 Nature; port of mullerlab/generalized-phase).

    `x` is wideband-filtered LFP [n, channels]. Phase of the analytic signal,
    with epochs of negative or sub-`lp` instantaneous frequency ("riding
    cycles") and the following (nwin-1) x their length replaced by
    shape-preserving cubic interpolation of the unwrapped phase.
    """
    z = analytic(x)
    ph = np.angle(z)
    inst_f = np.zeros_like(ph)
    inst_f[:-1] = np.angle(z[1:] * np.conj(z[:-1])) * fs / (2 * np.pi)
    if np.nanmean(inst_f) < 0:
        ph, inst_f = -ph, -inst_f
    out = np.empty_like(ph)
    t = np.arange(ph.shape[0])
    for c in range(ph.shape[1]):
        bad = inst_f[:, c] < lp
        bad[0] = False
        edges = np.diff(np.r_[0, bad.astype(np.int8), 0])
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1) - 1
        mask = bad.copy()
        for s, e in zip(starts, ends):
            mask[s:min(len(mask), s + (e - s) * nwin + 1)] = True
        p = np.unwrap(ph[:, c])
        if mask.any() and (~mask).sum() > 3:
            p[mask] = interpolate.PchipInterpolator(t[~mask], p[~mask], extrapolate=True)(t[mask])
        out[:, c] = np.angle(np.exp(1j * p))
    return out


# --------------------------------------------------------------------- spike-phase
@dataclass
class SpikePhase:
    spi: np.ndarray    # [spike channel, LFP channel] mean resultant length
    theta: np.ndarray  # preferred phase (rad)
    ppc: np.ndarray    # pairwise phase consistency (Vinck et al. 2010), unbiased by spike count
    n: np.ndarray      # spikes used per spike channel


def spike_phase(spike_ticks: list[np.ndarray], phase: np.ndarray, phase_ticks: np.ndarray,
                max_spikes: int = 20000, seed: int = 0) -> SpikePhase:
    """Spike-phase coupling between every spike channel and every LFP-phase channel."""
    rng = np.random.default_rng(seed)
    t0, per = phase_ticks[0], phase_ticks[1] - phase_ticks[0]
    n_s, n_l = len(spike_ticks), phase.shape[1]
    spi, theta, ppc = (np.full((n_s, n_l), np.nan) for _ in range(3))
    n = np.zeros(n_s, int)
    for i, st in enumerate(spike_ticks):
        idx = np.round((np.asarray(st) - t0) / per).astype(np.int64)
        idx = idx[(idx >= 0) & (idx < phase.shape[0])]
        if idx.size > max_spikes:
            idx = rng.choice(idx, max_spikes, replace=False)
        n[i] = idx.size
        if idx.size < 50:
            continue
        z = np.exp(1j * phase[idx, :]).sum(axis=0)
        spi[i] = np.abs(z) / idx.size
        theta[i] = np.angle(z)
        ppc[i] = (np.abs(z) ** 2 - idx.size) / (idx.size * (idx.size - 1))
    return SpikePhase(spi, theta, ppc, n)


# ------------------------------------------------------------------ spectrolaminar
def relative_power(lfp: np.ndarray, fs: float, nperseg: int = 1024) -> tuple[np.ndarray, np.ndarray]:
    """Welch PSD per channel, each frequency divided by its maximum across channels.

    Returns (freqs, relpow [freq, channel]). The alpha-beta (10-19 Hz) and gamma
    (75-150 Hz) profiles cross near layer 4 (Mendoza-Halliday et al. 2024).
    """
    good = ~np.isnan(lfp).any(axis=1)
    f, p = signal.welch(lfp[good], fs=fs, nperseg=nperseg, axis=0)
    return f, p / p.max(axis=1, keepdims=True)


def band_profile(f: np.ndarray, relpow: np.ndarray, band: tuple[float, float]) -> np.ndarray:
    sel = (f >= band[0]) & (f <= band[1])
    return relpow[sel].mean(axis=0)


# ----------------------------------------------------------------------------- CSD
def csd(lfp: np.ndarray, spacing_mm: float = 0.1, conductivity: float = 0.4) -> np.ndarray:
    """Standard CSD, -sigma d2V/dz2, along the last axis (channels in contact order).

    lfp in uV, spacing in mm, sigma in S/m -> nA/mm^3. Edge channels are NaN.
    Negative values are sinks.
    """
    v = np.asarray(lfp, float)
    out = np.full_like(v, np.nan)
    # (S/m) x (uV/mm^2) = 1e-6 A/(m mm^2) = 1 nA/mm^3
    out[..., 1:-1] = -conductivity * (v[..., 2:] - 2 * v[..., 1:-1] + v[..., :-2]) / spacing_mm ** 2
    return out


def earliest_sink(csd_avg: np.ndarray, t_ms: np.ndarray, window_ms: tuple[float, float] = (20, 100),
                  min_rel_amp: float = 0.5, onset_frac: float = 0.25) -> tuple[np.ndarray, int]:
    """Onset of the evoked sink per channel and the channel whose sink starts first.

    Amplitude-based, so it does not depend on the (tiny) baseline noise of a heavily
    averaged CSD. Only channels whose sink in `window_ms` reaches `min_rel_amp` of the
    strongest sink count; their onset is the first time the CSD falls below
    `onset_frac` of their own peak sink. `csd_avg` is [time, channel], sinks negative.
    Returns (onset_ms per channel, NaN where not a qualifying sink; earliest channel or -1).
    """
    win = np.flatnonzero((t_ms >= window_ms[0]) & (t_ms <= window_ms[1]))
    peak = np.nanmin(csd_avg[win], axis=0)
    strongest = np.nanmin(peak)
    onset = np.full(csd_avg.shape[1], np.nan)
    if not np.isfinite(strongest) or strongest >= 0:
        return onset, -1
    for c in np.flatnonzero(peak <= min_rel_amp * strongest):
        hit = np.flatnonzero(csd_avg[win, c] <= onset_frac * peak[c])
        onset[c] = t_ms[win[hit[0]]]
    return onset, int(np.nanargmin(onset)) if np.isfinite(onset).any() else -1
