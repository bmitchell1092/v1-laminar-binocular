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


# ------------------------------------------------------------- template fitting
_TEMPLATE_BOUNDARY = 32.5  # reversal of the Davis et al. average phase profile (between indices 32 and 33)


def davis_templates() -> dict[str, np.ndarray]:
    """Average laminar profiles from Davis et al. 2023 (index increases with depth); see data/ provenance."""
    from importlib.resources import files
    with np.load(files("v1lb") / "data" / "davis2023_templates.npz") as d:
        return {k: d[k] for k in d.files}


def phase_profile(sp: SpikePhase, min_spikes: int = 500) -> np.ndarray:
    """Preferred phase per LFP channel: SPI-weighted circular mean over spike channels with enough spikes."""
    good = sp.n >= min_spikes
    if not good.any():
        return np.full(sp.theta.shape[1], np.nan)
    return np.angle(np.nansum(np.exp(1j * sp.theta[good]) * sp.spi[good], axis=0))


@dataclass
class TemplateFit:
    orientation: str      # 'ascending' = contact number increases with depth
    shift: int
    score: float          # best score (mean cos for phase, Pearson r for power)
    score_other: float    # best score with the opposite orientation
    boundary: float       # contact number (1-based) of the template boundary; may lie off the probe
    n_overlap: int

    @property
    def margin(self) -> float:
        return self.score - self.score_other


def fit_template(profile: np.ndarray, template: np.ndarray, kind: str = "phase",
                 weight: np.ndarray | None = None, boundary_index: float = _TEMPLATE_BOUNDARY,
                 min_overlap: int = 12, edge_weight: float = 0.5) -> TemplateFit:
    """Slide a depth template along a probe profile in both orientations.

    kind='phase': score = weighted mean cos(profile - template) (absolute phase convention,
    validated against Davis et al. 2023). kind='linear': Pearson r (e.g. gamma - alpha power).
    Every finite profile channel must land on the template; the template's end plateaus are
    extended (weight x `edge_weight`) so the boundary may fall anywhere along the probe.
    """
    profile = np.asarray(profile, float)
    n = profile.size
    w_t = np.ones(template.size) if weight is None else np.asarray(weight, float)
    good = np.isfinite(template) & (w_t > 0)
    lo, hi = np.flatnonzero(good)[[0, -1]]
    # extend the template's end plateaus (at reduced weight) so the boundary can sit anywhere on the probe
    core_t = np.where(good[lo:hi + 1], template[lo:hi + 1], 0.0)
    core_w = np.where(good[lo:hi + 1], w_t[lo:hi + 1], 0.0)
    template = np.r_[np.full(n, core_t[0]), core_t, np.full(n, core_t[-1])]
    w_t = np.r_[np.full(n, core_w[0] * edge_weight), core_w, np.full(n, core_w[-1] * edge_weight)]
    boundary_index = boundary_index - lo + n
    m = template.size
    usable = w_t > 0
    c = np.flatnonzero(np.isfinite(profile))   # every usable probe channel must land on the template
    if c.size < min_overlap:
        raise ValueError("profile has too few finite channels for template fitting")
    best = {}
    for orient in ("ascending", "descending"):
        for s in range(-n, m + n):
            k = s + c if orient == "ascending" else s - c
            if k.min() < 0 or k.max() >= m or not usable[k].all():
                continue
            ok = c
            p, t, w = profile[ok], template[k], w_t[k]
            if kind == "phase":
                score = float(np.sum(w * np.cos(p - t)) / np.sum(w))
            else:
                if np.std(p) == 0 or np.std(t) == 0:
                    continue
                score = float(np.corrcoef(p, t)[0, 1])
            if orient not in best or score > best[orient][0]:
                best[orient] = (score, s, int(ok.size))
    if not best:
        raise ValueError("profile too short or too sparse for template fitting")
    orient = max(best, key=lambda o: best[o][0])
    other = "descending" if orient == "ascending" else "ascending"
    score, s, n_ok = best[orient]
    boundary = (boundary_index - s + 1) if orient == "ascending" else (s - boundary_index + 1)
    return TemplateFit(orient, s, score, best[other][0] if other in best else float("nan"), boundary, n_ok)


@dataclass
class SpectralCrossing:
    orientation: str   # 'ascending' if the gamma peak lies at a lower contact number than the alpha-beta peak
    gamma_peak: int    # contact numbers (1-based)
    alpha_peak: int
    crossing: float    # contact position where (gamma - alpha) changes sign between the peaks (NaN if none)


def spectral_crossing(alpha_beta: np.ndarray, gamma: np.ndarray, smooth: int = 3) -> SpectralCrossing:
    """Spectrolaminar motif (Mendoza-Halliday et al. 2024): gamma relative power peaks superficially,
    alpha-beta deep; their crossing marks layer 4. Peak-based, so it does not depend on LFP polarity."""
    k = np.ones(smooth) / smooth
    a = np.convolve(np.asarray(alpha_beta, float), k, mode="same")
    g = np.convolve(np.asarray(gamma, float), k, mode="same")
    ga, aa = int(np.argmax(g)), int(np.argmax(a))
    orient = "ascending" if ga < aa else "descending"
    lo, hi = sorted((ga, aa))
    d = g[lo:hi + 1] - a[lo:hi + 1]
    s = np.flatnonzero(np.sign(d[:-1]) != np.sign(d[1:]))
    if s.size:
        i = s[np.argmin(np.abs(s + 0.5 - (hi - lo) / 2))]          # sign change nearest the midpoint
        crossing = lo + i + d[i] / (d[i] - d[i + 1]) + 1             # linear interpolation, 1-based
    else:
        crossing = float("nan")
    return SpectralCrossing(orient, ga + 1, aa + 1, float(crossing))
