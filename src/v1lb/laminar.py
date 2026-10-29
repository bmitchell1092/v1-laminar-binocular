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
HAMMING3 = np.array([0.23, 0.54, 0.23])


def hamming_smooth(x: np.ndarray) -> np.ndarray:
    """3-point Hamming smoothing along the last axis (channels), edges replicated.

    Standard before a laminar CSD (Ulbert et al. 2001): it passes 8% of a pattern that
    alternates between neighbouring contacts, the signature of contact gain mismatch."""
    v = np.asarray(x, float)
    p = np.concatenate([v[..., :1], v, v[..., -1:]], axis=-1)
    return HAMMING3[0] * p[..., :-2] + HAMMING3[1] * p[..., 1:-1] + HAMMING3[2] * p[..., 2:]


def bad_channels(evoked: np.ndarray, lo: float = 0.6, hi: float = 1.6, half_width: int = 2) -> np.ndarray:
    """Contacts whose evoked LFP amplitude is out of line with their neighbours.

    `evoked` is a mean evoked LFP [time, channel] (channels in contact order). A contact is bad
    when its SD over time is outside [lo, hi] x the median SD of the contacts within
    `half_width` on either side. The LFP is smooth over a few hundred um, so a lone weak or
    strong contact is a contact or gain problem, not a laminar feature."""
    amp = np.nanstd(np.asarray(evoked, float), axis=0)
    n = amp.size
    bad = np.zeros(n, bool)
    for k in range(n):
        nb = [j for j in range(k - half_width, k + half_width + 1) if j != k and 0 <= j < n]
        ref = np.nanmedian(amp[nb])
        bad[k] = not (lo * ref <= amp[k] <= hi * ref)
    return bad


def interpolate_channels(x: np.ndarray, bad: np.ndarray) -> np.ndarray:
    """Replace bad channels (last axis) by linear interpolation between the nearest good ones;
    bad channels at either end copy the nearest good channel."""
    v = np.array(x, float, copy=True)
    good = np.flatnonzero(~np.asarray(bad, bool))
    if good.size == 0 or good.size == v.shape[-1]:
        return v
    for k in np.flatnonzero(bad):
        below, above = good[good < k], good[good > k]
        if below.size and above.size:
            g0, g1 = below[-1], above[0]
            v[..., k] = v[..., g0] + (v[..., g1] - v[..., g0]) * (k - g0) / (g1 - g0)
        else:
            v[..., k] = v[..., below[-1] if below.size else above[0]]
    return v


def csd(lfp: np.ndarray, spacing_mm: float = 0.1, conductivity: float = 0.4, smooth: bool = False,
        bad: np.ndarray | None = None) -> np.ndarray:
    """Standard CSD, -sigma d2V/dz2, along the last axis (channels in contact order).

    lfp in uV, spacing in mm, sigma in S/m -> nA/mm^3. Edge channels are NaN.
    Negative values are sinks. `bad` contacts are interpolated and `smooth` applies
    3-point Hamming smoothing before the second difference.
    """
    v = np.asarray(lfp, float)
    if bad is not None and np.any(bad):
        v = interpolate_channels(v, bad)
    if smooth:
        v = hamming_smooth(v)
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


# --------------------------------------------------------- probe-level decisions
def superficial_order(n: int, orientation: str) -> np.ndarray:
    """Contact indices (0-based) ordered from the superficial end to the deep end."""
    return np.arange(n) if orientation == "ascending" else np.arange(n)[::-1]


def top_of_cortex(rate: np.ndarray, evoked: np.ndarray, orientation: str, rate_floor: float = 2.5,
                  evoked_thresh: float = 6.0, run: int = 3) -> tuple[int, bool, str]:
    """First contact (1-based), from the superficial end, that starts `run` consecutive active
    contacts. Spontaneous spiking (threshold-crossing rate above `rate_floor` Hz; chance at -4 SD
    is ~1 Hz) is local, so it is used first. Evoked analog MUA spreads a few hundred um by volume
    conduction, so it is only a fallback (e.g. probes whose spikes are too small to threshold).
    Returns (contact or -1, top_contact_already_in_cortex, criterion used)."""
    order = superficial_order(len(rate), orientation)
    for name, active in (("spontaneous MUA", np.asarray(rate)[order] > rate_floor),
                         ("evoked MUA", np.asarray(evoked)[order] > evoked_thresh)):
        for i in range(len(order) - run + 1):
            if active[i:i + run].all():
                return int(order[i]) + 1, i == 0, name
    return -1, False, "none"


def step_kind(signed: float, ambiguous_near_pi: float = 2.8, min_turn: float = 0.3) -> str:
    """'canonical' (within-column direction going deeper), 'reversed' (a fold entered from its
    deep side), or 'ambiguous' (|step| near pi or tiny, so its direction is undefined)."""
    if min_turn < signed < ambiguous_near_pi:
        return "canonical"
    if -ambiguous_near_pi < signed < -min_turn:
        return "reversed"
    return "ambiguous"


def phase_steps(profile: np.ndarray, orientation: str, start_contact: int = 1, stop_contact: int | None = None,
                w: int = 3, min_step: float = 0.4, min_sep: int = 3) -> list[tuple[float, float, float]]:
    """Abrupt changes of preferred phase along the probe, from the superficial end down.

    For each split between neighbouring contacts, the step is the circular difference between
    the mean phase of `w` contacts below and `w` above. Local maxima of |step| above `min_step`
    rad, at least `min_sep` contacts apart, between `start_contact` (top of cortex) and
    `stop_contact` (deepest active contact). Returns [(boundary position as contact number,
    e.g. 23.5; |step| rad; signed step rad)], superficial to deep. See `step_kind` for how the
    sign is read (Davis et al.: +-pi superficial/input -> ~-1.2 deep is a positive step).
    """
    n = len(profile)
    order = superficial_order(n, orientation)
    q = np.exp(1j * np.asarray(profile, float)[order])
    first = int(np.flatnonzero(order == start_contact - 1)[0]) if start_contact > 0 else 0
    last = int(np.flatnonzero(order == stop_contact - 1)[0]) if stop_contact else n - 1
    signed = np.zeros(n - 1)
    for i in range(max(first, w - 1), min(n - w, last - w + 1)):
        a, b = np.nanmean(q[i - w + 1:i + 1]), np.nanmean(q[i + 1:i + w + 1])
        signed[i] = np.angle(b * np.conj(a))
    mag = np.abs(signed)
    steps = []
    for i in np.argsort(mag)[::-1]:
        if mag[i] < min_step:
            break
        if all(abs(i - j) >= min_sep for j, _ in steps):
            steps.append((i, mag[i]))
    steps.sort()
    return [((order[i] + order[i + 1]) / 2 + 1, float(m), float(signed[i])) for i, m in steps]


def l45_from_steps(steps: list[tuple[float, float, float]]) -> float:
    """Blake's rule (2026-10-28): among non-reversed steps below the top of cortex, the 2nd if
    there are two or more, otherwise the only one; NaN if none. Whether two steps are
    supra->input + input->deep in one column, or two folds of cortex, is judged separately
    (see `laminar_report`)."""
    canon = [s for s in steps if step_kind(s[2]) != "reversed"]
    if not canon:
        return float("nan")
    return canon[1][0] if len(canon) >= 2 else canon[0][0]


def mua_gap(rate: np.ndarray, a: float, b: float, rate_floor: float = 2.5, min_len: int = 2) -> bool:
    """True if >= `min_len` consecutive contacts strictly between boundary positions a and b
    spike at or below the noise floor (e.g. white matter between two folds)."""
    lo, hi = sorted((a, b))
    idx = [c for c in range(int(np.ceil(lo)), int(np.floor(hi)) + 1) if lo < c < hi]
    quiet = [rate[c - 1] <= rate_floor for c in idx]
    run = best = 0
    for q in quiet:
        run = run + 1 if q else 0
        best = max(best, run)
    return best >= min_len
