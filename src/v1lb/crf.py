"""Contrast-response functions (CRFs) and supersaturation.

Responses are window means of a [presentation, time, channel] signal (MUA envelope, spike counts) minus
the mean over blank (0/0) presentations in the same window, so the spontaneous level is the blank response
and each CRF starts at 0. Monocular CRFs (contra or ipsi alone) and the dioptic CRF share the blank point.

Supersaturation index (SI) = (max_c R(c) - R(c_max)) / max_c R(c): 0 for a CRF that rises or saturates,
up to 1 when the response at the top contrast falls back to blank. Its uncertainty comes from resampling
presentations within each contrast. The model is the Naka-Rushton function with Peirce's (2007)
supersaturation exponent, R(c) = Rmax c^n / (c^(n s) + c50^(n s)): s = 1 saturates, s > 1 turns down at
high contrast. Fitted on single presentations, compared with s = 1 by an F-test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

CONTRASTS = (0.0, 0.055, 0.11, 0.225, 0.45, 0.9)
WINDOWS = {"transient": (40, 140), "sustained": (150, 250)}
CONDITIONS = {"mono_contra": ("contrast_contra", "contrast_ipsi"), "mono_ipsi": ("contrast_ipsi", "contrast_contra"),
              "dioptic": ("contrast_contra", "contrast_ipsi")}


def window_mean(x: np.ndarray, t_ms: np.ndarray, window) -> np.ndarray:
    """Mean over `window` (ms, [start, stop)) of x [presentation, time, ...] -> [presentation, ...]."""
    w = (t_ms >= window[0]) & (t_ms < window[1])
    return np.nanmean(x[:, w], axis=1)


def crf_points(resp: np.ndarray, P: pd.DataFrame, condition: str) -> tuple[np.ndarray, np.ndarray]:
    """Presentations of one CRF: (contrast per presentation, responses [n, ...]), blank (0/0) included.

    mono_contra: ipsi eye at 0; mono_ipsi: contra eye at 0; dioptic: both eyes equal."""
    a, b = CONDITIONS[condition]
    ca, cb = P[a].to_numpy(float), P[b].to_numpy(float)
    sel = np.isclose(ca, cb) if condition == "dioptic" else np.isclose(cb, 0)
    return ca[sel], resp[sel]


def means_by_contrast(c: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Mean response at each of CONTRASTS [len(CONTRASTS), ...]."""
    return np.stack([np.nanmean(r[np.isclose(c, k)], axis=0) for k in CONTRASTS])


def supersaturation_index(means: np.ndarray) -> np.ndarray:
    """SI along axis 0 of [contrast, ...] means (blank-subtracted; NaN where the peak is not above 0)."""
    peak = np.nanmax(means, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        si = (peak - means[-1]) / peak
    return np.where(peak > 0, si, np.nan)


def bootstrap_si(c: np.ndarray, r: np.ndarray, n_boot: int = 1000, seed: int = 0) -> np.ndarray:
    """Bootstrap SIs [n_boot, ...]: presentations resampled with replacement within each contrast."""
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(np.isclose(c, k)) for k in CONTRASTS]
    out = []
    for _ in range(n_boot):
        m = np.stack([np.nanmean(r[rng.choice(g, len(g))], axis=0) for g in groups])
        out.append(supersaturation_index(m))
    return np.array(out)


def cv_supersaturation_index(c: np.ndarray, r: np.ndarray, n_splits: int = 200, seed: int = 0) -> np.ndarray:
    """Cross-validated SI per site [...]: the plain SI is biased upward by noise (the largest of noisy means
    overshoots the true peak), so a monotonic CRF still gets SI > 0. Here presentations are split in random
    halves within each contrast; half A picks the peak contrast c*, half B measures the drop
    R_B(c*) - R_B(c_max). The mean drop over splits is divided by the peak of the full-data means (a fixed,
    positive scale for a driven site; normalising by half A's peak blows up when that is near 0). Its
    expectation is <= 0 for any monotonic CRF and > 0 only for a real downturn; 0 when A peaks at the top."""
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(np.isclose(c, k)) for k in CONTRASTS]
    scale = np.nanmax(means_by_contrast(c, r)[1:], axis=0)
    acc = []
    for _ in range(n_splits):
        a_idx, b_idx = [], []
        for g in groups:
            p = rng.permutation(g)
            a_idx.append(p[: len(p) // 2])
            b_idx.append(p[len(p) // 2:])
        mA = np.stack([np.nanmean(r[i], axis=0) for i in a_idx])
        mB = np.stack([np.nanmean(r[i], axis=0) for i in b_idx])
        k = np.argmax(mA[1:], axis=0) + 1                                  # peak among non-blank contrasts
        acc.append(np.take_along_axis(mB, k[None], axis=0)[0] - mB[-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(scale > 0, np.nanmean(np.array(acc), axis=0) / scale, np.nan)


def downturn_p(c: np.ndarray, r: np.ndarray, n_splits: int = 50, seed: int = 0) -> np.ndarray:
    """Held-out one-sided test that the response falls from the peak contrast to the top contrast, per site.

    Each split: half A picks the peak non-blank contrast c*; half B tests R(c*) > R(c_max) with a one-sided
    Welch t-test (p = 1 when c* is the top contrast). Selection and test use disjoint presentations, so each
    split's p is valid; splits are combined as min(1, 2 x median p) (Meinshausen et al. 2009), valid whatever
    the dependence between splits. Under a non-decreasing CRF at most 5% of sites reach p < 0.05."""
    rng = np.random.default_rng(seed)
    r = r.reshape(len(r), -1)
    groups = [np.flatnonzero(np.isclose(c, k)) for k in CONTRASTS]
    ps = np.ones((n_splits, r.shape[1]))
    for i in range(n_splits):
        halves = [rng.permutation(g) for g in groups]
        a_idx = [h[: len(h) // 2] for h in halves]
        b_idx = [h[len(h) // 2:] for h in halves]
        mA = np.stack([np.nanmean(r[a], axis=0) for a in a_idx])
        k = np.argmax(mA[1:], axis=0) + 1
        top = r[b_idx[-1]]
        for j in np.flatnonzero(k < len(CONTRASTS) - 1):
            x, y = r[b_idx[k[j]], j], top[:, j]
            if np.nanstd(x) > 0 or np.nanstd(y) > 0:
                ps[i, j] = stats.ttest_ind(x, y, equal_var=False, alternative="greater", nan_policy="omit").pvalue
    return np.minimum(1.0, 2 * np.nanmedian(ps, axis=0))


def naka_rushton(c, rmax, c50, n, s=1.0):
    c = np.asarray(c, float)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        return rmax * c ** n / (c ** (n * s) + c50 ** (n * s))


@dataclass
class CRFFit:
    rmax: float
    c50: float
    n: float
    s: float
    sse: float
    sse_s1: float
    f: float
    p: float          # F-test, s free vs s = 1 (one extra parameter)
    peak_contrast: float


def fit_crf(c: np.ndarray, r: np.ndarray) -> CRFFit:
    """Fit R(c) = Rmax c^n / (c^(n s) + c50^(n s)) to single presentations (blank-subtracted)."""
    c, r = np.asarray(c, float), np.asarray(r, float)
    ok = np.isfinite(r)
    c, r = c[ok], r[ok]
    top = max(np.nanmax(means_by_contrast(c, r)), 1e-6)

    def fit(free_s: bool):
        best = None
        for c50 in (0.05, 0.15, 0.4):
            for n in (1.5, 3.0):
                x0 = [top * 1.5, c50, n] + ([1.5] if free_s else [])
                lo = [0, 0.01, 0.5] + ([0.5] if free_s else [])
                hi = [top * 20, 1.5, 8] + ([6.0] if free_s else [])
                f = lambda p: naka_rushton(c, p[0], p[1], p[2], p[3] if free_s else 1.0) - r
                try:
                    res = optimize.least_squares(f, x0, bounds=(lo, hi))
                except ValueError:
                    continue
                if best is None or res.cost < best.cost:
                    best = res
        return best

    full, nested = fit(True), fit(False)
    sse, sse1 = 2 * full.cost, 2 * nested.cost
    dof = len(r) - 4
    f = max((sse1 - sse) / (sse / dof), 0.0) if dof > 0 and sse > 0 else np.nan
    p = float(stats.f.sf(f, 1, dof)) if np.isfinite(f) else np.nan
    rmax, c50, n, s = full.x
    grid = np.linspace(0.001, 0.9, 900)
    peak_c = float(grid[np.argmax(naka_rushton(grid, rmax, c50, n, s))])
    return CRFFit(float(rmax), float(c50), float(n), float(s), sse, sse1, float(f), p, peak_c)
