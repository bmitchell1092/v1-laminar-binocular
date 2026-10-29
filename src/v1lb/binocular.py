"""Binocular combination measures that do not depend on the output nonlinearity or on response latency.

Eye weight a: the contrast scale that maps the non-dominant eye's monocular CRF onto the dominant eye's,
R_N(c) = R_D(a c) (0 < a <= 1). Under linear summation of eye-weighted contrasts followed by any monotonic
output nonlinearity, the binocular response to (c, c) equals R_D((1 + a) c).

Equivalent contrast k(t): at each time t, the factor by which the dominant eye's monocular contrast must be raised
to match the binocular response, R_D(k c, t) = R_bin(c, c, t). Because it compares responses at matched level
rather than at matched time, a latency difference between conditions (higher effective contrast starts earlier)
does not change k: if the binocular response is the monocular response at contrast k0 c, latency included,
k(t) = k0 at every t.

Linearity index L = (k - 1) / a: 1 for linear summation of the two eyes, 0 for ocularity invariance (binocular =
dominant eye alone); binocular normalization with exponent n gives 2^(1/n) - 1 (about 0.3 for n = 2.5) for a
balanced site (derived; report Q2).
"""

from __future__ import annotations

import numpy as np

CONTRASTS = np.array([0.0, 0.055, 0.11, 0.225, 0.45, 0.9])


def _monotone(r: np.ndarray) -> np.ndarray:
    """Running maximum along the contrast axis (axis 0): a CRF is non-decreasing up to noise."""
    return np.maximum.accumulate(r, axis=0)


def invert_crf(r_mono: np.ndarray, target: np.ndarray, contrasts: np.ndarray = CONTRASTS) -> tuple[np.ndarray, np.ndarray]:
    """Contrast at which a monocular CRF reaches `target`.

    r_mono [contrast, ...] (blank-subtracted, contrasts[0] = 0), target [...]. Interpolates linearly between 0 and
    the lowest contrast and in log contrast above it, on the running-maximum envelope. Returns (contrast, censored),
    censored = +1 where target exceeds the response at the top contrast (contrast set to NaN), -1 where target <= 0.
    """
    r = _monotone(np.asarray(r_mono, float))
    target = np.asarray(target, float)
    out = np.full(target.shape, np.nan)
    cens = np.zeros(target.shape, int)
    lc = np.log(contrasts[1:])
    flat_r = r.reshape(len(contrasts), -1)
    flat_t = target.reshape(-1)
    flat_o = out.reshape(-1)
    flat_c = cens.reshape(-1)
    for i in range(flat_t.size):
        ri, ti = flat_r[:, i], flat_t[i]
        if not np.isfinite(ti) or not np.all(np.isfinite(ri)):
            continue
        if ti <= 0:
            flat_c[i] = -1
            continue
        if ti > ri[-1]:
            flat_c[i] = 1
            continue
        if ti <= ri[1]:                                         # between blank (0) and the lowest contrast
            flat_o[i] = contrasts[1] * ti / ri[1] if ri[1] > 0 else np.nan
            continue
        j = np.searchsorted(ri[1:], ti)                          # ri[1:][j-1] < ti <= ri[1:][j]
        r0, r1 = ri[1:][j - 1], ri[1:][j]
        f = (ti - r0) / (r1 - r0) if r1 > r0 else 0.0
        flat_o[i] = np.exp(lc[j - 1] + f * (lc[j] - lc[j - 1]))
    return out, cens


def eye_weight(r_dom: np.ndarray, r_non: np.ndarray, contrasts: np.ndarray = CONTRASTS,
               grid: np.ndarray | None = None) -> float:
    """Contrast scale a in (0, 1] minimising sum_c (R_N(c) - R_D(a c))^2 over the non-blank contrasts.

    r_dom, r_non: monocular CRFs [contrast] (blank-subtracted). R_D(a c) is interpolated as in `invert_crf`."""
    grid = np.geomspace(0.05, 1.0, 120) if grid is None else grid
    rd = _monotone(np.asarray(r_dom, float))
    rn = np.asarray(r_non, float)
    lc = np.log(contrasts[1:])

    def r_at(c):
        if c <= contrasts[1]:
            return rd[1] * c / contrasts[1]
        return np.interp(np.log(c), lc, rd[1:])

    sse = [np.nansum([(rn[k] - r_at(a * contrasts[k])) ** 2 for k in range(1, len(contrasts))]) for a in grid]
    return float(grid[int(np.argmin(sse))])


def equivalent_contrast(r_dom_t: np.ndarray, r_bin_t: np.ndarray, c: float,
                        contrasts: np.ndarray = CONTRASTS) -> tuple[np.ndarray, np.ndarray]:
    """k(t) = (contrast at which R_D(., t) equals R_bin(c, c, t)) / c, per time point.

    r_dom_t [contrast, time] monocular dominant-eye responses; r_bin_t [time] binocular (c, c) response.
    Returns (k [time], censored [time]); censored +1 means k > contrasts[-1] / c (top-censored)."""
    cstar, cens = invert_crf(r_dom_t, r_bin_t, contrasts)
    return cstar / c, cens
