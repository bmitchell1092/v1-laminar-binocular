import numpy as np
import pandas as pd

from v1lb.crf import (CONTRASTS, bootstrap_si, crf_points, fit_crf, means_by_contrast, naka_rushton,
                      supersaturation_index)


def _trials(params, n=30, noise=1.0, seed=0):
    rng = np.random.default_rng(seed)
    c = np.repeat(CONTRASTS, n)
    return c, naka_rushton(c, *params) + rng.normal(0, noise, c.size)


def test_si_and_fit_separate_saturating_from_supersaturating():
    c, r = _trials((20.0, 0.1, 2.5, 1.0))                     # saturating
    m = means_by_contrast(c, r)
    assert supersaturation_index(m) < 0.1
    f = fit_crf(c, r)
    assert f.p > 0.01 and f.peak_contrast > 0.8
    c, r = _trials((20.0, 0.1, 2.5, 2.0))                     # supersaturating: peaks near 0.1, falls by 0.9
    m = means_by_contrast(c, r)
    assert supersaturation_index(m) > 0.3
    f = fit_crf(c, r)
    assert f.p < 0.001 and 1.6 < f.s < 2.5 and f.peak_contrast < 0.3
    b = bootstrap_si(c, r, n_boot=200)
    assert np.percentile(b, 2.5) > 0.1


def test_crf_points_select_conditions():
    P = pd.DataFrame(dict(contrast_contra=[0, .11, 0, .45, .9, .22], contrast_ipsi=[0, 0, .11, .45, .11, .22]))
    resp = np.arange(6.0)
    assert crf_points(resp, P, "mono_contra")[1].tolist() == [0, 1]
    assert crf_points(resp, P, "mono_ipsi")[1].tolist() == [0, 2]
    assert crf_points(resp, P, "dioptic")[1].tolist() == [0, 3, 5]


def test_cv_si_is_unbiased_for_monotonic_crfs_and_detects_downturns():
    from v1lb.crf import cv_supersaturation_index
    sites_sat, sites_sup = [], []
    for seed in range(40):                                           # noisy, saturating: plain SI > 0, CV SI ~ <= 0
        c, r = _trials((10.0, 0.1, 2.5, 1.0), noise=4.0, seed=seed)
        sites_sat.append(r)
        c, r = _trials((10.0, 0.1, 2.5, 2.0), noise=4.0, seed=100 + seed)
        sites_sup.append(r)
    R_sat, R_sup = np.stack(sites_sat, axis=1), np.stack(sites_sup, axis=1)
    plain = supersaturation_index(means_by_contrast(c, R_sat))
    cv_sat = cv_supersaturation_index(c, R_sat, n_splits=100)
    cv_sup = cv_supersaturation_index(c, R_sup, n_splits=100)
    assert np.nanmean(plain) > 0.03                                  # the bias the CV index removes
    assert np.nanmean(cv_sat) < 0.01 and np.nanmean(cv_sup) > 0.2


def test_downturn_p_is_calibrated_and_sensitive():
    from v1lb.crf import downturn_p
    sat, sup = [], []
    for seed in range(60):
        c, r = _trials((10.0, 0.1, 2.5, 1.0), noise=4.0, seed=seed)       # saturating (null)
        sat.append(r)
        c, r = _trials((10.0, 0.1, 2.5, 2.0), noise=4.0, seed=500 + seed)  # supersaturating
        sup.append(r)
    p_null = downturn_p(c, np.stack(sat, axis=1), n_splits=30)
    p_alt = downturn_p(c, np.stack(sup, axis=1), n_splits=30)
    assert np.mean(p_null < 0.05) <= 0.05 and np.mean(p_alt < 0.05) > 0.8
