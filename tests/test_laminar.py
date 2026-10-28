import numpy as np
import pytest

from v1lb.laminar import analytic, csd, earliest_sink, generalized_phase, spike_phase

FS = 1000.0


def test_ppc_matches_pairwise_definition():
    rng = np.random.default_rng(0)
    ph = rng.vonmises(0.5, 2.0, size=300)
    phase = ph[:, None]
    ticks = np.arange(300) * 30
    sp = spike_phase([ticks], phase, ticks)
    d = ph[:, None] - ph[None, :]
    iu = np.triu_indices(300, 1)
    assert sp.ppc[0, 0] == pytest.approx(np.cos(d[iu]).mean(), rel=1e-9)
    assert sp.theta[0, 0] == pytest.approx(np.angle(np.exp(1j * ph).sum()))


def test_spike_locking_recovered_per_lfp_channel():
    t = np.arange(60000) / FS
    lfp = np.stack([np.sin(2 * np.pi * 10 * t), np.sin(2 * np.pi * 10 * t + np.pi)], axis=1)
    phase = np.angle(analytic(lfp))
    ticks = np.arange(60000) * 30
    locked = ticks[np.flatnonzero(np.isclose(np.angle(np.exp(1j * phase[:, 0])), 0.0, atol=0.05))]
    sp = spike_phase([locked], phase, ticks)
    assert sp.spi[0, 0] > 0.95
    assert abs(np.angle(np.exp(1j * (sp.theta[0, 1] - np.pi)))) < 0.1   # reversed channel: opposite phase


def test_generalized_phase_equals_hilbert_for_clean_oscillation():
    t = np.arange(20000) / FS
    x = np.sin(2 * np.pi * 12 * t)[:, None]
    gp = generalized_phase(x, FS, lp=5)
    hp = np.angle(analytic(x))
    core = slice(500, -500)
    assert np.allclose(np.angle(np.exp(1j * (gp[core] - hp[core]))), 0, atol=1e-6)


def test_csd_of_quadratic_potential_is_constant():
    z = np.arange(10) * 0.1                      # mm
    v = 50 * z ** 2                              # uV, d2V/dz2 = 100 uV/mm^2
    c = csd(v[None, :])
    assert np.isnan(c[0, 0]) and np.isnan(c[0, -1])
    assert np.allclose(c[0, 1:-1], -0.4 * 100)   # nA/mm^3


def test_earliest_sink():
    t = np.arange(-50, 150).astype(float)
    x = np.random.default_rng(1).normal(0, 0.1, (t.size, 5))
    x[t >= 35, 2] -= 5                            # channel 2 sinks first
    x[t >= 50, 3] -= 5
    x[t >= 25, 4] -= 1                            # earlier but weak: not a qualifying sink
    onset, first = earliest_sink(x, t)
    assert first == 2 and onset[2] == 35 and onset[3] == 50 and np.isnan(onset[4])


def test_template_fit_recovers_orientation_and_boundary():
    from v1lb.laminar import davis_templates, fit_template
    T = davis_templates()
    tmpl, w = T["phase"], (T["n_phase"] >= 10).astype(float)
    # an ascending probe: contact 1 sits at template index 18 -> boundary (32.5) at contact 15.5
    prof = tmpl[18:50].copy()
    fit = fit_template(prof, tmpl, "phase", w)
    assert fit.orientation == "ascending" and fit.boundary == pytest.approx(15.5) and fit.margin > 0.2
    # the same probe wired the other way round
    fit = fit_template(prof[::-1], tmpl, "phase", w)
    assert fit.orientation == "descending" and fit.boundary == pytest.approx(32 - 15.5 + 1)
    # boundary near the deep end of the probe (contact 27.5): needs the extended plateau
    padded = np.r_[np.full(20, tmpl[14]), tmpl[14:52]]   # template index i -> padded index i + 6
    prof = padded[12:44]                                  # boundary (padded 38.5) at contact 27.5
    fit = fit_template(prof, tmpl, "phase", w)
    # only 5 contacts below the boundary, so orientation confidence (margin) is lower but still positive
    assert fit.orientation == "ascending" and fit.boundary == pytest.approx(27.5) and fit.margin > 0.1
    # spectrolaminar template, polarity-free
    g = T["gamma_minus_alpha"]
    fit = fit_template(g[10:42][::-1], g, "linear")
    assert fit.orientation == "descending"


def test_spectral_crossing():
    from v1lb.laminar import spectral_crossing
    z = np.arange(32)
    gamma = np.exp(-((z - 10) / 5) ** 2)
    alpha = np.exp(-((z - 26) / 5) ** 2)
    c = spectral_crossing(alpha, gamma)
    assert c.orientation == "ascending" and c.gamma_peak == 11 and c.alpha_peak == 27
    assert c.crossing == pytest.approx(19.0, abs=0.6)
    assert spectral_crossing(alpha[::-1], gamma[::-1]).orientation == "descending"
