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
