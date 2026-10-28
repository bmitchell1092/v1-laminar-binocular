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


def test_top_of_cortex_and_steps_follow_orientation():
    from v1lb.laminar import top_of_cortex, phase_steps, l45_from_steps
    rate = np.r_[np.full(4, 1.0), np.full(28, 20.0)]           # contacts 1-4 silent (above cortex)
    evoked = np.full(32, 20.0)                                   # volume-conducted: active everywhere
    assert top_of_cortex(rate, evoked, "ascending") == (5, False, "spontaneous MUA")
    assert top_of_cortex(rate[::-1], evoked, "descending") == (28, False, "spontaneous MUA")
    assert top_of_cortex(np.full(32, 20.0), evoked, "ascending") == (1, True, "spontaneous MUA")
    tiny = np.full(32, 1.0)                                      # spikes too small: fall back to evoked
    assert top_of_cortex(tiny, np.r_[np.zeros(6), np.full(26, 20.0)], "ascending") == (7, False, "evoked MUA")
    # two canonical steps below the top of cortex: +2.5 -> -2.0 at 8/9, -2.0 -> -0.6 at 23/24 (211008-like)
    prof = np.r_[np.full(8, 2.5), np.full(15, -2.0), np.full(9, -0.6)]
    steps = phase_steps(prof, "ascending", start_contact=3)
    assert [s[0] for s in steps] == [8.5, 23.5] and all(s[2] > 0 for s in steps)
    assert l45_from_steps(steps) == 23.5
    steps = phase_steps(prof[::-1], "descending", start_contact=30)
    assert [s[0] for s in steps] == [24.5, 9.5] and l45_from_steps(steps) == 9.5
    # a second fold entered from its deep side gives a reversed (negative) step
    fold = np.r_[np.full(10, 3.0), np.full(10, -1.2), np.full(12, 3.0)]
    steps = phase_steps(fold, "ascending")
    assert [s[2] > 0 for s in steps] == [True, False]
    assert l45_from_steps(steps) == 10.5
    # a ~pi step has no defined direction (211012-like): ambiguous, still an L4/5 candidate
    from v1lb.laminar import step_kind
    pi_step = np.r_[np.full(16, 2.6), np.full(16, -0.5)]
    steps = phase_steps(pi_step, "ascending")
    assert len(steps) == 1 and step_kind(steps[0][2]) == "ambiguous" and l45_from_steps(steps) == 16.5
    # steps among silent contacts below the deepest active contact are ignored
    assert phase_steps(prof, "ascending", start_contact=3, stop_contact=20) == phase_steps(prof, "ascending", 3)[:1]


def test_mua_gap():
    from v1lb.laminar import mua_gap
    rate = np.r_[np.full(15, 20.0), [1.0, 1.2, 0.9], np.full(14, 20.0)]
    assert mua_gap(rate, 10.5, 22.5) and not mua_gap(rate, 2.5, 12.5)


def test_initial_sink_extent():
    from v1lb.evp import initial_sink_extent
    t = np.arange(-50, 150).astype(float)
    c = np.zeros((t.size, 12))
    c[(t >= 35) & (t < 60), 5:9] = -4.0            # sink spanning channels 5-8
    c[(t >= 35) & (t < 60), 9] = -1.0              # weak tail: not part of the sink
    assert initial_sink_extent(c, t, 35.0, seed=6) == (5, 8)


def test_registration_recovers_shift():
    from v1lb.registration import register
    rng = np.random.default_rng(3)
    z = np.arange(48)
    tissue = {"log_rate": np.where((z > 8) & (z < 34), 1.2, 0.0) + 0.3 * np.sin(z / 3),
              "gamma": np.exp(-((z - 16) / 5) ** 2), "alpha_beta": np.exp(-((z - 28) / 6) ** 2)}
    ref = {k: v[8:40] + rng.normal(0, 0.02, 32) for k, v in tissue.items()}        # probe at final depth
    mov = {k: v[2:34] + rng.normal(0, 0.02, 32) for k, v in tissue.items()}        # 6 contacts shallower earlier
    r = register(mov, ref)
    assert r.shift == -6 and r.score > 0.9 and r.margin > 0.1


def test_dominant_early_sink_ignores_weak_early_blip():
    from v1lb.evp import dominant_early_sink
    t = np.arange(-50, 150).astype(float)
    c = np.zeros((t.size, 16))
    c[(t >= 22) & (t < 30), 2] = -1.0               # weak early blip (an onset detector picks this)
    c[(t >= 38) & (t < 70), 9:12] = -5.0            # dominant layer-4 sink
    peak, top, bottom, _ = dominant_early_sink(c, t)
    assert (top, bottom) == (9, 11) and peak in (9, 10, 11)
