import numpy as np

from v1lb.dynorm import fit, predict

C = np.array([0.0, 0.055, 0.11, 0.225, 0.45, 0.9])


def grid():
    cd, cn = np.meshgrid(C, C, indexing="ij")
    keep = (cd + cn) > 0
    return cd[keep], cn[keep]


def test_parameter_recovery_and_model_comparison():
    t = np.arange(0, 220.0)
    cd, cn = grid()
    true = dict(g=2.0, sigma=0.12, n=2.2, delta=30.0, tau1=4.0, tau2=35.0, a=0.6, w=1.2)
    rng = np.random.default_rng(0)
    y = predict(true, t, cd, cn) + rng.normal(0, 0.02, (len(cd), len(t)))
    f = fit(t, cd, cn, y, "full", starts=4)
    for k in ("n", "delta", "tau2", "a", "w"):
        assert abs(f.params[k] - true[k]) / true[k] < 0.15, (k, f.params[k])
    fs = fit(t, cd, cn, y, "static", starts=4)
    fw = fit(t, cd, cn, y, "w=a", starts=4)
    assert f.sse < fs.sse and f.sse < fw.sse


def test_binocular_below_dominant_eye_needs_w_above_a():
    t = np.arange(0, 220.0)
    base = dict(g=2.0, sigma=0.1, n=2.0, delta=30.0, tau1=4.0, tau2=35.0)
    late = slice(150, 220)
    for a, w, below in ((0.6, 1.2, True), (0.6, 0.4, False)):
        r = predict(base | dict(a=a, w=w), t, np.array([0.9, 0.9]), np.array([0.0, 0.9]))
        assert (r[1, late].mean() < r[0, late].mean()) == below
    early = predict(base | dict(a=0.6, w=1.2), t, np.array([0.9, 0.9]), np.array([0.0, 0.9]))[:, 33:38]
    assert (early[1] > early[0]).all()                      # at onset the second eye adds


def test_contrast_dependent_pool_time_constant_is_identifiable():
    t = np.arange(0, 220.0)
    cd, cn = grid()
    true = dict(g=2.0, sigma=0.12, n=1.5, delta=25.0, tau1=4.0, tau2=60.0, a=0.8, w=1.1, beta=8.0)
    rng = np.random.default_rng(2)
    y = predict(true, t, cd, cn) + rng.normal(0, 0.02, (len(cd), len(t)))
    fc = fit(t, cd, cn, y, "tau_c", starts=3)
    ff = fit(t, cd, cn, y, "full", starts=3)
    assert fc.sse < ff.sse and 4 < fc.params["beta"] < 16


def test_two_stage_model_selected_when_true():
    from v1lb.dynorm import predict_two_stage
    t = np.arange(0, 220.0)
    cd, cn = grid()
    true = dict(g=1.5, sigma=0.3, n=1.3, delta=22.0, tau1=4.0, tau2=30.0, a=0.85, S=0.08, omega=1.2, tau_s=10.0)
    rng = np.random.default_rng(3)
    y = predict_two_stage(true, t, cd, cn) + rng.normal(0, 0.01, (len(cd), len(t)))
    f2 = fit(t, cd, cn, y, "2stage", starts=3)
    f1 = fit(t, cd, cn, y, "full", starts=3)
    assert f2.sse < f1.sse
    assert np.allclose(predict(f2.params, t, cd, cn), predict_two_stage(f2.params, t, cd, cn))


def test_full2_nests_single_stage_and_seeded_fit_is_never_worse():
    t = np.arange(0, 220.0)
    cd, cn = grid()
    p1 = dict(g=2.0, sigma=0.1, n=1.3, delta=22.0, tau1=4.0, tau2=30.0, a=0.8, w=1.1)
    S = 400.0
    p2 = p1 | dict(S=S, sigma=p1["sigma"] / S, omega=1.0, tau_s=20.0)
    r1, r2 = predict(p1, t, cd, cn), predict(p2, t, cd, cn)
    assert np.nanmax(np.abs(r1 - r2)) < 0.02 * np.nanmax(r1)          # large-S limit = single-stage 'full'
    rng = np.random.default_rng(4)
    y = r1 + rng.normal(0, 0.02, r1.shape)
    f1 = fit(t, cd, cn, y, "full", starts=2)
    init = f1.params | dict(S=S, sigma=f1.params["sigma"] / S, omega=1.0, tau_s=20.0)
    f2 = fit(t, cd, cn, y, "full2", starts=1, init=init)
    assert f2.sse <= f1.sse * 1.01


def test_ding_sperling_identities_and_averaging():
    from v1lb.dynorm import predict_ds
    t = np.arange(0, 400.0)
    p = dict(g=1.0, sigma=1e-6, n=1.0, delta=20.0, tau1=4.0, tau2=1.0, a=1.0, rho=10.0, gamma=1.2, tau_g=50.0)
    # with a near-identity output stage (sigma tiny, n = 1, fast pool) R ~ Gamma / lowpass(Gamma) -> use Gamma directly
    from v1lb.dynorm import step_responses, lowpass
    u, _ = step_responses(t, 20.0, 4.0, 1.0, static=True)
    x = 0.45 * u
    e = (10.0 * lowpass(x[None], 50.0)[0]) ** 1.2
    g_bin = (x * (1 + e) + x * (1 + e)) / (1 + 2 * e)
    assert np.allclose(g_bin[300:], x[300:] * (2 + 2 * e[300:]) / (1 + 2 * e[300:]))
    # late, high contrast, unequal eyes: binocular below the dominant eye alone; at onset it adds
    q = dict(g=2.0, sigma=0.1, n=1.5, delta=20.0, tau1=4.0, tau2=30.0, a=0.6, rho=20.0, gamma=1.5, tau_g=40.0)
    r = predict_ds(q, t, np.array([0.45, 0.45]), np.array([0.0, 0.45]))
    assert r[1, 30:34].mean() > r[0, 30:34].mean() and r[1, 250:].mean() < r[0, 250:].mean()


def test_interaction_loss_runs_and_prefers_truth():
    from v1lb.dynorm import loss_value
    t = np.arange(0, 220.0)
    cd, cn = grid()
    true = dict(g=2.0, sigma=0.1, n=1.5, delta=20.0, tau1=4.0, tau2=30.0, a=0.6, rho=20.0, gamma=1.5, tau_g=40.0)
    rng = np.random.default_rng(5)
    y = predict(true, t, cd, cn) + rng.normal(0, 0.01, (len(cd), len(t)))
    f = fit(t, cd, cn, y, "ds", starts=2, loss="interaction")
    assert loss_value(f.params, t, cd, cn, y) <= loss_value(true, t, cd, cn, y) * 1.05


def test_org_steady_state_and_contrast_dependent_rise():
    from v1lb.dynorm import predict_org
    t = np.arange(0, 400.0)
    p = dict(g=1.0, b=1.0, sigma=0.1, n=2.0, delta=10.0, tau1=2.0, tau_y=2.0, tau_a=5.0, a=0.0, w=0.0)
    r = predict_org(p, t, np.array([0.03, 0.3]), np.array([0.0, 0.0]))
    z = np.array([0.03, 0.3])
    assert np.allclose(r[:, -1], z ** 2 / (0.01 + z ** 2), rtol=1e-3)      # Naka-Rushton steady state
    half = [t[np.flatnonzero(x >= 0.5 * x[-1])[0]] for x in r]
    assert half[0] > half[1] + 5                                           # weak input rises later and slower


def test_cond_nests_full_and_slows_weak_inputs():
    from v1lb.dynorm import predict_cond
    t = np.arange(0, 220.0)
    cd, cn = grid()
    p = dict(g=2.0, sigma=0.12, n=2.0, delta=30.0, tau1=4.0, tau2=35.0, a=0.6, w=1.2)
    assert np.allclose(predict_cond(p | {"tau_m": 0.05}, t, cd, cn), predict(p, t, cd, cn), atol=1e-6)
    assert np.allclose(predict(p | {"tau_m": 0.05}, t, cd, cn), predict(p, t, cd, cn), atol=1e-6)
    c = np.array([0.055, 0.9])
    slow = predict_cond(p | {"tau_m": 20.0}, t, c, np.zeros(2))
    fast = predict(p, t, c, np.zeros(2))
    lag = [t[np.flatnonzero(s >= 0.5 * s.max())[0]] - t[np.flatnonzero(f >= 0.5 * f.max())[0]]
           for s, f in zip(slow, fast)]
    assert lag[0] > lag[1] + 3                                             # integration slows weak inputs more
