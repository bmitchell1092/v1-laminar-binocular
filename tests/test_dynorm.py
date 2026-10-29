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
