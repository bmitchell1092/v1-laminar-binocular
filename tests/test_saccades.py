import numpy as np

from v1lb.saccades import detect


def _step(n, at, amp, dur=20):
    """Smooth position step (raised cosine) of `amp` deg starting at sample `at`, lasting `dur` samples."""
    p = np.zeros(n)
    k = np.arange(dur)
    p[at:at + dur] = amp * (1 - np.cos(np.pi * k / dur)) / 2
    p[at + dur:] = amp
    return p


def test_detects_saccade_and_microsaccade_and_ignores_drift():
    rng = np.random.default_rng(0)
    n = 4000
    x = _step(n, 1000, 5.0) + _step(n, 2500, 0.3, dur=12) + np.cumsum(rng.normal(0, 0.002, n))   # + drift
    y = np.cumsum(rng.normal(0, 0.002, n))
    x += rng.normal(0, 0.005, n); y += rng.normal(0, 0.005, n)
    s = detect(x, y)
    assert len(s) == 2
    assert abs(s.onset.iloc[0] - 1000) <= 5 and 4.5 < s.amplitude.iloc[0] < 5.5 and s.dx.iloc[0] > 4.5
    assert abs(s.onset.iloc[1] - 2500) <= 5 and 0.2 < s.amplitude.iloc[1] < 0.4


def test_invalid_samples_are_never_saccades():
    rng = np.random.default_rng(1)
    n = 3000
    x = rng.normal(0, 0.005, n); y = rng.normal(0, 0.005, n)
    x[1200:1350] = 40.0                                      # blink artefact
    valid = np.ones(n, bool); valid[1190:1360] = False
    assert len(detect(x, y, valid=valid)) == 0


def test_around_onsets():
    import pandas as pd
    from v1lb.saccades import around_onsets
    s = pd.DataFrame(dict(onset_tick=[30000, 60000, 90000], offset_tick=[30900, 60600, 90900], amplitude=[5.0, 0.4, 3.0]))
    a = around_onsets(s, [45000, 59000], window_ms=(-100, 300))
    assert a.n_saccades.tolist() == [0, 1] and a.n_micro.tolist() == [0, 1]
    assert a.since_saccade_ms.iloc[0] == (45000 - 30900) / 30 and a.to_saccade_ms.iloc[1] == 1000 / 30
