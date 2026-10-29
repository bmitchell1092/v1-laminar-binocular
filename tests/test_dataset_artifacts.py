import numpy as np

from v1lb.dataset import artifact_presentations


def test_artifact_flags_only_extreme_presentations():
    rng = np.random.default_rng(0)
    t = np.arange(-200, 500)
    mua = rng.normal(0, 1, (300, t.size, 4))
    mua[::3, 250:300, :] += 4.0                      # a visual response in a third of presentations
    mua[17, 400, 2] += 400.0                         # one artifact on one channel
    flagged = artifact_presentations(mua, t, np.array([True, True, True, True]))
    assert flagged.sum() == 1 and flagged[17]
    # the artifact is ignored if its channel is outside cortex
    assert not artifact_presentations(mua, t, np.array([True, True, False, True])).any()
