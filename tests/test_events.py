import numpy as np

from _synth import write_nev
from v1lb.events import ml_codes, presentations, split_trials, trial_bracket_problems
from v1lb.io import NevFile

# recording starts mid-trial (no 9-triplet), then: a trial broken during stimulus 1,
# a complete 3-stimulus trial, a trial broken during stimulus 2.
CODES = ([35, 11, 8, 23, 24, 18, 18, 18]
         + [9, 9, 9, 117, 117, 117, 35, 11, 8, 23, 97, 36, 18, 18, 18]
         + [9, 9, 9, 117, 117, 118, 35, 11, 8, 23, 24, 25, 26, 27, 32, 36, 96, 96, 18, 18, 18]
         + [9, 9, 9, 117, 117, 119, 35, 11, 8, 23, 24, 25, 97, 36, 18, 18, 18])


def test_nev_roundtrip_and_codes(tmp_path):
    ticks = (np.arange(len(CODES)) * 300 + 1000).tolist()
    f = NevFile(write_nev(tmp_path / "x.nev", [(t, c + 128) for t, c in zip(ticks, CODES)],
                          spikes=[(1500, 3, 0), (2500, 5, 1)]))
    t, v, reason = f.digital_events()
    assert t.tolist() == ticks and ml_codes(v).tolist() == CODES and (reason == 1).all()
    st, eid, unit = f.spikes()
    assert st.tolist() == [1500, 2500] and eid.tolist() == [3, 5] and unit.tolist() == [0, 1]


def test_split_trials_skips_partial_first_trial():
    codes = np.array(CODES)
    ticks = np.arange(len(codes)) * 300
    tr = split_trials(codes, ticks)
    assert len(tr) == 3
    assert tr["trial_digit"].tolist() == [1, 2, 3]
    assert tr["rewarded"].tolist() == [False, True, False]
    assert tr["broke_fixation"].tolist() == [True, False, True]
    assert trial_bracket_problems(codes) == []


def test_presentations_completion():
    codes = np.array(CODES)
    ticks = np.arange(len(codes)) * 300
    p = presentations(codes, ticks, split_trials(codes, ticks))
    assert list(zip(p.nev_trial, p.pres, p.completed)) == [
        (0, 1, False),
        (1, 1, True), (1, 2, True), (1, 3, True),
        (2, 1, True), (2, 2, False)]
    third = p[(p.nev_trial == 1) & (p.pres == 3)].iloc[0]
    assert third.off_code == 32                        # bminteroc ends stimulus 3 with 32, not 28
