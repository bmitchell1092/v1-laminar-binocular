"""Checks against the real recordings and an independent reader (neo). Marked `raw`."""

import numpy as np
import pytest

from v1lb.config import recording_files
from v1lb.events import ml_codes, split_trials
from v1lb.io import NevFile, NsxFile, read_stim_text
from v1lb.trials import match_trials

pytestmark = pytest.mark.raw


def _neo_chunk(base, stream_name, i_start, i_stop, labels):
    from neo.rawio import BlackrockRawIO
    r = BlackrockRawIO(filename=str(base), nsx_to_load=[2, 6])
    r.parse_header()
    streams = r.header["signal_streams"]
    si = [s["name"] for s in streams].index(stream_name)
    chans = r.header["signal_channels"]
    names = [c["name"] for c in chans if c["stream_id"] == streams[si]["id"]]
    raw = r.get_analogsignal_chunk(0, 0, i_start, i_stop, stream_index=si)
    return raw[:, [names.index(l) for l in labels]], r.get_signal_t_start(0, 0, si)


@pytest.mark.parametrize("ext,stream,t0,period", [(".ns6", "nsx6", 102, 1), (".ns2", "nsx2", 0, 30)])
def test_nsx_matches_neo_on_nsp_clock(ext, stream, t0, period):
    files = recording_files("211008_B", "bminteroc003")
    f = NsxFile(files[ext])
    assert f.t_start == t0 and f.period == period
    labels = ["eA01", "eA17", "ainp1"]
    i0, n = 30000, 600
    ours, ticks = f.read_raw(t0 + i0 * period, t0 + (i0 + n) * period, labels)
    theirs, neo_t0 = _neo_chunk(files[ext].with_suffix(""), stream, i0, i0 + n, labels)
    assert neo_t0 == pytest.approx(t0 / 30000)
    np.testing.assert_array_equal(ours, theirs)


def test_211008_trials_and_text_align_by_index():
    files = recording_files("211008_B", "bminteroc003")
    nev = NevFile(files[".nev"])
    t, v, _ = nev.digital_events()
    codes = ml_codes(v)
    assert len(codes) == 2773
    tr = split_trials(codes, t)
    assert len(tr) == 182
    m = match_trials(tr, nev.origin_utc, read_stim_text(files[".gBMINTEROCGrating_di"]))
    assert m.shift == 0 and round(m.offset_s / 3600) == -5      # CDT


def test_211009_blackrock_missed_ml_trial_1():
    files = recording_files("211009_B", "bminteroc001")
    nev = NevFile(files[".nev"])
    t, v, _ = nev.digital_events()
    tr = split_trials(ml_codes(v), t)
    m = match_trials(tr, nev.origin_utc, read_stim_text(files[".gBMINTEROCGrating_di"]))
    assert m.shift == 1 and m.max_dev_s < 0.1
    assert m.table["ml_trial"].iloc[0] == 2


def test_photodiode_211008():
    from v1lb.events import presentations
    from v1lb.photodiode import detect
    files = recording_files("211008_B", "bminteroc003")
    nev = NevFile(files[".nev"])
    t, v, _ = nev.digital_events()
    codes = ml_codes(v)
    p = presentations(codes, t, split_trials(codes, t))
    res, lv = detect(NsxFile(files[".ns6"]), p)
    done = res[res.completed]
    assert len(done) == 398 and (done.pd_flag == "").all()
    assert lv.frame_ticks == 250 and lv.polarity == -1
    assert done.pd_n_frames.mode()[0] == 26                       # 216.7 ms stimuli on 211008
    assert abs(done.pd_on_lat_ms.median() - 20.9) < 0.2
