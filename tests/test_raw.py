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


def _legacy_vs_new(pen, rec):
    from v1lb.bminteroc import presentation_table
    from v1lb.legacy import read_legacy_stim
    from v1lb.config import raw_root
    L, meta = read_legacy_stim(raw_root() / "STIM 2024 2" / f"{pen}.mat")
    P, m = presentation_table(pen[:8], rec)
    P = P[P.completed].copy()
    P["trl"] = P.nev_trial + 1
    M = L.merge(P, left_on=["trl", "prez"], right_on=["trl", "pres"], how="outer", indicator=True)
    return M, m, meta


def test_bminteroc_table_reproduces_legacy_211008():
    M, m, meta = _legacy_vs_new("211008_B_eA", "bminteroc003")
    assert (M._merge == "both").all() and len(M) == 398
    assert (M.tp_sp_1 == M.on_tick).all()                                  # same event codes
    assert np.allclose(M.contrast_1, M.contrast_L) and np.allclose(M.contrast_2, M.contrast_R)
    assert meta["V1"] == m["hemisphere"] == "RV1"
    d_ms = (M.tp_pt_1 - M.pd_on_tick) / 30                                # legacy t0 ~12 ms early
    assert abs(np.median(d_ms) + 12.1) < 0.5


def test_bminteroc_table_fixes_211009_labels():
    M, m, _ = _legacy_vs_new("211009_B_eA", "bminteroc001")
    assert (M._merge == "both").all() and m["trial_shift"] == 1
    agree = np.isclose(M.contrast_1, M.contrast_L) & np.isclose(M.contrast_2, M.contrast_R)
    assert agree.mean() < 0.05                                             # legacy labels were wrong


def test_dataset_layout_and_depths():
    from v1lb.dataset import build_penetration
    ds = build_penetration("211008_B_eA", max_presentations=6, write=False, alignment_file="laminar_alignment_v1.csv")
    assert ds.lfp.dims == ("presentation", "time", "channel") and ds.lfp.shape == (6, 700, 32)
    assert ds.time.values[0] == -200 and ds.time.values[-1] == 499
    assert ds.channel.values[0] == "eA01" and ds.depth_l45.values[0] == 22.5      # ascending, L4/5 at 23.5
    assert ds.depth_l45.values[-1] == -8.5 and not ds.in_cortex.values[1] and ds.in_cortex.values[2]
    assert ds.depth_below_toc.values[2] == 0 and ds.depth_below_toc.values[-1] == 29     # top of cortex = eA03
    assert ds.bad_lfp.dtype == bool
    assert np.isnan(ds.csd.values[:, :, 0]).all() and np.isfinite(ds.csd.values[:, :, 1]).all()
    assert ds.spikes.sum() > 0 and ds.attrs["orientation"] == "ascending"
    assert set(ds.condition.values) <= {"blank", "mono_contra", "mono_ipsi", "dioptic", "dichoptic"}
    fix = ds.eye_x.sel(time=slice(-50, 0)).values[ds.completed.values]          # fixating before onset
    assert ds.eye_x.attrs["source"] == "ainp3" and np.nanmedian(np.abs(fix)) < 3   # near the fixation point (deg)
    assert np.isfinite(ds.pupil.values).all()


def test_dataset_descending_probe_is_reordered():
    from v1lb.dataset import build_penetration
    ds = build_penetration("211012_B_eA", max_presentations=3, write=False, alignment_file="laminar_alignment_v1.csv")
    assert ds.channel.values[0] == "eA32" and ds.channel.values[-1] == "eA01"          # superficial first
    assert ds.depth_l45.values[0] == 6.5 and (np.diff(ds.depth_l45.values) == -1).all()


def test_eye_channels_identified_and_calibrated():
    from v1lb.eye import calibrate
    c = calibrate("211008_B", "bminteroc003")
    assert (c.x_channel, c.y_channel) == ("ainp3", "ainp4") and min(c.r2) > 0.95


def test_kilosort_offset_two_sort_generations():
    from v1lb.config import recording_files
    from v1lb.io import NsxFile
    from v1lb.units import load_sort, measure_offset
    for session, rec, expected in [("211008_B", "bminteroc003", -204.5), ("221206_J", "bminteroc001", 0)]:
        rep = measure_offset(load_sort(session, rec), NsxFile(recording_files(session, rec)[".ns6"]), n_units=3)
        assert rep["consistent"] and abs(rep["offset"] - expected) <= 1.5
