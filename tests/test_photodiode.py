import numpy as np
import pandas as pd
import pytest

from _synth import write_nsx
from v1lb.io import NsxFile
from v1lb.photodiode import detect

DARK, WHITE, GRAY = 980, -334, 880       # mV, as measured on the rig
FRAME, PULSE = 250, 105                  # ticks


def trace(n_ticks, stimuli, gray_until=0):
    """Photodiode trace on ticks 0..n_ticks-1 with strobed flashes for each (first_flash, n_frames)."""
    x = np.full(n_ticks, DARK, float)
    for k in range(0, gray_until, FRAME):        # gray background pulses (fixation scene)
        x[k:k + PULSE] = GRAY
    for first, n in stimuli:
        for k in range(n):
            s = first + k * FRAME
            x[s:s + PULSE] = WHITE
    return x


def nsx_with_trace(tmp_path, x, t0=0):
    mv_per_bit = 10000 / 65534
    data = np.round(x / mv_per_bit).astype(np.int16)[:, None]
    return NsxFile(write_nsx(tmp_path / "pd.ns6", ["ainp1"], [(t0, data)],
                             analog=(-5000, 5000), digital=(-32767, 32767)))


def test_onset_frames_offset(tmp_path):
    codes_on = [20000, 40000, 60000]
    lat = [625, 375, 625]                         # 20.8 ms, 12.5 ms (one frame early), 20.8 ms
    frames = [38, 38, 26]
    x = trace(80000, [(c + l, n) for c, l, n in zip(codes_on, lat, frames)], gray_until=19000)
    nsx = nsx_with_trace(tmp_path, x)
    pres = pd.DataFrame({"on_tick": codes_on, "off_tick": [c + 9500 for c in codes_on]})
    res, lv = detect(nsx, pres)
    assert lv.polarity == -1 and lv.frame_ticks == pytest.approx(FRAME)
    assert res["pd_on_tick"].tolist() == [c + l for c, l in zip(codes_on, lat)]
    assert res["pd_n_frames"].tolist() == frames
    assert (res["pd_off_tick"] - res["pd_on_tick"]).tolist() == [n * FRAME for n in frames]
    assert res["pd_flag"].tolist() == ["", "", ""]


def test_missing_flash_is_flagged(tmp_path):
    x = trace(40000, [(20625, 38)])
    nsx = nsx_with_trace(tmp_path, x)
    pres = pd.DataFrame({"on_tick": [5000, 20000], "off_tick": [-1, 29500]})
    res, _ = detect(nsx, pres)
    assert res["pd_flag"].tolist() == ["no_flash", ""]
    assert np.isnan(res["pd_on_lat_ms"].iloc[0])
