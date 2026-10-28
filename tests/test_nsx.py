import numpy as np
import pytest

from _synth import write_nsx
from v1lb.io import NsxFile


def ramp(n, nch):
    return (np.arange(n * nch, dtype=np.int16) % 1000).reshape(n, nch)


def test_header_and_clock(tmp_path):
    f = NsxFile(write_nsx(tmp_path / "a.ns6", ["eA01", "ainp1"], [(102, ramp(50, 2))]))
    assert f.fs == 30000 and f.t_start == 102 and f.t_stop == 152
    assert f.labels == ["eA01", "ainp1"]
    assert f.channels[0].scale == pytest.approx(16382 / 65528)  # 0.25 uV/bit
    assert f.origin_utc.isoformat() == "2021-10-08T23:25:45+00:00"


def test_read_is_on_nsp_clock(tmp_path):
    data = ramp(50, 2)
    f = NsxFile(write_nsx(tmp_path / "a.ns6", ["eA01", "ainp1"], [(102, data)]))
    x, t = f.read(100, 106, ["ainp1"])
    assert t.tolist() == [100, 101, 102, 103, 104, 105]
    assert np.isnan(x[:2]).all()                       # before the first sample
    np.testing.assert_allclose(x[2:, 0], data[:4, 1] * f.channels[1].scale)
    raw, _ = f.read_raw(102, 106, ["ainp1"])
    np.testing.assert_array_equal(raw[:, 0], data[:4, 1])


def test_gap_between_packets_is_nan(tmp_path):
    a, b = ramp(10, 1), ramp(10, 1) + 7
    f = NsxFile(write_nsx(tmp_path / "p.ns6", ["eA01"], [(102, a), (122, b)]))
    x, t = f.read(102, 132)
    assert np.isnan(x[10:20]).all()                    # ticks 112..121 not recorded
    np.testing.assert_allclose(x[20:, 0], b[:, 0] * f.channels[0].scale)
    with pytest.raises(ValueError):
        f.read_raw(108, 125)


def test_slow_stream_grid(tmp_path):
    f = NsxFile(write_nsx(tmp_path / "a.ns2", ["eA01"], [(0, ramp(20, 1))], period=30))
    assert f.fs == 1000
    _, t = f.read(15, 300)
    assert t.tolist() == list(range(30, 300, 30))       # samples sit on multiples of 30 ticks
