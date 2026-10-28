import datetime as dt

import pytest

from v1lb.io import read_stim_text
from v1lb.io.stimtext import datenum_to_datetime


def test_datenum_conversion():
    # MATLAB: datestr(738437.768098) -> '08-Oct-2021 18:26:03'
    d = datenum_to_datetime(738437.768098)
    assert d.replace(microsecond=0) == dt.datetime(2021, 10, 8, 18, 26, 3)


def test_read_adds_trial_and_presentation(tmp_path):
    p = tmp_path / "x.gBMINTEROCGrating_di"
    rows = ["trial\tgrating_contrast\tgrating_fixedc\ttimestamp"]
    for tr in (1, 2):
        for c in (0.0, 0.45, 0.9):
            rows.append(f"{tr:04d}\t{c}\t0.11\t{738437.768098 + tr / 86400}")
    p.write_text("\r\n".join(rows) + "\r\n")
    df = read_stim_text(p)
    assert df["trial"].tolist() == [1, 1, 1, 2, 2, 2]
    assert df["pres"].tolist() == [1, 2, 3, 1, 2, 3]
    dt_s = (df["ml_time_local"].iloc[3] - df["ml_time_local"].iloc[0]).total_seconds()
    assert dt_s == pytest.approx(1.0, abs=1e-3)
