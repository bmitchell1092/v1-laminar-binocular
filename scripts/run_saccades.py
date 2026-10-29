"""Saccades for every built dataset: <pen>_saccades.parquet (continuous) and <pen>_presentation_eye.parquet
(per presentation: saccades / microsaccades in -200...500 ms, time since the last saccade, to the next).

Usage: .venv\Scripts\python scripts\run_saccades.py [penetration ...]
"""

import sys

import pandas as pd
import xarray as xr

from v1lb.config import out_root
from v1lb.saccades import around_onsets, detect_continuous

if __name__ == "__main__":
    d = out_root() / "datasets"
    pens = [a for a in sys.argv[1:] if not a.startswith("--")] or sorted(
        p.name.removesuffix("_eye.nc") for p in d.glob("*_eye.nc"))
    for pen in pens:
        e = xr.open_dataset(d / f"{pen}_eye.nc", engine="h5netcdf")
        s = detect_continuous(e.tick.values, e.eye_x.values.astype(float), e.eye_y.values.astype(float),
                              e.pupil.values.astype(float))
        s.to_parquet(d / f"{pen}_saccades.parquet")
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        a = around_onsets(s, P.pd_on_tick.to_numpy())
        a.insert(0, "presentation", range(len(P)))
        a.to_parquet(d / f"{pen}_presentation_eye.parquet")
        minutes = (e.tick.values[-1] - e.tick.values[0]) / 30000 / 60
        print(f"{pen}: {len(s)} saccades ({len(s) / minutes / 60:.2f}/s, {(s.amplitude < 1).mean():.0%} micro); "
              f"presentations with any saccade in -200..500 ms: {(a.n_saccades > 0).mean():.0%}, "
              f"with a large one: {(a.n_large > 0).mean():.0%}", flush=True)
