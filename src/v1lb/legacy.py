"""Read legacy BM-ephys STIM files (MAT v7.3) for side-by-side validation.

Legacy conventions (see docs/01): `tp_sp` = event-code ticks, `tp_pt` =
"photodiode" ticks on a clock shifted by the .ns6 start (-102) and minus one
refresh; `contrast` = [fixedc contrast] = [left right]; `trl` = 1-based NEV
trial index (the legacy pipeline matched stimulus-file trials by index).
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def _char(ds) -> str:
    return "".join(map(chr, np.asarray(ds[()]).ravel()))


def _cellstr(f: h5py.File, ds) -> list[str]:
    return [_char(f[r]) for r in np.asarray(ds[()]).ravel()]


def read_legacy_stim(path: str | Path) -> tuple[pd.DataFrame, dict]:
    """Per-presentation table (one column per trial-length STIM field) and metadata."""
    with h5py.File(path, "r") as f:
        S = f["STIM"]
        n = S["trl"].shape[-1]
        cols = {}
        for key, ds in S.items():
            if not isinstance(ds, h5py.Dataset) or ds.dtype.kind == "O" or ds.attrs.get("MATLAB_class") == b"char":
                continue
            a = np.asarray(ds[()]).T
            if a.ndim != 2 or a.shape[0] != n:
                continue
            if a.shape[1] == 1:
                cols[key] = a[:, 0]
            else:
                for j in range(a.shape[1]):
                    cols[f"{key}_{j + 1}"] = a[:, j]
        meta = {}
        for key in ("penetration", "V1", "sortdirection", "header"):
            if key in S:
                meta[key] = _char(S[key])
        for key in ("el_labels", "filelist", "paradigm"):
            if key in S:
                meta[key] = _cellstr(f, S[key])
        for key in ("v1lim", "depths"):
            if key in S:
                meta[key] = np.asarray(S[key][()]).T
    df = pd.DataFrame(cols)
    for key in ("trl", "prez", "filen"):
        if key in df:
            df[key] = df[key].astype(int)
    return df, meta
