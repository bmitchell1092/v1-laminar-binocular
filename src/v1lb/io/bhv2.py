"""MonkeyLogic 2 `.bhv2` reader (read-only).

Layout (written recursively by MonkeyLogic's mlbhv2): every variable is
  uint64 name length, name (1 byte/char), uint64 type length, type (1 byte/char),
  uint64 ndims, uint64 dims[ndims], then the content:
  * struct: uint64 number of fields, then for each element (column-major) each field as a variable
  * cell: each element as a variable (empty name)
  * char: 1 byte per character; logical: 1 byte per value; numeric: raw little-endian values
  * function_handle: header 1x1, content = one nested (unnamed) char variable holding func2str
Top-level variables include MLConfig, TrialRecord and Trial1, Trial2, ... (AnalogData.Eye in degrees).
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

_NUMERIC = {"double": "<f8", "single": "<f4", "int8": "<i1", "uint8": "<u1", "int16": "<i2", "uint16": "<u2",
            "int32": "<i4", "uint32": "<u4", "int64": "<i8", "uint64": "<u8", "logical": "<u1"}


class _Reader:
    def __init__(self, data: bytes):
        self.b = data
        self.p = 0

    def u64(self) -> int:
        v = struct.unpack_from("<Q", self.b, self.p)[0]
        self.p += 8
        return v

    def text(self, n: int) -> str:
        s = self.b[self.p:self.p + n].decode("latin-1")
        self.p += n
        return s

    def variable(self):
        name = self.text(self.u64())
        typ = self.text(self.u64())
        dims = tuple(self.u64() for _ in range(self.u64()))
        count = int(np.prod(dims)) if dims else 0
        if typ == "struct":
            nfield = self.u64()
            elems = []
            for _ in range(count):
                e = {}
                for _ in range(nfield):
                    k, v = self.variable()
                    e[k] = v
                elems.append(e)
            return name, elems[0] if count == 1 else elems
        if typ == "cell":
            return name, [self.variable()[1] for _ in range(count)]
        if typ == "char":
            return name, self.text(count)
        if typ == "function_handle":                  # content is a nested char variable (func2str)
            return name, self.variable()[1]
        if typ in _NUMERIC:
            dt = np.dtype(_NUMERIC[typ])
            a = np.frombuffer(self.b, dt, count, self.p).copy()
            self.p += count * dt.itemsize
            a = a.reshape(dims[::-1]).T if len(dims) > 1 else a       # MATLAB column-major
            if typ == "logical":
                a = a.astype(bool)
            return name, (a.item() if a.size == 1 else a)
        raise ValueError(f"unsupported bhv2 type {typ!r} for {name!r} at byte {self.p}")


def read_bhv2(path: str | Path, names: set[str] | None = None) -> dict:
    """All top-level variables (or only `names`) of a .bhv2 file, as nested dicts/lists/arrays."""
    r = _Reader(Path(path).read_bytes())
    out = {}
    while r.p < len(r.b):
        name, value = r.variable()
        if names is None or name in names:
            out[name] = value
    return out
