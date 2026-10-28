"""Blackrock NEV reader (NEURALEV, file spec 2.2 / 2.3): digital events and spikes.

Timestamps are NSP clock ticks (30 kHz), the same clock as NSx packet timestamps.
"""

from __future__ import annotations

import datetime as dt
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._common import cstr, systemtime_utc

_BASIC_HEADER_BYTES = 336
_EXT_HEADER_BYTES = 32


@dataclass(frozen=True)
class NevElectrode:
    electrode_id: int
    label: str
    bank: int
    pin: int
    nv_per_bit: int
    low_threshold: int
    high_threshold: int
    spike_width: int  # samples


class NevFile:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        with open(self.path, "rb") as f:
            head = f.read(_BASIC_HEADER_BYTES)
            if head[:8] != b"NEURALEV":
                raise ValueError(f"{self.path.name}: not a NEURALEV file")
            self.version = (head[8], head[9])
            (self.flags, self.n_header_bytes, self.packet_bytes,
             self.time_resolution, self.sample_resolution) = struct.unpack_from("<HIIII", head, 10)
            self.origin_utc: dt.datetime = systemtime_utc(head[28:44])
            self.application = cstr(head[44:76])
            self.comment = cstr(head[76:332])
            n_ext = struct.unpack_from("<I", head, 332)[0]
            ext = f.read(_EXT_HEADER_BYTES * n_ext)
        self.electrodes, self.digital_labels = self._parse_ext(ext, n_ext)
        n_bytes = self.path.stat().st_size - self.n_header_bytes
        if n_bytes % self.packet_bytes:
            raise ValueError(f"{self.path.name}: data section is not a whole number of packets")
        self.n_packets = n_bytes // self.packet_bytes

    @staticmethod
    def _parse_ext(ext: bytes, n: int) -> tuple[dict[int, NevElectrode], list[tuple[str, int]]]:
        wav, lbl, dig = {}, {}, []
        for i in range(n):
            b = ext[i * 32:(i + 1) * 32]
            kind = b[:8]
            if kind == b"NEUEVWAV":
                eid, bank, pin, nv, _energy, hi, lo, _nunits, _bps, width = struct.unpack_from("<HBBHHhhBBH", b, 8)
                wav[eid] = (bank, pin, nv, lo, hi, width)
            elif kind == b"NEUEVLBL":
                eid = struct.unpack_from("<H", b, 8)[0]
                lbl[eid] = cstr(b[10:26])
            elif kind == b"DIGLABEL":
                dig.append((cstr(b[8:24]), b[24]))
        electrodes = {eid: NevElectrode(eid, lbl.get(eid, ""), *w[:3], w[3], w[4], w[5]) for eid, w in wav.items()}
        return electrodes, dig

    def _packets(self) -> np.memmap:
        return np.memmap(self.path, dtype=np.uint8, mode="r", offset=self.n_header_bytes,
                         shape=(self.n_packets, self.packet_bytes))

    @staticmethod
    def _u(a: np.ndarray, lo: int, hi: int, dtype: str) -> np.ndarray:
        return np.ascontiguousarray(a[:, lo:hi]).view(dtype).ravel()

    def digital_events(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(ticks int64, digital value uint16, insertion reason uint8) for packet ID 0."""
        a = self._packets()
        pid = self._u(a, 4, 6, "<u2")
        d = np.asarray(a[pid == 0])
        return (self._u(d, 0, 4, "<u4").astype(np.int64),
                self._u(d, 8, 10, "<u2").copy(),
                d[:, 6].copy())

    def spikes(self, waveforms: bool = False):
        """(ticks int64, electrode id uint16, unit uint8[, waveforms int16 n x samples])."""
        a = self._packets()
        pid = self._u(a, 4, 6, "<u2")
        s = np.asarray(a[(pid >= 1) & (pid <= 2048)])
        out = (self._u(s, 0, 4, "<u4").astype(np.int64), self._u(s, 4, 6, "<u2").copy(), s[:, 6].copy())
        if waveforms:
            out += (np.ascontiguousarray(s[:, 8:]).view("<i2").reshape(len(s), -1),)
        return out
