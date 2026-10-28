"""Blackrock NSx reader (NEURALCD, file spec 2.2 / 2.3).

Data packet timestamps are in NSP clock ticks (the header's time resolution,
30 kHz). Sample i of a packet with timestamp T sits at tick T + i * period.
In this dataset every .ns6 starts at tick 102 and every .ns2 at tick 0, so
reading by tick (not by file-relative sample index, as the legacy pipeline
did through NPMK) is what puts all streams on one clock.
"""

from __future__ import annotations

import datetime as dt
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from ._common import cstr, systemtime_utc

_BASIC_HEADER_BYTES = 314
_EXT_HEADER_BYTES = 66
_PACKET_HEADER_BYTES = 9


@dataclass(frozen=True)
class NsxChannel:
    electrode_id: int
    label: str
    bank: int
    pin: int
    min_digital: int
    max_digital: int
    min_analog: int
    max_analog: int
    units: str
    highpass_hz: float
    highpass_order: int
    lowpass_hz: float
    lowpass_order: int

    @property
    def scale(self) -> float:
        """Physical units per bit (uV for electrodes, mV for analog inputs)."""
        return (self.max_analog - self.min_analog) / (self.max_digital - self.min_digital)


@dataclass(frozen=True)
class NsxPacket:
    timestamp: int  # NSP tick of the first sample
    n_samples: int
    offset: int  # byte offset of the first sample in the file


class NsxFile:
    """Header + packet table of one NSx file; data are memory-mapped on demand."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        with open(self.path, "rb") as f:
            head = f.read(_BASIC_HEADER_BYTES)
            if head[:8] != b"NEURALCD":
                raise ValueError(f"{self.path.name}: not a NEURALCD (2.2/2.3) NSx file")
            self.version = (head[8], head[9])
            n_header_bytes = struct.unpack_from("<I", head, 10)[0]
            self.label = cstr(head[14:30])
            self.comment = cstr(head[30:286])
            self.period, self.time_resolution = struct.unpack_from("<II", head, 286)
            self.origin_utc: dt.datetime = systemtime_utc(head[294:310])
            n_channels = struct.unpack_from("<I", head, 310)[0]
            ext = f.read(_EXT_HEADER_BYTES * n_channels)
        self.channels = [self._parse_channel(ext[i * 66:(i + 1) * 66]) for i in range(n_channels)]
        self.packets = self._scan_packets(n_header_bytes)

    # ------------------------------------------------------------------ header
    @staticmethod
    def _parse_channel(b: bytes) -> NsxChannel:
        if b[:2] != b"CC":
            raise ValueError("bad NSx extended header")
        eid = struct.unpack_from("<H", b, 2)[0]
        bank, pin = b[20], b[21]
        dmin, dmax, amin, amax = struct.unpack_from("<4h", b, 22)
        hp_hz, hp_ord, _hp_type, lp_hz, lp_ord, _lp_type = struct.unpack_from("<IIHIIH", b, 46)
        return NsxChannel(eid, cstr(b[4:20]), bank, pin, dmin, dmax, amin, amax, cstr(b[30:46]),
                          hp_hz / 1000, hp_ord, lp_hz / 1000, lp_ord)

    def _scan_packets(self, offset: int) -> list[NsxPacket]:
        size = self.path.stat().st_size
        frame = 2 * len(self.channels)
        packets = []
        with open(self.path, "rb") as f:
            while offset < size:
                f.seek(offset)
                ph = f.read(_PACKET_HEADER_BYTES)
                if len(ph) < _PACKET_HEADER_BYTES or ph[0] != 1:
                    raise ValueError(f"{self.path.name}: bad data packet header at byte {offset}")
                ts, n = struct.unpack_from("<II", ph, 1)
                packets.append(NsxPacket(ts, n, offset + _PACKET_HEADER_BYTES))
                offset += _PACKET_HEADER_BYTES + n * frame
        if offset != size:
            raise ValueError(f"{self.path.name}: {size - offset} trailing bytes")
        return packets

    # -------------------------------------------------------------- properties
    @property
    def fs(self) -> float:
        return self.time_resolution / self.period

    @property
    def labels(self) -> list[str]:
        return [c.label for c in self.channels]

    @property
    def t_start(self) -> int:
        """NSP tick of the first sample."""
        return self.packets[0].timestamp

    @property
    def t_stop(self) -> int:
        """NSP tick just after the last sample."""
        p = self.packets[-1]
        return p.timestamp + p.n_samples * self.period

    def channel_index(self, channels: Sequence[str | int] | None) -> np.ndarray:
        if channels is None:
            return np.arange(len(self.channels))
        labels = self.labels
        idx = [labels.index(c) if isinstance(c, str) else int(c) for c in channels]
        return np.asarray(idx)

    def _memmap(self, p: NsxPacket) -> np.memmap:
        return np.memmap(self.path, dtype="<i2", mode="r", offset=p.offset,
                         shape=(p.n_samples, len(self.channels)))

    # -------------------------------------------------------------------- read
    def sample_ticks(self, start: int, stop: int) -> np.ndarray:
        """Sample grid (NSP ticks) covering [start, stop), phase-locked to the first packet."""
        t0, per = self.t_start, self.period
        k0 = -((t0 - start) // per)  # ceil((start - t0) / per)
        k1 = -((t0 - stop) // per)
        return t0 + per * np.arange(k0, k1, dtype=np.int64)

    def read(self, start: int, stop: int, channels: Sequence[str | int] | None = None
             ) -> tuple[np.ndarray, np.ndarray]:
        """Samples with NSP tick in [start, stop), in physical units (float32).

        Returns (data [n_samples, n_channels], ticks [n_samples]). Samples that
        fall outside the recording or in gaps between packets are NaN.
        """
        ci = self.channel_index(channels)
        ticks = self.sample_ticks(start, stop)
        out = np.full((ticks.size, ci.size), np.nan, dtype=np.float32)
        if ticks.size == 0:
            return out, ticks
        for p in self.packets:
            # grid positions of this packet's samples (rounded if a packet is off-phase)
            first = int(round((p.timestamp - ticks[0]) / self.period))
            lo, hi = max(first, 0), min(first + p.n_samples, ticks.size)
            if lo >= hi:
                continue
            block = self._memmap(p)[lo - first:hi - first][:, ci]
            out[lo:hi] = block
        scale = np.array([self.channels[i].scale for i in ci], dtype=np.float32)
        return out * scale, ticks

    def read_raw(self, start: int, stop: int, channels: Sequence[str | int] | None = None
                 ) -> tuple[np.ndarray, np.ndarray]:
        """Like `read` but int16 digital values; the range must lie inside one packet."""
        ci = self.channel_index(channels)
        ticks = self.sample_ticks(start, stop)
        if ticks.size == 0:
            return np.empty((0, ci.size), dtype="<i2"), ticks
        for p in self.packets:
            first = (ticks[0] - p.timestamp) / self.period
            if first >= 0 and first == int(first) and int(first) + ticks.size <= p.n_samples:
                i = int(first)
                return np.asarray(self._memmap(p)[i:i + ticks.size][:, ci]), ticks
        raise ValueError(f"[{start}, {stop}) is not inside a single data packet")
