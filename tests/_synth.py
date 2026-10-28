"""Writers for tiny synthetic Blackrock files (spec 2.3), used only by tests."""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

SYSTEMTIME_2021_10_08 = struct.pack("<8H", 2021, 10, 5, 8, 23, 25, 45, 0)


def write_nsx(path: Path, labels: list[str], packets: list[tuple[int, np.ndarray]],
              period: int = 1, analog: tuple[int, int] = (-8191, 8191),
              digital: tuple[int, int] = (-32764, 32764)) -> Path:
    n = len(labels)
    header_bytes = 314 + 66 * n
    head = bytearray(314)
    head[0:8] = b"NEURALCD"
    head[8:10] = bytes([2, 3])
    struct.pack_into("<I", head, 10, header_bytes)
    head[14:30] = b"synthetic".ljust(16, b"\x00")
    struct.pack_into("<II", head, 286, period, 30000)
    head[294:310] = SYSTEMTIME_2021_10_08
    struct.pack_into("<I", head, 310, n)
    ext = bytearray()
    for i, lab in enumerate(labels):
        e = bytearray(66)
        e[0:2] = b"CC"
        struct.pack_into("<H", e, 2, i + 1)
        e[4:20] = lab.encode().ljust(16, b"\x00")
        e[20], e[21] = 1, i + 1
        struct.pack_into("<4h", e, 22, digital[0], digital[1], analog[0], analog[1])
        e[30:46] = b"uV".ljust(16, b"\x00")
        struct.pack_into("<IIHIIH", e, 46, 300, 1, 1, 7500000, 3, 1)
        ext += e
    body = bytearray()
    for ts, data in packets:
        data = np.asarray(data, dtype="<i2").reshape(-1, n)
        body += struct.pack("<BII", 1, ts, data.shape[0]) + data.tobytes()
    path.write_bytes(bytes(head) + bytes(ext) + bytes(body))
    return path


def write_nev(path: Path, digital: list[tuple[int, int]], spikes: list[tuple[int, int, int]] = (),
              packet_bytes: int = 104) -> Path:
    head = bytearray(336)
    head[0:8] = b"NEURALEV"
    head[8:10] = bytes([2, 3])
    struct.pack_into("<HIIII", head, 10, 1, 336, packet_bytes, 30000, 30000)
    head[28:44] = SYSTEMTIME_2021_10_08
    struct.pack_into("<I", head, 332, 0)
    rows = []
    for ts, value in digital:
        p = bytearray(packet_bytes)
        struct.pack_into("<IHBBH", p, 0, ts, 0, 1, 0, value)
        rows.append((ts, p))
    for ts, eid, unit in spikes:
        p = bytearray(packet_bytes)
        struct.pack_into("<IHBB", p, 0, ts, eid, unit, 0)
        p[8:8 + 4] = struct.pack("<2h", -100, 50)
        rows.append((ts, p))
    rows.sort(key=lambda r: r[0])
    path.write_bytes(bytes(head) + b"".join(bytes(p) for _, p in rows))
    return path
