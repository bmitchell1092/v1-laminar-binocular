from __future__ import annotations

import datetime as dt
import struct


def cstr(b: bytes) -> str:
    """Decode a fixed-width, NUL-padded header string."""
    return b.split(b"\x00", 1)[0].decode("latin-1").strip()


def systemtime_utc(b: bytes) -> dt.datetime:
    """Windows SYSTEMTIME (8 x uint16) -> timezone-aware UTC datetime.

    Blackrock writes the recording origin in UTC (verified: MonkeyLogic local
    timestamps sit exactly -5 h / -6 h from it, CDT / CST).
    """
    year, month, _dow, day, hour, minute, second, ms = struct.unpack("<8H", b)
    return dt.datetime(year, month, day, hour, minute, second, ms * 1000, tzinfo=dt.timezone.utc)
