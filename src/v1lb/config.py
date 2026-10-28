"""Data locations and raw-file lookup.

The raw root can be overridden with the environment variable `V1LB_RAW_ROOT`.
Raw data are read-only: nothing in this package writes under `raw_root()`.
"""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_RAW_ROOT = r"D:\NEUROPHYS DATA\2021"


def raw_root() -> Path:
    return Path(os.environ.get("V1LB_RAW_ROOT", DEFAULT_RAW_ROOT))


def session_dir(session: str) -> Path:
    """Raw folder of a session, e.g. '211008_B' -> <root>/Raw Data B/211008_B."""
    monkey = session.split("_")[1]
    return raw_root() / f"Raw Data {monkey}" / session


def recording_files(session: str, recording: str) -> dict[str, Path]:
    """All files of one recording, keyed by extension (e.g. '.ns6', '.gBMINTEROCGrating_di').

    `recording` is the part after the session prefix, e.g. 'bminteroc001'.
    Matching is case-insensitive and tolerates the 'datafile_' prefix used for
    some 211210_B recordings. Raises if two files claim the same extension.
    """
    folder = session_dir(session)
    stem = f"{session}_{recording}".lower()
    found: dict[str, Path] = {}
    for p in folder.iterdir():
        if not p.is_file():
            continue
        name = p.name.lower().removeprefix("datafile_")
        if not name.startswith(stem + "."):
            continue
        ext = p.name[p.name.lower().find(stem) + len(stem):]
        if ext in found:
            raise ValueError(f"ambiguous {ext} for {session}_{recording}: {found[ext].name}, {p.name}")
        found[ext] = p
    if not found:
        raise FileNotFoundError(f"no files for {session}_{recording} in {folder}")
    return found
