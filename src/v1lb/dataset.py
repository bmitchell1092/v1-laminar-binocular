"""Per-penetration data files: one xarray Dataset (NetCDF/HDF5) per penetration.

Layout of `<out_root>/datasets/<penetration>.nc`:
  dims        presentation, time (ms, 1 kHz, 0 = photodiode flash), channel (superficial -> deep)
  data vars   lfp [uV], csd [nA/mm^3], mua [uV, analog envelope], spikes [threshold-crossing counts/ms]
  coords      presentation: every column of `bminteroc.presentation_table` (contrast_L/R/ipsi/contra,
              condition, NSP ticks, frames, previous presentation, outcome, ...)
              channel: contact (1-based label number), depth_l45 (contacts, + = superficial),
              depth_um, in_cortex, beyond_deep_limit
  attrs       provenance and definitions
Only presentations with a detected photodiode flash are included; filter on `completed` as needed.
"""

from __future__ import annotations

import datetime as dt
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import __version__
from .bminteroc import presentation_table
from .config import metadata_dir, out_root, recording_files
from .io import NsxFile
from .laminar import csd as compute_csd
from .signals import analog_mua_epochs, probe_labels, threshold_spikes

WINDOW_MS = (-200, 500)
ALIGNMENT_FILE = "laminar_alignment_v1.csv"


def _git_hash() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                              cwd=Path(__file__).resolve().parents[2], timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def channel_table(n: int, orientation: str, l45: float, top_of_cortex: int, deep_limit: float) -> pd.DataFrame:
    """Contacts ordered superficial -> deep with depth relative to L4/5 (positive = superficial)."""
    contacts = np.arange(1, n + 1) if orientation == "ascending" else np.arange(n, 0, -1)
    sign = 1 if orientation == "ascending" else -1
    depth = sign * (l45 - contacts)
    in_cortex = sign * (contacts - top_of_cortex) >= 0 if top_of_cortex > 0 else np.ones(n, bool)
    beyond = sign * (contacts - deep_limit) > 0 if np.isfinite(deep_limit) else np.zeros(n, bool)
    return pd.DataFrame({"contact": contacts, "depth_l45": depth.astype(float), "depth_um": depth * 100.0,
                         "in_cortex": in_cortex, "beyond_deep_limit": beyond})


def lfp_epochs(ns2: NsxFile, labels: list[str], onsets: np.ndarray, window_ms=WINDOW_MS) -> np.ndarray:
    """LFP [events, time, channels] (uV) on an exact 1 ms grid locked to each onset tick.

    The 1 kHz samples sit on multiples of the NS2 period; each epoch is linearly
    interpolated onto onset + k ms (shift < 1 sample; LFP is band-limited to 500 Hz)."""
    per = ns2.period
    t_rel = np.arange(window_ms[0], window_ms[1]) * per          # ticks relative to onset
    out = np.empty((len(onsets), t_rel.size, len(labels)), np.float32)
    for i, t in enumerate(np.asarray(onsets)):
        x, ticks = ns2.read(t + t_rel[0] - per, t + t_rel[-1] + 2 * per, labels)
        target = t + t_rel
        for c in range(x.shape[1]):
            out[i, :, c] = np.interp(target, ticks, x[:, c])
    return out


def bin_spikes(spike_ticks: list[np.ndarray], onsets: np.ndarray, window_ms=WINDOW_MS, ticks_per_ms: int = 30) -> np.ndarray:
    """Spike counts per 1 ms bin [events, time, channels] (uint8), bins [onset + k ms, onset + (k+1) ms)."""
    n_t = window_ms[1] - window_ms[0]
    out = np.zeros((len(onsets), n_t, len(spike_ticks)), np.uint8)
    for c, st in enumerate(spike_ticks):
        st = np.sort(np.asarray(st))
        for i, t in enumerate(np.asarray(onsets)):
            a = np.searchsorted(st, t + window_ms[0] * ticks_per_ms)
            b = np.searchsorted(st, t + window_ms[1] * ticks_per_ms)
            if b > a:
                k = (st[a:b] - t - window_ms[0] * ticks_per_ms) // ticks_per_ms
                np.add.at(out[i, :, c], k, 1)
    return out


def build_penetration(penetration: str, out_dir: Path | None = None, window_ms=WINDOW_MS,
                      max_presentations: int | None = None, write: bool = True) -> xr.Dataset:
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    align = pd.read_csv(metadata_dir() / ALIGNMENT_FILE).set_index("penetration")
    meta, al = pens.loc[penetration], align.loc[penetration]
    session, bank, rec = meta["session"], meta["bank"], meta["bminteroc"]

    files = recording_files(session, rec)
    ns2, ns6 = NsxFile(files[".ns2"]), NsxFile(files[".ns6"])
    labels = probe_labels(ns2, bank)
    chans = channel_table(len(labels), al["orientation"], float(al["l45"]), int(al["top_of_cortex"]),
                          float(al["deep_limit"]) if pd.notna(al["deep_limit"]) else np.nan)
    ordered = [labels[c - 1] for c in chans["contact"]]          # superficial -> deep

    P, pmeta = presentation_table(session, rec)
    P = P[P["pd_on_tick"] > 0].reset_index(drop=True)
    if max_presentations:
        P = P.iloc[:max_presentations]
    onsets = P["pd_on_tick"].to_numpy()

    lfp = lfp_epochs(ns2, ordered, onsets, window_ms)
    csd = compute_csd(lfp).astype(np.float32)                    # along channels (superficial -> deep)
    mua = analog_mua_epochs(ns6, ordered, onsets, pre_ms=-window_ms[0], post_ms=window_ms[1]).astype(np.float32)
    span = (int(onsets.min()) + window_ms[0] * 30 - 30000, int(onsets.max()) + window_ms[1] * 30 + 30000)
    spikes = bin_spikes(threshold_spikes(ns6, ordered, max(span[0], ns6.t_start), min(span[1], ns6.t_stop)),
                        onsets, window_ms)

    t_ms = np.arange(window_ms[0], window_ms[1]).astype(float)
    pres_coords = {}
    for col in P.columns:
        v = P[col]
        pres_coords[col] = ("presentation", v.astype(str).to_numpy() if v.dtype == object else v.to_numpy())
    ds = xr.Dataset(
        data_vars=dict(
            lfp=(("presentation", "time", "channel"), lfp, {"units": "uV", "source": "ns2 (0.3-500 Hz hardware band)"}),
            csd=(("presentation", "time", "channel"), csd,
                 {"units": "nA/mm^3", "definition": "-sigma d2V/dz2, sigma 0.4 S/m, 0.1 mm spacing; sinks negative; edge channels NaN"}),
            mua=(("presentation", "time", "channel"), mua,
                 {"units": "uV", "definition": "ns6: 500-5000 Hz band-pass, rectify, 200 Hz low-pass (4th-order, zero-phase), 1 kHz"}),
            spikes=(("presentation", "time", "channel"), spikes,
                    {"units": "count per 1 ms bin", "definition": "negative crossings of -4 SD (1 s windows) after 500 Hz high-pass; 1 ms dead time"}),
        ),
        coords=dict(time=("time", t_ms, {"units": "ms", "zero": "first backlight flash at the photodiode (NSP clock)"}),
                    channel=("channel", ordered),
                    presentation=("presentation", np.arange(len(P))),
                    **{k: ("channel", chans[k].to_numpy()) for k in chans.columns},
                    **pres_coords),
        attrs=dict(penetration=penetration, session=session, recording=rec, bank=bank,
                   hemisphere=meta["hemisphere"], orientation=al["orientation"], l45=float(al["l45"]),
                   top_of_cortex=int(al["top_of_cortex"]),
                   deep_limit=float(al["deep_limit"]) if pd.notna(al["deep_limit"]) else np.nan,
                   alignment=f"{ALIGNMENT_FILE} ({al['basis']})", contact_spacing_um=100.0,
                   time_zero=("first backlight flash at the photodiode (bottom-right of the screen), "
                              "detected per presentation; the possible 2.6-3.4 ms scan offset to the "
                              "stimulus row is NOT applied (docs/05 A2b)"),
                   scan_offset_note="ViewPixx backlight mode unknown",
                   trial_matching_shift=int(pmeta["trial_shift"]), ml_minus_nsp_s=float(pmeta["ml_minus_nsp_s"]),
                   pipeline="v1lb " + __version__, git=_git_hash(),
                   created=dt.datetime.now().isoformat(timespec="seconds")),
    )
    if write:
        out_dir = Path(out_dir or out_root() / "datasets")
        out_dir.mkdir(parents=True, exist_ok=True)
        enc = {v: {"zlib": True, "complevel": 4} for v in ("lfp", "csd", "mua", "spikes")}
        ds.to_netcdf(out_dir / f"{penetration}.nc", engine="h5netcdf", encoding=enc)
        P.to_parquet(out_dir / f"{penetration}_presentations.parquet")
    return ds


def open_penetration(penetration: str, data_dir: Path | None = None) -> xr.Dataset:
    """Open a built dataset lazily."""
    return xr.open_dataset(Path(data_dir or out_root() / "datasets") / f"{penetration}.nc", engine="h5netcdf")
