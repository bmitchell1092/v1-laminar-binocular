"""Kilosort 2 / phy single units (2023 curation), moved onto the NSP clock and re-qualified.

Sorts live in `<raw_root>/KiloSort-ed_phy/<session>_<recording>/`, one bminteroc file each: Kilosort ran per
file on a raw int16 copy of the .ns6 with every channel (`n_channels_dat`). `cluster_info.tsv` `ch` is the row
of that copy (= the .ns6 header order), mapped to electrode labels by `chanMap.mat` (`chanMap0ind` ->
`chanMapLabels`). On dual-probe days a sort covers both banks (x = 1000 / 2000 um), so units are assigned to a
penetration by the bank of their label.

Spike times: sample s of the copy is NSP tick T0 + s (T0 = first packet timestamp, 102), but in the sorts made
before 221206 the waveform trough sits 204 samples (6.8 ms) before T0 + s (05 A7). The offset is therefore
measured per sort from mean raw waveforms of its best units, separately in each half of the recording, and
corrected spike ticks put the trough at the spike time.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import signal
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter1d

from .config import raw_root, recording_files
from .io import NsxFile

KS_DIRNAME = "KiloSort-ed_phy"
REFRACTORY_S = 0.0015
SNIPPET = (-30, 60)               # samples around the spike (-1 ... +2 ms)

# well isolated: phy 'good' and all of these (Hill et al. 2011; Siegle et al. 2021 / Allen Brain Observatory)
MAX_ISI_CONTAMINATION = 0.1       # estimated fraction of spikes from other neurons, from refractory violations
MAX_AMPLITUDE_CUTOFF = 0.1        # estimated fraction of spikes lost below the detection threshold
MIN_PRESENCE_RATIO = 0.9          # fraction of 1-min bins with spikes
MIN_SNR = 3.0                     # mean waveform trough / noise SD on the peak channel


def sort_dir(session: str, recording: str) -> Path | None:
    d = raw_root() / KS_DIRNAME / f"{session}_{recording}"
    return d if d.is_dir() else None


@dataclass
class PhySort:
    path: Path
    row_labels: np.ndarray        # electrode label of each row of the sorted copy
    spike_samples: np.ndarray     # int64, as stored by Kilosort
    spike_clusters: np.ndarray
    amplitudes: np.ndarray        # Kilosort template scaling per spike
    info: pd.DataFrame            # cluster_info.tsv + 'label'


def load_sort(session: str, recording: str) -> PhySort:
    d = sort_dir(session, recording)
    if d is None:
        raise FileNotFoundError(f"no Kilosort/phy folder for {session}_{recording}")
    m = loadmat(d / "chanMap.mat", squeeze_me=True)
    rows = np.atleast_1d(m["chanMap0ind"]).astype(int)
    labels = np.full(rows.max() + 1, "", dtype=object)
    labels[rows] = [str(x) for x in np.atleast_1d(m["chanMapLabels"])]
    info = pd.read_csv(d / "cluster_info.tsv", sep="\t")
    info["label"] = labels[info["ch"].to_numpy()]
    return PhySort(d, labels, np.load(d / "spike_times.npy").ravel().astype(np.int64),
                   np.load(d / "spike_clusters.npy").ravel(), np.load(d / "amplitudes.npy").ravel(), info)


def _snippets(ns6: NsxFile, ticks: np.ndarray, labels: list[str], window=SNIPPET, pad: int = 90) -> np.ndarray:
    """High-passed (300 Hz) raw snippets [spike, time, channel] around each tick (uV)."""
    sos = signal.butter(3, 300, "highpass", fs=ns6.fs, output="sos")
    out = np.empty((len(ticks), window[1] - window[0], len(labels)), np.float32)
    for i, t in enumerate(ticks):
        x, _ = ns6.read(t + window[0] - pad, t + window[1] + pad, labels)
        out[i] = signal.sosfiltfilt(sos, np.nan_to_num(x), axis=0)[pad:pad + window[1] - window[0]]
    return out


def measure_offset(sort: PhySort, ns6: NsxFile, n_units: int = 5, n_spikes: int = 150,
                   search=(-400, 100), seed: int = 0) -> dict:
    """Samples to add to T0 + s so that the spike time falls on the waveform trough.

    Uses the `n_units` good units with most spikes; for each half of the recording, the trough of the mean
    high-passed waveform on the unit's channel within `search`. Returns the median lag over units with a clear
    waveform, the per-unit lags, and `consistent`: each unit's halves agree within 2 samples (no drift) and
    most units sit within 2 samples of the median (a unit whose trough is shaped differently may not)."""
    rng = np.random.default_rng(seed)
    t0 = ns6.t_start
    good = sort.info[sort.info.group == "good"].sort_values("n_spikes", ascending=False).head(n_units)
    mid = np.median(sort.spike_samples)
    lags = {}
    for _, u in good.iterrows():
        s = sort.spike_samples[sort.spike_clusters == u.id]
        s = s[(t0 + s + search[0] > ns6.t_start + 200) & (t0 + s + search[1] < ns6.t_stop - 200)]
        per_half = []
        for half in (s[s < mid], s[s >= mid]):
            if len(half) < 30:
                continue
            pick = rng.choice(half, min(n_spikes, len(half)), replace=False)
            w = _snippets(ns6, t0 + pick, [u.label], window=search).mean(axis=0)[:, 0]
            noise = np.std(w[: 50]) + 1e-9
            if -w.min() > 4 * noise:
                per_half.append(int(np.argmin(w)) + search[0])
        if per_half:
            lags[int(u.id)] = per_half
    flat = [l for v in lags.values() for l in v]
    if not flat:
        return dict(offset=np.nan, lags=lags, consistent=False)
    offset = int(np.median(flat))
    halves_agree = all(max(v) - min(v) <= 2 for v in lags.values())
    near = np.mean([abs(np.median(v) - offset) <= 2 for v in lags.values()])
    return dict(offset=offset, lags=lags, consistent=bool(halves_agree and near > 0.5))


def isi_contamination(spike_s: np.ndarray, duration_s: float, refractory_s: float = REFRACTORY_S) -> float:
    """Hill et al. (2011) false-positive rate from refractory-period violations (Allen implementation)."""
    n = len(spike_s)
    if n < 2:
        return np.nan
    viol = np.sum(np.diff(np.sort(spike_s)) < refractory_s)
    return float((viol / (2 * n * refractory_s)) / (n / duration_s))


def amplitude_cutoff(amplitudes: np.ndarray, n_bins: int = 500, smooth: float = 3) -> float:
    """Fraction of spikes missing below the detection threshold, from the truncated amplitude distribution
    (Siegle et al. 2021 / Allen implementation; capped at 0.5)."""
    if len(amplitudes) < 20:
        return np.nan
    pdf, edges = np.histogram(amplitudes, bins=n_bins, density=True)
    pdf = gaussian_filter1d(pdf, smooth)
    peak = int(np.argmax(pdf))
    g = int(np.argmin(np.abs(pdf[peak:] - pdf[0]))) + peak
    return float(min(np.sum(pdf[g:]) * np.mean(np.diff(edges)), 0.5))


def presence_ratio(spike_s: np.ndarray, start_s: float, stop_s: float, bin_s: float = 60.0) -> float:
    edges = np.arange(start_s, stop_s + bin_s, bin_s)
    counts, _ = np.histogram(spike_s, edges)
    return float(np.mean(counts > 0))


def noise_sd(ns6: NsxFile, labels: list[str], seconds: float = 10.0) -> np.ndarray:
    """Robust SD (MAD / 0.6745) of the 300 Hz high-passed signal, from the middle of the recording."""
    mid = (ns6.t_start + ns6.t_stop) // 2
    x, _ = ns6.read(mid, mid + int(seconds * ns6.fs) * ns6.period, labels)
    sos = signal.butter(3, 300, "highpass", fs=ns6.fs, output="sos")
    y = signal.sosfiltfilt(sos, np.nan_to_num(x), axis=0)
    return np.median(np.abs(y - np.median(y, axis=0)), axis=0) / 0.6745


def units(session: str, recording: str, bank: str | None = None, n_waveform: int = 300, seed: int = 0
          ) -> tuple[pd.DataFrame, dict[int, np.ndarray], dict]:
    """Curated units of one sorted recording on the NSP clock, with isolation metrics.

    Returns (table, spike ticks per cluster id, offset report). The table has one row per phy 'good' or
    'mua' cluster (noise dropped) of `bank` (all banks if None): phy group/KS label, electrode label and
    contact, peak channel and waveform from the raw data, firing rate, isi_contamination, amplitude_cutoff,
    presence_ratio, snr, trough_to_peak_ms, and `well_isolated`."""
    sort = load_sort(session, recording)
    ns6 = NsxFile(recording_files(session, recording)[".ns6"])
    rep = measure_offset(sort, ns6)
    if not rep["consistent"]:
        raise RuntimeError(f"{session}_{recording}: spike-time offset not consistent: {rep}")
    ticks_all = ns6.t_start + sort.spike_samples + rep["offset"]
    inside = (ticks_all >= ns6.t_start) & (ticks_all < ns6.t_stop)
    start_s, stop_s = ns6.t_start / 30000, ns6.t_stop / 30000
    info = sort.info[sort.info.group.isin(["good", "mua"])]
    if bank is not None:
        info = info[info.label.str.startswith(f"e{bank}")]
    bank_labels = sorted({lab for lab in sort.row_labels if isinstance(lab, str) and lab[:2] in
                          ({f"e{bank}"} if bank else {"eA", "eB"})}, key=lambda s: (s[:2], int(s[2:])))
    nsd = dict(zip(bank_labels, noise_sd(ns6, bank_labels)))
    rng = np.random.default_rng(seed)
    rows, spikes = [], {}
    for _, u in info.iterrows():
        sel = (sort.spike_clusters == u.id) & inside
        tk = np.sort(ticks_all[sel])
        spikes[int(u.id)] = tk
        s = tk / 30000
        dur = stop_s - start_s
        row = dict(cluster_id=int(u.id), group=u.group, ks_label=u.KSLabel, label=u.label,
                   bank=u.label[1], contact=int(u.label[2:]), n_spikes=len(tk), rate_hz=len(tk) / dur,
                   ks_contam_pct=u.ContamPct, isi_contamination=isi_contamination(s, dur),
                   amplitude_cutoff=amplitude_cutoff(sort.amplitudes[sel]),
                   presence_ratio=presence_ratio(s, start_s, stop_s))
        # waveform on the probe the unit belongs to, from the raw data
        probe = [lab for lab in bank_labels if lab[:2] == u.label[:2]]
        ok = tk[(tk > ns6.t_start + 200) & (tk < ns6.t_stop - 200)]
        pick = rng.choice(ok, min(n_waveform, len(ok)), replace=False) if len(ok) else ok
        if len(pick) >= 20:
            w = _snippets(ns6, pick, probe).mean(axis=0)                       # [time, channel]
            troughs = -w.min(axis=0)
            pk = int(np.argmax(troughs))
            wv = w[:, pk]
            it = int(np.argmin(wv))
            ip = it + int(np.argmax(wv[it:]))
            row.update(peak_label=probe[pk], peak_contact=int(probe[pk][2:]), trough_uv=float(wv[it]),
                       snr=float(troughs[pk] / nsd[probe[pk]]), trough_lag=it + SNIPPET[0],
                       trough_to_peak_ms=(ip - it) / 30.0)
        rows.append(row)
    tab = pd.DataFrame(rows)
    if len(tab):
        tab["well_isolated"] = ((tab.group == "good") & (tab.isi_contamination < MAX_ISI_CONTAMINATION)
                                & (tab.amplitude_cutoff < MAX_AMPLITUDE_CUTOFF)
                                & (tab.presence_ratio >= MIN_PRESENCE_RATIO) & (tab.snr >= MIN_SNR))
    return tab, spikes, rep
