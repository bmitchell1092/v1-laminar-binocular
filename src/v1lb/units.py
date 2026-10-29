"""Kilosort 2 / phy single units (2023 curation), moved onto the NSP clock and re-qualified.

Sorts live in `<raw_root>/KiloSort-ed_phy/<session>_<recording>/`, one bminteroc file each: Kilosort ran per
file on a raw int16 copy of the .ns6 with every channel (`n_channels_dat`). `cluster_info.tsv` `ch` is the row
of that copy (= the .ns6 header order), mapped to electrode labels by `chanMap.mat` (`chanMap0ind` ->
`chanMapLabels`). On dual-probe days a sort covers both banks (x = 1000 / 2000 um), so units are assigned to a
penetration by the bank of their label.

Spike times: sample s of the copy is NSP tick T0 + s (T0 = first packet timestamp, 102), but in the sorts made
before 221206 the waveform trough sits 204 samples = 2 x T0 (6.8 ms) before T0 + s (05 A7). The offset is
measured per sort from mean raw waveforms of its best units, separately in each half of the recording. Single
units scatter by several samples (trough timing depends on waveform shape), so for the 2023 sorts the measured
value is snapped to the generation it belongs to (-2 x T0 or 0). Kilosort 4 reports the trough sample itself
(offset 0; measured 0 in the pilot), so its sorts use 0 and the measurement is a check.

`well_isolated` is metric-based (any non-noise cluster); the curation label is kept in `group`.
"""

from __future__ import annotations

import json
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
# 2023 sorts not to use (Kilosort 4 instead)
BROKEN_SORTS = {"221118_J_bminteroc001": "truncated: spikes end 4.4 min into a 32 min recording; no real units"}
REFRACTORY_S = 0.0015
SNIPPET = (-30, 60)               # samples around the spike (-1 ... +2 ms)

# well isolated: all of these (Hill et al. 2011; Siegle et al. 2021 / Allen Brain Observatory)
MAX_ISI_CONTAMINATION = 0.1       # estimated fraction of spikes from other neurons, from refractory violations
MAX_AMPLITUDE_CUTOFF = 0.1        # estimated fraction of spikes lost below the detection threshold
                                  # (from raw per-spike trough amplitudes on the peak channel, so it is
                                  # comparable across sorters; Kilosort's own amplitudes: ks_amplitude_cutoff)
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
    info: pd.DataFrame            # id, group, KSLabel, ContamPct, ch (row), label
    generations: tuple | None = None   # known spike-time offsets as multiples of T0 (2023 sorts: (-2, 0))
    packets: list | None = None        # [(first tick, n samples)] when the sorted stream joined several packets


def base_ticks(sort: PhySort, ns6: NsxFile, samples: np.ndarray | None = None) -> np.ndarray:
    """NSP tick of each sorted sample (before the spike-time offset): T0 + s for a single packet; packet by
    packet for a stream that joined several."""
    s = sort.spike_samples if samples is None else np.asarray(samples, np.int64)
    if not sort.packets or len(sort.packets) == 1:
        return ns6.t_start + s
    starts = np.cumsum([0] + [n for _, n in sort.packets[:-1]])
    k = np.clip(np.searchsorted(starts, s, side="right") - 1, 0, len(starts) - 1)
    return np.array([t for t, _ in sort.packets], np.int64)[k] + (s - starts[k]) * ns6.period


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
                   np.load(d / "spike_clusters.npy").ravel(), np.load(d / "amplitudes.npy").ravel(), info,
                   generations=(-2, 0))


def load_ks4(results_dir: str | Path, ns6_labels: list[str]) -> PhySort:
    """Uncurated Kilosort 4 output (sorted straight from the .ns6, rows = .ns6 header order). Every cluster
    is kept (group = Kilosort's own label); the best channel is the template's largest peak-to-peak."""
    d = Path(results_dir)
    labels = np.array(ns6_labels, dtype=object)
    cmap = np.load(d / "channel_map.npy").ravel()
    tmpl = np.load(d / "templates.npy")                                   # [template, time, channel]
    best = cmap[np.ptp(tmpl, axis=1).argmax(axis=1)]
    ks = pd.read_csv(d / "cluster_KSLabel.tsv", sep="\t")
    contam = pd.read_csv(d / "cluster_ContamPct.tsv", sep="\t")
    info = ks.merge(contam, on="cluster_id").rename(columns={"cluster_id": "id"})
    clusters = np.load(d / "spike_clusters.npy").ravel()
    info["n_spikes"] = np.bincount(clusters, minlength=info.id.max() + 1)[info.id]
    info["ch"] = best[info.id.to_numpy()]                     # KS4: cluster id = template id
    info["group"] = info["KSLabel"]
    info["label"] = labels[info["ch"].to_numpy()]
    run = d / "run.json"
    packets = json.loads(run.read_text()).get("packets") if run.exists() else None
    return PhySort(d, labels, np.load(d / "spike_times.npy").ravel().astype(np.int64), clusters,
                   np.load(d / "amplitudes.npy").ravel(), info, generations=(0,), packets=packets)


def _snippets(ns6: NsxFile, ticks: np.ndarray, labels: list[str], window=SNIPPET, pad: int = 90,
              chunk: int = 500) -> np.ndarray:
    """High-passed (300 Hz) raw snippets [spike, time, channel] around each tick (uV).

    Single-packet files are indexed straight from the memory map (vectorised); otherwise one read per spike."""
    sos = signal.butter(3, 300, "highpass", fs=ns6.fs, output="sos")
    n_t = window[1] - window[0]
    ticks = np.asarray(ticks, np.int64)
    out = np.empty((len(ticks), n_t, len(labels)), np.float32)
    if len(ns6.packets) == 1:
        p = ns6.packets[0]
        mm = ns6._memmap(p)
        cols = ns6.channel_index(labels)
        scale = np.array([ns6.channels[i].scale for i in cols], np.float32)
        rel = np.arange(window[0] - pad, window[1] + pad)
        for a in range(0, len(ticks), chunk):
            idx = (ticks[a:a + chunk] - p.timestamp) // ns6.period
            rows = np.clip(idx[:, None] + rel[None, :], 0, p.n_samples - 1)
            x = mm[rows.ravel()][:, cols].reshape(len(idx), rel.size, len(cols)).astype(np.float32) * scale
            out[a:a + chunk] = signal.sosfiltfilt(sos, x, axis=1)[:, pad:pad + n_t]
        return out
    for i, t in enumerate(ticks):
        x, _ = ns6.read(t + window[0] - pad, t + window[1] + pad, labels)
        out[i] = signal.sosfiltfilt(sos, np.nan_to_num(x), axis=0)[pad:pad + n_t]
    return out


def measure_offset(sort: PhySort, ns6: NsxFile, n_units: int = 5, n_spikes: int = 150,
                   search=(-400, 100), seed: int = 0) -> dict:
    """Samples to add to T0 + s so that the spike time falls on the waveform trough.

    Uses the `n_units` clusters with most spikes and a clear waveform, 'good' first; for each half of the recording, the trough of the mean
    high-passed waveform on the unit's channel within `search`. Returns the median lag over units with a clear
    waveform, the per-unit lags, and `consistent`: most units' halves agree within 3 samples (no drift) and most
    units sit within 8 samples of the offset. With `sort.generations`, the offset is the nearest known value."""
    rng = np.random.default_rng(seed)
    t0 = ns6.t_start
    info = sort.info[sort.info.group.isin(["good", "mua"])]
    good = pd.concat([info[info.group == "good"].sort_values("n_spikes", ascending=False),
                      info[info.group != "good"].sort_values("n_spikes", ascending=False)]).head(2 * n_units)
    ticks = base_ticks(sort, ns6)
    mid = np.median(ticks)
    lags = {}
    for _, u in good.iterrows():
        s = ticks[sort.spike_clusters == u.id]
        s = s[(s + search[0] > ns6.t_start + 200) & (s + search[1] < ns6.t_stop - 200)]
        per_half = []
        for half in (s[s < mid], s[s >= mid]):
            if len(half) < 30:
                continue
            pick = rng.choice(half, min(n_spikes, len(half)), replace=False)
            w = _snippets(ns6, pick, [u.label], window=search).mean(axis=0)[:, 0]
            noise = np.std(w[: 50]) + 1e-9
            if -w.min() > 4 * noise:
                per_half.append(int(np.argmin(w)) + search[0])
        if per_half:
            lags[int(u.id)] = per_half
        if len(lags) == n_units:
            break
    if not lags:
        return dict(offset=np.nan, measured=np.nan, lags=lags, consistent=False)
    per_unit = [float(np.median(v)) for v in lags.values()]
    measured = float(np.median(per_unit))
    if sort.generations:
        offset = min((g * t0 for g in sort.generations), key=lambda o: abs(o - measured))
    else:
        offset = int(round(measured))
    stable = np.mean([max(v) - min(v) <= 3 for v in lags.values()]) > 0.5
    near = np.mean([abs(u - offset) <= 8 for u in per_unit]) > 0.5
    return dict(offset=int(offset), measured=measured, lags=lags, consistent=bool(stable and near))


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


def with_depth(tab: pd.DataFrame, alignment: pd.Series, n_contacts: int = 32) -> pd.DataFrame:
    """Add depth coordinates of each unit's peak contact from a laminar alignment row."""
    from .dataset import channel_table
    a = alignment
    ch = channel_table(n_contacts, a.orientation, float(a.l45), int(a.top_of_cortex),
                       float(a.deep_limit) if pd.notna(a.deep_limit) else np.nan).set_index("contact")
    tab = tab.copy()
    contact = tab.peak_contact.fillna(tab.contact).astype(int)
    for col in ("depth_l45", "depth_um", "depth_below_toc", "in_cortex", "beyond_deep_limit"):
        tab[col] = ch.loc[contact, col].to_numpy()
    return tab


def units(session: str, recording: str, bank: str | None = None, n_waveform: int = 300,
          n_amplitudes: int = 3000, seed: int = 0,
          sort: PhySort | None = None) -> tuple[pd.DataFrame, dict[int, np.ndarray], dict]:
    """Curated units of one sorted recording on the NSP clock, with isolation metrics.

    Returns (table, spike ticks per cluster id, offset report). The table has one row per phy 'good' or
    'mua' cluster (noise dropped) of `bank` (all banks if None): phy group/KS label, electrode label and
    contact, peak channel and waveform from the raw data, firing rate, isi_contamination, amplitude_cutoff,
    presence_ratio, snr, trough_to_peak_ms, and `well_isolated`."""
    sort = sort or load_sort(session, recording)
    ns6 = NsxFile(recording_files(session, recording)[".ns6"])
    rep = measure_offset(sort, ns6)
    if not rep["consistent"]:
        if not (sort.generations and len(sort.generations) == 1):
            raise RuntimeError(f"{session}_{recording}: spike-time offset not consistent: {rep}")
        rep["offset"] = sort.generations[0] * ns6.t_start      # single known convention (Kilosort 4: 0)
    ticks_all = base_ticks(sort, ns6) + rep["offset"]
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
                   ks_amplitude_cutoff=amplitude_cutoff(sort.amplitudes[sel]),
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
            # per-spike amplitude on the peak channel at the mean trough (sorter-independent amplitude cutoff)
            many = rng.choice(ok, min(n_amplitudes, len(ok)), replace=False)
            amps = -_snippets(ns6, many, [probe[pk]], window=(it + SNIPPET[0] - 2, it + SNIPPET[0] + 3)).min(axis=(1, 2))
            row.update(peak_label=probe[pk], peak_contact=int(probe[pk][2:]), trough_uv=float(wv[it]),
                       snr=float(troughs[pk] / nsd[probe[pk]]), trough_lag=it + SNIPPET[0],
                       trough_to_peak_ms=(ip - it) / 30.0, amplitude_cutoff=amplitude_cutoff(amps))
        rows.append(row)
    tab = pd.DataFrame(rows)
    if len(tab):
        tab["well_isolated"] = ((tab.isi_contamination < MAX_ISI_CONTAMINATION)
                                & (tab.amplitude_cutoff < MAX_AMPLITUDE_CUTOFF)
                                & (tab.presence_ratio >= MIN_PRESENCE_RATIO) & (tab.snr >= MIN_SNR))
    return tab, spikes, rep
