"""Depth registration between two recordings of the same probe (e.g. an evp file vs bminteroc).

The probe was often still being lowered while the evp files were recorded, so a landmark
measured in one file must be shifted into the other's contact numbers before comparing.
Profiles that depend on tissue rather than on the stimulus are lined up: log spontaneous
MUA rate and the gamma (75-150 Hz) and alpha-beta (10-19 Hz) relative-power profiles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import laminar as lam
from .config import recording_files
from .io import NsxFile
from .signals import probe_labels, read_lfp, threshold_spikes

TICKS = 30000


def depth_profiles(session: str, recording: str, bank: str, max_minutes: float = 5.0,
                   skip_s: float = 5.0) -> dict[str, np.ndarray]:
    """Per-contact profiles (contact-number order): log10(rate + 1), gamma and alpha-beta relative power."""
    f = recording_files(session, recording)
    ns2, ns6 = NsxFile(f[".ns2"]), NsxFile(f[".ns6"])
    labels = probe_labels(ns2, bank)
    start = ns6.t_start + int(skip_s * TICKS)
    stop = min(ns6.t_stop, start + int(max_minutes * 60 * TICKS))
    spikes = threshold_spikes(ns6, labels, start, stop)
    rate = np.array([s.size for s in spikes]) / ((stop - start) / TICKS)
    lfp, _ = read_lfp(ns2, labels, start, stop)
    fr, rp = lam.relative_power(lfp, 1000)
    return dict(log_rate=np.log10(rate + 1), gamma=lam.band_profile(fr, rp, (75, 150)),
                alpha_beta=lam.band_profile(fr, rp, (10, 19)))


@dataclass
class Registration:
    shift: int          # contact c in the moving recording sits at contact c + shift in the reference
    score: float        # mean Pearson r across profiles at the best shift
    margin: float       # best score minus the best score at least 3 contacts away
    n_overlap: int


def register(moving: dict[str, np.ndarray], reference: dict[str, np.ndarray], max_shift: int = 24,
             min_overlap: int = 12, keys=("log_rate", "gamma", "alpha_beta")) -> Registration:
    n = len(reference[keys[0]])
    scores = {}
    for s in range(-max_shift, max_shift + 1):
        lo, hi = max(0, -s), min(n, n - s)            # moving contacts c (0-based) with c + s inside the reference
        if hi - lo < min_overlap:
            continue
        rs = []
        for k in keys:
            a, b = moving[k][lo:hi], reference[k][lo + s:hi + s]
            if np.std(a) > 0 and np.std(b) > 0:
                rs.append(np.corrcoef(a, b)[0, 1])
        if rs:
            scores[s] = (float(np.mean(rs)), hi - lo)
    best = max(scores, key=lambda s: scores[s][0])
    far = [v[0] for s, v in scores.items() if abs(s - best) >= 3]
    return Registration(best, scores[best][0], scores[best][0] - (max(far) if far else np.nan), scores[best][1])
