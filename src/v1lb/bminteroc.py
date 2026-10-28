"""Per-presentation table for the bminteroc paradigm (6 x 6 left/right contrasts).

Stimulus-file column meanings (T_binocular_gratings.m, confirmed against the
explicit `_di_v2` files from 211208 on):
* grating_fixedc = LEFT-eye contrast, grating_contrast = RIGHT-eye contrast
  (ML eye code 3 = left, 2 = right).
* other_xpos / other_ypos = stimulus centre relative to fixation (dva) in every
  file version; `horzdva` and `grating_xpos` swapped meaning after 211008.
* Same orientation, SF and size in both eyes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import photodiode
from .config import recording_files
from .events import ml_codes, presentations, split_trials
from .io import NevFile, NsxFile, read_stim_text
from .trials import match_trials

# Every legacy STIM file has V1 = 'RV1' for monkey B and 'LV1' for monkey J.
SUBJECT_HEMISPHERE = {"B": "RV1", "J": "LV1"}

_STIM_COLUMNS = {
    "grating_fixedc": "contrast_L",
    "grating_contrast": "contrast_R",
    "grating_tilt": "tilt",
    "grating_sf": "sf",
    "grating_diameter": "diameter",
    "grating_phase": "phase",
    "other_xpos": "x",
    "other_ypos": "y",
}


def classify(contrast_l: np.ndarray, contrast_r: np.ndarray, hemisphere: str) -> np.ndarray:
    """blank / mono_contra / mono_ipsi / dioptic / dichoptic."""
    contra_is_left = hemisphere == "RV1"
    l, r = np.asarray(contrast_l), np.asarray(contrast_r)
    out = np.full(l.shape, "dichoptic", dtype=object)
    out[(l == 0) & (r == 0)] = "blank"
    out[(l > 0) & (r == 0)] = "mono_contra" if contra_is_left else "mono_ipsi"
    out[(l == 0) & (r > 0)] = "mono_ipsi" if contra_is_left else "mono_contra"
    out[(l > 0) & (l == r)] = "dioptic"
    out[np.isnan(l) | np.isnan(r)] = None
    return out


def presentation_table(session: str, recording: str) -> tuple[pd.DataFrame, dict]:
    """All stimulus presentations of one bminteroc recording, NSP clock.

    One row per onset code (completed or not), with event-code and photodiode
    times, the matched ML trial, stimulus parameters, condition labels and the
    previously *shown* presentation's contrasts.
    """
    files = recording_files(session, recording)
    nev = NevFile(files[".nev"])
    ticks, values, _ = nev.digital_events()
    codes = ml_codes(values)
    trials = split_trials(codes, ticks)
    pres = presentations(codes, ticks, trials)
    pres, levels = photodiode.detect(NsxFile(files[".ns6"]), pres)

    stim = read_stim_text(files[".gBMINTEROCGrating_di"])
    match = match_trials(trials, nev.origin_utc, stim)
    pres = pres.merge(match.table[["nev_trial", "ml_trial"]], on="nev_trial", how="left")
    s = stim.rename(columns={"trial": "ml_trial", **_STIM_COLUMNS})[["ml_trial", "pres", *_STIM_COLUMNS.values()]]
    pres = pres.merge(s, on=["ml_trial", "pres"], how="left", validate="one_to_one")
    pres = pres.merge(trials[["nev_trial", "rewarded", "broke_fixation"]], on="nev_trial", how="left")

    hemisphere = SUBJECT_HEMISPHERE[session.split("_")[1]]
    contra, ipsi = ("L", "R") if hemisphere == "RV1" else ("R", "L")
    pres["contrast_contra"] = pres[f"contrast_{contra}"]
    pres["contrast_ipsi"] = pres[f"contrast_{ipsi}"]
    pres["condition"] = classify(pres["contrast_L"].to_numpy(), pres["contrast_R"].to_numpy(), hemisphere)

    shown = pres[pres["pd_on_tick"] > 0].sort_values("pd_on_tick")
    prev = shown[["contrast_L", "contrast_R", "pd_off_tick"]].shift(1)
    pres.loc[shown.index, "prev_contrast_L"] = prev["contrast_L"].to_numpy()
    pres.loc[shown.index, "prev_contrast_R"] = prev["contrast_R"].to_numpy()
    pres.loc[shown.index, "ms_since_prev_off"] = (shown["pd_on_tick"] - prev["pd_off_tick"]).to_numpy() / 30

    meta = dict(session=session, recording=recording, hemisphere=hemisphere,
                origin_utc=nev.origin_utc.isoformat(), n_trials=len(trials),
                trial_shift=match.shift, ml_minus_nsp_s=match.offset_s, trial_match_max_dev_s=match.max_dev_s,
                pd_dark=levels.dark, pd_white=levels.white, pd_frame_ticks=levels.frame_ticks)
    return pres, meta
