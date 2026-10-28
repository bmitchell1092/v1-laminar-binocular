# 05 — Known issues, open questions, and the verification plan

Severity: **H** = could change a scientific conclusion; **M** = bias or
robustness problem; **L** = fragility or cosmetic.
Status: *measured* (checked on raw data 2026-10-28), *code* (follows from
reading code; not yet tested on data), *open* (unknown).

---

## A. Timing and alignment

| # | Issue | Sev | Status |
|---|---|---|---|
| A1 | `.ns6` data start at NSP timestamp 102 and `.ns2` at 0. NPMK windowed reads plus the "keep trailing samples" trim leave `tp_pt` 102 samples (3.4 ms) early. MUA and KLS cancel this; LFP/CSD, NEV and autosort spikes land **+3.4 ms late** relative to MUA/KLS. ([01 §6](01_preprocessing_pipeline.md)) | H for latency, L for magnitude | T0 *measured*; consequence *code* |
| A2 | Legacy onset = photodiode crossing − 1 full refresh (8.33 ms) − 102 samples (A1). **Measured** on `211008_B_bminteroc003` and `221206_J_bminteroc002` (`ainp1`, NSP clock): (1) the ViewPixx backlight is **strobed**: the white box gives a ~3.5 ms light pulse once per 8.33 ms frame. (2) The first pulse reaches the photodiode **20.0 ms after the event code** (50% rise; IQR 0.04 ms, B and J alike), with a ~1.9 ms 10–90% rise. (3) Stored legacy `tp_pt − tp_sp` = 265 / 255 samples, against 250 / 247 predicted from the PD time − 352. The mechanism is confirmed to within 0.5 ms. Net effect: legacy t = 0 is **≈ 8.3 ms before the photodiode pulse for MUA/KLS** and **≈ 11.7 ms before it for LFP/CSD/NEV**. Legacy latencies are therefore too long by that amount minus any scan offset (below). | **H** for latency | *measured* |
| A2b | Light at the stimulus location vs at the photodiode (box at bottom-right, stimulus near center): if the backlight **scans** top→bottom, the stimulus is lit earlier by Δ ≈ (row_PD − row_stim)/1125 × 8.33 ms ≈ **3.4 ms (B: PPD 40.9, ypos −1.3°)** and **≈ 2.6 ms (J 221206: PPD 48.6, ypos −3.47°)**. If it is a **global** strobe, Δ = 0. One photodiode can't tell these apart; need the ViewPixx backlight mode used on the rig. | M | *open* |
| A3 | Legacy detection (|z| > 1 + "energy" rule, applied to a pulsed signal) mis-times **12/398 (B) and 15/1450 (J)** presentations by 6–15 ms early or ~1 frame late. Offsets are ill-defined: legacy offset ≈ 4.5 ms *before* the off code, whereas the light actually stops ~20 ms after it. That makes legacy stimulus durations ~11–13 ms too short. Rebuild: detect the first pulse in a fixed 15–28 ms post-event window, flag deviations > 1 frame from the session median (real dropped frames), and define offset as the first missing pulse. | M | *measured* |
| A4 | 16 recordings have two data packets (7–30-sample gaps). NPMK treats them as paused, with no zero-padding, and windows can straddle segments. Includes `211027_B_bminteroc001`, `221103_J_bminteroc001`, `221128_J_bminteroc001`. | M | *measured* |
| A5 | Stimulus duration differs by session: event-code on→off is 216.7 ms (26 frames) on 211008 vs 316.7 ms (38 frames) later. The ISI is 316.7 ms. Both on and off flips land ~20 ms after their codes, so the true on-duration equals the event-code interval. SDF lengths and window clipping differ across sessions. | M | *measured* |
| A6 | 3 presentations per trial with about a 317 ms ISI. The baseline (−50…0 or −150…−50 ms) for presentations 2–3 sits in the ISI after the previous stimulus, so offset responses and contrast adaptation from the previous presentation are in play. Presentation history is not modelled anywhere. | M | *code* |

## B. Trials, stimulus and eye coding

| # | Issue | Sev | Status |
|---|---|---|---|
| B1 | The NEV ↔ text-file match is by index, and the `checkTrMatch` failure was demoted to a warning. Verified OK on 4 files (trial counts equal, first trial 1, 3 rows/trial); must be checked on all. | H if any file fails | partly *measured* |
| B2 | Trial-ID codes (116–125) are absent in early-B files such as 211008, so only the timestamp check is available there. | L | *measured* |
| B3 | Contrast pairs are drawn **with replacement** (`randi`), so per-condition trial counts are unbalanced. | M | *code* |
| B4 | Raw eye 3 = **left**, 2 = **right**. `diTP` swaps for LV1 so 2 = ipsi, 3 = contra. The OD project separately flips **B's** OD sign with hard-coded unit ranges. Possible double flip or inconsistency between projects. | H for OD/DE assignment | *code* |
| B5 | 211008 differs from later sessions: SF 4 cpd, 1° diameter, header `cinteroc`, and a different text-file column layout. | M (heterogeneity) | *measured* |

## C. Laminar alignment

| # | Issue | Sev | Status |
|---|---|---|---|
| C1 | ToC / SinkBtm / BoC were read **by eye** from CSD, with no algorithm. Between STIM generations, SinkBtm moved by up to 6 channels (600 µm): `211008_B_eA` went 24 → 22 → 18; `221204_J_eB` went 16 → 29 → 22. Layer assignment is the least reproducible element. | **H** | *measured* |
| C2 | The 2021–22 DataLog with ToC/BoC columns is missing. **Recovered**: every version from the 9 `ERR.mat` files plus the Feb-2022 `DataLog_2021.xlsx` from git (commit `574f625^`) is in [legacy_TuneList_all_versions.csv](legacy_TuneList_all_versions.csv) (143 rows). The J values were stable from Jan 2023 on; the B values changed in every generation. | M | *measured* |
| C3 | **Wrong channels in `STIM 2024 2` (= `LMNProject 2024 2`, the LMN input) for 211012, 211102 and 211217.** On 2024-03-07 they were re-marked `descending`/`Flipped` with limits in flipped numbering (1/5/17, 1/5/17, 1/11/20 ≈ raw ch 32→16). `runTuneList` doesn't remap, so it extracted raw `eA01–eA17`/`eA01–eA20`. One day earlier (`STIM 2024`) 211012 used `eA15–eA32`. So 3 of the 11 B penetrations in the LMN dataset are mostly the wrong contacts, with reversed depth labels. The true probe orientation for these sessions is still unknown. | **H** | *measured* |
| C4 | Channel pin order differs for 211007–211012. The label → physical-contact mapping needs checking via CSD. | M | *open* |
| C5 | `BadCh` (`STIM.rmch`) is stored but not applied for 2021 data. | M | *code* |
| C6 | At least 4 different supra/gran/infra rules, and plot labels say mm where the units are channels. | M | *code* |
| C7 | CSD: no spatial smoothing, no boundary padding (edge channels NaN); the bminteroc CSD RESP keeps sinks only. Spacing of 0.1 mm matches the Plexon 100 µm probes. | L | *code* |

## D. Units and spikes

| # | Issue | Sev | Status |
|---|---|---|---|
| D1 | `importPhy` drops good units that share a channel: it loops over unique depths but indexes clusters by the loop counter. | H for SUA results | *code* |
| D2 | In `diNeuralDat_usePhy`, trials from files that were not sorted get firing rate **0**, not NaN. | H if analyses mix files | *code* |
| D3 | Kilosort ran per file; several bminteroc files failed or were never sorted (`D:\KiloSort\KiloSort_log.xlsx`). The `.ns6.dat` inputs are gone, so re-sorting requires re-binarizing. | M | *measured* |
| D4 | The LMN "units" are **MUA channels** (100 µm apart, overlapping signal), and statistics treat them as independent. Pseudoreplication. | **H** for inference | *code* |
| D5 | The KS `chanMap.mat` for `221206_J_bminteroc002` references `bminteroc001.ns6`. | L | *measured* |

## E. Analysis code (details in [04 §6](04_LMN_analysis.md) and [03 §6](03_project_history_and_conventions.md))

| # | Issue | Sev |
|---|---|---|
| E1 | `getPARAMS_1D` (PARAMS / PARAMS_s used by the supersaturation figure) exists nowhere; several functions are defined under other file names. | H (not reproducible) |
| E2 | Saved workspaces were made by hand; `depth_array` columns 4–5 (used by v2 Figs 2D/3D/4A/4B) are created by no code. | H (not reproducible) |
| E3 | AIC is compared across fits with different N (6 vs 5 points); v1 and v2 fits are mixed in Fig 2C. | M |
| E4 | OD/DE is chosen from the sign of `occ_norm` on the same data that are then analysed, with no significance test (circular). | M |
| E5 | Inclusion is only a contrast-tuning ANOVA, with no visual-responsiveness or baseline test. | M |
| E6 | RESP and SDF are normalized by different constants; DI is rescaled separately in the CSD variant. | M |
| E7 | Indices averaged over contrasts include the blank row. | M |
| E8 | ORGaNICs simulation: Euler Δt = 1 ms against τ = 0.1–0.5 ms (unstable); `Ntheta` name mismatch; no saved fits. | M |
| E9 | Hard-coded unit IDs, row ranges, depth crops, 6-contrast ordering, paths. | L |
| E10 | Window index meaning depends on `resptype` (13 vs 50 windows); several scripts assume tw 12 = full. | M |

## F. Provenance and environment

| # | Issue |
|---|---|
| F1 | No working MATLAB on this machine (R2019b/R2021b folders lack `bin`). Python 3.14 + numpy/scipy/pandas/openpyxl are available; no h5py. |
| F2 | External dependencies (MLAnalysisOnline, KiloSortUtils, npy-matlab, ML2) are unversioned copies under `D:\MATLAB`; `setup.m` expects them under `Documents\MATLAB`. |
| F3 | LMN workspace (`OneDrive\...\workspaces\LMN Project`, ~25 GB) is mostly OneDrive cloud-only. `CSD_bminteroc.mat`, `MUA_bminteroc.mat` and `LFP_bminteroc*.mat` were downloaded on 2026-10-28 by a header read. |
| F4 | B sessions in 2022 (`220202_B`, OD project "211008–220302") are referenced but not on `D:`. |
| F5 | Two OneDrive roots in the code paths (`OneDrive - Vanderbilt` and `OneDrive\4. Vanderbilt`). |

---

## Verification plan (Phase 1). All read-only, scripted, one report per session.

**§A Clock and photodiode**
1. Read `ainp1` from `.ns6` on the **NSP clock** (file sample k ↔ 102 + k). Read event times from the NEV. Per presentation, measure event → photodiode rise and fall latency; plot the distribution per session.
2. Compare against stored `STIM.tp_pt` to confirm the predicted offset (−102 − 250 samples relative to the photodiode crossing) exactly.
3. From the photodiode waveform shape (ramp over a frame?), decide whether the panel scans and what correction to apply. Use the stimulus y position from the text file, the screen geometry from bhv2 `MLConfig`, and the photodiode box at the lower right.
4. Measure the true on-duration per presentation and session.
5. Characterize the two-packet files.

**§B Trials and stimulus**
1. For every bminteroc file: NEV trial count vs text-file trials, the 116-code trial digits where available, ML vs BR timestamps, and the bhv2 `TrialError` vs completed presentations.
2. Rebuild the per-presentation table independently of `diTP` and diff it against `STIM 2024 2`.
3. Tabulate per-condition trial counts (L × R contrast) per session.
4. Confirm the eye mapping from `.gBMINTEROCGrating_di_v2` (explicit `contrast_L`/`contrast_R` columns, from 211208) against `fixedc`/`contrast`.

**§C Laminar**
1. Recompute evp-flash and bminteroc CSD per penetration with a documented method.
2. Detect the initial sink automatically (onset/bottom) and compare with every historical v1lim generation.
3. Check sort direction and pin order for flagged sessions.

**§D Units**
1. Re-import phy `good` clusters for all 24 sorted recordings with a fixed importer; count units vs the old `STIM.units`.
2. Quantify the units lost to D1.
