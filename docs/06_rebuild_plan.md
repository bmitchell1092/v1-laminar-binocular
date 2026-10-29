# 06 — Rebuild plan

Goal: every session can be opened on its own, with an indexed trial table and
every signal (LFP, CSD, MUA, sorted units) on one clock. QC is built in, and
pooled analyses are derived by code, never from hand-saved workspaces. Old
results are corrected wherever they turn out to be wrong.

## Decisions (Blake, 2026-10-28)
| topic | decision |
|---|---|
| language | **Python** (changed from MATLAB on 2026-10-28: the institutional license couldn't be set up on this machine). Project venv `.venv`, `requirements.txt` |
| storage | per-penetration HDF5/NetCDF (xarray, labelled dims) + Parquet/CSV tables; still openable in MATLAB; NWB export at publication |
| workflow | commit straight to `main`, small commits with passing tests |
| clock | fix the NEV↔NSx offset: every stream is read by NSP timestamp |
| photodiode | t = 0 = first backlight flash at the photodiode, per presentation, NSP clock; the scan offset to the stimulus row is stored per session but not applied (backlight mode unknown; see 05 A2b) |
| layers | **spike-phase coupling**, following Davis et al. 2023 (SPI, i.e. mean resultant length, and preferred phase), plus PPC; replaces by-eye CSD calls and spectrolaminar crossings (CSD kept only as a cross-check) |
| units | keep the 2023 phy curation as the baseline and fix the importer; Kilosort 4 is a candidate (pilot first) |
| paper | redo all analyses later on the new data |
| missing 2022 B sessions | non-critical to-do |

## Phase 0 — Understand and document *(done)*
`docs/01–05`, plus the recovered TuneList (`legacy_TuneList_all_versions.csv`).

## Phase 1 — Verify raw-data facts (read-only)
Run 05 §A–§D. Deliverable: one QC page per session plus a summary table.
* §A: photodiode on the NSP clock (rise/fall latency, waveform shape, scan
  geometry from bhv2 `MLConfig` and stimulus y-position) gives the correction;
  measure true stimulus durations.
* §B: independent trial table for every bminteroc file, diffed against `STIM 2024 2`.

## Phase 2 — Session-level data model (Python)
Per penetration `YYMMDD_M_eX`:
* `trials` table, one row per presentation: file, ML trial, presentation index,
  `contrast_L/R`, raw eye, ipsi/contra, orientation, SF, size, x/y, NEV on/off,
  photodiode on/off (NSP samples), duration, previous presentation, outcome, QC flags.
* `channels` table: label, probe/bank, raw channel, physical contact index
  (orientation resolved), depth µm, depth relative to the laminar reference,
  layer, bad flag, alignment version.
* signals on the NSP clock: LFP (1 kHz), CSD, analog MUA (1 kHz, one documented
  filter chain), threshold-crossing MUA spike times, sorted-unit spike times
  plus metrics.
* sidecar paradigms: evp, darkrest, dotmapping, rfori, bmcBRFS.
* One version-controlled session table (CSV) replaces DataLog.xlsx: file
  numbers per paradigm, probe, bank, hemisphere, Rank, bad channels, notes.

## Phase 3 — Laminar alignment by spike-phase coupling
Method (Davis, Dotson, Franken, Muller & Reynolds, *eLife* 2023, 12:e84512;
code: github.com/zwdsalk/LaminarPhaseCoupling, generalized phase:
github.com/mullerlab/generalizedphase). Their code is MATLAB, so we port it
to Python and validate the port on their example data (Dryad
doi:10.5061/dryad.3r2280gmm) before using it on ours:
* MUA spikes: high-pass 500 Hz, threshold at 4 SD of a sliding 1 s window, per channel.
* LFP: low-pass 300 Hz, 1 kHz, wideband 5–50 Hz (4th-order Butterworth,
  zero-phase), generalized phase; narrowband Hilbert (4–8, 8–14, 15–30, 30–50 Hz)
  as a robustness check.
* For each spike channel i × LFP channel j: SPI = |mean e^{iφ_j(t_spk,i)}|,
  preferred phase θ_ij, and PPC (Vinck 2010), which is unbiased by spike count.
* Input/deep boundary = depth of the preferred-phase reversal. The paper picks
  it by eye; we will make it algorithmic (e.g. the split that maximizes
  circular separation of θ above vs below) and keep a manual review step.
* Data: ≥1 min continuous recording (darkrest, fixation, or full task files);
  one estimate per file gives a drift check within a session.
* **Marker hierarchy (Blake, 2026-10-28).** SPI as a function of depth is the
  final arbiter of **probe orientation** and the input/deep boundary. The
  bottom of L4 (L4/L5 boundary) is harder to read from PPC alone, so it is
  cross-checked with two supporting markers: the **spectrolaminar
  crossing** (alpha-beta vs gamma relative power; Mendoza-Halliday et al.
  2024) and the **earliest CSD sink** (stimulus-evoked).
* Orientation for the flagged sessions (05 C3, C4) is settled first, by SPI
  vs depth.
* The clock fix matters here: a 3.4 ms spike-vs-LFP offset is a 12° phase error
  at 10 Hz and 61° at 50 Hz.

## Phase 4 — Units (**mandatory for the supersaturation paper**, Blake 2026-10-28: MUA pools neurons with different c50s and can mask single-neuron supersaturation; well-isolated Kilosort units are required)
Fix `importPhy`: all good units, every sorted file, NaN for unsorted files.
**Kilosort 4 pilot** (optional, after Phase 1):
* For: only 24 recordings were sorted, several bminteroc files failed, and
  per-file sorting means units can't be followed across paradigms. KS4 on
  concatenated session files (bminteroc + rfori + dotmapping + darkrest) would
  give each unit its own tuning, RF and laminar position. Drift correction is better too.
* Against: re-binarizing about 96 h of `.ns6`; KS4 runs in Python (isolated
  environment, GPU present); re-curating in phy is the slow part.
* Plan: sort 2 sessions with KS4, compare yield and quality with the 2023
  curation, then decide.

## Phase 5 — Science (priorities set by Blake, 2026-10-28)
Not a full redo of the old LMN paper.
1. **Supersaturation as a function of depth and/or laminar compartment** (the
   paper). Start with the best days and assess power first (per-channel
   reliability of supersaturation, between-penetration variance, simulated
   power for a depth or compartment effect with penetration as a random effect).
   If laminar assignment leaves it underpowered, fall back to **depth from the
   top of cortex**, with V1 channels defined by dotmapping receptive fields (no
   layer labels).
2. **Gap-filling analyses** chosen after a dense, up-to-date literature sweep.
3. **Dual-probe days** (7 J sessions, eA + eB recorded together): response
   delays between columns, and traveling waves across V1.
Statistics: hierarchical, with penetration (and monkey) as random effects.

## To-do (non-critical)
* Locate the 2022 B sessions (`220202_B`, `220302_B`).
* Locate the final 2021–22 DataLog.xlsx (its contents are mostly recovered already).

## Progress log
* 2026-10-28: `v1lb` package started: NSx/NEV readers on the NSP clock (identical to
  neo, sample for sample), stimulus-file reader, event-code trial and
  presentation parsing, and timestamp-based trial matching (detects the
  211009_B one-trial shift). 16 tests pass, including raw-data cross-checks.
* 2026-10-28: `v1lb.photodiode`: per-presentation flash onset, frame count and
  offset (first missing flash). Validated on all 35 bminteroc files (37,167
  completed presentations, none missed). 19 tests pass.
* 2026-10-28: `v1lb.bminteroc.presentation_table` (events + photodiode +
  timestamp matching + stimulus parameters + ipsi/contra labels + previous
  presentation) and `v1lb.legacy` (reads legacy STIM v7.3). Rebuilt tables
  match all 30 legacy penetrations except the known errors (211009 labels;
  legacy t0 ~12 ms early, and 1-3 frames off on 7-29% of presentations in
  five B penetrations). 21 tests pass.
* 2026-10-28: `v1lb.signals` (LFP read, threshold-crossing MUA per Davis et al.)
  and `v1lb.laminar` (generalized phase port, SPI/preferred phase/PPC matrix,
  spectrolaminar relative power, CSD + amplitude-based earliest sink). Pilot on
  211008_B_eA, 211012_B_eA and 221206_J_eB (results in 05 C8). 26 tests pass.
* 2026-10-28: phase convention validated against Davis et al.'s published
  outputs (supra/input ~ ±π, deep ~ −1.2 rad); their average phase and
  gamma−alpha profiles saved as templates (`src/v1lb/data`, with provenance).
  Added `fit_template` (both orientations, full-probe overlap) and peak-based
  `spectral_crossing`. Their code calls generalized phase with a 0 Hz cutoff,
  which we match. V1 profiles can show two phase steps, so the boundary rule
  needs review (05 C8). 28 tests pass.
* 2026-10-28: laminar decision rules agreed with Blake and implemented
  (`laminar.top_of_cortex`, `phase_steps`, `step_kind`, `l45_from_steps`,
  `mua_gap`; `laminar_report` + `scripts/run_laminar_reports.py`):
  orientation from SPI (template fit) checked against spectrolaminar peaks;
  top of cortex from spontaneous MUA (evoked analog MUA only as fallback);
  L4/5 = 2nd non-reversed phase step below the top of cortex, else the only
  one; two steps may also be two folds (V1/V1 or V1/V2), so reversed steps
  (>= 0.8 rad) or MUA gaps between steps are flagged. Each penetration gets a
  complexity label and a confidence tier (A/B/C) from 8 explicit checks, for
  ranking and possible exclusion. `metadata/penetrations.csv` lists the 30
  penetrations. 30 tests pass.
* 2026-10-28: laminar reports for all 30 penetrations
  (`docs/qc/laminar_calls_2026-10-28.csv`): 3 tier A, 12 B, 15 C; markers agree
  on orientation in 22/30. In 5 penetrations both markers contradict the legacy
  sort direction (211027_B, 221109_J, 221128_J, 221209_J_eA, 221214_J_eA);
  these are likely second-fold cases that need manual review. Review page
  (private): https://claude.ai/artifact/3td44k5uBfgQtTznuB3YQN. Verdicts are
  stored in its db under `verdicts/<penetration>` ({verdict, note, updated}).
  Built by `scripts/build_laminar_review.py`.
* 2026-10-28: **laminar alignment v1 frozen** (`metadata/laminar_alignment_v1.csv`,
  built by `scripts/freeze_laminar_alignment.py`). Blake: keep tiers A/B, set C
  aside. Claude reviewed A/B and applied documented exceptions: 211208_B L4/5
  overridden to 13.5 (second fold beyond contact 18); set aside 221128_J
  (orientation conflicts with legacy, ambiguous depth), 211103_B (contradictory
  candidates), 221209_J_eB (cross-checks match neither step), 221214_J_eB (L4/5
  3.5 contacts below top of cortex). **Included: 11 penetrations (J 8, B 3).**
  Depth of contact c relative to L4/5 (contacts, + = superficial): ascending
  l45 − c, descending c − l45. `spiking_bottom` is advisory only (small-spike J
  probes under-threshold deep layers); `deep_limit` is set only for confirmed
  second folds. Verdicts mirrored on the review page (keep / unsure).
* 2026-10-28: `v1lb.dataset` builds one xarray/NetCDF file per penetration
  (`output/datasets/<pen>.nc` + `_presentations.parquet`): lfp, csd, mua
  (analog envelope), spikes (threshold-crossing counts), presentation × time
  (−200…499 ms, 0 = photodiode flash, 1 kHz) × channel (superficial → deep),
  with the full trial table and depth relative to L4/5 as coordinates.
  `open_penetration(pen)` loads lazily. Sanity check on 211008: MUA onset
  35–37 ms at +3…+12 contacts above L4/5, initial sink at +2…+4, CRFs show
  binocular summation at low contrast and suppression at high contrast.
  Built by `scripts/build_datasets.py` for the 11 included penetrations.
* 2026-10-28: all 11 datasets built (723–1,313 presentations each, ~2.4 GB
  total). Physiology check (`docs/qc/dataset_check_2026-10-28.csv`): MUA latency
  consistent; earliest-sink cross-check from bminteroc is unstable (05 C9).
  Next: evp flash CSD as the proper CSD check of alignment v1.
* 2026-10-28: flash-evoked CSD (`v1lb.evp`, `scripts/run_flash_csd.py`,
  `scripts/compare_flash_csd.py`) for every evp file. First comparison with
  alignment v1 was confounded: the evp files were recorded 1–3 h before
  bminteroc while the probe was still being lowered (221119: flash sink moves
  ~22 contacts between two evps 27 min apart; 211008 drifts +20 → +17 → +6 → +1
  over the day). Clean agreements: 211008 (+1), 221202_J_eB (+2). 211208 is stable
  across 4 evps and puts the sink bottom 5 contacts deeper than the override,
  matching the original 2022 value (18). Added `v1lb.registration` (log MUA rate +
  gamma/alpha relative-power profiles) to register each evp to the bminteroc
  recording before comparing (`scripts/register_flash_csd.py`).
  `events.ml_events` drops non-strobe port lines (stray 2048 events in
  211008_B_evp002).
* 2026-10-28: eye data added (`v1lb.io.bhv2` MonkeyLogic 2 reader, `v1lb.eye`): `ainp3`/`ainp4` = gaze X/Y calibrated against MonkeyLogic, `ainp5` = pupil (light reflex); per-presentation `eye_x`, `eye_y`, `pupil` in the datasets plus continuous `<pen>_eye.nc` (05 C11). Next for eyes: saccade/microsaccade detection and fixation-vs-post-saccade epochs. 38 tests pass.
* 2026-10-28: laminar alignment v2 frozen (05 C12, C13). All 70 evps registered to bminteroc (probe static). Raw flash CSDs stripe from contact gain mismatch; outlier-contact interpolation + 3-point Hamming smoothing stabilises the landmark. Flash and SPI rarely agree, so L4/5 needs corroboration (anatomical depth below the top of cortex or SPI agreement). 22 penetrations with layer labels (1 A, 21 B), 26 depth-usable. Datasets now built for all depth-usable penetrations, with `depth_below_toc`, `bad_lfp` and the smoothed CSD. 40 tests pass.
