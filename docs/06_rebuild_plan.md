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
* Probe orientation for flagged sessions (05 C3, C4) is settled first, from
  the MUA/LFP profile at the cortex edges and from the CSD sink–source order.
* The clock fix matters here: a 3.4 ms spike-vs-LFP offset is a 12° phase error
  at 10 Hz and 61° at 50 Hz.

## Phase 4 — Units
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

## Phase 5 — Science (later)
Redo the LMN paper on the new data: laminar CRFs, NRM vs LSM, NDE weight,
supersaturation by layer as a test of SSN predictions, laminar timing.
Hierarchical statistics with penetration as a random effect; model comparison
with equal N and held-out data.

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
