# v1-laminar-binocular — notes for Claude

Blake Mitchell's rebuild of his PhD pipeline (laminar V1, binocular contrast).
**Read `docs/README.md` first.** `docs/05` is the live issue register and
`docs/06` is the plan; keep both current as things change.

## Decisions
- **Python** for all pipeline code (2026-10-28; MATLAB could not be licensed on this machine).
- Storage: per-penetration files readable without this code. Proposed: labelled arrays as HDF5/NetCDF via xarray (MATLAB can still open them with `ncread`/`h5read`), plus Parquet/CSV tables. *Awaiting Blake's confirmation.*
- Read every stream by **NSP timestamp**. `.ns6` data start at NSP sample 102 and `.ns2` at 0 (neo agrees). The legacy reads were offset by this.
- Trials are matched to stimulus-file rows **by timestamp**, never by index (see 05 B1: 211009_B).
- Photodiode: t = 0 is the first backlight flash at the photodiode, detected per presentation on the NSP clock (~20.8 ms after the event code, but sometimes one frame earlier or later; see 05 A2c). The ViewPixx backlight mode is unknown, so the scan offset to the stimulus row (≈2.6–3.4 ms if the backlight scans) is stored per session and **not** applied by default.
- Laminar alignment: spike-phase coupling (Davis et al. 2023: SPI and preferred phase, plus PPC). CSD is only a cross-check.
- Units: start from the 2023 phy curation with a fixed importer. A Kilosort 4 pilot is optional.
- Correct legacy errors; don't reproduce old numbers for their own sake.

## Data (read-only — never write, move or rename anything under `D:\NEUROPHYS DATA`)
- Raw: `D:\NEUROPHYS DATA\2021\Raw Data B\`, `Raw Data J\`. Dual-probe J days have 72 channels (eA01–32 + eB01–32) in the same files.
- Legacy processed data: `STIM 2024 2\` (wrong channels for 211012/211102/211217, wrong labels for 211009; see 05). MAT v7.3, readable with h5py. Phy curation in `KiloSort-ed_phy\`.
- Legacy session metadata: `docs/legacy_TuneList_all_versions.csv`.
- Task code: `D:\MATLAB\NIMH-MonkeyLogic\33 binocular gratings\`. Legacy MATLAB helpers (reference only): `D:\MATLAB\MLAnalysisOnline`, `KiloSortUtils-master`.
- New processed outputs go to a separate output root, not inside the raw tree (location still to be decided).

## Code conventions
- Package code in `src/` (installable package); tests in `tests/` (pytest), with small fixtures and an optional marker for tests that need `D:`.
- Own small, explicit readers for NSx/NEV/stimulus text/bhv2, cross-checked against neo in tests. Every time value carries its clock (NSP samples) in the name or docstring.
- No dependencies on the legacy repo or MATLAB code. Port needed logic (`parsEventCodesML`, `readgGrating`, …) with attribution and tests.
- No hard-coded user paths. One config module returns the roots.
- Session-level functions take a penetration ID and return tables/arrays; no hidden global state.

## Environment
- Project venv: `.venv\Scripts\python.exe` (Python 3.14, packages in `requirements.txt`). Use it for all repo code and tests.
- The system Python (3.14, user site-packages) is fine for throwaway checks; `python -I` hides its packages.
- MATLAB R2026b sits unactivated at `C:\Users\bmitc\AppData\Local\Programs\MATLAB\R2026b` and is not used.
- Scratch scripts go in the session scratchpad, never in this repo or on D:.

## Git
- Private repo, default branch `main`. Commit and push routine doc updates without asking. Ask before pushing code changes until a code workflow is agreed. Never force-push.
