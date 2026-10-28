# v1-laminar-binocular — notes for Claude

Blake Mitchell's rebuild of his PhD pipeline (laminar V1, binocular contrast).
**Read `docs/README.md` first.** `docs/05` is the live issue register and
`docs/06` is the plan; keep both current as things change.

## Decisions (2026-10-28)
- MATLAB for all pipeline code. Per-session outputs are `.mat` files (v7.3 when large; use `matfile` for partial loads).
- Read every stream by **NSP timestamp**. The legacy NSx reads were offset by the `.ns6` start timestamp (102 samples).
- Photodiode: t = 0 is the first backlight flash at the photodiode, detected per presentation on the NSP clock (~20.0 ms after the event code). The ViewPixx backlight mode is unknown, so the scan offset to the stimulus row (≈2.6–3.4 ms if the backlight scans) is stored as a separate per-session field and **not** applied by default.
- Laminar alignment: spike-phase coupling (Davis et al. 2023: SPI and preferred phase, plus PPC). CSD is only a cross-check.
- Units: start from the 2023 phy curation with a fixed importer. A Kilosort 4 pilot is optional.
- Correct legacy errors; don't reproduce old numbers for their own sake.

## Data (read-only — never write, move or rename anything under `D:\NEUROPHYS DATA`)
- Raw: `D:\NEUROPHYS DATA\2021\Raw Data B\`, `Raw Data J\`. Dual-probe J days have 72 channels (eA01–32 + eB01–32) in the same files.
- Legacy processed data: `STIM 2024 2\` (has wrong channels for 211012/211102/211217, see 05 C3); phy curation in `KiloSort-ed_phy\`.
- Legacy session metadata: `docs/legacy_TuneList_all_versions.csv`.
- Task code: `D:\MATLAB\NIMH-MonkeyLogic\33 binocular gratings\`. Legacy helpers: `D:\MATLAB\MLAnalysisOnline`, `KiloSortUtils-master`.
- New processed outputs go to a separate output root, not inside the raw tree (location still to be decided).

## Code conventions
- No dependencies on the legacy repo, `MLAnalysisOnline` or `archive/`. Port the few functions that are needed (`parsEventCodesML`, `readgGrating`, …) with attribution and tests.
- No hard-coded user paths. One config function returns the roots.
- Keep functions small and testable. Session-level functions take a penetration ID and return structs/tables; nothing depends on the base workspace.

## Environment
- MATLAB R2026b is installed at `C:\Users\bmitc\AppData\Local\Programs\MATLAB\R2026b` (via mpm, with Signal, Stats/ML, Curve Fitting, Optimization, Image Processing). It is **not activated yet**: it needs Salk's network license (`licenses\network.lic` or `MLM_LICENSE_FILE`). Once active, run headless with `"<root>\bin\matlab.exe" -batch "..."`.
- Python 3.14 (numpy/scipy/pandas/openpyxl, user site, so don't use `python -I`) is fine for throwaway read-only checks.
- Scratch scripts go in the session scratchpad, never in this repo or on D:.

## Git
- Private repo, default branch `main`. Commit and push routine doc updates without asking. Ask before pushing code changes until a code workflow is agreed. Never force-push.
