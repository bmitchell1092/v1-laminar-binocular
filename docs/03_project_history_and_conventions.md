# 03 — Project history and conventions

Status: written 2026-10-28 from `3. visualization/` (all but 2024_LaminarCRFs,
which is in [04](04_LMN_analysis.md)), `utils/dichoptic/`, `archive/`.
[I] = inference.

## 1. Chronology

| # | Era | Question | Data | Code |
|---|---|---|---|---|
| 0 | 2012–2016 lab metadata | — | helga, bohr, einstein (E), ingrid (I) | `archive/jnm_sessioninfo.m` (J. Westerberg; pDep = L4c/5 boundary in channels, pSpc = 100 µm) |
| 1 | ~2016–2019, M. Cox era | monocular vs binocular vs dCOS, adaptation/SOA, attention (RSVP) | E, I (`diSTIM_*`, `IDX_*`) | `utils/dichoptic/diIDX*.m`, `diCRF`, `biCRF`, `diATTN` |
| 2 | 2020–2022, **iScience 2022** (Mitchell et al.) | dioptic vs monocular across contrast & time: sublinear summation, dynamic binocular facilitation, response- vs contrast-gain | `bincontrast_12101_08-Oct-2020-250ms`, datatype `auto` (dMUA); 234 E + 80 I = 314 units (Figs 1–2), 248 with ≥4 contrasts (Figs 3–5); E contrasts [0 .225 .45 .9] | `3. visualization/2021_iScience/`, `GCProjectinit.m`; windows = 16 × 100 ms sliding 10 ms + full (0–250) + baseline |
| 3 | 11/2020–5/2021 | binocular normalization as CRF shifts | iScience data | `2022_OD Project/Binocular_normalization/` (`modelcrf*`) |
| 4 | VSS 2022 | interocular normalization for balanced/unbalanced contrasts | early B data [I] | `2022_VSS/` |
| 5 | SfN 2022 | dichoptic SDF grids by layer; LSM/QDS/NLS/NRM over time | B [I] | `2022_SfN/` (`fitMODEL2` missing) |
| 6 | 2023, **Current Biology** OD project (submission + resubmission 1) | ocular dominance as CRF shift; dichoptic asymmetry; OD-weighted binocular models | B + J, 25 penetrations, 502 MUA channels (B 1–98, J 99–502); 75 SUA for S1 | `2022_OD Project/current_biology_resubmission_1/` (latest), `ODproject_init.m`; published code `D:\MATLAB\MitchellEtAl_CB_2023` |
| 7 | SfN 2023 | normalization across condition/depth/time | B + J | `2023_SfN/` |
| 8 | 2024–2025, **LMN** | laminar binocular normalization + supersaturation | B + J (11 + 19 penetrations) | `2024_LaminarCRFs/` → [04](04_LMN_analysis.md) |
| — | 2022–2024 tools | CSD/LFP alignment, RF & orientation QC | | `LFP and CSD viewer/`, `RF_Tuning/` |

## 2. Conventions

**Penetration**: `YYMMDD_M_eX` (11 chars). Banks eA–eD; front-panel bank 'E' = analog (excluded).

**Eye codes** — two conventions are in circulation; resolve explicitly in the rebuild:
* README / `diTP` comment: `STIM.eye` 1 = both, 2 = ipsi, 3 = contra (after the LV1 swap in `diTP`).
* `diUnitTuning2021` comments: 2 = RE, 3 = LE ("right eye for B").
* For B (RV1) both are true at once (RE = ipsi). For J (LV1) `diTP` swaps so 2 = LE = ipsi.
  The OD project nonetheless flips **B's** OD sign with hard-coded unit ranges
  (`resub/F2:49` 1:98; `initial/F2:48` 397:end; `version 2/...:111` 414:end) — check.
* Rebuild: store **raw eye (L/R)** per trial plus hemisphere; derive ipsi/contra and DE/NDE.

**Contrast**: E [.225 .45 .9]; I [.05 .15 .30 .50 1] or [.2 .4 .8]; B/J bminteroc [0 .055 .11 .225 .45 .9].
Some code uses % (0–100, blank → 1 %), some fractions (blank → 0.01).

**Conditions**: `UNIT.DE/NDE/BIN/DI`; Cox era `MON.{DE,NDE}_{PS,NS}`, `BIN.{PS,NS}`, dCOS;
`diKEY` = [cond, DE c, NDE c, soa, adaptnull, nTrials]. Standard trial filter:
`~adapted & ~blank & rns==0 & cued==0 & motion==0`.

**Unit types**: auto = autosort dMUA (.ppnev); nev = online sorted; kls = Kilosort; mua = analog MUA (ns6);
lfp, csd (ns2). (`diIDX_diKey_v1.m:48` loads `_AUTO` under the name 'nev'.)

**Depth**: `STIM.depths(:,2)` = contacts relative to SinkBtm (+ = superficial), 100 µm/contact.
SinkBtm/ToC/BoC were **entered by hand** from CSD plots (no algorithm in code).

## 3. Signal recipes found across the repo (they differ)

| signal | variants |
|---|---|
| LFP | ns2 as recorded (/4 → µV), no extra filtering; archive: ns6 → LP 500 Hz → ↓1 kHz |
| CSD | `calcCSD` (MLAnalysisOnline): second spatial difference, 0.1 mm spacing, edges NaN; ×0.4 (→ nA/mm³, σ = 0.4 S/m) in some scripts only; baseline-subtracted in EVPanalysis only; `filterCSD` smoothing default / 0.05 / 0.1; triggers code 23 only vs all onsets |
| analog MUA | HP 750 Hz (4th-order Butterworth, filtfilt) → rectify → then: LP 250 + LP 50 (`diNeuralDat`, 'roelfsema+'); LP 250 + LP 125 (EVPanalysis); LP 200 (archive); window mean only (getRF); PSP-kernel conv (ex_Alignment) |
| SDF | `jnm_kernel('psp')`: (1−e^{−t/τg})e^{−t/τd}, τd = 20 ms, τg = 1 sample (rise time depends on Fs), causal, unit area |
| laminar aids | LFP power profiles (alpha/beta vs gamma "cross"), MUA latency (first z > 3 in 50–100 ms), normalized EVP/MUA images |

## 4. RF mapping and orientation tuning
* **Dot mapping** (`RF_Tuning/analyDotMappingSPK.m`): response on→off code, spikes/s;
  0.05° grid, each dot fills its disc; RF = largest connected component of map > 0.5·max
  (or right-tailed t-test α = .001). `_singleColumn` variant uses tail 'left' (finds suppression — bug).
  `utils/dichoptic/getRF + meanRF + fitRF`: window on+50 ms→off, baseline −50→0, z ≥ 1,
  `fRF = [cx cy wx wy mx my]`; V1 RF limits = channels whose centroid lies within `dthresh` (0.5–1°) of the median.
* **Orientation** (`RF_Tuning/analyRFTuningSPK.m`): spike count per presentation, per eye;
  spline fit after duplicating at +180°; `anovan` p < .05; peak = argmax.
  Cox/bincontrast: `jnm_tune` (gauss1 over circular shifts), `gaussianFit`, `gauss2`.
* Preferred eye/orientation rules differ: ANOVA (diIDX), argmax (diUnitTuning), OD sign (2021).

## 5. External dependencies (not in repo)
* **MLAnalysisOnline** (`D:\MATLAB\MLAnalysisOnline`): `parsEventCodesML, calcCSD, filterCSD,
  getLFP, trigData, readgGrating, readgDotsXY, readgDRFTGrating, readBRFS, concatBHV,
  pEvtPhoto, photoReTrigger, getNeuralData, getMLcodes, f_ShadedLinePlotbyDepth, nanunique, …`
* **MonkeyLogic 2** `mlread` (bhv2), **npy-matlab** `readNPY`, `uCalcTilts0to179`,
  Jake Westerberg's `jnm_psd`, `jnm_reorder`, `jnm_ncorr`.
* **KiloSortUtils** (`D:\MATLAB\KiloSortUtils-master`): `KiloSort2SpikeStruct(_BM)`,
  `setupDataForKiloSort` (NSxToHL), `offlineBRAutoSort_2021`.
* Missing BM functions (called, defined nowhere or under another file name):
  `getPARAMS`, `getPARAMS_1D`, `getPARAMS_2D` (inside `fitSigmoid_2D.m`),
  `getPARAMS_1D_fullTW` (inside `fitSigmoid_1D_fullTW.m`), `fitMODEL2`, `newmodelcrf`, `modelAnalysis`.
* File Exchange: `distributionPlot`, `superbar`, `cprintf`, `boxplot2`, `violinplot`.
* MATLAB toolboxes: Signal, Statistics, Curve Fitting, Image Processing, Mapping, Deep Learning (`combvec`).

## 6. Cross-project inconsistencies to resolve once
1. Eye identity (above).
2. OD index: sign, normalization, balancing, and which field (`occ(1)` is an index in one struct, a p-value in another).
3. Preferred-eye criterion (ANOVA / argmax / sign) and post-hoc `occSwap`.
4. Response windows and baselines (≥6 variants; "tw" means a different window per dataset).
5. Layer boundaries (supra ≥ 3, ≥ 4, ≥ 5; one archived version never labels infragranular).
6. aMUA filter chain; CSD scaling/smoothing/baseline.
7. CRF fitting: two different `BMfitCRdata.m` on the path; % vs fraction; bounds; R² definitions.
8. Additivity index: BIN/(DE+NDE) vs (SUM−BIN)/(SUM+BIN).
9. Archive files shadow active files of the same name if `archive/` is on the path.
