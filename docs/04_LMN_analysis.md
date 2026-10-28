# 04 — LMN project (2024_LaminarCRFs): analysis chain, definitions, figures

Status: written 2026-10-28 from a read of `2. postprocess/`,
`3. visualization/2024_LaminarCRFs/` and `3. visualization/2023_SfN/`.
Nothing was re-run. **[inferred]** marks interpretation; everything else is what
the code literally does.

---

## 1. Data flow

```
<pen>_MUA.mat (STIM + SDF, from diRunDir)          C:\Users\<USER>\OneDrive - Vanderbilt\Data Analysis\MATLAB\processed data\LMNProject 2024 2\
   │
   ├─ 2. postprocess/bincontrast_2021.m  (resptype 5, dataform 6, animals J+B)  → UNIT, IDX, info   [current]
   │       └─ saved by hand → MUA_bminteroc_4_seperated/MUA_bminteroc_4_core.mat   [inferred]
   │
   ├─ LaminarCRFS_init.m → getDepth, getPOP, fitSigmoid_1D(_v2), fitLine_1D(_v2), fitNDE, fitModel
   │       └─ saved by hand → _FIT, _FIT2, _NDE, _MODEL, _depth                     [inferred]
   │
   └─ Manuscript_focused_v2/Figure*.m  (scripts on the base workspace; save .svg/.pdf to figDir)

<pen>_CSD.mat → bincontrast_CSD_Spiking.m → PEN{…} → CSD_visualization.m   (CSD_bminteroc.mat, hand-saved) [inferred]
.ns2 → LFP and CSD viewer/EVPanalysis_loop_improved.m → PEN(1:33) → LFP_workspace_initiation.m / LFP_analysis.m (LFP_bminteroc*.mat)
```

No code writes any `*_bminteroc.mat` filename; all workspaces were saved
interactively. `MUA_bminteroc_4_depth.mat` (2 KB) probably holds the 5-column
`depth_array` the v2 figures need — **no code in the repo creates columns 4–5**.

Workspace file formats seen: `CSD_bminteroc.mat` = MAT v7.3 (2024-03-07);
`LFP_bminteroc*.mat` = MAT v5 (2024-03-06/07).

### Current vs superseded
| Status | Files |
|---|---|
| Current | `bincontrast_2021.m`; `fitSigmoid_1D(_v2)`, `fitLine_1D(_v2)`, `fitMeanCRF`, `fitMeanSigmoid`, `fitMeanLine`, `fitMODEL.m` (defines `fitModel`), `getModels`, `fitNDE`, `refitBIN`, `getDepth`, `getPOP`; `Manuscript_focused_v2/` + `Simulation/` |
| Predecessor | `Manuscript_focused/` (v1, Fig 1–12; Fig 6 still edited Jul 2025), `LMNProject_*.m` |
| Superseded / broken | `bincontrast.m` (2022, monkey I), `bincontrast_CSD.m` (broken WIP), `fitBIN`, `fitGain`, `fitMODELxLAY`, `*_fullTW`, `fitSigmoid_2D`, `genMODEL`, `getCRF`, `getLAY*`, `getOCC*`, `getDepth2`, `2. postprocess/simulateV1Response.m` (differs from the `Simulation/` copy) |

---

## 2. Data structures (bincontrast_2021, resptype 5)

**Windows** (`info.respWin`, 50 rows; SDF sample k ↔ (k−151) ms, −150…200 ms):
w1–11 = 50 ms wide, 10 ms step (w1 = 50–100 … w11 = 150–200);
w12–48 = 20 ms wide, 5 ms step (w12 = 0–20, w21 = 45–65, w30 = 90–110, w33 = 105–125, w48 = 180–200);
**w49 = 50–200 ms "full"**; **w50 = −150…−50 ms baseline**.
(resptype 4 — used by the CSD/SfN code — has only 13 rows; w12 = full there.)

**`UNIT(u)`** — u is one V1 *channel* (MUA):
* `.DE / .NDE / .BIN`: `RESP_avg|_sd|_error` (6 contrasts × 50 windows, row 1 = blank),
  `RESP_trls{c}` (50 × trials), `SDF_avg…` (6 × 351), `SDF_trls{c}`.
* `.DI`: same with 36 rows, `diLevels = combvec(cLevels,cLevels)'`; col 1 = DE, col 2 = NDE,
  **DE varies fastest** (row r = i + 6(j−1)). r1 blank; r2–6 DE-mono; r7,13,19,25,31 NDE-mono;
  r8,15,22,29,36 dioptic; other 20 dichoptic.

**`IDX(u)`**: `penetration, N, depth (= STIM.depths(e,:)'), rate, tonic, layer, v1lim, X, DE, NDE,
occ, occana, cLevels, diLevels, trial counts…`; `X` = `diUnitTuning2021` output (`occ, occ_norm,
contrast_tuning, dicontrasts, prefeye, nulleye…`).

**`info`**: dataset, animal, datatype, N, sdfWin, respWin, cLevels, diLevels, resptype, dataform, normtype…

**Fits**: `SigmoidFit / LineFit .{DE,NDE,BIN}.{a,k,n,s | slope}` (window × unit), `.rsquare, .aic, .gof`;
`NDE.w.{value,rsquare,fitobject}`; `MODEL.{LSM,NRM}.unit.{a,k,n,s,rsquare(adj),AIC,AICc}` and `.pop`.

**Depth**: `depth_array` rows superficial→deep: col1 = depth+2, col2 = n units, col3 = layer 1–3
(v2 figures also use col4 = layer label 1–6, col5 = mm labels — only from the saved workspace);
`depth_idx{row}`, `layer_idx{1..3}`.

---

## 3. Definitions

**Response.** Each window mean minus w50 baseline; SDF minus mean of samples 1–101. Per-unit
normalization (dataform 6): divide by the trial-mean of the 95th percentile across windows (RESP)
or samples (SDF) — *separate constants for RESP and SDF*; negatives and NaN → 0. Trial rejected if
any value > 100 (raw), all < 1e-4, or Inf. Windows used: full = w49; coarse transient = w1;
time-resolved = w21:48; "transient" = w21:30; "sustained" = w31:48.

**Inclusion.** Channel inside ToC–BoC (preprocessing); `diUnitTuning2021` succeeds; contrast-tuned
(one-way ANOVA of monocular full-window responses across contrast, both eyes pooled, p ≤ 0.05);
≥1 trial in every non-blank DE/NDE/BIN condition; CRF fit succeeds (via missing `getPARAMS_1D`);
depth analyses keep depth(2) ∈ [−12, 13]. **No** orientation, visual-responsiveness, OD-significance
or rate criterion.

**Depth / layers** (units = channels from SinkBtm; positive = superficial; 100 µm/channel assumed):

| Definition | Supra | Gran | Infra | Used by |
|---|---|---|---|---|
| `IDX.layer` (bincontrast_2021) | d ≥ 5 | 0…4 | d < 0 | stored per unit |
| `getDepth` (+2 shift) | d ≥ 3 | −2…2 | ≤ −3 | v2 manuscript |
| `getLAY2` | d ≥ 4 | 0…3 | < 0 | SfN |
| `depth_array(:,4)` (workspace) | labels 1–3 | 4 | 5–6 | v2 only |

**Ocular dominance** (`diUnitTuning2021`): full-window, baseline-subtracted, *unnormalized* monocular
responses; `occ = (R2 − R3)/(R2 + R3)` (2 = ipsi, 3 = contra); `occ_norm` = mean across contrasts of
per-contrast index; **DE = sign(occ_norm)** — no significance test, and computed on the same data
later analysed (circular).

**CRFs** (Curve Fitting Toolbox `fit`, NonlinearLeastSquares, unweighted, 6 trial-mean points, blank at x = 0.01):
* Naka–Rushton `R = a·xⁿ/(xⁿ + kⁿ) + B`, B fixed at DE blank; supersaturating `R = a·xⁿ/(x^{sn} + k^{sn}) + B`.
  Bounds a ∈ [R(.9), 1.2R(.9)], k ∈ [0.01, 0.8], n ∈ [0.3, 4], s ∈ [1, 2]. AIC = 2k + N ln(SSE/N).
* Linear `y = m·x + B`. `_v2` versions intended to drop the blank (sigmoid_v2 doesn't — bug).
* `fitNDE 'w'`: NDE `R = A·w·x^N/(x^N + K^N) + B` with DE A,K,N,B fixed → NDE gain weight w ∈ [0.01, 1.5].
* Population: `fitMeanCRF`, `fitMeanSigmoid` (lsqcurvefit + bootstrap over units), `fitMeanLine` (robustfit on log10 c).

**Supersaturation**: SI = (Rmax − R(0.9))/(Rmax − R0) on BIN means (v1, SI > 1 → NaN; violin at tw 33 =
105–125 ms); v2 uses MI = 1 − SI. Parameter `s` from the supersaturating fit; examples use s > 1.2 at tw 30.

**Binocular indices** (dioptic BIN vs monocular at same c): BM = (BIN − DE)/(BIN + DE);
OD = (DE − NDE)/(DE + NDE); SUM = BIN/(DE + NDE); QD = BIN/√(DE² + NDE²).
Dichoptic summation: R(i,j)/(R_DE(i,0) + R_NDE(0,j)) for the 25 pairs with both > 0.

**Binocular models** (`fitModel`, x = DE, y = NDE contrast; 0 → 0.01; fit to all 36 DI means,
trial-count weighted):
* LSM (linear summation of monocular NR): `R = a·xⁿ/(xⁿ+kⁿ) + a·w·yⁿ/(yⁿ+kⁿ)`
* NRM (binocular normalization): `R = a(xⁿ + w·yⁿ)/(xⁿ + yⁿ + kⁿ)`; supersaturating version with `s·n` in the denominator.
* w fixed per unit/window from `NDE.w` (or 1 − ocularity). No baseline term.
* "Normalization strength" v2 = ΔR²(NRM − LSM); v1 = (NRM − LSM)/(NRM + LSM).

**Dynamic normalization** (`Figure_12`): conditions ordered by summed contrast; time-to-peak, peak
magnitude change (exp1 fits) and SD-after-50 ms (poly1) vs summed contrast, per unit, by depth.

**Simulation** (`Simulation/simulateV1Response.m`): Heeger ORGaNICs
(`toolboxes/ORGaNIC_Normalization`), scalar contrast step input, 2 orientation channels; cost
functions RMSE / L1 / correlation / xcorr / DTW; optimizers bayesopt (100 evals), fminsearch, random
search, manual. No fitted parameters saved — exploratory/unfinished. The model has no binocular structure.

---

## 4. Stimulus conditions as analysed
Contrasts {0, 0.055, 0.11, 0.225, 0.45, 0.90} (`cLevels = [0 X.dicontrasts]`, per session). DE =
monocular eye == DE; NDE = monocular other eye; BIN = dioptic; DI = all 36 pairs; blank = `STIM.blank`.
Orientation matched across eyes, never used for selection. Monkeys B and J.

---

## 5. Manuscript_focused_v2 figure list
| Script | Shows | Known problems |
|---|---|---|
| Figure1D_MUA_responses | 6×6 grid of depth×time DI SDF maps | — |
| Figure2A_contrast_response_functions | population + per-layer DE CRFs (w1) | bootstrap N = 5; 6×SEM error bars; crashes at AIC printout (`uSigmoid` field) |
| Figure2B_scatter_saturation | per-unit −AIC sigmoid vs line (w49) | AIC across different N (6 vs 5 points) |
| Figure2C_laminarHistograms | per-layer KDE of R² and −AIC | mixes v1 and v2 fits; overwrites `colors` |
| Figure2D_saturation_across_depth_v2 | R² / ΔR² vs depth | needs `depth_array(:,5)`; layer lines inconsistent (0/4 vs 1/5) |
| Figure3B_integration_scatter | NRM vs LSM adj R² per unit | — |
| Figure3C_integration_laminarHistograms | per-layer KDE NRM vs LSM | — |
| Figure3D_integration_across_depth_v2 | population surface, ΔR² vs depth | legend swapped; needs col 5 |
| Figure4A_integration_spatiotemporal | depth × time ΔR² (NRM − LSM) | — |
| Figure4B_saturation_spatiotemporal | depth × time AIC-weight sigmoid − linear | — |
| Figure5_Laminar_SDFs | per-layer DI SDFs, peak magnitude/latency | indexes time where it means condition; undefined vars |
| Figure5_observed_data_* / simulated / simulation_vs_observed | ORGaNICs vs observed SDFs | wrong blank subtracted; undefined `index_vector_2` |
| Figure6_supersaturation | SI distributions, SI/s/R² over depth × time | needs `PARAMS`, `PARAMS_s`, `BIN2` (missing `getPARAMS_1D`); MI row offset bug |
| unit_stats, Supplementary/* | counts, CRF params vs depth, movie | — |

---

## 6. Fragilities for the rebuild (selected; consolidated in 05)
1. `getPARAMS_1D` is missing (called by bincontrast_2021, Figure6, NRMproject_init); PARAMS ≈ SigmoidFit [inferred].
2. File/function name mismatches (`fitMODEL.m` defines `fitModel`, etc.); `LaminarCRFS_init.m:73` calls `fitModel` with wrong args.
3. Three-plus layer definitions in use; y-axes say mm but units are channels.
4. Window numbering depends on resptype (13 vs 50 windows); several scripts assume tw 12 = full.
5. AIC compared across different numbers of points; v1/v2 fits mixed.
6. `IDX.tonic` ≈ 0 (computed after baseline subtraction).
7. RESP and SDF normalized by different constants; CSD_Spiking rescales DI to [0,1] separately.
8. Indices averaged over contrast include the blank row (OD = 0, SUM = 0.5 there).
9. Blank trials not filtered by the same task mask as other conditions.
10. Hard-coded 6-contrast ordering, unit IDs, depth crops, row ranges, paths (`OneDrive - Vanderbilt` vs `OneDrive\4. Vanderbilt`).
11. `fitNDE` runs windows where PARAMS is zero (w1–20); `refitBIN` uses a linear index for B; `fitModel` population fit reuses last unit's start point.
12. Ties in summed-contrast ordering (e.g. 0.45+0.45 vs 0.9+0) resolved arbitrarily.
13. ORGaNICs simulation: Euler step 1 ms vs τ 0.1–0.5 ms (unstable); `Ntheta` name mismatch; bayesopt searches unused params.
