# 01 — Preprocessing pipeline: raw Blackrock/MonkeyLogic → STIM → neural data

Status: written 2026-10-28 from a full read of `1. preprocess/` plus the external
code it calls (`D:\MATLAB\MLAnalysisOnline`, `D:\MATLAB\KiloSortUtils-master`,
`D:\MATLAB\NIMH-MonkeyLogic\33 binocular gratings`). Nothing here has yet been
verified by re-running code on data; items marked **VERIFY** are hypotheses
derived from reading code and must be checked empirically (see
[05_known_issues_and_QC.md](05_known_issues_and_QC.md)).

---

## 0. Big picture

```
DataLog.xlsx (TuneList sheet) ──► importTuneList ──► runTuneList ('BMC_BM_2021')
                                                        │
             raw files per recording (RIGDIR/<YYMMDD_M>/<YYMMDD_M_task###>.*)
             .nev  .ns2  .ns6  .bhv2  .g<TASK>Grating_di
                                                        │
                       diTP        → trial table, event-code sample times (tp_sp)
                       diCheck     → sanity checks on task content
                       diPT_2021   → photodiode re-trigger (tp_pt)
                       V1 limits   → from DataLog ToC / SinkBtm / BoC columns
                                                        │
                                   STIMDIR/<YYMMDD_M_eX>.mat  (STIM struct, -v7.3)
                                                        │
                diRunDir (GUI: choose datatype) → diNeuralDat / diNeuralDat_usePhy
                                                        │
               STIMDIR/<YYMMDD_M_eX>_<DATATYPE>.mat   (STIM, RESP, SDF, sdftm, PSTH, psthtm[, SUA])
                                                        │
                                       2. postprocess / 3. visualization
```

One STIM file per **penetration** = per probe. Penetration name is
`YYMMDD_<Monkey>_e<Bank>`, e.g. `211008_B_eA`. Dual-probe days produce two
penetrations (`_eA`, `_eB`) from the *same* raw files; the bank letter selects
the electrode labels (`eA01…`, `eB01…`).

Paths come from globals set by `setup.m` (`setup('BMC_BM')` is the 2021 config):
`RIGDIR = D:\NEUROPHYS DATA\2021\Raw Data B\`, `AUTODIR = …\AutoSort-ed\`,
`SORTDIR = …\KiloSort-ed_phy\`, `ALIGNDIR = …\V1 Limits\`,
`STIMDIR = SAVEDIR = …\STIM 2024 2\`. `setup.m` also expects the repo and
`MLAnalysisOnline` under `C:\Users\<user>\Documents\MATLAB\` (not where they
are on this machine).

---

## 1. Session metadata: `DataLog.xlsx` / TuneList — [importTuneList.m](../1.%20preprocess/importTuneList.m)

* Reads sheet `TuneList`, first 29 columns. Columns used: `Date, Monkey, Bank,
  Structure (RV1/LV1), Rig, SortDirection, SinkBtm, Penetration, Rank,
  BadCh/BadBtmCh`, plus one column per paradigm whose cell is a MATLAB
  expression of file numbers, e.g. `bminteroc = [1]`, `evp = [1,2,3,4]`
  (`cellfun(@eval,…)`, [importTuneList.m:57-59](../1.%20preprocess/importTuneList.m#L57-L59)).
* Rows with `Rank == 0` are dropped ([L66](../1.%20preprocess/importTuneList.m#L66)).
* `importTuneList(2)` keeps only penetrations that have ≥1 `bminteroc` file.
* The 2021 V1-limit branch of `runTuneList` re-reads the workbook and needs
  columns `ToC`, `BoC`, `SinkBtm`, and the penetration name in **column 13**
  ([runTuneList.m:243-254](../1.%20preprocess/runTuneList.m#L243-L254)).

**Gap:** the workbook is git-ignored and is *not* in the repo. The only copy
found on this machine is `D:\MATLAB\ephys-analysis\stim\DataLog.xlsx`
(Jul 2022). Its TuneList has 2021 rows only for monkey **B, bank A**
(`211007 … 211217`), and **has no `ToC`/`BoC` columns**. The workbook the
current code expects (with ToC/BoC, eB rows, and monkey J 2022 sessions) has
not been located. Its contents are recoverable: every STIM folder's `ERR.mat`
saves the full `TuneList` struct (incl. `ToC, SinkBtm, BoC, BadCh, Rank,
SortDirection, Flipped`, per-paradigm file numbers), and each STIM file stores
`v1lim = [ToC SinkBtm BoC]`. Note the `ERR.mat` in `STIM 2024 2` holds only the
B penetrations (overwritten by the B run); J rows must come from other folders.
All versions are collected in [legacy_TuneList_all_versions.csv](legacy_TuneList_all_versions.csv)
(9 `ERR.mat` files + `DataLog_2021.xlsx` recovered from git commit `574f625^`).
`runTuneList` never remaps labels for `SortDirection = 'descending'`; see 05 C3.
**The V1 limits differ between STIM generations** (e.g. `211008_B_eA`:
[8 24 31] → [7 22 29] → [2 18 29] in `STIM 2024 2`), see 05.

---

## 2. Trial table and event codes — [diTP.m](../1.%20preprocess/diTP.m)

### 2.1 File ordering and paradigm
* Files sorted by NS2 header `DateTime` ([diTP.m:39-50](../1.%20preprocess/diTP.m#L39-L50)).
* Paradigm = `BRdatafile(10:end-3)` (e.g. `bminteroc`); must have a matching
  stimulus text file, e.g. `.gBMINTEROCGrating_di` ([L5-17](../1.%20preprocess/diTP.m#L5-L17)).

### 2.2 Event codes
* `NEV.Data.SerialDigitalIO.UnparsedData - 128` (strobe bit) with
  `TimeStamp` in 30 kHz samples ([L86-88](../1.%20preprocess/diTP.m#L86-L88)).
* `parsEventCodesML` (MLAnalysisOnline) splits into trials: a trial starts at
  the first `9` of a `9 9 9` triplet (detected as `diff(diff(find(code==9)))==0`)
  and ends at `18` + the next two codes. If the last code isn't 18 the final
  trial is dropped.
* **bminteroc event codes** (from `T_binocular_gratings.m`, the ML2 timing file,
  `D:\MATLAB\NIMH-MonkeyLogic\33 binocular gratings\`):

  | code | meaning |
  |---|---|
  | 9 / 18 | trial start / end (ML) |
  | 116+block, 116+condition, 116+(trial mod 10) | trial ID triplet, sent first |
  | 35, 11 | fix spot on, start wait-for-fixation |
  | 8 | fixation acquired |
  | 23 / 25 / 27 | stimulus 1 / 2 / 3 ON |
  | 24 / 26 | ISI scene start (= stim 1 / 2 OFF) |
  | 32, 36 | after stim 3: end scene (32) + fix cross off (36) — **there is no 28** |
  | 97, 36 | fixation broken |
  | 12 | end wait fixation (no fixation acquired) |
  | 96 | reward |

* Per trial, the number of *completed* presentations = count of
  {24, 26, 28, 30, 32} ([L709](../1.%20preprocess/diTP.m#L709)). Onset of
  presentation p = code `21+2p`, offset = `22+2p`; for p = 3 the offset is
  `28 | 32` — a hard-coded fix for exactly the "no 28" behaviour above
  ([L728-730](../1.%20preprocess/diTP.m#L728-L730)).
* `STIM.tp_ec = [21 22]+2p`, `STIM.tp_sp = [onset offset]` in 30 kHz NEV samples.

### 2.3 Trial ↔ stimulus-text-file matching (important)
* Stimulus parameters come from the text file (`readgGrating`, MLAnalysisOnline).
  The ML2 timing file appends **one row per planned presentation (3 per trial) for
  every trial**, including aborted trials, with `TrialRecord.CurrentTrialNumber`.
* `diTP` matches by index: NEV-parsed trial `t` ↔ text rows with
  `grating.trial == t`, and completed presentation p ↔ row p of that trial
  ([L708, L741](../1.%20preprocess/diTP.m#L708)). This is only correct if the
  Blackrock recording started before ML trial 1 and no trial was lost/split.
* `checkTrMatch` tests (1) sequential text-file trials, (2) sequential
  event-code trial digits (116-series), (3) ML timestamp vs NEV start time.
  **In Oct 2021 the failure was demoted to a warning**
  ([diTP.m:144-150](../1.%20preprocess/diTP.m#L144-L150)), with the comment
  "2021 problems arise here when we don't have the event code sequencing
  triplets OR ML and BR timestamp misalignment". → **VERIFY** for every file.

### 2.4 bminteroc stimulus coding ([L357-392](../1.%20preprocess/diTP.m#L357-L392))
* `contrast = [fixedc contrast]`. In the timing file `fixedc` = **left eye**,
  `contrast` = **right eye**. `eyes = [3 2]` for every row, so raw eye code
  **3 = left eye, 2 = right eye**.
* `eye`: 2 if only column 2 > 0, 3 if only column 1 > 0, 1 if both > 0; NaN for blank.
* `tilt = [tilt tilt]` (same orientation both eyes); `tilt(contrast==0) = NaN`.
* Flags: `blank`, `monocular`, `dioptic` (equal contrast & tilt),
  `diconflict` (both eyes, unequal contrast), `botheyes`, `tiltmatch`.
* **Hemisphere swap** ([L544-554](../1.%20preprocess/diTP.m#L544-L554)): for
  `LV1`, codes 2↔3 are swapped in `eye` and `eyes`. Result: **2 = ipsilateral,
  3 = contralateral** in all STIM files. (RV1: right eye is ipsi → raw 2 already
  ipsi; consistent with the timing file.) Contrast columns stay attached to
  `STIM.eyes(:,1)` / `STIM.eyes(:,2)`.

### 2.5 Stimulus parameters actually used for bminteroc
From `genGratingRecordML2.m` (`case 'bminteroc'`,
`D:\MATLAB\NIMH-MonkeyLogic\50 tuning and control\4. RF tuning\`):
* contrast per eye ∈ {0, 0.055, 0.11, 0.225, 0.45, 0.90} — full 6 × 6
  factorial (36 combos incl. 0/0 blank); `mintr = 30` per combo (1080
  presentations), **drawn with replacement** (`randi`), so per-condition counts
  vary around 30.
* 3 presentations / trial, **300 ms on, 300 ms ISI**, static (tf 0), phase 0,
  same orientation/SF/size/position in both eyes (set per session at the rig;
  defaults: rf [-2,-1] dva, 2° diameter, 0.9 cpd, 112°).
* Mirror stereoscope / split screen: left-eye image at `-0.25·screen width`,
  right-eye at `+0.25·screen width`; fixation tracked on the left-eye cross.
* `genBMinterocParams.m` (J. Rogatinsky, 2019; 0.05/0.10/0.22…) is an
  unused older draft — do not use its contrast values.

---

## 3. Photodiode re-trigger — [diPT_2021.m](../1.%20preprocess/diPT_2021.m), [photoReTriggerSTIM_2021.m](../1.%20preprocess/photoReTriggerSTIM_2021.m)

* Photodiode = analog input `ainp1` in the **.ns6** (30 kHz). The ML2 task draws
  a 3°×3° box in the **lower-right corner**: white during each stimulus scene,
  black during ISIs, gray (background) after stim 3.
* Requires `MLConfig.Resolution == '1920 x 1080 120 Hz'` (ViewPixx);
  `refresh = ceil(30000/120) = 250` samples (8.33 ms).
* For every onset and offset event time `tp`: read ±3 refreshes (±750 samples)
  of `ainp1`, take `|z-score|`, threshold **1** (was 3 in the original).
  If the post-event range exceeds the pre-event range by > 1, `adj` = first
  post-event sample above threshold; otherwise `adj` = last pre-event crossing
  + 1 refresh.
* **`tp_pt = tp + adj − refresh`** ([L157](../1.%20preprocess/photoReTriggerSTIM_2021.m#L157)).
  The original (CRT-era) code subtracted `refresh·yconv`,
  `yconv = 0.5 + ypos·PPD/ScreenYres` (≈ ½ frame for a central stimulus, to
  account for raster scan between stimulus row and photodiode row). The 2021
  version subtracts a **full** frame. With the photodiode at the *bottom* of a
  top-to-bottom-scanned panel, the stimulus row is reached before the photodiode
  row, so some negative correction is right, but a full frame likely
  over-corrects by ~3–4 ms for a near-central stimulus → **VERIFY**
  (latencies are directly affected).
* The NS6 read takes the *trailing* `diff(win)+1` samples
  ("typically a 103 sample offset", [L115-116](../1.%20preprocess/photoReTriggerSTIM_2021.m#L115-L116)) —
  see §6.
* Failures → `NaN` rows; whole-file failure aborts the penetration.

---

## 4. V1 limits and depth

**2021+ data** ([runTuneList.m:241-270](../1.%20preprocess/runTuneList.m#L241-L270)):
from DataLog `ToC` (top of cortex channel), `SinkBtm` (bottom of the initial
CSD sink = L4/L5-ish reference), `BoC` (bottom of cortex), read manually from
CSD profiles.
* `STIM.el_labels = eX<ToC…BoC>` (assumes channel number increases with depth;
  `ToC:1:BoC` is empty otherwise).
* `STIM.depths` (n × 3): `[contacts from top, contacts relative to SinkBtm
  (positive = above), contacts from bottom]` — units are **contacts**, not µm.
* `STIM.v1lim = [ToC SinkBtm BoC]`.
* `STIM.rmch` (BadCh) is stored but **not applied** for 2021 data.
* RF overlap (`rf_xyr`, `overlap`) is empty for 2021 data.

**Pre-2019 data** used `diV1Lim` + `V1 Limits` files (RF-mapping-derived limits
intersected with evoked-response limits). Not used for 2021.

---

## 5. Neural data extraction — [diRunDir.m](../1.%20preprocess/diRunDir.m), [diNeuralDat.m](../1.%20preprocess/diNeuralDat.m)

`diRunDir` (listdlg GUI) loops STIM files whose name is exactly 15 chars, and
writes `<pen>_<DATATYPE>.mat` with `STIM, RESP, win_ms, SDF, sdftm, PSTH,
psthtm` (+ `SUA` for KLS).

| datatype | source | processing |
|---|---|---|
| `lfp` | .ns2 (1 kHz) | `/4` (→ µV) |
| `csd` | .ns2 | `/4`, `calcCSD` ×0.4; RESP = mean |sink| (positive values zeroed) |
| `mua` | .ns6 (30 kHz) | `f_calcMUA2D 'roelfsema+'`: 4th-order Butterworth HP 750 Hz → rectify → LP 250 Hz → LP 50 Hz (all `filtfilt`), then downsample ×30 |
| `auto` | `AutoSort-ed/<file>.ppnev` | threshold-crossing "dMUA" from `offlineBRAutoSort_2021` |
| `nev` | .nev | online-sorted spikes |
| `kls` | `KiloSort-ed_phy/<file>/` via `importPhy` + `diNeuralDat_usePhy` | phy `good` clusters |

* `TP = STIM.tp_pt` (30 kHz). LFP/CSD use `round(tp_pt/30)` on the 1 kHz NS2.
* `win_ms = [50 100; 150 250; 50 250; -50 0]` → RESP(:, 1..4, trial) =
  transient, sustained, full, baseline. Windows are clipped at stimulus offset.
* `sdftm` = −300 ms … 150 ms + max stim duration, then trimmed to
  −150 ms … end−150 ms. With 300 ms ISI, the pre-stimulus part of
  presentations 2–3 contains the previous stimulus' offset response.
* Spikes → SDF by convolution with `jnm_kernel('psp', 20 ms)` (causal PSP
  kernel, padded with zeros so `conv(…,'same')` stays causal) ×1000 → spikes/s.
  PSTH: 10 ms bins.
* `calcCSD`: standard second spatial difference, **hard-coded 0.1 mm contact
  spacing**, no smoothing, no boundary (Vaknin) padding → first/last channel
  NaN. **VERIFY** spacing against the probes used.

### Kilosort / phy path ([importPhy.m](../1.%20preprocess/utilities/importPhy.m))
* Kilosort was run **per file** (not concatenated). Binary made with `NSxToHL`
  (raw samples, no header, no zero-padding) by `KiloSortUtils/setupDataForKiloSort`.
* `importPhy` uses **only one file per penetration** — the first file whose
  name contains the selected task.
* Only clusters labelled `good` in phy `cluster_info.group` are kept, mapped to
  channels via `chanMap.mat` (`chanMap0ind`, `chanMapLabels`), and restricted
  to `STIM.el_labels`.
* **Bug:** units are built by `for d = 1:length(unique(depth))` but index
  clusters by `d` ([L116-127](../1.%20preprocess/utilities/importPhy.m#L116-L127)),
  so whenever ≥2 good units share a channel, the *last* clusters in the list are
  silently dropped (count dropped = n_good − n_unique_channels).
* **Bug:** in `diNeuralDat_usePhy`, trials from files that were not sorted get
  `RESP = 0`, `SDF = 0` (not NaN), because `empty(e)` is never set when no unit
  matches the file.

---

## 6. Clock alignment between NEV and NSx — **VERIFY first**

Blackrock NSx data packets carry a start timestamp `T0` (NSP clock ticks). The
bundled NPMK (`openNSx` 6.2.0.0) handles windows like this:

* `openNSx(file,'t:a:b','sample')` reads **file-relative** samples `a…b`
  ([openNSx.m:814](../1.%20preprocess/NPMK/openNSx.m#L814)) and then, for a
  non-paused file, **prepends `T0` zeros** ([L869-873](../1.%20preprocess/NPMK/openNSx.m#L869-L873)).
* Event times (NEV) are on the NSP clock. File-relative sample `k` of an NSx
  file is NSP time `T0 + k`.
* The pipeline "fixes" the extra length by keeping the *trailing* samples
  (photodiode: `dx`; MUA: `trim`). That removes the zeros but leaves the data
  shifted by `T0` relative to the NEV clock.

**Measured (2026-10-28, headers of all 430 .ns2 / 429 .ns6):** first data-packet
timestamp is **102 in every .ns6** and **0 in every .ns2**. (Matches the code
comment "typically a 103 sample offset" = 102 zeros + 1-based index.) 16
recordings have two packets separated by 7–30 samples (e.g.
`211027_B_bminteroc001`, `221103_J_bminteroc001`, `221128_J_bminteroc001`);
NPMK treats these as "paused" files and does *not* zero-pad them.

Resulting offsets, derived from code + measured T0 (still to be confirmed by a
direct photodiode/event check, see 05 §A):

| signal | how read | offset vs. true photodiode onset |
|---|---|---|
| photodiode → `tp_pt` | NS6 window, trailing samples (file-relative) | `tp_pt` = true crossing − 102 samples (−3.4 ms), then − 1 refresh |
| MUA (NS6) | NS6 window, trailing samples | 102-sample error **cancels** → consistent with `tp_pt` |
| KLS spikes | `NSxToHL` binary = file-relative samples | **cancels** → consistent with `tp_pt` |
| LFP / CSD (NS2, T0 = 0, no padding) | NS2 window on NSP clock vs `round(tp_pt/30)` | responses appear **+3.4 ms late** relative to MUA/KLS |
| auto (.ppnev) | full `openNSx` read (zero-padded → NSP clock) | spikes **+3.4 ms late** |
| nev (online) | NSP clock | spikes **+3.4 ms late** |
| paused (2-packet) files | segment 2 shifted by the 7–30 sample gap | trials after the pause off by a further 0.2–1 ms |

So MUA and Kilosort units share one time base, and LFP/CSD/NEV/autosort share
another that is 3.4 ms later. That matters for any CSD-sink-vs-MUA latency
comparison and for laminar latency analyses; it is irrelevant for magnitudes.
Absolute latency also inherits the refresh-subtraction question in §3. The
rebuild must read every stream by **NSP timestamp**.

### 6.1 Event-code checks (measured 2026-10-28)
On `211008_B_bminteroc003`, `211217_B_bminteroc001`, `221109_J_bminteroc001`,
`221206_J_bminteroc002`:
* `9 9 9` / `18 18 18` triplets present on every trial → `parsEventCodesML` is valid.
* Trial-ID codes 116–125 are **absent on 211008** (present on 211217 and all J) →
  `checkTrMatch` test 2 cannot run for early-B sessions.
* Stimulus text file trial count = NEV trial count, first trial = 1, 3 rows per
  trial, on all four files → index matching is correct here (check all files in the rebuild).
* Event-code on→off interval: **216.7 ms (26 frames) on 211008**, **316.7 ms
  (38 frames) on the later sessions**; ISI 316.7 ms everywhere. 211008 therefore
  ran a shorter stimulus than the 300 ms in the current `genGratingRecordML2`.
  Event codes are frame-quantized; true on-durations must come from the photodiode.
* 211008 stimulus file `header` column = `cinteroc`, SF = 4 cpd, diameter 1°;
  column layout of the text file changed between 211008 and 211217
  (`gaborfilter_on, gabor_std` → `grating_phzdist, grating_posdist`).
