# 02 — The dataset: `D:\NEUROPHYS DATA\2021`

Status: inventory made 2026-10-28, read-only (file listing, NSx/NEV header
parsing, MAT headers, HDF5 field listing). Working tables are in the session
scratchpad (`inventory/listing.tsv`, `nsx_scan.tsv`, `raw_full.txt`).
Obs = observed, Inf = inferred.

21,783 files, ~1.84 TB. Sibling `D:\NEUROPHYS DATA\2016-2018` (monkeys E/I,
iScience-era data) exists but was not inventoried.

---

## 1. Subjects

| Monkey | Name in code/bhv2 | Hemisphere | Sessions on disk | Dates |
|---|---|---|---|---|
| B | "B52" | **RV1** (all penetrations) | 15 | 211007 – 211231 |
| J | "Jo" | **LV1** (all penetrations) | 18 | 221028 – 221215 |

* Older monkeys E ("einstein", 2015–16) and I ("ingrid", 2016) are the iScience
  2022 dataset; helga/bohr (2012–13) appear only in `archive/jnm_sessioninfo.m`.
* A monkey **T** appears in task test files (`211003_T_*`, `211014_T_rfori001`).
* Code and figures reference **B sessions in 2022** (`220202_B_evp004`,
  OD project B range "211008–220302") that are **not on this drive** — open question.

## 2. Hardware / file format (Obs from headers)

* Blackrock NSP, Central "File Dialog v7.0.6"; NSx `NEURALCD` v2.3, NEV `NEURALEV` v2.3.
* `.ns2` = 1 kHz (HP 0.3 Hz 1st order, LP 500 Hz 4th order);
  `.ns6` = 30 kHz raw (HP 0.3 Hz, LP 7.5 kHz). Electrodes in µV, ±8191 µV ↔ ±32764 → **0.25 µV/bit** (the `/4`).
  Analog inputs (`ainp*`) in mV, ±5000 mV ↔ ±32767.
* **First packet timestamp: 102 in .ns6, 0 in .ns2** (see 01 §6). 16 recordings
  have 2 packets with 7–30-sample gaps.
* Per recording: `.ccf` (Central config), `.nev`, `.ns2`, `.ns6`, and for ML tasks
  `.bhv2` + stimulus text file `.g<PARADIGM>Grating_di` (+ `_di_v2` from 211208) +
  `*_GRATINGRECORD####.mat`. Many recordings also have `<rec>.mat` = saved NPMK NEV struct.
  No `.nsv`, `.ns1/3/4/5`, or ML1 `.bhv` files exist.
* NEV for bminteroc in J (and many B) is event-only (~1 MB): no online spikes.
* Probes: Plexon 32-ch, 100 µm spacing (Kilosort `chanMap.mat`: probetype
  'Plexon', 100 µm y-steps). Bank 1 = eA, bank 2 = eB, bank 5 = analog.

### Channel configurations
| n ch | electrodes | analog | sessions |
|---|---|---|---|
| 39 | eA01–32 (early pin order) | ainp1,3,4,5,11,12,13 | 211007–211012 |
| 39 | eA01–32 (later pin order) | same | 211025–211105 |
| 40 | eA01–32 | + ainp15 | 211208, 211213, 211217, 211231, J 221028–221130 |
| 41 | eA01–32 | 9 inputs | 211210 |
| 32 | eA01–24 only | 8 | 211229 |
| **72** | **eA01–32 + eB01–32** | 8 | **J 221202, 204, 206, 209, 212, 214, 215 (dual probe)** |

`ainp1` = photodiode. Electrode labels in the files are always used to select
channels; the pin order differs between 211007–211012 and later (verify that
label → physical contact mapping is the same, via CSD).

## 3. Dual-probe days
Seven J sessions, every recording 72 channels (both probes in **one** set of
files — eA/eB is only a channel-label prefix, never in filenames). In processing:
* eA penetration uses `bminteroc001` / `rfori001`, eB uses `bminteroc002` /
  `rfori002` (221202, 206, 209, 214, 215) — Inf: stimulus placed on each probe's RF in turn.
  Each file still contains the other probe's responses (stimulus off its RF).
* 221212 processed for eA only; 221204 for eB only.

## 4. Paradigms

| paradigm (file token) | stim file | B sess / recs | J sess / recs | what it is |
|---|---|---|---|---|
| **bminteroc** | `.gBMINTEROCGrating_di[_v2]` | 11 / 12 | 18 / 23 | 6×6 L×R contrast grid, static gratings, 3 per trial (main task; see 01 §2.5) |
| bmcBRFS | `.gBMCBRFSGrating_di` | 12 / 14 | 18 / 22 | binocular rivalry flash suppression (Inf: Brock Carlson's "BMC" task) |
| dotmapping | `.gDotsXY_di` | 15 / 22 | 18 / 23 | RF mapping with dots |
| rfori | `.gRFORIGrating_di[_v2]` | 15 / 19 | 18 / 24 | orientation tuning, 16 orientations × eye {both, R, L}, 400 ms on / 200 ms ISI |
| rforiDRFT | | 6 / 7 | – | drifting orientation tuning |
| rfsf, rfsize, rfphase, rfcone/cone | | 9, 9, 2, 6 | cone4lesions 1 | SF, size, phase, cone-isolating tuning |
| natstereo | `.gImageXY_di` | 4 / 4 | 1 / 1 | natural stereo images |
| phzdisparity / posdisparity | | 2 / 1 | – | disparity |
| saliency* (5 variants) | `.lineSeg*_di` | 5 / 6 | – | line-segment saliency |
| **evp** | none | 15 / 60 | 17 / 25 | flash evoked potentials (CSD alignment) |
| **darkrest** | none | 12 / 13 | 8 / 8 | rest in darkness (spectral laminar markers) |
| calibrationPSD / recPSD / psdREC | none | 1 / 1 | 16 / 40 | 60 s PSD recordings |
| mcosinteroc | | 1 / 1 | – | 211012 only |

Total .ns6 recording time: B 44.6 h, J 51.8 h.

### Sessions
B (RV1, all eA): 211007, 008, 009, 012, 025, 027, 102, 103, 105, 208, 210, 213, 217, 229, 231.
**Processed** (STIM 2024 2): 211008, 009, 012, 025, 027, 102, 103, 208, 210, 213, 217 (11).
Not processed: 211007 (no .bhv2 at all), 211105 (no bminteroc), 211229 (24 ch), 211231.

J (LV1): 221028, 102, 103, 104, 109, 118, 119, 121, 123, 128, 130, 202*, 204*, 206*, 209*, 212*, 214*, 215* (*dual).
**Processed** (STIM 2024 2): eA 221109, 118, 119, 121, 123, 128, 130, 202, 206, 209, 212, 214, 215; eB 221202, 204, 206, 209, 214, 215 (19).
221028, 221103, 221104 exist only in older STIM folders; 221102 never processed.

## 5. Processed-data folders

| folder | what | notes |
|---|---|---|
| `STIM 2024 2` | **current** (`setup('BMC_BM')`): base, `_MUA`, `_CSD` for 11 B + 19 J | STIM has 52 fields (adds `sortdirection`) |
| `STIM 2024` | same set, 1 day earlier, different v1lim | |
| `STIM PHYSORTED MUA+GOOD` | base, `_KLS` (23), `_MUA`, `_CSD` | KLS = phy good units |
| `STIM PHYSORTED GOOD`, `STIM_KILOSORTED` | KLS variants | |
| `STIM 2023`, `STIM FULL`, `STIM J`, `STIM J Backup`, `STIM B`, `STIM`, `STIM_orig` | older generations | `STIM` also has `_AUTO`, `_NEV`; `STIM_orig` mixed rfori + bminteroc trials |
| `KiloSort-ed` / `KiloSort-ed_phy` | 24 recordings sorted per file (7 B, 17 J); `_phy` is curated | `.ns6.dat` inputs deleted; KS log at `D:\KiloSort\KiloSort_log.xlsx` |
| `AutoSort-ed` | 40 `.ppnev` (B only) | threshold-crossing dMUA |
| `V1 Limits`, `V1 Limits Auto` | RF + evoked-response limit files (`EV`, `RF`, `fRF`, `v1lim` 1×4) | 1.4–3.2 GB each; not used for 2021 V1 limits |
| `evp files`, `.ns2 files` | copies of raw evp / bminteroc recordings | redundant with raw |
| `Figures\EVPanalysis` | 27 SVGs (B) | |

**v1lim changes between generations** (Obs): `211008_B_eA` [8 24 31] (STIM_orig)
→ [7 22 29] (STIM B … 2023) → [2 18 29] (STIM 2024 2); `221204_J_eB`
[6 16 32] → [20 29 32] (2024) → [10 22 32] (2024 2).

### STIM struct (STIM 2024 2, e.g. `211008_B_eA`, 398 trials)
`header, penetration, V1, paradigm, filelist, filetime, fpath, runtime,
sortdirection, message` · per trial: `task, filen, trl, prez, reward, tp_ec,
tp_sp, tp_pt` · stimulus: `contrast(n×2), tilt(n×2), eyes(n×2), eye, sf, phase,
soa, motion, tf, xpos, ypos, diameter, gabor, timestamp` · flags: `monocular,
ditask, dioptic, tiltmatch, rns, blank, adapted, adapter, suppressor, diconflict,
botheyes, rsvpmask, cued` · laminar: `rmch, rf_xyr, overlap, el_labels,
depths(n_ch×3), v1lim, rank`. KLS-type STIM add `units, kls, phy`.
`_MUA` / `_CSD` add `RESP (ch×4×trials), SDF (ch×time×trials), sdftm, PSTH, psthtm, win_ms`.

## 6. Data-integrity flags
1. Nested duplicate `Raw Data J\221204_J\221204_J\` (47 GB); `bmInteroc001` vs `bminteroc001` case mix.
2. `datafile_` prefix on 9 recordings in 211210_B; unprefixed `211208_B\dotmapping001.*`; `211027_B_darkest001`.
3. 221206_J has no evp recording (TuneList says evp = 1); its `recPSD001` is 34.6 min.
4. Odd channel maps: stray eB02 in `211208_B_evp001.ns6`; 211229 only 24 electrodes; early-B pin order.
5. KS `chanMap.mat` for `221206_J_bminteroc002` points at `bminteroc001.ns6`.
6. Missing `.bhv2` for 211007_B (all), several bmcBRFS/dotmapping files.
7. Older processed files reference `D:\NEUROPHYS DATA\2021\Raw Data\` (renamed to `Raw Data B`).
8. 211012, 211102, 211217 have `SortDirection = 'descending'` / `Flipped = 1` in TuneList, but the 2021 depth code ignores sort direction → verify depth orientation for those.
