# Reviewer's read of draft V2 (2026-10-31)

Two reads of `draft_v2.docx` (text with tracked changes accepted, Figs 1–6 as embedded), done before revising into
V3:

* **Informed read** (Claude, with the analysis record): checks the text against everything learned since V2.
* **Cold read** (a separate reader given only the text and figures, as a journal referee): full report in
  [review_v2_cold_read.md](review_v2_cold_read.md). Paragraph numbers [n] below refer to it.

Where the two agree, the point is listed once. Blake's two new V2 comments and the parked Figure 5 plan are at the
end.

## Verdict

The core observation is strong and visible without any derived measure (Fig 2a): at 45% the response to both eyes
first tracks the dominant eye at 90%, then falls to the dominant eye at 45%. What V2 does not yet do is show that the
measure built on it is valid, that the delay it infers is specifically interocular, and that the laminar and model
claims follow. The writing problems mostly follow from those gaps: claims stated firmly in the abstract are hedged
in the Results.

In order of importance:

1. **Validity of the measure has not been shown** (cold read 1–4; informed concern 3). No fixed-rule null has been
   run through the pipeline; the dominant eye is selected on the same data; the 20–50 ms window sits on the
   transient's falling flank, where waveform-shape differences can mimic a collapse; MUA pools eye-specific neurons.
2. **"Add first, normalize later" is stronger than the evidence** (cold read 1). An index of 1 means the other eye
   acts as extra contrast, not that normalization is absent (any model with w = a gives 1 at all times). A *static*
   model with w > a also gives an index near 1 while the response rises and a deeper decline at high contrast. The
   static-vs-delayed comparison uses total error, which the monocular transient dominates, so it does not test
   whether the *interocular* part is delayed. The Discussion's "interocular slower than monocular" contradicts the
   one-delay model and Fig 2b/6a, which show the monocular response saturating more over time.
3. **Novelty is not stated against Cox et al. 2019 and Mitchell et al. 2022/2023** (both reads).
4. **Laminar and model claims over-reach** (both reads): "normalization appears first in the earliest sink" is
   uncorrected, at 15 ms resolution, and adds nothing beyond "the sink responds first" if each depth normalizes a
   fixed time after its own onset; Fig 4b's visible depth profile contradicts "it did not differ"; Fig 5e shows the
   single-stage model unbiased but imprecise, not failing. The slow-rise claim (Fig 6) is withdrawn (noise).
5. **The 6 × 6 design is under-used** (cold read): unequal-contrast cells are never analysed. Early-vs-late
   iso-response contours in the (dominant, other) contrast plane would connect directly to the psychophysics
   (Legge & Rubin 1981; Meese et al. 2006; Ding & Sperling 2006) and could become the centrepiece.

## Triage of the cold read

| Cold-read point | Assessment | Action |
|---|---|---|
| Fixed-rule null through the identical pipeline | **Agree, top priority.** The fitted model with w set to a is a ready-made null: binocular = monocular at (1 + a)c at every moment, with realistic transient shapes. | Run `windowed_L` on the w = a predictions of each penetration's fit, with the data's noise added; the index must stay ~1 in every window, especially 20–50 ms. |
| Dominant eye selected on the same data | **Agree** (also informed concern). Fig 1e's eyes diverge most at 90%, the selection contrast. | Select dominance and a on odd presentations, analyse even (and reverse); stratify by ocular dominance. |
| Static model with w > a may produce the sequence | **Agree; the key model question.** | Overlay static and delayed predictions of the index time course on Fig 2c; compare static vs delayed on the interaction loss alone; fit separate monocular and interocular pool delays. |
| MUA pools eye-specific neurons; units net-suppressed (Discussion only) | **Agree.** | Early/late index for isolated units (~58, small n) and for MUA by ocular-dominance bin. |
| Fixation (vergence) disparity sets the interocular phase; late suppression could be emerging disparity tuning (Menz & Freeman 2003) | **Real alternative; register #7.** Only 211008 has both eyes recorded, so it can barely be measured. Arguments against: 26/26 consistency across sites with different RFs, contrast dependence, MUA pools phase-invariant cells. | Blake: does any other session have the second eye? Discuss explicitly in the paper; test with 211008 if possible. |
| ±0.5–1 mm layer uncertainty "makes laminar statements uninterpretable"; depth axes reach 3 mm | **Partly.** The range is a loose summary; the paper should report the measured marker disagreement and why apparent thickness exceeds 2 mm (penetration angle). | Quantify from alignment v2 (docs/05 C12–C13). |
| Fig 3a, 4c use raw differences, the measure the paper calls confounded | **Agree.** The confound is 1–2 ms, but the inconsistency is real. | Move to the supplement or show them latency-aligned. |
| Same-day dual-probe penetrations are not independent | **Fair.** | Mixed model with session as a random effect, or session-level tests, as a check. |
| Late index "below zero" | **Agree** (informed concern 4); the population mean at 45% late is ~0.1, the penetration median −0.12. | "Matched the dominant eye alone." |
| Parameter recovery for w > a | **Agree.** | Simulate w = a with matched noise, fit, report the w − a distribution. |
| Fig 5b outlier = 211008 (217 ms stimuli, offset inside the fit window) | **Likely right.** | Check; fit 211008 to its own stimulus window or exclude it with the reason stated. |
| Rise is inside the fit window, not a prediction ([65]) | **Right.** | Moot (Fig 6 removed). |
| Orange means "other eye" in Figs 1–3 and "fixed-delay model" in Figs 5–6 | **Right.** | Recolour models (fixed-delay to a hue not used for conditions; revalidate). |
| Interocular timing: both halves on the same raster lines | **Right, and a strength.** | State it in Methods. |
| Internal IDs in figures (221212_J_eA, 211008), "cell" for a condition, protocol language, placeholders | **Right.** | Remove or rename throughout. |
| Missing citations: Smith et al. 1997 (binocular contrast combination in monkey V1); Macknik & Martinez-Conde 2004 (dichoptic masking); Legge & Rubin 1981; Baker, Meese & Georgeson 2007 | **Check and cite.** | Literature pass (with Blake's latency request). |

## Further concerns from the informed read

| # | Concern | Where | Action |
|---|---|---|---|
| 1 | Timing claims built on threshold crossings of a noisy index (Fig 3d half-transition; Fig 4d "first window with linearity index < 0.5"); groups differ in signal-to-noise | Abstract, [43], [50]–[51] | Noise check as for rise time; drop "normalization appeared first" unless it survives. |
| 2 | Two different "a": eye weight 0.70 ([28]) and model excitatory weight 0.88, some > 1 ([57], Fig 5c) | | Rename one; explain the difference. |
| 3 | "Several milliseconds" binocular latency advantage | Abstract, [20] | 1–2 ms, produced by every model; a methods motivation, not the hook. |
| 4 | Window reconciliation of the literature | Abstract, [74]–[75] | Checked internally, not supported; remove from the abstract, report plainly in the Discussion. |
| 5 | Registration language: "registered direction tests" ([110]), "specified but not yet run" ([68]), "the rule was fixed before…" ([100]) | | Remove (Blake's rule). |
| 6 | Stereoscope sentence ([25]) | | One monitor, two halves, septum. |

## Narrative

* **Question first, confound second.** Open with "normalization takes time" and the unresolved structure of
  binocular gain control, argued so far with time-averaged measures; state what is new against Cox et al. 2019.
  The latency issue becomes one correctly sized methods paragraph.
* **Title.** "Eyes do not normalize" (cold read). Candidates: "Binocular combination in macaque V1 shifts from
  summation to interocular suppression within tens of milliseconds"; or keep the add/normalize wording if the
  null and the interocular-delay test support it.
* **Results order** (both reads): (1) design and example; (2) the population sequence with the validated measure
  (null included); (3) contrast dependence and the 6 × 6 iso-response contours, early vs late; (4) model: what
  must be delayed, w > a, the conductance model; (5) depth and input, reduced to what holds (the sequence at all
  depths and in the input current).
* **Supplement:** Fig 3d, Fig 4c–d (unless the noise check supports d), Fig 5 two-stage and Ding–Sperling,
  offset transient, single units, controls.
* **Ending.** A concluding paragraph with a testable prediction (binocular summation larger for brief stimuli than
  for sustained ones; check whether duration dependence has been measured), not a limitations list.
* **Abstract.** ~330–370 words now; cut to one paragraph of ≤ 250: question → approach → result → model →
  implication. No laminar-order claim, no 1.5% model gain.

## Figures

| Figure | Keep | Change |
|---|---|---|
| 1 | b (relabel "cell" → "condition"; map ipsi/contra to dominant/other), d (onset-aligned like Fig 2a) | a, c, f → supplement or shrink; remove session IDs; units on e. |
| 2 | a, b (the paper's best panels) | c noisy at 22.5% and transient dip at 45% → fewer, wider windows, or model-null overlay; d: report clipped counts; explain 2c vs 2d late values. |
| 3 | merge b + c into Fig 2 or a contrast figure; relabel "independent recordings" ("held-out stimulus position runs", one monkey) | a → supplement (raw difference); d → supplement. |
| 4 | a (or move to Fig 1 instead of 1c), b with the mixed-model result stated honestly | c → supplement; d → supplement unless the noise check holds. |
| 5 (merged) | see the parked plan below | add a model-predicted index time course (static vs delayed) panel. |
| New | iso-response contours, early vs late | from the 6 × 6 grid, in contrast units per window. |

## Analyses before writing V3 (each specified in docs/09 before it is run)

1. **Validity:** fixed-rule null (w = a model + noise) through `windowed_L`; cross-validated dominance; MUA by
   ocular-dominance bin; isolated units.
2. **Interocular delay:** static vs delayed index time courses; static vs delayed on the interaction loss;
   separate monocular and interocular pool delays; parameter recovery for w > a.
3. **Iso-response contours** across the 6 × 6 grid, early vs late.
4. **Timing noise checks:** Fig 3d, Fig 4d.
5. **Fixation disparity:** 211008 (both eyes recorded) if usable; otherwise argue in the Discussion.
6. **Literature:** binocular latency advantage (Blake's comment); contrast–latency; duration dependence of
   binocular summation; the cold read's missing citations.

## Blake's new comments on V2

* **Authors:** Blake Mitchell, Brock Carlson, …, Alex Maier (last); Kacie Dougherty and/or Michele Cox possibly
  later, not yet.
* **Latency literature (abstract):** cite the contrast–latency work and check whether binocular vs monocular
  latency has been characterized (human VEP, cat, macaque); cite it if so; the argument holds either way.

## Parked: Figure 5 restructuring (agreed 2026-10-31, after the merged figure was built)

* Row 1: (a) simplified model schematic; (b) data vs conductance-model fit in the Fig 2a format (dominant 90%,
  dominant 45%, other 45%, both 45%), plotted first to check that the model does the paper's main effect; (c)
  held-out comparison of the three models with plain labels ("conductance better than fixed-delay",
  y = % less held-out error).
* Row 2: (d) static vs delayed, now as predicted index time courses as well as held-out error; (e) w vs a; (f)
  population time courses and rise times.
* Old panel g (late linearity index, model − data) → supplement with the model variants.
* Recolour the models so none uses orange ("other eye").
* Legend ≈150 words; noise-check explanation in the Results text; state the monkey split (gain mostly in J).
* Title more cautious: "A delayed, conductance-like normalization pool best predicts binocular responses".
