# 09 — Paper design: binocular combination over time and depth in V1 (draft, 2026-10-29)

Status: design draft for discussion. Builds on the literature sweep
([reports/V1 binocular contrast normalization.md](../reports/V1%20binocular%20contrast%20normalization.md)) and
first-pass feasibility ([08](08_binocular_feasibility.md)). Nothing here is a result until the figure's
"pans out if" criterion is met on the final analysis.

## The paper in one sentence

In awake macaque V1, the two eyes' signals **add almost linearly at response onset at every depth**, and
**interocular normalization takes over tens of milliseconds later**; the depth profile and the current-source
density say where that normalization comes from.

Why it matters: the field's disagreement about binocular gain control (monocular-first gain control in cat simple
cells; a fully binocular normalization pool in human fMRI/EEG; a two-stage model in psychophysics) may be a
disagreement about **when** and **where** one looks, not about mechanism. No study has combined millisecond
timing, laminar resolution and a full dichoptic contrast matrix. This dataset has all three, on a corrected clock
(05 A2, A7).

Working titles: *"The two eyes add first and normalize later in primate V1"*; *"Binocular summation precedes
interocular normalization across the layers of macaque V1"*.

## Core figures (essential: 1–3; strengthening: 4–5; optional: 6)

Each figure answers one question. The criteria are set now, before the final analysis, so that a null is
interpretable rather than explained away.

### Figure 1 — Design, alignment and data quality (setup)

* **Panels:** (a) stereoscope and 6×6 dichoptic contrast design; (b) example penetration: flash CSD, MUA
  profile, top of cortex, L4/5 landmark, depth axis; (c) photodiode-locked timing (stimulus onset error removed);
  (d) example channel: monocular and dioptic responses over time; (e) population monocular CRFs per eye;
  (f) ocular dominance distribution (MUA, units).
* **Message:** the measurement is clean enough to resolve milliseconds and depth.
* **Status:** all inputs exist (datasets, laminar alignment v2, photodiode).

### Figure 2 — The eyes add at onset, then normalize (the central figure)

* **Question:** does the rule by which the two eyes combine change over the course of the response?
* **Key measure: equivalent contrast over time, k(t).** For each channel and millisecond, k(t) is the factor by
  which the dominant eye's monocular contrast must be raised to match the binocular response:
  R_mono(k·c, t) = R_bin(c, c, t). Comparing at matched response, not matched time, makes it immune to the
  latency confound (binocular input is effectively higher contrast and so starts earlier; an earlier rising edge
  alone makes a transient "facilitation"). Predictions: linear summation k = 1 + a (a = non-dominant/dominant eye
  weight, ≈ 2 for balanced sites); normalization before summation k = 2^(1/n) ≈ 1.4; ocularity invariance k → 1.
* **Panels:** (a) population monocular vs binocular time courses; (b) binocular − monocular difference by
  second-eye contrast; (c) k(t) for each contrast with the three model bands; (d) the time at which k(t) leaves
  the linear-summation band (normalization latency), relative to response onset; (e) additivity index (binocular /
  sum of monocular) in onset-aligned early vs late windows across contrast.
* **Pans out if:** k(t) starts within the linear band (penetration-level CI overlapping 1 + a) and falls by ≥ 0.3
  toward the normalization or invariance bands within ~50 ms of onset.
* **If not:** k(t) constant over time means the transient facilitation was latency alone and a single static
  normalization rule holds. That is still a clean, publishable answer to the summation-rule question (report Q2),
  just a different paper.
* **First look:** binocular − monocular difference is positive from 20–37 ms (at the monocular onset, half-max
  30–39 ms), peaks near 45 ms, and is gone by ~70 ms; additivity index 0.85 → 0.60 with contrast in the early
  window, ~0.6 throughout late. Consistent with the claim, but not yet latency-proof.

### Figure 3 — Where normalization is strongest: depth

* **Question:** does interocular normalization (its latency, its strength w, the late net suppression) vary with
  depth?
* **Panels:** (a) k(t) heatmap, depth × time; (b) normalization latency vs depth; (c) late-epoch interocular
  weight w (Moradi–Heeger form with eye-specific gains, free w, fitted per channel) vs depth; (d) late net
  suppression vs depth.
* **Axis:** depth below the top of cortex as the primary axis (robust, 26 penetrations); landmark-based
  compartments (22 penetrations, ±0.5–1 mm) only as a secondary, robustness panel.
* **Pans out if:** a depth effect with a penetration-level CI excluding zero, in either latency or w.
* **If not:** depth-invariant timing and strength means binocular normalization is a property of the column, not of
  a laminar stage. That answer is still informative, and simpler.
* **First look:** late net suppression only at 20–32 contacts below the top of cortex (−6 to −8% of peak), none
  above; MUA saturates slightly less with depth (07).

### Figure 4 — Input or output: current-source density (strengthening)

* **Question:** is the normalization already present in the earliest synaptic input (current sink), or does it
  appear only in spiking output?
* **Panels:** (a) earliest-sink amplitude for monocular vs binocular stimulation; (b) additivity index and k for
  the sink vs MUA at the same depth; (c) same over time.
* **Pans out if:** the earliest sink adds near-linearly while MUA at the same depth normalizes (output-stage,
  intracortical normalization); or the sink is already normalized (inherited, feedforward). Either outcome is a
  mechanism statement; only an unidentifiable sink fails.
* **Risk:** an identifiable earliest sink in ≤ 22 penetrations (05 C13); after ~10 ms the CSD includes
  intracortical currents.

### Figure 5 — Single neurons and ocular dominance (strengthening)

* **Question:** does the onset-sum/late-normalize sequence hold in isolated neurons, and does ocular dominance
  decide who is suppressed late? (This reconciles Dougherty et al. 2019, facilitation of monocular neurons, with
  Zhang et al. 2024, suppression of eye-preferring neurons, as an epoch difference.)
* **Panels:** (a) k or additivity index early vs late in isolated units; (b) binocular modulation vs |ODI|,
  early vs late.
* **Pans out if:** units reproduce the Fig 2 direction (pooled across units; per-unit fits are underpowered), and
  late modulation becomes more negative with |ODI|.
* **Risk:** 68 strict / 136 lenient units, ~15–50 driven per condition (05 C14). Pooled, hierarchical
  statistics only.

### Figure 6 — One model for all of it (optional)

* A delayed binocular normalization model: the eyes sum linearly; a divisive pool with both eyes (weight w) and a
  delay/time constant τ normalizes the response. Fitted per depth bin, it should reproduce Figs 2–4 with few
  parameters (w, τ, n, per-depth gain).
* **Pans out if:** cross-validated, it beats a static normalization model (same parameters, no delay) on the
  time courses.

## Decisions to lock before running (analysis plan)

* **Presentations:** completed, no saccade ≥ 1° from −200 to 500 ms (05 C15); report microsaccade-free
  presentations as a control.
* **Time windows:** aligned to each channel's own response onset (laminar onsets differ by 15–25 ms), plus fixed
  windows as a check.
* **Eye weights:** dominant eye and the eye-weight ratio a fitted from the monocular CRFs per site. Every
  combination statistic is compared with the prediction for that site's a, never with a fixed "2".
* **Statistics:** the penetration is the unit of replication (neighbouring MUA channels share signal);
  hierarchical models with depth as a covariate; monkey as a fixed effect; penetration-level bootstrap CIs.
* **Depth:** depth below the top of cortex primary; landmark compartments secondary; re-run the key results under
  each laminar marker (report question P).
* **Inclusion:** driven channels (ANOVA across monocular contrasts, p < .01); all 26 depth-usable penetrations;
  211008 (median 8 presentations per cell) reported but checked by exclusion.

## Controls

* **Latency:** handled by k(t) (Fig 2); also show onset latency vs effective contrast.
* **Vergence:** zero-disparity gratings make vergence error into interocular phase noise. Estimate its effect on
  the B days with both eyes on Blackrock; report by ocularity class.
* **Adaptation:** three presentations per trial, 317 ms apart. Regress out the previous presentation's per-eye
  contrast (report Q15) or show the result with first presentations only.
* **Pupil and arousal:** pupil size per presentation as a covariate.
* **Monkey:** B and J separately in a supplementary figure.

## Out of scope (other papers or supplementary)

* **Supersaturation:** retired (07). One sentence and a supplementary panel: no supersaturation at any depth.
* **Dual-probe inter-column timing:** companion paper (priority 3).
* **Not testable with this stimulus set:** orientation conflict and rivalry; disparity and phase; temporal
  frequency; monoptic vs dichoptic masking. These are the limitations section.

## Novelty check (Blake to confirm)

The dataset underlies Mitchell et al. 2022 (iScience) and 2023 (Curr Biol), which used pooled MUA over 50–200 ms
(ocular-dominance-weighted interocular normalization, response-gain fits). The claims above depend on
timing, depth, CSD and units, which those papers did not resolve. Also to check: the laminar chapters of the 2024
dissertation; the legacy `2024_LaminarCRFs` manuscript (BM-ephys); Cox et al. 2019 (early facilitation, late
suppression, one contrast pair); Dougherty et al. 2019 (binocular modulation of monocular L4 neurons); Carlson et al.
2025 preprint (interlaminar coherence).

## Order of work

1. Fig 2 core: k(t) with eye-weight correction; normalization latency. This decides the paper.
2. Fig 3: depth profiles of k(t), latency, w.
3. Fig 4: per-condition CSD, earliest sink.
4. Fig 5: units.
5. Fig 1 assembly; Fig 6 model if 1–3 hold.
6. Controls and supplementary.

| Figure | Script (to write unless noted) | Inputs |
|---|---|---|
| 1 | `fig1_setup.py` | datasets, laminar alignment v2, flash CSD |
| 2 | `equivalent_contrast.py` (k(t), latency) | datasets (MUA), `run_dichoptic.py` (exists) |
| 3 | extends `equivalent_contrast.py`; `fit_interocular.py` (w) | as Fig 2 |
| 4 | `dichoptic_csd.py` | datasets (CSD), laminar alignment v2 |
| 5 | `units_binocular.py` | units, `run_dichoptic.py` |
| 6 | `delayed_normalization.py` | Figs 2–3 outputs |
