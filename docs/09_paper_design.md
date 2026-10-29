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

## Narrative arc: the paper as a sequence of figures

A reader should be able to get the paper from the figure titles alone. Each figure makes one claim and sets up
the question the next one answers. Claims are written as we hope they turn out; the figure sections below say what
replaces each claim if it does not.

| | Figure claim (title) | Question it answers | Sets up |
|---|---|---|---|
| 1 | We can measure binocular combination with millisecond and laminar precision | can we see it? | when do the eyes combine? |
| 2 | The two eyes add at response onset, then normalize within tens of milliseconds | **when** | is that the same at every depth? |
| 3 | The timing is shared across depth; the late interocular suppression is strongest in deep layers | **where** | is normalization in the input or made locally? |
| 4 | The earliest synaptic input adds the eyes linearly; normalization appears in the cortical output | **how** | which neurons carry it? |
| 5 | In single neurons, ocular dominance decides who is suppressed late | **who** | can one mechanism explain it all? |
| 6 | A binocular normalization pool acting with a delay accounts for time, contrast and depth | **synthesis** | discussion: reconciling cat, human and macaque results |
| S1 | Pitfalls and confounds: none of the above is explained by latency, eye imbalance, fixation, vergence, adaptation, layer assignment or monkey | **robustness** | — |

Supplementary beyond S1: S2 per monkey; S3 no supersaturation at any depth (07); S4 the main results under each
laminar marker; S5 unit isolation and spike-time verification (05 A7, C14).

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
* **First pass, latency-proof (2026-10-29; `scripts/equivalent_contrast.py`, `v1lb.binocular`,
  [qc/equivalent_contrast_windows.csv](qc/equivalent_contrast_windows.csv)).**
  * **Method:** k computed on window means aligned to each channel's onset. Per-millisecond inversion was too noisy
    where CRFs saturate.
  * **Early:** L in the first 20 ms after onset is ~1 (linear) at every contrast (medians 0.8–1.3, CIs include 1).
  * **Late (50–200 ms after onset):**

    | Contrast | Early → late L | p (26 penetrations, Wilcoxon) | When L drops (after onset) |
    |---|---|---|---|
    | .11 | 1.18 → 0.74 | .016 | ~50–75 ms |
    | .225 | 1.20 → 0.49 | < .001 | 35–50 ms |
    | .45 | 0.97 → −0.12 | < .001 | 20–35 ms |
    | .055 | 1.02 → 0.50 | .26, ns | — |

    Normalization arrives sooner and stronger at higher contrast.
  * **Controls:** eye weights refitted per window are stable (median a 0.64–0.70), and the drops are unchanged
    (register #25). Monocular onset shortens with contrast (46 → 35 ms); binocular onset is 1–3 ms earlier than
    monocular at the same contrast, as expected from higher effective contrast.
  * **Criterion:** passes the pre-set "pans out if" at .11, .225 and .45 (drop ≥ 0.3; ≤ ~50 ms at .225 and .45,
    ~50–75 ms at .11).
  * **Weaknesses:** per-window CIs are wide, and between-penetration variance is large. Top-censoring reaches 43%
    at .45 in the first window (L is a lower bound there).
* **Depth, first pass (Fig 3).**
  * **Drop:** similar at all depths (median 0.5–0.6).
  * **Late L:** middle depths (7–21 contacts below the top) stay nearer linear (0.92) than upper depths (−0.15;
    p = .005, 25 penetrations within-penetration). Deep (0.21) vs middle is not significant (p = .23, 15
    penetrations).
  * **Mixed model** for late L ~ depth + depth²: inverted U, peaking ~1.7 mm below the top of cortex (both terms
    p < .001). This is a candidate "middle (input-dominated) depths stay linear, superficial and deep normalize"
    pattern. It needs the landmark-based compartments and CSD (Fig 4) to interpret.
* **Compartments and per monkey (2026-10-29; `scripts/fig3_depth.py`, [qc/fig3_depth.csv](qc/fig3_depth.csv)).**
  * **Compartment medians** (L early / transition / late):

    | Compartment | Early | Transition | Late |
    |---|---|---|---|
    | Supragranular | 0.70 | 0.60 | 0.35 |
    | Granular | 1.27 | 1.09 | 0.81 |
    | Infragranular | 1.15 | 0.82 | 0.75 |

    Supragranular sites are the least linear already at onset.
  * **Within penetrations: not significant.** Granular is more linear than supragranular in 63% of 19
    penetrations (p = .09–.40); granular vs infragranular is inconsistent.
  * **The inverted U is J only:** J p = .001 (18 penetrations, peak ~2.1 mm below the top of cortex); B p = .29
    (8 penetrations). This fails the both-monkeys bar.
  * **Verdict:** Fig 3's pre-set criterion is not met robustly. The data favour the fallback: the early–late
    drop is similar at all depths, so the timing and size of interocular normalization look columnar.
    Superficial sites start less linear, consistent with receiving already-combined signals (hypothesis).

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
* **First pass (2026-10-29; `scripts/dichoptic_csd.py`, [qc/dichoptic_csd_sinks.csv](qc/dichoptic_csd_sinks.csv)).**
  * **Sink selection:** an earliest grating sink was found in all 26 penetrations. "Earliest" means the
    earliest-onset channel among strong sinks; the strongest sink is often a later, intracortical one.
    * Mostly 1–13 contacts below the top of cortex and 0.5–8.5 contacts above the L4/5 landmark. Two are deep,
      consistent with L6 input.
    * Onsets mostly 23–44 ms after the flash.
  * **Eye weights:** early sink a ≈ 0.77, the same as MUA at the same channels (0.70). The input is not monocular
    at this resolution.
  * **Combination:** at onset the sink adds the eyes roughly linearly. k at 0–20 ms after sink onset is 1.4–2.2,
    against ~1.8 predicted for linear summation. It then normalizes, like spiking: k 0.4–0.9 from 50 ms, and
    0.14–0.24 at .45.
  * **Verdict so far:** no evidence that normalization is present in the earliest input. That argues against
    geniculate inheritance at onset. Whether the sink normalizes later than spiking is unresolved: it is tentative
    at .45, and the windows are timed from different onsets.
  * **Method notes:** use k, not L = (k − 1)/a, for the CSD, because the sink's eye weight collapses late (register
    #30). Per-penetration CSD inversions are noisy; grand averages fail because of the late eye-weight collapse.
  * **Same clock (2026-10-29):** MUA at the sink channels starts a median 2 ms after the sink (IQR −2 to +6).
    Timed from the sink onset, k for the sink and for MUA do not differ consistently in any window at .225 or .45
    (sink more linear in 40–65% of penetrations, all p > .05). The largest gap is at .45, 35–50 ms: 1.03 vs 0.37,
    p = .053.
  * **Verdict:** between the pre-set outcomes. The input adds linearly at onset (not inherited), then
    normalization appears in currents and spikes together. That points to recurrent intracortical normalization:
    recurrent currents join the input within ~20–35 ms. It is a weaker statement than "output only", and the figure
    could be folded into Fig 2 or 3 rather than standing alone.
  * **Still open:** a sink-depth breakdown (granular vs L6 sinks).

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
* **First look (2026-10-29; `scripts/fig5_units_od.py`, [qc/fig5_units_od.csv](qc/fig5_units_od.csv)).**
  * **Units:** 62–65 lenient units (20 penetrations) and 23–25 strict (13–14) driven. Units are mostly net
    *suppressed* by binocular stimulation (median BMI −0.05 to −0.15; 29–45% facilitated). MUA is mostly
    facilitated (median +0.01 to +0.06; 53–69%). This is the aggregation effect Zhang et al. 2024 describe.
  * **BMI vs |ODI|:** in MUA, BMI falls with |ODI| in *both* epochs (transient 40–140 ms and sustained 150–250 ms;
    ρ −0.09 to −0.34, mixed-model slopes p < .05, strongest at .9). So the dependence on ocular dominance is not
    late-specific in these windows. In units the relationship is underpowered (mostly ns).
  * **Verdict:** the Dougherty vs Zhang "epoch reconciliation" is not supported in fixed windows. An onset-aligned
    first-20-ms epoch is the remaining test. Fig 5 is at best supporting material.

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

## Pitfalls and confounds register (basis of Supplementary Figure S1)

A living list: add every new suspect as it comes up. "Threatens" names the figure whose claim it could fake or
hide. Each item gets a test, and the S1 column says where the answer is shown.

| # | Pitfall / confound | How it could fake or hide the result | Threatens | Test or control | Shown in |
|---|---|---|---|---|---|
| 1 | **Latency shift with effective contrast** | Binocular input is effectively higher contrast, so it starts earlier; an earlier rising edge alone gives a transient "facilitation" with no change in the combination rule | 2, 3 | k(t) compares at matched response, not matched time; onset latency vs monocular contrast and vs binocular condition; onset-aligned windows | S1a |
| 2 | **Unequal eye weights in "balanced" MUA** | \|ODI\| < 0.3 allows eye weights ~2:1; with an accelerating CRF this alone pushes summation statistics below linear | 2, 3 | per-site eye-weight fit (contrast scale a); every statistic compared with that site's prediction | S1b |
| 3 | **Output nonlinearity** | An expansive or saturating CRF changes the additivity index without any change in how inputs combine | 2 | k(t) is invariant to the output nonlinearity; additivity index reported only alongside k | S1b |
| 4 | **Ceiling and the contrast range** | k·c above .9 cannot be measured, so k is censored at high contrast; interpolation between 6 levels adds error | 2, 3 | report k only where k·c ≤ .9; treat censoring explicitly; interpolation in log contrast, checked against Naka–Rushton inversion | S1b |
| 5 | **Noise bias in ratios and maxima** | Ratios of small responses near threshold and maxima of noisy means are biased (cf. supersaturation, 07) | 2, 3 | cross-validated estimates; bootstrap over presentations and penetrations; report low-contrast cells separately | S1b |
| 6 | **Fixation instability and microsaccades** | Saccade-locked modulation differs between conditions or epochs; B is a less steady fixator (9–24% of presentations with a large saccade) | 2, 3, 5 | large saccades excluded; results repeated on microsaccade-free presentations; saccade rate by condition | S1c |
| 7 | **Vergence error** | Zero-disparity gratings turn vergence error into interocular phase noise, which phase-sensitive cells read as facilitation or suppression | 2, 5 | B days with both eyes on Blackrock: vergence vs binocular modulation; results by ocularity class | S1d |
| 8 | **Adaptation from the previous presentation** | Three presentations per trial ~317 ms apart; the previous stimulus's per-eye contrast changes the current response | 2, 3 | regress out previous per-eye contrast; repeat with first presentations only | S1e |
| 9 | **Interocular crosstalk in the stereoscope** | Light from one half reaching the other eye makes "monocular" stimuli partly binocular and shrinks apparent summation | 2, 5 | strictly monocular units (expected in L4) should give ~0 to the other eye; check the septum geometry | S1f |
| 10 | **Interocular onset asymmetry** | If one eye's grating appeared earlier (display scan), binocular timing effects would be partly an eye delay | 2 | gratings sit on the same screen rows in both halves, so they share scan time (verify from stimulus positions); monocular latencies per eye | S1a |
| 11 | **Eye labelling (contra/ipsi)** | A swapped eye code would invert ocular dominance analyses | 5 | verified mapping (raw 3 = left, 2 = right; hemisphere swap); granular OD should favour the contralateral eye across penetrations | S1f |
| 12 | **Layer assignment** | Labels uncertain by ±0.5–1 mm (05 C13); a depth effect could be a mislabelled boundary | 3, 4 | depth below the top of cortex as the primary axis; rerun under each laminar marker | S4 |
| 13 | **Top-of-cortex error** | Shifts the whole depth axis per penetration and blurs depth effects | 3 | sensitivity analysis with ±2-contact shifts; penetration-level random effects | S4 |
| 14 | **CSD artefacts** | Contact gain mismatch makes sink/source stripes (05 C12); volume conduction | 4 | bad-contact interpolation and Hamming smoothing; results with and without smoothing | S1g |
| 15 | **MUA ocularity dilution** | Analog MUA pools neighbouring columns, so OD-specific effects are diluted | 5 | single units; threshold-crossing MUA as a check | Fig 5 |
| 16 | **Unit isolation and spike timing** | Contaminated clusters mix neurons; legacy spike times were offset (05 A7) | 5 | strict vs lenient units; corrected spike times verified on raw waveforms | S5 |
| 17 | **Pseudo-replication** | Neighbouring channels share signal; treating channels as independent inflates significance | all | the penetration is the unit of replication; hierarchical models; penetration bootstrap | all |
| 18 | **Monkey and session differences** | Species (M. radiata vs M. mulatta), hemisphere, stimulus duration (217 vs 317 ms, 05 A5), fixation quality | all | per monkey; duration-limited windows for early B; monkey as a fixed effect | S2 |
| 19 | **RF position and size relative to the grating** | Grating centred on one probe's RF; off-centre or larger-RF sites get surround suppression, which could look like depth-dependent normalization | 3 | RF maps (dotmapping) as a covariate; restrict to sites with RF inside the grating | S1h |
| 20 | **Orientation and SF tuning** | One grating per session; poorly tuned sites (more common in some layers) respond differently | 3 | tuning from the rfori runs as a covariate where available | S1h |
| 21 | **Trial selection** | Completed presentations only; if fixation breaks depend on condition, the surviving sample is biased | 2, 3 | break rate by condition; analysis including broken presentations up to the break | S1c |
| 22 | **Absolute timing reference** | The photodiode sits at the bottom right; any scan offset to the stimulus row is constant (≤ 3.4 ms, 05 A2b) | 1, 2 | affects absolute latencies only, not condition differences; state it | Methods |
| 23 | **Arousal and drift within session** | Slow changes in responsiveness or pupil interact with condition order | 2, 3 | randomised order; time in session and pupil as covariates | S1e |
| 24 | **Analysis flexibility** | Many windows, thresholds and depth bins | all | decisions locked in this document before the final run; cross-validated model comparison | Methods |
| 25 | **Time-varying eye weights** | If the non-dominant eye's monocular response fades faster, L = (k − 1)/a drops late even under linear summation | 2, 3 | eye weight refitted in each window. **Tested: a is stable (0.64–0.70) and the early–late drops are unchanged (L_w)** | S1b |
| 26 | **Small eye weights amplify L** | Strongly monocular channels (small a) make L = (k − 1)/a noisy and heavy-tailed | 2, 3 | medians over channels and penetrations; repeat with a ≥ 0.3; report k alongside L | S1b |
| 27 | **Top-censoring at the rising edge** | The binocular response can exceed even the .9 monocular response early (it starts earlier and rises faster); k is then only a lower bound (up to 43% at .45, 0–10 ms) | 2 | treat censored values as lower bounds (medians); report the censored fraction per window; primary contrasts .11 and .225 | S1a |
| 28 | **Inverting near-flat CRFs** | Where the monocular CRF saturates, small response changes map to large changes in k (noisy, biased late values) | 2, 3 | window means rather than per-ms values; primary contrasts below saturation; CIs from penetration bootstrap | S1b |
| 29 | **Strongest sink ≠ earliest sink** | The strongest grating sink is often a later intracortical one; testing it as "input" would misplace the input stage | 4 | earliest-onset sink among strong sinks; sink onset must not lag the MUA onset at the same channels | S1g |
| 30 | **Late collapse of the sink's eye weight** | L = (k − 1)/a divides by a small late a and explodes for CSD | 4 | use k for CSD; report a(t) | S1g |

Supplementary Figure S1 layout (one panel per row group): **a** latency and onset symmetry; **b** eye weights,
nonlinearity, censoring and noise; **c** fixation, microsaccades and trial selection; **d** vergence; **e**
adaptation and drift; **f** monocular crosstalk and eye labelling; **g** CSD processing; **h** RF and tuning
covariates.

## Out of scope (other papers or supplementary)

* **Supersaturation:** retired (07). One sentence and a supplementary panel: no supersaturation at any depth.
* **Dual-probe inter-column timing:** companion paper (priority 3).
* **Not testable with this stimulus set:** orientation conflict and rivalry; disparity and phase; temporal
  frequency; monoptic vs dichoptic masking. These are the limitations section.

## Competitive landscape (checked 2026-10-29)

* **Risk:** Neuropixels labs in macaque V1 (20 µm spacing, hundreds of units, sublaminar resolution, e.g. the
  [laminar identification study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10849622/), the Carr et al.
  [cell-class study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12132580/), anesthetized
  [orientation across layers](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11456443/)) would beat this dataset on
  unit counts and on 4A/4B/4Cα/4Cβ resolution.
* **No dichoptic Neuropixels work found.** The search found no Neuropixels study with dichoptic or stereoscopic
  contrast control in macaque V1 through October 2026. The nearest interocular-normalization work is population
  imaging: Reynaud, Roux, Chemla, Chavane & Hess ([AES abstract](https://aes.amegroups.org/article/view/5220/html):
  VSDI and frequency tagging, superficial, mechanism explicitly open), and Zhang et al. 2024 (two-photon,
  superficial, one contrast).
* **Edge:** a full dichoptic contrast matrix, laminar timing on a corrected clock, and CSD, across 26
  penetrations in two animals.
* **Implications for the design:**
  * Make the claims rest on population timing, depth and CSD (Figs 2–4), where 100 µm MUA/CSD is sufficient.
  * Keep single units as support (Fig 5), not as the evidence base.
  * Make no sublaminar claims (no 4Cα vs 4Cβ).
  * Treat Figs 1–3 as the minimum publishable unit, so a short-format paper is possible if speed matters.

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
