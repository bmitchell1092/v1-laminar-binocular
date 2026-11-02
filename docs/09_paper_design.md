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

### Revised arc after the first passes (2026-10-29; approved by Blake)

The first passes moved the weight of the paper from depth to **time × contrast**: Fig 2 holds, Fig 3 favours
the columnar fallback, Fig 4 lands between its outcomes, Fig 5 is supporting at best. A revised arc that uses
only what held:

| | Figure claim (title) | Evidence so far |
|---|---|---|
| 1 | We can measure binocular combination with millisecond and laminar precision | setup |
| 2 | The two eyes add at response onset, then normalize within tens of milliseconds | latency-proof k(t); linear at onset at all contrasts; early→late drops p ≤ .016 |
| 3 | Contrast sets how far and how fast binocular summation collapses | falls with contrast in 26/26 penetrations, both monkeys, odd/even; visible in raw data; held-out test 6/7 (p = .063, failed narrowly); L < 0 at .45 (net interocular suppression beyond invariance) |
| 4 | The same sequence runs at every depth and in the earliest synaptic input | drop similar across depth; earliest sink linear at onset; currents and spikes normalize together |
| 5 | A contrast-driven, delayed binocular normalization pool explains it | to do (model) |

Supplementary: units (neurons net-suppressed while MUA is facilitated; OD dependence early and late), depth
compartments and the J-only inverted U, per monkey, no supersaturation, the confounds register.

**Decisions (Blake, 2026-10-29):**
* **Claim 3** is presented, with the failed held-out test stated plainly (6/7 negative, p = .063, J only).
* **Depth** is reframed as columnar.
* **Fig 5** is the dynamic normalization model.

Supplementary beyond S1: S2 per monkey; S3 no supersaturation at any depth (07); S4 the main results under each
laminar marker; S5 unit isolation and spike-time verification (05 A7, C14).

**Sensitivity: every result after MUA artifact rejection (2026-10-29; register #34, 05 C16).** All analyses were
rerun with `clean_presentations`. Pre-registered results stay as run; these reruns are sensitivity analyses. The
figures (`scripts/figures/`) use the cleaned outputs.

| Result | As first reported | After artifact rejection |
|---|---|---|
| Claim 2, early → late L (.11 / .225 / .45) | 1.18 → 0.74 (p = .016) / 1.20 → 0.49 / 0.97 → −0.12 | 1.19 → 0.61 (p = .001) / 1.20 → 0.47 / 0.97 → −0.12 (both p < .001) |
| Claim 3, transition slope < 0 | 26/26, median −0.92 | 26/26, median −0.95; odd/even 25/25 and 23/24 |
| Claim 3, held-out (registered) | 6/7, p = .063, **failed** | 7/8, p = .035 (221209_J_eB now has driven channels); sensitivity only |
| "Sooner at higher contrast" | not supported (p = .23) | not supported (p = .25; Friedman p = .19) |
| Laminar timing, normalization above / below the sink | +5.5 ms (p = .026) / +6.5 ms (p = .006) | +6 ms (p = .016) / +7 ms (p = .003) |
| Delayed vs static; w > a | 23/25; 24/25 | 23/25; 24/25 |
| Two-stage vs single-stage, raw error (registered) | PASS, 21/24, ~1.1% | PASS, 23/24, ~1.5% (538 channels) |
| Two-stage vs single-stage, interaction loss | p = 7 × 10⁻⁶ | p = 6 × 10⁻⁸ |
| Capture of the late collapse (registered) | no model | no model |
| Rise time .11 → .9, data vs fitted model | 48 → 20.5 vs 20 → 11 ms | 46.5 → 20 vs 20 → 11 ms |

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
* **Laminar timing check (2026-10-29; flagged by the circuit literature, where intracortical suppression reaches
  L4C last).** Channels are grouped by position relative to each penetration's earliest current sink (±2
  contacts), using .225 and .45 and window-specific eye weights; 20–25 penetrations.
  * **The sink is the input:** its channels respond first, at 31 ms after the flash, vs 35.5 ms above (+4 ms,
    p = .02) and 36 ms below (+2.5 ms, p = .008).
  * **Each depth normalizes a fixed ~27 ms after its own onset:** first window with L < 0.5 at 27.5 ms after
    onset, both at the sink and above it (no difference). Below the sink it takes longer (42.5 ms).
  * **On the flash clock:** normalization arrives first at the input sink (~63 ms), then superficial channels
    ~5.5 ms later (p = .026), then deep channels (~74.5 ms).
  * **Interpretation:** interocular normalization follows the feedforward activation sequence with a fixed
    local lag, starting in the input layer. That is unlike intracortical near-surround suppression, which reaches
    L4C last (Bijanzadeh 2018).
  * **Confidence:** modest. The windows are coarse (15 ms) and p = .026 for the upward step. A finer,
    pre-registered timing analysis is needed before claiming it. A candidate for Figs 3–4.

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
* **First pass (2026-10-29; Fig 5 in the revised arc; `v1lb.dynorm`, `scripts/fig5_model.py`,
  [qc/fig5_model.csv](qc/fig5_model.csv)).**
  * **Model:** after Zhou, Benson, Kay & Winawer 2019. Excitation (c_D + a·c_N)·u(t); pool (c_D + w·c_N)·v(t),
    with v a low-passed u (τ₂); R = g·E^n / (σ^n + P^n). Fitted per penetration to the channel-averaged 6×6 time
    courses (0–219 ms). Cross-validated across odd/even presentations.
  * **Passes:**
    * Beats static normalization on held-out halves in 23/25 penetrations (B 5/7, J 18/18); held-out error
      2.95× lower; p < .001.
    * Beats w = a in 21/25 (median 2.4% lower error, p < .001).
    * **w > a in 24/25** (B 6/7, J 18/18): median w = 1.08 vs a = 0.88. Ocular dominance weights excitation;
      the pool is close to eye-balanced. This is the quantitative test of the "OD sets the numerator, the pool
      stays balanced" reading (Mitchell 2023, report Q5), and the condition for binocular < dominant eye at high
      contrast.
    * Variance explained median 0.88. Normalization delay τ₂ median 30 ms (IQR 23–37), matching the observed
      transition. n ≈ 1.1, σ ≈ 0.09, latency δ ≈ 21 ms.
  * **Shortfall:** the model reproduces the direction and contrast ordering of L(t) (ρ 0.23–0.56, all p ≤ .002)
    but is too slow and too shallow. At .45 the data reach L ≈ −0.4 by 40 ms; the model reaches ~0 only by ~90 ms.
    At .11 the model stays near 1 while the data fall to ~0.5. It also underpredicts the early dioptic peak at high
    contrast. A fixed pool delay explains *that* normalization is delayed, not that it is faster at high contrast.
  * **Tested and rejected: a pool time constant that shortens with contrast** ('tau_c', τ₂ / (1 + β·pool
    contrast); Bair & Movshon 2004). It beats the fixed delay in only 16/25 penetrations, with no gain in held-out
    error (median ratio 1.000, p = .045). β ≈ 0 (median 0.02), and L(t) is reproduced no better. So the faster,
    deeper collapse at high contrast is not a contrast-dependent pool speed.
  * **Where the misfit points:** at .11 the data fall to L ≈ 0.5 while the model stays ~1, because its pool is
    below σ there. That suggests interocular suppression already acting at low contrast: a separate early
    eye-vs-eye gain-control stage, as in the psychophysical two-stage model (Meese, Georgeson & Baker 2006). Next
    model: a dynamic two-stage model. Stage 1 divides each eye by S + both eyes' contrasts, with its own delay;
    stage 2 is the delayed binocular pool. Compare with the single-stage model on held-out halves. This is report
    Q2 (the summation rule) in dynamic form.
  * **Narrative (Blake, 2026-10-29):** if it holds, the paper can make the point that the single-stage model fails
    in a specific way and a *dynamic two-stage* model is the solution.
  * **Pre-specified comparison (written before fitting the two-stage model).**
    * **Single-stage ('full'):** E = (c_D + a·c_N)·u; P = (c_D + w·c_N)·v; R = g·E^n / (σ^n + P^n). 8 parameters.
    * **Two-stage ('2stage'):**
      * Stage 1, per eye: b_D = c_D·u / (S + v_s·(c_D + ω·c_N)) and b_N = a·c_N·u / (S + v_s·(c_N + ω·c_D)),
        where v_s is u low-passed with τ_s.
      * Stage 2: B = b_D + b_N; R = g·B^n / (σ^n + (low-passed B, τ₂)^n).
      * 10 parameters. It nests the single-stage model as S → ∞, except that the pool weight comes from stage 1
        instead of a free w.
    * **Success criterion (all three must hold):**
      1. Two-stage has lower cross-validated (odd/even) error in ≥ 70% of penetrations, Wilcoxon p < .05.
      2. It holds in both monkeys (B ≥ 4/7).
      3. It reproduces L(t) better: lower median |L_data − L_model| at .11 and at .45.
    * **If it fails:** the paper reports the single-stage delayed model with its stated shortfall, and two-stage is
      not claimed.
    * **Result (2026-10-29): FAILS the pre-specified criterion.**

      | Criterion | Result | Verdict |
      |---|---|---|
      | 1. CV wins in ≥ 70% of penetrations, p < .05 | 68% (17/25), p = .039; median CV error only 1% lower | fails, narrowly |
      | 2. Both monkeys | B 5/7, J 12/18 | passes |
      | 3. Better L(t) at .11 and .45 (median abs error) | .45: 0.56 → 0.43 better; .11: 0.48 → 0.50 worse | fails |

      Two-stage is therefore **not claimed**; the paper reports the single-stage delayed model with its shortfall.
    * **Exploratory only (not claimable):**
      * Rank agreement of L(t) improves at .11 (ρ 0.24 → 0.54) and .225 (0.38 → 0.55).
      * Fitted stage-1 interocular suppression weight ω ≈ 0.95 (IQR 0.76–1.14), close to the weight ~1 that
        psychophysics requires.
      * Stage-1 semi-saturation S ≈ 0.13; stage-1 delay τ_s ≈ 36 ms; stage-2 delay τ₂ ≈ 34 ms.
    * **What a fair second test would need:** an analysis fixed in advance and powered for it, e.g.
      channel-level hierarchical fits rather than channel averages, or new data with denser low contrasts. Not a
      re-run of this comparison with a relaxed criterion.
    * **Flaws in the first pre-registration (identified 2026-10-29, after the result; Blake raised the question).**
      1. **Not nested, and handicapped.** The single-stage model had a free pool weight w, strongly supported
         (w > a in 24/25). The two-stage model had none, so it nested only the w = a single-stage model. It fit
         even the *training* data worse than single-stage 'full' in 3/25 penetrations. The optimizer was not the
         problem: training SSE ≤ its nested special case in 25/25.
      2. **Criterion 3 used L.** L is a derived, heavy-tailed ratio that the models were not fitted to, and it is
         very noisy per penetration × window. It largely measures noise.
      3. **The 70% threshold was arbitrary** (17/25 vs 17.5/25) and redundant with the significance test.

      The first result stands as reported (failed as specified).
  * **Second, corrected pre-registration: channel-level two-stage test (registered 2026-10-29, before implementing
    or fitting).** Blake approved trying it.
    * **Models:**
      * Single-stage 'full' (as before, 8 parameters).
      * Two-stage 'full2' (11 parameters). Stage 1 is unchanged: b_D = c_D·u / (S + v_s·(c_D + ω·c_N)) and
        b'_N = c_N·u / (S + v_s·(c_N + ω·c_D)). Stage 2 has numerator B = b_D + a·b'_N and pool
        P = lowpass(b_D + w·b'_N, τ₂); R = g·B^n / (σ^n + P^n).
      * 'full2' nests 'full' exactly as S → ∞.
    * **Optimizer:** each 'full2' fit includes a start at the fitted single-stage solution embedded at the S bound,
      plus random starts. More starts than the first test for both models.
    * **Data:**
      * Each driven in-cortex MUA channel separately, under the same driven rule as before.
      * Oriented to its dominant eye and normalised by its dominant .9 peak.
      * 0–219 ms, all 35 non-blank cells.
      * Excluded a priori: 211008 (~4 presentations per cell per half; register #32) and 211210 (2 driven
        channels).
    * **Cross-validation:** fit odd presentations and score even ones, and the reverse; CV = summed held-out SSE.
      Per channel, r = log(CV_full / CV_full2), so positive favours two-stage.
    * **Unit of replication: the penetration** (median r over its channels).
    * **Pass criterion (both must hold):**
      1. One-sided Wilcoxon signed-rank test of penetration medians > 0, at **α = .025** (Bonferroni for a second
         attempt at the same question).
      2. The median of penetration medians is > 0 in both monkeys separately.
    * **Reported regardless (no pass/fail):** median % change in held-out error; fraction of channels favouring
      two-stage; fitted ω, S and τ_s; L(t) agreement.
    * **Reporting:** both tests appear in the paper, the first with its flaws stated. If this one fails, two-stage
      is not claimed and there is no third attempt on these data.
    * **Result (2026-10-29): PASS** (`scripts/fig5_channel_test.py`,
      [qc/fig5_channel_test.csv](qc/fig5_channel_test.csv)).
      * **Primary test:** 533 channels in 24 penetrations. Penetration medians of r > 0: one-sided Wilcoxon
        **p = 3 × 10⁻⁶** (α = .025).
      * **Each monkey:** median of penetration medians B +0.006 (5/6 penetrations favour two-stage), J +0.015
        (16/18).
      * **Breadth:** two-stage favoured in 21/24 penetrations and 76% of channels.
      * **Seeding worked:** training error of 'full2' ≤ 'full' in 100% of channels.
      * **Effect size, stated plainly: consistent but small.** Median held-out error is 1.1% lower.
      * **Two-stage parameters (per channel, median [IQR]):**
        * interocular suppression weight ω 0.86 [0.52, 1.33], near the ~1 psychophysics requires;
        * stage-1 semi-saturation S 0.18 [0.09, 0.37];
        * stage-1 delay τ_s 31 ms [17, 65]; stage-2 delay τ₂ 36 ms [19, 74];
        * pool weight w 1.12 vs excitatory weight a 0.92; w > a in 75% of channels.
    * **Claim for the paper:**

      > A dynamic two-stage model, with early interocular gain control (suppression weight ~0.9, ~30 ms delay)
      > before delayed binocular normalization, predicts held-out responses better than single-stage delayed
      > normalization, consistently across channels, penetrations and both monkeys, though by a small margin
      > (~1%).

      State both tests: the first (channel averages, non-nested, L-based criterion) failed narrowly, and its
      flaws are given above; this corrected one passed.
    * **Descriptive check: neither model reproduces the late binocular collapse** (`scripts/fig5_L_check.py`,
      [qc/fig5_L_check.csv](qc/fig5_L_check.csv); per-channel fits to all presentations, penetration medians).
      * **At onset** both match the data (L ≈ 1); two-stage is slightly closer.
      * **Late, both predict too linear a combination, and two-stage is worse:**

        | Late (≥ 50 ms) | Data | Single-stage | Two-stage |
        |---|---|---|---|
        | .11 | 0.62–0.64 | 0.97–0.99 | 1.15–1.19 |
        | .225 | 0.28–0.53 | 0.64–0.74 | ~1.02 |
        | .45 | −0.50 to 0.22 | 0.17–0.29 | 0.30–0.68 |

        Signed error late at .11: +0.29 single-stage, +0.53 two-stage.
      * **Interpretation:** the 1% held-out gain comes from time-course shape, not from the interocular
        interaction the paper is about. L is sensitive where CRFs flatten (a few % response error becomes a
        large L error), but the bias is systematic.
      * **Consequence for the paper:** "two-stage solves it" is only half-supported.
    * **Model claims that survive:**
      1. Delayed normalization is required (static is ~3× worse on held-out data).
      2. The pool weights the non-dominant eye more than excitation does (w > a).
      3. An early interocular stage improves held-out prediction consistently but slightly.
      4. No fitted model reproduces the strength of the late binocular collapse; it is stated as the open
         problem.
    * **Possible next model test (not run; would need its own pre-registration):** fit with a loss that weights
      the binocular interaction (binocular vs matched-monocular responses per window) rather than raw time
      courses.
  * **Third model test: interaction-targeted fits (registered 2026-10-29, before implementing; Blake approved).**
    * **Question:** can a delayed normalization model capture the *late binocular collapse* (L(t)), and with
      which architecture?
    * **Loss (all models).** Two parts, each normalised by its own data variance and number of points, then
      summed:
      * the monocular time courses (10 cells, one eye at 0);
      * the binocular interaction I(c_D, c_N, t) = R(c_D, c_N, t) − R(c_D, 0, t) for the 25 cells with both eyes
        > 0.

      This replaces raw time-course SSE, whose fit is dominated by the large onset transients.
    * **Models:**
      1. Single-stage 'full' (8 parameters).
      2. Two-stage 'full2' (11 parameters).
      3. **Dynamic Ding–Sperling 'ds'** (10 parameters), after Ding & Sperling 2006:
         * x_D = c_D·u(t) and x_N = a·c_N·u(t);
         * gain-control energy ε_e(t) = (ρ·lowpass(x_e, τ_g))^γ;
         * Γ(t) = Σ_e x_e·(1 + ε_e) / (1 + ε_D + ε_N);
         * R = g·Γ^n / (σ^n + lowpass(Γ, τ₂)^n).

         Γ(I, 0) = I and Γ(I, I) = I(2 + 2ε)/(1 + 2ε) match the identities stated in Ding & Sperling. The full
         equation form is to be verified against the PDF before reporting.
    * **Data:** per channel, as in the channel-level test (driven channels; 211008 and 211210 excluded);
      odd/even cross-validation; all-presentation fits for L.
    * **Primary criterion (capture), per model.** Late L (≥ 50 ms after the channel's onset) is computed from
      the model's all-presentation fit and from the data, with the same windows. Take the penetration median of
      the signed error at .11, .225 and .45. A model **captures** the late collapse iff, at all three contrasts,
      the 90% bootstrap CI (over penetrations) of the median signed error lies within ±0.25 (equivalence test).
    * **Secondary (architecture; a third look at the comparison question).** Held-out interaction loss is
      compared pairwise (ds vs full, full2 vs full, ds vs full2) with one-sided Wilcoxon tests on penetration
      medians of log loss ratios. A difference is claimed only if p < .005 and it holds in both monkeys.
    * **Reporting:** all three model tests are reported.
    * **Result (2026-10-29; `scripts/fig5_interaction_test.py`,
      [qc/fig5_interaction_test.csv](qc/fig5_interaction_test.csv); 533 channels, 24 penetrations).**
      * **PRIMARY: no model captures the late collapse** by the equivalence criterion. Median late-L signed
        error (model − data) at .11 / .225 / .45, with 90% CIs:

        | Model | .11 | .225 | .45 |
        |---|---|---|---|
        | Single-stage | −0.00 [−0.39, 0.41] | −0.19 [−0.35, 0.10] | −0.09 [−0.32, 0.18] |
        | Two-stage | +0.36 [0.10, 0.56] | +0.39 [0.00, 0.60] | +0.22 [−0.02, 0.35] |
        | Ding–Sperling | −0.19 [−0.38, 0.06] | −0.13 [−0.39, 0.15] | +0.34 [−0.08, 0.44] |

        * **The interaction loss removed the single-stage bias** (medians near 0, against +0.29 with the raw-SSE
          fit), but the CIs are too wide for equivalence.
        * **Two-stage stays too linear late.**
        * **Design weakness, stated (not used to rescue):** penetration-level L is so noisy that a ±0.25 margin
          was nearly unreachable for any model.
      * **SECONDARY:**
        * **Two-stage beats single-stage on held-out interaction loss: p = 7 × 10⁻⁶**, positive in B and J (median
          log ratio +0.008 / +0.011). Claimed under the stringent α. That makes two independent loss functions
          favouring two-stage (raw SSE p = 3 × 10⁻⁶; interaction p = 7 × 10⁻⁶), each by ~1%.
        * **Ding–Sperling** does not beat single-stage (p = .49), and two-stage beats it (p ≈ 1 for ds > full2).
          Its parameters sit at bounds (γ = 4; τ_g up to 400 ms), so it is poorly constrained by these data.
      * **What the model figure can claim:**
        1. Delayed normalization is required.
        2. The pool weights the non-dominant eye more than excitation does.
        3. Early interocular gain control (two-stage) predicts held-out responses better than single-stage,
           consistently under two loss functions and in both monkeys, by a small margin.
        4. Two-stage's advantage lies in the early interaction; it does not reproduce the late collapse, which
           remains incompletely explained.
        5. The single-stage model fitted to the interaction has no median bias, but its precision is too low to
           call it a capture.
      * **Reading:** the early binocular interaction favours an early interocular stage; the strength of the late
        collapse is not pinned down by any fitted model. No further model tests on these data without a new
        reason and a new pre-registration.
  * **Controls:** aggregation bias (register #31) is negligible. 211008 is unstable under cross-validation (#32).
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
| 31 | **Aggregation bias in model weights** | The model is fitted to the average of channels with different eye weights; averaging could itself make the pool look more balanced than excitation (w > a) | 5 (model) | simulation: channels with w = a and heterogeneous a, averaged and fitted. **Tested: bias +0.003 to +0.025** even for a spread of 0.1–1.0, against an observed median w − a of +0.22 | S1i |
| 32 | **Small trial counts destabilize cross-validation** | With ~4 presentations per cell per half (211008), held-out error of flexible models explodes | 5 | report per-penetration CV; exclude or flag penetrations with < 10 presentations per cell; results with and without 211008 | S1i |
| 33 | **Offset-decay estimates are method-sensitive; input may carry the slow decay; an OFF transient dominates** | Half-decay depends on baseline (blank vs pre-stimulus), channel set and eye pooling (MUA 27–59 ms). A slow decay already present in the input currents would not reflect cortical recurrence. The response to stimulus offset peaks at ~1.8× the sustained level ~7 ms after the response offset (88–100% of penetrations), so the half-decay mostly measures the OFF transient; model comparisons must use the same metric (the toy scan of the second report did not) | 5 (model) | fixed measurement protocol pre-registered; sink vs MUA comparison; model the input stage; separate the OFF transient (decay after the OFF peak, or ON/OFF drive) | S1i |
| 34 | **Unrejected MUA artifacts** (found 2026-10-29 while building Fig 1; 05 C16) | ~0.6% of presentations carry MUA excursions of robust z 20–2000, almost all in monkey J (up to 62 in one penetration). One such presentation can dominate a ~30-presentation cell mean, add heavy tails to k and L, and distort single-penetration panels | all | rule fixed before seeing its effect: drop presentations whose peak MUA (−50…450 ms) on any in-cortex channel exceeds robust z 20 (`v1lb.dataset.clean_presentations`); every analysis rerun; pre-registered results kept as run and the cleaned reruns reported as sensitivity analyses | S1c |

**S1 results, first pass (2026-10-29; `scripts/fig2_robustness.py`,
[qc/fig2_robustness.csv](qc/fig2_robustness.csv)).** The core statistics are recomputed under each control subset.

* **Claim 3 survives every control** (transition-window linearity falls with contrast):

  | Subset | Penetrations with a negative slope |
  |---|---|
  | No saccade of any size (#6) | 92% |
  | Pupil within its IQR (#23) | 92% |
  | Previous presentation absent or ≤ .11 (#8) | 100% |
  | First / second half of session (#23) | 100% / 88% |
  | Monkey B / J (#18) | 8/8 / 18/18 |
  | Without 211008 (#32); channels with a ≥ 0.3 (#26) | 100% / 100% |

  All p ≤ .0003, except monkey B alone (p = .008).
* **Claim 2** (early → late collapse):
  * **.45:** robust, p ≤ .02 in every subset.
  * **.11:** significant in most subsets, but not after a low-contrast previous presentation (p = .72), in
    the second session half (p = .33), or in monkey B alone (p = .84).
  * **.225: sensitive to windowing (#24).** With a single late window (50–200 ms) the baseline gives p = .07;
    with the median over late sub-windows (equivalent_contrast.py) it gives p < .001. The paper must state
    the window choice and show both.
* **Implication:** claim 3 is the most robust result in the paper. Claim 2 should be anchored at high contrast,
  with the lower-contrast collapse described as weaker and windowing-dependent.
* **Still to do:** vergence (#7; only 211008 has both eyes on Blackrock), RF position (#19, needs RF maps),
  monocular crosstalk (#9), trial selection (#21), CSD processing (#14).

Supplementary Figure S1 layout (one panel per row group): **a** latency and onset symmetry; **b** eye weights,
nonlinearity, censoring and noise; **c** fixation, microsaccades and trial selection; **d** vergence; **e**
adaptation and drift; **f** monocular crosstalk and eye labelling; **g** CSD processing; **h** RF and tuning
covariates.

## Out of scope (other papers or supplementary)

* **Supersaturation:** retired (07). One sentence and a supplementary panel: no supersaturation at any depth.
* **Dual-probe inter-column timing:** companion paper (priority 3).
* **Not testable with this stimulus set:** orientation conflict and rivalry; disparity and phase; temporal
  frequency; monoptic vs dichoptic masking. These are the limitations section.

## Framing: dynamic and recurrent normalization (Blake, 2026-10-29)

* **Heeger's remark.** When Blake showed David Heeger preliminary data, Heeger said: "Hmm, it seems like
  normalization takes time." That is our central result (claims 2 and 3), and it calls for a dynamic or recurrent
  normalization framing, and for relating the results to recurrent feedback in cortex and to the SSN, including
  SSN models across the column.
* **Why it may explain the open problem.** Our fitted models put normalization behind a *fixed* delay: a
  feedforward pool behind a low-pass filter. Recurrent and shunting accounts make the *speed* of normalization
  input-dependent:
  * In conductance-based normalization (Carandini & Heeger 1994; Carandini, Heeger & Movshon 1997), more input
    means more shunting, a shorter effective time constant and faster normalization.
  * ORGaNICs (Heeger & Mackey 2019; Heeger & Zemlianova 2020) produce input-dependent dynamics as well.
  * Claim 3 (faster and further collapse at higher contrast) could therefore be a signature of the mechanism.
    The ad hoc contrast-scaled pool time constant ('tau_c') was a crude stand-in and failed.
* **Prior work.** Blake's 2024 manuscript (BM-ephys `3. visualization/2024_LaminarCRFs/Manuscript_focused_v2/
  Simulation/`) fitted the Heeger & Zemlianova ORGaNICs code (`toolboxes/ORGaNIC_Normalization`) to monocular
  CRF dynamics: `simulateV1Response.m`, Bayesian-optimisation and latency-aware cost functions. It is MATLAB, so
  it would be ported to Python if used.
* **In progress:** a second literature sweep (`research_notes/Recurrent normalization and laminar SSN/`). It
  covers dynamic and recurrent normalization, SSN dynamics with laminar/columnar SSN models, and circuit evidence
  for recurrent vs feedforward normalization.
* **Discriminating diagnostic, stated before computing (2026-10-29): decay after stimulus offset.**
  * Recurrent or shunting normalization (ORGaNICs: τ ≈ 60 ms with no input; Carandini, Heeger & Movshon 1997)
    predicts a *slow* decay that does not depend on contrast. The input-dependent speed-up vanishes once the input
    is gone.
  * A fixed-delay feedforward pool predicts a *fast* decay set by the input filter (τ₁ ≈ 4 ms in our fits), which
    can undershoot below baseline while the delayed pool is still high.
  * The same models predict onset rise times that shorten strongly with contrast (recurrent), against a
    contrast-independent rise shape (feedforward with a fixed filter).
  * Measure: per driven channel and monocular dominant-eye contrast (.11–.9), the half-decay time after
    response offset and the onset rise time (10→90%); penetration medians; test the contrast dependence.
    Descriptive, not pass/fail.
* **Diagnostic result (2026-10-29; `scripts/dynamics_diagnostics.py`,
  [qc/dynamics_diagnostics.csv](qc/dynamics_diagnostics.csv)): both signatures favour recurrent/shunting
  normalization.**
  * **Rise time (10→90%)**, dominant eye alone:

    | | .11 | .225 | .45 | .9 | Slope per log contrast |
    |---|---|---|---|---|---|
    | Data, channel-averaged | 36 ms | 26 ms | 20 ms | 18 ms | −8.6 ms (84% of penetrations; B −27.9, J −5.8) |
    | Data, per channel | 48 ms | 39 ms | 26 ms | 20.5 ms | — |
    | Fitted fixed-delay model | 20 ms | 16 ms | 13 ms | 11 ms | −4.3 ms |

    The fixed-delay model's rise also shortens with contrast, so the trend alone does not discriminate. But it is
    too fast everywhere (by ~16 ms at .11) and half as contrast-dependent. Slow integration at low input is the
    recurrent signature.
  * **Half-decay after stimulus offset.** The fitted fixed-delay model (fitted to 0–219 ms only, so the decay is
    an out-of-sample prediction) predicts **5–12 ms**. The **data show 27–32 ms**, 2.5–5× slower and nearly
    contrast-independent (Friedman p = .008, slight lengthening with contrast). This is the slow,
    input-independent decay predicted by recurrent/shunting normalization (ORGaNICs: τ ≈ 60 ms with no input).
  * **Undershoot after offset:** not interpretable. The metric is biased by noise at low sustained levels.
  * **Implication:** the fixed-delay feedforward models (single-stage and two-stage) fail an out-of-sample
    dynamic test. A recurrent normalization model is the principled next fit: Carandini–Heeger shunting, or
    ORGaNICs (Blake has legacy MATLAB fits), fitted to onset *and* offset and binocular. To be pre-registered.
* **Inheritance control (recommended by the second literature report).** If the input (earliest current sink)
  also decays slowly after offset, the slow offset is inherited and stops discriminating the models.
  * **Inconclusive.** The sink's offset decay is unstable across methods: 97–100 ms with per-trial baselines and
    both eyes pooled; 60 / 75 / 26 / 19 ms at .11–.9 with blank subtraction and the dominant eye. The input-layer
    CSD after offset is too noisy to decide.
  * **MUA half-decay at the sink channels** is 39–59 ms, depending on method, against 27–32 ms over all driven
    channels. Decay estimates are method-sensitive (register #33).
  * **Rise times**, at the same channels:

    | | Rise at .11 → .9 | Slope per log contrast | Penetrations with a negative slope |
    |---|---|---|---|
    | Spiking (MUA) | 25.5 → 16 ms | −4.4 ms | 73% |
    | Input current (sink) | 46 → 32.5 ms | −2.7 ms | 55% (no consistent trend) |

    This weakly suggests the contrast speed-up arises at or after the input rather than being inherited. Not
    conclusive.
  * **Consequence:** the offset prediction remains a valid out-of-sample test of a fitted MUA model. As
    *mechanistic* evidence for cortical recurrence it is limited unless the input stage is modelled, as the second
    literature report recommends.
* **Correction (2026-10-29, found while building Fig 6): the half-decay statistic mostly measures an offset
  transient, so the "slow decay favours recurrence" reading is withdrawn.**
  * **What the figure showed.** Aligned to the response offset (stimulus duration + onset latency) and scaled to
    the sustained level, the population response does not simply decay. It rises to a peak ~7 ms later: median
    1.75–1.9× the sustained level at every contrast, above 1.2× in 88–100% of 26 penetrations, in both monkeys.
    This is an OFF response to the grating's disappearance. The 27–32 ms "half-decay" is mostly its duration.
  * **The toy scan was not like-for-like.** The second literature report's scan (ORGaNICs, recurrent shunting,
    SSN "matching" a flat ~30 ms decay) measured half-decay from the post-offset *peak*. Measured as the data are
    (from duration + t50), the toy ORGaNICs gives 0–25 ms half-decays that lengthen with contrast. No model has
    been shown to reproduce the data's offset behaviour under the data's own metric.
  * **What stands:** the fitted fixed-delay model, measured with the same metric, rises too fast at low contrast
    (20 vs 48 ms at .11, per channel) and has no offset transient. A model with no OFF pathway cannot produce one,
    whether it is feedforward or recurrent.
  * **Fig 6** therefore shows rise times and the offset transient, not the half-decay. A future decay test must
    separate the OFF transient: measure decay after the OFF peak, or model ON and OFF drive.
* **Next:** a theory-driven fourth model test (recurrent/shunting normalization vs the fixed-delay model). It
  needs its own pre-registration after the sweep. The condition set earlier ("no more model tests without a new
  reason") is met by this theoretical reason.
* **Fourth model test: binocular ORGaNICs (registered internally 2026-10-30, before implementing or fitting; Blake
  asked for a model of the slow low-contrast rise and the high-contrast collapse).**
  * **Why:** Heeger & Zemlianova 2020 (Figs 3–4, verified in the PDF) show the effective time constant long at low
    contrast and short at high, giving slow, gradual rises and a phase advance (shorter latency) with contrast: the
    feature the fixed-delay model lacks (Fig 6). An exploratory hand-set simulation (scratch, not reported) of two
    eye-preferring populations with pooled modulators gave rise times 41.5 → 11.8 ms (11 → 90%), a binocular
    latency advantage of ~4.5 ms, and add-then-normalize with a faster, deeper collapse at higher contrast; it
    normalized too little at 11%.
  * **Model 'org' (10 parameters):** populations D (recorded; drive z_D = (c_D + a·c_N)·u(t)) and N (drive
    z_N = (c_N + a·c_D)·u(t)); u = gamma-filtered step (δ, τ1) as in 'full'. Rawat form:
    τ_y dy_i/dt = −√a_i·y_i + b·z_i;  τ_a da_i/dt = −a_i + σ² + y_i²·a_i + w·y_j²·a_j.  R = g·max(y_D, 0)^n.
    Parameters g, b, σ, n, δ, τ1, τ_y, τ_a, a, w. Exponential (exact for piecewise-constant coefficients)
    integration at 1 ms.
  * **Data and fitting:** identical to 'full' in fig5_model.py: channel-averaged, dominant-oriented 6 × 6 time
    courses per penetration, 0–219 ms, 35 non-blank cells, odd/even cross-validation, cleaned data; same starts
    policy (base + random), least squares.
  * **Primary (architecture):** per penetration r = log(CV_full / CV_org); one-sided Wilcoxon over penetrations,
    org better if p < .05 and the median r > 0 in both monkeys.
  * **Secondary (dynamics; the reason for the test):** from the all-presentation fits, 10–90% rise time of the
    dominant eye alone at 11% and 90% (same metric as dynamics_diagnostics.py, on the model prediction).
    Data rise is measured identically on the same channel-averaged time course the model was fitted to (amended
    before fitting: per-channel medians are not comparable with a channel-average fit). 'org' passes if the median
    over penetrations of (model − data) rise at 11% lies within ±10 ms and the median ratio of model rise
    11%/90% is at least 1.5; the fixed-delay fits are scored the same way for comparison.
  * **Tertiary (collapse):** late linearity index (windows starting 50, 75, 100 ms after onset, channel-average
    readout as in fig5_model_L) model minus data at 11, 22.5, 45%; reported with 90% bootstrap CIs and the ±0.25
    equivalence margin, as for the other models.
  * **Reported regardless:** parameters (τ_y, τ_a, σ, w vs a), binocular-vs-monocular latency in the model, and
    the offset (no OFF drive is modelled, so no offset transient is expected).
  * **Result (2026-10-30; `scripts/fig5_org_model.py`, [qc/fig5_org.csv](qc/fig5_org.csv); 25 penetrations).**

    | Criterion | Result | Verdict |
    |---|---|---|
    | Primary: held-out error, org vs full | org better in 10/25; median log ratio −0.026; p = .75 (B +0.031, J −0.040) | **fails**: no better overall |
    | Secondary: rise at 11%, model − data; rise ratio 11/90% | org −5 ms, ratio 2.06 (data 33 ms, ratio 2.13); full −13 ms, ratio 1.88 | **org passes**, full fails |
    | Tertiary: late linearity index within ±0.25 | neither model at any contrast (both ~+0.4 at 11%) | fails, both |

    * **Parameters:** τ_y 2.1 ms, τ_a 82 ms (IQR 59–96), σ 0.22, n 1.02; pool weight w 1.40 > excitatory a 0.83
      (the same w > a as the delayed model).
    * **Binocular latency advantage, a prediction rather than a fitted parameter:** binocular − dominant half-max
      latency at 45% is −2 ms in the model and −1 ms in the data (−1 to −2 ms across 11–45%, earlier in 69–81% of
      penetrations). *(Corrected in the fifth test: the fixed-delay model gives the same −2 ms, so this is not
      specific to ORGaNICs.)*
  * **Check after the result (descriptive, not a rescue): the late 11% index is mostly noise.** Adding the data's
    own noise sample ((odd − even)/2) to each model's prediction moves the late 11% index by roughly ±0.4–0.7
    depending only on the sign of the noise (full: 1.46 vs 0.71; org: 1.24 vs 0.81), against a model − data gap of
    ~0.4. At 22.5 and 45% the noise moves it by ±0.1–0.2, and both models lie within that range of the data. The
    late linearity index at low contrast cannot discriminate models; the significant early-to-late drop at 11% stands.
  * **Reading:** recurrent (ORGaNICs) normalization reproduces the slow rise at low contrast better than the
    fixed-delay model. *(Corrected in the fifth test: the latency advance with contrast and the earlier binocular
    onset, listed here originally, are produced by the fixed-delay model too.)* It fits the full time courses no
    better overall. Candidate next model (needs its own registration):
    delayed normalization whose integration time shortens with the pooled input (conductance-based, Carandini,
    Heeger & Movshon 1997), which nests 'full'.
* **Fifth model test: conductance-based delayed normalization, reconciling 'full' and 'org' (registered internally
  2026-10-30, before implementing or fitting; Blake asked for a reconciliation of the two models as the paper's
  ending rather than stopping at the fourth test).**
  * **Why:** the ORGaNICs response equation, τ_y dy/dt = −√a·y + b·z, is a membrane equation: a drive integrated
    against a conductance (√a) that both divides the steady state and shortens the effective time constant. The
    delayed-normalization model 'full' is the same structure with instantaneous integration. Giving 'full' a finite
    integration time that shortens with the normalization pool (shunting normalization: Carandini & Heeger 1994;
    Carandini, Heeger & Movshon 1997) should give ORGaNICs' timing while keeping the fixed-delay model's fit.
    Distinct from the rejected 'tau_c' model, which made the *pool's* delay contrast-dependent; here the pool delay
    stays fixed and the *response* integrates against the pool.
  * **Model 'cond' (9 parameters: the 8 of 'full' + τ_m):** E = (c_D + a·c_N)·u(t), P = (c_D + w·c_N)·v(t) as in
    'full'; τ_m·σⁿ dR/dt = g·Eⁿ − (σⁿ + Pⁿ)·R. R relaxes toward the 'full' prediction g·Eⁿ/(σⁿ + Pⁿ) with
    effective time constant τ_m / (1 + (P/σ)ⁿ): long at low contrast, short once the pool is strong. τ_m bounds
    0.05–60 ms; at 0.05 ms the step factor is e^−20, so 'cond' reduces exactly to 'full'. Exact exponential
    integration at 1 ms (coefficients held over each step).
  * **Data and fitting:** identical to the fourth test (channel-averaged, dominant-oriented 6 × 6 time courses,
    0–219 ms, odd/even cross-validation, cleaned data, 4 starts), plus one extra start at the 'full' fit of the
    same split with τ_m at its lower bound, so 'cond' is never worse than 'full' in-sample.
  * **Primary (does it fit better than 'full'):** per penetration r = log(CV_full / CV_cond); one-sided Wilcoxon,
    'cond' better if p < .05 and the median r > 0 in both monkeys.
  * **Secondary (does it keep ORGaNICs' timing):** exactly as the fourth test: median (model − data) 10–90% rise
    at 11% within ±10 ms and median model rise ratio 11%/90% ≥ 1.5.
  * **Tertiary (collapse):** late linearity index (windows starting 50, 75, 100 ms) model − data at 22.5% and 45%
    within ±0.25 (90% bootstrap CI of the median). The 11% value is reported but not tested, because the fourth
    test's noise check (above, done before this registration) showed it moves ±0.4–0.7 with the data's own noise.
  * **Reported regardless:** 'cond' vs 'org' held-out error (descriptive; the 'org' fits from the fourth test on
    identical data and splits); τ_m and the effective time constant at 11% and 90%; half-max latency at 11% minus
    90% (the latency shift with contrast) for data and all three models; binocular − dominant latency at 45%.
  * **If the primary fails:** the ending is descriptive (which model reproduces which signature), not a winning
    model; no further variants on these data without a new theoretical reason.
  * **Result (2026-10-30; `scripts/fig5_cond_model.py`, [qc/fig5_cond.csv](qc/fig5_cond.csv); 25 penetrations).**

    | Criterion | Result | Verdict |
    |---|---|---|
    | Primary: held-out error, cond vs full | cond better in 21/25 (2 ties at τ_m floor); median log ratio 0.015 (1.5% less error; 3.3% where τ_m > 1 ms); p = .0002; J median 0.022 (15/18), B median 0.0002 (6/7) | **passes**; the gain is small and mostly in J |
    | Secondary: rise at 11%, model − data; ratio 11/90% | cond −10 ms, ratio 1.87 (full −13, 1.88; org −5, 2.06; data 33 ms, 2.13) | **passes at the margin** (−10 is the bound) |
    | Tertiary: late index at 22.5, 45% within ±0.25 | medians −0.21, +0.10; 90% CIs −0.33 to +0.29, −0.33 to +0.34 | **fails** (imprecise, not biased) |

    * **Descriptive, cond vs org:** cond has less held-out error in 18/25 (median log ratio 0.047, p = .002); J 16/18
      favour cond, B 5/7 favour org.
    * **Parameters:** τ_m 13 ms (IQR 0.4–32; at the floor, i.e. identical to 'full', in 5/25); pool delay τ2 38 ms;
      w > a in 24/25 (1.07 vs 0.88, p < 10⁻⁷), as in 'full'. Effective time constant at 11% vs 90%: 13 vs 8 ms at
      onset (15 ms after δ), 6.7 vs 1.2 ms once the pool is up (120 ms).
    * **Rise time (10–90%) by contrast, medians in ms, data / full / org / cond:** 5.5% 41/22/31/27; 11% 33/20/28/22;
      22.5% 28/16/24/19; 45% 20/14/18/16; 90% 18/11/14/12. Every model rises too fast at every contrast; ORGaNICs is
      closest throughout; cond slows the 11% rise by 3 ms over 'full' (paired: slower in 19/25, never faster).
  * **Correction to the fourth test's reading (found in this run):** the latency advance with contrast and the
    earlier binocular onset are *not* ORGaNICs-specific. Half-max latency 11% − 90% is 5 ms in the data and 6, 7 and
    6 ms in full, org and cond; binocular − dominant latency at 45% is −1 ms in the data and −2 ms in all three
    models. Every normalization model fitted here produces both, because a stronger drive reaches half its peak
    sooner. They do not discriminate models and should not be presented as a recurrent prediction. The rise time
    is the timing signature that discriminates.
  * **Reading:** delayed normalization acting as a conductance (the pool both divides the response and shortens its
    integration time) fits the time courses better than the fixed-delay model and than ORGaNICs, keeps the delayed
    pool weighted toward the other eye, and shares ORGaNICs' core mechanism (the normalization signal sets the
    response's speed). It moves the slow low-contrast rise only part of the way; ORGaNICs, whose normalization
    signal grows with the response itself, gets closer on rise but fits worse. Discussion candidate (no further
    model test on these data): a pool with a recurrent, response-driven component would combine the two.
    *(The rise-time part of this reading is superseded by the correction below.)*
* **Correction (2026-10-31, found while building the merged model figure): the slow low-contrast rise is mostly
  measurement noise.** Peak-normalized traces at 11% averaged only ~0.7 of their own peak in the data but ~0.95 in
  the models, the signature of a noise-inflated peak. The 10–90% rise is measured against that peak, so it was
  checked like the late linearity index (`scripts/rise_noise_check.py`, [qc/rise_noise.csv](qc/rise_noise.csv)).
  * **Per-penetration medians (the metric of the fourth and fifth tests' secondary criteria):** adding the data's
    own noise sample ((odd − even)/2) to a model's prediction lengthens its rise at 11% from 20 to 33–44 ms (full),
    22 to 36–46 (cond), 28 to 45–50 (org); data 33 ms. At 5.5% the noisy models exceed the data (52–61 vs 41 ms).
    At 22.5–90% the noise effect is 0–7 ms. The low-contrast rise is noise-dominated, so the secondary criteria of
    both tests (rise at 11% within ±10 ms) could not discriminate models: the "org passes, full fails" and
    "cond passes at the margin" verdicts are uninformative, not evidence for either model.
  * **Population-mean responses (noise ~1/5; descriptive, chosen after the per-penetration metric was shown to be
    noise-dominated):** data rise 22, 21, 21, 19, 17 ms at 5.5–90%, nearly flat; full 27/23/20/17/13, cond
    32/25/21/18/15, org 33/29/23/19/16. Model − data (95% penetration-bootstrap CI): every model inside the data's CI
    at 5.5–22.5%; at 45% full −4 [−6, −1], cond −2 [−4, 0], org −1 [−3, +2]; at 90% full −3 [−5, −1], cond −2
    [−4, −1], org −2 [−4, 0].
  * **What is withdrawn:** that the data rise slowly at low contrast and that the fixed-delay model "rises ~2.3×
    too fast at low contrast" (Fig 6 of drafts V1–V2; the per-channel 46 → 20 ms figures of the dynamics
    diagnostics are noisier still); that ORGaNICs' advantage is reproducing a slow low-contrast rise.
  * **What stands:** the fixed-delay model rises 3–4 ms too fast at high contrast; integration (conductance or
    recurrent) brings it within ~2 ms; the conductance model's held-out advantage (a whole-time-course measure,
    unaffected by this bias) stands.
  * **To check next (same bias, not yet tested):** every threshold-crossing time measured per penetration or per
    channel group at low contrast, i.e. the normalization half-transition times (Fig 3d) and the laminar
    normalization times ("first window with linearity index < 0.5", Fig 4), should get the same noise check.

## V3 analysis plan (registered internally 2026-10-31, before running; approved by Blake)

Motivated by the two reviews of draft V2 ([manuscript/review_v2.md](../manuscript/review_v2.md)). **Blake's
guidance:** keep a clear, clean main figure set that carries the central narrative (the eyes add, then normalize
over time; contrast sets how far; the model). These analyses go to the supplement and are cited in one sentence of
text unless one changes the story. Fixation (vergence) disparity is argued in the Discussion, not measured (only
211008 recorded both eyes). Each analysis reports its outcome plainly, pass or fail; no post-hoc rescue.

Common definitions: "channel-average readout" = `fig5_model.population` + `windowed_L` (the model-test readout);
"per-channel readout" = the Fig 2d pipeline (`equivalent_contrast.py`, penetration medians over channels); early =
0–20 ms after onset, late = 50–200 ms (per-channel) or windows starting 50, 75, 100 ms (channel-average); noise
sample = (odd − even)/2 of the same data.

**A. Is the measure valid?**
* **A1 Fixed-rule null.** Per penetration, the all-presentation 'full' fit with w set to a (the pool weighted like
  the drive): every binocular response then equals the monocular response at (1 + a)c at every moment, with the
  fitted transient shapes. Read out with the channel-average readout, noise-free and with the noise sample added
  and subtracted. *Valid if* (i) noise-free, the median index over penetrations stays within 1 ± 0.25 in every
  window 0–150 ms at 11, 22.5 and 45%; (ii) with noise, the null's early-minus-late drop at 45% is smaller than the
  data's (paired one-sided Wilcoxon over penetrations, p < .05, using each sign of the noise); (iii) the null's
  20–50 ms slope against log contrast is negative in fewer than 20 of 25 penetrations (data: 26/26 per channel).
* **A2 Cross-validated dominance.** Per channel, dominant eye and eye weight a from odd presentations; equivalent
  contrast and index from even presentations; and the reverse; averaged. Compared with the same-data pipeline run
  on the same halves. *Effect stands if* the early-to-late drop at 45% remains significant (Wilcoxon over
  penetrations, p < .05) and the cross-validated late index at 45% is within 0.25 of the same-data value. The shift
  (cross-validated − same-data) is the size of the selection bias.
* **A3 Ocular-dominance strata.** Channels split into terciles of |ocular dominance index| (from odd presentations;
  index from even). Early and late index at 22.5 and 45% per tercile. Descriptive: does the sequence hold in
  strongly monocular and in balanced channels?
* **A4 Isolated units.** Early and late index at 22.5 and 45% for driven units (lenient set), units as the sample,
  with the net-suppression observation (Discussion of V2) set beside it. Descriptive, small n.

**B. Is the delay specifically interocular?**
* **B1 Model index time courses.** Static, w = a and delayed ('full') fits read out with the channel-average
  readout in every window, against the data (median over penetrations). Descriptive; candidate panel for the model
  figure.
* **B2 Static vs delayed on the binocular interaction.** Held-out (odd/even) comparison of static and delayed
  models fitted and scored with the interaction loss (`dynorm.interaction_residuals`: monocular time courses and the
  binocular interaction weighted equally). *Delay needed for the interaction if* delayed is better with one-sided
  Wilcoxon p < .05 and median log-ratio > 0 in both monkeys.
* **B3 Separate pool delays.** Model 'sep' = 'full' with the other eye's pool term delayed by its own time constant:
  P(t) = c_D·v(t; τ2) + w·c_N·v(t; τ2N). Nests 'full' (τ2N = τ2; extra start at the 'full' fit). *Interocular
  normalization is slower if* median (τ2N − τ2) > 0, one-sided Wilcoxon p < .05, in both monkeys; held-out error
  vs 'full' reported. If not, the paper says monocular and interocular normalization share one delay.
* **B4 Parameter recovery for w > a.** Synthetic data from each penetration's w = a fit plus noise samples
  (independent random sign flips per condition for the odd and even halves); fit 'full'; distribution of fitted
  w − a. *w > a is not a fitting artefact if* the data's median w − a exceeds the 95th percentile of the synthetic
  medians (100 synthetic sets).

**C. Iso-response contours (6 × 6 grid).** Per penetration (channel-average) and window (early, 20–50 ms, late):
for target levels equal to the dominant eye's response at c* = 11, 22.5, 45%, find for each other-eye contrast the
dominant-eye contrast giving that response (monotone, log-contrast interpolation, as for equivalent contrast).
Summarize each contour by the Minkowski exponent m in (c_D/c*)^m + (a·c_N/c*)^m = 1 (m = 1 linear summation,
large m = dominant eye alone; m bounded 0.5–20), and by the fraction of contour points beyond the dominant-eye-alone
line (suppression). *Contours change from linear to max-like if* late m > early m (Wilcoxon over penetrations,
p < .05) with early median m within 0.7–1.5. Main-figure candidate (Blake decides); population contours with
bootstrap CIs.

**D. Timing measures and noise.** Half-transition times (Fig 3d) and laminar normalization times ("first window
with index < 0.5", Fig 4d): (i) split-half replication (odd and even presentations each give the same sign of the
group differences); (ii) noise sensitivity (each measure on all presentations vs on halves, per group); (iii)
SNR-matched comparison (subsample channels or add the noise sample so groups have equal signal-to-noise). *The
laminar order stands if* it survives (iii) with p < .05 in the same direction and replicates in both halves;
otherwise the paper states only that every depth shows the sequence.

**E. Literature.** Binocular vs monocular response latency (human VEP, cat, macaque); contrast–latency; duration
dependence of binocular summation and contrast matching in psychophysics; the citations flagged by the cold read
(Smith et al. 1997; Macknik & Martinez-Conde 2004; Legge & Rubin 1981; Baker, Meese & Georgeson 2007; Menz &
Freeman 2003).

### V3 analysis results

* **D, timing measures and noise (2026-10-31; `scripts/v3_timing_noise.py`, [qc/v3_D.csv](qc/v3_D.csv)).**
  **Laminar-order criterion fails.** Fig 4d's numbers reproduce exactly. Halves recomputed with
  `equivalent_contrast.penetration` on odd and even presentations (dominance, a and onset re-estimated per half);
  SNR matched by adding λ × the noise sample (both signs) to higher-SNR groups (group SNR: above 2.17, sink 2.77,
  below 1.51).
  * Normalization on the stimulus clock (group − sink): all +6 / +7 ms (p = .016 / .003); odd half +3 / +3.5
    (p = .19 / .069); even half +4 / +1 (p = .75 / .31); SNR-matched above vs sink p = .33 (+noise) / .055
    (−noise), below vs sink p = .026 / .004. Signs replicate; significance does not.
  * Time after each group's own onset: median difference 0 for both groups (ties in 16/21 and 11/23
    penetrations). The stimulus-clock order is the onset order plus 15 ms window quantization; noise shifts the
    timings by only 0 to −1 ms.
  * Onset order: sink channels first by 2–4 ms; same sign in both halves (odd p = .014 / .032, even p = .07 / .13).
  * Fig 3d half-transition times: 35 / 30 / 15 ms (Friedman p = .19) do not replicate (halves 25/35/25 and
    25/22.5/25) and do not survive SNR matching (contrast SNR 0.73 / 1.20 / 2.80; p = .66 / .19).
  * **For the paper:** keep "channels at the earliest sink respond first" and "every depth normalizes a similar
    time after its own onset"; drop "normalization reaches the input layer first" (abstract, Fig 4 title,
    Discussion). Fig 3d and Fig 4d go to the supplement. Deviations: noise added to group-median index curves, not
    raw responses; sink depth from all presentations in both halves.
* **C, iso-response contours (2026-10-31; `scripts/v3_iso_contours.py`, [qc/v3_C.csv](qc/v3_C.csv), diagnostic
  `<out_root>/figures/diagnostics/v3_iso_contours.png`; 25 penetrations). Criterion passes at every c*.**
  * Median Minkowski exponent m early / 20–50 / late: 11% 1.15 / 1.01 / 1.38 (late > early in 60%, p = .009);
    22.5% 1.00 / 1.13 / 2.18 (84%, p < .001); 45% 0.90 / 1.33 / 2.39 (76%, p < .001); early medians within
    0.7–1.5 at every c*; both monkeys' median changes positive.
  * Population contours (Fig 2b readout, penetration bootstrap): early m 1.12 / 1.00 / 0.89 with CIs including 1;
    late 1.41 [1.08, 2.12], 2.30 [1.53, 20], 1.81 [0.96, 20]. Late contours move partway toward the
    dominant-eye-alone line, not to an L shape, with wide CIs.
  * Suppression (contour points needing more dominant-eye contrast than c*): none early; late at low other-eye
    contrast (up to 1.36× at 22.5%), per-penetration fractions 0 / 0.4 / 0.4.
  * Undefined points: "below axis" (other eye alone reaches the target) 151 / 127 / 98 of 375, almost all where
    a·c_N ≥ c*; censored (target beyond the dominant eye at 90%) 0 / 31 / 52, set to 90% as lower bounds.
  * Deviations: `fig5_model.T` extended to 0–280 ms inside the script so the late window is complete; each c*
    tested separately; rows with c_N > 0 inverted from their own value at c_D = 0; population contours from
    population-average responses. Caveat: the Minkowski form cannot represent suppression, so suppressed and
    censored points push m to its bound; m mixes "closer to the dominant eye alone" with "suppressive".
  * **For the paper:** a strong main-panel candidate beside Fig 2a–b (within-window, so latency-proof; uses the
    dichoptic cells; the language of the psychophysics). Show 45% (or 22.5 and 45%), early vs late; give m in the
    text.
* **A2–A4, dominance, ocular-dominance strata, units (2026-10-31; `scripts/v3_dominance_cv.py`,
  `scripts/v3_units_index.py`, [qc/v3_A2.csv](qc/v3_A2.csv), [v3_A3](qc/v3_A3.csv), [v3_A4](qc/v3_A4.csv)).**
  The all-presentation variant reproduces Fig 2d exactly.
  * **A2 passes.** Cross-validated (dominance, a and onset from one half, index from the other, both directions
    averaged; 25 penetrations, 211008 drops out for too few presentations per half): 45% early 0.99 → late −0.32
    (p < .001); same-data on halves 0.90 → −0.25. Selection-bias shift (late, paired) +0.015 [−0.026, 0.040] at
    45%, +0.02 to +0.023 elsewhere (CIs include 0). Dominance flips between halves in 26% of channel-halves (43% in
    the most balanced tercile) yet the index barely moves.
  * **A3:** the sequence holds in every |ODI| tercile at 45% (balanced 0.99 → 0.40; middle 1.03 → −0.10;
    monocular 0.78 → −0.56; all p < .001) and deepens with monocularity; at 22.5% the balanced tercile does not
    collapse significantly (1.37 → 1.15, p = .08). Small eye weights in monocular channels magnify the index.
  * **A4 (descriptive):** lenient units (55–61 in 19 penetrations): 22.5% 1.50 → 0.19 (p < .001); 45% 0.19 → −0.78
    (p = .001), with very noisy early values (IQR −1.2 to 1.2). Strict and cross-validated variants weaker
    (p = .051–.08 at 45%). Net binocular modulation at 45%: transient −0.05, sustained −0.01, both n.s. (41–44%
    facilitated): **V2's "units mostly net-suppressed" overstates; correct it.**
  * **Noted:** halving the data lowers the late 45% index (−0.25 on halves vs −0.12 on all presentations): noise
    pushes the per-channel late index down (running maximum, censoring). With the channel-average (−0.04) and
    population (+0.1) readouts, "below the dominant eye alone" is partly noise: write "matched".
* **A1, B1, B4, null, model index time courses, parameter recovery (2026-10-31; `scripts/v3_null_model_index.py`,
  `scripts/v3_param_recovery.py`, [qc/v3_A1.csv](qc/v3_A1.csv), [v3_B1](qc/v3_B1.csv), [v3_B4](qc/v3_B4.csv)).**
  * **A1 passes all three criteria: the measure is valid.** Noise-free null (w = a fits): median index 0.96–1.03
    in every window at every contrast (largest deviation 0.035). Early-minus-late drop at 45%: data 1.25, null +
    noise 0.11, null − noise 0.23 (paired one-sided Wilcoxon p = .00017 / .0023). Null 20–50 ms slope negative in
    13/24 and 16/25 (noise-free 7/25) vs data 25/25.
  * **B1:** only the delayed model declines over time (45%: 0.97 → −0.08; 22.5%: 1.00 → 0.57; 11%: 1.00 → 0.93).
    The static model is set within ~10 ms and stays flat (45% ≈ 0.50, 22.5% ≈ 0.82, 11% ≈ 1.02), but **it
    reproduces the contrast dependence** (20–50 ms slope negative in 25/25; delayed 24/25; w = a 3/25). So the
    contrast dependence (Fig 3b–c) is not evidence for delay; the decline over time is. Answers the cold read's
    main model question.
  * **B4 passes:** synthetic w = a data refitted with 'full' (100 sets × 25 penetrations) give median w − a 0.000
    (5–95% of set medians −0.006 to 0.006; w > a in 12/25 per set) vs data 0.191 (24/25). Rules out an optimizer or
    noise artefact, not model misspecification.
  * **Observed (B1 readout, not registered):** at 45% the channel-average index undershoots (−0.61 at windows
    starting 35 ms, −0.28 at 50) and partly recovers (0.43 at 75–150 ms); at 22.5% 0.31 at 35–75 ms, 0.79 at
    100–150. The delayed model declines monotonically and misses it. The per-channel late window (50–200 ms)
    averages over the undershoot and the recovery. Logged in docs/10 as a candidate; needs a noise check before
    any claim.
  * Deviations: early / 20–50 / late = channel-average windows starting 0+10 / 20+35 / 50+75+100 ms; both noise
    signs required to pass; static and w = a refitted here (starts = 4); B4 noise added to all-presentation
    predictions with random sign per condition.
* **B2, B3, interocular delay (2026-10-31; `scripts/v3_interocular_delay.py`, `dynorm` model 'sep',
  [qc/v3_B2.csv](qc/v3_B2.csv), [v3_B3](qc/v3_B3.csv)).**
  * **B2 passes: the binocular interaction itself needs the delay.** Fitted and scored on the interaction loss,
    delayed beats static held out in 23/24 penetrations (B 5/6, J 18/18), median log ratio 0.72 (~52% less
    error; B 0.53, J 0.74), p = 1.2 × 10⁻⁷; in-sample delayed ≤ static in 24/24; τ2 under this loss 24 ms (IQR
    20–39).
  * **B3 passes its criterion on a negligible effect.** τ2N − τ2 = +2.1 ms (positive in 20/25; B +3.8, J +1.6;
    p = 9.4 × 10⁻⁵); τ2 28.7 ms, τ2N 31.0 ms; odd/even sign agreement 76%; no bound hits. Held out, 'sep' does not
    beat 'full' (13/25, p = .33). **For the paper:** monocular and interocular normalization share essentially one
    delay (the interocular term lags by ~2 ms at most); drop "interocular normalization is slower than monocular".
* **Correction found by B2 (applies to every held-out model comparison): 211008's even half is empty** (no
  condition has the 2 presentations `population()` requires), so its even-split fits never touch data and its
  "held-out error" is that of an unfitted model (cv_full 67 278 vs cv_org 1 985). This, not the 217 ms stimulus,
  is the Fig 5b outlier. Rechecked without 211008: cond vs full 20/24, median +0.016, p = .0002 (B +0.004, J
  +0.022); org vs full 9/24, p = .86; cond vs org 18/24, p = .0002; delayed vs static 23/24, p = 1.2 × 10⁻⁷.
  Conclusions unchanged; the paper reports held-out comparisons on 24 penetrations and says why.

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

### Closest prior paper: Cox et al. 2019 (J Vis 19(12):13, [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477/))

It motivated Mitchell et al. 2023. Awake *M. radiata*, laminar arrays (100 µm), SUA/MUA/CSD, 500 ms gratings.

**Design:**
* Conditions: monocular, binocular congruent, and binocular incongruent (orthogonal).
* **One fixed contrast pair per penetration set** (dominant eye about half the non-dominant: .45/.9, .5/1, .4/1).
* Fixed windows: transient 50–100 ms, sustained 150–250 ms.

**Findings that match ours:**
* Congruent binocular stimulation boosts the transient response (+13.7 spikes/s; mean modulation index 0.21) and
  the boost vanishes in the sustained window. This is our "add at onset, then normalize" (Fig 2).
* The initial granular sink is larger for binocular stimulation, and the congruent sink returns to the monocular
  level by 150–250 ms. This is our "earliest input adds, then currents normalize" (Fig 4).
* Equiocular granular neurons stay facilitated late; our granular sites are the most linear late (ns within
  penetrations).

**Possible difference:** they found ocularity correlates with *sustained* but not transient modulation; we see
the OD dependence in MUA in both fixed windows. Their units vs our MUA, and our 40–140 ms window already contains
the collapse.

**What they did not do and we add:**
* contrast as a parametric variable (full 6×6 matrix);
* a latency-proof, onset-aligned, time-resolved measure of the combination rule (k(t), L(t)) with eye weights;
* the normalization framing with model tests (linear summation vs normalization vs two-stage; L < 0 at high
  contrast);
* corrected timing, and depth across 26 penetrations;
* a dynamic normalization model.

They did not invoke divisive normalization, and called the early summation both "additive" (abstract) and
"subadditive" (discussion). Our k(t) settles that: linear at onset for the eye weights.

**Their incongruent (dichoptic cross-orientation) results are outside our design.**

## Novelty check (resolved by Blake, 2026-11-01)

* **Mitchell et al. 2022 (iScience) is a separate dataset**: mostly the same contrast in the two eyes, at different
  contrast levels. The earlier note here, "the dataset underlies Mitchell et al. 2022", was wrong. Its result
  (summation 50–100 ms, contrast-dependent suppression by 150–250 ms) is prior evidence for the sequence, not
  overlap: cite it as the starting point.
* **The 2024 dissertation's third chapter is markedly different** from this paper. No overlap concerns with either.
* **Framing for V3:** the first millisecond-resolution account of the transition, with the full dichoptic 6 × 6
  design: linear summation in contrast units for the first ~20 ms; normalization taking hold 20–35 ms after
  onset; iso-response contours; a model showing the delay is needed for the binocular interaction itself. Cite
  Mitchell et al. 2022, Cox et al. 2019, Dougherty et al. 2019 and Tanabe & Cumming 2014 for the coarser sequence
  and the delays already reported.
* Still to cite or check: the relation of these recordings to Mitchell et al. 2023 (Curr Biol) was not discussed;
  Carlson et al. 2025 preprint (interlaminar coherence, same sessions; see the Brock coherence paper notes).

## Order of work

1. Fig 2 core: k(t) with eye-weight correction; normalization latency. This decides the paper.
2. Fig 3: depth profiles of k(t), latency, w.
3. Fig 4: per-condition CSD, earliest sink.
4. Fig 5: units.
5. Fig 1 assembly; Fig 6 model if 1–3 hold.
6. Controls and supplementary.

**Figures as built (2026-10-29 to 2026-10-31; superseded by the figure set for V3 below, whose scripts replaced
these on 2026-11-01; the figure packet script was retired, as the Word manuscript carries the legends).** Each figure
script only plots saved analysis outputs and prints the statistics it shows.

| Figure (as built) | Figure script | Analysis scripts it reads |
|---|---|---|
| 1 Setup | `fig1_setup.py` | `figure_data.py`, flash CSD + alignment v2, `run_crf.py`, `run_dichoptic.py` |
| 2 Add, then normalize | `fig2_add_then_normalize.py` | `figure_data.py`, `equivalent_contrast.py` |
| 3 Contrast sets the collapse | `fig3_contrast.py` | `figure_data.py`, `normalization_timing.py` (incl. `--tag cross`) |
| 4 Depth and earliest sink | `fig4_depth_input.py` | `sink_data.py`, `equivalent_contrast.py`, `dichoptic_csd.py`, `laminar_timing.py` |
| 5 Normalization as a delayed conductance (2026-10-31: merges old Figs 5 and 6) | `fig5_models.py` | `fig5_model.py`, `fig5_org_model.py`, `fig5_cond_model.py`, `rise_noise_check.py` |
| S Model variants (old Fig 5d–e: two-stage, Ding–Sperling) | `figS_model_variants.py` | `fig5_channel_test.py`, `fig5_interaction_test.py` |
| S Offset transient (old Fig 6c) | `figS_offset.py` | `figure_data.py`, `fig5_model.py` |

Units (old Fig 5: units net-suppressed while MUA is facilitated) go to the supplement. The old Fig 6b (per-channel
rise times vs the fixed-delay fit) is dropped: its low-contrast rise times are noise-dominated (correction of
2026-10-31 above); the new Fig 5e shows rise times of the population-mean response instead.

### Figure set for V3 (approved by Blake 2026-11-01)

Four main figures, each with one claim; everything else is supplementary and cited in the text in a sentence.
Rules for every figure:
* no session IDs;
* "condition", never "cell";
* no registration language;
* the other eye's weight in the measure is the "eye weight", and in the model the "drive weight", beside the "pool
  weight";
* held-out model comparisons use 24 penetrations (211008 has no usable even half), with the reason in the legend.

| Figure | Claim | Panels | Script |
|---|---|---|---|
| 1 | Design and an example | (a) 6 × 6 dichoptic contrast design and trial sequence; (b) example penetration, onset-aligned, in the Fig 2a format; (c) monocular contrast response, dominant vs other eye; (d) ocular dominance index, MUA and units | `fig1_setup.py` |
| 2 | The eyes add, then normalize | (a) population responses at 45%, onset-aligned; (b) equivalent-contrast readout; (c) iso-response contour for the 45% target, early vs late, in percent contrast (population, penetration-bootstrap CIs; per-penetration exponent test in the note; 22.5% and the fitted contours in S3); (d) early vs late linearity index per penetration | `fig2_add_then_normalize.py` |
| 3 | Contrast sets how far it collapses, at every depth and in the input current | (a) transition-window index per penetration vs contrast; (b) slope on log contrast per penetration and in the held-out stimulus-position runs; (c) early and late index by depth; (d) population CSD aligned to the earliest sink; (e) sink current at 45% in the Fig 2a format | `fig3_contrast_depth.py` |
| 4 | Normalization is delayed, weights the other eye more, and acts like a conductance | (a) model schematic and colour key; (b) data vs conductance fit, Fig 2a format; (c) linearity index over time: data, static, delayed; (d) held-out interaction error, delayed vs static; (e) pool weight vs drive weight of the other eye; (f) held-out comparison of the three dynamic models | `fig4_model.py` |
| S1 | Apparatus and laminar alignment | (a) stereoscope; (b) example flash CSD with alignment | `figS1_methods.py` |
| S2 | The linearity index is valid | (a) fixed-rule null through the pipeline (noise-free and with the data's noise) vs data; (b) cross-validated dominance; (c) ocular-dominance terciles; (d) isolated units | `figS2_validity.py` |
| S3 | Iso-response contours at every target contrast | population contours at 11, 22.5, 45% early / 20–50 ms / late; per-penetration exponents | `figS3_iso_contours.py` |
| S4 | Timing measures and noise | (a) raw binocular gain (old 3a); (b) half-transition times with halves (old 3d); (c) laminar onset and normalization times, stimulus clock and after own onset (old 4d); (d) binocular gain of the sink current and MUA (old 4c) | `figS4_timing.py` |
| S5 | Model checks | (a) two-stage gain; (b) capture test; (c) late index, model − data (old 5g); (d) parameter recovery of w − a; (e) separate interocular pool delay | `figS5_model_checks.py` |
| S6 | Response dynamics | (a) dominant-eye time courses with fits (old 5d); (b) rise time of the population mean (old 5e); (c) offset transient | `figS6_dynamics.py` |
| Table S1 | Main statistics per monkey | | `tableS1_monkeys.py` |

**Built 2026-11-01** (all scripts in `scripts/figures/`; legends in `manuscript/draft_v3.docx`). What the build
changed in the claims:
* **Fig 3e, sink current in the Fig 2a format:** the input current adds at onset (both eyes > dominant eye alone at
  45% in the first 20 ms after sink onset, 22/26, p = .001; close to the dominant eye at 90%), but the later
  combination cannot be tested, because the dominant eye's sink currents at 45% and 90% no longer differ after 50 ms
  (p = .45). Claim narrowed to "the earliest input current sums the two eyes at onset". Display scaling (each
  penetration by its mean monocular sink current, 0–150 ms) was chosen after seeing that peak scaling let two
  penetrations dominate the mean; the per-penetration sign counts do not depend on it.
* **Fig 3c, depth:** at 45% the index falls at every depth; at 11 + 22.5% the drop does not depend on depth (mixed
  model p = .60); the inverted U of the late index (peak ~1.7–1.9 mm, p = .0002) is monkey J only (B p = .34).
* **Fig 4b–c, model:** the conductance fit reproduces the Fig 2a sequence but undershoots transient peaks (~8%);
  the fixed-delay fit misses the decline at 11% (stays ~0.93; data 0.48 [0.18, 0.67]) and the 45% undershoot.
  Stated in the text.
* **Monkey split (Table S1):** monkey B shows no early-to-late decline at 11% (4/8) and its 45% test is p = .098
  (7/8); its iso-contour changes are underpowered. The conductance gain is essentially monkey J's.

**Polished for readability 2026-11-02** (Blake: "Figure 2c reads as quite messy … not sure if this is the best way
to represent the data"). Display only: every script's printed statistics are unchanged except counts of points
drawn at an axis edge.
* **Fig 2c redrawn.** One target (45%) instead of two; axes in percent contrast ("dominant-eye contrast needed (%)"
  against "other-eye contrast × eye weight (%)"), not fractions of the target; "dominant eye alone" is a horizontal
  line at the target, not an L-shaped rule; no fitted Minkowski contours (they fit the late points poorly, exponent
  CI 0.96 to ≥ 20); points joined per window and labelled directly; the statistic is one note (19/25, p < 0.001).
  The 22.5% target and the fitted contours stay in S3, now drawn in the same percent coordinates.
* **Fig 2b:** each window on its own response axis, so the late readout is visible. **Fig 2d:** condition headers
  above each pair, windows as tick labels ("0–20", "50–200"), off-axis counts moved to the legend (S2 b–d likewise).
* **Fig 3:** (a) values beyond the axis drawn as open triangles (8), window in the title; (b) one tick label per
  group; (c) key moved off the reference line; (d) CSD drawn as a bilinear image (the Gouraud mesh drew saw-tooth
  seams between channel rows); (e) shorter note, lighter bands.
* **Fig 4:** (a) one key describing every model the same way; (b) key reads "solid: data / dashed: conductance fit"
  (the 90% trace stays, because the text's comparison with 90% cites this panel); (c) a key for data, static and
  fixed-delay replaces model names stacked beside the reference labels.
* **Fig 1c:** analysis window as the panel title, not in the axis label; (d) "fraction" (MUA and units).
* **S5:** monkey key inside panel e, empty third slot removed, panel c key and note moved off the error bars.

Carry into V4 (draft_v3.docx was open in Word, so the text was not touched):
* Results, contour paragraph: "the contour is instead a vertical line at the dominant eye's contrast" → a
  horizontal line at the target contrast (true of the V3 figure too); cite Fig 2c for 45% and Supplementary Fig S3
  for 11 and 22.5%.
* Legends: Fig 2 (b separate response axes; c rewritten for one target in percent contrast, no fitted curves; d
  off-axis counts 2, 1, 3, 3 at 5.5–45%); Fig 3a (8 values drawn at the axis edge); Fig 4b (solid data, dashed fit);
  S2 b–d (off-axis counts: b 5, 2, 4; c 0, 2, 7; d 46, 32); S3 (percent axes; horizontal "dominant eye alone").

Model colours (no model uses orange, the other eye's colour): data ink, static dark grey dashed, fixed-delay aqua,
conductance violet, recurrent magenta (aqua and magenta below 3:1 on white and close for deutans, so always
direct-labelled). Inside a Fig 2a-format panel a model is drawn as a dashed line in the condition's colour.

**Manuscript (from 2026-10-30):** the text lives in `manuscript/draft_vN.docx` (Word; Blake comments, Claude revises with
tracked changes). Blake's rule: pre-registration is internal discipline only. The paper reports every analysis on
the same (artifact-cleaned) data without registration language; the as-run registered outcomes stay in docs/09-10.
