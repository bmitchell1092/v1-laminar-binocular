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

## V4 check: dip and partial recovery (registered internally 2026-11-02, before running)

**Why.** The accessibility redesign of the figures (Blake, 2026-11-02: "very important that we make this paper
accessible") shows the result as *the response the other eye adds* (both eyes − dominant eye alone) against what
linear summation of contrast predicts (dominant eye at c(1 + a) − dominant eye at c). In that view the population
means at 45% look like a dip and a partial recovery, not a step: the other eye adds about what summation predicts in
the first 20 ms, nothing from ~30 to ~110 ms, and about two thirds of the prediction by 150–200 ms (the B1 undershoot
logged above, in response units). **Not blind:** the windows below were chosen after seeing those population means
in the mockups, so this check can confirm or reject the shape but is not an independent discovery. It also serves as
the noise check the B1 undershoot was waiting for.

**F. Measure.** Per penetration, onset-aligned channel-mean responses as `figure_data.py` (driven MUA channels,
dominant eye, scaling, onset and eye weight a from all presentations), for all presentations and for odd and even
presentations. In each window after response onset: added A = mean of R(c, c) − R(c, 0); predicted P = mean of
R_D(c(1 + a)) − R(c, 0), with R_D interpolated in log contrast between the dominant eye's measured contrasts at each
time point. Windows: early 0–20 ms, dip 30–100 ms, recovery 130–200 ms. Contrast 45% is primary; 11 and 22.5% are
reported.

**Tests (45%, Wilcoxon signed-rank over penetrations, one-sided in the stated direction):**
* F1 dip: A − P < 0 in 30–100 ms.
* F2 recovery: A(130–200) − A(30–100) > 0.
* F3 partial: A − P < 0 in 130–200 ms.
* F4 undershoot: A < 0 in 30–100 ms (the other eye's grating lowers the response).

**Criteria.** A test holds if p < .05 on all presentations *and* the median difference has the same sign in the odd
and in the even half. The paper may then say: F1, "adds nothing from ~30 to ~100 ms"; F2 + F3, "partly recovers but
stays below summation"; F4, "the other eye briefly suppresses". If F2 fails, the text keeps one late window and says
nothing about recovery; if F4 fails, "matched" stays (A2 note above). Script: `scripts/v4_dip_recovery.py`, output
[qc/v4_F.csv](qc/v4_F.csv).

**F results (2026-11-02; 26 penetrations; per-penetration window means in
[qc/v4_F_windows.csv](qc/v4_F_windows.csv)). F1 and F2 hold; F3 and F4 fail.** Population means at 45% (added /
predicted, dominant-eye peak = 1): early 0.199 / 0.164, dip −0.004 / 0.037, recovery 0.030 / 0.041.
* **F1 holds:** at 45% the other eye adds less than summation predicts in 30–100 ms (23/26, p = 4 × 10⁻⁵; odd
  p = .002, even p = .0001).
* **F2 holds:** its contribution returns, 130–200 vs 30–100 ms (21/26, p = .0004; odd p = .002, even p = .002).
* **F3 fails:** by 130–200 ms the contribution is not reliably below summation (16/26, median −0.011, p = .087; odd
  p = .20, even p = .051). In response units, at 45%, normalization is clear only from ~30 to ~100 ms.
* **F4 fails:** no undershoot (14/26 below zero, p = .24; halves disagree in sign). The other eye adds nothing in the
  dip; it does not suppress. "Matched" stays.
* **Reported, 11%:** the opposite pattern: below summation in the dip (18/26, p = .017; odd half p = .09) *and*
  late (19/26, p = .001, both halves), with no recovery (13/26, p = .81); the other eye adds about 60% of the prediction from ~50 ms on.
  **22.5%:** predictions are small after 30 ms (≤ 0.02) and nothing reaches p < .05.
* **For the paper:** at 45% write "adds nothing from ~30 to ~100 ms, then partly returns"; do not write that the
  binocular response stays at the dominant eye alone through 200 ms at 45%. The per-channel late index (Fig 2d,
  50–200 ms) averages the dip and the return and weights them through the flat late contrast response, so it reads
  lower than the response-unit view; the text should say both. Blake to decide how the abstract and title phrase
  "normalize later" given that the sustained partial normalization is clearest at 11% and the strong one at 45% is
  transient.

## V4 check G: onset latency by eye, contrast and depth (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): Fig 1 should make the point about monocular vs binocular latencies, and show MUA and
CSD latency as a function of depth, binocular vs monocular, "to highlight the importance of accounting for latency.
That's what this paper adds: 1 ms resolution of binocular summation in V1 and laminar organization of
normalization." Known so far, at integer-ms resolution only: binocular − dominant-eye half-max onset −1 to −2 ms at
11–45% (earlier in 69–81% of penetrations; fifth model test), dominant-eye onset 11% − 90% = 5 ms, and onsets 4 ms
(above) and 2 ms (below) later than at the earliest sink (S4c, p = .02, .012). The new question is whether the
binocular onset advance is the size that summing the two eyes' contrasts predicts, at every depth. **Not blind** to
those integer-ms numbers; the summation comparison (G2), the depth split of the advance (G3) and the other eye's
extra delay (G5) have not been computed.

**G. Measure.** Per driven in-cortex MUA channel (the `figure_data.py` rule), condition means over clean
presentations, blank subtracted, Gaussian σ = 2 ms; dominant eye as in `figure_data.py`; eye weight a per channel
(40–250 ms). For all presentations and for odd and even presentations (dominance, a and inclusion from all).
* **Onset latency** of a condition: the first time in 10–250 ms at which the response reaches half its own peak
  (10–250 ms), linearly interpolated between the samples on either side of the crossing (sub-millisecond). Undefined
  if that peak is below 0.2 × the dominant eye's peak at 90%.
* **Conditions:** dominant eye alone and other eye alone at 11, 22.5, 45, 90%; both eyes at 11, 22.5, 45%.
* **Summation prediction for the binocular onset:** the dominant eye's onset at contrast c(1 + a), interpolated
  linearly in log contrast between its onsets at the measured contrasts (11–90%).
* **Depth groups** relative to the earliest sink, as in S4c: above (> 2 contacts above), at (within 2), below
  (> 2 below).
* Per channel, then the median over channels per penetration (per penetration and group for depth), then over
  penetrations; 95% CIs by bootstrap over penetrations (10,000 resamples).

**Tests.**
* **G1, both eyes start earlier:** Δ = onset(both eyes, c) − onset(dominant eye, c) < 0, one-sided Wilcoxon over
  penetrations, at 11, 22.5 and 45%.
* **G2, as early as summation predicts:** r = onset(both eyes, c) − onset(dominant eye, c(1 + a)). Holds at a
  contrast if the 95% CI of the median r lies within ±1 ms (the resolution the paper claims) and both halves'
  medians lie within ±1 ms. A two-sided Wilcoxon of r against 0 is reported with it.
* **G3, at every depth:** at 45%, G1 and G2 within each depth group. A difference between groups (above − at,
  below − at, two-sided Wilcoxon) may be described only if p < .05 and both halves agree in sign.
* **G4, the laminar onset sequence (re-measured with this estimator; consistency with S4c, not new):** onset of the
  dominant eye at 90% later above and below than at the sink, one-sided Wilcoxon.
* **G5, exploratory (reported, no claim without its own confirmation):** the other eye's extra delay,
  onset(other eye, c) − onset(dominant eye, a·c) at 45 and 90% where a·c ≥ 11%; 0 means the other eye's signal is a
  weaker copy of the dominant eye's with the same timing; > 0 means it arrives later than its strength explains.
  Two-sided Wilcoxon, overall and by depth group; logged in docs/10 if it meets the bar there.
* **Reported:** onset vs contrast for each eye; the binocular advance at 90% (no prediction: c(1 + a) > 90%); the
  same advance in the sink current (`sink.npz`, sink-onset clock).

**Criteria.** G1, G3 (within a group) and G4 hold if p < .05 on all presentations and both halves' medians have the
predicted sign. The paper may then say: G1 + G2 at a contrast, "both eyes started X ms earlier than the dominant eye
alone, as early as summing the two eyes' contrasts predicts (within 1 ms)"; G1 without G2, "earlier, but by
less/more than summation predicts" (direction from the median r); G3, "at every depth" only if it holds in all three
groups. Script: `scripts/v4_latency.py`, output [qc/v4_G.csv](qc/v4_G.csv).

**G results (2026-11-03; 26 penetrations, 570 channels; per-penetration medians in
[qc/v4_G_penetrations.csv](qc/v4_G_penetrations.csv)). G1, G3 (the advance, in every group) and G4 hold; G2 fails
at every contrast; G5 finds no extra delay.**
* **Onsets** relative to the dominant eye at 90% (median of penetration medians; absolute 32.8 ms):

  | | 11% | 22.5% | 45% | 90% |
  |---|---|---|---|---|
  | dominant eye alone | +6.8 | +3.7 | +1.3 | 0 |
  | other eye alone | +7.9 | +4.0 | +1.9 | +0.7 |
  | both eyes | +4.5 | +1.1 | −0.7 | −1.8 |
  | if the contrasts add (dominant eye at c(1 + a)) | +4.9 | +2.1 | +0.1 | — |

* **G1 holds at every contrast:** both eyes start earlier than the dominant eye alone by 1.9 ms at 11% (23/26),
  1.6 at 22.5% (25/26) and 1.5 at 45% (23/26); all p < 10⁻⁴, both halves. Reported, 90%: 1.8 ms (25/26).
* **G2 fails at every contrast as registered** (95% CI of the median not within ±1 ms):
  * 11%: both eyes − summation prediction +0.13 ms (CI −1.10 to +0.86; 13/26 below zero; p = .67);
  * 22.5%: −0.54 ms (CI −1.15 to +0.12; 16/26; p = .14);
  * 45%: −0.79 ms (CI −1.59 to −0.31; 20/26; p = .001; odd half −1.01): both eyes start *earlier* than summation
    predicts (advance 1.5 ms vs predicted 0.8 ms). Both eyes at 45% start before the dominant eye alone at 90% in
    17/26 penetrations.
  * At 11 and 22.5% the advance is close to the prediction and not different from it, but the data do not resolve
    it to within 1 ms. The paper may not say "within 1 ms".
* **G3:** at 45% the advance holds in every depth group (above −2.0 ms, 21/23; at the sink −1.4, 24/25; below
  −1.9, 21/25; p ≤ .002, both halves), so "at every depth" is allowed for the advance. Equivalence fails in every
  group (both eyes − prediction −1.2, −1.0, −0.8 ms). No laminar difference (advance above − at −0.5 ms, p = .29;
  below − at −0.3 ms, p = .29).
* **G4 holds:** the dominant eye's onset at 90% is 4.0 ms later above the sink (16/22, p = .0095) and 1.4 ms later
  below (18/25, p = .015), both halves (S4c: 4 and 2 ms).
* **G5, exploratory:** no extra delay of the other eye: 45% −0.5 ms (p = .008, slightly *earlier* than its strength
  predicts; 8/25 above zero), 90% −0.1 ms (p = .45). The other eye's signal behaves as a weaker copy of the
  dominant eye's. Nothing to log in docs/10.
* **Reported, sink:** the current at the earliest sink also starts earlier with both eyes (−2.3, −1.8, −1.0 ms at
  11, 22.5, 45%; 16/22, 19/25, 16/25; p = .03, .045, .007), as does MUA on the sink channels (−2.2, −1.5, −1.2 ms;
  p ≤ .002).
* **For the paper:** "Raising the contrast in one eye from 11 to 90% shortened the onset by ~7 ms. Adding the other
  eye shortened it by a further 1.5–1.9 ms at every contrast and every depth, and in the current at the earliest
  sink. At 11 and 22.5% the advance was close to what summing the two eyes' contrasts predicts (differences 0.1 and
  0.5 ms, not significant); at 45% it exceeded the prediction by 0.8 ms: two eyes at 45% started sooner than one eye
  at 90%." Descriptive reading (not tested): within one eye the onset gains little above 45% (1.3 ms to 90%), but a
  second eye still gains 1.5–1.8 ms, as if the single-eye limit lay before the cortex and the second eye added drive
  the first could not.

## V4 check H: does onset latency follow drive across the whole 6 × 6 design? (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): "our eye conditions have a rank ordered drive, and it all lines up. If we rank order
the conditions by effective drive, then the latency to response onset also match that effective drive. This means
it's important to look at the response in a manner that compares their peaks, not just windows or bins." **Not
blind** to G (dioptic cells: both eyes start earlier than the dominant eye alone, and at 45% earlier than contrast
summation predicts); the dichoptic cells, the summed-response comparison and the peak readout are new.

**H. Measure.** G's channels, presentations, smoothing, scaling (dominant-eye peak at 90% = 1), splits and onset
estimator, for all 35 non-blank cells (dominant-eye contrast × other-eye contrast).
* **Drive** of a cell = its peak response in 10–250 ms (the peak the onset estimator uses).
* **Monocular drive–latency curve,** per channel: least-squares line of onset on log peak through the eight
  monocular cells (each eye alone at 11–90%) with defined onsets (at least 5). Read at a binocular peak only up to
  1.5 × the largest monocular peak.
* **If the responses added:** for a binocular cell, the onset of D(c_D, t) + O(c_O, t), the sum of the two measured
  monocular responses of the channel, with the same estimator.
* Binocular cells for H1–H2: both contrasts ≥ 11% (16 cells). Per channel, the median over cells with defined
  values; then the median over channels per penetration; tests over penetrations.

**Tests.**
* **H1, earlier than adding the responses:** onset(both eyes) − onset(D + O) < 0, one-sided Wilcoxon.
* **H2, on the monocular drive curve:** onset(both eyes) − curve(peak(both eyes)); holds if the 95% CI of the median
  lies within ±1 ms and both halves' medians do (two-sided Wilcoxon reported).
* **H3, rank order:** per penetration, Spearman ρ between peak and onset over all cells with defined onsets; holds
  if ρ < 0 in ≥ 80% of penetrations, sign test p < .05, and both halves' medians are negative.
* **Reported (descriptive):** a peak-based linearity index at 11, 22.5 and 45% in each eye: the binocular peak placed
  on the dominant eye's peak-vs-contrast curve (`invert_crf`), with a peak-based eye weight; L_peak = (k − 1)/a,
  censored values (k·c > 90%) kept as lower bounds and counted. And the binocular peak minus the dominant eye's peak
  at c(1 + a).

**Criteria.** H1 and H3 hold if p < .05 on all presentations and both halves agree in sign. The paper may then
say: H1, "both eyes started earlier than the sum of the two monocular responses would have"; H1 + H2, "the
binocular response started when a single eye driving the same peak would have"; H3, "onset followed drive across
all 35 conditions". Script: `scripts/v4_drive_latency.py`, output [qc/v4_H.csv](qc/v4_H.csv).

**H results (2026-11-03; 26 penetrations, 570 channels; population medians per cell in
[qc/v4_H_cells.csv](qc/v4_H_cells.csv)). H1 and H3 hold; H2 fails.**
* **H1 holds:** both eyes start 1.5 ms earlier than the sum of the two monocular responses would (median over the 16
  cells; 25/26 penetrations; p < 10⁻⁴; odd −1.70, even −1.31 ms). Adding the responses would make the response
  larger (peak 1.80 at 90/90, against 1.23 measured) but not earlier.
* **H2 fails:** the binocular onset is 0.93 ms earlier than a single eye driving the same peak (95% CI −1.22 to
  −0.43; p = .002; 13/26 within ±1 ms; odd half −1.58). The gap is clearest at high drive: within one eye the onset
  gains only 1.3 ms from 45 to 90%, while the binocular cells keep getting earlier (−1.8 ms at 90/90).
* **H3 holds:** per penetration, the rank correlation between peak and onset has median −0.64 (24/26 negative; sign
  test p < 10⁻⁵; halves −0.51 and −0.55). Over the 35 population medians it is −0.98 (descriptive).
* **Reported, peaks:** the binocular peak does not differ from the dominant eye's peak at c(1 + a) (−0.028, −0.008,
  −0.017 at 11, 22.5, 45%; 95% CIs within −0.04 to +0.03; p ≥ .17). Peak-based linearity index 0.70 (CI
  0.55–1.03), 0.99 (0.60–1.15) and 0.91 (0.74–1.00; 38% censored, so lower bounds).
* **For the paper:** "Across all 35 conditions the response started sooner the stronger it was. The response to
  both eyes started 1.5 ms earlier than the sum of the two monocular responses would have, and about 1 ms earlier
  than a single eye driving the same peak. At their peaks, binocular responses did not differ from what summing the
  two eyes' contrasts predicts." Not: "on the same curve as a single eye" (H2 failed).
* **Figure:** Fig 1, last panel (Blake: "a very powerful visual to set up the entire paper"); mock in the
  scratchpad (`figpolish/v5/drive_mock.py`).

## V4 check I: does the other eye suppress where it barely drives? (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): Fig 3 is the laminar organization of normalization, and "we need to think carefully
about what Figure 3 could add to really show 'normalization' is occurring here", without repeating the model
figure. A binocular response below summation can come from normalization, from output saturation, or from response
decay; the summation prediction (the site's own dominant eye at c(1 + a), same window) already controls the last
two. The defining signature of normalization is suppression by an input that drives the site little (a suppressive
pool broader than the excitatory drive; Fig 4e's w > a is the model version). Here: at strongly monocular sites, the
other eye's grating should *lower* the response to the dominant eye (added < 0). Linear summation cannot do that,
and neither can saturation, because the contrast response does not supersaturate in this dataset (07). **Not
blind** to A3 (the late linearity index is most negative in the most monocular tercile, −0.56 at 45%; index units,
magnified by small eye weights) or to F4 (no undershoot over all sites). The response-unit test by ocular dominance
and depth is new.

**I. Measure.** Channels, presentations, smoothing, dominance, scaling and onset as `v4_dip_recovery.py`
(onset-aligned, dominant-eye peak at 90% = 1). Added A = R(c, c) − R(c, 0) per channel, in the F windows after
onset: early 0–20, dip 30–100, late 130–200 ms. Ocular dominance |ODI| (90%, 40–250 ms) from one half of the
presentations; A from the other half; both directions (odd → even, even → odd), with global tercile edges per
direction. Per penetration and tercile, the mean over channels; the two directions averaged; tests over
penetrations. 45% primary; 22.5 and 11% reported.

**Tests (one-sided Wilcoxon over penetrations).**
* **I1, suppression without drive:** most monocular tercile, A < 0 in 30–100 ms.
* **I2, delayed, not immediate:** most monocular tercile, A > 0 in 0–20 ms. (An interocular phase or disparity
  mismatch would act from onset; normalization is delayed. Caveat: disparity tuning also sharpens over tens of ms.)
* **I3, graded with the other eye's drive:** A(most binocular) − A(most monocular) > 0 in 30–100 ms.
* **I4, by depth (reported; "at every depth" only if all three hold):** I1 within depth groups relative to the
  earliest sink (above / at / below), terciles as above.
* **Reported:** the same in 130–200 ms; 11 and 22.5%; all-presentation values for display.

**Criteria.** I1–I3 hold if p < .05 and the two cross-validation directions agree in sign. The paper may then say:
I1, "at sites the other eye barely drives, its grating lowered the response to the dominant eye"; I1 + I2, "…but
only after the first ~20 ms"; I3, "the more monocular the site, the stronger the suppression". If I1 fails, Fig 3
cannot claim suppression without drive and falls back to contrast dependence (the transition index, 26/26). Script:
`scripts/v4_suppression.py`, output [qc/v4_I.csv](qc/v4_I.csv).

## V4 check J: does the current at the earliest sink start and peak before the MUA there? (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): "Shouldn't the CSD peak before the MUA? It does, right?", and whether the LFP filtering
smears it in time. **Filter chains (read from the code and the NSx headers, all 36 bminteroc files):** the LFP is
the amplifier's own 1 kHz stream (0.3 Hz first-order high-pass, 500 Hz fourth-order low-pass, both causal; no
software filtering), and the CSD is computed from it in space only (bad contacts interpolated, 3-point Hamming
across contacts). The MUA envelope is software, zero-phase: 500–5000 Hz band-pass, rectification, 200 Hz low-pass
(fourth order, `sosfiltfilt`), from the 30 kHz stream (0.3 Hz–7.5 kHz hardware). The figure caches then add a
Gaussian σ = 2 ms (zero-phase) to both. Both streams are read on the NSP tick clock, so the legacy 3.4 ms LFP/CSD lag
(05 A1) is gone. The two chains bias a CSD-vs-MUA comparison in opposite directions (causal delay on the CSD,
acausal spread on the MUA), so the comparison needs the bias measured. **Not blind:** `dichoptic_csd.py` recorded
"MUA onset minus sink onset", but with different thresholds (MUA half-max, sink 25% of peak), which favours the sink.

**J. Measure.** `sink.npz`: the sink current (−CSD over the sink channel and its two neighbours) and MUA on the same
channels, per penetration, on the sink-onset clock, for the dominant eye at 90% and both eyes at 45%. The same
estimators for both signals: half-max onset of its own peak (−5 to 150 ms, interpolated, as G) and peak time
(0–80 ms). Per penetration, MUA − current. **Filter bias:** both chains applied to the same test responses
(rises of 4–12 ms, 10–90%), the expected MUA − current difference read with the same estimators; reported results
are shown raw and with this bias removed.

**Tests (one-sided Wilcoxon over penetrations, raw and bias-corrected).** J1: the MUA half-max onset is later than
the current's. J2: the MUA peak is later than the current's. **Criteria:** a test holds if p < .05 raw *and* bias-
corrected. The paper may then say "the current at the earliest sink started (J1) / peaked (J2) X ms before spiking
on the same channels". Script: `scripts/v4_csd_timing.py`, output [qc/v4_J.csv](qc/v4_J.csv).

**I results (2026-11-03; 26 penetrations, 570 channels; |ODI| tercile edges 0.048 and 0.165;
[qc/v4_I.csv](qc/v4_I.csv), means in [qc/v4_I_summary.csv](qc/v4_I_summary.csv)). I1, I3 and I4 fail; I2 holds.**
* **I1 fails:** at the most monocular sites the other eye's grating does not reliably lower the response in
  30–100 ms at 45% (median −0.021 of the dominant-eye peak, 15/24 below zero, p = .20; both directions negative).
  Reported: not at 22.5% (+0.007) or 11% (+0.053) either.
* **I2 holds:** at the same sites the other eye adds at onset (0–20 ms, 45%: +0.095, 18/24, p = .002, both
  directions; 22.5% +0.080; 11% +0.059).
* **I3 fails** in 30–100 ms (+0.010, 13/22, p = .29). Reported, 130–200 ms: +0.024 (15/22, p = .018, both
  directions).
* **I4 fails** in every depth group (above −0.031, at −0.034, below −0.015; p = .13–.52).
* **Caveat:** MUA sites are seldom monocular. The most monocular tercile starts at |ODI| 0.17 (about 1.4 : 1).
* **For the paper:** no claim of interocular suppression without drive in MUA. Fig 3 falls back to contrast
  dependence, as registered.

**J results (2026-11-03; 26 penetrations; [qc/v4_J.csv](qc/v4_J.csv),
[qc/v4_J_penetrations.csv](qc/v4_J_penetrations.csv)). J1 and J2 fail: spiking on the sink channels starts before
the sink current.**
* **Clock check (QC; `scripts/qc_nsx_clock.py`, two sessions per monkey):** the 1 kHz LFP equals the 30 kHz signal
  passed through the amplifier's causal 500 Hz filter at lag 0 (r = 1.000). Against a zero-phase version it lags by
  0.9 ms. The two streams share one clock; the LFP and CSD carry the 0.9 ms delay of the causal filter.
* **Filter chains, simulated (rises 4–12 ms, 10–90%):** current half-max +0.6 ms late, peak +1 to +3 ms; MUA
  half-max −0.1 to −0.2 ms, peak 0 to +2 ms. Expected MUA − current bias −0.8 ms (half-max) and −1 ms (peak).
* **Half-max onset, MUA − current:** dominant eye 90% −3.3 ms raw (MUA later in 4/26), −2.5 ms corrected (95% CI
  −5.1 to −0.5); both eyes 45% −3.0 raw, −2.2 corrected (CI −3.7 to −0.3).
* **Peak, MUA − current:** dominant eye 90% −4.5 raw, −3.5 corrected (CI −6 to +1); both eyes 45% −3.0 / −2.0
  (CI −10 to +2).
* **Descriptive (not registered):** at every crossing threshold from 10 to 75% of its own peak the MUA leads (−3.8,
  −3.5, −3.3, −3.5 ms at 10, 25, 50, 75%, dominant eye 90%).
* **Reading (hypothesis, not tested):** analog MUA at the input depth includes thalamic afferent (presynaptic)
  spikes, and the sink grows as synaptic input accumulates; the three-contact average and spatial smoothing blend
  in neighbouring, later currents.
* **Sensitivity (descriptive, not registered; `scripts/v4_csd_sensitivity.py`,
  [qc/v4_J_sensitivity.csv](qc/v4_J_sensitivity.csv)).** The current corrected for the 0.9 ms delay; CSD with or
  without Hamming smoothing, at the sink contact alone or over three contacts, with or without the 2 ms Gaussian.
  The MUA half-max lead shrinks with the most local CSD (single contact, no smoothing at all) to 0.5 ms (dominant eye
  90%, 13/26 MUA earlier), 1.4 ms (both eyes 45%, 17/26) and 2.6 ms (both eyes 90%, 16/26); over all variants 0.5–2.7
  ms. Peaks range from 5 ms earlier to 1 ms later. Blending neighbouring contacts accounts for part of the lead; MUA
  still reaches half its peak first in most variants. Literature: reports/CSD sink versus spiking onset order.md.
* **For the paper:** do not say the current precedes spiking at the sink. "Earliest sink" marks the input depth, not
  a time reference.

## V4 check M: do single units at the sink follow the sink current while MUA leads it? (registered internally 2026-11-03, before running)

**Why.** The literature review ([reports/CSD sink versus spiking onset order.md](../reports/CSD%20sink%20versus%20spiking%20onset%20order.md))
found analog-MUA studies reporting sink and spiking as coincident, and the one single-unit study (Westerberg et al.
2019) reporting the sink first. Single-axon work puts thalamic terminal spikes ~0.7 ms before their sink, and
cortical spikes after it. Blake (2026-11-03): run the single-unit test. **Prediction of the afferent account:**
sorted units (mostly cortical cells) start no earlier than the sink current, while analog MUA on the same contacts
starts before it; the MUA lead over units is largest at the input depth. **Not blind** to J (MUA half-max leads
the current by 0.5–2.7 ms depending on processing); unit timing has not been looked at.

**M. Measure.**
* **Units:** the lenient set (as `v3_units_index.py`), in cortex, with a response: peak of the smoothed rate in
  10–150 ms above the baseline mean (−50 to 0 ms) + 5 baseline SD. Strict (well-isolated) subset reported.
* **Presentations:** clean, either eye ≥ 45% (the sink-selection set; a strong-drive pool), blank presentations
  subtracted. Reported: dominant eye at 90%, both eyes at 45%.
* **Cross-validated sink contact:** the `dichoptic_csd.py` rule (earliest among strong sinks reaching 25% of its own
  peak, dataset CSD) applied to odd presentations, all signals measured on even ones, and the reverse. Units within
  2 contacts of that contact are "at the sink"; the others above or below.
* **Signals (1 ms, Gaussian σ = 2 ms):** the unit's rate; analog MUA on the unit's peak contact; the sink current
  (−CSD) at the selected contact from unsmoothed LFP (primary; the most local estimate), advanced by the measured
  0.9 ms filter delay. Reported: the three-contact, Hamming-smoothed current of J.
* **Estimators:** half-max of each signal's own peak (10–150 ms, interpolated; primary), 25% crossing and peak time
  (reported). Per unit, the mean of the two directions.

**Tests (units as the sample; per-penetration medians reported).**
* **M1, the sink leads units:** unit − current > 0 for units at the sink (one-sided Wilcoxon).
* **M2, MUA leads units on the same contact:** unit − MUA > 0 for units at the sink (one-sided).
* **M3, reported:** MUA − current on the same contacts (J's comparison on these contacts).
* **M4, the MUA lead is largest at the input:** (unit − MUA) at the sink > (unit − MUA) above and below (one-sided
  Mann–Whitney, units as the sample).
* **Criteria:** a test holds if p < .05 and both directions agree in sign; with fewer than 10 units at the sink the
  result is inconclusive, not a failure. The paper may then say, for M1 + M2: "isolated neurons at the input depth
  started after the sink current, while multiunit activity on the same contacts started with or before it",
  consistent with thalamic afferent spikes in the MUA. Script: `scripts/v4_units_sink.py`, output
  [qc/v4_M.csv](qc/v4_M.csv).

**M results (2026-11-03; responsive lenient units at the cross-validated sink: 27 in 13 penetrations, 10 strict;
above 40, below 42; per-unit values in [qc/v4_M_units.csv](qc/v4_M_units.csv)). M2 holds; M1 and M4 fail.**
* **M1 fails (primary: single-contact, unsmoothed current):** units at the sink start +3.1 ms after the current
  (15/27, p = .41), and the two directions disagree (+5.1, −1.4 ms). Reported, the three-contact Hamming current:
  units +4.6 ms after it (20/27, p = .009, both directions positive; per penetration 9/13). Reported, other
  conditions (both eyes 45%, the unit's dominant eye 90%): units 7–20 ms after the current (p ≤ .02), with fewer
  presentations and noisier unit onsets.
* **M2 holds:** at the sink, units start 4.1 ms after the analog MUA on their own contact (20/27, p = .003; both
  directions +2.0, +2.3 ms; per penetration 10/13). 25% crossings: +2.6 ms (17/27, p = .014). Peaks do not differ
  (+1.0 ms, p = .07).
* **M3, reported:** MUA − current −3.1 ms on these contacts (20/27, p = .008), but the two directions disagree.
* **M4 fails:** units lag the MUA by as much away from the sink as at it (above +4.8 ms, 32/40; below +3.5 ms,
  31/42; at +4.1 ms; Mann–Whitney p = .58). The MUA's lead over single units is not specific to the input depth.
* **Strict subset:** 9–10 units at the sink, inconclusive by the registered rule.
* **Reading.** Single units at the input depth start after the sink current in most variants (as Westerberg et al.
  2019 found), but the registered single-contact test is not reliable. Analog MUA reaches half its peak ~4 ms before
  the typical isolated neuron at every depth, so its lead over the sink is better explained by MUA being a
  population signal dominated by the earliest cells than by thalamic afferent spikes at the input (M4). For the
  paper: "analog MUA reached half its peak about 4 ms before isolated neurons on the same contacts, at every depth;
  isolated neurons at the input depth started after the sink current (three-contact estimate; the single-contact
  estimate was too noisy to resolve the order)". Timing claims compare a signal with itself, as before. Methods: state both filter chains and the 0.9 ms LFP delay. Legacy CSD-vs-MUA latencies carry
  the 3.4 ms artefact (05 A1) and must not be reused.

## V4 checks K and L: normalization and the model across depth (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): "The same strength at every depth is perfectly fine"; the model deserves its own
figure, "but we want to show the modeling effects across time and depth, together. We need to make better use of
spatiotemporal ways of showing the data." **Not blind (L):** a descriptive depth × time map (scratchpad
`figpolish/v5/fig3_mock.py`, `depth_share.py`) shows, at 45% and 30–100 ms after onset, both eyes − "if the contrasts
add" of −0.02 to −0.06 at and above the earliest sink (the other eye's contribution slightly negative) and about 0
from 0.6 mm below it (the other eye keeps adding about the prediction), with similar predictions (0.03–0.05) and eye
weights (0.64–0.72) at all depths. The V3 index-based depth tests found no depth dependence of the drop (mixed model
p = .60) and an inverted U in monkey J only.

**L. Measure.** Per channel, onset-aligned as `v4_dip_recovery.py`, all / odd / even presentations: added A and
predicted P (dominant eye at c(1 + a) − dominant eye at c, the channel's own a) in 0–20, 30–100 and 130–200 ms.
Per penetration and depth group (above / at / below the earliest sink, as S4c), the mean over channels.
* **L1, normalization in every group:** A − P < 0 in 30–100 ms at 45%, one-sided Wilcoxon per group.
* **L2, weaker below the input:** (A − P)below − (A − P)at > 0 in 30–100 ms at 45%, one-sided, paired over
  penetrations; (A − P)above − (A − P)at reported two-sided.
* **Reported:** 22.5 and 11%; 0–20 and 130–200 ms.
* **Criteria:** a test holds if p < .05 on all presentations and both halves' medians have the predicted sign.
  L2 holding would let the paper say "below the input layer the other eye kept adding what summation predicts while
  at and above it its contribution was cancelled"; if L2 fails, "normalization at every depth" stands.

**K. Measure.** Per driven channel, stimulus clock 0–219 ms (the `fig5_model.py` data, one channel at a time, scaled
by its own dominant-eye peak at 90%): the no-delay ('static') and fixed-delay ('full') models fitted to odd and to
even presentations (each scored on the other half; held-out error = summed held-out SSE) and to all presentations
(parameters, predictions), four starts each. Per penetration and depth group, the median over channels.
* **K1, the delayed pool fits better at every depth:** log(held-out error static / full) > 0 in each group,
  one-sided Wilcoxon; "at every depth" only if all three hold.
* **K2, the pool weights the other eye more than the drive at every depth:** w − a > 0 in each group.
* **K3, reported:** depth differences (above − at, below − at; two-sided Wilcoxon) in the pool delay τ2, w − a and
  the latency δ.
* **Reported (display):** depth × time maps on the response clock of what the other eye adds and of both eyes −
  "if the contrasts add", for data, static and full predictions, and their residuals (Fig 4 in the space-time
  figure set).
* **Criteria:** as for L (p < .05, both cross-validation directions or halves agree). Scripts:
  `scripts/v4_depth_norm.py` (L, [qc/v4_L.csv](qc/v4_L.csv)) and `scripts/v4_model_depth.py`
  (K, [qc/v4_K.csv](qc/v4_K.csv)).

**L results (2026-11-03; 26 penetrations; group means in [qc/v4_L_groups.csv](qc/v4_L_groups.csv)). L2 holds:
normalization is weaker below the input layer.** At 45%, 30–100 ms after onset (both eyes − "if the contrasts add",
dominant-eye peak = 1):
* **L1:** holds above the sink (−0.073, 21/23, p = .0001) and at it (−0.065, 22/25, p = .0001), both halves; fails
  below it (−0.011, 18/25, p = .082; odd half +0.015).
* **L2 holds:** below − at = +0.049 (18/25, p = .0003; halves +0.060, +0.043). Above − at: −0.009 (p = .57).
* **What the other eye adds** (medians of group means): above −0.025, at −0.029, below +0.007, against predictions
  of 0.048, 0.025 and 0.037. At and above the sink its contribution is cancelled and slightly negative; below it the
  other eye still adds a little, short of the prediction (not tested as a net suppression).
* **Reported:** 22.5% the same pattern (L2 +0.037, p = .006, both halves; 130–200 ms +0.021, p = .045); 11% no depth
  difference (L2 p = .12). At onset (0–20 ms, 45%) sites below the sink add more than predicted (+0.033; L2 +0.044,
  p = .040).
* **Not blind:** the pattern was seen in the population depth map before the test (above). The earlier index-based
  "no depth dependence" (V3) and this response-unit result differ because the index divides by the prediction and
  pools windows; the text should use the response-unit result, state the windows, and keep the index result in the
  supplement.
* **For the paper:** "From 30 to 100 ms after onset, at and above the earliest sink the other eye's contribution was
  cancelled; below it, the other eye kept adding (at 45% and 22.5%; the difference, 18/25 penetrations)." This
  replaces "the same strength at every depth".

**K results (2026-11-03; 24 penetrations, 520 channels; 211008 and one other lack two presentations per cell in a
half; group medians in [qc/v4_K_groups.csv](qc/v4_K_groups.csv)). K1 and K2 hold at every depth.**
* **K1:** the delayed pool predicts held-out data better than no delay above the sink (log error ratio 0.63, 19/21),
  at it (1.04, 22/23) and below it (0.55, 21/23); all p < 10⁻⁴.
* **K2:** the pool weights the other eye more than the drive does above (w − a = 0.24, 20/21), at (0.20, 21/23) and
  below (0.09, 21/23) the sink; all p < 10⁻⁴.
* **Deviation from the registered criterion:** the held-out error sums both cross-validation directions, so the
  "both directions agree" check was not computed for K1, and K2 uses the all-presentation fits. Both are reported as
  p < .05 only; a rerun can save each direction.
* **K3 (reported):** below the sink the excess pool weight is smaller (w − a below − at −0.06, 7/23 above zero,
  p = .030) and the delay helps less (log ratio −0.36, 3/23, p < 10⁻⁴), consistent with L. Pool delay τ2 shows no
  depth difference (medians 20, 29, 28 ms; above − at p = .35, below − at p = .67). Latency δ +2.5 ms above (p = .07).
* **Maps (mock, scratchpad `figpolish/v5/fig4_mock.py`):** at 45% the no-delay fits add a small constant amount
  everywhere and miss the onset burst; the delayed-pool fits reproduce the burst at every depth, slightly
  underestimate the 30–60 ms dip at and above the sink, and miss the late return (strongest deep).

## V4 check N: the four models by depth group, both cross-validation directions (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): on the model, "whatever you recommend … it's worth being thorough". K fitted only the
no-delay and fixed-delay models per channel and summed the two held-out directions. N fits all four models of the
paper per depth group, where channel averages are less noisy, and keeps each direction; K is rerun with the
directions kept (K1′). **Not blind** to K (delay and w > a at every depth, smaller excess below the sink) or to the
per-penetration model results (fifth model test: conductance better than fixed delay overall).

**N. Measure.** Per penetration and depth group (above / at / below the earliest sink, ≥ 3 driven channels): the
`fig5_model.py` population restricted to the group's channels (channel mean, 0–219 ms, σ = 2 ms, each channel
oriented to its dominant eye and scaled by its dominant-eye peak at 90%), all / odd / even presentations. Models:
'static' (no delay), 'full' (fixed delay), 'w=a' (pool weighted like the drive), 'cond' (conductance; extra start at
the full fit, as `fig5_cond_model.py`); four starts each. Held-out error per direction: fit odd, score even
(o→e), and the reverse (e→o). Parameters from the all-presentation fits.

**Tests (one-sided Wilcoxon over penetrations within each group; a test holds if p < .05 on the summed held-out
error and the median has the predicted sign in both directions; "at every depth" only if all three groups hold).**
* **N1, the delay:** log(CV static / CV full) > 0.
* **N2, the pool weights the other eye differently from the drive:** log(CV w=a / CV full) > 0.
* **N3, conductance over fixed delay:** log(CV full / CV cond) > 0.
* **N4, reported:** depth differences (above − at, below − at; two-sided Wilcoxon) in τ2, w − a, δ and τ_m (cond
  fits), and in the held-out gains.
* **K1′:** K1 with each direction kept (same channels and fits; the criterion as registered for K).
* Script: `scripts/v4_model_groups.py`, output [qc/v4_N.csv](qc/v4_N.csv); `scripts/v4_model_depth.py` rerun for K1′.

**N results (2026-11-03; 24 penetrations, 60 groups: above 14, at 23, below 23). N1, N2 and N3 hold in every
group; K1′ holds.**
* **N1, the delay:** log(CV static / CV full) 1.10 above (13/14), 1.25 at (22/23), 0.91 below (21/23); all
  p ≤ .0001; both directions.
* **N2, pool weighting:** log(CV w=a / CV full) 0.055 (12/14), 0.017 (21/23), 0.009 (19/23); p ≤ .0003; both
  directions. Small gains, consistent everywhere.
* **N3, conductance:** log(CV full / CV cond) 0.004 (11/14, p = .012), 0.008 (19/23, p = .0003), 0.024 (19/23,
  p < 10⁻⁴). Above the sink one direction's median is barely above zero (< 0.0001): the conductance gain there is
  marginal.
* **N4, by depth (reported):** the excess pool weight w − a is smaller below the sink (conductance fits: 0.23, 0.19,
  0.13; below − at −0.095, 5/22 above zero, p = .003), and the delay helps less there (p = .007), both as in K and L.
  No depth difference in the pool delay τ2 (33, 37, 39 ms; p = .30, .59) or the conductance time constant τ_m
  (3.3, 8.2, 19.8 ms; below − at p = .16). The input delay δ is 3.8 ms longer above (p = .049) and 1.5 ms shorter
  below (p = .023) than at the sink.
* **K1′ holds:** the per-channel delay gain is positive in both directions in every group (above 0.65 / 0.62, at
  1.02 / 1.03, below 0.51 / 0.56). The K deviation is resolved.
* **For the paper:** "At every depth, a delayed normalization pool that weights the other eye more than the drive
  does predicted held-out responses better than no delay or an equally weighted pool; acting as a conductance
  improved it slightly further. The extra weight on the other eye was smallest below the input layer, where the
  other eye kept adding (Fig 3). The pool delay, 33–39 ms, did not differ with depth."
* **Display:** the conductance model per channel (all presentations, no tests; `v4_model_depth.py 6 cond`) for the
  Fig 4 maps.
* **What the model misses (descriptive, Fig 4b, c, h):** on the response clock at 45%, data − conductance fits is
  +0.016 in 0–20 ms (the onset burst is captured at every depth), −0.003 in 30–100 ms averaged over depth but grey
  at and above the sink (the data dip deeper), and +0.031 in 130–150 ms (the late return, strongest deep, is
  missed). By depth group, both eyes − "if the contrasts add" (30–100 ms, medians): data −0.056 / −0.061 / −0.020
  above / at / below; conductance fits −0.039 / −0.020 / −0.017. The fits normalize about equally at every depth:
  they underestimate normalization where it is strongest and do not reproduce the laminar difference of check L,
  although their extra pool weight is smaller below. For the paper: a limit of a single-stage delayed pool fitted
  per depth; candidate next model (needs its own registration): a stronger or faster pool at and above the input.

## V4 checks O and P: summation and normalization in the synaptic current; which depth behaves like which model (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): "can we ask our questions of additivity and normalization at sub-spiking-threshold
CSD? Does it hold?"; and on the model, no laminar model: "report the effects we see across layers and which layer
behaves like which model". **Nulls for a current.** Spiking has a threshold and saturates, so the spiking null is
summation of the eyes' *contrasts*. Synaptic current is closer to linear, and at the input the two eyes arrive on
separate afferents, so the natural null for the current is that the *currents* add: both eyes = dominant eye alone +
other eye alone (the additivity index, Blake's legacy measure). **Not blind** to the earlier sink results (the sink
current adds at onset, Fig 3e of V4, 22/26; its late contrast response is flat between 45 and 90%; MUA on the sink
channels leads it), or to N (gains by depth).

**O. Measure.** Per penetration with a usable earliest sink: the sink current (−CSD, dataset CSD, over the sink
channel and its two neighbours, as `sink_data.py`) and MUA on the same channels; condition means over clean
presentations, blank subtracted, σ = 2 ms, on the sink-onset clock (`dichoptic_csd_sinks.csv`); the dominant eye =
larger current at 90% in 0–50 ms; all / odd / even presentations (dominance and scale from all). Windows after sink
onset: early 0–20, dip 30–100, late 130–200 ms. Per window, the additivity index of window means,
B / (D + O) (1 = the currents add; < 1 = less than the sum); and both eyes − "if the contrasts add" (the current's own
eye weight over 0–50 ms, reported). 45% primary; 22.5 and 11% reported.
* **O1, the input currents add at onset:** current additivity index in 0–20 ms; holds if the 95% bootstrap CI of the
  median lies within 0.8–1.2 and both halves' medians do.
* **O2, the currents normalize later:** additivity index 30–100 ms < 0–20 ms (paired, one-sided Wilcoxon) and
  30–100 ms < 1 (one-sided).
* **O3, spiking is less additive than the current feeding it:** MUA additivity index < current additivity index in
  0–20 ms (paired, one-sided); 30–100 ms reported.
* **Reported (display):** depth × time maps of the current's additivity (both eyes − (dominant + other), signed so
  that negative = less current than the sum) across all channels, depth relative to the sink.
* **Criteria:** O2 and O3 hold if p < .05 and both halves agree in sign. Script `scripts/v4_csd_additivity.py`,
  output [qc/v4_O.csv](qc/v4_O.csv).

**O results (2026-11-03; 26 penetrations; per penetration in [qc/v4_O_penetrations.csv](qc/v4_O_penetrations.csv);
n is lower where the monocular currents in a window sum to ≤ 0). O1 and O3 fail; O2 holds against 1 but not
against onset.**
* **O1 fails:** the currents do not add at onset. At 45%, in the first 20 ms after sink onset, the binocular current
  is 0.65 of the sum of the two monocular currents (95% CI 0.56–0.76; n = 23); 0.66 at 11% and 0.72 at 22.5%.
* **In contrast units the current follows the spiking rule at onset:** both eyes − "if the contrasts add" (the
  current's own contrast response and eye weight) is +0.02 at 45% (n.s.), +0.001 at 22.5%, −0.016 at 11%: the
  binocular current is what the dominant eye's current would be at the summed contrast.
* **O2:** later the current is further below the sum (30–100 ms: 0.49, 16/18 below 1, p < 10⁻⁴; 130–200 ms 0.43), but
  the drop from onset is not significant at 45% (−0.08, 11/17, p = .10) or 22.5% (p = .13); at 11% it is (−0.13,
  13/18, p = .007, both halves). In contrast units the late current is 0.03 below summation at 45% (16/26,
  p = .046), small because the current's late contrast response is flat between 45 and 90%.
* **O3 fails:** MUA on the same channels is not less additive than the current (0–20 ms: 0.72 vs 0.65; MUA − current
  +0.10, 6/23 below zero).
* **Maps (descriptive; [csd_additivity.npz] in the figure data, mock `figpolish/v5/csd_add_mock.py`):** at every
  depth, sinks and sources alike, the binocular current is well below the sum of the monocular currents from the
  first milliseconds after sink onset.
* **Reading.** The compression that makes "contrasts add" the right null is already in the current at the input,
  not added by the spike threshold: the eyes' currents combine like contrasts, not like independent inputs. Likely
  reasons (not tested): the sink current is mostly intracortical even in its first 20 ms (Happel et al. 2010: ~88%
  of the granular sink), and summed conductances do not add linearly. The time-dependent normalization (the
  onset-to-dip drop) is clear in spiking and weak in the current, except at 11%.
* **For the paper:** "Below the spiking threshold the rule was the same: the current at the earliest sink combined
  the eyes' contrasts, not their currents (the binocular current was 65% of the sum of the monocular currents in the
  first 20 ms, and matched the dominant eye's current at the summed contrast). Its later fall below summation was
  small and significant only at 11%, so the normalization dynamics are clearest in spiking."

**P. Measure.** N's depth groups, data and cross-validation, adding the two remaining models of the paper:
recurrent normalization ('org', ORGaNICs) and the two-stage model ('2stage'). For each group and penetration, the
held-out error of the six models (static, full, w=a, cond, org, 2stage).
* **P1, reported per group:** which model has the lowest held-out error (counts of penetrations), and each model's
  % less error than 'full'.
* **P2:** org vs cond, log(CV cond / CV org) by group (one-sided > 0 per group, Wilcoxon): does recurrent
  normalization describe any depth better than the conductance form? And below − at (two-sided): does the better
  model change with depth?
* **P3:** 2stage vs full, log(CV full / CV 2stage) by group (one-sided > 0).
* **Criteria:** as N (p < .05 on the summed held-out error, both directions in the predicted sign). The paper may
  say "depth X behaved like model Y" only for a held test. Script: `scripts/v4_model_groups.py` with the two models
  added (output [qc/v4_P.csv](qc/v4_P.csv)).

**P results (2026-11-03; 24 penetrations, 60 groups). P3 holds at and below the sink; P2 fails everywhere.**
* **P1, best model per group (lowest held-out error, counts):** above (14): two-stage 6, conductance 5, static,
  fixed delay and recurrent 1 each; at the sink (23): two-stage 15, recurrent 5, conductance 3; below (23): two-stage
  9, conductance 7, recurrent 4, static 2, fixed delay 1. % less error than the fixed delay (medians): two-stage
  +1.8 / +2.6 / +1.4, conductance +0.4 / +0.8 / +2.4, recurrent +1.1 / −5.1 / −0.8 (above / at / below).
* **P2 fails:** recurrent normalization does not beat the conductance form at any depth (log ratio −0.008, −0.082,
  −0.027; 5–6 penetrations in its favour). It is relatively less bad below than at the sink (p = .04) but never
  better.
* **P3:** the two-stage model beats the fixed delay at the sink (+0.026, 19/23, p = .0001, both directions) and below
  it (+0.014, 18/23, p = .0007, both directions); above, not reliably (+0.019, 9/14, p = .059).
* **Not tested:** two-stage vs conductance directly (P1 counts only).
* **For the paper ("which depth behaves like which model"):** at and below the earliest sink, the combination is
  best described by a two-stage model, in which each eye's signal is first gain-controlled by the other eye before
  the eyes combine and a delayed binocular pool normalizes; below the sink the conductance form helps as much. Above
  the sink no form beats the fixed-delay pool reliably, and what matters most there is the pool weighting the other
  eye separately from the drive (N2, 5.3%). Recurrent normalization is not the best description at any depth. All
  these gains (1–3%) are small beside the delay itself (60–71%).
* **Figure:** Fig 4e now shows all four comparisons by depth; best-model counts in the legend.

## Strength and novelty, re-assessed (2026-11-03; Blake asked to step back)

An independent, skeptical review (subagent given the figures, the qc tables and the literature notes, not our
conclusions) and our own read agree on most points.

* **Robust:**
  * onset follows drive, and the second eye counts as drive (G, H; both monkeys in sign; expected under any
    normalization model, so not diagnostic of one);
  * the eyes' contrasts sum for ~20 ms after onset and then fall short of summation (F1, Fig 2a), clearest at 45%;
  * the transition deepens with contrast (26/26);
  * the pool weights the other eye more than the drive (24/25; K2, N4).
* **Moderate to weak:**
  * the laminar difference (L2 holds, but it was seen before it was tested, is ~0.05 of peak, B is 5/7, and the
    index-based depth test was null);
  * normalization time by depth (Fig 3d, descriptive; the earlier timing measures failed split-half reliability);
  * conductance, two-stage and recurrent gains (0.4–3%);
  * all CSD results.
* **Known, so incremental:**
  * the summation-then-normalization sequence (Cox et al. 2019; Dougherty et al. 2019; Mitchell et al. 2022);
  * delayed normalization (Carandini, Heeger & Movshon 1997 and later);
  * interocular normalization weighted by ocular dominance, **time-averaged, in these same recordings** (Mitchell
    et al. 2023, Curr Biol: B + J, 25 penetrations, 502 MUA channels; docs/03 item 6);
  * a 2023 SfN abstract on normalization across condition, depth and time (docs/03 item 7).
  * The reviewer reports the 2024 dissertation says normalization is "weakest in the geniculate input layers and
    strengthened through the cortical circuit", which needs reconciling with L (strong at the sink). Blake to check.
* **Undersold (reviewer, post hoc, re-run and reproduced here; not registered):** in the unequal-contrast cells
  with a strong other eye, the shortfall from summation is larger and sustained. Both eyes − "if the contrasts add",
  130–200 ms: dominant 11% + other 45% −0.033 (23/26), + other 90% −0.048 (20/23); dominant 22.5% + other 90%
  −0.046 (16/18). On the equal-contrast diagonal the effect is transient at 45% (check F). This uses the full 6 × 6
  design and extends Mitchell et al. 2023 into time.
* **Main objections:** same data as Mitchell et al. 2023; two monkeys with B underpowered (8 penetrations) and
  weaker in several tests; MUA sites nearly all binocular, so onset summation could be across neurons; analysis
  flexibility (index vs response units disagree on depth); zero-disparity gratings with one eye tracked (disparity
  suppression lags excitation by ~7 ms, Tanabe & Cumming 2014).
* **Framing:** "1 ms resolution of binocular summation" overstates: onsets are resolved to ~1 ms, the combination
  rule to tens of ms.
* **Figure problems flagged:**
  * Fig 2c–d grey reads as suppression, which F4 rejected;
  * Fig 2b's sink/source bands above the sink may be contact-gain stripes, and the onset overlay invites the
    current-before-spikes reading that J rejected;
  * Fig 2a bands are wide, and a negative index is not a "share added";
  * Fig 3d has no test behind it;
  * Fig 4's lead comparison is against a no-delay model with no transient (a straw man);
  * Fig 4a is drawn with exponent 2 (fitted ~1);
  * Fig 1f's grey band extrapolates beyond the monocular range.
* **Tier:** J Neurosci or eLife as is; Curr Biol a stretch.
* **Proposed action list (awaiting Blake):**
  1. one hierarchical model for the main claims (monkey, session, penetration);
  2. the 6 × 6 time course as the centrepiece (registered, with a held-out test on the dual-probe cross-recordings);
  3. depth robustness for L (L4/5-landmark alignment, continuous depth, deep-sink exclusion);
  4. units hierarchically by ocular dominance;
  5. the figure fixes;
  6. an explicit novelty statement against Mitchell et al. 2023, the SfN abstract and the dissertation.
* Reviewer scripts (post hoc, not in the repo): scratchpad `offdiag.py`, `shift_null.py`. Shifting the summation
  prediction 1 ms earlier halves the early surplus (+0.031 → +0.017) and leaves the 30–100 ms shortfall (−0.037,
  23/26).

## V4 checks Q–U: the action list (registered internally 2026-11-03, before running; Blake: "Please proceed")

**Q. The full 6 × 6 design over time.** *Not blind:* the reviewer's descriptive cell medians (scratchpad
`offdiag.py`, above). Per penetration, onset-aligned channel means as check F (`v4_dip_recovery.penetration`; all /
odd / even presentations; dominance, scale, onset and the median eye weight a from all). For every binocular cell
with both contrasts ≥ 11% and c_D + a·c_O ≤ 90%: shortfall S = window mean of R(c_D, c_O) − R_D(c_D + a·c_O)
(`response_at`; dominant-eye peak = 1) in 0–20, 30–100 and 130–200 ms after onset. Linear mixed model (statsmodels,
REML): S ~ window × (x_O + x_D) + monkey, x = log2(c / 22.5%), random intercept per session with a penetration
variance component.
* **Q1, graded with the other eye's contrast after onset:** the x_O slope < 0 in 30–100 ms and in 130–200 ms.
* **Q2, develops after onset:** the x_O slope more negative in 30–100 and in 130–200 ms than in 0–20 ms
  (window × x_O < 0).
* **Q3, sustained when the other eye is stronger:** in 130–200 ms, S < 0 over the cells with c_O > c_D (the same
  model's random structure, intercept only).
* **Criteria:** one-sided p < .05, and the same sign refitted on odd and on even presentations and in each monkey
  alone.

**R. Held-out confirmation on the dual-probe cross-recordings** (`output/datasets_cross`: nine monkey J runs in which
the grating sat on the other probe's receptive field; used before only for the candidate-1 statistic, docs/10).
Driven channels by that test's rule (one-way ANOVA across the dominant eye's monocular contrasts, 40–250 ms,
p < .01); recordings with ≥ 5 driven channels. Same S and model, with a recording random intercept. **R1:** Q1 in
130–200 ms (x_O slope < 0); **R2:** Q3. A test holds if one-sided p < .05; with fewer than 5 recordings, inconclusive.

**S. Depth robustness of L2** (below − at the sink, 45%, 30–100 ms; per channel, `v4_L_channels.parquet`).
* **S1:** continuous depth relative to the sink: per-channel shortfall ~ depth, penetration random intercept;
  slope > 0.
* **S2:** the same with depth relative to the L4/5 landmark (laminar alignment v2).
* **S3:** L2 excluding penetrations whose earliest sink lies > 0.5 mm below L4/5 or > 1.5 mm below the top of cortex
  (deep sinks).
* **S4:** L2 with the sink chosen on one half of the presentations and the shortfall measured on the other (both
  directions averaged).
* **Criteria:** each variant holds if one-sided p < .05; "robust to the depth reference" only if all four hold.

**T. Single units by ocular dominance, hierarchically.** Lenient driven units (as `v3_units_index.py`): rates,
condition means (≥ 3), blank subtracted, σ = 3 ms, the unit's dominant eye, scaled by its dominant-eye peak at 90%,
half-max onset, eye weight a; at 45% in each eye: added A and predicted P in 0–20, 30–100 ms after the unit's onset;
|ODI| (90%, 40–250 ms). Mixed model with unit and penetration nested in session, monkey fixed.
* **T1:** the other eye adds at onset, A(0–20) > 0.
* **T2:** normalization in single neurons, A − P (30–100 ms) < 0.
* **T3:** T1 in the more monocular half of units (|ODI| above the median): summation within neurons, not only across
  them.
* **Criteria:** one-sided p < .05, same sign in both halves of the presentations.

**U. The main claims in one hierarchical model** (penetration as the observation, random intercept per session,
monkey as a fixed, sum-coded effect): U1 the 45% shortfall in 30–100 ms < 0 (F); U2 the slope of the transition
index on log contrast < 0 (per penetration and contrast, penetration nested in session); U3 the binocular onset
advance at 45% < 0 (G1); U4 the per-penetration rank correlation of onset and peak < 0 (H3); U5 below − at > 0
(L2); U6 the excess pool weight w − a > 0 at each depth (N, conductance fits). Each holds if the intercept is
significant one-sided at p < .05; the monkey effect is reported two-sided (does B differ from J?).

Scripts: `scripts/v4_full_design.py` (Q, R), `scripts/v4_depth_robust.py` (S), `scripts/v4_units_od.py` (T),
`scripts/v4_hierarchical.py` (U); outputs `docs/qc/v4_Q.csv` … `v4_U.csv`.

**Q–U results (2026-11-03).**
* **Q holds in full** (26 penetrations, 22 sessions, 843 cell × window observations).
  * The shortfall from summation grows with the other eye's contrast in 30–100 ms (slope −0.013 per doubling,
    p = .0002) and in 130–200 ms (−0.011, p = .0009).
  * Both slopes are more negative than at onset (−0.024 and −0.022, p < 10⁻⁴). At onset the slope is positive
    (+0.011): a stronger other eye adds slightly *more* than summation.
  * In 130–200 ms, with the other eye stronger than the dominant eye, both eyes fall 0.026 below summation
    (p < 10⁻⁵); on the matched diagonal, −0.010.
  * Every estimate has the same sign in the odd and even refits and in each monkey alone (late, other eye stronger:
    B −0.016, J −0.036).
* **R, held out (8 cross-recordings with ≥ 5 driven channels; 221209_J_eA had 3):**
  * **R2 holds:** with the other eye stronger, 130–200 ms, −0.024 (p = .032; halves −0.011, −0.037).
  * **R1 fails:** the late slope on the other eye's contrast is −0.006 (p = .20); the 30–100 ms slope is −0.011
    (p = .054).
  * These recordings had the grating on the other probe's receptive field, so their responses are weaker. The
    sustained shortfall replicates; its grading does not reach significance.
* **S holds in all four variants:** shortfall ~ depth below the sink +0.025 per mm (p < 10⁻⁴); relative to L4/5
  +0.025 per mm (p < 10⁻⁴, 22 penetrations); without the deep-sink penetrations (211027_B, 211208_B, 221206_J_eB)
  +0.053 (16/22, p = .0007); with a cross-validated sink +0.048 (17/24, p = .001). The laminar difference does not
  depend on the depth reference.
* **T (65 lenient units, 20 penetrations; B 30, J 35):**
  * **T2 holds:** single neurons fall short of summation in 30–100 ms (−0.076, p < 10⁻⁴, both halves), also in the
    more monocular half (−0.077, p = .003).
  * **T1 fails:** the other eye's addition at onset is +0.072 but not reliable (p = .20).
  * **T3 fails:** −0.007 in the more monocular half.
  * Early A − P −0.083 (n.s.).
  * Normalization is present within neurons; onset summation is shown in MUA but not resolved in single units.
* **U holds for every claim** (penetrations as observations, 22 sessions): U1 −0.041 (p = .0004); U2 slope −0.89
  per doubling (p < 10⁻⁵); U3 −1.55 ms (p = .013); U4 ρ −0.48 (p < 10⁻⁵); U5 +0.051 (p = .0006); U6 w − a 0.29,
  0.21, 0.11 above / at / below (p ≤ .034). Monkey B differs in U2 (shallower contrast slope, +0.79, p = .0001) and
  U4 (weaker rank correlation, +0.19, p = .0005), not in U1, U3, U5 or U6.
* **For the paper:**
  * The unequal-contrast time course is the strongest new result: suppression graded with the other eye's
    contrast, absent at onset, sustained when the other eye is stronger. It is confirmed in held-out recordings for
    the sustained part only.
  * The laminar difference is robust to the depth reference.
  * Single neurons normalize; onset summation within neurons remains unshown.
  * Report the hierarchical estimates with monkey effects in the main text and the per-penetration counts in the
    legends.

## Novelty statement (draft for the Introduction and Discussion, 2026-11-03; Blake to check)

**Prior work, and what each showed.**
* **Cox et al. 2019** (J Vis): transient binocular enhancement, absent in the sustained response; one contrast pair.
* **Dougherty et al. 2019:** facilitation at onset, then suppression ~12 ms later, in monocular L4 units.
* **Mitchell et al. 2022** (iScience; a separate dataset): sublinear summation at 50–100 ms, then contrast-dependent
  suppression by 150–250 ms; equal contrasts, 100 ms windows.
* **Mitchell et al. 2023** (Curr Biol; **these recordings**, B + J, 25 penetrations): ocular dominance shapes how the
  non-dominant eye's contrast is combined; time-averaged responses; models weighted by ocular dominance.
* **The 2024 dissertation** (public; "Neuronal Basis of V1 Binocular Normalization"):
  * "a single, dynamical process of divisive normalization operating both within and between the eyes";
  * this process "was weakest in the geniculate input layers and intensified as signals propagated through the
    microcircuit".
* **SfN 2023:** "normalization across condition/depth/time" (docs/03 item 7; content to be supplied by Blake).

**What this paper adds** (each with its check):
1. **The full dichoptic design resolved in time, on each site's own response clock (Q, R; Fig 2).** In the first
   ~20 ms the eyes' contrasts add in every condition. Afterwards the binocular response falls short of summation in
   proportion to the other eye's contrast, and the shortfall is sustained when the other eye is the stronger one.
   The sustained part replicates in held-out recordings. Earlier work had equal contrasts or time-averaged
   responses.
2. **Onset latency follows drive across all 35 conditions, and the second eye counts as drive (G, H; Fig 1f).** It
   shortens the onset by 1.5–1.9 ms at every contrast and depth. This is the first laminar measurement in macaque;
   the human VEP shows ~3 ms, and a macaque single-unit study found no significant advance.
3. **The laminar profile of the shortfall (L, S; Fig 3).** It is weaker below the earliest sink than at and above
   it, for every depth reference tested.
   * **To reconcile with the dissertation:** it reports normalization weakest in the input layers. Here normalization
     is as strong at the earliest sink as above it, and weakest below. The difference may lie in the depth reference
     (the sink versus the L4/5 landmark) or in the measure (an index versus the response-unit shortfall). Blake to
     say which measure and alignment the dissertation used.
4. **The dynamics of the pool weighting (K, N, P; Fig 4).**
   * A delayed pool that weights the other eye more than the drive does, at every depth. This is the dynamic version
     of the 2023 ocular-dominance-weighted result.
   * *(Revised on v3, check Y, 2026-11-11.)* A two-stage form improves on it slightly at every depth, and recurrent
     normalization does not.
   * The extra weight is no longer reliably smallest deep (p = .076).
   * The pool delay is similar at the input and deep (44–46 ms). The shorter superficial delay is unconfirmed.
   * The model fits neither the laminar difference nor the late return.
5. **In supplement:**
   * single neurons normalize (T2);
   * the input current combines the eyes' contrasts, not their currents (O);
   * analog MUA leads isolated neurons by ~4 ms at every depth (M).

**What it does not claim:**
* onset summation within single neurons (T1 and T3 failed);
* a ±1 ms equivalence between the binocular advance and contrast summation (G2 failed);
* interocular suppression by a non-driving eye (I failed);
* a current-before-spikes order at the sink (J failed).

## V4 check V: the additivity index as the displayed measure (registered internally 2026-11-03, before running)

**Why.** Blake (2026-11-03): the "what the other eye adds" framing is a transformation that needs long labels; "this
is why additivity index is so powerful"; figures should be simple, one point each (J Neurosci). The figures will show
the additivity index AI = both eyes / (dominant eye + other eye), with a reference line for what contrast summation
predicts in the same units, AI_sum = R_D(c_D + a·c_O) / (dominant + other). AI < AI_sum has the same sign as the
registered shortfall per observation (the denominator is positive), but the index is a ratio of window means, so
its tests are run anew. **Not blind** to F, L, Q (same data, response units).

**V. Measure.** Window means after each channel's response onset (0–20, 30–100, 130–200 ms):
* per penetration from `population.npz` (channel means, the penetration's median eye weight);
* per penetration and depth group from `spacetime.npz` (channel means within the group, each channel's own eye
  weight for AI_sum).
AI = mean(B) / mean(D + O), AI_sum = mean(R_D(c_D + a·c_O)) / mean(D + O), where D + O > 0.12 of the dominant-eye
peak.
* **V1:** matched 45%, 30–100 ms, AI < AI_sum.
* **V2:** depth, (AI − AI_sum) below − at > 0 (45%, 30–100 ms).
* **V3:** 130–200 ms, cells with the other eye stronger (mean over cells per penetration), AI < AI_sum.
* **V4, reported:** 0–20 ms, AI − AI_sum (two-sided; summation at onset predicts ≈ 0 or above).
* **Criteria:** one-sided Wilcoxon over penetrations, p < .05. The halves were checked in the response-unit versions
  (F, L, Q). Script `scripts/v4_ai_display.py`, output [qc/v4_V.csv](qc/v4_V.csv).

**V results (2026-11-03). V1–V3 hold.**
* **V1:** matched 45%, 30–100 ms: AI 0.54 against 0.62 if the contrasts add (23/25, p < 10⁻⁴).
* **V2:** depth: AI − AI_sum below − at +0.09 (21/24, p < 10⁻⁴). By group: above −0.135 (20/21), at −0.113 (22/24),
  below −0.020 (18/24, p = .018).
* **V3:** other eye stronger, 130–200 ms: 0.59 against 0.65 (24/25, p < 10⁻⁵).
* **V4, reported:** at onset, 0–20 ms, matched 45%: 0.73 against 0.71 (n.s.).
* The figures may therefore show the additivity index with the summation reference line; the statistics on the
  panels are V's.

### Simplified figure set (Blake, 2026-11-03: one point per figure, the additivity index, J Neurosci)

| Figure | The one point | Panels |
|---|---|---|
| 1 | The experiment | stereoscope and probe, trial, 6 × 6 design with swatches, example penetration (raw traces). The contrast responses and ocular dominance go to the supplement |
| 2 | Response onset follows drive, and the second eye is drive, at every depth | (a) onset vs response strength, 35 conditions; (b) laminar onset, 0–50 ms: MUA by depth for both eyes and the dominant eye alone, with each one's onset at each depth |
| 3 | Summation gives way to normalization within ~30 ms, more at and above the input layer, and more when the other eye is stronger | (a) additivity-index ridgeline by depth (45%), with the summation reference; (b) the index by depth group, 30–100 ms, data vs reference; (c, d) the index against the other eye's contrast at 0–20 and 130–200 ms |
| 4 | A delayed pool, weighted toward the other eye | (a) sketch; (b) held-out gain of the delay by depth; (c) the other eye's drive weight vs pool weight by depth |

The richer figures built earlier today (`fig2_full_design.py`, `fig3_laminar_normalization.py`,
`fig4_model_spacetime.py`, S7, S8) become supplementary candidates.

### Revised figure set after the strength review (2026-11-03)

| Figure | Claim | Script |
|---|---|---|
| 1 | The experiment; onset follows drive across all 35 conditions (1f: trend band only over the monocular range) | `fig1_setup.py` |
| 2 | **The full 6 × 6 design over time**: matched vs unequal contrasts; the shortfall against the other eye's contrast in three windows; the sustained shortfall per monkey and in the held-out recordings | `fig2_full_design.py` |
| 3 | Space and time: what the other eye adds on the stimulus and response clocks (negative side faint, `figstyle.ADDED_SOFT`); ridgeline at 45%; strength by depth group with the four depth-reference variants | `fig3_laminar_normalization.py` |
| 4 | The model across depth: sketch with the fitted exponent (n ≈ 1.1); data / no delay / conductance maps; residuals; **d**, each model vs the fixed-delay pool (leads); **e**, the delay itself for comparison; weights, delay and data vs model by depth | `fig4_model_spacetime.py` |
| S7 | The synaptic current at the input (J, M, O) | `figS7_input_current.py` |
| S8 | The CSD depth landmark with the MUA onset by depth; shortfall maps at three contrasts; descriptive normalization time | `figS8_spacetime_maps.py` |

The linearity index over time (formerly Fig 2a) stays in S2a. `fig2_spacetime.py` is retired: its maps moved to
Fig 3, its CSD landmark to S8. The V4-era `fig2_add_then_normalize.py`, `fig3_contrast_depth.py` and `fig4_model.py`
stay for the imports and until the set is approved.

### Space-time figure set: status (2026-11-03)

**Decided (Blake, 2026-11-03):** the CSD analyses go to the supplement (Supplementary Fig S7: current additivity,
contrast summation in the current, sink vs MUA vs unit timing); the CSD map stays in Fig 2 only as the landmark that
defines the depth axis. The two-stage result is reported in Fig 4e and one Results sentence, not in the title or
abstract (gains 1–3%, best-model counts descriptive, no direct test against conductance).

| Figure | Claim | Status | Script |
|---|---|---|---|
| 1 | The experiment; onset follows drive | **built** (a–e as before; f, onset vs response strength for all 35 conditions, check H) | `fig1_setup.py` |
| 2 | Space-time summation | **built** (a, linearity index over time at 11, 22.5, 45%; b, CSD with the spiking onset at each depth; c, d, what the other eye adds at 45% on the stimulus and response clocks) | `fig2_spacetime.py` |
| S7 | The synaptic current at the input | **built** (checks J, M, O) | `figS7_input_current.py` |
| 3 | Laminar organization of normalization | **built** (maps at three contrasts, ridgeline at 45%, strength by depth group, normalization time) | `fig3_laminar_normalization.py` |
| 4 | Delayed normalization across time and depth | **built** (sketch; data / no delay / conductance maps; residuals; held-out gains, pool weight and delay by depth; data vs model normalization by depth) | `fig4_model_spacetime.py` |

The V4-era `fig2_add_then_normalize.py`, `fig3_contrast_depth.py` and `fig4_model.py` stay until the set is
approved (S2 imports from `fig3_contrast_depth.py`; the new Fig 4 imports the sketch from `fig4_model.py`).

## Laminar audit W (registered internally 2026-11-03, before running; Blake: "Yes … align every session based on where the PPC flips")

**Why.** The depth results are not clean, and the input markers disagree.
* **The earliest MUA onset and the earliest sink agree in few penetrations.** Their depths are within ±2 contacts in
  6/26 penetrations; the median offset is +4 contacts, with the MUA deeper.
* **The sink and the L4/5 call agree in 6/22.**
* **The earliest-MUA depth is unreliable.** It is reproducible within 2 contacts across split halves in 14/24
  penetrations at 90%, and in 9/24 at 45%.
* **The in-cortex spans are too long for one sheet.** They run 20–32 contacts (2.0–3.2 mm), longer than V1 (~1.8 mm).
* **Blake:** a deep probe crossing a second fold of cortex "would explain why our deep layers appear to have
  response latencies as fast as the input layer above it."
* **Blake asked for three things:**
  * bring back the spectrolaminar motif and the SPI;
  * align every session on the preferred-phase ("PPC") flip as a programmatic depth;
  * check that the SPI descends with depth, and that alpha/beta power is higher in the deeper layers.

**Blind.** Every call uses laminar markers or monocular responses only, never a binocular response or a
combination measure:
* from `output/laminar/{pen}_laminar.npz` (30 min of the recording): preferred phase, SPI and PPC, spontaneous
  threshold-crossing rate and relative power;
* clean monocular presentations, one eye at 45 or 90% and the other eye at 0%, for the grating CSD and MUA.

The cards may be inspected and the rules revised before Blake reviews them. Nothing binocular is computed on the new
alignment until Blake has made the calls and the rerun is registered (step 3 below).

**W. Draft rules** (`scripts/laminar_audit.py`; one card per penetration; table [qc/laminar_audit.csv](qc/laminar_audit.csv)):
1. **Orientation.** Majority of three votes: the Davis et al. 2023 preferred-phase template fit; the
   spectrolaminar peak order (gamma above alpha/beta); and the SPI falling with contact order. The SPI here is each
   spike channel's mean over LFP contacts, using channels with ≥ 500 spikes and above the rate floor. The vote
   count is reported, together with the legacy call.
2. **Top of cortex.** The existing rule: the first of 3 consecutive contacts from the superficial end with a
   spontaneous rate above 2.5 Hz (evoked MUA as the fallback).
3. **Phase reversal: the primary depth reference.** The Davis template is slid along the preferred-phase profile
   with the orientation fixed, and the reversal is the template boundary (the input/deep flip). Local phase steps
   are shown on the card as a cross-check.
4. **Bottom of the first sheet.** Below the reversal, the first run of ≥ 2 contacts at or below the rate floor
   (white matter); the first sheet ends at the contact above it. A **second sheet** is ≥ 3 consecutive active
   contacts below that gap, or a reversed phase step (|step| ≥ 0.8 rad) below the reversal.
5. **Input layer, within the first sheet.**
   * The earliest monocular sink: the earliest sink reaching 25% of its own peak, among sinks ≥ 40% of the
     strongest (the dichoptic-CSD rule).
   * The earliest MUA onset: the half-max onset of each contact's monocular MUA (both eyes pooled), where the
     evoked t > 6, taking the minimum of a 3-contact running median.
6. **The checks Blake asked for, within the first sheet:**
   * the SPI falls with depth below the top of cortex (Spearman ρ < 0; PPC reported);
   * alpha/beta minus gamma relative power rises with depth (ρ > 0);
   * the depth of the spectrolaminar crossing relative to the reversal is reported.

**Expected if the rules work.**
* Top of cortex to the reversal ~1.0–1.3 mm (L1 to the bottom of L4C).
* The reversal to the bottom of the sheet ~0.5–0.8 mm.
* The first sheet ≤ ~2.0 mm.
* The earliest sink and the earliest MUA 0–0.5 mm above the reversal (L4C), and within 0.2 mm of each other.
* The SPI and alpha/beta checks in the expected direction in most penetrations.

Each is counted, not tested; the audit is descriptive.

**Then (registered when Blake has made the calls, before running).** Once Blake has reviewed the cards and fixed
the calls (`metadata/laminar_alignment_v3.csv`), the following are rerun once on that alignment, restricted to the
first sheet:
* the laminar onset order (G4);
* the depth tests (V2 and L2 with the depth groups, S);
* the ridgeline.

Depth is measured from the phase reversal and, as a sensitivity check, from the audited earliest sink. The new
groups and predictions are registered before any of them is computed.

**Implementation notes, made before Blake's review and still blind to binocular outcomes.**
* The datasets store contacts from superficial to deep in the v2 orientation, so they are re-indexed to raw contact
  order before the audit orientation is applied.
* The sink onset is found walking back from each sink's peak. The first-crossing version locked onto noise at the
  20 ms window edge where the CSD was weak.
* Two rules use the monocular MUA amplitude (3-contact median; 1 = the probe's 90th percentile), because the
  spontaneous rate is too low and patchy on several probes for a 2.5 Hz floor:
  * the top of cortex is also called from the evoked MUA (≥ 0.2 of that amplitude and t ≥ 6 on 3 contacts), and
    both tops are shown;
  * the white-matter dip is < 0.35 of that amplitude on 2 contacts.
* A second early MUA zone must lie > 0.2 mm below the first sheet, start within 3 ms of the earliest onset, and be
  separated from it by contacts ≥ 5 ms later.

**W results (2026-11-03; descriptive; cards in `<out_root>/laminar_audit/`, summary `laminar_audit_summary.png`).**
* **Blake's checks.**
  * Alpha/beta − gamma power rises with depth in the first sheet in 24/26 penetrations (13/14 with a clear
    reversal). The exceptions are 221128_J_eA and 221214_J_eA.
  * The SPI falls with depth in 18/26 (10/14). It does not fall in 211012_B_eA, 211217_B_eA, 221202_J_eA, 221202_J_eB,
    221206_J_eA, 221214_J_eA, 221214_J_eB or 221215_J_eA.
* **Orientation.** The three votes agree in 16/26 penetrations; 10 go 2-of-3. The call differs from the laminar
  report's in 2 (221204_J_eB, 221215_J_eA) and from v2 (which built the datasets) in 1 (221206_J_eB). *Corrected
  2026-11-03: this line first said "differs from v2"; the card's "legacy" field is the laminar-report call.*
* **Phase reversal.**
  * Clear (≥ 1 rad, on the probe) in 14/26.
  * Small (0.5–0.8 rad) in 6: 221212_J_eA, 211008_B_eA, 221202_J_eB, 221206_J_eA, 221209_J_eB, 221214_J_eB.
  * Off the probe (no flip) in 6: 211210_B_eA, 221123_J_eA, 221204_J_eB, 221209_J_eA, 221214_J_eA, 221206_J_eB.
* **With a clear reversal** (14 penetrations; medians, IQR):
  * top of cortex 1.25 mm above the reversal (0.8–1.45), as expected;
  * earliest monocular sink 0.30 mm above (0.25–0.55); between 0.05 and 0.95 mm above in 11/14, so L4C sits above the
    input/deep flip;
  * earliest MUA at the reversal (+0.05 mm; −0.32 to +0.22), deeper than the sink by a median 0.15 mm; within
    0.2 mm of the sink in 8/14;
  * spectrolaminar crossing 0.11 mm above (0.74 above to 0.11 below).
  * The three exceptions: 211025_B_eA and 221215_J_eA, with the sink 0.65–0.75 mm below the reversal; and
    221109_J_eA, with the reversal at contact 2.5, near the deep end, and the sink 2.45 mm above it, so its flip
    is probably not L4/5.
* **Bottom of the first sheet: the rule does not work.** The dip lies a median 0.30 mm below the reversal
  (0.12–0.68), shallower than V1's deep layers (expected 0.5–0.8 mm). In 221130_J_eA the dip sits right under the
  reversal and strong MUA returns below it, more like layer 5 than white matter. The MUA-amplitude dip is not a
  white-matter marker on its own; the bottom needs Blake's eye on the cards.
* **Second sheet: little evidence.**
  * A separate early MUA zone appears in 2 penetrations.
    * 211217_B_eA: 1.75 mm below the reversal, at 33.9 ms, against 53.7 ms in the upper sheet. This is a candidate
      second fold.
    * 221130_J_eA: 0.35 mm below, at 30.0 vs 29.4 ms. More likely its own deep layers.
  * Activity returns after a dip in 6 (221212_J_eA, 211103_B_eA, 211208_B_eA, 221118_J_eA, 221119_J_eA,
    221130_J_eA).
  * A reversed phase step below the reversal appears in 1 (221215_J_eA).
  * So a second fold does not generally explain the fast deep onsets. Where the reversal is clear, the earliest
    spiking sits at the input/deep flip itself, 0.15 mm below the earliest sink.
* **The two tops disagree by > 0.3 mm in 11/26;** the spontaneous-rate top is the deeper one in all 11.

**For Blake to decide (cards and summary):**
* whether the phase reversal is the depth reference for the 14 clear penetrations;
* what to do with the 12 without a clear flip:
  * the earliest sink plus 0.3 mm,
  * hand calls,
  * or leaving them out of the depth tests;
* which top-of-cortex rule to use;
* the bottom of V1 and any second sheet, by eye (candidate: 211217_B_eA);
* 221109_J_eA's reversal.

## Laminar alignment v3 and the depth rerun X (registered internally 2026-11-03, before running; Blake: "Okay yes. Go ahead … On the ones that need my eye, I would just make a decision that you think is best.")

**Why.** W showed that the phase flip alone gives a depth reference in 14/26 penetrations and is wrong in 2 of them.
Blake accepted the recommendations: a cross-checked flip as the anchor, the evoked-MUA top of cortex, a fixed deep
cap instead of a white-matter call, and the rule's outcome on the cards flagged for his eye. **Not blind:** the depth
results L, V2, G4 and S on the sink-based groups are known (all hold). Nothing binocular has been computed on v3.

**Alignment v3** (`scripts/freeze_laminar_alignment_v3.py` → `metadata/laminar_alignment_v3.csv`; laminar markers
only, frozen before this registration):
* **Three estimates of the input/deep boundary:**
  * the phase flip, if on the probe with a step ≥ 0.5 rad;
  * the earliest monocular sink + 0.3 mm;
  * the alpha/beta–gamma crossing + 0.1 mm.

  The offsets are W's medians over the 14 clear flips. The earliest-MUA depth is not used: it is noisy, and using
  it would make the onset-by-depth tests circular.
* **Anchor:**
  * the flip, if the sink or the crossing lies within 0.4 mm of it;
  * otherwise the mean of the sink and the crossing, if they lie within 0.4 mm of each other;
  * otherwise none, and the penetration is left out of the depth analyses.

  The anchor is rounded to the nearest boundary between contacts.
* **Top of cortex:** the evoked-MUA top. With it, the top lies 1.25 mm above the anchor (IQR 0.80–1.52), as expected
  for the pia to the bottom of layer 4C. With the spontaneous-rate top the distance is 0.85 mm, too short.
* **Deep cap:** contacts more than 0.7 mm below the anchor are excluded, which covers layers 5 and 6 plus slack for
  oblique penetrations. There is no white-matter call; any second sheet falls beyond the cap.
* **Depth groups:**
  * superficial: from the top of cortex to 0.5 mm above the anchor;
  * input: 0.5 mm above the anchor to the anchor;
  * deep: the anchor to 0.7 mm below.

  With the anchor between contacts, the input group has 5 contacts and the deep group 7.
* **Freeze result:**
  * 22/26 penetrations anchored: 15 on the flip (11 clear flips and 4 small steps confirmed by the sink or
    crossing), 7 on the sink and crossing.
  * Left out: 211210_B_eA, 221118_J_eA, 221204_J_eB, 221206_J_eB.
  * The orientation matches v2 in all 22 (it differs only in 221206_J_eB, which is left out), so the datasets'
    channel order is unaffected. Depth is computed from contact numbers (`dataset.anchor_depth`).
* **Cards flagged for Blake's eye** (delegated; each keeps the rule's outcome):
  * 221215_J_eA: orientation as v2; flip, crossing and earliest MUA agree. The sink (0.75 mm below) and the
    reversed step lie beyond the cap.
  * 221214_J_eA: the sink and crossing agree with an on-probe phase step (π → ~0 rad) that the template fit
    missed. Its alpha/beta inversion comes from contacts beyond the cap.
  * 221118_J_eA: left out; the flip, sink and crossing spread over 1.1 mm.
  * 211217_B_eA: no second sheet. The deep early zone is low-amplitude and patchy, with no reversed phase step,
    and it lies beyond the cap.
  * 211025_B_eA and 221109_J_eA, anchored on the sink and crossing: the flip lies 0.25 mm below the top of cortex
    (211025), or in unresponsive tissue at the deep end of the probe (221109).

**Channels.** MUA channels in the first sheet (below the evoked-MUA top, no deeper than the cap) that are driven by
the `run_crf.py` rule. That rule is a one-way ANOVA across the six monocular contrasts, p < .01, with a positive
peak, in either eye and either window. It is recomputed here for every contact, because v2 placed some superficial
contacts above the cortex.

**Measures, as registered for L, V and G** (per channel, clean presentations, condition means blank subtracted,
Gaussian σ = 2 ms, dominant eye and scale from all presentations):
* **L:** added A and predicted P at 45%, 30–100 ms after the channel's onset; all, odd and even presentations.
* **V:** AI and AI_sum per penetration and group, as the ratio of channel means; all presentations.
* **G:** the dominant eye's interpolated half-max onset at 90%, and the binocular advance at 45%. Per channel, then
  the median per penetration and group.

**Tests** (one-sided Wilcoxon over penetrations; a test holds if p < .05 and, where halves exist, both halves'
medians lie in the predicted direction):
* **X1 (L2 on v3):** (A − P)deep − (A − P)input > 0. Reported: superficial − input (two-sided), and L1 per group.
* **X2 (V2 on v3):** (AI − AI_sum)deep − (AI − AI_sum)input > 0.
* **X3 (G4 on v3):** the dominant eye's onset at 90% is later in the superficial group than in the input group, and
  later in the deep group than in the input group (each > 0). The deep comparison is the less certain of the two:
  layer 6 receives geniculate input, and W found the earliest MUA at the flip.
* **X4 (G3 on v3):** the binocular advance at 45% is < 0 in each group.
* **X5 (S1 on v3):** per-channel shortfall A − P ~ depth from the anchor (mm), with a penetration random intercept;
  slope > 0.
* **Sensitivity of X1 and X2:**
  * (a) deep cap 0.5 mm;
  * (b) the sink anchor (input group = the 5 contacts centred on the earliest monocular sink) in the same 22
    penetrations;
  * (c) only the 15 penetrations anchored on the flip.

  "Robust to the depth reference" only if X1 holds in all three.
* **Ridgeline (display only):** Fig 3a redrawn on v3 depth, in 0.2 mm bins from the anchor, with ≥ 10 penetrations
  per bin.

**Commitment.** v3 is the alignment for the paper whichever way X comes out. The figures and text take X's
numbers. If X1 fails, the paper reports normalization at every depth and drops "weaker below the input layer";
the sink-based result moves to the supplement with its dependence on the depth reference stated. Script:
`scripts/x_depth_v3.py`; outputs `docs/qc/x_depth_v3.csv`, `x_groups.csv`.

**X results (2026-11-03). X1, X2, X4 and X5 hold and are robust to the depth reference; X3 holds for the
superficial group and fails for the deep group.**
* **Data:** 22 penetrations; 424 first-sheet channels.
  * 92 of these are superficial channels that v2's spontaneous-rate top had placed above the cortex.
  * For the 341 channels shared with L, the per-channel values are identical to L's (differences < 10⁻¹⁵), and all
    341 were also driven in `run_crf.py`.
  * Group medians per penetration: superficial 8 channels, input 5, deep 7.
  * Group values in [qc/x_groups.csv](qc/x_groups.csv).
* **X1 holds:** (A − P) deep − input = +0.052 (18/22, p < 10⁻⁴; halves +0.040, +0.052).
  * L1 by group: superficial −0.078 (20/21, p = .0001), input −0.073 (20/22, p = .0001), deep −0.029 (15/22,
    p = .053, so it fails below the boundary).
  * Superficial − input: −0.006 (p = .52).
* **X2 holds:** (AI − AI_sum) deep − input = +0.11 (18/21, p = .0001).
  * The index from 30 to 100 ms, against summation: superficial 0.51 vs 0.67, input 0.49 vs 0.66, deep 0.55 vs
    0.60. Deep is still below summation (15/21, p = .010).
* **X3:**
  * **Superficial holds:** the superficial onset is 5.3 ms later than the input onset (18/21, p = .0001; halves
    5.7, 5.2).
  * **Deep fails:** deep − input = −1.35 ms (8/22 later; one-sided p = .88). The deep group starts no later than
    the input group, and earlier in 14/22 (two-sided p = .26, descriptive).
  * **By 0.2 mm bin** (median over penetrations, relative to each penetration's median onset; descriptive): the
    earliest onsets lie 0.1–0.6 mm below the anchor (−1.8 to −2.9 ms). Onsets rise to 0 ms at 0.7 mm below and
    to +3 to +8 ms above 0.8 mm above the anchor.
  * **Why G4 differed:** G4's "below the sink" group (+1.4 ms later) had no deep limit, so it included the later
    contacts beyond the cap.
  * **Reading (not tested):** fast deep onsets are real and sit in the 0.6 mm below the input/deep boundary. That
    fits direct geniculate input to layer 6, possibly with thalamic axon spikes in the MUA. It does not point to
    a second fold or a misalignment.
* **X4 holds in every group:** the binocular advance at 45% is −2.0 ms superficial (20/21), −1.5 ms input (20/22)
  and −1.1 ms deep (17/22, p = .005), with both halves in every group.
* **X5 holds:** the shortfall rises by +0.030 per mm of depth from the anchor (p < 10⁻⁴; 424 channels).
* **Sensitivity: X1 and X2 hold in all three variants**, so "robust to the depth reference" is allowed.

  | Variant | X1 | X2 |
  |---|---|---|
  | (a) 0.5 mm cap | +0.041 (17/22, p = .0001) | +0.090 (p = .0001) |
  | (b) sink anchor | +0.031 (17/22, p = .003) | +0.067 (p = .0003) |
  | (c) 15 flip-anchored penetrations | +0.035 (11/15, p = .006) | +0.061 (p = .003) |
* **Ridgeline:** Fig 3a and 3b redrawn on v3 (`fig3_normalization_ai.py`; the sink version is kept as
  `fig3_normalization_ai_sink`). From 30 to 100 ms the index falls furthest below summation 0.1–0.7 mm above the
  boundary:

  | Depth from the anchor | 0.7 above | 0.5 above | 0.3 above | 0.1 above | 0.1 below | 0.3 below | 0.5 and 0.7 below |
  |---|---|---|---|---|---|---|---|
  | Shortfall below summation | 0.14 | 0.17 | 0.14 | 0.11 | 0.07 | 0.04 | 0 |

  Superficial bins fall 0.06–0.09 short.
* **For the paper:** "Normalization was strongest in and above the input layer and faded below the input/deep
  boundary. From 30 to 100 ms after onset, the additivity index fell 0.17 below the summation prediction within
  0.5 mm above the boundary, and 0.05 below it in the 0.7 mm beneath (18 of 21 penetrations). The result held for
  three alternative depth references. Responses began 5 ms later in the superficial layers than in the input
  layer; deep sites began no later than the input layer. Adding the other eye advanced the onset at every depth."
* **Still on the sink groups** (rerun on v3 before they enter the figures): the model by depth (K3, N), U5 in the
  hierarchical model, and the laminar latency panel planned for Fig 1. *(Done 2026-11-11: check Y below.)*

## Check Y: the remaining depth results on alignment v3 (registered internally 2026-11-11, before running; Blake: "Let's rerun the remaining depth analyses on v3. Then, let's re-evaluate.")

**Why.** X moved L, V2, G4 and S to alignment v3. The model by depth (K, N, P), U5 and U6 of the hierarchical
model, Fig 2b (onset by depth) and Fig 4b–c still use the groups above / at / below the earliest sink. **Not blind:**
N, K and P on the sink groups are known (delay and w > a at every depth, the excess weight and the delay gain
smallest below the sink, no depth difference in τ2, two-stage best at and below the sink). So are X's groups
(normalization weaker in the deep group). Nothing model-based has been computed on v3.

**Groups and channels.** X's channels and groups (`x_channels.parquet`: first sheet, driven by the recomputed
`run_crf.py` rule; superficial / input / deep from the v3 anchor). Group fits need ≥ 3 channels in a group with a
positive dominant-eye peak, as N.

**Y1 (N on v3).** `v4_model_groups.py`'s group populations, models ('static', 'full', 'w=a', 'cond'), starts and
held-out errors, on the v3 groups. N1–N3 per group with N's criterion (p < .05 one-sided on the summed held-out
error, and both directions' medians in the predicted sign); "at every depth" only if all three groups hold.
N4, reported two-sided: deep − input and superficial − input in τ2, w − a, δ, τ_m (conductance fits) and the
held-out gains.

**Y2 (P on v3).** The recurrent and two-stage models on the same groups; P1 (best-model counts, % less error than
'full'), P2 (org vs cond) and P3 (2stage vs full) per group, with P's criteria; the depth differences against the
input group (two-sided).

**Y3 (K on v3, no refit).** The existing per-channel fits (`v4_K_channels.parquet`) for channels in X's set,
grouped on v3, median per penetration and group. K1′ and K2 per group as registered; K3 deep − input and
superficial − input (two-sided) in τ2, w − a, δ and the delay gain. Reported only: X's 92 newly included
superficial channels have no per-channel fits, so the superficial group is smaller here than in Y1.

**Y4 (U on v3).** The U model (penetration as the observation, random intercept per session, monkey sum-coded).
* **U5′:** (A − P) deep − input > 0 (45%, 30–100 ms; `x_groups.csv`). Reported: superficial − input (two-sided).
* **U6′:** w − a > 0 in each v3 group (Y1 conductance fits).
* U1–U4 do not depend on depth and are not rerun.

**Y5 (display, no tests).**
* **Fig 2b on v3 depth:** 0.2 mm bins from the anchor, first sheet, ≥ 10 penetrations per bin. Onsets are X's
  interpolated per-channel onsets at 45%. The other eye alone is added to `x_depth_v3.py`, and that script now also
  saves the stimulus-clock arrays. X is rerun for these, and its test table must come out unchanged.
* **Fig 4a–c on Y1's fits and groups.**

**Wording allowed in the paper:**
* "a delayed pool at every depth" only if Y1 N1 holds in all three groups;
* "the pool weights the other eye more than the drive at every depth" only if U6′ holds in all three;
* "the extra weight is smallest below the input/deep boundary" only if Y1 N4 (w − a, deep − input, conductance
  fits) is < 0 at two-sided p < .05;
* the pool delay by depth only as "no difference detected", with the medians;
* "depth X behaves like model Y" only for a P test that holds in Y2.

**Commitment.** v3 numbers replace the sink-group numbers in the figures and text whichever way Y comes out. The
sink-group results stay here as run. Script: `scripts/y_model_v3.py`; outputs `docs/qc/y_N.csv`, `y_P.csv`,
`y_K.csv`, `y_U.csv`.

**Y results (2026-11-11). Y1 N1–N3, U5′, U6′ and P3 hold in every group. The depth difference in the extra pool
weight no longer reaches significance on v3, and a shorter superficial pool delay appears (reported, not registered
as a claim).**
* **Data:**
  * Y1/Y2: 21 penetrations, 57 groups (superficial 16, input 20, deep 21). 211008_B has too few presentations per
    cell in a half, as in N.
  * Y3: 20 penetrations, 293 channels with per-channel fits (superficial 69, input 95, deep 129).
  * X's rerun with the added outputs reproduced `x_depth_v3.csv` and `x_groups.csv` byte for byte.
  * The refactored test functions reproduce `v4_N.csv`, `v4_P.csv` and `v4_K.csv` exactly from the saved fits.
* **Y1 (N on v3):** N1, N2 and N3 hold in all three groups, in both cross-validation directions.
  * **N1, the delay:** log(CV static / CV full) 0.69 superficial (16/16), 1.27 input (19/20), 0.97 deep (20/21);
    all p < 10⁻⁴.
  * **N2, pool weighting:** 0.070 (14/16), 0.037 (18/20), 0.010 (18/21); p ≤ .0003.
  * **N3, conductance:** 0.021 (14/16), 0.005 (14/20), 0.0002 (13/21, p = .018). In the deep group the gain is
    negligible: both directions' medians are ≤ 0.0001.
  * **N4, reported (two-sided, against the input group):**
    * **Excess pool weight:** w − a is 0.23 / 0.28 / 0.13 (superficial / input / deep). Deep − input is −0.127
      (6/20), p = .076, so it does not reach p < .05 (on the sink groups it was p = .003).
    * **Pool delay:** τ2 is 27 / 46 / 44 ms. Superficial − input is −20.5 ms (2/16 longer, p = .004); deep − input
      is +0.3 ms (p = .96). The fixed-delay fits agree (16 / 36 / 37 ms, p = .002). The conductance time constant
      τ_m is longer superficially (18 vs 7.5 ms, p = .006).
    * **Delay gain:** smaller superficially (0.69 vs 1.27, p = .0006) and slightly smaller deep (p = .040).
    * **Input delay δ:** no depth difference (19.4 / 20.4 / 18.0 ms).
* **Y2 (P on v3):**
  * **P3 holds in every group:** two-stage beats the fixed delay superficially (+0.046, 12/16, p = .009), at the
    input (+0.034, 16/20, p = .002) and deep (+0.015, 17/21, p = .0004), in both directions. On the sink groups it
    failed above the sink.
  * **P2 fails in every group:** recurrent normalization never beats the conductance form (superficial +0.019,
    p = .30; input −0.089; deep −0.065). It is relatively less bad superficially than at the input (p = .001).
  * **P1, best model (counts):**
    * superficial: recurrent 7, two-stage 7, conductance 2;
    * input: two-stage 14, conductance 4, fixed delay 1, w = a 1;
    * deep: two-stage 12, conductance 6, recurrent 2, static 1.
  * **% less error than the fixed delay:** two-stage 4.5 / 3.4 / 1.5; conductance 2.0 / 0.5 / 0.0; recurrent
    2.5 / −8.9 / −5.2.
* **Y3 (K on v3, no refit):**
  * **K1′ and K2 hold in every group:** delay gain 0.71 / 0.88 / 0.59, both directions; w − a 0.29 / 0.30 / 0.15;
    all p ≤ .0001.
  * **K3:**
    * w − a deep − input −0.09 (p = .13);
    * the delay gain is smaller deep (−0.29, 1/20, p < 10⁻⁴);
    * τ2 superficial − input −4.8 ms (3/16, p = .013; medians 22 / 37 / 35 ms);
    * δ deep − input −0.9 ms (p = .044).
* **Y4 (U on v3):**
  * **U5′ holds:** deep − input +0.048 (22 penetrations, 19 sessions, p = .0004; monkey B effect −0.009, p = .52).
    Superficial − input −0.011 (p = .52, two-sided).
  * **U6′ holds in every group:** w − a 0.18 superficial, 0.27 input, 0.19 deep (p ≤ .0001). Monkey B is lower
    superficially (−0.11, p = .009), not at the input (p = .16) or deep (p = .38).
* **Y5 (display):**
  * **Fig 2b on v3:** the earliest onsets are 0.1–0.5 mm below the boundary (29–30 ms after the flash for both
    eyes at 45%). Onsets rise to 37–39 ms 0.9–1.3 mm above it. Both eyes lead the dominant eye alone in nearly
    every bin.
  * **Fig 4 on Y1:** the delayed pool reduces held-out error by 50 / 72 / 63% (median). The other eye's weight is
    0.77 / 0.86 / 0.86 in the drive and 1.01 / 1.06 / 1.05 in the pool.
  * The sink versions are kept as `fig2_latency_drive_sink` and `fig4_delayed_pool_sink`.
* **Reading the new superficial pool delay (not tested).**
  * The group fits run on the stimulus clock, over channel means. Superficial onsets are later and more spread
    out (X3), so averaging smears the superficial transient. That alone would shorten the fitted delay and lower
    the delay gain.
  * The per-channel fits, which have no such smearing, show a much smaller difference (−4.8 ms).
  * Treat the shorter superficial delay as unconfirmed. It needs a response-clock test before any claim.
* **For the paper (wording allowed by the registration):**
  * "At every depth, a delayed normalization pool predicted held-out responses better than no delay (50–72% less
    error), and weighted the other eye more than the drive did (w − a 0.18–0.27). A two-stage form, in which each
    eye is first gain-controlled by the other, improved the fit slightly further at every depth (1.5–4.5%).
    Recurrent normalization did not improve on the conductance form at any depth."
  * Dropped: "the extra weight is smallest below the input layer" (p = .076 on v3; the medians still point that
    way, 0.13 vs 0.28).
  * Dropped: "the pool delay did not differ with depth". Superficially it is shorter in the group fits, but see the
    reading above. The text says the delay was similar at the input and deep (44–46 ms), and states the superficial
    estimate with its caveat or leaves it to the supplement.
  * Unchanged: the fitted model does not reproduce the laminar difference of X1. That comparison (N, "what the
    model misses") is still on the sink groups. On v3 the pool weights show no reliable depth difference, so the
    gap stands. The paper's laminar claim rests on the data (X1, U5′), not on the model.

All depth results are now on alignment v3.

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
  *(Superseded in part by the rebuild below: Fig 2b and 2c left the main figure.)*

### Figure set for V4: rebuilt for readers outside the field (2026-11-02)

**Why.** Blake: "if you were someone who was not as familiar with the literature, would you be able to understand
what each figure is showing? … I think it's very important that we make this paper accessible." A cold read found
four barriers: no picture of the two predictions before the data (the summation prediction hid in a legend
parenthesis, "dominant eye, 90% (≈ linear sum)"); "linear sum" reads as "the responses add"; the late effect sat in
small gaps between overlapping traces; and measures were built on measures (equivalent contrast → linearity index →
its slope on contrast; iso-response contours; ratios of held-out errors). Rules in docs/11 ("Written for readers
outside the field").

| Figure | Claim | Panels | Script |
|---|---|---|---|
| 1 | The experiment | (a) stereoscope with each eye's grating, laminar probe in V1, trial sequence (schematics); (b) 6 × 6 design with a grating swatch at each contrast; (c) example penetration: both eyes, each eye alone, "if the contrasts add" (its own eye weight); (d) monocular contrast responses; (e) ocular dominance | `fig1_setup.py` |
| 2 | The eyes add, then normalize | (a) what to look for: sketches of "contrasts add" vs "the eyes normalize each other"; (b) population at 45%: both eyes, dominant eye alone, "if the contrasts add" (45% × (1 + a), each penetration's a); (c) what the other eye adds vs what summation predicts it adds; (d) early vs late linearity index per penetration (unchanged statistics) | `fig2_add_then_normalize.py` |
| 3 | Contrast sets how far, at every depth; the input current adds at onset | (a) what the other eye adds at 11, 22.5, 45%; (b) transition-window index vs contrast per penetration, per-monkey slope counts in the note (8/8, 18/18); (c) depth profiles; (d) population CSD; (e) current added at the earliest sink vs summation (the current's own eye weight, median 0.48) | `fig3_contrast_depth.py` |
| 4 | Normalization is delayed and weights the other eye more | (a) timing sketch: drive (both eyes add at once), pool (lagging, median fitted 30 ms), response; (b) data vs conductance fit, four conditions (keeps 90% for the text); (c) what the other eye adds at 11 and 45%: data vs no-delay (static) and delayed-pool (fixed-delay) fits; (d) held-out interaction error; (e) pool vs drive weight | `fig4_model.py` |
| S1 | Reading the index; laminar alignment | (a) equivalent-contrast readout (former Fig 2b); (b) example flash CSD | `figS1_methods.py` |
| S2 | The index is valid | (a) fixed-rule null; (b) slope on contrast by monkey and held-out runs (former Fig 3b); (c) cross-validated dominance; (d) ocular-dominance terciles; (e) units | `figS2_validity.py` |
| S3 | Iso-response contours | as before (the former Fig 2c lives here) | `figS3_iso_contours.py` |
| S5 | Model checks | (a)–(e) as before; (f) held-out comparison of the three dynamic models (former Fig 4f) | `figS5_model_checks.py` |

S4, S6 and Table S1 are unchanged (Table S1's panel references updated). Shared code: `v1lb.binocular.response_at`
and `other_eye_added`, `figstyle.added_panel`, `SUM`, `LINEARITY`, `smooth` (5 ms running mean, display only).

**Statistics.** Every statistic carried over is unchanged (printouts diffed against the saved ones: Fig 2d, the
slopes now in S2b, Fig 3c–e tests, Fig 4b/d/e, S5f, S1 readout and landmarks). New numbers, all descriptive except the
registered check F above:
* Fig 2c / 3a window means (added / predicted, dominant-eye peak = 1): 45% 0.199 / 0.164 (0–20 ms), −0.004 / 0.037
  (30–100), 0.030 / 0.041 (130–200); 11% 0.119 / 0.115, 0.042 / 0.061, 0.030 / 0.051; 22.5% 0.167 / 0.136,
  0.008 / 0.018, 0.015 / 0.014.
* Fig 3e: the sink current's eye weight is lower than the MUA's (median 0.48, IQR 0.20–0.91, vs 0.70), so its
  summation prediction is smaller than the old 90% stand-in; the measured added current exceeds it at onset.
* Fig 4c (stimulus clock, added: data / static / fixed-delay): 45% 30–60 ms 0.179 / 0.031 / 0.156, 70–140 ms
  0.007 / 0.029 / 0.009, 150–219 ms 0.022 / 0.029 / −0.002; 11% 0.088 / 0.056 / 0.097, 0.047 / 0.057 / 0.052,
  0.028 / 0.057 / 0.036. The static model adds a small constant amount from the start and misses the onset burst
  at every contrast; the delayed model gets the burst and the drop but not the late return at 45% (check F).

Carry into V4 (on top of the list above):
* **Figure references in the text:** Fig 2b → Supplementary Fig S1a; Fig 2c (contours) → Supplementary Fig S3; old
  Fig 3a → Fig 3b; old Fig 3b → Fig 3b note and Supplementary Fig S2b; old Fig 3e → Fig 3e (now "added" format); old
  Fig 4a model schematic → Fig 4a timing sketch and the Methods; old Fig 4c (index over time with models) → Fig 4c
  ("added" format); old Fig 4f → Supplementary Fig S5f; S2 b–d → S2 c–e. Fig 1 letters: example c, CRF d, ODI e.
* **Text from check F:** at 45% the other eye adds about what summation predicts in the first 20 ms, nothing from
  ~30 to ~100 ms, and its contribution partly returns by 130–200 ms, where it no longer differs reliably from
  summation (F3 failed); no undershoot (F4 failed). At 11% it stays below summation late (about 60%) without a
  return. Remove "the 45% undershoot" from the model paragraph and the Fig 4 legend; say the delayed model misses the
  late return at 45%.
* **Legends:** rewrite all four main legends for the new panels (illustrations labelled; "if the contrasts add" =
  45% × (1 + a) with each penetration's eye weight; "added" = both eyes − dominant eye alone; 5 ms running mean for
  display); S1, S2, S5 legends for the moved panels.
* **Words:** "linear sum" → "if the contrasts add" / "linear summation of contrast" on first use; define the
  linearity index as the share of the other eye's contrast added.
* **Blake to decide:** how the title, abstract and Fig 2 heading phrase "normalize later", given that the strong
  normalization at 45% is transient in response units and the sustained partial normalization is clearest at 11%.

Model colours (no model uses orange, the other eye's colour): data ink, static dark grey dashed, fixed-delay aqua,
conductance violet, recurrent magenta (aqua and magenta below 3:1 on white and close for deutans, so always
direct-labelled). Inside a panel of condition time courses (Fig 4b) a model is drawn as a dashed line in the
condition's colour.

### Proposed figure set: latency as the thread, space-time maps (2026-11-03; awaiting Blake)

**Why.** Blake, reading the paper (2026-11-03): Fig 1 should make the monocular vs binocular latency point and
show MUA and CSD latency by depth, to set this paper apart from Mitchell et al. 2022 and later work. "That's what
this paper adds: 1 ms resolution of binocular summation in V1 and laminar organization of normalization." Fig 2
should make the temporal point a spatiotemporal one: additivity maps of depth × time aligned to the response onset,
beside the same maps aligned to the stimulus onset; its first panel should show that latency shortens with contrast
and with adding the other eye's input ("drive also reduces response latency"). Data for the latency panels: check G
above. Maps mocked in the scratchpad (`figpolish/v5/fig2_mock.py`, `st_mock.py`); descriptive, no new test.

| Figure | Claim | Panels |
|---|---|---|
| 1 | The experiment; we can tell when each depth starts | (a)–(e) as now; (f) population CSD aligned on the earliest sink (moved from Fig 3d) with the MUA onset at each depth overlaid: the current arrives first at the sink, spiking starts there, 4 ms later above and 1.4 ms later below (G4) |
| 2 | Drive shortens latency, and the second eye is drive; the eyes add at every depth on each site's own clock | (a) onset vs contrast: each eye alone, both eyes, "if the contrasts add" (G1, G2); (b) the same at 45% by depth group (G3); (c, d) what the other eye adds, depth × time at 45%, on the stimulus clock (with each depth's onset drawn) and on each site's response clock |
| 3 | The eyes add first, then normalize | the current Fig 2: what to look for; population at 45%; added vs the summation prediction (the depth average of 2d, with its prediction); early vs late per penetration |
| 4 | Contrast sets how far; laminar organization | maps (response clock) at 11 and 22.5%, or the current line panels; transition index vs contrast; depth profiles; current added at the earliest sink |
| 5 | Delayed normalization | the current Fig 4 |
| S | Additivity index maps | both eyes / (dominant + other), both clocks, 11–45% (the legacy measure, for comparison with Mitchell et al. 2022) |

**Map index (recommendation: "added by the other eye", the additivity index in the supplement).**
* **The additivity index has no fixed summation reference.** With saturating responses, contrasts that add exactly
  still give both / (dominant + other) well below 1 (late, ~0.55–0.65 everywhere in the mock), so a reader takes
  "0.6" as sublinear when it can be what summation predicts.
* **Its early peak is mostly latency.** The index is near 1 only in the first ~5 ms, on the rising edge, where a
  1.5 ms binocular lead inflates the ratio: the confound the paper removes.
* **The added map shows the result directly** (orange = what the other eye adds, white = nothing), in the same
  terms as Figs 3–5: a burst for the first ~25 ms at every depth, then nothing in 30–100 ms, then a weak return. On
  the stimulus clock the burst is smeared and tilted (superficial bins peak at ~55 ms, against ~35–45 ms at and below
  the sink).
* **Stated with the existing S4c tests, not new ones:** normalization is reached 6–7 ms later above and below the
  sink on the stimulus clock (p = .016, .003), and at the same time after each group's own onset above (median
  difference 0 ms, p = .50); below, the median difference is also 0 ms, but p = .045.

**Blake to decide:** the map index; five main figures or Fig 3 folded into Fig 2; moving the CSD map into Fig 1;
which "Mitchell et al. 2024" is meant (the text cites Mitchell et al. 2022, iScience, and 2023, Curr Biol; docs/10
also lists the 2024 dissertation); and, still open, how to phrase "normalize later".

**Revised the same day after Blake's next notes (supersedes the table above).** Blake: the linearity index across
time and contrast (S2a) "is a strong and easily readable result, and effectively all we would need to show for the
add first, normalize afterwards … it can literally be one single figure panel. But I think Figure 2 really needs to
get to the spatiotemporal profile fast, so as to not step on the toes of the iScience paper, current bio paper, and
journal of vision paper from Cox." Fig 1 can make the latency points: conditions rank by drive and onsets follow
("it's important to look at the response in a manner that compares their peaks, not just windows or bins"); the
drive-by-latency panel is the last panel of Fig 1. Check H supplies it.

| Figure | Claim | Panels |
|---|---|---|
| 1 | The experiment; onset follows drive, and both eyes are drive | (a)–(e) as now; (f) onset vs peak response for all 35 conditions, with what adding the two monocular responses would give (H) |
| 2 | Space-time profile, fast | (a) linearity index over time at 11, 22.5, 45% (the S2a data, without the null lines): the whole temporal claim in one panel; (b) population CSD with the MUA onset at each depth (G4; moved from Fig 3d); (c, d) what the other eye adds, depth × time at 45%, stimulus clock vs each site's response clock |
| 3 | Laminar organization of normalization | response-clock maps at 11 and 22.5% (or depth profiles); transition index vs contrast; current added at the earliest sink |
| 4 | Delayed normalization | the current Fig 4 |
| S | Moved out of the main text | the current Fig 2 (sketch, time courses, added vs summation, early vs late per penetration); G's onset-vs-contrast and by-depth panels; additivity-index maps; S2a's fixed-rule null |

Open for Blake: whether the one temporal panel opens Fig 2 (recommended, as the bridge to the maps) or closes Fig 1;
the map index; which Mitchell et al. 2024; "normalize later". S2a still says "linear sum": relabel "contrasts add"
when it moves.

**Decided (Blake, 2026-11-03):** the one temporal panel opens Fig 2; the maps show "added by the other eye", with
the additivity index (both / (dominant + other)), which readers know, in the supplement and its relation stated in
the legend; "Mitchell et al. 2024" was Mitchell et al. 2023 (Curr Biol). Fig 2 is the space-time resolution of
binocular summation; Fig 3 the laminar organization of normalization. Open: how Fig 3 differs from the model figure
and what shows that normalization occurs (checks I and J above; proposal in the reply of 2026-11-03).

**Manuscript (from 2026-10-30):** the text lives in `manuscript/draft_vN.docx` (Word; Blake comments, Claude revises with
tracked changes). Blake's rule: pre-registration is internal discipline only. The paper reports every analysis on
the same (artifact-cleaned) data without registration language; the as-run registered outcomes stay in docs/09-10.

### Figure set for draft V4 (Blake, 2026-11-11)

Blake agreed that the full dichoptic design over time (checks Q, R) is the strongest result and should lead as Fig 2,
but was not satisfied with how it had been shown (`fig2_full_design.py`). Mock-ups compared: one matrix per time window;
the measured responses cell by cell; one line per other-eye contrast. The measured traces alone make the effect look
like nothing (the shortfall is 3–6% of the dominant-eye peak), and the difference alone is abstract. Blake chose the
combination: one enlarged cell for scale, the window means for the statistics, the matrices for the design. Onset
follows drive stays as Fig 1's final panel (Blake).

| Figure | Claim | Panels | Script |
|---|---|---|---|
| 1 | The experiment; onset follows drive | (a) setup, (b) 6 × 6 design, (c) example penetration, (d) onset vs response strength, 35 conditions (H, G1, U4) | `fig1_setup.py` |
| 2 | The eyes add at onset; afterwards both eyes fall short of summation in proportion to the other eye's contrast, and stay short when the other eye is stronger | (a) one cell (dominant 22.5%, other 90%): both eyes, "if the contrasts add", dominant eye alone; (b) window means vs the other eye's contrast with Q's slopes and the held-out R2; (c) the 3 × 4 design as a matrix per window | `fig2_dichoptic_time.py` |
| 3 | The laminar profile: normalization in and above the input layer, fading below; onset first just below the boundary, and the second eye advances it at every depth | (a) additivity-index ridgeline by depth, 45%; (b) the index by depth group, 30–100 ms (X2); (c) onset at each depth (X3, X4; formerly simplified Fig 2b) | `fig3_normalization_ai.py` |
| 4 | A delayed pool, weighted toward the other eye, at every depth | (a) sketch; (b) held-out gain of the delay (Y1); (c) drive vs pool weight (U6′) | `fig4_delayed_pool.py` |
| S | Each eye alone | (a) monocular contrast responses, (b) ocular dominance (moved from Fig 1) | `figS_monocular.py` |

Fig 2 is in response units (% of the dominant-eye peak), the units of checks Q and R; Fig 3 keeps the additivity
index (check V). Fig 3c–d of the simplified set (the index against the other eye's contrast) is superseded by Fig 2.
`fig2_latency_drive.py` stays as the module for the laminar onset panel; `fig2_full_design.py` is retired from the
main set. Legend caveat for Fig 2: part of the 0–20 ms surplus reflects the binocular response starting ~1.5 ms
earlier (shifting the prediction 1 ms earlier halves it), so the text says "at or above summation at onset", not
"supra-additive". The older supplementary figures (S1–S8) still use the linearity index and alignment v2 in places
and need rebuilding before submission.

### Blake's comments on draft V4 (2026-11-11) and check Z

Six comments; all handled in `manuscript/draft_v5.docx`, each answered in its thread.
* **"Summation" was ambiguous** (Blake: "do you mean linear summation?"). It meant the eyes' weighted contrasts add
  before the compressive contrast-response function, R_D(c_D + a·c_O), not the sum of the two monocular responses
  (which predicts AI = 1 against 0.62 at 45%). The text now says **contrast summation** throughout and defines it
  once against linear summation of responses; the figure label `figstyle.SUM` is now "contrast summation".
* **Caption bold:** only the "Figure N. Title." sentence is bold (the V4 build had bolded whole captions).
* **Fig 3a busy:** replaced by three traces in bands of five contacts each (0.5–1.0 mm above the boundary; the
  0.5 mm above it, holding the earliest sink; the 0.5 mm below), on Blake's suggestion of above / at / below the
  sink matched in contacts. Bands from the boundary, not the sink, so panel a uses the same reference as b and c;
  the input band holds the earliest sink in 17/22 penetrations. 30–100 ms: AI 0.49 / 0.52 / 0.55 against contrast
  summation 0.63 / 0.66 / 0.60 (19 / 22 / 21 penetrations). The ridgeline moved to Supplementary Fig S2
  (`figS_ridgeline.py`). Fig 3c marks the median depth of the earliest sink.
* **CSD beyond the earliest sink:** nothing in the paper depends on it. The sink-current analysis (old Supplementary
  Fig S2, `figS7_input_current.py`; "spiking led the current by 2–3 ms") is dropped from the paper; the earliest sink
  remains a marker for the depth reference.
* **"for every depth reference we tried"** dropped from the Discussion (the limitations already say it).
* **"Responses began earliest 0.1–0.5 mm below the boundary"** (Blake: reviewers will expect the input layers first;
  does the earliest sink still mark the initial input, with spiking fastest below it?) → check Z.

**Check Z (descriptive, not registered; run after the question; `scripts/z_onset_sink.py`, docs/qc/z_onset_sink.csv,
z_onset_profile.csv).** Onsets of check X (first sheet, v3), relative to the earliest monocular sink.
* **Z1:** the earliest sink lies a median 0.20 mm above the boundary (IQR 0.28 above to 0.00), inside the input
  group in 17/22 penetrations. By construction the boundary sits near sink + 0.3 mm where the sink is used.
* **Z2 (within penetrations, minus onset at the sink ± 0.15 mm):**
  * 0.2–0.4 mm below: dominant eye 90% +0.19 ms (earlier in 10/21, p = .84); both eyes 45% +0.51 ms (p = .92).
  * 0.5–0.8 mm below: −1.1 ms (14/20, p = .14); both eyes 45% −1.7 ms (13/20, p = .31); dominant 45% −2.1 ms
    (p = .33); same sign in both halves, not significant.
  * 0.2–0.4 mm above: −0.1 ms (p = .59). 0.5–0.8 mm above: +4.0 ms (3/19 earlier, p < .001); both eyes +4.4 ms.
* **Reading:** spiking begins as early at the earliest sink as anywhere; onsets form a plateau from the sink to
  ~0.7 mm below the boundary and are later above the sink. The "earliest 0.1–0.5 mm below the boundary" of V4 was a
  1–2 ms dip in population medians that the within-penetration comparison does not support. The text now reports
  the plateau, with deep channels responding as early as (not earlier than) the input layers.

## Checks AA and AB: the model's laminar difference and the pupil control on the current measures (registered internally 2026-11-12, before running; Blake: "Go ahead and re-run the depth-by-depth model comparison and the pupil-size control on the new measures")

**Why.** Two statements in the paper still rest on earlier analyses. (1) "Fitted at each depth, the model did not
reproduce the weaker normalization of the deep group": N's descriptive comparison, on the sink groups. (2) The pupil
control (#23) was run on the linearity index (S1, `fig2_robustness.py`), not on the shortfall from contrast summation.
**Not blind:**
* AA: on the sink groups, the conductance fits normalized about equally at every depth (−0.039 / −0.020 / −0.017 against
  data −0.056 / −0.061 / −0.020). On v3, the excess pool weight w − a is 0.23 / 0.28 / 0.13, deep − input p = .076.
* AB: on the old index, every claim held for pupils within the penetration's IQR. Nothing is computed on the new
  measures.

**AA. Measure.** Y1's group populations and fits on v3 (`y_N_pred.npz`; two-stage predictions from the
`y_P_groups.parquet` parameters). Per penetration and group, the shortfall at 45% in 30–100 ms after the group's onset:
both eyes minus the dominant eye at c(1 + a), with a, onset and the prediction computed from the same array, for the data and for each fit
(conductance, the Fig 4 model; fixed delay; two-stage), as N's `group_deficit`.
* **AA1 (reported):** medians per group, data and each model.
* **AA2, the model's own depth difference:** model (deep − input) > 0, one-sided Wilcoxon over penetrations, per model.
* **AA3, the gap:** data (deep − input) − model (deep − input) > 0, one-sided Wilcoxon, per model.
* **Wording:**
  * "the model did not reproduce the weaker deep normalization" only if, for the conductance model, AA3 holds and
    AA2 fails;
  * "reproduced it only in part" if both hold;
  * "reproduced it" if AA2 holds and AA3 fails.
  * The same wording applies to the two-stage model, reported alongside.

**AB. Subsets.** Clean presentations, split per penetration at the median of the pre-stimulus pupil signal (median
over −200…0 ms; the median is robust to blinks). The analog pupil signal has a per-monkey offset (negative in monkey
B), so only the within-penetration split is used. Subsets: "lower" and "higher" halves; "IQR" (as #23), reported.
Dominant eye, scale, onset and eye weight come from all clean presentations, as for the odd/even halves; conditions
need ≥ 2 presentations in a subset.

**AB tests, per subset:**
* **AB1 (F1):** at 45% in each eye, the shortfall A − P < 0 in 30–100 ms. Population channel mean per penetration,
  as `v4_dip_recovery.py`; one-sided Wilcoxon.
* **AB2 (Q1, Q3):** the x_O slope of the shortfall < 0 in 30–100 and in 130–200 ms, and the 130–200 ms shortfall
  where c_O > c_D < 0. Mixed model as `v4_full_design.py`; one-sided.
* **AB3 (X1):** (A − P) deep − input > 0 at 45%, 30–100 ms. X's channels and v3 groups; one-sided Wilcoxon.
* **AB4 (X4):** the binocular onset at 45% minus the dominant eye's at 45% < 0. Penetration median over X's
  channels; one-sided Wilcoxon.

**AB criteria.**
* **The control passes** if AB1–AB4 hold (p < .05, predicted sign) in both the lower and the higher half.
* **Reported:** the IQR subset; higher − lower for AB1, AB3 and AB4 (two-sided Wilcoxon); for AB2, the
  window:x_O:half interaction in one mixed model over both halves (two-sided).
* **Wording:**
  * If the control passes: "the results held for presentations with smaller and with larger pupils".
  * If a test fails in one half: say which, and give the higher − lower difference.

**Scripts:** `scripts/aa_model_misses_v3.py` → `docs/qc/aa_model_misses.csv`, `aa_groups.csv`;
`scripts/ab_pupil.py` → `docs/qc/ab_pupil.csv`.

**AA results (2026-11-12; 21 penetrations, 57 groups: superficial 16, input 20, deep 21; 20 penetrations with both
input and deep). The conductance and fixed-delay pools do not reproduce the weaker deep normalization; the two-stage
model does.**
* **AA1, shortfall at 45%, 30–100 ms after the group's onset (medians; data / conductance / fixed delay / two-stage):**

  | Group | Data | Conductance | Fixed delay | Two-stage |
  |---|---|---|---|---|
  | superficial | −0.059 | −0.046 | −0.046 | −0.045 |
  | input | −0.060 | −0.028 | −0.029 | −0.031 |
  | deep | −0.026 | −0.018 | −0.018 | +0.002 |

  Data deep − input: +0.045 (15/20, p = .0002).
* **Conductance (the Fig 4 model):**
  * AA2 fails: the model's own deep − input is +0.008 (13/20, p = .083).
  * AA3 holds: the data exceed it by +0.019 (14/20, p = .002).
  * Wording allowed: "did not reproduce the weaker deep normalization".
* **Fixed delay:** the same pattern. AA2 +0.009 (p = .072); AA3 +0.020 (p = .002).
* **Two-stage:**
  * AA2 holds: +0.034 (15/20, p = .002).
  * AA3 fails narrowly: the remaining gap is +0.005 (14/20, p = .053).
  * Wording allowed: "reproduced it".
* **Reading.**
  * The single-stage pools underestimate the shortfall at the input, about half the data's value, so they flatten
    the laminar profile.
  * The two-stage form, in which each eye is first gain-controlled by the other, leaves the deep group at contrast
    summation and so recovers most of the difference. It is also the best-fitting model at every depth (Y2, P3).
  * The remaining gap is borderline (p = .053), so "reproduced most of it" is the safe phrasing.
  * Not a mechanism claim: the two-stage model's extra stage is also its extra flexibility.
* **For the paper (Results, model paragraph):** "Fitted at each depth, the single-stage pools did not reproduce the
  weaker normalization of the deep group (data, deep − input +0.045; conductance fits +0.008, p = 0.08), whereas the
  two-stage form reproduced most of it (+0.034, p = 0.002; remaining gap p = 0.053)." The sentence "the laminar
  difference of Fig 3 therefore rests on the data, not on the model" is replaced.

**AB results (2026-11-12; 26 penetrations for AB1 and AB2, 22 for AB3 and AB4; a median of 545 presentations per
half). The control passes: AB1–AB4 hold in the lower and in the higher half, and in the IQR subset.** The "all" subset
reproduces F1 (−0.039, 23/26), the Q slopes (−0.0126, −0.0111), Q3 (−0.026) and X1 (+0.052, 18/22) exactly.

| Test | Lower half | Higher half | IQR |
|---|---|---|---|
| AB1 A − P, 45%, 30–100 ms | −0.040 (20/26, p = .0004) | −0.051 (23/26, p = .0001) | −0.045 (23/26) |
| AB2 x_O slope, 30–100 ms | −0.0138 (p = .0007) | −0.0139 (p = .0002) | −0.0129 |
| AB2 x_O slope, 130–200 ms | −0.0104 (p = .009) | −0.0119 (p = .001) | −0.0093 (p = .011) |
| AB2 Q3, c_O > c_D, 130–200 ms | −0.023 (p = .0008) | −0.026 (p < .0001) | −0.031 |
| AB3 deep − input | +0.055 (17/22, p = .0003) | +0.047 (17/22, p = .0003) | +0.038 (15/22, p = .002) |
| AB4 binocular advance, 45% | −1.8 ms (19/22) | −1.3 ms (17/22, p = .0008) | −1.1 ms (18/22) |

* **Higher − lower (reported, two-sided):**
  * AB1: −0.015 (18/26, p = .038). The matched-contrast shortfall is slightly larger with a larger pupil signal.
  * Nothing else differs: AB3 +0.003 (p = .82); AB4 +0.4 ms (p = .28); AB2 slope interactions +0.0008 and −0.0010
    (p ≥ .87); AB2 shortfall at x = 0 −0.0075 and +0.0027 (p ≥ .25).
* **Caveat:** "larger" is the higher analog pupil signal within a penetration. The signal's polarity was not
  calibrated (monkey B's is offset negative).
* **For the paper (Alternative explanations):** "every result held separately for presentations with smaller and with
  larger pupils before stimulus onset; only the shortfall at matched 45% contrasts differed between them, slightly
  larger with larger pupils (by 1.5% of the dominant-eye peak, p = 0.038)". The Methods describe the split.
* Both statements are in `manuscript/draft_v5_clean.docx` (Results model paragraph, Discussion laminar and
  alternatives paragraphs, Methods eye position).

## Fig 4 rebuilt (2026-11-12, descriptive; Blake's comment on the V5 clean copy)

The held-out-error panel is replaced by panels that show what each model feature does, all from Y1's stored fits
(`y_N_pred.npz`; the two-stage values from AA). Each array (data or fit) is taken on its own onset, eye weight and
contrast-summation prediction; per penetration the mean over its depth groups (`scripts/figures/fig4_delayed_pool.py`).

* **(b) Delay**, 45% in each eye, both eyes − contrast summation (% of dominant-eye peak), medians over 21 penetrations:

  | | 0–20 ms | 30–100 ms | 130–180 ms |
  |---|---|---|---|
  | data | +4.8 | −5.5 | −3.1 |
  | conductance | +1.2 | −2.7 | −2.2 |
  | fixed delay | +1.0 | −2.7 | −2.3 |
  | static pool | −2.2 | −2.8 | −2.8 |
  | w = a | −0.7 | +0.1 | +0.1 |

* **(c) Weight**, 30–100 ms, mean over c_D = 11, 22.5, 45% (where summation is predictable), by c_O = 11 / 22.5 / 45 /
  90%: data −1.9 / −3.0 / −4.7 / −7.0; conductance −0.5 / −1.2 / −2.6 / −5.8; w = a +0.1 at every c_O.
* **(e) Depth**: AA's values (data, conductance, two-stage).
* **Reading:** the delayed, other-eye-weighted pool reproduces the sign, timing and contrast dependence of the sequence,
  at about half its size; a pool without a delay never exceeds contrast summation; equal weights never fall short of it.


## Check AC: deep onsets penetration by penetration (2026-11-12, descriptive, post hoc; Blake's comment on V6)

**Why.** Blake: the deep group responding as early as the input group, while normalizing less, could come from
penetrations whose boundary sits too high (the earliest sink placed too high), or from a probe crossing a second fold
of cortex, so that a second input layer is counted as deep. Is the deep onset a per-penetration regularity, or a few
sessions? **Not registered:** X3 (deep − input −1.35 ms, 14/22, p = .26) and Z were known. Script
`scripts/ac_deep_onset.py`; tables [qc/ac_penetrations.csv](qc/ac_penetrations.csv),
[qc/ac_tests.csv](qc/ac_tests.csv), [qc/ac_profiles.csv](qc/ac_profiles.csv); per-penetration panels
`<out_root>/figures/ac_deep_onset_cards.png` (onset against depth, with the laminar audit's whole-probe monocular onset
beyond the 0.7 mm cap). "Markers agree": the earliest monocular sink lies inside the input band (0–0.5 mm above the
v3 boundary).

* **AC1. Not a per-penetration regularity.** Deep − input onset (dominant eye, 90%): earlier in 14/22, later in 8
  (median −1.35 ms, p = .26). Where the odd and even halves agree in sign: 11 earlier, 5 later. The later ones
  include large delays in monkey B (211008 +8.5, 211012 +14.6 ms). 211025_B_eA has 3 input and 3 deep channels, and
  its halves differ by 30 ms.
* **AC3. It is a monkey J pattern and follows the markers' disagreement.**
  * Monkey J: 11/15 earlier, median −2.2 ms (p = .018); monkey B: 3/7 earlier, median +1.4 ms (p = .38).
  * Markers disagree (7, all J and all flip-anchored: the sink 0.55–1.05 mm above the flip in 6, 0.75 mm below it in
    221215): 7/7 earlier, median −2.7 ms (p = .016). Markers agree (15): median +0.5 ms, 7 earlier, 8 later (p = .85).
    Four of the 7 have the sink one contact above the input band (0.55 mm), so the split depends partly on the band
    edge. As a continuous measure, the higher the sink above the boundary, the earlier deep is than input
    (Spearman ρ = −0.38, p = .085, n = 22).
  * With the sink-based groups (input = the 5 contacts at the sink), deep − input is +0.75 ms (15/22 later, p = .39).
* **AC2 and the cards. No sign of a second fold.** No penetration has a second early zone separated from the first by
  later contacts within the deep group. In J the onset falls steadily down the probe, is earliest 0.3–0.7 mm below the
  boundary (lower deep − input −1.2 ms, 15/22; upper deep 0–0.3 mm −0.2 ms, 12/22) and becomes later beyond the cap.
  The audit's one candidate fold (211217_B_eA, an early zone 1.75 mm below) lies beyond the cap, and its deep group is
  later than its input group (+3.9 ms). Geometry also argues against a fold: a second sheet's layer 4C would lie
  beneath this sheet's layers 5–6, the white matter and the second sheet's 6 and 5, well beyond 0.7 mm.
* **AC4. The early deep channels do not normalize like the input group.** Deep channels with onsets at or before
  their input group's median fall short of contrast summation by less than the input group (+0.035, 16/21, p = .002;
  later deep channels +0.042, 15/17). Across penetrations, deep − input onset and deep − input shortfall are
  unrelated (ρ = 0.01). In the 7 where the markers disagree, the contacts 0.3–0.7 mm below the boundary are the
  earliest (−3.7 ms relative to the penetration) and barely fall short (A − P ≈ 0), while the input band falls short
  (−0.05). The normalization change sits at the boundary in both subsets.
* **The normalization result does not depend on these penetrations.**

  | deep − input | markers agree (15) | markers disagree (7) | B (7) | J (15) |
  |---|---|---|---|---|
  | X1 shortfall | +0.056, 14/15, p = .0002 | +0.048, 4/7, p = .11 | +0.035, 6/7, p = .023 | +0.056, 12/15, p = .0008 |
  | X2 AI − AI_sum | +0.127, 13/14, p = .0006 | +0.061, 5/7, p = .039 | +0.095, 5/6, p = .078 | +0.110, 13/15, p = .0004 |
  | X3 onset (ms) | +0.5, 8/15 later, p = .85 | −2.7, 0/7 later, p = .016 | +1.4, 4/7, p = .38 | −2.2, 4/15, p = .018 |

* **Reading.**
  * "Deep channels respond as early as the input layer" is not a general finding. Deep onsets do not differ from input
    onsets where the laminar markers agree, are later in several B penetrations, and are earlier in monkey J, most
    where the earliest sink sits well above the phase flip.
  * In those J penetrations the earliest spiking lies 0.3–0.7 mm below the flip and 0.6–1.5 mm below the earliest
    sink (audit W: sink and MUA disagree in J). These contacts could be layer 6 with direct geniculate input, or the
    bottom of a layer 4C that both markers place too high. Latency alone cannot tell which. Either way, they do not
    normalize like the input band above them.
  * The deeper normalization result holds where the markers agree (the stronger subset) and in each monkey, so it does
    not come from misplaced boundaries. A misplaced input layer counted as deep would, if anything, shrink the
    difference.
* **For the paper (proposed):** drop "deep channels respond as early" from the Fig 3 title, the Results and the
  Discussion. Report the superficial delay (5 ms), and that deep onsets did not differ consistently from input onsets
  (earlier in monkey J, not in B). Add the markers-agree subset as a sensitivity line for X1/X2.
* **Decision (Blake, 2026-11-12):** drop "deep channels respond as early as the input layer" from the paper. A
  literature review of where the initial sink, the phase flip and the spectral crossing lie relative to the anatomy
  (`reports/Macaque V1 laminar markers anatomy.md`) comes first, because it may change the anchor and the compartments.

**Laminar literature review (2026-11-12; `reports/Macaque V1 laminar markers anatomy.md`, notes in
`research_notes/Macaque V1 laminar markers anatomy/`).** No study has measured, against histology, how far the 4C/5
border lies below the initial sink in macaque V1. The offsets in use are conventions. The usable constraint is the
layer budget: about 1.1–1.4 mm of cortex above the 4C/5 border (4C ≈ 0.28 mm) and about 0.5 mm below it (layers 5–6
≈ 0.49 mm). 4Cα and layer 6 both spike earliest, so early spiking 0.3–0.7 mm below a correct border is expected.
* **First gate applied (evoked-MUA top to anchor, expected 1.1–1.4 mm; `top_to_anchor_mm` in alignment v3):**
  * The 7 J penetrations where the markers disagree: 1.25–1.55 mm (221214_J_eB 2.75). The flip passes. Their sinks
    (0.55–1.05 mm above the flip) lie 0.2–0.9 mm below the top of cortex, too high for layer 4C, so in these
    penetrations the sink is the marker that failed. The early deep contacts lie where layer 6 should be.
  * Short (< 0.9 mm): 221209_J_eA, 211012, 211027, 221109 and 221214_J_eA have the top of cortex at or near the end of
    a descending probe (top censored, so the distance is a lower bound); 211103_B_eA (0.75 mm, top contact 5) is a
    real flag.
  * Long (> 1.5 mm): 221202_J_eA 2.25, 221214_J_eB 2.75, 211008_B_eA 1.75, 221209_J_eB, 221206_J_eA and 221212_J_eA
    1.55 (oblique probe, a high top, or a deep anchor).
* **Recommendations for Blake (not yet adopted):** keep the 4C/5 border as the target and the v3 anchor, and describe
  the input band as 4C plus 4B. Add the thickness gates as flags. Add sensitivity analyses: rigid anchor shifts of
  ±0.1 and ±0.2 mm; input 0–0.3 mm; deep 0–0.5 mm, split into L5 (0 to −0.25) and L6 (−0.25 to −0.5); a guard band;
  and the gate-passing subset.
