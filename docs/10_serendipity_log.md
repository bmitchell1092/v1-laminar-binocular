# 10 — Serendipity log: robust effects found along the way

Standing practice (Blake, 2026-10-29). Any effect that turns up during the analyses, even off-target, is logged here
if it looks robust, for example a monotonic relationship between two variables, or a covariate that explains part of
the laminar response latencies and has not been described. A serendipitous finding can become its own paper, but
only if it can be nailed down beyond doubt.

## The bar ("beyond a shadow of a doubt")

A candidate is promoted to **confirmed** only when it meets all of these:

1. **Visible in raw data:** in plots of the raw responses, not only in a fitted statistic.
2. **Replicates:**
   * same sign in ≥ 80% of penetrations;
   * same sign in both monkeys separately;
   * monotonic (dose–response) where a graded variable is involved.
3. **Survives the confounds register** (docs/09): latency, eye weights, fixation, vergence, adaptation, depth
   assignment, RF position, monkey.
4. **Survives a held-out confirmation.** A discovery/confirmation split (e.g. alternate penetrations, or first vs
   second half of trials) is fixed before the confirmatory test, and the effect holds in the held-out part.
5. **Not already described:** checked against the literature, and against Mitchell et al. 2022/2023, the 2024
   dissertation and the legacy manuscripts.
6. **Mechanistically interpretable,** or at least a clear statement of what it would mean.

Status levels: **noted** (seen once) → **candidate** (meets 1–2) → **confirmed** (meets 1–5) → **paper** (6, with a
design of its own).

## Screen built into the analyses

The S1 confound controls double as a screen. Each of these covariates is estimated against response amplitude and
onset latency at each depth, so a robust off-target relationship surfaces without a separate fishing expedition:

* time since the last saccade or microsaccade (fixation age);
* pupil size before onset;
* position of the presentation within the trial (1st, 2nd, 3rd) and the previous presentation's contrast;
* time in session;
* RF eccentricity, grating size and spatial frequency (across penetrations);
* monkey.

Multiple comparisons are acknowledged: anything found this way is a discovery that needs the held-out confirmation
(criterion 4).

## Pre-registered confirmations

### Candidate 1: interocular normalization arrives sooner at higher contrast (registered 2026-10-29, before looking)

* **Held-out data:** the dual-probe "cross-recordings", never analyzed before. On each dual-probe day, each
  probe was also recorded during the bminteroc run whose grating sat on the *other* probe's RF. Pairs are
  221202, 221206, 221209 and 221214 eA/eB, and 221215 eA, as allowed by alignment v2. These were not used in the
  discovery: the discovery used each penetration's own run. Only channels driven by the grating count (one-way
  ANOVA across monocular contrasts, p < .01, as in the discovery). A cross-recording is included if it has ≥ 5
  such channels.
* **Statistic (identical to the discovery):** per channel, the linearity index L with window-specific eye
  weights, in the 20–50 ms window after that channel's onset, for dioptic .11, .225, .45. Per cross-recording:
  the slope of the channel-median L on log contrast.
* **Prediction and success criterion:**
  * the slope is negative in ≥ 80% of included cross-recordings;
  * a one-sided sign test is p < .05;
  * if fewer than 5 cross-recordings qualify, the test is inconclusive (not a failure).
* **Also reported, labelled as reliability only (not independent):** odd vs even presentations within the
  discovery penetrations.
* **Result (2026-10-29): FAIL by the pre-registered criterion, narrowly.**
  * **Sample:** 8 cross-recordings had driven channels; 7 had ≥ 5 and were included. 221209 eB's cross-recording
    had none: the grating missed its RF.
  * **Outcome:** the slope is negative in 6/7 (86%, meeting the ≥ 80% bar), but the one-sided sign test gives
    p = .063 (> .05).
  * **Design flaw:** with 7 recordings only 7/7 could pass (p = .008); the test was underpowered by design.
    Lesson: size held-out tests before registering them.
  * **Direction and size replicate:** median slope −0.87 (discovery −0.92); early–late drops 0.32 / 0.69 / 0.93
    at .11 / .225 / .45 (discovery 0.49 / 0.59 / 0.85).
  * **Limits:** the held-out data are J only, from 4 days, and eA/eB pairs share a day.
  * **Reliability:** odd/even splits 25/25 and 22/24 negative.
  * **Status stays "candidate", no post-hoc rescue.** A future confirmation needs new data, or a registered
    test that can pass (e.g. ≥ 10 independent recordings).

## Log

| # | Effect | Where seen | Size and consistency | Known? | Status | Next step |
|---|---|---|---|---|---|---|
| 1 | **Binocular summation collapses further and faster at higher contrast** (first framed as "normalization arrives sooner at higher contrast") | Fig 2 (`equivalent_contrast.py`, `normalization_timing.py`) | **Robust:** L falls with contrast in the 20–50 ms window in 26/26 penetrations (both monkeys). Replicates in odd and even presentations separately (25/25, 22/24; reliability, not independence). Visible in raw data: binocular gain at onset is larger and collapses faster at higher contrast. At .45 it peaks at ~0.25 of the monocular peak and is gone by ~30 ms; at .225 it peaks at 0.19 and is gone by ~40 ms; at .11 it peaks at 0.14, decays slowly and stays positive. Late strength rises with contrast (early–late drop 0.49 / 0.59 / 0.85). **Not robust:** "sooner" independent of strength. Amplitude-normalized half-transition times (median 25 / 30 / 15 ms) put .45 first in only 50% of penetrations (p = .23); odd/even reliability r = 0.39. | strength dependence: a contrast-dependent late k is the two-stage / Ding–Sperling signature, and L < 0 at .45 (binocular below the dominant eye alone) exceeds ocularity invariance. Timing dependence: not found in the sweep | **candidate** (strength and population-level speed); the timing-only claim is not supported. **Held-out test failed narrowly** (6/7 negative, p = .063; see above). Survives every S1 control (no saccades, pupil, adaptation, session halves, per monkey: 88–100% of penetrations, p ≤ .008; docs/09) | pre-registered held-out test on the cross-recordings (above); if it holds, model which normalization architectures give L < 0 at high contrast |
| 2 | **Monkey difference in MUA onset latency:** J 26–31 ms vs B 34–37 ms | dataset check (2026-10-28) | ~8 ms, all penetrations | latency varies with eccentricity, SF, size and species; B and J differ in all of these | noted, probably confounded | regress latency on eccentricity, SF, size; only interesting if a residual monkey effect remains |
| 3 | **Late binocular linearity is an inverted U across depth:** middle depths stay linear, superficial and deep normalize | Fig 3 first pass | mixed model p < .001; upper vs middle p = .005 | not found in the sweep | on-target (Fig 3) | landmark compartments, CSD (Fig 4), per-monkey |
| 4 | **Flash CSD stripes from contact gain mismatch** on 100 µm Plexon probes; outlier-contact interpolation + Hamming smoothing fixes them | 05 C12 | systematic, many penetrations | known in principle (CSD sensitivity to gain) | methods note | cite in methods; not a paper |
| 5 | **Kilosort spike times offset by exactly 2 × T0 (204 samples)** in one 2023 sort batch | 05 A7 | exact, all affected sorts | pipeline-specific | methods note | cite in methods |
| 6 | **Non-dominant-eye current at the earliest sink is transient.** The sink's eye weight a falls from 0.69 (0–35 ms after sink onset) to 0.40 (50–150 ms); MUA at the same channels shows a weaker version (0.68 → 0.59) | Fig 4 (`dichoptic_csd.py`) | sink: lower late in 65% of 26 penetrations (B 88%, J 56%), p = .004; MUA 69%, p = .02 | not checked | noted (fails the ≥ 80% bar) | could be a noise floor in the small late non-dominant sink; recheck with a noise-corrected eye weight and per-monkey; revisit if it strengthens |
