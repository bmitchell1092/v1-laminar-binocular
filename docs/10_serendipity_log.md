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

## Log

| # | Effect | Where seen | Size and consistency | Known? | Status | Next step |
|---|---|---|---|---|---|---|
| 1 | **Interocular normalization arrives sooner at higher contrast.** Early-linear → late-normalized transition 20–35 ms after onset at .45, 35–50 at .225, ~50–75 at .11 | Fig 2 (`equivalent_contrast.py`) | **Per penetration:** L is flat across contrast at onset (0–20 ms: 0.99/1.02/0.93; falls with contrast in 69% of penetrations, ns), falls with contrast in the transition window (20–50 ms: 0.86/0.58/−0.22; **96% of 26 penetrations, B 8/8, J 17/18, p < 10⁻⁵**; 35–50 ms alone: 26/26), and noisier late (69%). Transition time monotonic in 84% of penetrations (Friedman p = .007; window resolution causes ties). B transitions later than J. Medians are robust to top-censoring. | contrast-dependent dynamics of *interocular* normalization not found in the sweep (Albrecht 2002: time-invariant c50 for monocular normalization, anesthetized) | **candidate** (meets bar 2). Timing vs strength not yet separated: a contrast-dependent late strength would itself be the two-stage / Ding–Sperling signature | raw-data figure; finer time resolution (sliding windows or a fitted transition time per penetration × contrast); held-out split; model: does a contrast-driven pool with fixed delay predict earlier transitions? |
| 2 | **Monkey difference in MUA onset latency:** J 26–31 ms vs B 34–37 ms | dataset check (2026-10-28) | ~8 ms, all penetrations | latency varies with eccentricity, SF, size and species; B and J differ in all of these | noted, probably confounded | regress latency on eccentricity, SF, size; only interesting if a residual monkey effect remains |
| 3 | **Late binocular linearity is an inverted U across depth:** middle depths stay linear, superficial and deep normalize | Fig 3 first pass | mixed model p < .001; upper vs middle p = .005 | not found in the sweep | on-target (Fig 3) | landmark compartments, CSD (Fig 4), per-monkey |
| 4 | **Flash CSD stripes from contact gain mismatch** on 100 µm Plexon probes; outlier-contact interpolation + Hamming smoothing fixes them | 05 C12 | systematic, many penetrations | known in principle (CSD sensitivity to gain) | methods note | cite in methods; not a paper |
| 5 | **Kilosort spike times offset by exactly 2 × T0 (204 samples)** in one 2023 sort batch | 05 A7 | exact, all affected sorts | pipeline-specific | methods note | cite in methods |
| 6 | **Non-dominant-eye current at the earliest sink is transient.** The sink's eye weight a falls from 0.69 (0–35 ms after sink onset) to 0.40 (50–150 ms); MUA at the same channels shows a weaker version (0.68 → 0.59) | Fig 4 (`dichoptic_csd.py`) | sink: lower late in 65% of 26 penetrations (B 88%, J 56%), p = .004; MUA 69%, p = .02 | not checked | noted (fails the ≥ 80% bar) | could be a noise floor in the small late non-dominant sink; recheck with a noise-corrected eye weight and per-monkey; revisit if it strengthens |
