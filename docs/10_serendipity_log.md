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

## Log

| # | Effect | Where seen | Size and consistency | Known? | Status | Next step |
|---|---|---|---|---|---|---|
| 1 | **Interocular normalization arrives sooner at higher contrast.** The early-linear → late-normalized transition happens 20–35 ms after onset at .45, 35–50 at .225, ~50–75 at .11 | Fig 2 first pass (`equivalent_contrast.py`) | monotonic over 3 contrasts; p ≤ .016 at each; population level only so far | contrast-dependent dynamics of *interocular* normalization not found in the sweep (Albrecht 2002 reports time-invariant c50 for monocular normalization, anesthetized) | noted (on-target-adjacent: part of Fig 2, could stand alone) | per-penetration transition time vs contrast; both monkeys separately; split-half |
| 2 | **Monkey difference in MUA onset latency:** J 26–31 ms vs B 34–37 ms | dataset check (2026-10-28) | ~8 ms, all penetrations | latency varies with eccentricity, SF, size and species; B and J differ in all of these | noted, probably confounded | regress latency on eccentricity, SF, size; only interesting if a residual monkey effect remains |
| 3 | **Late binocular linearity is an inverted U across depth:** middle depths stay linear, superficial and deep normalize | Fig 3 first pass | mixed model p < .001; upper vs middle p = .005 | not found in the sweep | on-target (Fig 3) | landmark compartments, CSD (Fig 4), per-monkey |
| 4 | **Flash CSD stripes from contact gain mismatch** on 100 µm Plexon probes; outlier-contact interpolation + Hamming smoothing fixes them | 05 C12 | systematic, many penetrations | known in principle (CSD sensitivity to gain) | methods note | cite in methods; not a paper |
| 5 | **Kilosort spike times offset by exactly 2 × T0 (204 samples)** in one 2023 sort batch | 05 A7 | exact, all affected sorts | pipeline-specific | methods note | cite in methods |
