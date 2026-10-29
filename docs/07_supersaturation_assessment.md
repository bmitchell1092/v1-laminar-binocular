# 07 — Supersaturation vs depth: honest assessment (2026-10-28)

**Question (Blake):** does contrast supersaturation (responses that fall at high contrast) vary with cortical
depth or laminar compartment in V1, as a probe of the stabilized supralinear network (SSN)? Is this dataset
powered to answer it, and if layer labels are too uncertain, does depth below the top of cortex suffice?

**Short answer: no detectable supersaturation, at any depth.** Multi-unit activity is well powered and shows
none. Single units show a downturn at the rate expected by chance, but too few units are both isolated and
driven to rule out a minority of supersaturating neurons. Uncertain layer labels are not what limits this: the
effect itself is absent where we have power.

Scripts: `scripts/run_saccades.py`, `scripts/run_crf.py`, `scripts/assess_supersaturation.py`,
`scripts/plot_supersaturation.py`. Tables: [qc/supersaturation_assessment_*.csv](qc/). Figure:
`output/crf/supersaturation_assessment.png` (not in git).

## What went in

| | |
|---|---|
| Penetrations | 26 depth-usable (B 8, J 18); 22 with an L4/5 landmark (1 tier A, 21 tier B; 05 C13) |
| Stimulus | one grating per penetration at the probe's RF (1–4°, 1–4 cpd), 6×6 contrast design {0, .055, .11, .225, .45, .9} |
| Presentations | completed, no saccade ≥ 1° from −200 to 500 ms (05 C15); median ~30 per contrast cell (211008: 10; 211012: 19) |
| Responses | blank-subtracted window means, transient 40–140 ms and sustained 150–250 ms after the photodiode flash |
| CRFs | monocular contra (ipsi eye at 0), monocular ipsi, dioptic (both eyes equal) |
| MUA | 581 in-cortex channels driven by the grating (one-way ANOVA across contrasts p < .01) |
| Single units | 68 strictly isolated (contamination < 10%, 05 C14) / 136 lenient (< 50%); driven per condition: 14–18 strict, 34–51 lenient, from 7–19 penetrations, concentrated in 211012, 211025, 211008, 221123, 211217 |

## Statistics built not to fool us

* The plain supersaturation index (peak − response at 0.9) / peak is **biased upward by noise**: the largest of
  six noisy means overshoots the true peak, so even a rising CRF scores > 0. We use a **cross-validated index**
  (`crf.cv_supersaturation_index`: peak contrast picked on half the presentations, drop measured on the other
  half), whose expectation is ≤ 0 for any monotonic CRF.
* Per site, a **held-out downturn test** (`crf.downturn_p`: same split; one-sided Welch test that the response at
  the held-out peak exceeds the response at 0.9; splits combined as 2 × median p). Calibrated by simulation:
  ≤ 5% of saturating sites reach p < .05; > 80% of supersaturating ones do. The fraction of sites with p < .05
  is compared with 5% (binomial).
* Compartment effects two ways: a site-level mixed model with a random intercept per penetration, and a
  conservative per-penetration comparison of compartment means (neighbouring MUA channels share signal, so
  sites are not independent). Power from the per-penetration differences.

## Results

**1. MUA: no supersaturation.** Fraction of driven sites with a held-out downturn (5% expected by chance):

| window | mono contra | mono ipsi | dioptic |
|---|---|---|---|
| transient | 0.0% (0/529) | 0.0% (0/539) | 0.6% (3/534) |
| sustained | 2.8% (14/503) | 1.4% (7/510) | 3.9% (18/463) |

Median cross-validated index −0.03 to 0.00 in every condition; the median ratio R(0.9)/peak is 1.0 (MUA CRFs
peak at the top contrast). The few downturns (sustained, mostly dioptic) peak at .225–.45 and fall to ~75% of
peak at .9. Fits are physiological: MUA c50 median 0.07–0.10, exponent ~2.5.

**2. Single units: at chance, with wide bounds.** Lenient units 5.6–9.8% (2–4 units per condition), binomial
p 0.15–0.54; strict units 0–14% (0–2 of 14–18). Upper 95% bound on the fraction of driven units with a
downturn detectable at these trial counts ≈ 23% (4/41, dioptic sustained). Unit medians of the
cross-validated index are negative. A few units fall *below* spontaneous at high contrast (index > 1); that is
suppression, not a supersaturating CRF, and is worth a separate look.

**3. Depth and compartment.** MUA is powered for this: with 17–21 penetrations per comparison, 80% power for a
0.10 difference in the index needs 4–18 penetrations (most conditions ≤ 8). There is no supersaturation to
compare, but there is a small, consistent trend: **deeper sites saturate slightly less** (index slope −0.01 to
−0.02 per mm below the top of cortex, p < .05 in 5 of 6 condition × window cells; infragranular vs granular
−0.03 to −0.07, mixed-model p ≤ .007 in 4 of 6). Single units cannot be split by compartment (≤ 8
penetrations contribute units to any compartment).

**4. Binocular.** Dioptic CRFs are not more supersaturating than the better monocular CRF at the same site
(per-penetration mean difference −0.027, p = .03, 22 penetrations): if anything, slightly less.

## What could hide supersaturation here

* Contrast sampling is coarse at the top (only .45 and .9 above .225); a downturn between .5 and .9 shows only
  as R(.9) < R(.45).
* ~30 presentations per contrast cell; low single-unit rates (median 3 Hz).
* Isolation (05 C14): 100 µm contact spacing and small J spikes limit single-unit yield; Kilosort 4 did not
  improve it.
* One stimulus size/position per penetration, chosen for one probe's RF, so many units sit off-centre.
* Layer labels are uncertain (05 C13), but that does not drive the null: the depth-only analysis is equally flat.

## Verdict and options

Supersaturation as a function of depth is **not supported by these data**: absent in well-powered MUA at every
depth, and at chance in the single units we can isolate (too few to exclude a minority population). Options:

1. **Report as a bounded null** inside a broader paper (MUA: none at any depth; units: < ~23% at these trial
   counts).
2. **Stay within the SSN/normalization question at the level MUA is powered for**: the laminar profile of
   contrast saturation (c50, exponent, compression) and binocular combination across depth. The depth trend in
   item 3 is the kind of result this design can support. It overlaps the earlier LMN paper, so scope it narrowly.
3. **Move to priorities 2 and 3**: the literature sweep, and dual-probe inter-column timing (7 J days; the spike
   and stimulus timing fixes in 05 A2/A7 were prerequisites).
4. **If supersaturation stays the goal**, it needs new data: denser high-contrast sampling (.45–1.0), more
   presentations, stimuli optimised per unit, and a high-density probe for single-unit yield.
