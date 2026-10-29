# 08 — Which open binocular questions can this dataset answer? (first pass, 2026-10-29)

Follows the literature sweep ([reports/V1 binocular contrast normalization.md](../reports/V1%20binocular%20contrast%20normalization.md),
question numbers below refer to its ranked table). First-pass checks only: they size the effects and the
precision, they are not the final analyses. Scripts: `scripts/run_dichoptic.py` (6×6 matrices per site and
window), `scripts/binocular_feasibility.py` ([qc/binocular_feasibility.csv](qc/binocular_feasibility.csv)),
`scripts/binocular_timecourse.py` (figure `output/crf/interocular_timecourse.png`).

Novelty bar: this dataset underlies Mitchell et al. 2022 (iScience) and 2023 (Curr Biol), which used pooled MUA
over 50–200 ms. New work needs depth, epochs/onsets on the corrected clock (05 A2, A7), isolated units, the two
simultaneous probes, or LFP/CSD. Check overlap with the laminar chapters of Blake's 2024 dissertation and the
legacy `2024_LaminarCRFs` manuscript before claiming novelty.

## Sample

565 driven in-cortex MUA channels in 26 penetrations; 83% are ocularly balanced (|ODI| < 0.3 at .9 contrast; MUA
pools neighbouring columns), median |ODI| 0.10. Penetration-level bootstrap CIs (resampling penetrations, so
neighbouring channels do not inflate precision) are about ±0.05 for population medians and ±0.1 within a
5-contact (0.5 mm) depth bin. Isolated units: 68 strict / 136 lenient, ~15–50 driven per condition.

## First looks

**Q7 ocularity invariance (additivity index AI = R(c,c) / (R(c,0) + R(0,c)); 1 = linear, 0.5 = invariant for
balanced sites).** Transient (40–140 ms): 0.85 [0.75, 0.99] at .055, 0.70 at .11, 0.63 at .225, 0.60 [0.57, 0.62]
at .45. Sustained (150–250 ms): 0.65 at .055, 0.58–0.62 above. Combination is near-linear at low contrast early and
approaches ocularity invariance at high contrast and late. Caveat: AI depends on the output nonlinearity, so the
low-contrast value needs a model (Q2) before it can be read as linear summation.

**Q2 summation rule (octave ratio O(c) = R(c,c) / R_D(2c); early summation predicts 1 on balanced sites).**
O = 0.85 (.055), 0.89 (.11), 0.97 (.225), 0.95 (.45) transient; 0.71, 0.86, 1.01, 0.92 sustained. O < 1 at low
contrast is partly an artefact: |ODI| < 0.3 still allows eye weights down to ~0.5, and with an accelerating CRF
at low contrast unequal weights alone push O below 1. **Next:** fit each eye's gain (non-dominant CRF mapped onto
the dominant one by a contrast scale a), then test R(c,c) against the early-summation prediction R_D((1 + a)c), and
fit Moradi–Heeger Eq. 10 (eye-specific gains, free w) per channel and epoch.

**Q1 timing of interocular modulation (R(D=.9, N) − R(D=.9, 0), normalised to the monocular peak, by depth).**
Clear at the population level. Adding the second eye *facilitates* the response from 20–37 ms after the flash at
every depth, i.e. at the monocular response onset (half-max 30–39 ms), peaking near 45 ms (+4% to +8% of peak at
N = .9, 40–100 ms). The facilitation is gone by ~70 ms; from 150–250 ms the binocular response equals the
monocular one at 0–20 contacts below the top of cortex (difference −0.03 to +0.01) and is *suppressed* in the
deepest bin (20–32 contacts: −6% to −8%). So: summation at onset, normalization arriving ~20–30 ms later, net
suppression only deep. Population onsets resolve to a few ms; single-channel onsets will be noisy. **Next:**
onset of the facilitation collapse (normalization latency) per depth and per N; same for CSD sinks (Q4).

## Verdicts

| # | Question | Verdict | Why |
|---|---|---|---|
| 1 | Onset of interocular suppression/facilitation by depth | **Go** | population time courses resolve onsets to a few ms; effect +4–8% early, −6–8% late deep |
| 2 | Summation rule (k, octave test) by epoch | **Go, with eye-gain correction** | precision ±0.05; needs eye-weight fits |
| 3 | Depth profile of interocular weight w | **Go (depth axis)** | ±0.1 per 0.5 mm bin; named layers only as secondary |
| 4 | Normalization in the earliest current sink | **Go on ~22 penetrations** | requires per-condition CSD; not yet checked |
| 7 | Ocularity invariance by signal and depth | **Go** | AI already shows the contrast and epoch dependence |
| 5, 8, 9 | OD and the pool, gain mode, dichoptic dip | **Limited** | need isolated units (15–50 driven) or pooled MUA, where OD is diluted |
| 13 | Inter-column timing (dual probes) | **Priority 3** | 7 J days; first check that one grating drives both RFs |
