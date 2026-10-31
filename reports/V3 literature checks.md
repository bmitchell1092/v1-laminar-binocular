# V3 literature checks: latency, time course and vergence

Checked 2026-10-31 against PubMed records, PMC full text, Crossref metadata and publisher pages. Preparations are tagged inline: **mac-W** (awake, fixating macaque), **mac-A** (anesthetized macaque), **cat-A** (anesthetized cat), **human**, **mouse**. A plain "macaque" or "cat" means the abstract did not state the preparation. Verification status for every reference is in the list at the end. "Verified" means the record or full text was read at the source (PubMed abstract, PMC article, publisher or Crossref record). "Cited only" means the claim was seen in another paper and not checked at its source.

## Novelty: what is already published

"Summation first, suppression later" is **not new as a qualitative sequence**. Four Maier-lab papers already report it in awake macaque V1, two of them with the same rig and stimulus family.

- **Cox et al. 2019** (J Vis) describe "an initial additive combination of the two eyes' signals" followed by differentiation of concordant and discordant input. They used fixed windows (transient 50–100 ms, sustained 150–250 ms), one contrast pair (~0.45/0.9) and no linearity test. Their Discussion calls the early phase "(subadditive) binocular summation", which contradicts "additive" in the abstract.
- **Dougherty et al. 2019** (Curr Biol) show in monocular neurons that binocular **facilitation occurs at response onset** in layer 4, whereas **suppression begins "several milliseconds" later**. The median modulation onset was 62 ms against a median response onset of 50 ms, so suppression arrived about 12 ms after onset. All stimuli were high contrast (≥0.8) and at zero disparity.
- **Mitchell et al. 2022** (iScience) are the closest precedent. They report an early phase (50–150 ms) in which binocular facilitation tracked contrast, then contrast-dependent control: by 150–250 ms facilitation persisted at low contrast and became **suppression at medium and high contrast**. Their Introduction calls this "initial, rapid (50–100 ms) summation that was more contrast-invariant followed by slower, contrast-dependent processing". The paper already contains the qualitative claim "more strongly at higher contrast". Over 0–250 ms, binocular responses were ~56% of the summed monocular responses.
- **Carlson et al. 2023** (iScience) confirm dioptic facilitation in the transient (50–150 ms) window only.

Two non-Maier papers give **millisecond-scale delays between binocular excitation and suppression**:

- **Tanabe & Cumming 2014** (mac-W): suppressive binocular (disparity) components lag excitation by a median of **7 ms**.
- **Dougherty et al. 2019**: the ~12 ms figure above.

What remains new, as far as these searches reach:

1. An explicit, time-resolved test that the binocular response equals L + R for about 20 ms after onset, across a full 6×6 dichoptic contrast matrix.
2. A millisecond onset for interocular normalization (20–35 ms after response onset) and its dependence on contrast. Earlier work used windows of 50–100 ms or more.
3. A binocular onset latency 1–2 ms shorter than monocular in macaque V1 (none of the Maier-lab papers report this).

Two framing risks:

- **Reconciling the linear early window with the published sublinear average.** Mitchell 2022's 56% sublinearity over 0–250 ms has to be reconciled with "linear in the first 20 ms". Both can be true, but the paper should say so explicitly.
- **Albrecht et al. 2002** state that monocular contrast normalization is "fully expressed within the first 10 ms after the onset of the response" (Q2, Q4). Reviewers will ask why interocular normalization is slower. Truchard et al. 2000 (monocular gain control) and Sengpiel & Vorobyov 2005 (cortical interocular suppression) together support the answer: monocular gain control is early and interocular normalization cortical and later.

Blake Mitchell's Vanderbilt dissertation, *Neuronal Basis of V1 Binocular Normalization* (ir.vanderbilt.edu/handle/1803/19176), turned up in search. I did not check its content. If it contains the time-course result, it counts as prior disclosure.

---

## 1. Binocular vs monocular response latency

**Short answer.** Yes, but mostly in human VEPs, where binocular P100 and N75 peaks come a few ms earlier than monocular ones. The best-quantified value is **2.9 ms** in focus, rising to **8.8 ms** with +2 D blur. The effect is not universal: two VEP studies found no significant difference. Binocular reaction times are also faster. In cat, a 1960 study reported shorter binocular latencies to diffuse flashes. In macaque V1, the only direct single-unit comparison found is Romero et al. 2007: **70.7 ms binocular vs 78.2 ms dominant eye** (7.5 ms, not significant). None of the Maier-lab papers, including Cox et al. 2019, report binocular vs monocular onset latency. No MEG study comparing binocular and monocular latency was found.

**Human VEP**

- **Plainis et al. 2011** (human, n = 13; 10′ checks, 100% contrast, 4 reversals/s): binocular P100 was earlier than monocular in every case. The mean difference was **2.9 ms in focus** and grew with defocus to **8.8 ms at +2.00 D**. Binocular amplitude enhancement went from 2.1 to 3.1 over the same range. This is the clearest number. It also shows that the **binocular latency advantage grows as monocular signals weaken**, the same direction as a contrast-latency account predicts.
- **Knierim, Emser & Schimrigk 1985** (human, n = 100): P100 latencies were "significantly shorter" for binocular stimulation. The abstract gives no numbers.
- **McCulloch & Skarf 1991** (human infants, n = 161): "Monocular P-VEPs have slightly longer latencies than the binocular P-VEPs". The difference did not change with age and was larger with larger checks. The abstract gives no number.
- **di Summa et al. 1997** (human, n = 14): N75 and P100 both peaked earlier binocularly. N75 shortened more than P100 with larger checks, and N145 did not change. The authors argue that the latency shortening is "not consistent with" summing two separate monocular signals.
- **Null results.** Tobimatsu & Kato 1996 (n = 13, sinusoidal gratings at 90% contrast) found "no significant effects of binocular stimulation on latencies of N70 or P100". di Summa et al. 1999 (n = 10, mild refractive errors) found no difference between binocular and *better-eye* monocular latencies. **Whether the binocular latency is compared with the faster eye or with the average of the eyes matters.**
- A search snippet stated that "monocular P100 peak times are typically 0–4 ms later than binocular" (QUT eprint 49198). I could not open the source, so this range is unverified.

**Behaviour.** Blake, Martens & Di Gianfilippo 1980 (human) found binocular reaction times to gratings "consistently faster than monocular RT" across contrast, including where RT had reached its asymptote. The advantage exceeded probability summation in observers with good stereopsis.

**Cat**

- Burns, Heron & Grafstein 1960 (cat; surface potentials and units; diffuse light) report that latency was shortened both by binocular stimulation and by brighter stimuli. The Crossref abstract gives no ms values.
- Minke & Auerbach 1972 (*Exp Brain Res* 14:409–422) compared single-unit and VEP latencies after monocular and binocular stimulation in cat striate cortex. Only the bibliographic record was accessible. **Content not verified.**
- Olianezhad et al. 2024 (cat, n = 719 sites) report inter-eye latency mismatches. Contralateral latency was 35.9 ± 7.6 ms and ipsilateral 40.0 ± 8.8 ms; the absolute mismatch per site was 5.9 ± 5.7 ms (4–9 ms across samples). They did not compare binocular with monocular latency. The paper is relevant for analysis because **a binocular onset 1–2 ms earlier should be tested against the faster eye's onset**, not the mean of the two eyes.

**Macaque**

- **Romero et al. 2007** (macaque V1, n = 35, dynamic random-dot stimuli; the abstract does not state the preparation): mean latency was 97.9 ms for the non-dominant eye, 78.2 ms for the dominant eye and **70.7 ms for binocular stimulation**. The dominant-vs-binocular difference was **not significant**. The authors conclude that latency is "largely determined by the dominant eye".
- **Maier lab.**
  - Cox et al. 2019 and Mitchell et al. 2022 report no binocular vs monocular onset latency.
  - Mitchell et al. 2022's "onset of facilitation" is the time of peak facilitation, not a latency.
  - Mitchell et al. 2023 analysed 50–200 ms averages only.
  - Dougherty et al. 2019 computed onset from dominant-eye responses only. They do report modulation onset relative to visual onset, by layer: L4 42.1 vs 50.1 ms, L2/3 45.2 vs 47.6 ms, L5/6 46.4 vs 54.7 ms.
- **LGN.** Schroeder et al. 1990 (mac-W, flash MUA) found binocular suppression (15–70%) and facilitation in the LGN within the retinally driven window (15–65 ms), sometimes before excitation. Latency was not compared.
- **MEG.** Not found.

---

## 2. Contrast-dependent latency in V1 and earlier

**Short answer.** For "stronger stimuli have shorter latencies", cite **Gawne, Kjaer & Richmond 1996** (awake macaque V1: latency mainly tracks contrast) and **Albrecht et al. 2002** (macaque and cat V1, brief *stationary* gratings like ours: the response is scaled and shifted earlier with contrast). Add **Carandini, Heeger & Movshon 1997** for the mechanism: normalization accounts for "the decrease in response gain and latency observed at high contrasts". Add **Shapley & Victor 1978** for the retinal origin of the contrast-dependent phase advance. Carandini & Heeger 1994 is a weaker choice for latency because its abstract does not mention latency or phase.

| Reference | Prep | Finding | Magnitude |
|---|---|---|---|
| Shapley & Victor 1978 | cat retina (X, Y ganglion cells) | Contrast gain control: higher contrast speeds the dynamics, with "a marked phase advance" of the components above 4 Hz; larger in Y cells | qualitative in abstract |
| Dean & Tolhurst 1986 | cat simple cells | Response phase advances with contrast | 20–35° for a 5× contrast increase |
| Albrecht 1995 | monkey and cat simple cells | Phase advance with contrast, saturating like amplitude; from changes in both latency and temporal receptive-field shape | ≈45 ms (≈¼ cycle near 5 Hz) |
| Gawne, Kjaer & Richmond 1996 | mac-W V1 complex cells (n = 37) | Response strength codes orientation; "latency is more a function of the stimulus contrast" | ms range not accessible (full text blocked) |
| Carandini, Heeger & Movshon 1997 | mac-A simple cells | Normalization (divisive shunting) explains the lower gain and **latency at high contrast** | model fits |
| Reich, Mechler & Victor 2001 | macaque V1, static stimuli | The temporally coded contrast information is carried most strongly by **response-onset latency** | information measure |
| Albrecht et al. 2002 | macaque and cat V1; 200 ms stationary gratings | ~95% of the PSTH variance across contrast is captured by **scaling and shifting**; the CRF (c50, n) is constant from the first 10 ms | the 65.3 ms shift across the contrast range is cited by Hu et al. 2011, not verified at source |
| Sit et al. 2009 | mac-W, VSDI | "response onset depends on stimulus contrast"; peak spatial profile and offset do not | population |
| Hu, Wang & Wang 2011 | cat | Peak latency increases as contrast falls; the CRF changes over the response (contrast gain early, response gain later) | ~10 ms shift with transient stimuli |
| Zhou, Benson, Kay & Winawer 2019 | human ECoG, fMRI, macaque spikes (model) | A delayed divisive gain-control model reproduces responses that "are slower at low stimulus contrast" | model |

Not recommended for the latency claim:

- **Sclar, Maunsell & Lennie 1990**: the abstract covers contrast sensitivity only.
- **Bair et al. 2002**: onset vs offset latency, not contrast.
- **Maunsell et al. 1999**: magno vs parvo timing; the abstract states no latency–intensity relation.

**Scale for the manuscript.** Contrast-latency slopes in V1 run to tens of ms across the contrast range. A binocular advantage of 1–2 ms therefore corresponds to a modest equivalent increase in monocular contrast. Stating it that way, as the monocular contrast that gives the same latency, would tie the latency result to the "stronger stimulus" argument.

---

## 3. Duration dependence of binocular combination (psychophysics and human neural)

**Short answer.** Yes, and in the direction the paper needs.

- Binocular summation is largest for brief, low-contrast stimuli and disappears by about 100 ms in an orientation task (Bearse & Freeman 1994).
- Binocular combination is "more linear" at low contrast and short duration (117 ms). At high contrast and long duration (1 s), perceived contrast is constant across interocular ratios, which is ocularity invariance (Ding, Klein & Levi 2013).
- Orthogonal dichoptic gratings appear fused below about 150 ms and rival above it (Wolfe 1983).
- Interocular masking grows with duration and is temporally low-pass. Monocular suppression is the transient, fast one (Petrov & McKee 2009; Zhou et al. 2014; Huang, Baker & Hess 2012).

I found no psychophysical study that measures the time course of **ocularity invariance** directly as a function of duration. Ding et al. 2013 is the nearest, with two durations.

**Binocular summation vs duration**

- **Bearse & Freeman 1994** (human; orientation discrimination with 1-D DOG stimuli). With brief exposures, summation was **66% at 8% contrast** and fell to about zero above 15%. At low contrast it was greatest at **50 ms**, smaller at longer durations and **not significant at ≥100 ms**. The authors attribute the loss to saturation after binocular combination.
- **Baker et al. 2018** (meta-analysis, 65 studies). Threshold summation exceeds √2 and depends on spatial and temporal frequency: slower stimulus speeds give more summation. Duration is not singled out in the abstract.

**Contrast matching and phase combination vs duration**

- **Ding, Klein & Levi 2013** (human). Under "normal viewing conditions (high contrast and long stimulus duration), perceived contrast is constant, independent of the interocular contrast ratio". However, "at low contrasts and short stimulus durations, binocular combination is more linear and contrast summation is phase-dependent". The two durations were **1 s and 117 ms**, and the authors describe 117 ms as too brief for vergence eye movements. This is the closest psychophysical analogue to "linear first, normalized later, more so at high contrast".
- **Ding & Levi 2017** (human, brightness). Weighting between the eyes differed between 200 ms and 1000 ms presentations, with more unequal weighting at 1000 ms. Levelt's (1965) more asymmetric results are attributed to his long adjustment-task durations (Levelt is cited only).
- **Legge & Rubin 1981** and **Baker, Meese & Georgeson 2007** are the standard contrast-matching references. Their records were verified, but I did not check their stimulus durations.

**Interocular suppression vs duration**

- **Wolfe 1983** (human). Dichoptic orthogonal gratings "appear fused if flashed for less than 150 ms". Pairs of flashes separated by less than ~150 ms produce rivalry, and intervals over 150 ms restore abnormal fusion.
- **Carter & Cavanagh 2007** add a caveat. With brief intermittent presentations, rivalry shows strong, location-specific onset biases rather than fusion. Balanced sustained rivalry appears only when presentations exceed about 1 s.
- **Petrov & McKee 2009** (human surround masking). There are two components: a **strong, transient, largely monocular** one, which peaked at ~20 ms durations, began within one frame (12 ms) and fell sharply beyond ~100 ms; and a **weaker, sustained, binocular (interocular)** one, about 3× weaker, roughly constant with duration and dominant beyond ~150 ms. This is surround suppression, so it is indirect evidence for overlay interocular normalization.
- **Zhou et al. 2014** (normals and amblyopes). Dichoptic masking is "low-pass temporally": longer target presentations gave more threshold elevation, and masking was greatest at the longest durations.
- **Huang, Baker & Hess 2012**. Dichoptic masking measured with a probed-sinewave paradigm has low-pass spatiotemporal properties.
- **Meese & Baker 2009**. Monocular cross-orientation suppression depends on speed (TF/SF). Interocular suppression is speed-invariant, which fits a slower, cortical route.
- **van Boxtel et al. 2007** report that dichoptic masking and rivalry share perceptual dynamics. **Kaunitz et al. 2014** report that target duration and onset timing govern interocular suppression of transient targets.
- **Macknik & Martinez-Conde 2004** (macaque V1 units and human psychophysics, "standing wave of invisibility"). V1 binocular neurons are "dichoptic in excitation but monoptic in inhibition": a monocular target's response is suppressed strongly only by masks in the *same* eye, while interocular masking reaches a small subset of cells (as summarized by Petrov & McKee 2009). The abstract gives no SOA functions, and I could not reach the full text, so **the dichoptic-masking time-course numbers are not verified**. For this paper, it shows that strong interocular masking is not an early-V1 property under that paradigm.

**Human neural time course**

- **Moradi & Heeger 2009** (fMRI). Linear summation fails at ≥10% contrast, and interocular normalization fits. A 500 ms interocular onset asynchrony "had no impact". fMRI cannot resolve milliseconds.
- **Chadnova et al. 2014** (MEG, conference abstract, n = 4). Dichoptic suppression is present "as early as the primary visual cortex". No latencies are given.
- I found no EEG or MEG study that resolves when binocular summation gives way to interocular suppression on a millisecond scale.

---

## 4. Neurophysiology of the time course of binocular combination and normalization

**Short answer.** Beyond Cox et al. 2019, the time course has been addressed by **Dougherty et al. 2019** (facilitation at onset, suppression ~12 ms later in monocular neurons), **Mitchell et al. 2022** (early facilitation, then contrast-dependent suppression at medium and high contrast; **closest to this paper**), **Carlson et al. 2023** (dioptic facilitation only in 50–150 ms; dichoptic suppression at 150–250 ms) and **Tanabe & Cumming 2014** (disparity suppression lags excitation by a median of 7 ms). See "Novelty" at the top. No study tests B = L + R at millisecond resolution or times the onset of interocular normalization as a function of contrast.

**Macaque V1 binocular time course**

- **Cox et al. 2019** (mac-W, laminar arrays; ~0.45/0.9 contrast; 500 ms stimuli; zero disparity). Transient (50–100 ms) facilitation for both congruent and incongruent stimuli: congruent Δ = 13.7 spikes/s, incongruent Δ = 5.9 spikes/s. In the sustained window (150–250 ms), congruent facilitation is gone and incongruent stimuli are suppressed (Δ = −5.6 spikes/s). The 100–150 ms window was not analysed.
- **Dougherty et al. 2019** (mac-W, monocular neurons, all layers). In L4, facilitated units had a 0.25-quantile modulation onset of 49 ms against 54 ms for suppressed units. The authors attribute facilitation to "additive synaptic activation" by convergent thalamocortical input and suppression to intracortical inhibition. They never use the word "normalization".
- **Mitchell et al. 2022** (mac-W, arrays; contrast groups ~0.2, 0.45, 0.9; equal contrasts in the two eyes; 250 ms). Analyses used 100 ms windows sliding in 10 ms steps, plus fixed 50–150, 100–200 and 150–250 ms windows. Late suppression fitted response-gain better than contrast-gain control. No latency analysis.
- **Mitchell et al. 2023** (mac-W, 6 contrasts 0–0.9 per eye, 250 ms). Interocular normalization weighted by ocular dominance fits best. Analyses used 50–200 ms averages, with no time-resolved analysis.
- **Carlson et al. 2023** (mac-W). Dioptic facilitation in 50–150 ms, not significant in 150–250 ms. They state that "the response suppression from dichoptic stimulation is typically observed 150–250 ms after stimulus onset."
- **Tanabe & Cumming 2014** (mac-W, binocular subspace mapping). "The suppression was delayed relative to the excitation, with a median difference in latency of 7 ms." The suppression concerns the anti-preferred interocular phase at lower spatial frequencies. This is a disparity mechanism, not contrast normalization, but it is the shortest published binocular excitation-to-suppression delay.
- **Samonds, Potetz & Lee 2009** (macaque V1). Competing disparities decorrelate "soon after response onset", and disparity tuning improves over time.
- **Zhang et al. 2024** (mac-W two-photon, superficial layers, ~10,000 neurons). Binocular combination depends on ocular dominance (OD-dependent divisive interocular inhibition plus summation). Calcium kinetics cannot resolve these timescales.

**Cat and mouse**

- **Menz & Freeman 2003, 2004** (cat). Disparity tuning sharpens during the response, coarse to fine (see Q5).
- **Sengpiel, Blakemore & Harrad 1995; Sengpiel & Blakemore 1994** (cat). Interocular suppression by an orthogonal grating appears only if the cell is already responding to the other eye.
- **Sengpiel et al. 1998** (cat). Interocular suppression acts mainly as **response-gain** control, with a contrast-gain component. Cox et al. 2019 cite this paper for suppression "after ~150 ms" in more than half of neurons, but the abstract does not mention timing (**cited only**).
- **Sengpiel & Vorobyov 2005** (cat; species per Meese & Baker 2009). Interocular suppression is cortical and GABAergic: it cannot be evoked by fast stimuli, adaptation abolishes it, and bicuculline reduces it.
- **Truchard, Ohzawa & Freeman 2000** (cat). Contrast-gain reductions occur "primarily at a monocular site, before convergence", while binocular gain changes are modest.
- **Longordo et al. 2013** (mouse, whole-cell L2/3). Binocular integration is sublinear, amplitude-dependent and needs balanced inhibition. I could not access the time course.
- **Perumal et al. 2026** (mouse, bioRxiv v2, 24 Sep 2026; Stuart lab). Pyramidal and PV binocular summation was sublinear "throughout the visually evoked response". SST summation was supralinear early (<200 ms) and late (>250 ms), with sublinear summation in between. This is **not** linear-then-sublinear in pyramidal cells and does not threaten novelty.

**General normalization dynamics (non-binocular), for the Discussion**

- **Albrecht et al. 2002**. Contrast normalization is "present in full strength, virtually instantaneously, at the onset of the response". This is the main counterpoint.
- **Hu et al. 2011** (cat). The CRF does change over the response.
- **Smith, Bair & Movshon 2006** (macaque). Cross-orientation suppression is "very fast, starting even before the response to optimal excitatory stimuli".
- **Bair, Cavanaugh & Movshon 2003** (macaque). Surround suppression is delayed, and its latency depends on strength: strong suppression arrives ~30 ms earlier than weak.
- **Models with a delayed or low-pass normalization signal:** Zhou et al. 2019; Heeger & Zemlianova 2020; Gulati & Ray 2025 (preprint).
- **Zhou, Whitmire, Chen & Seidemann 2023** (bioRxiv, mac-W VSDI and GCaMP) found V1 responses "surprisingly near additive in time". This is a possible counterpoint for population signals.

---

## 5. Fixation disparity in trained macaques, interocular phase and disparity dynamics

**Short answer.**

- **Size of vergence noise.** Coil-based estimates in fixating macaques put the horizontal eye-position SD at **6–8 arcmin per eye** (Motter & Poggio 1984). The eyes were misaligned by more than 7 arcmin for ~60% of fixation time. However, Read & Cumming 2003 showed that search-coil error (~0.1° SD) accounts for this apparent vergence noise. From V1 responses they estimate the **true vergence SD at ≤1 arcmin**, like the 1–2 arcmin reported in humans.
- **Systematic offset.** I found **no measurement of mean fixation disparity in macaques viewing through a mirror stereoscope** (not found).
- **Effect on interocular phase.** A horizontal vergence error d (deg) changes interocular phase by Δφ = 360°·f·d·cos α for a grating of spatial frequency f (cpd) whose stripes are tilted α from vertical. At the default 0.9 cpd this is small for plausible errors. At 4 cpd (session 211008) it can reach the half-cycle range where phase-tuned V1 cells switch from facilitation to suppression.
- **Disparity dynamics.** Disparity selectivity sharpens during the first tens of ms (coarse to fine), and disparity suppression lags excitation by ~7 ms. A vergence-induced phase offset could therefore produce a reduction that also grows over the response. Contrast dependence argues against this explanation: Smith et al. 1997 found that phase tuning is **contrast-invariant**.

**Vergence and fixation disparity**

- **Motter & Poggio 1984** (mac-W, coils). Eye-position SD was 6–8 arcmin horizontally and 7–13 arcmin vertically during binocular fixation. "For nearly 60% of eye pause time the eyes were misaligned on the fixation target by more than 7 min arc."
- **Read & Cumming 2003** (mac-W). Coil error is "approximately the same magnitude as the small eye movements which occur during fixation: SD approximately 0.1°", "large enough to explain the SD of measured vergence". From disparity-tuning variability they estimate "the SD of the vergence angle in these two monkeys is ≤ 1 min arc" (per-neuron mean 0.011°, ≈0.7 arcmin).
- **Cumming & Parker 1999** (mac-W). Human fixation-disparity SDs are 1–2 arcmin, and Motter & Poggio's values "were at least fourfold greater". A 1% eye-calibration gain error produces an apparent fixation disparity above 0.025° (1.5 arcmin). This is relevant because the paper's eye tracking cannot settle the question either.
- **Short-latency vergence.** Busettini, Fitzgibbon & Miles 2001 measured disparity vergence at 80.9 ± 3.9 ms in humans and about 25 ms earlier in monkeys (~55 ms); see also Masson, Busettini & Miles 1997. A grating could therefore trigger corrective vergence only from ~55 ms after stimulus onset, with negligible movement in the first ~20 ms of the V1 response. The vergence state that matters is the one present at stimulus onset.

**Interocular phase from vergence error** (derived; Δφ = 360°·f·d for vertical stripes, scaled by cos α otherwise)

| Vergence error | 0.9 cpd (default) | 4 cpd (211008) |
|---|---|---|
| 1 arcmin (Read & Cumming SD) | 5° | 24° |
| 5 arcmin | 27° | 120° |
| 10 arcmin (Motter & Poggio range) | 54° | 240° (= −120°) |

**Phase tuning in V1**

- **Smith, Chino, Ni, Ridder & Crawford 1997** (mac-A, drifting dichoptic gratings). About 70% of simple cells and 40% of complex cells were phase-selective. Simple cells "typically demonstrated binocular facilitation at the optimal phase disparity and binocular suppression for disparities 180 degrees away." Crucially for the Discussion, "stimulus contrast had relatively little effect on a cell's disparity tuning, because response saturation occurred at the same contrast level for all relative interocular phase disparities." A phase offset caused by vergence should therefore scale responses by a contrast-independent factor. It should not produce a suppression that strengthens with contrast.
- **Smith, Chino, Ni & Cheng 1997** (mac-A). Near threshold, the two eyes' contrasts combine linearly before rectification for cooperative and antagonistic phases alike.
- Cox et al. 2019, Dougherty et al. 2019 and Zhang et al. 2024 all used zero or near-zero disparity. eLife reviewers raised the near-zero disparity of Zhang et al. as a limitation.

**Time course of disparity selectivity**

- **Menz & Freeman 2003** (cat). "The disparity tuning of cortical cells generally sharpened during the time course of response". Cells tuned to large spatial scales had shorter latencies than those tuned to fine scales. The ms values were not accessible (paywall).
- **Menz & Freeman 2004** (cat). Over a short time span, disparity frequency increases and disparity range decreases while optimal disparity stays constant (coarse to fine). LGN and simple-cell RF dynamics do not fully predict this.
- **Samonds, Potetz & Lee 2009** (macaque V1). Disparity tuning sharpens over time through local competition and distant cooperation.
- **Tanabe & Cumming 2014** (mac-W). Suppression lags excitation by a median of 7 ms (see Q4).

---

## 6. Citation verification

| Requested | Status | Correct details |
|---|---|---|
| Smith, Chino, Ni, Cheng, Crawford & Harwerth 1997, *J Neurophysiol*, binocular combination of contrast signals | **Author list wrong** | The paper has **four** authors: Smith EL 3rd, Chino Y, Ni J, Cheng H (1997). *Binocular combination of contrast signals by striate cortical neurons in the monkey.* J Neurophysiol 78(1):366–382. doi:10.1152/jn.1997.78.1.366 (PubMed 9242286). The six-author list belongs to a different paper: Smith, Chino, Ni, Cheng, Crawford & Harwerth (1997). *Residual binocular interactions in the striate cortex of monkeys reared with abnormal binocular vision.* J Neurophysiol 78(3):1353–1362. doi:10.1152/jn.1997.78.3.1353 (Crossref). |
| Smith et al. 1997, binocular spatial phase tuning | Verified | Smith EL 3rd, Chino YM, Ni J, Ridder WH 3rd, Crawford MLJ (1997). *Binocular spatial phase tuning characteristics of neurons in the macaque striate cortex.* J Neurophysiol 78(1):351–365. doi:10.1152/jn.1997.78.1.351 |
| Macknik & Martinez-Conde 2004 | Verified | J Cogn Neurosci 16(6):1049–1059. doi:10.1162/0898929041502788 |
| Legge & Rubin 1981 | Verified | Percept Psychophys 30(1):49–61. doi:10.3758/BF03206136 |
| Baker, Meese & Georgeson 2007 | Verified | Spat Vis 20(5):397–413. doi:10.1163/156856807781503622. PubMed gives the order Baker, Meese, Georgeson. **Crossref metadata lists it reversed** (Georgeson, Meese, Baker), so reference managers may import it wrongly. |
| Menz & Freeman 2003 | Verified | Nat Neurosci 6(1):59–65. doi:10.1038/nn986. Published online Dec 2002, so some databases show 2002. |
| Davis et al. 2023 eLife | Verified | Davis ZW, Dotson NM, Franken TP, Muller L, Reynolds JH. *Spike-phase coupling patterns reveal laminar identity in primate cortex.* eLife 12:e84512. doi:10.7554/eLife.84512 |
| Vinck et al. 2010 | Verified | Vinck M, van Wingerden M, Womelsdorf T, Fries P, Pennartz CMA. *The pairwise phase consistency: a bias-free measure of rhythmic neuronal synchronization.* NeuroImage 51(1):112–122. doi:10.1016/j.neuroimage.2010.01.073 |
| Zhang et al. 2024 eLife | Verified | Zhang SH, Zhao XN, Jiang DQ, Tang SM, Yu C. *Ocular dominance-dependent binocular combination of monocular neuronal responses in macaque V1.* eLife 13:RP92839. doi:10.7554/eLife.92839. Two-photon, three awake macaques, superficial layers. A search summary of the eLife review page said a reviewer flagged Mitchell et al. 2023 as having an "identical purpose" and a "highly similar conclusion", and that the authors then cited it. I did not read the review page directly. |
| Mitchell et al. 2023 Curr Biol | Verified | Mitchell BA, Carlson BM, Westerberg JA, Cox MA, Maier A. *A role for ocular dominance in binocular integration.* Curr Biol 33(18):3884–3895.e5. doi:10.1016/j.cub.2023.08.019 |
| Mitchell et al. 2022 iScience | Verified | Mitchell BA, Dougherty K, Westerberg JA, Carlson BM, Daumail L, Maier A, Cox MA. *Stimulating both eyes with matching stimuli enhances V1 responses.* iScience 25(5):104182. doi:10.1016/j.isci.2022.104182 |
| Cox et al. 2019 J Vis | Verified | Cox MA, Dougherty K, Westerberg JA, Schall MS, Maier A. *Temporal dynamics of binocular integration in primary visual cortex.* J Vis 19(12):13. doi:10.1167/19.12.13 |

Smaller inconsistencies noted inside papers (worth knowing if you quote them):

- **Cox et al. 2019.** Table 1 lists 30 penetrations, but Fig 1D gives n = 41. The Fig 3 caption gives −0.03 where the text gives −0.06.
- **Mitchell et al. 2023.** The Fig 1 legend lists the lowest nonzero contrast as 0.05, and the Methods list 0.055.
- **Zhou et al. 2014.** PubMed spells the second author McNeal and Crossref spells it McNeill.

---

## References

Format: Authors (Year). Title. *Journal* Volume(Issue):pages. DOI. [status]. "PubMed" or "PMC" means the abstract or full text was read at that source. "Crossref" means bibliographic metadata only.

1. Albrecht DG (1995). Visual cortex neurons in monkey and cat: effect of contrast on the spatial and temporal phase transfer functions. *Vis Neurosci* 12(6):1191–1210. doi:10.1017/S0952523800006817. [verified, PubMed]
2. Albrecht DG, Geisler WS, Frazor RA, Crane AM (2002). Visual cortex neurons of monkeys and cats: temporal dynamics of the contrast response function. *J Neurophysiol* 88(2):888–913. doi:10.1152/jn.2002.88.2.888. [verified, PubMed; the 65.3 ms latency shift is cited only, via Hu et al. 2011]
3. Bair W, Cavanaugh JR, Movshon JA (2003). Time course and time-distance relationships for surround suppression in macaque V1 neurons. *J Neurosci* 23(20):7690–7701. doi:10.1523/JNEUROSCI.23-20-07690.2003. [verified, PubMed]
4. Baker DH, Lygo FA, Meese TS, Georgeson MA (2018). Binocular summation revisited: beyond √2. *Psychol Bull* 144(11):1186–1199. doi:10.1037/bul0000163. [verified, PubMed]
5. Baker DH, Meese TS, Georgeson MA (2007). Binocular interaction: contrast matching and contrast discrimination are predicted by the same model. *Spat Vis* 20(5):397–413. doi:10.1163/156856807781503622. [verified, PubMed]
6. Bearse MA Jr, Freeman RD (1994). Binocular summation in orientation discrimination depends on stimulus contrast and duration. *Vision Res* 34(1):19–29. doi:10.1016/0042-6989(94)90253-4. [verified, PubMed]
7. Blake R, Martens W, Di Gianfilippo A (1980). Reaction time as a measure of binocular interaction in human vision. *Invest Ophthalmol Vis Sci* 19(8):930–941. No DOI (PMID 7409986). [verified, PubMed]
8. Burns BD, Heron W, Grafstein B (1960). Responses of cerebral cortex to diffuse monocular and binocular stimulation. *Am J Physiol* 198(1):200–204. doi:10.1152/ajplegacy.1960.198.1.200. [verified, Crossref abstract]
9. Busettini C, Fitzgibbon EJ, Miles FA (2001). Short-latency disparity vergence in humans. *J Neurophysiol* 85(3):1129–1152. doi:10.1152/jn.2001.85.3.1129. [verified, PubMed]
10. Carandini M, Heeger DJ (1994). Summation and division by neurons in primate visual cortex. *Science* 264(5163):1333–1336. doi:10.1126/science.8191289. [verified, PubMed]
11. Carandini M, Heeger DJ, Movshon JA (1997). Linearity and normalization in simple cells of the macaque primary visual cortex. *J Neurosci* 17(21):8621–8644. doi:10.1523/JNEUROSCI.17-21-08621.1997. [verified, PubMed]
12. Carlson BM, Mitchell BA, Dougherty K, Westerberg JA, Cox MA, Maier A (2023). Does V1 response suppression initiate binocular rivalry? *iScience* 26(8):107359. doi:10.1016/j.isci.2023.107359. [verified, PubMed + PMC]
13. Carter O, Cavanagh P (2007). Onset rivalry: brief presentation isolates an early independent phase of perceptual competition. *PLoS One* 2(4):e343. doi:10.1371/journal.pone.0000343. [verified, PubMed]
14. Chadnova E, Reynaud A, Clavagnier S, Baker DH, Baillet S, Hess RF (2014). Dynamics of dichoptic masking in the primary visual cortex. *BMC Neurosci* 15(Suppl 1):P145. doi:10.1186/1471-2202-15-S1-P145. [verified, conference abstract]
15. Cox MA, Dougherty K, Westerberg JA, Schall MS, Maier A (2019). Temporal dynamics of binocular integration in primary visual cortex. *J Vis* 19(12):13. doi:10.1167/19.12.13. [verified, PubMed + PMC]
16. Cumming BG, Parker AJ (1999). Binocular neurons in V1 of awake monkeys are selective for absolute, not relative, disparity. *J Neurosci* 19(13):5602–5618. doi:10.1523/JNEUROSCI.19-13-05602.1999. [verified, PMC]
17. Davis ZW, Dotson NM, Franken TP, Muller L, Reynolds JH (2023). Spike-phase coupling patterns reveal laminar identity in primate cortex. *eLife* 12:e84512. doi:10.7554/eLife.84512. [verified, PubMed]
18. Dean AF, Tolhurst DJ (1986). Factors influencing the temporal phase of response to bar and grating stimuli for simple cells in the cat striate cortex. *Exp Brain Res* 62(1):143–151. doi:10.1007/BF00237410. [verified, PubMed]
19. di Summa A, Polo A, Tinazzi M, Zanette G, Bertolasi L, Bongiovanni LG, Fiaschi A (1997). Binocular interaction in normal vision studied by pattern-reversal visual evoked potential (PR-VEPS). *Ital J Neurol Sci* 18(2):81–86. doi:10.1007/BF01999567. [verified, PubMed]
20. di Summa A, Fusina S, Bertolasi L, Vicentini S, Perlini S, Bongiovanni LG, Polo A (1999). Mechanism of binocular interaction in refraction errors: study using pattern-reversal visual evoked potentials. *Doc Ophthalmol* 98(2):139–151. doi:10.1023/A:1002190127573. [verified, PubMed]
21. Ding J, Klein SA, Levi DM (2013). Binocular combination of phase and contrast explained by a gain-control and gain-enhancement model. *J Vis* 13(2):13. doi:10.1167/13.2.13. [verified, PubMed + PMC]
22. Ding J, Levi DM (2017). Binocular combination of luminance profiles. *J Vis* 17(13):4. doi:10.1167/17.13.4. [verified, PubMed + PMC]
23. Dougherty K, Cox MA, Westerberg JA, Maier A (2019). Binocular modulation of monocular V1 neurons. *Curr Biol* 29(3):381–391.e4. doi:10.1016/j.cub.2018.12.004. [verified, PubMed + PMC]
24. Gawne TJ, Kjaer TW, Richmond BJ (1996). Latency: another potential code for feature binding in striate cortex. *J Neurophysiol* 76(2):1356–1360. doi:10.1152/jn.1996.76.2.1356. [verified, PubMed; full text not accessible]
25. Gulati D, Ray S (2025). Slow-varying normalization explains diverse temporal frequency masking interactions across scales in the macaque visual cortex. *bioRxiv* doi:10.1101/2025.03.11.642541 (v2). [verified, bioRxiv API; preprint]
26. Heeger DJ, Zemlianova KO (2020). A recurrent circuit implements normalization, simulating the dynamics of V1 activity. *Proc Natl Acad Sci USA* 117(36):22494–22505. doi:10.1073/pnas.2005417117. [verified, PubMed]
27. Hu M, Wang Y, Wang Y (2011). Rapid dynamics of contrast responses in the cat primary visual cortex. *PLoS One* 6(10):e25410. doi:10.1371/journal.pone.0025410. [verified, PMC]
28. Huang PC, Baker DH, Hess RF (2012). Interocular suppression in normal and amblyopic vision: spatio-temporal properties. *J Vis* 12(11):29. doi:10.1167/12.11.29. [verified, PubMed]
29. Kaunitz LN, Fracasso A, Skujevskis M, Melcher D (2014). Waves of visibility: probing the depth of inter-ocular suppression with transient and sustained targets. *Front Psychol* 5:804. doi:10.3389/fpsyg.2014.00804. [verified, PubMed]
30. Knierim A, Emser W, Schimrigk K (1985). [Latency differences of visual evoked potentials in binocular and monocular checkerboard stimulation as well as age and sex dependence] (in German). *EEG EMG Z Elektroenzephalogr Elektromyogr Verwandte Geb* 16(4):212–214. PMID 3935416; DOI 10.1055/s-2008-1060985 taken from the Thieme URL, not from PubMed. [verified, PubMed]
31. Legge GE, Rubin GS (1981). Binocular interactions in suprathreshold contrast perception. *Percept Psychophys* 30(1):49–61. doi:10.3758/BF03206136. [verified, PubMed]
32. Longordo F, To MS, Ikeda K, Stuart GJ (2013). Sublinear integration underlies binocular processing in primary visual cortex. *Nat Neurosci* 16(6):714–723. doi:10.1038/nn.3394. [verified, PubMed]
33. Macknik SL, Martinez-Conde S (2004). Dichoptic visual masking reveals that early binocular neurons exhibit weak interocular suppression: implications for binocular vision and visual awareness. *J Cogn Neurosci* 16(6):1049–1059. doi:10.1162/0898929041502788. [verified, PubMed; full text not accessible]
34. Masson GS, Busettini C, Miles FA (1997). Vergence eye movements in response to binocular disparity without depth perception. *Nature* 389(6648):283–286. doi:10.1038/38496. [verified, PubMed]
35. McCulloch DL, Skarf B (1991). Development of the human visual system: monocular and binocular pattern VEP latency. *Invest Ophthalmol Vis Sci* 32(8):2372–2381. No DOI (PMID 2071348). [verified, PubMed]
36. Meese TS, Baker DH (2009). Cross-orientation masking is speed invariant between ocular pathways but speed dependent within them. *J Vis* 9(5):2. doi:10.1167/9.5.2. [verified, PubMed]
37. Menz MD, Freeman RD (2003). Stereoscopic depth processing in the visual cortex: a coarse-to-fine mechanism. *Nat Neurosci* 6(1):59–65. doi:10.1038/nn986. [verified, PubMed]
38. Menz MD, Freeman RD (2004). Temporal dynamics of binocular disparity processing in the central visual pathway. *J Neurophysiol* 91(4):1782–1793. doi:10.1152/jn.00571.2003. [verified, PubMed]
39. Minke B, Auerbach E (1972). Latencies and correlation in single units and visual evoked potentials in the cat striate cortex following monocular and binocular stimulations. *Exp Brain Res* 14(4):409–422. doi:10.1007/BF00235036. [record only, PubMed + Crossref; content not seen]
40. Mitchell BA, Dougherty K, Westerberg JA, Carlson BM, Daumail L, Maier A, Cox MA (2022). Stimulating both eyes with matching stimuli enhances V1 responses. *iScience* 25(5):104182. doi:10.1016/j.isci.2022.104182. [verified, PubMed + PMC]
41. Mitchell BA, Carlson BM, Westerberg JA, Cox MA, Maier A (2023). A role for ocular dominance in binocular integration. *Curr Biol* 33(18):3884–3895.e5. doi:10.1016/j.cub.2023.08.019. [verified, PubMed + PMC]
42. Moradi F, Heeger DJ (2009). Inter-ocular contrast normalization in human visual cortex. *J Vis* 9(3):13. doi:10.1167/9.3.13. [verified, PubMed]
43. Motter BC, Poggio GF (1984). Binocular fixation in the rhesus monkey: spatial and temporal characteristics. *Exp Brain Res* 54(2):304–314. doi:10.1007/BF00236231. [verified, PubMed]
44. Olianezhad F, Jin J, Najafian S, Pons C, Mazade R, Kremkow J, Alonso JM (2024). Binocular receptive-field construction in the primary visual cortex. *Curr Biol* 34(11):2474–2486.e5. doi:10.1016/j.cub.2024.04.058. [verified, PMC]
45. Perumal MB, Bushnell H, Dopp D, Gharaei S, Arabzadeh E, Nair SS, Stuart GJ (2026). Cell-type specific inhibitory dynamics shape binocular response normalization in visual cortex. *bioRxiv* doi:10.64898/2026.09.20.753031 (v2). [verified, bioRxiv; preprint]
46. Petrov Y, McKee SP (2009). The time course of contrast masking reveals two distinct mechanisms of human surround suppression. *J Vis* 9(1):21. doi:10.1167/9.1.21. [verified, PubMed + PMC]
47. Plainis S, Petratou D, Giannakopoulou T, Atchison DA, Tsilimbaris MK (2011). Binocular summation improves performance to defocus-induced blur. *Invest Ophthalmol Vis Sci* 52(5):2784–2789. doi:10.1167/iovs.10-6545. [verified, PubMed]
48. Read JCA, Cumming BG (2003). Measuring V1 receptive fields despite eye movements in awake monkeys. *J Neurophysiol* 90(2):946–960. doi:10.1152/jn.01068.2002. [verified, PubMed + PMC]
49. Reich DS, Mechler F, Victor JD (2001). Temporal coding of contrast in primary visual cortex: when, what, and why. *J Neurophysiol* 85(3):1039–1050. doi:10.1152/jn.2001.85.3.1039. [verified, PubMed]
50. Romero MC, Castro AF, Bermudez MA, Perez R, Gonzalez F (2007). Eye dominance and response latency in area V1 of the monkey. *Vis Neurosci* 24(5):757–761. doi:10.1017/S0952523807070642. [verified, PubMed]
51. Samonds JM, Potetz BR, Lee TS (2009). Cooperative and competitive interactions facilitate stereo computations in macaque primary visual cortex. *J Neurosci* 29(50):15780–15795. doi:10.1523/JNEUROSCI.2305-09.2009. [verified, PubMed]
52. Schroeder CE, Tenke CE, Arezzo JC, Vaughan HG Jr (1990). Binocularity in the lateral geniculate nucleus of the alert macaque. *Brain Res* 521(1–2):303–310. doi:10.1016/0006-8993(90)91556-V. [verified, PubMed]
53. Sengpiel F, Blakemore C (1994). Interocular control of neuronal responsiveness in cat visual cortex. *Nature* 368(6474):847–850. doi:10.1038/368847a0. [verified, PubMed]
54. Sengpiel F, Blakemore C, Harrad R (1995). Interocular suppression in the primary visual cortex: a possible neural basis of binocular rivalry. *Vision Res* 35(2):179–195. doi:10.1016/0042-6989(94)00125-6. [verified, PubMed]
55. Sengpiel F, Baddeley RJ, Freeman TC, Harrad R, Blakemore C (1998). Different mechanisms underlie three inhibitory phenomena in cat area 17. *Vision Res* 38(14):2067–2080. doi:10.1016/S0042-6989(97)00413-6. [verified, PubMed; the ~150 ms timing is cited only, via Cox et al. 2019]
56. Sengpiel F, Vorobyov V (2005). Intracortical origins of interocular suppression in the visual cortex. *J Neurosci* 25(27):6394–6400. doi:10.1523/JNEUROSCI.0862-05.2005. [verified, PubMed]
57. Shapley RM, Victor JD (1978). The effect of contrast on the transfer properties of cat retinal ganglion cells. *J Physiol* 285:275–298. doi:10.1113/jphysiol.1978.sp012571. [verified, PubMed]
58. Sit YF, Chen Y, Geisler WS, Miikkulainen R, Seidemann E (2009). Complex dynamics of V1 population responses explained by a simple gain-control model. *Neuron* 64(6):943–956. doi:10.1016/j.neuron.2009.08.041. [verified, PubMed]
59. Smith EL 3rd, Chino Y, Ni J, Cheng H (1997). Binocular combination of contrast signals by striate cortical neurons in the monkey. *J Neurophysiol* 78(1):366–382. doi:10.1152/jn.1997.78.1.366. [verified, PubMed]
60. Smith EL 3rd, Chino YM, Ni J, Ridder WH 3rd, Crawford MLJ (1997). Binocular spatial phase tuning characteristics of neurons in the macaque striate cortex. *J Neurophysiol* 78(1):351–365. doi:10.1152/jn.1997.78.1.351. [verified, PubMed]
61. Smith EL 3rd, Chino YM, Ni J, Cheng H, Crawford MLJ, Harwerth RS (1997). Residual binocular interactions in the striate cortex of monkeys reared with abnormal binocular vision. *J Neurophysiol* 78(3):1353–1362. doi:10.1152/jn.1997.78.3.1353. [Crossref only]
62. Smith MA, Bair W, Movshon JA (2006). Dynamics of suppression in macaque primary visual cortex. *J Neurosci* 26(18):4826–4834. doi:10.1523/JNEUROSCI.5542-06.2006. [verified, PubMed]
63. Tanabe S, Cumming BG (2014). Delayed suppression shapes disparity selective responses in monkey V1. *J Neurophysiol* 111(9):1759–1769. doi:10.1152/jn.00426.2013. [verified, PubMed]
64. Tobimatsu S, Kato M (1996). The effect of binocular stimulation on each component of transient and steady-state VEPs. *Electroencephalogr Clin Neurophysiol* 100(3):177–183. doi:10.1016/0168-5597(95)00273-1. [verified, PubMed]
65. Truchard AM, Ohzawa I, Freeman RD (2000). Contrast gain control in the visual cortex: monocular versus binocular mechanisms. *J Neurosci* 20(8):3017–3032. doi:10.1523/JNEUROSCI.20-08-03017.2000. [verified, PubMed]
66. van Boxtel JJA, van Ee R, Erkelens CJ (2007). Dichoptic masking and binocular rivalry share common perceptual dynamics. *J Vis* 7(14):3. doi:10.1167/7.14.3. [verified, PubMed]
67. Vinck M, van Wingerden M, Womelsdorf T, Fries P, Pennartz CMA (2010). The pairwise phase consistency: a bias-free measure of rhythmic neuronal synchronization. *NeuroImage* 51(1):112–122. doi:10.1016/j.neuroimage.2010.01.073. [verified, PubMed]
68. Wolfe JM (1983). Influence of spatial frequency, luminance, and duration on binocular rivalry and abnormal fusion of briefly presented dichoptic stimuli. *Perception* 12(4):447–456. doi:10.1068/p120447. [verified, PubMed]
69. Zhang SH, Zhao XN, Jiang DQ, Tang SM, Yu C (2024). Ocular dominance-dependent binocular combination of monocular neuronal responses in macaque V1. *eLife* 13:RP92839. doi:10.7554/eLife.92839. [verified, PubMed + eLife]
70. Zhou J, Benson NC, Kay K, Winawer J (2019). Predicting neuronal dynamics with a delayed gain control model. *PLoS Comput Biol* 15(11):e1007484. doi:10.1371/journal.pcbi.1007484. [verified, PubMed]
71. Zhou J, McNeal S, Babu RJ, Baker DH, Bobier WR, Hess RF (2014). Time course of dichoptic masking in normals and suppression in amblyopes. *Invest Ophthalmol Vis Sci* 55(7):4098–4104. doi:10.1167/iovs.14-13969. [verified, PubMed]
72. Zhou J, Whitmire M, Chen Y, Seidemann E (2023). Delayed normalization model captures disparate nonlinear neural dynamics measured with different techniques in macaque and human V1. *bioRxiv* doi:10.1101/2023.01.30.525700. [verified, bioRxiv API; preprint, no journal version listed]

**Cited only (not checked at source):**

- Levelt WJM (1965). *Br J Psychol* 56:1–13, as cited in Ding & Levi 2017.
- Ding J, Sperling G (2007). Binocular combination: measurements and a model. In Harris & Jenkin (eds), *Computational Vision in Neural and Machine Systems*, pp. 257–305, Cambridge University Press; as cited in Yehezkel et al. 2016, *R Soc Open Sci* 3:160534.
- The "0–4 ms" VEP range from QUT eprint 49198 (source not identified).

**Not found:**

- An MEG binocular vs monocular latency comparison.
- A measured mean fixation disparity for macaques viewing through a mirror stereoscope.
- A single-unit study testing B = L + R at millisecond resolution.
- A psychophysical time course of ocularity invariance beyond the two durations in Ding et al. 2013.
