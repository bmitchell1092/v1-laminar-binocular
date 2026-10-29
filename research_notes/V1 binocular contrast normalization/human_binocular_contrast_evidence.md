# Human evidence on binocular contrast combination, summation, dichoptic masking and interocular normalization: what V1 neurons should (and might not) do

Scope: human psychophysics, EEG/SSVEP, MEG and (laminar) fMRI, foundational work plus 2015–2026. The aim is to pull out quantitative benchmarks and the places where human data and V1 neurophysiology disagree or leave gaps. Throughout, "the team's design" means awake macaque V1 laminar recordings with a 6×6 L×R contrast matrix {0, .055, .11, .225, .45, .9} using static, flashed, same-orientation gratings (no frequency tagging).

Verification note: items marked *(abstract/summary only)* come from abstracts, repository records or search summaries, not from a full-text read. Full text was read (directly or by text extraction) for Meese et al. 2006, Moradi & Heeger 2009, Baker & Wade 2017, Baker et al. 2018, Hou et al. 2021 and Carter et al. 2025.

---

## 1. Binocular summation at and above threshold, contrast matching / ocularity invariance, and dichoptic vs monoptic masking (psychophysics)

### Takeaway
At detection threshold, two eyes beat one by about 1.5–1.7× (3.5–4.8 dB). That is reliably more than the classical √2 and less than linear summation (2×). Above threshold the binocular advantage disappears for discrimination and perceived contrast ("ocularity invariance"), and a dichoptic mask suppresses much more strongly than a monocular or binocular mask. The model that fits all of this best is the Meese–Georgeson–Baker two-stage gain-control model. It has near-linear monocular transduction (m ≈ 1.3), strong interocular suppression *before* binocular summation (other-eye weight = 1), and a second, steeply accelerating gain-control stage *after* summation.

### Cited Findings
**Threshold summation (foundational)**
- Legge (1984a), "Binocular contrast summation—I. Detection and discrimination," *Vision Res* 24(4):373–383: for simple detection, binocular d′ was at least twice monocular d′, and binocular summation fell as background (pedestal) contrast rose. — [Legge 1984 I PDF](https://legge.psych.umn.edu/sites/legge.psych.umn.edu/files/files/media/legge84_binocular_contrast_summation._i._detection_and_discrimination.pdf) *(abstract/summary only)*
- Legge (1984b), "Binocular contrast summation—II. Quadratic summation," *Vision Res* 24(4):385–394, proposed a binocular energy detector: monocular signals are squared before summation, which predicts a √2 threshold ratio with the same psychometric steepness in both cases. — [Legge 1984 II PDF](https://legge.psych.umn.edu/sites/legge.psych.umn.edu/files/files/media/legge84_binocular_contrast_summation._ii._quadratic_summation.pdf) *(abstract/summary only)*
- Meese, Georgeson & Baker (2006), "Binocular contrast vision at and above threshold," *J Vis* 6(11):1224–1243, doi:10.1167/6.11.7:
  - Average binocular summation ratio at threshold was 1.70 (4.6 dB) across four stimulus durations. It was always above √2 (3 dB) and below linear (2, 6 dB), so it is "inconsistent with the quadratic summation model of Legge and the ideal [observer]."
  - Threshold summation was 4.5 dB in Exp 2 and 4.8 dB in Exp 3.
  - An ideal observer combining noisy monocular signals predicts only about 1.2 dB with accelerating transduction (3 dB if transduction is linear). Probability summation gives at most about 1.5 dB. Neither can explain 4.5 dB, so the summation must be physiological (neural). — [Meese et al. 2006 (Aston eprint PDF)](https://publications.aston.ac.uk/id/eprint/1422/1/jov_6_11_7.pdf)
- Baker, Lygo, Meese & Georgeson (2018), "Binocular summation revisited: Beyond √2," *Psychol Bull* 144(11):1186–1199, doi:10.1037/bul0000163. Meta-analysis of 65 studies (N = 716) plus a new spatiotemporal experiment (N = 41) — [PMC6195301](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301):
  - Pooled ratio: 1.53 (3.72 dB) equal-weighted, 1.50 (3.54 dB) sample-size-weighted, 1.47 (3.35 dB) inverse-variance-weighted (lower bound 1.43). All are significantly above √2. The new experiment averaged 1.65 (4.33 dB).
  - Model: resp = L^m + R^m, so ratio = 2^(1/m). Best-fit m = 1.75 in the meta-analysis (true summation ≈ 1.49) and m = 1.26 in the experiment (≈ 1.74).
  - Summation falls with stimulus speed (TF/SF): m = 1.14 + 0.28·log10(TF/SF), R² = 0.70. A tenfold TF increase cuts summation by about 12%; a tenfold SF increase raises it by about 12% (≈1 dB).
  - Interocular sensitivity imbalance lowers measured summation by about 0.3–0.5 dB per 1 dB of imbalance.
  - The nonlinearity that sets the ratio is placed before binocular combination; "later nonlinearities do not change the summation ratio."
  - The authors conclude that summation ranges from about √2 to 2 with conditions rather than having one canonical value.

**Dipper / masking (dichoptic vs monoptic vs binocular)**
- Meese et al. (2006) — [PDF](https://publications.aston.ac.uk/id/eprint/1422/1/jov_6_11_7.pdf):
  - Dichoptic masking (mask in one eye, test in the other) is "much more severe" than monocular or binocular masking.
  - Peak facilitation (dipper depth) was 4.1 dB for dichoptic vs 9.9 dB (monocular) and 8.0 dB (binocular).
  - Monocular and binocular discrimination functions converge at pedestals around 3% and above (replicating Legge 1984b). Near the dip, however, binocular summation is ≥ 4.5 dB.
  - In a half-binocular condition (pedestal in both eyes, test in one), summation persisted across the full range of mask contrasts (mean 5.1 dB; two-stage model predicted 5.2 dB). The authors conclude that binocular summation is not abolished above threshold; it is masked by the measurement.
  - Dichoptic masking is much more severe than half-binocular masking even though the half-binocular condition contains more total mask energy.
- Meese et al. (2006) model fits (same PDF):
  - A Legge-type model (power transducer, then summation, then gain control) fit reasonably overall (RMS 1.75 dB) but had "a complete absence of dichoptic facilitation."
  - The two-stage model fit best. Stage 1 (each eye): resp_L = C_L^m / (S + C_L + C_R). Binocular sum: B = resp_L + resp_R. Stage 2: resp = B^p / (Z + B^q).
  - Fitted parameters (Exp 2, RMS 0.873 dB): m = 1.28, S = 0.985, p = 7.99, q = 6.59, Z = 0.077. Exp 3 refit: m = 1.24, S = 1.42, p = 7.67, q = 6.10, Z = 0.062.
  - The authors note that the low m gives high binocular summation, "but dichoptic masking remains severe due to the inclusion of interocular suppression." With interocular suppression removed, the model reverts to the failing Legge-type behaviour.
  - Equation and parameter values were read from a text extraction of the PDF's tables; check them against the published Table 2/3 before reuse.
- Meese et al. (2006) also cite cat physiology (Truchard, Ohzawa & Freeman 2000), which concluded that "suppressive gain control occurs initially at a purely monocular stage before linear binocular summation," followed by weaker suppression after summation. — [PDF](https://publications.aston.ac.uk/id/eprint/1422/1/jov_6_11_7.pdf)

**Contrast matching / ocularity invariance**
- Legge & Rubin (1981), "Binocular interactions in suprathreshold contrast perception," *Percept Psychophys* 30(1):49–61: with 1 and 8 c/deg gratings, binocular contrast-matching functions departed from simple averaging, with individual differences in a form of ocular dominance. Later work cites it for the finding that perceived suprathreshold contrast is almost the same with one eye or two. — [Legge & Rubin 1981 PDF](https://legge.dl8.umn.edu/sites/legge.psych.umn.edu/files/2020-08/Legge%20Rubin%201981%20Binocular%20interactions%20in%20suprathreshold%20contrast%20perception%20.pdf) *(abstract/summary only)*
- Anderson & Movshon (1989), "Binocular combination of contrast signals," *Vision Res* 29:1115–1132. Secondary sources report that added noise and contrast adaptation changed the shape of binocular contrast contours, which a single-channel model could not explain. — [secondary summary via search; e.g., PMC4521337](https://pmc.ncbi.nlm.nih.gov/articles/PMC4521337) *(not read)*
- Baker, Meese & Georgeson (2007), "Binocular interaction: Contrast matching and contrast discrimination are predicted by the same model," *Spatial Vision* 20(5):397–413, doi:10.1163/156856807781503622 — [Aston record](https://publications.aston.ac.uk/id/eprint/1426) *(abstract only)*:
  - Design: 25 stimuli crossing five left-eye contrasts (0.3–32%) with five right-eye contrasts (0.3–32%) of 1 c/deg horizontal gratings.
  - Tasks: increment discrimination (increment in one eye) and contrast matching (monocular or binocular test vs standard).
  - "With all model parameters fixed by previous work," the two-stage model predicted both data sets and beat linear and quadratic summation models.
- Ding & Sperling (2006), "A gain-control theory of binocular combination," *PNAS* 103(4):1141–1146, doi:10.1073/pnas.0509629103 — [PMC1347993](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993):
  - Each eye attenuates the other eye's signal in proportion to its own contrast energy, and each eye also gain-controls the other eye's gain control.
  - Very low-contrast inputs add linearly; "the intrinsic nonlinearity manifests itself only as contrast increases."
  - Tested with a phase-combination paradigm: 0.68 c/deg gratings with interocular contrast and phase differences. A single estimated parameter explained about 95% of variance across 48 conditions *(the 95% figure is from a secondary copy of the PNAS page)*.
  - Later DSKL extensions add interocular gain *enhancement* and a fusion stage. — [Huang et al. 2010, PMC3000330](https://pmc.ncbi.nlm.nih.gov/articles/PMC3000330)

### Inferences
- **Benchmarks for V1.** The human benchmarks a V1 population read-out has to reproduce are:
  - (i) threshold-level binocular advantage of about 1.5–1.7, implying an effective pre-summation exponent m ≈ 1.2–1.75;
  - (ii) near-zero binocular advantage for suprathreshold discrimination and appearance;
  - (iii) dichoptic masking much stronger than monoptic masking, with weak dichoptic facilitation (about 4 dB vs 8–10 dB);
  - (iv) an interocular suppression weight around 1, i.e., the other eye's contrast enters each eye's denominator as strongly as its own.
- **Prediction for monocular neurons.** The two-stage model predicts that responses to the dominant eye should be divisively suppressed by contrast in the *non-driving* eye with roughly equal weight, and at a stage before binocular summation. This contrasts with the cat physiology quoted by Meese et al. (Truchard et al. 2000), where gain control is mainly monocular. For layer-resolved data the test is: do monocular (granular, L4C-like) units or CSD sinks show divisive suppression from the other eye? The 6×6 matrix tests this directly. For a strictly monocular unit, rows/columns where the non-driving eye goes 0 → .9 are a dichoptic-masking function.
- **Expected summation ratio for the team's stimuli.** The speed dependence (m = 1.14 + 0.28·log10(TF/SF)) means slow or static, flashed gratings should sit at the high-summation end, with m close to 1 and a ratio near 1.7–1.8. So near-threshold contrasts (.055, .11) are where binocular facilitation of V1 population responses should be largest. This is an inference; it assumes neural CRF exponents follow the psychophysical transducer.
- **Design analogues.**
  - The 6×6 matrix is close to the Baker, Meese & Georgeson (2007) 5×5 L×R design (0.3–32%). That allows a direct neural analogue of their discrimination and matching surfaces: iso-response contours in (C_L, C_R) space, fitted with a Minkowski exponent to give linear (1), quadratic (2) or max-like contours.
  - The "half-binocular" condition maps onto cells where both eyes have a pedestal and contrast changes in one eye.
  - Dichoptic masking maps onto rows/columns for units dominated by one eye.
- **Normalize by each eye's sensitivity first.** Baker et al. (2018) show that imbalance lowers apparent summation by 0.3–0.5 dB per dB of imbalance. Neural summation ratios should therefore be computed after equating each unit's monocular sensitivity (e.g., per-eye c50), or relative to the better eye, to compare with psychophysics.

### Gaps
- I could not open the full text of Baker & Meese (2007), "Binocular contrast interactions: dichoptic masking is not a single process" (*Vision Res*; [White Rose PDF](https://eprints.whiterose.ac.uk/109662/1/BakerMeese2007R.pdf); text extraction failed), or Baker, Meese & Summers (2007, *Neuroscience*, "two routes to suppression before binocular summation"). Their tuning results (orientation, SF and TF selectivity of dichoptic vs monocular suppression) are not verified here.
- I did not verify exact contrast-matching numbers from Legge & Rubin (1981), Anderson & Movshon (1989) or Baker et al. (2007): how much higher monocular contrast must be to match binocular contrast, and how that varies with contrast.
- I did not obtain a reliable source for "Fechner's paradox" magnitudes.

---

## 2. SSVEP/EEG (frequency tagging, intermodulation), MEG and fMRI (including laminar 7T): stage and strength of interocular suppression

### Takeaway
Eye-tagged SSVEP and MEG show that a single interocular normalization rule can describe cortical population responses: excitation is summed across eyes, and the divisive pool contains both eyes. The fits give exponents around 2.2–2.4 and ocularity invariance at suprathreshold contrasts. Dichoptic masking acts mainly as contrast gain at the first harmonic (5 Hz), though one data set looked more like response gain. 3T fMRI shows almost no BOLD increment in V1 from adding a matched grating to the second eye: the additivity index at 10%/10% was 0.01. Orthogonal (rivalrous) dichoptic gratings produce only about half the interocular suppression that matched gratings do. Laminar 7T work places eye-specific rivalry effects in superficial V1, but GE-BOLD has a superficial vascular bias, and only CBV/VASO peaks at the expected middle depth. The laminar locus of interocular normalization is therefore unresolved in humans.

### Cited Findings
**SSVEP/EEG**
- Baker & Wade (2017), "Evidence for an optimal algorithm underlying signal combination in human visual cortex," *Cereb Cortex* 27(1):254–264, doi:10.1093/cercor/bhw395 — [PMC5903417](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417):
  - Design ("eye" experiment): 1 c/deg horizontal grating patches (2° diameter) in a 5×5 grid (14° field), shown through shutter goggles. Eyes tagged at 5 and 7.5 Hz; contrasts 4–64%; mask fixed at 32%; N = 12; 64-channel EEG analyzed at Oz/POz.
  - Of six models, "late summation" fit best: resp = (A^p + B^p) / (Z^q + A^q + B^q).
    - Fit for the eye experiment: p = 2.22, q = 2.22, Z = 9.48 (% contrast), R_max = 0.71; SSE 0.01; 92% of variance explained.
    - Fit for the space experiment: p = 2.43, q = 2.18, Z = 7.46.
    - Early summation, linear summation and independent transducers fit worse (eye SSE 0.04, 1.11 and 0.24).
  - Summing two in-phase components produced "very little increase" in the response (ocularity invariance), and only late summation predicted this.
  - With no free parameters except one R_max change, the suppression (mask) conditions were predicted with 73% of variance explained. The eye fit was weaker because the curve was shallower than predicted, which the authors read as a response-gain effect.
  - Across 1–30 Hz (harmonics and intermodulation, e.g., SNR > 2 at 3, 4, 8, 9, 16, 18 and 19 Hz) the model explained 69% of variance.
  - "Optimal" refers to the form of the model, which resembles a Kalman filter weighting each input by inverse variance. This is speculative.
  - The authors present it as an algorithm-level account and expect "single-cell neurophysiology" to reveal the implementation. They make no single-neuron predictions and do not localize suppression to V1 specifically.
- Baker, Vilidaite & Wade (2021), "Steady-state measures of visual suppression," *PLoS Comput Biol* 17(10):e1009507 — [PMC8544832](https://pmc.ncbi.nlm.nih.gov/articles/PMC8544832) *(abstract/summary only)*:
  - A computational meta-analysis of 16 SSVEP studies found the most evidence for contrast gain with overlaid masks.
  - New experiments tested monocular and dichoptic overlay masks plus aligned and orthogonal surrounds on 5 Hz gratings. At the occipital pole and 1F (5 Hz), there was strong evidence for contrast gain for all four mask types.
  - At 2F (10 Hz), suppression was stronger and involved both contrast and response gain.
  - Suppression got stronger at more lateral electrodes.
  - Interpretation: gain control serves a common purpose across anatomically distinct routes, possibly "to suppress less reliable inputs."
- Hou, Nicholas & Verghese (2020), "Contrast normalization accounts for binocular interactions in human striate and extra-striate visual cortex," *J Neurosci* 40(13):2753–2763, doi:10.1523/JNEUROSCI.2043-19.2020. Source-imaged SSVEP of dichoptic, eye-tagged stimuli. Both forms of binocular interaction (eye-specific self-term suppression and intermodulation) are explained by a shared gain-control mechanism in striate and extrastriate cortex. — [SKI record](https://www.ski.org/publication/contrast-normalization-accounts-for-binocular-interactions-in-human-striate-and-extra-striate-visual-cortex/) *(abstract only; full text not read)*
- Hou, Tyson, Uner, Nicholas & Verghese (2021), *J Neurosci* 41(41):8632–8643, doi:10.1523/JNEUROSCI.0268-21.2021 — [PMC8513700](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513700):
  - Design: dichoptic parallel 2 c/deg gratings. Non-dominant eye tagged at 6.07 Hz, swept 1.7–40% in 10 log steps; dominant eye at 8.5 Hz, fixed at 20%. Self-terms 2F1/2F2; dominant IM term F1+F2 (14.57 Hz). ROIs: V1, V3a, hV4, hMT+, LOC.
  - The model sums the eyes linearly before divisive normalization. The authors note this is equivalent to pre-combination interocular suppression when monocular exponents = 1.
  - Controls (from the 2020 data): excitatory exponent p ≈ 2.4, suppressive q ≈ 2.4, target weight ≈ 1.45 relative to a mask weight of 1. V1 response gain was about twice that of extrastriate ROIs.
  - IM terms "result only from interaction of the two eyes' inputs."
- Carter, Baker, Morland, Lawton & Wade (2025), "Interocular suppression in chromoluminance patterns measured with SSVEP," *J Vis* 25(4):6, doi:10.1167/jov.25.4.6 — [PMC12007675](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007675):
  - Design: eyes tagged at 5 and 7 Hz; luminance, L–M and S-cone disks and 0.5 c/deg gratings; N = 12.
  - Raw suppression depended on chromaticity (p = 0.002) and SF (p = 0.005). Suppression as a *proportion* of the monocular response did not (p = 0.82 and 0.56).
  - First-order IM terms (2 and 12 Hz) were significant and second-order IM terms were not. The authors conclude binocular combination occurs in neurons carrying half-wave-rectified (simple-cell-like) signals that respond at the fundamental.
  - The data do not support strict segregation of weakly binocular low-SF chromatic "blob" pathways.
- Segala, Bruno, Martin, Aung, Wade & Baker (2023), "Different rules for binocular combination of luminance flicker in cortical and subcortical pathways," *eLife* 12:RP87048 — [eLife](https://elifesciences.org/articles/87048); [PMC10522334](https://pmc.ncbi.nlm.nih.gov/articles/PMC10522334) *(abstract/reviews only)*:
  - EEG and binocular pupillometry were recorded simultaneously with luminance flicker.
  - Binocular facilitation was substantial in both pathways and strongest in cortex. Interocular suppression was stronger in the subcortical (pupil) pathway than in cortex.
  - Perceptual luminance matching showed near-linear binocular additivity rather than ocularity invariance.
  - A reviewer noted these flicker results differ markedly from grating results, and that the two pathways were measured with different methods.
- Katyal, Vergeer, He, He & Engel (2018), "Conflict-sensitive neurons gate interocular suppression in human visual cortex," *Sci Rep* — [PMC5775389](https://pmc.ncbi.nlm.nih.gov/articles/PMC5775389) *(abstract/summary only)*:
  - SSVEP intermodulation amplitude grew with interocular conflict (orientation difference).
  - Adapting to conflict both reduced this IM effect and raised the conflict needed for perceptual suppression. The authors propose that suppression is gated by conflict-detecting neurons rather than only by eye-vs-eye competition.
  - A later paper notes direct neural evidence is still lacking.

**MEG**
- Chadnova, Reynaud, Clavagnier, Baker, Baillet & Hess (2018), "Interocular interaction of contrast and luminance signals in human primary visual cortex," *NeuroImage* (pp. 23–30), doi:10.1016/j.neuroimage.2017.10.035 — [White Rose record](https://eprints.whiterose.ac.uk/123813) *(abstract only)*:
  - Design: binary-noise stimuli tagged at 4 and 6 Hz, 0–32% contrast; MEG with individual fMRI retinotopy.
  - V1 SSVER magnitude rose monotonically with contrast. Dichoptic masking reduced power for both eyes, well explained by a normalization model.
  - A 1.5 ND filter (30× less light, same contrast) delayed that eye's V1 processing by about 38 ms and weakened its masking of the other eye. The model was extended with a temporal-filter stage.
  - An early 4-subject version concluded that suppression occurs "as early as" V1. — [BMC Neurosci abstract](https://bmcneurosci.biomedcentral.com/articles/10.1186/1471-2202-15-S1-P145)

**fMRI (3T)**
- Moradi & Heeger (2009), "Inter-ocular contrast normalization in human visual cortex," *J Vis* 9(3):13, 1–22 — [Heeger lab PDF](https://www.cns.nyu.edu/heegerlab/content/publications/Moradi-JoV2009.pdf):
  - Design: 3T, V1/V2/V3, 2-back task at fixation to divert attention; 1.5 s gratings. Conditions: monocular gratings, binocular gratings (same orientation/SF/phase, independent contrasts per eye), dichoptic plaids, monocular plaids.
  - Additivity index: AI = [r(c_ref, c) − r(0, c)] / (0.5·[r(c_ref, 0) + r(0, c_ref)]). Linear summation gives AI = 1.
  - Measured AI for 10% in both eyes was 0.01 (95% CI −0.23 to 0.25). A power-law hemodynamic nonlinearity alone would predict 0.46–0.57 (k = 0.47–0.65), so V1 BOLD shows near-complete ocularity invariance that the hemodynamics cannot explain.
  - Binocular normalization (eyes normalize each other) fit better than linear summation.
  - At 10% contrast: AI = 0.15 for binocular gratings, 0.31 for dichoptic plaids, 0.47 for monocular plaids. Fits required orthogonal orientations to contribute "only about half as much suppression as compatible orientations."
  - Responses to dichoptic plaids were significantly larger than to monocular gratings. There was no V1 correlate of flash suppression with attention diverted. V2 and V3 behaved like V1.
  - The authors flag a conflict with single-cell data suggesting that normalization and cross-orientation suppression are "mainly monocular" (Truchard et al. 2000; Walker et al. 1998). Walker et al. found interocular suppression for orthogonal stimuli in almost all cells, but it was generally weak.

**Laminar / high-field fMRI**
- Zaretskaya et al. (2020), "Eye-selective fMRI activity in human primary visual cortex: Comparison between 3 T and 9.4 T, and effects across cortical depth," *NeuroImage* 220:117078 — [DOAJ record](https://doaj.org/article/0d45249dd3de444aa9970c8f6cae2073); [Graz summary](https://neurovision.uni-graz.at/en/news/untersuchung-von-augendominanz-signalen-im-visuellen-kortex-mittels-hochaufloesendem-fmri) *(abstract only)*:
  - 9.4T gave more eye-selective voxels and about double the run-to-run consistency of eye preference.
  - GE-EPI eye selectivity was strongest in *superficial* depths, against the expected middle-layer peak. The authors attribute this mainly to draining veins and ascending venules.
- Decoding eye-of-origin across depth with BOLD vs CBV at 7T (bioRxiv 2023.09.28.560016; peer-reviewed 2025): only CBV-fMRI gave the highest decoding accuracy at the depth expected from neurophysiology (middle, L4C). Macrovascular contributions affect all methods' laminar profiles. — [PMC12407697](https://pmc.ncbi.nlm.nih.gov/articles/PMC12407697); [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.28.560016.full.pdf) *(abstract only; authorship likely Haenelt et al., not verified here)*
- Qian, Zhang, Chen, de Hollander, Knapen, He & Zhang (2025), "Mesoscale cortical mechanisms of perceptual conflict resolution in binocular rivalry," *Nat Hum Behav* (published 13 Nov 2025; listed as 10(2):403–416), doi:10.1038/s41562-025-02320-4 — [CAS news release](https://english.cas.cn/newsroom/research_news/life/202511/t20251124_1133129.shtml); [VU record](https://research.vu.nl/en/publications/mesoscale-cortical-mechanisms-of-perceptual-conflict-resolution-i/) *(summary only)*:
  - 7T laminar/columnar fMRI. Rivalry arises from lateral (interocular) inhibition in the superficial layers of V1 ocular dominance columns, not in LGN eye-specific layers.
  - Eye-specific feedback from IPS modulated and synchronized competition in superficial and deep V1 layers.
  - The authors claim to reconcile awake-monkey electrophysiology with earlier 3T fMRI.
  - In the 2021 ISMRM abstract, eye-specific modulation peaked in the middle layer during physical replay but was stronger in the superficial layer during rivalry, and the LGN showed minimal rivalry modulation. — [ISMRM 2021 abstract](https://archive.ismrm.org/2021/0636.html)
- A 7T VASO study of ocular dominance across depth (Akbari et al., ISMRM 2023) reported that the unstimulated eye's response is smaller than the stimulated eye's under monocular stimulation. — [ISMRM 2023](https://archive.ismrm.org/2023/0622.html) *(abstract only)*

### Inferences
- **Common ground.** Across EEG, MEG and fMRI, human population signals in V1 are consistent with interocular normalization: excitation is summed and the denominator pools both eyes. Matched suprathreshold gratings in the two eyes then give little or no extra response (AI ≈ 0.01 in BOLD; "very little increase" in SSVEP).
  - For V1 spiking or MUA, the matching benchmark is the diagonal of the 6×6 (c, c) compared with the monocular edges (c, 0) and (0, c).
  - Prediction: binocular/monocular ratio ≈ 1.5–2 at .055/.11, falling to ≈ 1.0–1.2 at .45/.9. A purely linear binocular neuron would give 2.
  - The Baker & Wade late-summation form with p = q ≈ 2.2 makes exactly this prediction: binocular ≈ 2× monocular at low contrast, ≈ 1× at saturation.
- **Additivity index without tagging.** Moradi & Heeger's AI can be computed directly from the 6×6 for every c_ref and every laminar compartment. The 2×2 interaction contrast r(cL, cR) − r(cL, 0) − r(0, cR) + r(0, 0) is the static analogue of SSVEP intermodulation: it isolates the nonlinear interocular interaction without frequency tagging.
- **Fixed-other-eye analogues.** The fixed-mask SSVEP designs map onto single rows/columns of the 6×6 with the other eye at .225 or .45: Hou et al. (target 1.7–40%, other eye fixed at 20%), Baker & Wade (4–64%, mask 32%) and Chadnova et al. (0–32%). Fit the same models (late summation; linear-sum-then-normalize with target weight about 1.45; two-stage) to unit, MUA and CSD data and compare exponents (≈ 2.2–2.4 in EEG) and suppression weights.
- **Contrast gain or response gain?** Baker, Vilidaite & Wade (2021) find contrast gain for dichoptic masks at 1F, while Baker & Wade (2017) saw a shallower, response-gain-like curve. For monocular-dominant units, the 6×6 can decide whether non-driving-eye contrast shifts the CRF horizontally (contrast gain) or scales it (response gain). The answer may differ by layer.
- **Laminar locus.** Human laminar fMRI cannot yet say where interocular normalization happens:
  - GE-BOLD is biased to the surface (Zaretskaya et al.);
  - CBV peaks at mid-depth for eye-of-origin;
  - Qian et al.'s superficial-layer rivalry locus comes with that caveat.
  - Macaque CSD and depth-resolved MUA in L4C vs L2/3 vs L5/6 can test the competing predictions directly. The two-stage model predicts interocular suppression already in monocular (granular) signals. A purely "superficial lateral inhibition between ODCs" account predicts interocular suppression emerging mainly in supragranular layers, with a time lag.
- **Phase-sensitive combination.** Carter et al.'s 1F intermodulation argument (combination in half-wave-rectified neurons) predicts that interocular interactions should be strong in simple-like, phase-sensitive V1 cells. With static, same-phase gratings the team cannot test phase dependence, but can ask whether simple-like units (by F1/F0 if drifting-grating data exist, or by response transients) show stronger interocular normalization than complex-like units.
- **Same-orientation stimuli measure the strongest case.** Moradi & Heeger found same-orientation interocular suppression to be about twice as strong as orthogonal. The team's same-orientation gratings therefore probe the condition with the strongest predicted interocular normalization. Comparing it with any orthogonal-dichoptic data would test their surprising result: the opposite of what rivalry intuition and Walker et al.'s cat data suggest.

### Gaps
- Exact numeric results of Hou et al. (2020) (V1 vs extrastriate weights) and Chadnova et al. (2018) (normalization weights, volume number) were not read in full.
- Segala et al. (2023) effect sizes (binocular/monocular ratios for EEG vs pupil) were not retrieved.
- I found no human laminar fMRI study that measures binocular *contrast* combination (summation or dichoptic masking across a contrast range) as a function of cortical depth. All the laminar work found concerns eye-of-origin selectivity or rivalry. This appears to be an open niche that macaque laminar data can fill.
- The Norcia lab's foundational dichoptic IM work (Brown, Candy & Norcia 1999, infants; Baitch & Levi 1988; Norcia et al. 2000) and Tsai, Wade & Norcia (2012) were not re-verified this session.

---

## 3. Population-level phenomena: low-contrast facilitation vs high-contrast suppression, contrast-dependent ocular dominance, eye-specific adaptation and plasticity, rivalry onset, amblyopia

### Takeaway
Human data consistently show the same pattern: near-linear binocular addition at low contrast, and normalization at high contrast that removes the binocular advantage and makes the higher-contrast eye dominate. Eye balance is plastic: 2 h of monocular patching boosts the deprived eye's V1 BOLD response, and the effect predicts the shift in rivalry dominance. Amblyopia studies indicate that the binocular summation mechanism survives. What changes is monocular attenuation and noise in the weak eye, and excitatory (intermodulation-producing) interactions are reduced. Interocular suppression stays intact or even weakens, rather than becoming hyper-strong.

### Cited Findings
- **Low vs high contrast.**
  - Ding & Sperling (2006): low-contrast left/right inputs "add linearly," and the gain-control nonlinearity appears "only as contrast increases." — [PMC1347993](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993)
  - Legge (1984a): binocular summation decreases as pedestal contrast rises. — [Legge 1984 I](https://legge.psych.umn.edu/sites/legge.psych.umn.edu/files/files/media/legge84_binocular_contrast_summation._i._detection_and_discrimination.pdf)
  - Meese et al. (2006): monocular and binocular dippers converge above about 3% pedestal contrast. — [PDF](https://publications.aston.ac.uk/id/eprint/1422/1/jov_6_11_7.pdf)
- **Interocular imbalance as a contrast-ratio effect.** Baker et al. (2018) give summation relative to the better eye as about 3.2 dB (meta subset) or 4.9 dB (experiment) at equal sensitivity, falling by 0.3–0.5 dB per dB of imbalance. — [PMC6195301](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301)
- **Luminance-induced delay changes interocular weight.** Reducing one eye's luminance delayed its V1 MEG response by about 38 ms and weakened its masking of the other eye (Chadnova et al. 2018). — [White Rose](https://eprints.whiterose.ac.uk/123813)
- **Short-term monocular deprivation.**
  - Binda, Kurzawski, Lunghi, Biagi, Tosetti & Morrone (2018), "Response to short-term deprivation of the human adult visual cortex measured with 7T BOLD," *eLife* 7:e40014: 2 h of patching boosted the deprived eye's BOLD response and shifted V1 ocular dominance toward it, consistent with homeostatic plasticity.
  - The V1 effect correlated across participants with the increase in deprived-eye dominance in binocular rivalry. It was strongest in V1, present in V2–V4, and implicated the parvocellular pathway.
  - — [eLife 40014](https://elifesciences.org/articles/40014v1)
  - A follow-up found a deprived-eye boost in ventral pulvinar but not LGN (Kurzawski et al., eLife). — [eLife 74565](https://elifesciences.org/articles/74565)
- **Rivalry onset / conflict.** Moradi & Heeger (2009) found no V1 correlate of flash suppression under diverted attention. Interocular suppression between incompatible (orthogonal) gratings was about half that between compatible ones and did not depend on eye of origin. — [PDF](https://www.cns.nyu.edu/heegerlab/content/publications/Moradi-JoV2009.pdf). See also Katyal et al. 2018 and Qian et al. 2025 in Section 2.
- **Amblyopia.**
  - Baker, Meese, Mansouri & Hess (2007), "Binocular summation of contrast remains intact in strabismic amblyopia," *IOVS* 48(11):5332–5338, doi:10.1167/iovs.07-0194. N = 9 strabismic amblyopes; 3 or 9 c/deg, 200 ms gratings. The conventional summation ratio was near 1. After scaling contrast to equate monocular sensitivities it rose to about √2 or higher, similar to controls, and the deficit could be simulated with an ND filter. — [Southampton record](https://eprints.soton.ac.uk/49281)
  - Baker, Meese & Hess (2008), *Vision Res*, "Contrast masking in strabismic amblyopia: attenuation, noise, interocular suppression and binocular summation": the best model had signal attenuation and added noise in the amblyopic eye, with *intact* interocular suppression and binocular summation stages. — [Aston record](https://research.aston.ac.uk/en/publications/contrast-masking-in-strabismic-amblyopia-attenuation-noise-intero)
  - Hou et al. (2021): in strabismic amblyopia, mean excitatory exponent p fell to 0.8 (vs 2.4 in controls; anisometropes 1.9) while the suppressive exponent q was similar (2.7 vs 2.4). The amblyopic-eye target weight dropped to about 0.33 in strabismics with fusion (controls ≈ 1.45). IM responses were significantly reduced in all ROIs (absent without fusion). σ rose and R_m fell, with hMT+ relatively spared. The authors frame this as electrophysiological support for preferential loss of *excitatory* binocular connections (Sengpiel & Blakemore 1996). — [PMC8513700](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513700)
  - Chadnova et al. (2017, *Vision Res*, MEG): amblyopes showed about a 20 ms interocular processing delay and *weaker* suppressive interactions than normal observers. — [search record, PMID 28842260](https://neuinfo.org/literature/pmid:28842260?rpKey=on) *(abstract/summary only)*

### Inferences
- **Contrast-dependent ocular dominance.** Ding & Sperling and the two-stage model predict that the effective weight of each eye depends on the interocular contrast ratio *and* on absolute contrast: linear at low contrast, increasingly winner-take-more at high contrast. In the 6×6, compute an ocular dominance index (or fitted eye weight) separately at each total-contrast level. The human prediction is that population OD should not be a fixed property of a recording site: dominance by the higher-contrast eye should sharpen as contrast rises.
- **Amblyopia logic carried over to normal V1.** The amblyopia results separate "excitatory binocular summation" from "interocular suppression" as dissociable components. In normal V1 the matching question is whether units (or layers) differ in the balance of the two. For example, are some units binocularly facilitated at low contrast with little suppression, while others are strongly suppressed by the non-driving eye? The Maier-lab finding that monocular V1 neurons are binocularly modulated is the relevant animal analogue (not re-verified here).
- **Eye-specific adaptation.** Eye-specific plasticity (deprived-eye boost) and the luminance/latency effects imply that interocular weights adapt on the scale of minutes to hours and depend on timing. The team cannot manipulate this, but trial-to-trial or session-to-session drift in eye weights, and latency differences between eyes (e.g., the 3.4 ms NSx/NEV clock offset must be corrected first), are confounds and observables worth checking.

### Gaps
- No human study found here directly maps "binocular facilitation at low contrast vs suppression at high contrast" onto a single population CRF surface across the full L×R matrix in V1 BOLD or EEG. Moradi & Heeger and Baker & Wade sample only parts of it.
- The Lunghi, Burr & Morrone (2011, *Curr Biol*) psychophysical monocular-deprivation result and its EEG follow-ups (Lunghi et al. 2015, *J Physiol*) were not re-verified this session.
- Models of rivalry onset as interocular normalization (Said & Heeger 2013, *PLoS Comput Biol*; Li, Rankin, Rinzel, Carrasco & Heeger 2017, *PNAS*) and fMRI rivalry correlates (Polonsky et al. 2000; Tong & Engel 2001; Haynes, Deichmann & Rees 2005; Lee, Blake & Heeger 2005) are relevant background but were not verified in this session.

---

## 4. Recent controversies (2015–2026): monocular vs binocular gain-control stages, psychophysics vs SSVEP disagreements, V1 vs extrastriate locus, and explicit calls for neurophysiology

### Takeaway
The live disagreements are concrete and testable in macaque V1:
- (a) whether interocular suppression sits at a monocular stage before summation with a weight of about 1 (psychophysics), or whether normalization is mainly monocular (cat single units);
- (b) whether dichoptic suppression is contrast gain or response gain;
- (c) whether ocularity invariance holds for population signals generally, or only for spatial patterns (it fails for luminance flicker);
- (d) whether interocular inhibition is strongest for matched or mismatched orientations;
- (e) which cortical layer carries interocular suppression.

Human authors repeatedly defer the implementation question to single-cell and laminar neurophysiology.

### Cited Findings
- **(a) Stage of gain control.**
  - Psychophysics needs strong pre-summation interocular suppression (two-stage model, other-eye weight in the stage-1 denominator) plus post-summation gain control. Without the interocular term, dichoptic masking cannot be fit. — [Meese et al. 2006](https://publications.aston.ac.uk/id/eprint/1422/1/jov_6_11_7.pdf)
  - Cat V1 single-unit work (Truchard, Ohzawa & Freeman 2000) concluded gain control is "initially at a purely monocular stage before linear binocular summation," with weaker post-summation suppression, as cited by Meese et al. (2006). Moradi & Heeger (2009) state that their fMRI results "might seem at odds with single cell recordings that suggest contrast normalization and cross-orientation suppression are mainly monocular." — [Moradi & Heeger PDF](https://www.cns.nyu.edu/heegerlab/content/publications/Moradi-JoV2009.pdf)
  - EEG models differ: Hou et al. sum the eyes linearly before normalization and note equivalence to pre-combination suppression only when monocular exponents = 1 — [PMC8513700](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513700). Baker & Wade's late-summation model normalizes each eye by both eyes before summation — [PMC5903417](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417).
- **(b) Contrast vs response gain.**
  - SSVEP meta-analysis and new data: contrast gain dominates for dichoptic and monocular overlay masks at 1F, with both types at 2F. — [Baker, Vilidaite & Wade 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8544832)
  - The Baker & Wade (2017) eye data showed a shallower-than-predicted, response-gain-like curve. — [PMC5903417](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417)
- **(c) Stimulus-dependent rules.**
  - For luminance flicker, perception shows near-linear additivity and cortical interocular suppression is weaker than subcortical. Reviewers note this differs markedly from grating results. — [Segala et al. 2023](https://elifesciences.org/articles/87048)
  - Summation ratios vary with speed (TF/SF) and are not a fixed √2. — [Baker et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301)
  - Proportional interocular suppression is similar across luminance, L–M and S-cone gratings. — [Carter et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007675)
- **(d) Orientation tuning of interocular suppression.**
  - Same-orientation interocular suppression is about twice orthogonal in human V1 BOLD (AI 0.15 vs 0.31). — [Moradi & Heeger 2009](https://www.cns.nyu.edu/heegerlab/content/publications/Moradi-JoV2009.pdf)
  - SSVEP IM grows with interocular conflict, interpreted as conflict-sensitive neurons gating suppression. — [Katyal et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5775389)
- **(e) Locus.**
  - MEG localizes dichoptic masking and normalization to V1. — [Chadnova et al. 2018](https://eprints.whiterose.ac.uk/123813)
  - Source-imaged SSVEP finds a shared gain-control mechanism in V1 and extrastriate cortex, with V1 response gain about 2× extrastriate in controls. — [Hou et al. 2020 record](https://www.ski.org/publication/contrast-normalization-accounts-for-binocular-interactions-in-human-striate-and-extra-striate-visual-cortex/); [Hou et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513700)
  - fMRI finds V1, V2 and V3 similar. — [Moradi & Heeger 2009](https://www.cns.nyu.edu/heegerlab/content/publications/Moradi-JoV2009.pdf)
  - 7T places rivalry-related interocular inhibition in superficial V1 ODCs, not LGN — [Qian et al. 2025](https://english.cas.cn/newsroom/research_news/life/202511/t20251124_1133129.shtml) — but GE-BOLD depth profiles are vascularly biased toward the surface. — [Zaretskaya et al. 2020](https://doaj.org/article/0d45249dd3de444aa9970c8f6cae2073); [CBV vs BOLD decoding](https://pmc.ncbi.nlm.nih.gov/articles/PMC12407697)
- **Explicit deferrals to neurophysiology.**
  - Baker & Wade (2017) call their account algorithm-level (Marr) and expect "single-cell neurophysiology" to reveal how it is implemented. — [PMC5903417](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417)
  - Qian et al. (2025) frame their 7T result as reconciling awake-monkey electrophysiology with 3T fMRI. — [CAS](https://english.cas.cn/newsroom/research_news/life/202511/t20251124_1133129.shtml)
  - Later work on Katyal's conflict-neuron proposal says more direct neural evidence is lacking. — [PMC5775389 search summary](https://pmc.ncbi.nlm.nih.gov/articles/PMC5775389)
  - Hou et al. (2021) cite macaque work (Smith et al. 1997; Hallum et al. 2017) as the physiological basis for their excitatory-loss interpretation. — [PMC8513700](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513700)

### Inferences
Prioritized human-derived predictions testable with the team's 6×6 laminar data:

1. **Pre-summation interocular suppression (two-stage model).**
   - Prediction: monocular-dominant units, especially in granular L4C, are divisively suppressed by the non-driving eye with a weight comparable to their own eye (w ≈ 1).
   - Test: fit stage-1 R = C_d^m / (S + C_d + w·C_nd) on rows/columns with the non-driving eye at 0 → .9.
   - A finding of w ≈ 0 in L4C with w > 0 only in L2/3 would favour the cat "monocular-first" view and the superficial-layer locus (Qian et al.). It would challenge the psychophysical architecture, or relocate its "stage 1" to supragranular monocular cells.
2. **Ocularity invariance at high contrast, facilitation at low.**
   - Prediction: the diagonal-to-edge ratio is about 1.5–2 at .055/.11 and about 1 at .45/.9. The Moradi AI is about 0 at moderate contrast.
   - Population read-outs (MUA, CSD, LFP gamma) may differ from individual units. Human BOLD and SSVEP likely weight superficial pyramidal activity (Baker & Wade note SSVEP mostly reflects superficial pyramidal cells).
3. **Iso-response contour shape.**
   - Prediction: Minkowski exponent m ≈ 1.2–1.75 near threshold; contours that are more max-like at high contrast (Ding & Sperling).
   - Testable per unit and per layer.
4. **Contrast gain vs response gain of dichoptic suppression.**
   - Prediction: mostly contrast gain (Baker, Vilidaite & Wade 2021).
   - Laminar differences would explain why EEG harmonics disagree (1F contrast gain; 2F mixed).
5. **Timing.**
   - Prediction: interocular suppression arrives later than monocular excitation, and a delayed eye suppresses less (Chadnova et al. 2018; Richard, Chadnova & Baker 2018, "Binocular vision adaptively suppresses delayed monocular signals," *NeuroImage* — [White Rose PDF](https://eprints.whiterose.ac.uk/127595/1/Richardetal2018repositoryversion.pdf), not read).
   - Flashed ~317 ms stimuli give onset transients for latency analysis by layer. This partly substitutes for the phase information SSVEP provides.
6. **Not testable with this design:** frequency-tagged intermodulation spectra (use the static 2×2 interaction contrast instead), orientation dependence (same-orientation only; Moradi's matched > orthogonal claim and Katyal's conflict neurons need orthogonal dichoptic conditions), phase-combination (Ding–Sperling) and interocular-delay manipulations.

### Gaps
- I found no 2024–2026 human paper that explicitly proposes a single-neuron or laminar test of the two-stage model. The calls for neurophysiology are general rather than specific.
- An abstract titled "Interocular normalization in monkey primary visual cortex (V1)" exists ([AES abstract AB015](https://aes.amegroups.com/article/view/5220/html)) but was not read. It may be directly relevant prior animal work for the coordinator's other researchers.
- Whether human psychophysical "stage 1" maps onto monocular L4C cells or onto monocular-biased supragranular cells is not addressed in any human source found. This is the central open question the team's data could answer.
