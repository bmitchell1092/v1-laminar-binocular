# Single-neuron and laminar neurophysiology of binocular vs monocular responses in primate (and cat) V1, with emphasis on contrast

*Tags used throughout: [cat-A] = anesthetized, paralyzed cat; [mac-A] = anesthetized macaque; [mac-W] = awake, fixating macaque; [human], [mouse], [model]. Sources were retrieved via Europe PMC/PMC abstracts and open full texts (Oct 2026). Links are DOIs. Numbers come from the papers' abstracts, results text, tables or figure legends unless marked otherwise.*

## 1. Monocular vs binocular contrast response functions (CRFs): summation, facilitation, binocular and interocular suppression, and changes in c50 / Rmax / exponent

### Takeaway
Near threshold, the two eyes' contrast signals combine roughly linearly in single V1 neurons: in cat and anesthetized macaque, binocular thresholds are lower and the iso-response contours in (CL, CR) space are linear. Above threshold, dioptic responses are strongly sublinear. In awake macaque, multiunit binocular responses are about 56% of the summed monocular responses. In superficial-layer 2-photon imaging, binocular responses are about the same as monocular responses, because neurons dominated by one eye are suppressed while balanced neurons are enhanced. The available CRF-parameter evidence points to response-gain (Rmax) effects for interocular imbalance, ocular dominance and interocular suppression, and to contrast-gain (c50) effects for monoptic cross-orientation masking. Per-neuron monocular-vs-dioptic changes in the Naka-Rushton exponent are essentially unreported.

### Cited Findings
**Cat foundations (Freeman/Ohzawa lab)**
- [cat-A] In simple cells, dichoptic drifting gratings at varying interocular phase give phase-specific binocular interaction. At the optimal phase, the response is larger than either monocular response and often larger than their sum. At 180° from optimal it is "strongly inhibited and often completely suppressed." The degree of binocular interaction "does not depend critically on the ocular dominance." Most cells that are strongly dominated by one eye, and even cells that appear monocular, show phase-dependent changes. The data fit a linear binocular-summation model followed by a threshold. A small proportion of apparently monocular simple cells show "a purely inhibitory influence from the silent eye." — [Ohzawa & Freeman 1986a, J Neurophysiol 56:221](https://doi.org/10.1152/jn.1986.56.1.221)
- [cat-A] In complex cells, about 40% are binocularly phase-specific (consistent with linear convergence at subunits) and about 40% are non-phase-specific. 8% of complex cells that look monocular under alternate-eye tests show purely inhibitory, mostly phase-independent influence from the silent eye. — [Ohzawa & Freeman 1986b, J Neurophysiol 56:243](https://doi.org/10.1152/jn.1986.56.1.243)
- [cat-A] Most binocular convergence is linear summation, "however, there is a nonlinear mechanism that maintains binocular interaction regardless of large differences in stimulus strength between the eyes". On this view, a cell dominated by one eye under monocular tests "may respond equally under binocular conditions," so binocular processing could be "independent of the ocular dominance values determined monocularly." — [Freeman & Ohzawa 1990, Vision Res 30:1661](https://doi.org/10.1016/0042-6989(90)90151-a)
- [cat-A] Monocular and binocular CRFs analyzed with ROC neurometrics: binocular contrast thresholds were generally lower and binocular neurometric slopes steeper. The mean single-cell binocular advantage ratio was "slightly higher than the typical ratios reported in human psychophysics." "No single rule appeared to account for the various degrees of binocular summation seen in individual cells." — [Anzai, Bearse, Freeman & Cai 1995, Vis Neurosci 12:77](https://doi.org/10.1017/s0952523800007331)
- [cat-A] Truchard, Ohzawa & Freeman (2000) fit binocular simple cells with a monocular-gain stage per eye, a binocular summation stage with its own gain, and a static half-squaring output. They held one eye's contrast fixed at high (48–51%) or low (4–10%) contrast and varied the other eye's contrast over at least a tenfold range at 8 or more interocular phases. Results: "steep reductions in monocular gain, whereas binocular gain exhibits modest and variable changes," so gain reductions "occur primarily at a monocular site, before convergence." — [Truchard et al. 2000, J Neurosci 20:3017 (abstract)](https://doi.org/10.1523/jneurosci.20-08-03017.2000)
- [cat-A] Truchard 2000 details from the PMC full text:
  - Sample: n = 23 high-fixed-contrast and n = 13 low-fixed-contrast experiments.
  - Gain slopes, high-contrast group: monocular −0.64 ± 0.15 vs binocular −0.05 ± 0.11. Low-contrast group: −0.48 ± 0.32 vs −0.10 ± 0.25.
  - Over a tenfold contrast change, monocular gain fell by about 5× on average; the binocular factor was about 1.2.
  - Ocular dominance did not correlate significantly with the gain slopes.
  - The site of the monocular gain reduction "might arise at a subcortical site," though cortical involvement could not be ruled out. Transient dynamics were not studied.
  — [Truchard et al. 2000, PMC6772209](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- [cat-A] Sengpiel et al. (1998) fit hyperbolic-ratio CRFs for the dominant eye under three kinds of suppression. Cross-orientation inhibition (monoptic orthogonal grating) acted mainly as a rightward shift ("contrast-gain control"). Surround inhibition and interocular suppression (orthogonal grating in the other eye) acted "primarily" as downward shifts ("response-gain control"), with "a component of contrast-gain control." — [Sengpiel, Baddeley, Freeman, Harrad & Blakemore 1998, Vision Res 38:2067](https://doi.org/10.1016/s0042-6989(97)00413-6)
- [cat-A] Sengpiel & Blakemore (1994): a sudden "inappropriate stimulus in one eye can interocularly suppress the activity of cortical neurons if they are already responding to an optimally oriented stimulus in the other eye." In strabismic cats, which lack binocular facilitation, even similarly oriented contours trigger this suppression. — [Sengpiel & Blakemore 1994, Nature 368:847](https://doi.org/10.1038/368847a0)
- [cat-A] Sengpiel, Blakemore & Harrad (1995):
  - Almost half of LGN neurons showed modest interocular inhibition that was independent of orientation.
  - Monocularly driven layer-4 units behaved similarly.
  - In most binocular cortical cells, a very different orientation in the other eye suppressed the response over a wide range of spatial frequency, independent of relative phase.
  - The suppression is nonlinear: it is significant only if the cell is already responding.
  — [Sengpiel, Blakemore & Harrad 1995, Vision Res 35:179](https://doi.org/10.1016/0042-6989(94)00125-6)
- [cat-A] In strabismic cats, only 1 of 50 apparently monocular cells showed subliminal excitation from the silent eye. Most cells showed nonspecific interocular suppression that depended on temporal sequence: it was significant only if the dominant eye was already being driven, and not with simultaneous onset. — [Sengpiel, Blakemore, Kind & Harrad 1994, J Neurosci 14:6855](https://doi.org/10.1523/jneurosci.14-11-06855.1994)

**Macaque, anesthetized (Smith/Chino)**
- [mac-A] Smith et al. (1997), binocular contrast summation: CRFs were measured for optimal dichoptic grating pairs at left:right contrast ratios of 0.1–10. For all cell classes, and for both cooperative and antagonistic phase relations, there was "a linear relationship between the left- and right-eye contrast components required to produce a threshold response." Signals combine linearly at the subunit level before rectification, and "the eye of origin of a light-evoked signal has no specific consequence." — [Smith, Chino, Ni & Cheng 1997, J Neurophysiol 78:366](https://doi.org/10.1152/jn.1997.78.1.366)
- [mac-A] Smith et al. (1997), phase tuning: about 70% of simple cells and about 40% of complex cells were sensitive to interocular phase disparity, with facilitation at the optimal phase and suppression 180° away. "Stimulus contrast had relatively little effect on a cell's disparity tuning, because response saturation occurred at the same contrast level for all relative interocular phase disparities." In other words, disparity tuning was contrast-invariant. — [Smith, Chino, Ni, Ridder & Crawford 1997, J Neurophysiol 78:351](https://doi.org/10.1152/jn.1997.78.1.351)
- [mac-A] Monkeys reared with anisometropic lenses or prisms: many apparently monocular neurons still showed clear binocular interactions under dichoptic stimulation. Contrast signals from the two eyes were still combined in "a normal linear fashion." Prism-reared monkeys had an unusually high proportion of complex cells showing binocular suppression. — [Smith et al. 1997, J Neurophysiol 78:1353](https://doi.org/10.1152/jn.1997.78.3.1353)
- [mac-A] After only 3 days of optically imposed strabismus during the critical period, there was "a striking increase in the prevalence of V1 neurons that exhibited binocular suppression" (binocular response weaker than monocular), while interocular phase-disparity sensitivity was unchanged. The authors argue this suppression originates "beyond where binocular signals are initially combined." — [Zhang et al. 2005, PNAS 102:9026](https://doi.org/10.1073/pnas.0500280102)

**Macaque, awake: dioptic and dichoptic contrast**
- [mac-W] Concordant, zero-disparity dioptic gratings at Michelson contrasts 0, 0.20–0.225, 0.40–0.45 and 0.80–0.90 were recorded with multielectrode arrays.
  - Binocular responses were about 56% of the summed monocular responses at every contrast (regression slopes 0.55, 0.56 and 0.57 for low, medium and high).
  - 219 units (69.7%) were facilitated overall relative to monocular viewing.
  - Facilitation "scaled most proportionately with contrast during the initial transient." Later in the response, facilitation persisted longer at low contrast and turned into binocular suppression at high contrast.
  - Response-gain fits beat contrast-gain fits early (50–150 ms: R² 0.90 vs 0.87). They were comparable at 100–200 ms (0.86 vs 0.84), and response gain was still better late (150–250 ms: 0.81 vs 0.73).
  - The late binocular suppression of the CRF was "better captured by response-gain."
  - The sample was biased toward one of the two monkeys.
  — [Mitchell, Dougherty, Westerberg, Carlson, Daumail, Maier & Cox 2022, iScience 25:104182](https://doi.org/10.1016/j.isci.2022.104182)
- [mac-W] Mitchell et al. (2023) used a 6×6 matrix of contrasts {0, 0.05, 0.11, 0.225, 0.45, 0.90} per eye (36 conditions), with n = 502 units over 25 sessions in 2 monkeys.
  - Contrast in a unit's dominant eye was weighted more than contrast in the non-dominant eye. The comparison held total contrast energy fixed, e.g. 11%/45% vs 45%/11%.
  - Absolute OD predicted binocular response asymmetry (adj. R² = 0.46) and binocular summation asymmetry (adj. R² = 0.16).
  - Non-dominant-eye CRFs were best explained by letting only Rmax differ from the dominant-eye CRF (mean adj. R² 0.92), versus only C50 (0.82) or only the exponent n (0.79).
  - An equivalent contrast weight w in the numerator correlated with OD (adj. R² = 0.56).
  — [Mitchell, Carlson, Westerberg, Cox & Maier 2023, Curr Biol 33:3884](https://doi.org/10.1016/j.cub.2023.08.019)
- [mac-W] Zhang et al. (2024) imaged 10,168 superficial-layer V1 neurons with 2-photon calcium imaging in 3 awake macaques; 9,390 (92.3%) were tuned.
  - Under monocular stimulation, neurons preferring the stimulated eye responded more strongly than binocular neurons.
  - Under binocular stimulation, eye-preferring neurons were on average suppressed and binocular neurons enhanced. This gave "similar neuronal responses irrespective of their eye preferences, and an overall response level similar to that with monocular viewing."
  — [Zhang, Zhao, Jiang, Tang & Yu 2024, eLife 13:RP92839](https://doi.org/10.7554/elife.92839)
- [mac-W] Zhang et al. (2024) summarize earlier electrode work as showing that binocular responses of macaque V1 neurons are "lower than the arithmetic sum of their respective monocular responses" but stronger than the preferred-eye response (citing Prince et al. 2002 and Mitchell et al. 2022). — [Zhang et al. 2024 full text](https://doi.org/10.7554/elife.92839)
- [mac-W] In a large awake V1 random-dot-stereogram sample (n = 787), "energy models make several predictions about the relationship between the response rates to monocular and binocular dot patterns. Few of the predictions were fulfilled exactly," though simple modifications reconcile them. Disparity sensitivity did not correlate with ocular dominance. — [Prince, Pointon, Cumming & Parker 2002, J Neurophysiol 87:191](https://doi.org/10.1152/jn.00465.2000)

**Precortical contributions**
- [mac-W] In the awake macaque LGN, single units of all types (P, M, K) "were significantly suppressed when stimuli were presented at low contrast to the dominant eye and at high contrast to the non-dominant eye." The two eyes interacted antagonistically: suppression weakened with high dominant-eye contrast or low non-dominant-eye contrast. In the paper's discussion, the largest proportion of suppressed neurons was about 25%. — [Dougherty, Carlson, Cox, Westerberg, Zinke, Schmid, Martin & Maier 2021, eNeuro 8:ENEURO.0364-20.2020](https://doi.org/10.1523/eneuro.0364-20.2020)
- [mac-W] Flash-evoked LGN multiunit activity in alert macaque showed binocular suppression as a 15–70% reduction of peak MUA in all laminae except 4. Non-dominant-eye stimulation alone reduced MUA by 10–45% below spontaneous. Some effects occurred at 15–65 ms, before the earliest V1 layer-6 response (30–35 ms), arguing against corticofugal origin for the early effects. — [Schroeder, Tenke, Arezzo & Vaughan 1990, Brain Res 521:303](https://doi.org/10.1016/0006-8993(90)91556-v)
- [marmoset, anesth] In marmoset LGN K layers, binocular flashed stimuli reduced peak rate relative to the dominant eye in 61% (14/23) of cells; 17% showed facilitation and 22% partial summation. "Binocular summation in the marmoset LGN is weak and predominantly sublinear." — [Belluccini, Zeater, Pietersen, Eiber & Martin 2019, Vis Neurosci 36:E012](https://doi.org/10.1017/s0952523819000099)

**Human population-level comparators (for model predictions)**
- [human fMRI] Dichoptic same-orientation gratings evoked more activity than a monocular grating only when contrast was below 10%. Linear summation failed at 10% contrast and above, and "binocular contrast normalization, such that activity from each eye reduced the gain for the other eye," fit well. Cross-orientation suppression was weaker than interocular same-orientation suppression, and onset asynchrony had no effect when attention was diverted. — [Moradi & Heeger 2009, J Vis 9(3):13](https://doi.org/10.1167/9.3.13)
- [human SSVEP] Interocular suppression was strongest in V1 and V3a (about 50% reduction) and weakest in LOC (about 20%). Self-terms and intermodulation terms shared a common divisive gain-control nonlinearity. — [Hou, Nicholas & Verghese 2020, J Neurosci 40:2753](https://doi.org/10.1523/jneurosci.2043-19.2020)
- [human psychophysics] The two-stage model is favored: nearly linear contrast transduction before binocular summation, plus suppressive interocular interactions in the contrast gain control. — [Meese, Georgeson & Baker 2006, J Vis 6:1224](https://doi.org/10.1167/6.11.7)
- [human psychophysics] Gain-control theory: each eye exerts gain control on the other in proportion to its contrast energy, and also on the other eye's gain control. Very low contrasts add linearly, and the nonlinearity appears only as contrast rises. — [Ding & Sperling 2006, PNAS 103:1141](https://doi.org/10.1073/pnas.0509629103)

### Inferences
- Near-threshold linearity (Smith 1997; Anzai 1995) and strong suprathreshold sublinearity (Mitchell 2022's ~0.56 slope; Zhang 2024's binocular ≈ monocular) are both consistent with a normalization stage whose strength grows with total binocular contrast. This matches the low-contrast-only facilitation in human fMRI (Moradi & Heeger) and the Ding–Sperling low-contrast linear regime.
- The CRF-parameter evidence is split by suppression type:
  - Response gain (Rmax): interocular imbalance or OD (Mitchell 2023), orthogonal-in-other-eye interocular suppression in cat (Sengpiel 1998), and the dioptic facilitation and late suppression in awake macaque (Mitchell 2022).
  - Contrast gain (c50): monoptic cross-orientation masking (Sengpiel 1998), and Truchard's monocular-stage gain reduction, which behaves like saturation within each eye's channel.
- The model equations in Mitchell 2023 make diagnostic predictions for balanced dioptic vs monocular CRFs (my algebra from those equations):
  - "Averaging" (eq. 1): no change.
  - Linear sum of monocular CRFs (eq. 2): Rmax doubles, c50 unchanged.
  - Interocular normalization (eq. 4, R = (cLⁿ+cRⁿ)/(cLⁿ+cRⁿ+c50ⁿ)): Rmax unchanged and a leftward c50 shift by exactly 2^(1/n) (e.g., ×0.71 for n = 2).
  - Mitchell 2022's finding that response gain fits the dioptic transient better therefore argues against pure eq. 4 in the transient epoch. This is a sharp test that the 6×6 matrix supports.
- **Testability for the team's dataset:** the 6×6 matrix contains the monocular CRF for each eye (the 0-contrast row and column) and the dioptic diagonal for every unit. That allows a per-unit, per-layer, per-epoch fit of c50, Rmax and n for monocular vs dioptic stimulation, plus dichoptic iso-response contours in (CL, CR). With the matrix's lowest non-zero contrast at 0.055, iso-response contours near threshold may be poorly constrained; Smith 1997 used contrast ratios of 0.1–10, which is comparable. The contours discriminate:
  - linear summation (straight contour; Smith 1997),
  - quadratic/Minkowski summation (convex),
  - winner-take-all (rectangular).

### Gaps
- I did not retrieve Kato, Bishop & Orban (1981) beyond its title ("Binocular interaction on monocularly discharged lateral geniculate and striate neurons in the cat," J Neurophysiol 46:932, [doi](https://doi.org/10.1152/jn.1981.46.5.932)). The publisher page returned 403 and PubMed has no abstract. Later papers cite it to say that only a small minority of V1 neurons respond exclusively to one eye ([Cox et al. 2019](https://doi.org/10.1167/19.12.13)), and as precedent for modeling binocular responses as interocular suppression plus summation ([Chen et al. 2026](https://doi.org/10.7554/elife.107518)). Its quantitative content is unverified here.
- I found no awake-macaque study reporting the per-neuron change in the Naka-Rushton exponent between monocular and dioptic stimulation, and none reporting monocular-vs-dioptic c50/Rmax changes broken down by layer.
- Read & Cumming did not, as far as I found, publish a study varying interocular contrast in V1 single units. Their contrast-relevant contribution is the requirement for a monocular nonlinearity before binocular combination (see Q5).
- Bishop-lab binocular interaction-field papers (e.g., Bishop, Henry & Smith 1971) and Poggio & Fischer 1977 ([doi](https://doi.org/10.1152/jn.1977.40.6.1392); no abstract retrieved) were not examined in detail.

## 2. Maier lab and collaborators (Vanderbilt; Dougherty, Cox, Westerberg, Carlson, Mitchell, Schall et al.), 2018–2026

### Takeaway
Since 2018, Maier and collaborators have built an awake-macaque, laminar (CSD-aligned) case for four claims:
1. Binocular interactions start earlier than the classic model holds, including in monocular L4C neurons and in the LGN.
2. Binocular integration unfolds in two phases: a transient additive phase (about 50–150 ms), then a sustained phase in which concordance and discordance are differentiated and suppression is widespread.
3. Dioptic facilitation is contrast-dependent and turns into suppression at high contrast late in the response.
4. Ocular dominance acts as a response-gain weight on the non-dominant eye within interocular normalization.

The 2023 Current Biology paper used a 6×6 contrast design identical to the team's dataset.

### Cited Findings
- **Dougherty, Cox, Westerberg & Maier (2019), "Binocular Modulation of Monocular V1 Neurons," Curr Biol 29(3):381–391.e4.** [doi:10.1016/j.cub.2018.12.004](https://doi.org/10.1016/j.cub.2018.12.004); bioRxiv [10.1101/320218](https://doi.org/10.1101/320218). [mac-W], 2 *M. radiata*, mirror stereoscope, linear arrays (0.1 mm spacing), layers assigned by CSD relative to the L4/L5 boundary. Key results:
  - Of 290 single units, 78% were binocular (n = 226) and 22% monocular (n = 64). 997 multiunits were also recorded.
  - Monocular neurons' median response to dioptic stimulation was smaller than to the dominant eye alone (n = 33, p = 0.029). The result held under a stricter definition (n = 8, p = 0.0078) and after controlling for fixational eye movements.
  - 11/17 monocular L4 neurons (64.7%) and 11 monocular neurons outside L4 (68.8%) showed significant binocular modulation.
  - Most monocular neurons were suppressed by binocular stimulation; more than 80% (19/23) of these individually significant.
  - Of 10 facilitated monocular neurons, 7 were in L4. Facilitation in L4 "tended to occur at the onset of the visual response."
  - Suppression "occurred several milliseconds following the onset," which the authors attribute to cortical inhibition.
  - A power analysis indicated 80% power with 11 monocular neurons.
  — [Dougherty et al. 2019](https://doi.org/10.1016/j.cub.2018.12.004)
- **Cox, Dougherty, Westerberg, Schall & Maier (2019), "Temporal dynamics of binocular integration in primary visual cortex," J Vis 19(12):13.** [doi:10.1167/19.12.13](https://doi.org/10.1167/19.12.13); preprint "Binocular Integration Manifests as a Transient Spiking Increase Followed by Selective Suppression in Primary Visual Cortex," [10.1101/318667](https://doi.org/10.1101/318667). [mac-W], 2 *M. radiata*, 30 sessions. Key results:
  - Transient (50–100 ms) single-unit responses exceeded monocular responses for both congruent dioptic stimuli (ΔM = +13.73 spikes/s) and incongruent dichoptic stimuli (ΔM = +5.85).
  - In the sustained phase (150–250 ms), congruent facilitation vanished (ΔM = +0.76, n.s.) and incongruent stimuli were suppressed (ΔM = −5.57 spikes/s).
  - Sustained responses to congruent stimuli were enhanced or suppressed depending on ocularity and laminar location.
  - Over 50% of neurons showed suppression after about 150 ms, at about 10 spikes/s, which the authors note matches cat dCOS (Sengpiel et al. 1998).
  — [Cox et al. 2019](https://doi.org/10.1167/19.12.13)
- **Dougherty, Schmid & Maier (2019), "Binocular response modulation in the lateral geniculate nucleus," J Comp Neurol 527:522–534** (review). Binocular modulation of LGN responses is "predominantly suppressive." It might reflect intrinsic LGN circuitry or input from binocular structures (cortex, colliculus). — [doi:10.1002/cne.24417](https://doi.org/10.1002/cne.24417)
- **Dougherty, Carlson, Cox, Westerberg, Zinke, Schmid, Martin & Maier (2021), "Binocular Suppression in the Macaque Lateral Geniculate Nucleus Reveals Early Competitive Interactions between the Eyes," eNeuro 8(2).** [mac-W] single LGN units. The methods list display settings for 66 units, and a reviewer comment quotes N = 56. Findings:
  - P, M and K cells were all suppressed at low dominant-eye / high non-dominant-eye contrast, with an antagonistic two-eye interaction.
  - Onset latency of modulation could not be analyzed because of fixation variability. "More work will be needed to determine the source of binocular modulation in the LGN."
  — [doi:10.1523/eneuro.0364-20.2020](https://doi.org/10.1523/eneuro.0364-20.2020)
- **Mitchell, Dougherty, Westerberg, Carlson, Daumail, Maier & Cox (2022), "Stimulating both eyes with matching stimuli enhances V1 responses," iScience 25(5):104182.** [mac-W], multielectrode arrays, monkeys "E" and "I", 4 contrast levels. Findings:
  - Binocular responses were about 56% of the summed monocular responses.
  - 69.7% (219) of units were facilitated.
  - The transient tracked contrast proportionally. Later, facilitation was sustained at low contrast and became suppression at high contrast.
  - Response-gain models were favored (details in Q1).
  — [doi:10.1016/j.isci.2022.104182](https://doi.org/10.1016/j.isci.2022.104182)
- **Maier, Cox, Westerberg & Dougherty (2022), "Binocular Integration in the Primate Primary Visual Cortex," Annu Rev Vis Sci 8:345–360** (review). "Long-standing consensus held that the two eyes' signals remain largely segregated until they are combined by neurons in the upper layers of the primary visual cortex. However, new experimental data challenge this classic model." The review relates the physiology to psychophysical models of binocular combination. — [doi:10.1146/annurev-vision-100720-112922](https://doi.org/10.1146/annurev-vision-100720-112922)
- **Mitchell, Carlson, Westerberg, Cox & Maier (2023), "A role for ocular dominance in binocular integration," Curr Biol 33(18):3884–3895.e5.** [mac-W]: monkey J (male *M. mulatta*) and monkey B (male *M. radiata*), 25 sessions, 24- or 32-channel 0.1 mm linear probes, mirror stereoscope, 6×6 contrasts {0, .05, .11, .225, .45, .90}.
  - The main analyses use n = 502 single-contact multiunits (404 from monkey J), trial-averaged over 50–200 ms. A small sample of Kilosort/Phy-curated single units appears only in a supplementary figure.
  - Four binocular models were compared: averaging, linear summation, quadratic summation and interocular normalization (R = (cLⁿ+cRⁿ)/(cLⁿ+cRⁿ+c50ⁿ)).
  - Adding a unit-specific OD weight w, derived from monocular data only, improved every model (Wilcoxon p < 0.00001).
  - The mean AIC gain was largest for interocular normalization (+21.5%) and quadratic summation (+19.7%), smaller for summation (+11.7%) and averaging (+9.79%).
  - About half of units (49.2%) preferred the contralateral eye.
  - Limitations as stated: the data "do not directly address the question of why ocular dominance exists," and "it remains to be seen whether ocular dominance columns have a function related to their structure and spatial organization."
  — [doi:10.1016/j.cub.2023.08.019](https://doi.org/10.1016/j.cub.2023.08.019)
- **Carlson, Mitchell, Dougherty, Westerberg, Cox & Maier (2023), "Does V1 response suppression initiate binocular rivalry?," iScience 26(8):107359.** [mac-W], 2 *M. radiata* (E-48, I-34), 14 sessions, binocular rivalry flash suppression (BRFS).
  - BRFS produced "response facilitation rather than response suppression across V1 neurons." A subset of neurons was still reduced, which the authors attribute to adaptation from asynchronous presentation.
  - They propose this selective suppression, rather than generic dCOS, as an alternative initiator of rivalry.
  — [doi:10.1016/j.isci.2023.107359](https://doi.org/10.1016/j.isci.2023.107359)
- **Carlson, Mitchell, Westerberg, Poggi, Woodman & Maier (2025), "V1 interlaminar coherence decreases with interocular conflict," bioRxiv** (preprint, not peer reviewed). [mac-W] laminar MUA and LFP.
  - Interocular conflict "induces modest suppression of single neuron spiking" and "reliably reduces interlaminar coherence, particularly between granular and infragranular layers."
  - During BRFS, coherence stayed reduced even when firing rates were unchanged, and was higher for perceptually dominant stimuli.
  — [doi:10.1101/2025.10.04.680461](https://doi.org/10.1101/2025.10.04.680461)
- Other lab context, author sets only: Daumail et al. 2023, "Rapid adaptation of primate LGN neurons to drifting grating stimulation," J Neurophysiol 129:1447 ([doi](https://doi.org/10.1152/jn.00058.2022)), with Mitchell, Carlson, Martin, Tong, Maier and Dougherty. Earlier lab rivalry work: Maier et al. 2008, "Divergence of fMRI and neural signals in V1 during perceptual suppression in the awake monkey," Nat Neurosci 11:1193 ([doi](https://doi.org/10.1038/nn.2173)). Neither abstract was examined in detail.
- Outside the lab, Zhang et al. (2024) explicitly discuss a discrepancy with Dougherty et al. (2019): "they found that monocular and binocular neurons exhibit similar responses under monocular stimulation, with only monocular neurons, but not binocular neurons, being significantly suppressed by binocular stimulation". In contrast, Zhang et al. found much stronger preferred-eye responses in monocular than in binocular neurons, and enhancement of binocular neurons. They suggest the difference may reflect temporal resolution (spikes at 50–200 ms vs calcium signals integrated up to 1000 ms) or small, biased electrode samples. — [Zhang et al. 2024](https://doi.org/10.7554/elife.92839)
- An eLife reviewer of Zhang et al. (2024) wrote: "Future work would do well to combine the spatial power of 2P imaging with the temporal power of electrophysiology to assess ocular dominance-dependent binocular combination across the V1 laminar microcircuit." — [Zhang et al. 2024 reviews](https://doi.org/10.7554/elife.92839)

### Inferences
- **The team's dataset appears to be the Mitchell et al. 2023 dataset, or a superset of it.** The match is on several points:
  - Same monkey designations and species (J = male *M. mulatta*, B = male *M. radiata*).
  - Identical contrast set {0, .05/.055, .11, .225, .45, .90} in a 6×6 design.
  - Mirror stereoscope and 24/32-channel probes at 0.1 mm spacing.
  - 25 sessions vs the team's roughly 30 penetrations.

  Therefore "does OD weight the non-dominant eye, and does interocular normalization fit best (time-averaged, pooled MUA)" is already published. New contributions would need to add one or more of:
  - laminar resolution, i.e. CSD-assigned L4C vs supra- vs infragranular;
  - time resolution, i.e. transient vs sustained epochs within the ~317 ms stimulus;
  - isolated-single-unit analyses (about 70 units);
  - full CRF parameter comparisons, monocular vs dioptic;
  - two-probe, between-column analyses;
  - LFP/CSD/coherence measures.
- The earlier Maier-lab papers (Dougherty 2019; Cox 2019; Carlson 2023) used different monkeys (E48/I34) and mostly single high-contrast stimuli. Mitchell 2022 used only 4 contrast levels and no dichoptic imbalance. None of these crossed laminar position with a full dichoptic contrast matrix, which is the cell the team's dataset fills.
- Given the "correct over consistent" policy, re-deriving the Mitchell 2023 results with the rebuilt pipeline is a validation step and may legitimately change the numbers. Known timing issue: the 3.4 ms NSx offset noted in project docs.

### Gaps
- I did not find any 2024–2026 peer-reviewed Maier-lab paper specifically on laminar binocular CRFs or laminar normalization. The 2025 coherence preprint is the most recent lab item indexed in Europe PMC. An unindexed preprint by B. Mitchell on laminar CRFs may exist; I found none.
- I could not determine from the abstract whether Maier et al. 2022 (Annual Review) gives quantitative laminar summaries beyond the primary papers; the full text was not open.
- I did not obtain exact single-unit sample sizes per layer for Mitchell 2022.

## 3. Laminar differences: monocularity of L4C, where binocular combination first occurs, supra- vs infragranular summation and normalization, ocular dominance, contrast gain, CSD evidence

### Takeaway
Classic anatomy and physiology make L4C (4Cα magno-recipient, 4Cβ parvo-recipient) and 4A mostly monocular and untuned, with binocular cells dominating layers 2/3, 4B, 5 and 6. Awake laminar work now shows that roughly two-thirds of "monocular" L4 neurons are modulated by the other eye. Early facilitation in L4 is consistent with converging thalamocortical input, while the later, more common suppression is consistent with intracortical inhibition. Binocular stimulation also enlarges the initial granular CSD sink. Laminar evidence on how normalization differs by layer is thin. In cat, Sengpiel et al. (1998) found response-gain interocular suppression in L4 vs contrast-gain suppression in L2/3, 5 and 6. In macaque, Cox et al. (2019) report layer-dependent sustained facilitation, but with very small per-layer samples (n = 4–7).

### Cited Findings
- [mac, lightly anesth] Hubel & Wiesel (1968):
  - "In layer IV most cells are driven by one eye only," forming an eye-alternating mosaic.
  - Layers II–III and V–VI contain mostly binocular complex and hypercomplex cells.
  - Ocular dominance columns run surface to white matter, independent of orientation columns.
  — [Hubel & Wiesel 1968, J Physiol 195:215](https://doi.org/10.1113/jphysiol.1968.sp008455)
- [mac-A] Blasdel & Fitzpatrick (1984): layers 4A, 4Cα and 4Cβ are "populated mostly by monocularly driven cells having small receptive fields and lacking orientation selectivity." 4B cells are orientation- and often direction-selective. 4Cα cells have lower contrast thresholds than 4Cβ cells; 4A has a mixture. — [Blasdel & Fitzpatrick 1984, J Neurosci 4:880](https://doi.org/10.1523/jneurosci.04-03-00880.1984)
- [mac-A] Neurons in IVcα (magno-recipient) have high contrast sensitivity and IVcβ (parvo-recipient) low. IVa and IVb contain mixtures. Non-oriented cells occur only in parvo-recipient IVa and IVcβ. — [Hawken & Parker 1984, Exp Brain Res 54:367](https://doi.org/10.1007/bf00236238)
- [mac-A] Contrast sensitivity functions become progressively steeper from LGN to V1 to MT. — [Sclar, Maunsell & Lennie 1990, Vision Res 30:1](https://doi.org/10.1016/0042-6989(90)90123-3)
- [mac-W] In alert macaques, zones with ongoing activity correspond to the cytochrome-oxidase-rich geniculorecipient layers 4A, 4C and 6. Layers 2/3, 4B and 5 have little ongoing activity, and about 90% of cells there are orientation-selective. — [Snodderly & Gur 1995, J Neurophysiol 74:2100](https://doi.org/10.1152/jn.1995.74.5.2100)
- [mac-W] Dougherty et al. (2019), laminar results:
  - Monocular neurons in all layers, "including the input layers," distinguish dominant-eye from binocular stimulation (11/17 in L4).
  - 7 of the 10 facilitated monocular neurons were in L4.
  - Visual onset order was L4 (42.1 ms), then L2/3 (45.2 ms), then L5/6 (46.4 ms).
  - Methods: the L4/L5 boundary was set from the CSD sink to a white flash (33 penetrations) or grating (45 penetrations). The supragranular/granular boundary was "set to 0.5 mm above" the granular/infragranular boundary because it is hard to define from CSD.
  — [Dougherty et al. 2019](https://doi.org/10.1016/j.cub.2018.12.004)
- [mac-W] Cox et al. (2019), CSD: the initial granular current sink (50–100 ms) was larger for both congruent (ΔM ≈ −480 nA/mm²) and incongruent (≈ −495) binocular stimuli than for monocular stimuli. In 150–250 ms, the congruent sink approached the monocular level (n.s.), while the incongruent sink remained significantly larger (−337). — [Cox et al. 2019](https://doi.org/10.1167/19.12.13)
- [mac-W] Cox et al. (2019), equiocular tuned single units by layer (upper ≥0.6 mm, middle 0.1–0.5 mm, deep ≤−0.1 mm relative to the L4/5 boundary; n = 7, 5 and 4 per figure legend):
  - Transient facilitation: upper ΔM = 9.04 (p = 0.039), middle 28.26 (p = 0.030), deep 36.23 (p = 0.088).
  - Sustained facilitation: absent in upper layers (0.19, n.s.), persisting in middle layers (16.68, p = 0.046).
  — [Cox et al. 2019 Table 2](https://doi.org/10.1167/19.12.13)
- [cat-A] "Response-gain control prevailed in layer 4, whereas cells in layers 2/3, 5 and 6 mainly showed contrast-gain control" (pooled across the three suppression types tested, including interocular suppression). — [Sengpiel et al. 1998](https://doi.org/10.1016/s0042-6989(97)00413-6)
- [cat-A] Monocularly driven layer-4 units in area 17 showed modest, orientation-independent interocular inhibition, like LGN cells. Binocular cells showed orientation-dependent interocular suppression. — [Sengpiel, Blakemore & Harrad 1995](https://doi.org/10.1016/0042-6989(94)00125-6)
- [cat-A] Contrast gain control (adaptation) showed "no obvious differences regarding cell type (simple vs. complex) or laminar distribution." It is almost certainly cortical, because LGN cells show minimal effects. — [Ohzawa, Sclar & Freeman 1985, J Neurophysiol 54:651](https://doi.org/10.1152/jn.1985.54.3.651)
- [cat-A] Interocular transfer of adaptation was about 55% at the population level, present in binocular and monocular neurons alike. A cell's position relative to OD column borders "had no significant effect" on IOT strength. — [Howarth, Vorobyov & Sengpiel 2009, Cereb Cortex 19:1835](https://doi.org/10.1093/cercor/bhn211)
- [mac-W] Laminar MUA and LFP to drifting gratings: sustained MUA is roughly constant across depth. Gamma-band LFP peaks in output layers 2, 3 and 4B, which account for 67% of V1 gamma power and 56% of evoked spikes. — [Xing, Yeh, Burns & Shapley 2012, PNAS 109:13871](https://doi.org/10.1073/pnas.1201478109)
- [mac-W] Interocular conflict reduces interlaminar LFP coherence, "particularly between granular and infragranular layers." — [Carlson et al. 2025 bioRxiv](https://doi.org/10.1101/2025.10.04.680461)
- [mac-W, 2P] Zhang et al. (2024) is limited to superficial layers. The authors note that "binocular neurons also exist in deeper layers, and even neurons in the input layer are affected by feedback from downstream binocular neurons" (citing Dougherty 2019). — [Zhang et al. 2024](https://doi.org/10.7554/elife.92839)
- [model, macaque 4Cα] A multiscale network model reproduces "narrow binocular strips along borders of ocular dominance columns." The data are best matched when 10–30% of interactions near OD-column boundaries are cross-columnar, and L6 feedback is "largely monocular." (Preprint.) — [Xiao, Lin & Young 2026 bioRxiv](https://doi.org/10.64898/2026.06.19.733440)
- [human 7T fMRI, layer/column-resolved] "Eye-specific rivalry arises from interocular inhibition between adjacent ocular dominance columns in the superficial layers of the primary visual cortex, but not between ocular layers of the lateral geniculate nucleus." — [Qian et al. 2026, Nat Hum Behav 10:403](https://doi.org/10.1038/s41562-025-02320-4)
- [human 7T + EEG, amblyopia] Monocular activity was attenuated and delayed in the thalamic-input layers of V1, followed by "imbalanced binocular suppression and weakened binocular integration in the superficial layers." — [Wang et al. 2025, Imaging Neurosci 3:imag_a_00561](https://doi.org/10.1162/imag_a_00561)
- [mouse, awake] Binocular integration is widespread across layers, in both regular- and fast-spiking units. Ocular dominance, interocular matching and disparity selectivity are mutually independent. — [Fu, Tanabe & Cang 2023, J Neurosci 43:6495](https://doi.org/10.1523/jneurosci.0925-23.2023)

### Inferences
- Under the classic model, L4C neurons get eye-specific LGN input and binocular combination first happens in 4B, 2/3 and the deep layers. The current evidence instead implies combination at three levels:
  1. Weak interaction in the LGN: suppression in about 25% of awake LGN units, and 15–70% MUA reduction.
  2. Interaction within L4C: early facilitation in a minority of cells and later suppression in a majority. The facilitation fits the binocular strips at OD-column borders and cross-columnar 4C dendrites (model), or thalamocortical convergence.
  3. Full excitatory combination in 2/3, 4B and 5/6.
- Laminar contrast gain differs by sublayer (4Cα vs 4Cβ), so a 32-channel penetration spanning L4C can test whether magno- and parvo-recipient monocular populations differ in binocular modulation. Low-c50 units in upper L4C would be the magno-recipient candidates. Dougherty 2021 found LGN suppression across P, M and K cells, which predicts modulation in both sublayers.
- A larger granular CSD sink under binocular stimulation (Cox 2019) is a direct, high-SNR laminar readout. Across the 6×6 matrix it can show whether granular synaptic input already reflects interocular normalization: whether the sink's contrast dependence is better fit by eq. 4 than by summation.
- **Testability:** good for MUA and CSD across layers. Single-unit power in L4C is the limiting factor: truly monocular L4C single units will be rare, since only 17 L4 monocular single units across Dougherty's sample were enough for significance, and the team has about 70 single units total.

### Gaps
- I found no primate laminar study comparing monocular vs dioptic CRF parameters (c50, Rmax, n) by layer. The macaque laminar evidence on binocular normalization is limited to Cox 2019 (n = 4–7 per layer) and Dougherty 2019 (monocular units only).
- I found no macaque evidence directly testing Sengpiel 1998's cat result (response gain in L4 vs contrast gain elsewhere for interocular suppression).
- I found no single-unit macaque study separating 4Cα and 4Cβ for binocular modulation.
- How binocular modulation is distributed in L6 (corticogeniculate output) vs L5 has not been reported separately in macaque.
- Koniocellular input to layers 1–3 (CO blobs) and its binocularity were not covered by any retrieved V1 study.

## 4. Timing: latency of binocular facilitation vs suppression; transient vs sustained epochs

### Takeaway
Across preparations, binocular facilitation (or additive combination) appears at or near response onset, and suppression is delayed. The delay is about 7 ms for disparity-specific suppressive subunits (Tanabe & Cumming 2014) and about 12 ms median for binocular modulation of monocular neurons (Dougherty 2019: modulation onset 62 ms vs response onset 50 ms). At the population level, dichoptic suppression and the distinction between concordant and discordant stimuli emerge only in a sustained phase starting about 150 ms after onset (Cox 2019). For dioptic gratings, high contrast speeds the transition from facilitation to suppression (Mitchell 2022). Cat work adds that interocular suppression can depend on stimulus sequence: it is strongest when the suppressed eye is already driving the cell. This was not replicated in awake strabismic macaques.

### Cited Findings
- [mac-W] Dougherty et al. (2019), latencies:
  - Median visual response onset was 50 ms; median binocular modulation onset was 62 ms (p = 2.01 × 10⁻¹⁷).
  - Modulation onset by layer: L2/3 47.6 ms, L4 50.1 ms, L5/6 54.7 ms.
  - In L4, the 0.25-quantile modulation onset was 49 ms for facilitated and 54 ms for suppressed neurons.
  - The authors interpret L4 facilitation at response onset as converging thalamocortical input, and the later suppression as "cortical inhibition."
  — [Dougherty et al. 2019](https://doi.org/10.1016/j.cub.2018.12.004)
- [mac-W] Cox et al. (2019) describe a transient phase "lasting 50–150 ms from stimulus onset" with enhanced responses to both congruent and incongruent binocular stimuli. A sustained phase follows, with suppression for incongruent stimuli and ocularity- and layer-dependent facilitation or suppression for congruent stimuli. The authors note that rivalry requires at least about 150 ms of exposure, which matches this timing. — [Cox et al. 2019](https://doi.org/10.1167/19.12.13); [preprint](https://doi.org/10.1101/318667)
- [mac-W] Mitchell et al. (2022): binocular facilitation "scaled most proportionately with contrast during the initial transient." Later, "longer maintenance of facilitation for lower contrast and binocular suppression at high contrast." This supports "initial response enhancement followed by contrast-dependent control of excitation." — [Mitchell et al. 2022](https://doi.org/10.1016/j.isci.2022.104182)
- [mac-W] With binocular noise and subspace mapping, suppression to the opposite interocular phase (at lower spatial frequencies) "was delayed relative to the excitation, with a median difference in latency of 7 ms." — [Tanabe & Cumming 2014, J Neurophysiol 111:1759](https://doi.org/10.1152/jn.00426.2013)
- [mac-W] Spike-triggered covariance on binocular white noise reveals excitatory and suppressive binocular subunits arranged "in a specific push-pull manner for disparity," which reduces responses to false matches. — [Tanabe, Haefner & Cumming 2011, J Neurosci 31:8295](https://doi.org/10.1523/jneurosci.5000-10.2011)
- [cat-A] Interocular suppression "exhibits a remarkable non-linearity: a grating of non-preferred orientation in one eye causes significant interocular suppression only if the neuron is already responding to an appropriate stimulus in the other eye." — [Sengpiel, Blakemore & Harrad 1995](https://doi.org/10.1016/0042-6989(94)00125-6). In strabismic cats this suppression appeared only with sequential onset, not simultaneous. — [Sengpiel et al. 1994](https://doi.org/10.1523/jneurosci.14-11-06855.1994)
- [mac-W, strabismic] Simultaneous dichoptic presentation in 500 single neurons from 2 exotropic monkeys produced "no net interocular suppression or facilitation." — [Economides, Adams & Horton 2021, J Neurosci 41:5522](https://doi.org/10.1523/jneurosci.0044-21.2021). Staggered presentation (5 s offsets) gave mean interaction indices of 0.48 ± 0.05 and 0.46 ± 0.08 (0.50 = no suppression; n = 82), possibly explained by adaptation. — [Economides et al. 2021, J Neurophysiol 126:1101](https://doi.org/10.1152/jn.00275.2021)
- [mac-W] BRFS (about 1 s asynchrony) produced facilitation, not suppression, across V1 neurons. A subset was reduced through adaptation to the asynchronous presentation. — [Carlson et al. 2023](https://doi.org/10.1016/j.isci.2023.107359)
- [mac-W, LGN MUA] In lamina 3, the onset of non-dominant-eye suppression "often preceded that of the excitatory response to stimulation of the dominant eye." Early effects (15–65 ms) precede V1 L6 responses (30–35 ms). — [Schroeder et al. 1990](https://doi.org/10.1016/0006-8993(90)91556-v)
- [mac-W, VSD] With dichoptic patterns, brief local monocular contrast increments triggered traveling waves across V1 at about 30 mm/s. — [Yang, Heeger, Blake & Seidemann 2015, J Neurophysiol 113:277](https://doi.org/10.1152/jn.00610.2013)
- [mouse] Interocular phase disparity reduced early firing (40–80 ms), while orientation conflict increased late firing (100–200 ms) together with suppression of SST interneurons. The authors infer feedforward vs recurrent mechanisms. — [Montgomery et al. 2026, Front Syst Neurosci 20:1786396](https://doi.org/10.3389/fnsys.2026.1786396)
- [cat-A] Truchard et al. (2000) characterized only time-averaged (4 s) gain control and state that its "transient dynamics" were not investigated. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)

### Inferences
- The two timescales may reflect different mechanisms:
  - Suppression delayed by about 7–12 ms relative to excitation fits a single extra synapse, i.e. local feedforward inhibition. That is the expected cost of disynaptic inhibition.
  - The ~150 ms sustained-phase differentiation fits recurrent or lateral (horizontal or traveling-wave) interactions, or feedback.
- The team's ~317 ms stimuli allow the same split as Cox 2019: transient about 40–150 ms, sustained about 150–317 ms. The sustained window is about 170 ms, comparable to Cox's 150–250 ms.
- Across the 6×6 matrix, the analysis can ask:
  1. Does suppression latency shorten with contrast in the non-dominant eye, as expected for a normalization pool driven by that eye?
  2. Does it differ by layer, e.g. earliest in L2/3 as for modulation onset in Dougherty 2019?
  3. Is the transient fit by summation (eq. 2) and the sustained by interocular normalization (eq. 4)?
- **Caution from project documentation:** NSx (.ns6) windowed reads are offset about 3.4 ms from the NEV event clock in this pipeline (docs/01 §6). That is about half of Tanabe & Cumming's 7 ms excitation–suppression delay and about 30% of Dougherty's 12 ms median. Latencies from continuous data (aMUA, CSD) and NEV-aligned spikes must be put on the same clock before comparing facilitation and suppression onsets.
- Simultaneous onset in both eyes means the team cannot test the Sengpiel "already responding" sequence effect directly. One exception: if successive presentations within a trial (3 per trial) change one eye's contrast while the other persists, sequential effects are partly accessible. Whether the task's sequences do this is not verified.

### Gaps
- I found no awake-macaque study measuring how the latency of interocular suppression depends on contrast across a full dichoptic matrix.
- No macaque study has used CSD to locate where the delayed suppression arrives (granular vs supragranular sinks and sources).
- How the delay relates to cell type (narrow- vs broad-spiking) is untested in macaque.

## 5. Eye-specific gain control: monocular (pre-combination) vs binocular (post-combination) normalization pools; dichoptic masking in single units

### Takeaway
The evidence supports at least two suppressive stages:
- A fast, broadband monocular or precortical stage. It mediates monoptic cross-orientation suppression and Truchard's monocular gain reduction, works at high temporal frequency, resists adaptation, and is largely explained by LGN saturation and rectification plus threshold.
- A cortical, GABA-dependent, adaptable interocular stage. It mediates dichoptic cross-orientation and interocular suppression, prefers low temporal frequency, and is abolished by adaptation and reduced by bicuculline.

In awake macaque, interocular normalization weighted by ocular dominance fits V1 contrast data best (Mitchell 2023; Zhang 2024). Georgeson & Sengpiel's re-analysis of cat data places interocular-suppression signals before OD attenuation and the contrast-gain-control pool after it, i.e. a two-stage model. One awake study (Macknik 2004) dissents, finding inhibition in V1 to be mostly monoptic.

### Cited Findings
- [cat-A] Binocular cross-orientation suppression was "as strong as monocular COS, but disparity changes are of no consequence." Interocular (dichoptic) effects were mostly weak but "present in almost every cell," and onset asynchrony did not matter. The authors concluded COS is "primarily monocular and acts prior to the convergence of the two monocular streams." — [Walker, Ohzawa & Freeman 1998, J Neurophysiol 79:227](https://doi.org/10.1152/jn.1998.79.1.227)
- [cat-A] Cross-orientation suppression originates within the classical receptive field, occurs at all orientations, and is largely phase-independent. Its dichoptic component was examined in 3 binocular cells; that result is truncated in the retrieved abstract. — [DeAngelis, Robson, Ohzawa & Freeman 1992, J Neurophysiol 68:144](https://doi.org/10.1152/jn.1992.68.1.144)
- [cat-A] Monoptic COS survives masks drifting too fast to drive cortex, and is immune to adaptation-induced hyperpolarization of mask-driven neurons. It is therefore likely subcortical, through thalamic suppression and thalamocortical depression. — [Freeman, Durand, Kiper & Carandini 2002, Neuron 35:759](https://doi.org/10.1016/s0896-6273(02)00819-x)
- [cat-A, intracellular] Cross-oriented stimuli suppressed both excitation and inhibition. LGN contrast saturation and rectification in a feedforward model predict F1 suppression, and spike-threshold amplification nearly doubles the suppression. — [Priebe & Ferster 2006, Nat Neurosci 9:552](https://doi.org/10.1038/nn1660)
- [cat-A] "Dichoptic COS is best elicited with lower temporal frequencies and is substantially reduced after prolonged adaptation to a mask grating," whereas monoptic COS prefers high temporal frequency and resists adaptation. Conclusion: monoptic COS is subcortical and dichoptic COS intracortical. — [Li, Peterson, Thompson, Duong & Freeman 2005, J Neurophysiol 94:1645](https://doi.org/10.1152/jn.00203.2005)
- [cat-A] Interocular suppression "cannot be evoked by very fast-moving stimuli," is "essentially" eliminated by adaptation to the orthogonal suppressor, and "is greatly reduced in the presence of the GABA antagonist bicuculline." It is "mediated by inhibitory circuitry within the visual cortex." — [Sengpiel & Vorobyov 2005, J Neurosci 25:6394](https://doi.org/10.1523/jneurosci.0862-05.2005)
- [cat-A re-analysis, model] Georgeson & Sengpiel (2021) compared six two-stage variants, with divisive contrast gain control (CGC) before and after binocular summation plus interocular suppression (IOS). In the best model, "signals driving IOS arise 'early', before attenuation due to OD, while suppressive CGC signals are 'late' and so affected by OD." Late CGC "yields partial compensation for OD," which reduces ocular imbalance at the input to binocular summation. — [Georgeson & Sengpiel 2021, Vision Res 185:29](https://doi.org/10.1016/j.visres.2021.03.004)
- [cat-A] Truchard et al. (2000) found contrast-gain reductions mostly in each eye's monocular channel, with binocular-stage gain changing only modestly. The authors could not tell whether the binocular-stage site was the exponent, an offset or the binocular gain. — [Truchard et al. 2000](https://doi.org/10.1523/jneurosci.20-08-03017.2000); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- [cat-A] Interocular suppression acts on the CRF mostly as response gain, while monoptic cross-orientation suppression acts as contrast gain. — [Sengpiel et al. 1998](https://doi.org/10.1016/s0042-6989(97)00413-6)
- [mac-W] Single units in awake monkey V1 show dichoptic masking in two separate circuits, one binocular and one monocular, "possibly even subcortical." V1 neurons "are dichoptic in terms of excitation, but monoptic in terms of inhibition": excitation to monocular targets was "inhibited strongly only by masks presented to the same eye, and not by masks presented to the opposite eye." — [Macknik & Martinez-Conde 2004, J Cogn Neurosci 16:1049](https://doi.org/10.1162/0898929041502788)
- [mac-W] Dichoptic (incongruent, orthogonal) stimulation suppressed sustained responses in over 50% of neurons by about 10 spikes/s, after a transient enhancement. — [Cox et al. 2019](https://doi.org/10.1167/19.12.13)
- [mac-W] Including an OD weight on the non-dominant eye improved interocular normalization most (AIC +21.5%), and this was the best-fitting model overall. — [Mitchell et al. 2023](https://doi.org/10.1016/j.cub.2023.08.019)
- [mac-W, 2P] Zhang et al. (2024) proposed a model in which each eye's monocular response is normalized by OD-weighted divisive interocular inhibition, then the two are summed (wc = (ODI+1)/2, wi = 1 − wc, with a divisive constant k). — [Zhang et al. 2024](https://doi.org/10.7554/elife.92839)
- [mac-W, 2P] Continuous flash suppression (3,564 neurons; 3,004 orientation-tuned) "strongly suppressed V1 orientation responses in an ocular-dominance-dependent manner." It "nearly abolish[ed] responses in neurons preferring the masker eye or both eyes," and "significantly" reduced responses in grating-eye-preferring neurons. — [Chen, Wang, Jiang, Tang & Yu 2026, eLife 14:RP107518](https://doi.org/10.7554/elife.107518)
- [mac-A, Utah arrays] In normal and amblyopic *M. nemestrina*, a noise mask in the fellow eye attenuated responses to amblyopic-eye gratings, "as in a normal control animal." The reverse masking failed in amblyopes, an "imbalance of suppression that favors the dominant eye." — [Shooner et al. 2017, J Neurosci 37:8734](https://doi.org/10.1523/jneurosci.1760-17.2017)
- [mac-A] OD was computed separately for receptive-field excitation and suppression. Excitatory drive was dominated by the fellow eye, but "suppression from both the FE and AEs was prevalent in all animals." The authors describe "balanced binocular suppression that is unaltered in amblyopia." — [Hallum et al. 2017, J Neurosci 37:8216](https://doi.org/10.1523/jneurosci.0449-17.2017); [Shooner et al. 2017](https://doi.org/10.1523/jneurosci.1760-17.2017)
- [mac-W] The energy model failed for anticorrelated stereograms, where response amplitude is reduced rather than simply inverted ([Cumming & Parker 1997, Nature 389:280](https://doi.org/10.1038/38487)). The failure is fixed by passing each eye's input "through a monocular simple cell, incorporating a threshold, before being combined at a binocular simple cell." This is evidence for a monocular nonlinearity before binocular combination. — [Read, Parker & Cumming 2002, Vis Neurosci 19:735](https://doi.org/10.1017/s0952523802196052)
- [mac-W] Ocular dominance "predicts neither strength nor class of disparity selectivity" in awake V1 with random-dot stimuli. — [Read & Cumming 2004, J Neurophysiol 91:1271](https://doi.org/10.1152/jn.00588.2003)
- [human fMRI] Mutual suppression between orthogonal orientations "does not depend on the eye of origin, nor on the onset times," and cross-orientation suppression is weaker than interocular same-orientation suppression. — [Moradi & Heeger 2009](https://doi.org/10.1167/9.3.13)
- [mouse, preprint] PV interneurons' binocular responses track pyramidal cells with sublinear summation (inhibition-stabilized-network-like). SST interneurons show linear or supralinear summation and fast dynamics, so normalization in binocular V1 depends on cell type. — [Perumal et al. 2026 bioRxiv](https://doi.org/10.64898/2026.09.20.753031)

### Inferences
- The cat literature consistently splits suppression into a monocular, fast, untuned, adaptation-resistant component (possibly subcortical) and an interocular, cortical, GABAergic, adaptable component. In a contrast-only paradigm with iso-oriented gratings, these predict different signatures:
  - Monocular-stage normalization acts within each eye's channel before summation, as in Truchard and the Meese two-stage model. It should appear as c50 shifts of each eye's own CRF, independent of the other eye.
  - Interocular normalization means the other eye's contrast enters the denominator. It should show up as a non-dominant-eye contrast that reduces the response to dominant-eye contrast, most strongly when the non-dominant eye is at high contrast and the dominant eye at low contrast (the same corner where LGN suppression appears in Dougherty 2021).
- Macknik 2004's "monoptic inhibition" conflicts with Cox 2019, Zhang 2024, Chen 2026 and Shooner 2017 in normal controls, all of which find substantial dichoptic suppression in macaque V1. The differences may be stimulus-driven (brief flashed targets and masks in Macknik vs sustained gratings, noise or CFS elsewhere). This is unresolved.
- Because gratings in the team's dataset are presumably identical in orientation across eyes, any interocular suppression they show is "same-orientation" interocular normalization, i.e. the Moradi & Heeger type. It cannot be compared directly with dCOS literature that uses orthogonal masks.

### Gaps
- No primate study I found used pharmacology (e.g., GABA antagonists) or adaptation to separate monocular from interocular normalization in V1 single units. The decisive manipulations are cat-only (Sengpiel & Vorobyov 2005; Li 2005; Freeman 2002).
- No primate study has tested the temporal-frequency dependence of dichoptic vs monoptic masking.
- The site of Truchard's "monocular" gain reduction (LGN, thalamocortical synapse, or monocular L4C neurons) remains unresolved.
- No study has directly compared the normalization pool's ocular composition with the neuron's excitatory OD in macaque L4C vs other layers. Hallum 2017 did so for excitation vs suppression with Utah arrays, which have no laminar resolution.

## 6. Controversies and explicitly stated open questions (with discriminating measurements and testability for the team's dataset)

### Takeaway
The live controversies are:
- where the eyes first interact (LGN vs L4C vs upper layers);
- whether contrast gain control is monocular (Truchard) or interocular/binocular (Moradi & Heeger; Mitchell 2023), or both (two-stage);
- whether ocular dominance is functionally meaningful for binocular combination (Freeman & Ohzawa 1990; Read & Cumming 2004 vs Mitchell 2023; Zhang 2024);
- whether binocular stimulation suppresses or enhances monocular and eye-preferring neurons (Dougherty 2019 vs Zhang 2024);
- whether interocular suppression in V1 is substantial and early (Cox 2019; CFS 2026) or weak and monoptic (Macknik 2004);
- whether V1 suppression has a perceptual role in rivalry and strabismus (Economides 2021; Carlson 2023; Keliris 2010).

Several of these can be settled with a laminar, time-resolved analysis of the full dichoptic contrast matrix.

### Cited Findings
- Stated open question: "Future work may determine the degree of binocular modulation in primate LGN as well as the role of corticofugal feedback in its implementation." — [Dougherty et al. 2019](https://doi.org/10.1016/j.cub.2018.12.004)
- Stated open question: "More work will be needed to determine the source of binocular modulation in the LGN." The modulation could come from intrinsic LGN circuits, corticogeniculate feedback, or superior colliculus. Fixation variability prevented latency analysis. — [Dougherty et al. 2021](https://doi.org/10.1523/eneuro.0364-20.2020)
- Stated controversy: "Controversy surrounds whether ocular dominance... is functionally relevant." Open question: "It remains to be seen whether ocular dominance columns have a function related to their structure and spatial organization." Future work is to understand "the interplay between OD/binocularity, contrast gain-control, and interocular interactions." — [Mitchell et al. 2023](https://doi.org/10.1016/j.cub.2023.08.019)
- Contrary positions on OD: "binocular processing for all cortical cells could be uniform and independent of the ocular dominance values determined monocularly" ([Freeman & Ohzawa 1990](https://doi.org/10.1016/0042-6989(90)90151-a)); the degree of binocular interaction "does not depend critically on the ocular dominance" ([Ohzawa & Freeman 1986a](https://doi.org/10.1152/jn.1986.56.1.221)); OD predicts neither strength nor class of disparity tuning ([Read & Cumming 2004](https://doi.org/10.1152/jn.00588.2003)); OD is uncorrelated with disparity sensitivity ([Prince et al. 2002](https://doi.org/10.1152/jn.00465.2000)); OD is uncorrelated with gain slopes ([Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)).
- Supporting positions on OD: OD predicts the contrast weight ([Mitchell et al. 2023](https://doi.org/10.1016/j.cub.2023.08.019)), the sign of binocular modulation in 2P ([Zhang et al. 2024](https://doi.org/10.7554/elife.92839)), and suppression under CFS ([Chen et al. 2026](https://doi.org/10.7554/elife.107518)).
- Stated open questions (Zhang et al. 2024):
  - "It remains unclear whether the discrepancies observed [vs Dougherty 2019] are caused by differences in temporal resolutions." "We cannot rule out the possibility that neuronal response changes beyond 200ms may contribute."
  - Future work should test "uneven effective stimulus contrasts... monocular adaptation... short-term monocular deprivation."
  - 2P is limited to superficial layers.
  — [Zhang et al. 2024](https://doi.org/10.7554/elife.92839)
- Stated caveats (Truchard et al. 2000): the site of monocular gain reduction "might arise at a subcortical site"; transient dynamics were not investigated; the binocular-stage site could not be resolved among exponent, offset and binocular gain. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- Stated question (Cox et al. 2018 preprint): "One open question is whether binocular integration occurs at a single stage of sensory processing or in a sequence of computational steps." Their answer was at least two steps. — [Cox et al. 2018 bioRxiv](https://doi.org/10.1101/318667)
- Stated contingency (Mitchell 2022): "the contrast dependency of this dynamic remains unknown" (the authors' framing before their study). Limitation: data "biased towards one" subject. — [Mitchell et al. 2022](https://doi.org/10.1016/j.isci.2022.104182)
- Rivalry and awareness:
  - About 20% of V1 neurons reflect perception in flash suppression, equally often monocular and binocular cells. "Perception may be strongly reflected in more elaborate aspects of activity in V1 circuits (e.g., specific neuronal subtypes)... (e.g., synchronization)." — [Keliris, Logothetis & Tolias 2010, J Neurosci 30:12353](https://doi.org/10.1523/jneurosci.0677-10.2010)
  - Rivalry is proposed to arise through interactions between binocular neurons at several levels, not reciprocal inhibition between monocular V1 neurons. — [Leopold & Logothetis 1996, Nature 379:549](https://doi.org/10.1038/379549a0)
  - Perception-related modulation in V1 shows up in LFP power and coherence below 30 Hz, but not in MUA or gamma. — [Gail, Brinksmeyer & Eckhorn 2004, Cereb Cortex 14:300](https://doi.org/10.1093/cercor/bhg129)
  - Rivalry-like counterphase alternation of OD columns occurs under anesthesia ([Xu et al. 2016, J Neurosci 36:3231](https://doi.org/10.1523/jneurosci.3660-15.2016)). In anesthetized 2P imaging, the strength of a neuron's fluctuation correlated with its OD and orientation selectivity ([Wang et al. 2025, Prog Neurobiol 254:102820](https://doi.org/10.1016/j.pneurobio.2025.102820)).
  - Interlaminar coherence tracks interocular conflict and perceptual dominance even when firing rates do not. — [Carlson et al. 2025 bioRxiv](https://doi.org/10.1101/2025.10.04.680461)
- Strabismus: V1 neurons in awake exotropic monkeys "do not appear to reflect strabismic suppression"; elimination of diplopia "is likely to be mediated at a higher cortical level." — [Economides et al. 2021](https://doi.org/10.1523/jneurosci.0044-21.2021). This contrasts with cat strabismic interocular suppression ([Sengpiel et al. 1994](https://doi.org/10.1523/jneurosci.14-11-06855.1994)) and with macaque prism-rearing results ([Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.3.1353); [Zhang et al. 2005](https://doi.org/10.1073/pnas.0500280102)).
- Column-level mechanism: Sengpiel et al. (1995) proposed "intracortical inhibition between adjacent ocular dominance columns" ([doi](https://doi.org/10.1016/0042-6989(94)00125-6)), supported by human 7T fMRI ([Qian et al. 2026](https://doi.org/10.1038/s41562-025-02320-4)). In cat, interocular transfer of adaptation did not depend on position relative to OD-column borders ([Howarth et al. 2009](https://doi.org/10.1093/cercor/bhn211)).
- Excitation vs inhibition in awake V1: inhibition is "monoptic" ([Macknik & Martinez-Conde 2004](https://doi.org/10.1162/0898929041502788)), vs widespread dichoptic suppression ([Cox et al. 2019](https://doi.org/10.1167/19.12.13); [Chen et al. 2026](https://doi.org/10.7554/elife.107518)).

### Inferences
Each open question below gives the competing hypotheses, a discriminating measurement, and how testable it is with the team's dataset: awake macaque, about 30 laminar penetrations, 32-channel 0.1 mm probes, 7 dual-probe days, 6×6 contrast matrix, about 30 trials per condition, about 317 ms stimuli, LFP/CSD/aMUA/~70 single units, eye tracking.

1. **Where is binocular combination first expressed: LGN-inherited, thalamocortical/L4C, or intracortical upper layers?**
   - **Measurement:** in CSD-defined L4C, compare when binocular modulation begins with the onset of the earliest granular sink and MUA, separately for facilitation and suppression and across non-dominant-eye contrast. Modulation present within the first ~5–10 ms of the granular response, and scaling with non-dominant-eye contrast, favors thalamocortical or LGN origin. Modulation lagging by more than ~10 ms and appearing first in L2/3 (as Dougherty 2019's L2/3 47.6 ms vs L4 50.1 ms hints) favors intracortical origin.
   - **Testability:** high for MUA/CSD. Low for monocular single units in L4C because of their rarity. Correct the 3.4 ms NSx/NEV offset first.
2. **Monocular vs interocular vs two-stage contrast gain control.**
   - **Measurement:** fit the full 6×6 matrix per unit or contact with nested models: monocular CGC then summation (Truchard/Meese), interocular normalization (Mitchell eq. 4; Moradi–Heeger), Ding–Sperling gain-control-of-gain-control, and the Georgeson–Sengpiel two-stage model with early IOS and late CGC. Use cross-validated likelihood or AIC, compared by layer and by epoch. Diagnostic patterns (my derivation):
     - Interocular normalization predicts a pure leftward c50 shift of 2^(1/n) for dioptic vs monocular stimulation, with Rmax unchanged.
     - Linear summation of normalized monocular responses predicts an Rmax increase.
     - In the off-diagonal corner (non-dominant eye high, dominant eye low), only interocular normalization predicts suppression below the dominant-eye monocular response.
   - **Testability:** high. This is the dataset's main strength. Lowest contrast is 0.055, so threshold-regime linearity is weakly constrained.
3. **Does the type of gain change differ by layer and epoch (response vs contrast gain; Sengpiel 1998 in cat: L4 response gain, L2/3, 5, 6 contrast gain)?**
   - **Measurement:** layer-specific parameter shifts (ΔRmax vs Δc50 vs Δn) of the dominant-eye CRF as non-dominant-eye contrast increases, in transient and sustained windows.
   - **Testability:** high for MUA. For single units, about 30 trials × 36 conditions supports 3–4-parameter fits but per-layer single-unit counts will be small.
4. **Is OD functionally relevant, or an artifact of subthreshold input (Freeman & Ohzawa 1990; Read & Cumming 2004 vs Mitchell 2023; Zhang 2024)?**
   - **Measurement:** test whether the OD weight w predicts binocular responses out of sample, i.e. with w fit on monocular rows and columns and binocular cells held out. Also test whether this holds equally in L4C, where OD is extreme, and in upper and deep layers. Compare OD from spikes with OD from the laminar CSD sink at the same contact, to address subthreshold input.
   - **Testability:** high. This extends Mitchell 2023 to layers and epochs and is non-redundant.
5. **Dougherty 2019 vs Zhang 2024: under binocular viewing, are eye-preferring neurons suppressed and balanced neurons enhanced, and does this depend on epoch?**
   - **Measurement:** time-resolved binocular modulation index as a function of |ODI| in 40–150 ms vs 150–317 ms windows, using single units and MUA. Zhang's calcium readout integrates up to 1 s, so convergence on the Zhang pattern in the sustained window would reconcile the two.
   - **Testability:** high for MUA; moderate for single units. The 317 ms duration cannot probe beyond about 300 ms.
6. **Is dioptic facilitation additive in the transient and normalized in the sustained epoch, with a contrast-dependent switch (Cox 2019; Mitchell 2022)?**
   - **Measurement:** epoch-specific model comparison on the diagonal and off-diagonal of the matrix, plus the time at which binocular minus monocular changes sign as a function of total contrast.
   - **Testability:** high.
7. **Is interocular suppression delayed relative to excitation, by how much, and does the delay depend on non-dominant-eye contrast or layer?**
   - **Measurement:** onset latency of significant binocular vs monocular divergence, e.g. ROC or sliding-window with 1 ms steps as in Dougherty 2019, per condition and layer. Expected delays are about 7 ms (Tanabe & Cumming 2014) or about 12 ms median (Dougherty 2019).
   - **Testability:** moderate. Use MUA or CSD for SNR and correct clock offsets.
8. **Column-level interocular inhibition (adjacent OD columns; Sengpiel 1995; Qian 2026) vs no column-border dependence (Howarth 2009).**
   - **Measurement:** on the 7 dual-probe days, classify each probe's column by OD from MUA or CSD. Then test whether stimulating probe B's preferred eye suppresses probe A's response to its own eye beyond what A's local OD predicts. Also test whether between-probe LFP coherence or noise correlations fall under dichoptically imbalanced contrasts.
   - **Testability:** partial. With only 7 pairs, results are descriptive or population-pooled across contacts. Probe separation and retinotopic overlap must be checked.
9. **Do granular–infragranular coherence changes track interocular imbalance even without rate changes (Carlson 2025)?**
   - **Measurement:** interlaminar LFP coherence and Granger causality as a function of the (CL, CR) matrix, with balanced vs imbalanced conditions matched for total contrast.
   - **Testability:** high, using LFP across 32 channels.
10. **Monoptic vs dichoptic masking mechanisms (Li 2005; Sengpiel & Vorobyov 2005; Macknik 2004).**
    - **Measurement:** requires orthogonal masks, varied temporal frequency, adaptation, or pharmacology.
    - **Testability:** not testable with an iso-oriented contrast matrix (verify grating orientation and phase across eyes in the task code). The dataset can test only same-orientation interocular normalization.
11. **Interocular phase or vergence confound.** For simple-like cells, binocular modulation of gratings depends on interocular phase ([Ohzawa & Freeman 1986a](https://doi.org/10.1152/jn.1986.56.1.221); [Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.1.351)). Trial-to-trial vergence errors under a mirror stereoscope would therefore add variable facilitation or suppression.
    - **Measurement:** if eye tracking includes both eyes, regress binocular modulation on vergence error. If only one eye is tracked, use the F1/F0 ratio to separate simple from complex units and check whether binocular-modulation variance is larger in simple-like units. Dougherty 2019 reported controlling only for fixational eye movements.
    - **Testability:** depends on whether eye tracking was monocular or binocular (not verified).
12. **Cell-type-specific normalization (PV vs SST in mouse; Perumal 2026; Montgomery 2026).**
    - **Measurement:** compare narrow- vs broad-spiking single units' summation across the matrix. Narrow-spiking (putative PV) units should show sublinear summation that tracks broad-spiking units.
    - **Testability:** low to moderate with about 70 single units.
13. **LGN contribution (Dougherty 2021; Schroeder 1990).** Not directly testable without LGN recordings. A partial test is whether L4C modulation is strongest in the corner where LGN suppression is strongest (low dominant-eye, high non-dominant-eye contrast).

### Gaps
- I did not find explicit discussion-section open questions in Maier et al. 2022 (Annual Review); the full text was not open.
- No study reconciles Macknik & Martinez-Conde 2004 with later dichoptic-suppression findings in awake macaque.
- No primate study tests interocular-suppression mechanisms (temporal frequency, adaptation, GABA) at the single-unit level. Cat results may not transfer, because eye-specific geniculocortical terminations are more segregated in monkeys than in cats ([Dougherty et al. 2019](https://doi.org/10.1016/j.cub.2018.12.004)).
- Some details of the team's task are unverified here and determine items 10 and 11: identical vs differing grating orientation and phase across eyes, the presentation order within trials, and whether vergence was recorded.
- Recent (2024–2026) non-Maier-lab laminar electrophysiology on binocular contrast in macaque seems absent from the indexed literature. The closest work is 2P superficial-layer imaging (Tang/Yu lab: Zhang 2024; Chen 2026; Wang et al. 2026 [Prog Neurobiol 258:102874](https://doi.org/10.1016/j.pneurobio.2025.102874) on OD population geometry) and human 7T laminar fMRI. A Neuropixels laminar study of monkey V1 ([Carr et al. 2026, PNAS](https://doi.org/10.1073/pnas.2521556123)) was found by title only; whether it includes binocular conditions was not checked.
