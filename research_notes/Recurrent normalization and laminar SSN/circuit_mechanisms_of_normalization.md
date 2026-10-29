# Circuit mechanisms of normalization in V1: feedforward vs recurrent vs feedback, timing, layers, and interocular suppression

*Research notes compiled 2026-10-29. Each finding is tagged with species and preparation (AW = awake, AN = anesthetized). The team's results are referred to as F1–F4: (F1) the two eyes combine near-linearly in the first ~20 ms, including in the earliest current sink; (F2) interocular normalization arrives ~20–35 ms after onset, in CSD and spiking together, and is faster and stronger at higher contrast; (F3) the normalization pool weights the non-dominant eye (NDE) more than excitation does; (F4) timing and strength are similar across depth. "Prep. unverified" means the preparation is the lab's usual one but was not stated in the abstract or text I retrieved.*

---

## Q1. Evidence for feedforward vs recurrent (vs feedback) mechanisms of contrast normalization, cross-orientation suppression, surround suppression, and interocular suppression

### Takeaway
Monocular cross-orientation suppression (COS) in V1 is largely feedforward. It comes from LGN contrast saturation and rectification, amplified by spike threshold, with no added cortical inhibition. Dichoptic (interocular) suppression has the opposite signatures in every test so far: a cortical temporal-frequency cutoff, it is abolished by adaptation, and it is GABA-A dependent. That points to an intracortical, recurrent or lateral origin. Surround suppression has at least three parts: an LGN-inherited fast component, an intracortical component (horizontal connections and SOM cells in mouse, consistent with an inhibition-stabilized network), and a feedback component. In cortex, normalization often works by *withdrawing* recurrent excitation (and inhibition) rather than by adding shunting inhibition.

### Cited Findings

#### Monocular cross-orientation suppression: feedforward/thalamic
- **Priebe & Ferster 2006** (cat, simple cells, in vivo intracellular; prep. unverified, AN): cross-oriented stimuli *suppressed rather than enhanced* synaptic inhibition while also suppressing excitation. LGN relay-cell contrast saturation and rectification, fed into a linear feedforward model, predicted COS of the F1 component of excitation. Spike suppression was about twice the suppression of synaptic input, so threshold amplifies it. Priebe NJ, Ferster D. "Mechanisms underlying cross-orientation suppression in cat visual cortex." *Nat Neurosci* 9:552–561. **Paywalled; no PMC copy.** — [doi:10.1038/nn1660](https://doi.org/10.1038/nn1660)
- **Freeman, Durand, Kiper & Carandini 2002** (cat, AN): masks drifting too fast to drive cortex still suppress. The suppression high-cut was 19.4 ± 3.1 Hz (median 19.6 Hz, n = 44), versus 7.5 ± 1.0 Hz (median 7.7 Hz) for V1 responses. Adapting to the mask barely changed suppression: maximal response fell by a factor of 1.1 ± 1.1 and c50 rose by 1.5 ± 1.1. They proposed a thalamic origin plus thalamocortical synaptic depression. The paper explicitly says this account does *not* cover dichoptic COS, "as its effects can take hundreds of milliseconds to reach peak strength" (attributed to Sengpiel, personal communication), and that "dichoptic suppression must be due to intracortical effects, perhaps to lateral inhibition." "Suppression without inhibition in visual cortex." *Neuron* 35:759–771. — [doi:10.1016/S0896-6273(02)00819-X](https://doi.org/10.1016/S0896-6273(02)00819-X); PDF at [Kiper lab](https://tilde.ini.uzh.ch/~kiper/my_pdfs/freeman.pdf)
- **Carandini, Heeger & Senn 2002** (model): thalamocortical depression predicts that suppression is monocular, immune to cortical adaptation, and evoked by masks that drive thalamus but not cortex. It also predicts response saturation with contrast-invariant tuning. "A synaptic explanation of suppression in visual cortex." *J Neurosci* 22:10053–10065. — [doi:10.1523/JNEUROSCI.22-22-10053.2002](https://doi.org/10.1523/JNEUROSCI.22-22-10053.2002)
- **Li, Thompson, Duong, Peterson & Freeman 2006** (cat, prep. unverified, AN): COS recovers "too rapid[ly] to be explained by synaptic depression." LGN contrast saturation, combined with LGN rectification and cortical spike threshold, "can account for most of the COS." "Origins of cross-orientation suppression in the visual cortex." *J Neurophysiol* 96:1755–1764. — [doi:10.1152/jn.00425.2006](https://doi.org/10.1152/jn.00425.2006)
- **Smith, Bair & Movshon 2006** (macaque, AN, 112 complex cells): recovery from suppression came within ≤20 ms of response onset, whereas synaptic depression recovers in 60–600 ms in vitro, so depression is unlikely. The data fit Priebe & Ferster's LGN-nonlinearity account but cannot rule out fast local inhibition. Extrastriate feedback was judged implausible for COS. "Dynamics of suppression in macaque primary visual cortex." *J Neurosci* 26:4826–4834. — [PMC6674161](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161)
- **Bonin, Mante & Carandini 2005** (cat LGN, prep. unverified): LGN neurons already show contrast saturation, size tuning and masking. One "suppressive field" model, inherited from retina, explains all three. That field is untuned for orientation and responds to high temporal frequencies, which is "hardly consistent with feedback from cortex." *J Neurosci* 25:10844–10856. — [doi:10.1523/JNEUROSCI.3562-05.2005](https://doi.org/10.1523/JNEUROSCI.3562-05.2005)

#### Interocular (dichoptic) suppression: intracortical
- **Li, Peterson, Thompson, Duong & Freeman 2005** (cat, AN): dichoptic COS is strongest at *low* temporal frequencies and is markedly weakened by adaptation to the mask. Monoptic COS peaks at much higher temporal frequencies and resists adaptation. So monoptic COS is subcortical and dichoptic COS is intracortical. "Cross-orientation suppression: monoptic and dichoptic mechanisms are different." *J Neurophysiol* 94:1645–1650. — [doi:10.1152/jn.00203.2005](https://doi.org/10.1152/jn.00203.2005)
- **Sengpiel & Vorobyov 2005** (10 young cats, AN and paralyzed; 74 cells for interocular, 73 for monocular COS):
  - Temporal frequency: interocular suppression peaked at 4 Hz (28.6%), was not significant at 16 Hz, and was absent at 32 Hz. Monocular COS was maximal at 8 Hz (44.6%) and still significant at 32 Hz (8.9%).
  - Adaptation to the orthogonal grating cut interocular suppression from 33.4 ± 2.9% to 5.6 ± 8.7% (24/27 cells).
  - Bicuculline (GABA-A antagonist) cut it from 46.1 ± 4.0% to 18.1 ± 4.9% (25/28 cells).
  - No latency data.
  - "Intracortical origins of interocular suppression in the visual cortex." *J Neurosci* 25:6394–6400. — [doi:10.1523/JNEUROSCI.0862-05.2005](https://doi.org/10.1523/JNEUROSCI.0862-05.2005); [PMC6725266](https://pmc.ncbi.nlm.nih.gov/articles/PMC6725266)
- **Sengpiel, Blakemore & Harrad 1995** (cat LGN and area 17): almost half of LGN neurons showed modest interocular inhibition that did not depend on orientation, and monocular layer-4 units behaved similarly. In most binocular cells, an orthogonal grating in the other eye suppressed the response to the optimal grating, independent of phase and over a wide SF range. Suppression appeared "only if the neuron is already responding" to the other eye. They proposed intracortical inhibition between adjacent ocular dominance columns. *Vision Res* 35:179–195. **Paywalled.** — [doi:10.1016/0042-6989(94)00125-6](https://doi.org/10.1016/0042-6989(94)00125-6)
- **Sengpiel & Blakemore 1994** (cat): the sudden appearance of an inappropriate stimulus in one eye suppresses neurons already responding to the other eye. In strabismic cats even similar orientations trigger this. *Nature* 368:847–850. **Paywalled.** — [doi:10.1038/368847a0](https://doi.org/10.1038/368847a0)
- **Walker, Ohzawa & Freeman 1998** (cat, prep. unverified): binocular COS was as strong as monocular COS and unaffected by disparity. Interocular COS was weak but present in almost every cell. Temporal offsets between test and mask did not change COS. They concluded that COS is "primarily monocular and acts prior to the convergence." *J Neurophysiol* 79:227–239. — [doi:10.1152/jn.1998.79.1.227](https://doi.org/10.1152/jn.1998.79.1.227)
- **Truchard, Ohzawa & Freeman 2000** (cat binocular simple cells, AN, n = 23): raising contrast from 5 to 50% reduced monocular gain about 5-fold (relative gain 5.0 ± 1.9). Binocular gain changed only 1.2 ± 0.3-fold. A 50% grating in the other eye did not consistently suppress (relative gain 1.2 ± 0.53, n.s.). Conclusion: gain reductions occur "primarily at a monocular site, before convergence." These used drifting gratings averaged over 4 s. "Contrast gain control in the visual cortex: monocular versus binocular mechanisms." *J Neurosci* 20:3017–3032. — [PMC6772209](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209)
- **Sengpiel, Baddeley, Freeman, Harrad & Blakemore 1998** (cat): COS acted mainly as contrast gain (rightward CRF shift). Surround and interocular suppression acted mainly as response gain (downward shift), with some contrast-gain component. *Vision Res* 38:2067–2080. **Paywalled.** — [doi:10.1016/S0042-6989(97)00413-6](https://doi.org/10.1016/S0042-6989(97)00413-6)
- **Macknik & Martinez-Conde 2004** (macaque, AW, dichoptic masking), a dissenting result: V1 binocular neurons are "dichoptic in terms of excitation, but monoptic in terms of inhibition." Targets were strongly inhibited only by masks in the same eye. *J Cogn Neurosci* 16:1049–1059. **Paywalled.** — [doi:10.1162/0898929041502788](https://doi.org/10.1162/0898929041502788)
- **Dougherty et al. 2021** (macaque LGN, AW, 56 units: 20 P, 30 M, 6 K): across all conditions the median binocular effect was −2.65% (n.s.). With dominant-eye (DE) contrast at 0.12 and NDE contrast at 0.9, the median was −9.95% (p = 0.0024), and at most ~25% of units were suppressed. There was no significant modulation at medium or high DE contrast. Onset latency of the modulation could not be analyzed. "Binocular suppression in the macaque lateral geniculate nucleus reveals early competitive interactions between the eyes." *eNeuro* 8(2). — [doi:10.1523/ENEURO.0364-20.2020](https://doi.org/10.1523/ENEURO.0364-20.2020); [PMC8035044](https://pmc.ncbi.nlm.nih.gov/articles/PMC8035044/)
- **Dougherty, Schmid & Maier 2019** (review): binocular modulation in LGN is predominantly suppressive. It could be intrinsic (local LGN or TRN circuits) or arrive indirectly from binocular structures such as V1 feedback. *J Comp Neurol* 527:522–534. — [doi:10.1002/cne.24417](https://doi.org/10.1002/cne.24417)

#### Surround suppression and other normalization: feedforward plus recurrent plus feedback
- **Alitto & Usrey 2008** (macaque LGN, AW and AN): extraclassical suppression is stronger in M than P cells. It "arises too quickly for involvement from cortical feedback," and retinal suppression matches LGN suppression. So cortical surround suppression likely includes a retinal component. *Neuron* 57:135–146. — [doi:10.1016/j.neuron.2007.11.019](https://doi.org/10.1016/j.neuron.2007.11.019)
- **Sceniak, Chatterjee & Callaway 2006** (macaque, AN; geniculocortical afferents recorded with cortex silenced by muscimol): more than 80% of afferents in 4C and 4A had strong surround suppression (SI > 0.5). *J Neurophysiol* 96:3474–3484. — [doi:10.1152/jn.00734.2006](https://doi.org/10.1152/jn.00734.2006)
- **Webb, Dhruv, Solomon, Tailby & Lennie 2005** (macaque V1; prep. unverified, AN): there are two surround mechanisms.
  - One is prominent with high-contrast CRF drive. It is orientation selective, has relatively sharp spatiotemporal tuning, is *binocularly driven* (substantial interocular transfer), and is strongly desensitized by adaptation.
  - The other is more prominent at low CRF contrast. It is very broadband, monocular and adaptation-immune, and is attributed to the "input layers of primary visual cortex, or earlier."
  - *J Neurosci* 25:11666–11675. — [doi:10.1523/JNEUROSCI.3414-05.2005](https://doi.org/10.1523/JNEUROSCI.3414-05.2005); [PDF](https://www.cnbc.cmu.edu/braingroup/papers/webb_etal_2005.pdf)
- **Ozeki, Finn, Schaffer, Miller & Ferster 2009** (cat, AN, intracellular): iso-oriented surrounds cut excitatory conductance by 51% and inhibitory conductance by 38% on average. In 6 cells, surround onset caused a transient rise in inhibition lasting 30–50 ms before both conductances settled lower. LGN surround suppression (~20% at iso-orientation) was too weak and too weakly tuned to explain cortical suppression. Conclusion: V1 is an inhibition-stabilized network (ISN). *Neuron* 62:578–592. — [PMC2691725](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725)
- **Haider et al. 2010** (cat, prep. unverified), a partial contrast to Ozeki: combined CRF + nCRF stimulation increased inhibitory barrages onto regular-spiking pyramidal cells, and fast-spiking interneurons fired more. *Neuron* 65:107–121. — [doi:10.1016/j.neuron.2009.12.005](https://doi.org/10.1016/j.neuron.2009.12.005)
- **Sato, Haider, Häusser & Carandini 2016** (mouse V1 L2/3 pyramidal cells, lightly isoflurane-anesthetized, whole-cell; normalization driven by ChR2 activation of callosal input from a site ~60° away):
  - At 100% contrast, distal activation *reduced* both IPSC conductance (−0.91 ± 0.25 nS) and EPSC conductance (−0.34 ± 0.08 nS), measured 150–300 ms after onset.
  - At 0% contrast it added both.
  - No transient IPSC increase was found (80–120 ms window, 0/10 cells).
  - Its effect was mainly additive from 0–150 ms and divisive from 150–300 ms.
  - Normalization "operate[s] largely by modulating the level of excitation."
  - *Nat Neurosci* 19:568–570. — [doi:10.1038/nn.4249](https://doi.org/10.1038/nn.4249); [PMC4817833](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817833)
- **Adesnik, Bruns, Taniguchi, Huang & Scanziani 2012** (mouse, prep. unverified): superficial-layer SOM cells *increase* their firing with surround stimulation because horizontal axons preferentially excite them. Perturbing SOM cells altered pyramidal surround suppression. *Nature* 490:226–231. — [doi:10.1038/nature11526](https://doi.org/10.1038/nature11526)
- **Adesnik 2017** (mouse L2/3, AW, whole-cell): the E/I ratio falls systematically with increasing contrast or size. Suppressing SOM cells enhanced the E and I that underlie size tuning. Presented as support for the stabilized supralinear network (SSN). *Neuron* 95:1147–1159. — [doi:10.1016/j.neuron.2017.08.014](https://doi.org/10.1016/j.neuron.2017.08.014)
- **Wilson, Runyan, Wang & Sur 2012** (mouse, prep. unverified): PV cells mainly *divide* responses and preserve selectivity; SOM cells mainly *subtract* and sharpen selectivity. *Nature* 488:343–348. — [doi:10.1038/nature11347](https://doi.org/10.1038/nature11347)
- **Atallah, Bruns, Carandini & Scanziani 2012** (mouse, prep. unverified), which conflicts with a purely divisive PV role: PV effects on L2/3 pyramidal cells were "a linear transformation, both additive and multiplicative, with a threshold." *Neuron* 73:159–170. — [doi:10.1016/j.neuron.2011.12.013](https://doi.org/10.1016/j.neuron.2011.12.013)
- **Rockwell & MacLean 2025** (mouse, calcium imaging of E and PV cells; preprint, not peer reviewed): suppression in PV cells was *transiently reduced* by the dynamics of subcortical input. Heterogeneous normalization strength was linked to suppression of both subcortical and local cortical inputs. — [doi:10.1101/2025.08.19.671076](https://doi.org/10.1101/2025.08.19.671076)
- **Nassi, Avery, Cetin, Roe & Reynolds 2015** (macaque, AW; C1V1 under CaMKIIα in excitatory cells, 249 units):
  - Light and visual responses summed sub-additively (median additive index 0.67 ± 0.02, n = 93).
  - A normalization model explained 82% of variance (C1V1) and 86% (eArch3.0).
  - Light facilitated responses at low contrast and suppressed them at high contrast.
  - Direct activation was likely supragranular because light attenuates with depth; there were no laminar recordings.
  - Effects were judged "likely, though not definitively," local to V1.
  - *Neuron* 86:1504–1517. — [doi:10.1016/j.neuron.2015.05.040](https://doi.org/10.1016/j.neuron.2015.05.040); [PMC4534089](https://pmc.ncbi.nlm.nih.gov/articles/PMC4534089)
- **Bair, Cavanaugh & Movshon 2003** (macaque, AN): surround suppression propagated across 6–8 mm of cortex at ~1 m/s (>1 m/s in ~40% of cells), faster than horizontal axons (~0.1–0.2 m/s). Feedback was favored, plus a "second source of suppression that is local, strong, and fast." *J Neurosci* 23:7690–7701. — [PMC6740744](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744)
- **Nurminen, Merlin, Bijanzadeh, Federer & Angelucci 2018** (marmoset; ArchT silencing of V2→V1 feedback terminals; prep. not verified in my sources): reducing feedback enlarged V1 RFs, reduced surround suppression and responses to RF-confined stimuli, and increased responses to stimuli reaching into the proximal surround. Stronger reduction lowered response gain further. *Nat Commun* 9:2281. — [doi:10.1038/s41467-018-04500-5](https://doi.org/10.1038/s41467-018-04500-5)
- **Vangeneugden et al. 2019** (mouse, AW vs AN): silencing lateral visual areas (LM) reduced surround suppression in awake but not anesthetized mice. *Curr Biol* 29:4268–4275. **Paywalled.** — [doi:10.1016/j.cub.2019.10.037](https://doi.org/10.1016/j.cub.2019.10.037)
- **Reynaud, Masson & Chavane 2012** (macaque, AW, VSDI): V1 normalization is "a timely push–pull between local and lateral networks." Local recurrent loops increase contrast gain over time, while centripetal lateral input (0.1–0.3 m/s, in the horizontal-axon range) clamps it. *J Neurosci* 32:12558–12569. — [PMC6621242](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242)
- **Lien & Scanziani 2013** (mouse L4, optogenetic cortical silencing; prep. unverified): thalamic excitation is about one third of total excitation onto L4 neurons. "A principal role of cortical circuits is to amplify tuned thalamic excitation." *Nat Neurosci* 16:1315–1323. — [doi:10.1038/nn.3488](https://doi.org/10.1038/nn.3488)

#### Theory
- **Rubin, Van Hooser & Miller 2015 (SSN)**: surround suppression and sublinear normalization emerge from supralinear input/output functions, strong recurrent excitation and feedback inhibition. For weak drive, summation becomes linear or supralinear and contextual effects facilitate. E and I cells behave similarly. *Neuron* 85:402–417. — [doi:10.1016/j.neuron.2014.12.026](https://doi.org/10.1016/j.neuron.2014.12.026)
- **Shushruth et al. 2012** (macaque, prep. unverified, plus model): the orientation tuning of surround modulation arises "only if V1 operates in a regime of strong local recurrence," through "maximal withdrawal of recurrent excitation." *J Neurosci* 32:308–321. — [doi:10.1523/JNEUROSCI.3789-11.2012](https://doi.org/10.1523/JNEUROSCI.3789-11.2012)
- **Sanzeni et al. 2020** (mouse V1/S1/M1, AW and light AN): paradoxical ISN effects appear even without sensory stimulation. *eLife* 9:e54875. — [doi:10.7554/eLife.54875](https://doi.org/10.7554/eLife.54875)
- **Heeger & Zemlianova 2020 (ORGaNICs)**: a family of recurrent circuits implements weighted normalization through recurrent amplification and reproduces V1 dynamics. *PNAS* 117:22494–22505. — [doi:10.1073/pnas.2005417117](https://doi.org/10.1073/pnas.2005417117)
- **Carandini & Heeger 2012**: review of normalization as a canonical computation. *Nat Rev Neurosci* 13:51–62. — [doi:10.1038/nrn3136](https://doi.org/10.1038/nrn3136)

### Inferences
- **Two routes to division.** The literature consistently separates a fast, monocular, broadband, adaptation-immune route (precortical: LGN saturation and rectification plus threshold) from a slower, eye-crossing, adaptable, GABA-dependent route (cortical). Interocular normalization belongs to the second class in every test that separated them: Li 2005, Sengpiel & Vorobyov 2005, Webb 2005, and Baker/Meese psychophysics (Q4).
- **Bearing on F1.** A feedforward, monocular-stage gain control (Truchard 2000; Priebe & Ferster 2006) would already be present in each eye's signal at onset. It would not create *interocular* sublinearity, so near-linear binocular summation at onset is expected even if each eye's signal is saturated. LGN interocular suppression is small and appears only when DE contrast is low (Dougherty 2021). That fits the absence of interocular normalization in the earliest L4C sink, at least for medium and high DE contrasts. One caveat: at low DE contrast with high NDE contrast, an LGN-inherited suppression of up to ~10% could be present at onset.
- **Bearing on F2.** Delayed interocular suppression appearing in synaptic currents (CSD) fits a cortical (lateral or recurrent) mechanism. It does not fit thalamocortical depression, which predicts suppression from the first spikes and immunity to cortical state. The cat evidence favors an intracortical GABA-A pathway, while Ozeki 2009 and Sato 2016 show cortical normalization can act through *withdrawal* of recurrent excitation. These two are not mutually exclusive in an ISN: inhibition drives the change, but the steady state shows lower E and I.
- **Bearing on F3.** Truchard 2000 found weak binocular gain control in cat; the team's data favor a substantial interocular term. The difference may reflect species (cat vs macaque), stimulus (drifting gratings averaged over 4 s vs static flashed gratings) or analysis (F1 simple-cell responses vs MUA/CSD).
- **Bearing on F4.** Bijanzadeh 2018 (Q3) shows LGN-inherited suppression appearing first in L4C and intracortical suppression appearing last in L4C. A depth-invariant interocular normalization argues against both a purely precortical origin and a strictly superficial-first horizontal origin. It fits a distributed recurrent mechanism or one that spreads interlaminarly faster than the measurement's resolution.

### Gaps
- No study I found measured *intracellular* E and I conductances during interocular suppression in any species. The mouse binocular whole-cell studies (Longordo 2013; Zhao 2013; Q4) measure sublinear binocular summation, not dichoptic orthogonal suppression time courses.
- Priebe & Ferster 2006 (abstract only; full text paywalled) gives no latency numbers in the abstract.
- The Freeman 2002 claim that dichoptic suppression takes "hundreds of milliseconds to reach peak" rests on a personal communication, not published latencies.
- The preparation for Nurminen 2018, Adesnik 2012, Wilson 2012, Atallah 2012 and Lien & Scanziani 2013 was not confirmed in the text I retrieved.
- PDFs the team may need to supply: Priebe & Ferster 2006 (Nat Neurosci); Sengpiel et al. 1995 and 1998 (Vision Res); Sengpiel & Blakemore 1994 (Nature); Macknik & Martinez-Conde 2004 (J Cogn Neurosci); Vangeneugden 2019 (Curr Biol).

---

## Q2. Timing: latency of normalization or suppression relative to excitation, its contrast dependence, and conductance time courses

### Takeaway
Latencies differ by suppression type:

| Suppression type | Delay after excitatory onset | Source |
|---|---|---|
| Monocular COS | ≈0 ms (suppression latency ~6–7.5 ms *shorter* than excitatory onset latency: paired mean 6.4 ms, difference of means 7.5 ms) | Smith 2006 |
| Near-surround suppression | ~5–10 ms on average | Bair 2003; Smith 2006 |
| Awake laminar "fast suppression" (flashed gratings) | ~11–21 ms | Wang 2020 |
| Binocular suppression of monocular neurons (awake macaque) | ~12 ms between pooled medians | Dougherty 2019 |
| Disparity-tuned suppression | ~7 ms | Tanabe & Cumming 2014 |
| Population contrast gain (VSDI) | develops with τ ≈ 31 ms | Reynaud 2012 |
| "Slow suppression" | ~45–61 ms | Wang 2020 |
| Feedback-like signals | ≈30+ ms later still | Bijanzadeh 2018; Self 2013 |

The team's 20–35 ms falls between the feedforward/near-surround group and the slow/feedback group. It best matches local recurrent or lateral intracortical estimates. Stronger suppression tends to arrive earlier (Bair 2003; Smith 2006 for surround). Conductance studies show that visually evoked conductance increases are orientation-tuned rather than untuned shunting (Anderson 2000). Cortical normalization can show up as a delayed fall in both E and I (Ozeki 2009; Sato 2016).

### Cited Findings

#### Monocular CRF suppression and normalization latencies
- **Smith, Bair & Movshon 2006** (macaque, AN; continuous random sequences where the mask switched on and off while the optimal grating stayed on):
  - Latencies: response onset 50.0 ± 15.8 ms; response offset 30.1 ± 9.3 ms; COS 42.5 ± 11.1 ms; COS release 40.9 ± 12.4 ms. So COS began ~6–7 ms *before* the onset latency of excitation and 13.5 ms after the offset latency.
  - Within cells, surround suppression and release were ~12 ms later than COS.
  - Grand-PSTH surround suppression "did not stabilize until ∼50 ms after the onset response began."
  - Mean suppression was 25% for COS and 42% for surround.
  - Surround suppression strength was anticorrelated with latency (r = −0.48, n = 87). For COS the correlation was weaker and n.s. (r = −0.25).
  - Contrast was fixed at 50% per component.
  - — [PMC6674161](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161)
- **Bair, Cavanaugh & Movshon 2003** (macaque, AN, 92 cells): surround suppression was delayed by 9 ms on average (SD 15) relative to CRF onset (~52 vs ~61 ms), and some cells were suppressed before CRF onset. "Strong suppression arrived on average 30 msec earlier than weak suppression." — [PMC6740744](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744)
- **Albrecht, Geisler, Frazor & Crane 2002** (monkey and cat V1, stationary gratings flashed for 200 ms; prep. unverified, AN):
  - Contrast-set gain control and response expansion are "fully expressed within the first 10 ms after the onset of the response."
  - One scaled Naka–Rushton function, with the same c50 and exponent, fits CRFs from the first through the last 10 ms.
  - Scaling and shifting one PSTH explains ~95% of the contrast-dependent variation.
  - *J Neurophysiol* 88:888–913. — [doi:10.1152/jn.2002.88.2.888](https://doi.org/10.1152/jn.2002.88.2.888)
- **Reynaud, Masson & Chavane 2012** (macaque, AW, VSDI; 2° drifting grating, 2.5–80% contrast), which **conflicts with Albrecht 2002 at the population level**:
  - c50 fell from ~30% to ~5% between 40 and 100 ms, with a time constant ~31 ms, so contrast gain increases over time.
  - Response latency shortened from 87 to 40 ms as contrast rose.
  - Center-driven activity arrived in the center ROI at ~33 ms. Surround-driven activity arrived at 51, 62 and 92 ms for proximal, intermediate and distal surrounds, and clamped the center CRF at ~64, 82 and 127 ms.
  - — [PMC6621242](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242)
- **Sit, Chen, Geisler, Miikkulainen & Seidemann 2009** (macaque, AW, VSDI, briefly flashed Gabor): the complex V1 population dynamics are explained by a *feedforward* population gain-control model in which population activity sets local dynamics and gain. *Neuron* 64:943–956. — [doi:10.1016/j.neuron.2009.08.041](https://doi.org/10.1016/j.neuron.2009.08.041)
- **Zhou, Benson, Kay & Winawer 2019** (human ECoG/fMRI and macaque V1 spikes; model): a gain-control signal that is *slower than* the linear drive is critical for sub-additive temporal summation and for slower responses at low contrast. Later visual areas show faster and stronger gain control. *PLoS Comput Biol* 15:e1007484. — [doi:10.1371/journal.pcbi.1007484](https://doi.org/10.1371/journal.pcbi.1007484)
- **Wang et al. 2020** (macaque, AW; 24-channel laminar U-probe; flashed gratings; MUA split into excitation, fast suppression (FS) and slow suppression (SS)):

  | Layer | FS latency after excitation (ms) | SS latency after excitation (ms) |
  |---|---|---|
  | L2/3 | 20.5 | 44.7 |
  | L4B | 18.1 | 46.3 |
  | L4Cα | 11.4 | 60.3 |
  | L4Cβ | 21.3 | 61.0 |
  | L5 | 18.3 | 46.8 |
  | L6 | 14.9 | 58.8 |

  The authors call FS "likely … a local process" (feedforward or recurrent) and attribute SS to feedback. They note that the "10 ms" FS figure in their own Discussion only roughly matches the 11–21 ms in Table 1. — [doi:10.1523/JNEUROSCI.1129-20.2020](https://doi.org/10.1523/JNEUROSCI.1129-20.2020); [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/)
- **Henry, Jazayeri, Shapley & Hawken 2020** (macaque, AN, 106 neurons, 81 laminar-assigned): eCRF modulation *peak* times were untuned suppression 55 ± 2 ms (duration 23 ms), facilitation 60 ± 3 ms (16 ms), and tuned suppression 76 ± 3 ms (36 ms). Tuned suppression peaked at 60–75 ms in 4Cα and ~100 ms in 4Cβ and L2/3. These are peaks, not onsets. *eLife* 9:e54264. — [PMC7253173](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173)
- **Ringach, Hawken & Shapley 2003; Shapley, Hawken & Ringach 2003** (macaque; prep. unverified, AN): orientation dynamics show "early global enhancement followed by global and tuned suppression." Global suppression is plausibly cortical inhibition. *J Neurophysiol* 90:342–352 — [doi:10.1152/jn.01018.2002](https://doi.org/10.1152/jn.01018.2002); *Neuron* 38:689–699 — [doi:10.1016/S0896-6273(03)00332-5](https://doi.org/10.1016/S0896-6273(03)00332-5). **No ms values in the abstracts.**

#### Binocular and interocular timing
- **Dougherty, Cox, Westerberg & Maier 2019** (macaque, AW; 290 neurons, 64 monocular, 33 tested binocularly: 23 suppressed, 10 facilitated; contrast ≥0.8):
  - Visual response onsets: L4 42.1 ms, L2/3 45.2 ms, L5/6 46.4 ms.
  - Median response onset was 50 ms and median binocular-modulation onset 62 ms (p = 2 × 10⁻¹⁷). The ~12 ms gap is computed from pooled medians, not a figure the authors report.
  - Modulation onset by layer: L2/3 47.6, L4 50.1, L5/6 54.7 ms.
  - In L4, the 0.25-quantile onset was 49 ms for facilitated and 54 ms for suppressed neurons (n.s., p = 0.16).
  - Abstract: facilitation in input layers "tended to occur at the onset of the visual response" (explainable by converging thalamocortical inputs), whereas suppression "occurred several milliseconds following the onset … suggesting that the bulk of binocular modulation involves cortical inhibition."
  - Stronger eye preference predicted stronger suppression (r = −0.275 for neurons, −0.348 for MUA).
  - *Curr Biol* 29:381–391. — [doi:10.1016/j.cub.2018.12.004](https://doi.org/10.1016/j.cub.2018.12.004); [PMC6363852](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852)
- **Cox, Dougherty, Westerberg, Schall & Maier 2019** (macaque, AW; 152 units; DE contrast typically 0.45, NDE 0.9):
  - Transient window (50–100 ms): binocular > monocular for both congruent (+13.73 sp/s) and incongruent (+5.85 sp/s) stimuli, and binocular stimulation enlarged the initial granular CSD sink.
  - Sustained window (150–250 ms): congruent facilitation disappeared and incongruent stimuli were suppressed (−5.57 sp/s). Ocular-biased neurons showed transient facilitation and then suppression (−9.50 sp/s incongruent).
  - The transient is interpreted as feedforward and the sustained suppression as cortico-cortical, possibly horizontal links between neighbouring OD columns.
  - The coarse windows limit latency resolution.
  - *J Vis* 19(12):13. — [doi:10.1167/19.12.13](https://doi.org/10.1167/19.12.13); [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- **Mitchell, Dougherty, Westerberg, Carlson, Daumail, Maier & Cox 2022** (macaque, AW; dioptic stimuli): binocular facilitation "scaled most proportionately with contrast during the initial transient." Later came "longer maintenance of facilitation for lower contrast and binocular suppression at high contrast." *iScience* 25:104182. — [doi:10.1016/j.isci.2022.104182](https://doi.org/10.1016/j.isci.2022.104182)
- **Tanabe & Cumming 2014** (macaque, AW; binocular RF subspace mapping): disparity-tuned suppression lags excitation by a median 7 ms, consistently across the population, so false-match responses are progressively attenuated. *J Neurophysiol* 111:1759–1769. — [doi:10.1152/jn.00426.2013](https://doi.org/10.1152/jn.00426.2013)
- **Menz & Freeman 2003** (cat; prep. unverified): disparity tuning generally *sharpened* over the response time course (coarse-to-fine). *Nat Neurosci* 6:59–65. **Paywalled.** — [doi:10.1038/nn986](https://doi.org/10.1038/nn986)
- **Baker, Meese & Summers 2007** (human psychophysics): the weight of interocular cross-orientation suppression *increases* with stimulus duration, whereas the ipsiocular weight *decreases*. The ipsiocular route is attributed to LGN or retina and the interocular route to cortex. *Neuroscience* 146:435–448. — [doi:10.1016/j.neuroscience.2007.01.030](https://doi.org/10.1016/j.neuroscience.2007.01.030); [PDF](https://eprints.soton.ac.uk/47919/1/Baker_et_al_2007.pdf)
- **Walker, Ohzawa & Freeman 1998** (cat): asynchronous test and mask presentation did not alter COS strength. — [doi:10.1152/jn.1998.79.1.227](https://doi.org/10.1152/jn.1998.79.1.227)
- **Moradi & Heeger 2009** (human fMRI): a 500 ms interocular onset asynchrony had no effect on dichoptic-plaid suppression when attention was diverted. — [doi:10.1167/9.3.13](https://doi.org/10.1167/9.3.13)

#### Conductance time courses and shunting
- **Anderson, Carandini & Ferster 2000** (cat, intracellular; prep. unverified, AN):
  - Input conductance rose with contrast by 20–300% (simple cells 40–300%, complex 20–120%) and was "invariably maximal" at the preferred orientation. This contradicts both the untuned-shunting normalization prediction and the veto model.
  - In simple cells, E and I were often push-pull.
  - E and I had similar preferred orientations and tuning widths.
  - *J Neurophysiol* 84:909–926. — [doi:10.1152/jn.2000.84.2.909](https://doi.org/10.1152/jn.2000.84.2.909)
- **Borg-Graham, Monier & Frégnac 1998** (cat, in vivo whole-cell voltage clamp; prep. unverified): optimally oriented flashed bars raised somatic input conductance to more than 3× rest. The peak came at short latency and looked GABA-A dominated ("transient and strong shunting inhibition"). *Nature* 393:369–373. **Paywalled.** — [doi:10.1038/30735](https://doi.org/10.1038/30735)
- **Carandini & Heeger 1994**: the original shunting-inhibition normalization model. *Science* 264:1333–1336. — [doi:10.1126/science.8191289](https://doi.org/10.1126/science.8191289)
- **Ozeki et al. 2009** (cat, AN): at surround onset there was a 30–50 ms transient rise in inhibition, after which steady-state E fell 51% and I fell 38%. The paper cites cortical surround suppression as "delayed by at least 10 ms" after center onset. — [PMC2691725](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725)
- **Sato et al. 2016** (mouse, light AN): no transient IPSC increase at 80–120 ms; the divisive effect dominated at 150–300 ms. — [PMC4817833](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817833)
- **Haider, Häusser & Carandini 2013** (mouse, AW vs AN, whole-cell): under anesthesia inhibition tracked excitation. Awake, "inhibition was much stronger than excitation and had extremely broad spatial selectivity," and responses were briefer. *Nature* 493:97–100. — [doi:10.1038/nature11665](https://doi.org/10.1038/nature11665)
- **Nassi et al. 2015** (macaque, AW): single-unit latency to light was 10.4 ± 1.9 ms (17/55 units ≤2 ms). Pulsed 5 ms stimulation also produced suppression (post-stimulation window 10–210 ms after offset). — [PMC4534089](https://pmc.ncbi.nlm.nih.gov/articles/PMC4534089)

#### Reference latencies for laminar interpretation
- **Maunsell & Gibson 1992** (macaque; prep. unverified):
  - Most V1 latencies were 30–50 ms; the shortest per animal were 20–31 ms.
  - Magnocellular LGN lesions delayed the earliest V1 responses by 7–10 ms, so "the first 7–10 ms of cortical activation" are magnocellular.
  - The shortest-latency units were in 4C or 6.
  - "The most superficial neurons [were] not becoming active until 15 ms after responses were observed in layer 4C."
  - *J Neurophysiol* 68:1332–1344. — [doi:10.1152/jn.1992.68.4.1332](https://doi.org/10.1152/jn.1992.68.4.1332)
- **Self, van Kerkoerle, Supèr & Roelfsema 2013** (macaque, AW, figure-ground task): responses started at 40 ms in L4 and L6. Boundary signals (attributed to horizontal connections) appeared at ~70 ms in upper L4 and superficial layers. Feedback-recipient inputs to L1, 2 and 5 arrived ~30 ms later still. *Curr Biol* 23:2121–2129. **Paywalled.** — [doi:10.1016/j.cub.2013.09.013](https://doi.org/10.1016/j.cub.2013.09.013)

### Inferences
- **Placing F2 on the latency scale.** Awake V1 responses start at ~40–50 ms (Self 2013; Dougherty 2019), so the team's interocular normalization at 20–35 ms after onset corresponds to ~60–85 ms after stimulus onset.
  - That is later than monocular COS (≈0 ms) and near-surround suppression (~5–10 ms), which have feedforward or LGN-inherited signatures.
  - It overlaps awake fast suppression in L2/3 and 4Cβ (~20 ms; Wang 2020), the horizontal boundary signal (~30 ms after onset; Self 2013) and the VSDI contrast-gain time constant (~31 ms; Reynaud 2012).
  - It precedes awake slow or feedback suppression (45–61 ms after excitation; Wang 2020) and feedback-layer filling-in (~60 ms after onset; Self 2013).
  - This puts interocular normalization in the local recurrent or lateral (OD-column-to-column) regime rather than in feedforward or long-range feedback regimes.
- **Comparison with Dougherty 2019.** The ~12 ms median gap there for binocular modulation of monocular neurons is shorter than the team's 20–35 ms. Their measure pooled facilitation and suppression, used only high contrast and measured *any* binocular modulation onset. The team's normalization-specific estimate need not match it.
- **Albrecht vs F2.** Albrecht 2002's "gain control fully expressed in the first 10 ms" concerns *monocular* self-normalization, which is plausibly inherited feedforward saturation (Priebe & Ferster 2006; Li 2006). It is compatible with F1 and F2 if interocular division is a separate, delayed cortical stage. Reynaud 2012's evolving c50 shows that population-level gain control can also be dynamic. The two studies conflict in measure (spikes vs VSDI) and stimulus (flashed vs drifting).
- **Contrast dependence (F2, faster and stronger at high contrast).** This matches:
  - Bair 2003: strong surround suppression arrives ~30 ms earlier than weak;
  - Smith 2006: surround strength anticorrelated with latency;
  - Nassi 2015: optogenetic drive switches from facilitating to suppressive as contrast rises;
  - SSN theory (Rubin 2015): sublinearity grows with drive, so stronger drive should recruit the pool sooner.
  - A purely feedforward depression account predicts suppression from the first spikes at all contrasts, which is not what F1 and F2 show.
- **Conductance predictions.** Anderson 2000 and Ozeki 2009 argue against untuned shunting as the main divisive mechanism.
  - For CSD, an ISN-style withdrawal of recurrent excitation predicts that normalization reduces *net sinks* in the same layers where spiking drops, at the same time. That is F2's "synaptic currents and spiking together."
  - A pure shunting mechanism could in principle reduce spiking with less change in net transmembrane current. CSD cannot separate E and I, though, so this is weak evidence either way.
- **Magnocellular-first input.** Maunsell & Gibson 1992 show the first ~7–10 ms of cortical activation is magnocellular. The near-linear first ~20 ms (F1) may therefore partly cover a period before parvocellular-driven and intracortical pool activity has built up. This is a hypothesis to test by checking whether onset linearity depends on contrast, which shifts M vs P weighting.

### Gaps
- I found no published *intracellular* or CSD latency for interocular suppression specifically. The cat interocular papers (Sengpiel & Vorobyov 2005; Li 2005) report no latencies.
- The Smith 2006 surround suppression and release latencies are relative to COS (~12 ms later), not reported as absolute values in the text I retrieved.
- Ringach 2003 and Shapley 2003 give no ms values in their abstracts.
- Wang 2020 has an internal inconsistency between the FS delay in its Discussion ("10 ms") and Table 1 (11–21 ms).
- The awake vs anesthetized difference is large: Haider 2013 and Vangeneugden 2019 show anesthesia changes inhibition and feedback. Latencies from AN cat and macaque studies (Bair, Smith, Henry, Bijanzadeh) may underestimate awake cortical contributions.
- PDFs to supply: Borg-Graham 1998; Menz & Freeman 2003; Self 2013. Albrecht 2002 and Anderson 2000 are J Neurophysiol and likely free after the APS embargo (not checked).

---

## Q3. Laminar evidence on where normalization and suppression are computed

### Takeaway
Simultaneous laminar recordings show that LGN-inherited suppression appears first in L4C, at response onset. Intracortically generated suppression appears first *outside* L4C (upper L4, L2/3, upper-deep layers for near surround; deep and superficial layers almost together for far surround) and last in L4C. Fast suppression is strongest in L4Cα and L6; slow, feedback-like suppression is strongest in L2/3 and L5. Binocular suppression of monocular neurons is present in all layers, including L4, and appears slightly earlier in L2/3. The team's depth-uniform interocular normalization (F4) does not reproduce any of these surround-suppression laminar sequences.

### Cited Findings
- **Bijanzadeh, Nurminen, Merlin, Clark & Angelucci 2018** (macaque, AN with sufentanil; linear arrays; MUA, LFP and CSD; latencies given as L, relative to the earliest contact in the same penetration):
  - **RF-only stimulation:** earliest spiking in L4C (L = 3.5 ± 1.35 ms).
  - **Full-field 20° grating:** earliest suppression in L4C (4 ± 1.15 ms), "at the same time as the onset of visual responses." Interpreted as "inherited from LGN as reduced feedforward excitation." Other layers: upper-L4 13.8 ± 2.48 ms, Deep-U 16.81 ± 2.97 ms, L2/3 22.71 ± 2.07 ms, Deep-L 32.17 ± 4.06 ms.
  - **Iso-oriented near-surround annuli not abutting the RF:** earliest suppression in upper-L4 (8.84 ± 2 ms) and Deep-U (19.55 ± 5.6 ms), then L2/3 (25.46 ± 3 ms). It was delayed in L4C (34.21 ± 7.1 ms) and Deep-L (28.22 ± 8.05 ms). Interpreted as "an intracortical origin, outside L4C."
  - **Tuned near-surround suppression:** first in upper-L4 (12.01 ms), then L2/3 (37.15 ms). It is "already orientation-tuned when it first emerges."
  - **Far surround:** Deep-U 8.23 ms, upper-L4 14.18 ms, L2/3 14.53 ms, Deep-L 15.76 ms; L4C 42.48 ± 12.66 ms, with many L4C contacts not suppressed at all.
  - **CSD:** near-surround squares evoked the earliest sinks in L2/3 and Deep-U, 13–16 ms after the earliest RF-evoked L4C sink. Far annuli evoked the earliest sinks in L1 and Deep-L, 22–28 ms after the RF-evoked L4C sink.
  - The authors caution that anesthesia effects on surround suppression under sufentanil are unknown.
  - *Neuron* 100:259–274. — [doi:10.1016/j.neuron.2018.08.020](https://doi.org/10.1016/j.neuron.2018.08.020); [preprint PDF](https://www.biorxiv.org/content/10.1101/171793.full.pdf)
- **Wang et al. 2020** (macaque, AW, laminar MUA):

  | Layer | FS index | SS index |
  |---|---|---|
  | L2/3 | 0.23 | 0.28 |
  | L4B | 0.46 | 0.17 |
  | L4Cα | 0.72 | 0.07 |
  | L4Cβ | 0.47 | 0.10 |
  | L5 | 0.39 | 0.23 |
  | L6 | 0.59 | 0.12 |

  - FS is "strongest in input layers (4C and 6)"; SS is "3 times stronger in output layers (2/3 and 5)."
  - SS was earliest near L1, interpreted as feedback.
  - FS in L4C correlated with orientation selectivity (r = −0.65). SS in L2/3 correlated with surround suppression (r = 0.61).
  - CSD was used only for laminar alignment.
  - — [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/)
- **Henry et al. 2020** (macaque, AN):
  - Untuned suppression was pronounced in 4Cα and 4B, the magnocellular stream, and attributed partly to M-LGN feedforward suppression plus local lateral inhibition.
  - Tuned suppression peaked earlier in 4Cα (60–75 ms) than in 4Cβ and L2/3 (~100 ms).
  - L6 early suppression resembled 4Cα; L5 showed delayed orthogonal facilitation.
  - — [PMC7253173](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173)
- **Sengpiel et al. 1998** (cat): response-gain control (typical of surround and interocular suppression) "prevailed in layer 4," whereas cells in layers 2/3, 5 and 6 "mainly showed contrast-gain control." — [doi:10.1016/S0042-6989(97)00413-6](https://doi.org/10.1016/S0042-6989(97)00413-6)
- **Sengpiel, Blakemore & Harrad 1995** (cat): monocular layer-4 units showed only modest, orientation-independent interocular inhibition like that in LGN. Strong orientation-dependent interocular suppression was a property of binocular cells. — [doi:10.1016/0042-6989(94)00125-6](https://doi.org/10.1016/0042-6989(94)00125-6)
- **Dougherty et al. 2019** (macaque, AW):
  - Monocular neurons were concentrated in L4, between 0.0 and 0.5 mm above the L4/L5 border.
  - Binocular modulation affected monocular neurons in *all* layers; 11/17 L4 monocular neurons were significantly modulated.
  - 7/10 facilitated monocular neurons were in L4, whereas suppression (more than two-thirds of monocular neurons) occurred across layers.
  - Modulation appeared first in L2/3 (47.6 ms), then L4 (50.1 ms), then L5/6 (54.7 ms).
  - — [PMC6363852](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852)
- **Cox et al. 2019** (macaque, AW): incongruent transient enhancement was confined to middle and upper layers, while deep layers showed transient suppression. Equiocular granular neurons kept congruent facilitation in the sustained window (+16.68 sp/s). In the sustained window the incongruent granular CSD sink stayed larger than the monocular one, while the congruent sink returned toward monocular levels. — [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- **Westerberg, Cox, Dougherty & Maier 2019** (macaque, AW; methodological template): repetition suppression was strongest outside granular layers. "Stimulus repetition first affected supragranular synaptic currents [CSD], which translated into a reduction of stimulus-evoked spiking across layers." *J Neurophysiol* 121:1938–1952. — [doi:10.1152/jn.00113.2019](https://doi.org/10.1152/jn.00113.2019)
- **Carlson, Mitchell, Westerberg, Poggi, Woodman & Maier 2025** (macaque, laminar MUA and LFP; preprint, not peer reviewed): interocular conflict "reliably reduces interlaminar coherence, particularly between granular and infragranular layers." During binocular rivalry flash suppression, coherence stayed reduced even when firing rates were unchanged. "V1 interlaminar coherence decreases with interocular conflict." bioRxiv. — [doi:10.1101/2025.10.04.680461](https://doi.org/10.1101/2025.10.04.680461)
- **Sceniak et al. 2006** (macaque, AN): geniculocortical afferents with strong surround suppression were found throughout 4C and 4A, so some suppression is already present in L4 *input*. — [doi:10.1152/jn.00734.2006](https://doi.org/10.1152/jn.00734.2006)
- **Maunsell & Gibson 1992**: superficial neurons became active ~15 ms after L4C. — [doi:10.1152/jn.1992.68.4.1332](https://doi.org/10.1152/jn.1992.68.4.1332)
- **Rodent laminar context:**
  - **Adesnik 2012** (mouse): the SOM-mediated surround circuit is in superficial layers. — [doi:10.1038/nature11526](https://doi.org/10.1038/nature11526)
  - **Keller, Roth & Scanziani 2020** (mouse; awake implied, since effects are "reduced by anaesthesia"): a feedback-generated second RF is "slower and … delayed" relative to the feedforward RF. It is found in feedback-recipient layers and "absent in the main input layer." Large uniform stimuli suppress it, and SOM cells are driven by them. *Nature* 582:545–549. — [doi:10.1038/s41586-020-2319-4](https://doi.org/10.1038/s41586-020-2319-4)
  - **Lien & Scanziani 2013** (mouse L4): thalamic excitation is ~1/3 of total L4 excitation; the rest is cortical amplification. — [doi:10.1038/nn.3488](https://doi.org/10.1038/nn.3488)
- **Nassi et al. 2015** (macaque, AW): optogenetically triggered normalization probably engaged mainly supragranular cells (light attenuation); there was no laminar analysis. — [PMC4534089](https://pmc.ncbi.nlm.nih.gov/articles/PMC4534089)

### Inferences
- **Laminar order as a diagnostic.** Bijanzadeh 2018 shows the laminar *order* of suppression onset can separate inherited from intracortical suppression. Inherited suppression appears in L4C at response onset; intracortical suppression appears outside L4C first and arrives in L4C 25–40 ms later; far/feedback suppression starts in deep layers and L1 (CSD).
  - F1 (no interocular normalization in the earliest L4C sink) rules out the inherited pattern for the interocular term.
  - F4 (similar timing across depth) matches neither the horizontal sequence (upper L4/L2/3 first, L4C last) nor the feedback sequence (L1/deep first).
  - The closest published analogue is far-surround suppression emerging "almost simultaneously in superficial and deep layers" except L4C. A key check: is interocular normalization in L4C equally fast in the team's data? If it is, that differs from every surround-suppression result.
- **Why depth uniformity could arise in L4.** Dougherty 2019's binocular suppression of L4 monocular neurons, and Cat-L4 response-gain-type suppression (Sengpiel 1998), show the input layer receives interocular suppression too. A depth-uniform normalization could come from:
  - L4 lateral or interneuronal connections between OD columns (as Dougherty 2019 proposes);
  - a network-wide ISN-like state in which a local perturbation changes E and I in all layers within one recurrent time constant (Ozeki 2009; Sanzeni 2020);
  - fast interlaminar loops, since L4C to superficial delays are only ~5–15 ms (Maunsell & Gibson 1992; Bijanzadeh 2018), close to the resolution of onset estimates in noisy MUA or CSD.
- **Cross-study contrast in fast suppression.** In Wang 2020 fast suppression is earliest (11 ms) and strongest in 4Cα, plausibly an M-pathway feedforward or local component. The team's depth-uniform interocular normalization differs, which again argues interocular normalization is not simply the monocular fast-suppression process.

### Gaps
- No laminar *optogenetic* manipulation of normalization in macaque V1 was found. Nassi 2015 was not laminar.
- No study reports CSD-resolved onset of *interocular* suppression by layer. Cox 2019 reports CSD only in coarse 50–100 ms and 150–250 ms windows. Carlson 2025 reports coherence, not CSD onset; its full text was not retrievable (rate-limited), so penetration counts and timing are unknown.
- Bijanzadeh 2018 and Henry 2020 are anesthetized; whether laminar latency orders hold in awake animals is untested except for Wang 2020 (awake, monocular, no surround laminar ordering in the same format).
- Bijanzadeh 2018 absolute latencies sit in Supplementary Tables S1–S6, which I did not read.

---

## Q4. Binocular-specific circuits: OD-column interactions, interocular suppression mechanisms, rivalry dynamics, and recent laminar and two-photon work

### Takeaway
Interocular suppression in V1 is cortical and inhibitory in cat (low temporal-frequency cutoff, adaptable, GABA-A dependent). It is gated by ongoing activity in the other eye, is strongest in neurons dominated by one eye, and is weak or contrast-conditional in LGN. In awake macaque, binocular modulation of even monocular L4 neurons is delayed relative to onset and depends on ocular dominance. Two-photon imaging shows OD-dependent suppression vs enhancement in superficial layers. The team's lab's earlier spiking work found the NDE's influence enters mainly as response gain on its excitatory drive. Human psychophysics and fMRI support interocular normalization (each eye divides the other) plus a separate monocular route. Rivalry models place fast, eye-of-origin mutual inhibition early and slow attentional or feature competition later. Opponency neurons reconcile rivalry with normalization.

### Cited Findings

#### Macaque laminar and two-photon
- **Zhang, Zhao, Jiang, Tang & Yu 2024** (macaque, AW, two-photon GCaMP5G; superficial layers at ~200 µm and one other depth; 10,168 neurons, 9,390 tuned; 0.9-contrast Gabors; 1 s on):
  - On average, neurons preferring either eye were *suppressed* by binocular stimulation, while binocular neurons were *enhanced*. Binocular modulation vs |ODI|: y = −0.36x + 0.14 (p < 0.001).
  - Model: each eye's response is divided by an OD-dependent interocular suppression factor (larger for more monocular neurons), then summed.
  - Calcium integrates up to ~1 s, so there is no timing information.
  - "Ocular dominance-dependent binocular combination of monocular neuronal responses in macaque V1." *eLife* 13:RP92839. — [doi:10.7554/eLife.92839](https://doi.org/10.7554/eLife.92839)
- **Mitchell, Carlson, Westerberg, Cox & Maier 2023** (macaque, AW; 502 units, 25 penetrations; 6×6 contrasts {0, .05/.055, .11, .225, .45, .9}; 50–200 ms window):
  - Contrast in the DE was weighted more than contrast in the NDE, and the asymmetry covaried with OD.
  - The NDE effect was best captured as response gain: an Rmax-only fit gave adjusted R² 0.92, vs 0.82 for c50-only and 0.79 for n-only.
  - The authors place the OD weight on the NDE's *excitatory (numerator) drive*, "not on the normalizing suppressive pool." Fitted weights correlated with ODI (adj. R² = 0.56).
  - Interocular normalization with OD gave the best fit (median adjusted R² 0.82; AIC improved 21.5% by adding OD).
  - "A role for ocular dominance in binocular integration." *Curr Biol* 33:3884–3895. — [doi:10.1016/j.cub.2023.08.019](https://doi.org/10.1016/j.cub.2023.08.019); [PMC10530424](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530424)
- **Dougherty et al. 2019**: the authors propose "lateral connections within L4 (including interneurons crossing ocular dominance columns) and interlaminar feedback from L2/3 or L5/6" as the substrate for delayed binocular suppression of monocular neurons. — [PMC6363852](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852)
- **Cox et al. 2019**: the authors propose "short-range horizontal connections between neighboring ocular dominance columns" for the sustained suppression. Deep-layer suppression might come from retinogeniculate collaterals bypassing L4C or from interlaminar spread. — [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- **Reynaud, Roux, Chemla, Chavane & Hess 2019** (macaque, AN; VSDI with 6 Hz / 10 Hz frequency tagging of the two eyes; abstract only): stimulating one eye suppressed population activity driven by the other. An interocular normalization model fit well. The mechanisms (local, feedforward, feedback or long-range) are undetermined, and no timing is reported. *Ann Eye Sci* 4:AB015. — [doi:10.21037/aes.2019.AB015](https://aes.amegroups.org/article/view/5220/html)

#### Cat and rodent
- **Sengpiel & Vorobyov 2005; Li et al. 2005; Sengpiel et al. 1995; Sengpiel & Blakemore 1994**: interocular suppression is cortical, GABA-A dependent, adaptable and gated by activity in the other eye. The 1995 paper proposes inhibition between adjacent OD columns. See Q1 for numbers.
- **Longordo, To, Ikeda & Stuart 2013** (mouse L2/3 pyramidal, in vivo whole-cell; prep. not in abstract):
  - Binocular inputs sum *sublinearly* in an amplitude-dependent way, most strongly at preferred orientation and disparity and at highest contrast.
  - The sublinearity is postsynaptic and "requires balanced recruitment of inhibition."
  - It acts as "a divisive form of gain control."
  - *Nat Neurosci* 16:714–723. **Paywalled.** — [doi:10.1038/nn.3394](https://doi.org/10.1038/nn.3394)
- **Zhao, Liu & Cang 2013** (mouse simple cells, AN and AW, whole-cell): binocular Vm integration is sublinear, which keeps non-preferred-orientation depolarizations below threshold. Conductance-based simulations, plus genetic and pharmacological tests, show "inhibition promotes sublinear binocular integration." *Nat Commun* 4:2088. — [doi:10.1038/ncomms3088](https://doi.org/10.1038/ncomms3088)
- **Cang, Fu & Tanabe 2023** (review, mostly mouse, ferret, tree shrew): most OD studies use only monocular stimulation, "which could lead to a mischaracterization of binocularity." The circuit basis of interocular interactions remains largely unknown. *Front Neural Circuits* 17:1084027. — [doi:10.3389/fncir.2023.1084027](https://doi.org/10.3389/fncir.2023.1084027)

#### Human psychophysics and fMRI
- **Moradi & Heeger 2009** (human fMRI, attention diverted):
  - Linear summation failed at ≥10% contrast. A model in which "activity from each eye reduced the gain for the other eye" fit very well.
  - Dichoptic plaids evoked the same activity as monocular plaids.
  - Cross-orientation suppression was weaker than same-orientation interocular suppression.
  - *J Vis* 9(3):13. — [doi:10.1167/9.3.13](https://doi.org/10.1167/9.3.13)
- **Ding & Sperling 2006** (human psychophysics): each eye exerts gain control on the other eye's signal in proportion to its own contrast energy, *and* on the other eye's gain control. Very low contrasts add linearly; the nonlinearity appears only as contrast rises. A one-parameter fit explained 95% of variance. *PNAS* 103:1141–1146. — [doi:10.1073/pnas.0509629103](https://doi.org/10.1073/pnas.0509629103)
- **Meese, Georgeson & Baker 2006** (human; two-stage model): nearly linear contrast transduction before binocular summation, plus suppressive ocular interactions in gain control. *J Vis* 6:1224–1243. — [doi:10.1167/6.11.7](https://doi.org/10.1167/6.11.7)
- **Meese & Baker 2009** (human): monocular cross-orientation masking is a power function of speed (TF/SF), consistent with a "fast-acting, early locus." Dichoptic masking is speed-invariant, implying different mechanisms. *J Vis* 9(5):2. — [doi:10.1167/9.5.2](https://doi.org/10.1167/9.5.2)
- **Baker, Meese & Summers 2007** (human): the interocular route is spatially tuned and adaptable, and its weight increases with duration. The ipsiocular route is broadband, non-adaptable, and its weight decreases with duration. The authors propose ipsiocular suppression "precedes" interocular suppression. — [doi:10.1016/j.neuroscience.2007.01.030](https://doi.org/10.1016/j.neuroscience.2007.01.030)

#### Rivalry models
- **Wilson 2003**: rivalry data "require at least two hierarchic rivalry stages." Competitive inhibition at the first (monocular or eye-based) stage can be removed by stimulus dynamics, exposing a later stage. *PNAS* 100:14499–14503. — [doi:10.1073/pnas.2333622100](https://doi.org/10.1073/pnas.2333622100)
- **Said & Heeger 2013**: conventional rivalry models based on strong competition between orientation-selective neurons wrongly predict rivalry for plaids. Adding ocular-opponency neurons gives rivalry only for dichoptic gratings and "reconciles models of binocular rivalry with the divisive normalization framework." *PLoS Comput Biol* 9:e1002991. — [doi:10.1371/journal.pcbi.1002991](https://doi.org/10.1371/journal.pcbi.1002991)
- **Li, Rankin, Rinzel, Carrasco & Heeger 2017**: rivalry is driven by attentional modulation, which is feature-selective and "relatively slow," and by mutual inhibition, which is eye-of-origin selective and "relatively fast." *PNAS* 114:E6192–E6201. — [doi:10.1073/pnas.1620475114](https://doi.org/10.1073/pnas.1620475114)

### Inferences
- **Bearing on F3.** Zhang 2024 and Dougherty 2019 both find that monocularly dominated neurons are *more* suppressed binocularly, while their excitation is driven mostly by one eye. So the NDE must contribute more to suppression than to excitation in those neurons, which is qualitatively F3. Mitchell 2023 placed the OD weight on the NDE numerator (down-weighting its excitation) with an untouched pool. That is algebraically the same asymmetry: the pool weights the NDE relatively more than excitation does. Whether the pool is eye-balanced or NDE-biased is not separable without dynamic or laminar constraints like F2 and F4.
- **Circuit candidates for an NDE-heavy pool.**
  - Lateral inhibition from the neighbouring OD column (Sengpiel 1995; Cox 2019; Dougherty 2019). A neuron in a DE column would be suppressed by interneurons driven from the adjacent NDE column.
  - Broadly binocular inhibitory interneurons. L4 FS cells' spike tuning is only moderately broader than RS cells' (orientation and SF, not ocularity), per Cardin, Palmer & Contreras 2007, *J Neurosci* 27:10333 — [doi:10.1523/JNEUROSCI.1692-07.2007](https://doi.org/10.1523/JNEUROSCI.1692-07.2007). The species is not stated in the abstract; the lab normally uses cat. I found no macaque data on interneuron ocularity.
- **Bearing on F2.** Mutual inhibition between eyes is "relatively fast" in rivalry models (Li 2017), and interocular psychophysical suppression *grows* with duration (Baker 2007). Both fit an interocular divisive term that builds up over tens of ms after onset rather than being present at onset.
- **Bearing on F1.** Ding & Sperling 2006 predict linear binocular addition at low contrast and normalization only at higher contrast. An onset period of near-linear summation that becomes sublinear as the pool builds is the time-domain analogue: low effective pool drive early behaves like low contrast.

### Gaps
- No direct anatomical or physiological identification of the interneuron class (PV, SOM or other) mediating interocular suppression in primate V1 was found. Ocularity of macaque interneurons is not reported in any source I retrieved.
- No published measurement of interocular suppression latency in LGN (Dougherty 2021 could not measure it) or of its contrast dependence in time.
- Longordo 2013 (Nat Neurosci) is paywalled; its preparation and any timing data need the PDF.
- Carlson 2025 full text could not be retrieved (bioRxiv rate-limited).
- Reynaud 2019 is an abstract only.

---

## Q5. Which empirical signatures in time-resolved, laminar, dichoptic contrast data would separate feedforward from recurrent interocular normalization?

### Takeaway
The literature offers several practical discriminating signatures: (i) whether interocular suppression appears in the earliest L4C sink; (ii) the laminar order of suppression onset; (iii) the latency of suppression relative to each eye's own excitatory onset and how it scales with contrast; (iv) whether suppression appears as withdrawal of net sinks together with spiking (recurrent/ISN) or as spiking loss without matching sink loss; (v) dependence on ocular dominance and interocular contrast asymmetry; (vi) the transient-then-steady-state structure predicted by ISN. Temporal-frequency and adaptation tests are the classic cat-literature discriminators but are not available in static-grating data. F1, F2 and F4 together most closely match a local recurrent, cortical interocular normalization, not an LGN-inherited or long-range feedback mechanism.

### Cited Findings
- **Signature 1: suppression at response onset in L4C means inherited; suppression delayed in L4C means intracortical.** Full-field surround suppression appeared in L4C "at the same time as the onset of visual responses" (inherited from LGN), whereas near and far surround suppression was delayed in L4C by tens of ms (Bijanzadeh 2018, macaque, AN — [doi:10.1016/j.neuron.2018.08.020](https://doi.org/10.1016/j.neuron.2018.08.020)). Monocular COS arrives with a latency shorter than excitatory onset latency (Smith 2006, macaque, AN — [PMC6674161](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161)).
- **Signature 2: laminar order of onset.** Horizontal near-surround suppression appears first in upper-L4, Deep-U and L2/3. Feedback-driven far-surround sinks start in L1 and Deep-L (Bijanzadeh 2018). Feedback-recipient layers L1, 2 and 5 receive figure-ground input ~30 ms after horizontal boundary signals (Self 2013, AW — [doi:10.1016/j.cub.2013.09.013](https://doi.org/10.1016/j.cub.2013.09.013)). Repetition suppression appears first in supragranular CSD and then in spiking across layers (Westerberg 2019, AW — [doi:10.1152/jn.00113.2019](https://doi.org/10.1152/jn.00113.2019)).
- **Signature 3: temporal-frequency cutoff and adaptation (precortical vs cortical).** Precortical suppression extends to ≥20–32 Hz and resists adaptation. Cortical interocular suppression cuts off at ~8–16 Hz, is abolished by adaptation (33.4% to 5.6%) and is GABA-A dependent (Freeman 2002 — [doi:10.1016/S0896-6273(02)00819-X](https://doi.org/10.1016/S0896-6273(02)00819-X); Sengpiel & Vorobyov 2005 — [PMC6725266](https://pmc.ncbi.nlm.nih.gov/articles/PMC6725266); Li 2005 — [doi:10.1152/jn.00203.2005](https://doi.org/10.1152/jn.00203.2005); Webb 2005 — [doi:10.1523/JNEUROSCI.3414-05.2005](https://doi.org/10.1523/JNEUROSCI.3414-05.2005)).
- **Signature 4: E and I withdrawal, plus a transient inhibitory overshoot (ISN).** Under cortical (recurrent) suppression, both E and I fall at steady state (E −51%, I −38%), with a 30–50 ms transient inhibitory increase at onset (Ozeki 2009, cat, AN — [PMC2691725](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725)). Normalization by distal cortical drive reduces both EPSC and IPSC (Sato 2016, mouse, light AN — [PMC4817833](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817833)). Feedforward COS also lowers both E and I (Priebe & Ferster 2006 — [doi:10.1038/nn1660](https://doi.org/10.1038/nn1660)), so withdrawal alone does not discriminate. Timing relative to the first sink does.
- **Signature 5: drive dependence (SSN).** In the SSN, weak drive gives linear or supralinear summation and strong drive gives sublinear summation and stronger suppression, with E and I cells behaving alike (Rubin 2015 — [doi:10.1016/j.neuron.2014.12.026](https://doi.org/10.1016/j.neuron.2014.12.026)). Awake mouse L2/3 shows the E/I ratio falling with contrast (Adesnik 2017 — [doi:10.1016/j.neuron.2017.08.014](https://doi.org/10.1016/j.neuron.2017.08.014)). Stronger surround suppression arrives earlier (Bair 2003 — [PMC6740744](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744)). Optogenetic drive in alert macaque switches from facilitation to suppression as contrast rises (Nassi 2015 — [PMC4534089](https://pmc.ncbi.nlm.nih.gov/articles/PMC4534089)).
- **Signature 6: ocular dominance.** Monocular-stage gain control (Truchard 2000 — [PMC6772209](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209)) acts within each eye. Cortical interocular suppression grows with eye preference (Dougherty 2019 — [PMC6363852](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852); Zhang 2024 — [doi:10.7554/eLife.92839](https://doi.org/10.7554/eLife.92839)). In LGN, suppression appears only when the driving eye's contrast is low and the other eye's is high (Dougherty 2021 — [PMC8035044](https://pmc.ncbi.nlm.nih.gov/articles/PMC8035044/)).
- **Signature 7: delayed vs simultaneous gain control in models.** Models in which gain control is slower than the linear drive reproduce sub-additive dynamics and contrast-dependent latency (Zhou 2019 — [doi:10.1371/journal.pcbi.1007484](https://doi.org/10.1371/journal.pcbi.1007484)). Recurrent normalization circuits reproduce V1 dynamics (Heeger & Zemlianova 2020 — [doi:10.1073/pnas.2005417117](https://doi.org/10.1073/pnas.2005417117)). Feedforward population gain control can also produce complex dynamics (Sit 2009 — [doi:10.1016/j.neuron.2009.08.041](https://doi.org/10.1016/j.neuron.2009.08.041)). So dynamics alone do not uniquely identify a recurrent mechanism.

### Inferences
Concrete tests the team can run on existing 6×6 dichoptic, laminar, 1 ms-resolution data:
1. **Earliest-sink test (feedforward vs cortical).** Fit the interocular term in the first L4C sink, in 1–2 ms bins from sink onset, separately for low-DE/high-NDE cells. LGN-inherited suppression would be expected only there (Dougherty 2021) and should be present at sink onset. A null across all contrast pairs (F1) argues against a precortical interocular component.
2. **Latency vs each eye's onset.** Measure normalization onset relative to the *NDE's* own excitatory onset, not just the binocular onset. An intracortical, NDE-driven pool should lag NDE-driven excitation by one or more synaptic or recurrent delays. A feedforward pool would not lag it. Contrast-dependent shifts in the lag would also be expected. If normalization latency tracks NDE onset latency across NDE contrasts (lower contrast meaning later onset; Reynaud 2012 found 87 vs 40 ms over the contrast range), the pool is driven by the NDE signal arriving in cortex.
3. **Laminar ordering at high resolution.** Estimate normalization onset per layer with bootstrap CIs, and compare with the Bijanzadeh 2018 templates (near-horizontal: upper-L4/L2/3 first, L4C last; far/feedback: L1 and deep first; inherited: L4C first). True depth simultaneity (F4) within CIs of a few ms points to a distributed recurrent state change, ISN-like, rather than a laminar-specific pathway. Any L1 or deep-first onset would implicate feedback.
4. **CSD sinks vs sources.** If normalization appears as reduced sinks in L4C and supragranular layers at the same time as reduced MUA (F2), that matches withdrawal of recurrent excitation (Ozeki 2009; Sato 2016). A transient *increase* in sources or sinks just before suppression in some layer would be a candidate ISN inhibitory overshoot (Ozeki's 30–50 ms transient).
5. **Contrast scaling of onset and strength (F2).** "Faster and stronger at higher contrast" is predicted by SSN drive dependence (Rubin 2015) and by Bair 2003's strength–latency relation. A pure thalamocortical-depression model predicts suppression locked to the first spikes at all contrasts. Fitting a delayed-gain-control model (Zhou 2019 style) per layer gives a pool time constant to compare with Reynaud 2012's ~31 ms.
6. **OD-dependent pool weights over time.** Fit time-resolved weights of each eye in the numerator and the pool. Monocular-stage gain control (Truchard 2000) predicts pool weight ∝ own-eye drive. Lateral OD-column inhibition predicts an NDE-heavy pool (F3) that emerges *with* the delayed normalization, not at onset.
7. **Monocular-neuron test.** In L4 monocular units (Dougherty 2019), the NDE produces no excitation, so any NDE-driven suppression must be cortical or LGN-mediated. Its latency relative to the unit's onset is a clean estimate of the interocular pool delay without excitatory confounds.

How F1–F4 map onto mechanisms in this literature:

| Finding | Interpretation |
|---|---|
| F1 | No precortical or feedforward interocular normalization at onset |
| F2 | Delay in the local recurrent or lateral range (≈Wang 2020 FS in L2/3 and 4Cβ; ≈Reynaud 2012 τ), shorter than feedback (Wang 2020 SS; Bijanzadeh 2018 far surround; Self 2013); contrast-dependent as SSN or ISN predict |
| F3 | Matches OD-dependent interocular suppression (Zhang 2024; Dougherty 2019) and the OD-numerator weighting of Mitchell 2023 |
| F4 | Atypical compared with surround-suppression laminar sequences; points to a column-wide recurrent mechanism with fast interlaminar spread |

The strongest single falsifier of "recurrent" would be interocular normalization present in the L4C sink at onset for any contrast pair. The strongest falsifier of "local" (pointing to feedback) would be an L1- or deep-first onset.

### Gaps
- No published dataset combines dichoptic contrast manipulation with layer-resolved CSD onset analysis at ms resolution. The team's data appear to be the first. Cox 2019 and Dougherty 2019 used coarse windows or high contrast only.
- Temporal-frequency and adaptation discriminators (Sengpiel & Vorobyov 2005; Li 2005) have not been applied in awake macaque. Static flashed gratings cannot test them without new experiments.
- I found no direct primate test of ISN signatures (paradoxical inhibitory responses) during interocular suppression. Nassi 2015 tested sub-additivity but not ISN paradoxical effects or laminar specificity.
- Interneuron-type contributions (PV vs SOM) to interocular normalization have not been tested in any species in the sources found. The mouse binocular studies (Longordo 2013; Zhao 2013) implicate inhibition generally.
