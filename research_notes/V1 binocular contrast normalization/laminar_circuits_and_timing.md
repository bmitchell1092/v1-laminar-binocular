# Laminar circuit dynamics and timing of contrast responses and normalization in primate V1, including inter-columnar interactions

Scope notes for the report writer:
- Species is macaque unless marked otherwise. Cat, marmoset and mouse data are labelled.
- **Anesth.** means sufentanil-anesthetized and paralyzed. **Awake** means fixating.
- Numbers are as reported by the source. Where I read a number from a secondary summary rather than the primary paper, I say so.
- **Attribution corrections to the brief:**
  - The 2014 *Nat Commun* travelling-wave paper is by Muller, **Reynaud, Chavane & Destexhe**, not Muller/Reynolds/Sejnowski.
  - Davis et al. 2020 *Nature* recorded **MT in marmosets**, not V1.
  - Davis et al. 2023 *eLife* recorded **no V1** (V4, MT, PFC).
  - I could not find the "Li, Reynolds et al." or "Chen, Dhruv, Shapley" laminar contrast papers named in the brief (see Gaps).
- **Design context**, stated once and reused in the "Design fit" lines:
  - The team recorded with 32-channel, 100 µm linear probes in awake macaques, about 30 penetrations.
  - On 7 days there were two probes in separate V1 columns, both written to the same acquisition files, so inter-probe timing is sample-locked.
  - Stimuli were flashed static gratings of about 317 ms in a 6×6 dichoptic contrast matrix {0, .055, .11, .225, .45, .9}.
  - Signals are LFP/CSD, analog MUA and about 70 SUs. Photodiode timing is good to about 1 ms.
  - Per the legacy-repo notes, the NSx-vs-NEV clock offset of 3.4 ms must be corrected before any latency comparison between continuous signals and event codes.

---

## 1. Laminar differences in contrast response functions (c50, exponent, saturation; contrast gain vs response gain)

### Takeaway
No modern study reports Naka-Rushton c50, n and Rmax layer by layer across all V1 layers in awake macaque. The laminar picture is still built from three things:
- old 4C sub-lamina contrast-sensitivity and 2-DG data: 4Cα/4B/6 are high-sensitivity (magno-dominated) and 4Cβ is low-sensitivity (parvo-dominated);
- area-level CRF comparisons;
- recent model predictions.

The team's per-layer CRFs across a 6×6 dichoptic matrix would fill a real gap, especially the contrast-gain vs response-gain classification of *interocular* normalization by layer.

### Cited Findings
**Layer 4 sub-laminae (old physiology and 2-DG)**
- Hawken & Parker (1984, *Exp Brain Res* 54:367–372, doi:10.1007/BF00236238), with layers assigned by histology:
  - The magno-recipient IVcα held high-contrast-sensitivity neurons and the parvo-recipient IVcβ held low-sensitivity neurons.
  - IVa and IVb held a mix.
  - Non-oriented cells appeared only in the parvo-recipient sub-laminae (IVa, IVcβ).
  - [ORA Oxford record](https://ora.ox.ac.uk/objects/uuid:0efd6468-18fe-4aec-a3b9-6ec8bc9400e5)
- Vanni, Hokkanen, Werner & Angelucci (2020, *Cereb Cortex* 30:3483–3517, doi:10.1093/cercor/bhz322) compiled the physiology:
  - Contrast threshold rises up to 3.5-fold from the top of 4Cα to the bottom of 4Cβ. This is a threshold measure, not c50.
  - Citing Sclar et al. 1990, V1 has lower contrast sensitivity than V2 and V5/MT, and the difference is carried by c50 rather than the exponent.
  - I saw both points only in a search summary, not in the review's physiology section. [Vanni et al. 2020](https://academic.oup.com/cercor/article/30/6/3483/5691251)
- Tootell, Hamilton & Switkes (1988, *J Neurosci* 8:1594–1609), 2-deoxyglucose with achromatic gratings at several contrasts:
  - Uptake is essentially saturated at ≥50% contrast.
  - At 8% contrast, 4Cα showed clear uptake and 4Cβ showed none.
  - Layers 4B, 4Cα and 6 were strongly active, while 2, 3, 4A and 5 were weak.
  - [PMC6569196](https://pmc.ncbi.nlm.nih.gov/articles/PMC6569196); [MGH abstract](https://nmr.mgh.harvard.edu/publications/journal_articles/6108)

**CRF dynamics and saturation**
- Albrecht, Geisler, Frazor & Crane (2002, *J Neurophysiol* 88:888–913), anesth. monkey and cat, 200 ms stationary gratings:
  - PSTH shape is roughly constant across contrast within a cell. Scaling plus shifting the PSTH explains about 95% of the variance, which implies a contrast-dependent latency shift, but the abstract gives no ms values.
  - Contrast normalization and expansive half-squaring are "fully expressed within the first 10 ms" of the response.
  - One Naka-Rushton fit, with the same c50 and exponent, fits from the first to the last 10 ms.
  - [PubMed 12163540](https://pubmed.ncbi.nlm.nih.gov/12163540/)
- Yoshida, Chen, Geisler & Seidemann (bioRxiv, June 2026, doi:10.64898/2026.06.15.732419), fixating macaques, electrophysiology plus calcium imaging, Gaussian increments up to 900% Weber contrast:
  - Average responses keep rising well above 100% Weber contrast.
  - "Only a small minority of neurons saturate below 100% Weber contrast."
  - No laminar breakdown in the abstract.
  - [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.06.15.732419v1)
- Network model of 4Cα plus L6 feedback (Xiao, Lin & Young, arXiv:2606.21785, June 2026):
  - L4 populations reach half-max at about 10% contrast (cardinal-preferring) or 10–25% (mixed-input), then saturate.
  - Only L4 is modelled; L2/3 is not simulated.
  - [arXiv](https://arxiv.org/html/2606.21785)
  - An earlier comprehensive V1 model addresses the same issue: "Contrast response in a comprehensive network model of macaque V1" (Chariker, Shapley, Hawken & Young, *J Vis* 2020). I saw only the title. [JOV](https://jov.arvojournals.org/article.aspx?articleid=2765264)

**Binocular facilitation and interocular normalization (same lab, not laminar-resolved in abstracts)**
- Mitchell, Dougherty, Westerberg, Carlson, Daumail, Maier & Cox (2022, *iScience* 25:104182, doi:10.1016/j.isci.2022.104182), N = 314 units, 2 monkeys, 250 ms stimuli, contrasts 0 / ~.2 / ~.45 / ~.85:
  - Binocular (concordant) stimuli gave mean modulation index 0.057 (SD 0.067).
  - Facilitation scaled most with contrast in the early transient (50–150 ms).
  - Later (150–250 ms), facilitation lasted longer at low contrast and turned into binocular suppression at medium and high contrast.
  - A one-parameter response-gain model best matched the early CRF shift, although contrast-gain variants performed comparably.
  - [PubMed 35494250](https://pubmed.ncbi.nlm.nih.gov/35494250/)
- Mitchell, Carlson, Westerberg, Cox & Maier (2023, *Curr Biol* 33:3884–3895.e5, doi:10.1016/j.cub.2023.08.019), with monocular, equal-contrast binocular and unequal-contrast dichoptic gratings:
  - Contrast in a neuron's dominant eye was weighted more strongly than contrast in its non-dominant eye, and the asymmetry tracked ocular dominance.
  - Divisive normalization with OD weighting fit best, with *interocular* normalization performing best.
  - [PubMed 37657450](https://pubmed.ncbi.nlm.nih.gov/37657450/); title and journal also confirmed via the [eLife 92839 review thread](https://elifesciences.org/reviewed-preprints/92839v2/reviews)
- Zhang, Zhao, Tang & Yu (2024, *eLife*, doi:10.7554/eLife.92839), two-photon GCaMP5 of superficial V1, 3 awake macaques, 9,390 tuned neurons, 0.9-contrast gratings:
  - Binocular viewing suppresses monocular (strongly OD-biased) neurons and enhances binocular neurons, so responses converge across ocular dominance.
  - An OD-dependent divisive interocular suppression plus binocular summation model fits.
  - Dichoptic unequal contrast was *not* tested.
  - [eLife RP 92839](https://elifesciences.org/reviewed-preprints/92839v2)

**Surround context and contrast across layers**
- Sceniak, Ringach, Hawken & Shapley (1999, *Nat Neurosci* 2:733–739):
  - Spatial summation extent is on average 2.3× larger at low contrast.
  - This holds "in cells throughout V1" and is independent of surround inhibition.
  - [Nature](https://www.nature.com/articles/nn0899_733)
- Geniculocortical afferents have contrast-invariant summation, so the contrast dependence arises intracortically. [PubMed 16928793](https://pubmed.ncbi.nlm.nih.gov/16928793/)
- Surround suppression of firing rate is strongest in supragranular (SG) layers (Nurminen, Bijanzadeh & Angelucci 2026, *eLife* 15:e86334, anesth., 100% contrast). [PMC12931924](https://pmc.ncbi.nlm.nih.gov/articles/PMC12931924/)
- Surround suppression is more sharply orientation-tuned above L4C (Henry et al. 2013; Shushruth et al. 2013, as cited by Bijanzadeh et al. 2018). [Bijanzadeh et al. author manuscript](https://par.nsf.gov/biblio/10110109)

### Inferences
- **Two streams.** The robust laminar CRF prediction from the literature is two-stream. Magno-dominated compartments (4Cα, 4B, 6) should have low c50 and early saturation. Parvo-dominated compartments (4Cβ, and partly 2/3 and 5) should have higher c50 and less saturation at 0.9 Michelson.
  - "Deep layers saturate less/more" is therefore not a single prediction: L6 is magno-like (Henry et al. 2020, §2; Tootell 1988), while L5 is mixed.
- **Contrast range.** The 2026 Seidemann-lab result implies that apparent "saturation" at 0.9 Michelson grating contrast may be partly a ceiling of stimulus range. For gratings, Michelson 0.9 is near the physical maximum, so Rmax is estimable, but a "supersaturation" decline should be checked per layer.
- **Albrecht invariance.** Albrecht et al. 2002 imply that c50 and n are time-invariant, which is testable per layer. Bijanzadeh, Wang and Henry (§2) show suppression components with 10–40 ms delays, which predicts that layers with strong slow suppression (output layers 2/3, 4B, 5) should show a time-varying CRF.
- **Open question: does interocular normalization act as contrast gain (c50 shift) or response gain (Rmax scaling), and does this differ by layer?**
  - Discriminating measurement: fit each unit's or layer's monocular CRF with the other eye held at each of the 6 contrasts, then test whether c50 or Rmax shifts as a function of fellow-eye contrast and layer.
  - Design fit: **Yes.** The full 6×6 matrix is exactly this design. Single-unit counts per layer (about 70 SUs total) will be thin, so analog MUA/CSD per layer is the main workhorse.
- **Open question: are 4Cα and 4Cβ CRFs separable in awake data?**
  - Discriminating measurement: CSD/MUA at ≤50 µm resolution, or histology.
  - Design fit: **Partial.** At 100 µm spacing, 4Cα vs 4Cβ is at most 2–3 contacts. Cox et al. 2019 cite ±50 µm alignment error with the same probe class. Treat 4C as one or two compartments, and use latency (4Cα earliest; see §2) as a sub-lamina marker.

### Gaps
- I found no source with per-layer Naka-Rushton c50, n or Rmax for all layers in macaque V1 (awake or anesthetized). The closest are 4C sub-laminar sensitivity (1984), 2-DG (1988) and models (2020, 2026).
- Sclar, Maunsell & Lennie (1990, *Vision Res* 30:1–10) values for c50 and exponents in LGN, V1, V2 and MT were not retrievable; only the citation was confirmed ([ordb record](https://ordb.biotech.ttu.edu/odordb/Data/44537)). Do not quote numbers from memory.
- Hansen et al. 2012 (*Neuron*, laminar correlated variability) and Henry, Joshi, Xing, Shapley & Hawken 2013 were not retrieved directly. Xing, Yeh, Burns & Shapley 2012 PNAS used only >90% contrast and has no CRF ([PMC3427063](https://pmc.ncbi.nlm.nih.gov/articles/PMC3427063)).
- Nandy et al. 2017 laminar attention is, to my knowledge, a V4 study; this was not verified this session.
- The Mitchell 2024 Vanderbilt dissertation, "Neuronal Basis of V1 Binocular Normalization", proposes a single divisive normalization operating within and between eyes ([ir.vanderbilt.edu](https://ir.vanderbilt.edu/handle/1803/19176)). Its laminar results were not visible in the abstract.

---

## 2. Timing: laminar latency sequence, contrast-dependent latency, timing of suppression/normalization vs excitation, and timing of binocular interactions

### Takeaway
The sequence is consistent across methods. 4C (especially 4Cα) and L6 are first, then the extragranular layers. Superficial layers lag 4C by about 15–25 ms, and 4Cβ lags 4Cα by about 20 ms (anesth.).

Suppression has at least two clocks:
- A fast component at or near response onset (0–10 ms after excitation). It is untuned or cross-orientation, strongest in input layers, and plausibly LGN-inherited or fast local inhibition.
- A slow component 15–60 ms later. It is tuned, surround or feedback-like, and strongest in output layers (2/3, 4B, 5).

Binocular interactions follow the same two-stage logic. Early (50–100/150 ms) facilitation is followed by later (150–250 ms) contrast- and conflict-dependent suppression.

### Cited Findings
**Laminar latency sequence**
- Maunsell & Gibson (1992, *J Neurophysiol* 68:1332–1344):
  - Most V1 units have 30–50 ms latencies.
  - The shortest-latency units were all in layers 4C or 6 and were among the most transient.
  - Latencies increase monotonically toward the surface, with the most superficial neurons active about 15 ms after 4C.
  - [PubMed 1432087](https://pubmed.ncbi.nlm.nih.gov/1432087/)
- Nowak, Munk, Girard & Bullier (1995, *Vis Neurosci*), anesth.: the shortest V1 latencies are in 4Cα and 4B, and 4Cβ latencies are on average about 20 ms longer. I read this from a search summary. [Cambridge](https://www.cambridge.org/core/product/2E43515431CBCFCC2451939F734E3427)
- Schroeder, Mehta & Givre (1998, *Cereb Cortex* 8:575–592), awake, CSD plus MUA: V1 and V2 show a feedforward profile, with excitation beginning in Lamina 4 and then the extragranular laminae. [abstract](https://pesquisa.bvsalud.org/sms/resource/pt/mdl-9823479)
- Bijanzadeh, Nurminen, Merlin, Clark & Angelucci (2018, *Neuron* 100:259–274), anesth., 4 macaques, 10 perpendicular penetrations. Relative CSD sink onset after RF stimulation, versus the earliest layer:
  - L4C: 0 ms
  - upper-L4: 16.9 ± 7.0 ms (n.s. from L4C)
  - Deep-upper: 14.3 ± 2.0 ms
  - L2/3: 23.4 ± 3.1 ms
  - L1: 23.4 ± 1.6 ms
  - Deep-lower: 61 ± 6.7 ms
  - Absolute L4C onset was about 39 ms for CSD but about 23 ms for spiking, so CSD and spike latencies are not directly comparable.
  - [author manuscript, par.nsf.gov](https://par.nsf.gov/biblio/10110109) (text extracted from the PDF)
- Gieselmann & Thiele (2022, *eLife* 11:e62949), awake, 16-contact laminar:
  - Granular is fastest, then infragranular, with supragranular generally lagging.
  - Stimulus size altered latencies in the granular and infragranular layers.
  - [eLife](https://elifesciences.org/articles/62949)
- Wang et al. (Xing lab; 2020, *J Neurosci* 40:7436–7450), awake, 24-channel 100 µm U-probe, 4Cα located at the earliest CSD sink: responses start in the input layers (4C, 6) and then reach the output layers (2/3, 4B, 5). [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183)
- Self, van Kerkoerle, Supèr & Roelfsema (2013, *Curr Biol* 23:2121–2129, doi:10.1016/j.cub.2013.09.013), figure-ground task:
  - Boundary signals appear in upper L4 and superficial layers about 70 ms after onset, attributed to horizontal connections.
  - Figure enhancement is strongest in feedback-recipient layers (1, 2, 5). This is the characterization in the follow-up Poort et al. 2016.
  - [VU record](https://research.vu.nl/en/publications/distinct-roles-of-the-cortical-layers-of-area-v1-in-figure-ground/); [Poort et al. 2016, PMC5028009](https://pmc.ncbi.nlm.nih.gov/articles/PMC5028009)

**Timing of suppression relative to excitation**
- Bair, Cavanaugh & Movshon (2003, *J Neurosci* 23:7690–7701), anesth., 92 complex cells:
  - CRF onset 52 ± 13 ms; offset 35 ± 10 ms.
  - Surround-suppression onset 61 ± 17 ms, a mean delay of 9 ± 15 ms after CRF onset. Strong suppression often arrives at CRF onset, and weak suppression is about 30 ms later (r = −0.48).
  - Suppression is mostly sustained for near surrounds and mostly transient for far surrounds.
  - The time-distance relationship is weak. The implied speed is often >1 m/s (about 40% of cells), and only 1/31 cells were <0.2 m/s.
  - Conclusion: feedback, plus a fast local source.
  - [PMC6740744](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744)
- Smith, Bair & Movshon (2006, *J Neurosci* 26:4826–4834), anesth., 112 complex cells:
  - Cross-orientation suppression latency is 42.5 ± 11.1 ms, versus excitation onset of 50.0 ± 15.8 ms. Suppression thus starts about 6–7 ms before response onset.
  - Surround suppression comes about 12 ms later than cross-orientation suppression within the same cells.
  - The two strengths are uncorrelated (r = 0.11).
  - Conclusion: cross-orientation suppression is feedforward (LGN nonlinearity) or fast local inhibition; surround suppression is more feedback-like. Monocular viewing only.
  - [PMC6674161](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161)
- Henry, Jazayeri, Shapley & Hawken (2020, *eLife* 9:e54264, doi:10.7554/eLife.54264), anesth., 106 neurons, 15 macaques, histological layer assignment, reverse correlation of eCRF:
  - Untuned suppression peaks at 55 ± 2 ms (55% of cells).
  - Tuned facilitation peaks at 60 ± 3 ms (25%).
  - Tuned suppression peaks at 76 ± 3 ms (79%).
  - By layer: 4Cα has tuned suppression peaking at 60–75 ms with untuned suppression present. 4Cβ has tuned suppression peaking at about 100 ms with no average untuned suppression. L2/3 is dominated by tuned suppression at about 100 ms. 4B and L6 are 4Cα-like (magno). L5 is mixed (n = 8).
  - [PMC7253173](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)
- Wang et al. (2020, *J Neurosci*), awake MUA:
  - Fast suppression (FS) delay relative to excitation: 4Cα 11.4, L6 14.9, 4B 18.1, L5 18.3, L2/3 20.5, 4Cβ 21.3 ms.
  - Slow suppression (SS) delay: L2/3 44.7, 4B 46.3, L5 46.8, L6 58.8, 4Cα 60.3, 4Cβ 61.0 ms.
  - FS is strongest in input layers (4C, 6), is attributed to GABA-A and is local.
  - SS is about 3× stronger in output layers (2/3, 4B, 5) and is attributed to long-range or feedback input, possibly with a GABA-B contribution.
  - Gratings were at 90% contrast only.
  - [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183)
- Bijanzadeh et al. (2018), MUA suppression onset relative to the earliest contact:
  - Full-field surround suppression appears first in L4C (4 ± 1.2 ms), coincident with response onset, so it is inherited from the LGN.
  - Near-surround suppression appears first in upper-L4 (8.8 ± 2 ms), then Deep-U (19.6 ms), L2/3 (25.5 ms) and L4C (34.2 ms). This is intracortical, about 14 ms after the L4C RF response.
  - Far-surround suppression appears first in the deep layers. In absolute terms this is about 16–19 ms after the L4C RF response. It is weak or absent in L4C (42.5 ± 12.7 ms where present), consistent with feedback.
  - Orientation-tuned near suppression is already tuned when it first emerges.
  - [author manuscript](https://par.nsf.gov/biblio/10110109)

**Timing of binocular facilitation and interocular suppression across layers**
- Cox, Dougherty, Westerberg, Schall & Maier (2019, *J Vis* 19(12):13, doi:10.1167/19.12.13), 2 awake *M. radiata*, 30 penetrations, 24/32-channel, 100 µm probes, 152 SUs, dominant eye about 0.45 vs non-dominant about 0.9.
  - Transient window (50–100 ms): binocular > monocular for congruent (+13.7 sp/s) and incongruent (+5.9 sp/s) stimuli.
  - Sustained window (150–250 ms): congruent facilitation disappears, and incongruent stimuli suppress (−5.6 sp/s).
  - In the sustained window, stronger ocular bias predicts more suppression (congruent r = −0.40).
  - By layer:
    - 4C equiocular neurons keep congruent facilitation into the sustained window.
    - Upper- and deep-layer equiocular neurons facilitate only transiently.
    - Deep layers show unexpected early (transient) suppression for incongruent stimuli.
  - CSD: the granular sink is larger for both binocular conditions at 50–100 ms. At 150–250 ms the congruent sink returns to monocular levels, while the incongruent sink stays larger.
  - [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- Dougherty, Cox, Westerberg & Maier (2019, *Curr Biol*), laminar arrays:
  - Monocular neurons concentrate in L4C, while binocular neurons dominate the layers above and below.
  - Monocular neurons are nonetheless modulated by binocular stimulation, more often suppressed than facilitated (as characterized by Zhang et al. 2024).
  - Spike analysis window was about 50–200 ms.
  - [bioRxiv 320218](https://www.biorxiv.org/content/10.1101/320218v2.full); [eLife 92839 reviews](https://elifesciences.org/reviewed-preprints/92839v2/reviews)
- Mitchell et al. (2022): early facilitation (50–150 ms) scales with contrast, and late suppression (150–250 ms) appears at medium and high contrast. [PubMed 35494250](https://pubmed.ncbi.nlm.nih.gov/35494250/)
- Carlson, Mitchell, Dougherty, Westerberg, Cox & Maier (2023, *iScience* 26:107359, doi:10.1016/j.isci.2023.107359) asked whether V1 response suppression initiates binocular rivalry. [PubMed 37520732](https://pubmed.ncbi.nlm.nih.gov/37520732/)
  - The departmental summary of Carlson's dissertation says that at the onset of binocular rivalry flash suppression (BRFS), response facilitation, not suppression, dominates.
  - The same summary says interocular conflict reduces interlaminar coherence, especially between deep and superficial layers. [Vanderbilt news](https://as.vanderbilt.edu/psychology/2025/03/06/brock-carlson-successfully-defends-dissertation/)
- Precortical interocular suppression (Dougherty et al. 2021, *eNeuro* 8(2), doi:10.1523/ENEURO.0364-20.2020), awake single units in LGN:
  - Monocular LGN neurons are suppressed by high contrast in the other eye.
  - Suppression is smaller with high dominant-eye contrast or low non-dominant-eye contrast.
  - [eNeuro](https://www.eneuro.org/content/8/2/ENEURO.0364-20.2020)
- Sengpiel & Vorobyov (2005, *J Neurosci* 25:6394–6400, cat, anesth.):
  - Interocular suppression cuts off at 8–16 Hz temporal frequency, whereas monocular cross-orientation suppression extends to 32 Hz (LGN-like).
  - Bicuculline reduced interocular suppression from 46.1 ± 4.0% to 18.1 ± 4.9%, implicating intracortical inhibition.
  - These numbers come from the authors' conference abstract. [physoc abstract](https://www.physoc.org/abstracts/interocular-cross-orientation-suppression-in-cat-striate-cortex-involves-intracortical-inhibition/); [Cardiff ORCA](https://orca.cardiff.ac.uk/id/eprint/1153/)
- Dichoptic masking: monoptic masking is found in all LGN and V1 neurons, but dichoptic masking only in a subset of binocular V1 neurons. V1 binocular neurons are "dichoptic for excitation but monoptic for inhibition", i.e. interocular suppression is weak in early binocular neurons. My recollection is that this is Macknik & Martinez-Conde, monkey. [SUNY record](https://researchconnect.suny.edu/en/publications/dichoptic-visual-masking-reveals-that-early-binocular-neurons-exh/)
- Continuous flash suppression, two-photon (Chen, Wang, Jiang, Tang & Yu, 2025/2026, *eLife*, doi:10.7554/eLife.107518), awake, superficial V1:
  - Orientation-tuning amplitude falls 84% and 61% in the two monkeys.
  - Tuning is abolished in masker-eye neurons, and grating-eye neurons fall 78% and 42%.
  - Fit by an OD-dependent gain-control model (R² 0.98–0.99).
  - Calcium kinetics cannot resolve suppression onset.
  - [eLife RP 107518](https://elifesciences.org/reviewed-preprints/107518v3)

### Inferences
- **What the latency differences require.** The differences that discriminate mechanisms are about 5–60 ms:
  - fast vs slow suppression, about 10 vs 40 ms;
  - 4Cα vs 4Cβ, about 20 ms;
  - L4C vs L2/3, about 15–25 ms;
  - near vs far surround, about 5–15 ms.
  
  Photodiode-locked 1 ms timing is more than sufficient, provided the 3.4 ms NSx offset is corrected. Fine onset differences of 2–5 ms (e.g., Smith 2006's −6 ms) need many trials. Analog MUA averaged across trials per contact is the right signal for onset.
- **Flashed vs drifting stimuli.** The team's flashed static gratings are closer to Bair/Smith/Albrecht and Wang 2020 (flashed) than to drifting-grating studies. Onset transients in flashed stimuli exaggerate the early facilitation window, which is consistent with the Maier-lab transient facilitation findings.
- **Open question: is interocular suppression a fast, input-layer, LGN-inherited component or a slow, output-layer, intracortical or feedback component? Or both, as for surround suppression?**
  - Discriminating measurement: per layer, the latency at which the dichoptic (fixed-eye contrast + other-eye contrast) response diverges from the monocular response, as a function of other-eye contrast. Compare against excitation onset in the same layer.
    - LGN-inherited suppression should appear in 4C at response onset (Bijanzadeh logic) and scale with other-eye contrast per Dougherty 2021.
    - Intracortical suppression should appear first outside 4C, or in 4C with a delay of 10 ms or more.
  - Design fit: **Yes.** This is directly computable from the 6×6 matrix with analog MUA and CSD per contact. The 317 ms duration covers both the transient and the 150–250 ms sustained windows used by Cox 2019 and Mitchell 2022.
- **Open question: does the latency of binocular facilitation or suppression shift with contrast in the same way as excitation?**
  - Discriminating measurement: compare the contrast-dependent onset shift of the monocular response with the onset of the binocular-minus-monocular difference.
  - Design fit: **Yes,** with 5 nonzero contrasts per eye.
- **Open question: why do deep layers show early suppression for incongruent stimuli (Cox 2019)?**
  - Candidate explanations are direct LGN input to L6 or propagation from upper layers.
  - Discriminating measurement: the relative onset of suppression in L6 vs L2/3 vs 4C, and whether L6 suppression precedes the 4C sink.
  - Design fit: **Partial.** The team's stimuli are dichoptic contrast with the same orientation, not orthogonal conflict, so the incongruent-orientation condition is absent. Contrast imbalance does not substitute for orientation conflict.

### Gaps
- No primary-source ms values were found for the contrast-dependent latency shift of macaque V1 (Albrecht 2002, Gawne 1996 and Bair 2002 were not retrieved). Only the "scaling and shifting" description is verified.
- Exact Schroeder 1998 laminar ms values and Self 2013 per-layer latencies were not retrieved.
- Most precise laminar suppression timing (Bijanzadeh, Henry, Bair, Smith) is from anesthetized animals. Only Wang 2020, Gieselmann & Thiele 2022 and the Maier-lab binocular studies are awake.
- I found no published per-layer *onset latency* of interocular suppression in macaque V1. Cox 2019 used fixed windows and reported no onset latency for facilitation or suppression. This looks like a genuine gap the team could fill.

---

## 3. Mechanisms of normalization: feedforward/LGN vs recurrent intracortical vs feedback; inhibitory cell types; laminar CSD under high vs low contrast

### Takeaway
Primate laminar evidence supports a multi-source normalization.
- **Feedforward/LGN-inherited:** at response onset, in L4C. Covers cross-orientation suppression and full-field surround suppression, and part of interocular suppression, since it is present in LGN.
- **Local intracortical:** about 10 ms after onset, strongest in input layers. Covers near-surround suppression and fast suppression.
- **Feedback or long-range horizontal:** 15–60 ms after onset, in feedback-recipient layers 1 and 5/6 and the output layers. Covers far-surround and slow suppression.

Primate cell-type causal evidence is thin. One macaque study shows normalization persisting under optogenetic drive of excitatory cells. A marmoset preprint shows PV cells give linear gain control in granular layers and tuned suppression in extragranular layers.

### Cited Findings
**Optogenetics and cell types**
- Nassi, Avery, Cetin, Roe & Reynolds (2015, *Neuron* 86:1504–1517, doi:10.1016/j.neuron.2015.05.040), alert macaque V1, optogenetic depolarization of excitatory neurons:
  - Stimulation both raised and lowered baseline activity, consistent with indirect recruitment of inhibition.
  - Combined optogenetic and visual responses were sub-additive and depended lawfully on stimulation intensity and luminance contrast.
  - "Normalization computation persists" under artificial drive.
  - [PMC4534089](https://pmc.ncbi.nlm.nih.gov/articles/PMC4534089)
  - The follow-up model (Avery, Nassi & Reynolds 2018, *PLoS One*) attributes the suppression to slow excitatory and inhibitory conductances. [PMC6203252](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6203252/)
- Marmoset V1 PV+ (Vafaei, Clark, Khadir, Dai, Federer & Angelucci, bioRxiv Dec 2025, doi:10.64898/2025.12.23.696300), bidirectional optogenetics plus laminar recordings. Not peer reviewed.
  - In granular layers, PV+ cells implement divisive/multiplicative *linear* gain control.
  - In extragranular layers, they provide tuned nonlinear suppression that sharpens orientation tuning.
  - Models attribute gain control to threshold shifts and tuning to input noise.
  - [bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.23.696300v1)
- Macaque GABAergic targeting tools exist but have not been used for normalization:
  - Dlx-enhancer optogenetic inactivation in macaque V1, with 468/543 mCherry+ neurons PV+ (De, El-Shamayleh & Horwitz 2020, *eLife*). [eLife 52658](https://elifesciences.org/articles/52658)
  - Enhancer-driven PV/SST access with Neuropixels in development at Salk (FENS 2026 abstract). [FENS abstract](https://world-wide.org/fens-26/cell-type-specific-optogenetic-access-ec6945e1)

**Laminar origin of suppression components**
- Bijanzadeh et al. 2018: full-field surround suppression in L4C coincides with response onset, so it is LGN-inherited. Near-surround suppression is intracortical, outside L4C. Far-surround suppression arises from feedback via L1 and deep-lower layers. [author manuscript](https://par.nsf.gov/biblio/10110109)
- Henry et al. 2020 on untuned suppression:
  - It likely has multiple sources: LGN magno input plus cortical lateral inhibition.
  - The delayed component in parvo-dominated neurons is probably cortical.
  - Whether eCRF modulation is "inherited from the LGN versus generated within the cortical microcircuit" remains open.
  - [PMC7253173](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)
- Wang et al. 2020: two laminar subnetworks.
  - Input network (4C, 6) has fast GABA-A-like suppression.
  - Output network (2/3, 4B, 5) has slow suppression about 3× stronger, attributed to long-range or feedback input.
  - [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183)
- Smith, Bair & Movshon 2006: cross-orientation suppression is too fast for extrastriate feedback, and synaptic depression is unlikely given its slow recovery. They favour LGN nonlinearities or fast local interneurons. [PMC6674161](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161)

**Population (VSDI) dynamics**
- Reynaud, Masson & Chavane (2012, *J Neurosci* 32:12558–12569), VSDI in 2 behaving macaques:
  - Generalized normalization explains the dynamic CRF through recurrent polysynaptic intracortical loops. The recurrent pool "swells" at 0.007–0.014 m/s.
  - It explains surround suppression through monosynaptic horizontal spread at 0.1–0.3 m/s.
  - [PMC6621242](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242)
- Sit, Chen, Geisler, Miikkulainen & Seidemann (2009, *Neuron* 64:943–956), VSDI of superficial V1 over about 1 cm² in fixating macaques with flashed Gabors: complex dynamics are explained by a *feedforward* population gain-control model. [PMC2807412](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412)

**Interocular mechanisms**
- Interocular suppression is partly precortical, present in LGN (Dougherty et al. 2021). [eNeuro](https://www.eneuro.org/content/8/2/ENEURO.0364-20.2020)
- In cat, interocular suppression is GABA-A sensitive and low-pass in temporal frequency, i.e. cortical (Sengpiel & Vorobyov 2005). [physoc abstract](https://www.physoc.org/abstracts/interocular-cross-orientation-suppression-in-cat-striate-cortex-involves-intracortical-inhibition/)

**CSD and LFP under binocular, size and adaptation manipulations (laminar CSD as input measure)**
- Cox et al. 2019: the granular sink grows with binocular input in the transient. It then returns to monocular level for congruent input, i.e. "normalizes", but not for incongruent input. [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- Westerberg, Cox, Dougherty & Maier (2019, *J Neurophysiol* 121:1938–1952): repetition adaptation transfers across the eyes and is localized mainly to supragranular intracortical processing. Current-sink changes were compared with spiking output. [PMC6589708](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708)
- Xing et al. 2012 (anesth., >90% contrast): sustained MUA is constant across depth (P = 0.52), while sustained γ-LFP is concentrated in 2/3 and 4B (67% of γ). [PMC3427063](https://pmc.ncbi.nlm.nih.gov/articles/PMC3427063)

### Inferences
- **What CSD vs spikes can separate.** Normalization that is already present in the granular current sink implies reduced feedforward drive (LGN or thalamocortical). Normalization that appears only in spikes, or only in extragranular sinks, implies intracortical or feedback inhibition. The Cox-2019 CSD approach applied across the full contrast matrix is a direct, unexploited discriminator.
- **Open question: is interocular normalization implemented as reduced input (smaller 4C sink at high other-eye contrast) or as output inhibition (unchanged sink, suppressed spikes)?**
  - Discriminating measurement: compare L4C sink amplitude and onset against L4C MUA as a function of other-eye contrast, and do the same in L2/3 and L5/6.
  - Design fit: **Yes,** for CSD and MUA. Monocular L4C units are the key population.
- **Open question: does the PV-like "linear gain in granular, tuned suppression extragranular" split (marmoset) predict layer-specific contrast-gain vs response-gain signatures of normalization in macaque?**
  - Discriminating measurement: CRF parameter shifts by layer (see §1). Optionally separate putative fast-spiking (narrow-waveform) from regular-spiking units.
  - Design fit: **Partial.** Correlational only, with no cell-type manipulation. About 70 SUs is small for an FS/RS-by-layer split.
- **Open question: does feedback contribute to interocular normalization?**
  - Discriminating measurement: the earliest suppression appearing in L1 and deep-lower sinks (Bijanzadeh signature) vs in 4C.
  - Design fit: **Partial.** The team's gratings are RF-centred with no surround manipulation, so the feedback signature would have to come from timing and laminar sink pattern alone.

### Gaps
- I found no primary study comparing laminar CSD sink amplitude or timing for **high vs low luminance contrast** in macaque V1. Tootell 1988 2-DG is the closest laminar contrast-dependence evidence. Treat this as a gap the team could fill.
- I found no peer-reviewed macaque study isolating PV, SOM or VIP contributions to contrast normalization (searched 2023–2026).
- The "Li, Reynolds et al." optogenetic normalization paper named in the brief was not located.

---

## 4. Inter-columnar interactions: travelling waves, propagation speed, contrast dependence, dual-probe/array latency differences, ocular-dominance relationships

### Takeaway
Stimulus-evoked and spontaneous activity in primate visual cortex propagates horizontally at about 0.1–0.6 m/s, the unmyelinated horizontal-axon range. VSDI single-trial estimates run higher (median 0.57 m/s).

Lateral coupling *shrinks* with contrast. The spike-triggered-LFP space constant falls about 3-fold from blank to 100% contrast in macaque, and summation fields are 2.3× larger at low contrast. Surround suppression travels inward at about 0.1–0.3 m/s (VSDI) but behaves as if >1 m/s in single units, which points to feedback.

No study I found measured contrast-dependent onset-latency differences between two *laminar* probes in separate V1 columns, or related cross-column coupling to ocular-dominance column identity. This is a gap the 7 dual-probe days could partly address.

### Cited Findings
**Propagation speeds and contrast dependence**
- Muller, Reynaud, Chavane & Destexhe (2014, *Nat Commun* 5:3675, doi:10.1038/ncomms4675), VSDI in 2 awake macaques, 0.5°-s.d. Gaussian blob, 50 ms:
  - The evoked response propagates as a travelling wave on single trials.
  - Speeds were 0.25–1.35 m/s (median 0.57, MAD 0.18 m/s).
  - Waves were detected on 80% of 50 ms trials, 52.6% of 100 ms trials and 20% of 600 ms trials.
  - Spread was an ellipse of about 7.5 × 3.7 mm (about 2° × 2.7°).
  - Attributed to the superficial horizontal fibre network. Contrast dependence was not tested.
  - [PMC4015334](https://pmc.ncbi.nlm.nih.gov/articles/PMC4015334); [Nature](https://www.nature.com/articles/ncomms4675)
- Nauhaus, Busse, Carandini & Ringach (2009, *Nat Neurosci* 12:70–76), anesth., Utah 10×10 at 400 µm (macaque V1 and cat area 18):
  - Spike-triggered LFP waves travel at 0.24 ± 0.20 m/s (macaque) and 0.31 ± 0.23 m/s (cat).
  - The spontaneous space constant λ is 2.1 ± 0.6 mm (macaque) and 5.8 ± 1.9 mm (cat).
  - From blank to 100% contrast, λ falls by 3.0× (macaque) and 5.5× (cat), and amplitude falls by 1.7× and 2.7×. 25% contrast is intermediate.
  - Coupling is stronger between similar-orientation sites.
  - Conclusion: "lateral connections dominate at low contrast."
  - Recordings were likely superficial; laminar depth was not confirmed.
  - [PMC2610236](https://pmc.ncbi.nlm.nih.gov/articles/PMC2610236)
- Caveat (Ray & Maunsell 2011, *J Neurosci* 31:12674–12682): the apparent outward-travelling STA wave is a superposition of two parts.
  - One is a brief spike-locked negativity that decays quickly with distance.
  - The other is a slow network rhythm peaking about 25 ms after the spike that decays more gradually.
  - Network activity must therefore be accounted for before reading STAs as connectivity.
  - [PMC3488382](https://pmc.ncbi.nlm.nih.gov/articles/PMC3488382)
- Reynaud et al. 2012 (VSDI, awake): surround suppression spreads inward at 0.1–0.3 m/s. [PMC6621242](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242)
- Chemla et al. (2019, *J Neurosci* 39:4282–4298), VSDI in awake macaque with long-range apparent motion: *suppressive* travelling waves produce a spatiotemporal normalization of the stimulus representation. The model generates them from recurrent intracortical gain control. [PMC6538863](https://pmc.ncbi.nlm.nih.gov/articles/PMC6538863)
- Bair et al. 2003: horizontal connections are cited as about 0.1–0.2 m/s with a 3.5–4.5 mm reach. Surround suppression reached a median of 7.5 mm (range 2–15 mm) at an implied speed of about 1 m/s, which favours feedback. [PMC6740744](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744)
- Zhu, Oh, Trepka, Chen & Moore (2026, *eLife* 13:RP103255), anesth., Neuropixels, 4 penetrations, 677 neurons:
  - They cite horizontal conduction at 0.1–0.4 mm/ms vs feedback at 20–60 mm/ms. These are literature values.
  - Under border-ownership (nCRF) stimulation, layers 2/3 and 5/6 lead the input layer more than under CRF stimulation.
  - [PMC12774416](https://pmc.ncbi.nlm.nih.gov/articles/PMC12774416/)
- Grinvald, Lieke, Frostig & Hildesheim (1994, *J Neurosci* 14:2545–2568), VSDI in macaque V1:
  - A listing reports spread at 100–250 µm/ms (0.1–0.25 m/s), with space constants of about 2.7 mm parallel and 1.5 mm perpendicular to the V1/V2 border.
  - The listing's journal label conflicts with the paper, so treat these numbers as unverified.
  - [IBM Research listing](https://research.ibm.com/publications/cortical-point-spread-function-and-long-range-lateral-interactions-revealed-by-real-time-optical-imaging-of-macaque-monkey-primary-visual-cortex)
- Benucci, Frazor & Carandini (2007, *Neuron* 55:103–117), **cat** V1 VSDI: activity is a *standing* wave in orientation space and a *travelling* wave in retinotopic space at 0.2–0.5 m/s. [PMC2171365](https://pmc.ncbi.nlm.nih.gov/articles/PMC2171365)
- Sato, Nauhaus & Carandini (2012, *Neuron* 75:218) review V1 travelling waves and their likely horizontal-connection substrate. [Monash record](https://research.monash.edu/en/publications/traveling-waves-in-visual-cortex/)
- Davis, Muller, Martinez-Trujillo, Sejnowski & Reynolds (2020, *Nature* 587:432–436), **marmoset MT**, Utah arrays:
  - Spontaneous travelling waves at 0.1–0.6 m/s, consistent with unmyelinated horizontal axons.
  - Their phase at target onset predicts target-evoked response and detection.
  - [PMC7677221](https://pmc.ncbi.nlm.nih.gov/articles/PMC7677221)
- Davis et al. (2021, *Nat Commun*): spontaneous waves emerge from horizontal-fibre time delays and travel through locally asynchronous-irregular states. [DOAJ](https://doaj.org/article/03f4f65aa15d40559abd129057a7ae9b)

**Low-contrast spatial integration**
- Sceniak et al. 1999: summation is 2.3× larger at low contrast. [Nature](https://www.nature.com/articles/nn0899_733)

**Ocular-dominance relationships**
- Model of 4Cα plus L6 (Xiao, Lin & Young 2026):
  - Binocular strips at ODC borders emerge in the model.
  - The authors infer that about 10–30% of connections from border neurons cross into the opposite-eye column.
  - L6 projections cross ODC borders freely (citing Wiser & Callaway 1996).
  - They predict L6 feedback to L4 is largely monocular.
  - [arXiv](https://arxiv.org/html/2606.21785)
- Yoshioka, Blasdel, Levitt & Lund (1996, *Cereb Cortex* 6:297–310) examined lateral connectivity relative to ODCs and CO blobs. Only the citation was retrieved; the findings were not verified. [ordb record](https://ordb.biotech.ttu.edu/odordb/Data/77760)

**Laminar inter-column context (recent)**
- Xu, Morton, Hudson, Nandy & Jadi (2024 preprint): the spatial configuration of flankers differentially modulates inter-laminar communication efficacy and the feedforward/feedback balance. [OUCI record](https://ouci.dntb.gov.ua/works/7PXLDM07)

### Inferences
- **Expected inter-probe lags.** If the two probes are d mm apart on cortex (d estimated from RF separation × magnification factor), horizontal propagation at 0.1–0.6 m/s predicts onset or peak lags of d/0.6 to d/0.1 ms.
  - At d = 2 mm this is about 3–20 ms; at d = 5 mm, about 8–50 ms.
  - Feedback at about 2–6 m/s predicts lags of only about 1–2.5 ms.
  - These ranges are separable with sample-locked dual-probe data and 1 ms timing, *if* d is known and RFs do not overlap (overlap would add direct feedforward drive).
- **Contrast prediction (Nauhaus).** Cross-column coupling should be strongest at low contrast and in spontaneous or blank periods: larger amplitude and longer reach. Evoked cross-column latency differences may also depend on contrast, because low-contrast RFs expand (Sceniak) and recruit lateral input.
- **Open question: is the influence of one column on another carried by horizontal connections or by feedback, and does contrast change the balance?**
  - Discriminating measurements:
    - the lag of the cross-probe spike-triggered CSD or LFP (or CCG peak) and its laminar locus — horizontal targets L2/3 and 5, feedback targets L1 and 5/6;
    - whether the lag scales with inter-probe distance;
    - how coupling amplitude and reach change with contrast (0 / 0.055 … 0.9).
  - Ray & Maunsell require controlling for γ/network rhythms, e.g. by trial shuffling or spike-time jitter.
  - Design fit: **Partial.** 7 probe pairs give 7 distances, a small N for a speed estimate. The same-file 72-channel recordings remove clock alignment error. The 0-contrast cells of the 6×6 matrix provide within-session spontaneous baselines.
- **Open question: does cross-column coupling depend on whether the two columns share ocular dominance?** (Relevant to whether interocular normalization pools across ODCs horizontally.)
  - Discriminating measurement: classify each probe's column by monocular dominance, then compare coupling for same-eye vs opposite-eye column pairs, and for stimulation of the eye dominant at the source vs the target column.
  - Design fit: **Weak.** With 7 pairs, the number of same-eye vs opposite-eye pairs will be tiny. Report as exploratory or case-study.
- **Open question: do evoked responses propagate across columns as waves on single trials (Muller 2014), and is the wave suppressive at high contrast (Chemla 2019)?**
  - Design fit: **Limited.** Two sampling points cannot define a 2-D wave. Only a pairwise lag consistent with a wave is testable.

### Gaps
- I found no study with dual *laminar* probes in two V1 columns reporting inter-column onset-latency differences by layer, or their contrast dependence.
- I found no primary-source quantitative link between horizontal connectivity and ODC identity in macaque. The Yoshioka 1996 findings were not retrieved.
- Zanos et al. 2015 (saccade-related waves, *Neuron*) was not retrieved.
- Kajikawa & Schroeder 2011 ("How local is the LFP?") and the Seidemann-lab VSDI–spiking relationship (Chen et al. 2012) were not retrieved.

---

## 5. Spike–LFP phase coupling and spectral/CSD markers for laminar alignment, and their caveats

### Takeaway
Three alignment markers are in use:
- the early CSD sink, the classic "gold standard", which in practice is identifiable in only about half of probes;
- the spectrolaminar α/β–γ crossover (Mendoza-Halliday et al. 2024);
- spike–LFP phase reversal (Davis et al. 2023).

Each has documented failure rates and biases. V1 specifically is problematic for the α/β marker because of volume conduction from V2. Davis et al. 2023 contains no V1 data. Best practice is to require agreement of at least two markers per penetration.

### Cited Findings
**Spike–LFP phase coupling (Davis et al. 2023)**
- Davis, Dotson, Franken, Muller & Reynolds (2023, *eLife* 12:e84512) recorded V4 (2 macaques), MT and PFC (marmosets) with 32-channel linear arrays. **No V1.**
- Method:
  - Spiking phase-locks near ±π to superficial and input-layer LFP and near 0 to deep-layer LFP, using 5–50 Hz generalized phase.
  - The reversal marks the input-layer/deep boundary, with the input layer defined as the earliest CSD sink.
- Agreement with CSD:
  - r = 0.64 with CSD depth (34 sessions).
  - Reported mean difference 271 ± 49 µm (SEM) with no significant difference between estimates (p = 0.98).
  - 14/34 sessions differed by ≥200 µm; in those, the phase-coupling estimate better matched firing-rate and LFP-power patterns.
- Practical limits:
  - Works with about 1 min of data.
  - Estimates drifted about 132 µm deeper over the first minutes as tissue settled.
  - The reversal channel was chosen by eye, and only sessions with a recognizable CSD were included.
  - Mechanism untested; generality to anesthesia or other areas unknown.
- Sources: [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800); [eLife](https://elifesciences.org/articles/84512)

**Spectrolaminar motif (Mendoza-Halliday et al. 2024)**
- Mendoza-Halliday, Major, Lee, Lichtenfeld, Carlson, Mitchell, Meng, Xiong, Westerberg, Jia, Johnston, Selvanayagam, Everling, Maier, Desimone, Miller & Bastos (2024, *Nat Neurosci* 27:547–560, doi:10.1038/s41593-023-01554-7). The team's own lab members are co-authors.
- Scale: 5 macaques, 14 areas, 810–942 probes.
- Motif: γ peaks in L3 and α/β in L6, with the crossover at L4. Crossover distance from L4 across areas is 46 ± 51 µm (SEM).
- Identification rates:
  - Spectrolaminar pattern: 61% manual, 64% FLIP, 81% vFLIP.
  - Early CSD sink: only 51% of 810 probes.
- **V1 exception:** the highest α/β power in V1 was in white-matter channels, likely volume-conducted from V2. V1 data came from one bonnet macaque (10 sessions), and the V1 histology n is inconsistent (3 vs 5).
- Source: [PMC10917659](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659)

**Critique (Mackey et al.)**
- Mackey, Duecker, Neymotin, Dura-Bernal, Haegens, Barczak, O'Connell, Jones, Ding, Ghuman & Schroeder (bioRxiv 2024, v4 Dec 2025, doi:10.1101/2024.09.18.613490), using A1, belt cortex and V1:
  - The gradient was present in only 61–64% of recordings ("common but by no means universal").
  - The crossover identified L4 in only 29–33%.
  - FLIP/vFLIP are variable and spuriously detect inversions.
  - Cited biophysical models do not reproduce the LFP trends.
  - They call for community best practices on L4 identification.
- Source: [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.18.613490v4)

**LFP spatial reach**
- Xing, Yeh & Shapley (2009, *J Neurosci* 29:11540–11549): LFP integrates over about 250 µm radius in macaque V1, with laminar variation. [PMC2910581](https://pmc.ncbi.nlm.nih.gov/articles/PMC2910581)
- A 2016 *J Neurophysiol* study reports that LFP spread in V1 is band-pass, largest in high-γ (60–150 Hz). [PMC5082390](https://pmc.ncbi.nlm.nih.gov/articles/PMC5082390)

**Alignment practice in related laminar V1 studies**
- Bijanzadeh 2018 and Nurminen 2026: earliest CSD sink plus histology.
- Wang 2020: earliest sink = mid-4Cα at relative depth 0.49.
- Zhu 2026: bottom of the major sink = 4C/5 border, plus DiI histology.
- Cox 2019: ±50 µm alignment error at 100 µm spacing.
- Sources: [par.nsf.gov](https://par.nsf.gov/biblio/10110109); [PMC12931924](https://pmc.ncbi.nlm.nih.gov/articles/PMC12931924/); [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183); [PMC12774416](https://pmc.ncbi.nlm.nih.gov/articles/PMC12774416/); [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)

### Inferences
- **Which markers to trust in V1.** The CSD early sink (driven by high-contrast gratings) is the most V1-appropriate marker. The α/β deep peak is suspect in V1 (V2 volume conduction). The phase-reversal method is unvalidated in V1, but is cheap to compute from inter-trial intervals.
- **Open question: do the three markers agree in awake macaque V1 at 100 µm spacing, and does marker choice change laminar conclusions?**
  - Discriminating measurement: compute all three per penetration, quantify pairwise offsets, and re-run the key laminar analyses (CRF parameters, suppression latencies) under each alignment.
  - Design fit: **Yes.** About 30 penetrations with CSD plus spontaneous periods give a first V1 validation of the Davis phase-reversal marker. This is of independent methodological value.
- **Contrast as a lever on CSD reliability.** The team's contrast matrix also enables a CSD-sink-onset-vs-contrast analysis. If the earliest sink shifts depth with contrast (e.g., 4Cα dominant at low contrast per Tootell 1988), that is both a laminar marker for 4Cα vs 4Cβ and a confound for alignment that uses a single contrast.

### Gaps
- No V1 validation of the Davis et al. phase-reversal marker was found.
- Whether the laminar position of the CSD early sink depends on stimulus contrast in macaque V1 was not found in the literature.

---

## 6. Explicitly stated open questions in recent papers (with design-fit notes)

### Takeaway
Recent laminar papers repeatedly leave four questions open:
- LGN-inherited vs cortically generated suppression;
- whether L2/3 selectivity is inherited from L4C or generated recurrently;
- the circuit origin of early deep-layer binocular suppression;
- how to identify L4 reliably.

Binocular two-photon papers explicitly list *unequal interocular contrast* and *fast dynamics* as future work, which is the team's design niche.

### Cited Findings
- Henry et al. 2020: unclear whether eCRF modulation arises from one or several mechanisms, and "to what extent it is inherited from the LGN versus generated within the cortical microcircuit". Layer 5 orthogonal facilitation and 4B tuned suppression rest on small samples. [PMC7253173](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)
- Wang et al. 2020 open questions:
  - whether L2/3 selectivity inherits fast suppression from L4C or recurrent excitation dominates;
  - the origin of L2/3 excitation (4C feedforward vs L2/3 recurrent vs L2/3–L5 loop);
  - why weak slow suppression in input layers does not track drifting-grating surround suppression.
  - [PMC7511183](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183)
- Cox et al. 2019:
  - The mechanism of early deep-layer suppression is unresolved: direct retinogeniculate input to L6 vs propagation from upper layers.
  - Whether sustained incongruent suppression initiates rivalry is untested.
  - Results are limited to zero disparity; small per-layer samples; 100 µm spacing (Neuropixels suggested).
  - [PMC6797477](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)
- Zhang et al. 2024 (two-photon):
  - Uneven interocular contrast is future work.
  - Calcium "cannot reveal the fast dynamics", which may explain discrepancies with electrode data (Dougherty 2019).
  - Two-photon imaging samples mainly L2/3.
  - [eLife RP 92839](https://elifesciences.org/reviewed-preprints/92839v2); [reviews](https://elifesciences.org/reviewed-preprints/92839v2/reviews)
- Zhu et al. 2026: how much L4 input comes from feedback/horizontal laminae is unknown; feedback origin of border ownership is untested; the role of anesthesia is unresolved. [PMC12774416](https://pmc.ncbi.nlm.nih.gov/articles/PMC12774416/)
- Mendoza-Halliday et al. 2024: the generating mechanisms of the γ and α/β gradients are unknown. [PMC10917659](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659)
- Mackey et al.: a reliable (online) L4 identification method is needed, as are community best practices. [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.18.613490v4)
- Davis et al. 2023: generality to non-sensory areas, anesthesia and sleep is unknown, and the dipole mechanism is untested. [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Xiao, Lin & Young 2026: explicit testable predictions are 10–30% cross-ODC connections from border neurons, largely monocular L6→L4 feedback, and L6 only slightly more binocular than L4. [arXiv](https://arxiv.org/html/2606.21785)
- Nurminen et al. 2026: variability is size-tuned in a layer-dependent way.
  - The Fano factor drop for an RF-sized stimulus is −34% in SG, −24% in G and −14% in IG.
  - A far surround *raises* SG Fano factor (+17%).
  - So "surround suppression does not always quench variability". Contrast effects on laminar variability were not tested.
  - [PMC12931924](https://pmc.ncbi.nlm.nih.gov/articles/PMC12931924/)
- Carlson dissertation summary: facilitation, not suppression, dominates BRFS onset, and conflict reduces deep–superficial coherence. [Vanderbilt news](https://as.vanderbilt.edu/psychology/2025/03/06/brock-carlson-successfully-defends-dissertation/)

### Inferences
Each item gives the open question, the discriminating measurement, and the design fit.

1. **LGN-inherited vs cortical interocular normalization.**
   - Measurement: suppression onset in 4C relative to response onset, and 4C sink amplitude vs other-eye contrast. Compare with the LGN interocular-suppression contrast dependence (Dougherty 2021).
   - Design fit: **Yes.**
2. **Whether L2/3 interocular normalization is inherited from L4C or generated recurrently.**
   - Measurement: the laminar order of suppression onset, plus CRF parameter changes from 4C to 2/3 under dichoptic masking. An L2/3 c50 shift with no 4C shift implies a recurrent or local origin.
   - Design fit: **Yes,** at MUA/CSD level.
3. **Whether early deep-layer suppression comes from L6 thalamic input or descends from L2/3.**
   - Measurement: whether L6 suppression onset precedes or follows L2/3 suppression onset.
   - Design fit: **Partial.** Feasible for contrast-imbalance suppression; the orientation-conflict condition is absent.
4. **Whether the time course of interocular normalization differs between OD-biased and binocular neurons and between layers, at ms resolution.**
   - Measurement: time-resolved (10 ms bin) dichoptic CRF fits by OD class and layer.
   - Design fit: **Yes,** though SU numbers per OD class and layer will be small.
5. **Whether the contrast dependence of cross-column coupling in awake animals matches anesthetized Utah-array data (Nauhaus).**
   - Measurement: cross-probe STA or CCG reach and amplitude vs contrast, controlling for γ.
   - Design fit: **Partial** (7 pairs).
6. **Laminar-marker agreement in V1.**
   - Measurement: CSD vs spike–LFP phase reversal vs spectral crossover per penetration.
   - Design fit: **Yes.**
7. **Whether layer-specific variability or noise-correlation changes accompany interocular normalization** (an extension of Nurminen 2026 to the contrast/eye dimension).
   - Measurement: Fano factor and shared variance by layer across the 6×6 matrix.
   - Design fit: **Yes** for MUA and SUs, given enough repeats per cell of the matrix. Check trial counts per condition first.

### Gaps
- Explicit "open question" statements from Self 2013/2019, Xing 2012, Gieselmann & Thiele 2022 and Bijanzadeh 2018 Discussion were not fully extracted.
- I could not verify any paper explicitly posing the question of *laminar timing of interocular contrast normalization in awake macaque with a full dichoptic contrast matrix*. Its absence from the sources reviewed suggests the question is open, but absence of evidence across a limited search is not proof.
