# LFP-based laminar markers in primate V1: spectrolaminar crossing, spike–LFP phase flip, and the CSD initial sink

Scope: where each marker sits relative to histological layers in macaque V1 (in µm from layer 4 / the 4C–5 boundary or the initial sink), how reliable it is, and published critiques. Researched 2026-11-12. Tool budget was limited, so several primary papers were read only in part (noted in Gaps).

## Q1. Spectrolaminar motif (Mendoza-Halliday, Major et al. 2024, Nat Neurosci): what does the alpha-beta/gamma crossing mark in V1, and how was it validated?

### Takeaway
The paper says the crossing marks layer 4 ("layer 4 or the layer 4-to-5 boundary" in examples). That claim rests on histology (electrolytic lesions) from only about 24 probes pooled across areas, with 3 to 5 in V1 (the text and the figure legend disagree), and gives no V1-specific distance. The pooled mean crossover-to-L4-centre distance is about 46 µm (±51 SEM) for unipolar LFP. The authors themselves flag V1 as an exception: alpha-beta power peaks in white matter there, probably by volume conduction from V2.

### Cited Findings
- Definition: the crossover is the channel where relative (per-frequency, normalised to the maximum channel) alpha-beta and gamma power are equal. Gamma (50–150 Hz) peaks in L2/3 and alpha-beta (10–30 Hz) in L5/6. The authors state "the crossover between relative gamma and alpha-beta power corresponded to layer 4" (Results, Fig. 4). In the example probes the crossover was "in layer 4 or the layer 4-to-5 boundary." — [Mendoza-Halliday et al. 2024, Nat Neurosci (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- Histological validation used electrolytic lesions plus Nissl in LIP (n=8), LPFC (n=10), MST (n=2), V1 (n=3) and PMd (n=1). The text says V1 n=3 but the Extended Data Fig. 4 legend reportedly says n=5 (an internal discrepancy). Distances were measured "from the center of cortical layer 4." The V1 lesions were in one monkey ("Bo"), made at −10 to −20 µA for 20 s. — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- Pooled over LIP, LPFC, MST and V1, the median laminar positions (layer numbers) were: gamma peak 3 (95% CI 2.63–3.2), crossover 4 (3.82–4.45), alpha-beta peak 6 (5.43–6.05). The CSD sink median was 4.0 (3.6–4.67, n=14) against n=24 for the crossover (Fig. 4e). — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- Distance figures: the mean crossover-to-L4 distance was 58 µm in LIP (CI −19 to 135 µm) and 27 µm in LPFC (CI −183 to 237 µm). Extended Data Fig. 9 gives 46 ± 51 µm (SEM) for unipolar LFP across areas, against 96 µm for bipolar and 125 µm for CSD-derived "drop-off" distances. **No V1-specific distance is given in the text.** — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- V1 caveat from the authors: "V1 had the highest alpha-beta power in white matter channels, likely due to volume conduction" from V2. — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- FLIP regresses mean relative power in the alpha-beta (10–19 Hz) and gamma (75–150 Hz) bands over a contiguous channel range r, maximising a goodness-of-fit G built from two regressions of opposite slope. A negative G is read as an inverted probe. FLIP identified 64% of probes and manual inspection 61%. vFLIP searches band pairs up to 150 Hz and identified 81% (766/942) of macaque probes. The probes in the paper were Neuropixels with 20 µm spacing, interpolated to 24 channels spanning cortex for FLIP. — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/); code/data on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.9w0vt4bnp)
- Identifiability against CSD: the CSD early sink was identifiable in 51% of 810 probes, against 61% (manual), 64% (FLIP) and 81% (vFLIP) for the motif. CSD identifiability fell for high-angle probes (P=0.002) and the motif's did not. Image similarity between perpendicular and angled probes was 0.24 for CSD against 0.92 for the motif (Fig. 6f–i). CSD sinks were more identifiable in V4, MST and MT than in LIP and LPFC (Extended Data Table 1). — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- The authors say the generating mechanisms "remain unknown," and that deep alpha-beta may partly reflect pyramidal-cell geometry and volume conduction. — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)

### Inferences
- In V1 specifically, the histological support for "crossing = L4" is thin: a handful of lesioned probes, all from one animal, and no reported V1 distance. White-matter alpha-beta in V1 would steepen the deep end of the alpha-beta gradient and could push the crossing downward, towards the 4C/5 boundary or below.
- The ~46 µm pooled figure is an SEM-level summary across mostly non-V1 areas, measured from the **centre** of L4. In macaque V1, L4 (4A to 4C) spans several hundred µm, so "within L4" is not the same as "at the 4C/5 boundary."

### Gaps
- Per-area V1 identifiability percentages (Extended Data Fig. 8b) and the V1 lesion results (Extended Data Fig. 4) were not read directly. The n=3 vs n=5 discrepancy for V1 is unresolved.
- No V1-specific crossover-to-sink distance in µm was found in the text.

## Q2. Critiques and replications of the spectrolaminar motif, including V1 data

### Takeaway
An independent V1/A1/belt dataset (Mackey, Schroeder and colleagues) found the FLIP crossing in L4 in only **8% of fitted V1 penetrations**, with 71% of all crossings infragranular, using a CSD/MUA-defined L4. The original authors' rebuttal finds the motif present in 67% of Mackey's V1 probes and reports their revised vFLIP2 crossover 80 µm from histological L4 centre (CI −92 to 260 µm). Both pieces were published as Nature Neuroscience Matters Arising on 12 Dec 2025. The dispute partly reduces to which ground truth is used: CSD sink or histology.

### Cited Findings
- Mackey CA et al., "Is there a ubiquitous spectrolaminar motif of local field potential power across primate neocortex?" (bioRxiv 2024.09.18.613490, v3, Sep 2024; later a Nat Neurosci Matters Arising). They used A1, belt and V1 with 2 monkeys per area. L4 was defined as the earliest current sink within 100 µm of the peak-MUA site. — [Mackey et al. 2024 bioRxiv (PMC11429918)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11429918)
- Mackey et al., FLIP results:
  - FLIP fitted 64% of penetrations overall, but only **50% in V1** (73% in A1, 71% in belt).
  - The crossing fell in L4 in only 29% of fitted cases (A1 32%, belt 40%, **V1 8%**), and was **infragranular in 71%**.
  - FLIP returned negative G ("inverted") on 2 penetrations each in A1 and V1, which they call anatomically implausible for pial-surface penetrations of flat cortex ("motif inversion errors").
  - vFLIP identified L4 correctly in 42% of fitted V1 penetrations.
  - vFLIP and FLIP agreed on only 49% of penetrations in the authors' shared dataset (56% at a 500 µm tolerance).
  - Referential LFP is subject to volume conduction, which inflates false positives, and false positives drop with bipolar or CSD derivations.
  - Distances are reported only as layer assignments, not µm.
  
  — [Mackey et al. 2024 (PMC11429918)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11429918)
- The Transmitter (12 Dec 2025) reports that both Matters Arising pieces appeared that day in Nat Neurosci: the critique [doi:10.1038/s41593-025-02167-y](https://doi.org/10.1038/s41593-025-02167-y) and the reply [doi:10.1038/s41593-025-02168-x](https://doi.org/10.1038/s41593-025-02168-x). It summarises the critique as "~60% of recordings" showing the motif, with L4 found in 40%/32%/8% (belt/A1/visual), shifted deep. An independent expert (Morales-Gregorio) calls the 2024 conclusions "way too optimistic" about universality, while noting that CSD depends on sensory responses and is sensitive to noise and electrode angle. A second expert (van Bree) says the motif "holds up to scrutiny." — [The Transmitter, López Lloreda 2025](https://www.thetransmitter.org/brain-waves/dispute-erupts-over-universal-cortical-brain-wave-claim/)
- Rebuttal preprint: Major, Abdaltawab, ... Maier, Desimone, Miller, Bastos & Mendoza-Halliday, "A ubiquitous spectrolaminar motif across independent studies, including Mackey et al.'s own data" (bioRxiv 10.1101/2025.05.07.652644, posted 8 May 2025; the peer-reviewed form is presumably the Matters Arising reply).
  - Manual inspection found the motif in 66% (42/64) of Mackey probes: A1 65% (17/26), **V1 67% (16/24)**, belt 64% (9/14).
  - vFLIP2 found it in 78% of Mackey probes and 80% of their own (N=923).
  - The **vFLIP2 crossover was 0.08 mm from the centre of histologically verified L4 (N=18 probes; CI −0.092 to 0.26 mm)**.
  - They argue that using the first CSD sink as ground truth is "circular" and "less accurate" than histology: "Anatomy remains the only definitive ground truth."
  - No new V1 dataset is presented beyond Mackey's.
  
  — [Major et al. 2025 bioRxiv (PMC12247776)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12247776/)
- Earlier V1 laminar LFP work does not cleanly show "alpha-beta deep" in V1. Maier, Adams, Aura & Leopold 2010 (Front Syst Neurosci 4:31) used 24 contacts at 100 µm, 13 sessions and 2 monkeys, with depth zero at the centre of the initial sink (attributed to 4Cα):
  - coherence segregated into superficial and deep domains with a **sharp transition "near the bottom of layer 4"**, at about one contact (~100 µm) below the sink centre (Fig. 4B, 5B);
  - gamma (30–100 Hz) power was ~50% higher in granular and supragranular than infragranular sites;
  - **5–20 Hz power showed no laminar difference**;
  - inter-session alignment was estimated as good to within ~200 µm, and there was no histology.
  
  — [Maier et al. 2010](https://www.frontiersin.org/journals/systems-neuroscience/articles/10.3389/fnsys.2010.00031/full)
- van Kerkoerle et al. 2014 (PNAS 111:14332) found that in V1, alpha is directed from superficial and deep layers toward L4 and gamma the opposite way. Microstimulating V1 evoked gamma in V4, and microstimulating V4 evoked alpha in V1. — [van Kerkoerle et al. 2014 (PMC4210002)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002) (details from search summaries, not a full-text read)

### Inferences
- For V1 specifically, the only independent test (Mackey) suggests that the FLIP crossing usually sits **below** the CSD-defined L4, in the infragranular layers. This fits the V1 white-matter alpha-beta caveat in the original paper and the absence of a deep alpha-beta peak in Maier et al. 2010. Mackey's ground truth is the sink (an upper/mid-4C marker in V1; see Q4), so part of the "infragranular" displacement could be the distance between the sink and the bottom of 4C. That distance is probably ~100–250 µm, an estimate from Maier's ±200–250 µm granular zone.
- The motif is more robust as a coarse superficial-to-deep gradient (presence 60–80%) than as a precise L4 point marker. As a sub-layer anchor in V1, its error bars (vFLIP2 CI up to +260 µm; CSD-based V1 L4 hit rate 8–42%) are larger than the thickness of 4C.

### Gaps
- The final published Matters Arising texts (Dec 2025) were not read, and their numbers may differ from the preprints.
- I did not find an independent histology-based V1 test of the crossing (lesions plus Nissl in V1 other than the original paper's 3 to 5 probes).
- Spaak et al. 2012, Bollimunta et al. 2008/2011, Haegens et al. 2015, Bastos et al. 2015/2018 and Sotero were not fetched. Their laminar alpha/gamma claims for V1 are unverified here. A search snippet attributed a layer-5 origin of alpha in macaque V1 to a CSD-plus-coupling study, which I could not confirm as Spaak or Bollimunta.
- I found no Gieselmann/Thiele or Haegens commentary specifically on the 2024 motif.

## Q3. Spike–LFP phase coupling (Davis, Dotson, Franken, Muller & Reynolds 2023, eLife): where is the preferred-phase flip, and has it been used in V1?

### Takeaway
The preferred-phase flip marks the boundary between the **input layer (earliest CSD sink) and the deep layers (L5/6)**. It was validated **against CSD only, not histology**, in macaque V4 and marmoset MT/PFC, **never V1**. It agreed with CSD on average (r=0.64) but differed by at least 200 µm in 14 of 34 sessions, where it tracked the LFP power inversion better than CSD did.

### Cited Findings
- Paper: Davis ZW, Dotson NM, Franken TP, Muller L, Reynolds JH (2023), eLife 12:e84512, peer-reviewed. Method name: laminar-phase coupling (LPC). — [Davis et al. 2023 eLife](https://elifesciences.org/articles/84512); [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Location: the preferred spike phase reverses at the input/deep boundary (Fig. 3c–e). "Input" is tied to the earliest CSD sink, and "deep" to L5/6. Layer 4 is not named and there is no histology of their own. — [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Data: macaque V4 (2 monkeys, 6 sessions each), marmoset MT (8 sessions) and PFC (14 sessions), 34 sessions in all, recorded with 32-channel linear probes inserted perpendicular to the surface. **No V1.** Inclusion required a recognisable CSD source/sink pattern. — [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Agreement with CSD: across 34 sessions, LPC depth correlated with CSD depth (r=0.64, p=0.00005), with no significant mean difference (the extracted figure is "271 ± 49 µm SEM; p=0.98"; this may be the mean absolute difference and should be checked in Fig. 3f). 14/34 sessions disagreed by ≥200 µm, up to ~600 µm (Fig. 5d). In those sessions LPC matched the 8–30 vs 65–100 Hz power inversion better (r=0.87 vs 0.65 for CSD), and only LPC boundaries gave significant depth differences in firing rate (Fig. 6). — [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Signals and metrics:
  - LFP was filtered 5–50 Hz and phase was taken by "generalized phase." The reversal also appears in theta, alpha, beta and low gamma with Hilbert phase (Fig. 3–S2).
  - SPI is the mean resultant length: 0.162 superficial, 0.122 input, 0.073 deep (Fig. 2).
  - Preferred phase is about −2.7 to −2.9 rad above the input layer against −1.44 rad in deep layers, i.e. ~1.3–1.4 rad (~75–80°) apart, not a full π inversion (Fig. 2f).
  - The paper does not use PPC; any PPC step is an addition in our pipeline.
  
  — [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- Practicalities and caveats:
  - About 1 minute of data is enough. The result is robust to task state (task, fixation, dark), referencing (Fig. 3–S5) and filtering, and holds for single units (Fig. 3–S3).
  - Hundreds of MUA spikes per channel are needed, and thousands are better.
  - The boundary was chosen by eye from cross-channel plots, blind to CSD.
  - The reversal drifted 132 ± 26 µm deeper between the first and last 2 minutes (Fig. 3–S6), which they attribute to probe drift.
  - The mechanism is unresolved.
  - Generality beyond these areas and states is untested.
  - Channel spacing was not stated in the text that was extracted.
  
  — [PMC10162800](https://pmc.ncbi.nlm.nih.gov/articles/PMC10162800)
- I found no 2024–2026 paper applying LPC to macaque V1, and no independent critique of LPC. — searches summarised in [eLife 84512](https://elifesciences.org/articles/84512)

### Inferences
- If LPC transfers to V1, it should mark roughly the **bottom of 4C / top of L5**, not the initial-sink centre. The sink centre in V1 is usually placed in 4Cα (Q4), so the flip would be expected **~100–300 µm below the sink centre**. This estimate rests on Maier et al. 2010's coherence transition one contact (~100 µm) below the sink and on 4C thickness. It is not a measured V1 value.
- The Maier et al. 2010 V1 coherence transition "near the bottom of layer 4" is the closest V1 analogue, and supports the existence of a sharp superficial/deep LFP boundary at about 4C/5.
- The ~1.4 rad phase difference and ≥200 µm disagreements in 41% of sessions suggest that LPC should be cross-checked per penetration, not taken as ground truth.

### Gaps
- There is no histology for LPC in any area. V1 behaviour (with its very thick, dense L4 and strong thalamic sink) is untested.
- Fig. 3f statistics need direct reading to confirm what "271 ± 49 µm" refers to.

## Q4. Where the CSD initial sink sits, and other automated or multi-feature laminar markers in V1

### Takeaway
In macaque V1 the earliest sink is conventionally placed in 4C, specifically 4Cα. The highest-resolution test (Neuropixels at 20 µm; Zhang, Li & Callaway 2024) finds that CSD is low-resolution and inconsistent across penetrations and stimuli, with sink-to-source transitions that often do not match laminar borders. They identify layers instead from spike waveform, unit density, AP-band power/coherence and latency. No paper found compares the phase flip, the spectral crossing and the sink against histology in the same V1 penetrations.

### Cited Findings
- Maier et al. 2010 placed depth zero at the initial-sink centre, attributed to 4Cα following Mitzdorf & Singer 1979 and Schroeder et al. 1991, which used microlesions and histology. They took L4 as approximately ±200 µm (Methods) or ±250 µm (Results, Fig. 2A) around the sink. — [Maier et al. 2010](https://www.frontiersin.org/journals/systems-neuroscience/articles/10.3389/fnsys.2010.00031/full)
- Zhang LA, Li P, Callaway EM, "High-resolution laminar identification in macaque primary visual cortex using Neuropixels probes." bioRxiv 2024.01.23.576944 (v3 Sep 2024), also an eLife reviewed preprint (RP97290). Methods: anaesthetised, paralysed monkeys (2), 35 penetrations (31 analysed), NP1.0 at 20 µm.
  - The shortest-latency sink "roughly corresponds" to L4C, but the strong sink spans only ~200 µm, narrower than L4C.
  - The sink-to-source criterion would put the L4C bottom at an implausible 1800 µm in one example (Fig. 2E).
  - Profiles vary across penetrations (Figs. 2E–H, 3). Black stimuli to the dominant eye gave the most consistent L4C sink. White stimuli sometimes gave early L4C sources, and L/M cone-opponent stimuli often lacked a source below L4C.
  - Gamma (30–100 Hz) LFP power and coherence profiles gave "less precise boundaries" than AP-band signals (Suppl. Fig. S8).
  - Layers are instead identified from unit density, spike spread and duration, peak–trough ratio, AP-band power, local coherence, ΔP/P and latency (separating 4Cα from 4Cβ). They did no histology of their own, comparing instead against published anatomy.
  - Cortex was ~1.5–2.0 mm thick, L4A <100 µm, and L1 ~120–150 µm.
  
  — [Zhang, Li & Callaway 2024 (PMC10849622)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849622/); [eLife reviews](https://elifesciences.org/reviewed-preprints/97290/reviews)
- eLife reviewers on Zhang et al.: CSD "produced inconsistent results across array penetrations and across visual stimulus conditions," and in some cases predicted boundaries inconsistent with other measures. One reviewer calls the paper a compelling argument for action potentials and stimulus-driven responses over CSD for V1 layer assignment. — [eLife RP97290 reviews](https://elifesciences.org/reviewed-preprints/97290/reviews) (seen through search snippets, because the page would not render)
- Mendoza-Halliday et al. 2024 found CSD sink identifiability 51% across areas, falling sharply with probe angle. V1 and V4 sinks were more identifiable than prefrontal ones. — [Mendoza-Halliday et al. 2024 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659/)
- Davis et al. 2023 note that CSD "can be susceptible to noise on individual channels resulting in errors in assigning laminar boundaries." — [Davis et al. 2023](https://elifesciences.org/articles/84512)

### Inferences
- The three markers target different anatomical points in V1:
  - the initial sink marks mid/upper 4C (4Cα);
  - the spike-phase flip marks the input/deep boundary (≈4C/5, by its definition);
  - the spectral crossing nominally marks "L4" (centre ±~50–250 µm), but in V1 data it tends to land deeper.
  
  Anchoring on any one of them implies different zero points, which may be 100–300 µm apart. This is an estimate from the numbers above, not a measured V1 offset.
- For 32-channel, ~100 µm-pitch V1 penetrations with strong visual drive (the situation in this project), the initial sink is likely to be well defined: V1 has a robust thalamic sink, and the Zhang et al. inconsistency mostly concerns sub-layer borders and source placement at 20 µm. The sink is therefore a defensible primary anchor. LPC is a sensible independent per-penetration check that targets 4C/5. The spectral crossing has the weakest V1-specific support and is best used as a coarse superficial/deep sanity check.

### Gaps
- No study was found (2023–2026) comparing CSD sink, LPC flip and spectral crossing against histology within the same macaque V1 penetrations. The proposed "Westerberg/Maier laminar alignment" multi-marker comparison and a "LaminarToolbox" could not be located, so they may not exist under those names.
- The final eLife version-of-record status of Zhang et al. could not be confirmed, because the eLife pages failed to load.
- Precise histological 4C thickness and the sink-centre to 4C/5 distance in awake macaque V1 were not established here.
