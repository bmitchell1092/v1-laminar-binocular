# Timing of the granular (L4C) CSD sink versus multiunit activity in primate V1

Research notes compiled 2026-11-03. Scope: macaque V1 first, then other primate visual and auditory cortex. "Verified" means I read the statement in the paper's full text (or the abstract where noted). Anything I could not check against the source is marked **unverified**. Full texts that returned HTTP 403 (J Neurosci, PNAS publisher site, Cell Press) were read through PMC where a PMC copy existed.

---

## Q1. Schroeder, Mehta & Givre (1998) and the earlier Schroeder/Givre/Tenke/Mehta papers: L4C sink onset versus MUA onset

### Takeaway
Schroeder et al. (1998) measured onset latency on the CSD only, using a rectified, channel-averaged CSD (AVREC) once per penetration. They report lamina 4C onsets of about 26–29 ms after diffuse flashes but give **no separate MUA onset values** and **no explicit statement** that MUA leads or lags the 4C sink. They describe the two as co-occurring ("current sinks … and the associated increase in multiunit activity") and chose CSD over MUA because it is "a more direct index of afferent input". Givre, Arezzo & Schroeder (1995) put the onset of 4C transmembrane current at about 23 ms after full-field flashes; that figure is arithmetic from their abstract, not a stated value.

### Cited Findings
**Schroeder CE, Mehta AD, Givre SJ (1998). A spatiotemporal profile of visual system activation revealed by current source density analysis in the awake macaque. *Cereb Cortex* 8(7):575–592. DOI 10.1093/cercor/8.7.575.** Full text read from the PDF.
- Stimulus: diffuse light flashes "10 µs in duration and 7.8 × 10^5 lux intensity", from a Grass PS-2 photostimulator, presented at 2/s from a diffuser subtending 20° directly in front of the animal. Pattern flash and pattern reversal (3 c/deg square-wave grating, 60% contrast) were used only in limited comparisons. — [Schroeder et al. 1998 PDF](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf); [DOI](https://doi.org/10.1093/cercor/8.7.575)
- Why CSD and not MUA was the primary measure: "the main analyses of spatial (laminar) and temporal response patterns used the CSD measure in preference to multiunit activity because it is a more direct index of afferent input. That is, the transmembrane current flow measures index the PSP, which is the first-order response to synaptic input." MUA was used "to help in defining co-located CSD features as indices of excitation and inhibition." — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- How onset latency was defined: measured "once for each electrode penetration, using the … AVREC measures". AVREC means every CSD trace in the profile was baseline-corrected, full-wave rectified and averaged across channels. Onset was "the earliest significant (>2 SD units) deviation from the baseline", and the deviation had to stay at or above 2 SD "for 8 out of 10 of the following sampling points (each point corresponding to 0.5 ms at the 2 kHz sampling rate)". The baseline SD came from a prestimulus period of at least 20 ms. The window from 5 ms before to 10 ms after the stimulus was excluded to avoid stimulus artifacts. — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Signal processing:
  - Field potentials were digitized at 2 kHz.
  - MUA was "band pass filtered (0.5–2 kHz, roll-off 24 dB/octave), full-wave rectified, further amplified (gain = 8) and digitized at 4 kHz" (analog MUA envelope, not thresholded spikes).
  - CSD used the three-point second-derivative formula (Freeman & Nicholson 1975), with the grid equal to the intercontact spacing (100 or 150 µm).
  - [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Figure 3 (V1, 100 µm spacing, 14 contacts) has "Drop lines in the CSD and MUA columns [that] indicate the onset latency as determined by the statistical criteria used in this study". The figure therefore marks onsets for both signals, but the text gives no MUA-versus-CSD numbers. I could not read the figure values from the PDF text. — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Co-occurrence statement: "in both granular and supragranular laminae, the initial response is excitatory, as indicated by the current sinks (signaling net local depolarization) and the associated increase in multiunit activity." In summary: "there is a clear correspondence between postsynaptic potential and action potential indices (CSD and MUA respectively) of local activation." — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Latencies (CSD/AVREC-based):
  - "Lamina 4C response latencies for the different eccentricities of stimulation … ranged from 26 to 29 ms" (fixation offsets 0–32° from the diffuser).
  - Lamina 4C of foveal V1: "mean = 27.9 ms, SD = 3.5".
  - V1 overall: "mean 26.3 ms, range 18–37 ms".
  - Supragranular V1: 35.5 ms.
  - Foveal/parafoveal versus peripheral V1: 26.4 ms (SD 5.1) versus 19.3 ms (SD 2.9).
  - LGN grand-mean AVREC onset at 15 ms; LGN-to-V1 step 10.9 ms.
  - [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Comparison with spike data: their V1 latencies "align very well with those from … unanesthetized monkeys (Maunsell and Gibson, 1992), whose Lamina 4 latencies range from 20 to 31 ms". They also cite an 8–15 ms (elsewhere 9–15 ms) lag between thalamorecipient and upper laminae (Maunsell & Gibson 1992; Givre et al. 1994). — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- MUA can occur without a CSD signature: "Active white matter regions can be characterized … by multiunit activity in the presence of a near-zero second derivative". At the grey/white boundary the low-frequency CSD drops to near zero "though multiunit activity often extends considerably beyond this point". — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- The 4C subdivisions: "the stimulus sensitivities of the upper and lower Lamina 4C sinks show that they receive input from M- and P-type thalamic afferents respectively (Givre et al., 1995)". No separate 4Cα and 4Cβ onset values are given in the text. — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)

**Givre SJ, Arezzo JC, Schroeder CE (1995). Effects of wavelength on the timing and laminar distribution of illuminance-evoked activity in macaque V1. *Vis Neurosci* 12(2):229–239. DOI 10.1017/S0952523800007914.** Abstract only.
- Stimuli: full-field red, blue and green flashes versus illuminance-matched white flashes, in two alert macaques. Recordings in V1, optic radiations and LGN. — [Givre et al. 1995 (Cambridge)](https://www.cambridge.org/core/journals/visual-neuroscience/article/effects-of-wavelength-on-the-timing-and-laminar-distribution-of-illuminanceevoked-activity-in-macaque-v1/A6F43C4DE063E98F41A3A3526EF1849D)
- "In layer 4c, the red effect reached significance, on average, at 47 ms … or ≈24 ms after the onset of transmembrane current flow". In the supragranular layers: 55 ms, about 14 ms after current-flow onset. — [Givre et al. 1995](https://www.cambridge.org/core/journals/visual-neuroscience/article/effects-of-wavelength-on-the-timing-and-laminar-distribution-of-illuminanceevoked-activity-in-macaque-v1/A6F43C4DE063E98F41A3A3526EF1849D)
- Implied onset of 4c current flow: about 23 ms (47 − 24); supragranular about 41 ms (55 − 14). These are derived figures, not stated in the paper. The abstract gives no MUA latencies and does not separate 4Cα from 4Cβ. — [Givre et al. 1995](https://www.cambridge.org/core/journals/visual-neuroscience/article/effects-of-wavelength-on-the-timing-and-laminar-distribution-of-illuminanceevoked-activity-in-macaque-v1/A6F43C4DE063E98F41A3A3526EF1849D)

**Mehta AD, Ulbert I, Schroeder CE (2000). Intermodal selective attention in monkeys. II: physiological mechanisms of modulation. *Cereb Cortex* 10(4):359–370. DOI 10.1093/cercor/10.4.359.** Abstract only.
- Laminar ERP, CSD and MUA profiles in V1, V2, V4 and LGN. Attention "did not modulate either the response in lateral geniculate nucleus (LGN) or the initial excitation in lamina 4C of V1". The abstract gives no onset values. — [Mehta et al. 2000 (DOI)](https://doi.org/10.1093/cercor/10.4.359)

### Inferences
- In this lineage, the textual answer to "which is first in 4C" is that they co-occur. The quantitative onset rule was applied only to the CSD (AVREC), so these papers cannot be cited for a measured MUA-versus-sink order in 4C. Any claim that "Schroeder 1998 showed MUA precedes (or follows) the sink" would have to come from reading Fig. 3's drop lines, which I could not do.
- Because AVREC averages rectified CSD over all channels, its onset is the earliest current anywhere in the array, sink or source. It is not strictly the 4C sink onset.
- The note that MUA appears in white matter with a flat CSD shows the authors knew MUA can pick up afferent axonal activity with no local CSD. That is the mechanism by which granular MUA could lead the granular sink.

### Gaps
- Schroeder, Tenke, Givre, Arezzo & Vaughan (1991, *Vision Res* 31:1143–1157, pattern-reversal VEP), Givre, Schroeder & Arezzo (1994, *Vision Res*, V4/flash VEP), and Mehta, Ulbert & Schroeder (2000) Part I (*Cereb Cortex* 10:343–358) were not accessed. They may contain per-lamina MUA onset values. **Unverified.**
- The Fig. 3 and Fig. 14 values in Schroeder et al. (1998), including the MUA drop-line times, could not be extracted from the PDF text layer.
- Full text of Givre et al. (1995), which may contain separate 4Cα and 4Cβ sink and MUA onsets, was not accessed.

---

## Q2. Later macaque laminar studies (Maier, Self/van Kerkoerle/Roelfsema, Westerberg/Cox/Dougherty/Maier, Nandy, Bastos, Thiele, Angelucci, Kajikawa/Schroeder, Lakatos): explicit sink-versus-MUA timing statements in L4C

### Takeaway
Few papers report numbers for this comparison. The explicit statements are:
- **Coincidence:**
  - van Kerkoerle et al. (2014): the sinks and sources "coincided with the onset of the visual response in the MUA".
  - Gieselmann & Thiele (2022): the "early current sink coincides with the earliest activation from thalamo-cortical afferents"; the MUAe RF emerged at 30–50 ms.
- **Sink first:**
  - Westerberg et al. (2019): "the onset latency of current sinks … generally preceded that of spiking activity" (normalized granular and supragranular compartments, no ms values given).
- **Spikes nominally earlier, with a warning:** Bijanzadeh et al. (2018) give mean L4C onsets of about 39 ms for CSD and about 23 ms for spiking. They state these absolute values "cannot be directly compared" because the signals and methods differ.
- **During gamma oscillations, not stimulus onset:** van Kerkoerle et al. (2014) found the L4 CSD sink led L4 MUA by 17 ± 4° of phase.

### Cited Findings
**Maier A, Adams GK, Aura C, Leopold DA (2010). Distinct superficial and deep laminar domains of activity in the visual cortex during rest and stimulation. *Front Syst Neurosci* 4:31. DOI 10.3389/fnsys.2010.00031.** Full text via PMC.
- "This analysis revealed a robust short-latency current sink in the middle layers for each session". "Previous studies have shown that this sink in V1 corresponds most closely to layer 4Cα (Mitzdorf and Singer, 1979; Schroeder et al., 1991)". The initial sink "originates in layer 4C, possibly with its peak in layer 4Cα". — [Maier et al. 2010 (PMC2928665)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2928665/)
- The sink centre was the alignment "zero point", and ±200 µm around it was taken as the approximate extent of layer 4. CSD was computed with a standard second-derivative method and the iCSD toolbox spline methods (Pettersen et al. 2006). — [Maier et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2928665/)
- I found no statement comparing sink onset with MUA onset in this paper.

**Maier A, Aura CJ, Leopold DA (2011). Infragranular sources of sustained local field potential responses in macaque primary visual cortex. *J Neurosci* 31(6):1971–1980. DOI 10.1523/JNEUROSCI.5300-09.2011.** Abstract only; full text returned 403.
- "In agreement with previous studies, we found a short-latency phase-locked current sink, thought to correspond to thalamocortical input to layer 4C." A sustained, non-phase-locked CSD component sat about 500 µm deeper than the initial sink. — [Maier et al. 2011 (DOI)](https://doi.org/10.1523/JNEUROSCI.5300-09.2011)
- An earlier search-engine summary credited "a layer 4C sink that coincided with the onset of the visual response in the multiunit activity" to "one study". The wording I verified is in van Kerkoerle et al. (2014), below. I did not confirm it for Maier et al. (2011). **Unverified for Maier 2011.**

**van Kerkoerle T, Self MW, Dagnino B, Gariel-Mathis MA, Poort J, van der Togt C, Roelfsema PR (2014). Alpha and gamma oscillations characterize feedback and feedforward processing in monkey visual cortex. *PNAS* 111(40):14332–14341. DOI 10.1073/pnas.1402773111.** Full text via PMC.
- Laminar U-probes at 100 µm spacing. Depth was set with "the CSD response triggered by the appearance of a checkerboard stimulus". "The stimulus evoked a characteristic CSD profile with a current sink in layer 4C and a current source in the deep layers … These sinks and sources coincided with the onset of the visual response in the MUA." The paper gives no ms values for this. — [van Kerkoerle et al. 2014 (PMC4210002)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/)
- Gamma, with troughs detected after a 55–65 Hz filter (not stimulus onset): "The first coherence peak reflected the coincidence of the MUA with a current sink in layer 4 [mean phase advance of CSD sink relative to MUA: 17 ± 4° (SEM); n = 24]". Also, "Spiking activity in layer 4 coincided with the layer 4 current sink, and was delayed by 1 ms in the deep layers … and by 2 ms in the superficial layers". — [van Kerkoerle et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/)
- Alpha (8–12 Hz): MUA was in phase with sinks in L1/2 and L5 but "with sources in layers 4 and 6 (mean phase in L4, −174 ± 7°…)". "MUA in layer 4 lagged MUA in the deep and superficial layers by 12 and 10 ms". This concerns alpha-trough alignment, not stimulus-onset latency. — [van Kerkoerle et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/)

**Self MW, van Kerkoerle T, Supèr H, Roelfsema PR (2013). Distinct roles of the cortical layers of area V1 in figure-ground segregation. *Curr Biol* 23(21):2121–2129. DOI 10.1016/j.cub.2013.09.013.** Abstract only; Cell Press returned 403 and there is no PMC copy.
- "We found that the visual response started 40 ms after stimulus presentation in layers 4 and 6". Boundary-related synaptic currents and MUA appeared in upper L4 and the superficial layers at about 70 ms, and feedback-layer (1, 2, 5) inputs about 30 ms later. The abstract does not say whether the 40 ms is CSD or MUA. — [Self et al. 2013 (DOI)](https://doi.org/10.1016/j.cub.2013.09.013)
- Gieselmann & Thiele (2022) state that their MUAe latency procedure was "adapted from Self et al., 2013". Self et al. therefore used a curve fit for MUA latency; details are under Q4. — [Gieselmann & Thiele 2022 (PMC8916775)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/)

**Westerberg JA, Cox MA, Dougherty K, Maier A (2019). V1 microcircuit dynamics: altered signal propagation suggests intracortical origins for adaptation in response to visual repetition. *J Neurophysiol* 121(5):1938–1952. DOI 10.1152/jn.00113.2019.** Full text via PMC.
- Granular contacts were "determined … on the basis of the initial sink of a flash-evoked CSD". "The initial current sink in V1 granular layer 4C following onset of visual stimulation is indicative of retinogeniculate synaptic activation". — [Westerberg et al. 2019 (PMC6589708)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/)
- "In line with the causal relationship between synaptic depolarization and spiking activity, we found that the onset latency of current sinks (SiCSD) generally preceded that of spiking activity (Fig. 9A)." Fig. 9A shows the "initial transient spiking (solid lines) and current source density sink (SiCSD; dashed lines)" for the granular and supragranular compartments. Responses were averaged over 61 sessions and "normalized between the baseline and maximum response for comparison". The text gives no ms values. — [Westerberg et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/)
- The spiking signal here is single-unit activity, with the abstract and results referring to "single-unit spiking". Stimuli were static gratings presented repeatedly, 200 ms on (500 ms in a control). — [Westerberg et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/)

**Cox MA, Dougherty K, Westerberg JA, Schall MS, Maier A (2019). Temporal dynamics of binocular integration in primary visual cortex. *J Vis* 19(12):13. DOI 10.1167/19.12.13.** Full text via Europe PMC. I found no sentence comparing CSD sink timing with MUA or spike timing. — [Cox et al. 2019 (PMC6797477)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477/)

**Gieselmann MA, Thiele A (2022). Stimulus dependence of directed information exchange between cortical layers in macaque V1. *eLife* 11:e62949. DOI 10.7554/eLife.62949.** Full text via Europe PMC. The earlier search result attributed this paper to Ferro et al.; Europe PMC lists Gieselmann & Thiele.
- In an example session, the alignment channel ("putative layer 4c at zero depth", the early sink) also showed "the early emergence of a visual RF in the same channel between 30 and 50ms" in the MUAe RF maps, "demonstrating that the early current sink coincides with the earliest activation from thalamo-cortical afferents". — [Gieselmann & Thiele 2022 (PMC8916775)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/)
- MUAe latency by depth: "The shortest latencies occurred in granular layers, followed by infragranular layers, with supragranular layers generally lagging behind". Larger gratings lengthened granular latency by up to about 3 ms and shortened infragranular latency by up to about 6 ms, so "the lag between granular and infragranular layers was reduced to ~2ms (from ~7ms for small stimuli)". Stimuli were 0.5° to 15° gratings. — [Gieselmann & Thiele 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/)

**Ferro D, van Kempen J, Boyd M, Panzeri S, Thiele A (2021). Directed information exchange between cortical layers in macaque V1 and V4 and its modulation by selective attention. *PNAS* 118(12):e2022097118. DOI 10.1073/pnas.2022097118.** Full text via Europe PMC.
- "The earliest current sink of the CSD and the shortest multiunit response latency (multiunit activity envelope [MUA_E]) identified input layer IV". The two are used together as layer markers. No order between them is reported in the main text. — [Ferro et al. 2021 (PMC8000025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8000025/)

**Bijanzadeh M, Nurminen L, Merlin S, Clark AM, Angelucci A (2018). Distinct laminar processing of local and global context in primate primary visual cortex. *Neuron* 100(1):259–274. DOI 10.1016/j.neuron.2018.08.020.** Full text via PMC.
- "Stimuli inside the RF … evoked the earliest current sinks in L4C (mean ΔL ± s.e.m. = 0 ± 0 ms)". Relative CSD onset lags were: L2/3 23.44 ± 3.09 ms, L1 23.44 ± 1.55 ms, upper deep 14.29 ± 1.96 ms, lower deep 61 ± 6.68 ms. Upper L4 (16.91 ± 6.96 ms) did not differ significantly from L4C. — [Bijanzadeh et al. 2018 (PMC6279245)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- "RF stimulation evoked the earliest spiking in L4C (mean ΔL ± s.e.m = 3.5 ± 1.35 ms …)". — [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- Their caveat: "the absolute latencies of CSD signals … cannot be directly compared to the absolute latencies of spiking responses … due to the fundamentally different nature of CSD signals and spikes, and the different methods we used to compute onset latency of the two signals. For example, in L4C the average onset latency of responses to RF stimulation was about 39 ms, for CSD signals (Table S1), but 23 ms for spiking responses (Table S6)." — [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- The stimuli also differed. CSD latencies came from "a black square of 0.5° side systematically flashed" (plus annular gratings for surround conditions). The spiking latency in Fig. 7A is "onset latency of spiking responses to gratings in the RF" (6 penetrations, 79 contacts). — [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- Their rationale for using sinks: sink onsets "reflect the net post-synaptic potentials of local neurons" and are "better suited to localize laminar activation than measuring onset latency of spiking responses", because dendrites "rarely co-localize with the soma". — [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)

**Zhang LA, Li P, Callaway EM (2024). High-resolution laminar identification in macaque primary visual cortex using Neuropixels probes. bioRxiv 2024.01.23.576944 (v3), also an eLife reviewed preprint. DOI 10.1101/2024.01.23.576944.** Full text via PMC.
- "The shortest latency current sink in cortical area V1 CSD profiles evoked from flashed, full-field visual stimuli" is taken to roughly correspond to L4C. They state that at high resolution the CSD profiles are fragmented and the short-latency sink spans only about 200 µm, narrower than L4C. They also warn that a "zero latency sink corresponds to residual signal from the prior white to black transition". — [Zhang, Li & Callaway (PMC10849622)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849622/); [eLife reviewed preprint](https://elifesciences.org/reviewed-preprints/97290)
- "responses in L4Cα are typically faster than L4Cβ", measured as AP-band (300–3000 Hz) power change, not MUA onset. No sink-versus-spike onset comparison is reported. — [Zhang, Li & Callaway](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849622/)

**Other papers checked, with no sink-versus-MUA timing statement found:**
- Nandy AS, Nassi JJ, Reynolds JH (2017). *Neuron* 93(1):235–246, V4 (PMC full text). — [PMC5217483](https://pmc.ncbi.nlm.nih.gov/articles/PMC5217483/)
- Xing D, Yeh CI, Burns S, Shapley RM (2012). *PNAS* 109:13871–13876. MUA and gamma-LFP laminar dynamics to drifting gratings; no CSD-versus-MUA onset statement found. — [PMC3427063](https://pmc.ncbi.nlm.nih.gov/articles/PMC3427063/)
- Ninomiya T, Dougherty K, Godlove DC, Schall JD, Maier A (2015). *J Neurophysiol* 113:3242–3255. — [PMC4440241](https://pmc.ncbi.nlm.nih.gov/articles/PMC4440241/)
- Westerberg JA, Schall MS, Maier A, Woodman GF, Schall JD (2022). *eLife* 11:e72139, V4. — [PMC8846592](https://pmc.ncbi.nlm.nih.gov/articles/PMC8846592/)
- van Kerkoerle T, Self MW, Roelfsema PR (2017). *Nat Commun* 8:13804. — [DOI](https://doi.org/10.1038/ncomms13804)
- Kajikawa Y, Schroeder CE (2011). *Neuron* 72:847–858. A paper on the spatial reach of the LFP; no onset-order statement in the Europe PMC text. — [DOI](https://doi.org/10.1016/j.neuron.2011.09.029)

### Inferences
- Read together, the literature does not establish a consistent, quantified lead for either signal at stimulus onset in macaque L4C. The two statements that use the same timescale are close to the noise level of the methods. van Kerkoerle et al. (2014) report coincidence at stimulus onset. Westerberg et al. (2019) report a sink lead on normalized, compartment-averaged traces, with single units, a granular compartment spanning several contacts, and no ms value.
- The 17° gamma-phase advance of the L4 sink over L4 MUA in van Kerkoerle et al. (2014) is about 0.8 ms at 60 Hz (my arithmetic; the paper does not state the frequency for this conversion). That fits a sub-millisecond to roughly 1 ms "sink before spikes" order within the ongoing gamma cycle, but it is not an onset-latency measure.
- Bijanzadeh et al. (2018) give a nominal 39 ms (CSD) versus 23 ms (spikes) in L4C, which would put spikes about 16 ms first. The authors explicitly say this difference reflects methods and stimuli, not physiology, so it should not be cited as "MUA precedes sink".

### Gaps
- Self et al. (2013) and Self et al. (2019, *Curr Biol* 29:1019–1029, DOI 10.1016/j.cub.2019.02.016) full texts were not accessible (403, no PMC copy). Their L4 CSD-versus-MUA latencies, if reported in figures or supplements, remain **unverified**.
- Maier et al. (2011) full text returned 403. Dougherty et al. (2017, *Cereb Cortex*, DOI 10.1093/cercor/bhv304) and Dougherty et al. (2019, *Curr Biol*, binocular modulation) were not read in full.
- I found no Bastos et al. paper in macaque V1 with an explicit sink-versus-MUA onset comparison. The Europe PMC query returned nothing relevant.
- The supplementary tables in Bijanzadeh et al. (Tables S1, S6) and the van Kerkoerle et al. SI (Fig. S1, MUA filtering) were not read.

---

## Q3. Reports that granular MUA precedes the granular sink (thalamic afferent activity), or that the sink precedes MUA (postsynaptic depolarization before spikes)

### Takeaway
I found no primate V1 paper that explicitly reports granular MUA onset **preceding** the granular CSD sink, measured with matched methods. The only "spikes earlier" numbers (Bijanzadeh 2018) are disowned by their authors as a methods artifact. The explicit directional claim in primate V1 runs the other way: sink onset generally precedes spiking (Westerberg 2019), plus a small gamma-phase sink lead (van Kerkoerle 2014). Several papers describe the two as coinciding. The idea that MUA can lead the sink rests on indirect evidence. MUA is present in white matter where the CSD is flat (Schroeder 1998), and in primate A1 the MUA of the initial lamina 4 response extends into the white matter (Steinschneider 1992).

### Cited Findings
- Sink first: "the onset latency of current sinks (SiCSD) generally preceded that of spiking activity (Fig. 9A)", explained by "the causal relationship between synaptic depolarization and spiking activity". The comparison is macaque V1, granular and supragranular compartments, single units, gratings. — [Westerberg et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/)
- Sink slightly first within the gamma cycle: "mean phase advance of CSD sink relative to MUA: 17 ± 4°" in layer 4. — [van Kerkoerle et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/)
- Coincident at stimulus onset: the checkerboard-evoked L4C sink and deep source "coincided with the onset of the visual response in the MUA". — [van Kerkoerle et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/)
- Coincident at stimulus onset: "the early current sink coincides with the earliest activation from thalamo-cortical afferents" (MUAe RF at 30–50 ms in the sink channel). — [Gieselmann & Thiele 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/)
- Co-occurrence without an order: "current sinks (signaling net local depolarization) and the associated increase in multiunit activity". — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Spikes nominally first, but disclaimed: L4C onset about 39 ms for CSD versus 23 ms for spiking, which the authors say "cannot be directly compared". — [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- Afferent axonal MUA without CSD: "Active white matter regions can be characterized … by multiunit activity in the presence of a near-zero second derivative". — [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- Primate auditory cortex (secondary): Steinschneider, Tenke, Schroeder, Javitt et al. (1992, *Electroencephalogr Clin Neurophysiol* 84:196–200). Per the abstract as summarized by the search engine, the click-evoked N8 peaks at 8–9 ms and is generated by a lamina 4 current sink with a deeper source. Simultaneous MUA spans lower lamina 3 to the subjacent white matter, and N8 is attributed to thalamocortical afferents (possibly with lamina 4 stellate cells). I did not read the abstract text directly, and it does not state an MUA-versus-sink onset order. — [Einstein research portal entry](https://einstein.elsevierpure.com/en/publications/cellular-generators-of-the-cortical-auditory-evoked-potential-ini-2/) (**partially verified**)
- Primate auditory cortex (secondary): a 2018 *eNeuro* study of awake macaque A1 and belt (ENEURO.0420-17.2018) reportedly found that "LFP onset latencies predominantly reflecting synaptic input arrival typically preceded spike onset latencies". This comes from the search-engine summary; authors and values are not verified. — [eNeuro PDF](https://www.eneuro.org/content/eneuro/5/3/ENEURO.0420-17.2018.full.pdf) (**unverified**)
- Lakatos P, Chen CM, O'Connell MN, Mills A, Schroeder CE (2007), *Neuron* 53:279–292 (A1). The abstract describes laminar CSD and MUA but makes no onset-order statement. — [DOI](https://doi.org/10.1016/j.neuron.2006.12.011)

### Inferences
- Biophysically, both orders are plausible and depend on what the MUA electrode picks up. Thalamocortical terminal and axonal spikes in 4C should precede the postsynaptic sink they generate, by roughly the synaptic delay plus rise time (about 1–2 ms). Somatic spiking of 4C neurons should follow the sink. A 0.5–2 kHz or 0.6–9 kHz analog MUA envelope in 4C mixes both sources. Which one dominates the first detectable MUA deflection is not measured in any primate V1 paper I found.
- The two published directional statements, coincidence and a small sink lead, are consistent with a 4C MUA onset dominated by postsynaptic somatic spiking, not afferent spikes. But neither used a per-contact, matched-criterion onset analysis with sub-millisecond precision.

### Gaps
- No primate V1 study found reports per-contact 4C MUA onset minus 4C sink onset (in ms) with identical onset criteria and filtering on both signals. This is the measurement needed to settle the order.
- The cat and rodent literature on afferent volleys preceding layer 4 sinks (e.g., Mitzdorf & Singer 1979, cat) was outside scope and not reviewed.

---

## Q4. Onset criteria and filtering used

### Takeaway
Onset criteria vary widely: a 2 SD threshold on channel-averaged rectified CSD (Schroeder 1998), 3–5 SD (or 7 SD for averaged rectified CSD) for three consecutive bins on per-depth CSD (Bijanzadeh 2018), curve-fit models of MUAe (Self 2013, Gieselmann & Thiele 2022), and visual comparison of normalized traces (Westerberg 2019). MUA filtering ranges from 0.5–2 kHz (Schroeder 1998) to 0.6–9 kHz high-pass, rectified and low-passed at 200 Hz (Roelfsema/Thiele MUAe). These differences alone can shift onsets by many milliseconds, as Bijanzadeh et al. state explicitly.

### Cited Findings
- **Schroeder et al. 1998:**
  - Onset rule: on AVREC, the first point >2 SD above a baseline of ≥20 ms, sustained for 8 of the next 10 samples (0.5 ms each), with −5 to +10 ms excluded.
  - MUA: 0.5–2 kHz bandpass (24 dB/oct), full-wave rectified, digitized at 4 kHz.
  - LFP digitized at 2 kHz.
  - CSD: three-point second derivative at 100 or 150 µm spacing.
  - [Schroeder et al. 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf)
- **Bijanzadeh et al. 2018:**
  - CSD was z-scored against a 200 ms prestimulus baseline and computed every 0.5 ms.
  - Sink onset was "the earliest time bin after stimulus onset in which the CSD amplitude was 3–5SD below baseline for three consecutive bins", taken at each interpolated depth. Grand averages of rectified CSD used 7 SD.
  - Laminar ordering used relative latency ΔL, the onset minus the shortest onset across layers in the same penetration.
  - [Bijanzadeh et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/)
- **Gieselmann & Thiele 2022 (Thiele lab, adapted from Self et al. 2013):**
  - Filtering: LFP low-pass "0.75–300 Hz" (3rd-order IIR Butterworth), downsampled to about 1 kHz. MUAe used a "high-pass filter (IIR-Butterworth, 3rd order, 0.6–9 kHz)", full-wave rectification, then low-pass at 200 Hz and downsampling to about 1 kHz (citing Supèr & Roelfsema 2005).
  - Latency method: a model of the transient from 0–200 ms was fitted as the sum of two Gaussians plus a cumulative Gaussian, by nonlinear least squares, to MUAe z-scored against 0–300 ms of spontaneous activity.
  - "different from Self et al., 2013, we used an absolute criterion to determine the response latency". This implies Self et al. (2013) used a relative criterion, such as a fraction of the fitted maximum. The exact fraction is **unverified**.
  - [Gieselmann & Thiele 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/)
- **Ferro et al. 2021:**
  - LFP band-pass 0.5–300 Hz, downsampled to 1,107.375 Hz.
  - Spiking band-pass 600–9,000 Hz, then MUA_E or sorted units.
  - Baseline from 200 ms before to 30 ms after stimulus onset.
  - [Ferro et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8000025/)
- **Westerberg et al. 2019:**
  - Latency comparison: Fig. 9A compares traces "normalized between the baseline and maximum response". No onset threshold is stated for that panel.
  - Adaptation-difference onsets (Fig. 9B): "a t-test at each data point with a P value <0.001 for more than 10 consecutive samples".
  - CSD smoothing: CSD maps were smoothed "with a two-dimensional Gaussian filter (s = 0.1 mm and 15 ms)", but this is stated as "For display purposes".
  - [Westerberg et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/)
- **Zhang, Li & Callaway (Neuropixels):**
  - LFP: 60 Hz notch, 1–100 Hz band-pass, downsampled to 1 kHz.
  - CSD: inverse-CSD with a Gaussian spatial smoothing of σ = 30 µm.
  - Windows: baseline −40 to 10 ms; response 10 to 150 ms.
  - Spiking proxy: AP-band 300–3000 Hz power.
  - [Zhang, Li & Callaway](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849622/)
- **Maier et al. 2010:** a standard second-derivative CSD plus iCSD spline methods (Pettersen et al. 2006). The centre of the initial sink was the alignment zero, with ±200 µm taken as L4. — [Maier et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2928665/)
- **Stimuli across studies:**
  - Diffuse 10 µs flashes at 2/s (Schroeder 1998).
  - Full-field colored and white flashes (Givre 1995).
  - Checkerboard onset (van Kerkoerle 2014).
  - Flashed stimuli for CSD alignment, gratings for spiking (Westerberg 2019).
  - A flashed 0.5° black square for CSD versus RF gratings for spikes (Bijanzadeh 2018).
  - 0.5–15° gratings (Gieselmann & Thiele 2022).
  - Full-field black/white flashes at 1.5 s per phase (Zhang et al. 2024).
  - [Schroeder 1998](https://pages.ucsd.edu/~msereno/569/readings/CSD-visualareas.pdf); [Givre 1995](https://www.cambridge.org/core/journals/visual-neuroscience/article/effects-of-wavelength-on-the-timing-and-laminar-distribution-of-illuminanceevoked-activity-in-macaque-v1/A6F43C4DE063E98F41A3A3526EF1849D); [van Kerkoerle 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210002/); [Westerberg 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589708/); [Bijanzadeh 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279245/); [Gieselmann & Thiele 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8916775/); [Zhang et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849622/)

### Inferences
- Several method choices bias the apparent order by more than the 1–2 ms physiological difference in question:
  - **Threshold on rectified versus signed CSD.** Rectified AVREC registers sources as well as sinks.
  - **Threshold level.** 2 SD versus 3–7 SD.
  - **Persistence rule.** 8 of 10 samples versus 3 consecutive bins.
  - **Temporal smoothing of the CSD.** A 15 ms Gaussian would pull onsets earlier, though Westerberg 2019 says it was for display only.
  - **MUA envelope filtering.** IIR filters and 200 Hz low-pass add group delay unless applied zero-phase; whether filtering was zero-phase is not stated in the sources I read.
  - **Relative (fraction-of-peak) versus absolute criteria.** Relative criteria interact with response amplitude and rise time, which differ between CSD and MUA.
- A fair comparison needs the same criterion on both signals (e.g., the same k·SD rule or the same fraction-of-peak on signed sink and MUA), zero-phase filtering, the same stimulus, and the same or adjacent contacts.
- For the BM-ephys/v1-laminar-binocular data specifically, the repo notes record a 3.4 ms offset between `.ns6` windowed reads and the NEV clock (`.ns6` first timestamp = 102, `.ns2` = 0). If LFP/CSD and MUA come from differently offset NSx streams, that offset is larger than the sub-millisecond to roughly 1 ms orders reported here and would have to be corrected before drawing any conclusion about order.

### Gaps
- The exact latency criterion in Self et al. (2013) (fraction of maximum or other) and in van Kerkoerle et al. (2014) (SI methods) is **unverified**.
- None of the accessed papers states whether MUA band-pass or low-pass filters were applied zero-phase (forward-backward) or causally.
- I did not check how Bijanzadeh et al. (2018) defined spiking onset; their STAR Methods for the MUA PSTH onset were not extracted.
