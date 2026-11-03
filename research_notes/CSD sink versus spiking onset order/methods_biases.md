# Methodological biases in CSD-vs-MUA latency comparisons (filtering, envelope, CSD computation, onset criteria)

Scope: how processing choices move the measured onset (half-max, threshold) and peak times of laminar CSD and analog MUA, judged against the claim that "MUA reaches half its peak ~3 ms before the sink current" (project result, `docs/09_paper_design.md`, V4 check J: MUA − current half-max −3.3 ms raw / −2.5 ms after a −0.8 ms filter-bias correction, dominant eye 90%; −3.0 / −2.2 ms at 45% binocular; 26 penetrations).

Source-quality legend: **[full text]** = read in the article itself this session; **[abstract/summary]** = only abstract, repository listing or search-engine summary seen, so treat details as unverified; **[project]** = this repo's own measured docs; **[own sim]** = my simulation (in Inferences only, never as a cited fact). Simulation script: `C:\Users\bmitc\AppData\Local\Temp\claude\c--Users-bmitc-Documents-GitHub-BM-ephys\e97c5ad9-a6a5-4557-a778-4a7637ac879a\scratchpad\sim\latbias.py` and `latbias2.py` (session scratchpad; copy them out if they need to be kept). Simulations ran at 30 kHz. Test responses were a raised-cosine rise of R ms, with 10–90% rise ≈ 0.59 R, then an exponential decay with τ = 25–40 ms. Filters were scipy Butterworth SOS.

---

## Q1. Zero-phase (acausal) vs causal filtering: how do they shift onsets?

### Takeaway
Causal filters delay a response by roughly their low-frequency group delay; for the amplifier's 4th-order 500 Hz low-pass that is about 0.83 ms, which matches the measured 0.9 ms. Zero-phase filters add no delay at half-max for a symmetric smoother, but they smear energy backward, so low thresholds and change-point onsets come out early. The ERP literature puts that error at tens of ms for 30 Hz cutoffs. Here the dominant zero-phase smoother is the Gaussian σ = 2 ms (≈ 66 Hz −3 dB), not the 200 Hz envelope low-pass. The net filter bias on the CSD-vs-MUA half-max comparison is small (≈ −0.8 ms) and is already measured and corrected in the project.

### Cited Findings
- **VanRullen 2011.** Typical 30–40 Hz low-pass filters, being zero-phase, are "mostly harmless" for peak amplitude or latency analyses but not for onsets. A neural event gets "smeared out in time for several tens or even hundreds of milliseconds before and after the event." Simulation: 50 trials of a step whose onset was uniform in 150–180 ms (noise SD 0.05), filtered at 30 Hz with EEGLAB `eegfilt` defaults; the filtered ERP departed from zero "as early as 100 ms". His rule of thumb is that the problem applies to cutoffs below ~100 Hz. Built-in amplifier filters at ≥1 kHz "could theoretically distort neuronal onset latencies by 1 or 2 ms." Remedy: analyze unfiltered data or restrict interpretation to peaks. **[full text]** VanRullen R (2011) *Front Psychol* 2:365, "Four common conceptual fallacies in mapping the time course of recognition." [doi:10.3389/fpsyg.2011.00365](https://frontiersin.org/articles/10.3389/fpsyg.2011.00365/full). Note: the search engine listed Kanai and Kreiman as co-authors, but the fetched article page gives VanRullen as sole author.
- **Rousselet 2012** (reply to VanRullen; real face-ERP data from 30 sessions, electrode B8, 500 Hz sampling) **[full text; tables garbled in extraction, verify]**:
  - In a 10 Hz sine-onset demo, non-causal filtering advanced the onset. A causal high-pass preserved the onset. A causal low-pass delayed it by 10 samples (20 ms), and a shift can undo that delay.
  - Non-causal FIR low-pass at 20/30/40/50 Hz moved median onsets earlier by only **1.7/2.0/2.5/1.5 ms** relative to raw data.
  - Non-causal high-pass at 1 and 2 Hz advanced onsets by **88 and 106 ms**.
  - Causal elliptic high-pass at 0.5–5 Hz moved median onsets by **−0.3 to −2.6 ms**.
  - A 4th-order non-causal Butterworth distorted more than FIR at low high-pass cutoffs (0.3–0.4 Hz).
  - Recommendation: "causal high-pass filtering provides an excellent solution to study ERP onsets."
  - Rousselet GA (2012) *Front Psychol* 3:131, "Does filtering preclude us from studying ERP time-courses?" [doi:10.3389/fpsyg.2012.00131](https://www.frontiersin.org/articles/10.3389/fpsyg.2012.00131/full)
- **Widmann & Schröger 2012** (commentary on both papers above; simulation at 500 Hz sampling, 30 Hz cutoff, 49-point FIR) **[full text; numbers as extracted, verify]**:
  - The EEGLAB `firls` default underestimated onset by up to ~42 ms. Better non-causal filters cut the underestimate to about **4 ms (windowed sinc) to 12 ms (discrete Gaussian)**.
  - Causal filtering overestimated onset by **4 ms** (16 ms with jack-knifing).
  - At −14 dB single-trial SNR, onsets were overestimated by 4–6 ms (non-causal) and 18 ms (causal).
  - Only linear-phase FIR filters can be made exactly zero-phase by shifting by the group delay. Minimum-phase causal filters avoid pre-onset response but add a delay that cannot be removed exactly, because their phase is non-linear.
  - Advice: use filters "as short as possible"; avoid default settings; keep cutoffs and transition bands away from the signal band.
  - Widmann A, Schröger E (2012) *Front Psychol* 3:233, "Filter effects and filter artifacts in the analysis of electrophysiological data." [doi:10.3389/fpsyg.2012.00233](https://www.frontiersin.org/articles/10.3389/fpsyg.2012.00233/full)
- **Acunzo, MacKenzie & van Rossum 2012.** Acausal high-pass filtering "can generate a systematic bias easily leading to misinterpretations of neural activity", with the early C1 component as the case study. Of 185 ERP/ERF papers surveyed, 80 used high-pass cutoffs above 0.1 Hz. A corrigendum appeared in *J Neurosci Methods* 211(2):309; its content was not checked. **[abstract/summary]** *J Neurosci Methods* 209:212–218. [doi:10.1016/j.jneumeth.2012.06.011](https://eprints.nottingham.ac.uk/49640/)
- **Widmann, Schröger & Maess 2015.** Practical guidance on FIR vs IIR filters, cutoff, order, roll-off, ripple, delay and causality, aimed at protecting SNR while limiting distortion. Secondary sources report that it favours FIR filters and discourages software defaults. **[abstract/summary only; no onset numbers verified]** *J Neurosci Methods* 250:34–46. [abstract](https://home.uni-leipzig.de/biocog/archived-20210701/index.php@page=abstract&signature=widmann_a2015jneuroscimethods250_34.html)
- **de Cheveigné & Nelken 2019.** Review of what filters do, the problems they cause and how to avoid filtering. It recommends reporting filter characteristics, including an inset plot of the impulse or step response. **[abstract/summary only; full text 403]** *Neuron* 102:280–293. [doi:10.1016/j.neuron.2019.02.039](https://cris.huji.ac.il/en/publications/filters-when-why-and-how-not-to-use-them/); author manuscript on [UCL Discovery](https://discovery-pp.ucl.ac.uk/id/eprint/10067236)
- **Project measurements [project]:**
  - The 1 kHz LFP equals the 30 kHz signal passed through the amplifier's causal 500 Hz filter at lag 0 (r = 1.000), and lags a zero-phase version by **0.9 ms**. CSD therefore carries that 0.9 ms delay; the MUA chain is zero-phase.
  - Project simulation over rises of 4–12 ms (10–90%): current half-max +0.6 ms late and peak +1 to +3 ms late; MUA half-max −0.1 to −0.2 ms and peak 0 to +2 ms. Expected MUA − current bias: **−0.8 ms** (half-max), −1 ms (peak).
  - Source: [docs/09_paper_design.md, V4 check J](../../docs/09_paper_design.md)

### Inferences
- **[own sim] Causal 4th-order Butterworth, 500 Hz (the LFP/CSD chain).**
  - Analytic DC group delay = 2.613/(2π·500) = **0.832 ms**, consistent with the measured 0.9 ms.
  - Simulated half-max delay: **+0.83 to +0.84 ms** for rises of 4–16 ms. Peak delay: +0.93 / +0.80 / +0.60 ms for rises of 4 / 8 / 16 ms.
  - A zero-phase version of the same filter shifts half-max by 0.00 ms.
  - So the CSD is late by a near-rigid ~0.83 ms, independent of shape. Spatial operations (second difference, Hamming, 3-channel mean) are linear in space and add no temporal delay.
- **[own sim] Zero-phase 4th-order 200 Hz Butterworth envelope low-pass (`sosfiltfilt`).**
  - On a unit step: 50% crossing at exactly 0 ms (by symmetry), 10% at **−1.07 ms**, 5% at **−1.28 ms**, plus a pre-ringing undershoot of −6.8% at −2.47 ms.
  - On responses with 2 / 4 / 8 ms rises: half-max shift +0.03 / +0.01 / −0.01 ms; 10% shift −0.59 / −0.24 / −0.03 ms.
  - A causal version would have delayed MUA by **+2.27 ms** at half-max (+1.27 ms at 10%). A pipeline that used a causal envelope filter would therefore bias the comparison by ~2 ms toward the MUA being *late*.
- **[own sim + analytic] Gaussian σ = 2 ms.**
  - Its −3 dB point is at √(ln 2)/(2πσ) ≈ **66 Hz** (amplitude ½ at ≈ 94 Hz). That is inside VanRullen's "<100 Hz" danger zone and makes it the strongest temporal smoother in either chain.
  - On a step it crosses 10% / 5% / 1% at −2.56 / −3.29 / −4.65 ms (−1.28σ, −1.64σ, −2.33σ).
  - On transients with 2 / 4 / 8 ms rises: half-max −0.21 / −0.21 / −0.20 ms (nearly unbiased); 10% −2.14 / −1.72 / −1.16 ms; peak **+3.30 / +2.53 / +1.47 ms**. Peaks are delayed because the fast-rise, slow-decay transient is asymmetric.
  - Because the Gaussian is applied to both signals, it biases the comparison only through shape differences. A sharper MUA transient has its peak delayed more and its 10% crossing advanced more than a broader sink.
  - Net: low-threshold onsets favour the sharper signal (MUA). Peak comparisons are biased against it, so the observed MUA peak lead is conservative on this count. Half-max is nearly neutral.
- The causal first-order 0.3 Hz high-pass (τ ≈ 531 ms) lets a step decay by only ~1.9% in 10 ms. It is irrelevant to ms-scale onsets. Rousselet's finding that causal high-pass filters preserve onsets supports this, and Acunzo's acausal high-pass problem does not apply because the CSD chain's high-pass is causal.
- The zero-phase 500–5000 Hz band-pass acts before rectification. Its envelope response is about the step response of a low-pass with half its bandwidth (general filter theory, [Analog Devices, Basic Linear Design ch. 8](https://www.analog.com/media/en/training-seminars/design-handbooks/basic-linear-design/chapter8.pdf) [summary]), i.e. sub-ms. Each spike's energy is spread symmetrically by ≲1 ms. The full-chain simulation in Q2 shows no half-max bias from it.
- The 30 kHz stream's hardware low-pass at 7.5 kHz (stated in the project docs as 0.3 Hz–7.5 kHz hardware; order not verified here) would contribute ~0.04 ms of causal delay if 3rd-order Butterworth (2/(2π·7500)). That is negligible.
- **Bottom line for Q1:** all filter effects together predict an MUA − current half-max bias of about **−0.8 ms**, almost entirely the causal 500 Hz LFP filter. That is well short of the observed −3 ms, so filtering alone cannot produce the lead. This matches the project's own estimate.

### Gaps
- The full texts of de Cheveigné & Nelken 2019 and Widmann et al. 2015 were not retrieved (403 / paywall), so neither paper's specific quantitative statements are verified here. Acunzo 2012 details and the corrigendum's content are also unverified.
- The ERP filtering literature works at ~30 Hz cutoffs and 500 Hz sampling. No source scales its errors to the 200 Hz and σ = 2 ms smoothers used here, so the transfer is by simulation only.
- One step should be confirmed in code: that MUA downsampling from 30 kHz to 1 kHz is zero-phase (scipy `decimate(zero_phase=True)`, `resample_poly`, or picking samples after the 200 Hz low-pass). This was not checked here.

---

## Q2. Analog MUA definitions and how envelope choices affect MUA latency

### Takeaway
The standard analog MUA is band-pass (~0.5–5 kHz), full-wave rectification, then low-pass (~200 Hz). Full-wave rectification keeps volume-conducted contamination low. With a zero-phase envelope filter, the 200 Hz cutoff has essentially no effect on half-max timing. A causal envelope filter would have delayed the MUA by ~2 ms. The rectified envelope is sub-linear in spike rate when spikes overlap densely, which can move MUA half-max ≲1 ms earlier than the underlying rate. That is too small to explain a 3 ms lead.

### Cited Findings
- **Legatt, Arezzo & Vaughan 1980.** Volume-conducted, time-locked activity can contaminate averaged MUA. What matters is its magnitude relative to the locally generated MUA, not its absolute size. With **full-wave** rectification, volume-conducted activity that is a sizeable fraction of the MUA level does not much distort the average. **Half-wave** rectification is discouraged because small far-field potentials then have a much larger effect. **[abstract/summary; the snippet did not show author names, attribution is via the standard citation]** *J Neurosci Methods* 2:203–217 (1980), "Averaged multiple unit activity as an estimate of phasic changes in local neuronal activity: effects of volume-conducted potentials." [PubMed 6771471](https://pubmed.ncbi.nlm.nih.gov/6771471)
- **Supèr & Roelfsema.** Multiunit recording is reliable and useful when neurons near the electrode have similar response properties, and it is advantageous for cross-correlation analysis. The chapter's exact envelope parameters were **not verified**. The Roelfsema-lab convention of 500–5000 Hz band-pass, rectification and 200 Hz low-pass is from memory, not this session. **[abstract/summary]** *Prog Brain Res* 147:263–282 (vol. 2005; online 2004), "Chronic multiunit recordings in behaving animals: advantages and limitations." [doi:10.1016/S0079-6123(04)47020-4](https://eurekamag.com/research/048/521/048521212.php); [UvA record](https://dare.uva.nl/id/fe35d97f-a374-456f-a80e-804ac605c606)
- **Kajikawa & Schroeder 2011.** In macaque A1, frequency-tuning bandwidths show that the signals' "listening areas" order as **MUA < CSD < LFP**, i.e. MUA is the most local. LFPs receive non-local contributions; a PDF copy reads "local", which is probably an extraction artifact. Their Fig. 1 shows the initial current sink and the peak MUA both in layer 4. **[abstract + figure legends]** *Neuron* 72:847–858. [PMC3240862](https://pmc.ncbi.nlm.nih.gov/articles/PMC3240862); [doi:10.1016/j.neuron.2011.09.029](https://doi.org/10.1016/j.neuron.2011.09.029)
- **Project pipeline [project]:**
  - MUA: 30 kHz stream (0.3 Hz–7.5 kHz hardware), 500–5000 Hz band-pass, rectification, 200 Hz low-pass (4th order, `sosfiltfilt`), all zero-phase, then Gaussian σ = 2 ms.
  - Clock: both streams are read on the NSP tick clock, so the legacy 3.4 ms LFP/CSD lag (docs/05 A1) is gone.
  - Source: [docs/09_paper_design.md](../../docs/09_paper_design.md); [docs/05_known_issues_and_QC.md A1](../../docs/05_known_issues_and_QC.md)

### Inferences
- **[own sim] Envelope low-pass cutoff (zero-phase).** At 200 Hz, half-max moves by ≤ 0.03 ms for rises of 2–8 ms (Q1). Pre-onset spread shows only at low thresholds (10% crossing −0.6 ms for a 2 ms rise, ≈0 for an 8 ms rise). A lower cutoff behaves like a wider symmetric smoother: still ~0 at half-max, with growing early leakage at low thresholds. With a causal envelope filter, the cutoff sets a delay of roughly the filter's group delay: +2.27 ms at half-max for a 4th-order 200 Hz filter on a step. Cutoff matters for latency mainly when the envelope filter is causal.
- **[own sim] Rectification nonlinearity.**
  - Model: Poisson spike sums of a biphasic ~1 ms waveform with 30% amplitude jitter, plus white noise, through the full MUA chain plus σ = 2 ms, 400 trials. Compared against the true rate smoothed by the same Gaussian.
  - aMUA − rate half-max: **−0.03 to +0.10 ms** (rise 4 ms) and **−0.15 to −0.54 ms** (rise 8 ms) in sparse or noise-dominated regimes. Dense overlapping spikes with near-zero baseline: **−0.41 ms** (rise 4) and **−0.94 ms** (rise 8).
  - Analytic check with envelope ∝ √(baseline + rate): baseline/peak 0 → −0.67 / −1.33 ms; 0.1 → −0.35 / −0.69 ms; 1.0 → −0.11 / −0.22 ms, for 4 / 8 ms rises.
  - Monte Carlo scatter is ~0.3–0.5 ms, judged from the sparse regime where theory predicts ~0. Peak shifts are +0.1 to +0.6 ms.
  - Conclusion: compressive envelope statistics bias MUA half-max early by **≲1 ms**, largest for dense firing from a low baseline (as at a strongly driven input layer). The bias does not scale to 3 ms.
- Full-wave rectification plus the 500 Hz high-pass keeps residual LFP-band (sink) energy out of the MUA (Legatt et al.). LFP leakage into MUA is therefore unlikely to have pulled the MUA's timing toward the sink's.
- Because MUA is the most local signal (Kajikawa & Schroeder), its timing at the sink channels is less spatially blended than the CSD's (see Q3).

### Gaps
- The full texts of the specific MUA-definition papers named in the brief were not retrieved. These include the Supèr & Roelfsema envelope parameters, Xing et al. (2012 PNAS, [pdf link](https://www.pnas.org/doi/pdf/10.1073/pnas.1201478109), 403), and Choi et al. (likely Choi, Koenig, Jia & Thakor 2010, IEEE TBME, "Quantifying time-varying multiunit neural activity using entropy based measures", [PubMed 20460201](https://pubmed.ncbi.nlm.nih.gov/20460201/); not read).
- No source was found that measures MUA onset latency as a function of envelope low-pass cutoff (targeted search returned nothing).
- No published measurement of how the analog MUA envelope scales with spike rate (√ vs linear regime) in macaque V1 was found. The simulation regimes are assumptions.

---

## Q3. CSD estimation (standard vs iCSD, spatial smoothing, volume conduction) and timing of the earliest sink

### Takeaway
None of standard CSD, iCSD or spatial smoothing introduces a *temporal* delay; all are linear in space. They do mix currents across depth. The 3-contact average at the sink, the Hamming smoothing and the second difference itself blend neighbouring currents into the "sink current". If those currents start later, the sink's half-max moves later by about +1 to +4.5 ms in simple simulations. This is the one processing factor found that is large enough to account for most of a ~3 ms MUA lead.

### Cited Findings
- **Pettersen et al. 2006 (iCSD).** Inverse CSD inverts the electrostatic forward solution and works with various array geometries, in three variants: δ-source thin discs, step-wise constant, and cubic-spline smooth. The standard CSD gives significant errors when the lateral source extent is comparable to a cortical column (≲ ~1 mm). Conductivity discontinuities can cause sizable errors even for wider sources. **[abstract]** Pettersen KH, Devor A, Ulbert I, Dale AM, Einevoll GT (2006) *J Neurosci Methods* 154:116–133. [record](https://nmr.mgh.harvard.edu/publications/journal_articles/3942). The Elephant `icsd.py` implementation (derived from Pettersen's code) offers Hamming spatial filtering, e.g. `num = ss.hamming(self.f_order)`. That is a property of the software and does not show the 2006 paper used Hamming. [Elephant icsd.py](https://gin.g-node.org/liangzai951/multielectrode_grasp/src/master/code/elephant/elephant/current_source_density_src/icsd.py)
- **Kajikawa & Schroeder 2011.** LFP reflects non-local, volume-conducted contributions, whereas CSD is more local and MUA most local (listening areas MUA < CSD < LFP, macaque A1). **[abstract]** [PMC3240862](https://pmc.ncbi.nlm.nih.gov/articles/PMC3240862)
- **Mitzdorf 1985.** The canonical review of CSD's physical basis and its interpretation in cat cortex: sinks reflect net inward transmembrane current, dominated by synaptic (EPSC) currents. The search did **not** surface the passage on presynaptic/fibre contributions, so that point is unverified. *Physiol Rev* 65:37–100. [doi:10.1152/physrev.1985.65.1.37](https://docs.backyardbrains.com/assets/files/Mitzdorf-1985-62c16807b94586a7ef1ee976c48b39ec.pdf) [bibliographic only]
- **Project pipeline [project]:** CSD is computed from the causally filtered 1 kHz LFP in space only (bad contacts interpolated, 3-point Hamming across contacts). The "sink current" is −CSD averaged over the sink channel and its two neighbours. Source: [docs/09_paper_design.md](../../docs/09_paper_design.md)

### Inferences
- **No temporal delay from spatial operations.** Second difference, Hamming smoothing, the 3-channel mean and iCSD are all linear and time-invariant across contacts. The CSD inherits exactly the LFP's temporal filtering (the causal 0.83–0.9 ms) and nothing more.
- **Spatial footprint (computed from the stated operations; the project's actual weights should be checked in code).**
  - scipy/MATLAB `hamming(3)` = [0.08, 1, 0.08], normalized [0.069, 0.862, 0.069].
  - Convolving with the 3-contact mean gives effective CSD-contact weights [0.023, 0.310, 0.333, 0.310, 0.023]. That is essentially an equal-weight average over the sink contact ±100 µm, plus ~2% from ±200 µm.
  - The second difference reaches one more contact on each side in LFP space, so the sink current depends on LFP contacts within ±300 µm.
  - If the implementation instead used classic Hamming-shaped weights such as [0.23, 0.54, 0.23], the outer weights would be roughly 10× larger.
- **[own sim] Depth-mixing bias.**
  - Model: centre current with true onset 0 and an 8 ms rise; both neighbours have the same shape but start d ms later, with amplitude w relative to the centre. The 3-contact average's half-max shifts by:
    - w = 0.5: **+0.94 / +1.83 / +2.54 ms** for d = 2 / 4 / 6 ms.
    - w = 1.0: **+1.31 / +2.71 / +4.50 ms** for d = 2 / 4 / 6 ms.
  - So averaging the earliest sink with contacts that activate 2–6 ms later can by itself make the current appear 1–4.5 ms later than at the true earliest point, i.e. most of the observed −3 ms MUA lead. The MUA from the same three channels is also averaged, but MUA is more local (Kajikawa & Schroeder) and its laminar onset gradient may differ.
  - **Test:** recompute the MUA − current difference on the single earliest-sink contact, without the 3-channel mean and with/without Hamming. Also try iCSD (spline) and compare.
- **Standard vs iCSD.** Pettersen's errors for laterally limited sources are mainly spatial (amplitude, smearing in depth, spurious sources at the edges). By the reasoning above they could move timing only through depth mixing, not as a direct delay. No source quantifies the effect on sink onset times.
- **Volume conduction.** Far-field LFP is largely removed by the second spatial derivative (Kajikawa & Schroeder). Any residual common far-field component mixes timing from distant generators into the sink estimate, by the same mechanism as depth mixing. Its size is unknown here.

### Gaps
- No paper found that quantifies how CSD method (standard vs iCSD) or spatial smoothing changes the onset latency of the earliest sink.
- The anatomical thickness of macaque L4Cα/β relative to the ±100–300 µm footprint was not sourced here.
- The Mitzdorf 1985 statements on presynaptic (afferent fibre) contributions to the earliest sink were not verified.

---

## Q4. Onset criteria (half-max, threshold, change-point) and rise-time confounds in "which comes first"

### Takeaway
Half-max of each signal's own peak confounds onset with rise time. Two signals with the same true onset but different rises differ at half-max by about half the difference in their full rise times: ≈ 2 / 3 / 4 ms for rise differences of 4 / 6 / 8 ms. A rise-time confound predicts a lead that grows with the threshold fraction. A lead that stays constant across thresholds (the project's MUA lead is −3.8 to −3.3 ms at 10–75%) points to a rigid time shift instead. Model-based estimators (maximum likelihood, change-point) are preferred in the methods literature.

### Cited Findings
- **Friedman & Priebe 1998.** Derived maximum-likelihood and least-squares latency estimators and compared them with those common in the literature. MLE performed best by mean-squared error and was recommended for automation. **[abstract]** *J Neurosci Methods* 83:185–194. [PubMed 9765132](https://pubmed.ncbi.nlm.nih.gov/9765132); [doi:10.1016/S0165-0270(98)00075-2](https://scholars.duke.edu/publication/669158). Follow-up: Friedman & Priebe 1999, "Smoothing bandwidth selection for response latency estimation," *J Neurosci Methods* 87:13–16, which addresses the smoothing step of the half-height technique. [pdf](https://www.cis.jhu.edu/~parky/CEP-Publications/FP-JNM1999.pdf) [abstract only; 403 on pdf]
- **Lee, Williford & Maunsell 2007** (*J Neurosci* 27:9632–9637, not *J Neurophysiol*) **[full text]**:
  - Latency was "time to half peak": PSTH smoothed with a Gaussian of **SD 8 ms**; peak = first maximum exceeding baseline by 3.72 SEM; latency = first time the filtered PSTH exceeded half the peak plus baseline; bootstrap 95% CIs.
  - Rationale: "a vertical offset does not lead to large differences in latency." Time-to-first-significant-change and time-to-peak would also be insensitive to the offset. A fixed-rate criterion would show a modest change.
  - Half-peak latency depended strongly on contrast (V4 population, attended): **62 ms (100%), 73 ms (25%), 96 ms (6.2%), 209 ms (1.6%)**.
  - Attention, which changed magnitude by a similar amount, did not change latency significantly (p > 0.13). So half-peak latency is not simply a function of response magnitude.
  - [PMC6672969](https://pmc.ncbi.nlm.nih.gov/articles/PMC6672969)
- **Levakova et al. 2015.** Review of neuronal response-latency estimators, covering threshold, change-point and model-based methods. **[title/listing only]** Levakova M, Tamborrino M, Ditlevsen S, Lansky P (2015) *BioSystems* 136:23–34. [listing](https://test01.ku.dk/:obvius/pureproxy/3974268/da/publications/a-review-of-the-methods-for-neuronal-response-latency-estimation)
- **Project result [project]:**
  - The MUA leads at every crossing threshold from 10 to 75% of its own peak: **−3.8, −3.5, −3.3, −3.5 ms** at 10, 25, 50, 75% (dominant eye 90%).
  - The legacy `dichoptic_csd.py` used MUA half-max vs sink 25% of peak, which favoured the sink.
  - Source: [docs/09_paper_design.md](../../docs/09_paper_design.md)

### Inferences
- **[own sim] Pure rise-time confound** (identical true onsets, raised-cosine rises, both Gaussian σ = 2 ms):

  | Rise A vs B | Half-max diff | 10% diff | Peak diff |
  |---|---|---|---|
  | 4 vs 8 ms | −2.00 ms | −1.38 ms | −2.93 ms |
  | 6 vs 12 ms | −3.00 ms | −1.79 ms | −4.77 ms |
  | 4 vs 12 ms | −4.00 ms | −2.52 ms | −6.20 ms |

  - A ~3 ms half-max "lead" can therefore arise with zero true onset difference if the current's full rise is ~6 ms longer than the MUA's (10–90% ≈ 3.5 ms longer).
  - **Diagnostic:** under a rise-time confound, the lead shrinks toward low thresholds (about 60–70% of the half-max value at 10%) and grows toward the peak.
- **Threshold invariance argues against rise time.** The project's lead is nearly constant from 10% to 75% (−3.8 to −3.3 ms), with a somewhat larger peak lead (−4.5 ms). The near-invariance across thresholds argues against a rise-time difference as the main cause and fits a near-rigid shift of the rising edge. Candidate rigid shifts:
  - The causal LFP filter (0.83–0.9 ms; already corrected).
  - A clock offset (the legacy 3.4 ms; QC says it is gone, r = 1.000 at lag 0).
  - A true biological lead.
  - Depth mixing (Q3). Its threshold profile was not computed here.
- **Baseline and amplitude.** Half-max of each signal's own peak is invariant to scale and baseline offset (Lee et al.), so the CSD's and MUA's different units and SNRs do not by themselves bias it. Fixed-threshold criteria (e.g., baseline + k·SD) do depend on SNR: the noisier signal crosses later, and the bias grows with lower SNR. Widmann & Schröger 2012 also showed onset errors grow at low SNR. Such criteria would bias "which comes first" toward the higher-SNR signal.
- **Mixed criteria.** Using different fractions for the two signals biases the comparison. For a raised-cosine rise of R ms, the 25% crossing comes R/6 ms before the 50% crossing (crossings at R·acos(1−2f)/π). So MUA at 50% vs sink at 25% favours the sink by ≈ 0.7–2 ms for full rises of 4–12 ms, before smoothing. This is the legacy bias already flagged.
- **Change-point / MLE estimators** (Friedman & Priebe) estimate the departure from baseline, not the mid-rise. They would separate "onset" from "rise". They are sensitive to backward smearing by zero-phase smoothers, though: the Gaussian σ = 2 ms makes a step depart from baseline ~2.5–4.6 ms early (Q1). Change-point methods should therefore be run on data without the σ = 2 ms Gaussian, or with a causal smoother whose delay is corrected.

### Gaps
- The Commenges & Seal change-point papers (1985/1986) were not retrieved, so their properties are not documented here.
- The full text of Friedman & Priebe 1998 was blocked (403), so its numerical bias/variance tables are not reported.
- No source was found that gives a closed-form bias of half-max latency as a function of rise time. The table above is from simulation only.

---

## Q5. Papers that explicitly compared CSD and MUA onset latencies (ideally under different processing pipelines)

### Takeaway
No paper was found that compares CSD-sink and MUA onset latencies across processing pipelines. The classic macaque and A1 laminar work describes the earliest sink and the MUA as starting together at the thalamorecipient layer. That work also attributes part of the earliest response to the thalamocortical afferent (presynaptic) volley, which is a plausible biological reason for MUA at the input layer to lead the summed synaptic current.

### Cited Findings
- **Schroeder, Mehta & Givre 1998.** In awake macaque, laminar CSD indexed postsynaptic potentials and MUA indexed action potentials. V1 and V2 showed a "feedforward" profile, with excitation beginning at the depth of lamina 4 and spreading to extragranular layers. **[abstract]** *Cereb Cortex* 8:575–592, "A spatiotemporal profile of visual system activation revealed by current source density analysis in the awake macaque." [PubMed 9823479](https://pubmed.ncbi.nlm.nih.gov/9823479)
- **Givre, Arezzo & Schroeder** (*Vis Neurosci*, wavelength study, macaque V1, CSD + MUA). "In layer 4c, the red effect reached significance, on average, at 47 ms, or ≈24 ms after the onset of transmembrane current flow." That implies a 4C current-onset latency of ≈23 ms for their flash stimuli (inference by subtraction). **[abstract/summary]** [Cambridge abstract](https://www.cambridge.org/core/journals/visual-neuroscience/article/abs/effects-of-wavelength-on-the-timing-and-laminar-distribution-of-illuminanceevoked-activity-in-macaque-v1/A6F43C4DE063E98F41A3A3526EF1849D). A secondary review notes that Schroeder et al. and Givre et al. define onset latency differently, which affects cross-area latency comparisons. [Ales et al., NeuroImage, pdf](https://www.st-andrews.ac.uk/~jma23/papers/alesInPressNeuroimage_OnDeterminingSources.pdf) [summary; not verified]
- **A1 laminar work (Steinschneider/Schroeder lineage).** "A comparison of concurrent AEP, CSD and multiunit activity (MUA) implicates direct contributions to the initial negativity from the thalamocortical afferent volley, immediately followed by postsynaptic activation." The early lamina IV sinks correspond "in timing and location to multiple high amplitude MUA." **[search excerpt shown with Tenke & Kayser 2012, Clin Neurophysiol, "Generator localization by current source density (CSD): implications of volume conduction and field closure at intracranial and scalp resolutions"; which primary paper the sentence summarizes was not verified]** [PMC3498576](https://pmc.ncbi.nlm.nih.gov/articles/PMC3498576/)
- **Kajikawa & Schroeder 2011.** For a best-frequency tone in macaque A1, the initial current sink and the peak MUA are both in layer 4 (Fig. 1). **[figure legend]** [PMC3240862](https://pmc.ncbi.nlm.nih.gov/articles/PMC3240862/figure/F1)
- **Rat V1.** A search summary stated that a deep-layer CSD sink, "highlighting the arrival of synaptic input", was "always followed after 10 ms by a MUA onset at same depths". **[unverified; the source page appears to concern spontaneous Up-state onsets in anesthetized rats, not visual responses]** [PMC3704866](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3704866/)
- **Project [project]:**
  - Under one pipeline, with identical half-max and peak estimators, spiking on the sink channels starts before the sink current: half-max MUA − current −3.3 ms raw / −2.5 ms bias-corrected (95% CI −5.1 to −0.5), dominant eye 90%; −3.0 / −2.2 ms (CI −3.7 to −0.3), both eyes 45%.
  - Peaks: −4.5 / −3.5 ms (CI −6 to +1) and −3.0 / −2.0 ms (CI −10 to +2).
  - The project's working hypothesis: analog MUA at the input depth includes thalamic afferent spikes; the sink grows as input accumulates; and the 3-contact average and smoothing blend in later currents.
  - Legacy CSD-vs-MUA latencies carry the 3.4 ms clock artefact.
  - Source: [docs/09_paper_design.md](../../docs/09_paper_design.md); [docs/05 A1, A7](../../docs/05_known_issues_and_QC.md)

### Inferences
- **Bias budget for "MUA half-max ~3 ms before the sink current", from the sources and simulations above.** Negative means it pushes MUA earlier than the current.

  | Source | Bias (MUA − current) | Status |
  |---|---|---|
  | Causal 500 Hz LFP filter | ≈ −0.83 ms | Measured 0.9 ms; corrected |
  | Zero-phase envelope 200 Hz | ≈ 0 ms | |
  | Gaussian σ = 2 ms | ≈ 0 ms at half-max for both | Matters for 10% / change-point / peak |
  | Envelope compressive nonlinearity | 0 to ≈ −1 ms | Regime-dependent |
  | Depth mixing in 3-contact current | 0 to ≈ −4.5 ms | Unmeasured; plausibly the largest |
  | Rise-time confound | 0 to −4 ms | Disfavoured by threshold invariance |
  | Clock | 0 | QC r = 1.000 at lag 0; legacy would have added −3.4 ms |

  After the filter correction, the residual −2.2 to −2.5 ms is within reach of depth mixing alone, plus at most ~1 ms of envelope compression. It is also consistent with a genuine afferent-spike lead. The claim "spiking precedes the sink current" is therefore defensible only as "under this pipeline". The critical control is the single-contact, unsmoothed (or iCSD) sink.
- **Physiological prior.** Postsynaptic spikes should follow the synaptic current that drives them. A real MUA lead at the input layer would therefore imply presynaptic (thalamocortical axon/terminal) spikes in the MUA, or early-firing neurons outside the averaged current's footprint. The A1 statements about the afferent volley support the first option qualitatively. No source found here quantifies the afferent-volley share of L4C analog MUA in macaque V1.
- **Recommended robustness checks, each mapped to a factor above:**
  1. Single earliest-sink contact vs the 3-contact mean; Hamming on vs off; iCSD (spline) vs standard CSD.
  2. Repeat without the σ = 2 ms Gaussian, using raw 1 kHz signals with trial bootstrap.
  3. A causal 500 Hz filter applied to the MUA as well, so both chains match, instead of subtracting a simulated bias.
  4. Thresholds from 5 to 90% of peak plus a change-point/MLE onset, checking whether the lead is rigid or grows with threshold.
  5. MUA with a "spike-rate-linear" envelope (e.g., power rather than |x|), or thresholded spike counts, to bound the rectification nonlinearity.
  6. Report the filter step/impulse responses, as de Cheveigné & Nelken recommend.

### Gaps
- No study was found that varied the processing pipeline and reported how CSD-vs-MUA onset order changed. This is a genuine gap in the located literature, not just an access failure.
- The primary Steinschneider et al. (1992/1994) and Givre et al. (1994/1995) full texts with exact sink-vs-MUA onset numbers were not retrieved.
- Self et al. 2013, van Kerkoerle et al. 2014, Maier et al. 2011 and Xing et al. 2012 report laminar CSD and MUA timing in macaque V1. They were not read this session and may contain directly comparable onset numbers.
- No source was found on the fraction of L4C analog MUA that comes from LGN afferent terminals.
