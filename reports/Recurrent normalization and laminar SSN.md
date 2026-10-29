# Offset kinetics point beyond feedforward delayed normalization

> **Correction (2026-10-29).** Two findings after this report weaken its offset argument (docs/09, "Correction",
> and register #33). First, the data's 27–32 ms "half-decay" mostly measures an OFF transient: after the response
> offset the population response rises to ~1.8× its sustained level before falling. Second, the toy scan measured
> half-decay from the post-offset peak, not as the data were measured; under the data's metric the toy ORGaNICs does
> not reproduce the data. The rise-time result (the fitted DN rises too fast at low contrast) stands. Read the
> offset claims below with this in mind.

Normalization in awake macaque V1 does take time. In the team's data it starts about 20–35 ms after each channel's own response onset. What the new dynamics diagnostics add is that normalization changes **how fast the response integrates**, which a delayed divisor applied to a fast response does not do. Rise times halve with contrast (36 → 18 ms, 10–90%). Offsets decay with a nearly contrast-independent half-decay of **27–32 ms**. The fitted feedforward delayed-normalization (DN) model predicts 20 → 11 ms rise times and 5–12 ms offsets. Conductance and recurrent theories predict this rise/decay asymmetry, and only they do. Carandini, Heeger & Movshon measured membrane time constants falling from 25 to 4.9 ms with contrast. Sit et al. found a falling edge that did not depend on contrast in macaque voltage-sensitive-dye imaging. Heeger & Zemlianova's ORGaNICs circuit predicts an effective time constant falling from ~60 to ~6 ms and an offset time constant of ~60 ms "regardless of stimulus". A toy parameter scan run for this report (derived) found that ORGaNICs, recurrent shunting and a two-population stabilized supralinear network (SSN) can each reproduce the team's rise and decay numbers when given a saturating front end, and DN cannot. Two further derivations narrow the field. First, any circuit driven only by the summed eye signal gives an equivalent-contrast factor of exactly 2 at every moment. The binocular collapse therefore needs eye-specific units whose pool weights the other eye more than their excitation does, which is the team's w > a. Second, timing that is uniform across depth when measured from each channel's own onset is compatible with local normalization in every layer. It does not match any published surround-suppression laminar sequence, so stimulus-locked timing must be checked. In the circuit literature, the team's 20–35 ms sits between fast feedforward suppression (~0–10 ms) and feedback-like suppression (45–61 ms), in the range of local recurrent or lateral inhibition. That fits cat evidence that interocular suppression is intracortical and depends on GABA-A. The next model to fit is a binocular ORGaNICs circuit with eye-specific drive, a cross-eye pool weight and a loop delay, and with a front end fixed from the earliest current sink. It should be nested with the current DN and pre-registered to predict offsets it was not fitted to. Heeger's remark has no written source. The closest published statement is that "the recurrence takes time (multiples of the time constant)".

Throughout, **(derived)** marks algebra or simulation done for this report or the underlying notes; these are not published results. **(unverified)** marks numbers the notes could not confirm against a primary source. Preparations are tagged mac-W (awake macaque), mac-A (anesthetized macaque), cat-A, mouse or human. "Team" marks the unpublished findings given in the brief, numbered F1–F6.

## 1. Every candidate theory speeds up with contrast, but they differ at offset

Four families of dynamic normalization compete. They differ in where the contrast-dependent speed sits. In conductance (shunting) models it sits in the membrane time constant. In ORGaNICs it sits in recurrent gain. In delayed normalization it sits only in when a separate pool begins to divide. In the SSN it sits in the effective coupling of an E/I network. All four normalize sooner at high contrast, so the faster collapse at high contrast (F2) does not separate them. They part ways on whether the response's own integration time changes and on how the response decays after the stimulus ends.

### Conductance models make the time constant itself the normalizer

In the original account, simple cells sum LGN input linearly, and that sum is divided by pooled cortical activity. A closed-form version predicted macaque V1 responses across contrast, orientation and spatiotemporal frequency ([Carandini & Heeger 1994](https://doi.org/10.1126/science.8191289)). Zhou et al. later classed that closed form as **instantaneous** gain control that "misses the transient-decay shape" ([Zhou et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/)).

The explicit dynamic model is a membrane equation, `C dV/dt + g(t)·V = I_d(t)`, with output `R ∝ max(0,V)^n`. The conductance depends on the pool's output responses, `g = g0/(1 − kΣR)` (reconstructed from extraction), so it is feedback normalization. The authors themselves call this form "completely ad hoc" ([Carandini, Heeger & Movshon 1997 · mac-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)). Because "the increase in membrane conductance decreases the time constant", the fitted population medians fall from **τ0 = 25 ms at rest to τ1 = 4.9 ms at unit contrast for gratings**. The corresponding values are 51 → 5.4 ms for plaids and 150 → 7.5 ms for noise masking, with n ≈ 2.4–2.6. One example cell (37 → 9 ms) "responds ∼20 msec sooner" at unit contrast than at 0.12. The authors judged the 18-fold conductance increase that noise masking would need "unlikely to be possible in real cells".

Sit et al. applied the same structure to macaque population responses: `C dV/dt = A − gV` with `g = g0(1 + B)` and `B = b·I⊗H`. The model has two stages, n = 2 and a 20 ms delay between the stages ([Sit et al. 2009 · mac-W VSD](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/)). In the data, rising-edge latency fell and rising slope grew with contrast, while "the falling edge is largely independent of contrast" (p > 0.15) and slower than the rise (t10 = 65 ms). The model assumes conductance follows the input instantaneously, and a 10 ms conductance time constant made "no qualitative difference".

One consequence of the membrane equation is an early linear window (derived). Starting from rest, the initial slope of `C dV/dt = I − gV` is I/C, which does not depend on g. The two eyes' currents therefore add linearly until roughly t ≈ C/g. At the CHM medians, that window shrinks from about 25 ms near zero contrast to about 5 ms at unit contrast.

### ORGaNICs puts the contrast-dependent speed into recurrent gain

Heeger framed normalization from the start as striate cells that "mutually inhibit one another" ([Heeger 1992](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1017/S0952523800009640%22&resultType=core&format=json)). The oscillatory recurrent gated neural integrator circuits (ORGaNICs) make this explicit through modulators ([Heeger & Mackey 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6842604/)). A recurrent modulator a and an input modulator b "determine the state of the circuit by controlling the recurrent gain and effective time constant". When both modulators are large, responses follow the input drive and have a short effective time constant. When they are small, responses are dominated by recurrence and the time constant is long. Oscillations "emerged because of the inherent delay" in the multiplicative normalization loop.

The V1 normalization circuit is the key paper ([Heeger & Zemlianova 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)). In single-unit form its equations are:

```
τ_v dv/dt = −v + [b0/(1+b0)]·z + ŷ/(1+a)        y = [v]₊²,  ŷ = √y
τ_a da/dt = −a + √u + a·√u
τ_u du/dt = −u + Σ_k w_jk y_k u_k + (σ·b0/(1+b0))²
steady state:  y = z² / (σ² + W z²)        Eq. 10:  τ = τ_v · [(1+b0)/b0] · √g,  g = 1/(σ² + c_t² + β c_m²)
```

These equations are verified for the single-unit case. They match Heeger's own MATLAB implementation, `BM-ephys/toolboxes/ORGaNIC_Normalization/organicNormalization3.m` (DJH 2018), which sets the floor `uMin = (σ·b0/(1+b0))²`. That is the constant the web extraction had garbled. The published parameters are b0 = 0.2, σ = 0.1, τ_v = 1 ms, τ_a = 2 ms and τ_u = 1 ms (10 ms for the onset-transient simulation). The authors call these intrinsic constants "∼10 times shorter than experimental measurements", and state that the simulations "have not (yet) been optimized to replicate published data".

The published predictions are specific. The time constant is ≈ 60 ms at zero contrast and ≈ 6 ms at high contrast, so "high-contrast stimuli are integrated over much briefer periods of time (by a factor of ∼10×)", while responses are amplified about 100× at low contrast and only about 1× at high contrast. The **offset time constant is ≈60 ms regardless of stimulus**, so responses rise faster than they decay. Slow modulators produce onset transients; with fast modulators the fixed point can become unstable and the circuit oscillates. Gamma at about 40–50 Hz appears only above ~50% contrast and its frequency rises with contrast, although the simulated rise of gamma frequency with stimulus size contradicts data, as the authors acknowledge.

The paper also contains Heeger's most direct published version of "normalization takes time": "the recurrence takes time (multiples of the time constant)". A related Heeger framework finds that recurrent dynamics delay responses by about 100 ms between successive levels of its hierarchy even with τ = 5 ms ([Heeger 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5338385/)).

Two derived facts follow from these equations. First, linearizing with the modulators held fixed gives τ_eff = τ_v/(B·√(σ² + Σ_k w_jk z_k²)), where B = b0/(1+b0). The speed therefore depends on **total drive to the pool, which in a binocular circuit includes the other eye**. Second, Rawat et al. rewrote the circuit as `τ_y ẏ = −y + b z + (1 − √a)·y` and `τ_a ȧ = −a + b0²σ² + W(y² ⊙ a)` ([Rawat, Heeger & Martiniani 2024](https://arxiv.org/html/2409.18946v3)). They proved unconditional local stability for identity recurrent weights, and in two dimensions for 0 < w_r ≤ 1. In this form τ_eff = τ_y/√a. Holding y fixed, the normalization signal relaxes with a time constant that grows with contrast, and it is self-exciting. A toy simulation in the notes produced **transient over-normalization** that the DN and shunting models did not (derived). At c = 0.45 the binocular summation ratio was 0.25 at +10 ms against a sustained 0.52. This depends on parameters and has not been checked against the coupled eigenvalues in Rawat et al.'s appendix.

A 2025 follow-up proposes that a breakdown of normalization is "an early warning signal for the loss of stability" ([Morone et al. 2025, preprint](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12139785/)).

Blake's legacy ORGaNICs fits are in `BM-ephys/3. visualization/2024_LaminarCRFs/Manuscript_focused_v2/Simulation/`. They used Bayesian optimization over monocular 0–200 ms responses, with bounds σ 0.05–0.1, b0 0.1–0.3, τ_v 0.3–2 ms, τ_a 0.4–5 ms and τ_u 7–15 ms. Under Eq. 10, with input linear in contrast, σ ≤ 0.1 forces τ_eff(0.11)/τ_eff(0.9) = √((σ² + 0.81)/(σ² + 0.0121)) ≥ 6.1 (derived). The team's rise-time ratio over the same contrasts is 2. The legacy bounds therefore cannot reach the data unless σ is freed or the input saturates.

### Delayed normalization divides a fast response by a slow pool

The main feedforward competitor divides the linear response by a low-passed copy of itself ([Zhou et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/)): `R_L = S ∗ h1`, `R_DN = |R_L|^n / (σ^n + (|R_L| ∗ h2)^n)`, with `h2 = e^(−t/τ2)`. One parameter set fit all 10 contrasts of each cell, so **τ2 does not change with contrast**. Lower, later peaks at low contrast come from normalization building up later when σ dominates. Fitted τ1 was "~40 ms for the single cells". Scaled and shifted copies of the 70% response explained more than 98% of the variance at other contrasts. The model beat instantaneous normalization, a two-channel model and Carandini & Heeger 1994. In human ECoG, "delayed normalization is critical to predict the slower dynamics associated with changes in stimulus contrast", and peak latency "sharply decreases with increasing contrast" ([Groen et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9546476/)).

Independent fits agree that normalization is "rapid, but not instantaneous". SSVEP masking gave a short-memory pool with τ = **26 ms** (95% CI 16–36) ([Tsai, Wade & Norcia 2012 · human](https://pmc.ncbi.nlm.nih.gov/articles/PMC3337145/)). In voltage-sensitive-dye imaging, population contrast gain rose with τ ≈ 30 ms: c50 fell from 30% to 5% over 40–100 ms (τ = 31 ms), and the exponent fell from ~3 to below 1 (τ = 35 ms) ([Reynaud, Masson & Chavane 2012 · mac-W](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242/)). Across spikes, LFP, ECoG and EEG, normalization integrated over about 14–23 ms for parallel masks and 5–10 ms for orthogonal masks, and the authors call a single filter "unlikely" to suffice ([Gulati & Ray 2025 · bonnet macaque, preprint](https://www.biorxiv.org/content/10.1101/2025.03.11.642541)). The team's ~30 ms pool sits inside this range.

Three things follow for DN (derived). It collapses sooner at high contrast because the pool crosses σ at t* ≈ τ2·σ/c. Its kinetics are still fixed by τ2, and its depth is capped at the static asymptote. Its offset is limited by the numerator, so it decays about as fast as the front end and may undershoot.

One objection to any delay is that contrast-set gain control is "fully expressed within the first 10 ms" once responses are aligned for latency ([Albrecht et al. 2002 · mac-A, cat-A](https://pubmed.ncbi.nlm.nih.gov/12163540/)). That does not refute DN (derived). Aligning "scaled and shifted" responses removes exactly the contrast-dependent latency that DN generates, and Albrecht's result concerns monocular self-normalization, which is plausibly inherited from LGN saturation.

Contrast-dependent integration is also seen directly. In direction-selective V1 cells, the mean temporal cutoff was **10 Hz at low contrast against 21 Hz at high contrast** ([Bair & Movshon 2004 · mac-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC6729763/)).

### Rate models with an instantaneous target cannot open a linear window

The temporal-attention normalization model relaxes toward an algebraically normalized target: `τ dr/dt = −r + e/(s + σ^n)`, `e = a(w·x)^n`, `s = Σe`, with n = 1.5 and τ = 52 ms. τ only "determines how long the response takes to rise to steady state … and to return to zero" ([Denison, Carrasco & Heeger 2021](https://www.cns.nyu.edu/heegerlab/content/publications/Denison-NatHumBehav2021.pdf)). Its spatiotemporal successor adds exponential excitation and suppression windows and concludes that "suppression always acts later and across longer time intervals than excitation" ([Chapman & Denison 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12680323/)).

The Heeger-lab binocular models have the same structure. The rivalry model uses a dynamical normalization with τ = 50 ms and subtractive ocular-opponency units, simulated only at 0.5 contrast ([Said & Heeger 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3610603/)). The attention model of rivalry is already two-stage. A monocular layer with τ_s = 10 ms is normalized by **both eyes' drives**. It feeds a binocular summation layer, with opponency (τ = 20 ms), attention (τ = 150 ms) and adaptation (τ = 2 s) on top ([Li et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544261/)).

Because the pool is computed from instantaneous drives, these models apply interocular suppression in full proportion from the first moment (derived). In the notes' toy simulation the binocular summation ratio at +2 ms was already 0.66 at high contrast, against a sustained 0.50. To match F1, the pool must be computed from delayed signals or from outputs.

### The SSN normalizes by gaining coupling, and no laminar SSN predicts timing

The SSN is an E/I rate network with power-law units. Rubin et al. use `τ_E ṙ_E = −r_E + k[W_EE r_E − W_EI r_I + h_E]₊^n` with k = 0.04 and n = 2 ([Rubin, Van Hooser & Miller 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)). An L2/3 model uses n = 2.2, τ_E = 10 ms and τ_I = 6.67 ms ([Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025)). Responses to weak inputs sum supralinearly and responses to strong inputs sublinearly. Both E and I cells sum sublinearly, as "an emergent outcome of network dynamics". Equal-strength pairs give "a sublinear multiple (~0.5 to 0.7) of the sum", and unequal pairs approach winner-take-all.

Supersaturation occurs if and only if Ω_E = J_II·g_E − J_EI·g_I < 0 ([Ahmadian, Rubin & Miller 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)). The strongest sublinearity occurs for supersaturating parameters, yet most parameter sets leave "a broad dynamic range between the onset of normalization and of supersaturation". A ring-model criterion is closely related (Persi et al. 2011, cited within Ahmadian et al.).

Dynamically, the effective gain G ∝ r^((n−1)/n) grows with drive (derived). Stronger drive therefore makes the E–I eigenvalues more negative (faster) and adds ringing at higher frequencies. Simulated timescales were "nearly as fast as the single-neuron time constant (20 ms in these simulations)", against 4–15× longer in chaotic networks. At offset, activity "rapidly decayed back to spontaneous levels with the cellular time constant" ([Hennequin et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5971207/)). Loose balance can give effective time constants "as small as a few milliseconds" ([Ahmadian & Miller 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9129095/)).

Several network signatures scale with contrast. Gamma peak frequency rises with contrast while spontaneous activity "generically does not oscillate", and only a columnar model, with inhibition confined to ~0.4 mm, reproduces both the contrast dependence of gamma and surround suppression ([Holt, Miller & Ahmadian 2024](https://doi.org/10.1371/journal.pcbi.1012190)). A sampling-trained E/I network shows onset transients and gamma frequency that "scaled with stimulus contrast" ([Echeveste et al. 2020](https://doi.org/10.1038/s41593-020-0671-1)). An inhibition-stabilized network (ISN) shows a transient rise in inhibition of **30–50 ms** before suppression settles ([Ozeki et al. 2009 · cat-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725/)). In the L2/3 SSN, by contrast, E activity decays with τ ≈ 10 ms after thalamic input is removed, "regardless of stimulus size or contrast" ([Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025)).

No published SSN has interacting L4, L2/3 and L5/6 populations with laminar timing predictions. The canonical model treats "a single layer, e.g. layer 2/3 … ignoring interlaminar processing" ([Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)). A four-cell-type L2/3 SSN driven by measured L4 input switches from a positive PC→VIP⊣SST⊣PC loop at low contrast to a negative PC→SST⊣PC loop at high contrast ([Mossing et al. 2021, preprint](https://doi.org/10.1101/2021.03.31.437953)). A laminar neural-mass model built from Allen connectivity puts ISN operation in L4 and L5 but not L2/3 ([Eskikand et al. 2022, preprint](https://doi.org/10.1101/2022.12.09.519705)). Mouse optogenetics finds ISN operation in "both superficial and deep layers" ([Sanzeni et al. 2020](https://doi.org/10.7554/eLife.54875)), and a 2026 analysis warns that in coupled networks a paradoxical response "is neither necessary nor sufficient for local inhibition stabilization" ([Wu & Miller 2026, preprint](https://doi.org/10.64898/2026.09.09.750483)).

On the binocular side there is one dedicated binocular SSN. It adds ocular-dominance columns and reproduces interocular transfer of surround suppression ([Selby & Tripp 2017](https://www.sciencedirect.com/science/article/abs/pii/S0893608017300552)), but in excitatory and not inhibitory units ([Selby 2016 thesis](http://hdl.handle.net/10012/10264)). In mouse binocular V1, PV binocular summation is sublinear and tracks pyramidal cells, while SST summation is "linear or supralinear" with rapid dynamics ([Perumal et al. 2026, preprint](https://doi.org/10.64898/2026.09.20.753031)).

The SSN and every other family face one strong constraint, which comes from the team's metric (derived). Suppose both eyes feed the same units through identical temporal filters, so that the network is driven by h(t) = (c_L + c_R)·g(t). Then R_bin(c,c)(t) = R_mono(2c)(t) at every t, whatever the nonlinearity, recurrence or delays. **The equivalent-contrast factor is then exactly 2 at all times and contrasts.** Such a network produces "sublinear summation" of rates, but never the fall of k below 2 that defines F1 → F2. A binocular collapse in that metric requires eye-specific units whose suppression weights the eyes differently from their excitation.

The team's DN parametrizes exactly this (derived). Early on, before the pool builds, k = 1 + a, which is ≈ 2 for balanced channels. Once the pool builds, the collapse requires w > a, and w > a held in 24 of 25 penetrations. For the SSN, this points to the ocular-dominance-column architecture rather than the unstructured two-population loop.

### What the team's own kinetics imply

The team's diagnostics can be converted to first-order time constants (derived). Here 10–90% rise = 2.2τ and half-decay = 0.69τ, with the front-end filter ignored. Rises of 36 → 18 ms correspond to τ_on ≈ 16 → 8 ms. Half-decays of 27–32 ms correspond to τ_off ≈ 39–46 ms. **Decay is 2.4–5.6× slower than rise.** That is the asymmetry predicted by shunting and ORGaNICs. It places τ_off near the CHM resting τ0 (25–51 ms for gratings and plaids) and below the H&Z zero-contrast 60 ms. It is above the 10–20 ms cellular τ used in SSN papers. Real responses are not first-order and include the front-end filter, so these are upper bounds on intrinsic τ_on.

To check whether any model class can hit these numbers, a coarse toy scan was run for this report (derived). It is a single unit with a saturating per-eye front end f(c) = c/(c + c50_in), two low-pass input stages and a stimulus offset at 200 ms. The scripts are `research_notes/Recurrent normalization and laminar SSN/sims/rise_decay.py` and `rise_decay_scan.py`. The scan is illustrative, it was not fitted to data, and it did not test the binocular collapse jointly.

| Model class (best on a coarse grid) | Rise 10–90% at c = .11 / .45 / .9 (ms) | Half-decay at c = .11 / .45 / .9 (ms) |
|---|---|---|
| **Team data** | **36 / — / 18** | **27–32, ~flat** |
| Team's fitted DN | 20 / — / 11 | 5–12 |
| DN (τ2 = 15 ms, σ = .3, front end 2×15 ms) | 32.1 / 19.2 / 17.0 | 17.6 / 19.5 / 19.9 |
| Recurrent shunting (τ0 = 60 ms, conductance driven by output) | 34.5 / 22.1 / 19.8 | 29.0 / 25.3 / 24.5 |
| ORGaNICs, Rawat form (τ_y = 8, τ_a = 2 ms, σ = .2) | 38.0 / 19.3 / 16.7 | 30.0 / 30.7 / 32.2 |
| Two-population SSN (τ_E = 10, τ_I = 5 ms; illustrative weights) | 31.2 / 20.5 / 18.9 | 29.0 / 29.1 / 29.1 |

DN reached even 17–20 ms of half-decay only with the slowest front end tested. In DN, rise and decay are tied to the same front-end filter, and the pool shortens the rise without slowing the decay. The three models that act on the integration time constant each matched both diagnostics. **Finding F4 therefore separates "divide a fast response" from "change the integration time", but it does not choose among ORGaNICs, recurrent shunting and the SSN.** The SSN match also needed a fairly slow front end, and the scan did not check whether that parameter set is in the inhibition-stabilized, normalizing regime.

| Family | Where normalization acts | Speed vs contrast | Early linear window | Offset | Key published numbers |
|---|---|---|---|---|---|
| Static normalization | Algebraic division | Input filter only | Only while drive < σ (derived) | As fast as input | [CH 1994](https://doi.org/10.1126/science.8191289) |
| Instantaneous-target rate | Division inside a relaxing target | One fixed τ | None (toy) | Same fixed τ | τ 52 ms ([Denison 2021](https://www.cns.nyu.edu/heegerlab/content/publications/Denison-NatHumBehav2021.pdf)); 10 ms monocular layer ([Li 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544261/)) |
| Feedforward DN | Division by low-passed pool | Fixed τ1, τ2; onset ~τ2σ/c (derived) | Yes | Numerator-limited, fast | τ 26 ms ([Tsai 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3337145/)); ~31 ms ([Reynaud 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242/)) |
| Conductance / shunting | Membrane τ = C/g | 25 → 4.9 ms (gratings) | ~C/g (derived) | Slow, contrast-independent | [CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/); [Sit 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/) |
| ORGaNICs | Recurrent gain via output-driven modulators | 60 → 6 ms, ∝ 1/√(σ²+Σwz²) | First few τ | ≈60 ms regardless of stimulus | [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/) |
| SSN / ISN | Effective coupling ∝ r^((n−1)/n) | Faster with drive; ~cellular τ or less | Yes, low gain at onset | Cellular τ (~10 ms), contrast-independent | [Hennequin 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5971207/); [Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025) |

## 2. Only integration-time models with eye-specific units fit all six findings

The table maps each model against F1–F6. **Explains** means the model reproduces the finding for a stated reason. **Partly** means it gets the sign or only part of the effect. **Contradicts** means it predicts the opposite or a number the data exclude. **Silent** means it makes no prediction. The model-fit outcomes in column F3 are the team's. Every other entry is an inference from the cited equations unless it carries a source.

| Model | F1 early linear, then normalized at 20–35 ms | F2 further/faster collapse at high c | F3 fit outcomes | F4 contrast-dependent rise, slow flat offset | F5 depth-uniform timing; superficial less linear at start | F6 no supersaturation |
|---|---|---|---|---|---|---|
| Static normalization ([CH 1994](https://doi.org/10.1126/science.8191289); [MH 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)) | **Contradicts:** division is proportional from the first spike | **Partly:** with w > a, k falls ~1/c once saturated (derived), but nothing time-dependent | **Contradicts:** loses to DN in 23/25 | **Contradicts:** kinetics belong to the input filter | Silent | **Explains:** saturation is built in |
| Feedforward DN, fixed pool ([Zhou 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/); team) | **Explains:** pool starts at zero; early k = 1 + a (derived) | **Partly:** sooner (τ2σ/c), but kinetics fixed and depth capped | The reference model: beats static, needs w > a, too slow and shallow at high c | **Contradicts:** predicts 20 → 11 ms rise and 5–12 ms offset | **Explains** if the pool is shared | **Explains** |
| DN with contrast-scaled pool τ (team) | Explains | Intended to | **Failed** (team) | **Contradicts:** moves when division starts, not the response's own rise or decay (derived) | Silent | Explains |
| Dynamic two-stage, early interocular gain ([Meese et al. 2006](https://doi.org/10.1167/6.11.7); [Li 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544261/)) | **Explains** with a ~31 ms delay | **Partly:** static two-stage k falls ~1.19 → 1.02 (derived) | Beats single-stage by only ~1% held-out error, so architecture is not the bottleneck | **Contradicts** if both stages divide fast responses | Silent | Explains |
| Dynamic Ding–Sperling ([Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)) | **Explains:** combination of 50 ms exposures is "well approximated by … linear addition" | **Partly:** static k falls ~1.31 → 1.02 (derived) | **Fails:** parameters hit bounds; its psychophysical gain filter of ~110 ms is far slower than a 20–35 ms collapse | **Contradicts:** division of a fast drive | Silent | Explains |
| Instantaneous-target rate models ([Said & Heeger 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3610603/); [Li 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544261/); [Denison 2021](https://www.cns.nyu.edu/heegerlab/content/publications/Denison-NatHumBehav2021.pdf)) | **Contradicts:** suppression proportional from onset (toy) | Partly | Static plus low-pass, so no better than static (derived) | **Contradicts:** one τ, rise and decay symmetric and contrast-independent | Silent | Explains |
| Feedforward conductance ([Sit 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/); CH94 conductance form) | **Partly:** linear for ~C/g, which is only ~5 ms at high c unless conductance is delayed | **Explains:** g grows with both eyes' drive | **Explains** why scaling the pool τ fails | **Partly:** Sit saw a flat, slow falling edge, but in the toy decay lengthened with contrast | Shared pool OK | Explains |
| Recurrent conductance ([CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)) | **Explains:** the pool needs output first | **Explains:** τ 25 → 4.9 ms | Explains | **Explains:** toy rise 34 → 20 ms, decay 25–29 ms | **Explains** if shared | Explains; authors call g(ΣR) "ad hoc" |
| ORGaNICs with eye-specific units ([H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/); [Rawat 2024](https://arxiv.org/html/2409.18946v3)) | **Explains,** but published τ (1–2 ms) give only a few ms at high c, so it needs a loop delay or slower τ_y (derived) | **Explains:** τ_eff ∝ 1/√(σ²+Σwz²); transient over-normalization could supply the missing late collapse (toy) | **Explains** why scaling the pool τ fails; steady state nests the team's static model (derived) | **Explains:** offset ≈60 ms regardless of stimulus; toy 38 → 17 ms rise, 30–32 ms decay | **Explains:** shared modulator | **Explains:** steady state is normalization |
| Generic SSN, eye-summed input ([Rubin 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/); [Hennequin 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5971207/)) | **Explains:** low gain near spontaneous rates | **Contradicts in the k metric:** an eye-summed drive gives k ≡ 2 (derived) | Explains the absence of a fixed delay | **Partly:** speeds with drive, but decays at cellular τ (~10 ms) unless the front end is slow | Silent (single layer) | **Constrains:** needs Ω_E > 0 or a peak beyond 90% |
| L2/3-only SSN with L4 as input ([Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025); [Mossing 2021](https://doi.org/10.1101/2021.03.31.437953)); laminar neural mass ([Eskikand 2022](https://doi.org/10.1101/2022.12.09.519705)) | **Partly:** linear L4 fits the earliest sink | Within L2/3 only | Not addressed | **Contradicts:** offset ≈10 ms regardless of contrast | **Contradicts:** predicts L4 near-linear while L2/3 normalizes, or ISN in L4/L5 but not L2/3 | Uses a saturating Naka–Rushton input |
| Ocular-dominance-column binocular SSN ([Selby & Tripp 2017](https://www.sciencedirect.com/science/article/abs/pii/S0893608017300552)) | Silent (static) | Could give k < 2 through cross-column structure (inference) | Not fitted | Silent | Silent | Could give interocular normalization without supersaturation (untested inference) |

Three filters do most of the work. F4 removes every model that divides a fast response: static normalization, DN and its contrast-scaled and two-stage variants, Ding–Sperling, and instantaneous-target rate models. F2, read in the equivalent-contrast metric, removes any circuit driven only by the summed eye signal. F5 removes architectures that confine normalization to L2/3. Three classes survive on paper: ORGaNICs with eye-specific units, recurrent conductance normalization, and an ocular-dominance-column SSN. **None has been fitted to time-resolved V1 data of any kind**, binocular or monocular, as Heeger & Zemlianova state for ORGaNICs.

## 3. Interocular normalization at 20–35 ms falls in the local recurrent range

### Monocular suppression is fast and feedforward; interocular suppression is slow and cortical

The circuit literature separates two routes to division.

**Monocular cross-orientation suppression is largely feedforward.** In cat simple cells, cross-oriented stimuli suppressed rather than enhanced synaptic inhibition. LGN saturation and rectification, fed into a linear model, predicted the suppression, and the spike threshold roughly doubled it ([Priebe & Ferster 2006 · cat-A](https://doi.org/10.1038/nn1660)). Masks that drift too fast to drive cortex still suppress: the suppression cut off at a median 19.6 Hz, against 7.7 Hz for V1 responses, and was barely touched by adaptation ([Freeman et al. 2002 · cat-A](https://doi.org/10.1016/S0896-6273(02)00819-X)). LGN saturation, rectification and threshold "can account for most of the COS" ([Li et al. 2006](https://doi.org/10.1152/jn.00425.2006)). LGN neurons already have an untuned, high-temporal-frequency suppressive field ([Bonin, Mante & Carandini 2005](https://doi.org/10.1523/JNEUROSCI.3562-05.2005)). In macaque, cross-orientation suppression began ~6–7 ms *before* the excitatory onset latency (42.5 vs 50.0 ms), and the authors "confidently rule out feedback" ([Smith, Bair & Movshon 2006 · mac-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161/)).

**Interocular suppression has the opposite signature in every cat test.** It peaked at 4 Hz (28.6%), was not significant at 16 Hz and was absent at 32 Hz. Adaptation cut it from **33.4% to 5.6%**, and the GABA-A antagonist bicuculline cut it from **46.1% to 18.1%** ([Sengpiel & Vorobyov 2005 · cat-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC6725266)). Dichoptic cross-orientation suppression is strongest at low temporal frequency and is weakened by adaptation ([Li et al. 2005 · cat-A](https://doi.org/10.1152/jn.00203.2005)). It appears "only if the neuron is already responding", which suggests inhibition between neighbouring ocular-dominance columns ([Sengpiel, Blakemore & Harrad 1995 · cat-A](https://doi.org/10.1016/0042-6989(94)00125-6)). Human psychophysics matches: the interocular suppression weight *grows* with stimulus duration while the ipsiocular weight shrinks ([Baker, Meese & Summers 2007 · human](https://doi.org/10.1016/j.neuroscience.2007.01.030)).

**Primate LGN interocular suppression is weak and conditional.** The median binocular effect was −2.65% (n.s.), reaching −9.95% only with dominant-eye contrast 0.12 and non-dominant-eye contrast 0.9 ([Dougherty et al. 2021 · mac-W LGN](https://pmc.ncbi.nlm.nih.gov/articles/PMC8035044/)). So an inherited component at onset could exist only in that corner of the matrix, and the team's finding that the earliest sink combines near-linearly should be checked there specifically.

### The team's delay matches local recurrent and lateral estimates

| Suppressive signal | Delay after excitatory onset | Prep. | Source |
|---|---|---|---|
| Monocular cross-orientation suppression | ≈0 (latency 42.5 vs 50.0 ms) | mac-A | [Smith 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161/) |
| Surround suppression | 9 ms mean (SD 15); strong suppression ~30 ms earlier than weak | mac-A | [Bair 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744/) |
| Disparity-tuned suppression | 7 ms median | mac-W | [Tanabe & Cumming 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4044366) |
| Binocular modulation of monocular neurons | ~12 ms between pooled medians (50 vs 62 ms) | mac-W | [Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/) |
| Fast suppression, laminar MUA | 11.4 (4Cα) to 21.3 ms (4Cβ) | mac-W | [Wang 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/) |
| **Team: interocular normalization** | **~20–35 ms after each channel's onset** | mac-W | team |
| Population contrast-gain build-up | τ ≈ 31 ms | mac-W | [Reynaud 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242/) |
| Slow suppression, laminar MUA | 44.7–61.0 ms | mac-W | [Wang 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/) |
| Tuned extra-classical suppression (peak) | 76 ms peak, vs 55 ms untuned | mac-A | [Henry 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/) |
| Feedback-recipient figure-ground signal | ~30 ms after horizontal boundary signals | mac-W | [Self 2013](https://doi.org/10.1016/j.cub.2013.09.013) |

Awake V1 responses begin at about 40–50 ms (L4 42.1 ms in [Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/); 40 ms in L4 and L6 in [Self 2013](https://doi.org/10.1016/j.cub.2013.09.013)). The team's normalization therefore falls at roughly 60–85 ms after stimulus onset (derived). That is later than the feedforward, LGN-inherited group. It overlaps fast suppression in L2/3 (20.5 ms) and 4Cβ (21.3 ms) and the ~31 ms contrast-gain build-up. It comes before slow, feedback-attributed suppression, which is about 3× stronger in output layers ([Wang 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/)).

Dougherty et al. reach the same conclusion about mechanism. Binocular suppression "occurred several milliseconds following the onset of the visual response", which suggests "cortical inhibition". Significant modulation needed dominant-eye contrast ≥ 0.8, which matches F2's stronger normalization at high contrast. The strength–latency relation for surround suppression also fits F2: strong suppression arrived "on average 30 msec earlier than weak suppression" ([Bair, Cavanaugh & Movshon 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744/)). In alert macaque, optogenetic drive of excitatory cells facilitated at low contrast and suppressed at high contrast. Light and visual responses summed sub-additively, with a median additive index of 0.67 ([Nassi et al. 2015 · mac-W](https://doi.org/10.1016/j.neuron.2015.05.040)).

Human MEG adds a timing constraint. A neutral-density filter that "delayed monocular processing by approximately 38 ms in V1" reduced that eye's masking ability ([Chadnova et al. 2018 · human](https://eprints.whiterose.ac.uk/123813)). A pool integrating over much more than ~30 ms would have absorbed a 38 ms lag, so this argues for a fast or coincidence-sensitive interocular pool (inference, untested).

The team's w > a (F3) has a circuit counterpart. Monocularly dominated neurons are *more* suppressed binocularly: in two-photon imaging binocular modulation varies with |ODI| as y = −0.36x + 0.14 ([Zhang et al. 2024 · mac-W](https://doi.org/10.7554/eLife.92839)), and in electrophysiology stronger eye preference predicts more suppression ([Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)). Mitchell et al. placed the ocular-dominance weight on the non-dominant eye's excitatory drive, with an Rmax-only fit at adjusted R² 0.92 ([Mitchell et al. 2023 · mac-W, this dataset](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530424)). Algebraically that is the same asymmetry: the pool counts the non-dominant eye relatively more than excitation does (inference). The proposed substrates are interneurons that cross ocular-dominance columns within L4, plus interlaminar feedback ([Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)), or short-range horizontal links between neighbouring columns ([Cox et al. 2019 · mac-W](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477)).

### Conductance data say "shunting" must be read as an effective time constant

Literal untuned shunting is not what cortex shows. Input conductance rose 20–300% with contrast but was "invariably maximal" at the preferred orientation, which contradicts untuned-shunting normalization ([Anderson, Carandini & Ferster 2000 · cat-A](https://doi.org/10.1152/jn.2000.84.2.909)). Optimal flashed bars raised somatic conductance to more than 3× rest at short latency ([Borg-Graham et al. 1998](https://doi.org/10.1038/30735)). Cortical suppression usually works by **withdrawing** recurrent drive. Iso-oriented surrounds cut excitatory conductance by 51% and inhibitory conductance by 38%, after a 30–50 ms transient rise in inhibition ([Ozeki et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725)). Callosally driven normalization at 100% contrast reduced IPSC (−0.91 nS) and EPSC (−0.34 nS) conductances and showed no transient IPSC increase. It was additive over 0–150 ms and divisive over 150–300 ms ([Sato et al. 2016 · mouse](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817833)). In awake mice, inhibition "was much stronger than excitation and had extremely broad spatial selectivity" ([Haider et al. 2013](https://doi.org/10.1038/nature11665)).

If interocular normalization also withdraws excitation and inhibition, membrane time constants would lengthen during suppression rather than shorten (derived). The speed-up with contrast would then be a property of the network, set by effective recurrent gain as in ORGaNICs or the SSN, not of literal membrane shunting. That is one reason to prefer ORGaNICs over a literal conductance model.

ISN-style withdrawal also predicts that net current sinks shrink at the same time and depth as spiking. That is F1's "synaptic currents and spiking together" (inference). CSD cannot separate E from I, so this is weak evidence either way. No study has measured intracellular E and I conductances during interocular suppression in any species.

## 4. Depth-uniform timing is measured from each layer's onset, and that changes what it means

### The team's pattern matches none of the published laminar sequences

Every published laminar sequence of suppression differs from the team's. In anesthetized macaque, full-field suppression inherited from LGN appears in L4C "at the same time as the onset of visual responses" (4 ± 1.15 ms) ([Bijanzadeh et al. 2018 · mac-A](https://doi.org/10.1016/j.neuron.2018.08.020)). Near-surround suppression appears first in upper L4 (8.8 ms) and the upper deep layers (19.6 ms), then L2/3 (25.5 ms), and last in L4C (34.2 ms). Far-surround suppression appears almost together in the upper deep layers, upper L4, L2/3 and the lower deep layers (8–16 ms) but reaches L4C only at 42.5 ms. Near-surround current sinks begin in L2/3 and the upper deep layers 13–16 ms after the earliest L4C sink, and far-surround sinks in L1 and the lower deep layers 22–28 ms after it.

Awake laminar data agree that suppression differs by layer. Fast suppression follows excitation by 11.4 ms in 4Cα, against 14.9–21.3 ms elsewhere, and is strongest in the input layers (index 0.72 in 4Cα). Slow suppression is about three times stronger in L2/3 and L5 ([Wang et al. 2020 · mac-W](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/)). Tuned suppression peaks at 60–75 ms in 4Cα but ~100 ms in 4Cβ and L2/3 ([Henry et al. 2020 · mac-A](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)). Binocular modulation of monocular neurons begins first in L2/3 (47.6 ms), then L4 (50.1 ms) and L5/6 (54.7 ms) ([Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)). The overall median of 62 ms in the same paper implies different subsets, so these numbers are unverified as a set.

Theory does not settle which layer should normalize more. Single-layer SSNs put the computation in L2/3 and treat L4 as linear input. The laminar neural-mass model puts ISN operation in L4 and L5. Mouse data show ISN behaviour superficially and deep. In marmoset, granular PV cells "implement divisive/multiplicative linear gain control" while extragranular PV cells "exert tuned nonlinear suppression" ([Vafaei et al. 2025, preprint](https://www.biorxiv.org/content/10.64898/2025.12.23.696300.full.pdf)). So similar timing across layers may hide different mechanisms in different layers.

### Three hypotheses, separated by stimulus-locked timing

The key subtlety is that the team measures normalization relative to **each channel's own response onset**. Onsets themselves are staggered: superficial neurons lag L4C by ~15 ms ([Maunsell & Gibson 1992](https://doi.org/10.1152/jn.1992.68.4.1332)). That gives three distinguishable hypotheses (derived). Under **H1, local normalization in every layer**, each layer has its own recurrent loop with similar time constants, so onset-locked delays are uniform with depth and stimulus-locked normalization times are staggered like the onsets. Under **H2, a single column-wide pool or ISN state change**, a local perturbation changes E and I in every layer within one recurrent time constant ([Ozeki 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725); [Sanzeni 2020](https://doi.org/10.7554/eLife.54875)). Stimulus-locked normalization time is then uniform and onset-locked delays shrink in late-onset layers, and superficial channels begin "less linear" because the pool is already active when they start, which matches F5's second clause. Under **H3, feedforward inheritance from L4**, onset-locked delays are again uniform, but superficial normalization should be predictable from L4 output lagged by the onset difference, and superficial channels can never be less normalized than their L4 input.

F5's onset-locked uniformity favours H1 or H3. It does not exclude H2 if the onset stagger in this dataset is small compared with the resolution of the delay estimates.

The checks follow directly, and the first is to **report onset latency by depth**: if onsets are nearly simultaneous across the probe, the onset-locked versus stimulus-locked distinction cannot be resolved and "uniform" says little about mechanism. Normalization onset should then be **estimated both ways**, per channel with bootstrap CIs, with uniformity tested against pre-set equivalence bounds (for example ±5 ms) rather than inferred from a failure to reject differences. The depth axis should be **distance from the earliest sink in each penetration**. With ±0.5–1 mm label uncertainty, a named "L4C" bin can mix 4C with upper L4 or deep L2/3 and blur the Bijanzadeh-type distinctions of 10–30 ms, and an early sink is identifiable in only 51% of 810 probes in a multi-lab sample ([Mendoza-Halliday et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917659)).

The **CSD sequence of the normalization** is the next readout: where do net sinks shrink first? Repetition suppression in this lab's data first affected supragranular synaptic currents and then spiking across layers ([Westerberg et al. 2019 · mac-W](https://doi.org/10.1152/jn.00113.2019)), which is the analysis template. An onset in L1 or the deep layers first would implicate feedback, and suppression in the L4C sink at onset would implicate inheritance from LGN. The most publishable laminar question is **whether interocular normalization in the earliest-sink band is as fast as elsewhere**; if it is, the result differs from every surround-suppression sequence. The team should also **state which way superficial channels are "less linear"**. In an SSN, binocular units at weak drive sum supralinearly in rates, whereas inherited normalization gives k < 2, so the direction decides between the two accounts.

Two cleaner estimates are available. In **monocular L4 units** the non-dominant eye drives no excitation, so the latency of its suppression is a direct estimate of the pool delay ([Dougherty 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)). And because the first 7–10 ms of cortical activation are magnocellular ([Maunsell & Gibson 1992](https://doi.org/10.1152/jn.1992.68.4.1332)), the team should **test whether early linearity depends on contrast**; if it does, part of F1 may reflect a magnocellular-first input window rather than a pool delay (hypothesis).

## 5. Fit a delayed binocular ORGaNICs circuit that must predict the offset

### The model and how it nests the current fits

The recommended model is ORGaNICs in the Rawat et al. form. It adds eye-specific excitation, a shared output-driven modulator fed by a pool unit with its own eye weights, a loop delay, and a front end fixed from the earliest current sink. It is written so that the team's static model is its exact steady state and the team's DN is approximately its feedforward, instantaneous-response limit (both derived). The equations are proposed for this report and are not published:

```
Front end, per eye e ∈ {D, N}; fixed per penetration from the earliest current sink
  x_e(t) = f(c_e) · [h_in ∗ s](t − δ_e)       f(c) = c^q / (c^q + c_in^q)       h_in = m-stage low-pass, τ_in
Recorded unit and pool unit (eye weights as in the team's DN)
  z(t)   = x_D(t) + a · x_N(t)                  z_p(t) = x_D(t) + w · x_N(t)
Responses (ORGaNICs, Rawat et al. 2024 form)
  τ_y dy/dt   = b · z(t)   − √A(t) · y(t)       r   = [y]₊²
  τ_y dy_p/dt = b · z_p(t) − √A(t) · y_p(t)     r_p = [y_p]₊²
Shared modulator (recurrent, output-driven, loop delay Δ)
  τ_A dA/dt = −A + σ² + r_p(t − Δ) · A(t − Δ)
Steady state (derived):  r = b² (x_D + a x_N)² / (σ² + b² (x_D + w x_N)²)
Recorded channel:        R(t) = β · r(t) + baseline
```

Before the modulator rises above σ², the response follows the summed excitatory drive, so k = 1 + a, the linear value for a balanced channel (derived). With w > a, binocular stimulation drives the pool harder than monocular stimulation at matched excitation, and the collapse follows. The response integrates with τ_y/√A. That gives a zero-contrast time constant τ_y/σ, which sets the offset, and a shorter one at high pool drive, which sets the rise. Two features of the dynamics need care. The delay Δ supplies the ~20 ms linear window at high contrast that published ORGaNICs parameters cannot produce. Rawat et al.'s stability proof does not cover an explicit delay, so stability must be checked numerically. An oscillation at high contrast would be a gamma prediction, not a bug.

| Variant | Response | Modulator drive | What it nests |
|---|---|---|---|
| M0 DN | Instantaneous (τ_y → 0): r = b²z²/A | Feedforward: A = σ² + low-pass_τ2[b² z_p²] | Team's single-stage DN (approximately, at n = 2) |
| M1 feedforward conductance | τ_y free | Feedforward: A = σ² + b² z_p²(t − Δ), optionally low-passed | [Sit 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/); CH94 conductance form |
| M2 recurrent ORGaNICs | τ_y free | Recurrent: r_p(t − Δ)·A(t − Δ) with τ_A | [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/); CHM97 class |
| M3 = M2 + early interocular stage | as M2 | as M2, plus x_e → x_e/(1 + w1·low-pass[x_e'](t − Δ1)) | Team's dynamic two-stage (w1 ≈ 0.9, Δ1 ≈ 31 ms) |

The free parameters per channel are β, baseline, a, w, σ, τ_y, τ_A and Δ (τ2 for M0; w1 and Δ1 for M3). The output exponent stays at 2. The equivalent-contrast factor is unaffected by any monotone output nonlinearity, so this choice matters for absolute fits but not for k (derived).

The front end (τ_in, m, δ_D, δ_N, c_in, q) should be estimated from the earliest current sink of each penetration and then held fixed. Without that, the slow offset could be absorbed into an input filter that is slow at offset, and the test would lose its power.

The proposed bounds are derived from the sources above. σ should range from 0.02 to 0.6 in front-end units, much wider than the legacy 0.05–0.1. τ_y should be chosen so that the zero-contrast time constant τ_y/σ spans 15–150 ms, a range that covers the CHM τ0 of 25–150 ms, the H&Z 60 ms and the ~39–46 ms implied by the team's offsets. τ_A should run from 1 to 30 ms. Δ should run from 0 to 40 ms, bracketing Dougherty's ~12 ms lag, Wang's 11–21 ms fast suppression and the team's 31 ms two-stage delay. a and w should run from 0 to 1.5, initialized at the team's DN estimates.

Fit hierarchically across channels, with depth relative to the earliest sink as a covariate on Δ, τ_y and w.

ORGaNICs is preferred over a CHM-style recurrent conductance model for four reasons. Its equations are published and stability-analysed. Heeger's code and Blake's fitting scaffolding already exist (MATLAB, to be ported). Its steady state is exactly the team's static model. Its binocular extension is just a choice of eye weights. The CHM conductance form is one its authors call ad hoc.

The SSN should not be fitted next. MUA and CSD cannot identify separate E and I parameters, and its distinguishing predictions (P6 and P7 below) need low-contrast and cell-class data first. Instead, simulate an ocular-dominance-column SSN as a forward model. Give it eye-specific E/I pairs, cross-column inhibition and the same front end, and constrain it to Ω_E > 0 so it cannot supersaturate below 90% ([Ahmadian 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)). Then compute its P1–P6 predictions alongside M0–M3.

### Data and fitting protocol

Use trial-averaged aMUA per channel for all 36 conditions, in 1 ms bins on the corrected photodiode clock. Each window runs from −50 ms through the ~317 ms on period and the following off period up to the next onset. Cover the 25 fitted penetrations, with CSD for the front end.

Fit M0–M3 on the 0–219 ms onset window only, as for the earlier DN fit. Then, without refitting, predict three things: the post-offset decay curve and its contrast slope, the held-out contrast pairs (leaving one cell of the matrix out at a time), and the seven dual-probe J penetrations. Only after these predictions are scored should onset and offset be fitted jointly, for final parameter estimates. Compare variants by held-out error under both loss functions already in use.

### Pre-registerable discriminating tests and expected outcomes

Every entry in this table is a derived prediction unless it carries a citation. The laminar test (stimulus-locked versus onset-locked timing) belongs here too and is specified in section 4.

| Test (criterion) | M0 DN | M1 feedforward conductance | M2 ORGaNICs / recurrent | Ocular-dominance-column SSN | Inherited from front end |
|---|---|---|---|---|---|
| **P1.** Onset-only fit predicts offset half-decay within ±5 ms; predicted contrast slope inside the observed CI | **Fails:** 5–12 ms predicted vs 27–32 ms observed | Slow; decay may lengthen with contrast (toy) | **Passes:** slow, ≈ τ_y/σ, ~flat (H&Z; toy 30–32 ms) | ≈ cellular τ (~10 ms) unless front end slow | MUA decay equals sink decay |
| **P2.** Earliest-sink half-decay ÷ MUA half-decay at the same depth (cortical if < 0.5; inherited if ≈ 1) | ≈ 1, both fast | < 0.5 | < 0.5 | < 0.5 | ≈ 1, both slow; P1 then stops discriminating |
| **P3.** Rise-phase τ across all 36 conditions, plotted against pool drive Σ = (c_D + w c_N)² | No single curve; only collapse onset moves (~τ2σ/c) | Single curve, τ ∝ 1/(1 + kΣ) | Single curve, τ ∝ 1/√(σ² + Σ) | Single, continuous shrinkage with total drive | Rise set by front end |
| **P4.** k(t) or summation ratio dips below its sustained value 5–15 ms after collapse onset at c ≥ .45 | No (monotone) | No | Yes, if modulators lag outputs (toy; depends on parameters) | Onset overshoot growing with contrast ([Echeveste 2020](https://doi.org/10.1038/s41593-020-0671-1)) | No |
| **P5.** LFP/CSD gamma peak frequency vs contrast; binocular vs monocular at matched per-eye contrast | No intrinsic gamma | None | Only above ~50% contrast, frequency rising ([H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)) | Frequency rises with total drive, so binocular > monocular; absent at 5.5% ([Holt 2024](https://doi.org/10.1371/journal.pcbi.1012190); derived) | — |
| **P6.** Early-window additivity at (5.5, 5.5) and (11, 11), against a threshold-only feedforward control | ≤ linear | ≤ linear | ≤ linear | ≥ linear; supralinear for binocular units ([Rubin 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)) | — |
| **P7.** Interocular normalization in narrow- vs broad-spiking units | Pool units sum linearly and are not normalized | — | Both normalized | Both normalized; narrow-spiking grow faster ([Hennequin 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5971207/)); mouse PV/SST split ([Perumal 2026](https://doi.org/10.64898/2026.09.20.753031)) | — |

P1 is the decisive test, and P2 decides how to read it. If the earliest sink decays as slowly as MUA, the offset tells nothing about cortical normalization. The front end would then need its own asymmetric model. If the sink decays fast, a recurrent circuit that predicts a ~30 ms offset from onset data alone would be strong evidence. P3 tests the functional form that distinguishes ORGaNICs (a √ law) from conductance models (a linear law). In the notes' toy simulations, a suppressor in the other eye shortened t10, t50 and t_peak in every model, partly because thresholds were set relative to each response's own maximum (derived). P3 should therefore use absolute divergence times rather than relative latencies.

P4 is the only test that could explain why every fitted model underestimates the late collapse. If no dip appears, a second, slower suppressive component is the next candidate. The awake slow suppression at 45–61 ms is one ([Wang 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7511183/)), and it is consistent with the view that multiple normalization filters operate ([Gulati & Ray 2025](https://www.biorxiv.org/content/10.1101/2025.03.11.642541); [Henry 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)).

Two further analyses decompose effects rather than discriminate between models. Comparing the contrast slope of k between the 20–50 ms and the sustained windows separates the **dynamic part of F2** (pool time τ2σ/c, or τ_eff(c)) from its **static part** (w > a in saturation), and every surviving model must reproduce both. **Carry-over across the three presentations** constrains slow components: with τ2 ≈ 30 ms, DN predicts none across the ~317 ms gap, so any carry-over points to a slow pool component.

### Statistical design lessons from the narrow miss

The dual-probe test failed at 6/7 with p = .063. With n = 7, a one-sided sign test has only one passing outcome, 7/7 (p = 1/128 ≈ .008). The 6/7 result was the best possible near-miss: P(≥ 6 of 7) = 8/128 = .0625 (derived). Future pre-registrations should use continuous, hierarchical effect estimates, for example the predicted-minus-observed offset half-decay with penetration random effects. They should also pre-declare equivalence bounds for "contrast-independent" and "depth-uniform", and a minimum meaningful held-out improvement. The ~1% advantage of the two-stage model shows that small improvements will appear and need a bar set in advance.

## 6. Gaps, conflicts, unverified items, and PDFs for the literature/ folder

### What the literature has not done

No study has fitted ORGaNICs and DN head-to-head on the same time-resolved V1 dataset. Nor has any study fitted a dynamic model to millisecond-resolved neural binocular combination across a dichoptic contrast matrix. The temporal binocular models found are psychophysical, SSVEP or MEG. No SSN models dichoptic onset dynamics, and none has interacting laminar populations. No study reports CSD-resolved onset of interocular suppression by layer. No study reports the pool time constant as an explicit function of contrast in V1 spiking. Interocular E and I conductances, and the class of interneuron that mediates interocular suppression in primates, are unmeasured. No primate study has run a paradoxical-response ISN test, laminar or otherwise. The dataset can fill the first, fourth and fifth of these gaps.

### Conflicts and unverified items

| Item | Detail | Status |
|---|---|---|
| Fast versus delayed normalization | Gain control is complete within 10 ms in anesthetized spikes ([Albrecht 2002](https://pubmed.ncbi.nlm.nih.gov/12163540/)); c50 evolves with τ 31 ms in awake VSD ([Reynaud 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242/)); team sees delayed interocular normalization | Reconcilable if Albrecht's effect is monocular and inherited; untested |
| Dichoptic suppression "hundreds of ms" to peak ([Freeman 2002](https://doi.org/10.1016/S0896-6273(02)00819-X)) | Rests on a personal communication | Conflicts with the team's 20–35 ms and Dougherty's ~12 ms |
| Weak binocular gain control in cat ([Truchard 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209)) | Monocular gain fell ~5×, binocular 1.2×; drifting gratings averaged over 4 s | Conflicts with F3; species, stimulus or analysis |
| ISN inhibitory transient | 30–50 ms rise in inhibition ([Ozeki 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725)) vs none at 80–120 ms ([Sato 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817833)) | Unresolved; different preparations |
| Wang 2020 fast-suppression delay | "10 ms" in the Discussion vs 11–21 ms in Table 1 | Internal inconsistency |
| Dougherty 2019 latencies | Overall median 62 ms vs layer medians 47.6–54.7 ms | Different subsets; verify |
| SSN offset | Cellular τ (~10 ms) predicted; toy SSN matched 29 ms only with a slow front end, regime unchecked | Open |
| ORGaNICs equations | Single-unit H&Z equations verified against Heeger's MATLAB; Heeger & Mackey Eq. 1 input-gain term garbled; SI not read | Partly verified |
| Transient over-normalization | Toy only; Rawat et al. Appendix C eigenvalues not read | Unverified |
| τ conversions from team data | Assume first-order kinetics; front end included | Derived; upper bounds |
| Toy scan | Coarse grid, single unit, illustrative weights; binocular collapse not tested jointly | Derived, illustrative |
| Heeger's "normalization takes time" | No written source found | Unverified recollection |
| Sengpiel & Vorobyov bicuculline numbers | This folder's notes cite the J Neurosci paper; an earlier note had them from a conference abstract | Verify in the PDF |
| Citations partly from memory | Tanabe & Cumming 2014 title; Chadnova 2018 author list; Richard, Chadnova & Baker 2018 title | Verify |
| Details not extracted | Echeveste 2020 unit nonlinearity; Rubin 2015 τ_E and τ_I; Ahmadian 2013 β; Holt 2024 macaque gamma data source | Verify |
| Preprints | Gulati & Ray 2025; Mossing 2021; Eskikand 2022; Perumal 2026; Vafaei 2025; Wu & Miller 2026; Morone 2025; Carlson 2025 | Not peer reviewed |
| Team metric details | The exact definitions of "binocular linearity" and the half-decay reference point were taken as given | Confirm before modelling |
| Mitchell 2024 dissertation | Reportedly proposes normalization that grows through the microcircuit (search summary) | Check overlap before claiming novelty |

### PDFs the PI should place in literature/

| PDF | Why it is needed | Access status from the notes |
|---|---|---|
| Heeger 1992, *Vis Neurosci* 9:181–197 | Original mutual-inhibition equations | Paywalled |
| Heeger 1993, *J Neurophysiol* 70:1885–1898 | Stated contrast-dependent dynamics | APS returned 403; may be free on site |
| Carandini & Heeger 1994, *Science* 264:1333–1336 | Shunting equations | Paywalled |
| Albrecht et al. 2002, *J Neurophysiol* 88:888–913 | Latency shift vs contrast | APS returned 403 |
| Heeger & Zemlianova 2020 SI (pnas.2005417117.sapp.pdf) | Derivations and V1 comparisons | Not retrieved |
| Heeger & Mackey 2019 SI | Input-gain term, derivations | Not retrieved |
| Rawat et al. 2024, Appendix C (arXiv 2409.18946) | Coupled (y, a) eigenvalues | Free; not retrieved |
| Denison et al. 2021, *Nat Hum Behav* | Parameter table | Author PDF on the Heeger lab site; save a copy |
| Chadnova et al. 2018, *NeuroImage* 167:23–30 | Temporal-filter parameters | Repository version exists |
| Selby & Tripp 2017, *Neural Networks* 90:29–41 | Only binocular SSN; whether dynamics were modelled | Paywalled (403) |
| Goris et al. 2024, *Nat Rev Neurosci* | Sub-additivity and variability review | Paywalled |
| Sanzeni et al. 2023, *Neuron* | Mouse–monkey operating regime | Possibly paywalled |
| Palmigiano et al., bioRxiv 2020.11.11.378729 | Optogenetic regime of cortex | Free; fetch failed (HTTP 429) |
| Priebe & Ferster 2006, *Nat Neurosci* 9:552–561 | Cross-orientation suppression mechanism and latencies | Paywalled |
| Sengpiel, Blakemore & Harrad 1995, *Vision Res* 35:179–195 | Ocular-dominance-column interocular inhibition | Paywalled |
| Sengpiel et al. 1998, *Vision Res* 38:2067–2080 | Laminar gain modes | Paywalled |
| Sengpiel & Blakemore 1994, *Nature* 368:847–850 | Interocular control | Paywalled |
| Macknik & Martinez-Conde 2004, *J Cogn Neurosci* 16:1049–1059 | Dissenting awake result | Paywalled |
| Vangeneugden et al. 2019, *Curr Biol* | Feedback and awake surround suppression | Paywalled |
| Borg-Graham et al. 1998, *Nature* 393:369–373 | Conductance time course | Paywalled |
| Menz & Freeman 2003, *Nat Neurosci* 6:59–65 | Coarse-to-fine disparity dynamics | Paywalled |
| Self et al. 2013, *Curr Biol* 23:2121–2129 | Awake laminar timing | Paywalled |
| Longordo et al. 2013, *Nat Neurosci* 16:714–723 | Binocular sublinearity; preparation and timing | Paywalled |
| Bijanzadeh et al. 2018, Supplementary Tables S1–S6 | Absolute laminar latencies | Not read |
| Krishnakumaran et al. 2025, *J Neurosci* | Delayed inhibitory accumulation and gamma | Title only |
| Gulati & Ray, peer-reviewed version | Per-monkey filter values | Check status |
| Anderson et al. 2000, *J Neurophysiol* 84:909–926 | Conductance tuning | Likely free after embargo |
| Carlson et al. 2025, bioRxiv full text | Interlaminar coherence timing | Rate-limited |
| Tsodyks et al. 1997, *J Neurosci* 17:4382–4388 | Origin of the paradoxical effect | Prior-knowledge citation |

## Conclusion

Asking whether "normalization takes time" now has an answer, and it was the less informative question. The more informative one is what normalization changes. If it changes only the response amplitude through a delayed divisor, offsets are fast and rise times are fixed by the front end. If it changes the response's integration time, rise times shrink with pool drive and offsets stay slow and flat. The team's onset and offset kinetics favour the second. That shifts the model contest from "how long is the pool delay?" to "is the integration time set by recurrent gain or by input-driven conductance, and is the slow offset cortical at all?" Existing data can answer both. The model must be fitted to onsets only and required to predict offsets, and the earliest current sink can serve as a proxy for the thalamic input.

The equivalent-contrast metric also carries a theoretical consequence the field has not exploited. Any circuit that sees only the summed eye signal, whether recurrent, delayed, supralinear or inhibition-stabilized, is locked at k = 2 for all time. The fall of k from linear to normalized therefore proves that interocular normalization runs through eye-specific units whose suppression counts the other eye more than their excitation does. That makes ocular-dominance-column structure a required ingredient, not an optional one, for any SSN account. It also makes the team's w > a the parametric form of their main result, rather than a separate finding.

On the laminar side, the uniform onset-locked delays should be re-expressed in stimulus-locked time before they are interpreted. A column-wide pool and local normalization in every layer predict opposite patterns there. Whichever holds, interocular normalization that reaches the earliest-sink band as fast as other layers would set it apart from every published surround-suppression laminar sequence.
