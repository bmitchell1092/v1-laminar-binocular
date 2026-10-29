# Dynamic and recurrent divisive normalization in V1: the Heeger lineage, delayed-pool models, and whether "normalization takes time"

**Notation.** Contrast is c, or c_L and c_R for the two eyes. z is the linear input drive. E is the excitatory drive and S is the suppressive (pool) drive. σ is the semisaturation constant and n is the exponent. τ is a time constant. BSR(t) = R_bin(t) / [R_L(t) + R_R(t)].

**Labels.** Equations come from the cited source unless they are marked **[derived]**. Equations rebuilt from a garbled web extraction are marked **[reconstructed; verify vs PDF]**.

**Context.** The team's findings (1)–(5) are as stated in the brief. I did not research them.

**Paywalled sources the team should supply as PDFs:**
- Heeger 1992 (Vis Neurosci, Cambridge).
- Carandini & Heeger 1994 (Science).
- Heeger 1993 and Albrecht et al. 2002 (J Neurophysiol). The APS site returned 403 to automated fetch and they are not in PMC. They may be free on journals.physiology.org.
- Denison et al. 2021 (Nat Hum Behav). An author PDF is on the Heeger lab site.
- Chadnova et al. 2018 (NeuroImage). A repository version exists.
- SI appendices of Heeger & Mackey 2019 and Heeger & Zemlianova 2020 (pnas.2005417117.sapp.pdf). These were not retrieved, and they contain the derivations and V1 comparisons.

---

## Q1. Heeger-lineage models: the equations, the parameters that set the speed of normalization, and what each predicts for time course vs contrast

### Takeaway
Every Heeger-lineage mechanism puts normalization inside the dynamics, and every one predicts that **the effective time constant shrinks as contrast rises**. Heeger & Zemlianova 2020 quantify this as τ_eff ≈ 60 ms at zero contrast falling to ≈ 6 ms at high contrast, or τ_eff ∝ 1/√(σ² + Σ w c²). This holds for:
- **Shunting / conductance models** (Carandini & Heeger 1994; Carandini, Heeger & Movshon 1997): τ = C/g, and g grows with the pool.
- **ORGaNICs recurrent circuits** (Heeger & Mackey 2019; Heeger & Zemlianova 2020; Rawat et al. 2024).

So normalization arrives faster and is stronger at high contrast. Low-contrast inputs are integrated slowly and amplified about 100×. Responses rise faster than they decay. Slow modulators produce onset transients, and fast modulators produce gamma oscillations at high contrast. Heeger & Zemlianova say this explicitly: "the recurrence takes time (multiples of the time constant)". That is the most direct published form of "normalization takes time" from Heeger.

### Cited Findings

**Heeger 1992: normalization by mutual inhibition**
- Citation: Heeger DJ (1992). Normalization of cell responses in cat striate cortex. *Vis Neurosci* 9(2):181–197. doi:10.1017/S0952523800009640. Paywalled.
- The model modifies the linear/energy model so that striate cells "mutually inhibit one another", normalizing responses for contrast. This is a recurrent framing from the start. — [Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1017/S0952523800009640%22&resultType=core&format=json)
- The paper's dynamic equations were not retrieved; see Gaps.

**Heeger 1993: half-squaring and normalization**
- Citation: Heeger DJ (1993). Modeling simple-cell direction selectivity with normalized, half-squared, linear operators. *J Neurophysiol* 70(5):1885–1898. doi:10.1152/jn.1993.70.5.1885.
- Structure: linear stage → half-squaring → divisive normalization.
- Per the Europe PMC abstract record, the model accounts for nonlinear direction-selectivity data and for contrast-dependent response dynamics. I saw only a paraphrase of the abstract, so get the PDF. — [doi](https://doi.org/10.1152/jn.1993.70.5.1885)

**Carandini & Heeger 1994: summation and division (shunting)**
- Citation: Carandini M, Heeger DJ (1994). Summation and division by neurons in primate visual cortex. *Science* 264(5163):1333–1336. doi:10.1126/science.8191289. Paywalled.
- Simple cells linearly sum LGN input, and that sum is divided by pooled cortical activity. A closed-form version predicted macaque V1 responses across contrast, orientation and spatiotemporal frequency. — [doi](https://doi.org/10.1126/science.8191289)
- Zhou et al. (2019) classify Carandini & Heeger 1994 as **instantaneous** gain control that "misses the transient-decay shape". — [Zhou et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/)

**Carandini, Heeger & Movshon 1997: the explicit dynamic shunting model**
- Citation: Carandini M, Heeger DJ, Movshon JA (1997). Linearity and normalization in simple cells of the macaque primary visual cortex. *J Neurosci* 17(21):8621–8644. doi:10.1523/JNEUROSCI.17-21-08621.1997. Free in PMC.
- Equations:
  - Eq. 1, linear drive: I_d(t) = ∭ W(x,y,T)·c(x,y,t−T) dx dy dT.
  - Eq. 2, membrane: C dV/dt + g(t)V = I_d.
  - Eq. 3: R ∝ max(0,V)^n.
  - Eq. 4: g = g0 / (1 − kΣR) **[reconstructed]**. The conductance depends on the pool's *output* responses ΣR, so this is feedback normalization. The authors call this form "completely ad hoc".
  - Eq. 6, phase advance: arctan(2πfτ0) − arctan(2πfτ1).
  - Here τ0 is the membrane time constant at rest and τ1 is the time constant at unit contrast. "The increase in membrane conductance decreases the time constant."
  — [CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)
- Fitted time constants:
  - Population medians, τ0 / τ1: 25 / 4.9 ms (gratings), 51 / 5.4 ms (plaids), 150 / 7.5 ms (noise masking).
  - Median n ≈ 2.4–2.6.
  - Example cell: τ0 = 37 ms, τ1 = 9 ms. It "responds ∼20 msec sooner" at unit contrast than at 0.12.
  - The authors say the noise-masking ratio would require an 18-fold conductance increase and call this "unlikely to be possible in real cells".
  — [CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)
- Predicted and observed effects of contrast:
  - Raising contrast increases responsivity at high temporal frequencies.
  - Phase advance grows with temporal frequency (e.g., 29.5° at 1.6 Hz to 51.9° at 6.5 Hz in one cell), i.e., integration time shortens.
  - "temporal frequencies <3 Hz yielded strong saturation but little phase advance, whereas temporal frequencies much >6 Hz showed large phase advances but little saturation".
  — [CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)

**Heeger 2017: theory of cortical function**
- Citation: Heeger DJ (2017). Theory of cortical function. *PNAS* 114(8):1773–1782. doi:10.1073/pnas.1619788114.
- Energy function: E = Σ_i α^(i) [ λ^(i) Σ_j (f_j^(i))² + (1−λ^(i)) Σ_j (p_j^(i))² ], with:
  - f = y − z (feedforward), p = y − ŷ (prior), z_j = (v_j)², v_j = Σ_k w_jk y_k^(i−1).
- Dynamics: τ dy/dt = −dE/dy.
- The energy function contains **no explicit normalization**. The paper says quadratic nonlinearities were used "only for the sake of simplicity" and that "the normalization model" is better.
- Relevant timing result: recurrent dynamics delayed responses "by ∼100 ms in layer 2 relative to layer 3" even with τ = 5 ms. Changing the state λ changes response dynamics.
- No prediction about latency vs contrast.
— [Heeger 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5338385/)

**Heeger & Mackey 2019: ORGaNICs**
- Citation: Heeger DJ, Mackey WE (2019). Oscillatory recurrent gated neural integrator circuits (ORGaNICs), a unifying theoretical framework for neural dynamics. *PNAS* 116(45):22783–22794. doi:10.1073/pnas.1911633116.
- Equations **[reconstructed; verify]**:
  - Eq. 1: τ_y dy_j/dt = −y_j + (b_j⁺/(1+b_j⁺)) z_j + (1/(1+a_j⁺)) ŷ_j. The extraction garbled the input-gain term; Heeger & Zemlianova 2020 use b0/(1+b0).
  - Eq. 2: τ_a da/dt = −a + W_ax x + f(y) and τ_b db/dt = −b + W_bx x. The **recurrent modulator a is driven by outputs y** to implement normalization.
  - Eq. 3, steady state: |y_j|² = |z_j|² / (σ² + Σ_k |z_k|²).
- Speed is set by the modulators:
  - When a and b are both large, "the response time courses are dominated by the input drive, so the responses exhibit a short effective time constant".
  - When they are small, the circuit is "dominated by the recurrent drive", giving a long effective time constant.
  - The modulators "determine the state of the circuit by controlling the recurrent gain and effective time constant".
  - Oscillations "emerged because of the inherent delay" in the multiplicative normalization loop.
- No explicit prediction about time course vs contrast.
— [Heeger & Mackey 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6842604/)

**Heeger & Zemlianova 2020: the V1 ORGaNICs normalization circuit (the key paper)**
- Citation: Heeger DJ, Zemlianova KO (2020). A recurrent circuit implements normalization, simulating the dynamics of V1 activity. *PNAS* 117(36):22494–22505. doi:10.1073/pnas.2005417117. Free in PMC.
- Equations **[reconstructed; verify]**:
  - Eq. 1: τ_v dv_j/dt = −v_j + (b0/(1+b0)) z_j + (1/(1+a_j)) ŷ_j.
  - Eq. 2: y_j = [v_j]₊².
  - Eq. 3: z = W_zx x.
  - Eq. 4: ŷ = W_ŷy √y.
  - Eq. 5: τ_a da_j/dt = −a_j + √u_j + a_j √u_j.
  - Eq. 6: τ_u du_j/dt = −u_j + Σ_k w_jk y_k u_k + (σ b0/(1+b0))².
  - Eq. 7: y = [z]₊² / (σ² + W z²).
  - Eq. 9: g = 1/(σ² + c_t² + β c_m²).
  - **Eq. 10: τ = τ_v · ((1+b0)/b0) · √g.**
  — [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)
  - The constant term in Eq. 6 was garbled in extraction. The form (σ b0/(1+b0))² is the one that reproduces Eq. 7 at steady state; see Inferences.
- Parameters:
  - b0 = 0.2, σ = 0.1.
  - τ_v = 1 ms, τ_a = 2 ms, τ_u = 1 ms. τ_u = 10 ms in the onset-transient simulation (Fig. 3E).
  - Δt = 0.1 ms; 250 ms stimuli at 5–100% contrast.
  - The intrinsic constants (1–2 ms) are "∼10 times shorter than experimental measurements".
  - The simulations "have not (yet) been optimized to replicate published data".
  — [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)
- Contrast predictions (published):
  - g = 100 and τ = 60 ms at zero contrast; g ≈ 1 and τ ≈ 6 ms at high contrast. "high-contrast stimuli are integrated over much briefer periods of time (by a factor of ∼10×)".
  - "Responses were amplified by 100× when contrast was low but by only ∼1× when contrast was high."
  - "high-contrast stimuli evoked rapid increases in activity, whereas low-contrast stimuli evoked much slower and more gradual increases."
  - The **offset time constant is ≈60 ms regardless of stimulus**, so rise and decay are asymmetric.
  - Phase advances with contrast.
  - With slow modulators (τ_u = 10 ms), high-contrast responses show strong onset transients. With fast modulators, "the fixed point may be unstable and the responses may exhibit oscillations".
  - The "effective time constant should increase with the square root of effective gain".
  - "This is a challenge for any model that relies heavily on recurrent amplification because the recurrence takes time (multiples of the time constant)."
  — [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)
- Other predictions:
  - Gamma at about 40–50 Hz appears only at high contrast (>50%), and its frequency rises with contrast.
  - Simulated gamma frequency rises with stimulus size, which contradicts data. The authors acknowledge this.
  - An orthogonal 50% mask suppresses a 50% grating by about a factor of 2.
  - Comparisons with dynamics data (Sit et al. 2009; Bair et al. 2003; Albrecht 1995) are by citation only, not by fits.
  — [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)

**Rawat, Heeger & Martiniani 2024: unconditional stability of ORGaNICs**
- Citation: Rawat S, Heeger DJ, Martiniani S (2024). Unconditional stability of a recurrent neural circuit implementing divisive normalization. NeurIPS 2024; arXiv:2409.18946.
- Equation 1:
  - τ_y ⊙ ẏ = −y + b⊙z + (1 − a⁺)⊙(W_r(√y⁺ − √y⁻)).
  - τ_a ⊙ ȧ = −a + b0²⊙σ² + W((y⁺ + y⁻)⊙a⁺²).
  - y⁺ = ⌊y⌋², a⁺ = √⌊a⌋.
- Steady state for W_r = I: y_s = b⊙z / √(b0²σ² + W(b²⊙z²)) and a_s = b0²σ² + W(b²z²).
- Stability:
  - Local stability is unconditional when W_r = I.
  - In 2-D it is proven for 0 < w_r ≤ 1. For w_r > 1 there are 1 or 3 fixed points, at least one of them stable.
  - The ratio τ_y/τ_a sets how much a unit weights normalization against matching its input.
- No fits to V1 data.
— [Rawat et al. 2024](https://arxiv.org/html/2409.18946v3); also in [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11469413/)

**Morone, Rawat, Heeger & Martiniani 2025 (bioRxiv)**
- Citation: Stabilization of recurrent neural networks through divisive normalization. doi:10.1101/2025.05.16.654567.
- ORGaNICs with GOE recurrent weights stay stable even when the spectral radius exceeds 1. "Breakdown of normalization is an early warning signal for the loss of stability."
- I read only the abstract and introduction. — [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12139785/)

**Denison, Carrasco & Heeger 2021: dynamic normalization model of temporal attention**
- Citation: Denison RN, Carrasco M, Heeger DJ (2021). A dynamic normalization model of temporal attention. *Nat Hum Behav* 5:1674–1685. doi:10.1038/s41562-021-01129-1. Paywalled; author PDF on the Heeger lab site.
- Equations:
  - Eq. 1: τ dr_i/dt = −r_i + e_i/(s_i + σ^n).
  - Eq. 2: e_i = a_i (w_i·x)^n.
  - Eq. 3: s_i = Σ_j e_j.
- The pool is computed from **instantaneous** drives, so normalization lives inside the steady-state target and τ only "determines how long the response takes to rise to steady state … and to return to zero".
- Parameter table: n = 1.5, τ_S1 = 52 ms, τ_S2 = 100 ms.
- The paper cites the recurrent circuit (H&Z 2020) as one mechanism with an equivalent steady state.
— [Denison et al. 2021, author PDF](https://www.cns.nyu.edu/heegerlab/content/publications/Denison-NatHumBehav2021.pdf)

**Chapman & Denison 2025: dynamic spatiotemporal normalization (follow-up)**
- Citation: Chapman AF, Denison RN (2025). A dynamic spatiotemporal normalization model captures perceptual and neural effects of spatial and temporal context. *PLoS Biol* 23(12):e3003546.
- Equations:
  - e_i(t) = Σ_T d_i(T)^n · w^E_T, where w^E is exponential with time constant τ_E.
  - s(t) = Σ_i Σ_T e_i(T) · w^S_T, where w^S is exponential with time constant τ_S.
  - τ_R dr_i/dt = −r_i + e_i/(s + σ^n).
- Fixed parameters: τ_R = 52 and n = 1.5 "from prior work".
- "suppression always acts later and across longer time intervals than excitation".
— [Chapman & Denison 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12680323/)

### Inferences

**[derived] Checking the H&Z Eq. 6 constant and Eq. 10.**
- Steady state of Eq. 5: a = √u/(1−√u), so (1+a)/a = 1/√u.
- With ŷ = √y = v, Eq. 1 gives v = B z/√u, where B = b0/(1+b0).
- Eq. 6 with uniform u gives u = K/(1 − Wy).
- Choosing K = (σB)² gives y = z²/(σ² + Wz²), which is Eq. 7. The alternative rendering K = σ²B does not.
- Linearizing Eq. 1 with the modulators held fixed: τ_eff = τ_v(1+a)/a = τ_v/√u = τ_v/(B√(σ² + Wz²)). This is Eq. 10.
- Conclusion: **the ORGaNICs effective time constant scales as 1/√(σ² + Σ_k w_jk z_k²)**. It depends on the total drive in the normalization pool, which in a binocular circuit includes the other eye.

**[derived] Same result in the Rawat et al. form.**
- With a fixed, ẏ = (bz − √a·y)/τ_y, so τ_eff = τ_y/√a_s = τ_y/√(b0²σ² + W b²z²).
- Holding y fixed, the a-equation relaxes with time constant τ_a·a_s/(b0²σ²). This *grows* with contrast, and a is self-exciting (the a⁺²·y² term).
- So in ORGaNICs the output speeds up with contrast while the normalization signal itself can transiently run away and overshoot at high contrast.
- This needs checking against Appendix C of Rawat et al. (the coupled eigenvalues), which I did not retrieve.

**[derived] How each Heeger-lineage model accounts for the team's findings (1)–(5).**
- **Static normalization (closed-form Carandini & Heeger 1994; Moradi & Heeger 2009)**
  - (1) Only weakly. A rising input drive passes through σ early, so a brief quasi-linear phase exists, but no true unnormalized window.
  - (2) Partly: drive reaches σ sooner at high contrast.
  - (3) The team found it fits worse.
  - (4) and (5) Silent.
- **Shunting with membrane dynamics (Carandini, Heeger & Movshon 1997)**
  - (1) Yes. Starting from rest, the initial slope of C dV/dt = I − gV is I/C, independent of g. Inputs from the two eyes therefore add linearly until t ~ C/g.
  - (2) Yes. g grows with pool drive from both eyes, so the collapse is earlier and deeper at high contrast. With CHM's medians, the time constant runs from about 25 ms near zero contrast to about 5 ms at unit contrast.
  - (3) Consistent. The quantity that scales with contrast is the integration time constant of the *response*, not a low-pass on the pool.
  - (4) An eye-specific conductance stage before binocular summation would be the two-stage version.
  - (5) A shared pool driving local conductances gives timing that does not depend on depth.
- **ORGaNICs (H&Z 2020)**
  - (1) Yes. The feedforward term (b0/(1+b0))z dominates the first few τ_v, and the binocular drive z_L + z_R adds before rectification. Normalization appears as recurrent amplification develops, more for weak total drive than strong.
  - (2) Yes. τ_eff ∝ 1/√(σ² + Σwc²) is about 10× faster at high contrast, and self-exciting modulators can over-normalize transiently.
  - (3) Explains why scaling the pool τ with contrast fails. In ORGaNICs the contrast-dependent speed lives in the output's recurrent gain, and the normalization signal is driven by outputs, not by a low-passed copy of the input.
  - (4) Natural implementation: the input modulator b (an input gain), driven by the other eye, or a stacked ORGaNICs stage.
  - (5) A modulator pool shared across depths predicts similar timing at all depths.
- **Dynamic rate model with instantaneous normalization target (Denison 2021; also Said & Heeger 2013 and Li et al. 2017, Q3)**
  - (1) No. Normalization is present from the first moment in proportion, unless the pool S is computed from delayed outputs.
  - (2) Only through the drive crossing σ.
  - (3) Essentially static normalization plus a low-pass filter, so it should do no better than static.
  - (4) Li et al. 2017 already has a monocular interocular-normalization layer followed by a binocular layer.
  - (5) Silent.

### Gaps
- I did not read the full text of Heeger 1992, Heeger 1993 or Carandini & Heeger 1994; all are paywalled. Their exact feedback and shunting equations, and Heeger 1993's stated contrast-dependent dynamics (phase advance, transients), need the PDFs.
- The H&Z 2020 and Heeger & Mackey 2019 equations were extracted from web text with radicals and fractions lost. I verified them by re-deriving Eqs. 7 and 10, but they should still be checked against the PDF and SI.
- H&Z 2020 and ORGaNICs have **not been quantitatively fitted to V1 contrast-dependent time courses**; the authors say so explicitly. I found no published fit of ORGaNICs to single-unit, MUA or laminar dynamics.
- The coupled (y, a) eigenvalues that set the speed of the normalization signal in ORGaNICs are in an appendix of Rawat et al. 2024 that I did not retrieve.
- I found no Heeger-authored paper that applies ORGaNICs dynamics to binocular or interocular normalization.
- I found no written source for Heeger's remark that "normalization takes time". The closest published statement is the H&Z quote that "the recurrence takes time".

---

## Q2. Other dynamic normalization models, and the evidence for delayed (or recurrent) versus instantaneous normalization

### Takeaway
**Feedforward delayed normalization (DN)** divides by a low-passed copy of the drive. It is the main competitor to the recurrent models and fits:
- human ECoG and fMRI (Zhou et al. 2019; Groen et al. 2022);
- single units in Albrecht's dataset (Zhou et al. 2019);
- macaque spikes, LFP, ECoG and EEG (Gulati & Ray 2025).

The fitted pool time constants are 14–30 ms, close to the team's ~30 ms. Several studies independently conclude that normalization is **"rapid but not instantaneous"**:
- Reynaud et al. 2012: τ ≈ 30 ms for contrast gain.
- Tsai et al. 2012: τ ≈ 26 ms.
- Gulati & Ray 2025: 14–23 ms for parallel masks, 5–10 ms for orthogonal masks.

Against delay:
- Albrecht et al. (2002) report gain control "fully expressed within the first 10 ms" once responses are shifted to align latency.
- Smith, Bair & Movshon (2006) find cross-orientation suppression starting about 6 ms *before* the excitatory response onset.
- Sit et al. (2009) explain VSD dynamics with instantaneous conductance.

Surround and orientation-tuned suppression clearly lag, by about 9–25 ms, and weaker suppression lags more. Different suppressive components therefore have different delays.

### Cited Findings

**Zhou, Benson, Kay & Winawer 2019: the delayed normalization (DN) model**
- Citation: Zhou J, Benson NC, Kay K, Winawer J (2019). Predicting neuronal dynamics with a delayed gain control model. *PLoS Comput Biol* 15(11):e1007484. doi:10.1371/journal.pcbi.1007484.
- Equations **[reconstructed; matches Groen 2022]**:
  - h1(t) = t·e^(−t/τ1) − w·t·e^(−t/1.5τ1).
  - R_L = S ∗ h1.
  - R_DN = |R_L|^n / (σ^n + (|R_L| ∗ h2)^n), with h2 = e^(−t/τ2).
  - BOLD = g·ΣR_DN.
- Search bounds: τ1 and τ2 in 0.07–1 s, n in 1–6, σ in 0.01–0.5; shifts were under 10 ms.
- Fitted τ1 was "~40 ms for the single cells, but ~100 ms for broadband ECoG V1".
- "Because delayed gain control is essential for the desired model behavior, we refer to the model as a delayed normalization (DN) model."
- "We fit one set of DN parameters to all 10 response time courses for each cell", so **τ2 is contrast-invariant**.
- At low contrast, responses peak lower and later. The model attributes this to normalization building up later when σ dominates, and fits Albrecht's cells with 96–97% accuracy.
- Scaled and shifted copies of the 70% response explain more than 98% of the variance at other contrasts.
- Rejected alternatives:
  - compressive temporal summation with *instantaneous* normalization, R = R_L^n/(σ^n + R_L^n), which lacks the onset transient;
  - two temporal channels, which have no gain control and miss the slower low-contrast dynamics;
  - Carandini & Heeger 1994.
- A cascaded two-stage DN qualitatively reproduces the pattern across visual areas (S4B Fig).
— [Zhou et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/)

**Groen et al. 2022: human ECoG**
- Citation: Groen IIA, Piantoni G, Montenegro S, Flinker A, Devore S, Devinsky O, Doyle W, Dugan P, Friedman D, Ramsey NF, Petridou N, Winawer J (2022). Temporal dynamics of neural responses in human visual cortex. *J Neurosci* 42(40):7562–7580. doi:10.1523/JNEUROSCI.1812-21.2022.
- Uses the same DN model, with the pool as "a delayed version of the numerator".
- "peak latency sharply decreases with increasing contrast". The transient is less pronounced at low contrast.
- "delayed normalization is critical to predict the slower dynamics associated with changes in stimulus contrast". A larger normalization term captures both contrast-related and adaptation-related reductions.
- Recovery is near-full by an inter-stimulus interval of 533 ms.
- Fitted values are on GitHub, not in the text.
— [Groen et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9546476/)

**Tsai, Wade & Norcia 2012: human SSVEP, dynamics of normalization in masking**
- Citation: Tsai JJ, Wade AR, Norcia AM (2012). Dynamics of normalization underlying masking in human visual cortex. *J Neurosci* 32(8):2783–2789.
- Three models were compared:
  - memory-less pool: 94% variance explained;
  - long-memory pool: 93%;
  - short-memory pool with h(t) = (1/τ)e^(−t/τ): **τ = 26 ms (95% CI 16–36)**, 96%, the best χ².
- "normalization is rapid, but not instantaneous". The abstract gives about 30 ms.
— [Tsai et al. 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3337145/)

**Gulati & Ray 2025 (bioRxiv): macaque V1 across scales**
- Citation: Gulati D, Ray S (2025). Slow-varying normalization explains diverse temporal frequency masking interactions across scales in the macaque visual cortex. bioRxiv doi:10.1101/2025.03.11.642541 (version posted 29 Sep 2025). Not peer reviewed.
- Recordings: spikes, LFP and ECoG from bonnet macaque V1, plus EEG, with counterphasing plaids.
- "adding a simple low-pass filter in the normalization signal sufficiently explained these diverse effects."
- Filter cutoffs were about 7–11 Hz (parallel plaids) and 15–30 Hz (orthogonal). These correspond to normalization integration times of **about 14–23 ms (parallel) and 5–10 ms (orthogonal)**.
- The authors say a single filter is "unlikely" to be sufficient and that multiple filters probably operate.
- They link the slower, parallel-plaid normalization to SST interneurons, whose onset is about 20–25 ms later than pyramidal and PV cells. This is a speculative link.
— [Gulati & Ray 2025](https://www.biorxiv.org/content/10.1101/2025.03.11.642541)

**Krishnakumaran, Pavuluri & Ray 2025**
- Citation: Delayed accumulation of inhibitory input explains gamma frequency variation with changing contrast in an inhibition stabilized network. *J Neurosci* 45. doi:10.1523/JNEUROSCI.1279-24.2024.
- I have only the title, from Gulati & Ray's reference list. It is relevant because it ties delayed inhibitory accumulation to contrast-dependent gamma. — [cited in Gulati & Ray 2025](https://www.biorxiv.org/content/10.1101/2025.03.11.642541)

**Reynaud, Masson & Chavane 2012: macaque VSD**
- Citation: Reynaud A, Masson GS, Chavane F (2012). Dynamics of local input normalization result from balanced short- and long-range intracortical interactions in area V1. *J Neurosci* 32(36):12558–12569. doi:10.1523/JNEUROSCI.1618-12.2012.
- "contrast gain of the population increases over time with a time constant of ∼30 ms".
- c50 fell from 30% to 5% over 40–100 ms (τ = 31 ms), and the CRF slope n fell from about 3 to below 1 (τ = 35 ms).
- The normalization pool swells by about 0.35–0.49 mm over the first 30–90 ms, slowly enough (0.007–0.014 m/s) to suggest polysynaptic recruitment.
- Conclusion: normalization is "rapid but not instantaneous" and "should be considered as a dynamic process".
— [Reynaud et al. 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6621242/)

**Sit, Chen, Geisler & Seidemann 2009: macaque VSD, gain-control model**
- Citation: Sit YF, Chen Y, Geisler WS, Seidemann E (2009). Complex dynamics of V1 population responses explained by a simple gain-control model. *Neuron* 64(6):943–956. doi:10.1016/j.neuron.2009.08.041.
- Equations:
  - C dV/dt = A − gV.
  - g = g0(1 + B).
  - B = b·I ⊗ H.
  - Stage 2 input = V₁^n, with n = 2, a 20 ms input delay, and fitted C1 = 3.19, b1 = 1521.
- Findings:
  - Rising-edge latency decreases and rising slope increases with contrast.
  - "the falling edge is largely independent of contrast" (p > 0.15), and is slower than the rise (t10 = 65 ms).
  - The model "assumes that conductance changes instantaneously with the input". With a 10 ms conductance time constant there was "no qualitative difference".
— [Sit et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/)

**Albrecht, Geisler, Frazor & Crane 2002: single units, argues "normalization is fast"**
- Citation: Albrecht DG, Geisler WS, Frazor RA, Crane AM (2002). Visual cortex neurons of monkeys and cats: temporal dynamics of the contrast response function. *J Neurophysiol* 88(2):888–913. doi:10.1152/jn.2002.88.2.888.
- 200 ms stationary gratings, with the CRF measured every 2 ms.
- Response shapes are "simply scaled and shifted" across contrast, which explains about 95% of the variance.
- Contrast-set gain control and response expansion are "fully expressed within the first 10 ms after the onset of the response". The same c50 and exponent fit from the first to the last 10 ms.
— [PubMed abstract](https://pubmed.ncbi.nlm.nih.gov/12163540/)

**Bair & Movshon 2004: contrast-dependent temporal integration**
- Citation: Bair W, Movshon JA (2004). Adaptive temporal integration of motion in direction-selective neurons in macaque visual cortex. *J Neurosci* 24(33):7305–7323. doi:10.1523/JNEUROSCI.0554-04.2004.
- Lower contrast gives wider, shorter spike-triggered-average peaks.
- Mean cutoff frequency was **10 Hz at low contrast versus 21 Hz at high contrast in V1** (MT: 11 vs 22 Hz).
- The authors argue that signal strength at the soma does not control integration, and they do not frame the effect as normalization.
— [Bair & Movshon 2004](https://pmc.ncbi.nlm.nih.gov/articles/PMC6729763/)

**Smith, Bair & Movshon 2006: cross-orientation suppression is fast, surround is later**
- Citation: Smith MA, Bair W, Movshon JA (2006). Dynamics of suppression in macaque primary visual cortex. *J Neurosci* 26(18):4826–4834.
- Latencies: onset 50.0 ± 15.8 ms, offset 30.1 ± 9.3 ms, cross-orientation suppression 42.5 ± 11.1 ms.
- "suppression started ∼6.4 ms before response onset". It lags the offset response by about 13.5 ms.
- Surround suppression and release came about 12 ms later than cross-orientation suppression.
- "we can confidently rule out feedback as the source of cross-orientation suppression". The authors favour direct feedforward input or fast local inhibition.
- They argue a single shared normalization pool is unlikely because "their response properties should match, but this does not appear to be the case".
— [Smith et al. 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674161/)

**Bair, Cavanaugh & Movshon 2003: surround suppression delay**
- Citation: Bair W, Cavanaugh JR, Movshon JA (2003). Time course and time-distance relationships for surround suppression in macaque V1 neurons. *J Neurosci* 23(20):7690–7701.
- Surround suppression delay averaged 9 ms (SD 15). Mean latencies: CRF onset 52 ms, suppression 61 ms.
- "Strong suppression often occurred as early as CRF onset, whereas weak suppression was delayed by ∼30 msec". Delay correlated with strength (r = −0.48).
- Delay depended little on distance, which points to feedback or fast sources rather than slow horizontal connections.
— [Bair et al. 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740744/)

**Henry, Jazayeri, Shapley & Hawken 2020: laminar timing of extra-classical suppression**
- Citation: Henry CA, Jazayeri M, Shapley RM, Hawken MJ (2020). Distinct spatiotemporal mechanisms underlie extra-classical receptive field modulation in macaque V1 microcircuits. *eLife* 9:e54264. doi:10.7554/eLife.54264.
- Peak times: facilitation 60 ± 3 ms, untuned suppression 55 ± 2 ms, tuned suppression 76 ± 3 ms.
- By layer:
  - 4Cα shows early untuned suppression that is absent in 4Cβ.
  - Tuned suppression peaks at about 60–75 ms in 4Cα and about 100 ms in 4Cβ and L2/3.
- "it is incorrect to equate one cortical area with one canonical normalization process common to all neurons."
— [Henry et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7253173/)

### Inferences

**[derived] Albrecht's result does not refute delayed normalization.**
- "Scaled and shifted" removes the contrast-dependent latency, and Zhou et al. showed that DN itself produces scaled and shifted copies at more than 98% variance explained.
- Albrecht's result is therefore compatible with DN. It is also compatible with shunting or ORGaNICs, where the time constant is a few ms at high contrast, so normalization looks "instantaneous" at the contrasts that dominate CRF fits.

**[derived] Reynaud's falling c50 is what recurrent amplification predicts.**
- A c50 that falls over time (rising contrast gain) means low-contrast responses grow relatively more over time.
- That is the signature of slow recurrent amplification of weak inputs in H&Z 2020, where amplification is about 100× at low contrast.
- A pure static model with a rising input, c50(t) = σ/h(t), would also show it, but only during the rise.

**The fitted pool time constants agree with the team's DN.**
- The team's ~30 ms pool delay sits within the published values: Tsai 26 ms; Reynaud about 30 ms; Gulati & Ray 14–23 ms.
- The orthogonal or untuned suppressive components are much faster: 5–10 ms (Gulati & Ray), and cross-orientation suppression is roughly simultaneous with onset (Smith et al.).
- A single-τ pool is probably a mixture of a fast untuned component and a slower tuned or lateral component. This fits Henry et al.'s separate untuned and tuned suppression peaks.

**[derived] Per-model line for findings (1)–(5).**
- **DN with a fixed low-pass pool (Zhou 2019; Groen 2022; Tsai 2012; Gulati 2025)**
  - (1) Yes. The pool starts at zero, so the early response is unnormalized.
  - (2) Yes qualitatively. The pool reaches σ sooner at high contrast, at t* ≈ τ2·σ/c to first order. But the kinetics are fixed by τ2, and depth is capped at the static asymptote. This matches the team's "too slow/shallow".
  - (4) Zhou's cascaded DN is a two-stage option.
  - (5) A single input-driven pool gives the same timing at all depths.
- **Sit-style shunting with instantaneous conductance**
  - (1) and (2) Yes, through RC charging; see Q1.
  - Also predicts contrast-independent falling edges, as Sit et al. observed.
- **Reynaud's dynamic pool recruitment**
  - (1) Yes.
  - (2) Partly. Their time constants (15–34 ms) were not reported as contrast-dependent.
  - (5) Predicts a spatial (pool-size) dynamic that a single-site laminar probe sees only indirectly.

### Gaps
- I found no study that reports the pool or normalization **time constant as an explicit function of contrast** in V1 spiking data. Zhou fixed τ2 across contrasts, and Reynaud reported τ per location, not per contrast.
- I could not retrieve Albrecht et al. 2002's numbers for latency shift versus contrast (APS returned 403).
- Gulati & Ray 2025 is a preprint; the peer-reviewed version and its fitted per-monkey values should be checked.
- I did not retrieve 2015–2026 single-unit or laminar studies that directly measure the onset of normalization relative to excitation for **same-orientation, overlapping** stimuli (the case closest to dichoptic same-orientation gratings). Candidate leads I did not verify: Freeman et al. 2002 (*Neuron*, "Suppression without inhibition"), Priebe & Ferster 2006, and Krishnakumaran et al. 2025.

---

## Q3. Interocular and binocular dynamics: dynamic models of binocular combination or rivalry, and timing evidence for interocular suppression

### Takeaway
The Heeger-lab binocular models (Said & Heeger 2013; Li et al. 2017) are **dynamic rate models with instantaneous normalization inside a single-τ relaxation**. Li et al. 2017 is already two-stage:
- a monocular layer normalized by both eyes' drives (τ = 10 ms);
- then a binocular summation layer with its own normalization;
- plus opponency (τ = 20 ms), attention (τ = 150 ms) and slow adaptation (τ = 2 s).

These models were built for rivalry on second-to-minute timescales and make no claims about the first 50 ms. The physiology agrees that interocular suppression lags the response:
- macaque V1 binocular modulation of monocular neurons begins about 12 ms after response onset (median 62 vs 50 ms) and needs high dominant-eye contrast (Dougherty et al. 2019);
- disparity-related suppression lags excitation by a median 7 ms (Tanabe & Cumming 2014).

Human MEG shows that an interocular timing difference (a luminance-induced delay of about 38 ms) weakens the delayed eye's ability to mask the other.

### Cited Findings

**Said & Heeger 2013: rivalry and cross-orientation suppression model**
- Citation: Said CP, Heeger DJ (2013). A model of binocular rivalry and cross-orientation suppression. *PLoS Comput Biol* 9(3):e1002991. doi:10.1371/journal.pcbi.1002991.
- Eq. 1 is "a dynamical variant of the normalization equation". The half-squared drive is divided by the pool's half-squared drives plus a semisaturation constant.
- Time constant 50 ms; semisaturation 0.5 (monocular and summation units) and 0.9 (opponency units).
- Opponency neurons, R−L or L−R, "linearly inhibited both monocular neurons in the opposite eye", i.e., subtractive feedback.
- Adaptation time constant 80 s; Euler step 2 ms.
- Simulations used fixed 0.5 contrast, so the paper makes no contrast-dependent timing predictions.
— [Said & Heeger 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3610603/)

**Li, Rankin, Rinzel, Carrasco & Heeger 2017: attention model of binocular rivalry**
- Citation: Li HH, Rankin J, Rinzel J, Carrasco M, Heeger DJ (2017). Attention model of binocular rivalry. *PNAS* 114(30):E6192–E6201. doi:10.1073/pnas.1620475114.
- Monocular layer **[reconstructed]**:
  - τ_s dR_l/dt = −R_l + α E_l / (S_m + H_l^n + σ^n).
  - E_l = [D_l^n − w_o O_r]₊ · [1 + w_a R_a].
  - **S_m = Σ_k E_lk + Σ_k E_rk**, i.e., interocular normalization at the monocular stage.
  - τ_h dH/dt = −H + w_h R.
- Binocular summation layer: E_b = (R_l + R_r)^n with S_b = E_b.
- Attention layer: τ_a = 150 ms. Opponency layer: τ_o = 20 ms.
- Parameters: n = 1 (monocular) or 2 (other layers); σ = 0.5 (σ_a = 0.2); α = 2; τ_s = 10 ms; τ_h = 2000 ms; w_a = 0.6; w_o = 0.65; w_h = 2.
- Normalization has no time constant of its own: each layer is a first-order relaxation towards an algebraic normalized target.
— [Li et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544261/)

**Dougherty, Cox, Westerberg & Maier 2019: binocular modulation timing in macaque V1 (Maier lab)**
- Citation: Dougherty K, Cox MA, Westerberg JA, Maier A (2019). Binocular modulation of monocular V1 neurons. *Curr Biol* 29(3):381–391.e4. doi:10.1016/j.cub.2018.12.004.
- "median response onset: 50 ms, median binocular modulation onset: 62 ms" (p = 2.0 × 10⁻¹⁷).
- "binocular suppression occurred several milliseconds following the onset of the visual response". The authors read this as evidence that suppression relies on intracortical processing.
- In L4, facilitation tended to come before suppression (0.25-quantiles 49 vs 54 ms), but the difference was not significant (p = 0.16).
- Binocular modulation onset by layer: L2/3 47.6 ms, L4 50.1 ms, L5/6 54.7 ms. Visual onsets: L4 42.1, L2/3 45.2, L5/6 46.4 ms.
- "only high contrast stimulation of the dominant eye resulted in significant binocular modulation" (≥0.8 contrast).
— [Dougherty et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)

**Tanabe & Cumming 2014: delayed suppression in disparity coding (monkey V1)**
- Citation as I recall it; title and journal not confirmed from the page: Tanabe S, Cumming BG (2014). Delayed suppression shapes disparity selective responses in monkey V1. *J Neurophysiol*.
- "The suppression was delayed relative to the excitation, with a median difference in latency of 7 ms". Responses to false stereo matches weakened as the response developed.
— [PMC4044366](https://pmc.ncbi.nlm.nih.gov/articles/PMC4044366)

**Chadnova et al. 2018: human MEG interocular normalization and timing**
- Citation: Chadnova E, Reynaud A, Clavagnier S, Baker DH, Baillet S, Hess RF (2018). Interocular interaction of contrast and luminance signals in human primary visual cortex. *NeuroImage* 167:23–30. The author list is from my memory of the paper, not from the retrieved record.
- In V1, SSVER magnitude rose monotonically from 0 to 32% contrast. Dichoptic masking reduced SSVER power, and a normalization model explained the reduction.
- A neutral-density filter (30× lower luminance, same contrast) "delayed monocular processing by approximately 38 ms in V1" and reduced the filtered eye's masking ability.
- The authors added a temporal filtering stage before binocular combination to the normalization model.
— [White Rose repository record](https://eprints.whiterose.ac.uk/123813)

**Richard, Chadnova & Baker 2018: suppression of delayed monocular signals**
- The title I recall, "Binocular vision adaptively suppresses delayed monocular signals" (*NeuroImage*), is from memory; verify it.
- "interocular interactions adaptively suppress signals from the filtered eye and minimize its contribution to the binocular" signal.
— [White Rose repository PDF](https://eprints.whiterose.ac.uk/127595/1/Richardetal2018repositoryversion.pdf)

**Static models for reference (from the sibling notes file)**
- Moradi & Heeger 2009: a 500 ms interocular onset asynchrony did not change the fMRI results.
- Baker & Wade 2017: the static late-summation normalization model reproduced SSVEP intermodulation terms.
— see [binocular_normalization_models.md](../V1%20binocular%20contrast%20normalization/binocular_normalization_models.md); [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/); [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)

### Inferences
- **The binocular models give no constraint on the first 50 ms.** Li et al. 2017 and Said & Heeger 2013 were built and fitted for rivalry alternations over seconds. Their τ values (10–50 ms) are plausible but unconstrained for onset dynamics.
- **[derived] Li 2017 predicts no early linear window.** Because S_m is computed from instantaneous drives, the monocular layer's interocular suppression is proportionally present from onset. This conflicts with the team's finding (1). To get (1), S_m must be computed from delayed or recurrent signals: replace S_m with a low-passed copy (DN-style) or with outputs (ORGaNICs-style).
- **Li et al.'s architecture already matches finding (4).** It has monocular interocular normalization followed by binocular summation with its own normalization.
- **Dougherty et al. match the team's timescale.** Their roughly 12 ms lag of binocular modulation after onset, and suppression arriving later than facilitation, fit the team's early-linear-then-normalized pattern ((1), ~20 ms).
- **Their laminar order is not a simple feedforward cascade.** Modulation appears first in L2/3, then L4, then L5/6. That is weakly in favour of a recurrent or shared source, and matches the team's similar timing across depth (5).
- **Their contrast condition matches (2).** Modulation reached significance only at high dominant-eye contrast, consistent with faster and stronger normalization at high contrast.
- **[derived] The Chadnova and Richard results imply a fast pool.** If an eye whose signal is delayed by about 38 ms loses masking power, the interocular pool cannot be a long low-pass integrator. A pool with τ ≈ 30 ms would still largely overlap a 38 ms-delayed input. This favours a fast or coincidence-sensitive interocular gain control, consistent with an early interocular stage (finding 4). This is an inference and has not been tested quantitatively.

### Gaps
- I did not retrieve Wilson's rivalry models (e.g., Wilson 2003 PNAS, "Computational evidence for a rivalry hierarchy in vision", doi:10.1073/pnas.2333622100). From memory they use Naka–Rushton rate units with τ of about 15–20 ms and slow adaptation; verify.
- I found **no published dynamic model fitted to millisecond-resolved neural binocular-combination data** (single units, MUA or laminar) over a dichoptic contrast matrix. The temporal models of binocular combination I found are psychophysical or SSVEP/MEG.
- The full text of Chadnova et al. 2018 (temporal filter parameters) and of the 2014 conference abstract "Dynamics of dichoptic masking in the primary visual cortex" was not retrieved.
- I did not retrieve "Baker/Wade temporal models" beyond static SSVEP fits; I found no explicit time-domain Baker–Meese–Georgeson binocular model.

---

## Q4. Which predictions would discriminate recurrent normalization (ORGaNICs or shunting) from feedforward delayed normalization (fixed low-pass pool) in time-resolved data?

### Takeaway
Faster collapse at high contrast does **not** by itself discriminate between the models. A fixed-τ2 DN also normalizes sooner at high contrast, because its pool crosses σ sooner (t* ≈ τ2σ/c). The discriminating signatures are listed below. Most are directly testable on a 6×6 dichoptic matrix with 1 ms MUA/CSD.
- **(a) Offset kinetics.** DN responses fall about as fast as the numerator's linear filter. Recurrent and shunting models decay slowly with a contrast-independent time constant (H&Z: about 60 ms; Sit: falling edge independent of contrast and slower than the rise).
- **(b) Rise/decay asymmetry and contrast-dependent integration time of the response itself.** Shunting/recurrent: τ_eff ∝ 1/√(σ² + Σwc²) or C/g. DN: fixed τ1, τ2.
- **(c) Intrinsically generated contrast-dependent latency.** Large in ORGaNICs, small in DN unless the input latency is made contrast-dependent.
- **(d) Transient over-normalization.** In ORGaNICs, a non-monotonic BSR(t) can dip below its sustained value at high contrast.
- **(e) Gamma oscillations.** ORGaNICs predicts gamma whose frequency rises with contrast; DN predicts none.
- **(f) Collapse speed set by total binocular energy.** A single collapse curve should hold when conditions are plotted against total binocular energy √(σ² + c_L² + w c_R²).

### Cited Findings
- ORGaNICs: τ_eff ≈ 60 ms at zero contrast and ≈ 6 ms at high contrast. The offset time constant is ≈60 ms "regardless of stimulus". Gamma appears only at contrasts above 50% and its frequency rises with contrast. Slow modulators give onset transients. — [H&Z 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7486719/)
- Shunting: the time constant shortens with conductance (τ0 ≈ 25 ms to τ1 ≈ 5 ms), producing phase advance. — [CHM 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573724/)
- VSD data: rising latency and slope depend on contrast, and the falling edge does not (p > 0.15). — [Sit et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807412/)
- DN with fixed τ2 fits all contrasts, and lower, later peaks at low contrast come from the pool. Zhou et al. note that DN "can predict convergence of falling edges" but that this is "not a consistent feature across the three cells". — [Zhou et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6892546/)
- Integration bandwidth roughly doubles at high contrast: cutoff 10 Hz to 21 Hz. — [Bair & Movshon 2004](https://pmc.ncbi.nlm.nih.gov/articles/PMC6729763/)
- A pool filter is needed to stop the divisive signal oscillating at the stimulus frequency. The authors say multiple filters probably act at different levels. — [Gulati & Ray 2025](https://www.biorxiv.org/content/10.1101/2025.03.11.642541)

### Inferences

**[derived] Toy simulations.** These are my own; they are not published and are not fits to data.
- Script: `C:\Users\bmitc\AppData\Local\Temp\claude\c--Users-bmitc-Documents-GitHub-BM-ephys\e97c5ad9-a6a5-4557-a778-4a7637ac879a\scratchpad\sim\toy_dyn.py` and `toy_dyn2.py`. The scratchpad is temporary; copy the scripts if needed.
- Common input: c × a step passed through two 5 ms low-pass stages. Binocular drive is z = c_L + c_R; σ = 0.1; output exponent 2.
- The models:
  - **DN30:** τ2 = 30 ms.
  - **ORGaNICs** (Rawat et al. form, W_r = I): τ_y = τ_a = 2 ms.
  - **Shunting:** τ0 = 25 ms, τ = 5 ms at unit input, conductance ∝ z².
  - **Instantaneous-target rate model:** τ = 10 ms.
- Results at equal contrasts c_L = c_R = c:
  - **Early BSR is near the expansive "linear" value (about 2 for n = 2) in DN, ORGaNICs and shunting, but not in the rate model.**
    - At high contrast, the rate model's BSR at +2 ms is 0.66 versus a sustained 0.50.
    - In other words, the instantaneous-target model has no early linear window.
  - **Time to half-collapse shortens with contrast in every model.** Contrasts 0.11, 0.45, 0.9:

    | Model | Time to half-collapse (ms) |
    |---|---|
    | DN30 | 19.9, 7.8, 5.3 |
    | ORGaNICs | 14.5, 5.7, 4.0 |
    | Shunting | 26.7, 9.6, 5.3 |

    Faster collapse at high contrast therefore does not discriminate.
  - **Only ORGaNICs gave a non-monotonic BSR (transient over-normalization).**
    - At c = 0.45, BSR at +10 ms was 0.25 against a sustained 0.52; at c = 0.9 it was 0.35 at +10 ms.
    - Mechanism: the binocular response has the shorter τ_eff, so it overshoots and settles while the monocular responses are still rising.
    - DN and shunting collapsed monotonically to their asymptotes.
    - A DN that is "too shallow at high contrast" (team finding 3) is what this missing over-normalization would produce.
  - **Offset (stimulus off at 200 ms; time to fall to 37%):**

    | Model | Offset time (ms) | Depends on contrast? |
    |---|---|---|
    | DN30 | about 7 | no |
    | ORGaNICs | about 17–20 | roughly no |
    | Shunting | about 23–36 | rises with suppressor contrast |

    So DN falls fast and recurrent or shunting models fall slowly.
  - **Contrast-dependent latency generated by the model itself:**

    | Model | t10 at c = 0.11 (ms) | t10 at c = 0.45 (ms) | Peak at c = 0.11 (ms) | Peak at c = 0.45 (ms) |
    |---|---|---|---|---|
    | ORGaNICs | 13.7 | 8.8 | 52 | 19 |
    | DN30 | 5.2 | 4.2 | 24.5 | 14.1 |

  - **A suppressor in the other eye shortened t10, t50 and t_peak of the suppressed unit in every model.** In DN this is partly an artifact of thresholds set relative to each response's own maximum. Use absolute divergence times, not relative latency, when testing "the suppressor speeds the suppressed".
- Caveats:
  - ORGaNICs' over-normalization depends on τ_a/τ_y. H&Z report oscillations with fast modulators and onset transients with τ_u = 10 ms.
  - The parameters are illustrative only.
  - The early BSR ≈ 2 comes from the n = 2 output. With a team metric where "linear" means BSR = 1, compare time courses, not absolute values.

**[derived] Concrete tests on the 6×6 dichoptic matrix (~317 ms on/off, 3 presentations per trial).**
1. **Offset decay time constant versus contrast cell.**
   - DN predicts a fast, numerator-limited fall, and possibly a post-offset undershoot after high-contrast stimuli because the pool lingers.
   - ORGaNICs and shunting predict a slow, contrast-independent decay (H&Z about 60 ms; Sit's falling edge) and rise faster than they decay at high contrast.
2. **Rise-phase time constant versus √(σ² + c_L² + w c_R²).**
   - Fit an exponential to each cell's rise, or to the deviation from the linear (sum-of-monocular) prediction.
   - ORGaNICs predicts τ ∝ 1/√(σ² + Σwc²). Shunting predicts τ ∝ 1/(1 + kΣwc²).
   - DN predicts a fixed τ1 rise, with only the start of the collapse moving earlier.
   - The 36 conditions should collapse onto one curve against total pool drive if the pool is binocular. If they collapse against each eye's own drive instead, that supports a monocular (two-stage) pool.
3. **Non-monotonic BSR(t).** Look for an over-collapse of BSR at about 5–15 ms after onset at medium and high contrast, below the sustained value. ORGaNICs-type models predict it; DN and pure shunting do not.
4. **Gamma in LFP/CSD.**
   - ORGaNICs predicts gamma only above about 50% contrast, with frequency rising with contrast and power changing with binocular versus monocular drive.
   - DN predicts no oscillation.
   - Caveat: H&Z's size dependence is opposite to data.
5. **Carry-over between the 3 presentations.** DN with τ2 ≈ 30 ms predicts negligible carry-over after 317 ms. Groen et al. found recovery by 533 ms with a longer human τ. Any adaptation-like carry-over constrains slow pool components.
6. **Phase or latency advance versus contrast for the monocular conditions.** CHM 1997: about 20 ms earlier at full contrast versus 12%. ORGaNICs predicts strong intrinsic latency shifts. DN needs an explicitly contrast-dependent input latency to match.

**[derived] Bottom-line ranking against findings (1)–(5).**
- The best qualitative match is a recurrent or conductance-based normalization whose **response time constant**, not pool time constant, scales with total binocular contrast energy (ORGaNICs, or CHM97-style shunting with feedback g(ΣR)).
- Add an early, fast, eye-specific gain-control stage for (4), as in Li 2017's monocular layer or a cascaded DN.
- Treat the "~30 ms delay" as the low-contrast or tuned component, plus a faster untuned component (5–10 ms; Gulati & Ray, Smith et al.).
- A fixed low-pass DN captures (1) and the sign of (2), but its kinetics are fixed by τ2 and its depth is capped at the static asymptote. That explains (3).

### Gaps
- No published study has fitted ORGaNICs and DN head to head on the same time-resolved V1 dataset; this would be a novel contribution.
- I found no published laminar or depth-resolved test of whether normalization time constants differ by layer for same-orientation stimuli. The evidence is mixed:
  - Henry et al. 2020 find layer-specific suppression timing for extra-classical stimuli.
  - Dougherty et al. 2019 find a small laminar order for binocular modulation.
- The toy-simulation over-normalization in ORGaNICs depends on parameters and has not been checked against Rawat et al.'s eigenvalue analysis; treat it as a hypothesis to test by simulation with the team's fitted parameters.
