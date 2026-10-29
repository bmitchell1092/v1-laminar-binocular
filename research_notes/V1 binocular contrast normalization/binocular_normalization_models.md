# Computational models of binocular contrast combination and interocular normalization: equations, competing predictions, and what a 6x6 dichoptic contrast matrix can adjudicate

Notation used throughout: L, R = contrast in left/right eye (or D, N = dominant / non-dominant eye of a neuron); grid = {0, 5.5, 11, 22.5, 45, 90}% in each eye (36 cells). "Dioptic" = same grating, same contrast in both eyes; "dichoptic" = unequal contrasts (or different patterns) in the two eyes. BSR = binocular summation ratio. Unless stated otherwise, equations are copied or reconstructed from the cited source, and anything I derived myself is labelled as such under **Inferences**.

---

## Q1. How has the normalization model been applied to binocular combination? (monocular vs binocular pools, interocular weights, where the exponents go)

### Takeaway
Most binocular normalization models are one family. The neuron (or population) response is a ratio. Each eye's contrast is raised to an exponent in the numerator, and the denominator contains each eye's contrast weighted by an **interocular weight w**. The models differ in three ways: (i) whether the exponent is applied before or after the eyes are summed; (ii) whether w is about 0 (normalization within each eye, followed by linear summation), about 1 (fully binocular normalization, which gives "ocularity invariance"), or somewhere between; and (iii) whether a second gain-control stage follows the binocular sum. Human fMRI and EEG with gratings favour w ≈ 1 with the exponent applied before summation. Cat simple cells (Truchard et al. 2000) favour mostly monocular gain control. Cortical EEG with luminance flicker gives w ≈ 0.02. **The evidence on this question conflicts.**

### Cited Findings
**Moradi & Heeger (2009): explicit monocular vs binocular normalization models, fitted to human V1 BOLD**
- Monocular CRF (Eq. 2): r(c) = a·c^n/(c^n + σ^n). Fit: a = 0.98, σ = 3.6% contrast, n = 1.7 (R² = 0.93). — [Moradi & Heeger 2009, J Vis 9(3):13, doi:10.1167/9.3.13](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Monocular normalization followed by linear binocular summation (Eq. 3/4): r(cL,cR) = r(cL) + r(cR). Each eye is normalized only by its own contrast. This predicts an additivity index of 1. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Nonlinear binocular summation (Eq. 5): r = a·f^n/(f^n + σ^n), where f(cL,cR) is the sum, root-sum-square or max of the two eyes' contrasts, taken before normalization. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Binocular normalization (Eq. 6): r(cL,cR) = a·(cL^n + cR^n)/(cL^n + cR^n + σ^n). The exponent is applied to each eye's contrast before pooling, and each eye's contrast divides the other eye's gain. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Orientation-weighted binocular normalization (Eq. 10): r = a·cL^n/(cL^n + w·cR^n + σ^n) + a·cR^n/(cR^n + w·cL^n + σ^n). Fit: w = 1 for same-orientation gratings and w = 0.52 for orthogonal gratings (R² = 0.83). — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Power-law combination (Eq. 9): r = [r(cL)^(1/k) + r(cR)^(1/k)]^k, with k = 0.65 (range 0.47–0.65). The power law and linear summation both failed. The power law predicted an additivity index of about 0.46–0.57, but the observed index was about 0.01 at 10% contrast. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Conclusion: the data favour binocular normalization ("each eye's activity reduces the gain for the other"). Monocular normalization followed by linear summation did not fit. Orthogonal suppression "does not depend on the eye of origin" and was weaker than same-orientation interocular suppression. Additivity indices were 0.15 for binocular gratings, 0.31 for dichoptic plaids and 0.47 for monocular plaids. A 500 ms onset asynchrony did not change the results. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Stimuli: monocular contrasts 1.5–48%. In the binocular conditions, the left eye was fixed at 2% (Exp. 2a) or 10% (Exp. 2b) while the right eye varied from 0 to 50%. This is a partial dichoptic matrix (two rows), not a full one. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)

**Baker & Wade (2017): six combination rules compared with SSVEP (eye and space)**
- Models (A, B = component contrasts):
  - Linear: A + B.
  - Independent transducers.
  - Early summation (Eq. 1): (A+B)^p/(Z^q + (A+B)^q).
  - Linear numerator (Eq. 2): (A+B)^p/(Z^q + A^q + B^q).
  - Linear denominator (Eq. 3): (A^p + B^p)/(Z^q + (A+B)^q).
  - **Late summation (Eq. 4, the winner):** (A^p + B^p)/(Z^q + A^q + B^q).
  - A masker at a different flicker frequency enters only the denominator (Eq. 5): A^p/(Z^q + A^q + B^q).
  — [Baker & Wade 2017, Cereb Cortex 27:254–264, doi:10.1093/cercor/bhw395](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
- Fits, eye experiment: p = 2.22, q = 2.22, Z = 9.48, Rmax = 0.71. Space experiment: p = 2.43, q = 2.18, Z = 7.46. Late summation had the lowest SSE and AIC; in the eye experiment it explained 92% of the variance and 73% of held-out conditions. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
- The authors call the model's prediction "ocularity invariance": adding a second component at the same frequency barely increases the response. They "speculate that this model implements the optimal algorithm for combining multiple noisy inputs" and link it to a Kalman-filter gain, which strictly holds only in special cases. The eye-experiment fit was worse "mostly due to a shallower-than-predicted contrast response function", which they describe as resembling "a response gain effect... rather than a contrast gain effect". Intermodulation terms (e.g., 2 and 12 Hz) appeared in the data and the model reproduced them. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)

**Meese, Georgeson & Baker (2006) / Baker, Meese & Georgeson (2007): the two-stage model**
- Stage 1 (gain control before summation, includes interocular suppression): binsum = L^m/(S + L + R) + R^m/(S + R + L) (Eq. 1). Stage 2 (after summation): resp = binsum^p/(Z + binsum^q) (Eq. 2). Performance: resp(mask+test) − resp(mask) = σ·k for a criterion d′ (Eq. 3). — [Baker, Meese & Georgeson 2007, Spatial Vision 20:397–413, doi:10.1163/156856807781503622 (Aston post-print)](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- Parameters fitted to Meese et al. (2006) data (Table 1; contrast in % units):
  - Observer DHB: m = 1.19, S = 0.65, p = 8.61, q = 6.75, Z = 0.16, k = 0.14.
  - Second observer: m = 1.28, S = 0.89, p = 6.45, q = 5.15, Z = 0.19, k = 0.09.
  - Direct fits to the 2007 data: m = 1.50 / 1.44, S = 0.20 / 0.60, p = 7.50 / 8.18, q = 6.57 / 7.07, Z = 1.29 / 0.21.
  — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- "p is quite large (about 6 or 7)". If p > q + 1, "the second stage imposes expansion (not compression)". How to interpret stage 2 is "presently unclear". — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- Rewritten as weighted averaging: binsum = w_L·L^(m−1) + w_R·R^(m−1), with w_L = L/(S+L+R) and w_R = R/(S+L+R). S is small (about 1%), so above threshold w_L ≈ L/(L+R) and w_R ≈ R/(L+R). With equal contrasts each eye gets weight 0.5 and the signals are averaged (with compression exponent m − 1 ≈ 0.3). With one eye only, that eye gets weight 1. This is "Hering's law of complementary shares" applied at suprathreshold contrast. Response to monocular stimulation is "(almost) the same as for two eyes, except at low contrasts". In the model, monocular and binocular responses differ markedly only at or below about 1% contrast. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- MGB (2006) abstract: three earlier models were rejected. A development of one rejected model ("twin summation") and the new two-stage model both fit well, and the authors favour the two-stage model. Both have "nearly linear contrast transduction before binocular summation" and "suppressive interocular interactions that contribute to contrast gain control". — [Meese, Georgeson & Baker 2006, J Vis 6(11):1224–1243, doi:10.1167/6.11.7](https://doi.org/10.1167/6.11.7)
- Threshold BSR in this framework: for monocular or binocular stimulation, response ∝ C^m at very low contrast. The threshold BSR is therefore about 2^(1/m) ≈ 1.7 with m ≈ 1.3. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)

**Truchard, Ohzawa & Freeman (2000): monocular vs binocular gain control in cat simple cells**
- Design: one eye's contrast was varied over at least a 10-fold range (e.g., 2.5–50%). The other eye was fixed high (48–51%, n = 23 cells) or low (4–10%, n = 13). Interocular phase was varied over 360° in at least 8 steps. Gratings drifted at 2 Hz for 4 s, and the F1 response was analysed. — [Truchard et al. 2000, J Neurosci 20:3017–3032, doi:10.1523/JNEUROSCI.20-08-03017.2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- Model: each eye has a linear filter (gain G_V or G_F, phase φ); the two are summed linearly; then a binocular gain G_B and an output nonlinearity f(x) = [x − T]_+^n with n = 2. Gains were allowed to vary with contrast. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- Results:
  - Mean monocular gain slope −0.64 (SD 0.15). Mean binocular gain slope −0.05 (SD 0.11), positive in 7 of 23 cells.
  - Relative monocular gain at 5% vs 50% contrast: 5.0×. Relative binocular gain: 1.2×.
  - Residual error: full model 4.9%, monocular-only 6.3%, binocular-only 8.2%. Monocular-only beat binocular-only (p < 0.001).
  - Conclusion: "contrast gain reductions occur primarily at a monocular site". ODI did not correlate with gain slopes.
  - Caveat from the authors: the binocular site could not be separated from the output exponent.
  — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)

**Segala, Bruno, Wade & Baker (2023): the interocular weight w varies with pathway**
- Single-stage simplification of MGB stage 1 ("The second gain control stage was omitted"). The parameters are Z, w, Rmax and n, and the numerator exponent is fixed at 2. The extracted text reads "RespL = L²/(Z + L) + w·R". This is probably a rendering error for the MGB-style form L²/(Z + L + wR); check the original. — [Segala et al. 2023, eLife 12:RP87048, doi:10.7554/eLife.87048](https://elifesciences.org/articles/87048)
- Group medians for w: pupillometry (subcortical) 0.61; EEG experiments 0.02 and 0.02; perceptual matching 0.09. All are below the w = 1 value expected from grating studies. Cortical EEG combination of luminance flicker is near-linear with "no evidence of interocular suppression": "Ocularity invariance is, therefore, not a ubiquitous feature of visual processing". The authors argue the brain repurposes "a generic architecture... by adjusting parameter values". — [Segala et al. 2023](https://elifesciences.org/articles/87048)

**Ocular-dominance-dependent normalization in macaque V1 (recent)**
- Zhang, Zhao, Tang & Yu (2024, eLife), two-photon imaging of superficial V1 in awake macaques (9,390 tuned neurons; Gabors at 0.9 contrast only). Monocular responses were modelled with ODI-derived weights (w_c = (ODI+1)/2, w_i = 1 − w_c) and divisive gain control. Binocular responses add an ODI-dependent interocular suppression factor with exponent b. Results:
  - Monocular neurons are suppressed by binocular stimulation, and binocular neurons are enhanced.
  - Binocular modulation scales with |ODI| (slope −0.41).
  - The binocular response is nearly flat across ODI.
  - No contrast-response or dichoptic-contrast data were collected; the authors list these as future work.
  — [Zhang et al. 2024, eLife reviewed preprint 92839](https://elifesciences.org/reviewed-preprints/92839v2)
- Mitchell et al. (2023): in macaque V1, contrast in a neuron's dominant eye carried more weight, and including ocular dominance in divisive normalization improved fits. — [Mitchell et al. 2023, Curr Biol 33:3884–3895, doi:10.1016/j.cub.2023.08.019](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530424/)
- Georgeson & Sengpiel (2021) re-analysed cat V1 adaptation and interocular-transfer data and compared six two-stage gain-control models. In the best model, "suppression signals arise early, before ocular-dominance attenuation, while suppressive gain control arises late and is affected by ocular dominance". Late gain control partly compensates for ocular imbalance. When the dominant eye was tested, adapting that eye reduced gain more than adapting the other eye. When the non-dominant eye was tested, adapting either eye reduced gain about equally. (Abstract-level information only.) — [Georgeson & Sengpiel 2021, Vision Res 185:29–49, doi:10.1016/j.visres.2021.03.004](https://doi.org/10.1016/j.visres.2021.03.004)

**General normalization equation (background)**
- Canonical form: R_j = γ·D_j^n / (σ^n + Σ_k w_jk·D_k^n). The weights w_jk define the normalization pool. In a binocular version, the pool can include same-eye and other-eye drives. [Citation from prior knowledge; not re-fetched this session: Heeger 1992, Vis Neurosci 9:181–197, doi:10.1017/S0952523800009640; Carandini & Heeger 2012, Nat Rev Neurosci 13:51–62, doi:10.1038/nrn3136](https://doi.org/10.1038/nrn3136)

### Inferences
- **One family covers most models.** Write R(L,R) = F( [a_L·L^p + a_R·R^p] / [σ + b_LL·L^q + b_LR·R^q ...] ) and its per-eye variant (MH Eq. 10 / MGB stage 1). Each published model is then a set of constraints on that template:

  | Model | Exponent placement | Interocular weight | Constraints |
  |---|---|---|---|
  | MH Eq. 3 (monocular normalization, linear sum) | before | w = 0 | — |
  | MH Eq. 6 = BW Eq. 4 (late summation) | before | w = 1 | p = q in MH; p ≈ q in BW |
  | MH Eq. 5 / BW Eq. 1 (early summation) | after linear sum | (implicitly 1) | — |
  | MGB/BMG stage 1 | numerator exponent m ≈ 1.2–1.5 | 1 | denominator linear in contrast; stage 2 = monotone expansive/compressive transducer |
  | Segala et al. | numerator 2, denominator 1 | free (≈ 0.02–0.61) | — |
  | Zhang/Mitchell ocular-dominance variants | — | — | a_L ≠ a_R, and possibly b's ∝ ocular dominance |

  Moradi & Heeger Eq. 6 and Baker & Wade Eq. 4 are algebraically the same model when p = q. Two independent groups, using fMRI and SSVEP, therefore converged on "exponent before summation, binocular pool, w ≈ 1" for same-orientation gratings.
- Where the exponent goes changes the dioptic-vs-monocular relation in a way that does not depend on parameters (derivation, mine):
  - **Exponent before summation, w = 1:** R(c,c) = R_mono(2^(1/n)·c). The binocular CRF is the monocular CRF shifted left by a constant factor 2^(1/n) (≈ 1.5 for n = 1.7; ≈ 1.37 for q = 2.22), a pure contrast-gain shift with the same Rmax.
  - **Linear summation before the nonlinearity (early summation, or Smith et al. 1997 "linear combination"):** R(c,c) = R_mono(2c) for any monotone output nonlinearity.
  - **w = 0:** R(c,c) = 2·R_mono(c), so the binocular response exceeds the monocular ceiling.
  - **Two-stage model or Ding–Sperling-type contrast weighting:** the equivalent-contrast factor falls toward 1 as contrast increases. This is ocularity invariance.
- The grid is nearly octave-spaced (5.5 → 11 → 22.5 → 45 → 90). The cells (c, c) and (2c, 0) are therefore both measured directly. This gives a parameter-free test of "linear binocular summation, then any monotone nonlinearity": the test is R(c,c) = R(2c,0) for equal eye weights. The test is invariant to the neuron's output nonlinearity, so it does not depend on c50 or spike threshold.

### Gaps
- I could not read the original MGB (2006) JOV equations or table (403 from the publisher). The parameter values above come from the BMG (2007) post-print Table 1, which re-fits the MGB (2006) data per observer. Commonly quoted "canonical" MGB values (e.g., m = 1.28, S = 0.99, p = 7.99, q = 6.59, Z = 0.0089) were **not verified**.
- I did not read the exact equations of Mitchell et al. (2023), by instruction (team context). Zhang et al.'s equations are summarized from the eLife page text, not checked against typeset equations.
- Heeger (1992) and Carandini & Heeger (2012) were not re-fetched. Whether Heeger (1992) discussed binocular pools explicitly is unverified.

---

## Q2. Model families: core equations, fitted values, neural status

### Takeaway
There are five families: (1) psychophysical gain-control models (Ding–Sperling, DSKL, MGB two-stage, Georgeson 2016 "max", Baker–Wade late summation); (2) normalization models fitted to pooled neural signals (Moradi–Heeger fMRI, Busse et al. population, Chadnova MEG, Baker EEG); (3) single-neuron binocular combination models (Truchard, Smith/Chino, Anzai, Hallum/Shooner, Zhang, Georgeson–Sengpiel); (4) the disparity energy model, in which contrast enters quadratically and phase matters; and (5) dynamical and recurrent circuits (Said–Heeger opponency, Wilson, SSN). The suppression mechanisms differ: monoptic cross-orientation suppression appears to be fast and subcortical, whereas interocular suppression appears to be cortical, GABAergic, slower and adaptable.

### Cited Findings
**Ding & Sperling (2006) gain-control theory (PNAS 103:1141–1146)**
- Correct DOI: **10.1073/pnas.0509629103** (not ...0508695103). — [Ding & Sperling 2006, PMC1347993](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- Mechanism: "each eye (i) exerts gain control on the other eye in proportion to the strength of its own input and (ii) exerts gain control on the other eye's gain control". Each eye's total contrast energy (TCE) is computed by channel filtering, full-wave rectification with exponent γ, and spatiotemporal summation. For identical images, Γ(I,0) = Γ(0,I) = I and Γ(I,I) = I(2 + 2ε)/(1 + 2ε). — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- Behaviour at the extremes:
  - At very low contrast "the left- and right-eye add linearly". At 50 ms exposure, combination "is well approximated by... linear addition".
  - At high contrast the output depends only on the interocular contrast ratio ("Contrast-Weighted Summation"). Equal combination is the rule (a monocular image matches a dioptic one) only asymptotically.
  - Fitted γ = 1.18 for the observer shown. One parameter explained 95% of the variance across 48 phase × contrast conditions (0.68 c/deg).
  - The gain-control temporal filter has a time constant of about 110 ms.
  — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- Model 2 (simple gain control without gain-control-of-gain-control) was rejected: "the perceived cyclopean sine wave becomes increasingly weaker relative to a monocular sin wave". — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)

**Ding, Klein & Levi (2013) DSKL gain-control + gain-enhancement model**
- The DS model with a fusion mechanism could fit perceived phase or perceived contrast separately, but not both. DSKL adds interocular gain enhancement:
  - Each eye's signal is multiplied by (1 + E*), where E* is the other eye's contrast energy with enhancement threshold g_e.
  - Each eye's TCE exerts gain control on the other eye's signal, gain-control and gain-enhancement layers, with efficiencies 1, α and β.
  - The output exponent is γ. A disparity-energy fusion stage (g_f = 0.038, γ_f = 1) follows.
  — [Ding, Klein & Levi 2013, J Vis 13(2):13, doi:10.1167/13.2.13](https://pmc.ncbi.nlm.nih.gov/articles/PMC4521337/)
- At high contrast, perceived contrast is nearly constant across interocular contrast ratios. The authors argue this is not winner-take-all; it comes from balanced interocular suppression and enhancement, while perceived phase shifts smoothly. Fitted γ ranged from about 1.2 to 3.4 across models and observers. The Table 2 values could not be extracted reliably. — [Ding, Klein & Levi 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4521337/)

**Two-stage model (MGB 2006; BMG 2007); dichoptic masking has several components**
- See the equations under Q1. With all parameters fixed from MGB (2006), the model predicted both discrimination and matching over a **5 × 5 L × R contrast grid (0.3–32%, 1 c/deg)** and beat linear and quadratic (Legge) summation.
  - nMS errors for matching: two-stage 0.78–0.98 dB; linear 2.54–3.52 dB; quadratic 1.09–1.69 dB.
  - When one eye's standard contrast reached 32%, perceived contrast became "essentially independent" of the other eye's contrast.
  — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- Matching is governed by only two parameters (m and S). If responses are equal at the output of the binocular summation stage, they are equal at the final output. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- Varying mask–target phase shows that "interocular suppression and binocular summation are separable". Suppression has tuned and untuned components, and binocular summation is more narrowly tuned for orientation and spatial frequency. — [Baker & Meese 2007, Vision Res 47:3096–3107, doi:10.1016/j.visres.2007.08.013](https://doi.org/10.1016/j.visres.2007.08.013)
- Two suppressive routes act before binocular summation:
  - Ipsiocular: spatially broadband, resists adaptation, weight falls with duration; proposed to arise in retina or LGN.
  - Interocular: spatially tuned, weakens with adaptation, weight rises with duration; proposed to be cortical.
  - The authors predict interocular suppression in cortical cells driven mainly by one eye.
  — [Baker, Meese & Summers 2007, Neuroscience 146:435–448, doi:10.1016/j.neuroscience.2007.01.030](https://doi.org/10.1016/j.neuroscience.2007.01.030)
- The weight of monocular cross-orientation suppression is a power function of speed (TF/SF), whereas dichoptic masking functions overlap across speeds, i.e., they are scale-invariant. — [Meese & Baker 2009, J Vis 9(5):2, doi:10.1167/9.5.2](https://doi.org/10.1167/9.5.2)

**Georgeson, Wallis, Meese & Baker (2016): "Contrast and lustre"**
- Eleven forms of binocular contrast discrimination imply separate binocular and monocular channels, with "only the strongest response... available to perception" (a "max" model). Adding binocular lustre explains out-of-phase pedestals. (Abstract-level; equations not retrieved.) — [Georgeson et al. 2016, Vision Res 129:98–118, doi:10.1016/j.visres.2016.08.001](https://doi.org/10.1016/j.visres.2016.08.001)

**Binocular summation at threshold (meta-analysis)**
- Linear summation L + R predicts BSR = 2. With resp = L^m + R^m, BSR = 2^(1/m), so m = 2 gives √2. The mean BSR across 65 studies is 1.53 unweighted (1.47–1.50 weighted), significantly above √2. Summation declines with speed (TF/SF; R² = 0.70), fitted as m = 1.14 + 0.28·log10(TF/SF). Summation falls with interocular sensitivity imbalance. "There is no single value for binocular summation." Nonlinearities after the summation stage do not change the threshold BSR. — [Baker, Lygo, Meese & Georgeson 2018, Psychol Bull 144:1186–1199, doi:10.1037/bul0000163](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301/)

**Pooled-neural normalization fits**
- Moradi & Heeger 2009 (fMRI): see Q1.
- Chadnova et al. 2018 (MEG source-localized V1; binary noise tagged at 4 and 6 Hz, 0–32% contrast): dichoptic masking reduced responses for both eyes and was well fit by an interocular normalization model. A 1.5 ND filter over one eye delayed that eye's V1 processing by about 38 ms and weakened its masking, which required a temporal filter before binocular combination. — [Chadnova et al. 2018, NeuroImage 167:23–30, doi:10.1016/j.neuroimage.2017.10.035](https://doi.org/10.1016/j.neuroimage.2017.10.035)
- Short-term monocular patching changes dichoptic inhibition in opposite directions for the two eyes (MEG): the unpatched eye becomes more inhibited and the patched eye is released. The authors place this before binocular combination. — [Chadnova et al. 2017, Sci Rep 7:41747](https://pmc.ncbi.nlm.nih.gov/articles/PMC5288724/)
- Busse, Wade & Carandini 2009 (cat V1 arrays; two superimposed gratings, i.e., plaids): population responses are weighted sums. Similar contrasts were summed roughly equally, and very different contrasts gave winner-take-all behaviour. One normalization model fitted spikes and LFP and predicted human V1 source currents. **The abstract reports no dichoptic manipulation**, so this is a monocular-pattern (cross-orientation) normalization result, not an interocular one. — [Busse, Wade & Carandini 2009, Neuron 64:931–942](https://pmc.ncbi.nlm.nih.gov/articles/PMC2807406/)

**Single-neuron binocular combination**
- Smith, Chino, Ni & Cheng (1997), macaque V1: "for all neuron classes, contrast signals from the two eyes combine linearly", meaning a reduction in one eye's contrast is offset by an increase in the other's. In simple cells and phase-specific complex cells, combination happens before nonlinear rectification. (Abstract-level.) — [Smith et al. 1997, J Neurophysiol 78:366–382, doi:10.1152/jn.1997.78.1.366](https://doi.org/10.1152/jn.1997.78.1.366)
- Anzai, Bearse, Freeman & Cai (1995), cat V1 neurometrics: binocular thresholds were generally lower and slopes steeper than monocular. Summation varied across cells "without a single governing rule", and the binocular advantage ratio was slightly higher than typical human psychophysics. — [Anzai et al. 1995, Vis Neurosci 12:77–93, doi:10.1017/s0952523800007331](https://doi.org/10.1017/s0952523800007331)
- Truchard et al. (2000): gain control is mostly monocular in cat simple cells (see Q1).
- Li, Peterson, Thompson, Duong & Freeman (2005), actual title "Cross-orientation suppression: **monoptic and dichoptic mechanisms are different**". Dichoptic suppression is strongest at lower temporal frequencies and is reduced by prolonged adaptation. The authors conclude that monoptic suppression is subcortical and dichoptic suppression is intracortical. — [Li et al. 2005, J Neurophysiol 94:1645–1650, doi:10.1152/jn.00203.2005](https://doi.org/10.1152/jn.00203.2005)
- Freeman, Durand, Kiper & Carandini (2002): (monoptic) cross-orientation suppression survives masks drifting too fast to drive cortex and is unaffected by adapting mask-driven cortical cells. The authors propose a subcortical or thalamocortical origin. — [Freeman et al. 2002, Neuron 35:759–771, doi:10.1016/s0896-6273(02)00819-x](https://doi.org/10.1016/s0896-6273(02)00819-x)
- Sengpiel & Vorobyov (2005), cat: interocular suppression depends on the suppressor's drift rate (fast stimuli cannot trigger it), is almost eliminated by adaptation to the orthogonal stimulus, and is greatly reduced by bicuculline. The authors conclude it is mediated by intracortical (GABAergic) inhibition. — [Sengpiel & Vorobyov 2005, J Neurosci 25:6394–6400](https://doi.org/10.1523/jneurosci.0862-05.2005)
- Sengpiel, Baddeley, Freeman, Harrad & Blakemore (1998), cat area 17: cross-orientation inhibition mainly shifted CRFs (contrast gain), while surround and **interocular suppression mainly lowered response gain**. The modes differed across cortical layers. — [Sengpiel et al. 1998, Vision Res 38:2067–2080, doi:10.1016/s0042-6989(97)00413-6](https://doi.org/10.1016/s0042-6989(97)00413-6). Baker et al. (2021) summarize the layer finding as layer 4 showing response gain and other layers being more consistent with contrast gain (secondary citation). — [Baker, Vilidaite & Wade 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8544832/)

**Amblyopia and interocular imbalance models (neural)**
- Hallum et al. (2017), amblyopic macaque V1/V2 arrays: dichoptic drifting carriers whose contrast was "sinusoidally modulated between 0% and 40%". The modulators drifted at 0.5 and 1.67 Hz, so the two eyes' contrasts varied independently. A difference-of-Gaussians envelope receptive field was fitted for each eye. Findings:
  - Excitation from the fellow eye dominated, especially in severe amblyopes.
  - "suppression from both the FE and AEs was prevalent in all animals"; amblyopia barely affected how common suppression was.
  - The authors propose that interocular inhibitory connections are preserved while excitatory drive from the amblyopic eye is lost.
  — [Hallum et al. 2017, J Neurosci 37:8216–8226](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566869/)
- Shooner et al. (2017): with dichoptic noise masks, amblyopic-eye responses were suppressed by fellow-eye masks, but fellow-eye responses could not be masked from the amblyopic eye. This is asymmetric dichoptic masking. — [Shooner et al. 2017, J Neurosci 37:8734–8741](https://pmc.ncbi.nlm.nih.gov/articles/PMC5588465/)
- Baker et al. (2018) argue that apparent loss of binocular summation in amblyopia reflects monocular sensitivity imbalance, and that "neural binocular mechanisms remain intact". — [Baker et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301/)

**Disparity energy model**
- Ohzawa, DeAngelis & Freeman (1990) proposed a hierarchical (simple → complex) binocular model whose complex cells are "ideally suited as disparity detectors". — [Ohzawa et al. 1990, Science 249:1037–1041, doi:10.1126/science.2396096](https://doi.org/10.1126/science.2396096)
- Read & Cumming (2004), awake macaque V1 with random-dot stereograms: ocular dominance predicted neither the strength nor the class of disparity selectivity, contrary to energy-model expectations. — [Read & Cumming 2004, J Neurophysiol 91:1271–1281](https://pmc.ncbi.nlm.nih.gov/articles/PMC1410815/)
- Tanabe & Cumming (2014), awake macaque V1: strong excitation to one interocular phase difference was often paired with suppression to the opposite phase at lower spatial frequencies. The suppression **lagged excitation by a median of 7 ms** and explains a mismatch between monocular and binocular signals. — [Tanabe & Cumming 2014, J Neurophysiol 111:1759–1769](https://pmc.ncbi.nlm.nih.gov/articles/PMC4044366/)

**Dynamical and rivalry models**
- Said & Heeger (2013): monocular neurons, binocular-summation neurons and **ocular-opponency neurons** (rectified L − R or R − L).
  - Each neuron's rate is its half-squared drive divided by its pool (exponent 2). Pools: all four monocular units; the two summation units; the opponency units by sign.
  - Opponency neurons "linearly inhibited (i.e., via subtraction)" the monocular neurons of the eye supplying their inhibitory input. This is subtractive feedback on top of divisive normalization.
  - Parameters: τ = 50 ms; semisaturation 0.5 for monocular and summation units, 0.9 for opponency units; all weights 1; adaptation τ = 80 s.
  - Rivalry for dichoptic gratings was more than three times as strong as for plaids. Opponency neurons are silent for binocular plaids.
  — [Said & Heeger 2013, PLoS Comput Biol 9:e1002991](https://pmc.ncbi.nlm.nih.gov/articles/PMC3610603/)
- Wilson (2017): a nonlinear dynamical model with two classes of inhibitory interneuron accounts for binocular contrast combination, fusion of similar orientations, and rivalry of dissimilar ones. It predicted reduced perceived tilt with interocular contrast differences, and the author tested this. — [Wilson 2017, Vision Res 140:89–95, doi:10.1016/j.visres.2017.07.016](https://doi.org/10.1016/j.visres.2017.07.016)
- **SSN, binocular extension:** the only binocular SSN found was a 2016 University of Waterloo master's thesis (Selby). It reportedly found SSN parameters giving interocular transfer of suppression in excitatory but not inhibitory units, and suggested an alternative stabilization scheme might be more realistic (search-result summary; not read in full). — [Selby, UWSpace 10012/10264](https://uwspace.uwaterloo.ca/handle/10012/10264)
- SSN analysis reference: Ahmadian, Rubin & Miller 2013, Neural Comput 25:1994–2037 (cited via search summary). — [search result referencing Ahmadian et al. 2013](https://uwspace.uwaterloo.ca/handle/10012/10264)

### Inferences
Equations below are reconstructed or derived by me where marked.

- **DS 2006, effective weights (reconstruction, mine):**
  - With ε_i = (m_i/g)^γ, the output is Out = I_L·(1+ε_L)/(1+ε_L+ε_R) + I_R·(1+ε_R)/(1+ε_L+ε_R). This reproduces the quoted Γ(I,I) = I(2+2ε)/(1+2ε) and Γ(I,0) = I exactly.
  - High-contrast limit: the weights are ∝ m_i^γ/(m_L^γ + m_R^γ), "contrast-weighted summation".
  - For phase-shifted gratings (±θ/2, contrast ratio δ = m_R/m_L) the cyclopean phase follows tan φ = [(1 − δ^(γ+1))/(1 + δ^(γ+1))]·tan(θ/2). I believe this matches DS Eq. 11, but the equation image was unreadable, so verify before quoting.
  - Dioptic BSR (in cyclopean-contrast units) = 2(1+ε)/(1+2ε). This is 2 at threshold, 1.5 at m = g, and tends to 1 at high contrast.
- **MGB and DS make the same suprathreshold prediction, Hering-style weighted averaging:** DS weights ∝ contrast^γ with γ ≈ 1.2; MGB weights ∝ contrast^1 (L/(S+L+R)). Busse et al.'s population "weighted sum, winner-take-all when contrasts differ" for monocular plaids has the same algebra. These are the same computation applied across eyes versus across orientations.
- **Disparity energy model (standard textbook form, from prior knowledge):** complex response C = Σ_{j=1,2} (s_Lj + s_Rj)², where s_Ej = ⟨f_Ej, I_E⟩ for a quadrature pair. For gratings with contrasts c_L, c_R, interocular stimulus phase Δ and receptive-field phase disparity Δφ: C ∝ a_L²c_L² + a_R²c_R² + 2a_La_Rc_Lc_R·cos(Δ − Δφ).
  - Contrast enters quadratically, and the cross-term sign depends on disparity.
  - At preferred phase, the dioptic response is 4× monocular for a balanced cell, with iso-response contours that are straight lines (c_L + c_R = const).
  - At anti-preferred phase, the response is (a_Lc_L − a_Rc_R)², which is zero on the balanced diagonal.
  - The pure model has no normalization. Fleet/Heeger-style or Read–Cumming threshold extensions were not verified this session.
- **SSN (from prior knowledge, flag):**
  - Dynamics: τ dr/dt = −r + k[W·r + h]_+^n with n ≈ 2. Strong recurrent inhibition stabilizes the network.
  - Generic SSN prediction (Rubin et al. 2015): weak inputs sum supralinearly and strong inputs sum sublinearly, with a transition from facilitation to suppression as input strength rises.
  - Applied to two eyes with excitatory cross-eye input (my extrapolation, not a published binocular SSN): dichoptic facilitation at low contrast, binocular response > sum of monocular responses (additivity index > 1) at low contrast, and normalization-like suppression at high contrast.
- Psychophysics (Baker–Meese–Summers 2007; Meese–Baker 2009) and cat physiology (Li 2005; Freeman 2002; Sengpiel–Vorobyov 2005) agree that monoptic suppression is fast, broadband, adaptation-resistant and probably precortical, while interocular suppression is cortical, GABAergic, adaptable and slower or low-TF. That predicts interocular normalization should be weaker or absent at the earliest geniculate-recipient stage and should strengthen after the first intracortical synapse (see Q5).

### Gaps
- DSKL exact equations and Table 2 parameters could not be extracted (equations are images in PMC).
- Georgeson et al. (2016) "max" model equations were not retrieved (paywalled).
- Not researched or verified this session:
  - Kingdom's binocular sum/difference-channel work (e.g., May & Zhaoping 2016; Li & Atick 1994).
  - Read, Parker & Cumming (2002) monocular-threshold energy-model variant.
  - Legge & Rubin (1981) Fechner's paradox for contrast (cited inside BMG 2007 only).
  - Anderson & Movshon (1989).
  - Cogan (1987).
  - Rubin, Van Hooser & Miller (2015) itself.
  - ORGaNICs (Heeger & Zemlianova 2020).
- No peer-reviewed binocular or interocular SSN was found. The Selby thesis claim comes only from a search summary.
- Hallum (2017) model equations are images; the text calls the model linear (difference-of-Gaussians envelope).

---

## Q3. What does each model predict for (a) monocular vs dioptic CRF and BSR(contrast), (b) the dichoptic surface, (c) ocularity invariance, (d) ocular-dominance dependence, (e) temporal dynamics? Where do they conflict?

### Takeaway
The cleanest conflicts on a same-orientation, zero-disparity 6 × 6 grid are:
- **Equivalent-contrast factor k** (how much monocular contrast matches a dioptic stimulus): constant 2 for linear early summation, constant 2^(1/n) for binocular normalization with the exponent before summation, and declining toward 1 for two-stage or Ding–Sperling contrast weighting. Monocular normalization with linear summation predicts that the binocular response exceeds the monocular ceiling, so k is undefined.
- **Monotonicity of the dichoptic surface**: two-stage and DS predict a Fechner-like dip (a weak second-eye contrast lowers the response to a strong first-eye contrast by about 5–15%); normalization models with matched numerator and denominator exponents predict no dip.
- **Iso-response contour geometry**: lines for linear or subtractive models, superellipses for exponent-before-summation, rays for ratio weighting.
- **Timing**: summation in the transient window, suppression in the sustained window. Several macaque data sets already show early facilitation followed by later suppression.

### Cited Findings
- **(a)/(c)** Two-stage model: threshold BSR ≈ 2^(1/m) ≈ 1.7. Above about 1% contrast, monocular and binocular responses are nearly the same (ocularity invariance) and differ "markedly" only at or below about 1% contrast. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- **(a)/(c)** Late summation: adding a same-frequency component to the other eye "barely increases" the response ("ocularity invariance"), and the data showed this. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
- **(a)/(c)** Human V1 BOLD: additivity index ≈ 0.01 at 10% contrast. Two identical gratings gave a response barely above one grating. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- **(c) counter-evidence:** cortical SSVEP to luminance flicker shows substantial binocular facilitation and near-linear combination (w ≈ 0.02). Pupil responses: binocular/monocular ratio near √2 at low-to-intermediate contrast and near 1 at high contrast. — [Segala et al. 2023](https://elifesciences.org/articles/87048)
- **(a)** Ding–Sperling: linear at very low contrast or brief (50 ms) presentations; contrast-weighted at high contrast; monocular = dioptic only asymptotically. — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- **(b)** Perceived contrast with unbalanced interocular contrasts: once one eye reached 32%, perception was "essentially independent" of the other eye's contrast. Linear and quadratic summation underestimated the binocular contrasts needed for a match when one eye's contrast was high. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- **(b)** DSKL: perceived contrast nearly constant across interocular contrast ratios at high contrast. — [Ding, Klein & Levi 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4521337/)
- **(b)** Dichoptic pedestal masking: a dichoptic mask shifts the dipper to higher pedestals while preserving facilitation. With a 32% mask, facilitation reached about 18 dB (8×). — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- **(b)** Macaque V1 (neurons): contrast signals combine linearly, so a decrease in one eye is offset by an increase in the other. — [Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.1.366)
- **(d)** Monocular neurons are suppressed and binocular neurons enhanced by binocular stimulation at 0.9 contrast; the binocular response is nearly flat across ODI (macaque superficial layers). — [Zhang et al. 2024](https://elifesciences.org/reviewed-preprints/92839v2)
- **(d)** Dominant-eye contrast weighted more; ocular dominance in normalization improves fits. — [Mitchell et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530424/)
- **(d)** Late gain control is affected by ocular dominance and partly compensates for imbalance. Adapting the dominant eye reduces dominant-eye gain more than adapting the non-dominant eye does, whereas non-dominant-eye gain is reduced equally by adapting either eye. — [Georgeson & Sengpiel 2021](https://doi.org/10.1016/j.visres.2021.03.004)
- **(d)** No correlation between ODI and monocular or binocular gain slopes in cat simple cells. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- **(d)** Amblyopic V1: excitation dominated by the fellow eye while suppression comes from both eyes. Suppression is therefore not ocular-dominance-weighted the way excitation is. — [Hallum et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566869/)
- **(d)** Asymmetric masking in amblyopic macaques: the fellow eye masks the amblyopic eye, but not the reverse. — [Shooner et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5588465/)
- **(e)** Macaque V1: monocular neurons, including in the input layers, show early facilitation and later suppression under binocular stimulation, attributed largely to cortical inhibition. — [Dougherty et al. 2019, Curr Biol 29:381–391](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)
- **(e)** Macaque V1: an early transient enhancement, followed by a sustained phase that is suppressed for incongruent stimuli and either suppressed or enhanced for congruent ones, depending on ocularity and laminar location. — [Cox et al. 2019, J Vis 19(12):13](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477/)
- **(e)** Disparity-related suppression in macaque V1 lags excitation by a median of 7 ms. — [Tanabe & Cumming 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4044366/)
- **(e)** Interocular suppression weight rises with stimulus duration, while ipsiocular weight falls. — [Baker, Meese & Summers 2007](https://doi.org/10.1016/j.neuroscience.2007.01.030)
- **(e)** DS gain-control filter τ ≈ 110 ms; linear summation at 50 ms. — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- **(e)** An interocular luminance difference delays one eye by about 38 ms in V1 and weakens that eye's masking. — [Chadnova et al. 2018](https://doi.org/10.1016/j.neuroimage.2017.10.035)
- **(e)** Interocular suppression fails for fast-drifting suppressors and is low-TF tuned. — [Sengpiel & Vorobyov 2005](https://doi.org/10.1523/jneurosci.0862-05.2005); [Li et al. 2005](https://doi.org/10.1152/jn.00203.2005)

### Inferences
**Computed predictions on the team's grid.** I computed these in scratch Python using the published parameter values above. Contrast is in % and both eyes are equally weighted. k is the equivalent-contrast factor, defined by R_mono(k·c) = R(c,c). It is invariant to any monotone output nonlinearity applied after combination, so it is comparable across cells with different c50s.

| Model (params) | R(c,c)/R(c,0) at c = 5.5 / 11 / 22.5 / 45 / 90 | k(c) | Additivity index at 5.5% → 90% | R(90,y)/R(90,0), y = 5.5…45 | Dichoptic dip? |
|---|---|---|---|---|---|
| Early (linear) summation + Naka-Rushton (BW Eq. 1, p = q = 2.22, Z = 9.48) | 2.53 / 1.49 / 1.11 / 1.02 / 1.01 | 2.00 constant | 1.53 → 0.01 (super-additive at low c) | 1.00–1.01 | No |
| Late summation = MH binocular normalization (BW Eq. 4, same params) | 1.63 / 1.26 / 1.07 / 1.02 / 1.00 | 1.37 constant (= 2^(1/q)) | 0.63 → 0.00 | 1.00 | No |
| MH Eq. 6 (n = 1.7, σ = 3.6%) | 1.20 / 1.07 / 1.02 / 1.01 / 1.00 | 1.50 constant (= 2^(1/n)) | 0.20 → 0.00 | 1.00 | No |
| Monocular normalization + linear sum (MH Eq. 3) | 2.0 at all c | undefined (binocular > monocular max) | 1.0 at all c | 1.68–1.99 (facilitation) | No; strong facilitation |
| MH Eq. 10, w = 0.5 | 1.50 / 1.39 / 1.35 / 1.34 / 1.34 (→ 2/(1+w)) | undefined at high c | 0.50 → 0.34 | 1.01–1.25 | No |
| Two-stage (BMG 2007 obs. 2: m 1.28, S 0.89, p 6.45, q 5.15, Z 0.19) | 1.11 / 1.05 / 1.03 / 1.01 / 1.01 | 1.19 → 1.02 | 0.11 → 0.01 | 0.96 / 0.94 / 0.92 / 0.93 | **Yes**, 4–8% (and R(45,y)/R(45,0) = 0.94 / 0.92 / 0.93 for y = 5.5–22.5) |
| Ding–Sperling (γ = 1.18, g = 5%) | 1.31 / 1.16 / 1.08 / 1.04 / 1.02 | 1.31 → 1.02 | 0.31 → 0.02 | 0.97 / 0.94 / 0.89 / 0.86 | **Yes**, 3–14% |
| Ding–Sperling (γ = 1.18, g = 20%) | 1.70 / 1.50 / 1.30 / 1.16 / 1.08 | 1.70 → 1.08 | 0.70 → 0.08 | 0.98 / 0.96 / 0.92 / 0.92 | **Yes**, 2–8% |
| Energy model, matched phase, no normalization | 4.0 at all c | 2.0 | 3.0 | 1.13 → 2.25 | No; strong facilitation |

(The DS output is a "cyclopean contrast". k and the sign of the dip are unchanged by any monotone CRF applied afterwards, so these columns still apply to neurons. The g values are illustrative because the DSKL thresholds were not reliably extractable.)

**What each model predicts, by question**
- **(a) BSR vs contrast.** Every normalizing model predicts the response-ratio BSR falls toward 1 at high contrast. The ratio itself therefore discriminates poorly at 45–90%; most of the information is in the 5.5–22.5% cells. The equivalent-contrast factor k discriminates better: constant (early or late summation) versus declining toward 1 (two-stage, DS). Note that neural c50s are typically far above MH's fMRI σ of 3.6%, so neurons should show larger response-ratio BSRs across more of the grid than the BOLD fits imply.
- **(b) Dichoptic surface.**
  - Normalization models with matched exponents (MH Eq. 6, BW Eq. 4) are monotone: adding contrast to the other eye never lowers the response.
  - The two-stage and DS models have a **non-monotone** surface. The derivative of MGB stage 1 with respect to the weak eye at zero contrast is −R^m/(S+R)², which is negative (my derivation). This is the neural analogue of Fechner's paradox, and it is a sign test that does not depend on the output nonlinearity.
  - The w = 0 and energy models predict facilitation along every row.
- **(c) Ocularity invariance.** w = 1 models give it at saturation only through Rmax equality. Two-stage and DS give it at all suprathreshold contrasts (k ≈ 1). w < 1 models and the energy model never give it.
- **(d) Ocular dominance.** The models make conflicting predictions:
  - (i) **Excitatory (numerator) weights ∝ ocular dominance, with an OD-independent binocular pool** (Hallum-like; Zhang-like). The non-dominant eye of a monocular neuron is purely suppressive: off-diagonal cells fall below the dominant-eye monocular row. The two eyes' monocular CRFs differ by response gain (Rmax), not c50.
  - (ii) **Numerator and pool weights both ∝ ocular dominance** (fully OD-weighted normalization). A monocular neuron is unaffected by the other eye, and the non-dominant-eye CRF is a contrast-gain (c50) shift of the dominant-eye CRF.
  - (iii) **Georgeson–Sengpiel two-stage**: suppression signals arise before ocular-dominance attenuation (so they are OD-independent, as in i), while late gain control is OD-weighted and compensatory.
  - The monocular CRF pair (row L = 0 and column R = 0) already separates (i) from (ii): equal c50 with different Rmax means (i); shifted c50 with equal Rmax means (ii).
- **(e) Dynamics.** The static models (MGB, BW, MH) are silent on timing. DS (gain control integrated over about 110 ms; linear at 50 ms), Baker–Meese–Summers (interocular weight grows with duration), recurrent normalization (Heeger-type feedback; SSN; ORGaNICs) and the macaque data (Dougherty 2019; Cox 2019) all predict that the transient window behaves more like linear summation (larger k, larger additivity index, smaller w) and the sustained window more like ocularity-invariant normalization or a dip.

**Conflicts and whether the team's 6 × 6 neural matrix can decide them** (about 70 single units plus MUA/CSD; about 30 trials per condition; transient and sustained windows)

| Model contrast | Discriminating signature | Can a 6 × 6 neural matrix decide it? | Analysis |
|---|---|---|---|
| 1. Monocular normalization + linear sum (w ≈ 0; Truchard-like, Segala-cortex-like) vs binocular normalization (w ≈ 1; MH/BW) | High-contrast off-diagonal facilitation: about 2× vs about 1× | **Yes**, easily; this is the largest effect | Fit MH Eq. 10 (with eye-specific Rmax, c50 and n) per unit, with w free. Use hierarchical Bayesian estimation of w (Segala/Baker 2021 approach) by layer and ocular dominance. |
| 2. Exponent before vs after binocular summation (late vs early) | k = 2^(1/n) vs k = 2; additivity index at 5.5–11% below 1 vs above 1 | **Yes**, especially in units with c50 > 10% | Model-free check: compare cells (c, c) with (2c, 0) and (0, 2c), which the octave grid gives directly. Then fit Minkowski B = (L^k + R^k)^(1/k) inside a Naka-Rushton function; k = 1 is linear, k = n is the superellipse. |
| 3. Contrast-weighted averaging (two-stage, DS) vs constant-k binocular normalization | k(c) declining to 1 vs flat; dip R(90,y) < R(90,0) | **Probably, at population level.** The dip is only about 5–15%, so a single unit with 30 trials is marginal. | Pool normalized responses across units (and MUA). Test the sign of R(90, 5.5–45) − R(90, 0) and R(45, 5.5–22.5) − R(45, 0), and the slope of log k against log c. |
| 4. Two-stage (MGB) vs Ding–Sperling | Weights ∝ c^1 vs ∝ c^γ (γ ≈ 1.2); low-contrast details | **Probably not** from contrast alone; they differ mainly in phase predictions | Would need interocular phase or disparity variation (DS was built on phase). |
| 5. Interocular suppression as contrast gain vs response gain vs subtraction | Lateral shift vs vertical scaling vs downward offset of the eye-A CRF family across eye-B levels | **Partially.** Best in monocular-driven units, where eye B adds little excitation. The literature is itself inconclusive (Baker 2021). | Fit nested models (B in the denominator; B dividing Rmax; B subtracted) with model comparison, pooled hierarchically. Report per layer, because Sengpiel 1998 found layer differences. |
| 6. Divisive interocular normalization vs ocular-opponency or difference channels (Said & Heeger; subtractive) | Opponency: response ≈ f([L − R]_+), near zero on and above the diagonal; iso-response lines parallel to the diagonal. Divisive: iso-response rays (L/R = const) at high contrast, non-zero diagonal. | **Yes**, for classifying units | Fit iso-response contours, or compare a model in L − R with a model in L/R per unit. Look for units with a near-zero diagonal. |
| 7. OD-independent pool vs OD-weighted pool | c50 vs Rmax difference between the eyes' monocular CRFs; size of suppression by the non-dominant eye in monocular units | **Yes** | Regress fitted pool weights (b_N/b_D) on excitatory weights (a_N/a_D) across units. (i) predicts slope 0 with b_N/b_D ≈ 1; (ii) predicts slope 1. |
| 8. Disparity energy model vs contrast-only normalization | Phase/disparity-dependent cross-term; dioptic suppression in tuned-inhibitory cells | **No.** A zero-disparity matrix confounds disparity tuning with interocular normalization. | Needs each unit's disparity or interocular-phase tuning. Without it, treat "binocular suppression at the diagonal" as possibly reflecting a non-preferred zero-disparity phase. |
| 9. Feedforward vs recurrent normalization | Window dependence of w, k and the dip; latency of the suppression onset | **Partially**, using the transient vs sustained windows and PSTHs | Fit the same model separately to transient and sustained windows and test Δw and Δk. Estimate the onset of R(L,R) − R(L,0) differences relative to response onset (cf. 7 ms in Tanabe & Cumming; early facilitation then suppression in Dougherty/Cox). |
| 10. SSN (supralinear → sublinear) vs normalization | Additivity index > 1 at low c changing to < 1 at high c | **Weakly.** Early summation + Naka-Rushton predicts the same static crossover. | Needs dynamics or perturbation (e.g., paradoxical inhibitory responses), or a laminar CSD sink/source balance. The 6 × 6 alone can only show the crossover. |
| 11. Single-stage vs two-stage gain control | Stage-2 parameters (p, q, Z) act like a generic monotone output nonlinearity | **No** for a single unit (stage 2 cannot be told apart from spike-output nonlinearity). Iso-response contours depend only on stage 1 (m, S, w), as BMG 2007 note for matching. | Use contour shape to estimate stage 1. Compare layers: a "second stage" would appear as a change in output nonlinearity between input and output layers. |

**Rough power note (my estimate).** With mean counts around 20 spikes per window and Fano factor about 1–1.5, the SEM over 30 trials is about 5% of the mean. A condition difference therefore has an SE of about 6.5% per unit. That puts 5–8% dips at or below single-unit detectability, while w ≈ 0 vs 1 (up to 100% differences) and k = 2 vs 1.4 differences are detectable per unit. MUA and CSD will have better SNR for the subtler tests.

### Gaps
- No paper found that reports a full dichoptic contrast matrix (more than 4 × 4 levels) for isolated single units in awake macaque V1 with model comparison across these families. Hallum et al. (2017) come closest by modulating contrast continuously in both eyes, but in amblyopic animals with a linear-RF model, and their equations were unreadable.
- Smith et al. (1997): the exact contrast design and whether any saturation or normalization was fitted are not known from the abstract alone.
- The grid predictions use psychophysical or BOLD parameters. Neuron-specific parameters (c50, n) will change the response-ratio BSR values, though not k or the sign of the dip, for the families noted.

---

## Q4. Which parameters are identifiable from a full 6 × 6 dichoptic contrast matrix, which need other stimuli, and who has fitted full dichoptic surfaces?

### Takeaway
A same-orientation, zero-disparity 6 × 6 matrix identifies the following:
- each eye's monocular CRF (Rmax, c50, n);
- the interocular weight w;
- whether the exponent is applied before or after summation (via k);
- the monotonicity of the dichoptic surface (dip);
- the iso-response contour geometry (linear, superellipse, ratio or difference);
- for each neuron, whether ocular dominance enters the excitatory weights, the pool weights, or both.

It cannot separate the following, each of which needs extra stimulus dimensions:
- tuned binocular summation from untuned interocular suppression (needs orthogonal or anti-phase dichoptic masks);
- monoptic from interocular suppression routes (needs monoptic plaids or masks);
- disparity and phase effects (energy model);
- temporal-frequency or speed dependence;
- the role of adaptation;
- stage-2 nonlinearities from the neuron's output nonlinearity.

### Cited Findings
- **Full or partial dichoptic contrast surfaces already fitted:**
  - Human psychophysics, 5 × 5 L × R grid (0.3–32%), discrimination and matching, two-stage model with zero free parameters. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
  - Human psychophysics, binocular, monocular and dichoptic masking functions (dipper functions); six-parameter two-stage model. — [Meese, Georgeson & Baker 2006](https://doi.org/10.1167/6.11.7)
  - Human psychophysics, 48 phase × contrast-ratio conditions, one free parameter (γ). — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
  - Human psychophysics, perceived phase and perceived contrast across interocular contrast ratios. — [Ding, Klein & Levi 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4521337/)
  - Human SSVEP, 30 conditions with 4–64% contrast in each eye, including different-frequency (dichoptic masking) pairs. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
  - Human fMRI, two rows of the matrix (left fixed at 2% or 10%; right 0–50%), plus dichoptic and monocular plaids. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
  - Human EEG, pupillometry and matching for luminance flicker with a 48% dichoptic mask; hierarchical Bayesian fit of w. — [Segala et al. 2023](https://elifesciences.org/articles/87048)
  - Human MEG, noise stimuli at 0–32% contrast with dichoptic masking; interocular normalization fit. — [Chadnova et al. 2018](https://doi.org/10.1016/j.neuroimage.2017.10.035)
  - Cat simple cells, one eye fixed (high or low), the other varied over at least a 10-fold range, × 8+ interocular phases. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
  - Macaque V1/V2 arrays (amblyopic), both eyes' contrasts modulated independently between 0 and 40% (modulators at 0.5 and 1.67 Hz). — [Hallum et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566869/)
  - Macaque V1 neurons, binocular contrast combination ("combine linearly"). — [Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.1.366)
  - Cat V1 neurometric monocular vs binocular detection. — [Anzai et al. 1995](https://doi.org/10.1017/s0952523800007331)
  - Macaque superficial V1 (2-photon), only 0.9 contrast, monocular vs binocular. — [Zhang et al. 2024](https://elifesciences.org/reviewed-preprints/92839v2)
- **Parameters that need other stimuli:**
  - Tuned vs untuned suppression and summation need mask–target phase and orientation manipulations. — [Baker & Meese 2007](https://doi.org/10.1016/j.visres.2007.08.013)
  - Orthogonal dichoptic suppression has its own weight (w = 0.52 vs 1 for same orientation). — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
  - Monoptic vs dichoptic suppression differ in TF tuning and adaptation. — [Li et al. 2005](https://doi.org/10.1152/jn.00203.2005); [Sengpiel & Vorobyov 2005](https://doi.org/10.1523/jneurosci.0862-05.2005)
  - Monocular suppression weight depends on speed while interocular suppression does not. — [Meese & Baker 2009](https://doi.org/10.1167/9.5.2)
  - Duration dependence of the two suppression routes. — [Baker, Meese & Summers 2007](https://doi.org/10.1016/j.neuroscience.2007.01.030)
  - Interocular phase is needed to locate combination before vs after rectification. — [Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.1.366); [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
  - Adaptation and interocular transfer constrain where ocular dominance enters. — [Georgeson & Sengpiel 2021](https://doi.org/10.1016/j.visres.2021.03.004)
  - Interocular luminance differences need a pre-combination temporal filter. — [Chadnova et al. 2018](https://doi.org/10.1016/j.neuroimage.2017.10.035)
- **Identifiability within these models:**
  - In BMG 2007, matching depends only on stage-1 parameters (m, S), while discrimination depends on all six parameters. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
  - Truchard et al. could not separate the binocular gain site from the output exponent. — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
  - In Segala et al., Rmax was omitted for matching because "it has no effect in this paradigm". — [Segala et al. 2023](https://elifesciences.org/articles/87048)

### Inferences
- **Identifiable from the 6 × 6 matrix alone, per unit or pooled:**
  1. **Per-eye monocular CRF** from the L = 0 and R = 0 edges (5 nonzero points each): Rmax_i, c50_i, n_i, and from them the relation between ocular dominance and the CRFs.
  2. **Interocular weight(s)** w_{L→R} and w_{R→L} from the 25 interior cells. These can be asymmetric (cf. Shooner's asymmetry in amblyopia).
  3. **Exponent placement** from k and the additivity index on the diagonal relative to the edges.
  4. **Combination geometry** (Minkowski exponent; ratio vs difference) from iso-response contours. Under any model of the form R = F(B(L,R)) with monotone F, contours depend only on B. Contours are therefore the right summary to compare with psychophysical matching (BMG 2007) and are robust to unit-specific output nonlinearities.
  5. **Monotonicity (dip)**: a sign test.
  6. **Time-window dependence** of all of the above.
- **Degeneracies to expect:**
  - Rmax, the numerator weight and σ trade off (fix σ or normalize Rmax per unit).
  - Stage-2 (p, q, Z) cannot be separated from the output nonlinearity.
  - With only 5 nonzero levels per eye, p and q are weakly identified separately (super-saturation needs dense high-contrast sampling).
  - Contrast gain vs response gain for interocular suppression is identifiable only with good SNR in units where eye B drives little excitation.
- **Need other stimuli:**
  - (i) orthogonal or anti-phase dichoptic gratings, to split tuned summation from untuned suppression;
  - (ii) monoptic plaids or masks, to split ipsiocular from interocular normalization;
  - (iii) disparity or interocular phase, for the energy-model cross-term and the DS/DSKL phase predictions;
  - (iv) TF or speed;
  - (v) adaptation;
  - (vi) causal or pharmacological manipulations, to distinguish recurrent SSN regimes.
- **Specific to this design:** with the same grating presented to both eyes at zero disparity, the diagonal and off-diagonal cells mix binocular summation with disparity tuning. Units whose preferred disparity is far from zero, or which are tuned-inhibitory, can show "binocular suppression" for reasons unrelated to interocular normalization (energy-model cross-term). A conservative analysis restricts the strongest claims to units that are monocular or have documented near-zero disparity preference, or else reports results by ocularity class.

### Gaps
- I found no study fitting MGB, DS and normalization families head-to-head on single-neuron dichoptic contrast surfaces in awake macaque.
- The Smith et al. (1997) design details are unknown here.
- Whether the Baker-lab SSVEP chromatic binocular normalization study (PMC12007675, surfaced in search) adds fitted surfaces was not checked.

---

## Q5. Open controversies: monocular vs binocular gain control; divisive vs subtractive interocular suppression; linear vs nonlinear summation; mapping to V1 layers

### Takeaway
None of these controversies is settled. The conflicting results line up along measurement scale (single simple cells vs pooled BOLD/EEG), stimulus (gratings vs flicker vs noise), and time window. Evidence on cortical circuitry (interocular suppression is intracortical and GABAergic; monoptic suppression is precortical; laminar differences in gain mode; ocular-dominance effects in macaque) points to interocular normalization that emerges and strengthens beyond the geniculate input stage. Laminar recordings across the full dichoptic grid are positioned to test that directly.

### Cited Findings
**Monocular vs binocular gain-control stage**
- Cat simple cells: gain control is "primarily at a monocular site" (binocular gain change about 1.2×). — [Truchard et al. 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772209/)
- Human V1 BOLD: binocular normalization is required; monocular normalization followed by linear summation fails. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Human SSVEP: late-summation (binocular pool) model wins. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
- Gain control both before and after summation (two-stage). — [Meese, Georgeson & Baker 2006](https://doi.org/10.1167/6.11.7)
- Cortical EEG to flicker: w ≈ 0.02 (near-monocular normalization); pupil w ≈ 0.61. — [Segala et al. 2023](https://elifesciences.org/articles/87048)
- Early suppression signals and late, OD-affected gain control in cat V1. — [Georgeson & Sengpiel 2021](https://doi.org/10.1016/j.visres.2021.03.004)
- Perception reads the maximum of separate monocular and binocular channels. — [Georgeson et al. 2016](https://doi.org/10.1016/j.visres.2016.08.001)

**Divisive vs subtractive; contrast gain vs response gain**
- Interocular suppression lowers response gain in cat area 17 (cross-orientation shifts contrast gain), with layer differences. — [Sengpiel et al. 1998](https://doi.org/10.1016/s0042-6989(97)00413-6)
- Bayesian meta-analysis of 16 studies: contrast gain is best supported for overlay masks, but group HDIs for all mask types overlap 1, so the literature is inconclusive. New EEG at 5 Hz: monocular and dichoptic overlay masks both show contrast gain at Oz. At 10 Hz, both contrast and response gain appear. The authors suggest that all suppression routes may serve to suppress less reliable inputs. — [Baker, Vilidaite & Wade 2021, PLoS Comput Biol 17:e1009507](https://pmc.ncbi.nlm.nih.gov/articles/PMC8544832/)
- SSVEP eye-experiment misfit resembles response gain. — [Baker & Wade 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903417/)
- Subtractive feedback from ocular-opponency neurons on top of divisive normalization explains rivalry vs plaid differences. — [Said & Heeger 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3610603/)
- Two inhibitory interneuron classes in a dynamical model unify contrast combination, fusion and rivalry. — [Wilson 2017](https://doi.org/10.1016/j.visres.2017.07.016)

**Linear vs nonlinear summation**
- Linear combination in macaque V1 neurons, before rectification in simple and phase-specific complex cells. — [Smith et al. 1997](https://doi.org/10.1152/jn.1997.78.1.366)
- Near-linear transduction before summation (m ≈ 1.2–1.5); linear and quadratic summation rejected for matching. — [Baker, Meese & Georgeson 2007](https://publications.aston.ac.uk/id/eprint/1426/1/AURA1426.pdf)
- Threshold summation m varies with speed (BSR ≈ 1.5 on average, above √2). — [Baker et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6195301/)
- Linear at low contrast or brief presentations, contrast-weighted at high contrast. — [Ding & Sperling 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347993/)
- Power-law and linear rules rejected in BOLD. — [Moradi & Heeger 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2747122/)
- Near-linear cortical combination of luminance flicker. — [Segala et al. 2023](https://elifesciences.org/articles/87048)

**Mapping onto V1 layers and circuitry**
- Monoptic cross-orientation suppression is subcortical or thalamocortical. — [Freeman et al. 2002](https://doi.org/10.1016/s0896-6273(02)00819-x); [Li et al. 2005](https://doi.org/10.1152/jn.00203.2005)
- Interocular suppression is intracortical and GABA-dependent. — [Sengpiel & Vorobyov 2005](https://doi.org/10.1523/jneurosci.0862-05.2005)
- Ipsiocular suppression may arise in retina or LGN; interocular suppression is cortical and predicted to act in cortical cells driven mainly by one eye. — [Baker, Meese & Summers 2007](https://doi.org/10.1016/j.neuroscience.2007.01.030)
- Monocular neurons in the input layers are binocularly modulated (early facilitation, later suppression). — [Dougherty et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6363852/)
- Sustained-phase binocular modulation depends on laminar location and ocularity. — [Cox et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6797477/)
- Layer-dependent gain modes: layer 4 response gain, other layers contrast gain, as summarized secondarily. — [Baker, Vilidaite & Wade 2021, citing Sengpiel et al. 1998](https://pmc.ncbi.nlm.nih.gov/articles/PMC8544832/)
- A 2024 Vanderbilt dissertation (Mitchell) proposes a single divisive normalization acting within and between the eyes, whose strength grows through the cortical microcircuit, with ocular dominance weighting the effectiveness of non-dominant-eye contrast (search-result summary). — [Mitchell 2024, Vanderbilt IR](https://ir.vanderbilt.edu/handle/1803/19176)
- Superficial-layer monocular neurons are suppressed and binocular neurons enhanced by binocular stimulation. — [Zhang et al. 2024](https://elifesciences.org/reviewed-preprints/92839v2)

### Inferences
- **Monocular vs binocular gain control may depend on scale and cell class.** Truchard's simple cells (largely layer 4 in cat) showing mostly monocular gain control is compatible with BOLD and EEG showing binocular normalization if interocular normalization strengthens after the input stage. Pooled signals are dominated by superficial and output layers.
  - Prediction for laminar data: w (and ocularity invariance, k → 1) increases from granular to supragranular and infragranular layers, and from the transient to the sustained window.
  - The 6 × 6 matrix plus CSD-defined layers can test this directly. Fit w per layer and per window.
- **Divisive vs subtractive.** These are distinguishable on the 6 × 6 grid by contour geometry: ratio-like rays vs difference-like diagonals. Contrast gain vs response gain is distinguishable by how eye-A CRFs change across eye-B levels. Baker 2021's inconclusive meta-analysis suggests per-unit tests will be underpowered, so use hierarchical pooling and pre-registered nested comparisons. Sengpiel 1998 found layer differences, so analyse by layer rather than pooling across layers.
- **Linear vs nonlinear summation are not mutually exclusive across the response time course.** Ding–Sperling (linear at 50 ms) and the macaque early-facilitation data suggest summation is close to linear (k ≈ 2, additivity index ≥ 1) in the onset transient and becomes contrast-weighted or normalized (k → 1, possible dip) in the sustained response. A window-specific fit of the same model is the most direct test available with this dataset.
- **Where the competing claims point.** Interocular suppression appears to be cortical, GABAergic, low-TF, adaptable and slower to build up (psychophysics and cat physiology). Monoptic suppression appears to be fast, broadband and precortical. That predicts:
  - (i) off-diagonal suppression with a later onset than monocular responses, in the transient-vs-sustained difference;
  - (ii) weaker interocular suppression in the earliest granular responses;
  - (iii) suppression by the non-dominant eye in monocular units (Baker–Meese–Summers prediction; Dougherty and Zhang findings).

### Gaps
- No direct laminar single-unit test of w, k or the dichoptic dip across the full dichoptic contrast matrix was found in the literature (outside the team's own work, which I did not research).
- The layer-specific gain-mode claim from Sengpiel et al. (1998) rests on the abstract plus Baker et al.'s (2021) secondary description. The primary layer-by-layer numbers were not read.
- No source was found that directly tests divisive vs subtractive interocular suppression with full contrast surfaces in neurons.
- Binocular SSN and recurrent-normalization (ORGaNICs) predictions for dichoptic contrast surfaces appear unpublished in peer-reviewed form. Only a master's thesis was found for SSN.
