# SSN / ISN theory, dynamics, and laminar and binocular models: notes for relating laminar dichoptic V1 data to SSN predictions

Compiled 2026-10-29. Sources were read as full text through Europe PMC / PMC XML where available (Ahmadian 2013, Rubin 2015, Hennequin 2018, Echeveste 2020, Ahmadian & Miller 2021, Ozeki 2009, Sanzeni 2020, Kraynyukova 2018, Nassi 2015, De 2020, Obeid & Miller 2025, Holt 2024). For the rest, only the abstract or a summary was available, as noted. Tags: **[preprint]** = not peer reviewed as of the record found; **[paywalled]** = the team may need to supply the PDF; **[prior-knowledge citation]** = standard reference not re-fetched in this session, so verify the details before quoting.

The team's data features referred to throughout as (1)–(6):
(1) near-linear interocular combination for the first ~20 ms, then normalization within ~20–35 ms;
(2) binocular summation collapses further and faster at higher contrast;
(3) a feedforward delayed-normalization model with a fixed pool delay is too slow and too shallow at high contrast;
(4) normalization timing and strength are similar across depth, but superficial sites start out less linear;
(5) no contrast supersaturation at any depth;
(6) responses at 5.5% contrast may be near threshold.

---

## Q1. Core SSN/ISN theory: equations, supralinear vs sublinear summation, and supersaturation

### Takeaway
The SSN is an E/I rate network whose units have an expansive power-law input/output function. The effective coupling therefore grows with activity. Weak input gives a weakly coupled, feedforward-like network in which responses sum supralinearly. Strong input gives an inhibition-stabilized, "loosely balanced" network in which E and I responses both sum sublinearly (normalization-like). In the two-population analysis, excitatory rates eventually *decrease* with input (supersaturation) if and only if $\Omega_E = J_{II}g_E - J_{EI}g_I < 0$. The absence of supersaturation therefore constrains the ratio of feedforward drive to I versus E relative to the I→E and I→I weights.

### Cited Findings
**Model equations**
- The canonical rate model is
  $\tau_E \dot r_E = -r_E + k[W_{EE}r_E - W_{EI}r_I + h_E]_+^{\,n}$ and $\tau_I \dot r_I = -r_I + k[W_{IE}r_E - W_{II}r_I + h_I]_+^{\,n}$, with $n>1$. Rubin et al. use "k = 0.04, n = 2.0" for the power-law I/O function. Obeid & Miller's L2/3 model uses the same form with $n_E=n_I=2.2$, $\tau_E = 10$ ms, $\tau_I = 6.67$ ms and $K = 0.01$. — [Rubin, Van Hooser & Miller 2015, Neuron 85:402–417, doi:10.1016/j.neuron.2014.12.026 (PMC4344127)](https://doi.org/10.1016/j.neuron.2014.12.026); [Obeid & Miller 2025, eNeuro 12(5), ENEURO.0459-24.2025](https://doi.org/10.1523/ENEURO.0459-24.2025)
- Why a power law: it "arises when the mean input to neurons is subthreshold, so that neurons fire on input fluctuations about the mean". In V1, "neurons remain in the unsaturated, power-law region of the I/O function throughout the full range of firing induced by visual stimuli, with powers in the range 2–5". — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- With external input $h = c\,g$ (contrast $c$ times a fixed input-weight vector $g$), the steady state $r = k[Wr + cg]_+^n$ can be rewritten exactly as $r = -cW^{-1}g + k^{-1/n}W^{-1}r^{\cdot 1/n}$ (their Eq. 15). In dimensionless form it is $y = -J^{-1}g + \beta J^{-1} y^{\cdot 1/n}$ (Eq. 13), where $\beta \to 0$ in the strongly driven limit. The leading term $r_0 = -cW^{-1}g$ is the (loosely) balanced solution, which "is of course identical to the condition we have found for approximate dynamic cancellation in the SSN". — [Ahmadian, Rubin & Miller 2013, Neural Comput 25:1994–2037, doi:10.1162/NECO_a_00472 (PMC4026108)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- Condition for normalization: "the condition for normalization can alternatively be expressed as the requirement that the recurrent input Wr grow more slowly than linearly with c." — [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- Stability: "dynamic stabilization will occur whenever the determinant of the weight matrix is positive and the inhibitory time constant is sufficiently small." — [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)

**Supralinear-to-sublinear transition**
- "responses to weak inputs sum supralinearly while those to stronger inputs generically sum sublinearly". "For a wide range of network and stimulus parameters, this dynamic stabilization yields a transition from supralinear to sublinear summation." — [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- "when individual neurons have supralinear input/output functions, inhibitory stabilization drives a transition from weak coupling and supralinear response summation for weak inputs to ISN behavior and sublinear summation for strong inputs." "Both E- and I-cell neuronal responses sum sublinearly, an emergent outcome of network dynamics, as opposed to the more intuitive scenario that suppression in E cells results from increased I-cell firing." — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- Normalization phenomenology reproduced by the SSN: for two gratings of equal contrast, the population response is "a sublinear multiple (~0.5 to 0.7) of the sum". As contrasts become unequal, the response approaches "winner-take-all", and summation becomes "increasingly winner-take-all" with increasing difference in stimulus strengths. The model "for the first time provides a network explanation of normalizing and winner-take-all behavior of both E and I cells". From MT data: "average summation was linear when at least one stimulus had contrast below that which drove half-maximal response". — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- In the stochastic spiking SSN, "population mean Vm increased monotonically with growing external input, first linearly or supralinearly for small inputs, but strongly sublinearly for larger inputs, with V̄I growing faster than V̄E". The faster growth of I activity is "a robust outcome of dynamic stabilization by feedback inhibition" and has been observed in rodent S1 (Shao et al. 2013) and V1 (Adesnik 2017). — [Hennequin et al. 2018, Neuron 98:846–860, doi:10.1016/j.neuron.2018.04.017 (PMC5971207)](https://doi.org/10.1016/j.neuron.2018.04.017)
- Loosely balanced regime: a hallmark is that "the net input a neuron receives grows sublinearly as a function of its external input". The network "flow[s] to this loosely balanced solution, without any need for fine tuning of parameters". The authors argue the evidence "is most consistent with a loosely balanced regime" for cortex. — [Ahmadian & Miller 2021, Neuron 109:3373–3391, doi:10.1016/j.neuron.2021.07.031 (PMC9129095)](https://doi.org/10.1016/j.neuron.2021.07.031)

**ISN definition and paradoxical effect**
- An ISN is a network "in which excitatory recurrence alone is strong enough to destabilize visual responses but feedback inhibition maintains stability". In nullcline terms, the two nullclines "intersect on the positive-sloping portion of the E nullcline". — [Ozeki, Finn, Schaffer, Miller & Ferster 2009, Neuron 62:578–592, doi:10.1016/j.neuron.2009.03.028 (PMC2691725)](https://doi.org/10.1016/j.neuron.2009.03.028)
- The paradoxical effect in ISNs (extra excitatory drive to I cells lowers their steady-state rate) originates with Tsodyks et al. — [Tsodyks, Skaggs, Sejnowski & McNaughton 1997, J Neurosci 17:4382–4388, doi:10.1523/JNEUROSCI.17-11-04382.1997 [prior-knowledge citation]](https://doi.org/10.1523/JNEUROSCI.17-11-04382.1997)
- In the SSN the network passes from non-ISN to ISN as input grows: "corresponding to a transition from non-ISN to ISN behavior". "The SSN shows ISN-regime operation for strong visual stimuli, but predicts that as sensory stimuli decrease in strength, the network should transition into a non-ISN state." — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/); [Sanzeni et al. 2020, eLife 9:e54875](https://doi.org/10.7554/eLife.54875)

**Supersaturation conditions (two-population SSN)**
- Definitions: $\Omega_E \equiv \mathrm{Det}J\,(-J^{-1}g)_E = J_{II}g_E - J_{EI}g_I$ and $\Omega_I \equiv \mathrm{Det}J\,(-J^{-1}g)_I = J_{IE}g_E - J_{EE}g_I$. Here $J_{EI}, J_{II}$ are the magnitudes of the inhibitory weights, so $J = \begin{pmatrix}J_{EE}&-J_{EI}\\J_{IE}&-J_{II}\end{pmatrix}$ and $\mathrm{Det}J = J_{EI}J_{IE}-J_{EE}J_{II} > 0$. Three cases are possible: (1) both elements of $-J^{-1}g$ positive; (2) E-element negative and I-element positive; (3) both negative. The fourth combination is impossible for $g_E, g_I \ge 0$. — [Ahmadian et al. 2013, Eqs. 26–27](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- When the E element of $-J^{-1}g$ is negative ($\Omega_E<0$), "the sublinear term becomes dominant ... and network behavior becomes strongly sublinear". That is the supersaturation case ("5.2 The case (−J⁻¹g)E < 0 and supersaturation"). — [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- "The strongest sublinear behavior, and hence behavior most likely to underly biological observations in cerebral cortex, occurs for parameters that lead to supersaturation." Yet "the model shows a broad dynamic range between the onset of normalization and of supersaturation for most parameter choices we have explored, the only exception ... being the case ΩE < 0 and ΩI > 0 for small JEE." — [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- Related criterion from a ring model with LGN input widths: supersaturation when $J_{EI}\sigma_{I,LGN} - J_{II}\sigma_{E,LGN} > 0$, which Ahmadian et al. say is "closely related to our criterion ΩE < 0" (Persi, Hansel, Nowak, Barone & van Vreeswijk 2011, PLoS Comput Biol, doi:10.1371/journal.pcbi.1001078). — [cited within Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- Restated by Obeid & Miller: "for ΩE < 0, as external input grows from zero, excitatory responses first grow supralinearly, but their growth then becomes sublinear and saturates (and ultimately would supersaturate) for increasingly strong external input." Their L2/3 model drives the network through a *saturating* Naka–Rushton contrast input $f(C) = f_{max}C^{3.5}/(C_{50}^{3.5}+C^{3.5})$ with $f_{max}=50$ and $C_{50}=11$. — [Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025)

**Extensions**
- An SSN can also show bistability, persistent activity after input withdrawal (with an exact connectivity condition), and a supercritical Hopf bifurcation that produces global oscillations. — [Kraynyukova & Tchumatchenko 2018, PNAS 115:3464–3469, doi:10.1073/pnas.1700080115 (PMC5879648)](https://doi.org/10.1073/pnas.1700080115)
- Multiple interneuron types (PV, SOM, VIP): the theory gives the conditions under which such E/PV/SOM/VIP networks are stable and shows how perturbing each subtype propagates through the recurrence. It is applied to disinhibition, surround suppression, and subtractive vs divisive tuning modulation. — [Litwin-Kumar, Rosenbaum & Doiron 2016, J Neurophysiol 115:1399–1409, doi:10.1152/jn.00732.2015](https://doi.org/10.1152/jn.00732.2015)
- With a single inhibitory class, "network stabilization and gain control work in opposition", so high gain coincides with low stability. Multiple interneuron classes can untangle the two. — [Bos, Miehl, Oswald & Doiron 2025, eLife, doi:10.7554/eLife.99808](https://doi.org/10.7554/eLife.99808)
- Palmigiano, Miller et al. work: statistical models of V1 recurrent circuitry were fit to each cell type's visual responses, and the distributions of optogenetic responses were linked to circuit stability. "The drive to VIP cells ... is largely offset by the drive to SOM cells" during locomotion. "Small mean effects of optogenetic perturbation seen in monkey V1, compared with mouse, could be explained by stronger recurrent connectivity in the monkey." (Source is a 2022 seminar abstract; see Gaps.) — [Stanford seminar abstract, A. Palmigiano](https://neuroscience.stanford.edu/node/2842)

### Inferences
- **Effective gain grows with rate.** From $r = k u^n$ it follows that $dr/du = nku^{n-1} = n k^{1/n} r^{(n-1)/n}$. The effective (linearized) connectivity $G W$, with $G = \mathrm{diag}(n k^{1/n} r_i^{(n-1)/n})$, therefore strengthens as activity rises. This is the single mechanism behind every input-strength-dependent SSN effect: weak→strong coupling, non-ISN→ISN, supralinear→sublinear, and the faster dynamics discussed in Q2. (Derivation; consistent with Hennequin 2018 and Holt et al. 2024.)
- **Dichoptic mapping.** The natural SSN formulation of the 6×6 matrix is $h = c_L g_L + c_R g_R$ (or $f(c_L)g_L + f(c_R)g_R$ with a saturating LGN/L4 front end, as in Obeid & Miller). For a single binocular unit in the weakly coupled regime, $k(a+b)^n > ka^n + kb^n$ for $n>1$, so the SSN predicts *supralinear* (not linear) interocular summation for binocular cells at weak drive. Near-linear early MUA summation (1) is therefore most consistent with the SSN if early MUA is dominated by monocularly driven neurons (e.g., L4C), or if the earliest phase reflects input currents rather than the rate nonlinearity. This is an inference, not a published prediction.
- **(5) No supersaturation constrains parameters.** In the two-population SSN, the absence of supersaturation up to 90% contrast implies either
  (a) $\Omega_E > 0$, i.e. $g_E/g_I > J_{EI}/J_{II}$ (feedforward drive to E large relative to drive to I, scaled by the I→E/I→I weight ratio), or
  (b) $\Omega_E < 0$ with the response peak beyond the tested contrast range (possible, because the dynamic range between normalization onset and supersaturation is usually broad, except when $\Omega_E<0$, $\Omega_I>0$ and $J_{EE}$ is small).
  Under (a), E rates approach the balanced solution $r_E \approx c\,\Omega_E/\mathrm{Det}J$, which grows *linearly* in $c$. Saturating CRFs then have to come from saturating input (a Naka–Rushton front end, as in Obeid & Miller) or from structured cross-population terms in $W^{-1}$, not from the local E/I loop alone.
  There is a tension: Ahmadian et al. say the strongest sublinearity occurs for supersaturating parameters, while the team sees strong interocular normalization (2) without supersaturation (5). That combination favours (i) $\Omega_E$ near or above 0 plus input saturation, or (ii) normalization generated by structured (e.g., ocular-dominance-column) interactions rather than the unstructured 2-population loop.
- **(6) Near threshold.** At 5.5% the SSN predicts weak coupling, non-ISN behaviour, and supralinear or at most linear binocular summation. A binocular summation ratio ≥ linear at (5.5%, 5.5%) would be an SSN-consistent signature. Strictly sub-additive summation at near-threshold contrast would be harder to reconcile with the SSN, though a fixed feedforward divisive model allows it.
- Per-paper lines for (1)–(6):
  - **Ahmadian 2013 / Rubin 2015 (generic SSN).** (1) early phase ≈ feedforward, weakly coupled; (2) yes, sublinearity grows with input; (3) normalization is emergent and input-dependent, so no fixed delay; (4) silent on layers (single layer); (5) requires Ω_E > 0 or a late peak; (6) predicts supralinear summation at 5.5%.
  - **Hennequin 2018.** Same as above, plus I grows faster than E, which predicts that putative inhibitory units normalize earlier and more strongly.
  - **Kraynyukova 2018.** Relevant only if bistability or rivalry-like alternation appears. Not needed for (1)–(6).
  - **Litwin-Kumar 2016 / Bos 2025.** Multiple interneuron classes allow normalization (gain) and stability to be set separately. This could explain strong normalization without supersaturation (5), but it adds free parameters.

### Gaps
- Not verified: the exact definition of the dimensionless parameter β (and the related α, ψ) in Ahmadian et al. 2013, because the equations in the extracted text are garbled. Only "β → 0 in the strong-input limit" and the series expansion in β could be confirmed.
- Not verified: Rubin et al. 2015's $\tau_E, \tau_I$ values (only k = 0.04, n = 2.0 were extracted).
- The Palmigiano et al. preprint ("Structure and variability of optogenetic responses identify the operating regime of cortex", bioRxiv 2020.11.11.378729, Palmigiano, Fumarola, Mossing, Kraynyukova, Adesnik & Miller) could not be fetched (HTTP 429), and Europe PMC did not index it under that title. Its claims above rest on a seminar abstract. Its peer-reviewed descendant appears to be Sanzeni et al. 2023 (Q5). The team should obtain the preprint PDF.
- There is no published closed-form condition for supersaturation in multi-population (OD-column or orientation-ring) SSNs beyond the Persi et al. ring criterion.

---

## Q2. SSN dynamics: transients, speed of normalization vs input strength, contrast-dependent dynamics

### Takeaway
The SSN literature consistently predicts that stronger drive gives stronger effective inhibitory feedback, which in turn gives *faster* relaxation, larger inhibition-dominated onset transients, and higher-frequency damped (gamma-band) oscillations. Variability quenching happens on roughly a single-neuron time constant, and offset decay follows cellular time constants. So yes: an SSN normalizes faster at stronger input. However, no paper found computes "latency to sublinear summation vs contrast" directly, and recurrent shunting normalization (non-SSN) also predicts contrast-dependent speeding.

### Cited Findings
- "Given the supralinear input/output functions of cortical neurons, increased stimulus drive strengthens effective network connectivity. This shifts the balance from interactions that amplify variability to suppressive inhibitory feedback, quenching correlated variability around more strongly driven steady states." The SSN mechanism "uniquely accounts for the spatial patterns and fast temporal dynamics of variability suppression". — [Hennequin et al. 2018](https://doi.org/10.1016/j.neuron.2018.04.017)
- Timescales in Hennequin et al. 2018:
  - "In the SSN, the timescales on which both suppression and recovery of variability occurred were nearly as fast as the single-neuron time constant (20 ms in these simulations)". In chaotic networks these timescales were "several (4–15) times longer".
  - Strong input "reduced the LFP power at low frequencies, due to enhanced effects of feedback inhibition; increased it at intermediate frequencies, due to the faster timescales associated with relatively enhanced inhibition".
  - At stimulus offset "activity rapidly decayed back to spontaneous levels with the cellular time constant".
  — [Hennequin et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5971207/)
- Loose balance "can produce negative eigenvalues large enough to produce quite fast dynamics, with effective time constants on the order of single neurons' membrane time-constant, or even as small as a few milliseconds, depending on parameters". Tight balance gives O(K) negative eigenvalues, "approaching instantaneous population response". — [Ahmadian & Miller 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9129095/)
- Gamma frequency vs contrast:
  - "As the transfer function is supralinear, this translates to higher neural gains and stronger effective connections", so "as contrast increases, we expect the gamma peak in the LFP power-spectrum to move to higher frequencies".
  - The frequency is Im(λ)/2π of the linearized E–I Jacobian; Re(λ) sets damping.
  - Spontaneous activity sits near threshold where gains are low, so the network "generically does not oscillate".
  - Without noise, damped oscillations are transient. Latency vs contrast was not analyzed.
  - Only a *columnar* model (a sharp local E boost plus long-range E; inhibition confined to ~0.4 mm columns) reproduces both the local contrast dependence of gamma and surround suppression.
  — [Holt, Miller & Ahmadian 2024, PLoS Comput Biol 20(6):e1012190, doi:10.1371/journal.pcbi.1012190](https://doi.org/10.1371/journal.pcbi.1012190) (preprint: [bioRxiv 10.1101/2023.05.11.540442](https://doi.org/10.1101/2023.05.11.540442))
- An E/I network trained for sampling-based inference "displayed several key biological properties, including divisive normalization, as well as stimulus-modulated noise variability, inhibition-dominated transients at stimulus onset, and strong gamma oscillations". Other findings:
  - "Both the magnitude of transients and the frequency of gamma oscillations scaled with stimulus contrast."
  - Population rates "had marked contrast-dependent overshoots at stimulus onset".
  - "inhibition transiently dominated over excitation during stimulus presentation".
  - The prediction that "onset transients should be tuned to stimuli" was confirmed in published awake-monkey V1 recordings.
  - A control network without Dale's law lacked oscillations and transients.
  — [Echeveste, Aitchison, Hennequin & Lengyel 2020, Nat Neurosci 23:1138–1149, doi:10.1038/s41593-020-0671-1 (PMC7610392)](https://doi.org/10.1038/s41593-020-0671-1)
- ISN temporal signature: "prior to reaching its steady state, ISNs (and not non-ISNs) should show a transient increase in inhibition". The I rate "will therefore initially increase before it ultimately decreases toward the new center-plus-surround fixed point". Experimentally, when a stimulus abruptly expanded from center to center+surround, "each cell showed a transient increase in inhibition of 30–50 ms duration" (6 cat V1 cells with stable conductance measurements). — [Ozeki et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725/)
- In the L2/3 SSN, after thalamic input removal, E decay time constants are ≈10 ms and I ≈6.67 ms (equal to the cellular τ) "regardless of stimulus size or contrast". This reproduces the fast decay seen in mouse V1 when thalamic input is silenced. — [Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025)
- Spatial-frequency analogue: "We tested the prediction that network resonant frequencies increase with contrast" by measuring contrast-modulated preferred spatial frequencies in ferret V1. A nonzero I-unit resonant spatial frequency "requires an ISN". — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- Review linking sub-additivity and variability quenching: across modelling approaches, "response sub-additivity and variability quenching often co-occur". — [Goris, Coen-Cagli, Miller, Priebe & Lengyel 2024, Nat Rev Neurosci, doi:10.1038/s41583-024-00795-0 [paywalled]](https://doi.org/10.1038/s41583-024-00795-0)
- Recurrent (non-SSN) normalization models also speed up with contrast:
  - Carandini & Heeger's shunting model: the normalization conductance grows with contrast, shortening the membrane time constant. — [Carandini & Heeger 1994, Science 264:1333–1336, doi:10.1126/science.8191289 [prior-knowledge citation; paywalled]](https://doi.org/10.1126/science.8191289); [Carandini, Heeger & Movshon 1997, J Neurosci 17:8621–8644, doi:10.1523/JNEUROSCI.17-21-08621.1997 [prior-knowledge citation]](https://doi.org/10.1523/JNEUROSCI.17-21-08621.1997)
  - ORGaNICs is a recurrent circuit that implements normalization and simulates V1 response dynamics. — [Heeger & Zemlianova 2020, PNAS 117:22494–22505, doi:10.1073/pnas.2005417117 [prior-knowledge citation]](https://doi.org/10.1073/pnas.2005417117)
  - The feedforward "delayed gain control" model uses a low-pass-filtered normalization pool. — [Zhou, Benson, Kay & Winawer 2019, PLoS Comput Biol 15:e1007484, doi:10.1371/journal.pcbi.1007484 [prior-knowledge citation]](https://doi.org/10.1371/journal.pcbi.1007484)

### Inferences
- **Why stronger input normalizes faster.** Linearizing about the fixed point gives $\dot{\delta r} = T^{-1}(-I + GW)\,\delta r$ with $G \propto r^{(n-1)/n}$. As contrast rises, $G$ grows, and in the ISN regime the E–I "sum" mode gets more negative eigenvalues (faster decay) and a larger imaginary part (faster, higher-frequency ringing; Holt 2024). The time for inhibition to catch up with excitation, which is when summation turns sublinear, should therefore shrink with total drive. This matches (2), and it explains (3): a model with a fixed pool delay cannot shorten its delay at high contrast.
- **Early-linear-then-normalized (1).** At onset, rates start near spontaneous levels, where gains are low (Holt 2024). The earliest response is therefore a weakly coupled, essentially feedforward copy of the input, so near-linear interocular combination for the first ~20 ms is the expected SSN behaviour. Normalization then appears once recurrent gain has built up. With $\tau_E \approx 10$–20 ms and $\tau_I \approx 6$–10 ms, the time to the normalized state is plausibly 1–3 τ, i.e. ~15–40 ms. That brackets the observed ~20–35 ms and is in line with Ozeki's 30–50 ms inhibitory transient and Hennequin's quenching timescale of about one τ. (Order-of-magnitude inference, not a fitted prediction.)
- **Low contrast (6).** At 5.5%, low gain predicts slower settling, no or weak overshoot, no gamma resonance (spontaneous-like low gain), and linear-to-supralinear summation.
- **Speed-up alone is not diagnostic.** Shunting normalization models (Carandini–Heeger; ORGaNICs) also speed up with contrast. To tell them apart, use the signatures that depend specifically on an expansive nonlinearity plus inhibitory stabilization: supralinearity at low contrast, paradoxical or transient inhibitory dynamics, contrast-dependent ringing frequency, and both E and I summing sublinearly (Q6).
- Per-paper lines for (1)–(6):
  - **Hennequin 2018.** (1) and (2) yes, via gain growth; (3) predicts timescales of about τ, not a fixed delay; (4) no layers; (5) inherits the Ω_E constraint; (6) predicts that variability quenching is weakest at 5.5%.
  - **Echeveste 2020.** (1)–(2): an onset overshoot that grows with contrast. Test whether the binocular overshoot exceeds the monocular one at matched per-eye contrast. (4)–(6) are not addressed.
  - **Holt 2024.** Predicts gamma peak frequency rising with total (binocular) drive. Gamma should be absent or weak at 5.5% (6). The columnar structure is not laminar.
  - **Ozeki 2009.** A transient rise in inhibition before steady-state suppression. With simultaneous dichoptic onset this mixes with the onset transient, so it is testable only with onset-asynchronous or added-eye designs.
  - **Obeid & Miller 2025.** Predicts contrast-independent offset decay at ~τ. That is a cheap test in the existing ~317 ms on/off data.

### Gaps
- No published SSN analysis of how the *latency or time constant of the supralinear→sublinear transition* scales with input strength was found. The predictions above come from the eigenvalue argument plus the qualitative statements in Hennequin 2018 and Holt 2024, so a simulation by the team is needed.
- No SSN paper simulates dichoptic onset dynamics or the time course of interocular summation.
- Echeveste 2020's network is an E/I rate network with Dale's law trained for sampling. I did not re-verify in the extracted text that its units use a threshold-power-law (SSN) nonlinearity. I believe they do, but check the Methods.

---

## Q3. SSN/ISN models across the cortical column and layers

### Takeaway
No published SSN with explicit, interacting L4 ↔ L2/3 ↔ L5/6 populations and laminar timing predictions was found. The canonical SSN papers model a single layer (L2/3) and ignore interlaminar processing. The closest laminar SSN work is:
- Mossing et al. 2021 [preprint]: a four-cell-type L2/3 SSN driven by measured L4 input in mouse.
- Eskikand et al. 2022 [preprint]: a laminar E/I neural-mass model in which L4 and L5, but not L2/3, are ISN.
- Holt et al. 2024: a columnar (horizontal), not laminar, model.

The empirical laminar evidence is mixed: in mouse, ISN signatures appear in both superficial and deep layers, and in primates inhibitory control is layer-specific.

### Cited Findings
- Rubin et al. 2015: "We model interactions in a single layer, e.g. layer 2/3 (L2/3), ignoring interlaminar processing. We assume that the net effect of externally-driven input ... to this layer is excitatory." — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- **Obeid & Miller L2/3 SSN.**
  - Architecture: a 75×75 grid of E and I cells covering 8×8 mm with an orientation map. E cells make short- and long-range connections, I cells short-range only, with E→I spreading farther than E→E ($\sigma_{IE}=6$ vs $\sigma_{EE}=3$).
  - Surround suppression: E is suppressed more than I (high-contrast suppression index 0.79 vs 0.27). Suppression strengthens with contrast and summation fields shrink. Both excitatory and inhibitory input currents decrease under surround suppression, as in an ISN. The model predicts that E→I connectivity grows relative to E→E with distance.
  - L4 is not modelled. The authors note that feedback loops via higher areas must satisfy the same net connectivity conditions.
  — [Obeid & Miller 2025, eNeuro](https://doi.org/10.1523/ENEURO.0459-24.2025); preprint "Stabilized Supralinear Network: Model of Layer 2/3 of the Primary Visual Cortex" [bioRxiv 10.1101/2020.12.30.424892](https://doi.org/10.1101/2020.12.30.424892)
- **Mossing et al. [preprint].** "We combined calcium imaging and electrophysiology to constrain a stabilized supralinear network model that explains how the four principal cell types in layer 2/3 (L2/3) of mouse V1 — SST, PV, VIP interneurons, and PCs — transform inputs from layer 4 (L4) PCs to encode drifting gratings of varying size and contrast." "Recurrent amplification drives a transition from a positive PC→VIP⊣SST⊣PC feedback loop at small size and low contrast to a negative PC→SST⊣PC feedback loop at large size and high contrast". Key predictions were confirmed with bidirectional optogenetics. — [Mossing, Veit, Palmigiano, Miller & Adesnik 2021, bioRxiv 10.1101/2021.03.31.437953](https://doi.org/10.1101/2021.03.31.437953)
- **Eskikand et al. [preprint].** A "model of populations of neurons across different layers of cortex", with one E and one I population per layer and connectivity from the Allen Institute synaptic physiology database.
  - "superficial layers are more inhibitory dominated compared to deeper layers".
  - "layer 2/3 in the model does not operate in the ISN regime but layers 4 and 5 do operate in the ISN regime".
  - "the strength of the paradoxical response increases as the model moves closer to bifurcation points".
  It is a neural-mass model, not a power-law SSN. — [Eskikand, Soto-Breceda, Cook, Burkitt & Grayden 2022, bioRxiv 10.1101/2022.12.09.519705](https://doi.org/10.1101/2022.12.09.519705)
- **Holt et al. columnar SSN.** Horizontal E connections are "a sharp local (intra-columnar) term plus a smooth long-range term", inhibition is short-range, and columns are ~0.4 mm. Laminar structure is not modelled. — [Holt et al. 2024](https://doi.org/10.1371/journal.pcbi.1012190)
- **Mouse laminar ISN evidence.** "We find clear evidence for ISN operation even without sensory stimulation in multiple mouse cortical areas and both in superficial and deep layers." The main analyses were restricted to ≈400 µm (L2/3 and upper L4). Surface light also produced paradoxical suppression in deeper inhibitory neurons, although "at our highest stimulation intensities we did not always see clear mean increases in inhibitory firing rates" at deep sites. — [Sanzeni, Akitake, Goldbach, Leedy, Brunel & Histed 2020, eLife 9:e54875](https://doi.org/10.7554/eLife.54875)
- **Caveat for any layer- or area-specific ISN test.** "inter-regional coupling can generate a paradoxical response in a region that is not locally inhibition-stabilized or abolish it in a region that is. Thus, in an interconnected network, a paradoxical response is neither necessary nor sufficient for local inhibition stabilization." — [Wu & Miller 2026, bioRxiv doi:10.64898/2026.09.09.750483 [preprint]](https://doi.org/10.64898/2026.09.09.750483)
- **Primate laminar inhibitory control (marmoset).** Laminar recordings with bidirectional PV optogenetics show that "in granular layers, PV+ cells implement divisive/multiplicative linear gain control, whereas in extra-granular layers they exert tuned nonlinear suppression that enhances orientation tuning". In their modelling, PV controls gain by changing spike threshold and controls tuning by changing input noise, which reshapes the I/O function. — [Vafaei, Clark, Khadir, Dai, Federer & Angelucci 2025, bioRxiv doi:10.64898/2025.12.23.696300 [preprint]](https://www.biorxiv.org/content/10.64898/2025.12.23.696300.full.pdf) (also posted on Research Square, doi:10.21203/rs.3.rs-9117440/v1)
- **Macaque laminar variability.**
  - "In all layers, stimulation of the receptive field (RF) reduced both individual and shared variability relative to pre-stimulus baseline."
  - Near-surround stimulation further decreased variability in infragranular layers but "had little effect in granular and supragranular layers".
  - Far-surround stimuli increased variability in supragranular layers.
  - "surround suppression and variability are not governed by a single mechanism".
  — [Nurminen, Bijanzadeh & Angelucci 2026, eLife, doi:10.7554/eLife.86334](https://doi.org/10.7554/eLife.86334)
- **L4 input fraction (cat/rodent)**, used to argue for strong recurrent amplification:
  - Cooling cortex reduced L4 first-harmonic responses to about 1/3 of the intact value.
  - Transient L4 voltage responses were about 1/2 of intact with cortex silenced (Chung & Ferster 1998).
  - Stimulus-evoked excitatory conductance in L4 with cortex silenced was 30–40% of intact (Lien & Scanziani 2013; Li et al. 2013).
  — [Ahmadian & Miller 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9129095/)
- In L4, "recurrent connections strongly amplify feed-forward inputs so that most of the stimulus-evoked depolarization in layer 4 neurons originates intracortically" (Douglas et al. 1995 proposal, as summarized by Ozeki). — [Ozeki et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691725/)

### Inferences
- **Which layer is more strongly normalized?** The models disagree:
  - Rubin, Obeid and Mossing place the SSN in L2/3 and treat L4 as input.
  - Eskikand's connectivity-based column puts ISN operation in L4/L5, not L2/3.
  - Mouse data show ISN behaviour in both superficial and deep layers (Sanzeni).
  - Marmoset PV data show a qualitative laminar split (granular = linear divisive gain; extragranular = nonlinear tuned suppression).
  The theory therefore does not firmly predict that one layer normalizes more. That fits the team's observation (4) of similar normalization timing and strength across depth.
- **(4) Similar timing across depth** is most naturally explained if (i) each layer has its own fast local E/I stabilization with similar inhibitory time constants (I time constants of a few ms are generic in the SSN), or (ii) normalization arises within L4 and is inherited by L2/3 within a few ms. A purely L2/3-located SSN with linear L4 input (the Rubin/Obeid/Mossing assumption) would predict *later or stronger* normalization superficially than in L4, which the data do not show.
- **Superficial sites "start less linear" (4).** In the SSN this would mean superficial units sit at higher effective gain at onset. Possible reasons: higher spontaneous or early rates, already-summed binocular input from L4, or a larger share of binocular neurons, for which the weak-coupling summation is supralinear $(a+b)^n$. "Less linear" therefore needs a direction: supralinear early summation superficially would favour the SSN, while sublinear early summation would favour inheritance of an already normalized L4 signal.
- **Wu & Miller 2026 implies that interlaminar loops are "inter-regional" coupling.** A paradoxical response in one layer would not prove that the layer is locally ISN, so laminar ISN inferences from optogenetics need multi-layer modelling.
- **Nurminen 2026 cautions** against assuming a single SSN mechanism links suppression and variability quenching across layers. Laminar variability-quenching tests (Q6) should be read layer by layer.
- Per-paper lines for (1)–(6):
  - **Mossing 2021.** (1)–(2) via contrast-dependent recurrent amplification in L2/3. It predicts a qualitative change in which inhibitory loop dominates from low to high contrast (VIP-mediated facilitation at low contrast, SST-mediated suppression at high), so it allows supralinear facilitation at 5.5% (6). (4) predicts L2/3 transformation of L4 input, i.e. layer differences. (5) not addressed. Mouse, no binocular inputs.
  - **Eskikand 2022.** Predicts stronger ISN (and paradoxical) behaviour in L4/L5 than L2/3, which could match superficial sites being "different" at onset. No contrast or binocular dynamics. Neural-mass, not SSN.
  - **Obeid & Miller 2025.** L2/3 only. Contrast-dependent suppression via E→I spread. Contrast-independent offset decay. (5) handled by a Naka–Rushton input front end.

### Gaps
- I found no SSN or ISN model with explicit L4→L2/3→L5/6 (and feedback) populations that makes laminar *timing* predictions, and none with binocular input across layers. My searches covered Europe PMC with "supralinear/inhibition-stabilized" + "laminar/layer 4/layer 2/3", Miller/Ahmadian 2021–2026 publication records, and web search.
- Mossing et al. 2021 was found only as a preprint. Whether it has been published in a journal was not determined.
- No primate laminar optogenetic test of the paradoxical effect was found (see Q5).

---

## Q4. Binocular and interocular SSN models

### Takeaway
There is one dedicated binocular SSN: Selby & Tripp 2017 (from Selby's 2016 MASc thesis). It adds ocular dominance columns and reproduces interocular transfer of surround suppression in E but not I units. No SSN model of binocular rivalry, dichoptic contrast summation, or a dichoptic contrast matrix was found. A 2026 mouse preprint (Perumal et al.) directly tests ISN predictions for binocular summation by cell type: PV summation is ISN-like (sublinear, tracks pyramidal cells), while SST summation is linear or supralinear.

### Cited Findings
- **Selby & Tripp 2017.** "Extending the Stabilized Supralinear Network model for binocular image processing", Neural Networks 90:29–41. The authors built a binocular SSN incorporating ocular dominance columns, and it "could reproduce data on interocular transfer of surround suppression"; the paper describes itself as preliminary. The DOI was not retrieved (ScienceDirect PII S0893608017300552). **[paywalled]** — [ScienceDirect record](https://www.sciencedirect.com/science/article/abs/pii/S0893608017300552); bibliographic record via [DBLP](https://dblp.org/pid/155/6705)
- **Selby 2016 MASc thesis** (Univ. Waterloo, supervisor B. Tripp), "Development of an Integrated Model of Primary Visual Cortex": "The SSN was also expanded to use binocular stimuli. Using an optimization procedure, SSN parameters were found that produce interocular transfer of suppression in excitatory units, but not inhibitory ones. The lack of interocular transfer in inhibitory units may indicate that an alternate inhibition-stabilization scheme is more biophysically realistic." — [UWSpace hdl:10012/10264](http://hdl.handle.net/10012/10264)
- **Perumal et al. 2026 [preprint]** (mouse binocular V1). "A key prediction from ISNs is that inhibition tracks and scales with excitation, leading to sublinear response summation." "Binocular responses in PV interneurons tracked pyramidal neuron responses, exhibiting sublinear summation consistent with ISN models. In contrast, binocular responses [in] SST interneurons exhibited rapid temporal dynamics and linear or supralinear summation, deviating from ISN predictions." Both classes received callosal excitatory drive. Simulations implicate "class-specific differences in local inhibitory circuit motifs" rather than cellular or synaptic differences alone. — [Perumal, Bushnell, Dopp, Gharaei, Arabzadeh, Nair & Stuart 2026, bioRxiv doi:10.64898/2026.09.20.753031](https://doi.org/10.64898/2026.09.20.753031)
- **Non-SSN normalization model of rivalry and cross-orientation suppression**, given here for contrast with the SSN. — [Said & Heeger 2013, PLoS Comput Biol (PMC3610603)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3610603/)
- **Bistability in the SSN**, a possible substrate for rivalry; no rivalry application found. — [Kraynyukova & Tchumatchenko 2018](https://doi.org/10.1073/pnas.1700080115)
- **Rubin 2015's two-stimulus normalization**: about 0.5–0.7 of the sum at equal strength, approaching winner-take-all as strengths diverge, with linear average summation when one stimulus is below its half-maximal contrast (MT data). This is the template for two eyes treated as two inputs. — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)

### Inferences
- **Mapping Rubin's two-stimulus results to the 6×6 matrix.** Equal-contrast dichoptic pairs should give binocular/(L+R) ratios that fall with contrast toward ~0.5–0.7. Unequal pairs should approach the response to the stronger eye ("winner-take-all" interocular dominance). Pairs with one eye below its c50 should sum about linearly. This matches (2), and points to an analysis: plot the summation ratio against the *difference* and the *sum* of the eyes' contrasts separately.
- **Testing Selby–Tripp in the units.** Their finding (interocular suppression in E but not I units) and Perumal's PV/SST split both point to a testable prediction: interocular normalization in narrow-spiking (putative PV) units should exist and track broad-spiking units if the circuit is ISN-like. A narrow-spiking population that sums linearly would argue for a non-ISN (feedforward-pool) mechanism for that class.
- **Per-paper lines for (1)–(6):**
  - **Selby & Tripp 2017.** Static (steady-state) model, so silent on (1)–(3). OD columns provide the structured $W^{-1}$ cross-eye terms that can produce interocular normalization without supersaturation (5), but this has not been shown. (4) and (6) are not addressed.
  - **Perumal 2026.** (1)–(2): "temporal dynamics" of PV vs SST binocular summation. This is the closest empirical analogue to the team's time-resolved interocular normalization, in mouse. The mapping to laminar MUA is unknown.
- In an OD-column SSN, interocular normalization at the neuron level can arise from the *linear* balanced solution $-cW^{-1}g$ through negative cross-column elements of $W^{-1}$. That would explain strong interocular normalization (2) without supersaturation (5), unlike the 2-population mechanism (Q1). This is my inference from Ahmadian 2013's balanced solution; it is not a published result.

### Gaps
- No SSN or ISN model of binocular rivalry was found.
- No SSN model of dichoptic contrast summation (a contrast-pair matrix) in primate was found.
- The Selby & Tripp 2017 abstract, methods, and DOI could not be retrieved (HTTP 403), so it is unknown whether they modelled dynamics or contrast dependence. Obtain the PDF.
- No SSN-based treatment of interocular timing (onset asynchrony, monocular vs binocular latency) was found.

---

## Q5. Empirical tests of SSN/ISN in primate V1 and in layered recordings

### Takeaway
- The paradoxical-effect test of the ISN is established in cat (Ozeki 2009; anesthetized, intracellular) and mouse (Sanzeni 2020; several areas, superficial and deep layers).
- In macaque V1, no direct paradoxical-inhibition test was found. The available evidence is indirect: optogenetically evoked normalization (Nassi 2015), and "reshuffling" analyses implying strong coupling (Sanzeni et al. 2023). mDlx-ChR2 work in macaque demonstrates the tool (De et al. 2020) but does not test paradoxical responses.
- Laminar primate optogenetics is just appearing in marmoset (Vafaei et al. 2025 [preprint]).

### Cited Findings
- **Cat V1 (Ozeki 2009).** "surround stimuli decrease both the inhibition and excitation received by suppressed cells". These decreases "cannot be accounted for by suppression of input from the LGN", and "an ISN is the only solution to this apparent paradox". Two further ISN predictions were confirmed: (i) a transient increase in inhibition (30–50 ms) when the surround is added, and (ii) the source of synaptic input to suppressed cells. — [Ozeki et al. 2009](https://doi.org/10.1016/j.neuron.2009.03.028)
- **Mouse (Sanzeni 2020).** ISN operation was found "even without sensory stimulation" in V1, S1 and M1, "both in superficial and deep layers", using VGAT-ChR2 stimulation of all inhibitory cells. With PV-only stimulation, differences between viral and transgenic expression were explained by the fraction of PV cells stimulated. — [Sanzeni et al. 2020](https://doi.org/10.7554/eLife.54875)
- **Mouse V1 currents.** "in mouse V1, one study (Adesnik, 2017) found that excitatory and inhibitory currents scale as predicted in the SSN with increasing visual stimulation." — [Sanzeni et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7324160/)
- **Earlier mouse L2/3 results suggesting non-ISN behaviour, as discussed by Rubin et al.** In one 2012 study, PV suppression increased E responses with a non-paradoxical increase in inhibitory conductance, "suggesting a dearth of E→E coupling and non-ISN behavior". Another 2012 study found SOM cells surround-facilitated. The author names were cut off in the extracted text; these are probably Atallah et al. 2012 and Adesnik et al. 2012, so verify before citing. Rubin's "more robust prediction" is that "response suppression in E cells should be accompanied by a decrease in the I conductance they receive." — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/)
- **Alert macaque V1 (Nassi 2015).** "Optogenetic depolarization of excitatory neurons either facilitated or suppressed baseline activity, consistent with indirect recruitment of inhibitory networks. As predicted by the normalization model, neurons exhibited sub-additive responses to optogenetic and visual stimulation, which depended lawfully on stimulation intensity and luminance contrast." — [Nassi, Avery, Cetin, Roe & Reynolds 2015, Neuron 86:1504–1517, doi:10.1016/j.neuron.2015.05.040 (PMC4534089)](https://doi.org/10.1016/j.neuron.2015.05.040)
- **Mouse and monkey reshuffling (Sanzeni et al. 2023).** "optogenetic stimulation of excitatory neurons strongly modulated the activity of single neurons yet had weak or no effects on the distribution of firing rates across the population". "Key statistics of mouse and monkey responses lay on a continuum, with mice/monkeys occupying the low-/high-rate regions". Reshuffling "emerges generically in randomly connected excitatory/inhibitory networks, provided the coupling strength ... is sufficient that powerful inhibitory feedback cancels the mean optogenetic input". — [Sanzeni, Palmigiano, Nguyen, Luo, Nassi, Reynolds, Histed, Miller & Brunel 2023, Neuron, doi:10.1016/j.neuron.2023.09.018 [possibly paywalled]](https://doi.org/10.1016/j.neuron.2023.09.018)
- **Macaque V1 mDlx-ChR2 (De et al. 2020).** Units excited at short latency (<5 ms) "presumably expressed ChR2 and suppressed other neurons via synaptic inhibition". Expression in one monkey was concentrated superficially, and L4 "is difficult to transduce efficiently irrespective of the vector". No paradoxical-response analysis was found in the text. — [De, El-Shamayleh & Horwitz 2020, eLife 9:e52658](https://doi.org/10.7554/eLife.52658)
- **Marmoset V1 laminar PV manipulation.** Granular layers show linear divisive/multiplicative gain control; extragranular layers show tuned nonlinear suppression. — [Vafaei et al. 2025 [preprint]](https://www.biorxiv.org/content/10.64898/2025.12.23.696300.full.pdf)
- **Macaque V1 laminar variability quenching** is layer-specific with stimulus size. — [Nurminen et al. 2026](https://doi.org/10.7554/eLife.86334)
- **Contrast-dependent gamma in macaque V1.** The gamma peak frequency increases with grating contrast, and Holt et al. reproduce this from published macaque data and replicate a local-Gabor experiment (σ = 0.5°). — [Holt et al. 2024](https://doi.org/10.1371/journal.pcbi.1012190)
- **Monkey onset transients.** These are stimulus-tuned, as predicted by the sampling-SSN (re-analysis of published awake-monkey V1 data). — [Echeveste et al. 2020](https://doi.org/10.1038/s41593-020-0671-1)

### Inferences
- Primate support for the SSN/ISN is currently indirect: optogenetic normalization, reshuffling implying strong coupling at higher monkey rates, contrast-dependent gamma, and tuned onset transients. The strongest direct test, paradoxical I-cell responses or decreases in I conductance during suppression, has not been reported in macaque V1 in anything I found.
- Monkeys sitting in the "high-rate" region of the reshuffling continuum (Sanzeni 2023) fits the SSN expectation that awake macaque V1 under grating stimulation is strongly coupled, and so in the regime where (2) is expected.
- The laminar primate data (Vafaei; Nurminen) suggest the "normalization" signature may arise from different mechanisms in granular and extragranular layers even when its timing and strength look similar (4). That is a reason to check *how* normalization manifests by layer (divisive gain vs tuned suppression; variability quenching), not only its time course.

### Gaps
- No macaque V1 study testing the paradoxical inhibitory response (ISN) with cell-type-specific optogenetics was found, nor any layer-resolved ISN test in primate.
- The Holt et al. macaque gamma data source ("reference 1") was not identified by name in the retrieved text.

---

## Q6. Which time-resolved laminar dichoptic measurements would discriminate SSN dynamics from feedforward delayed normalization?

### Takeaway
The most discriminating SSN signatures are those that follow from an expansive nonlinearity plus inhibitory stabilization, not from divisive normalization per se:
(a) supralinear interocular summation at low (near-threshold) contrast that turns sublinear with contrast;
(b) a normalization latency or time constant that shrinks with *total* drive;
(c) contrast- and drive-dependent onset overshoot and gamma ringing frequency;
(d) sublinear summation in putative inhibitory units as well as excitatory ones, with I growing faster than E;
(e) contrast-independent offset decay at the cellular τ;
(f) variability quenching that co-occurs with sublinearity on a ~τ timescale.

The team's (2) and (3) already favour an input-strength-dependent recurrent mechanism over a fixed delay. They do not yet separate the SSN from recurrent shunting normalization (Carandini–Heeger, ORGaNICs).

### Cited Findings
- Supralinear→sublinear transition with input strength; both E and I sublinear; winner-take-all at unequal strengths. — [Rubin et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344127/); [Ahmadian et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026108/)
- Effective connectivity strengthens with drive, giving faster timescales and quenched variability; variability timescales ≈ single-neuron τ; offset decay at the cellular τ. — [Hennequin et al. 2018](https://doi.org/10.1016/j.neuron.2018.04.017)
- Gamma peak frequency rises with contrast via gain; spontaneous state is non-oscillatory. — [Holt et al. 2024](https://doi.org/10.1371/journal.pcbi.1012190)
- Onset transient magnitude scales with contrast and is inhibition-dominated. — [Echeveste et al. 2020](https://doi.org/10.1038/s41593-020-0671-1)
- ISN-specific transient increase in inhibition (30–50 ms) before steady-state suppression. — [Ozeki et al. 2009](https://doi.org/10.1016/j.neuron.2009.03.028)
- Offset decay set by the cellular τ, independent of stimulus size or contrast, in the L2/3 SSN. — [Obeid & Miller 2025](https://doi.org/10.1523/ENEURO.0459-24.2025)
- Sub-additivity and variability quenching co-occur across model classes. — [Goris et al. 2024 [paywalled]](https://doi.org/10.1038/s41583-024-00795-0)
- PV (ISN-like, sublinear) vs SST (linear or supralinear, fast) binocular summation in mouse. — [Perumal et al. 2026 [preprint]](https://doi.org/10.64898/2026.09.20.753031)
- Alternative models that also produce contrast-dependent dynamics:
  - shunting/recurrent normalization — [Carandini & Heeger 1994 [prior-knowledge citation]](https://doi.org/10.1126/science.8191289); [Carandini, Heeger & Movshon 1997 [prior-knowledge citation]](https://doi.org/10.1523/JNEUROSCI.17-21-08621.1997); [Heeger & Zemlianova 2020 [prior-knowledge citation]](https://doi.org/10.1073/pnas.2005417117)
  - delayed gain control with a low-pass normalization pool — [Zhou et al. 2019 [prior-knowledge citation]](https://doi.org/10.1371/journal.pcbi.1007484)

### Inferences
**Measurements, each with the SSN prediction and the prediction of a feedforward delayed normalization (FDN) model with a fixed pool delay:**
1. **Normalization onset latency and time constant vs drive.**
   - How: fit, per site and condition, the time at which the binocular/(L+R) ratio departs from its early value, and the time constant of approach to the sustained ratio.
   - SSN: both shrink monotonically with total drive. Test with $f(c_L)+f(c_R)$, or with the stronger eye's contrast, to see which governs. The relation should be continuous and close to a function of effective gain.
   - FDN: a fixed delay (pool latency plus filter τ); only the *depth* of normalization changes with contrast.
   - The team's (3) already rejects the strict FDN. The next step is to check whether the latency scales like an eigenvalue of an E/I loop (SSN-like) or like a conductance-shortened membrane τ (shunting-like). The second scales roughly as $\tau_m/(1+g_{shunt}/g_{leak})$ (standard membrane relation; my simplification, not taken from a specific paper).
2. **Low-contrast summation sign.**
   - How: compute summation at (5.5, 5.5), (5.5, 11) and (11, 11) in the sustained window and in the early window.
   - SSN: ≥ linear (supralinear for binocular units) at near-threshold drive, then a continuous transition to sublinear.
   - FDN / pure divisive: ≤ linear everywhere, except for threshold effects.
   - Caveat (6): a spiking threshold alone can produce apparent supralinearity, so compare against a feedforward model with a threshold.
3. **Onset overshoot and ringing.**
   - How: measure the peak-to-sustained ratio and the frequency of damped oscillations (MUA, CSD sinks, LFP gamma) for monocular vs binocular stimulation at matched per-eye contrast.
   - SSN: overshoot magnitude and gamma frequency increase with total drive, so binocular > monocular at equal per-eye contrast, with no gamma at 5.5%.
   - FDN: an overshoot (from the delayed divisor) whose timing is fixed, and no intrinsic resonance.
4. **Cell-class summation.**
   - How: split units into narrow- and broad-spiking.
   - SSN: both classes normalize interocularly, and narrow-spiking units grow faster with contrast (Hennequin; Rubin).
   - FDN: the pool (inhibitory) signal sums the eyes roughly linearly and is not itself normalized.
   - Perumal 2026 shows a class split in mouse; test for its primate analogue.
5. **Offset dynamics.**
   - How: measure the decay τ after the 24/26 offset codes as a function of contrast and binocularity.
   - SSN: decay at the cellular τ (~10–20 ms), independent of contrast (Obeid & Miller; Hennequin). Possibly an inhibitory-rebound undershoot proportional to drive.
   - Shunting normalization: contrast-dependent decay, faster after high contrast if the conductance persists.
   - FDN: the decay shape depends on the pool filter.
   - This is a cheap test with existing data.
6. **Variability quenching.**
   - How: time-resolved Fano factor or shared variability per site across the 6×6 matrix.
   - SSN: quenching scales with total drive and follows a ~τ time course matched to the onset of sublinearity.
   - FDN: no necessary link between normalization depth and variability.
   - Interpret layer by layer (Nurminen 2026).
7. **Laminar ordering.**
   - How: compare onset latency (CSD sink), normalization latency, and early-window summation sign across L4C, L2/3 and L5/6.
   - L4-located local SSN: normalization emerges in L4 as early as in L2/3, with L4 also showing a supralinear→sublinear transition.
   - L2/3-located SSN with linear L4 input (Rubin/Obeid/Mossing assumption): L4 stays near-linear while L2/3 normalizes.
   - Inherited feedforward normalization: L2/3 normalization latency equals L4's plus the conduction delay, with no extra contrast-dependent speed-up superficially.
   - The team's (4) (similar timing across depth) argues against the second option.
8. **Supersaturation constraint (5).**
   - Fit a 2-population SSN (or an OD-column SSN) with a Naka–Rushton front end to the sustained 6×6 matrix, enforcing $\Omega_E > 0$ (no supersaturation within 90%). Then check whether it can still produce the observed interocular collapse and its speed.
   - If it can do so only with $\Omega_E < 0$ and a peak just beyond 90%, the absence of supersaturation is a fine-tuning constraint. If an OD-column $W$ yields the collapse with $\Omega_E > 0$, the structured-SSN account is preferred.

**How each model or paper would account for (1)–(6) (inferences):**

| Model / paper | (1) early linear → norm. by 20–35 ms | (2) deeper and faster at high contrast | (3) fixed-delay FDN too slow at high c | (4) similar across depth; superficial less linear at start | (5) no supersaturation | (6) 5.5% near threshold |
|---|---|---|---|---|---|---|
| Generic SSN (Ahmadian 2013; Rubin 2015) | Yes: low initial gain gives feedforward-like start; recurrence builds within ~1–3 τ | Yes: gain ∝ r^{(n−1)/n} gives stronger, faster inhibition | Explained: delay is emergent and shrinks with drive | Single layer, so silent; superficial "less linear" would mean higher onset gain | Requires Ω_E > 0 (g_E/g_I > J_EI/J_II) or a peak beyond 90%; tension with "strongest sublinearity when Ω_E < 0" | Predicts supralinear or linear summation and slow dynamics |
| Hennequin 2018 (stochastic SSN) | Yes | Yes, plus variability quenching on a ~τ timescale | Yes | Silent | As generic | Weak quenching at 5.5% |
| Echeveste 2020 (sampling SSN) | Onset overshoot then settling | Overshoot and gamma frequency grow with contrast | Yes | Silent | Not addressed | Small or no overshoot |
| Holt 2024 (columnar SSN, gamma) | Not addressed (steady-state spectra) | Gamma frequency rises with drive | Implicitly | Columns horizontal, not laminar | Not addressed | No resonance near threshold |
| Ozeki 2009 (ISN) | Predicts a transient rise in I (30–50 ms) before suppression | Not contrast-resolved | n/a | Not laminar | n/a | Below ISN onset, non-paradoxical |
| Obeid & Miller 2025 (L2/3 SSN) | Not analyzed at onset | Suppression grows with contrast | Not analyzed | L2/3 only; L4 is input | Uses a saturating Naka–Rushton input; Ω_E < 0 would eventually supersaturate | Not addressed |
| Mossing 2021 (L4→L2/3 four-type SSN, mouse) [preprint] | Not analyzed in time | Positive VIP loop at low contrast, negative SST loop at high contrast | Not analyzed | Predicts an L2/3 transformation of L4 input (layer difference) | Not addressed | Facilitation possible at low contrast |
| Eskikand 2022 (laminar neural mass) [preprint] | Not addressed | Not addressed | n/a | L4/L5 ISN, L2/3 not, so predicts layer differences | n/a | n/a |
| Selby & Tripp 2017 (binocular SSN) | Static | Not reported | n/a | Not laminar | OD-column structure could give interocular normalization without supersaturation (untested) | n/a |
| Perumal 2026 (mouse bV1 PV/SST) [preprint] | PV dynamics track pyramidal cells; SST fast | PV sublinear (ISN-like) | n/a | Not laminar in the abstract | n/a | n/a |
| FDN with fixed pool delay (Zhou 2019-type) | Yes, by construction (delay) | Deeper but not faster | Fails (team's finding) | Can reproduce similar timing if the pool is shared | Saturation without supersaturation is easy | Sub-additive or threshold only |
| Recurrent shunting normalization (Carandini–Heeger 1994/1997; ORGaNICs 2020) | Yes | Yes: conductance shortens τ | Explained | Silent | Easy (divisive saturation) | No supralinearity unless a threshold is added |

### Gaps
- I found no published head-to-head comparison of SSN and delayed-normalization or shunting-normalization models on time-resolved data, laminar or otherwise. The discriminating predictions above are my synthesis and should be checked by simulating the candidate models on the 6×6 dichoptic design with 1-ms resolution.
- How "less linear at start" in superficial sites splits into supralinear vs sublinear is not specified in the context provided. That direction decides whether the SSN or an inherited-normalization account is favoured (Q3 Inferences).
- Whether MUA or CSD at each depth reflects E-dominated or mixed E/I output (needed for measurement 4) is a data question outside this literature search.
