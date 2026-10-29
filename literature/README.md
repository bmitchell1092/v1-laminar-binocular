# literature/

Full-text papers for close reading (model equations, parameter tables, methods). PDFs are **not** committed:
`*.pdf` is in `.gitignore`, so publisher PDFs stay local. Only this README and any notes (`*.md`) are versioned.
Scanned PDFs (no text layer) are read by rendering pages with PyMuPDF (`pymupdf` in the project venv).

## On disk (added by Blake, 2026-10-29)

| File | Paper | Used for |
|---|---|---|
| `ding_sperling.pdf` | Ding & Sperling 2006, PNAS 103:1141, "A gain-control theory of binocular combination" (scan) | Eq. 6 verified: Γ = Σ_e I_e (1 + ε_e)/(1 + ε_L + ε_R); `v1lb.dynorm.predict_ds` |
| `Ding-A unified model for binocular fusion and depth percepti.pdf` | Ding & Levi 2021, Vision Res, "A unified model for binocular fusion and depth perception" | binocular combination model family |
| `74235376.pdf` | Baker, Meese & Georgeson 2007, "Binocular interaction: contrast matching and contrast discrimination are predicted by the same model" (accepted version) | two-stage model equations and parameters |
| `jov-6-11-7.pdf` | Meese, Georgeson & Baker 2006, J Vis 6(11):7, "Binocular contrast vision at and above threshold" | Table 1 Eq. 3 verified: Stage1(L) = L^m/(s + L + R), resp = (ΣStage1)^p/(z + (ΣStage1)^q); Table 2 two-stage fit m 1.28, s 0.985, p 7.99, q 6.59, z 0.077, σ 0.194 (settles the parameter conflict in the report; "z 0.0089" was wrong); twin-summation model best (RMS 0.664 vs 0.873 dB) |
| `Heeger-A recurrent circuit implements normalization, simulat.pdf` | Heeger & Zemlianova 2020, PNAS 117:22494, "A recurrent circuit implements normalization, simulating the dynamics of V1 activity" | ORGaNICs Eqs. 1–7 verified; recurrent-normalization model test |
| `nihms650360.pdf` | Rubin, Van Hooser & Miller 2015, Neuron, "The stabilized supralinear network" (author manuscript) | SSN |
| `Pattadkal_2024_Neuron.pdf` | Pattadkal, Zemelman, Fiete & Priebe 2024, Neuron 112:661, "Primate neocortex performs balanced sensory amplification" (scan) | amplification circuits distinguished by dynamics |

## Still useful if available

* Carandini, Heeger & Movshon 1997 (J Neurosci): conductance-based normalization dynamics.
* Albrecht et al. 2002 (J Neurophysiol): the "normalization is fast" counter-evidence.
