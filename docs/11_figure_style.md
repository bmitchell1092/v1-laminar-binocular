# 11 — Figure style (house style for the paper figures)

All paper figures are built by `scripts/figures/*.py` on top of `v1lb.figstyle`. Each script starts with
`fs.apply()`, takes its colours from the module, and ends with `fs.save(fig, name)`. Change a convention here and in
`figstyle.py`, never in one script.

## Build

| Step | Command | Output |
|---|---|---|
| shared inputs | `scripts/figures/figure_data.py` | `<out_root>/figures/data/population.npz` (per-penetration 6×6 time courses, flash- and onset-aligned) |
| figures | `scripts/figures/fig1_setup.py` … `fig4_model.py`; supplementary `figS1_*.py` … `figS6_*.py` (figure set: docs/09) | `<out_root>/figures/<name>.pdf` (vector, editable text) and `.png` (300 dpi) |

Each figure script reads only saved analysis outputs (`<out_root>/binocular/*.parquet`, `docs/qc/*.csv`) or the
shared cache, and prints the statistics it displays so they can be checked against docs/09.

## Format

* **Widths:** journal columns: single 85 mm, 1.5 column 114 mm, double 174 mm (`fs.figure(width, height_mm)`).
  Heights are chosen per figure; most are double column, 60–120 mm tall.
* **Type:** Arial (Helvetica fallback). Tick labels 6.5 pt, axis labels 7 pt, in-panel notes 5.5–6 pt, panel letters
  10 pt bold, lower case, top-left outside the axes (`fs.panel_label`).
* **Axes:** left and bottom spines only, 0.6 pt, dark grey; ticks outward, 2.5 pt. No grid. No titles that restate
  the axis; short in-panel notes say what is plotted (contrast, number of penetrations, statistic).
* **Lines:** data 1.0–1.2 pt (the binocular "both eyes" trace 1.2 pt); individual penetrations 0.4–0.5 pt light grey.
* **Uncertainty:** 95% bootstrap CI over penetrations (the unit of replication), as a band for time courses and as
  error bars for windowed points. Medians where the statistic is heavy-tailed (L, k), means for responses. The
  legend says which.
* **Keys:** no boxed legends. A compact line key (`fs.line_key`) or direct labels. Text is always in ink, never in
  the series colour, except direct labels placed beside their line.
* **Statistics on the panel:** counts as "k/n penetrations", tests as "p = 0.02" or "p < 0.001"
  (`fs.fmt_p`). The test is named in the legend.
* **Clipping:** when a few penetrations fall outside sensible limits, they are drawn at the axis edge as open
  triangles, and the legend says so and gives the count (the script prints it). No silent clipping, and no count
  text on the panel.
* **One reading per panel** (2026-11-02, after Fig 2c): a panel should be readable at first glance without the
  legend. Prefer one condition shown well to two crammed side by side (the rest go to the supplement); put at most
  one short statistic on the panel; leave fitted curves out of a main panel when they fit the late data poorly.
  Check every label against every line, error bar and key in the rendered PNG: no text may cross a mark.
* **Images:** colour maps (CSD) are drawn with `imshow(..., interpolation="bilinear")`. A Gouraud-shaded
  `pcolormesh` leaves saw-tooth seams between channel rows at 300 dpi.

## Colour (validated with the dataviz palette validator on a white surface)

| Role | Encoding | Values |
|---|---|---|
| Contrast (ordinal) | one-hue blue ramp, light → dark | 5.5% `#86b6ef`, 11% `#3987e5`, 22.5% `#256abf`, 45% `#184f95`, 90% `#0d366b` |
| Condition at one contrast | categorical | dominant eye blue `#2a78d6`, other eye orange `#eb6834`, both eyes ink `#0b0b0b` |
| Linear-summation reference | muted grey, dashed | `#8a8985` |
| Data vs models | categorical | data ink; static model dark grey `#52514e`, dashed; fixed-delay aqua `#1baf7a`; conductance violet `#6250d6`; recurrent magenta `#e87ba4`. No model uses orange, the other eye's colour (2026-11-01). Aqua and magenta are below 3:1 on white and 6.1 apart for deutans, so always direct-labelled. Inside a panel in the Fig 2a format, a model is a dashed line in the condition's colour; every data trace there is solid (the 90% one too) and the key reads "solid: data / dashed: fit" |
| Model variants | marker shape and fill, same hue | single-stage filled circle, two-stage open square, Ding–Sperling open triangle |
| Monkey | marker shape, never colour | B circle, J square |
| Epoch (early vs late) | grey ramp | early `#8a8985`, late ink |
| CSD | diverging, grey midpoint | sink dark blue ↔ source dark red; colour bar labelled "sink / source" |
| Neutral | ink `#0b0b0b`, secondary `#52514e`, rules `#d9d8d3` | |

Rules: contrast always uses the blue ramp and nothing else uses it; the dominant/other eye pair keeps its colours in
every figure; a reference line for linear summation is always grey dashed at L = 1 (or k = 1 + a), and "binocular =
dominant eye alone" is the solid light rule at L = 0 (k = 1).

## Conventions of content

* **No shorthand in figures or text.** Measures are named in words on every axis ("linearity index", not "L"); a
  symbol may appear only in a formula in the Methods (Blake, 2026-10-30).
* **Population readouts for noisy ratios.** Where a per-channel ratio is heavy-tailed (the linearity index), time
  courses are read out from population-average responses with penetration-bootstrap bands, and per-penetration
  statistics are shown separately (Fig 4c vs Fig 2d).
* **CSD shown as population averages** aligned on a landmark (the earliest sink), not single noisy examples.

* **Time axes** say their zero: "time from stimulus onset" (photodiode), "time from response onset" (each channel's
  half-max onset, dominant eye at 90%), "time from sink onset", "time from response offset".
* **Contrast** is shown in percent: 5.5, 11, 22.5, 45, 90.
* **Iso-response contours** are drawn in percent contrast with equal axes: "dominant-eye contrast needed (%)"
  against "other-eye contrast × eye weight (%)". Linear summation is the grey dashed diagonal; "dominant eye alone"
  is the light horizontal rule at the target contrast (no L-shaped rule).
* **Early vs late per penetration** (paired panels): tick labels give the windows ("0–20", "50–200") over "window
  after onset (ms)"; the condition ("both eyes 45%") and its test head each pair.
* **Responses** are MUA scaled per channel to the dominant-eye peak at 90% (40–250 ms) unless the axis says
  otherwise.
* **Depth** is "depth below top of cortex (mm)", increasing downward. Laminar landmarks are drawn as dashed lines
  labelled "top of cortex" and "L4/5".
* **No registration language in figures** (Blake, 2026-10-30: registration is internal). A test that failed is
  shown as plainly as one that passed.
* **No session IDs, and "condition", never "cell"** (2026-11-01).
* **Illustrations** (population means used to explain a measure, hand-set simulations) are labelled as such in
  the panel, not only in the legend.
