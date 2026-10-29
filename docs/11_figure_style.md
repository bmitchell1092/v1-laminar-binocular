# 11 — Figure style (house style for the paper figures)

All paper figures are built by `scripts/figures/*.py` on top of `v1lb.figstyle`. Each script starts with
`fs.apply()`, takes its colours from the module, and ends with `fs.save(fig, name)`. Change a convention here and in
`figstyle.py`, never in one script.

## Build

| Step | Command | Output |
|---|---|---|
| shared inputs | `scripts/figures/figure_data.py` | `<out_root>/figures/data/population.npz` (per-penetration 6×6 time courses, flash- and onset-aligned) |
| figures | `scripts/figures/fig1_setup.py` … `fig6_dynamics.py` | `<out_root>/figures/<name>.pdf` (vector, editable text) and `.png` (300 dpi) |

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
  triangles, and the legend says so. No silent clipping.

## Colour (validated with the dataviz palette validator on a white surface)

| Role | Encoding | Values |
|---|---|---|
| Contrast (ordinal) | one-hue blue ramp, light → dark | 5.5% `#86b6ef`, 11% `#3987e5`, 22.5% `#256abf`, 45% `#184f95`, 90% `#0d366b` |
| Condition at one contrast | categorical | dominant eye blue `#2a78d6`, other eye orange `#eb6834`, both eyes ink `#0b0b0b` |
| Linear-summation reference | muted grey, dashed | `#8a8985` |
| Data vs models | categorical | data ink; fixed-delay (feedforward) models orange; recurrent models aqua `#1baf7a` (below 3:1 on white, so always labelled) |
| Model variants | marker shape and fill, same hue | single-stage filled circle, two-stage open square, Ding–Sperling open triangle |
| Monkey | marker shape, never colour | B circle, J square |
| Epoch (early vs late) | grey ramp | early `#8a8985`, late ink |
| CSD | diverging, grey midpoint | sink dark blue ↔ source dark red; colour bar labelled "sink / source" |
| Neutral | ink `#0b0b0b`, secondary `#52514e`, rules `#d9d8d3` | |

Rules: contrast always uses the blue ramp and nothing else uses it; the dominant/other eye pair keeps its colours in
every figure; a reference line for linear summation is always grey dashed at L = 1 (or k = 1 + a), and "binocular =
dominant eye alone" is the solid light rule at L = 0 (k = 1).

## Conventions of content

* **Time axes** say their zero: "time from stimulus onset" (photodiode), "time from response onset" (each channel's
  half-max onset, dominant eye at 90%), "time from sink onset", "time from response offset".
* **Contrast** is shown in percent: 5.5, 11, 22.5, 45, 90.
* **Responses** are MUA scaled per channel to the dominant-eye peak at 90% (40–250 ms) unless the axis says
  otherwise.
* **Depth** is "depth below top of cortex (mm)", increasing downward. Laminar landmarks are drawn as dashed lines
  labelled "top of cortex" and "L4/5".
* **Pre-registered tests** are marked as such on the panel, and a failed one is shown as plainly as a passed one
  (Fig 3c, Fig 5e).
* **Illustrations** (population means used to explain a measure, hand-set simulations) are labelled as such in
  the panel, not only in the legend.
