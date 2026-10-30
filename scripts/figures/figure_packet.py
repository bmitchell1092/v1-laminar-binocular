"""Figure packet: one PDF with a cover page and each figure with its title and legend (for discussion).

Usage: .venv\\Scripts\\python scripts\\figures\\figure_packet.py
Run after the figure scripts. Legends below are the reviewed text; numbers in them were checked against the figure
scripts' printed statistics on the date in COVER. Writes <out_root>/figures/figure_packet.pdf.
"""

import pymupdf

from v1lb.config import out_root

A4 = pymupdf.paper_rect("a4")
MARGIN = 50
FONTS = {"arial": "C:/Windows/Fonts/arial.ttf", "arialbd": "C:/Windows/Fonts/arialbd.ttf"}   # has en dash, minus, <=
FONT, BOLD = "arial", "arialbd"


def new_page(doc):
    page = doc.new_page(width=A4.width, height=A4.height)
    for name, path in FONTS.items():
        page.insert_font(fontname=name, fontfile=path)
    return page

COVER = [
    ("title", "Binocular combination over time in macaque V1"),
    ("sub", "Draft figures for discussion · Blake Mitchell · 2026-10-29"),
    ("p", "In awake macaque V1 the two eyes' signals add almost linearly at response onset, and interocular "
          "normalization takes over 20–35 ms later. It collapses binocular summation further at higher contrast, "
          "at every depth, starting in the earliest synaptic input. A delayed normalization pool explains the "
          "direction of this sequence but not its full strength or its time course."),
    ("h", "Data"),
    ("p", "Two macaques (B: 8 penetrations, J: 18), laminar recordings in V1 (32 contacts, 100 µm), static gratings "
          "on the receptive field through a mirror stereoscope, each eye at 0, 5.5, 11, 22.5, 45 or 90% contrast "
          "(6 × 6 conditions, ~30 presentations each). Multi-unit activity (MUA), current-source density (CSD) and "
          "single units, on a photodiode-corrected clock. Presentations with saccades ≥ 1° or MUA artifacts are "
          "excluded. The penetration is the unit of replication throughout."),
    ("h", "Where each claim stands"),
    ("li", "Fig 2, linear at onset, normalized later: robust at 45%, significant at 11–45% (p ≤ .001)."),
    ("li", "Fig 3, stronger collapse at higher contrast: 26/26 penetrations, survives every control so far. "
           "Pre-registered held-out test failed as run (6/7, p = .06); 7/8 (p = .035) after artifact rejection, "
           "a sensitivity analysis. Not demonstrably earlier at higher contrast."),
    ("li", "Fig 4, same sequence at every depth and in the earliest current sink: consistent, modest statistics."),
    ("li", "Fig 5, delayed pool with w > a: robust (23/25, 24/25). Two-stage gain ~1.5%, consistent. "
           "No fitted model captures the late collapse."),
    ("li", "Fig 6, the fitted model rises too fast at low contrast, and the offset transient is unexplained. "
           "A recurrent-normalization test is planned, not run."),
    ("h", "Not yet done"),
    ("p", "Vergence, RF-position and eye-crosstalk controls (Supplementary S1); the recurrent-model test; "
          "dual-probe (inter-column) timing; manuscript text."),
]

FIGURES = [
    ("fig1_setup", "Figure 1. Binocular combination measured with millisecond and laminar precision",
     "(a) Mirror stereoscope: each eye viewed its own monitor. (b) Gratings matched to the receptive field were shown "
     "at six contrasts to each eye independently (36 conditions, ~30 presentations each; three per fixation trial, "
     "317 ms on, 317 ms off). (c) Flash-evoked CSD in an example penetration (221212_J_eA). Depth is measured from "
     "the top of cortex (onset of spontaneous MUA). The L4/5 boundary lies just below the dominant early flash sink. "
     "The flash recording is registered to the grating recording. (d) Photodiode-measured stimulus onset differs "
     "from the event code by whole 120 Hz frames; every analysis uses photodiode onsets. (e) Example penetration, "
     "mean of 28 driven channels, each scaled to its dominant-eye peak at 90%. Shown: each eye alone at 45%, both "
     "eyes at 45%, and the dominant eye at 90%, which approximates linear summation of contrast ((1 + a)·45%, eye "
     "weight a ≈ 0.7). (f) Monocular contrast responses (40–250 ms) of the dominant and other eye: mean and 95% "
     "bootstrap CI over 26 penetrations. (g) Ocular dominance index at 90% for driven MUA channels and isolated "
     "single units."),
    ("fig2_add_then_normalize", "Figure 2. The two eyes add at response onset, then normalize within tens of "
     "milliseconds",
     "(a) Population responses at 45%, aligned to each channel's response onset (half-maximum of the dominant-eye "
     "response at 90%): mean and 95% CI over 26 penetrations. For ~20 ms (shaded) the binocular response follows the "
     "linear-summation reference; it then falls to the dominant eye's level. (b) Equivalent contrast. The binocular "
     "response (diamond) is placed on the dominant eye's contrast response measured in the same window; k is the "
     "equivalent contrast divided by 45%. Early, k sits near the linear-summation prediction (1 + a)·45% (triangle); "
     "late, it sits near 45%. Population means, as an illustration. Because k compares responses at the same "
     "time, it is not biased by the earlier onset of binocular responses. (c) Linearity index L = (k − 1)/a in "
     "onset-aligned windows (L = 1, linear summation; L = 0, binocular equals the dominant eye alone). Median over "
     "penetrations of channel medians, 95% bootstrap CI. (d) Early (0–20 ms after onset) vs late (50–200 ms) L per "
     "penetration (grey lines; a few clipped at the axis). Wilcoxon signed-rank tests, 26 penetrations."),
    ("fig3_contrast", "Figure 3. Contrast sets how far binocular summation collapses",
     "(a) Raw binocular gain, with no model or inversion: both eyes minus the dominant eye alone at the same "
     "contrast, onset-aligned. The early gain is larger, and is lost sooner, at higher contrast. (b) L 20–50 ms "
     "after onset per penetration (grey) and the median: L falls with contrast in every penetration. (c) Slope of L "
     "on log contrast for the discovery penetrations (monkeys B and J) and for held-out cross-recordings (each "
     "probe's recording during the other probe's grating run on dual-probe days; monkey J), with one-sided sign "
     "tests. The pre-registered held-out test failed as run (6/7, p = .06). After MUA artifact rejection, a "
     "cleaning rule adopted later for an unrelated reason, the same test gives the counts shown. This is reported as "
     "a sensitivity analysis, not a confirmation. (d) Amplitude-normalized half-transition time: the collapse is "
     "deeper at higher contrast, but not reliably earlier."),
    ("fig4_depth_input", "Figure 4. The same sequence runs at every depth and starts in the earliest synaptic input",
     "(a) Grating-evoked CSD in the example penetration, with the earliest current sink (arrow) and the L4/5 "
     "boundary. (b) L by depth below the top of cortex, early and late (11% and 22.5% pooled): the collapse occurs "
     "at all depths. (c) Equivalent contrast k of the earliest sink (CSD) and of MUA on the same channels, timed from "
     "sink onset. The input current also adds near-linearly at onset and then normalizes. (d–f) Channels above, at "
     "(± 2 contacts) and below the earliest sink. (d) The sink channels respond first. (e) Each depth normalizes a "
     "similar time after its own onset (first 15 ms window with L < 0.5). (f) On the stimulus clock, normalization "
     "therefore arrives first at the input layer. Wilcoxon tests against the sink group, uncorrected."),
    ("fig5_models", "Figure 5. A delayed binocular normalization pool explains the direction, but not the full "
     "strength, of the collapse",
     "(a) Single-stage delayed normalization (after Zhou et al. 2019). The drive E sums the eyes, with weight a for "
     "the other eye; the pool sums them with its own weight w and acts through a low-pass delay τ2. The two-stage "
     "variant first divides each eye by both eyes' contrasts (Meese, Georgeson & Baker 2006, made dynamic). Fitted "
     "per penetration to the full 6 × 6 time courses. (b) Delayed normalization predicts held-out presentations "
     "(odd/even) better than static normalization. (c) The pool weights the other eye more than excitation does "
     "(w > a). (d) The two-stage model predicts held-out responses consistently but only slightly better (~1.5%), "
     "under two loss functions; both tests were pre-registered at the channel level. (e) The pre-registered capture "
     "test: late L of each model's fit minus the data. A model would capture the collapse only if every 90% CI lay "
     "within the ±0.25 band; none does, and the two-stage model stays too linear late."),
    ("fig6_dynamics", "Figure 6. The fitted model rises too fast at low contrast, and responses carry an offset "
     "transient it cannot produce",
     "(a) Dominant-eye responses, each penetration scaled to its peak: the rise steepens with contrast. (b) 10–90% "
     "rise time: data (penetrations grey; median and 95% CI) against the fitted delayed-normalization model of Fig 5, "
     "which was fitted to 0–219 ms, so this is out of sample. The model rises about 2.3× too fast at low contrast. "
     "(c) After the response offset, the data show an offset transient (about 1.8× the sustained level, in 88–100% "
     "of penetrations) that the model cannot produce. Decay statistics taken over this period mostly measure the "
     "transient, so they are not used to compare mechanisms. Whether recurrent or shunting normalization explains the "
     "slow, contrast-dependent rise is the subject of a planned, pre-registered model test."),
]


def cover(doc):
    page = new_page(doc)
    y = MARGIN + 40
    w = A4.width - 2 * MARGIN
    for kind, text in COVER:
        size, font, gap, indent = {"title": (20, BOLD, 10, 0), "sub": (10, FONT, 24, 0), "h": (11, BOLD, 6, 0),
                                   "p": (9.5, FONT, 14, 0), "li": (9.5, FONT, 6, 12)}[kind]
        if kind == "li":
            page.insert_text((MARGIN + 2, y + size), "•", fontsize=size, fontname=FONT)
        rect = pymupdf.Rect(MARGIN + indent, y, MARGIN + w, A4.height - MARGIN)
        used = page.insert_textbox(rect, text, fontsize=size, fontname=font, lineheight=1.3)
        y = A4.height - MARGIN - used + gap


def figure_page(doc, name, title, legend):
    src = pymupdf.open(out_root() / "figures" / f"{name}.pdf")
    r = src[0].rect
    w = A4.width - 2 * MARGIN
    fw = min(r.width, w)                                   # native size: figures are drawn at journal widths
    fh = r.height * fw / r.width
    x0 = MARGIN + (w - fw) / 2
    page = new_page(doc)
    used = page.insert_textbox(pymupdf.Rect(MARGIN, MARGIN, MARGIN + w, MARGIN + 60), title, fontsize=11,
                               fontname=BOLD, lineheight=1.25)
    top = MARGIN + 60 - used + 10
    h = fh
    page.show_pdf_page(pymupdf.Rect(x0, top, x0 + fw, top + fh), src, 0)
    page.insert_textbox(pymupdf.Rect(MARGIN, top + h + 14, MARGIN + w, A4.height - MARGIN), legend, fontsize=8.5,
                        fontname=FONT, lineheight=1.3)


def main():
    doc = pymupdf.open()
    cover(doc)
    for name, title, legend in FIGURES:
        figure_page(doc, name, title, legend)
    out = out_root() / "figures" / "figure_packet.pdf"
    out.unlink(missing_ok=True)
    doc.save(out, garbage=4, deflate=True)
    print(out)


if __name__ == "__main__":
    main()
