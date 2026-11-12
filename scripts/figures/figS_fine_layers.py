"""Supplementary Figure S3: the compartments, and the laminar results with the finer, more error-prone labels (check AD; Blake, 2026-11-12:
coarse labels in the main figure, fine labels here).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_fine_layers.py
Main laminar set (tiers 1 + 2 of laminar alignment v4). Fine compartments, by height h above the layer 4C/5 anchor
scaled by the squish factor s: L2/3 h > 0.59 s; L4A/B (0.28 s, 0.59 s]; L4C (0, 0.28 s]; L5 [-0.22 s, 0);
L6 [-0.49 s, -0.22 s) (Zhu et al. 2026 thicknesses; at 100 um, two or three contacts per compartment below L2/3).
(a) monocular grating CSD averaged over monkey J's tier-1 penetrations, aligned to the 4C/5 anchor and scaled
    (figS1_methods.tier1_csd; Blake: examples from a clean J session or the tier-1 average of one monkey), with the
    coarse compartment boundaries;
(b) additivity index 30-100 ms after onset at 45% in each eye, measured and contrast summation, per penetration
    (docs/qc/ad_groups.csv), medians and interquartile ranges;
(c) the dominant eye's onset at 90% relative to L4C, per penetration (grey) and median; AD3.
Writes <out_root>/figures/figS_fine_layers.pdf and .png; prints the values shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from figS1_methods import tier1_csd  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402

warnings.filterwarnings("ignore")
QC = Path(__file__).resolve().parents[2] / "docs" / "qc"
FINE = ("L2/3", "L4A/B", "L4C", "L5", "L6")


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def main():
    fs.apply()
    G = pd.read_csv(QC / "ad_groups.csv")
    G = G[(G.scheme == "fine") & G.tier.isin((1, 2))]
    T = pd.read_csv(QC / "ad_tests.csv").set_index("test")
    fig = fs.figure("double", 62)
    fig.set_layout_engine("none")
    c0, ccb = axes_mm(fig, 20, 10, 40, 44), axes_mm(fig, 62, 10, 2, 44)
    info = tier1_csd(c0, ccb, "J")
    c0.set_title(f"monkey J, mean of {info['n']} tier-1 penetrations", fontsize=6, loc="left")
    a, b = axes_mm(fig, 90, 10, 34, 44), axes_mm(fig, 138, 10, 34, 44)
    rows = []
    for k, comp in enumerate(FINE):
        x = G[G.compartment == comp]
        for col, colr, dy, mfc in (("AI", fs.INK, -0.12, fs.INK), ("AI_sum", fs.MUTED, 0.12, "white")):
            v = x[col].dropna()
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            a.plot([q1, q3], [k + dy] * 2, color=colr, lw=0.9)
            a.plot([md], [k + dy], marker="o", ms=3.4, mfc=mfc, mec=colr, mew=0.9, zorder=3)
            rows.append(dict(compartment=comp, series=col, n=len(v), median=md))
    a.set_yticks(range(len(FINE)), FINE)
    a.set_ylim(len(FINE) - 0.4, -0.6)
    a.set_xlim(0.3, 0.95)
    a.tick_params(axis="y", length=0)
    a.set_xlabel("additivity index, 30–100 ms")
    a.set_title("measured (filled); summation (open)", fontsize=6, loc="left")
    on = G.pivot(index="penetration", columns="compartment", values="dom_90")
    rel = on.sub(on["L4C"], axis=0)
    for k, comp in enumerate(FINE):
        v = rel[comp].dropna()
        b.scatter(v, np.full(len(v), k) + np.linspace(-0.15, 0.15, len(v)), s=5, color=fs.MUTED, lw=0)
        b.plot([v.median()], [k], marker="|", ms=9, mew=1.6, color=fs.INK)
        rows.append(dict(compartment=comp, series="onset - L4C", n=len(v), median=v.median()))
    b.axvline(0, color=fs.RULE, lw=0.6)
    b.set_yticks(range(len(FINE)), FINE)
    b.set_ylim(len(FINE) - 0.4, -0.6)
    b.tick_params(axis="y", length=0)
    b.set_xlim(-12, 16)
    b.set_xlabel("onset relative to L4C (ms)\ndominant eye, 90%")
    t6 = T.loc["AD2 shortfall L6 - L4C [tiers 1+2]"]                      # stated in the caption
    print(f"AD2 L6 - L4C: {int(t6.n_in_direction)}/{int(t6.n)}, p = {t6.p:.3f}")
    for ax, letter in ((c0, "a"), (a, "b"), (b, "c")):
        fs.panel_label(ax, letter, x=-0.2)
    print(fs.save(fig, "figS_fine_layers"))
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
