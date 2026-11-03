"""Supplementary Figure S7: the synaptic current at the input (docs/09 checks J, M and O; supplementary by Blake's
decision, 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\figures\\figS7_input_current.py
(a) current-source density at 45% in each eye, depth relative to the earliest sink x time from sink onset: both eyes,
    the sum of the two monocular currents, and both eyes - that sum signed so that grey = less current than the sum
    (shown where the summed current exceeds 8% of the sink's peak); each penetration scaled by its sink's peak
    current, then averaged (v4_csd_additivity.py);
(b) additivity index B / (D + O) of the sink current (the sink channel and its two neighbours) and of MUA on the same
    channels, in windows after sink onset, at 45%; per penetration, medians and interquartile ranges (check O);
(c) the sink current, both eyes - "if the contrasts add" (the current's own contrast response and eye weight), same
    windows (check O, reported);
(d) timing at the earliest sink (half-max of each signal's own peak): MUA - current (check J, three-contact current
    corrected for the 0.9 ms LFP filter delay and the simulated filter bias; and the single-contact, unsmoothed
    current), isolated units - current (single-contact current, units within 0.2 mm of the cross-validated sink;
    check M) and isolated units - MUA on their own contact, at the sink and away from it.
Inputs: <out_root>/figures/data/csd_additivity.npz; docs/qc/v4_O_penetrations.csv, v4_J.csv, v4_J_penetrations.csv,
v4_J_sensitivity.csv, v4_M_units.csv.
Writes <out_root>/figures/figS7_input_current.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 118
WINDOWS = ("early", "dip", "late")
WIN_LABEL = {"early": "0–20", "dip": "30–100", "late": "130–200"}
POS = {"a1": (16, 70, 34, 36), "a2": (54, 70, 34, 36), "acb1": (90, 70, 2, 36), "a3": (116, 70, 34, 36),
       "acb3": (152, 70, 2, 36), "b": (18, 12, 38, 38), "c": (72, 12, 30, 38), "d": (126, 12, 44, 38)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def maps(axs, caxes, Z, c=0.45):
    rel, off = Z["rel"], Z["offsets"] * 0.1
    B, S = np.nanmean(Z[f"both_{c:g}"], axis=0), np.nanmean(Z[f"sum_{c:g}"], axis=0)
    n = np.isfinite(Z[f"both_{c:g}"][:, :, 0]).sum(axis=0)
    ok = n >= 10
    D = np.where(np.abs(S) > 0.08, (B - S) * np.sign(S), np.nan)
    sink_cmap = fs.CSD_CMAP.reversed()                            # the current is drawn sink-positive
    specs = ((B, sink_cmap, 1.2, f"both eyes {fs.contrast_label(c)}"), (S, sink_cmap, 1.2, "dominant + other alone"),
             (D, fs.ADDED_CMAP, 0.4, "both − sum"))
    ims = []
    for k, (ax, (M, cmap, lim, title)) in enumerate(zip(axs, specs)):
        im = ax.imshow(M[ok], aspect="auto", origin="lower", cmap=cmap, vmin=-lim, vmax=lim, interpolation="bilinear",
                       extent=(rel[0], rel[-1], off[ok][0] - 0.05, off[ok][-1] + 0.05), rasterized=True)
        ax.set_ylim(off[ok][-1] + 0.05, off[ok][0] - 0.05)
        ax.set_xlim(-10, 150)
        ax.set_xticks([0, 50, 100, 150])
        ax.axvline(0, color=fs.INK, lw=0.6)
        ax.axhline(0, color=fs.INK, lw=0.5, ls=(0, (2, 2)))
        ax.set_title(title, fontsize=6.5)
        ax.set_xlabel("time from sink onset (ms)")
        if k:
            ax.tick_params(labelleft=False)
        ims.append(im)
    axs[0].set_ylabel("depth relative to\nearliest sink (mm)")
    cb = axs[0].figure.colorbar(ims[0], cax=caxes[0], ticks=[-1, 0, 1])
    cb.ax.set_yticklabels(["source", "0", "sink"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    cb = axs[0].figure.colorbar(ims[2], cax=caxes[1], ticks=[-0.4, -0.2, 0, 0.2, 0.4])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    cb.set_label("grey: less current\nthan the sum", fontsize=6)
    sel = (rel >= 0) & (rel < 20)
    return dict(both_over_sum_at_sink_0_20=float(np.nanmean(B[off == 0][:, sel]) / np.nanmean(S[off == 0][:, sel])))


def windows_panel(ax, Pn, series, ref, ylabel, ylim):
    rng = np.random.default_rng(4)
    A = Pn[Pn.split == "all"]
    rows = []
    for k, w in enumerate(WINDOWS):
        for lab, (col, colr, dx) in series.items():
            v = A[col.format(w=w)].dropna()
            x = k + dx
            ax.scatter(x + rng.uniform(-0.05, 0.05, v.size), v.clip(*ylim), s=3, color=fs.RULE, lw=0, zorder=1)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([x, x], [q1, q3], color=colr, lw=0.9, zorder=2)
            ax.plot([x], [md], marker="o", ms=3.2, color=colr, zorder=3)
            rows.append(dict(series=lab, window=w, n=len(v), median=md))
    ax.axhline(ref, color=fs.MUTED, lw=0.8, ls=fs.SUM_LS)
    ax.set_xticks(range(3), [WIN_LABEL[w] for w in WINDOWS])
    ax.set_xlabel("ms after sink onset")
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.5, 2.5)
    return pd.DataFrame(rows)


def timing_panel(ax):
    J = pd.read_csv(QC / "v4_J_penetrations.csv")
    Jt = pd.read_csv(QC / "v4_J.csv")
    bias = float(Jt.loc[Jt.test.str.startswith("filter bias"), "bias_onset"].median())
    Js = pd.read_csv(QC / "v4_J_sensitivity.csv")
    one = Js[(Js.condition == "both 45%") & (~Js.hamming) & (Js.span == "1 contact") & (Js.sigma == 2)]
    U = pd.read_csv(QC / "v4_M_units.csv")
    at = U[U.grp == "at"].groupby(["penetration", "cluster_id"]).mean(numeric_only=True)
    away = U[U.grp != "at"].groupby(["penetration", "cluster_id"]).mean(numeric_only=True)
    rows = [("MUA − current\n(three contacts)", J["d_onset_both eyes 45%"] - bias, fs.INK),
            ("MUA − current\n(one contact)", one.d_onset, fs.INK),
            ("unit − current\n(one contact)", at.unit_minus_current_half, fs.INK2),
            ("unit − MUA,\nat the sink", at.unit_minus_mua_half, fs.INK2),
            ("unit − MUA,\nelsewhere", away.unit_minus_mua_half, fs.INK2)]
    rng = np.random.default_rng(5)
    out = []
    for k, (lab, v, col) in enumerate(rows):
        v = pd.Series(v).dropna()
        ax.scatter(v.clip(-15, 20), k + rng.uniform(-0.12, 0.12, v.size), s=3, color=fs.RULE, lw=0, zorder=1)
        q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
        ax.plot([q1, q3], [k, k], color=col, lw=0.9, zorder=2)
        ax.plot([md], [k], marker="o", ms=3.2, color=col, zorder=3)
        out.append(dict(comparison=lab.replace("\n", " ").replace("−", "-"), n=len(v), median=md,
                        n_positive=int((v > 0).sum()),
                        p_two_sided=stats.wilcoxon(v).pvalue))
    ax.axvline(0, color=fs.MUTED, lw=0.6)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=5.5)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(-15, 20)
    ax.set_xlabel("half-max onset difference (ms)\nboth eyes 45% (MUA) / strong drive (units)")
    ax.text(0.98, 0.98, "> 0: the first\nsignal is later", transform=ax.transAxes, ha="right", va="top",
            fontsize=5.5, color=fs.INK2, linespacing=1.0)
    return pd.DataFrame(out), bias


def main():
    fs.apply()
    Z = np.load(out_root() / "figures" / "data" / "csd_additivity.npz")
    Pn = pd.read_csv(QC / "v4_O_penetrations.csv")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    m = maps([ax["a1"], ax["a2"], ax["a3"]], [ax["acb1"], ax["acb3"]], Z)
    Rb = windows_panel(ax["b"], Pn, {"current": ("current_ai_0.45_{w}", fs.INK, -0.12),
                                     "MUA, same channels": ("mua_ai_0.45_{w}", fs.MUTED, 0.12)},
                       ref=1.0, ylabel="additivity index\nboth eyes / (dominant + other)", ylim=(-0.5, 1.6))
    ax["b"].text(2.45, 1.03, "1: the currents add", ha="right", va="bottom", fontsize=5.5, color=fs.INK2)
    ax["b"].text(0.03, 0.03, "current", transform=ax["b"].transAxes, fontsize=5.5, color=fs.INK)
    ax["b"].text(0.03, 0.11, "MUA, same channels", transform=ax["b"].transAxes, fontsize=5.5, color=fs.MUTED)
    Rc = windows_panel(ax["c"], Pn, {"current": ("current_deficit_0.45_{w}", fs.INK, 0.0)}, ref=0.0,
                       ylabel="current: both eyes −\nif the contrasts add", ylim=(-0.5, 0.5))
    Rd, bias = timing_panel(ax["d"])
    for k, letter in (("a1", "a"), ("b", "b"), ("c", "c"), ("d", "d")):
        fs.panel_label(ax[k], letter, x={"b": -0.32, "c": -0.4, "d": -0.55}.get(k, -0.02))
    print(fs.save(fig, "figS7_input_current"))
    pd.set_option("display.width", 200)
    print("(a)", {k: round(v, 3) for k, v in m.items()})
    print("(b) additivity index, medians:")
    print(Rb.round(3).to_string(index=False))
    print("(c) both eyes - if the contrasts add (current), medians:")
    print(Rc.round(3).to_string(index=False))
    print(f"(d) timing (simulated filter bias removed from the three-contact comparison: {bias:+.2f} ms):")
    print(Rd.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
