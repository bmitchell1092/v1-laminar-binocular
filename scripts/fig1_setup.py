"""Figure 1 draft: design, alignment and data quality (docs/09).

Usage: .venv\\Scripts\\python scripts\\fig1_setup.py [example penetration, default 221212_J_eA]
(a) the 6 x 6 dichoptic contrast design: completed presentations per cell for the example penetration;
(b) example penetration: flash-evoked CSD (bad contacts interpolated, Hamming-smoothed) with the top of cortex and
    the L4/5 landmark (laminar alignment v2);
(c) photodiode-measured stimulus onset relative to the event code, all presentations (why onset must be measured);
(d) example penetration (channel average, dominant eye first) at .45: each eye alone, both eyes, and the dominant
    eye at .9 (the linear-contrast-summation reference, (1 + a) c with a ~ 0.9);
(e) population monocular CRFs per eye (driven MUA, 40-250 ms, normalised to each channel's best response);
(f) ocular dominance (ODI at .9) for driven MUA channels and single units.
Writes <out_root>/figures/fig1_setup.png.
"""

import glob
import sys
import warnings

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v1lb import laminar as lam  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import metadata_dir, out_root  # noqa: E402
from v1lb.dataset import open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7)


def main():
    pen = sys.argv[1] if len(sys.argv) > 1 else "221212_J_eA"
    d = out_root() / "datasets"
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration").loc[pen]
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
    fig, axs = plt.subplots(2, 3, figsize=(15, 9), facecolor=SURF)

    # (a) design
    ax = axs[0, 0]
    n = pd.crosstab(P[keep].contrast_ipsi, P[keep].contrast_contra).reindex(index=CONTRASTS, columns=CONTRASTS)
    ax.imshow(n.values, cmap="Greys", vmin=0, vmax=n.values.max() * 2, origin="lower")
    for i in range(6):
        for j in range(6):
            ax.text(j, i, int(n.values[i, j]), ha="center", va="center", fontsize=7, color=INK)
    ax.set_xticks(range(6))
    ax.set_xticklabels([f"{c:g}" for c in CONTRASTS], fontsize=7)
    ax.set_yticks(range(6))
    ax.set_yticklabels([f"{c:g}" for c in CONTRASTS], fontsize=7)
    ax.set_xlabel("contralateral eye contrast", fontsize=8, color=INK2)
    ax.set_ylabel("ipsilateral eye contrast", fontsize=8, color=INK2)
    ax.set_title(f"a  6 x 6 dichoptic design ({pen}; completed presentations)", fontsize=9, color=INK, loc="left")

    # (b) flash CSD with landmarks
    ax = axs[0, 1]
    style(ax)
    f = sorted(glob.glob(str(out_root() / "flash_csd" / f"{pen}_evp*.npz")))[0]
    z = np.load(f)
    lfp, t, con = z["lfp"], z["t_ms"], z["contacts"]
    bad = lam.bad_channels(lfp[(t >= 20) & (t <= 150)])
    csd = lam.csd(lfp, smooth=True, bad=bad)
    order = np.argsort(con) if al.orientation == "ascending" else np.argsort(-con)
    C = csd[:, order].T
    C -= np.nanmean(C[:, (t >= -50) & (t < 0)], axis=1, keepdims=True)
    lim = np.nanpercentile(np.abs(C[:, (t > 20) & (t < 150)]), 98)
    ax.imshow(C, aspect="auto", cmap="RdBu", vmin=-lim, vmax=lim, extent=[t[0], t[-1], len(con) - 0.5, -0.5])
    ax.set_xlim(-20, 150)
    cc_sorted = con[order]
    row = lambda c: np.interp(c, cc_sorted, np.arange(len(con))) if al.orientation == "ascending" else np.interp(-c, -cc_sorted, np.arange(len(con)))  # noqa: E731
    ax.axhline(row(al.top_of_cortex) - 0.5, color=INK, lw=1.2, ls="--", label="top of cortex (spontaneous MUA)")
    ax.axhline(row(al.l45), color=AQUA, lw=1.6, label="L4/5 (flash sink, alignment v2)")
    ax.set_xlabel("ms from flash", fontsize=8, color=INK2)
    ax.set_ylabel("contact (superficial at top)", fontsize=8, color=INK2)
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    ax.set_title("b  flash-evoked CSD (sinks blue) and laminar landmarks", fontsize=9, color=INK, loc="left")

    # (c) photodiode latency
    ax = axs[0, 2]
    style(ax)
    lat = pd.concat([pd.read_parquet(f).pd_on_lat_ms for f in glob.glob(str(d / "*_presentations.parquet"))])
    lat = lat[np.isfinite(lat)]
    ax.hist(lat, bins=np.arange(0, 45, 0.5), color=BLUE)
    ax.set_xlabel("photodiode flash minus event code (ms)", fontsize=8, color=INK2)
    ax.set_ylabel("presentations", fontsize=8, color=INK2)
    ax.set_title(f"c  stimulus onset varies by whole frames (n = {len(lat):,})", fontsize=9, color=INK, loc="left")

    # (d) example penetration: channel average (each channel oriented to its dominant eye) at .45
    ax = axs[1, 0]
    style(ax)
    ds = open_penetration(pen)
    t = ds.time.values
    R0 = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    dch = set(R0[(R0.penetration == pen) & (R0.site == "mua") & R0.condition.str.startswith("mono") & R0.driven].site_id)
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and str(ch) in dch]
    mua = ds.mua.values[keep][:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy()[keep], P.contrast_ipsi.to_numpy()[keep]
    w = (t >= 40) & (t < 250)
    mean = lambda a, b: mua[np.isclose(cc, a) & np.isclose(ci, b)].mean(axis=0)  # noqa: E731
    blank = mean(0, 0)
    c_ = 0.45
    rc, ri, both = mean(c_, 0) - blank, mean(0, c_) - blank, mean(c_, c_) - blank
    dom_c = (mean(.9, 0) - blank)[w].mean(axis=0) >= (mean(0, .9) - blank)[w].mean(axis=0)
    dom = np.where(dom_c, rc, ri)
    non = np.where(dom_c, ri, rc)
    scale = np.where(dom_c, mean(.9, 0) - blank, mean(0, .9) - blank)[w].max(axis=0)
    avg = lambda x: (x / scale).mean(axis=1)  # noqa: E731
    dom9 = np.where(dom_c, mean(.9, 0) - blank, mean(0, .9) - blank)
    # linear summation of contrast predicts the dominant eye at (1 + a) c, about .9 here (a ~ 0.9); the sum of the
    # two monocular *responses* is not a linear prediction (compressive CRF)
    for y, col, ls, lab in ((avg(dom), BLUE, "-", "dominant eye alone at .45"), (avg(non), ORANGE, "-", "other eye alone at .45"),
                            (avg(dom9), INK2, "--", "dominant eye at .9 (~ linear contrast summation)"),
                            (avg(both), INK, "-", "both eyes at .45")):
        ax.plot(t, y, color=col, lw=1.4, ls=ls, label=lab)
    ax.set_xlim(-20, 250)
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_xlabel("ms from photodiode flash", fontsize=8, color=INK2)
    ax.set_ylabel("MUA / dominant-eye peak at .9", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"d  {pen}, {len(chans)} driven channels, .45 contrast", fontsize=9, color=INK, loc="left")

    # (e) population monocular CRFs and (f) ocular dominance
    D = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "dichoptic" / "*.parquet"))], ignore_index=True)
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[R.condition.str.startswith("mono") & R.driven].groupby(["penetration", "site", "site_id"]).size()
    D = D.join(drv.rename("nd"), on=["penetration", "site", "site_id"])
    D = D[D.nd.notna() & (D.window == "full")]
    isolated = (D.isi_contamination < 0.5) & (D.snr >= 3)                     # lenient isolation (05 C14)
    D = D[(D.site == "mua") | ((D.site == "unit") & isolated)]
    m = D[D.site == "mua"].set_index(["penetration", "site_id", "contrast_contra", "contrast_ipsi"])["mean"]
    curves = {"contralateral eye": [], "ipsilateral eye": []}
    odi = {"mua": [], "unit": []}
    for (p_, s_), g in D.groupby(["penetration", "site_id"]):
        mm = g.set_index(["contrast_contra", "contrast_ipsi"])["mean"]
        try:
            con_c = np.array([mm[(c, 0.0)] for c in CONTRASTS])
            ips_c = np.array([mm[(0.0, c)] for c in CONTRASTS])
        except KeyError:
            continue
        top = max(con_c.max(), ips_c.max())
        if top <= 0:
            continue
        site = g.site.iloc[0]
        if site == "mua":
            curves["contralateral eye"].append(con_c / top)
            curves["ipsilateral eye"].append(ips_c / top)
        odi[site].append((con_c[-1] - ips_c[-1]) / (abs(con_c[-1]) + abs(ips_c[-1])))
    ax = axs[1, 1]
    style(ax)
    for (lab, cs), col in zip(curves.items(), (BLUE, ORANGE)):
        A = np.array(cs)
        med = np.median(A, axis=0)
        lo, hi = np.percentile(A, [25, 75], axis=0)
        ax.plot(range(6), med, color=col, marker="o", ms=4, lw=1.5, label=f"{lab} (n = {len(A)})")
        ax.fill_between(range(6), lo, hi, color=col, alpha=0.15, lw=0)
    ax.set_xticks(range(6))
    ax.set_xticklabels([f"{c:g}" for c in CONTRASTS], fontsize=7)
    ax.set_xlabel("contrast (other eye 0)", fontsize=8, color=INK2)
    ax.set_ylabel("response / best monocular response", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("e  monocular contrast responses, driven MUA (median, IQR)", fontsize=9, color=INK, loc="left")
    ax = axs[1, 2]
    style(ax)
    bins = np.linspace(-1, 1, 21)
    ax.hist(odi["mua"], bins=bins, color=BLUE, alpha=0.8, label=f"MUA channels (n = {len(odi['mua'])})")
    ax.hist(odi["unit"], bins=bins, color=ORANGE, alpha=0.8, label=f"isolated single units (n = {len(odi['unit'])})")
    ax.set_xlabel("ocular dominance index (contra - ipsi) / sum at .9", fontsize=8, color=INK2)
    ax.set_ylabel("sites", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("f  ocular dominance", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    outd = out_root() / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    fn = outd / "fig1_setup.png"
    fig.savefig(fn, dpi=90, facecolor=SURF)
    print(fn)


if __name__ == "__main__":
    main()
