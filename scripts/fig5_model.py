"""Figure 5: a delayed binocular normalization model fitted per penetration (v1lb.dynorm; docs/09).

Usage: .venv\\Scripts\\python scripts\\fig5_model.py
Data: per penetration, driven in-cortex MUA channels (driven in a monocular condition, run_crf.py); condition means
for all 36 cells (completed presentations, no saccade >= 1 deg), Gaussian sigma 2 ms, blank-subtracted, oriented
so axis 0 is the channel's dominant eye (larger monocular response at .9, 40-250 ms), normalised by the dominant
.9 peak, averaged over channels. Time 0-219 ms after the photodiode flash.
Models: full, static (no pool delay), w = a (pool weighted like excitation). Cross-validation: fit odd
presentations and score even ones, and the reverse; CV error = summed held-out SSE. Parameters from fits to all
presentations. Model L(t) is computed from the fitted predictions with the same onset-aligned windows as the data.
Writes <out_root>/binocular/fig5_model.parquet, docs/qc/fig5_model.csv and <out_root>/binocular/fig5_model.png.
"""

import glob
import time
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
T = np.arange(0, 220)
MODELS = ["full", "static", "w=a"]
L_WINDOWS = [(0, 10), (10, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150)]
CD, CN = np.meshgrid(CONTRASTS, CONTRASTS, indexing="ij")
NONBLANK = (CD + CN) > 0


def population(pen: str, driven: set):
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if len(chans) < 3:
        return None
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    out = {}
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    orient = None
    for split, sel in (("all", keep), ("odd", keep & (idx % 2 == 1)), ("even", keep & (idx % 2 == 0))):
        M = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= 2:
                    M[i, j] = mua[s].mean(axis=0)
        M = gaussian_filter1d(M - M[0, 0], 2, axis=2)
        if orient is None:                                       # dominant eye and scale from all presentations
            rc, ri = M[-1, 0][w].mean(axis=0), M[0, -1][w].mean(axis=0)
            dom_contra = rc >= ri
            scale = np.where(dom_contra, M[-1, 0][w].max(axis=0), M[0, -1][w].max(axis=0))
            orient = (dom_contra, scale)
        dom_contra, scale = orient
        Mo = np.where(dom_contra[None, None, None, :], M, M.transpose(1, 0, 2, 3))
        good = scale > 0
        out[split] = np.nanmean(Mo[..., good] / scale[good], axis=3)[:, :, tsel]
    return out, int(np.sum(orient[1] > 0))


def windowed_L(Y: np.ndarray) -> dict:
    """Onset-aligned windowed L for dioptic .11/.225/.45 from a [6, 6, T] array (dominant eye on axis 0)."""
    rd, rn = Y[:, 0], Y[0, :]
    pk = np.nanmax(rd[-1])
    on_idx = np.flatnonzero(rd[-1] >= 0.5 * pk)
    if not on_idx.size:
        return {}
    on = T[on_idx[0]]
    out = {}
    for c in (0.11, 0.225, 0.45):
        i = list(CONTRASTS).index(c)
        for w0, w1 in L_WINDOWS:
            sel = (T >= on + w0) & (T < on + w1)
            if not sel.any():
                continue
            d, n, b = rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1), Y[i, i, sel].mean()
            a = eye_weight(d, n)
            k, cens = equivalent_contrast(d[:, None], np.array([b]), c)
            k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
            out[(c, w0)] = (k - 1) / a
    return out


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    rows, Lrows, keep_example = [], [], {}
    for pen in pens:
        t0 = time.time()
        pop = population(pen, driven)
        if pop is None:
            continue
        Y, n_ch = pop
        ys = {k: v[NONBLANK] for k, v in Y.items()}
        res = dict(penetration=pen, monkey=pen[7], n_channels=n_ch)
        for m in MODELS:
            f_odd = dynorm.fit(T, cd, cn, ys["odd"], m, starts=4)
            f_even = dynorm.fit(T, cd, cn, ys["even"], m, starts=4)
            cv = dynorm.sse(f_odd.params, T, cd, cn, ys["even"], m) + dynorm.sse(f_even.params, T, cd, cn, ys["odd"], m)
            f_all = dynorm.fit(T, cd, cn, ys["all"], m, starts=4)
            res[f"cv_{m}"] = cv
            res[f"sse_{m}"] = f_all.sse
            if m == "full":
                res |= {f"p_{k}": v for k, v in f_all.params.items()}
                pred = np.full((6, 6, T.size), np.nan)
                pred[NONBLANK] = dynorm.predict(f_all.params, T, cd, cn)
                pred[0, 0] = 0
                for (c, w0), v in windowed_L(pred).items():
                    Lrows.append(dict(penetration=pen, source="model", c=c, w0=w0, L=v))
                keep_example[pen] = (Y["all"], pred)
        for (c, w0), v in windowed_L(Y["all"]).items():
            Lrows.append(dict(penetration=pen, source="data", c=c, w0=w0, L=v))
        res["var_all"] = float(np.nansum((ys["all"] - np.nanmean(ys["all"])) ** 2))
        rows.append(res)
        print(f"{pen}: {n_ch} ch, CV full {res['cv_full']:.2f} static {res['cv_static']:.2f} w=a {res['cv_w=a']:.2f}; "
              f"a {res['p_a']:.2f} w {res['p_w']:.2f} tau2 {res['p_tau2']:.0f} n {res['p_n']:.2f} "
              f"({time.time() - t0:.0f}s)", flush=True)
    F = pd.DataFrame(rows)
    Ld = pd.DataFrame(Lrows)
    out = out_root() / "binocular"
    F.to_parquet(out / "fig5_model.parquet")
    Ld.to_parquet(out / "fig5_model_L.parquet")

    # ---- tests
    summ = []
    for other in ("static", "w=a"):
        lr = np.log(F[f"cv_{other}"] / F["cv_full"])
        mk = F.monkey
        summ.append(dict(test=f"CV: full beats {other}", n=len(F), frac=(lr > 0).mean(),
                         frac_B=(lr[mk == "B"] > 0).mean(), frac_J=(lr[mk == "J"] > 0).mean(),
                         median=np.exp(lr.median()), p=stats.wilcoxon(lr).pvalue))
    d = F.p_w - F.p_a
    summ.append(dict(test="w > a (full fit, all presentations)", n=len(F), frac=(d > 0).mean(),
                     frac_B=(d[F.monkey == "B"] > 0).mean(), frac_J=(d[F.monkey == "J"] > 0).mean(),
                     median=d.median(), p=stats.wilcoxon(F.p_w, F.p_a).pvalue))
    r2 = 1 - F.sse_full / F.var_all
    summ.append(dict(test="variance explained, full (all presentations)", n=len(F), frac=r2.median(),
                     frac_B=r2[F.monkey == "B"].median(), frac_J=r2[F.monkey == "J"].median(), median=r2.min(), p=np.nan))
    for k in ("tau2", "n", "sigma", "delta", "a", "w"):
        summ.append(dict(test=f"param {k}: median [IQR]", n=len(F), frac=F[f"p_{k}"].median(),
                         frac_B=F[f"p_{k}"].quantile(.25), frac_J=F[f"p_{k}"].quantile(.75), median=np.nan, p=np.nan))
    # data vs model L: per window and contrast, median over penetrations; correlation of penetration-level values
    for c in (0.11, 0.225, 0.45):
        piv = Ld[Ld.c == c].pivot_table(index=["penetration", "w0"], columns="source", values="L").dropna()
        piv = piv[(piv.abs() < 5).all(axis=1)]
        summ.append(dict(test=f"L data vs model, c={c}: Spearman over penetration x window", n=len(piv),
                         frac=stats.spearmanr(piv.data, piv.model).correlation, frac_B=np.nan, frac_J=np.nan,
                         median=np.nan, p=stats.spearmanr(piv.data, piv.model).pvalue))
    S = pd.DataFrame(summ)
    pd.set_option("display.width", 200)
    print(S.round(3).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "fig5_model.csv", index=False)

    # ---- figure
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2), facecolor="#fcfcfb")
    Fi = F.set_index("penetration")
    stable = Fi[Fi.cv_full <= Fi.cv_static]                            # skip unstable fits (e.g. 211008, ~4 per cell per half)
    ex = (1 - stable.sse_full / stable.var_all).idxmax()
    Yd, Yp = keep_example[ex]
    cols = {0.11: "#2a78d6", 0.225: "#eb6834", 0.45: "#1baf7a"}
    for c in (0.11, 0.225, 0.45):
        i = list(CONTRASTS).index(c)
        axs[0].plot(T, Yd[i, 0], color=cols[c], lw=1, alpha=0.6)
        axs[0].plot(T, Yp[i, 0], color=cols[c], lw=1.6, ls="--")
        axs[1].plot(T, Yd[i, i], color=cols[c], lw=1, alpha=0.6, label=f"c = {c:g}")
        axs[1].plot(T, Yp[i, i], color=cols[c], lw=1.6, ls="--")
    axs[0].set_title(f"{ex}: dominant eye alone (data solid, model dashed)", fontsize=8)
    axs[1].set_title("both eyes at c (dioptic)", fontsize=8)
    axs[1].legend(fontsize=7, frameon=False)
    for ax in axs[:2]:
        ax.set_xlabel("ms from photodiode flash", fontsize=8)
    for c in (0.11, 0.225, 0.45):
        for src, ls in (("data", "-"), ("model", "--")):
            g = Ld[(Ld.c == c) & (Ld.source == src)].groupby("w0").L.median()
            xs = [(a + b) / 2 for a, b in L_WINDOWS if a in g.index]
            axs[2].plot(xs, g.values, color=cols[c], ls=ls, marker="o" if src == "data" else None, ms=3, lw=1.4)
    axs[2].axhline(1, color="#0b0b0b", lw=0.8, ls=":")
    axs[2].axhline(0, color="#52514e", lw=0.8)
    axs[2].set_ylim(-1.2, 1.8)
    axs[2].set_title("linearity L(t): data (solid) vs model (dashed), median of penetrations", fontsize=8)
    axs[2].set_xlabel("ms from onset", fontsize=8)
    axs[3].scatter(F.p_a, F.p_w, s=18, c=np.where(F.monkey == "B", "#eb6834", "#2a78d6"))
    lim = [0, max(F.p_a.max(), F.p_w.max()) * 1.05]
    axs[3].plot(lim, lim, color="#52514e", lw=0.8)
    axs[3].set_xlabel("excitatory weight a (non-dominant eye)", fontsize=8)
    axs[3].set_ylabel("pool weight w (non-dominant eye)", fontsize=8)
    axs[3].set_title("pool vs excitation weight (B orange, J blue)", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig5_model.png", dpi=85, facecolor="#fcfcfb")
    print(out / "fig5_model.png")


if __name__ == "__main__":
    main()
