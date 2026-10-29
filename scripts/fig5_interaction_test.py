"""Third pre-registered model test (docs/09): interaction-targeted fits of 'full', 'full2' and dynamic
Ding-Sperling 'ds', per channel.

Usage: .venv\\Scripts\\python scripts\\fig5_interaction_test.py
As registered:
* loss = monocular time courses + binocular interaction R(c_D, c_N) - R(c_D, 0), each normalised
  (v1lb.dynorm.interaction_residuals); 'full2' seeded at the single-stage solution;
* per driven channel (as fig5_channel_test.py; 211008, 211210 excluded); odd/even CV of the interaction loss;
  all-presentation fits for L;
* PRIMARY (capture, per model): late L (windows starting >= 50 ms after the channel's onset, mean of 50-75,
  75-100, 100-150) from model and data; penetration median of the signed error (model - data); captures iff the
  90% bootstrap CI over penetrations lies within +-0.25 at .11, .225 and .45;
* SECONDARY: pairwise held-out loss (ds vs full, full2 vs full, ds vs full2), one-sided Wilcoxon on penetration
  medians of log ratios; claimed iff p < .005 and the median is positive in both monkeys.
Writes <out_root>/binocular/fig5_interaction_test.parquet and docs/qc/fig5_interaction_test.csv.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import glob  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
T = np.arange(0, 220).astype(float)
CD, CN = np.meshgrid(CONTRASTS, CONTRASTS, indexing="ij")
NONBLANK = (CD + CN) > 0
EXCLUDED = {"211008_B_eA", "211210_B_eA"}
S_EMBED = 400.0
WINDOWS = [(0, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150)]
LATE = (50, 75, 100)
MODELS = ["full", "full2", "ds"]
MARGIN = 0.25


def channel_data(pen, driven):
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return []
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    Ms = {}
    for split, sel in (("all", keep), ("odd", keep & (idx % 2 == 1)), ("even", keep & (idx % 2 == 0))):
        M = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= 2:
                    M[i, j] = mua[s].mean(axis=0)
        Ms[split] = gaussian_filter1d(M - M[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    out = []
    for q, k in enumerate(chans):
        A = Ms["all"][..., q]
        dom_contra = A[-1, 0, w].mean() >= A[0, -1, w].mean()
        orient = (lambda X: X) if dom_contra else (lambda X: X.transpose(1, 0, 2))
        scale = orient(A)[-1, 0, w].max()
        if not np.isfinite(scale) or scale <= 0:
            continue
        ys = {s: orient(Ms[s][..., q])[:, :, tsel] / scale for s in ("all", "odd", "even")}
        if np.isnan(ys["all"]).any():
            continue
        out.append((pen, str(ds.channel.values[k]), ys))
    return out


def windowed_L(Y, onset):
    rd, rn = Y[:, 0], Y[0, :]
    res = {}
    for c in (0.11, 0.225, 0.45):
        i = list(CONTRASTS).index(c)
        for w0, w1 in WINDOWS:
            sel = (T >= onset + w0) & (T < onset + w1)
            if not sel.any():
                continue
            d, n, b = rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1), Y[i, i, sel].mean()
            a = eye_weight(d, n)
            k, cens = equivalent_contrast(d[:, None], np.array([b]), c)
            k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
            res[(c, w0)] = (k - 1) / a
    return res


def fit_model(m, y, init_full=None):
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    if m == "full2":
        init = init_full | dict(S=S_EMBED, sigma=max(init_full["sigma"] / S_EMBED, 1.01e-5), omega=1.0, tau_s=20.0)
        return dynorm.fit(T, cd, cn, y, "full2", starts=5, init=init, loss="interaction")
    return dynorm.fit(T, cd, cn, y, m, starts=6, loss="interaction")


def work(task):
    pen, ch, ys = task
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    y = {s: ys[s][NONBLANK] for s in ys}
    res = dict(penetration=pen, monkey=pen[7], channel=ch)
    fits = {}
    for split in ("odd", "even", "all"):
        f_full = fit_model("full", y[split])
        fits[(split, "full")] = f_full
        fits[(split, "full2")] = fit_model("full2", y[split], f_full.params)
        fits[(split, "ds")] = fit_model("ds", y[split])
    for m in MODELS:
        res[f"cv_{m}"] = (dynorm.loss_value(fits[("odd", m)].params, T, cd, cn, y["even"])
                          + dynorm.loss_value(fits[("even", m)].params, T, cd, cn, y["odd"]))
    for k in ("a", "rho", "gamma", "tau_g", "tau2"):
        res[f"ds_{k}"] = fits[("all", "ds")].params[k]
    on_idx = np.flatnonzero(ys["all"][-1, 0] >= 0.5 * np.nanmax(ys["all"][-1, 0]))
    if not on_idx.size:
        return res, []
    onset = T[on_idx[0]]
    rows = []
    for src in ["data"] + MODELS:
        if src == "data":
            Z = ys["all"]
        else:
            Z = np.zeros((6, 6, T.size))
            Z[NONBLANK] = dynorm.predict(fits[("all", src)].params, T, cd, cn)
        for (c, w0), v in windowed_L(Z, onset).items():
            rows.append(dict(penetration=pen, monkey=pen[7], channel=ch, source=src, c=c, w0=w0, L=v))
    return res, rows


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    tasks = []
    for pen in sorted({p for p, _ in driven} - EXCLUDED):
        tasks += channel_data(pen, driven)
    print(f"{len(tasks)} channels", flush=True)
    t0 = time.time()
    resrows, Lrows = [], []
    with ProcessPoolExecutor(max_workers=20) as ex:
        for i, (r, lr) in enumerate(ex.map(work, tasks, chunksize=1)):
            resrows.append(r)
            Lrows += lr
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(tasks)} ({time.time() - t0:.0f}s)", flush=True)
    C, L = pd.DataFrame(resrows), pd.DataFrame(Lrows)
    out = out_root() / "binocular"
    C.to_parquet(out / "fig5_interaction_test.parquet")
    L.to_parquet(out / "fig5_interaction_test_L.parquet")

    rng = np.random.default_rng(0)
    summ = []
    # PRIMARY: capture of the late collapse
    late = L[L.w0.isin(LATE)].groupby(["penetration", "monkey", "channel", "source", "c"]).L.mean().unstack("source")
    for m in MODELS:
        captured = True
        for c in (0.11, 0.225, 0.45):
            g = late.xs(c, level="c")[["data", m]].dropna()
            per = (g[m] - g["data"]).groupby(level="penetration").median()
            boot = [np.median(rng.choice(per.values, len(per))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            ok = (lo > -MARGIN) and (hi < MARGIN)
            captured &= ok
            summ.append(dict(section="PRIMARY capture", model=m, c=c, n_pen=len(per), median_signed_error=per.median(),
                             ci90_lo=lo, ci90_hi=hi, within_margin=ok, p=np.nan))
        summ.append(dict(section="PRIMARY capture", model=m, c="all three", n_pen=np.nan, median_signed_error=np.nan,
                         ci90_lo=np.nan, ci90_hi=np.nan, within_margin=captured, p=np.nan))
    # SECONDARY: pairwise held-out interaction loss (second model better if log ratio > 0)
    for better, worse in (("ds", "full"), ("full2", "full"), ("ds", "full2")):
        r = np.log(C[f"cv_{worse}"] / C[f"cv_{better}"])
        per = r.groupby([C.penetration, C.monkey]).median()
        p = stats.wilcoxon(per.values, alternative="greater").pvalue
        mk = per.groupby(level=1).median()
        claim = (p < 0.005) and (mk > 0).all()
        summ.append(dict(section="SECONDARY architecture", model=f"{better} beats {worse}", c="", n_pen=len(per),
                         median_signed_error=per.median(), ci90_lo=mk.get("B", np.nan), ci90_hi=mk.get("J", np.nan),
                         within_margin=claim, p=p))
    S = pd.DataFrame(summ)
    pd.set_option("display.width", 220)
    print("(PRIMARY: median_signed_error = model - data late L; within_margin = 90% CI inside +-0.25. "
          "SECONDARY: median_signed_error = median log loss ratio, ci90_lo/hi = monkey B/J medians, within_margin = claimed)")
    print(S.to_string(index=False))
    print("ds parameters, median [IQR]:", {k: f"{C[f'ds_{k}'].median():.3g} [{C[f'ds_{k}'].quantile(.25):.3g}, {C[f'ds_{k}'].quantile(.75):.3g}]"
                                           for k in ("a", "rho", "gamma", "tau_g", "tau2")})
    S.to_csv(REPO / "docs" / "qc" / "fig5_interaction_test.csv", index=False)


if __name__ == "__main__":
    main()
