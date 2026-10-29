"""Descriptive (not pass/fail): does the two-stage model ('full2') reproduce the linearity time course L(t)
better than single-stage 'full', especially the low-contrast misfit? (docs/09, after the channel-level test)

Usage: .venv\\Scripts\\python scripts\\fig5_L_check.py
Per driven channel (same rule and exclusions as fig5_channel_test.py): both models fitted to all presentations
('full2' seeded at the single-stage solution), then onset-aligned windowed L (window-specific eye weights) from the
data and from each model's prediction, for dioptic .11, .225, .45. Penetration medians per window; compares
median |L_model - L_data| and the signed error late at .11.
Writes <out_root>/binocular/fig5_L_check.parquet and docs/qc/fig5_L_check.csv.
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
    mua = ds.mua.values[keep][:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy()[keep], P.contrast_ipsi.to_numpy()[keep]
    M = np.full((6, 6, t.size, len(chans)), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= 3:
                M[i, j] = mua[s].mean(axis=0)
    M = gaussian_filter1d(M - M[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    out = []
    for q, k in enumerate(chans):
        dom_contra = M[-1, 0, w, q].mean() >= M[0, -1, w, q].mean()
        Mq = M[..., q] if dom_contra else M[..., q].transpose(1, 0, 2)
        scale = Mq[-1, 0, w].max()
        if not np.isfinite(scale) or scale <= 0 or np.isnan(Mq[:, :, tsel]).any():
            continue
        out.append((pen, str(ds.channel.values[k]), Mq[:, :, tsel] / scale))
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


def work(task):
    pen, ch, Y = task
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    y = Y[NONBLANK]
    f1 = dynorm.fit(T, cd, cn, y, "full", starts=6)
    init = f1.params | dict(S=S_EMBED, sigma=max(f1.params["sigma"] / S_EMBED, 1.01e-5), omega=1.0, tau_s=20.0)
    f2 = dynorm.fit(T, cd, cn, y, "full2", starts=5, init=init)
    on_idx = np.flatnonzero(Y[-1, 0] >= 0.5 * np.nanmax(Y[-1, 0]))
    if not on_idx.size:
        return []
    onset = T[on_idx[0]]                                  # data onset, used for all three sources
    rows = []
    for src, Z in (("data", Y), ("full", None), ("full2", None)):
        if Z is None:
            Z = np.zeros((6, 6, T.size))
            Z[NONBLANK] = dynorm.predict((f1 if src == "full" else f2).params, T, cd, cn)
        for (c, w0), v in windowed_L(Z, onset).items():
            rows.append(dict(penetration=pen, monkey=pen[7], channel=ch, source=src, c=c, w0=w0, L=v))
    return rows


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    tasks = []
    for pen in sorted({p for p, _ in driven} - EXCLUDED):
        tasks += channel_data(pen, driven)
    print(f"{len(tasks)} channels", flush=True)
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=20) as ex:
        for i, r in enumerate(ex.map(work, tasks, chunksize=1)):
            rows += r
            if (i + 1) % 100 == 0:
                print(f"  {i + 1}/{len(tasks)} ({time.time() - t0:.0f}s)", flush=True)
    L = pd.DataFrame(rows)
    out = out_root() / "binocular"
    L.to_parquet(out / "fig5_L_check.parquet")
    per = L.groupby(["penetration", "monkey", "source", "c", "w0"]).L.median().unstack("source").dropna()
    summ = []
    for (c, w0), g in per.groupby(level=["c", "w0"]):
        summ.append(dict(c=c, w0=w0, n_pen=len(g), L_data=g.data.median(), L_full=g.full.median(), L_full2=g.full2.median(),
                         abs_err_full=(g.full - g.data).abs().median(), abs_err_full2=(g.full2 - g.data).abs().median(),
                         full2_closer=float(((g.full2 - g.data).abs() < (g.full - g.data).abs()).mean())))
    S = pd.DataFrame(summ)
    pd.set_option("display.width", 200)
    print(S.round(3).to_string(index=False))
    late11 = per.xs(0.11, level="c")
    late11 = late11[late11.index.get_level_values("w0") >= 50]
    print(f"\nlate (>= 50 ms) at .11: signed error full {np.median(late11.full - late11.data):+.3f}, "
          f"full2 {np.median(late11.full2 - late11.data):+.3f}")
    S.to_csv(REPO / "docs" / "qc" / "fig5_L_check.csv", index=False)


if __name__ == "__main__":
    main()
