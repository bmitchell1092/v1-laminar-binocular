"""Second, corrected pre-registered test (docs/09): single-stage 'full' vs two-stage 'full2', fitted per channel.

Usage: .venv\\Scripts\\python scripts\\fig5_channel_test.py
Pre-registered (docs/09, before implementing):
* models 'full' (8 parameters) and 'full2' (11, nests 'full' as S -> inf); each 'full2' fit starts at the
  single-stage solution embedded at a large S, plus random starts (6 starts per model and half);
* each driven in-cortex MUA channel separately (driven in a monocular condition, run_crf.py), oriented to its
  dominant eye, normalised by its dominant .9 peak, 0-219 ms, 35 non-blank cells; excluded a priori: 211008,
  211210;
* CV: fit odd presentations, score even ones, and the reverse; per channel r = log(CV_full / CV_full2);
* unit of replication = penetration (median r); PASS iff one-sided Wilcoxon of penetration medians > 0 at
  alpha = .025, and the median of penetration medians > 0 in each monkey.
Writes <out_root>/binocular/fig5_channel_test.parquet and docs/qc/fig5_channel_test.csv.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

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
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
T = np.arange(0, 220).astype(float)
CD, CN = np.meshgrid(CONTRASTS, CONTRASTS, indexing="ij")
NONBLANK = (CD + CN) > 0
EXCLUDED = {"211008_B_eA", "211210_B_eA"}
S_EMBED = 400.0
ALPHA = 0.025


def channel_data(pen: str, driven: set) -> list:
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
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    Ms = {}
    for split, sel in (("all", keep), ("odd", keep & (idx % 2 == 1)), ("even", keep & (idx % 2 == 0))):
        M = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= 2:
                    M[i, j] = mua[s].mean(axis=0)
        Ms[split] = gaussian_filter1d(M - M[0, 0], 2, axis=2)
    A = Ms["all"]
    rc, ri = A[-1, 0][w].mean(axis=0), A[0, -1][w].mean(axis=0)
    out = []
    for q, k in enumerate(chans):
        dom_contra = rc[q] >= ri[q]
        scale = (A[-1, 0, w, q] if dom_contra else A[0, -1, w, q]).max()
        if not np.isfinite(scale) or scale <= 0:
            continue
        ys = {}
        for split in ("odd", "even"):
            Mq = Ms[split][..., q] if dom_contra else Ms[split][..., q].transpose(1, 0, 2)
            ys[split] = (Mq[:, :, tsel] / scale)[NONBLANK]
        out.append((pen, str(ds.channel.values[k]), float(ds.depth_below_toc.values[k]), ys["odd"], ys["even"]))
    return out


def fit_channel(task):
    pen, ch, depth, y_odd, y_even = task
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    fits = {}
    for half, y in (("odd", y_odd), ("even", y_even)):
        f1 = dynorm.fit(T, cd, cn, y, "full", starts=6)
        init = f1.params | dict(S=S_EMBED, sigma=max(f1.params["sigma"] / S_EMBED, 1.01e-5), omega=1.0, tau_s=20.0)
        f2 = dynorm.fit(T, cd, cn, y, "full2", starts=5, init=init)
        fits[half] = (f1, f2)
    cv1 = dynorm.sse(fits["odd"][0].params, T, cd, cn, y_even, "full") + dynorm.sse(fits["even"][0].params, T, cd, cn, y_odd, "full")
    cv2 = dynorm.sse(fits["odd"][1].params, T, cd, cn, y_even, "full2") + dynorm.sse(fits["even"][1].params, T, cd, cn, y_odd, "full2")
    p2 = {k: float(np.mean([fits[h][1].params[k] for h in ("odd", "even")])) for k in ("S", "omega", "tau_s", "w", "a", "tau2")}
    return dict(penetration=pen, monkey=pen[7], channel=ch, depth=depth, cv_full=cv1, cv_full2=cv2,
                train_full=fits["odd"][0].sse + fits["even"][0].sse, train_full2=fits["odd"][1].sse + fits["even"][1].sse,
                **{f"p2_{k}": v for k, v in p2.items()})


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven} - EXCLUDED)
    tasks = []
    for pen in pens:
        tasks += channel_data(pen, driven)
    print(f"{len(tasks)} channels in {len(pens)} penetrations", flush=True)
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=20) as ex:
        for i, r in enumerate(ex.map(fit_channel, tasks, chunksize=1)):
            rows.append(r)
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(tasks)} channels ({time.time() - t0:.0f}s)", flush=True)
    C = pd.DataFrame(rows)
    out = out_root() / "binocular"
    C.to_parquet(out / "fig5_channel_test.parquet")

    C["r"] = np.log(C.cv_full / C.cv_full2)
    per = C.groupby(["penetration", "monkey"]).r.median()
    w = stats.wilcoxon(per.values, alternative="greater")
    med_mk = per.groupby(level="monkey").median()
    passed = (w.pvalue < ALPHA) and (med_mk > 0).all()
    rows = [dict(item="channels", value=len(C)), dict(item="penetrations", value=len(per)),
            dict(item="one-sided Wilcoxon p (penetration medians > 0)", value=w.pvalue),
            dict(item="median of penetration medians r", value=per.median()),
            dict(item="median r, monkey B", value=med_mk.get("B", np.nan)),
            dict(item="median r, monkey J", value=med_mk.get("J", np.nan)),
            dict(item="penetrations with median r > 0", value=(per > 0).mean()),
            dict(item="channels with r > 0", value=(C.r > 0).mean()),
            dict(item="median % change in held-out error (full2 vs full)", value=100 * (np.exp(-C.r.median()) - 1)),
            dict(item="training SSE full2 <= full (fraction of channels)", value=(C.train_full2 <= C.train_full * 1.001).mean()),
            dict(item="PRE-REGISTERED RESULT", value="PASS" if passed else "FAIL")]
    for k in ("omega", "S", "tau_s", "w", "a", "tau2"):
        rows.append(dict(item=f"full2 {k}: median [IQR]",
                         value=f"{C[f'p2_{k}'].median():.3g} [{C[f'p2_{k}'].quantile(.25):.3g}, {C[f'p2_{k}'].quantile(.75):.3g}]"))
    S = pd.DataFrame(rows)
    print(S.to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "fig5_channel_test.csv", index=False)


if __name__ == "__main__":
    main()
