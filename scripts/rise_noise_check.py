"""Noise check on the 10-90% rise time (docs/09, correction of 2026-10-31): is the slow low-contrast rise real?

Usage: .venv\\Scripts\\python scripts\\rise_noise_check.py [n_workers]
Per penetration (the channel-averaged time courses the models were fitted to, fig5_model.population): the dominant
eye alone at every contrast for the data, the data's own noise sample (odd - even) / 2 (same variance as the noise in
the all-presentation mean), and the all-presentation fits of 'full', 'cond' and 'org' (fig5_cond.parquet,
fig5_org.parquet).
* per penetration: 10-90% rise (as fig5_org_model.timing) of the data, of each model, and of each model with the
  noise sample added and subtracted;
* population mean: rise of the mean over penetrations (noise ~1/5 of one penetration's), data and models, with a
  2000-sample penetration bootstrap of each and of model - data.
Writes <out_root>/binocular/dominant_traces.npz (traces [penetration, source, contrast, time]) and
docs/qc/rise_noise.csv; prints both tables.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig5_model import CD, CN, NONBLANK, T, population  # noqa: E402
from fig5_org_model import full6  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
SOURCES = ("data", "noise", "full", "cond", "org")
MODELS = ("full", "cond", "org")


def fitted():
    b = out_root() / "binocular"
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    O = pd.read_parquet(b / "fig5_org.parquet").set_index("penetration")
    return F.join(O[[c for c in O.columns if c.startswith("org_p_")]])


def params(r, m):
    pre = f"{m}_p_"
    return {k[len(pre):]: float(v) for k, v in r.items() if k.startswith(pre) and np.isfinite(v)}


def rise(x):
    """10-90% of the peak within 150 ms, first crossings (fig5_org_model.timing)."""
    pk = np.nanmax(x[T < 150])
    if not np.isfinite(pk) or pk <= 0:
        return np.nan
    t10, t90 = (np.flatnonzero(x >= f * pk) for f in (0.1, 0.9))
    return float(T[t90[0]] - T[t10[0]]) if t10.size and t90.size else np.nan


def work(task):
    pen, r, driven = task
    Y = population(pen, driven)[0]
    out = {"data": Y["all"][:, 0], "noise": np.nan_to_num((Y["odd"] - Y["even"]) / 2)[:, 0]}
    for m in MODELS:
        out[m] = full6(dynorm.predict(params(r, m), T, CD[NONBLANK], CN[NONBLANK]))[:, 0]
    return np.stack([out[s] for s in SOURCES])


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    F = fitted()
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    with ProcessPoolExecutor(n_workers) as ex:
        A = np.stack(list(ex.map(work, [(p, F.loc[p], driven) for p in F.index])))   # [pen, source, contrast, T]
    np.savez_compressed(out_root() / "binocular" / "dominant_traces.npz", A=A, pens=np.array(F.index, dtype=str),
                        sources=np.array(SOURCES), t=T, contrasts=CONTRASTS)

    k = {s: i for i, s in enumerate(SOURCES)}
    rows = []
    for ci in range(1, 6):
        c = CONTRASTS[ci]
        per = {"data": [rise(a[k["data"], ci]) for a in A]}
        for m in MODELS:
            per[m] = [rise(a[k[m], ci]) for a in A]
            per[m + "+noise"] = [rise(a[k[m], ci] + a[k["noise"], ci]) for a in A]
            per[m + "-noise"] = [rise(a[k[m], ci] - a[k["noise"], ci]) for a in A]
        for s, v in per.items():
            rows.append(dict(measure="per-penetration median", c=c, source=s, rise=np.nanmedian(v), lo=np.nan,
                             hi=np.nan))
        rng = np.random.default_rng(ci)
        boots = np.array([[rise(A[i].mean(axis=0)[k[s], ci]) for s in ("data",) + MODELS]
                          for i in (rng.integers(0, len(A), len(A)) for _ in range(2000))])
        mean = A.mean(axis=0)
        for j, s in enumerate(("data",) + MODELS):
            lo, hi = np.nanpercentile(boots[:, j], [2.5, 97.5])
            rows.append(dict(measure="population mean", c=c, source=s, rise=rise(mean[k[s], ci]), lo=lo, hi=hi))
            if s != "data":
                d = boots[:, j] - boots[:, 0]
                lo, hi = np.nanpercentile(d, [2.5, 97.5])
                rows.append(dict(measure="population mean, model - data", c=c, source=s, rise=np.nanmedian(d), lo=lo,
                                 hi=hi))
    S = pd.DataFrame(rows)
    S.to_csv(REPO / "docs" / "qc" / "rise_noise.csv", index=False)
    pd.set_option("display.width", 220)
    for meas, g in S.groupby("measure", sort=False):
        print(meas)
        print(g.pivot_table(index="source", columns="c", values="rise", sort=False).round(1).to_string())
        if g.lo.notna().any():
            print(g.assign(ci=g.lo.round(0).astype(str) + " to " + g.hi.round(0).astype(str))
                  .pivot(index="source", columns="c", values="ci").to_string())


if __name__ == "__main__":
    main()
