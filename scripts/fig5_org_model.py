"""Fourth model test: binocular ORGaNICs ('org') vs delayed normalization ('full') (docs/09, registered 2026-10-30).

Usage: .venv\\Scripts\\python scripts\\fig5_org_model.py [n_workers]
Data, fitting and cross-validation exactly as for 'full' in fig5_model.py (channel-averaged, dominant-oriented 6 x 6
time courses per penetration, 0-219 ms, odd/even cross-validation, 4 starts). Per penetration:
* primary: r = log(CV_full / CV_org), with CV_full refitted here identically;
* secondary: 10-90% rise and half-max latency of the dominant eye alone at each contrast, on the data and on each
  model's all-presentation prediction (same channel-averaged time course, same metric);
* tertiary: late linearity index (windows starting 50, 75, 100 ms after onset) of data and models (windowed_L).
Writes <out_root>/binocular/fig5_org.parquet, fig5_org_L.parquet and docs/qc/fig5_org.csv; prints the tests.
"""

import glob
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig5_model import CD, CN, NONBLANK, T, population, windowed_L  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
MODELS = ["full", "org"]
LATE = (50, 75, 100)


def first_cross(x, level):
    idx = np.flatnonzero(x >= level)
    return float(T[idx[0]]) if idx.size else np.nan


def timing(Y):
    """Rise (10-90% of the peak within 150 ms) and half-max latency of the dominant eye alone; Y [6, 6, T]."""
    out = {}
    for i, c in enumerate(CONTRASTS):
        if c == 0:
            continue
        x = Y[i, 0]
        pk = np.nanmax(x[T < 150])
        if not np.isfinite(pk) or pk <= 0:
            continue
        out[f"rise_{c:g}"] = first_cross(x, 0.9 * pk) - first_cross(x, 0.1 * pk)
        out[f"lat_{c:g}"] = first_cross(x, 0.5 * pk)
    j = list(CONTRASTS).index(0.45)
    xb = Y[j, j]
    out["lat_both_0.45"] = first_cross(xb, 0.5 * np.nanmax(xb[T < 150]))
    return out


def full6(pred_nb):
    P = np.full((6, 6, T.size), np.nan)
    P[NONBLANK] = pred_nb
    P[0, 0] = 0
    return P


def work(task):
    pen, driven = task
    t0 = time.time()
    pop = population(pen, driven)
    if pop is None:
        return None
    Y, n_ch = pop
    ys = {k: v[NONBLANK] for k, v in Y.items()}
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    res = dict(penetration=pen, monkey=pen[7], n_channels=n_ch)
    res |= {f"data_{k}": v for k, v in timing(Y["all"]).items()}
    Lrows = [dict(penetration=pen, source="data", c=c, w0=w0, L=v) for (c, w0), v in windowed_L(Y["all"]).items()]
    for m in MODELS:
        f_odd = dynorm.fit(T, cd, cn, ys["odd"], m, starts=4)
        f_even = dynorm.fit(T, cd, cn, ys["even"], m, starts=4)
        res[f"cv_{m}"] = (dynorm.sse(f_odd.params, T, cd, cn, ys["even"], m)
                          + dynorm.sse(f_even.params, T, cd, cn, ys["odd"], m))
        f_all = dynorm.fit(T, cd, cn, ys["all"], m, starts=4)
        res[f"sse_{m}"] = f_all.sse
        res |= {f"{m}_p_{k}": v for k, v in f_all.params.items()}
        P = full6(dynorm.predict(f_all.params, T, cd, cn))
        res |= {f"{m}_{k}": v for k, v in timing(P).items()}
        Lrows += [dict(penetration=pen, source=m, c=c, w0=w0, L=v) for (c, w0), v in windowed_L(P).items()]
    print(f"{pen}: CV full {res['cv_full']:.2f} org {res['cv_org']:.2f}; rise .11 data {res.get('data_rise_0.11')} "
          f"full {res.get('full_rise_0.11')} org {res.get('org_rise_0.11')} ({time.time() - t0:.0f}s)", flush=True)
    return res, Lrows


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        done = [r for r in ex.map(work, [(p, driven) for p in pens]) if r is not None]
    F = pd.DataFrame([d[0] for d in done])
    Ld = pd.DataFrame([row for d in done for row in d[1]])
    out = out_root() / "binocular"
    F.to_parquet(out / "fig5_org.parquet")
    Ld.to_parquet(out / "fig5_org_L.parquet")

    rows = []
    r = np.log(F.cv_full / F.cv_org)
    p = stats.wilcoxon(r, alternative="greater").pvalue
    mk = r.groupby(F.monkey).median()
    rows.append(dict(test="PRIMARY org beats full (log CV ratio)", n=len(r), median=r.median(),
                     frac=(r > 0).mean(), B=mk.get("B", np.nan), J=mk.get("J", np.nan), p=p,
                     passed=(p < 0.05) and (mk > 0).all()))
    for m in MODELS:
        d11 = F[f"{m}_rise_0.11"] - F["data_rise_0.11"]
        ratio = F[f"{m}_rise_0.11"] / F[f"{m}_rise_0.9"]
        rows.append(dict(test=f"SECONDARY {m}: rise at 11% model - data (ms); rise ratio 11/90%", n=int(d11.notna().sum()),
                         median=d11.median(), frac=ratio.median(), B=np.nan, J=np.nan, p=np.nan,
                         passed=(abs(d11.median()) <= 10) and (ratio.median() >= 1.5)))
    rows.append(dict(test="data rise 11% / 90% (ms), ratio", n=len(F), median=F["data_rise_0.11"].median(),
                     frac=(F["data_rise_0.11"] / F["data_rise_0.9"]).median(), B=F["data_rise_0.9"].median(),
                     J=np.nan, p=np.nan, passed=np.nan))
    rng = np.random.default_rng(0)
    late = Ld[Ld.w0.isin(LATE)].groupby(["penetration", "source", "c"]).L.mean().unstack("source")
    for m in MODELS:
        for c in (0.11, 0.225, 0.45):
            g = late.xs(c, level="c")
            dlt = (g[m] - g["data"]).dropna()
            boot = [np.median(rng.choice(dlt.values, len(dlt))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            rows.append(dict(test=f"TERTIARY {m}: late linearity index model - data, c={c}", n=len(dlt),
                             median=dlt.median(), frac=np.nan, B=lo, J=hi, p=np.nan,
                             passed=(lo > -0.25) and (hi < 0.25)))
    for k in ("tau_y", "tau_a", "sigma", "b", "n", "a", "w", "delta", "tau1"):
        v = F[f"org_p_{k}"]
        rows.append(dict(test=f"org param {k}: median [IQR]", n=len(v), median=v.median(), frac=np.nan,
                         B=v.quantile(.25), J=v.quantile(.75), p=np.nan, passed=np.nan))
    lat = (F["org_lat_both_0.45"] - F["org_lat_0.45"])
    dlat = (F["data_lat_both_0.45"] - F["data_lat_0.45"])
    rows.append(dict(test="binocular minus dominant-eye latency at 45% (ms): org | data", n=len(F),
                     median=lat.median(), frac=dlat.median(), B=np.nan, J=np.nan, p=np.nan, passed=np.nan))
    S = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print("(TERTIARY rows: B, J = 90% CI of the median; rise rows: frac = ratio 11/90%)")
    print(S.round(3).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "fig5_org.csv", index=False)


if __name__ == "__main__":
    main()
