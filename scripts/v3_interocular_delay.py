"""V3 analysis plan B2 and B3 (docs/09, registered 2026-10-31): is the normalization delay specifically interocular?

Usage: .venv\\Scripts\\python scripts\\v3_interocular_delay.py [n_workers]
Data, splits and starts as in fig5_cond_model.py: channel-averaged, dominant-oriented 6 x 6 time courses per
penetration (fig5_model.population), 0-219 ms, odd/even cross-validation, 4 starts.
* B2: static vs delayed ('full') normalization, each fitted and scored with the interaction loss
  (dynorm.interaction_residuals: monocular time courses and the binocular interaction R(c_D, c_N) - R(c_D, 0),
  weighted equally). r = log(CV_static / CV_full); delay needed for the interaction if one-sided Wilcoxon p < .05
  and median r > 0 in both monkeys.
* B3: model 'sep' (pool P = c_D v(tau2) + w c_N v(tau2N)) vs 'full', least squares; 'sep' gets one extra start at
  the 'full' fit of the same split with tau2N = tau2. Interocular normalization slower if median (tau2N - tau2) of
  the all-presentation fits > 0, one-sided Wilcoxon p < .05, in both monkeys; held-out error vs 'full' reported.
Writes <out_root>/binocular/v3_B2.parquet, v3_B3.parquet and docs/qc/v3_B2.csv, v3_B3.csv; prints the tests.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import glob  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig5_model import CD, CN, NONBLANK, T, population  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
STARTS = 4


def interaction_loss(p, m, cd, cn, y):
    pred = dynorm.predict(p, T, cd, cn, static=(m == "static"))
    return float(np.sum(dynorm.interaction_residuals(pred, y, cd, cn) ** 2))


def work(task):
    pen, driven = task
    t0 = time.time()
    pop = population(pen, driven)
    if pop is None:
        return None
    Y, n_ch = pop
    ys = {k: v[NONBLANK] for k, v in Y.items()}
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    b2 = dict(penetration=pen, monkey=pen[7], n_channels=n_ch)
    empty_half = any(np.isnan(ys[s]).all() for s in ("odd", "even"))
    b3_empty_half = empty_half
    for m in ("static", "full") if not empty_half else ():   # interaction loss undefined on an empty half (211008)
        f = {s: dynorm.fit(T, cd, cn, ys[s], m, starts=STARTS, loss="interaction") for s in ("odd", "even")}
        b2[f"cv_{m}"] = (interaction_loss(f["odd"].params, m, cd, cn, ys["even"])
                         + interaction_loss(f["even"].params, m, cd, cn, ys["odd"]))
        b2[f"insample_{m}"] = f["odd"].sse + f["even"].sse
        b2[f"{m}_tau2_odd"], b2[f"{m}_tau2_even"] = f["odd"].params["tau2"], f["even"].params["tau2"]

    b3 = dict(penetration=pen, monkey=pen[7], n_channels=n_ch, empty_half=b3_empty_half)
    fits = {}
    for s in ("odd", "even", "all"):
        ff = dynorm.fit(T, cd, cn, ys[s], "full", starts=STARTS)
        fs_ = dynorm.fit(T, cd, cn, ys[s], "sep", starts=STARTS, init=ff.params | {"tau2N": ff.params["tau2"]})
        fits[s] = {"full": ff, "sep": fs_}
    for m in ("full", "sep"):
        b3[f"cv_{m}"] = (dynorm.sse(fits["odd"][m].params, T, cd, cn, ys["even"], m)
                         + dynorm.sse(fits["even"][m].params, T, cd, cn, ys["odd"], m))
        b3[f"sse_{m}"] = fits["all"][m].sse
        b3 |= {f"{m}_p_{k}": v for k, v in fits["all"][m].params.items()}
    for s in ("odd", "even"):
        b3[f"sep_dtau_{s}"] = fits[s]["sep"].params["tau2N"] - fits[s]["sep"].params["tau2"]
    print(f"{pen}: B2 CV static {b2.get('cv_static', np.nan):.3f} full {b2.get('cv_full', np.nan):.3f} | "
          f"B3 CV full {b3['cv_full']:.2f} "
          f"sep {b3['cv_sep']:.2f}, tau2 {b3['sep_p_tau2']:.1f} tau2N {b3['sep_p_tau2N']:.1f} "
          f"({time.time() - t0:.0f}s)", flush=True)
    return b2, b3


def test_row(name, r, F, extra=None):
    r = pd.Series(np.asarray(r, float), index=F.index)
    p = stats.wilcoxon(r.dropna(), alternative="greater").pvalue
    mk = r.groupby(F.monkey).median()
    nb = r.groupby(F.monkey).apply(lambda x: f"{int((x > 0).sum())}/{x.notna().sum()}")
    row = dict(test=name, n=int(r.notna().sum()), median=r.median(), n_positive=int((r > 0).sum()),
               B_median=mk.get("B", np.nan), J_median=mk.get("J", np.nan), B_positive=nb.get("B", ""),
               J_positive=nb.get("J", ""), p=p, passed=bool((p < 0.05) and (mk > 0).all()))
    return row | (extra or {})


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        done = [r for r in ex.map(work, [(p, driven) for p in pens]) if r is not None]
    B2 = pd.DataFrame([d[0] for d in done]).set_index("penetration").dropna(subset=["cv_full"])
    B3 = pd.DataFrame([d[1] for d in done]).set_index("penetration")
    out = out_root() / "binocular"
    B2.to_parquet(out / "v3_B2.parquet")
    B3.to_parquet(out / "v3_B3.parquet")
    pd.set_option("display.width", 240)
    pd.set_option("display.max_colwidth", 80)

    rows2 = [test_row("B2 delayed beats static on the interaction loss: log(CV_static / CV_full)",
                      np.log(B2.cv_static / B2.cv_full), B2,
                      dict(median_pct_less_error=100 * (1 - B2.cv_full / B2.cv_static).median()))]
    rows2.append(dict(test="B2 in-sample check: full loss <= static loss (fraction of penetrations)", n=len(B2),
                      median=float((B2.insample_full <= B2.insample_static * (1 + 1e-6)).mean())))
    rows2.append(dict(test="B2 fitted tau2 of the delayed model (odd/even fits): median [IQR]", n=2 * len(B2),
                      median=float(np.median(np.r_[B2.full_tau2_odd, B2.full_tau2_even])),
                      B_median=float(np.percentile(np.r_[B2.full_tau2_odd, B2.full_tau2_even], 25)),
                      J_median=float(np.percentile(np.r_[B2.full_tau2_odd, B2.full_tau2_even], 75))))
    S2 = pd.DataFrame(rows2)

    d = B3.sep_p_tau2N - B3.sep_p_tau2
    lo, hi = dynorm.BOUNDS["tau2N"]
    rows3 = [test_row("B3 interocular pool slower: tau2N - tau2 (all-presentation 'sep' fits, ms)", d, B3),
             test_row("B3 held-out: sep beats full, log(CV_full / CV_sep) (reported, not a criterion)",
                      np.log(B3.cv_full / B3.cv_sep), B3,
                      dict(median_pct_less_error=100 * (1 - B3.cv_sep / B3.cv_full).median()))]
    keep = ~B3.empty_half.astype(bool)
    rows3.append(test_row("B3 held-out, excluding penetrations with an empty half (211008)",
                          np.log(B3.cv_full / B3.cv_sep)[keep], B3[keep]))
    for k in ("tau2", "tau2N"):
        v = B3[f"sep_p_{k}"]
        rows3.append(dict(test=f"B3 sep {k} (ms): median [IQR in B_median, J_median columns]", n=len(v),
                          median=v.median(), B_median=v.quantile(.25), J_median=v.quantile(.75)))
    rows3.append(dict(test="B3 full tau2 (ms): median [IQR]", n=len(B3), median=B3.full_p_tau2.median(),
                      B_median=B3.full_p_tau2.quantile(.25), J_median=B3.full_p_tau2.quantile(.75)))
    rows3.append(dict(test="B3 tau2N at a bound (fraction)", n=len(B3),
                      median=float(((B3.sep_p_tau2N < lo * 1.05) | (B3.sep_p_tau2N > hi * 0.95)).mean())))
    agree = np.sign(B3.sep_dtau_odd) == np.sign(B3.sep_dtau_even)
    rows3.append(dict(test="B3 sign of tau2N - tau2 agrees between odd and even fits (fraction)", n=len(B3),
                      median=float(agree.mean())))
    for mk, g in d.groupby(B3.monkey):
        rows3.append(dict(test=f"B3 tau2N - tau2, monkey {mk}: median [IQR]", n=len(g), median=g.median(),
                          B_median=g.quantile(.25), J_median=g.quantile(.75)))
    S3 = pd.DataFrame(rows3)
    print(S2.round(4).to_string(index=False))
    print(S3.round(4).to_string(index=False))
    S2.to_csv(REPO / "docs" / "qc" / "v3_B2.csv", index=False)
    S3.to_csv(REPO / "docs" / "qc" / "v3_B3.csv", index=False)


if __name__ == "__main__":
    main()
