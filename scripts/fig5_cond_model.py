"""Fifth model test: conductance-based delayed normalization ('cond') vs 'full' and 'org' (docs/09, registered 2026-10-30).

Usage: .venv\\Scripts\\python scripts\\fig5_cond_model.py [n_workers]
Data, fitting, cross-validation and timing metrics exactly as in fig5_org_model.py; 'cond' gets one extra start at the
'full' fit of the same split with tau_m at its lower bound (so it nests 'full' in practice as well as in form).
* primary: r = log(CV_full / CV_cond);
* secondary: 10-90% rise at 11% (model - data) and rise ratio 11/90%;
* tertiary: late linearity index model - data at 22.5 and 45% (11% reported, not tested);
* descriptive: 'cond' vs 'org' held-out error (fig5_org.parquet), effective time constants, latency shift 11 - 90%,
  binocular - dominant latency at 45%.
Needs fig5_org.parquet / fig5_org_L.parquet from fig5_org_model.py. Writes <out_root>/binocular/fig5_cond.parquet,
fig5_cond_L.parquet and docs/qc/fig5_cond.csv; prints the tests.
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
from fig5_org_model import full6, timing  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
LATE = (50, 75, 100)
TAU_M_MIN = dynorm.BOUNDS["tau_m"][0]


def tau_eff(p, c, after_onset):
    """Effective time constant tau_m / (1 + (P / sigma)^n) for the dominant eye alone at contrast c."""
    _, v = dynorm.step_responses(T, p["delta"], p["tau1"], p["tau2"])
    i = min(int(p["delta"] + after_onset), T.size - 1)
    return p["tau_m"] / (1 + (c * v[i] / p["sigma"]) ** p["n"])


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
    fits = {}
    for split in ("odd", "even", "all"):
        f_full = dynorm.fit(T, cd, cn, ys[split], "full", starts=4)
        f_cond = dynorm.fit(T, cd, cn, ys[split], "cond", starts=4, init=f_full.params | {"tau_m": TAU_M_MIN})
        fits[split] = {"full": f_full, "cond": f_cond}
    for m in ("full", "cond"):
        res[f"cv_{m}"] = (dynorm.sse(fits["odd"][m].params, T, cd, cn, ys["even"], m)
                          + dynorm.sse(fits["even"][m].params, T, cd, cn, ys["odd"], m))
        f_all = fits["all"][m]
        res[f"sse_{m}"] = f_all.sse
        res |= {f"{m}_p_{k}": v for k, v in f_all.params.items()}
        P = full6(dynorm.predict(f_all.params, T, cd, cn))
        res |= {f"{m}_{k}": v for k, v in timing(P).items()}
        Lrows += [dict(penetration=pen, source=m, c=c, w0=w0, L=v) for (c, w0), v in windowed_L(P).items()]
    pc = fits["all"]["cond"].params
    for c in (0.11, 0.9):
        res[f"cond_taueff_{c:g}_early"] = tau_eff(pc, c, 15)
        res[f"cond_taueff_{c:g}_late"] = tau_eff(pc, c, 120)
    print(f"{pen}: CV full {res['cv_full']:.2f} cond {res['cv_cond']:.2f} (tau_m {pc['tau_m']:.1f}); rise .11 data "
          f"{res.get('data_rise_0.11')} full {res.get('full_rise_0.11')} cond {res.get('cond_rise_0.11')} "
          f"({time.time() - t0:.0f}s)", flush=True)
    return res, Lrows


def boot_ci(x, rng, n=5000):
    boot = [np.median(rng.choice(x, len(x))) for _ in range(n)]
    return np.percentile(boot, [5, 95])


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
    F.to_parquet(out / "fig5_cond.parquet")
    Ld.to_parquet(out / "fig5_cond_L.parquet")

    O = pd.read_parquet(out / "fig5_org.parquet").set_index("penetration")
    Lo = pd.read_parquet(out / "fig5_org_L.parquet")
    F = F.set_index("penetration")
    F["cv_org"] = O.cv_org
    for k in [c for c in O.columns if c.startswith("org_lat") or c.startswith("org_rise")]:
        F[k] = O[k]
    rows = []
    rng = np.random.default_rng(0)

    def compare(test, r, passed_rule):
        p = stats.wilcoxon(r.dropna(), alternative="greater").pvalue
        mk = r.groupby(F.monkey).median()
        rows.append(dict(test=test, n=int(r.notna().sum()), median=r.median(), frac=(r > 0).mean(),
                         B=mk.get("B", np.nan), J=mk.get("J", np.nan), p=p,
                         passed=passed_rule and (p < 0.05) and (mk > 0).all()))

    compare("PRIMARY cond beats full (log CV full/cond)", np.log(F.cv_full / F.cv_cond), True)
    compare("descriptive: cond vs org (log CV org/cond)", np.log(F.cv_org / F.cv_cond), False)
    compare("descriptive: full vs org (log CV org/full; 4th test reversed)", np.log(F.cv_org / F.cv_full), False)
    for m in ("full", "org", "cond"):
        d11 = F[f"{m}_rise_0.11"] - F["data_rise_0.11"]
        ratio = F[f"{m}_rise_0.11"] / F[f"{m}_rise_0.9"]
        rows.append(dict(test=f"SECONDARY {m}: rise at 11% model - data (ms); frac = rise ratio 11/90%",
                         n=int(d11.notna().sum()), median=d11.median(), frac=ratio.median(), B=np.nan, J=np.nan,
                         p=np.nan, passed=(abs(d11.median()) <= 10) and (ratio.median() >= 1.5) if m == "cond" else np.nan))
    rows.append(dict(test="data rise 11% (ms); frac = ratio 11/90%; B = rise 90%", n=len(F),
                     median=F["data_rise_0.11"].median(), frac=(F["data_rise_0.11"] / F["data_rise_0.9"]).median(),
                     B=F["data_rise_0.9"].median(), J=np.nan, p=np.nan, passed=np.nan))
    late = (pd.concat([Ld, Lo[Lo.source == "org"]])
            .loc[lambda d: d.w0.isin(LATE)].groupby(["penetration", "source", "c"]).L.mean().unstack("source"))
    for m in ("full", "org", "cond"):
        for c in (0.11, 0.225, 0.45):
            dlt = (late.xs(c, level="c")[m] - late.xs(c, level="c")["data"]).dropna()
            lo, hi = boot_ci(dlt.values, rng)
            tested = (m == "cond") and c in (0.225, 0.45)
            rows.append(dict(test=f"{'TERTIARY' if tested else 'descriptive'} {m}: late index model - data, c={c}",
                             n=len(dlt), median=dlt.median(), frac=np.nan, B=lo, J=hi, p=np.nan,
                             passed=((lo > -0.25) and (hi < 0.25)) if tested else np.nan))
    for src in ("data", "full", "org", "cond"):
        shift = F[f"{src}_lat_0.11"] - F[f"{src}_lat_0.9"]
        bino = F[f"{src}_lat_both_0.45"] - F[f"{src}_lat_0.45"]
        rows.append(dict(test=f"{src}: latency 11% - 90% (ms); frac = binocular - dominant latency at 45%",
                         n=int(shift.notna().sum()), median=shift.median(), frac=bino.median(), B=np.nan, J=np.nan,
                         p=np.nan, passed=np.nan))
    for k in ("tau_m", "tau2", "sigma", "n", "a", "w", "delta", "tau1", "g"):
        v = F[f"cond_p_{k}"]
        rows.append(dict(test=f"cond param {k}: median [IQR]", n=len(v), median=v.median(), frac=np.nan,
                         B=v.quantile(.25), J=v.quantile(.75), p=np.nan, passed=np.nan))
    rows.append(dict(test="cond tau_m at lower bound (fraction of penetrations)", n=len(F),
                     median=np.nan, frac=(F.cond_p_tau_m < 2 * TAU_M_MIN).mean(), B=np.nan, J=np.nan, p=np.nan,
                     passed=np.nan))
    for c in (0.11, 0.9):
        for when in ("early", "late"):
            v = F[f"cond_taueff_{c:g}_{when}"]
            rows.append(dict(test=f"cond effective tau at {c:g}, {when} (ms): median [IQR]", n=len(v),
                             median=v.median(), frac=np.nan, B=v.quantile(.25), J=v.quantile(.75), p=np.nan,
                             passed=np.nan))
    S = pd.DataFrame(rows)
    pd.set_option("display.width", 240)
    pd.set_option("display.max_colwidth", 90)
    print("(late-index rows: B, J = 90% CI of the median; comparison rows: B, J = monkey medians)")
    print(S.round(3).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "fig5_cond.csv", index=False)


if __name__ == "__main__":
    main()
