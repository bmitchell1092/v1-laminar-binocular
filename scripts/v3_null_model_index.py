"""V3 analysis plan A1 (fixed-rule null) and B1 (model index time courses) (docs/09, registered 2026-10-31).

Usage: .venv\\Scripts\\python scripts\\v3_null_model_index.py [n_workers]
Per penetration, the channel-average 6 x 6 time courses of fig5_model.population (0-219 ms) and the channel-average
readout windowed_L (onset-aligned windows 0-10 ... 100-150 ms; linearity index at 11, 22.5, 45%).
* A1: the all-presentation 'full' fit (fig5_cond.parquet, full_p_*) with w set to a, so that every binocular response
  equals the monocular response at (1 + a)c at every moment with the fitted transient shapes; read out noise-free and
  with the data's noise sample (odd - even) / 2 added and subtracted. Criteria: (i) noise-free, median index over
  penetrations within 1 +- 0.25 in every window at 11, 22.5, 45%; (ii) the data's early-minus-late drop at 45%
  (early = windows starting 0, 10 ms; late = 50, 75, 100 ms) exceeds the noisy null's, paired one-sided Wilcoxon
  p < .05 for each sign of the noise; (iii) the noisy null's 20-50 ms slope against log contrast (windows starting
  20, 35 ms) is negative in fewer than 20 of 25 penetrations, for each sign of the noise.
* B1: static, w = a (fitted here, dynorm.fit, starts=4, all presentations) and delayed ('full', stored) fits read
  out in every window against the data; median over penetrations. Descriptive.
Writes <out_root>/binocular/v3_A1_L.parquet (all sources, both analyses), v3_B1_params.parquet (static and w = a
fits), v3_A1_population.npz (data and noise sample per penetration, reused by v3_param_recovery.py), and
docs/qc/v3_A1.csv, docs/qc/v3_B1.csv; prints both.
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
from fig5_org_model import full6  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
DIOPTIC = (0.11, 0.225, 0.45)
EARLY, MID, LATE = (0, 10), (20, 35), (50, 75, 100)
NULL_SOURCES = ("null", "null+noise", "null-noise")
MODEL_SOURCES = ("static", "w=a", "full")


def full_params(r):
    return {k[len("full_p_"):]: float(v) for k, v in r.items() if k.startswith("full_p_") and np.isfinite(v)}


def rows_for(pen, source, Y):
    return [dict(penetration=pen, source=source, c=c, w0=w0, L=v) for (c, w0), v in windowed_L(Y).items()]


def work(task):
    pen, r, driven = task
    t0 = time.time()
    pop = population(pen, driven)
    if pop is None:
        return None
    Y = pop[0]
    e = np.nan_to_num((Y["odd"] - Y["even"]) / 2)
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    rows = rows_for(pen, "data", Y["all"])
    p = full_params(r)
    P_full = full6(dynorm.predict(p, T, cd, cn))
    P_null = full6(dynorm.predict(p | {"w": p["a"]}, T, cd, cn))
    rows += rows_for(pen, "null", P_null) + rows_for(pen, "null+noise", P_null + e)
    rows += rows_for(pen, "null-noise", P_null - e) + rows_for(pen, "full", P_full)
    par = dict(penetration=pen, monkey=pen[7])
    for m in ("static", "w=a"):
        f = dynorm.fit(T, cd, cn, Y["all"][NONBLANK], m, starts=4)
        par |= {f"{m}_p_{k}": v for k, v in f.params.items()}
        rows += rows_for(pen, m, full6(dynorm.predict(f.params, T, cd, cn, static=(m == "static"))))
    print(f"{pen} ({time.time() - t0:.0f}s)", flush=True)
    return rows, par, Y["all"], e


def window_mean(L, w0s):
    return L[L.w0.isin(w0s)].groupby(["penetration", "source", "c"]).L.mean().unstack("source")


def slope(g):
    """Slope of the 20-50 ms index against log contrast, per penetration and source; g indexed (penetration, c)."""
    out = {}
    for pen, h in g.groupby(level="penetration"):
        h = h.droplevel("penetration").reindex(list(DIOPTIC))
        x = np.log(np.array(DIOPTIC))
        out[pen] = {s: (np.polyfit(x, h[s].to_numpy(), 1)[0] if np.isfinite(h[s]).all() else np.nan)
                    for s in h.columns}
    return pd.DataFrame(out).T


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    b = out_root() / "binocular"
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    with ProcessPoolExecutor(n_workers) as ex:
        done = [d for d in ex.map(work, [(p, F.loc[p], driven) for p in F.index]) if d is not None]
    L = pd.DataFrame([row for d in done for row in d[0]])
    P = pd.DataFrame([d[1] for d in done])
    L.to_parquet(b / "v3_A1_L.parquet")
    P.to_parquet(b / "v3_B1_params.parquet")
    np.savez_compressed(b / "v3_A1_population.npz", Y=np.stack([d[2] for d in done]),
                        e=np.stack([d[3] for d in done]), pens=np.array([d[1]["penetration"] for d in done], dtype=str))

    # ---- A1
    rows = []
    med = L.groupby(["source", "c", "w0"]).L.median().unstack("w0")
    for c in DIOPTIC:
        v = med.loc[("null", c)]
        rows.append(dict(criterion=f"A1(i) noise-free null median index per window, c={c}", n=int(L[(L.source == "null") & (L.c == c)].penetration.nunique()),
                         value="; ".join(f"{int(w)}:{x:.2f}" for w, x in v.items()), stat=float((v - 1).abs().max()),
                         p=np.nan, passed=bool(((v - 1).abs() <= 0.25).all())))
    E, Lt = window_mean(L, EARLY), window_mean(L, LATE)
    drop = (E - Lt).xs(0.45, level="c")
    for s in ("null+noise", "null-noise"):
        d = (drop["data"] - drop[s]).dropna()
        p = stats.wilcoxon(d, alternative="greater").pvalue
        rows.append(dict(criterion=f"A1(ii) early-late drop at 45%: data - {s} (median), paired one-sided Wilcoxon",
                         n=len(d), value=f"data drop {drop['data'].median():.2f}; {s} drop {drop[s].median():.2f}; "
                         f"noise-free null drop {drop['null'].median():.2f}", stat=float(d.median()), p=p,
                         passed=bool(p < 0.05)))
    S = slope(window_mean(L, MID))
    for s in ("data",) + NULL_SOURCES + MODEL_SOURCES:
        v = S[s].dropna()
        crit = s in ("null+noise", "null-noise")
        rows.append(dict(criterion=f"A1(iii) 20-50 ms slope vs log contrast negative: {s}" + ("" if crit else " (reported)"),
                         n=len(v), value=f"{int((v < 0).sum())}/{len(v)} negative; median slope {v.median():.2f}",
                         stat=float((v < 0).sum()), p=np.nan, passed=bool((v < 0).sum() < 20) if crit else np.nan))
    A1 = pd.DataFrame(rows)
    A1.to_csv(REPO / "docs" / "qc" / "v3_A1.csv", index=False)

    # ---- B1
    B1 = (L[L.source.isin(("data",) + MODEL_SOURCES + ("null",))].groupby(["c", "w0", "source"]).L
          .agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75), "count"]))
    B1.columns = ["median", "q25", "q75", "n"]
    B1 = B1.reset_index()
    B1.to_csv(REPO / "docs" / "qc" / "v3_B1.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print(A1.to_string(index=False))
    print(B1.pivot_table(index=["c", "source"], columns="w0", values="median").round(2).to_string())
    for m in ("static", "w=a"):
        print(m, "params median:", P[[c for c in P.columns if c.startswith(f"{m}_p_")]].median().round(3).to_dict())


if __name__ == "__main__":
    main()
