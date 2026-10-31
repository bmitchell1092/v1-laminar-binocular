"""V3 analysis plan B4: is w > a a fitting artefact? (docs/09, registered 2026-10-31).

Usage: .venv\\Scripts\\python scripts\\v3_param_recovery.py [n_workers] [n_sets]
Per penetration, synthetic data = the prediction of its w = a fit (v3_B1_params.parquet, from
v3_null_model_index.py) plus the data's noise sample (odd - even) / 2 (v3_A1_population.npz) with an independent
random sign per condition; the 'full' model is fitted exactly as for the data (dynorm.fit, starts=4, all
presentations) and w - a recorded. Repeated for n_sets synthetic sets (default 100). Criterion: the data's median
w - a over penetrations (fig5_cond.parquet, full_p_*) exceeds the 95th percentile of the synthetic medians.
Writes <out_root>/binocular/v3_B4_fits.parquet and docs/qc/v3_B4.csv; prints the result.
"""

import sys
import warnings
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig5_model import CD, CN, NONBLANK, T  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
B = out_root() / "binocular"


def load():
    P = pd.read_parquet(B / "v3_B1_params.parquet").set_index("penetration")
    z = np.load(B / "v3_A1_population.npz")
    e = dict(zip(z["pens"], z["e"]))
    return P, e


def work(task):
    pen, s, p_wa, e = task
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    rng = np.random.default_rng([s, zlib.crc32(pen.encode())])
    sign = rng.choice([-1.0, 1.0], size=(int(NONBLANK.sum()), 1))
    y = dynorm.predict(p_wa, T, cd, cn) + sign * e[NONBLANK]
    f = dynorm.fit(T, cd, cn, y, "full", starts=4)
    return dict(penetration=pen, set=s, w=f.params["w"], a=f.params["a"], w_minus_a=f.params["w"] - f.params["a"])


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    n_sets = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    P, e = load()
    tasks = []
    for pen, r in P.iterrows():
        p_wa = {k[len("w=a_p_"):]: float(v) for k, v in r.items() if k.startswith("w=a_p_") and np.isfinite(v)}
        tasks += [(pen, s, p_wa, e[pen]) for s in range(n_sets)]
    with ProcessPoolExecutor(n_workers) as ex:
        D = pd.DataFrame(list(ex.map(work, tasks, chunksize=4)))
    D.to_parquet(B / "v3_B4_fits.parquet")

    F = pd.read_parquet(B / "fig5_cond.parquet")
    data_d = F.full_p_w - F.full_p_a
    med = D.groupby("set").w_minus_a.median()
    frac = D.assign(pos=D.w_minus_a > 0).groupby("set").pos.sum()
    p95 = float(np.percentile(med, 95))
    rows = [dict(quantity="data: median w - a over penetrations (full fit)", value=float(data_d.median()),
                 detail=f"w > a in {int((data_d > 0).sum())}/{len(data_d)}"),
            dict(quantity=f"synthetic (w = a truth): median of per-set medians of w - a ({n_sets} sets)",
                 value=float(med.median()), detail=f"5-95%: {np.percentile(med, 5):.3f} to {p95:.3f}"),
            dict(quantity="synthetic: penetrations with fitted w > a per set (median [5-95%])", value=float(frac.median()),
                 detail=f"{np.percentile(frac, 5):.0f} to {np.percentile(frac, 95):.0f} of {D.penetration.nunique()}"),
            dict(quantity="synthetic: median |w - a| per penetration fit (all fits)",
                 value=float(D.w_minus_a.abs().median()), detail=""),
            dict(quantity="CRITERION: data median w - a > 95th percentile of synthetic medians",
                 value=float(data_d.median() - p95), detail=f"passed={bool(data_d.median() > p95)}")]
    S = pd.DataFrame(rows)
    S.to_csv(REPO / "docs" / "qc" / "v3_B4.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 100)
    print(S.to_string(index=False))
    print("per-penetration synthetic median w - a vs data:")
    per = D.groupby("penetration").w_minus_a.median().rename("synthetic_median").to_frame()
    per["data"] = (F.set_index("penetration").full_p_w - F.set_index("penetration").full_p_a)
    print(per.round(3).to_string())


if __name__ == "__main__":
    main()
