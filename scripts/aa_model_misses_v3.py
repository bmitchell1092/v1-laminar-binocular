"""Check AA: does the model fitted per depth group reproduce the weaker deep normalization? (docs/09, registered
2026-11-12; N's descriptive comparison, "what the model misses", rerun on laminar alignment v3).

Usage: .venv\\Scripts\\python scripts\\aa_model_misses_v3.py
Y1's group populations and fits (y_model_v3.py: <out_root>/binocular/y_N_pred.npz; two-stage predictions rebuilt from
the parameters in y_P_groups.parquet). Per penetration and v3 group, the shortfall at 45% in each eye, 30-100 ms after
the group's onset: both eyes - the dominant eye at c (1 + a), with a, the onset and the prediction computed from the
same array (data or fit), as fig4_model_spacetime.group_deficit.
AA1 (reported) medians per group; AA2 model (deep - input) > 0; AA3 data (deep - input) - model (deep - input) > 0;
one-sided Wilcoxon over penetrations; conductance (the Fig 4 model), fixed delay and two-stage fits.
Writes docs/qc/aa_model_misses.csv (tests) and docs/qc/aa_groups.csv (per penetration and group).
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from v1lb import dynorm
from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import out_root

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
MODELS = {"conductance": "cond", "fixed delay": "full", "two-stage": "2stage"}


def deficit(Y, cd, cn, T, c=0.45, w=(30, 100)):
    """Both eyes - contrast summation at c, in w ms after the array's own onset (fig4_model_spacetime.group_deficit)."""
    full = np.zeros((6, 6, T.size))
    for k, (a_, b_) in enumerate(zip(cd, cn)):
        full[list(CONTRASTS).index(a_), list(CONTRASTS).index(b_)] = Y[k]
    dom, oth = full[:, 0], full[0, :]
    j = np.flatnonzero((T >= 10) & (dom[-1] >= 0.5 * dom[-1].max()))
    on = T[j[0]]
    wt = T >= 40
    a = eye_weight(dom[:, wt].mean(axis=1), oth[:, wt].mean(axis=1))
    i = list(CONTRASTS).index(c)
    pred = response_at(dom, c * (1 + a))
    sel = (T >= on + w[0]) & (T < on + w[1])
    return float(np.mean(full[i, i, sel] - pred[sel]))


def main():
    b = out_root() / "binocular"
    Z = np.load(b / "y_N_pred.npz", allow_pickle=True)
    cd, cn, T = Z["cd"], Z["cn"], Z["t"]
    Pg = pd.read_parquet(b / "y_P_groups.parquet").set_index(["penetration", "grp"])
    rows = []
    for q, (pen, g) in enumerate(zip(Z["penetration"], Z["grp"])):
        r = dict(penetration=pen, monkey=pen[7], grp=g, data=deficit(Z["data"][q], cd, cn, T))
        for name, key in MODELS.items():
            if key == "2stage":
                p = {k[len("2stage_p_"):]: v for k, v in Pg.loc[(pen, g)].items() if k.startswith("2stage_p_")}
                Y = dynorm.predict(p, T, cd, cn)
            else:
                Y = Z[key]
                Y = Y[q]
            r[name] = deficit(Y, cd, cn, T)
        rows.append(r)
    G = pd.DataFrame(rows)
    G.to_csv(QC / "aa_groups.csv", index=False)

    out = []
    for g in ("superficial", "input", "deep"):
        s = G[G.grp == g]
        out.append(dict(test=f"AA1 (reported) shortfall at 45%, 30-100 ms, {g}", n=len(s), data=s.data.median(),
                        **{m: s[m].median() for m in MODELS}))
    W = G.pivot(index="penetration", columns="grp")
    d_data = (W["data"]["deep"] - W["data"]["input"]).dropna()
    out.append(dict(test="(reported) data deep - input", n=len(d_data), median=d_data.median(),
                    n_in_direction=int((d_data > 0).sum()), p=stats.wilcoxon(d_data, alternative="greater").pvalue))
    for m in MODELS:
        d_mod = (W[m]["deep"] - W[m]["input"]).dropna()
        out.append(dict(test=f"AA2 {m}: model deep - input > 0", n=len(d_mod), median=d_mod.median(),
                        n_in_direction=int((d_mod > 0).sum()),
                        p=stats.wilcoxon(d_mod, alternative="greater").pvalue))
        gap = (d_data - d_mod).dropna()
        out.append(dict(test=f"AA3 {m}: data (deep - input) - model (deep - input) > 0", n=len(gap),
                        median=gap.median(), n_in_direction=int((gap > 0).sum()),
                        p=stats.wilcoxon(gap, alternative="greater").pvalue))
    R = pd.DataFrame(out)
    R["holds"] = np.where(R.test.str.startswith(("AA2", "AA3")), R.p < .05, np.nan)
    R.to_csv(QC / "aa_model_misses.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"{G.penetration.nunique()} penetrations, {len(G)} groups: {G.grp.value_counts().to_dict()}")
    print(R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
