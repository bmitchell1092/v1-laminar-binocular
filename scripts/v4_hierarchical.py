"""V4 check U: the main claims in one hierarchical model (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_hierarchical.py
Penetration as the observation, a random intercept per session (the dual-probe days put two penetrations in one
session), monkey as a fixed sum-coded effect. U1 the 45% shortfall from summation in 30-100 ms < 0 (check F,
v4_F_windows.csv); U2 the slope of the transition-window linearity index (20-50 ms) on log contrast < 0 (V4 Fig 3b
data, normalization_timing.parquet; penetration x contrast observations, a penetration variance component); U3 the
binocular onset advance at 45% < 0 (check G); U4 the per-penetration rank correlation of onset and peak < 0 (check
H); U5 below - at the sink > 0 (check L); U6 the excess pool weight w - a > 0 in each depth group (check N,
conductance fits). Each holds if the intercept (U2: the slope) is significant one-sided at p < .05; the monkey effect
(B relative to the grand mean) is reported two-sided.
Writes docs/qc/v4_U.csv.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent / "figures"))
from fig3_contrast_depth import per_recording, transition  # noqa: E402

from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"


def model(D: pd.DataFrame, formula="y ~ 1 + C(monkey, Sum)", vc=None):
    D = D.dropna(subset=["y"]).assign(monkey=lambda d: d.penetration.str[7], session=lambda d: d.penetration.str[:8])
    return smf.mixedlm(formula, D, groups=D.session, re_formula="1", vc_formula=vc).fit(reml=True), D


def row(name, D, want, term="Intercept", formula="y ~ 1 + C(monkey, Sum)", vc=None):
    m, D = model(D, formula, vc)
    est, p = float(m.params[term]), float(m.pvalues[term])
    return dict(test=name, n_penetrations=D.penetration.nunique(), n_sessions=D.session.nunique(), estimate=est,
                p_one_sided=p / 2 if np.sign(est) == want else 1 - p / 2,
                monkey_B_effect=float(m.params["C(monkey, Sum)[S.B]"]), p_monkey=float(m.pvalues["C(monkey, Sum)[S.B]"]))


def main():
    rows = []
    F = pd.read_csv(QC / "v4_F_windows.csv")
    f = F[(F.split == "all") & np.isclose(F.c, 0.45) & (F.window == "dip")]
    rows.append(row("U1 shortfall from summation, 45%, 30-100 ms < 0",
                    pd.DataFrame(dict(penetration=f.penetration, y=f.added - f.predicted)), -1))
    tr = transition(per_recording(pd.read_parquet(out_root() / "binocular" / "normalization_timing.parquet")))
    T = tr.stack().rename("y").reset_index().rename(columns={"name": "penetration"})
    T["x"] = np.log2(T.c / 0.225)
    rows.append(row("U2 transition-window linearity index ~ log2 contrast, slope < 0", T, -1, term="x",
                    formula="y ~ x + C(monkey, Sum)", vc={"pen": "0 + C(penetration)"}))
    G = pd.read_csv(QC / "v4_G_penetrations.csv")
    G = G[G.split == "all"] if "split" in G.columns else G
    rows.append(row("U3 binocular onset advance at 45% < 0", pd.DataFrame(dict(penetration=G.penetration,
                                                                              y=G["adv_0.45"])), -1))
    C = pd.read_parquet(out_root() / "binocular" / "v4_H_cells.parquet")
    A = C[C.split == "all"]
    rho = A.groupby(["penetration", "site_id"]).apply(
        lambda g: stats.spearmanr(g.peak, g.onset, nan_policy="omit").statistic if g.onset.notna().sum() >= 6
        else np.nan).groupby(level=0).median()
    rows.append(row("U4 rank correlation of onset and peak < 0", pd.DataFrame(dict(penetration=rho.index,
                                                                                  y=rho.values)), -1))
    L = pd.read_csv(QC / "v4_L_groups.csv")
    piv = L[(L.split == "all") & np.isclose(L.c, 0.45) & (L.window == "dip")].pivot(index="penetration", columns="grp",
                                                                                      values="deficit")
    d = (piv["below"] - piv["at"]).dropna()
    rows.append(row("U5 below - at the sink > 0 (45%, 30-100 ms)", pd.DataFrame(dict(penetration=d.index,
                                                                                    y=d.values)), 1))
    N = pd.read_parquet(out_root() / "binocular" / "v4_N_groups.parquet")
    for g in ("above", "at", "below"):
        s = N[N.grp == g]
        rows.append(row(f"U6 excess pool weight w - a > 0, {g} sink",
                        pd.DataFrame(dict(penetration=s.penetration, y=s.cond_p_w - s.cond_p_a)), 1))
    R = pd.DataFrame(rows)
    R["holds"] = R.p_one_sided < .05
    R.to_csv(QC / "v4_U.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
