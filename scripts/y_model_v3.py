"""Check Y: the model by depth (N, P, K) and U5/U6 on laminar alignment v3 (docs/09, registered 2026-11-11).

Usage: .venv\\Scripts\\python scripts\\y_model_v3.py [n_workers]
Groups and channels are check X's (x_channels.parquet: first-sheet MUA channels driven by the recomputed run_crf.py
rule; superficial / input / deep from the v3 anchor). Run x_depth_v3.py first.
Y1 (N on v3): v4_model_groups.py's group populations (>= 3 channels), models 'static', 'full', 'w=a', 'cond', four
starts, held-out error per cross-validation direction; N1-N3 per group (one-sided Wilcoxon, p < .05 on the summed
error and both directions' medians in the predicted sign); N4 deep - input and superficial - input (two-sided).
Y2 (P on v3): 'org' and '2stage' on the same groups; P1-P3 as registered for P, depth differences against input.
Y3 (K on v3, no refit): the per-channel fits of v4_model_depth.py (v4_K_channels.parquet) for channels in X's set,
grouped on v3; K1', K2 per group, K3 deep - input and superficial - input.
Y4 (U on v3): U5' (A - P) deep - input > 0 (x_groups.csv), superficial - input reported two-sided; U6' w - a > 0
per v3 group (Y1 conductance fits); v4_hierarchical.py's model.
Writes <out_root>/binocular/y_N_groups.parquet, y_P_groups.parquet, y_N_pred.npz; docs/qc/y_N.csv, y_P.csv, y_K.csv,
y_U.csv.
"""

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import v4_model_depth as K_
import v4_model_groups as N_
from v4_hierarchical import row as u_row

from v1lb.config import out_root

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
V3 = ("superficial", "input", "deep")
REF, SUFFIX = "input", ""


def channel_groups() -> pd.DataFrame:
    """X's channel set: one row per (penetration, site_id) with its v3 group."""
    C = pd.read_parquet(out_root() / "binocular" / "x_channels.parquet")
    C = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip") & (C.grp != "")]
    return C[["penetration", "site_id", "grp", "rel_anchor_mm"]].drop_duplicates(["penetration", "site_id"])


def work(task):
    pen, chan_grp = task
    ds, P, keep = N_.load(pen)
    labels = [str(ch) for ch in ds.channel.values]
    chans = [k for k, ch in enumerate(labels) if ch in chan_grp]
    grp = np.array([chan_grp[labels[k]] for k in chans])
    pops, n = N_.group_means(ds, P, keep, chans, grp, V3)
    return N_.fit_groups(pen, pops, n, False), N_.fit_groups(pen, pops, n, True)


def u_rows(N: pd.DataFrame) -> pd.DataFrame:
    rows = []
    X = pd.read_csv(QC / "x_groups.csv").pivot(index="penetration", columns="grp", values="A_minus_P")
    d = (X["deep"] - X["input"]).dropna()
    rows.append(u_row("U5' (A - P) deep - input > 0 (45%, 30-100 ms)", pd.DataFrame(dict(penetration=d.index,
                                                                                         y=d.values)), 1))
    d = (X["superficial"] - X["input"]).dropna()
    r = u_row("(reported) U5' superficial - input (two-sided)", pd.DataFrame(dict(penetration=d.index, y=d.values)), 1)
    r["p_two_sided"] = 2 * min(r["p_one_sided"], 1 - r["p_one_sided"])
    rows.append(r)
    for g in V3:
        s = N[N.grp == g]
        rows.append(u_row(f"U6' excess pool weight w - a > 0, {g}",
                          pd.DataFrame(dict(penetration=s.penetration, y=s.cond_p_w - s.cond_p_a)), 1))
    R = pd.DataFrame(rows)
    R["holds"] = np.where(R.test.str.startswith("(reported)"), np.nan, R.p_one_sided < .05)
    return R


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    CG = channel_groups()
    tasks = [(p, dict(zip(g.site_id, g.grp))) for p, g in CG.groupby("penetration")]
    with ProcessPoolExecutor(n_workers) as ex:
        done = list(ex.map(work, tasks))
    nrows = [r for d, _ in done for r in d]
    prows = [r for _, d in done for r in d]
    b = out_root() / "binocular"
    N = pd.DataFrame([d[0] for d in nrows])
    N.to_parquet(b / "y_N_groups.parquet")
    np.savez_compressed(b / "y_N_pred.npz", data=np.stack([d[1] for d in nrows]),
                        **{m.replace("=", "_eq_"): np.stack([d[2][m] for d in nrows]) for m in N_.MODELS},
                        penetration=N.penetration.to_numpy(), grp=N.grp.to_numpy(), cd=N_.CD[N_.NONBLANK],
                        cn=N_.CN[N_.NONBLANK], t=N_.T)
    Pg = pd.DataFrame([d[0] for d in prows])
    Pg.to_parquet(b / "y_P_groups.parquet")
    Tn = N_.tests(N, V3, REF, SUFFIX)
    Tn.to_csv(QC / "y_N.csv", index=False)
    A = N.merge(Pg[["penetration", "grp"] + [c for c in Pg.columns if c.startswith(("cv_", "org_", "2stage_"))]],
                on=["penetration", "grp"])
    Tp = N_.tests_p(A, V3, REF, SUFFIX)
    Tp.to_csv(QC / "y_P.csv", index=False)

    Kc = pd.read_parquet(b / "v4_K_channels.parquet").drop(columns="grp").merge(
        CG[["penetration", "site_id", "grp"]], on=["penetration", "site_id"])
    Tk = K_.tests(K_.group_medians(Kc), V3, REF, SUFFIX)
    Tk.to_csv(QC / "y_K.csv", index=False)

    U = u_rows(N)
    U.to_csv(QC / "y_U.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"Y1/Y2: {N.penetration.nunique()} penetrations, {len(N)} groups: {N.grp.value_counts().to_dict()}")
    print(Tn.round(4).to_string(index=False))
    print(Tp.round(4).to_string(index=False))
    print(f"Y3: {Kc.penetration.nunique()} penetrations, {len(Kc)} channels with per-channel fits: "
          f"{Kc.grp.value_counts().to_dict()}")
    print(Tk.round(4).to_string(index=False))
    print(U.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
