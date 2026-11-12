"""Check AE: the model by depth (Y1, Y2, AA) on the v4 coarse compartments, tiers 1 + 2 (docs/09, registered
2026-11-12).

Usage: .venv\\Scripts\\python scripts\\ae_model_v4.py [n_workers]
Channels: check X's set relabelled by ad_compartments_v4 (laminar alignment v4; scaled coarse compartments) in the
tier-1 and tier-2 penetrations. Fits as y_model_v3.py (v4_model_groups.fit_groups: 'static', 'full', 'w=a', 'cond',
then 'org' and '2stage'); the shortfall of each array as aa_model_misses_v3.deficit.
Writes <out_root>/binocular/ae_N_groups.parquet, ae_N_pred.npz, ae_P_groups.parquet and docs/qc/ae_N.csv,
ae_P.csv, ae_groups.csv, ae_model_misses.csv.
"""

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import v4_model_groups as N_
from aa_model_misses_v3 import MODELS, deficit
from ad_compartments_v4 import relabel

from v1lb import dynorm
from v1lb.config import metadata_dir, out_root

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
GROUPS = ("supragranular", "granular", "infragranular")
REF, SUFFIX = "granular", " (v4)"


def channel_groups():
    V = pd.read_csv(metadata_dir() / "laminar_alignment_v4.csv")
    V = V[V.tier.isin((1, 2))]
    C = pd.read_parquet(out_root() / "binocular" / "x_channels.parquet")
    C = relabel(C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip")], V)
    C = C[C.coarse != ""]
    return C[["penetration", "site_id", "coarse"]].drop_duplicates(["penetration", "site_id"])


def work(task):
    pen, chan_grp = task
    ds, P, keep = N_.load(pen)
    labels = [str(ch) for ch in ds.channel.values]
    chans = [k for k, ch in enumerate(labels) if ch in chan_grp]
    grp = np.array([chan_grp[labels[k]] for k in chans])
    pops, n = N_.group_means(ds, P, keep, chans, grp, GROUPS)
    return N_.fit_groups(pen, pops, n, False), N_.fit_groups(pen, pops, n, True)


def misses(Pg, Z):
    cd, cn, T = Z["cd"], Z["cn"], Z["t"]
    Pg = Pg.set_index(["penetration", "grp"])
    rows = []
    for q, (pen, g) in enumerate(zip(Z["penetration"], Z["grp"])):
        r = dict(penetration=pen, monkey=pen[7], grp=g, data=deficit(Z["data"][q], cd, cn, T))
        for name, key in MODELS.items():
            if key == "2stage":
                p = {k[len("2stage_p_"):]: v for k, v in Pg.loc[(pen, g)].items() if k.startswith("2stage_p_")}
                Y = dynorm.predict(p, T, cd, cn)
            else:
                Y = Z[key][q]
            r[name] = deficit(Y, cd, cn, T)
        rows.append(r)
    G = pd.DataFrame(rows)
    out = []
    for g in GROUPS:
        s = G[G.grp == g]
        out.append(dict(test=f"AE3 (reported) shortfall at 45%, 30-100 ms, {g}", n=len(s), data=s.data.median(),
                        **{m: s[m].median() for m in MODELS}))
    W = G.pivot(index="penetration", columns="grp")
    d_data = (W["data"]["infragranular"] - W["data"]["granular"]).dropna()
    out.append(dict(test="(reported) data infragranular - granular", n=len(d_data), median=d_data.median(),
                    n_in_direction=int((d_data > 0).sum()), p=stats.wilcoxon(d_data, alternative="greater").pvalue))
    for m in MODELS:
        d_mod = (W[m]["infragranular"] - W[m]["granular"]).dropna()
        out.append(dict(test=f"AE3 {m}: model infragranular - granular > 0", n=len(d_mod), median=d_mod.median(),
                        n_in_direction=int((d_mod > 0).sum()), p=stats.wilcoxon(d_mod, alternative="greater").pvalue))
        gap = (d_data - d_mod).dropna()
        out.append(dict(test=f"AE3 {m}: data (infra - granular) - model (infra - granular) > 0", n=len(gap),
                        median=gap.median(), n_in_direction=int((gap > 0).sum()),
                        p=stats.wilcoxon(gap, alternative="greater").pvalue))
    R = pd.DataFrame(out)
    R["holds"] = np.where(R.test.str.startswith("AE3 ") & R.test.str.endswith("> 0"), R.p < .05, np.nan)
    return G, R


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    CG = channel_groups()
    tasks = [(p, dict(zip(g.site_id, g.coarse))) for p, g in CG.groupby("penetration")]
    with ProcessPoolExecutor(n_workers) as ex:
        done = list(ex.map(work, tasks))
    b = out_root() / "binocular"
    nrows = [r for d, _ in done for r in d]
    N = pd.DataFrame([d[0] for d in nrows])
    N.to_parquet(b / "ae_N_groups.parquet")
    np.savez_compressed(b / "ae_N_pred.npz", data=np.stack([d[1] for d in nrows]),
                        **{m.replace("=", "_eq_"): np.stack([d[2][m] for d in nrows]) for m in N_.MODELS},
                        penetration=N.penetration.to_numpy(), grp=N.grp.to_numpy(), cd=N_.CD[N_.NONBLANK],
                        cn=N_.CN[N_.NONBLANK], t=N_.T)
    Pg = pd.DataFrame([d[0] for _, ds in done for d in ds])
    Pg.to_parquet(b / "ae_P_groups.parquet")
    Tn = N_.tests(N, GROUPS, REF, SUFFIX)
    Tn.to_csv(QC / "ae_N.csv", index=False)
    A = N.merge(Pg[["penetration", "grp"] + [c for c in Pg.columns if c.startswith(("cv_", "org_", "2stage_"))]],
                on=["penetration", "grp"])
    Tp = N_.tests_p(A, GROUPS, REF, SUFFIX)
    Tp.to_csv(QC / "ae_P.csv", index=False)
    G, R = misses(Pg, np.load(b / "ae_N_pred.npz", allow_pickle=True))
    G.to_csv(QC / "ae_groups.csv", index=False)
    R.to_csv(QC / "ae_model_misses.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"{N.penetration.nunique()} penetrations, {len(N)} groups: {N.grp.value_counts().to_dict()}")
    for T_ in (Tn, Tp, R):
        print(T_.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
