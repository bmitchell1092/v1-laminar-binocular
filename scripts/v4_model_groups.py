"""V4 check N: the four normalization models by depth group, both cross-validation directions (docs/09, registered
2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_model_groups.py [n_workers] [extra]
With "extra" (check P, registered 2026-11-03): the recurrent ('org') and two-stage ('2stage') models on the same
groups; with N's results, the best model per group and the org-vs-cond and 2stage-vs-full comparisons by depth
(docs/qc/v4_P.csv).
Per penetration and depth group (above / at, within 2 contacts / below the earliest sink; >= 3 driven channels): the
fig5_model.py population restricted to the group's channels (channel mean, 0-219 ms after the flash, Gaussian sigma
2 ms, each channel oriented to its dominant eye and scaled by its dominant-eye peak at 90%), all / odd / even
presentations. Models 'static', 'full', 'w=a', 'cond' (cond with an extra start at the full fit and tau_m at its lower
bound, as fig5_cond_model.py); four starts each. Held-out error per direction (fit odd, score even; and the reverse);
parameters from the all-presentation fits.
N1 log(CV static / CV full) > 0; N2 log(CV w=a / CV full) > 0; N3 log(CV full / CV cond) > 0, per group (one-sided
Wilcoxon over penetrations; holds if p < .05 on the summed error and both directions have the predicted sign);
N4 depth differences in tau2, w - a, delta, tau_m and the gains (two-sided).
Writes <out_root>/binocular/v4_N_groups.parquet, v4_N_pred.npz (group data and predictions), docs/qc/v4_N.csv.
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
from scipy.ndimage import gaussian_filter1d

from v1lb import dynorm
from v1lb.binocular import CONTRASTS
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
T = np.arange(0, 220)
CD, CN = np.meshgrid(CONTRASTS, CONTRASTS, indexing="ij")
NONBLANK = (CD + CN) > 0
MODELS = ("static", "full", "w=a", "cond")
EXTRA = ("org", "2stage")                     # check P: recurrent (ORGaNICs) and two-stage
ALL_MODELS = MODELS + EXTRA
GROUPS = ("above", "at", "below")
TAU_M_MIN = dynorm.BOUNDS["tau_m"][0]
GAINS = {"N1 delay": ("static", "full"), "N2 pool weighting": ("w=a", "full"), "N3 conductance": ("full", "cond")}


def populations(pen, driven, sink):
    """Group means {grp: {split: [6, 6, T]}} and the number of channels per group."""
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    rel = ds.depth_below_toc.values[chans] - sink
    grp = np.where(np.abs(rel) <= 2, "at", np.where(rel < -2, "above", "below"))
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    out, orient = {}, None
    for split, sel in (("all", keep), ("odd", keep & (idx % 2 == 1)), ("even", keep & (idx % 2 == 0))):
        M = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= 2:
                    M[i, j] = mua[s].mean(axis=0)
        M = gaussian_filter1d(M - M[0, 0], 2, axis=2)
        if orient is None:
            rc, ri = M[-1, 0][w].mean(axis=0), M[0, -1][w].mean(axis=0)
            dom_contra = rc >= ri
            scale = np.where(dom_contra, M[-1, 0][w].max(axis=0), M[0, -1][w].max(axis=0))
            orient = (dom_contra, scale)
        dom_contra, scale = orient
        Mo = np.where(dom_contra[None, None, None, :], M, M.transpose(1, 0, 2, 3))
        for g in GROUPS:
            use = (grp == g) & (scale > 0)
            if use.sum() >= 3:
                out.setdefault(g, {})[split] = np.nanmean(Mo[..., use] / scale[use], axis=3)[:, :, tsel]
    n = {g: int(((grp == g) & (orient[1] > 0)).sum()) for g in GROUPS}
    return out, n


def work(task):
    pen, driven, sink, extra = task
    t0 = time.time()
    if not np.isfinite(sink):
        return []
    pops, n = populations(pen, driven, sink)
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    models = EXTRA if extra else MODELS
    rows = []
    for g, Y in pops.items():
        if not all(s in Y for s in ("all", "odd", "even")):
            continue
        ys = {k: v[NONBLANK] for k, v in Y.items()}
        if any(np.isnan(v).all(axis=1).any() for v in ys.values()):
            continue
        r = dict(penetration=pen, monkey=pen[7], grp=g, n_channels=n[g])
        preds = {}
        for split in ("odd", "even", "all"):
            if extra:                                              # check P: the two remaining models of the paper
                fits = {m: dynorm.fit(T, cd, cn, ys[split], m, starts=4) for m in EXTRA}
            else:
                fits = {m: dynorm.fit(T, cd, cn, ys[split], m, starts=4) for m in ("static", "full", "w=a")}
                fits["cond"] = dynorm.fit(T, cd, cn, ys[split], "cond", starts=4,
                                          init=fits["full"].params | {"tau_m": TAU_M_MIN})
            if split == "all":
                for m, f in fits.items():
                    r |= {f"{m}_p_{k}": v for k, v in f.params.items()}
                    preds[m] = dynorm.predict(f.params, T, cd, cn, static=(m == "static")).astype(np.float32)
            else:
                other = "even" if split == "odd" else "odd"
                for m, f in fits.items():
                    r[f"cv_{m}_{split[0]}>{other[0]}"] = dynorm.sse(f.params, T, cd, cn, ys[other], m)
        for m in models:
            r[f"cv_{m}"] = r[f"cv_{m}_o>e"] + r[f"cv_{m}_e>o"]
        rows.append((r, ys["all"].astype(np.float32), preds))
    print(f"{pen}: {[x[0]['grp'] for x in rows]} ({time.time() - t0:.0f}s)", flush=True)
    return rows


def tests(N: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, (worse, better) in GAINS.items():
        for g in GROUPS:
            G = N[N.grp == g]
            v = np.log(G[f"cv_{worse}"] / G[f"cv_{better}"]).dropna()
            dirs = {d: np.log(G[f"cv_{worse}_{d}"] / G[f"cv_{better}_{d}"]).dropna() for d in ("o>e", "e>o")}
            p = stats.wilcoxon(v, alternative="greater").pvalue
            same = all(x.median() > 0 for x in dirs.values())
            rows.append(dict(test=f"{name}: log(CV {worse} / CV {better}) > 0, {g} sink", n=len(v),
                             median=float(v.median()), n_in_direction=int((v > 0).sum()), p=p,
                             median_o_to_e=float(dirs["o>e"].median()), median_e_to_o=float(dirs["e>o"].median()),
                             holds=bool(p < .05 and same)))
    N = N.assign(w_minus_a=N.cond_p_w - N.cond_p_a, w_minus_a_full=N.full_p_w - N.full_p_a,
                 gain_delay=np.log(N.cv_static / N.cv_full), gain_cond=np.log(N.cv_full / N.cv_cond),
                 gain_weight=np.log(N["cv_w=a"] / N.cv_full))
    for col in ("cond_p_tau2", "cond_p_tau_m", "cond_p_delta", "w_minus_a", "full_p_tau2", "w_minus_a_full",
                "gain_delay", "gain_weight", "gain_cond"):
        piv = N.pivot(index="penetration", columns="grp", values=col)
        for g in ("above", "below"):
            d = (piv[g] - piv["at"]).dropna()
            rows.append(dict(test=f"N4 {col}, {g} - at sink", n=len(d), median=float(d.median()),
                             n_in_direction=int((d > 0).sum()), p=stats.wilcoxon(d).pvalue if len(d) > 5 else np.nan))
        for g in GROUPS:
            rows.append(dict(test=f"(reported) median {col}, {g} sink", n=int(piv[g].notna().sum()),
                             median=float(piv[g].median())))
    return pd.DataFrame(rows)


def tests_p(A: pd.DataFrame) -> pd.DataFrame:
    """Check P: the best model per group, each model against 'full', org vs cond and 2stage vs full by depth."""
    rows = []
    for g in GROUPS:
        G = A[A.grp == g]
        best = G[[f"cv_{m}" for m in ALL_MODELS]].idxmin(axis=1).str[3:]
        rows.append(dict(test=f"P1 best model (lowest held-out error), {g} sink", n=len(G),
                         **{f"best_{m}": int((best == m).sum()) for m in ALL_MODELS}))
        for m in ALL_MODELS:
            if m != "full":
                v = 100 * (1 - G[f"cv_{m}"] / G.cv_full)
                rows.append(dict(test=f"(reported) % less error than full, {m}, {g} sink", n=len(v),
                                 median=float(v.median()), n_in_direction=int((v > 0).sum())))
    for name, (worse, better) in {"P2 org vs cond": ("cond", "org"), "P3 2stage vs full": ("full", "2stage")}.items():
        for g in GROUPS:
            G = A[A.grp == g]
            v = np.log(G[f"cv_{worse}"] / G[f"cv_{better}"]).dropna()
            dirs = {d: np.log(G[f"cv_{worse}_{d}"] / G[f"cv_{better}_{d}"]).dropna() for d in ("o>e", "e>o")}
            p = stats.wilcoxon(v, alternative="greater").pvalue
            rows.append(dict(test=f"{name}: log(CV {worse} / CV {better}) > 0, {g} sink", n=len(v),
                             median=float(v.median()), n_in_direction=int((v > 0).sum()), p=p,
                             median_o_to_e=float(dirs["o>e"].median()), median_e_to_o=float(dirs["e>o"].median()),
                             holds=bool(p < .05 and all(x.median() > 0 for x in dirs.values()))))
        piv = A.assign(v=np.log(A[f"cv_{worse}"] / A[f"cv_{better}"])).pivot(index="penetration", columns="grp",
                                                                            values="v")
        for g in ("above", "below"):
            d = (piv[g] - piv["at"]).dropna()
            rows.append(dict(test=f"{name}, {g} - at sink (two-sided)", n=len(d), median=float(d.median()),
                             n_in_direction=int((d > 0).sum()), p=stats.wilcoxon(d).pvalue if len(d) > 5 else np.nan))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    sink = {p: (float(S.sink_depth_below_toc[p]) if p in S.index and S.usable[p] else np.nan) for p, _ in driven}
    pens = sorted({p for p, _ in driven})
    extra = len(sys.argv) > 2 and sys.argv[2] == "extra"
    with ProcessPoolExecutor(n_workers) as ex:
        done = [r for res in ex.map(work, [(p, driven, sink[p], extra) for p in pens]) for r in res]
    if extra:
        P = pd.DataFrame([d[0] for d in done])
        P.to_parquet(out_root() / "binocular" / "v4_P_groups.parquet")
        np.savez_compressed(out_root() / "binocular" / "v4_P_pred.npz", data=np.stack([d[1] for d in done]),
                            **{m: np.stack([d[2][m] for d in done]) for m in EXTRA},
                            penetration=P.penetration.to_numpy(), grp=P.grp.to_numpy(), cd=CD[NONBLANK],
                            cn=CN[NONBLANK], t=T)
        N = pd.read_parquet(out_root() / "binocular" / "v4_N_groups.parquet")
        A = N.merge(P[["penetration", "grp"] + [c for c in P.columns if c.startswith(("cv_", "org_", "2stage_"))]],
                    on=["penetration", "grp"])
        Tp = tests_p(A)
        Tp.to_csv(REPO / "docs" / "qc" / "v4_P.csv", index=False)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_columns", 20)
        print(f"{A.penetration.nunique()} penetrations, {len(A)} groups")
        print(Tp.round(4).to_string(index=False))
        return
    N = pd.DataFrame([d[0] for d in done])
    N.to_parquet(out_root() / "binocular" / "v4_N_groups.parquet")
    np.savez_compressed(out_root() / "binocular" / "v4_N_pred.npz", data=np.stack([d[1] for d in done]),
                        **{m.replace("=", "_eq_"): np.stack([d[2][m] for d in done]) for m in MODELS},
                        penetration=N.penetration.to_numpy(), grp=N.grp.to_numpy(), cd=CD[NONBLANK], cn=CN[NONBLANK],
                        t=T)
    Tt = tests(N)
    Tt.to_csv(REPO / "docs" / "qc" / "v4_N.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"{N.penetration.nunique()} penetrations, {len(N)} groups: {N.grp.value_counts().to_dict()}")
    print(Tt.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
