"""V4 check K: the no-delay and fixed-delay normalization models across depth (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_model_depth.py [n_workers] [cond]
With "cond" (display only, no tests): the conductance model per channel, all presentations (extra start at the full
fit, as fig5_cond_model.py), written to v4_K_cond_channels.parquet and v4_K_cond_pred.npz for the space-time maps.
Per driven channel: the fig5_model.py data for one channel (condition means for all 36 cells, Gaussian sigma 2 ms,
blank-subtracted, axis 0 = the channel's dominant eye, scaled by its dominant-eye peak at 90%; 0-219 ms after the
photodiode flash), for all, odd and even presentations (>= 2 per cell). 'static' and 'full' (v1lb.dynorm) fitted to
odd and to even presentations, each scored on the other half (held-out error = summed held-out SSE), and to all
presentations (parameters, predictions), four starts each. Per penetration and depth group (above / at / below the
earliest sink), the median over channels.
K1 log(CV static / CV full) > 0 per group; K2 w - a > 0 per group (one-sided Wilcoxon); K3 depth differences in tau2,
w - a and delta (two-sided). Each channel's half-max onset (dominant eye at 90%) is kept for response-clock maps.
Writes <out_root>/binocular/v4_K_channels.parquet, v4_K_pred.npz (data and predictions per channel),
docs/qc/v4_K.csv and docs/qc/v4_K_groups.csv.
"""

import glob
import sys
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
MODELS = ("static", "full")
GROUPS = ("above", "at", "below")
TAU_M_MIN = dynorm.BOUNDS["tau_m"][0]


def penetration(task):
    pen, driven, sink, cond_only = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans or not np.isfinite(sink):
        return []
    mua = ds.mua.values[:, :, chans].astype(float)
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    w = (t >= 40) & (t < 250)
    tsel = (t >= T[0]) & (t <= T[-1])
    M = {}
    for split, sel in (("all", keep), ("odd", keep & (idx % 2 == 1)), ("even", keep & (idx % 2 == 0))):
        m = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= 2:
                    m[i, j] = mua[s].mean(axis=0)
        M[split] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    cd, cn = CD[NONBLANK], CN[NONBLANK]
    out = []
    for q, k in enumerate(chans):
        m = M["all"][..., q]
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        pk = (m[-1, 0] if dom_contra else m[0, -1])[w].max()
        if not np.isfinite(pk) or pk <= 0:
            continue
        Y = {}
        for s in M:
            x = M[s][..., q]
            Y[s] = ((x if dom_contra else x.transpose(1, 0, 2)) / pk)[:, :, tsel][NONBLANK]
        if any(np.isnan(Y[s]).all(axis=1).any() for s in Y):
            continue
        dom90 = Y["all"][np.flatnonzero((cd == 0.9) & (cn == 0))[0]]
        j = np.flatnonzero((T >= 10) & (dom90 >= 0.5))
        if not j.size:
            continue
        rel = float(ds.depth_below_toc.values[k]) - sink
        r = dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), rel_sink=rel,
                 grp="at" if abs(rel) <= 2 else "above" if rel < -2 else "below", onset=float(T[j[0]]))
        preds = {}
        if cond_only:                                              # display mode: conductance fits, all presentations
            f_full = dynorm.fit(T, cd, cn, Y["all"], "full", starts=4)
            f = dynorm.fit(T, cd, cn, Y["all"], "cond", starts=4, init=f_full.params | {"tau_m": TAU_M_MIN})
            r |= {f"cond_{p}": v for p, v in f.params.items()}
            preds["cond"] = dynorm.predict(f.params, T, cd, cn).astype(np.float32)
            out.append((r, Y["all"].astype(np.float32), preds))
            continue
        for mdl in MODELS:
            f_odd = dynorm.fit(T, cd, cn, Y["odd"], mdl, starts=4)
            f_even = dynorm.fit(T, cd, cn, Y["even"], mdl, starts=4)
            f_all = dynorm.fit(T, cd, cn, Y["all"], mdl, starts=4)
            r[f"cv_{mdl}_o>e"] = dynorm.sse(f_odd.params, T, cd, cn, Y["even"], mdl)
            r[f"cv_{mdl}_e>o"] = dynorm.sse(f_even.params, T, cd, cn, Y["odd"], mdl)
            r[f"cv_{mdl}"] = r[f"cv_{mdl}_o>e"] + r[f"cv_{mdl}_e>o"]
            r |= {f"{mdl}_{p}": v for p, v in f_all.params.items()}
            preds[mdl] = dynorm.predict(f_all.params, T, cd, cn, static=(mdl == "static")).astype(np.float32)
        out.append((r, Y["all"].astype(np.float32), preds))
    return out


def group_medians(K: pd.DataFrame) -> pd.DataFrame:
    K = K.assign(gain=np.log(K.cv_static / K.cv_full), w_minus_a=K.full_w - K.full_a,
                 gain_oe=np.log(K["cv_static_o>e"] / K["cv_full_o>e"]),
                 gain_eo=np.log(K["cv_static_e>o"] / K["cv_full_e>o"]))
    cols = ["gain", "gain_oe", "gain_eo", "w_minus_a", "full_tau2", "full_a", "full_w", "full_delta"]
    return K.groupby(["penetration", "grp"])[cols].median().reset_index()


def tests(G: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in GROUPS:
        for col, name in (("gain", "K1 log(CV static / CV full) > 0"), ("w_minus_a", "K2 w - a > 0")):
            v = G[G.grp == g][col].dropna()
            p = stats.wilcoxon(v, alternative="greater").pvalue
            row = dict(test=f"{name}, {g} sink", n=len(v), median=float(v.median()),
                       n_in_direction=int((v > 0).sum()), p=p, holds=bool(p < .05))
            if col == "gain":                                      # K1': each cross-validation direction kept
                d = {k: G[G.grp == g][k].dropna().median() for k in ("gain_oe", "gain_eo")}
                row |= dict(median_o_to_e=float(d["gain_oe"]), median_e_to_o=float(d["gain_eo"]),
                            holds=bool(p < .05 and all(x > 0 for x in d.values())))
            rows.append(row)
    for col in ("full_tau2", "w_minus_a", "full_delta", "gain"):
        piv = G.pivot(index="penetration", columns="grp", values=col)
        for g in ("above", "below"):
            d = (piv[g] - piv["at"]).dropna()
            rows.append(dict(test=f"K3 {col}, {g} - at sink", n=len(d), median=float(d.median()),
                             n_in_direction=int((d > 0).sum()), p=stats.wilcoxon(d).pvalue, holds=np.nan))
    for g in GROUPS:
        rows.append(dict(test=f"(reported) median tau2, {g} sink", n=int(G[G.grp == g].full_tau2.notna().sum()),
                         median=float(G[G.grp == g].full_tau2.median())))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    sink = {p: (float(S.sink_depth_below_toc[p]) if p in S.index and S.usable[p] else np.nan) for p, _ in driven}
    pens = sorted({p for p, _ in driven})
    cond_only = len(sys.argv) > 2 and sys.argv[2] == "cond"
    with ProcessPoolExecutor(n_workers) as ex:
        done = [r for res in ex.map(penetration, [(p, driven, sink[p], cond_only) for p in pens]) for r in res]
    K = pd.DataFrame([d[0] for d in done])
    if cond_only:
        K.to_parquet(out_root() / "binocular" / "v4_K_cond_channels.parquet")
        np.savez_compressed(out_root() / "binocular" / "v4_K_cond_pred.npz", data=np.stack([d[1] for d in done]),
                            cond=np.stack([d[2]["cond"] for d in done]), cd=CD[NONBLANK], cn=CN[NONBLANK], t=T)
        print(f"conductance fits: {K.penetration.nunique()} penetrations, {len(K)} channels")
        return
    K.to_parquet(out_root() / "binocular" / "v4_K_channels.parquet")
    np.savez_compressed(out_root() / "binocular" / "v4_K_pred.npz", data=np.stack([d[1] for d in done]),
                        static=np.stack([d[2]["static"] for d in done]), full=np.stack([d[2]["full"] for d in done]),
                        cd=CD[NONBLANK], cn=CN[NONBLANK], t=T)
    G = group_medians(K)
    G.to_csv(REPO / "docs" / "qc" / "v4_K_groups.csv", index=False)
    Tt = tests(G)
    Tt.to_csv(REPO / "docs" / "qc" / "v4_K.csv", index=False)
    pd.set_option("display.width", 250)
    print(f"{K.penetration.nunique()} penetrations, {len(K)} channels")
    print(Tt.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
