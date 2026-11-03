"""V4 check L: is interocular normalization weaker below the input layer? (docs/09, registered 2026-11-03)

Usage: .venv\\Scripts\\python scripts\\v4_depth_norm.py [n_workers]
Per driven channel, onset-aligned as v4_dip_recovery.py (dominant eye, scale, onset and eye weight a from all
presentations), for all, odd and even presentations: added A = R(c, c) - R(c, 0) and predicted P = R_D(c (1 + a)) -
R(c, 0) (the channel's own a, v1lb.binocular.response_at) in 0-20, 30-100 and 130-200 ms after onset. Per penetration
and depth group (above / at, within 2 contacts / below the earliest sink), the mean over channels.
L1 A - P < 0 in 30-100 ms at 45% in each group; L2 (A - P) below - at > 0 (one-sided), above - at reported (two-sided).
Writes <out_root>/binocular/v4_L_channels.parquet, docs/qc/v4_L.csv (tests) and docs/qc/v4_L_groups.csv.
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

from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
SPLITS = ("all", "odd", "even")
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
DIOPTIC = (0.11, 0.225, 0.45)
GROUPS = ("above", "at", "below")


def penetration(task):
    pen, driven, sink = task
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
    M = {}
    for split, sel, nmin in (("all", keep, 3), ("odd", keep & (idx % 2 == 1), 2), ("even", keep & (idx % 2 == 0), 2)):
        m = np.full((6, 6, t.size, len(chans)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= nmin:
                    m[i, j] = mua[s].mean(axis=0)
        M[split] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    rows = []
    for q, k in enumerate(chans):
        m = M["all"][..., q]
        dom_contra = m[-1, 0][w].mean() >= m[0, -1][w].mean()
        mo = m if dom_contra else m.transpose(1, 0, 2)
        pk = np.nanmax(mo[-1, 0][w])
        if not np.isfinite(pk) or pk <= 0 or np.isnan(mo[:, 0]).any() or np.isnan(mo[0, :]).any():
            continue
        mo = mo / pk
        j = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] >= 0.5))
        if not j.size:
            continue
        on = float(t[j[0]])
        a = eye_weight(mo[:, 0][:, w].mean(axis=1), mo[0, :][:, w].mean(axis=1))
        rel = float(ds.depth_below_toc.values[k]) - sink
        grp = "at" if abs(rel) <= 2 else "above" if rel < -2 else "below"
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            for c in DIOPTIC:
                i = list(CONTRASTS).index(c)
                added = x[i, i] - x[i, 0]
                pred = response_at(x[:, 0], c * (1 + a)) - x[i, 0]
                for name, (w0, w1) in WINDOWS.items():
                    sel = (t >= on + w0) & (t < on + w1)
                    rows.append(dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), split=s,
                                     rel_sink=rel, grp=grp, a=a, c=c, window=name,
                                     added=float(np.nanmean(added[sel])), predicted=float(np.nanmean(pred[sel]))))
    return rows


def groups(C: pd.DataFrame) -> pd.DataFrame:
    C = C.assign(deficit=C.added - C.predicted)
    return C.groupby(["split", "c", "window", "penetration", "grp"])[["added", "predicted", "deficit"]].mean().reset_index()


def tests(G: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in DIOPTIC:
        for win in WINDOWS:
            tag = "" if (c == 0.45 and win == "dip") else "(reported) "
            piv = {s: G[(G.split == s) & np.isclose(G.c, c) & (G.window == win)].pivot(
                index="penetration", columns="grp", values="deficit") for s in SPLITS}
            for g in GROUPS:
                v = {s: piv[s][g].dropna() for s in SPLITS}
                p = stats.wilcoxon(v["all"], alternative="less").pvalue
                same = all(v[s].median() < 0 for s in ("odd", "even"))
                rows.append(dict(test=f"{tag}L1 both - prediction < 0, {g} sink, {win}, {c:g}", n=len(v["all"]),
                                 median=float(v["all"].median()), n_in_direction=int((v["all"] < 0).sum()), p=p,
                                 median_odd=float(v["odd"].median()), median_even=float(v["even"].median()),
                                 holds=bool(p < .05 and same)))
            for g, alt in (("below", "greater"), ("above", "two-sided")):
                d = {s: (piv[s][g] - piv[s]["at"]).dropna() for s in SPLITS}
                p = stats.wilcoxon(d["all"], alternative=alt).pvalue
                want = 1 if alt == "greater" else np.sign(d["all"].median())
                same = all(np.sign(d[s].median()) == want for s in ("odd", "even"))
                rows.append(dict(test=f"{tag}L2 {g} - at sink, {win}, {c:g}", n=len(d["all"]),
                                 median=float(d["all"].median()), n_in_direction=int((d["all"] > 0).sum()), p=p,
                                 median_odd=float(d["odd"].median()), median_even=float(d["even"].median()),
                                 holds=bool(p < .05 and same)))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    sink = {p: (float(S.sink_depth_below_toc[p]) if p in S.index and S.usable[p] else np.nan) for p, _ in driven}
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        C = pd.DataFrame([r for done in ex.map(penetration, [(p, driven, sink[p]) for p in pens]) for r in done])
    C.to_parquet(out_root() / "binocular" / "v4_L_channels.parquet")
    G = groups(C)
    G.to_csv(REPO / "docs" / "qc" / "v4_L_groups.csv", index=False)
    T = tests(G)
    T.to_csv(REPO / "docs" / "qc" / "v4_L.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(T.round(4).to_string(index=False))
    A = G[G.split == "all"].groupby(["c", "window", "grp"])[["added", "predicted", "deficit"]].median().round(3)
    print(A.to_string())


if __name__ == "__main__":
    main()
