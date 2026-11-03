"""V4 check I: does the other eye suppress where it barely drives? (docs/09, registered 2026-11-03)

Usage: .venv\\Scripts\\python scripts\\v4_suppression.py [n_workers]
Channels, presentations, smoothing, dominance, scaling (dominant-eye peak at 90% = 1) and half-max onset as
v4_dip_recovery.py, onset-aligned. Added A = R(c, c) - R(c, 0) per channel in windows after onset: early 0-20, dip
30-100, late 130-200 ms. |ODI| (90%, 40-250 ms) from one half of the presentations, A from the other, both directions
(odd -> even, even -> odd), with global tercile edges per direction; per penetration and tercile the mean over
channels, the two directions averaged. I1 most monocular tercile, A < 0 in 30-100 ms (45%); I2 same tercile, A > 0 in
0-20 ms; I3 A(most binocular) - A(most monocular) > 0 in 30-100 ms; I4 I1 by depth group relative to the earliest
sink. One-sided Wilcoxon over penetrations; a test holds if p < .05 and both directions agree in sign.
Writes <out_root>/binocular/v4_I_channels.parquet, docs/qc/v4_I.csv (tests) and docs/qc/v4_I_summary.csv
(all-presentation means by tercile, window and contrast).
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

from v1lb.binocular import CONTRASTS
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
SPLITS = ("all", "odd", "even")
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
DIOPTIC = (0.11, 0.225, 0.45)
TERCILES = ("most binocular", "middle", "most monocular")
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
    if not chans:
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
        j = np.flatnonzero((t >= 10) & (t < 250) & (mo[-1, 0] / pk >= 0.5))
        if not j.size:
            continue
        on = float(t[j[0]])
        rel = float(ds.depth_below_toc.values[k]) - sink if np.isfinite(sink) else np.nan
        r = dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), rel_sink=rel,
                 grp=(np.nan if not np.isfinite(rel) else "at" if abs(rel) <= 2 else "above" if rel < -2 else "below"))
        for s in SPLITS:
            x = M[s][..., q]
            rc, ri = x[-1, 0][w].mean(), x[0, -1][w].mean()
            r[f"odi_{s}"] = abs(rc - ri) / (abs(rc) + abs(ri)) if abs(rc) + abs(ri) > 0 else np.nan
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            for c in DIOPTIC:
                i = list(CONTRASTS).index(c)
                added = x[i, i] - x[i, 0]
                for name, (w0, w1) in WINDOWS.items():
                    sel = (t >= on + w0) & (t < on + w1)
                    r[f"A_{s}_{c:g}_{name}"] = float(np.nanmean(added[sel]))
        rows.append(r)
    return rows


def terciles(v: pd.Series) -> pd.Series:
    e = v.quantile([1 / 3, 2 / 3]).to_numpy()
    return pd.Series(np.select([v <= e[0], v <= e[1]], TERCILES[:2], TERCILES[2]), index=v.index).where(v.notna())


def per_penetration(C: pd.DataFrame, c: float, win: str, by=("tercile",)) -> pd.DataFrame:
    """Mean A per penetration (and group), for each cross-validation direction and for all presentations."""
    out = {}
    for name, sel, meas in (("odd>even", "odi_odd", "even"), ("even>odd", "odi_even", "odd"), ("all", "odi_all", "all")):
        D = C.assign(tercile=terciles(C[sel]))
        out[name] = D.groupby(["penetration", *by])[f"A_{meas}_{c:g}_{win}"].mean()
    T = pd.DataFrame(out)
    T["cv"] = T[["odd>even", "even>odd"]].mean(axis=1)
    return T.reset_index()


def run_test(name, v_cv: pd.Series, v_dirs: dict, alt: str) -> dict:
    v = v_cv.dropna()
    p = stats.wilcoxon(v, alternative=alt).pvalue
    want = -1 if alt == "less" else 1
    same = all(np.sign(d.dropna().median()) == want for d in v_dirs.values())
    return dict(test=name, n=len(v), median=float(v.median()), mean=float(v.mean()),
                n_in_direction=int((np.sign(v) == want).sum()), p=p,
                median_odd_to_even=float(v_dirs["odd>even"].dropna().median()),
                median_even_to_odd=float(v_dirs["even>odd"].dropna().median()), holds=bool(p < .05 and same))


def tests(C: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in DIOPTIC:
        for win in ("dip", "early", "late"):
            T = per_penetration(C, c, win)
            mono = T[T.tercile == "most monocular"].set_index("penetration")
            bino = T[T.tercile == "most binocular"].set_index("penetration")
            tag = "" if (c == 0.45 and win in ("dip", "early")) else "(reported) "
            if win == "early":
                rows.append(run_test(f"{tag}I2 most monocular, A > 0, {win}, {c:g}", mono.cv,
                                     {k: mono[k] for k in ("odd>even", "even>odd")}, "greater"))
            else:
                rows.append(run_test(f"{tag if win == 'dip' else '(reported) '}I1 most monocular, A < 0, {win}, {c:g}",
                                     mono.cv, {k: mono[k] for k in ("odd>even", "even>odd")}, "less"))
                d = {k: (bino[k] - mono[k].reindex(bino.index)) for k in ("odd>even", "even>odd", "cv")}
                rows.append(run_test(f"{tag if win == 'dip' else '(reported) '}I3 binocular - monocular tercile, {win}, "
                                     f"{c:g}", d["cv"], {k: d[k] for k in ("odd>even", "even>odd")}, "greater"))
    G = per_penetration(C[C.grp.notna()], 0.45, "dip", by=("tercile", "grp"))
    for g in GROUPS:
        m = G[(G.tercile == "most monocular") & (G.grp == g)].set_index("penetration")
        rows.append(run_test(f"I4 most monocular, A < 0, dip, 0.45, {g} sink", m.cv,
                             {k: m[k] for k in ("odd>even", "even>odd")}, "less"))
    return pd.DataFrame(rows)


def summary(C: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in DIOPTIC:
        for win in WINDOWS:
            T = per_penetration(C, c, win)
            for tc, g in T.groupby("tercile"):
                rows.append(dict(c=c, window=win, tercile=tc, n=int(g["all"].notna().sum()),
                                 mean_all=float(g["all"].mean()), median_all=float(g["all"].median()),
                                 median_cv=float(g.cv.median())))
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
    C.to_parquet(out_root() / "binocular" / "v4_I_channels.parquet")
    T = tests(C)
    T.to_csv(REPO / "docs" / "qc" / "v4_I.csv", index=False)
    Sm = summary(C)
    Sm.to_csv(REPO / "docs" / "qc" / "v4_I_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(f"{C.penetration.nunique()} penetrations, {len(C)} channels; |ODI| tercile edges (all): "
          f"{C.odi_all.quantile([1 / 3, 2 / 3]).round(3).tolist()}")
    print(T.round(4).to_string(index=False))
    print(Sm.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
