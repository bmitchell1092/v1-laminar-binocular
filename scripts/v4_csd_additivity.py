"""V4 check O: summation and normalization in the synaptic current at the earliest sink (docs/09, registered
2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_csd_additivity.py
Per penetration with a usable earliest sink (dichoptic_csd_sinks.csv): the sink current (-CSD, dataset CSD, over the
sink channel and its two neighbours, as sink_data.py) and MUA on the same channels; condition means over clean
presentations (>= 3 for all, >= 2 for a half), per-presentation baseline -50...0 ms removed from the CSD, blank
subtracted, Gaussian sigma 2 ms, on the sink-onset clock; all / odd / even presentations. Windows after sink onset:
early 0-20, dip 30-100, late 130-200 ms. Additivity index of window means B / (D + O) (orientation-free). Reported:
both eyes - "if the contrasts add" for the current (dominant eye = larger current at 90% in 0-50 ms, eye weight from
the current's monocular responses in 0-50 ms, scaled by the dominant-eye peak at 90% in 0-50 ms).
O1 current additivity 0-20 ms: 95% bootstrap CI of the median within 0.8-1.2 (and both halves' medians); O2 30-100 <
0-20 ms and 30-100 ms < 1; O3 MUA < current, 0-20 ms (paired, one-sided). 45% primary; 22.5 and 11% reported.
Also depth x time maps (descriptive) of both eyes - (dominant + other) for the CSD at every channel, -8...+12 contacts
from the sink, each penetration scaled by its sink's peak current, signed so that negative = less current than the sum.
Writes docs/qc/v4_O.csv, docs/qc/v4_O_penetrations.csv and <out_root>/figures/data/csd_additivity.npz.
"""

import warnings
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
REL = np.arange(-10, 201)
OFFSETS = np.arange(-8, 13)
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
DIOPTIC = (0.11, 0.225, 0.45)
SPLITS = ("all", "odd", "even")


def cell_means(x, cc, ci, sel, nmin):
    M = np.full((6, 6) + x.shape[1:], np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            s = sel & np.isclose(cc, a) & np.isclose(ci, b)
            if s.sum() >= nmin:
                M[i, j] = x[s].mean(axis=0)
    return M - M[0, 0]


def boot_ci(v, n=10_000, seed=0):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    rng = np.random.default_rng(seed)
    b = np.median(v[rng.integers(0, v.size, (n, v.size))], axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def penetration(s):
    pen = s.penetration
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    csd = ds.csd.values.astype(float)
    csd -= np.nanmean(csd[:, (t >= -50) & (t < 0)], axis=1, keepdims=True)
    k = list(ds.channel.values.astype(str)).index(s.sink_channel)
    incx = ds.in_cortex.values.astype(bool)
    chans = [q for q in (k - 1, k, k + 1) if 0 <= q < csd.shape[2] and incx[q]]
    sigs = {"current": -np.nanmean(csd[:, :, chans], axis=2),
            "mua": ds.mua.values[:, :, chans].astype(float).mean(axis=2)}
    on = float(s.sink_onset_ms)
    rows = []
    for split, sel, nmin in (("all", keep, 3), ("odd", keep & (idx % 2 == 1), 2), ("even", keep & (idx % 2 == 0), 2)):
        r = dict(penetration=pen, monkey=pen[7], split=split)
        for name, x in sigs.items():
            M = cell_means(gaussian_filter1d(x, 2, axis=1), cc, ci, sel, nmin)
            M = np.stack([[np.interp(on + REL, t, M[i, j]) for j in range(6)] for i in range(6)])
            w50 = (REL >= 0) & (REL < 50)
            Mo = M if np.nanmean(M[-1, 0, w50]) >= np.nanmean(M[0, -1, w50]) else M.transpose(1, 0, 2)
            scale = np.nanmax(Mo[-1, 0, w50])
            a = eye_weight(Mo[:, 0][:, w50].mean(axis=1) / scale, Mo[0, :][:, w50].mean(axis=1) / scale)
            for c in DIOPTIC:
                i = list(CONTRASTS).index(c)
                B, D, O = M[i, i], M[i, 0], M[0, i]
                pred = response_at(Mo[:, 0], c * (1 + a))
                for wn, (w0, w1) in WINDOWS.items():
                    sw = (REL >= w0) & (REL < w1)
                    den = np.nanmean(D[sw] + O[sw])
                    r[f"{name}_ai_{c:g}_{wn}"] = float(np.nanmean(B[sw]) / den) if den > 0 else np.nan
                    if name == "current":
                        r[f"current_deficit_{c:g}_{wn}"] = float(np.nanmean(Mo[i, i, sw] - pred[sw]) / scale)
        rows.append(r)
    # depth map, all presentations: both eyes - (dominant + other) at every contrast, per channel offset
    Mall = cell_means(gaussian_filter1d(csd, 2, axis=1), cc, ci, keep, 3)            # [6, 6, t, ch]
    strong = keep & (np.maximum(cc, ci) >= 0.45)
    peak = abs(np.nanmin(np.nanmean(csd[strong], axis=0)[(t >= 20) & (t < 70), k]))
    maps = {}
    for c in DIOPTIC:
        i = list(CONTRASTS).index(c)
        B, S_ = Mall[i, i], Mall[i, 0] + Mall[0, i]
        mB, mS = np.full((OFFSETS.size, REL.size), np.nan), np.full((OFFSETS.size, REL.size), np.nan)
        for oi, o in enumerate(OFFSETS):
            q = k + o
            if 0 <= q < csd.shape[2] and incx[q]:
                mB[oi] = np.interp(on + REL, t, -B[:, q]) / peak                    # sink positive
                mS[oi] = np.interp(on + REL, t, -S_[:, q]) / peak
        maps[c] = (mB, mS)
    return rows, maps


def test_rows(Pn: pd.DataFrame) -> pd.DataFrame:
    rows = []
    per = {s: Pn[Pn.split == s].set_index("penetration") for s in SPLITS}
    for c in DIOPTIC:
        tag = "" if c == 0.45 else "(reported) "
        e, dp = f"current_ai_{c:g}_early", f"current_ai_{c:g}_dip"
        v = per["all"][e].dropna()
        lo, hi = boot_ci(v)
        halves = all(0.8 <= per[s][e].median() <= 1.2 for s in ("odd", "even"))
        rows.append(dict(test=f"{tag}O1 current additivity index 0-20 ms within 0.8-1.2, {c:g}", n=len(v),
                         median=float(v.median()), ci_lo=lo, ci_hi=hi, median_odd=float(per["odd"][e].median()),
                         median_even=float(per["even"][e].median()),
                         holds=bool(lo >= 0.8 and hi <= 1.2 and halves)))
        for name, a, b, alt in ((f"O2 current additivity 30-100 < 0-20 ms, {c:g}", dp, e, "less"),):
            d = {s: (per[s][a] - per[s][b]).dropna() for s in SPLITS}
            p = stats.wilcoxon(d["all"], alternative=alt).pvalue
            rows.append(dict(test=f"{tag}{name}", n=len(d["all"]), median=float(d["all"].median()),
                             n_in_direction=int((d["all"] < 0).sum()), p=p, median_odd=float(d["odd"].median()),
                             median_even=float(d["even"].median()),
                             holds=bool(p < .05 and all(d[s].median() < 0 for s in ("odd", "even")))))
        x = {s: (per[s][dp] - 1).dropna() for s in SPLITS}
        p = stats.wilcoxon(x["all"], alternative="less").pvalue
        rows.append(dict(test=f"{tag}O2 current additivity 30-100 ms < 1, {c:g}", n=len(x["all"]),
                         median=float(per["all"][dp].median()), n_in_direction=int((x["all"] < 0).sum()), p=p,
                         median_odd=float(per["odd"][dp].median()), median_even=float(per["even"][dp].median()),
                         holds=bool(p < .05 and all(x[s].median() < 0 for s in ("odd", "even")))))
        for wn in ("early", "dip"):
            d = {s: (per[s][f"mua_ai_{c:g}_{wn}"] - per[s][f"current_ai_{c:g}_{wn}"]).dropna() for s in SPLITS}
            p = stats.wilcoxon(d["all"], alternative="less").pvalue
            t3 = tag if wn == "early" else "(reported) "
            rows.append(dict(test=f"{t3}O3 MUA - current additivity < 0, {wn}, {c:g}", n=len(d["all"]),
                             median=float(d["all"].median()), n_in_direction=int((d["all"] < 0).sum()), p=p,
                             median_odd=float(d["odd"].median()), median_even=float(d["even"].median()),
                             holds=bool(p < .05 and all(d[s].median() < 0 for s in ("odd", "even")))))
        for wn in WINDOWS:
            for sig in ("current", "mua"):
                v = per["all"][f"{sig}_ai_{c:g}_{wn}"].dropna()
                rows.append(dict(test=f"(reported) {sig} additivity index, {wn}, {c:g}", n=len(v),
                                 median=float(v.median()), ci_lo=boot_ci(v)[0], ci_hi=boot_ci(v)[1]))
            v = per["all"][f"current_deficit_{c:g}_{wn}"].dropna()
            rows.append(dict(test=f"(reported) current: both eyes - if the contrasts add, {wn}, {c:g}", n=len(v),
                             median=float(v.median()), n_in_direction=int((v < 0).sum()),
                             p=stats.wilcoxon(v).pvalue))
    return pd.DataFrame(rows)


def main():
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    rows, maps = [], {c: ([], []) for c in DIOPTIC}
    for _, s in S.iterrows():
        r, m = penetration(s)
        rows += r
        for c in DIOPTIC:
            maps[c][0].append(m[c][0])
            maps[c][1].append(m[c][1])
        print(s.penetration, flush=True)
    Pn = pd.DataFrame(rows)
    Pn.to_csv(REPO / "docs" / "qc" / "v4_O_penetrations.csv", index=False)
    T = test_rows(Pn)
    T.to_csv(REPO / "docs" / "qc" / "v4_O.csv", index=False)
    np.savez_compressed(out_root() / "figures" / "data" / "csd_additivity.npz", rel=REL, offsets=OFFSETS,
                        **{f"both_{c:g}": np.array(maps[c][0]) for c in DIOPTIC},
                        **{f"sum_{c:g}": np.array(maps[c][1]) for c in DIOPTIC})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
