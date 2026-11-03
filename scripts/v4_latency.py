"""V4 check G: onset latency by eye, contrast and depth (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_latency.py [n_workers]
Per driven in-cortex MUA channel (figure_data.py rule): condition means over clean presentations, blank subtracted,
Gaussian sigma 2 ms, for all presentations and for odd and even presentations (>= 3 presentations per condition for
all, >= 2 for a half); dominant eye (larger response at 90%, 40-250 ms), eye weight a and the scale (dominant-eye peak
at 90%) from all presentations. Onset latency = first time in 10-250 ms at which a condition reaches half its own peak,
linearly interpolated between samples; undefined if that peak is < 0.2 x the dominant-eye peak at 90%.
G1 both eyes earlier: onset(both, c) - onset(dominant, c) < 0; G2 as early as summation predicts: onset(both, c) -
onset(dominant, c (1 + a)) within +-1 ms (95% bootstrap CI of the median); G3 G1 and G2 at 45% per depth group
relative to the earliest sink (above / at, within 2 contacts / below); G4 dominant-eye onset at 90% later above and
below the sink than at it; G5 (exploratory) onset(other eye, c) - onset(dominant, a c) at 45 and 90%.
Per channel, then the median over channels per penetration (and group), then Wilcoxon over penetrations.
Writes <out_root>/binocular/v4_G_channels.parquet, docs/qc/v4_G.csv (tests) and docs/qc/v4_G_penetrations.csv.
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

from v1lb.binocular import CONTRASTS, eye_weight
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
SPLITS = ("all", "odd", "even")
MONO = (0.11, 0.225, 0.45, 0.9)
DIOPTIC = (0.11, 0.225, 0.45)
FLOOR = 0.2
MARGIN = 1.0                     # ms, equivalence margin of G2
GROUPS = ("above", "at", "below")


def onset(x, t, lo=10, hi=250, floor=FLOOR):
    """Half-max onset of x (scaled so that the dominant-eye peak at 90% = 1), interpolated between samples."""
    w = (t >= lo) & (t < hi)
    pk = np.nanmax(x[w])
    if not np.isfinite(pk) or pk < floor:
        return np.nan
    j = np.flatnonzero(w & (x >= 0.5 * pk))[0]
    if j == 0 or not w[j - 1]:
        return float(t[j])
    x0, x1 = x[j - 1], x[j]
    return float(t[j - 1] + (0.5 * pk - x0) / (x1 - x0) * (t[j] - t[j - 1]))


def onset_at(lat: dict, c: float) -> float:
    """The dominant eye's onset at contrast c, linear in log contrast between the two measured contrasts around c."""
    cs = np.array(MONO)
    if c < cs[0] or c > cs[-1]:
        return np.nan
    j = int(np.clip(np.searchsorted(cs, c), 1, len(cs) - 1))
    lo, hi = lat[cs[j - 1]], lat[cs[j]]
    f = (np.log(c) - np.log(cs[j - 1])) / (np.log(cs[j]) - np.log(cs[j - 1]))
    return float(lo * (1 - f) + hi * f)


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
        a = eye_weight((mo[:, 0] / pk)[:, w].mean(axis=1), (mo[0, :] / pk)[:, w].mean(axis=1))
        depth = float(ds.depth_below_toc.values[k])
        rel = depth - sink if np.isfinite(sink) else np.nan
        grp = (np.nan if not np.isfinite(rel) else "at" if abs(rel) <= 2 else "above" if rel < -2 else "below")
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            r = dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), split=s, a=a, depth=depth,
                     rel_sink=rel, grp=grp)
            D = {c: onset(x[list(CONTRASTS).index(c), 0], t) for c in MONO}
            for c in MONO:
                i = list(CONTRASTS).index(c)
                r[f"dom_{c:g}"] = D[c]
                r[f"oth_{c:g}"] = onset(x[0, i], t)
                r[f"bin_{c:g}"] = onset(x[i, i], t)
            for c in DIOPTIC:
                r[f"sum_{c:g}"] = onset_at(D, c * (1 + a))             # summation prediction for the binocular onset
            for c in (0.45, 0.9):
                r[f"othpred_{c:g}"] = onset_at(D, a * c)               # the other eye as a weaker dominant eye
            rows.append(r)
    return rows


def derived(C: pd.DataFrame) -> pd.DataFrame:
    C = C.copy()
    for c in MONO:
        C[f"adv_{c:g}"] = C[f"bin_{c:g}"] - C[f"dom_{c:g}"]
    for c in DIOPTIC:
        C[f"res_{c:g}"] = C[f"bin_{c:g}"] - C[f"sum_{c:g}"]
        C[f"pred_{c:g}"] = C[f"sum_{c:g}"] - C[f"dom_{c:g}"]
    for c in (0.45, 0.9):
        C[f"oxd_{c:g}"] = C[f"oth_{c:g}"] - C[f"othpred_{c:g}"]
    for c in MONO:                                                    # relative to the dominant eye at 90%
        for e in ("dom", "oth", "bin"):
            C[f"{e}rel_{c:g}"] = C[f"{e}_{c:g}"] - C["dom_0.9"]
    for c in DIOPTIC:
        C[f"sumrel_{c:g}"] = C[f"sum_{c:g}"] - C["dom_0.9"]
    return C


def boot_ci(v, n=10_000, seed=0):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    rng = np.random.default_rng(seed)
    b = np.median(v[rng.integers(0, v.size, (n, v.size))], axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def pen_medians(C: pd.DataFrame, by=("penetration",)) -> pd.DataFrame:
    cols = [c for c in C.columns if c.split("_")[0] in {"adv", "res", "pred", "oxd", "domrel", "othrel", "binrel",
                                                         "sumrel", "dom", "oth", "bin", "sum"}]
    return C.groupby(["split", *by])[cols].median().reset_index()


def one(name, P, col, alt, want, groups=None, equiv=False):
    """Test `col` over penetrations for all presentations, with the two halves."""
    out = {}
    for s in SPLITS:
        v = P[P.split == s].set_index(list(groups or ["penetration"]))[col].dropna()
        out[s] = v
    v = out["all"]
    row = dict(test=name, n=len(v), median=float(v.median()), ci_lo=np.nan, ci_hi=np.nan,
               median_odd=float(out["odd"].median()), median_even=float(out["even"].median()))
    if equiv:
        lo, hi = boot_ci(v)
        p = stats.wilcoxon(v).pvalue
        halves = all(abs(out[s].median()) <= MARGIN for s in ("odd", "even"))
        row.update(ci_lo=lo, ci_hi=hi, p=p, n_in_direction=int((v.abs() <= MARGIN).sum()),
                   holds=bool(lo >= -MARGIN and hi <= MARGIN and halves))
    else:
        lo, hi = boot_ci(v)
        p = stats.wilcoxon(v, alternative=alt).pvalue
        same = all(np.sign(out[s].median()) == want for s in ("odd", "even"))
        row.update(ci_lo=lo, ci_hi=hi, p=p, n_in_direction=int((np.sign(v) == want).sum()), holds=bool(p < .05 and same))
    return row


def tests(C: pd.DataFrame) -> pd.DataFrame:
    P = pen_medians(C)
    rows = []
    for c in DIOPTIC:
        rows.append(one(f"G1 both eyes earlier, {c:g}", P, f"adv_{c:g}", "less", -1))
        rows.append(one(f"G2 both eyes - summation prediction, {c:g}", P, f"res_{c:g}", None, 0, equiv=True))
        rows.append(one(f"(reported) predicted advance, {c:g}", P, f"pred_{c:g}", "less", -1))
    rows.append(one("(reported) both eyes earlier, 0.9", P, "adv_0.9", "less", -1))
    G = pen_medians(C[C.grp.notna()], by=("penetration", "grp"))
    for g in GROUPS:
        Gg = G[G.grp == g]
        rows.append(one(f"G3 both eyes earlier, 0.45, {g} sink", Gg, "adv_0.45", "less", -1))
        rows.append(one(f"G3 both eyes - summation prediction, 0.45, {g} sink", Gg, "res_0.45", None, 0, equiv=True))
        rows.append(one(f"(reported) predicted advance, 0.45, {g} sink", Gg, "pred_0.45", "less", -1))
    for col, lab in (("adv_0.45", "advance"), ("res_0.45", "residual"), ("dom_0.9", "G4 onset dominant 90%")):
        for g in ("above", "below"):
            piv = {s: G[G.split == s].pivot(index="penetration", columns="grp", values=col) for s in SPLITS}
            d = {s: (piv[s][g] - piv[s]["at"]).dropna() for s in SPLITS}
            v = d["all"]
            g4 = lab.startswith("G4")
            p = stats.wilcoxon(v, alternative="greater" if g4 else "two-sided").pvalue
            same = all(np.sign(d[s].median()) == np.sign(v.median()) for s in ("odd", "even"))
            if g4:
                same = all(d[s].median() > 0 for s in ("odd", "even"))
            lo, hi = boot_ci(v)
            rows.append(dict(test=f"{'G4' if g4 else 'G3 difference'} {lab}, {g} - at sink", n=len(v),
                             median=float(v.median()), ci_lo=lo, ci_hi=hi, p=p,
                             n_in_direction=int((v > 0).sum()), median_odd=float(d["odd"].median()),
                             median_even=float(d["even"].median()), holds=bool(p < .05 and same)))
    for c in (0.45, 0.9):
        r = one(f"G5 other eye extra delay, {c:g}", P, f"oxd_{c:g}", "two-sided", 1)
        r["holds"] = np.nan
        rows.append(r)
        for g in GROUPS:
            r = one(f"G5 other eye extra delay, {c:g}, {g} sink", G[G.grp == g], f"oxd_{c:g}", "two-sided", 1)
            r["holds"] = np.nan
            rows.append(r)
    return pd.DataFrame(rows)


def sink_current():
    """Reported: onset advance of the sink current and of MUA at the sink channels (sink.npz, sink-onset clock)."""
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    rel = Z["rel"].astype(float)
    out = []
    for sig in ("csd", "mua"):
        X = Z[sig]
        for c in DIOPTIC:
            i = list(CONTRASTS).index(c)
            adv = [onset(X[p, i, i], rel, lo=-5, hi=150) - onset(X[p, i, 0], rel, lo=-5, hi=150)
                   for p in range(X.shape[0])]
            v = pd.Series(adv).dropna()
            out.append(dict(test=f"(reported) sink {sig}: both eyes earlier, {c:g}", n=len(v), median=float(v.median()),
                            ci_lo=boot_ci(v)[0], ci_hi=boot_ci(v)[1], p=stats.wilcoxon(v, alternative="less").pvalue,
                            n_in_direction=int((v < 0).sum())))
    return pd.DataFrame(out)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    sink = {p: (float(S.sink_depth_below_toc[p]) if p in S.index and S.usable[p] else np.nan) for p, _ in driven}
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(n_workers) as ex:
        rows = [r for done in ex.map(penetration, [(p, driven, sink[p]) for p in pens]) for r in done]
    C = derived(pd.DataFrame(rows))
    C.to_parquet(out_root() / "binocular" / "v4_G_channels.parquet")
    T = pd.concat([tests(C), sink_current()], ignore_index=True)
    T.to_csv(REPO / "docs" / "qc" / "v4_G.csv", index=False)
    pen_medians(C[C.split == "all"]).to_csv(REPO / "docs" / "qc" / "v4_G_penetrations.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    A = C[C.split == "all"]
    print(f"{A.penetration.nunique()} penetrations, {len(A)} channels; defined: " +
          ", ".join(f"{k} {A[k].notna().mean():.0%}" for k in ("dom_0.11", "oth_0.11", "oth_0.225", "bin_0.11",
                                                               "sum_0.11", "sum_0.45", "othpred_0.45")))
    P = pen_medians(A)
    print("population medians of penetration medians, onset relative to the dominant eye at 90% (ms):")
    for e in ("dom", "oth", "bin", "sum"):
        cs = DIOPTIC if e == "sum" else MONO
        print(f"  {e}: " + ", ".join(f"{c:g} {P[f'{e}rel_{c:g}'].median():+.2f}" for c in cs))
    print("  absolute dominant-eye onset at 90%: median", round(P["dom_0.9"].median(), 2), "ms")
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
