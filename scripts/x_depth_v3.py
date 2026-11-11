"""Check X: the depth results rerun on laminar alignment v3 (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\x_depth_v3.py [n_workers]
Depth relative to the v3 input/deep anchor (metadata/laminar_alignment_v3.csv, `dataset.anchor_depth`): superficial
(top of cortex to 0.5 mm above the anchor), input (0.5 mm above to the anchor), deep (the anchor to 0.7 mm below).
Channels: MUA in the first sheet, driven by the run_crf.py rule (one-way ANOVA across the six monocular contrasts,
p < .01, positive peak; either eye, either window), recomputed here for every contact. Per channel, as L, V and G:
condition means over clean presentations (>= 3 all, >= 2 per half), blank subtracted, Gaussian sigma 2 ms, dominant
eye and scale (dominant-eye peak at 90%, 40-250 ms) and eye weight a from all presentations; window means after the
channel's half-max onset of D, O, B and the summation prediction R_D(c (1 + a)); interpolated onsets (v4_latency).
X1 (A - P) deep - input > 0 (45%, 30-100 ms); X2 (AI - AI_sum) deep - input > 0; X3 dominant-eye onset at 90%
superficial - input > 0 and deep - input > 0; X4 binocular advance at 45% < 0 per group; X5 shortfall ~ depth from
the anchor, slope > 0. Sensitivity of X1, X2: deep cap 0.5 mm; the sink anchor; the 15 flip-anchored penetrations.
Writes <out_root>/binocular/x_channels.parquet, x_latency.parquet, <out_root>/figures/data/spacetime_v3.npz and
spacetime_v3_meta.parquet (onset and stimulus clocks), docs/qc/x_depth_v3.csv (tests) and docs/qc/x_groups.csv.
"""

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v4_latency import onset as interp_onset

from v1lb import crf
from v1lb.binocular import CONTRASTS, eye_weight, response_at
from v1lb.config import metadata_dir, out_root
from v1lb.dataset import ALIGNMENT_V3, anchor_depth, clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
SPLITS = ("all", "odd", "even")
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
DIOPTIC = (0.11, 0.225, 0.45)
GROUPS = ("superficial", "input", "deep")
FLOOR = 0.12                                   # V: D + O must exceed this fraction of the dominant-eye peak
REL = np.arange(-30, 251)                      # onset clock of the space-time arrays (spacetime_data.py)
T_STIM = np.arange(-20, 251)                   # stimulus clock of the space-time arrays (spacetime_data.py)
CS4 = (0.11, 0.225, 0.45, 0.9)


def driven_mask(mua_keep, t, Pk):
    """run_crf.py rule per channel: monocular contrast ANOVA p < .01 and a positive peak, either eye or window."""
    blank = np.isclose(Pk.contrast_contra, 0) & np.isclose(Pk.contrast_ipsi, 0)
    out = np.zeros(mua_keep.shape[2], bool)
    for w in crf.WINDOWS.values():
        R = crf.window_mean(mua_keep, t, w)
        R = R - np.nanmean(R[blank], axis=0)
        for cond in ("mono_contra", "mono_ipsi"):
            c, r = crf.crf_points(R, Pk, cond)
            m = crf.means_by_contrast(c, r)
            for k in range(R.shape[1]):
                groups = [r[np.isclose(c, cc), k] for cc in crf.CONTRASTS]
                if min(len(g) for g in groups) > 1:
                    out[k] |= bool(stats.f_oneway(*groups).pvalue < 0.01 and np.nanmax(m[:, k]) > 0)
    return out


def penetration(task):
    pen, al = task
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    contacts = ds.contact.values
    main = anchor_depth(contacts, al["orientation"], al["anchor_contact"], al["top_contact"])
    sink = anchor_depth(contacts, al["orientation"], al["sink_anchor_contact"], al["top_contact"])
    cand = np.flatnonzero((main.grp != "") | (sink.grp != ""))
    if not cand.size:
        return [], [], None
    mua_all = ds.mua.values[:, :, cand].astype(float)
    drv = driven_mask(mua_all[keep], t, P[keep].reset_index(drop=True))
    cand, mua = cand[drv], mua_all[:, :, drv]
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    M = {}
    for split, sel, nmin in (("all", keep, 3), ("odd", keep & (idx % 2 == 1), 2), ("even", keep & (idx % 2 == 0), 2)):
        m = np.full((6, 6, t.size, len(cand)), np.nan)
        for i, a in enumerate(CONTRASTS):
            for j, b in enumerate(CONTRASTS):
                s = sel & np.isclose(cc, a) & np.isclose(ci, b)
                if s.sum() >= nmin:
                    m[i, j] = mua[s].mean(axis=0)
        M[split] = gaussian_filter1d(m - m[0, 0], 2, axis=2)
    w = (t >= 40) & (t < 250)
    ii = [list(CONTRASTS).index(c) for c in CS4]
    rows, lat, st = [], [], dict(onset=[], stim=[], meta=[])
    for q, k in enumerate(cand):
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
        site = dict(penetration=pen, monkey=pen[7], site_id=str(ds.channel.values[k]), contact=int(contacts[k]),
                    rel_anchor_mm=float(main.rel_anchor_mm[k]), grp=main.grp[k],
                    rel_sink_mm=float(sink.rel_anchor_mm[k]), grp_sink=sink.grp[k], a=a, onset=on)
        X = np.concatenate([mo[:, 0], mo[0, ii], mo[ii, ii]])
        st["onset"].append(np.stack([np.interp(on + REL, t, x) for x in X]))
        st["stim"].append(np.stack([np.interp(T_STIM, t, x) for x in X]))
        st["meta"].append(site)
        for s in SPLITS:
            x = M[s][..., q]
            x = (x if dom_contra else x.transpose(1, 0, 2)) / pk
            i45, i90 = list(CONTRASTS).index(0.45), list(CONTRASTS).index(0.9)
            lat.append(site | dict(split=s, dom_90=interp_onset(x[i90, 0], t), dom_45=interp_onset(x[i45, 0], t),
                                   oth_45=interp_onset(x[0, i45], t),
                                   bin_45=interp_onset(x[i45, i45], t)))
            for c in DIOPTIC:
                i = list(CONTRASTS).index(c)
                pred = response_at(x[:, 0], c * (1 + a))
                for name, (w0, w1) in WINDOWS.items():
                    sel = (t >= on + w0) & (t < on + w1)
                    D, O, B, Pf = (float(np.nanmean(v[sel])) for v in (x[i, 0], x[0, i], x[i, i], pred))
                    rows.append(site | dict(split=s, c=c, window=name, D=D, O=O, B=B, P_sum=Pf,
                                            added=B - D, predicted=Pf - D))
    return rows, lat, st if st["meta"] else None


# ------------------------------------------------------------------------------------------------- tests
def wil(v, alt, name, halves=None, extra=None):
    v = pd.Series(v).dropna()
    want = {"less": -1, "greater": 1}.get(alt, np.sign(v.median()))
    p = stats.wilcoxon(v, alternative=alt).pvalue if len(v) >= 5 else np.nan
    row = dict(test=name, n=len(v), median=float(v.median()), n_in_direction=int((np.sign(v) == want).sum()), p=p)
    if halves is not None:
        row |= dict(median_odd=float(halves[0].dropna().median()), median_even=float(halves[1].dropna().median()))
        same = all(np.sign(h.dropna().median()) == want for h in halves)
    else:
        same = True
    row["holds"] = bool(p < .05 and same) if alt != "two-sided" else np.nan
    return row | (extra or {})


def shortfall_groups(C, grp="grp", split="all", cap=None, rel="rel_anchor_mm"):
    X = C[(C.split == split) & np.isclose(C.c, 0.45) & (C.window == "dip") & (C[grp] != "")]
    if cap is not None:
        X = X[X[rel] <= cap + 1e-9]
    return X.assign(S=X.added - X.predicted).groupby(["penetration", grp]).S.mean().unstack()


def ai_groups(C, grp="grp", cap=None, rel="rel_anchor_mm", window="dip"):
    X = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == window) & (C[grp] != "")]
    if cap is not None:
        X = X[X[rel] <= cap + 1e-9]
    g = X.assign(DO=X.D + X.O).groupby(["penetration", grp])[["B", "P_sum", "DO"]].mean()
    ok = g.DO > FLOOR
    return pd.DataFrame(dict(AI=np.where(ok, g.B / g.DO, np.nan), AI_sum=np.where(ok, g.P_sum / g.DO, np.nan)),
                        index=g.index)


def tests(C, L, A):
    rows = []
    S = {s: shortfall_groups(C, split=s) for s in SPLITS}
    d = {s: S[s]["deep"] - S[s]["input"] for s in SPLITS}
    rows.append(wil(d["all"], "greater", "X1 (A - P) deep - input > 0, 45%, 30-100 ms", (d["odd"], d["even"])))
    u = {s: S[s]["superficial"] - S[s]["input"] for s in SPLITS}
    rows.append(wil(u["all"], "two-sided", "(reported) X1 superficial - input", (u["odd"], u["even"])))
    for g in GROUPS:
        rows.append(wil(S["all"][g], "less", f"(reported) L1 A - P < 0, {g}", (S["odd"][g], S["even"][g])))
    G = ai_groups(C)
    gap = (G.AI - G.AI_sum).unstack()
    rows.append(wil(gap["deep"] - gap["input"], "greater", "X2 (AI - AI_sum) deep - input > 0, 45%, 30-100 ms"))
    rows.append(wil(gap["superficial"] - gap["input"], "two-sided", "(reported) X2 superficial - input"))
    for g in GROUPS:
        x = G.xs(g, level=1)
        rows.append(wil(x.AI - x.AI_sum, "less", f"(reported) AI - AI_sum < 0, {g}",
                        extra=dict(AI=float(x.AI.median()), AI_sum=float(x.AI_sum.median()))))
    Lg = {s: L[(L.split == s) & (L.grp != "")].groupby(["penetration", "grp"])[["dom_90", "adv_45"]].median()
          for s in SPLITS}
    on = {s: Lg[s].dom_90.unstack() for s in SPLITS}
    for g in ("superficial", "deep"):
        v = {s: on[s][g] - on[s]["input"] for s in SPLITS}
        rows.append(wil(v["all"], "greater", f"X3 onset dominant eye 90%, {g} - input > 0 (ms)", (v["odd"], v["even"])))
    adv = {s: Lg[s].adv_45.unstack() for s in SPLITS}
    for g in GROUPS:
        rows.append(wil(adv["all"][g], "less", f"X4 binocular advance at 45% < 0, {g} (ms)",
                        (adv["odd"][g], adv["even"][g])))
    X = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip") & (C.grp != "")]
    X = X.assign(S=X.added - X.predicted)
    fit = smf.mixedlm("S ~ rel_anchor_mm", X, groups=X.penetration).fit(reml=True)
    est, p = float(fit.params["rel_anchor_mm"]), float(fit.pvalues["rel_anchor_mm"])
    rows.append(dict(test="X5 shortfall ~ depth from the anchor (per mm), slope > 0", n=X.penetration.nunique(),
                     n_channels=len(X), median=est, p=p / 2 if est > 0 else 1 - p / 2,
                     holds=bool(est > 0 and p / 2 < .05)))
    flip = set(A.penetration[A.anchor_basis.str.startswith("flip")])
    variants = (("(a) deep cap 0.5 mm", dict(cap=0.5), None),
                ("(b) sink anchor", dict(grp="grp_sink", rel="rel_sink_mm"), None),
                ("(c) flip-anchored penetrations", {}, flip))
    for name, kw, keep in variants:
        S2 = {s: shortfall_groups(C, split=s, **kw) for s in SPLITS}
        if keep is not None:
            S2 = {s: v[v.index.isin(keep)] for s, v in S2.items()}
        d2 = {s: S2[s]["deep"] - S2[s]["input"] for s in SPLITS}
        rows.append(wil(d2["all"], "greater", f"X1 {name}", (d2["odd"], d2["even"])))
        G2 = ai_groups(C, **kw)
        if keep is not None:
            G2 = G2[G2.index.get_level_values(0).isin(keep)]
        g2 = (G2.AI - G2.AI_sum).unstack()
        rows.append(wil(g2["deep"] - g2["input"], "greater", f"X2 {name}"))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    A = pd.read_csv(metadata_dir() / ALIGNMENT_V3)
    A = A[A.depth_usable].reset_index(drop=True)
    tasks = [(r.penetration, r._asdict()) for r in A.itertuples(index=False)]
    with ProcessPoolExecutor(n_workers) as ex:
        done = list(ex.map(penetration, tasks))
    C = pd.DataFrame([r for rows, _, _ in done for r in rows])
    L = pd.DataFrame([r for _, lat, _ in done for r in lat]).assign(adv_45=lambda x: x.bin_45 - x.dom_45)
    C.to_parquet(out_root() / "binocular" / "x_channels.parquet")
    L.to_parquet(out_root() / "binocular" / "x_latency.parquet")
    st = [s for _, _, s in done if s]
    fd = out_root() / "figures" / "data"
    np.savez_compressed(fd / "spacetime_v3.npz", onset=np.concatenate([np.stack(s["onset"]) for s in st]).astype(np.float32),
                        stim=np.concatenate([np.stack(s["stim"]) for s in st]).astype(np.float32), rel=REL, t=T_STIM)
    pd.DataFrame([m for s in st for m in s["meta"]]).to_parquet(fd / "spacetime_v3_meta.parquet")

    T = tests(C, L, A)
    T.to_csv(QC / "x_depth_v3.csv", index=False)
    sf = shortfall_groups(C).stack().rename("A_minus_P")
    ai = ai_groups(C)
    on = L[(L.split == "all") & (L.grp != "")].groupby(["penetration", "grp"])[["dom_90", "adv_45"]].median()
    nch = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip") & (C.grp != "")].groupby(
        ["penetration", "grp"]).size().rename("n_channels")
    Gr = pd.concat([nch, sf, ai, on], axis=1).reset_index().rename(columns={"level_1": "grp"})
    Gr.to_csv(QC / "x_groups.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(T.round(4).to_string(index=False))
    print(Gr.groupby("grp")[["n_channels", "A_minus_P", "AI", "AI_sum", "dom_90", "adv_45"]].median().round(3))
    print("penetrations:", C.penetration.nunique(), " channels (all split, 45%, dip, first sheet):", int(nch.sum()))


if __name__ == "__main__":
    main()
