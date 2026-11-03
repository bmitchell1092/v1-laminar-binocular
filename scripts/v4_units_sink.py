"""V4 check M: do single units at the sink follow the sink current while MUA leads it? (docs/09, registered
2026-11-03)

Usage: .venv\\Scripts\\python scripts\\v4_units_sink.py [n_workers]
Units: lenient set (well isolated, or refractory contamination < 0.5 and SNR >= 3), in cortex, responding to the
strong-drive pool (peak of the smoothed rate in 10-150 ms > baseline mean + 5 SD, -50...0 ms). Presentations: clean,
either eye >= 45%, blank presentations subtracted (reported: the unit's dominant eye at 90%, both eyes at 45%).
Sink contact chosen on odd presentations by the dichoptic_csd.py rule and every signal measured on even ones, and the
reverse. Signals (1 ms, Gaussian sigma 2 ms): unit rate, analog MUA on the unit's peak contact, sink current (-CSD)
at the chosen contact from unsmoothed LFP, advanced by the 0.9 ms LFP filter delay (reported: the three-contact,
Hamming-smoothed current of check J). Half-max of each signal's own peak (10-150 ms; primary), 25% crossing and peak
time (reported). Units within 2 contacts of the chosen contact are at the sink.
M1 unit - current > 0 at the sink; M2 unit - MUA > 0 at the sink; M3 MUA - current (reported); M4 (unit - MUA) at the
sink > elsewhere (Mann-Whitney). Units are the sample; a test holds if p < .05 and both directions agree in sign.
Writes <out_root>/binocular/v4_M_units.parquet, docs/qc/v4_M.csv and docs/qc/v4_M_units.csv.
"""

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v4_latency import onset

from v1lb import laminar as lam
from v1lb.config import out_root
from v1lb.dataset import bin_spikes, clean_presentations, open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
HW_DELAY = 0.9
DIRECTIONS = (("odd", "even"), ("even", "odd"))


def earliest_sink(csd_mean, t, incx):
    """The dichoptic_csd.py rule on a mean CSD [time, channel] (negative = sink): among channels whose 20-70 ms sink
    is >= 40% of the strongest, the one whose sink reaches 25% of its own peak first (15-70 ms)."""
    w = (t >= 20) & (t < 70)
    mm = np.where(incx[None, :], csd_mean, np.nan)[w]
    peaks = np.nanmin(mm, axis=0)
    strong = np.flatnonzero(incx & (peaks <= 0.4 * np.nanmin(peaks)))
    onsets = {}
    for k in strong:
        idx = np.flatnonzero((t >= 15) & (t < 70) & (-csd_mean[:, k] >= 0.25 * abs(peaks[k])))
        if idx.size:
            onsets[k] = float(t[idx[0]])
    return min(onsets, key=lambda k: (onsets[k], peaks[k])) if onsets else int(np.nanargmin(peaks))


def cross(x, t, frac, lo=10, hi=150):
    """First crossing of frac x own peak, interpolated (frac = 0.5 is the half-max onset)."""
    w = (t >= lo) & (t < hi)
    pk = np.nanmax(x[w])
    if not np.isfinite(pk) or pk <= 0:
        return np.nan
    j = np.flatnonzero(w & (x >= frac * pk))[0]
    if j == 0 or not w[j - 1]:
        return float(t[j])
    return float(t[j - 1] + (frac * pk - x[j - 1]) / (x[j] - x[j - 1]))


def timing(x, t):
    w = (t >= 10) & (t < 100)
    if not np.isfinite(x[w]).any():                                # edge contacts have no CSD
        return dict(half=np.nan, q25=np.nan, peak=np.nan)
    return dict(half=onset(x / np.nanmax(x[(t >= 10) & (t < 150)]), t, lo=10, hi=150, floor=0),
                q25=cross(x, t, 0.25), peak=float(t[w][np.nanargmax(x[w])]))


def work(task):
    pen = task
    upath = out_root() / "units" / f"{pen}_units.parquet"
    if not upath.exists():
        return []
    U = pd.read_parquet(upath)
    U = U[U.in_cortex.fillna(False).astype(bool)]
    U = U[U.well_isolated.fillna(False).astype(bool) | ((U.isi_contamination < 0.5) & (U.snr >= 3))].reset_index(drop=True)
    if U.empty:
        return []
    Sp = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values.astype(float)
    labels = list(ds.channel.values.astype(str))
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    idx = np.arange(len(P))
    half = {"odd": idx % 2 == 1, "even": idx % 2 == 0}
    blank = keep & (cc == 0) & (ci == 0)
    conds = {"strong pool": keep & (np.maximum(cc, ci) >= 0.45),
             "both eyes 45%": keep & np.isclose(cc, .45) & np.isclose(ci, .45),
             "contra 90%": keep & np.isclose(cc, .9) & (ci == 0), "ipsi 90%": keep & (cc == 0) & np.isclose(ci, .9)}
    win = (int(t[0]), int(t[-1]) + 1)
    ticks = [Sp.loc[Sp.cluster_id == cid, "tick"].to_numpy() for cid in U.cluster_id]
    rates = bin_spikes(ticks, P.pd_on_tick.to_numpy(), window_ms=win).astype(float) * 1000.0     # [pres, t, unit]
    lfp, csd_h, mua = (ds[v].values.astype(float) for v in ("lfp", "csd", "mua"))
    bad = ds.bad_lfp.values
    incx = ds.in_cortex.values.astype(bool)
    base = (t >= -50) & (t < 0)
    # response criterion on the strong pool, all presentations
    r_all = gaussian_filter1d(rates[conds["strong pool"]].mean(axis=0), 2, axis=0)
    resp = [np.nanmax(r_all[(t >= 10) & (t < 150), q]) > r_all[base, q].mean() + 5 * r_all[base, q].std()
            for q in range(len(U))]
    # each unit's dominant eye at 90% (all presentations, 40-250 ms)
    w = (t >= 40) & (t < 250)
    dom = ["contra 90%" if rates[conds["contra 90%"]][:, w, q].mean() >= rates[conds["ipsi 90%"]][:, w, q].mean()
           else "ipsi 90%" for q in range(len(U))]
    rows = []
    for X, Y in DIRECTIONS:
        sx = conds["strong pool"] & half[X]
        cm = csd_h[sx].mean(axis=0)
        cm -= cm[base].mean(axis=0)
        ch = earliest_sink(cm, t, incx & np.isfinite(cm).all(axis=0))
        sink_depth = float(ds.depth_below_toc.values[ch])
        for cname in ("strong pool", "both eyes 45%", "dominant eye 90%"):
            for q, u in U.iterrows():
                if not resp[q]:
                    continue
                sel = conds[dom[q] if cname == "dominant eye 90%" else cname] & half[Y]
                by = blank & half[Y]
                if sel.sum() < 5 or by.sum() < 3:
                    continue
                L = lfp[sel].mean(axis=0) - lfp[by].mean(axis=0)
                cur = gaussian_filter1d(-lam.csd(L, smooth=False, bad=bad)[:, ch], 2)
                nb = [k for k in (ch - 1, ch, ch + 1) if 0 <= k < len(labels) and incx[k]]
                cur3 = gaussian_filter1d(-(csd_h[sel].mean(axis=0) - csd_h[by].mean(axis=0))[:, nb].mean(axis=1), 2)
                pc = labels.index(str(u.peak_label))
                m_ = gaussian_filter1d(mua[sel][:, :, pc].mean(axis=0) - mua[by][:, :, pc].mean(axis=0), 2)
                r_ = gaussian_filter1d(rates[sel][:, :, q].mean(axis=0) - rates[by][:, :, q].mean(axis=0), 2)
                T = {"unit": timing(r_, t), "mua": timing(m_, t), "current": timing(cur, t),
                     "current3": timing(cur3, t)}
                for k in ("current", "current3"):
                    T[k] = {m: v - HW_DELAY for m, v in T[k].items()}
                rel = float(u.depth_below_toc) - sink_depth
                rows.append(dict(penetration=pen, monkey=pen[7], cluster_id=int(u.cluster_id),
                                 strict=bool(u.well_isolated) if pd.notna(u.well_isolated) else False,
                                 direction=f"{X}>{Y}", condition=cname, sink_channel=labels[ch], rel_sink=rel,
                                 grp="at" if abs(rel) <= 2 else "above" if rel < -2 else "below",
                                 rate_hz=float(u.rate_hz),
                                 **{f"{s}_{m}": v for s, d_ in T.items() for m, v in d_.items()}))
    return rows


def differences(M: pd.DataFrame) -> pd.DataFrame:
    M = M.copy()
    for m in ("half", "q25", "peak"):
        M[f"unit_minus_current_{m}"] = M[f"unit_{m}"] - M[f"current_{m}"]
        M[f"unit_minus_mua_{m}"] = M[f"unit_{m}"] - M[f"mua_{m}"]
        M[f"mua_minus_current_{m}"] = M[f"mua_{m}"] - M[f"current_{m}"]
        M[f"unit_minus_current3_{m}"] = M[f"unit_{m}"] - M[f"current3_{m}"]
    return M


def per_unit(M: pd.DataFrame, cond: str, grp_sel) -> dict:
    """Per-direction and combined (mean over directions) values per unit, restricted by group."""
    g = M[(M.condition == cond) & grp_sel(M)]
    cols = [c for c in g.columns if "_minus_" in c]
    out = {dname: g[g.direction == dname].set_index(["penetration", "cluster_id"])[cols] for dname in ("odd>even",
                                                                                                         "even>odd")}
    out["combined"] = g.groupby(["penetration", "cluster_id"])[cols].mean()
    return out


def one_sided(name, D, col, alt="greater"):
    v = D["combined"][col].dropna()
    if len(v) < 3:
        return dict(test=name, n=len(v))
    p = stats.wilcoxon(v, alternative=alt).pvalue
    want = 1 if alt == "greater" else -1
    same = all(np.sign(D[k][col].dropna().median()) == want for k in ("odd>even", "even>odd"))
    pens = D["combined"][col].dropna().groupby(level=0).median()
    return dict(test=name, n=len(v), n_penetrations=int(pens.size), median=float(v.median()),
                n_in_direction=int((np.sign(v) == want).sum()), p=p,
                median_odd_to_even=float(D["odd>even"][col].dropna().median()),
                median_even_to_odd=float(D["even>odd"][col].dropna().median()),
                penetration_median=float(pens.median()), n_pen_in_direction=int((np.sign(pens) == want).sum()),
                holds=bool(p < .05 and same and len(v) >= 10), inconclusive=bool(len(v) < 10))


def tests(M: pd.DataFrame) -> pd.DataFrame:
    rows = []
    at = lambda m: m.grp == "at"                                   # noqa: E731
    away = lambda m: m.grp != "at"                                 # noqa: E731
    for cond in ("strong pool", "both eyes 45%", "dominant eye 90%"):
        tag = "" if cond == "strong pool" else "(reported) "
        for subset in ("lenient", "strict"):
            D = per_unit(M[M.strict] if subset == "strict" else M, cond, at)
            t2 = tag if subset == "lenient" else "(reported) "
            rows.append(one_sided(f"{t2}M1 unit - current > 0, half-max, at sink, {cond}, {subset}", D,
                                  "unit_minus_current_half"))
            rows.append(one_sided(f"{t2}M2 unit - MUA > 0, half-max, at sink, {cond}, {subset}", D,
                                  "unit_minus_mua_half"))
            rows.append(one_sided(f"(reported) M3 MUA - current < 0, half-max, at sink, {cond}, {subset}", D,
                                  "mua_minus_current_half", "less"))
            for m in ("q25", "peak"):
                rows.append(one_sided(f"(reported) M1 unit - current > 0, {m}, at sink, {cond}, {subset}", D,
                                      f"unit_minus_current_{m}"))
                rows.append(one_sided(f"(reported) M2 unit - MUA > 0, {m}, at sink, {cond}, {subset}", D,
                                      f"unit_minus_mua_{m}"))
            rows.append(one_sided(f"(reported) M1 unit - three-contact current > 0, half-max, at sink, {cond}, "
                                  f"{subset}", D, "unit_minus_current3_half"))
        Da, Dw = per_unit(M, cond, at), per_unit(M, cond, away)
        a, b = Da["combined"]["unit_minus_mua_half"].dropna(), Dw["combined"]["unit_minus_mua_half"].dropna()
        p = stats.mannwhitneyu(a, b, alternative="greater").pvalue if len(a) and len(b) else np.nan
        same = all((Da[k]["unit_minus_mua_half"].median() - Dw[k]["unit_minus_mua_half"].median()) > 0
                   for k in ("odd>even", "even>odd"))
        rows.append(dict(test=f"{tag}M4 (unit - MUA) at sink > elsewhere, half-max, {cond}", n=len(a), n_away=len(b),
                         median=float(a.median()), median_away=float(b.median()), p=p,
                         holds=bool(p < .05 and same and len(a) >= 10), inconclusive=bool(len(a) < 10)))
        for g in ("above", "below"):
            Dg = per_unit(M, cond, lambda m, g=g: m.grp == g)
            rows.append(one_sided(f"(reported) M2 unit - MUA > 0, half-max, {g} sink, {cond}", Dg, "unit_minus_mua_half"))
    return pd.DataFrame(rows)


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    pens = sorted(S[S.usable].penetration)
    with ProcessPoolExecutor(n_workers) as ex:
        M = pd.DataFrame([r for res in ex.map(work, pens) for r in res])
    M = differences(M)
    M.to_parquet(out_root() / "binocular" / "v4_M_units.parquet")
    T = tests(M)
    T.to_csv(REPO / "docs" / "qc" / "v4_M.csv", index=False)
    s = M[M.condition == "strong pool"].groupby(["penetration", "cluster_id", "grp", "strict"])[
        [c for c in M.columns if c.endswith("_half") and not c.startswith("unit_minus") or c.startswith("unit_minus")]
    ].mean().reset_index()
    s.to_csv(REPO / "docs" / "qc" / "v4_M_units.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    g = M[M.condition == "strong pool"].drop_duplicates(["penetration", "cluster_id", "direction"])
    print("responsive units per group and direction:", g.groupby(["direction", "grp"]).size().to_dict())
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
