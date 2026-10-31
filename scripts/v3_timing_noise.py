"""V3 analysis D (docs/09, registered 2026-10-31): are the timing measures of Fig 3d and Fig 4d robust to noise?

Usage: .venv\\Scripts\\python scripts\\v3_timing_noise.py [n_workers] [--reuse]
Two threshold-crossing timing measures, each with the source script's own definitions:
* Fig 3d half-transition time (normalization_timing.half_time on the channel-median sliding-window L(s) per
  penetration and contrast; normalization_timing.parquet already holds the all / odd / even splits);
* Fig 4d laminar normalization time (laminar_timing.py logic: penetration-median L_w per group above / at / below the
  earliest sink, dioptic .225 and .45 pooled; t_rel = centre of the first window starting >= 20 ms with L < 0.5;
  t_abs = t_rel + the group's median onset). The odd and even halves are recomputed with
  equivalent_contrast.penetration itself, restricted to odd or even presentations (dominant eye, eye weight and onset
  re-estimated within each half), via a wrapper around its clean_presentations.
For each measure:
(i) split-half replication: the group differences computed in each half separately;
(ii) noise sensitivity: each measure on all presentations vs on each half (noise x sqrt 2), per group;
(iii) SNR-matched comparison: per group, signal = median over penetrations of the early - late drop of the L curve,
     noise = median over penetrations of the RMS over windows of the noise sample e = (L_odd - L_even) / 2 (same
     variance as the noise in the all-presentation curve); SNR = signal / noise. Every group with a higher SNR than
     the lowest gets lambda * e added (both signs), lambda = sqrt((SNR_g / SNR_min)^2 - 1), so all groups share the
     lowest SNR; the measure and the tests are then recomputed.
Writes <out_root>/binocular/v3_D_kwindows_halves.parquet, v3_D_laminar.parquet, v3_D_halftime.parquet and
docs/qc/v3_D.csv; prints the tables.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import equivalent_contrast as ec  # noqa: E402
from laminar_timing import CENTRES  # noqa: E402
from normalization_timing import half_time  # noqa: E402

from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
OUT = out_root() / "binocular"
LAM_GROUPS = ("above sink", "at sink", "below sink")
HT_C = (0.11, 0.225, 0.45)
_ORIG_CLEAN = ec.clean_presentations
_PARITY = {"odd": 1, "even": 0}


# ----------------------------------------------------------------------------------------------- odd / even windows
def _kw_half(task):
    pen, driven, split = task
    parity = _PARITY[split]
    ec.clean_presentations = lambda ds, P, E: _ORIG_CLEAN(ds, P, E) & (np.arange(len(P)) % 2 == parity)
    ec.WROWS.clear()
    ec.penetration(pen, driven)
    return [dict(r, split=split) for r in ec.WROWS]


def kwindows_halves(n_workers):
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    tasks = [(p, driven, s) for s in ("odd", "even") for p in pens]
    with ProcessPoolExecutor(n_workers) as ex:
        rows = [r for part in ex.map(_kw_half, tasks) for r in part]
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------- Fig 4d: laminar timing
def laminar_curves(W, S):
    """laminar_timing.main, up to the per-group curves: group onsets [pen, grp] and median L_w [pen, grp, w0]."""
    W = W[W.c.isin([0.225, 0.45]) & W.penetration.isin(S.index[S.usable])].copy()
    W["rel"] = W.depth - W.penetration.map(S.sink_depth_below_toc)
    W["grp"] = np.select([W.rel.abs() <= 2, W.rel < -2, W.rel > 2], ["at sink", "above sink", "below sink"], default="")
    on = W.drop_duplicates(["penetration", "site_id"]).groupby(["penetration", "grp"]).onset.median()
    per = W.groupby(["penetration", "grp", "w0"]).L_w.median().unstack("w0")
    return on, per


def t_rel_of(r: pd.Series) -> float:
    below = [w for w in r.index if w >= 20 and np.isfinite(r[w]) and r[w] < 0.5]
    return CENTRES[below[0]] if below else np.nan


def laminar_table(on, per):
    rows = [dict(penetration=pen, grp=grp, onset=on.get((pen, grp), np.nan), t_rel=t_rel_of(r))
            for (pen, grp), r in per.iterrows()]
    T = pd.DataFrame(rows)
    T["t_abs"] = T.t_rel + T.onset
    return T


def laminar_tests(T, variant):
    out = []
    for var in ("onset", "t_rel", "t_abs"):
        piv = T.pivot_table(index="penetration", columns="grp", values=var)
        for g in ("above sink", "below sink"):
            if g not in piv or "at sink" not in piv:
                continue
            ok = piv[["at sink", g]].dropna()
            d = ok[g] - ok["at sink"]
            p = stats.wilcoxon(ok[g], ok["at sink"]).pvalue if len(ok) > 5 and (d != 0).any() else np.nan
            out.append(dict(measure="Fig 4d laminar", variant=variant, variable=var, comparison=f"{g} - at sink",
                            n=len(ok), median_ref=ok["at sink"].median(), median_other=ok[g].median(),
                            median_diff=d.median(), frac_later=(d > 0).mean(), p=p))
    return out


def snr_table(curves_all, e, early, late, group_level):
    """Per group: signal = median over penetrations of early - late; noise = median RMS of e over windows."""
    rows = []
    for g, A in curves_all.groupby(level=group_level):
        E = e.loc[A.index]
        sig = (A[early].mean(axis=1) - A[late].mean(axis=1)).median()
        noise = np.sqrt((E ** 2).mean(axis=1)).median()
        rows.append(dict(group=g, signal=sig, noise=noise, snr=sig / noise, n=len(A)))
    S = pd.DataFrame(rows).set_index("group")
    S["lam"] = np.sqrt(np.clip((S.snr / S.snr.min()) ** 2 - 1, 0, None))
    return S


# ------------------------------------------------------------------------------------- Fig 3d: half-transition
def halftimes(per_split: pd.Series) -> pd.DataFrame:
    """per_split: channel-median L indexed [name, c, s] -> half-transition time [name x c]."""
    H = pd.Series({k: half_time(v.droplevel([0, 1])) for k, v in per_split.groupby(level=["name", "c"])}).unstack()
    return H.reindex(columns=list(HT_C))


def halftime_tests(H, variant):
    out = []
    ok3 = H.dropna()
    p_f = stats.friedmanchisquare(*[ok3[c] for c in HT_C]).pvalue if len(ok3) >= 6 else np.nan
    out.append(dict(measure="Fig 3d half-transition", variant=variant, variable="t_half", comparison="Friedman 11/22.5/45",
                    n=len(ok3), median_ref=ok3[0.11].median() if len(ok3) else np.nan,
                    median_other=ok3[0.45].median() if len(ok3) else np.nan, median_diff=np.nan, frac_later=np.nan,
                    p=p_f, note="medians " + "/".join(f"{ok3[c].median():.1f}" for c in HT_C) if len(ok3) else ""))
    for c in (0.225, 0.45):
        ok = H[[0.11, c]].dropna()
        d = ok[c] - ok[0.11]
        p = stats.wilcoxon(ok[c], ok[0.11]).pvalue if len(ok) > 5 and (d != 0).any() else np.nan
        out.append(dict(measure="Fig 3d half-transition", variant=variant, variable="t_half",
                        comparison=f"{c:g} - 0.11", n=len(ok), median_ref=ok[0.11].median(), median_other=ok[c].median(),
                        median_diff=d.median(), frac_later=(d > 0).mean(), p=p))
    return out


def shift_rows(measure, var, groups, all_tab, half_tabs):
    """(ii) noise sensitivity: per group, median over penetrations of (half - all)."""
    out = []
    for g in groups:
        for split, tab in half_tabs.items():
            both = pd.concat([all_tab[g], tab[g]], axis=1, keys=["all", "half"]).dropna()
            d = both.half - both["all"]
            out.append(dict(measure=measure, variant=f"{split} - all", variable=var, comparison=str(g), n=len(both),
                            median_ref=both["all"].median(), median_other=both.half.median(), median_diff=d.median(),
                            frac_later=(d > 0).mean(),
                            p=stats.wilcoxon(both.half, both["all"]).pvalue if len(both) > 5 and (d != 0).any() else np.nan,
                            note=f"defined: all {all_tab[g].notna().sum()}, {split} {tab[g].notna().sum()}"))
    return out


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 12
    rows = []

    # ================================================================ Fig 4d
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv").set_index("penetration")
    W_all = pd.read_parquet(OUT / "k_windows.parquet")
    f_half = OUT / "v3_D_kwindows_halves.parquet"
    if "--reuse" in sys.argv and f_half.exists():
        Wh = pd.read_parquet(f_half)
    else:
        Wh = kwindows_halves(n_workers)
        Wh.to_parquet(f_half)
    on_a, per_a = laminar_curves(W_all, S)
    curves = {"all": (on_a, per_a)}
    for split in ("odd", "even"):
        curves[split] = laminar_curves(Wh[Wh.split == split], S)
    tabs = {k: laminar_table(*v) for k, v in curves.items()}
    lam_rows = []
    for k, T in tabs.items():
        rows += laminar_tests(T, k)                               # all = Fig 4d as published; odd / even = (i)
        lam_rows.append(T.assign(variant=k))
    for var in ("t_rel", "t_abs", "onset"):                       # (ii)
        piv = {k: T.pivot_table(index="penetration", columns="grp", values=var) for k, T in tabs.items()}
        rows += shift_rows("Fig 4d laminar", var, LAM_GROUPS, piv["all"], {s: piv[s] for s in ("odd", "even")})
    # (iii) SNR matching on the group curves
    windows = list(per_a.columns)
    common = per_a.index.intersection(curves["odd"][1].index).intersection(curves["even"][1].index)
    A = per_a.loc[common, windows]
    e = ((curves["odd"][1].loc[common, windows] - curves["even"][1].loc[common, windows]) / 2)
    early = [w for w in windows if w < 20]
    late = [w for w in windows if w >= 50]
    snr_l = snr_table(A[A.index.get_level_values("grp") != ""], e, early, late, "grp")
    print("Fig 4d group SNR:\n", snr_l.round(3).to_string())
    for g, r in snr_l.iterrows():
        rows.append(dict(measure="Fig 4d laminar", variant="SNR", variable="L curve", comparison=g, n=int(r.n),
                         median_ref=r.signal, median_other=r.noise, median_diff=r.snr, frac_later=r.lam, p=np.nan,
                         note="median_ref = signal (early - late), median_other = noise RMS, median_diff = SNR, "
                              "frac_later = lambda added"))
    for sgn, lab in ((1, "SNR-matched +e"), (-1, "SNR-matched -e")):
        lam = A.index.get_level_values("grp").map(lambda g: snr_l.lam.get(g, 0.0)).to_numpy()[:, None]
        M = A + sgn * lam * e.fillna(0).to_numpy()
        T = laminar_table(on_a, M)
        rows += laminar_tests(T, lab)
        lam_rows.append(T.assign(variant=lab))
    pd.concat(lam_rows, ignore_index=True).to_parquet(OUT / "v3_D_laminar.parquet")

    # ================================================================ Fig 3d
    R = pd.read_parquet(OUT / "normalization_timing.parquet")
    per = R.groupby(["split", "name", "c", "s"]).L.median()
    Hs = {split: halftimes(per.xs(split, level="split")) for split in ("all", "odd", "even")}
    ht_rows = []
    for k, H in Hs.items():
        rows += halftime_tests(H, k)
        ht_rows.append(H.stack().rename("t_half").reset_index().assign(variant=k))
    rows += shift_rows("Fig 3d half-transition", "t_half", HT_C, Hs["all"], {s: Hs[s] for s in ("odd", "even")})
    # (iii): curves [name, c] x s
    Ca = per.xs("all", level="split").unstack("s")
    Co, Ce = per.xs("odd", level="split").unstack("s"), per.xs("even", level="split").unstack("s")
    common = Ca.index.intersection(Co.index).intersection(Ce.index)
    Ca, Co, Ce = Ca.loc[common], Co.loc[common], Ce.loc[common]
    e = (Co - Ce) / 2
    s_all = list(Ca.columns)
    early = [s for s in s_all if s <= 10]
    late = [s for s in s_all if 75 <= s <= 150]
    snr_h = snr_table(Ca, e, early, late, "c")
    print("Fig 3d contrast SNR:\n", snr_h.round(3).to_string())
    for c, r in snr_h.iterrows():
        rows.append(dict(measure="Fig 3d half-transition", variant="SNR", variable="L curve", comparison=f"{c:g}",
                         n=int(r.n), median_ref=r.signal, median_other=r.noise, median_diff=r.snr, frac_later=r.lam,
                         p=np.nan, note="median_ref = signal, median_other = noise RMS, median_diff = SNR, "
                                        "frac_later = lambda added"))
    for sgn, lab in ((1, "SNR-matched +e"), (-1, "SNR-matched -e")):
        lam = Ca.index.get_level_values("c").map(lambda c: snr_h.lam.get(c, 0.0)).to_numpy()[:, None]
        M = (Ca + sgn * lam * e.fillna(0).to_numpy()).stack().rename("L")
        H = halftimes(M)
        rows += halftime_tests(H, lab)
        ht_rows.append(H.stack().rename("t_half").reset_index().assign(variant=lab))
    pd.concat(ht_rows, ignore_index=True).to_parquet(OUT / "v3_D_halftime.parquet")

    Q = pd.DataFrame(rows)
    Q.to_csv(REPO / "docs" / "qc" / "v3_D.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 40)
    print(Q.drop(columns=["note"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
