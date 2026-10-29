"""Candidate 1 (docs/10): does interocular normalization arrive *sooner* at higher contrast, or only end
*stronger*? Sliding-window linearity and amplitude-normalized half-transition times.

Usage: .venv\\Scripts\\python scripts\\normalization_timing.py [--data-dir datasets_cross] [--tag cross]
Per driven in-cortex MUA channel (one-way ANOVA across the dominant eye's monocular contrasts, 40-250 ms,
p < .01; same rule as the discovery), condition means (completed presentations, no saccade >= 1 deg), Gaussian
sigma 3 ms, blank-subtracted. Sliding 15 ms windows starting 0, 5, ..., 150 ms after the channel's onset: k and
L = (k - 1) / a with the eye weight refitted per window (v1lb.binocular), for dioptic .11, .225, .45.
Per penetration: the channel median L(s). Half-transition time per penetration x contrast: first s >= 10 ms at which
L(s) has fallen halfway from its early value (s in 0-10) to its late value (s in 75-150); defined only when the
early-late drop is >= 0.2. Odd/even presentation splits are computed for reliability.
Writes <out_root>/binocular/normalization_timing[_tag].parquet and .png; prints the tests.
"""

import sys
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import open_penetration  # noqa: E402
from v1lb.saccades import around_onsets, detect_continuous  # noqa: E402

warnings.filterwarnings("ignore")
DIOPTIC = [0.11, 0.225, 0.45]
STARTS = np.arange(0, 151, 5)
WIDTH = 15


def eye_table(d: Path, name: str, P: pd.DataFrame) -> pd.DataFrame:
    f = d / f"{name}_presentation_eye.parquet"
    if f.exists():
        return pd.read_parquet(f)
    import xarray as xr
    e = xr.open_dataset(d / f"{name}_eye.nc", engine="h5netcdf")
    s = detect_continuous(e.tick.values, e.eye_x.values.astype(float), e.eye_y.values.astype(float),
                          e.pupil.values.astype(float))
    return around_onsets(s, P.pd_on_tick.to_numpy())


def channel_curves(M, t, rd_idx_contra: bool):
    """Per-contrast sliding-window L for one channel; M [6, 6, time] condition means (blank-subtracted)."""
    rc, ri = M[:, 0], M[0, :]
    rd, rn = (rc, ri) if rd_idx_contra else (ri, rc)
    return rd, rn


def analyse(name: str, d: Path, split: str = "all"):
    ds = open_penetration(name, d)
    P = pd.read_parquet(d / f"{name}_presentations.parquet")
    E = eye_table(d, name, P)
    keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
    if split != "all":
        keep &= (np.arange(len(P)) % 2 == (1 if split == "odd" else 0))
    Pk = P[keep].reset_index(drop=True)
    t = ds.time.values
    cc, ci = Pk.contrast_contra.to_numpy(), Pk.contrast_ipsi.to_numpy()
    incx = np.flatnonzero(ds.in_cortex.values)
    mua = ds.mua.values[keep][:, :, incx].astype(float)
    w = (t >= 40) & (t < 250)
    win = mua[:, w].mean(axis=1)
    M = np.full((6, 6, t.size, len(incx)), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            sel = np.isclose(cc, a) & np.isclose(ci, b)
            if sel.sum() >= 3:
                M[i, j] = mua[sel].mean(axis=0)
    M = gaussian_filter1d(M - M[0, 0], 3, axis=2)
    rows, raw = [], []
    for q in range(len(incx)):
        rc, ri = M[:, 0, :, q], M[0, :, :, q]
        if np.isnan(rc).any() or np.isnan(ri).any():
            continue
        contra_dom = rc[-1, w].mean() >= ri[-1, w].mean()
        rd, rn = (rc, ri) if contra_dom else (ri, rc)
        # driven: ANOVA across the dominant eye's monocular contrasts (other eye at 0)
        groups = [win[(np.isclose(cc, c) & np.isclose(ci, 0)) if contra_dom else (np.isclose(ci, c) & np.isclose(cc, 0)), q]
                  for c in CONTRASTS]
        if min(len(g) for g in groups) < 3 or stats.f_oneway(*groups).pvalue >= 0.01:
            continue
        pk = np.nanmax(rd[-1][(t >= 10) & (t < 250)])
        on_idx = np.flatnonzero((t >= 10) & (rd[-1] >= 0.5 * pk))
        if not on_idx.size or pk <= 0:
            continue
        on = t[on_idx[0]]
        for c in DIOPTIC:
            i = list(CONTRASTS).index(c)
            rb = M[i, i, :, q]
            for s0 in STARTS:
                sel = (t >= on + s0) & (t < on + s0 + WIDTH)
                dd, nn, bb = rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1), rb[sel].mean()
                a = eye_weight(dd, nn)
                k, cens = equivalent_contrast(dd[:, None], np.array([bb]), c)
                k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
                rows.append(dict(name=name, channel=q, c=c, s=s0, L=(k - 1) / a))
            rel = np.arange(-10, 200)
            mono = np.interp(on + rel, t, rd[i])
            bino = np.interp(on + rel, t, rb)
            raw.append(dict(name=name, channel=q, c=c, rel=rel, gain=(bino - mono) / pk))
    return pd.DataFrame(rows), raw


def half_time(Ls: pd.Series) -> float:
    e = Ls[(Ls.index >= 0) & (Ls.index <= 10)].mean()
    late = Ls[(Ls.index >= 75) & (Ls.index <= 150)].mean()
    if not np.isfinite(e) or not np.isfinite(late) or e - late < 0.2:
        return np.nan
    mid = (e + late) / 2
    after = Ls[(Ls.index >= 10) & (Ls < mid)]
    return float(after.index[0]) if len(after) else np.nan


def main():
    args = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    d = out_root() / args.get("data-dir", "datasets")
    tag = args.get("tag", "")
    names = sorted(p.name.removesuffix("_presentations.parquet") for p in d.glob("*_presentations.parquet"))
    allrows, raws = [], []
    for split in (("all", "odd", "even") if not tag else ("all",)):
        for name in names:
            r, raw = analyse(name, d, split)
            if len(r):
                allrows.append(r.assign(split=split))
            if split == "all":
                raws += raw
        print(f"{split}: done", flush=True)
    R = pd.concat(allrows, ignore_index=True)
    out = out_root() / "binocular"
    R.to_parquet(out / f"normalization_timing{'_' + tag if tag else ''}.parquet")
    per = R.groupby(["split", "name", "c", "s"]).L.median()
    nch = R[R.split == "all"].groupby("name").channel.nunique()
    print("driven channels per recording:", nch.to_dict())

    # (1) slope of L on log contrast in the 20-50 ms window (the discovery / pre-registered statistic)
    def slope_tab(split):
        rows = []
        for name in per.xs(split, level="split").index.get_level_values("name").unique():
            p = per.xs((split, name), level=("split", "name")).unstack("c")
            v = p.loc[(p.index >= 20) & (p.index <= 35)].median()          # windows starting 20-35 cover 20-50 ms
            if v.reindex(DIOPTIC).notna().all():
                rows.append(dict(name=name, slope=np.polyfit(np.log(DIOPTIC), v.reindex(DIOPTIC).values, 1)[0]))
        return pd.DataFrame(rows)
    for split in sorted(R.split.unique()):
        st = slope_tab(split)
        if len(st):
            n_neg = int((st.slope < 0).sum())
            print(f"[{split}] 20-50 ms slope of L on log contrast: negative in {n_neg}/{len(st)} "
                  f"({n_neg / len(st):.0%}); one-sided sign test p {stats.binomtest(n_neg, len(st), .5, alternative='greater').pvalue:.4f}; "
                  f"median slope {st.slope.median():+.2f}")
    # (2) sooner vs stronger, from the 'all' split
    A = per.xs("all", level="split")
    ht, dr = {}, {}
    for (name, c), Ls in A.groupby(level=["name", "c"]):
        Ls = Ls.droplevel([0, 1])
        ht[(name, c)] = half_time(Ls)
        dr[(name, c)] = Ls[(Ls.index <= 10)].mean() - Ls[(Ls.index >= 75)].mean()
    H = pd.Series(ht).unstack()
    D = pd.Series(dr).unstack()
    print("\nhalf-transition time (ms after onset), median [n defined]:",
          {c: (round(H[c].median(), 1), int(H[c].notna().sum())) for c in DIOPTIC})
    print("early-late drop (strength), median:", {c: round(D[c].median(), 2) for c in DIOPTIC})
    ok = H[[0.11, 0.45]].dropna()
    if len(ok) >= 6:
        print(f"t_half(.45) < t_half(.11) in {(ok[0.45] < ok[0.11]).mean():.0%}, equal {(ok[0.45] == ok[0.11]).mean():.0%} "
              f"of {len(ok)}; Wilcoxon p {stats.wilcoxon(ok[0.45], ok[0.11]).pvalue:.4f}")
    ok3 = H.dropna()
    if len(ok3) >= 6:
        print(f"all three defined in {len(ok3)}; Friedman p {stats.friedmanchisquare(*[ok3[c] for c in DIOPTIC]).pvalue:.4f}")
    # odd/even reliability of the half-time
    if "odd" in R.split.unique():
        hs = {}
        for split in ("odd", "even"):
            Ps = per.xs(split, level="split")
            hs[split] = pd.Series({k: half_time(v.droplevel([0, 1])) for k, v in Ps.groupby(level=["name", "c"])})
        both = pd.concat(hs, axis=1).dropna()
        print(f"odd/even half-time reliability: Spearman r {stats.spearmanr(both.odd, both.even).correlation:.2f} (n {len(both)})")

    # raw-data figure: binocular gain (bino - mono at the same contrast) / monocular peak, onset-aligned
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), facecolor="#fcfcfb")
    cols = {0.11: "#2a78d6", 0.225: "#eb6834", 0.45: "#1baf7a"}
    rel = np.arange(-10, 200)
    rng = np.random.default_rng(0)
    for c in DIOPTIC:
        G = pd.DataFrame([r for r in raws if r["c"] == c])
        pen_mean = np.stack([np.mean(np.stack(g.gain.values), axis=0) for _, g in G.groupby("name")])
        m = pen_mean.mean(axis=0)
        B = np.stack([pen_mean[rng.integers(0, len(pen_mean), len(pen_mean))].mean(axis=0) for _ in range(1000)])
        lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
        axs[0].plot(rel, m, color=cols[c], lw=1.5, label=f"c = {c:g} (n {len(pen_mean)})")
        axs[0].fill_between(rel, lo, hi, color=cols[c], alpha=0.15, lw=0)
        mm = A.xs(c, level="c").groupby(level="s").median()
        axs[1].plot(mm.index + WIDTH / 2, mm.values, color=cols[c], marker="o", ms=3, lw=1.4, label=f"c = {c:g}")
    axs[0].axhline(0, color="#52514e", lw=0.8)
    axs[0].set_xlabel("ms from each channel's response onset", fontsize=8)
    axs[0].set_ylabel("(binocular - monocular) / monocular peak", fontsize=8)
    axs[0].set_title("raw binocular gain: dominant eye at c alone vs both eyes at c", fontsize=9)
    axs[0].legend(fontsize=7, frameon=False)
    axs[1].axhline(1, color="#0b0b0b", lw=0.8, ls="--")
    axs[1].axhline(0, color="#52514e", lw=0.8)
    axs[1].set_ylim(-1, 2)
    axs[1].set_xlabel("ms from onset (sliding 15 ms window centre)", fontsize=8)
    axs[1].set_ylabel("linearity index L (median of recordings)", fontsize=8)
    axs[1].set_title("linearity over time by contrast", fontsize=9)
    axs[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fn = out / f"normalization_timing{'_' + tag if tag else ''}.png"
    fig.savefig(fn, dpi=85, facecolor="#fcfcfb")
    print(fn)


if __name__ == "__main__":
    main()
