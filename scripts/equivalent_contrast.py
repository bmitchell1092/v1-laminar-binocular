"""Figure 2 core: equivalent contrast k(t) and linearity index L(t) of binocular combination (docs/09).

Usage: .venv\\Scripts\\python scripts\\equivalent_contrast.py
Per driven in-cortex MUA channel (driven in a monocular condition, run_crf.py): condition-mean time courses
(completed presentations, no saccade >= 1 deg), Gaussian-smoothed (sigma 3 ms), blank (0/0) subtracted.
Dominant eye D = larger monocular response at .9 (40-250 ms); eye weight a = contrast scale mapping the
non-dominant CRF onto D's (v1lb.binocular.eye_weight). Onset = first time the D response at .9 reaches half its
peak. For dioptic contrast c, k(t) = equivalent dominant-eye contrast / c and L(t) = (k - 1) / a (1 = linear
summation, 0 = ocularity invariance), computed where the D response at .9 exceeds 20% of its peak.
Population: per penetration the median over channels, then mean over penetrations with a bootstrap 95% CI over
penetrations. Windows relative to each channel's onset: early 0-20 ms, late 60-200 ms.
Writes <out_root>/binocular/k_t.parquet, docs/qc/equivalent_contrast_summary.csv and
<out_root>/binocular/equivalent_contrast.png. Also prints onset latency vs contrast (S1a latency control).
"""

import glob
import time
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
DIOPTIC = [0.055, 0.11, 0.225, 0.45]
T_KEEP = (-20, 320)
EARLY, LATE = (0, 20), (60, 200)
WINDOWS_REL = [(0, 10), (10, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150), (150, 200)]
WROWS: list = []                     # windowed estimates, filled by penetration()


def onset_halfmax(x, t, lo=10, hi=250):
    w = (t >= lo) & (t < hi)
    pk = np.nanmax(x[w])
    if not np.isfinite(pk) or pk <= 0:
        return np.nan
    idx = np.flatnonzero(w & (x >= 0.5 * pk))
    return float(t[idx[0]]) if idx.size else np.nan


def penetration(pen: str, driven: set) -> pd.DataFrame:
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
    Pk = P[keep].reset_index(drop=True)
    t = ds.time.values
    tw = (t >= T_KEEP[0]) & (t < T_KEEP[1])
    chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
    if not chans:
        return pd.DataFrame()
    mua = ds.mua.values[keep][:, :, chans].astype(float)
    cc, ci = Pk.contrast_contra.to_numpy(), Pk.contrast_ipsi.to_numpy()
    M = np.full((6, 6, t.size, len(chans)), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            sel = np.isclose(cc, a) & np.isclose(ci, b)
            if sel.sum() >= 3:
                M[i, j] = mua[sel].mean(axis=0)
    M = gaussian_filter1d(M - M[0, 0], 3, axis=2)
    w_full = (t >= 40) & (t < 250)
    rows = []
    for q, k in enumerate(chans):
        rc, ri = M[:, 0, :, q], M[0, :, :, q]                       # [contrast, time] per eye
        if np.isnan(rc).any() or np.isnan(ri).any():
            continue
        contra_dom = rc[-1, w_full].mean() >= ri[-1, w_full].mean()
        rd, rn = (rc, ri) if contra_dom else (ri, rc)
        a = eye_weight(rd[:, w_full].mean(axis=1), rn[:, w_full].mean(axis=1))
        odi = (rc[-1, w_full].mean() - ri[-1, w_full].mean()) / (abs(rc[-1, w_full].mean()) + abs(ri[-1, w_full].mean()))
        on = onset_halfmax(rd[-1], t)
        if not np.isfinite(on):
            continue
        present = rd[-1] >= 0.2 * np.nanmax(rd[-1, w_full])
        onsets_mono = {f"on_mono_{c:g}": onset_halfmax(rd[i], t) for i, c in enumerate(CONTRASTS) if c > 0}
        for c in DIOPTIC:
            i = list(CONTRASTS).index(c)
            rb = M[i, i, :, q]
            # window means relative to this channel's onset: the robust estimate (inverting per ms amplifies noise)
            for w0, w1 in WINDOWS_REL:
                sel = (t >= on + w0) & (t < on + w1)
                kw, cw = equivalent_contrast(rd[:, sel].mean(axis=1)[:, None], np.array([rb[sel].mean()]), c)
                kw = float(CONTRASTS[-1] / c) if cw[0] == 1 else float(kw[0])
                # eye weight refitted in this window: guards against a non-dominant eye whose monocular response
                # fades faster, which would lower L under purely linear summation (confound register #25)
                a_w = eye_weight(rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1))
                WROWS.append(dict(penetration=pen, site_id=str(ds.channel.values[k]), depth=ds.depth_below_toc.values[k],
                                  depth_l45=ds.depth_l45.values[k], a=a, a_w=a_w, odi=odi, onset=on, c=c, w0=w0, w1=w1,
                                  k=kw, cens=int(cw[0]), L=(kw - 1) / a if np.isfinite(kw) else np.nan,
                                  L_w=(kw - 1) / a_w if np.isfinite(kw) else np.nan))
            k_t, cens = equivalent_contrast(rd[:, tw], rb[tw], c)
            k_t = np.where(cens == 1, CONTRASTS[-1] / c, k_t)            # top-censored: lower bound
            ok = present[tw]
            rows.append(dict(penetration=pen, site_id=str(ds.channel.values[k]), depth=ds.depth_below_toc.values[k],
                             depth_l45=ds.depth_l45.values[k], a=a, odi=odi, onset=on, c=c,
                             on_bin=onset_halfmax(rb, t), **onsets_mono,
                             t=t[tw], k=np.where(ok, k_t, np.nan), cens=np.where(ok, cens, 0),
                             L=np.where(ok, (k_t - 1) / a, np.nan)))
    return pd.DataFrame(rows)


def pen_curve(df: pd.DataFrame, rel: bool, grid: np.ndarray) -> np.ndarray:
    """Per-penetration median L over channels on `grid` (stimulus time or time from each channel's onset)."""
    curves = {}
    for pen, g in df.groupby("penetration"):
        Ls = []
        for _, r in g.iterrows():
            x = r.t - (r.onset if rel else 0)
            Ls.append(np.interp(grid, x, r.L, left=np.nan, right=np.nan))
        curves[pen] = np.nanmedian(np.array(Ls), axis=0)
    return pd.DataFrame(curves, index=grid).T


def boot(m: np.ndarray, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    m = m[~np.all(np.isnan(m), axis=1)]
    b = np.array([np.nanmean(m[rng.integers(0, len(m), len(m))], axis=0) for _ in range(n)])
    return np.nanmean(m, axis=0), np.nanpercentile(b, 2.5, axis=0), np.nanpercentile(b, 97.5, axis=0)


def windowed_summary(W: pd.DataFrame, out: Path):
    """Primary statistic: per penetration, median L over channels in each onset-aligned window; across
    penetrations, median with a bootstrap CI. Plus the paired early (0-20) vs late (50-200) test."""
    rng = np.random.default_rng(0)
    rows = []
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3), facecolor="#fcfcfb")
    cols = {0.055: "#6250d6", 0.11: "#2a78d6", 0.225: "#eb6834", 0.45: "#1baf7a"}
    for c in DIOPTIC:
        g = W[W.c == c]
        per = g.groupby(["penetration", "w0"]).L.median().unstack()
        xs, ms, los, his = [], [], [], []
        for (w0, w1) in WINDOWS_REL:
            v = per[w0].dropna().to_numpy()
            b = [np.median(rng.choice(v, len(v))) for _ in range(2000)]
            xs.append((w0 + w1) / 2)
            ms.append(np.median(v))
            los.append(np.percentile(b, 2.5))
            his.append(np.percentile(b, 97.5))
            rows.append(dict(c=c, window=f"{w0}-{w1}", n_pen=len(v), L_median=np.median(v), lo=los[-1], hi=his[-1],
                             frac_top_censored=float((g[g.w0 == w0].cens == 1).mean())))
        axs[0].errorbar(xs, ms, [np.subtract(ms, los), np.subtract(his, ms)], color=cols[c], marker="o", ms=4,
                        lw=1.4, capsize=2, label=f"c = {c:g}")
        early = per[[0, 10]].mean(axis=1)
        late = per[[50, 75, 100, 150]].mean(axis=1)
        ok = early.notna() & late.notna()
        p = stats.wilcoxon(early[ok], late[ok]).pvalue
        rows.append(dict(c=c, window="early 0-20 vs late 50-200", n_pen=int(ok.sum()), L_median=np.median(early[ok]),
                         lo=np.median(late[ok]), hi=np.median(early[ok] - late[ok]), frac_top_censored=p))
    for ax in axs:
        ax.axhline(1, color="#0b0b0b", lw=0.8, ls="--")
        ax.axhline(2 ** (1 / 2.5) - 1, color="#52514e", lw=0.8, ls=":")
        ax.axhline(0, color="#52514e", lw=0.8)
    axs[0].set_xlabel("ms from each channel's response onset (window centre)", fontsize=8)
    axs[0].set_ylabel("linearity index L = (k - 1) / a  (median of penetrations)", fontsize=8)
    axs[0].set_title("binocular combination over time: 1 = linear sum, 0 = invariance, dotted = normalization n=2.5",
                     fontsize=8)
    axs[0].legend(fontsize=7, frameon=False)
    # depth, c = .225 and .11 pooled per penetration x depth bin
    bins = [(0, 7), (7, 14), (14, 21), (21, 32)]
    dcol = ["#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"]
    for (lo_d, hi_d), col in zip(bins, dcol):
        g = W[W.c.isin([0.11, 0.225]) & (W.depth >= lo_d) & (W.depth < hi_d)]
        per = g.groupby(["penetration", "w0"]).L.median().unstack()
        if len(per) < 5:
            continue
        ms = [np.nanmedian(per[w0]) for w0, _ in WINDOWS_REL]
        axs[1].plot([(a + b) / 2 for a, b in WINDOWS_REL], ms, color=col, marker="o", ms=4, lw=1.4,
                    label=f"{lo_d}-{hi_d} contacts below top (n pen {len(per)})")
        early = per[[0, 10]].mean(axis=1)
        late = per[[50, 75, 100, 150]].mean(axis=1)
        rows.append(dict(c=".11+.225", window=f"depth {lo_d}-{hi_d}: early vs late", n_pen=len(per),
                         L_median=np.nanmedian(early), lo=np.nanmedian(late), hi=np.nanmedian(early - late),
                         frac_top_censored=np.nan))
    axs[1].set_xlabel("ms from each channel's response onset", fontsize=8)
    axs[1].set_title("c = .11 and .225 by depth (median of penetrations)", fontsize=8)
    axs[1].legend(fontsize=6, frameon=False)
    for ax in axs:
        ax.set_ylim(-1, 2.5)
    fig.tight_layout()
    fig.savefig(out / "equivalent_contrast_windows.png", dpi=85, facecolor="#fcfcfb")
    S = pd.DataFrame(rows)
    S.to_csv(REPO / "docs" / "qc" / "equivalent_contrast_windows.csv", index=False)
    pd.set_option("display.width", 200)
    print("windowed (rows 'early vs late': L_median = early, lo = late, hi = early - late, last column = Wilcoxon p)")
    print(S.round(3).to_string(index=False))


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    parts = []
    for pen in pens:
        t0 = time.time()
        df = penetration(pen, driven)
        parts.append(df)
        print(f"{pen}: {df.site_id.nunique() if len(df) else 0} channels ({time.time() - t0:.0f}s)", flush=True)
    K = pd.concat(parts, ignore_index=True)
    out = out_root() / "binocular"
    out.mkdir(parents=True, exist_ok=True)
    K.to_parquet(out / "k_t.parquet")
    W = pd.DataFrame(WROWS)
    W.to_parquet(out / "k_windows.parquet")
    windowed_summary(W, out)

    # ---- latency control (S1a): onsets vs contrast, monocular and dioptic
    ch = K.drop_duplicates(["penetration", "site_id"])
    lat = {f"mono {c:g}": ch[f"on_mono_{c:g}"].median() for c in CONTRASTS[1:]}
    lat |= {f"dioptic {c:g}": K[K.c == c].on_bin.median() for c in DIOPTIC}
    print("median half-max onset (ms):", {k: round(v, 1) for k, v in lat.items()})
    print(f"eye weight a: median {ch.a.median():.2f} (IQR {ch.a.quantile(.25):.2f}-{ch.a.quantile(.75):.2f}); "
          f"channels {len(ch)} in {ch.penetration.nunique()} penetrations")

    # ---- population L(t), onset-aligned, and early vs late per penetration
    grid = np.arange(-10, 250)
    summary = []
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), facecolor="#fcfcfb")
    cols = {0.055: "#6250d6", 0.11: "#2a78d6", 0.225: "#eb6834", 0.45: "#1baf7a"}
    for c in DIOPTIC:
        g = K[K.c == c]
        C = pen_curve(g, True, grid)
        m, lo, hi = boot(C.values)
        axs[0].plot(grid, m, color=cols[c], lw=1.6, label=f"c = {c:g}")
        axs[0].fill_between(grid, lo, hi, color=cols[c], alpha=0.15, lw=0)
        e = C.loc[:, (grid >= EARLY[0]) & (grid < EARLY[1])].mean(axis=1)
        l = C.loc[:, (grid >= LATE[0]) & (grid < LATE[1])].mean(axis=1)
        ok = e.notna() & l.notna()
        p = stats.wilcoxon(e[ok], l[ok]).pvalue if ok.sum() >= 6 else np.nan
        # half-decay time: first time after the early window that the mean curve falls below (early + late) / 2
        mid = (np.nanmean(e) + np.nanmean(l)) / 2
        after = (grid >= EARLY[1]) & (m < mid)
        t_half = float(grid[np.argmax(after)]) if after.any() else np.nan
        frac_cens = float((g.explode("cens").cens.astype(float) == 1).mean())
        summary.append(dict(c=c, n_pen=int(ok.sum()), n_ch=g.site_id.groupby(g.penetration).nunique().sum(),
                            L_early=np.nanmean(e), L_early_lo=np.nanpercentile(e[ok], 25),
                            L_late=np.nanmean(l), L_late_lo=np.nanpercentile(l[ok], 25),
                            early_minus_late=np.nanmean(e[ok] - l[ok]), p_wilcoxon=p,
                            half_decay_ms_from_onset=t_half, frac_top_censored=frac_cens))
        if c in (0.11, 0.225):
            axs[1].scatter(np.full(ok.sum(), 0 if c == 0.11 else 2), e[ok], s=12, color=cols[c])
            axs[1].scatter(np.full(ok.sum(), 1 if c == 0.11 else 3), l[ok], s=12, color=cols[c])
            for x0, (ee, ll) in zip([0 if c == 0.11 else 2] * ok.sum(), zip(e[ok], l[ok])):
                axs[1].plot([x0, x0 + 1], [ee, ll], color=cols[c], lw=0.5, alpha=0.5)
    for ax in axs[:2]:
        ax.axhline(1, color="#0b0b0b", lw=0.8, ls="--")
        ax.axhline(0, color="#52514e", lw=0.8)
        ax.axhline(2 ** (1 / 2.5) - 1, color="#52514e", lw=0.8, ls=":")
    axs[0].set_xlabel("ms from each channel's response onset", fontsize=8)
    axs[0].set_ylabel("linearity index L = (k - 1) / a", fontsize=8)
    axs[0].set_title("binocular combination over time (1 = linear sum, 0 = invariance)", fontsize=9)
    axs[0].set_ylim(-0.5, 2)
    axs[0].legend(fontsize=7, frameon=False)
    axs[1].set_xticks([0, 1, 2, 3])
    axs[1].set_xticklabels(["early\n.11", "late\n.11", "early\n.225", "late\n.225"], fontsize=8)
    axs[1].set_title("per penetration (median over channels)", fontsize=9)
    axs[1].set_ylim(-0.5, 2.5)
    # depth: c = .225, L(t) by depth bin
    bins = [(0, 7), (7, 14), (14, 21), (21, 32)]
    dcol = ["#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"]
    for (lo_d, hi_d), col in zip(bins, dcol):
        g = K[(K.c == 0.225) & (K.depth >= lo_d) & (K.depth < hi_d)]
        if g.penetration.nunique() < 5:
            continue
        C = pen_curve(g, True, grid)
        m, lo, hi = boot(C.values)
        axs[2].plot(grid, m, color=col, lw=1.4, label=f"{lo_d}-{hi_d} contacts below top (n pen {len(C)})")
        axs[2].fill_between(grid, lo, hi, color=col, alpha=0.12, lw=0)
        e = C.loc[:, (grid >= EARLY[0]) & (grid < EARLY[1])].mean(axis=1)
        l = C.loc[:, (grid >= LATE[0]) & (grid < LATE[1])].mean(axis=1)
        summary.append(dict(c=0.225, depth=f"{lo_d}-{hi_d}", n_pen=len(C), L_early=e.mean(), L_late=l.mean(),
                            early_minus_late=(e - l).mean()))
    axs[2].axhline(1, color="#0b0b0b", lw=0.8, ls="--")
    axs[2].axhline(0, color="#52514e", lw=0.8)
    axs[2].set_ylim(-0.5, 2)
    axs[2].set_title("c = .225 by depth", fontsize=9)
    axs[2].legend(fontsize=6, frameon=False)
    axs[2].set_xlabel("ms from each channel's response onset", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "equivalent_contrast.png", dpi=85, facecolor="#fcfcfb")
    S = pd.DataFrame(summary)
    pd.set_option("display.width", 220)
    print(S.round(3).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "equivalent_contrast_summary.csv", index=False)
    print(out / "equivalent_contrast.png")


if __name__ == "__main__":
    main()
