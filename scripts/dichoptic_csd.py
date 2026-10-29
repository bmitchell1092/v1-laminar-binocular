"""Figure 4: binocular combination in the earliest current sink vs in spiking (MUA) at the same depth (docs/09).

Usage: .venv\\Scripts\\python scripts\\dichoptic_csd.py
Per penetration (completed presentations, no saccade >= 1 deg; CSD from the dataset: bad contacts interpolated,
Hamming-smoothed, 05 C12):
* Earliest grating sink: mean CSD over presentations with either eye >= .45, baseline (-50...0) removed; among
  in-cortex channels whose 20-70 ms sink peak is >= 40% of the strongest, the one whose sink reaches 25% of its own
  peak first (its onset). Sink channels = that channel and its two neighbours. Usable if peak / baseline SD >= SINK_Z.
* Responses per contrast cell: sink amplitude = -CSD averaged over the sink channels (blank-subtracted, Gaussian
  sigma 2 ms); MUA averaged over the same channels.
* Windows relative to the sink onset (CSD) and to the MUA half-max onset at the same channels (MUA).
* Equivalent contrast k and linearity index L = (k - 1) / a per window (v1lb.binocular), eye weights refitted per
  window, for dioptic .11, .225, .45.
Writes <out_root>/binocular/csd_windows.parquet and docs/qc/dichoptic_csd_summary.csv; prints the comparison.
"""

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v1lb.binocular import CONTRASTS, equivalent_contrast, eye_weight
from v1lb.config import out_root
from v1lb.dataset import open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
SINK_Z = 5.0
DIOPTIC = [0.11, 0.225, 0.45]
WINDOWS = [(0, 10), (10, 20), (20, 35), (35, 50), (50, 75), (75, 100), (100, 150)]


def cell_means(x, cc, ci):
    """[pres, time] -> [6, 6, time] condition means (NaN where < 3 presentations)."""
    M = np.full((6, 6, x.shape[1]), np.nan)
    for i, a in enumerate(CONTRASTS):
        for j, b in enumerate(CONTRASTS):
            sel = np.isclose(cc, a) & np.isclose(ci, b)
            if sel.sum() >= 3:
                M[i, j] = x[sel].mean(axis=0)
    return M


def windowed_L(M, t, onset, signal, pen, extra):
    rows = []
    w_full = (t >= onset) & (t < onset + 150)
    rc, ri = M[:, 0], M[0, :]
    contra_dom = np.nanmean(rc[-1, w_full]) >= np.nanmean(ri[-1, w_full])
    rd, rn = (rc, ri) if contra_dom else (ri, rc)
    for c in DIOPTIC:
        i = list(CONTRASTS).index(c)
        for w0, w1 in WINDOWS:
            sel = (t >= onset + w0) & (t < onset + w1)
            d, n, b = rd[:, sel].mean(axis=1), rn[:, sel].mean(axis=1), M[i, i, sel].mean()
            if not (np.isfinite(d).all() and np.isfinite(n).all() and np.isfinite(b)):
                continue
            a = eye_weight(d, n)
            k, cens = equivalent_contrast(d[:, None], np.array([b]), c)
            k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
            rows.append(dict(penetration=pen, signal=signal, c=c, w0=w0, a=a, k=k, cens=int(cens[0]),
                             L=(k - 1) / a if np.isfinite(k) else np.nan, **extra))
    return rows


GRAND: dict = {}          # signal -> penetration -> [6 (dominant eye), 6, time from sink onset], normalised
GA_WINDOWS = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 55), (55, 75), (75, 100), (100, 150)]


def ga_L(stack: np.ndarray, rel: np.ndarray, c: float):
    """L per window from a grand-average 6x6 response (dominant eye on axis 0)."""
    G = np.nanmean(stack, axis=0)
    out = []
    i = list(CONTRASTS).index(c)
    for w0, w1 in GA_WINDOWS:
        sel = (rel >= w0) & (rel < w1)
        d, n, b = G[:, 0, sel].mean(axis=1), G[0, :, sel].mean(axis=1), G[i, i, sel].mean()
        a = eye_weight(d, n)
        k, cens = equivalent_contrast(d[:, None], np.array([b]), c)
        k = float(CONTRASTS[-1] / c) if cens[0] == 1 else float(k[0])
        out.append((k - 1) / a)
    return np.array(out)


def grand_average(out: Path, n_boot: int = 1000):
    """Population-level L from grand-average responses (both signals timed from the sink onset), with a
    bootstrap over penetrations: far less noisy than inverting each penetration's CSD contrast response."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rel = np.arange(-10, 200)
    rng = np.random.default_rng(0)
    pens = sorted(set(GRAND.get("csd_sink", {})) & set(GRAND.get("mua_same_depth", {})))
    rows = []
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), facecolor="#fcfcfb")
    xs = [(a + b) / 2 for a, b in GA_WINDOWS]
    for ax, c in zip(axs, DIOPTIC):
        diffs = []
        for sig, col in (("csd_sink", "#eb6834"), ("mua_same_depth", "#2a78d6")):
            stack = np.stack([GRAND[sig][p] for p in pens])
            L0 = ga_L(stack, rel, c)
            B = np.array([ga_L(stack[rng.integers(0, len(pens), len(pens))], rel, c) for _ in range(n_boot)])
            lo, hi = np.nanpercentile(B, [2.5, 97.5], axis=0)
            ax.errorbar(xs, L0, [L0 - lo, hi - L0], color=col, marker="o", ms=4, lw=1.4, capsize=2,
                        label="earliest current sink" if sig == "csd_sink" else "MUA, same channels")
            diffs.append(B)
            for (w0, w1), v, l, h in zip(GA_WINDOWS, L0, lo, hi):
                rows.append(dict(c=c, signal=sig, window=f"{w0}-{w1}", L=v, lo=l, hi=h, n_pen=len(pens)))
        D = diffs[0] - diffs[1]                                            # sink - MUA, paired by resample
        for (w0, w1), dl, dh in zip(GA_WINDOWS, *np.nanpercentile(D, [2.5, 97.5], axis=0)):
            rows.append(dict(c=c, signal="sink - MUA", window=f"{w0}-{w1}", L=np.nanmedian(D, axis=0)[GA_WINDOWS.index((w0, w1))],
                             lo=dl, hi=dh, n_pen=len(pens)))
        ax.axhline(1, color="#0b0b0b", lw=0.8, ls="--")
        ax.axhline(0, color="#52514e", lw=0.8)
        ax.set_ylim(-1.5, 2.5)
        ax.set_title(f"dioptic {c:g}: grand average, time from sink onset", fontsize=9)
        ax.set_xlabel("ms from earliest-sink onset", fontsize=8)
        ax.legend(fontsize=7, frameon=False)
    axs[0].set_ylabel("linearity index L (1 = linear sum, 0 = invariance)", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "dichoptic_csd_grand_average.png", dpi=85, facecolor="#fcfcfb")
    T = pd.DataFrame(rows)
    T.to_csv(REPO / "docs" / "qc" / "dichoptic_csd_grand_average.csv", index=False)
    pd.set_option("display.width", 200)
    print(T.round(3).to_string(index=False))


def main():
    al = pd.read_csv(REPO / "metadata" / "laminar_alignment_v2.csv").set_index("penetration")
    d = out_root() / "datasets"
    pens = sorted(p.name.removesuffix("_presentation_eye.parquet") for p in d.glob("*_presentation_eye.parquet"))
    rows, sinks = [], []
    for pen in pens:
        t0 = time.time()
        ds = open_penetration(pen)
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
        Pk = P[keep].reset_index(drop=True)
        t = ds.time.values
        cc, ci = Pk.contrast_contra.to_numpy(), Pk.contrast_ipsi.to_numpy()
        csd = ds.csd.values[keep].astype(float)
        csd -= np.nanmean(csd[:, (t >= -50) & (t < 0)], axis=1, keepdims=True)
        strong = np.maximum(cc, ci) >= 0.45
        m = np.nanmean(csd[strong], axis=0)                                  # [time, channel]
        incx = ds.in_cortex.values & np.isfinite(m).all(axis=0)
        w = (t >= 20) & (t < 70)
        mm = np.where(incx[None, :], m, np.nan)[w]
        # earliest input sink: among channels whose sink peak is >= 40% of the strongest sink, the one whose own
        # sink reaches 25% of its peak first (the strongest sink is often a later, intracortical one)
        peaks = np.nanmin(mm, axis=0)                                     # most negative per channel
        strong_ch = np.flatnonzero(incx & (peaks <= 0.4 * np.nanmin(peaks)))
        onsets = {}
        for k in strong_ch:
            idx = np.flatnonzero((t >= 15) & (t < 70) & (-m[:, k] >= 0.25 * abs(peaks[k])))
            if idx.size:
                onsets[k] = float(t[idx[0]])
        ch = min(onsets, key=lambda k: (onsets[k], peaks[k])) if onsets else int(np.nanargmin(peaks))
        peak = peaks[ch]
        ti = int(np.nanargmin(mm[:, ch]))
        z = abs(peak) / np.nanstd(m[(t >= -50) & (t < 0), ch])
        sink_on = onsets.get(ch, np.nan)
        chans = [k for k in (ch - 1, ch, ch + 1) if 0 <= k < m.shape[1] and incx[k]]
        info = dict(penetration=pen, sink_channel=str(ds.channel.values[ch]), sink_contact=int(ds.contact.values[ch]),
                    sink_depth_below_toc=float(ds.depth_below_toc.values[ch]),
                    sink_depth_l45=float(ds.depth_l45.values[ch]), sink_peak_ms=float(t[w][ti]), sink_onset_ms=sink_on,
                    sink_z=round(float(z), 1), usable=bool(z >= SINK_Z and np.isfinite(sink_on)),
                    l45_basis=al.loc[pen, "basis"] if pen in al.index else "")
        sinks.append(info)
        if not info["usable"]:
            print(f"{pen}: sink z {z:.1f} -> skipped", flush=True)
            continue
        sink_amp = gaussian_filter1d(-np.nanmean(csd[:, :, chans], axis=2), 2, axis=1)
        M = cell_means(sink_amp, cc, ci)
        M = M - M[0, 0]
        rows += windowed_L(M, t, sink_on, "csd_sink", pen, dict(depth=info["sink_depth_below_toc"]))
        mua = gaussian_filter1d(ds.mua.values[keep][:, :, chans].astype(float).mean(axis=2), 2, axis=1)
        Mm = cell_means(mua, cc, ci)
        Mm = Mm - Mm[0, 0]
        mx = np.nanmax(Mm[-1, 0][(t >= 10) & (t < 250)])
        on_m = t[np.flatnonzero((t >= 10) & (Mm[-1, 0] >= 0.5 * mx))[0]] if np.isfinite(mx) and mx > 0 else np.nan
        if np.isfinite(on_m):
            rows += windowed_L(Mm, t, float(on_m), "mua_same_depth", pen, dict(depth=info["sink_depth_below_toc"]))
        # same clock as the sink (time from the sink onset), for a direct sink-vs-spiking timing comparison
        rows += windowed_L(Mm, t, sink_on, "mua_sink_clock", pen, dict(depth=info["sink_depth_below_toc"],
                                                                      mua_onset_minus_sink=float(on_m - sink_on)))
        # grand-average inputs: both signals on one clock (time from the sink onset), dominant eye first,
        # each normalised by its dominant-eye response at .9 in the first 50 ms after the sink onset
        rel = np.arange(-10, 200)
        for sig, MM in (("csd_sink", M), ("mua_same_depth", Mm)):
            Mi = np.stack([[np.interp(sink_on + rel, t, MM[i, j]) for j in range(6)] for i in range(6)])
            w50 = (rel >= 0) & (rel < 50)
            dom_contra = np.nanmean(Mi[-1, 0, w50]) >= np.nanmean(Mi[0, -1, w50])
            Mi = Mi if dom_contra else Mi.transpose(1, 0, 2)                  # axis 0 = dominant eye
            scale = np.nanmax(Mi[-1, 0, w50])
            if np.isfinite(scale) and scale > 0:
                GRAND.setdefault(sig, {})[pen] = Mi / scale
        print(f"{pen}: sink {info['sink_channel']} ({info['sink_depth_below_toc']:.0f} below top, onset "
              f"{sink_on:.0f} ms, z {z:.0f}); MUA onset there {on_m} ms ({time.time() - t0:.0f}s)", flush=True)
    S = pd.DataFrame(sinks)
    R = pd.DataFrame(rows)
    out = out_root() / "binocular"
    out.mkdir(parents=True, exist_ok=True)
    grand_average(out)
    R.to_parquet(out / "csd_windows.parquet")
    S.to_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv", index=False)
    print(f"\nusable sinks: {S.usable.sum()}/{len(S)}; sink onset median {S[S.usable].sink_onset_ms.median():.0f} ms; "
          f"sink depth below top median {S[S.usable].sink_depth_below_toc.median():.0f} contacts")
    summ = []
    for c in DIOPTIC:
        for sig in ("csd_sink", "mua_same_depth"):
            g = R[(R.c == c) & (R.signal == sig)]
            per = g.groupby(["penetration", "w0"]).L.median().unstack()
            for w0, w1 in WINDOWS:
                if w0 in per:
                    v = per[w0].dropna()
                    summ.append(dict(c=c, signal=sig, window=f"{w0}-{w1}", n_pen=len(v), L_median=v.median(),
                                     L_q25=v.quantile(.25), L_q75=v.quantile(.75)))
        # paired: early (0-20) vs late (50-150) within each signal; and sink vs MUA late
        pc = R[R.c == c].groupby(["penetration", "signal", "w0"]).L.median().unstack("w0")
        early, late = pc[[0, 10]].mean(axis=1).unstack("signal"), pc[[50, 75, 100]].mean(axis=1).unstack("signal")
        for sig in ("csd_sink", "mua_same_depth"):
            ok = early[sig].notna() & late[sig].notna()
            summ.append(dict(c=c, signal=sig, window="early-late drop", n_pen=int(ok.sum()),
                             L_median=float(np.median(early[sig][ok] - late[sig][ok])),
                             L_q25=stats.wilcoxon(early[sig][ok], late[sig][ok]).pvalue if ok.sum() >= 6 else np.nan,
                             L_q75=np.nan))
        ok = late.notna().all(axis=1)
        if ok.sum() >= 6:
            summ.append(dict(c=c, signal="sink - MUA, late", window="50-150", n_pen=int(ok.sum()),
                             L_median=float(np.median(late.csd_sink[ok] - late.mua_same_depth[ok])),
                             L_q25=stats.wilcoxon(late.csd_sink[ok], late.mua_same_depth[ok]).pvalue, L_q75=np.nan))
    T = pd.DataFrame(summ)
    pd.set_option("display.width", 200)
    print("(rows 'early-late drop' and 'sink - MUA': L_median = median difference, L_q25 = Wilcoxon p)")
    print(T.round(3).to_string(index=False))
    T.to_csv(REPO / "docs" / "qc" / "dichoptic_csd_summary.csv", index=False)


if __name__ == "__main__":
    main()
