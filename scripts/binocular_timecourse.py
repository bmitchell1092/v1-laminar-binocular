"""Population time course of interocular modulation by depth (first look for report Q1).

Usage: .venv\\Scripts\\python scripts\\binocular_timecourse.py
Per driven in-cortex MUA channel: dominant eye D (larger monocular response at .9, 40-250 ms); difference
R(D=.9, N) - R(D=.9, 0) normalised to the monocular peak, for N in {.225, .45, .9}; penetration means, then the
mean over penetrations with a bootstrap CI over penetrations; onsets = first 10 consecutive ms with the CI
excluding 0. Writes <out_root>/crf/interocular_timecourse.png and prints the onset table.
"""

import glob
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from v1lb.config import out_root
from v1lb.dataset import clean_presentations, open_penetration
warnings.filterwarnings("ignore")
R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
drv = R[(R.site == "mua") & R.condition.str.startswith("mono")].groupby(["penetration", "site_id"]).driven.any()
rows = []
for f in sorted(glob.glob(str(out_root() / "datasets" / "*_presentation_eye.parquet"))):
    pen = Path(f).name.replace("_presentation_eye.parquet", "")
    ds = open_penetration(pen)
    P = pd.read_parquet(out_root() / "datasets" / f"{pen}_presentations.parquet")
    E = pd.read_parquet(f)
    keep = clean_presentations(ds, P, E)
    Pk = P[keep].reset_index(drop=True)
    t = ds.time.values
    mua = ds.mua.values[keep].astype(float)
    mua -= mua[:, (t >= -50) & (t < 0)].mean(axis=1, keepdims=True)
    cc, ci = Pk.contrast_contra.to_numpy(), Pk.contrast_ipsi.to_numpy()
    def tc(a, b):
        sel = np.isclose(cc, a) & np.isclose(ci, b)
        return mua[sel].mean(0) if sel.sum() >= 3 else None
    for k, ch in enumerate(ds.channel.values):
        if not ds.in_cortex.values[k] or not drv.get((pen, str(ch)), False):
            continue
        a9, b9 = tc(0.9, 0), tc(0, 0.9)
        if a9 is None or b9 is None:
            continue
        rc, ri = a9[:, k], b9[:, k]
        w = (t >= 40) & (t < 250)
        dom = "c" if rc[w].mean() >= ri[w].mean() else "i"
        for D in (0.45, 0.9):
            base = tc(D, 0) if dom == "c" else tc(0, D)
            for N in (0.225, 0.45, 0.9):
                both = tc(D, N) if dom == "c" else tc(N, D)
                if base is None or both is None:
                    continue
                pk = base[w, k].max()
                rows.append(dict(penetration=pen, ch=str(ch), depth=ds.depth_below_toc.values[k], D=D, N=N,
                                 base=base[:, k] / pk, diff=(both[:, k] - base[:, k]) / pk))
    print(pen, len(rows), flush=True)
T = pd.DataFrame(rows)
t = np.arange(-200, 500)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), facecolor="#fcfcfb")
bins = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 32)]
cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
summary = []
for ax, N in zip(axs, (0.225, 0.45, 0.9)):
    for (lo, hi), col in zip(bins, cols):
        g = T[(T.D == 0.9) & (T.N == N) & (T.depth >= lo) & (T.depth < hi)]
        if len(g) < 10: continue
        # penetration-level means, then mean over penetrations; bootstrap over penetrations for CI
        per = g.groupby("penetration")["diff"].apply(lambda s: np.mean(np.stack(s.values), 0))
        M = np.stack(per.values)
        m = M.mean(0)
        rng = np.random.default_rng(0)
        B = np.stack([M[rng.integers(0, len(M), len(M))].mean(0) for _ in range(500)])
        lo_ci, hi_ci = np.percentile(B, [2.5, 97.5], axis=0)
        perb = g.groupby("penetration")["base"].apply(lambda s: np.mean(np.stack(s.values), 0))
        mb = np.stack(perb.values).mean(0)
        onset = t[np.argmax((t > 0) & (mb > 0.5 * mb[(t > 0) & (t < 250)].max()))]
        sig = (hi_ci < 0) & (t > 0)
        # first time of 10 consecutive significant ms
        run = np.convolve(sig, np.ones(10, int), "valid") == 10
        sup_on = t[np.argmax(run)] if run.any() else np.nan
        sigf = (lo_ci > 0) & (t > 0)
        runf = np.convolve(sigf, np.ones(10, int), "valid") == 10
        fac_on = t[np.argmax(runf)] if runf.any() else np.nan
        summary.append(dict(N=N, depth=f"{lo}-{hi}", n_pen=len(M), n_ch=len(g), half_max_onset=onset, facilitation_onset=fac_on, suppression_onset=sup_on,
                            mean_diff_40_100=m[(t >= 40) & (t < 100)].mean(), mean_diff_150_250=m[(t >= 150) & (t < 250)].mean()))
        ax.plot(t, m, color=col, lw=1.5, label=f"{lo}-{hi} below ToC (n pen {len(M)})")
        ax.fill_between(t, lo_ci, hi_ci, color=col, alpha=0.15, lw=0)
    ax.axhline(0, color="#52514e", lw=0.8); ax.set_xlim(-20, 350)
    ax.set_title(f"R(D=.9, N={N}) - R(D=.9, 0), normalised to monocular peak", fontsize=9)
    ax.set_xlabel("ms from photodiode flash", fontsize=8); ax.legend(fontsize=6, frameon=False)
fig.tight_layout(); fn = out_root() / "crf" / "interocular_timecourse.png"; fig.savefig(fn, dpi=80, facecolor="#fcfcfb")
pd.set_option("display.width", 200)
print(pd.DataFrame(summary).round(3).to_string(index=False)); print(fn)
