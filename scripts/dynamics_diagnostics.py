"""Recurrent vs fixed-delay normalization diagnostics (docs/09; predictions stated before computing).

Usage: .venv\\Scripts\\python scripts\\dynamics_diagnostics.py
Per driven in-cortex MUA channel, dominant eye alone at .11, .225, .45, .9 (completed presentations, no saccade >= 1
deg; Gaussian sigma 2 ms; blank-subtracted):
* rise time: 10% -> 90% of the response peak (peak within 150 ms of the flash);
* offset decay: response offset reference t_off = stimulus duration (photodiode, per penetration) + the channel's
  half-max onset latency at that contrast; sustained level S = mean over t_off - 30...t_off; half-decay time = first
  time after t_off at which the response falls to S / 2; undershoot = minimum over t_off...t_off + 100 relative to S
  (negative = below baseline). Decay is defined only where S > 10% of the peak.
Recurrent/shunting normalization predicts contrast-independent, slow decay and rise times that shorten with
contrast; a fixed-delay feedforward pool predicts fast decay (possibly undershooting) and contrast-independent rise.
Writes docs/qc/dynamics_diagnostics.csv.
"""

import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d

from v1lb.config import out_root
from v1lb.dataset import open_penetration

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
CS = [0.11, 0.225, 0.45, 0.9]


def first_cross(t, x, level, after, below=False):
    idx = np.flatnonzero((t >= after) & ((x <= level) if below else (x >= level)))
    return float(t[idx[0]]) if idx.size else np.nan


def main():
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    rows = []
    for pen in sorted({p for p, _ in driven}):
        d = out_root() / "datasets"
        ds = open_penetration(pen)
        P = pd.read_parquet(d / f"{pen}_presentations.parquet")
        E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
        keep = P.completed.to_numpy() & (E.n_large.to_numpy() == 0)
        dur = float(np.median((P.pd_off_tick - P.pd_on_tick)[keep & (P.pd_off_tick > 0)]) / 30.0)
        t = ds.time.values
        chans = [k for k, ch in enumerate(ds.channel.values) if ds.in_cortex.values[k] and (pen, str(ch)) in driven]
        if not chans:
            continue
        mua = ds.mua.values[keep][:, :, chans].astype(float)
        cc, ci = P.contrast_contra.to_numpy()[keep], P.contrast_ipsi.to_numpy()[keep]
        mean = lambda a, b: mua[np.isclose(cc, a) & np.isclose(ci, b)].mean(axis=0)  # noqa: E731
        blank = mean(0, 0)
        w = (t >= 40) & (t < 250)
        for q in range(len(chans)):
            rc9, ri9 = (mean(.9, 0) - blank)[:, q], (mean(0, .9) - blank)[:, q]
            dom_c = rc9[w].mean() >= ri9[w].mean()
            for c in CS:
                x = gaussian_filter1d(((mean(c, 0) if dom_c else mean(0, c)) - blank)[:, q], 2)
                pk_w = (t >= 0) & (t < 150)
                pk = np.nanmax(x[pk_w])
                if not np.isfinite(pk) or pk <= 0:
                    continue
                t10, t50, t90 = (first_cross(t, x, f * pk, 0) for f in (0.1, 0.5, 0.9))
                t_off = dur + (t50 if np.isfinite(t50) else np.nan)
                row = dict(penetration=pen, monkey=pen[7], channel=q, c=c, dur=dur, onset=t50, rise=t90 - t10)
                if np.isfinite(t_off) and t_off + 60 <= t[-1]:
                    S = x[(t >= t_off - 30) & (t < t_off)].mean()
                    if S > 0.1 * pk:
                        row["half_decay"] = first_cross(t, x, S / 2, t_off, below=True) - t_off
                        row["undershoot"] = x[(t >= t_off) & (t < t_off + 100)].min() / S
                rows.append(row)
        print(pen, flush=True)
    D = pd.DataFrame(rows)
    per = D.groupby(["penetration", "monkey", "c"])[["onset", "rise", "half_decay", "undershoot"]].median()
    out = []
    for var in ("onset", "rise", "half_decay", "undershoot"):
        piv = per[var].unstack("c").reindex(columns=CS)
        meds = piv.median()
        ok = piv.dropna()
        sl = ok.apply(lambda r: np.polyfit(np.log(CS), r.values, 1)[0], axis=1) if len(ok) else pd.Series(dtype=float)
        fr = stats.friedmanchisquare(*[ok[c] for c in CS]).pvalue if len(ok) >= 6 else np.nan
        mk = sl.groupby(level="monkey").median() if len(sl) else pd.Series(dtype=float)
        out.append(dict(variable=var, **{f"median_c{c:g}": meds[c] for c in CS}, n_pen=len(ok),
                        slope_median=sl.median() if len(sl) else np.nan, frac_slope_neg=(sl < 0).mean() if len(sl) else np.nan,
                        slope_B=mk.get("B", np.nan), slope_J=mk.get("J", np.nan), friedman_p=fr))
    S = pd.DataFrame(out)
    pd.set_option("display.width", 220)
    print("(slopes per log contrast; ms for onset/rise/half_decay; undershoot = min after offset / sustained level)")
    print(S.round(3).to_string(index=False))
    S.to_csv(REPO / "docs" / "qc" / "dynamics_diagnostics.csv", index=False)


if __name__ == "__main__":
    main()
