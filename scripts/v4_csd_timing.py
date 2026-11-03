"""V4 check J: does the current at the earliest sink start and peak before the MUA there? (docs/09, registered
2026-11-03)

Usage: .venv\\Scripts\\python scripts\\v4_csd_timing.py
Signals: sink.npz (sink current = -CSD over the sink channel and its two neighbours, and MUA on the same channels,
sink-onset clock, Gaussian sigma 2 ms, each scaled to its dominant-eye peak at 90%), dominant eye at 90% and both eyes
at 45%. Same estimators for both: half-max onset of the own peak (-5...150 ms, interpolated; v4_latency.onset) and peak
time (0-80 ms). Per penetration, MUA - current; one-sided Wilcoxon (> 0), raw and with the filter bias removed.
Filter bias: the two acquisition chains applied to the same test responses (rise 4, 8, 12 ms, 10-90%):
* current: causal 0.3 Hz first-order high-pass and 500 Hz fourth-order low-pass at 30 kHz (the amplifier's LFP
  stream, NSx header), sampled at 1 kHz, Gaussian sigma 2 ms;
* MUA: noise carrier modulated by the response at 30 kHz, zero-phase 500-5000 Hz band-pass, rectification, zero-phase
  200 Hz low-pass (signals.analog_mua_epochs), every 30th sample, Gaussian sigma 2 ms; mean over 200 noise draws.
Bias = (MUA chain - truth) - (current chain - truth); bias-corrected difference = raw - bias (median over rises).
Writes docs/qc/v4_J.csv (tests and bias) and docs/qc/v4_J_penetrations.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import signal, stats
from scipy.ndimage import gaussian_filter1d

from v4_latency import boot_ci, onset

from v1lb.binocular import CONTRASTS
from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]
FS = 30_000
CELLS = {"dominant eye 90%": (0.9, 0.0), "both eyes 45%": (0.45, 0.45)}


def peak_time(x, t, lo=0, hi=80):
    w = (t >= lo) & (t < hi)
    return float(t[w][np.nanargmax(x[w])])


def test_response(rise_ms: float, t_ms: np.ndarray, t0=30.0) -> np.ndarray:
    """sin^2 rise to 1 (10-90% in rise_ms), then decay to a sustained 0.4 (tau 30 ms)."""
    tr = rise_ms / 0.59
    s = np.clip((t_ms - t0) / tr, 0, 1)
    r = np.sin(np.pi * s / 2) ** 2
    after = t_ms >= t0 + tr
    r[after] = 0.4 + 0.6 * np.exp(-(t_ms[after] - t0 - tr) / 30.0)
    return r


def chains(rise_ms: float, n_draws: int = 200, seed: int = 0):
    t30 = np.arange(-100 * FS // 1000, 300 * FS // 1000) / FS * 1000            # ms at 30 kHz
    r = test_response(rise_ms, t30)
    t1 = t30[::30]
    truth = r[::30]
    hp = signal.butter(1, 0.3, "highpass", fs=FS, output="sos")
    lp = signal.butter(4, 500, "lowpass", fs=FS, output="sos")
    cur = gaussian_filter1d(signal.sosfilt(lp, signal.sosfilt(hp, r))[::30], 2)
    bp = signal.butter(4, (500, 5000), "bandpass", fs=FS, output="sos")
    lp2 = signal.butter(4, 200, "lowpass", fs=FS, output="sos")
    rng = np.random.default_rng(seed)
    env = np.zeros(t30.size)
    for _ in range(n_draws):
        x = rng.standard_normal(t30.size) * (0.2 + r)                            # background + driven activity
        env += signal.sosfiltfilt(lp2, np.abs(signal.sosfiltfilt(bp, x)))
    mua = env / n_draws
    mua = gaussian_filter1d(mua[::30] - mua[::30][t1 < 0].mean(), 2)
    truth_s = truth
    out = {}
    for name, x in (("current", cur), ("mua", mua)):
        out[name] = dict(onset=onset(x / np.nanmax(x[(t1 >= 10) & (t1 < 250)]), t1) - onset(truth_s, t1),
                         peak=peak_time(x, t1, 0, 120) - peak_time(truth_s, t1, 0, 120))
    return out


def main():
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    rel = Z["rel"].astype(float)
    rows = []
    for p, pen in enumerate(Z["pens"]):
        r = dict(penetration=str(pen))
        for lab, (cd, co) in CELLS.items():
            i, j = list(CONTRASTS).index(cd), list(CONTRASTS).index(co)
            cur, mua = Z["csd"][p, i, j], Z["mua"][p, i, j]
            for sig, x in (("current", cur), ("mua", mua)):
                pk = np.nanmax(x[(rel >= -5) & (rel < 150)])
                r[f"{sig}_onset_{lab}"] = onset(x / pk, rel, lo=-5, hi=150, floor=0) if pk > 0 else np.nan
                r[f"{sig}_peak_{lab}"] = peak_time(x, rel)
            r[f"d_onset_{lab}"] = r[f"mua_onset_{lab}"] - r[f"current_onset_{lab}"]
            r[f"d_peak_{lab}"] = r[f"mua_peak_{lab}"] - r[f"current_peak_{lab}"]
        rows.append(r)
    P = pd.DataFrame(rows)
    P.to_csv(REPO / "docs" / "qc" / "v4_J_penetrations.csv", index=False)
    B = []
    for rise in (4, 8, 12):
        ch = chains(rise)
        B.append(dict(rise_ms=rise, current_onset=ch["current"]["onset"], mua_onset=ch["mua"]["onset"],
                      current_peak=ch["current"]["peak"], mua_peak=ch["mua"]["peak"],
                      bias_onset=ch["mua"]["onset"] - ch["current"]["onset"],
                      bias_peak=ch["mua"]["peak"] - ch["current"]["peak"]))
    B = pd.DataFrame(B)
    bias = {"onset": float(B.bias_onset.median()), "peak": float(B.bias_peak.median())}
    out = []
    for lab in CELLS:
        for kind in ("onset", "peak"):
            v = P[f"d_{kind}_{lab}"].dropna()
            for corr in ("raw", "bias-corrected"):
                x = v - (bias[kind] if corr == "bias-corrected" else 0)
                lo, hi = boot_ci(x)
                out.append(dict(test=f"{'J1' if kind == 'onset' else 'J2'} MUA - current, {kind}, {lab}, {corr}",
                                n=len(x), median=float(x.median()), ci_lo=lo, ci_hi=hi,
                                n_in_direction=int((x > 0).sum()), p=stats.wilcoxon(x, alternative="greater").pvalue))
    T = pd.DataFrame(out)
    for lab in CELLS:
        for kind in ("onset", "peak"):
            sel = T.test.str.contains(f"{kind}, {lab}")
            T.loc[sel, "holds"] = bool((T[sel].p < .05).all())
    T = pd.concat([T, B.assign(test="filter bias (MUA - current), simulated")], ignore_index=True)
    T.to_csv(REPO / "docs" / "qc" / "v4_J.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print("filter chains, estimate - truth (ms):")
    print(B.round(2).to_string(index=False))
    print(T[T.test.str.startswith("J")].round(4).to_string(index=False))
    for lab in CELLS:
        print(f"{lab}: current onset median {P[f'current_onset_{lab}'].median():.1f}, MUA onset "
              f"{P[f'mua_onset_{lab}'].median():.1f}; current peak {P[f'current_peak_{lab}'].median():.0f}, MUA peak "
              f"{P[f'mua_peak_{lab}'].median():.0f} ms after sink onset")


if __name__ == "__main__":
    main()
