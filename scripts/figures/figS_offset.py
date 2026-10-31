"""Supplementary figure: responses carry an offset transient that no fitted model produces (docs/09 register #33;
formerly Fig 6c, moved when Figs 5 and 6 were merged into one model figure).

Usage: .venv\\Scripts\\python scripts\\figures\\figS_offset.py
The response after stimulus offset at 45%, dominant eye alone, aligned to the response offset (stimulus duration +
onset latency) and scaled to the sustained level before it: data (mean, 95% CI) and the fitted fixed-delay model
(median; out of sample, the model was fitted to 0-219 ms). Decay statistics over this period mostly measure the
transient, so they are not used to compare mechanisms.
Inputs: figure_data.py (population.npz, population_meta.parquet), <out_root>/binocular/fig5_model.parquet.
Writes <out_root>/figures/figS_offset.pdf and .png.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as st  # noqa: E402

from v1lb import dynorm  # noqa: E402
from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
TM = np.arange(0, 600, 1.0)
MODEL_KEYS = ("g", "sigma", "n", "delta", "tau1", "tau2")


def first_cross(t, x, level, after, below=False):
    idx = np.flatnonzero((t >= after) & ((x <= level) if below else (x >= level)))
    return float(t[idx[0]]) if idx.size else np.nan


def rise_decay(t, x, dur):
    """As dynamics_diagnostics.py: 10-90% rise of the peak within 150 ms; half-decay after t_off = dur + t50, from
    the sustained level S (mean over t_off - 30...t_off); returns (rise, half-decay, t_off, S)."""
    pk = np.nanmax(x[(t >= 0) & (t < 150)])
    t10, t50, t90 = (first_cross(t, x, f * pk, 0) for f in (0.1, 0.5, 0.9))
    t_off = dur + t50
    S = x[(t >= t_off - 30) & (t < t_off)].mean()
    hd = first_cross(t, x, S / 2, t_off, below=True) - t_off if S > 0.1 * pk else np.nan
    return t90 - t10, hd, t_off, S


def fixed_delay(p: dict, t, c, dur):
    """Fitted single-stage model (v1lb.dynorm 'full') for the dominant eye alone, with the stimulus switched off."""
    on = st.gamma.cdf(np.maximum(t - p["delta"], 0.0), dynorm.K_GAMMA, scale=p["tau1"])
    off = st.gamma.cdf(np.maximum(t - p["delta"] - dur, 0.0), dynorm.K_GAMMA, scale=p["tau1"])
    u = on - off
    k = np.exp(-np.arange(0, 6 * p["tau2"] + 1, 1.0) / p["tau2"])
    v = np.convolve(u, k / k.sum())[: len(t)]
    return p["g"] * (c * u) ** p["n"] / (p["sigma"] ** p["n"] + (c * v) ** p["n"])



def offset_panel(ax, pop, meta, M):
    t, F = pop["t"], pop["flash"]
    i45 = list(CONTRASTS).index(0.45)
    rel = np.arange(-30, 81)
    Yd, Ym = [], []
    for k, pen in enumerate(pop["pens"]):
        if pen not in M.index:
            continue
        dur = meta.loc[pen, "dur"]
        x = F[k, i45, 0]
        t_off = dur + first_cross(t, x, 0.5 * np.nanmax(x[(t >= 0) & (t < 150)]), 0)
        if not np.isfinite(t_off) or t_off + 80 > t[-1]:
            continue
        S = x[(t >= t_off - 30) & (t < t_off)].mean()
        if S <= 0.1 * np.nanmax(x[(t >= 0) & (t < 150)]):
            continue
        Yd.append(np.interp(t_off + rel, t, x) / S)
        y = fixed_delay({kk: M.loc[pen, f"p_{kk}"] for kk in MODEL_KEYS}, TM, 0.45, dur)
        _, _, toff_m, Sm = rise_decay(TM, y, dur)
        Ym.append(np.interp(toff_m + rel, TM, y) / Sm)
    fs.ci_line(ax, rel, np.array(Yd), fs.MODEL["data"], alpha=0.15)
    ax.plot(rel, np.nanmedian(Ym, axis=0), color=fs.MODEL["fixed-delay"], lw=1.2)
    fs.zero_line(ax)
    ax.axhline(0.5, color=fs.RULE, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_xlim(-30, 80)
    ax.set_xlabel("time from response offset (ms)")
    ax.set_ylabel("response (sustained level = 1)")
    ax.text(0.97, 0.97, f"45%, one eye\n{len(Yd)} penetrations", transform=ax.transAxes, ha="right", va="top",
            fontsize=5.5, color=fs.INK2, linespacing=1.0)
    fs.line_key(ax, [("data", fs.MODEL["data"], "-"), ("fixed-delay fit", fs.MODEL["fixed-delay"], "-")], x=0.5,
                y=0.72, dy=0.075, length=0.08)
    k = int(np.nanargmax(np.nanmean(Yd, axis=0)))
    ax.annotate("offset transient", xy=(rel[k], np.nanmean(Yd, axis=0)[k]), xytext=(12, 4), textcoords="offset points",
                fontsize=5.5, color=fs.INK, arrowprops=dict(arrowstyle="-", lw=0.5, color=fs.INK2))


def main():
    fs.apply()
    d = out_root() / "figures" / "data"
    pop = dict(np.load(d / "population.npz"))
    meta = pd.read_parquet(d / "population_meta.parquet").set_index("penetration")
    M = pd.read_parquet(out_root() / "binocular" / "fig5_model.parquet").set_index("penetration")
    fig = fs.figure("single", 58)
    ax = fig.subplots()
    offset_panel(ax, pop, meta, M)
    print(fs.save(fig, "figS_offset"))

if __name__ == "__main__":
    main()
