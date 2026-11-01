"""Supplementary Figure S6: response dynamics of the fitted models (docs/09, figure set for V3; absorbs panels d-e of
the former fig5_models.py and the former figS_offset.py).

Usage: .venv\\Scripts\\python scripts\\figures\\figS6_dynamics.py
(a) dominant eye alone at 11% and 90%, channel-averaged and in units of each channel's 90% peak: data (mean over
    penetrations, 95% bootstrap CI) and the mean of each model's all-presentation fit (rise_noise_check.py,
    dominant_traces.npz);
(b) 10-90% rise time of the population-mean response vs contrast: data (95% penetration-bootstrap CI; upper ends
    beyond the axis are drawn at its edge, open triangles) and model fits (rise_noise_check.py). Per-penetration rise
    times are not shown: at <= 11% they are dominated by noise (docs/09 correction of 2026-10-31);
(c) the response after stimulus offset at 45%, dominant eye alone, aligned to the response offset (stimulus duration
    + onset latency) and scaled to the sustained level before it: data (mean, 95% CI) and the fixed-delay fit
    (median; out of sample, the model was fitted to 0-219 ms). Decay statistics over this period mostly measure the
    offset transient, so they are not used to compare mechanisms (docs/09 register #33).
Inputs: <out_root>/binocular/dominant_traces.npz, fig5_model.parquet; docs/qc/rise_noise.csv; figure_data.py
(population.npz, population_meta.parquet).
Writes <out_root>/figures/figS6_dynamics.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

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
REPO = Path(__file__).resolve().parents[2]
RISE_C = [0.055, 0.11, 0.225, 0.45, 0.9]
RISE_MAX = 45.0                        # ms axis limit in (b)
TM = np.arange(0, 600, 1.0)
MODEL_KEYS = ("g", "sigma", "n", "delta", "tau1", "tau2")
MODELS = {"full": ("fixed-delay", fs.MODEL["fixed-delay"]), "cond": ("conductance", fs.MODEL["conductance"]),
          "org": ("recurrent", fs.MODEL["recurrent"])}


def time_courses(axs, tr):
    A, src, t = tr["A"], list(tr["sources"]), tr["t"]
    for ax, c in zip(axs, (0.11, 0.9)):
        i = list(CONTRASTS).index(c)
        for m, (_, col) in MODELS.items():
            ax.plot(t, A[:, src.index(m), i].mean(axis=0), color=col, lw=1.0, zorder=3)
        fs.ci_line(ax, t, A[:, src.index("data"), i], fs.MODEL["data"], lw=1.2, alpha=0.15, zorder=4)
        fs.zero_line(ax)
        ax.set_xlim(0, 150)
        ax.text(0.97, 0.94, f"{fs.contrast_label(c)}, one eye", transform=ax.transAxes, ha="right", va="top",
                fontsize=5.5, color=fs.INK2)
    axs[0].set_ylim(-0.03, 0.42)
    axs[1].set_ylim(-0.08, 1.1)
    axs[0].tick_params(labelbottom=False)
    axs[1].set_xlabel("time from stimulus onset (ms)")
    axs[0].set_ylabel("response (dominant-eye 90% peak = 1)", y=0.0)
    return len(A)


def rise(ax, Q):
    q = Q[Q.measure == "population mean"].set_index(["source", "c"])
    x = np.arange(len(RISE_C))
    for k, (lab, col) in MODELS.items():
        ax.plot(x, q.loc[k].reindex(RISE_C).rise.to_numpy(), color=col, marker="o", ms=2.8, lw=1.0, zorder=3)
    d = q.loc["data"].reindex(RISE_C)
    m, lo, hi = d.rise.to_numpy(), d.lo.to_numpy(), np.minimum(d.hi.to_numpy(), RISE_MAX)
    ax.errorbar(x, m, [m - lo, hi - m], color=fs.MODEL["data"], marker="o", ms=3, lw=1.2, elinewidth=0.7, capsize=0,
                zorder=4)
    off = d.hi.to_numpy() > RISE_MAX
    ax.scatter(x[off], np.full(off.sum(), RISE_MAX), marker="^", s=10, facecolor="white", edgecolor=fs.INK, lw=0.5,
               zorder=5, clip_on=False)
    fs.line_key(ax, [("data", fs.MODEL["data"], "-", "o")] + [(lab, col, "-", "o") for lab, col in MODELS.values()],
                x=0.52, y=0.97, dy=0.075, length=0.1)
    ax.set_xticks(x, ["5.5", "11", "22.5", "45", "90"])
    ax.set_xlim(-0.4, len(RISE_C) - 0.6)
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel("rise time of the mean\nresponse, 10–90% (ms)")
    ax.set_ylim(0, RISE_MAX)
    return q[["rise", "lo", "hi"]].unstack("c"), int(off.sum())


def first_cross(t, x, level, after, below=False):
    idx = np.flatnonzero((t >= after) & ((x <= level) if below else (x >= level)))
    return float(t[idx[0]]) if idx.size else np.nan


def fixed_delay(p: dict, t, c, dur):
    """Fitted fixed-delay model (v1lb.dynorm 'full') for the dominant eye alone, with the stimulus switched off."""
    on = st.gamma.cdf(np.maximum(t - p["delta"], 0.0), dynorm.K_GAMMA, scale=p["tau1"])
    off = st.gamma.cdf(np.maximum(t - p["delta"] - dur, 0.0), dynorm.K_GAMMA, scale=p["tau1"])
    u = on - off
    k = np.exp(-np.arange(0, 6 * p["tau2"] + 1, 1.0) / p["tau2"])
    v = np.convolve(u, k / k.sum())[: len(t)]
    return p["g"] * (c * u) ** p["n"] / (p["sigma"] ** p["n"] + (c * v) ** p["n"])


def offset_level(t, x, dur):
    """Response offset (stimulus duration + half-max onset) and the sustained level over the 30 ms before it."""
    pk = np.nanmax(x[(t >= 0) & (t < 150)])
    t_off = dur + first_cross(t, x, 0.5 * pk, 0)
    S = x[(t >= t_off - 30) & (t < t_off)].mean() if np.isfinite(t_off) else np.nan
    return t_off, S, pk


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
        t_off, S, pk = offset_level(t, x, dur)
        if not np.isfinite(t_off) or t_off + 80 > t[-1] or S <= 0.1 * pk:
            continue
        Yd.append(np.interp(t_off + rel, t, x) / S)
        y = fixed_delay({kk: M.loc[pen, f"p_{kk}"] for kk in MODEL_KEYS}, TM, 0.45, dur)
        toff_m, Sm, _ = offset_level(TM, y, dur)
        Ym.append(np.interp(toff_m + rel, TM, y) / Sm)
    Yd = np.array(Yd)
    fs.ci_line(ax, rel, Yd, fs.MODEL["data"], alpha=0.15)
    ax.plot(rel, np.nanmedian(Ym, axis=0), color=MODELS["full"][1], lw=1.2)
    fs.zero_line(ax)
    ax.axhline(0.5, color=fs.RULE, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_xlim(-30, 80)
    ax.set_xlabel("time from response offset (ms)")
    ax.set_ylabel("response (sustained level = 1)")
    ax.text(0.97, 0.97, f"45%, one eye\n{len(Yd)} penetrations", transform=ax.transAxes, ha="right", va="top",
            fontsize=5.5, color=fs.INK2, linespacing=1.0)
    fs.line_key(ax, [("data", fs.MODEL["data"], "-"), ("fixed-delay fit", MODELS["full"][1], "-")], x=0.5,
                y=0.72, dy=0.075, length=0.08)
    mean = np.nanmean(Yd, axis=0)
    k = int(np.nanargmax(mean))
    ax.annotate("offset transient", xy=(rel[k], mean[k]), xytext=(12, 4), textcoords="offset points",
                fontsize=5.5, color=fs.INK, arrowprops=dict(arrowstyle="-", lw=0.5, color=fs.INK2))
    return dict(n=len(Yd), peak_after_offset=float(mean[k]), peak_time=int(rel[k]),
                model_peak_after_offset=float(np.nanmax(np.nanmedian(Ym, axis=0)[rel > 0])))


def main():
    fs.apply()
    b = out_root() / "binocular"
    d = out_root() / "figures" / "data"
    tr = dict(np.load(b / "dominant_traces.npz"))
    Q = pd.read_csv(REPO / "docs" / "qc" / "rise_noise.csv")
    pop = dict(np.load(d / "population.npz"))
    meta = pd.read_parquet(d / "population_meta.parquet").set_index("penetration")
    M = pd.read_parquet(b / "fig5_model.parquet").set_index("penetration")

    fig = fs.figure("double", 58)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1, 1])
    ax_a = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    ax_b, ax_c = fig.add_subplot(gs[:, 1]), fig.add_subplot(gs[:, 2])
    n_tc = time_courses(ax_a, tr)
    R, n_off = rise(ax_b, Q)
    O = offset_panel(ax_c, pop, meta, M)
    for ax, letter in zip((ax_a[0], ax_b, ax_c), "abc"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS6_dynamics"))
    print(f"(a) time courses: n = {n_tc} penetrations")
    print(f"(b) rise of the population-mean response (ms; CI ends above {RISE_MAX:g} ms drawn at the edge: {n_off}):")
    print(R.round(1).to_string())
    print(f"(c) offset: {O}")


if __name__ == "__main__":
    main()
