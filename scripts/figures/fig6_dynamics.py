"""Figure 6: the fitted delayed-normalization model rises too fast at low contrast, and the responses carry an
offset transient it cannot produce (docs/09, dynamics diagnostics; out of sample: the model was fitted to 0-219 ms).

Usage: .venv\\Scripts\\python scripts\\figures\\fig6_dynamics.py
(a) dominant eye alone, population responses on the stimulus clock, each penetration scaled to its peak;
(b) 10-90% rise time vs contrast: data (penetration medians of channels, dynamics_diagnostics.py; grey lines =
    penetrations, black = median with 95% bootstrap CI) and the fitted fixed-delay model (fig5_model.py parameters,
    per penetration, median);
(c) the response after stimulus offset at 45% (aligned to the response offset = stimulus duration + onset latency,
    scaled to the sustained level before it): data (mean, 95% CI) and the fitted model (median);
The half-decay statistic (dynamics_diagnostics.py) is not shown: (c) shows that it mostly measures the offset
transient, not the decay of the stimulus response (docs/09 register #33). No simulation panel: measured as the data
are, the toy recurrent models of the second literature report do not reproduce them either.
Inputs: figure_data.py (population.npz, population_meta.parquet), <out_root>/binocular/dynamics_channels.parquet,
fig5_model.parquet.
Writes <out_root>/figures/fig6_dynamics.pdf and .png; prints the numbers shown.
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
CS = [0.11, 0.225, 0.45, 0.9]
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


def rise_panel(ax, pop):
    t, F = pop["t"], pop["flash"]
    sel = (t >= 0) & (t <= 120)
    for c in CS:
        i = list(CONTRASTS).index(c)
        Y = F[:, i, 0, :]
        Y = Y / np.nanmax(Y[:, (t >= 0) & (t < 150)], axis=1, keepdims=True)
        fs.ci_line(ax, t[sel], Y[:, sel], fs.CONTRAST[c], alpha=0.12)
    fs.zero_line(ax)
    ax.set_xlim(0, 120)
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("MUA (own peak = 1)")
    fs.line_key(ax, [(fs.contrast_label(c), fs.CONTRAST[c], "-") for c in CS], x=0.66, y=0.42, dy=0.075, length=0.08)
    ax.text(0.03, 0.97, "dominant eye\nalone", transform=ax.transAxes, va="top", fontsize=5.5, color=fs.INK2,
            linespacing=1.0)


def vs_contrast(ax, per, mod, var, ylabel):
    piv = per[var].unstack("c").reindex(columns=CS)
    x = np.arange(len(CS))
    for _, r in piv.iterrows():
        ax.plot(x, r.values, color=fs.RULE, lw=0.4, zorder=1)
    m, lo, hi = fs.boot_ci(piv.to_numpy(), stat=np.nanmedian)
    ax.errorbar(x, m, [m - lo, hi - m], color=fs.MODEL["data"], marker="o", ms=3, lw=1.2, elinewidth=0.7, capsize=0,
                zorder=3)
    mm = mod[var].unstack("c").reindex(columns=CS).median()
    ax.plot(x, mm.values, color=fs.MODEL["fixed-delay"], marker="o", ms=3, lw=1.2, zorder=3)
    ax.set_xticks(x, ["11", "22.5", "45", "90"])
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 75 if var == "rise" else None)
    return m, mm


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
    D = pd.read_parquet(out_root() / "binocular" / "dynamics_channels.parquet")
    M = pd.read_parquet(out_root() / "binocular" / "fig5_model.parquet").set_index("penetration")

    per = D.groupby(["penetration", "c"])[["rise", "half_decay"]].median()
    mod = []
    for pen, r in M.iterrows():
        p = {k: r[f"p_{k}"] for k in MODEL_KEYS}
        dur = meta.loc[pen, "dur"] if pen in meta.index else 317.0
        for c in CS:
            rise, hd, _, _ = rise_decay(TM, fixed_delay(p, TM, c, dur), dur)
            mod.append(dict(penetration=pen, c=c, rise=rise, half_decay=hd))
    mod = pd.DataFrame(mod).groupby(["penetration", "c"]).median()

    fig = fs.figure("onehalf", 58)
    ax_a, ax_b, ax_c = fig.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.1, 0.9, 1.1]))
    rise_panel(ax_a, pop)
    rise_d, rise_m = vs_contrast(ax_b, per, mod, "rise", "rise time, 10–90% (ms)")
    offset_panel(ax_c, pop, meta, M)
    hd_d = per["half_decay"].unstack("c").reindex(columns=CS).median().to_numpy()
    hd_m = mod["half_decay"].unstack("c").reindex(columns=CS).median()
    fs.line_key(ax_b, [("data", fs.MODEL["data"], "-", "o"), ("fixed-delay fit", fs.MODEL["fixed-delay"], "-", "o")],
                x=0.4, y=0.97, dy=0.075, length=0.1)
    for ax, letter in zip((ax_a, ax_b, ax_c), "abc"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig6_dynamics"))
    print("rise (ms): data", dict(zip(CS, np.round(rise_d, 1))), "fixed-delay", rise_m.round(1).to_dict())
    print("half-decay statistic (not shown; offset-transient dominated): data", dict(zip(CS, np.round(hd_d, 1))),
          "fixed-delay", hd_m.round(1).to_dict())


if __name__ == "__main__":
    main()
