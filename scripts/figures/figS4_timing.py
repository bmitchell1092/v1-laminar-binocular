"""Supplementary Figure S4: timing measures and noise (docs/09, figure set for V3; the former Fig 3a, 3d, 4c and 4d,
with the split-half checks of V3 analysis D).

Usage: .venv\\Scripts\\python scripts\\figures\\figS4_timing.py
(a) raw binocular gain: both eyes at c minus the dominant eye alone at c, in units of the dominant-eye peak, aligned
    to each channel's onset; penetration means, 95% bootstrap CI (no model, no inversion; former Fig 3a). Because
    binocular responses start 1-2 ms earlier, part of the early difference is latency, which is why the paper's
    measure compares responses within windows aligned to each channel's onset;
(b) amplitude-normalized half-transition time (normalization_timing.half_time: first window start >= 10 ms where the
    channel-median linearity index has fallen halfway from its early to its late value; defined only for drops
    >= 0.2), for all presentations and for the odd and even presentations separately (v3_timing_noise.py): medians
    and interquartile ranges over the penetrations defined at all three contrasts, and Friedman tests over them;
(c) the laminar sequence (laminar_timing.py, v3_timing_noise.py) for channels above, at (within 2 contacts) and below
    the earliest sink: left, response onset and the time normalization is reached (centre of the first 15 ms window
    starting >= 20 ms after the group's onset with linearity index < 0.5), on the stimulus clock; right, the time
    normalization is reached after each group's own onset; dots = penetrations, markers = medians with interquartile
    ranges; Wilcoxon tests of each group against the channels at the sink;
(d) binocular gain (both eyes minus the dominant eye alone, at the same contrast) of the current at the earliest sink
    and of MUA on the same channels, on the sink-onset clock, each scaled by its dominant-eye peak at 90% in the first
    50 ms (sink_data.py); medians and 95% bootstrap CI over penetrations (medians: with this scaling a few
    penetrations with small or late-growing sinks dominate the mean; Fig 3e rescales and shows means).
Inputs: figure_data.py (population.npz), sink_data.py (sink.npz), <out_root>/binocular/v3_D_halftime.parquet and
v3_D_laminar.parquet (v3_timing_noise.py).
Writes <out_root>/figures/figS4_timing.pdf and .png; prints the statistics shown.
"""

import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
DIOPTIC = [0.11, 0.225, 0.45]
SPLITS = {"all": ("all presentations", fs.INK, "o", fs.INK), "odd": ("odd presentations", fs.INK2, "o", "white"),
          "even": ("even presentations", fs.INK2, "s", "white")}
GROUPS = ["above sink", "at sink", "below sink"]
GROUP_LABEL = {"above sink": "above", "at sink": "at sink", "below sink": "below"}


def gain(ax, pop):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= -10) & (rel <= 120)
    for c in (0.11, 0.225, 0.45, 0.9):
        i = list(CONTRASTS).index(c)
        fs.ci_line(ax, rel[sel], (O[:, i, i] - O[:, i, 0])[:, sel], fs.CONTRAST[c], alpha=0.1)
    fs.zero_line(ax, color=fs.INK2)
    fs.line_key(ax, [(fs.contrast_label(c), fs.CONTRAST[c], "-") for c in (0.11, 0.225, 0.45, 0.9)], x=0.7, y=0.95,
                dy=0.075, length=0.07)
    ax.set_xlim(-10, 120)
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("both eyes − dominant eye\n(dominant-eye peak = 1)")
    ax.text(0.7, 0.62, "same contrast\nin each eye", transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=fs.INK2, linespacing=1.0)


def half_times(ax, H: pd.DataFrame):
    rows = []
    for k, (split, (lab, col, mk, mfc)) in enumerate(SPLITS.items()):
        P = H[H.variant == split].pivot_table(index="level_0", columns="level_1", values="t_half")
        P = P.reindex(columns=DIOPTIC)
        ok = P.dropna()                                       # penetrations defined at all three contrasts
        x = np.arange(3) + (k - 1) * 0.22
        q1, md, q3 = (ok.quantile(q) for q in (0.25, 0.5, 0.75))
        ax.errorbar(x, md, [md - q1, q3 - md], color=col, lw=0.9, elinewidth=0.7, marker=mk, ms=3.4, mfc=mfc,
                    mec=col, capsize=0)
        p = stats.friedmanchisquare(*[ok[c] for c in DIOPTIC]).pvalue
        rows.append(dict(split=split, n_defined={c: int(P[c].notna().sum()) for c in DIOPTIC}, n_friedman=len(ok),
                         medians_complete=md.round(1).to_dict(), friedman_p=p,
                         medians_all_defined=P.median().round(1).to_dict()))
        ax.plot([0.04], [0.97 - k * 0.075], marker=mk, ms=3.2, mfc=mfc, mec=col, color=col, transform=ax.transAxes,
                clip_on=False)
        ax.text(0.08, 0.97 - k * 0.075, f"{lab}: Friedman {fs.fmt_p(p)} (n = {len(ok)})", transform=ax.transAxes,
                va="center", fontsize=5.5, color=fs.INK)
    ax.set_xticks(range(3), ["11", "22.5", "45"])
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(0, 80)
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_ylabel("half-transition time (ms after onset)")
    return pd.DataFrame(rows)


def laminar_clock(ax, T):
    rng = np.random.default_rng(3)
    for row, g in enumerate(GROUPS):
        G = T[T.grp == g]
        med = {}
        for var, mfc in (("onset", "white"), ("t_abs", fs.INK)):
            v = G[var].dropna()
            ax.scatter(v.clip(15, 130), row + rng.uniform(-0.18, 0.18, len(v)), s=4, color=fs.RULE, zorder=1, lw=0)
            q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
            ax.plot([q1, q3], [row, row], color=fs.INK, lw=0.7, zorder=2)
            ax.plot([md], [row], marker="o", ms=4.5, mfc=mfc, mec=fs.INK, mew=0.8, zorder=4)
            med[var] = md
        ax.plot([med["onset"], med["t_abs"]], [row, row], color=fs.INK, lw=1.6, zorder=3)
    ax.set_yticks(range(3), [GROUP_LABEL[g] for g in GROUPS])
    ax.set_ylim(2.6, -0.9)
    ax.set_xlim(15, 130)
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.tick_params(axis="y", length=0)
    ax.set_ylabel("channels relative to earliest sink")
    for k, (mfc, lab) in enumerate((("white", "response onset"), (fs.INK, "normalized (index < 0.5)"))):
        ax.plot([0.42], [0.97 - k * 0.07], marker="o", ms=4, mfc=mfc, mec=fs.INK, mew=0.8, transform=ax.transAxes,
                clip_on=False)
        ax.text(0.45, 0.97 - k * 0.07, lab, transform=ax.transAxes, va="center", fontsize=5.5, color=fs.INK)


def laminar_own(ax, T):
    rng = np.random.default_rng(4)
    for row, g in enumerate(GROUPS):
        v = T[T.grp == g].t_rel.dropna()
        ax.scatter(v.clip(0, 110), row + rng.uniform(-0.18, 0.18, len(v)), s=4, color=fs.RULE, zorder=1, lw=0)
        q1, md, q3 = v.quantile([0.25, 0.5, 0.75])
        ax.plot([q1, q3], [row, row], color=fs.INK, lw=0.7, zorder=2)
        ax.plot([md], [row], marker="o", ms=4.5, mfc=fs.INK, mec=fs.INK, mew=0.8, zorder=4)
    ax.set_ylim(2.6, -0.9)
    ax.set_xlim(0, 110)
    ax.set_xlabel("time after the group's own onset (ms)")
    ax.tick_params(axis="y", length=0, labelleft=False)
    ax.text(0.97, 0.97, "normalized (index < 0.5);\nmedian difference from\nthe sink channels: 0 ms",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)


def laminar_tests(T) -> pd.DataFrame:
    rows = []
    for var in ("onset", "t_abs", "t_rel"):
        piv = T.pivot_table(index="penetration", columns="grp", values=var)
        for g in ("above sink", "below sink"):
            ok = piv[["at sink", g]].dropna()
            d = ok[g] - ok["at sink"]
            rows.append(dict(variable=var, comparison=f"{g} - at sink", n=len(ok), median_at=ok["at sink"].median(),
                             median_group=ok[g].median(), median_diff=d.median(), n_later=int((d > 0).sum()),
                             n_tied=int((d == 0).sum()), p=stats.wilcoxon(ok[g], ok["at sink"]).pvalue))
    return pd.DataFrame(rows)


def sink_gain(ax, Z, c):
    rel = Z["rel"]
    i = list(CONTRASTS).index(c)
    sel = (rel >= -10) & (rel <= 150)
    out = {}
    for sig, ls in (("csd", (0, (3, 1.5))), ("mua", "-")):
        X = Z[sig]
        out[sig] = fs.ci_line(ax, rel[sel], (X[:, i, i] - X[:, i, 0])[:, sel], fs.INK, ls=ls, alpha=0.12,
                              stat=np.nanmedian)
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xlim(-10, 150)
    ax.set_xlabel("time from sink onset (ms)")
    ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
    return rel[sel], out


def main():
    fs.apply()
    pd.set_option("display.width", 200)
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    H = pd.read_parquet(out_root() / "binocular" / "v3_D_halftime.parquet")
    L = pd.read_parquet(out_root() / "binocular" / "v3_D_laminar.parquet")

    fig = fs.figure("double", 112)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1], hspace=0.05)
    ta = top.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.2, 1]))
    gain(ta[0], pop)
    HT = half_times(ta[1], H)
    ba = bot.subplots(1, 4, gridspec_kw=dict(width_ratios=[1.15, 0.85, 0.8, 0.8]))
    T_all = L[L.variant == "all"]
    laminar_clock(ba[0], T_all)
    laminar_own(ba[1], T_all)
    sg = {c: sink_gain(ax, Z, c) for ax, c in ((ba[2], 0.225), (ba[3], 0.45))}
    ba[3].sharey(ba[2])
    ba[3].tick_params(labelleft=False)
    ba[2].set_ylabel("both eyes − dominant eye\n(dominant-eye peak = 1)")
    fs.line_key(ba[2], [("current at the sink", fs.INK, (0, (3, 1.5))), ("MUA, same channels", fs.INK, "-")],
                x=0.25, y=0.95, dy=0.08, length=0.15)
    for ax, letter in zip((ta[0], ta[1], ba[0], ba[2]), "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "figS4_timing"))

    rel, O = pop["rel"], pop["onset"]
    print("(a) raw binocular gain, penetration means, peak over 0-40 ms and mean over 50-200 ms after onset:")
    for c in (0.11, 0.225, 0.45, 0.9):
        i = list(CONTRASTS).index(c)
        g = np.nanmean(O[:, i, i] - O[:, i, 0], axis=0)
        e, late = (rel >= 0) & (rel <= 40), (rel >= 50) & (rel < 200)
        print(f"   {fs.contrast_label(c)}: peak {np.nanmax(g[e]):.3f} at {rel[e][np.nanargmax(g[e])]} ms; "
              f"late mean {np.nanmean(g[late]):+.3f}")
    print("(b) half-transition times:")
    print(HT.to_string(index=False))
    print("(c) laminar timing (all presentations), Wilcoxon against the channels at the sink:")
    print(laminar_tests(T_all).round(4).to_string(index=False))
    for split in ("odd", "even"):
        print(f"    {split} presentations:")
        print(laminar_tests(L[L.variant == split]).round(4).to_string(index=False))
    print("(d) binocular gain at the sink, medians over penetrations:")
    for c, (t, out) in sg.items():
        for sig, m in out.items():
            e, late = (t >= 0) & (t <= 40), (t >= 50) & (t <= 150)
            print(f"   {fs.contrast_label(c)} {sig}: peak {np.nanmax(m[e]):.3f} at {t[e][np.nanargmax(m[e])]} ms; "
                  f"mean 50-150 ms {np.nanmean(m[late]):+.3f}")


if __name__ == "__main__":
    main()
