"""Figure 3: contrast sets how far binocular summation collapses (docs/09 claim 3; docs/10 candidate 1).

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_contrast.py
(a) raw binocular gain: both eyes at c minus the dominant eye alone at c, in units of the dominant-eye peak,
    aligned to each channel's onset; penetration means, 95% bootstrap CI (no model, no inversion);
(b) linearity index in the transition window (sliding 15 ms windows starting 20-35 ms after each channel's onset)
    per penetration (grey) and the median (normalization_timing.py);
(c) slope of the transition linearity index on log contrast per recording: all penetrations by monkey, and the
    independent cross-recordings of the dual-probe days (monkey J; one-sided sign tests). The internal record of
    the held-out test (failed as first run, 6/7, before MUA artifact rejection) is in docs/10; the paper reports
    the cleaned data, like every other analysis (Blake, 2026-10-30: registration is internal);
(d) amplitude-normalized half-transition time for penetrations where the index drops by >= 0.2 at all three
    contrasts (lines; bars = medians): the collapse is deeper at higher contrast but not reliably earlier.
Inputs: figure_data.py (population.npz), <out_root>/binocular/normalization_timing.parquet and
normalization_timing_cross.parquet.
Writes <out_root>/figures/fig3_contrast.pdf and .png; prints the statistics shown.
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
MIN_CROSS_CHANNELS = 5                       # cross-recordings need at least this many driven channels (docs/10)


def per_recording(R: pd.DataFrame) -> pd.DataFrame:
    """Channel-median L per recording, contrast and window start."""
    return R[R.split == "all"].groupby(["name", "c", "s"]).L.median()


def transition(per: pd.Series) -> pd.DataFrame:
    """Per recording: median over window starts 20-35 ms (windows covering 20-50 ms), by contrast."""
    p = per.unstack("s")
    return p.loc[:, (p.columns >= 20) & (p.columns <= 35)].median(axis=1).unstack("c").reindex(columns=DIOPTIC)


def slopes(tr: pd.DataFrame) -> pd.Series:
    ok = tr.dropna()
    return ok.apply(lambda r: np.polyfit(np.log(DIOPTIC), r.values, 1)[0], axis=1)


def half_time(Ls: pd.Series) -> float:
    """As normalization_timing.half_time: first window start >= 10 ms where L has fallen halfway from its early
    (starts 0-10) to its late (75-150) value; defined only when the drop is >= 0.2."""
    e = Ls[(Ls.index >= 0) & (Ls.index <= 10)].mean()
    late = Ls[(Ls.index >= 75) & (Ls.index <= 150)].mean()
    if not np.isfinite(e) or not np.isfinite(late) or e - late < 0.2:
        return np.nan
    after = Ls[(Ls.index >= 10) & (Ls < (e + late) / 2)]
    return float(after.index[0]) if len(after) else np.nan


def gain(ax, pop):
    rel, O = pop["rel"], pop["onset"]
    sel = (rel >= -10) & (rel <= 120)
    for c in (0.11, 0.225, 0.45, 0.9):
        i = list(CONTRASTS).index(c)
        fs.ci_line(ax, rel[sel], (O[:, i, i] - O[:, i, 0])[:, sel], fs.CONTRAST[c], alpha=0.1)
    fs.zero_line(ax, color=fs.INK2)
    fs.line_key(ax, [(fs.contrast_label(c), fs.CONTRAST[c], "-") for c in (0.11, 0.225, 0.45, 0.9)], x=0.62, y=0.95,
                dy=0.075, length=0.07)
    ax.set_xlim(-10, 120)
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("both eyes − dominant eye\n(dominant-eye peak = 1)")


def transition_panel(ax, tr: pd.DataFrame, n_neg: int):
    for _, r in tr.dropna().iterrows():
        ax.plot(DIOPTIC, r.values.clip(-2, 3), color=fs.RULE, lw=0.5, zorder=1)
    med = tr.median()
    ax.plot(DIOPTIC, med.values, color=fs.INK, lw=1.2, zorder=2)
    for c in DIOPTIC:
        ax.plot([c], [med[c]], marker="o", ms=3.6, color=fs.CONTRAST[c], zorder=3)
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    fs.zero_line(ax)
    fs.log_contrast_axis(ax, DIOPTIC)
    ax.set_xticklabels(["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_xlim(0.09, 0.55)
    ax.set_ylim(-2.1, 3.1)
    ax.set_ylabel("linearity index, 20–50 ms after onset")
    ax.text(0.97, 0.97, f"slope < 0 in\n{n_neg}/{len(tr.dropna())} penetrations", transform=ax.transAxes, ha="right",
            va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)


def slope_panel(ax, sl: pd.Series, sl_x: pd.Series):
    rng = np.random.default_rng(1)
    groups = [("B", sl[sl.index.str[7] == "B"], 0, "o"), ("J", sl[sl.index.str[7] == "J"], 1.3, "s"),
              ("held-out\n(J)", sl_x, 2.75, "s")]
    for lab, v, x, mk in groups:
        xs = x + rng.uniform(-0.18, 0.18, len(v))
        ax.scatter(xs, v.clip(-6, 3), s=7, marker=mk, facecolor="white" if x > 2 else fs.INK2,
                   edgecolor=fs.INK, lw=0.5, zorder=3)
        ax.plot([x - 0.28, x + 0.28], [v.median()] * 2, color=fs.INK, lw=1.4, zorder=4)
        n_neg = int((v < 0).sum())
        p = stats.binomtest(n_neg, len(v), 0.5, alternative="greater").pvalue
        ax.text(x, 3.3, f"{n_neg}/{len(v)}\n{fs.fmt_p(p)}", ha="center", va="bottom", fontsize=5.5, color=fs.INK,
                linespacing=1.0)
    ax.axvline(2.05, color=fs.RULE, lw=0.6)
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks([0, 1.3, 2.75], ["B", "J", "J"], fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, 3.45)
    ax.set_ylim(-6.3, 3.2)
    ax.set_ylabel("slope of linearity index\non log contrast")
    for x, lab in ((0.65, "all penetrations"), (2.75, "independent\nrecordings")):
        ax.text(x, -7.6, lab, ha="center", va="top", fontsize=5.5, color=fs.INK2, linespacing=1.0)


def half_time_panel(ax, per: pd.Series):
    H = pd.Series({k: half_time(v.droplevel([0, 1])) for k, v in per.groupby(level=["name", "c"])}).unstack()
    H = H.reindex(columns=DIOPTIC)
    ok = H.dropna()
    for _, r in ok.iterrows():
        ax.plot(range(3), r.values, color=fs.RULE, lw=0.5, zorder=1)
    for k, c in enumerate(DIOPTIC):
        ax.plot([k - 0.25, k + 0.25], [ok[c].median()] * 2, color=fs.CONTRAST[c], lw=1.8, zorder=3)
    p = stats.friedmanchisquare(*[ok[c] for c in DIOPTIC]).pvalue
    ax.set_xticks(range(3), ["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_ylabel("half-transition time (ms after onset)")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(0, ok.to_numpy().max() * 1.35)
    ax.text(0.97, 0.97, f"{len(ok)} penetrations with a drop\nat all three contrasts\nFriedman {fs.fmt_p(p)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0)
    return H, p


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    R = pd.read_parquet(out_root() / "binocular" / "normalization_timing.parquet")
    X = pd.read_parquet(out_root() / "binocular" / "normalization_timing_cross.parquet")
    per, per_x = per_recording(R), per_recording(X)
    nch = X[X.split == "all"].groupby("name").channel.nunique()
    per_x = per_x[per_x.index.get_level_values("name").isin(nch.index[nch >= MIN_CROSS_CHANNELS])]
    tr, tr_x = transition(per), transition(per_x)
    sl, sl_x = slopes(tr), slopes(tr_x)

    fig = fs.figure("double", 62)
    axs = fig.subplots(1, 4, gridspec_kw=dict(width_ratios=[1.25, 0.95, 1.2, 0.85]))
    gain(axs[0], pop)
    transition_panel(axs[1], tr, int((sl < 0).sum()))
    slope_panel(axs[2], sl, sl_x)
    H, p_f = half_time_panel(axs[3], per)
    for ax, letter in zip(axs, "abcd"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig3_contrast"))
    print(f"transition slope: negative in {(sl < 0).sum()}/{len(sl)} (B {(sl[sl.index.str[7] == 'B'] < 0).sum()}/"
          f"{(sl.index.str[7] == 'B').sum()}, J {(sl[sl.index.str[7] == 'J'] < 0).sum()}/{(sl.index.str[7] == 'J').sum()}), "
          f"median {sl.median():+.2f}")
    print(f"held-out: negative in {(sl_x < 0).sum()}/{len(sl_x)}, median {sl_x.median():+.2f}, one-sided sign test p "
          f"{stats.binomtest(int((sl_x < 0).sum()), len(sl_x), 0.5, alternative='greater').pvalue:.3f}")
    print("half-transition medians:", H.median().round(1).to_dict(), "n defined:", H.notna().sum().to_dict(),
          f"Friedman p {p_f:.3f}")


if __name__ == "__main__":
    main()
