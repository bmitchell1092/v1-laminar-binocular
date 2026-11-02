"""Figure 3: contrast sets how far binocular summation collapses, and the sequence runs at every depth and in the
input current (docs/09, figure set for V3; combines the former Fig 3b-c and Fig 4a-b, plus a new sink-current panel).

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_contrast_depth.py
(a) linearity index in the transition window (sliding 15 ms windows starting 20-35 ms after each channel's onset,
    i.e. covering 20-50 ms) per penetration (grey) and the median (normalization_timing.py); values beyond the axis
    are drawn at its edge as open triangles;
(b) slope of the transition linearity index on log contrast per recording: all penetrations by monkey, and the
    held-out runs of the dual-probe days, in which the grating sat on the other probe's receptive field (monkey J,
    recordings with >= 5 driven channels; one-sided sign tests);
(c) linearity index by depth below the top of cortex (window-specific eye weights), early (0-20 ms after each
    channel's onset) and late (50-200 ms), for dioptic 11 and 22.5% pooled (as the former Fig 4b and fig3_depth.py)
    and for 45%; median of penetration medians per 0.4 mm bin, 95% bootstrap CI over penetrations. Per bin: early vs
    late (Wilcoxon over penetrations); late index and early-late drop on depth (mixed model, linear + quadratic,
    random intercept per penetration, as fig3_depth.py), overall and per monkey;
(d) population CSD aligned to each penetration's earliest current sink (presentations with either eye >= 45%; each
    penetration scaled by its sink's peak, then averaged; sink_data.py), with the median L4/5 boundary; drawn as a
    bilinearly interpolated image (Gouraud-shaded meshes leave saw-tooth seams between channel rows);
(e) current at the earliest sink (-CSD over the sink channel and its two neighbours) at 45% in the format of Fig 2a,
    on the sink-onset clock: dominant eye at 90%, dominant eye at 45%, other eye at 45%, both eyes at 45%. Each
    penetration is scaled by its mean absolute monocular sink current over 0-150 ms after sink onset (dominant 90%,
    dominant 45%, other 45%; a scale that does not involve the binocular condition), then all by one constant so the
    population mean of the dominant eye at 90% peaks at 1. Means with 95% penetration-bootstrap CI. (sink.npz scales
    each penetration by its early dominant-eye peak; that divides by small or unstable peaks where the sink reverses
    or grows late, so a few penetrations dominate the mean, which is why the former Fig 4c used medians. The mean is
    kept here because, unlike per-trace medians, it makes the traces' differences equal the mean per-penetration
    differences that the tests use.) Tests: per penetration, both eyes minus the dominant eye alone (and minus the
    dominant eye at 90%) in the windows of Fig 2 (0-20, 20-50, 50-200 ms after sink onset); Wilcoxon over
    penetrations.
Inputs: <out_root>/binocular/normalization_timing.parquet, normalization_timing_cross.parquet, k_windows.parquet;
<out_root>/figures/data/sink.npz (sink_data.py); docs/qc/dichoptic_csd_sinks.csv.
Writes <out_root>/figures/fig3_contrast_depth.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
DIOPTIC = [0.11, 0.225, 0.45]
MIN_CROSS_CHANNELS = 5                       # held-out runs need at least this many driven channels (docs/10)
EPOCHS = {"early": (0, 20), "late": (50, 200)}
EPOCH_COLOR = {"early": fs.MUTED, "late": fs.INK}
DEPTH_BINS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 32)]     # contacts below the top
DEPTH_SETS = {"11 and 22.5%": [0.11, 0.225], "45%": [0.45]}
SINK_WINDOWS = [(0, 20), (20, 50), (50, 200)]
TR_LIM = (-2.0, 3.0)                                                               # (a) axis limits                                      # ms after sink onset (Fig 2)


# ------------------------------------------------------------------------------------------------ (a, b) contrast
def per_recording(R: pd.DataFrame) -> pd.Series:
    """Channel-median L per recording, contrast and window start."""
    return R[R.split == "all"].groupby(["name", "c", "s"]).L.median()


def transition(per: pd.Series) -> pd.DataFrame:
    """Per recording: median over window starts 20-35 ms (windows covering 20-50 ms), by contrast."""
    p = per.unstack("s")
    return p.loc[:, (p.columns >= 20) & (p.columns <= 35)].median(axis=1).unstack("c").reindex(columns=DIOPTIC)


def slopes(tr: pd.DataFrame) -> pd.Series:
    ok = tr.dropna()
    return ok.apply(lambda r: np.polyfit(np.log(DIOPTIC), r.values, 1)[0], axis=1)


def transition_panel(ax, tr: pd.DataFrame):
    lo, hi = TR_LIM
    v = tr.dropna()
    for _, r in v.iterrows():
        ax.plot(DIOPTIC, r.values.clip(lo, hi), color=fs.RULE, lw=0.5, zorder=1)
    for out, edge, mk in ((v.to_numpy() > hi, hi, "^"), (v.to_numpy() < lo, lo, "v")):
        xs = np.array(DIOPTIC)[np.nonzero(out)[1]]
        ax.scatter(xs, np.full(len(xs), edge), marker=mk, s=8, facecolor="white", edgecolor=fs.INK2, lw=0.5,
                   zorder=2, clip_on=False)
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
    ax.set_ylim(TR_LIM[0] - 0.1, TR_LIM[1] + 0.1)
    ax.set_ylabel("linearity index")
    ax.set_title("20–50 ms after onset", fontsize=6.5)
    ax.text(0.97, 0.97, f"{len(v)} penetrations\ngrey: each; black: median", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.5, color=fs.INK, linespacing=1.0, zorder=4,
            bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    n_clip = int(((v < lo) | (v > hi)).to_numpy().sum())
    for y, lab in ((1.0, "linear sum"), (0.0, "dominant eye alone")):
        ax.text(0.093, y + 0.04, lab, ha="left", va="bottom", fontsize=5.5, color=fs.INK2, zorder=4,
                bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    return med, n_clip


def slope_panel(ax, sl: pd.Series, sl_x: pd.Series):
    rng = np.random.default_rng(1)
    groups = [("B", sl[sl.index.str[7] == "B"], 0, "o", fs.INK2), ("J", sl[sl.index.str[7] == "J"], 1.0, "s", fs.INK2),
              ("J", sl_x, 2.35, "s", "white")]
    out = {}
    for lab, v, x, mk, fc in groups:
        xs = x + rng.uniform(-0.16, 0.16, len(v))
        clipped = v < -6
        ax.scatter(xs[~clipped], v[~clipped], s=7, marker=mk, facecolor=fc, edgecolor=fs.INK, lw=0.5, zorder=3)
        ax.scatter(xs[clipped], np.full(clipped.sum(), -6.1), s=9, marker="v", facecolor="white", edgecolor=fs.INK,
                   lw=0.5, zorder=3, clip_on=False)
        ax.plot([x - 0.26, x + 0.26], [v.median()] * 2, color=fs.INK, lw=1.4, zorder=4)
        n_neg = int((v < 0).sum())
        p = stats.binomtest(n_neg, len(v), 0.5, alternative="greater").pvalue
        ax.text(x, 3.3, f"{n_neg}/{len(v)}\n{fs.fmt_p(p)}", ha="center", va="bottom", fontsize=5.5, color=fs.INK,
                linespacing=1.0)
        out[(lab, x)] = (n_neg, len(v), p, v.median(), int(clipped.sum()))
    ax.axvline(1.68, color=fs.RULE, lw=0.6)
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks([0, 1.0, 2.35], ["monkey B", "monkey J", "monkey J,\nheld-out\nruns"], fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 2.85)
    ax.set_ylim(-6.3, 3.2)
    ax.set_ylabel("slope of linearity index\non log contrast")
    ax.text(-0.42, -5.6, "below 0: deeper\ncollapse at higher\ncontrast", ha="left", va="bottom", fontsize=5.5,
            color=fs.INK2, linespacing=1.0)
    return out


# ------------------------------------------------------------------------------------------------------- (c) depth
def depth_table(W: pd.DataFrame, cs) -> pd.DataFrame:
    X = W[W.c.isin(cs)].copy()
    X["epoch"] = np.select([X.w0 < 20, X.w0 < 50, X.w0 >= 50], ["early", "transition", "late"], default="")
    X["bin"] = pd.cut(X.depth, [b[0] for b in DEPTH_BINS] + [DEPTH_BINS[-1][1]], right=False,
                      labels=range(len(DEPTH_BINS)))
    return X


def depth_profile(ax, X: pd.DataFrame, title: str):
    per = X.groupby(["penetration", "bin", "epoch"], observed=True).L_w.median().unstack(["epoch", "bin"])
    y = np.array([(a + b - 1) / 2 for a, b in DEPTH_BINS]) * 0.1
    for k, ep in enumerate(EPOCH_COLOR):
        M = per[ep].reindex(columns=range(len(DEPTH_BINS))).to_numpy()
        m, lo, hi = fs.boot_ci(M, stat=np.nanmedian)
        good = np.isfinite(M).sum(axis=0) >= 5
        ax.errorbar(m[good], y[good] + (k - 0.5) * 0.04, xerr=[m[good] - lo[good], hi[good] - m[good]],
                    color=EPOCH_COLOR[ep], lw=1.1, elinewidth=0.6, marker="o", ms=2.6, capsize=0)
    ax.axvline(1, color=fs.COND["linear"], lw=0.8, ls=(0, (3, 1.5)), zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_ylim(3.45, -0.1)                                   # room below the deepest bin for the key
    ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.set_xlim(-1.3, 2.0)
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_xlabel("linearity index")
    ax.set_title(f"both eyes {title}", fontsize=6.5)
    # early vs late per bin, over penetrations
    rows = []
    for b in range(len(DEPTH_BINS)):
        if ("early", b) not in per.columns or ("late", b) not in per.columns:
            continue
        d = (per[("early", b)] - per[("late", b)]).dropna()
        rows.append(dict(set=title, bin_mm=f"{DEPTH_BINS[b][0] / 10:.1f}-{DEPTH_BINS[b][1] / 10:.1f}", n=len(d),
                         early=per[("early", b)].median(), late=per[("late", b)].median(), drop=d.median(),
                         frac_drop=(d > 0).mean(), p=stats.wilcoxon(d).pvalue))
    return pd.DataFrame(rows)


def depth_models(X: pd.DataFrame, title: str) -> pd.DataFrame:
    """Mixed models (as fig3_depth.py): late L and the early - late drop on depth (mm) + depth^2."""
    X = X.assign(monkey=X.penetration.str[7])
    S = X.groupby(["penetration", "monkey", "site_id", "depth", "epoch"]).L_w.median().unstack("epoch").reset_index()
    S["d"] = S.depth / 10
    rows = []
    for var in ("late", "drop"):
        g = S.assign(y=S.late if var == "late" else S.early - S.late).dropna(subset=["y"])
        g = g.assign(y=g.y.clip(-3, 4) if var == "late" else g.y.clip(-5, 5))
        for lab, h in (("all", g), ("B", g[g.monkey == "B"]), ("J", g[g.monkey == "J"])):
            m = smf.mixedlm("y ~ d + I(d**2)", h, groups="penetration").fit()
            q = m.params["I(d ** 2)"]
            rows.append(dict(set=title, variable=var, monkey=lab, n_pen=h.penetration.nunique(),
                             linear=m.params["d"], p_linear=m.pvalues["d"], quadratic=q,
                             p_quadratic=m.pvalues["I(d ** 2)"],
                             peak_mm=-m.params["d"] / (2 * q) if q < 0 else np.nan))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------- (d, e) the sink
def sink_csd(ax, cax, Z, S):
    M = np.nanmean(Z["maps"], axis=0)
    off, t = Z["offsets"] * 0.1, Z["t_map"]
    lim = np.nanpercentile(np.abs(M[:, (t > 20) & (t < 150)]), 99)
    dt, dz = t[1] - t[0], off[1] - off[0]
    im = ax.imshow(M, cmap=fs.CSD_CMAP, vmin=-lim, vmax=lim, interpolation="bilinear", aspect="auto", origin="lower",
                   extent=(t[0] - dt / 2, t[-1] + dt / 2, off[0] - dz / 2, off[-1] + dz / 2), rasterized=True)
    ax.set_ylim(off.max(), off.min())
    ax.set_xlim(-10, 150)
    l45 = np.nanmedian(S.sink_depth_l45) * 0.1                # L4/5 lies this far below the sink (median)
    ax.axhline(l45, color=fs.INK, lw=0.7, ls=(0, (3, 2)))
    ax.text(147, l45 - 0.03, "L4/5 (median)", ha="right", va="bottom", fontsize=5.5, color=fs.INK,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.plot([-10], [0], marker=">", ms=4, color=fs.INK, clip_on=False, zorder=5)
    ax.text(-6, -0.05, "earliest sink", ha="left", va="bottom", fontsize=5.5, color=fs.INK,
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to earliest sink (mm)")
    ax.text(0.97, 0.03, f"{Z['maps'].shape[0]} penetrations", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.5, color=fs.INK2, bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    cb = ax.figure.colorbar(im, cax=cax, ticks=[-lim, 0, lim])
    cb.ax.set_yticklabels(["sink", "0", "source"])
    cb.ax.tick_params(length=0, labelsize=5.5)
    cb.outline.set_visible(False)
    return l45


def rescaled_sink(Z):
    """Sink current per penetration scaled by its mean absolute monocular current over 0-150 ms after sink onset,
    then by one constant so the population mean of the dominant eye at 90% peaks at 1."""
    X, rel = Z["csd"], Z["rel"]
    w = (rel >= 0) & (rel <= 150)
    mono = np.stack([X[:, 5, 0], X[:, 4, 0], X[:, 0, 4]], axis=1)
    sc = np.nanmean(np.abs(mono[:, :, w]), axis=(1, 2))
    Y = X / sc[:, None, None, None]
    return Y / np.nanmax(np.nanmean(Y[:, 5, 0], axis=0)[w])


def sink_tests(Y, rel) -> pd.DataFrame:
    rows = []
    for a, b in SINK_WINDOWS:
        s = (rel >= a) & (rel < b)
        m = {k: np.nanmean(Y[:, i, j][:, s], axis=1) for k, (i, j) in
             dict(d90=(5, 0), d45=(4, 0), o45=(0, 4), b45=(4, 4)).items()}
        for lab, x, ref in (("both - dominant 45%", "b45", "d45"), ("both - dominant 90%", "b45", "d90"),
                            ("dominant 90% - dominant 45%", "d90", "d45")):
            d = m[x] - m[ref]
            rows.append(dict(window=f"{a}-{b} ms", comparison=lab, n=len(d), n_pos=int((d > 0).sum()),
                             mean_diff=d.mean(), median_diff=np.median(d), p=stats.wilcoxon(d).pvalue))
    return pd.DataFrame(rows)


def sink_traces(ax, Y, rel, T):
    sel = (rel >= -10) & (rel <= 150)
    items = [("dominant eye, 90%  (≈ linear sum)", fs.COND["linear"], (0, (3, 1.5)), Y[:, 5, 0]),
             ("dominant eye, 45%", fs.COND["dominant"], "-", Y[:, 4, 0]),
             ("other eye, 45%", fs.COND["other"], "-", Y[:, 0, 4]),
             ("both eyes, 45% each", fs.COND["both"], "-", Y[:, 4, 4])]
    for lab, col, ls, A in items:
        fs.ci_line(ax, rel[sel], A[:, sel], col, ls=ls, lw=1.2 if col == fs.INK else 1.0, alpha=0.09)
    fs.zero_line(ax)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-10, 150)
    ax.set_ylim(-0.25, 1.5)
    ax.set_xlabel("time from sink onset (ms)")
    ax.set_ylabel("current at the earliest sink\n(dominant eye 90%, peak = 1)")
    fs.line_key(ax, [(lab, col, ls) for lab, col, ls, _ in items], x=0.45, y=0.95)
    ax.text(10, 1.47, "first\n20 ms", ha="center", va="top", fontsize=5.5, color=fs.INK2, linespacing=0.95)
    r = T[T.comparison == "both - dominant 45%"].set_index("window")
    lines = [f"{w.replace('-', '–')}: {int(r.loc[w, 'n_pos'])}/{int(r.loc[w, 'n'])} penetrations, "
             f"{fs.fmt_p(r.loc[w, 'p'])}" for w in r.index[:2]]
    ax.text(0.45, 0.62, "both eyes > dominant eye at 45%\n" + "\n".join(lines), transform=ax.transAxes, ha="left",
            va="top", fontsize=5.5, color=fs.INK, linespacing=1.1)


def main():
    fs.apply()
    pd.set_option("display.width", 200)
    R = pd.read_parquet(out_root() / "binocular" / "normalization_timing.parquet")
    Xr = pd.read_parquet(out_root() / "binocular" / "normalization_timing_cross.parquet")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))

    per, per_x = per_recording(R), per_recording(Xr)
    nch = Xr[Xr.split == "all"].groupby("name").channel.nunique()
    per_x = per_x[per_x.index.get_level_values("name").isin(nch.index[nch >= MIN_CROSS_CHANNELS])]
    tr, tr_x = transition(per), transition(per_x)
    sl, sl_x = slopes(tr), slopes(tr_x)
    Y = rescaled_sink(Z)
    T = sink_tests(Y, Z["rel"])

    fig = fs.figure("double", 118)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.0, 0.95], hspace=0.04)
    ta = top.subplots(1, 4, gridspec_kw=dict(width_ratios=[1.0, 1.05, 0.78, 0.78]))
    med, n_clip_a = transition_panel(ta[0], tr)
    sres = slope_panel(ta[1], sl, sl_x)
    depth_rows, model_rows = [], []
    for k, (title, cs) in enumerate(DEPTH_SETS.items()):
        X = depth_table(W, cs)
        depth_rows.append(depth_profile(ta[2 + k], X, title))
        model_rows.append(depth_models(X, title))
    ta[3].sharey(ta[2])
    ta[3].tick_params(labelleft=False)
    ta[2].set_ylabel("depth below top of cortex (mm)")
    fs.line_key(ta[2], [(f"{EPOCHS[ep][0]}–{EPOCHS[ep][1]} ms", EPOCH_COLOR[ep], "-", "o") for ep in EPOCH_COLOR],
                x=0.04, y=0.075, dy=0.055, length=0.12)
    ba = bot.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.15, 0.045, 1.35]))
    l45 = sink_csd(ba[0], ba[1], Z, S)
    sink_traces(ba[2], Y, Z["rel"], T)
    for ax, letter in zip((ta[0], ta[1], ta[2], ba[0], ba[2]), "abcde"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig3_contrast_depth"))

    print("\n(a) transition-window index, median over penetrations:", med.round(2).to_dict(),
          f"n = {len(tr.dropna())}; values beyond the axis: {n_clip_a}")
    print("(b) slopes (negative = deeper collapse at higher contrast):")
    for (lab, x), (n_neg, n, p, m, nclip) in sres.items():
        grp = "held-out runs" if x > 2 else f"penetrations, monkey {lab}"
        print(f"   {grp}: negative in {n_neg}/{n}, one-sided sign test p = {p:.4f}, median {m:+.2f}, "
              f"drawn at the axis edge {nclip}")
    allp = stats.binomtest(int((sl < 0).sum()), len(sl), 0.5, alternative="greater").pvalue
    print(f"   all penetrations: negative in {(sl < 0).sum()}/{len(sl)}, p = {allp:.2g}, median {sl.median():+.2f}")
    print(f"   held-out runs: {len(sl_x)} with >= {MIN_CROSS_CHANNELS} driven channels: {sorted(sl_x.index)}")
    print("   held-out transition index medians:", tr_x.median().round(2).to_dict())
    print("(c) early vs late per depth bin (penetration medians; Wilcoxon):")
    print(pd.concat(depth_rows).round(3).to_string(index=False))
    print("(c) mixed models on depth (mm):")
    print(pd.concat(model_rows).round(4).to_string(index=False))
    print(f"(d) L4/5 median {l45:.2f} mm below the sink; {len(S)} sinks; {Z['maps'].shape[0]} penetrations in the map")
    print("(e) sink current at 45%, per-penetration window means (rescaled; Wilcoxon over penetrations):")
    print(T.round(4).to_string(index=False))
    rel = Z["rel"]
    for a, b in SINK_WINDOWS:
        s = (rel >= a) & (rel < b)
        print(f"   {a}-{b} ms population means: " + ", ".join(
            f"{k} {np.nanmean(Y[:, i, j][:, s]):.2f}" for k, (i, j) in
            dict(dominant90=(5, 0), dominant45=(4, 0), other45=(0, 4), both45=(4, 4)).items()))


if __name__ == "__main__":
    main()
