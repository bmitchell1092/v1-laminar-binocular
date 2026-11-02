"""Figure 3: contrast sets how far binocular summation collapses, the sequence runs at every depth, and the input
current adds at onset (docs/09, figure set for V3; rebuilt 2026-11-02 for readers outside the field, docs/11).

Usage: .venv\\Scripts\\python scripts\\figures\\fig3_contrast_depth.py
(a) what the other eye adds (both eyes minus the dominant eye alone) at 11, 22.5 and 45%, against what linear
    summation of contrast predicts it adds (the dominant eye at c (1 + a) minus at c; each penetration's median eye
    weight a; v1lb.binocular.other_eye_added), onset-aligned population means with 95% penetration-bootstrap CIs,
    5 ms running mean for display (figure_data.py); window means printed with the dip-and-recovery check (docs/09
    F) at each contrast;
(b) linearity index in the transition window (sliding 15 ms windows starting 20-35 ms after each channel's onset,
    i.e. covering 20-50 ms) per penetration (grey) and the median (normalization_timing.py); values beyond the axis
    are drawn at its edge as open triangles; the note gives, per monkey, the penetrations whose index falls with
    contrast (slope on log contrast < 0; one-sided sign tests). The slopes, and the held-out runs of the dual-probe
    days, are drawn in Supplementary Fig S2 (formerly panel b here);
(c) linearity index by depth below the top of cortex (window-specific eye weights), early (0-20 ms after each
    channel's onset) and late (50-200 ms), for dioptic 11 and 22.5% pooled (as the former Fig 4b and fig3_depth.py)
    and for 45%; median of penetration medians per 0.4 mm bin, 95% bootstrap CI over penetrations. Per bin: early vs
    late (Wilcoxon over penetrations); late index and early-late drop on depth (mixed model, linear + quadratic,
    random intercept per penetration, as fig3_depth.py), overall and per monkey;
(d) population CSD aligned to each penetration's earliest current sink (presentations with either eye >= 45%; each
    penetration scaled by its sink's peak, then averaged; sink_data.py), with the median L4/5 boundary; drawn as a
    bilinearly interpolated image (Gouraud-shaded meshes leave saw-tooth seams between channel rows);
(e) what the other eye adds to the current at the earliest sink (-CSD over the sink channel and its two neighbours)
    at 45%, on the sink-onset clock, against what summation predicts, with the sink current's own eye weight
    (v1lb.binocular.eye_weight on its monocular responses 0-50 ms after sink onset, the interval sink_data.py uses
    for dominance). Each penetration is scaled by its mean absolute monocular sink current over 0-150 ms after sink
    onset (dominant 90%, dominant 45%, other 45%; a scale that does not involve the binocular condition), then all
    by one constant so the population mean of the dominant eye at 90% peaks at 1. Means with 95% penetration-
    bootstrap CI. Tests: per penetration, both eyes minus the dominant eye alone (and minus the dominant eye at 90%)
    in the windows of Fig 2 (0-20, 20-50, 50-200 ms after sink onset); Wilcoxon over penetrations.
Inputs: figure_data.py (population.npz, population_meta.parquet); <out_root>/binocular/normalization_timing.parquet,
k_windows.parquet; <out_root>/figures/data/sink.npz (sink_data.py); docs/qc/dichoptic_csd_sinks.csv, v4_F.csv.
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
from v1lb.binocular import CONTRASTS, eye_weight, other_eye_added  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
DIOPTIC = [0.11, 0.225, 0.45]
EPOCHS = {"early": (0, 20), "late": (50, 200)}
EPOCH_COLOR = {"early": fs.MUTED, "late": fs.INK}
DEPTH_BINS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 32)]     # contacts below the top
DEPTH_SETS = {"11 + 22.5%": [0.11, 0.225], "45%": [0.45]}
SINK_WINDOWS = [(0, 20), (20, 50), (50, 200)]                                      # ms after sink onset (Fig 2)
F_WINDOWS = {"0-20 ms": (0, 20), "30-100 ms": (30, 100), "130-200 ms": (130, 200)}
TR_LIM = (-2.0, 3.0)                                                               # (b) axis limits
H_MM = 128
POS = {"a1": (15, 80, 29, 37), "a2": (48, 80, 29, 37), "a3": (81, 80, 29, 37), "b": (128, 80, 43, 37),
       "c1": (15, 11, 20, 50), "c2": (38, 11, 20, 50), "d": (73, 11, 38, 50), "dcb": (112.5, 11, 2.2, 50),
       "e": (140, 11, 32, 50)}                                                     # (x, y, w, h) in mm


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


# ------------------------------------------------------------------------------------------------ (a, b) contrast
def contrast_added(axs, pop, a):
    O, rel = pop["onset"], pop["rel"]
    sel = (rel >= -20) & (rel <= 200)
    rows = []
    for k, (ax, c) in enumerate(zip(axs, DIOPTIC)):
        A, P = other_eye_added(O, a, c)
        fs.added_panel(ax, rel[sel], fs.smooth(A)[:, sel], fs.smooth(P)[:, sel], label_added="measured",
                       key=k == 0, key_xy=(0.36, 0.95), zero_label=k == 2)
        ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
        ax.set_xlim(-20, 200)
        ax.set_ylim(-0.1, 0.34)
        ax.set_xticks([0, 100, 200])
        ax.set_title(f"both eyes {fs.contrast_label(c)}", fontsize=6.5)
        if k:
            ax.tick_params(labelleft=False)
        for name, (w0, w1) in F_WINDOWS.items():
            s = (rel >= w0) & (rel < w1)
            rows.append(dict(c=c, window=name, added=float(A[:, s].mean()), predicted=float(P[:, s].mean())))
    axs[0].set_ylabel("added by the other eye (MUA)\n= both eyes − dominant eye alone")
    axs[1].set_xlabel("time from response onset (ms)")
    return pd.DataFrame(rows)


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


def slope_counts(sl: pd.Series) -> dict:
    """Per monkey: penetrations whose transition index falls with contrast, one-sided sign test."""
    out = {}
    for mk in ("B", "J"):
        v = sl[sl.index.str[7] == mk]
        n_neg = int((v < 0).sum())
        out[mk] = (n_neg, len(v), stats.binomtest(n_neg, len(v), 0.5, alternative="greater").pvalue, v.median())
    return out


def transition_panel(ax, tr: pd.DataFrame, counts: dict):
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
    ax.axhline(1, color=fs.COND["linear"], lw=0.8, ls=fs.SUM_LS, zorder=0)
    fs.zero_line(ax)
    fs.log_contrast_axis(ax, DIOPTIC)
    ax.set_xticklabels(["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_xlim(0.09, 0.55)
    ax.set_ylim(TR_LIM[0] - 0.1, TR_LIM[1] + 0.1)
    ax.set_ylabel(fs.LINEARITY)
    ax.set_title("20–50 ms after onset", fontsize=6.5)
    (nb, Nb, _, _), (nj, Nj, _, _) = counts["B"], counts["J"]
    ax.text(0.03, 0.03, f"lower at higher contrast in\n{nb + nj}/{Nb + Nj} penetrations\n(monkey B {nb}/{Nb}, "
            f"J {nj}/{Nj})", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5, color=fs.INK,
            linespacing=1.0, zorder=4, bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    n_clip = int(((v < lo) | (v > hi)).to_numpy().sum())
    for y, lab in ((1.0, "1: contrasts add"), (0.0, "0: dominant eye alone")):
        ax.text(0.093, y + 0.04, lab, ha="left", va="bottom", fontsize=5.5, color=fs.INK2, zorder=4,
                bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    return med, n_clip


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
    ax.axvline(1, color=fs.COND["linear"], lw=0.8, ls=fs.SUM_LS, zorder=0)
    ax.axvline(0, color=fs.RULE, lw=0.6, zorder=0)
    ax.set_ylim(3.45, -0.1)                                   # room below the deepest bin for the key
    ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.set_xlim(-1.3, 2.0)
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_title(f"both eyes\n{title}", fontsize=6.5, linespacing=1.0)
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
    ax.text(-6, -0.05, "earliest sink:\ninput arrives first", ha="left", va="bottom", fontsize=5.5, color=fs.INK,
            linespacing=1.0, bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.85))
    ax.set_xticks([0, 50, 100, 150])
    ax.set_xlabel("time from stimulus onset (ms)")
    ax.set_ylabel("depth relative to earliest sink (mm)")
    ax.set_title(f"current-source density, {Z['maps'].shape[0]} penetrations", fontsize=6.5)
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


def sink_weights(Y, rel) -> np.ndarray:
    """Eye weight of the sink current per penetration: monocular responses 0-50 ms after sink onset."""
    w = (rel >= 0) & (rel < 50)
    return np.array([eye_weight(np.nanmean(y[:, 0][:, w], axis=1), np.nanmean(y[0, :][:, w], axis=1)) for y in Y])


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


def sink_added(ax, Y, rel, a, T):
    sel = (rel >= -10) & (rel <= 150)
    A, P = other_eye_added(Y, a, 0.45)
    fs.added_panel(ax, rel[sel], fs.smooth(A)[:, sel], fs.smooth(P)[:, sel], label_added="measured",
                   key_xy=(0.30, 0.95), zero_label=False, alpha=0.12)
    ax.axvspan(0, 20, color=fs.RULE, alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(-10, 150)
    ax.set_xticks([0, 50, 100, 150])
    ax.set_ylim(-0.3, 0.6)
    ax.set_xlabel("time from sink onset (ms)")
    ax.set_ylabel("current added by the other eye\n(earliest sink, dominant eye 90% peak = 1)")
    ax.set_title("both eyes 45%", fontsize=6.5)
    r = T[(T.comparison == "both - dominant 45%") & (T.window == "0-20 ms")].iloc[0]
    ax.text(0.97, 0.03, f"added in the first 20 ms:\n{int(r.n_pos)}/{int(r.n)} penetrations, {fs.fmt_p(r.p)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.05)
    return A, P


def main():
    fs.apply()
    pd.set_option("display.width", 200)
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    a = meta.loc[[str(p) for p in pop["pens"]], "a"].to_numpy()
    R = pd.read_parquet(out_root() / "binocular" / "normalization_timing.parquet")
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    S = pd.read_csv(REPO / "docs" / "qc" / "dichoptic_csd_sinks.csv")
    S = S[S.usable]
    Z = dict(np.load(out_root() / "figures" / "data" / "sink.npz"))
    Fq = pd.read_csv(REPO / "docs" / "qc" / "v4_F.csv")

    tr = transition(per_recording(R))
    sl = slopes(tr)
    counts = slope_counts(sl)
    Y = rescaled_sink(Z)
    a_sink = sink_weights(Y, Z["rel"])
    T = sink_tests(Y, Z["rel"])

    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")                     # axes placed in mm: in-panel keys stay inside their panels
    ax_a = [axes_mm(fig, *POS[n]) for n in ("a1", "a2", "a3")]
    ax_b, ax_c1, ax_c2, ax_d, ax_dcb, ax_e = (axes_mm(fig, *POS[n]) for n in ("b", "c1", "c2", "d", "dcb", "e"))
    Wa = contrast_added(ax_a, pop, a)
    med, n_clip_b = transition_panel(ax_b, tr, counts)
    depth_rows, model_rows = [], []
    for ax, (title, cs) in zip((ax_c1, ax_c2), DEPTH_SETS.items()):
        X = depth_table(W, cs)
        depth_rows.append(depth_profile(ax, X, title))
        model_rows.append(depth_models(X, title))
    ax_c2.sharey(ax_c1)
    ax_c2.tick_params(labelleft=False)
    ax_c1.set_ylabel("depth below top of cortex (mm)")
    fig.text((POS["c1"][0] + POS["c2"][0] + POS["c2"][2]) / 2 / 174, 2.2 / H_MM,
             "linearity index (share of the\nother eye's contrast added)", ha="center", va="bottom", fontsize=7,
             color=fs.INK, linespacing=1.0)
    fs.line_key(ax_c1, [(f"{EPOCHS[ep][0]}–{EPOCHS[ep][1]} ms", EPOCH_COLOR[ep], "-", "o") for ep in EPOCH_COLOR],
                x=0.04, y=0.075, dy=0.055, length=0.16)
    l45 = sink_csd(ax_d, ax_dcb, Z, S)
    sink_added(ax_e, Y, Z["rel"], a_sink, T)
    for ax, letter in zip((ax_a[0], ax_b, ax_c1, ax_d, ax_e), "abcde"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig3_contrast_depth"))

    print("(a) population means of what the other eye adds and of the summation prediction:")
    print(Wa.round(4).to_string(index=False))
    print("    per-penetration check (docs/09 F) at 11 and 22.5%:")
    print(Fq[~np.isclose(Fq.c, 0.45)][["test", "c", "n", "median", "n_in_direction", "p", "p_odd", "p_even"]]
          .round(4).to_string(index=False))
    print("\n(b) transition-window index, median over penetrations:", med.round(2).to_dict(),
          f"n = {len(tr.dropna())}; values beyond the axis: {n_clip_b}")
    for mk, (n_neg, n, p, m) in counts.items():
        print(f"    monkey {mk}: slope on log contrast negative in {n_neg}/{n}, one-sided sign test p = {p:.4f}, "
              f"median {m:+.2f}")
    allp = stats.binomtest(int((sl < 0).sum()), len(sl), 0.5, alternative="greater").pvalue
    print(f"    all penetrations: negative in {(sl < 0).sum()}/{len(sl)}, p = {allp:.2g}, median {sl.median():+.2f}")
    print("(c) early vs late per depth bin (penetration medians; Wilcoxon):")
    print(pd.concat(depth_rows).round(3).to_string(index=False))
    print("(c) mixed models on depth (mm):")
    print(pd.concat(model_rows).round(4).to_string(index=False))
    print(f"(d) L4/5 median {l45:.2f} mm below the sink; {len(S)} sinks; {Z['maps'].shape[0]} penetrations in the map")
    print(f"(e) sink-current eye weight: median {np.median(a_sink):.2f} (IQR {np.percentile(a_sink, 25):.2f}-"
          f"{np.percentile(a_sink, 75):.2f}); MUA eye weight median {np.median(a):.2f}")
    print("(e) sink current at 45%, per-penetration window means (rescaled; Wilcoxon over penetrations):")
    print(T.round(4).to_string(index=False))
    rel = Z["rel"]
    for w0, w1 in SINK_WINDOWS:
        s = (rel >= w0) & (rel < w1)
        print(f"   {w0}-{w1} ms population means: " + ", ".join(
            f"{k} {np.nanmean(Y[:, i, j][:, s]):.2f}" for k, (i, j) in
            dict(dominant90=(5, 0), dominant45=(4, 0), other45=(0, 4), both45=(4, 4)).items()))


if __name__ == "__main__":
    main()
