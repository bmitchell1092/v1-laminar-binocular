"""Figure 5: binocular normalization behaves like a delayed conductance (docs/09, model tests 1, 4 and 5; merges the
former Figs 5 and 6, whose two-stage, Ding-Sperling and offset panels moved to figS_model_variants.py and
figS_offset.py).

Usage: .venv\\Scripts\\python scripts\\figures\\fig5_models.py
Three normalization models fitted per penetration to the channel-averaged 6 x 6 time courses (0-219 ms after
stimulus onset; v1lb.dynorm):
* fixed-delay ('full'): R = g E^n / (sigma^n + P^n), drive E = (c_D + a c_N) u(t), pool P = (c_D + w c_N) v(t)
  delayed by a low-pass tau2;
* conductance ('cond'): R relaxes toward the fixed-delay value with time constant tau_m / (1 + (P / sigma)^n), so the
  pool both divides the response and speeds it up; tau_m -> 0 gives the fixed-delay model;
* recurrent ('org'): binocular ORGaNICs, whose normalization signal is built recurrently from the responses.
(a) the conductance model and its relatives (also the colour key);
(b) held-out (odd/even) error of static vs fixed-delay normalization per penetration (fig5_model.py);
(c) pool weight w vs excitatory weight a of the other eye in the conductance model (fig5_cond_model.py);
(d) dominant eye alone at 11% and 90%, in units of each channel's 90% peak: data (mean over penetrations, 95%
    bootstrap CI) and the mean of each model's fit (rise_noise_check.py, dominant_traces.npz);
(e) 10-90% rise time of the population-mean response vs contrast: data (95% penetration-bootstrap CI; upper ends
    beyond the axis are drawn at its edge, open triangles) and model fits (rise_noise_check.py). Per-penetration
    rise times are not shown: at <= 11% they are dominated by noise (docs/09 correction of 2026-10-31);
(f) % less held-out error per penetration for three pairings; one-sided Wilcoxon on log error ratios; values
    beyond the axis are drawn at its edge (open triangles);
(g) late linearity index (windows starting 50-100 ms after onset) of each model's fit minus the data, median with
    90% bootstrap CI, against the +-0.25 equivalence band; at 11% the data's own noise moves the index by +-0.4-0.7
    (docs/09), so that contrast is drawn open (reported, not tested).
Inputs: <out_root>/binocular/fig5_model.parquet, fig5_cond.parquet, fig5_cond_L.parquet, fig5_org.parquet,
fig5_org_L.parquet, dominant_traces.npz; docs/qc/rise_noise.csv.
Writes <out_root>/figures/fig5_models.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
RISE_C = [0.055, 0.11, 0.225, 0.45, 0.9]
DIOPTIC = [0.11, 0.225, 0.45]
LATE = (50, 75, 100)                   # window starts, as registered
MARGIN = 0.25                          # registered equivalence margin
GAIN_LIM = (-15.0, 25.0)               # % axis limits in (f)
RISE_MAX = 45.0                        # ms axis limit in (e)
MODELS = {"full": ("fixed-delay", fs.MODEL["fixed-delay"]), "cond": ("conductance", fs.MODEL["conductance"]),
          "org": ("recurrent", fs.MODEL["recurrent"])}
ORDER = ("full", "cond", "org")


def schematic(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec=fs.INK2, lw=0.6)
    arrow = dict(arrowstyle="-|>", lw=0.6, color=fs.INK2, mutation_scale=6, shrinkA=1, shrinkB=1)

    def node(x, y, s, **kw):
        return ax.text(x, y, s, ha="center", va="center", fontsize=6, color=fs.INK, bbox=kw.pop("bbox", box), **kw)

    def arr(src, dst):
        """Arrow from a node centre (or a point) to a node centre, clipped at the boxes."""
        is_node = hasattr(src, "get_bbox_patch")
        ax.annotate("", xy=dst.get_position(), xytext=src.get_position() if is_node else src,
                    arrowprops=dict(arrow, patchA=src.get_bbox_patch() if is_node else None,
                                    patchB=dst.get_bbox_patch()))

    yt, yb = 0.88, 0.56
    ax.text(0.0, yt, "dominant\neye $c_D$", ha="left", va="center", fontsize=6, color=fs.COND["dominant"],
            linespacing=0.95)
    ax.text(0.0, yb, "other\neye $c_N$", ha="left", va="center", fontsize=6, color=fs.COND["other"],
            linespacing=0.95)
    drive = node(0.37, yt, "drive  $E = c_D + a\\,c_N$")
    pool = node(0.30, yb, "pool  $c_D + w\\,c_N$")
    delay = node(0.56, yb, "delay $\\tau_2$")
    resp = node(0.81, 0.72, "$R \\rightarrow E^n / (\\sigma^n + P^n)$\nwith time constant\n"
                "$\\tau_m / (1 + (P/\\sigma)^n)$", linespacing=1.3, bbox=dict(box, ec=fs.MODEL["conductance"], lw=0.9))
    for y0, dst in ((yt, drive), (yt - 0.04, pool), (yb, pool), (yb + 0.04, drive)):
        arr((0.155, y0), dst)
    arr(pool, delay)
    arr(drive, resp)
    arr(delay, resp)
    ax.text(0.64, yb - 0.1, "$P(t)$ divides the response\nand sets its speed", fontsize=5.5, color=fs.INK2,
            ha="left", va="top", linespacing=1.0)
    items = [("data", fs.MODEL["data"], "-"),
             ("fixed-delay: $\\tau_m \\rightarrow 0$, R equals its target at every moment", MODELS["full"][1], "-"),
             ("conductance: $\\tau_m$ fitted", MODELS["cond"][1], "-"),
             ("recurrent (ORGaNICs): the pool is built from the responses", MODELS["org"][1], "-")]
    fs.line_key(ax, items, x=0.0, y=0.27, dy=0.085, length=0.06)


def static_vs_delayed(ax, M):
    for mk, g in M.groupby("monkey"):
        ax.scatter(g.cv_static, g.cv_full, s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                   lw=0.4, zorder=3)
    lo = min(M.cv_full.min(), M.cv_static.min()) / 1.5
    hi = max(M.cv_full.max(), M.cv_static.max()) * 1.5
    ax.plot([lo, hi], [lo, hi], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("held-out error, static")
    ax.set_ylabel("held-out error, fixed-delay")
    worst = M.loc[(M.cv_full / M.cv_static).idxmax()]
    ax.annotate(f"{worst.penetration[:6]} (few trials)", xy=(worst.cv_static, worst.cv_full), xytext=(-4, 0),
                textcoords="offset points", ha="right", va="center", fontsize=5.5, color=fs.INK2)
    better = int((M.cv_full < M.cv_static).sum())
    p = stats.wilcoxon(M.cv_full, M.cv_static).pvalue
    ax.text(0.97, 0.06, f"delayed better\n{better}/{len(M)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.0)
    for k, (mk, lab) in enumerate((("o", "monkey B"), ("s", "monkey J"))):
        y = 0.42 - k * 0.08
        ax.scatter([0.62], [y], s=9, marker=mk, facecolor=fs.INK2, edgecolor=fs.INK, lw=0.4, transform=ax.transAxes)
        ax.text(0.67, y, lab, va="center", fontsize=5.5, color=fs.INK, transform=ax.transAxes)
    return better, p


def weights(ax, F):
    a, w = F.cond_p_a, F.cond_p_w
    for mk, idx in F.groupby("monkey").groups.items():
        ax.scatter(a[idx], w[idx], s=9, marker=fs.MONKEY_MARKER[mk], facecolor=fs.MODEL["conductance"],
                   edgecolor="white", lw=0.3, zorder=3)
    top = max(1.6, float(max(a.max(), w.max())) * 1.06)
    ax.plot([0, top], [0, top], color=fs.MUTED, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal")
    ax.set_xlabel("excitatory weight a")
    ax.set_ylabel("pool weight w")
    n = int((w > a).sum())
    p = stats.wilcoxon(w, a).pvalue
    ax.text(0.04, 0.96, f"w > a\n{n}/{len(F)}, {fs.fmt_p(p)}", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.5, color=fs.INK, linespacing=1.0)
    ax.text(0.96, 0.04, "other eye's weight\nin the conductance model", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=5.5, color=fs.INK2, linespacing=1.0)
    return n, p


def time_courses(axs, tr):
    A, src, t = tr["A"], list(tr["sources"]), tr["t"]
    for ax, c in zip(axs, (0.11, 0.9)):
        i = list(CONTRASTS).index(c)
        for m in ORDER:
            ax.plot(t, A[:, src.index(m), i].mean(axis=0), color=MODELS[m][1], lw=1.0, zorder=3)
        fs.ci_line(ax, t, A[:, src.index("data"), i], fs.MODEL["data"], lw=1.2, alpha=0.15, zorder=4)
        fs.zero_line(ax)
        ax.set_xlim(0, 150)
        ax.text(0.97, 0.94, f"{fs.contrast_label(c)}, one eye", transform=ax.transAxes, ha="right", va="top",
                fontsize=5.5, color=fs.INK2)
    axs[0].set_ylim(-0.03, 0.42)
    axs[1].set_ylim(-0.08, 1.1)
    axs[0].tick_params(labelbottom=False)
    axs[1].set_xlabel("time from stimulus onset (ms)")
    axs[0].set_ylabel("response (90% peak = 1)", y=0.0)
    return len(A)


def rise(ax, Q):
    q = Q[Q.measure == "population mean"].set_index(["source", "c"])
    x = np.arange(len(RISE_C))
    for k in ORDER:
        ax.plot(x, q.loc[k].reindex(RISE_C).rise.to_numpy(), color=MODELS[k][1], marker="o", ms=2.8, lw=1.0, zorder=3)
    d = q.loc["data"].reindex(RISE_C)
    m, lo, hi = d.rise.to_numpy(), d.lo.to_numpy(), np.minimum(d.hi.to_numpy(), RISE_MAX)
    ax.errorbar(x, m, [m - lo, hi - m], color=fs.MODEL["data"], marker="o", ms=3, lw=1.2, elinewidth=0.7, capsize=0,
                zorder=4)
    off = d.hi.to_numpy() > RISE_MAX
    ax.scatter(x[off], np.full(off.sum(), RISE_MAX), marker="^", s=10, facecolor="white", edgecolor=fs.INK, lw=0.5,
               zorder=5, clip_on=False)
    fs.line_key(ax, [("data", fs.MODEL["data"], "-", "o")] + [(MODELS[k][0], MODELS[k][1], "-", "o") for k in ORDER],
                x=0.5, y=0.97, dy=0.07, length=0.1)
    ax.set_xticks(x, ["5.5", "11", "22.5", "45", "90"])
    ax.set_xlim(-0.4, len(RISE_C) - 0.6)
    ax.set_xlabel("contrast, one eye (%)")
    ax.set_ylabel("rise time of the mean\nresponse, 10–90% (ms)")
    ax.set_ylim(0, RISE_MAX)
    return q.rise.unstack("c")


def held_out(ax, F):
    pairs = (("org", "full"), ("cond", "full"), ("cond", "org"))
    rows = []
    for x, (ch, ref) in enumerate(pairs):
        r = np.log(F[f"cv_{ref}"] / F[f"cv_{ch}"])
        pct = 100 * (1 - np.exp(-r))
        rng = np.random.default_rng(x)
        for mk, idx in F.groupby("monkey").groups.items():
            v = pct[idx].to_numpy()
            xs = x + rng.uniform(-0.15, 0.15, len(v))
            inside = (v >= GAIN_LIM[0]) & (v <= GAIN_LIM[1])
            ax.scatter(xs[inside], v[inside], s=8, marker=fs.MONKEY_MARKER[mk], facecolor=fs.INK2, edgecolor=fs.INK,
                       lw=0.4, zorder=3)
            for edge, mkr in ((GAIN_LIM[1], "^"), (GAIN_LIM[0], "v")):
                out = (v > edge) if edge > 0 else (v < edge)
                ax.scatter(xs[out], np.full(out.sum(), edge), s=10, marker=mkr, facecolor="white", edgecolor=fs.INK,
                           lw=0.5, zorder=3, clip_on=False)
        ax.plot([x - 0.25, x + 0.25], [pct.median()] * 2, color=MODELS[ch][1], lw=1.6, zorder=4)
        p = stats.wilcoxon(r, alternative="greater").pvalue
        ax.text(x, GAIN_LIM[1] + 2, f"{int((r > 0).sum())}/{len(r)}\n{fs.fmt_p(p)}", ha="center", va="bottom",
                fontsize=5.5, color=fs.INK, linespacing=1.0)
        rows.append(dict(pair=f"{ch} vs {ref}", n=len(r), better=int((r > 0).sum()), median_pct=pct.median(), p=p,
                         clipped=int(((pct < GAIN_LIM[0]) | (pct > GAIN_LIM[1])).sum())))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks(range(3), ["recurrent\nvs\nfixed-delay", "conductance\nvs\nfixed-delay", "conductance\nvs\nrecurrent"],
                  fontsize=5.2)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(*GAIN_LIM)
    ax.set_ylabel("% less held-out error\n(first model vs second)")
    return pd.DataFrame(rows)


def collapse(ax, L):
    late = L[L.w0.isin(LATE)].groupby(["penetration", "source", "c"]).L.mean().unstack("source")
    rng = np.random.default_rng(0)
    ax.axhspan(-MARGIN, MARGIN, color=fs.RULE, alpha=0.6, lw=0, zorder=0)
    out = []
    for k, c in enumerate(DIOPTIC):
        g = late.xs(c, level="c")
        for j, m in enumerate(ORDER):
            d = (g[m] - g["data"]).dropna()
            boot = [np.median(rng.choice(d.values, len(d))) for _ in range(5000)]
            lo, hi = np.percentile(boot, [5, 95])
            x = k + (j - 1) * 0.24
            tested = c != 0.11
            ax.errorbar([x], [d.median()], [[d.median() - lo], [hi - d.median()]], color=MODELS[m][1], marker="o",
                        mfc=MODELS[m][1] if tested else "white", ms=3.2, elinewidth=0.8, capsize=0, lw=0)
            out.append(dict(c=c, model=m, n=len(d), median=d.median(), lo=lo, hi=hi,
                            inside=(lo > -MARGIN) and (hi < MARGIN), tested=tested))
    fs.zero_line(ax, color=fs.INK2)
    ax.set_xticks(range(3), ["11", "22.5", "45"])
    ax.set_xlabel("contrast, both eyes (%)")
    ax.set_ylabel("late linearity index,\nmodel − data")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.8, 0.8)
    ax.text(0, -0.76, "noise-\nlimited", ha="center", va="bottom", fontsize=5.5, color=fs.INK2, linespacing=1.0)
    ax.text(-0.45, -MARGIN - 0.03, "equivalence\nmargin", ha="left", va="top", fontsize=5.5, color=fs.INK2,
            linespacing=1.0)
    ax.text(-0.45, 0.78, "model too linear", fontsize=5.5, color=fs.INK2, va="top")
    return pd.DataFrame(out)


def main():
    fs.apply()
    b = out_root() / "binocular"
    M = pd.read_parquet(b / "fig5_model.parquet")
    F = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    O = pd.read_parquet(b / "fig5_org.parquet").set_index("penetration")
    F = F.join(O[[c for c in O.columns if c.startswith(("cv_org", "org_"))]])
    L = pd.concat([pd.read_parquet(b / "fig5_cond_L.parquet"),
                   pd.read_parquet(b / "fig5_org_L.parquet").query("source == 'org'")])
    tr = dict(np.load(b / "dominant_traces.npz"))
    Q = pd.read_csv(REPO / "docs" / "qc" / "rise_noise.csv")

    fig = fs.figure("double", 122)
    top, bot = fig.subfigures(2, 1, height_ratios=[1, 1.1], hspace=0.04)
    ta = top.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.75, 1, 1]))
    schematic(ta[0])
    nb, pb = static_vs_delayed(ta[1], M)
    nw, pw = weights(ta[2], F)
    gs = bot.add_gridspec(2, 4, width_ratios=[1.0, 0.95, 1.35, 1.0])
    ax_d = [bot.add_subplot(gs[0, 0]), bot.add_subplot(gs[1, 0])]
    ax_e, ax_f, ax_g = bot.add_subplot(gs[:, 1]), bot.add_subplot(gs[:, 2]), bot.add_subplot(gs[:, 3])
    n_tc = time_courses(ax_d, tr)
    R = rise(ax_e, Q)
    H = held_out(ax_f, F)
    K = collapse(ax_g, L)
    for ax, letter in zip((ta[0], ta[1], ta[2], ax_d[0], ax_e, ax_f, ax_g), "abcdefg"):
        fs.panel_label(ax, letter)
    print(fs.save(fig, "fig5_models"))
    print(f"static vs delayed: delayed better {nb}/{len(M)} (p {pb:.2g}); conductance w > a {nw}/{len(F)} "
          f"(p {pw:.2g}), median w {F.cond_p_w.median():.2f}, a {F.cond_p_a.median():.2f}; time courses n = {n_tc}")
    print("rise of the population-mean response (ms):")
    print(R.round(1).to_string())
    print(H.round(4).to_string(index=False))
    print(K.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
