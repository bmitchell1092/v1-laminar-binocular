"""V3 analysis plan, item C: iso-response contours across the 6 x 6 grid, early vs late (docs/09, registered
2026-10-31, commit 9f8dcfc).

Usage: .venv\\Scripts\\python scripts\\v3_iso_contours.py [n_workers]
Per penetration, the channel-average readout (fig5_model.population, all presentations; dominant eye on axis 0) is
averaged in windows aligned to the response onset (half-max of the dominant-eye response at 90%, as in
fig5_model.windowed_L): early 0-20 ms, mid 20-50 ms, late 50-200 ms after onset. In each window:
* the eye weight a is fitted to the two monocular contrast responses (v1lb.binocular.eye_weight, as in windowed_L);
* for target levels equal to the dominant eye's response at c* = 11, 22.5, 45% (running-maximum envelope), each
  other-eye contrast c_N gives one contour point: the dominant-eye contrast c_D at which the row R(., c_N) reaches the
  target, on the row's running-maximum envelope, interpolated linearly between c_D = 0 and 5.5% and in log contrast
  above (as for equivalent contrast). If the target is at or below R(0, c_N), the contour has already met the
  other-eye axis (c_D = 0, "below axis"); if it exceeds R(90%, c_N), c_D is set to 90% as a lower bound ("censored");
* each contour is summarized by the Minkowski exponent m in (c_D/c*)^m + (a c_N/c*)^m = 1 (m = 1 linear summation,
  large m = dominant eye alone; m on a log grid over 0.5-20, least squares in c_D/c* over the points with
  0 < a c_N/c* < 1), and by the suppression fraction: of the five contour points with c_N > 0, the fraction with
  c_D > c* (beyond the dominant-eye-alone line).
Registered criterion, evaluated at each c*: late m > early m (one-sided Wilcoxon over penetrations, p < .05) with the
early median m within 0.7-1.5.
Population contours: from the population-average window responses (mean over penetrations of each 6 x 6 window
matrix, as in Fig 2b), with 2000-sample penetration-bootstrap 95% CIs.
Writes <out_root>/binocular/v3_C_penetration.parquet (one row per penetration x window x c*),
v3_C_points.parquet (contour points), v3_C_population.parquet, v3_C_population_boot.parquet (bootstrap contour
points, for the figure CIs), v3_C_windows.npz (window matrices [penetration, window, 6, 6]), docs/qc/v3_C.csv and
the diagnostic figure
<out_root>/figures/diagnostics/v3_iso_contours.png; prints the tests.
"""

import glob
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig5_model  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, eye_weight  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
T = np.arange(0, 281)
fig5_model.T = T          # population() takes its time range from the module-level T; extend it to cover the late window
WINDOWS = {"early": (0, 20), "mid": (20, 50), "late": (50, 200)}
CSTAR = (0.11, 0.225, 0.45)
M_GRID = np.geomspace(0.5, 20.0, 400)
N_BOOT = 2000


def monotone(r):
    return np.maximum.accumulate(np.asarray(r, float))


def invert_row(r, target):
    """Dominant-eye contrast at which the row R(c_D, c_N) reaches target; returns (c_D, status)."""
    r = monotone(r)
    if target <= r[0]:
        return 0.0, "below_axis"
    if target > r[-1]:
        return float(CONTRASTS[-1]), "censored"
    if target <= r[1]:
        return float(CONTRASTS[1] * (target - r[0]) / (r[1] - r[0])), "ok"
    lc = np.log(CONTRASTS[1:])
    j = int(np.searchsorted(r[1:], target))
    r0, r1 = r[1:][j - 1], r[1:][j]
    f = (target - r0) / (r1 - r0) if r1 > r0 else 0.0
    return float(np.exp(lc[j - 1] + f * (lc[j] - lc[j - 1]))), "ok"


def fit_m(x, y):
    """Minkowski exponent from contour points (x = c_D/c*, y = a c_N/c*), points with 0 < y < 1."""
    use = (y > 0) & (y < 1) & np.isfinite(x)
    if not use.any():
        return np.nan, 0
    xs, ys = x[use], y[use]
    pred = (1 - ys[None, :] ** M_GRID[:, None]) ** (1 / M_GRID[:, None])
    sse = ((pred - xs[None, :]) ** 2).sum(axis=1)
    return float(M_GRID[int(np.argmin(sse))]), int(use.sum())


def contours(M):
    """All contours of one 6 x 6 window matrix M[c_D, c_N]: (eye weight, {c*: points DataFrame})."""
    a = eye_weight(M[:, 0], M[0, :])
    D = monotone(M[:, 0])
    out = {}
    for cs in CSTAR:
        target = D[list(CONTRASTS).index(cs)]
        rows = []
        for j, cn in enumerate(CONTRASTS):
            cd, status = invert_row(M[:, j], target)
            rows.append(dict(c_N=cn, c_D=cd, x=cd / cs, y=a * cn / cs, status=status))
        out[cs] = pd.DataFrame(rows)
    return a, out


def summarize(P):
    """m, number of fitted points and suppression fraction of one contour (points DataFrame)."""
    m, n_fit = fit_m(P.x.to_numpy(), P.y.to_numpy())
    pos = P[P.c_N > 0]
    return dict(m=m, n_fit=n_fit, suppression=float((pos.x > 1).mean()),
                n_below_axis=int((pos.status == "below_axis").sum()), n_censored=int((pos.status == "censored").sum()),
                n_below_axis_fit=int(((pos.status == "below_axis") & (pos.y < 1)).sum()),
                n_censored_fit=int(((pos.status == "censored") & (pos.y < 1)).sum()))


def work(task):
    pen, driven = task
    pop = fig5_model.population(pen, driven)
    if pop is None:
        return None
    Y, n_ch = pop[0]["all"], pop[1]
    rd = Y[-1, 0]
    idx = np.flatnonzero(rd >= 0.5 * np.nanmax(rd))
    if not idx.size:
        return None
    on = float(T[idx[0]])
    mats = {}
    for w, (w0, w1) in WINDOWS.items():
        sel = (T >= on + w0) & (T < on + w1)
        mats[w] = np.nanmean(Y[:, :, sel], axis=2)
    return dict(penetration=pen, monkey=pen[7], n_channels=n_ch, onset=on,
                late_truncated=bool(on + WINDOWS["late"][1] > T[-1] + 1), mats=mats)


def population_contours(mats):
    """Contours and summaries from the population-average window matrix."""
    a, C = contours(np.mean(mats, axis=0))
    return a, C, {cs: summarize(C[cs]) for cs in CSTAR}


def figure(Pop, Boot):
    fs.apply()
    fig = fs.figure("double", 62)
    axs = fig.subplots(1, 3, sharey=True)
    shade = {"early": fs.MUTED, "mid": fs.INK2, "late": fs.INK}
    ymax = 2.2
    for ax, cs in zip(axs, CSTAR):
        xx = np.linspace(0, 1, 50)
        ax.plot(xx, 1 - xx, color=fs.MUTED, ls=(0, (3, 1.5)), lw=0.8, zorder=1)
        ax.axhline(1, color=fs.RULE, lw=0.8, zorder=0)
        for k, w in enumerate(WINDOWS):
            p = Pop[(Pop.window == w) & (Pop.cstar == cs)].sort_values("c_N")
            b = Boot[(Boot.window == w) & (Boot.cstar == cs)].groupby("c_N").x
            lo, hi = b.quantile(0.025).reindex(p.c_N).to_numpy(), b.quantile(0.975).reindex(p.c_N).to_numpy()
            show = p.y.to_numpy() <= 1.6
            x, y = p.y.to_numpy()[show], np.minimum(p.x.to_numpy()[show], ymax)
            ax.errorbar(x + (k - 1) * 0.015, y, [np.clip(y - np.minimum(lo[show], ymax), 0, None),
                                                 np.clip(np.minimum(hi[show], ymax) - y, 0, None)],
                        color=shade[w], marker="o", ms=3, lw=1.1, elinewidth=0.6, capsize=0, zorder=3,
                        mfc=shade[w])
            cens = (p.status.to_numpy()[show] == "censored") | (p.x.to_numpy()[show] > ymax)
            ax.scatter(x[cens], y[cens], marker="^", s=14, facecolor="white", edgecolor=shade[w], zorder=4)
        ax.set_xlim(0, 1.6)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("other eye, a·c$_N$ / c*")
        ax.set_title(f"c* = {fs.contrast_label(cs)}", loc="left")
    axs[0].set_ylabel("dominant eye, c$_D$ / c*")
    axs[0].text(0.03, 0.6, "linear sum", rotation=-38, fontsize=5.5, color=fs.INK2, transform=axs[0].transData)
    axs[0].text(1.05, 1.03, "dominant eye alone", fontsize=5.5, color=fs.INK2, ha="left", va="bottom")
    fs.line_key(axs[2], [("0–20 ms", fs.MUTED, "-", "o"), ("20–50 ms", fs.INK2, "-", "o"),
                         ("50–200 ms after onset", fs.INK, "-", "o")], x=0.35, y=0.97, dy=0.08, length=0.1)
    d = out_root() / "figures" / "diagnostics"
    d.mkdir(parents=True, exist_ok=True)
    return fs.save(fig, "diagnostics/v3_iso_contours")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    R = pd.concat([pd.read_parquet(f) for f in glob.glob(str(out_root() / "crf" / "*_crf.parquet"))], ignore_index=True)
    drv = R[(R.site == "mua") & R.condition.str.startswith("mono") & R.driven]
    driven = set(zip(drv.penetration, drv.site_id.astype(str)))
    pens = sorted({p for p, _ in driven})
    with ProcessPoolExecutor(min(n_workers, 12)) as ex:
        done = [r for r in ex.map(work, [(p, driven) for p in pens]) if r is not None]

    # per penetration
    rows, pts = [], []
    for d in done:
        for w, Mw in d["mats"].items():
            a, C = contours(Mw)
            for cs in CSTAR:
                rows.append(dict(penetration=d["penetration"], monkey=d["monkey"], window=w, cstar=cs, a=a,
                                 onset=d["onset"], **summarize(C[cs])))
                pts += [dict(penetration=d["penetration"], window=w, cstar=cs, **r) for r in C[cs].to_dict("records")]
    S = pd.DataFrame(rows)
    Pt = pd.DataFrame(pts)

    # population (mean of window matrices) with penetration bootstrap
    pop_rows, pop_pts, boot_pts, boot_m = [], [], [], []
    rng = np.random.default_rng(0)
    for w in WINDOWS:
        mats = np.stack([d["mats"][w] for d in done])
        a, C, summ = population_contours(mats)
        for cs in CSTAR:
            pop_rows.append(dict(window=w, cstar=cs, a=a, **summ[cs]))
            pop_pts += [dict(window=w, cstar=cs, **r) for r in C[cs].to_dict("records")]
        for b in range(N_BOOT):
            ab, Cb, sb = population_contours(mats[rng.integers(0, len(mats), len(mats))])
            for cs in CSTAR:
                boot_m.append(dict(window=w, cstar=cs, boot=b, m=sb[cs]["m"], suppression=sb[cs]["suppression"]))
                boot_pts += [dict(window=w, cstar=cs, boot=b, c_N=r["c_N"], x=r["x"]) for r in Cb[cs].to_dict("records")]
    Pop, PopPts = pd.DataFrame(pop_rows), pd.DataFrame(pop_pts)
    Boot, BootM = pd.DataFrame(boot_pts), pd.DataFrame(boot_m)
    ci = BootM.groupby(["window", "cstar"]).m.quantile([0.025, 0.975]).unstack()
    Pop = Pop.merge(ci.rename(columns={0.025: "m_lo", 0.975: "m_hi"}), left_on=["window", "cstar"], right_index=True)

    out = out_root() / "binocular"
    S.to_parquet(out / "v3_C_penetration.parquet")
    Pt.to_parquet(out / "v3_C_points.parquet")
    PopPts.merge(Pop, on=["window", "cstar"], suffixes=("", "_summary")).to_parquet(out / "v3_C_population.parquet")
    Boot.to_parquet(out / "v3_C_population_boot.parquet")
    np.savez(out / "v3_C_windows.npz", penetrations=np.array([d["penetration"] for d in done], dtype=str),
             windows=np.array(list(WINDOWS), dtype=str),
             mats=np.stack([np.stack([d["mats"][w] for w in WINDOWS]) for d in done]))

    # registered tests and summaries
    qc = []
    for cs in CSTAR:
        g = S[S.cstar == cs].pivot(index="penetration", columns="window", values="m")
        both = g[["early", "late"]].dropna()
        p_one = stats.wilcoxon(both.late, both.early, alternative="greater").pvalue
        p_two = stats.wilcoxon(both.late, both.early).pvalue
        early_med = both.early.median()
        mk = (both.late - both.early).groupby(both.index.str[7]).median()
        qc.append(dict(test=f"REGISTERED c*={cs:g}: late m > early m (one-sided Wilcoxon); early median in 0.7-1.5",
                       n=len(both), early=early_med, mid=g["mid"].median(), late=both.late.median(),
                       p=p_one, p_two_sided=p_two, frac_late_gt_early=float((both.late > both.early).mean()),
                       B=mk.get("B", np.nan), J=mk.get("J", np.nan),
                       passed=bool((p_one < 0.05) and (0.7 <= early_med <= 1.5))))
        mm = S[S.cstar == cs].groupby("window")[["m", "suppression", "a"]].median()
        n_at_bound = S[S.cstar == cs].groupby("window").m.apply(lambda v: f"{(v <= 0.5).sum()} at 0.5, {(v >= 20).sum()} at 20")
        for w in WINDOWS:
            r = S[(S.cstar == cs) & (S.window == w)]
            qc.append(dict(test=f"per-penetration c*={cs:g} {w}: median m / suppression fraction / a; m at bounds: "
                                f"{n_at_bound[w]}", n=len(r), early=mm.loc[w, "m"], mid=mm.loc[w, "suppression"],
                           late=mm.loc[w, "a"], p=np.nan, p_two_sided=np.nan,
                           frac_late_gt_early=np.nan, B=r.n_below_axis.sum(), J=r.n_censored.sum(), passed=np.nan))
        for w in WINDOWS:
            r = Pop[(Pop.window == w) & (Pop.cstar == cs)].iloc[0]
            qc.append(dict(test=f"population c*={cs:g} {w}: m [95% CI] / suppression", n=len(done), early=r.m,
                           mid=r.m_lo, late=r.m_hi, p=np.nan, p_two_sided=np.nan, frac_late_gt_early=r.suppression,
                           B=np.nan, J=np.nan, passed=np.nan))
    tot = S.assign(points=5).groupby("window")[["points", "n_below_axis", "n_censored", "n_below_axis_fit",
                                                 "n_censored_fit"]].sum()
    for w, r in tot.iterrows():
        qc.append(dict(test=f"contour points {w} (all c*, c_N > 0): total / below axis / censored "
                            f"(fitted range: below {r.n_below_axis_fit}, censored {r.n_censored_fit})",
                       n=int(r.points), early=r.n_below_axis, mid=r.n_censored, late=np.nan, p=np.nan,
                       p_two_sided=np.nan, frac_late_gt_early=np.nan, B=np.nan, J=np.nan, passed=np.nan))
    qc.append(dict(test="penetrations with the late window truncated at 280 ms", n=len(done),
                   early=sum(d["late_truncated"] for d in done), mid=np.nan, late=np.nan, p=np.nan,
                   p_two_sided=np.nan, frac_late_gt_early=np.nan, B=np.nan, J=np.nan, passed=np.nan))
    Q = pd.DataFrame(qc)
    Q.to_csv(REPO / "docs" / "qc" / "v3_C.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print("(registered rows: early/mid/late = median m; B, J = monkey medians of late - early m;\n"
          " per-penetration rows: early = m, mid = suppression, late = a, B = below-axis points, J = censored points;\n"
          " population rows: early = m, mid/late = 95% CI, frac = suppression)")
    print(Q.round(3).to_string(index=False))
    print(figure(PopPts, Boot))


if __name__ == "__main__":
    main()
