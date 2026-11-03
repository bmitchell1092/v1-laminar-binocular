"""Figure 2 (revised set, docs/09 checks Q, R and U): the full 6 x 6 design over time.

Usage: .venv\\Scripts\\python scripts\\figures\\fig2_full_design.py
Shortfall from summation = both eyes - "if the contrasts add" (the dominant eye's response at c_D + a c_O, the
penetration's median eye weight a, interpolated in log contrast), in units of the dominant-eye peak at 90%, with each
channel aligned to its own response onset (figure_data.py, population.npz); cells with both contrasts >= 11% and
c_D + a c_O <= 90%.
(a) the shortfall over time, averaged over the cells in which the two eyes' contrasts are matched (11, 22.5, 45%),
    the other eye's contrast is higher, or lower; penetration means, 95% bootstrap CI; 5 ms running mean for display;
(b-d) the shortfall against the other eye's contrast in 0-20, 30-100 and 130-200 ms after onset, one line per
    dominant-eye contrast; penetration means, 95% bootstrap CI (v4_full_design.py, v4_Q_cells.csv); the slope on the
    other eye's contrast from the hierarchical model (check Q) in each panel;
(e) the 130-200 ms shortfall with the other eye's contrast higher, per penetration, by monkey, and per held-out
    cross-recording (the grating on the other probe's receptive field; check R); medians.
Inputs: <out_root>/figures/data/population.npz, population_meta.parquet; docs/qc/v4_Q_cells.csv, v4_Q.csv,
v4_R_cells.csv, v4_R.csv.
Writes <out_root>/figures/fig2_full_design.pdf and .png; prints the statistics shown.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.binocular import CONTRASTS, response_at  # noqa: E402
from v1lb.config import out_root  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
QC = REPO / "docs" / "qc"
H_MM = 108
CS = (0.11, 0.225, 0.45, 0.9)
WIN = {"early": ("0–20 ms after onset", (0, 20)), "dip": ("30–100 ms", (30, 100)), "late": ("130–200 ms", (130, 200))}
CATS = {"matched": ("contrasts matched", fs.INK), "higher": ("other eye's contrast higher", fs.COND["other"]),
        "lower": ("other eye's contrast lower", fs.COND["dominant"])}
POS = {"a": (16, 62, 74, 38), "e": (122, 62, 48, 38), "b": (16, 12, 44, 34), "c": (68, 12, 44, 34),
       "d": (120, 12, 44, 34)}


def axes_mm(fig, x, y, w, h):
    W, H = fig.get_size_inches() / fs.MM
    return fig.add_axes((x / W, y / H, w / W, h / H))


def time_courses(ax, pop, meta):
    rel, O = pop["rel"], pop["onset"]
    curves = {k: [] for k in CATS}
    for k, p in enumerate(pop["pens"]):
        a = float(meta.loc[str(p), "a"])
        acc = {c: [] for c in CATS}
        for cd in CS:
            for co in CS:
                if cd + a * co > CONTRASTS[-1] or (cd == co and cd == 0.9):
                    continue
                iD, iO = list(CONTRASTS).index(cd), list(CONTRASTS).index(co)
                s = O[k, iD, iO] - response_at(O[k, :, 0], cd + a * co)
                acc["matched" if cd == co else "higher" if co > cd else "lower"].append(s)
        for c in CATS:
            if acc[c]:
                curves[c].append(np.nanmean(acc[c], axis=0))
    sel = (rel >= -10) & (rel <= 200)
    out = {}
    for c, (lab, col) in CATS.items():
        Y = fs.smooth(np.array(curves[c]))[:, sel]
        out[c] = fs.ci_line(ax, rel[sel], Y, col, alpha=0.14)
    fs.zero_line(ax, color=fs.MUTED)
    for w0, w1 in (v[1] for v in WIN.values()):
        ax.axvspan(w0, w1, color=fs.RULE, alpha=0.35, lw=0, zorder=0)
    ax.set_xlim(-10, 200)
    ax.set_xlabel("time from response onset (ms)")
    ax.set_ylabel("both eyes − if the contrasts add\n(dominant-eye peak = 1)")
    ax.text(102, 0.005, "0: the contrasts add", ha="left", va="bottom", fontsize=5.5, color=fs.INK2)
    ax.set_ylim(-0.085, 0.13)
    fs.line_key(ax, [(lab, col, "-") for lab, col in CATS.values()], x=0.5, y=0.97, dy=0.075, length=0.07)
    return rel[sel], out


def by_other_contrast(axs, Q, T):
    A = Q[Q.split == "all"]
    rows = []
    for ax, (w, (title, _)) in zip(axs, WIN.items()):
        for k, cd in enumerate((0.11, 0.225, 0.45)):
            g = A[np.isclose(A.c_dom, cd) & (A.window == w)].pivot_table(index="unit", columns="c_oth", values="S")
            g = g.reindex(columns=list(CS))
            m, lo, hi = fs.boot_ci(g.to_numpy())
            x = np.array(CS) * (1 + (k - 1) * 0.04)
            ok = np.isfinite(m) & (g.notna().sum().to_numpy() >= 10)
            ax.errorbar(x[ok], m[ok], [m[ok] - lo[ok], hi[ok] - m[ok]], color=fs.CONTRAST[cd], lw=1.0, marker="o",
                        ms=2.6, elinewidth=0.6, capsize=0)
            rows += [dict(window=w, c_dom=cd, c_oth=co, n=int(g[co].notna().sum()), mean=mm)
                     for co, mm, o in zip(CS, m, ok) if o]
        fs.zero_line(ax, color=fs.MUTED)
        fs.log_contrast_axis(ax, CS, label="")
        ax.set_xticklabels(["11", "22.5", "45", "90"])
        ax.set_xlim(0.095, 1.05)
        ax.set_ylim(-0.11, 0.07)
        ax.set_title(title, fontsize=6.5)
        r = T.set_index("test").loc[f"Q x_O slope, {w}"]
        ax.text(0.03, 0.04, f"slope {r.estimate:+.3f} per doubling", transform=ax.transAxes, fontsize=5.5,
                color=fs.INK2)
    axs[0].set_ylabel("both eyes − if the\ncontrasts add")
    axs[1].set_xlabel("other eye's contrast (%)")
    for ax in axs[1:]:
        ax.tick_params(labelleft=False)
    for k, cd in enumerate((0.11, 0.225, 0.45)):
        axs[0].text(0.97, 0.36 - k * 0.09, f"dominant eye {fs.contrast_label(cd)}", transform=axs[0].transAxes,
                    ha="right", va="top", fontsize=5.5, color=fs.CONTRAST[cd])
    return pd.DataFrame(rows)


def sustained(ax, Q, Rc, T, Tr):
    rng = np.random.default_rng(6)
    late = lambda D: D[(D.split == "all") & (D.window == "late") & (D.c_oth > D.c_dom)].groupby("unit").S.mean()  # noqa
    m = late(Q)
    groups = [("monkey B", m[m.index.str[7] == "B"], "o"), ("monkey J", m[m.index.str[7] == "J"], "s"),
              ("held-out\nrecordings (J)", late(Rc), "D")]
    out = []
    for k, (lab, v, mk) in enumerate(groups):
        ax.scatter(np.full(v.size, k) + rng.uniform(-0.12, 0.12, v.size), v, s=7, marker=mk, facecolor="white",
                   edgecolor=fs.INK2, lw=0.6, zorder=2)
        ax.plot([k - 0.25, k + 0.25], [v.median()] * 2, color=fs.INK, lw=1.4, zorder=3)
        out.append(dict(group=lab.replace("\n", " "), n=len(v), median=float(v.median()), n_below=int((v < 0).sum())))
    fs.zero_line(ax, color=fs.MUTED)
    ax.set_xticks(range(3), [g[0] for g in groups])
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("both eyes − if the contrasts add,\nother eye higher, 130–200 ms")
    q = T.set_index("test").loc["Q late shortfall, other eye stronger (c_O > c_D)"]
    r = Tr.set_index("test").loc["R late shortfall, other eye stronger (c_O > c_D)"]
    ax.text(0.0, 1.03, f"model estimate, all penetrations: {q.estimate:+.3f}, {fs.fmt_p(q.p_one_sided)}\n"
            f"held-out recordings: {r.estimate:+.3f}, {fs.fmt_p(r.p_one_sided)}", transform=ax.transAxes,
            va="bottom", fontsize=5.5, color=fs.INK, linespacing=1.1)
    return pd.DataFrame(out)


def main():
    fs.apply()
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    Q, T = pd.read_csv(QC / "v4_Q_cells.csv"), pd.read_csv(QC / "v4_Q.csv")
    Rc, Tr = pd.read_csv(QC / "v4_R_cells.csv"), pd.read_csv(QC / "v4_R.csv")
    fig = fs.figure("double", H_MM)
    fig.set_layout_engine("none")
    ax = {k: axes_mm(fig, *v) for k, v in POS.items()}
    t, curves = time_courses(ax["a"], pop, meta)
    B = by_other_contrast([ax["b"], ax["c"], ax["d"]], Q, T)
    E = sustained(ax["e"], Q, Rc, T, Tr)
    for k in "aebcd":
        fs.panel_label(ax[k], k, x=-0.25 if k in "ab" else -0.3 if k == "e" else -0.04)
    print(fs.save(fig, "fig2_full_design"))
    pd.set_option("display.width", 220)
    print("(a) window means of the time courses (penetration means):")
    for c, m in curves.items():
        print(f"   {c}: " + "; ".join(f"{w0}-{w1} ms {np.nanmean(m[(t >= w0) & (t < w1)]):+.3f}"
                                     for w0, w1 in (v[1] for v in WIN.values())))
    print("(b-d) cell means:")
    print(B.round(3).pivot_table(index=["window", "c_dom"], columns="c_oth", values="mean").to_string())
    print("(e)")
    print(E.round(4).to_string(index=False))
    print(T.round(4).to_string(index=False))
    print(Tr.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
