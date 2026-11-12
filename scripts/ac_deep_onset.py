"""Check AC (descriptive, post hoc; Blake's comment on draft V6, 2026-11-12): is "deep responds as early as input" a
per-penetration regularity or a few sessions, and could a misplaced boundary or a second fold of cortex explain it?

Usage: .venv\\Scripts\\python scripts\\ac_deep_onset.py
Not registered: run after X3 (deep - input onset -1.35 ms, 14/22 earlier, p = .26 two-sided) was known. Uses check X's
channels and onsets (alignment v3; first-sheet driven MUA; half-max onsets interpolated, dominant eye at 90%) and the
laminar audit's whole-probe monocular cache.
AC1 per penetration: deep - input onset (median channel per group) in all, odd and even presentations;
AC2 the deep group split into its upper (0-0.3 mm below the boundary) and lower (0.3-0.7 mm) part;
AC3 subsets: by monkey; whether the earliest monocular sink lies inside the input band (0-0.5 mm above the boundary,
    "markers agree"); X1-X3 rerun on those subsets from docs/qc/x_groups.csv; the sink's height above the boundary
    against deep - input onset (Spearman over penetrations);
AC4 do the early deep channels normalize like the input group? Per penetration, the shortfall (A - P, 45%, 30-100 ms)
    of deep channels with onsets at or before the input group's median, minus the input group's;
AC5 profiles of shortfall and onset in 0.2 mm bins from the boundary, by subset;
cards: onset against depth per penetration, with the whole-probe monocular onset (beyond the 0.7 mm cap) for second
    early zones.
Writes docs/qc/ac_penetrations.csv, ac_tests.csv, ac_profiles.csv and <out_root>/figures/ac_deep_onset_cards.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import figstyle as fs  # noqa: E402
from v1lb.config import out_root  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
SPLITS = ("all", "odd", "even")
INPUT_MM, CAP_MM, UPPER_MM = 0.5, 0.7, 0.3


def load():
    B = out_root() / "binocular"
    L = pd.read_parquet(B / "x_latency.parquet")
    X = pd.read_parquet(B / "x_channels.parquet")
    L = L[L.grp != ""]
    x = X[(X.split == "all") & (X.c == 0.45) & (X.window == "dip") & (X.grp != "")].copy()
    x["AmP"] = x.added - x.predicted
    al = pd.read_csv(REPO / "metadata" / "laminar_alignment_v3.csv").set_index("penetration")
    return L, x[["penetration", "contact", "AmP"]], al


def half_max_onset(m, t, lo=10, hi=120):
    """As laminar_audit.half_max_onset (stimulus clock, ms)."""
    w = (t >= lo) & (t < hi)
    pk = np.nanmax(m[w])
    j = np.flatnonzero(w & (m >= 0.5 * pk))[0]
    return float(t[j - 1] + (0.5 * pk - m[j - 1]) / (m[j] - m[j - 1])) if w[j - 1] else float(t[j])


def whole_probe(pens, al):
    """Monocular MUA onset (one eye at >= 45%, eyes pooled) on every driven contact, as in the laminar audit."""
    rows = []
    for p in pens:
        M = np.load(out_root() / "laminar_audit" / f"{p}_mono.npz")
        t, c = M["t"], M["contact"].astype(int)
        mua = M["mua_contra"] + M["mua_ipsi"]
        ev = np.maximum(M["evoked_t_contra"], M["evoked_t_ipsi"])
        amp = np.nanmax(mua[(t >= 30) & (t < 120)], axis=0)
        amp_n = amp / np.nanpercentile(amp, 90)
        a = al.loc[p]
        sign = 1 if a.orientation == "ascending" else -1
        for k in range(len(c)):
            rows.append(dict(penetration=p, contact=c[k], rel_anchor_mm=sign * (c[k] - a.anchor_contact) / 10,
                             in_cortex=sign * (c[k] - a.top_contact) >= 0,
                             onset=half_max_onset(mua[:, k], t) if (ev[k] > 6 and amp_n[k] >= 0.2) else np.nan))
    return pd.DataFrame(rows)


def per_penetration(L, x, al):
    T = pd.DataFrame(index=sorted(L.penetration.unique()))
    T.index.name = "penetration"
    T["monkey"] = T.index.str.split("_").str[1]
    for s in SPLITS:
        d = L[L.split == s].copy()
        d["band"] = np.where((d.grp == "deep") & (d.rel_anchor_mm <= UPPER_MM + 1e-9), "deep_upper",
                             np.where(d.grp == "deep", "deep_lower", d.grp))
        g = d.groupby(["penetration", "grp"]).dom_90.median().unstack()
        b = d.groupby(["penetration", "band"]).dom_90.median().unstack()
        T[f"deep_minus_input_{s}"] = g.deep - g.input
        T[f"sup_minus_input_{s}"] = g.superficial - g.input
        T[f"upper_deep_minus_input_{s}"] = b.deep_upper - b.input
        T[f"lower_deep_minus_input_{s}"] = b.deep_lower - b.input
    sign = np.where(al.orientation == "ascending", 1, -1)
    T["sink_height_mm"] = (-(al.sink_contact - al.anchor_contact) * sign / 10).reindex(T.index)
    T["markers_agree"] = (T.sink_height_mm >= 0) & (T.sink_height_mm <= INPUT_MM)
    T["anchor_basis"] = al.anchor_basis.reindex(T.index)
    # AC4: early deep channels' shortfall relative to the input group's
    A = L[L.split == "all"].merge(x, on=["penetration", "contact"], how="left")
    ap = A.groupby(["penetration", "grp"]).AmP.median().unstack()
    T["shortfall_deep_minus_input"] = ap.deep - ap.input
    for p, d in A.groupby("penetration"):
        inp, dp = d[d.grp == "input"], d[d.grp == "deep"].dropna(subset=["dom_90", "AmP"])
        ref = inp.dom_90.median()
        early, late = dp[dp.dom_90 <= ref], dp[dp.dom_90 > ref]
        T.loc[p, "early_deep_shortfall_minus_input"] = early.AmP.median() - inp.AmP.median() if len(early) else np.nan
        T.loc[p, "late_deep_shortfall_minus_input"] = late.AmP.median() - inp.AmP.median() if len(late) else np.nan
        T.loc[p, "n_early_deep"] = len(early)
    return T


def wil(v, label, alternative="two-sided", halves=None):
    v = v.dropna()
    r = dict(test=label, n=len(v), median=v.median(), n_negative=int((v < 0).sum()), n_positive=int((v > 0).sum()),
             p=stats.wilcoxon(v, alternative=alternative).pvalue if len(v) > 5 else np.nan, alternative=alternative)
    if halves is not None:
        r["median_odd"], r["median_even"] = halves[0].dropna().median(), halves[1].dropna().median()
    return r


def tests(T):
    R = []
    on = "deep_minus_input"
    agree_halves = np.sign(T[f"{on}_odd"]) == np.sign(T[f"{on}_even"])
    R.append(wil(T[f"{on}_all"], "AC1 deep - input onset, all 22 (ms)", halves=(T[f"{on}_odd"], T[f"{on}_even"])))
    R.append(wil(T[f"{on}_all"][agree_halves], "AC1 deep - input onset, penetrations whose halves agree in sign"))
    for band in ("upper", "lower"):
        c = f"{band}_deep_minus_input"
        R.append(wil(T[f"{c}_all"], f"AC2 {band} deep - input onset", halves=(T[f"{c}_odd"], T[f"{c}_even"])))
    G = pd.read_csv(QC / "x_groups.csv").set_index(["penetration", "grp"])
    meas = {"X1 shortfall (A - P)": (G.A_minus_P, "greater"), "X2 AI - AI_sum": (G.AI - G.AI_sum, "greater"),
            "X3 onset dominant eye 90%": (G.dom_90, "two-sided")}
    subsets = {"monkey B": T.index[T.monkey == "B"], "monkey J": T.index[T.monkey == "J"],
               "markers agree": T.index[T.markers_agree], "markers disagree": T.index[~T.markers_agree]}
    for name, (col, alt) in meas.items():
        u = col.unstack()
        for sub, idx in subsets.items():
            v = (u.reindex(idx).deep - u.reindex(idx).input)
            R.append(wil(v, f"AC3 {name}, deep - input, {sub}", alt))
    for lab, d in (("all", T), ("without 221215_J_eA (sink below the boundary)", T.drop("221215_J_eA", errors="ignore"))):
        rho, p = stats.spearmanr(d.sink_height_mm, d[f"{on}_all"])
        R.append(dict(test=f"AC3 sink height vs deep - input onset (Spearman), {lab}", n=len(d), median=rho, p=p))
        rho, p = stats.spearmanr(d.sink_height_mm, d.shortfall_deep_minus_input)
        R.append(dict(test=f"AC3 sink height vs deep - input shortfall (Spearman), {lab}", n=len(d), median=rho, p=p))
    rho, p = stats.spearmanr(T[f"{on}_all"], T.shortfall_deep_minus_input)
    R.append(dict(test="AC4 deep - input onset vs deep - input shortfall (Spearman)", n=len(T), median=rho, p=p))
    R.append(wil(T.early_deep_shortfall_minus_input, "AC4 early deep channels' shortfall - input's"))
    R.append(wil(T.late_deep_shortfall_minus_input, "AC4 later deep channels' shortfall - input's"))
    return pd.DataFrame(R)


def profiles(L, x, T):
    A = L[L.split == "all"].merge(x, on=["penetration", "contact"], how="left")
    A["height"] = -A.rel_anchor_mm
    A["bin"] = np.round(np.floor(A.height / 0.2 + 1e-9) * 0.2 + 0.1, 1)
    out = []
    for sub, idx in (("markers agree", T.index[T.markers_agree]), ("markers disagree", T.index[~T.markers_agree]),
                     ("monkey B", T.index[T.monkey == "B"]), ("monkey J", T.index[T.monkey == "J"])):
        per = A[A.penetration.isin(idx)].groupby(["penetration", "bin"])[["AmP", "dom_90"]].median().reset_index()
        per["onset_rel"] = per.dom_90 - per.groupby("penetration").dom_90.transform("median")
        t = per.groupby("bin").agg(n=("AmP", "size"), shortfall=("AmP", "median"), onset_rel=("onset_rel", "median"))
        out.append(t[t.n >= 4].reset_index().assign(subset=sub))
    return pd.concat(out)


def cards(L, W, T, al):
    order = T.sort_values(["markers_agree", "deep_minus_input_all"]).index.tolist()
    COL = {"superficial": fs.MUTED, "input": fs.BLUE, "deep": fs.ORANGE}
    fig, axs = plt.subplots(4, 6, figsize=(180 * fs.MM, 175 * fs.MM), sharey=True, sharex=True)
    for ax, p in zip(axs.flat, order):
        a, r = al.loc[p], T.loc[p]
        sign = 1 if a.orientation == "ascending" else -1
        ax.axhspan(0, INPUT_MM, color=fs.BLUE, alpha=0.08, lw=0)
        ax.axhspan(-CAP_MM, 0, color=fs.ORANGE, alpha=0.08, lw=0)
        ax.axhline(-CAP_MM, color=fs.MUTED, lw=0.5, ls=":")
        ax.axhline(0, color=fs.INK, lw=0.5)
        w = W[(W.penetration == p) & W.in_cortex].sort_values("rel_anchor_mm")
        ax.plot(w.onset, -w.rel_anchor_mm, color=fs.MUTED, lw=0.7, zorder=1)
        d = L[(L.penetration == p) & (L.split == "all")]
        for g, dd in d.groupby("grp"):
            ax.scatter(dd.dom_90, -dd.rel_anchor_mm, s=6, color=COL[g], lw=0, zorder=3)
        for g, y in (("input", (0, INPUT_MM)), ("deep", (-CAP_MM, 0))):
            m = d[d.grp == g].dom_90.median()
            ax.plot([m, m], y, color=COL[g], lw=1.3, zorder=4)
        ax.plot(72, r.sink_height_mm, marker="<", color=fs.INK, ms=3.5, clip_on=False, zorder=5)
        ax.set_xlim(15, 72)
        ax.set_ylim(-1.8, 2.2)
        ax.text(0.03, 0.97, f"{p[:6]} {p[7]} {p[9:]}", transform=ax.transAxes, fontsize=5.5, va="top",
                color=fs.INK if r.markers_agree else fs.MAGENTA)
        ax.text(0.03, 0.86, f"deep−input {r.deep_minus_input_all:+.1f} ms\nhalves {r.deep_minus_input_odd:+.1f} / "
                f"{r.deep_minus_input_even:+.1f}", transform=ax.transAxes, fontsize=4.8, va="top", color=fs.INK2)
        ax.tick_params(labelsize=5)
    for ax in axs.flat[len(order):]:
        ax.axis("off")
    axs.flat[-1].text(
        0, 1, "x: onset, dominant eye at 90% (ms after\nstimulus onset); y: mm above the v3\nboundary (blue band, "
        "input; orange,\ndeep to the 0.7 mm cap, dotted).\nDots: channels (grey superficial, blue\ninput, orange deep); "
        "bars: group\nmedians. Grey line: monocular onset on\nevery driven contact (laminar audit),\nincluding beyond the "
        "cap. Arrow:\nearliest monocular sink. Magenta name:\nthe sink lies outside the input band.\nOrder: markers "
        "disagree, then agree;\nwithin each, by deep − input.", fontsize=5, va="top", transform=axs.flat[-1].transAxes)
    for ax in axs[:, 0]:
        ax.set_ylabel("mm above boundary", fontsize=5.5)
    for ax in axs[-1, :]:
        ax.set_xlabel("onset (ms)", fontsize=5.5)
    fig.tight_layout(h_pad=0.5, w_pad=0.3)
    path = out_root() / "figures" / "ac_deep_onset_cards.png"
    fig.savefig(path, dpi=220)
    return path


def main():
    L, x, al = load()
    T = per_penetration(L, x, al)
    R = tests(T)
    P = profiles(L, x, T)
    W = whole_probe(T.index, al)
    T.to_csv(QC / "ac_penetrations.csv")
    R.to_csv(QC / "ac_tests.csv", index=False)
    P.to_csv(QC / "ac_profiles.csv", index=False)
    print(cards(L, W, T, al))
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 20)
    print(T[["monkey", "sink_height_mm", "markers_agree", "deep_minus_input_all", "deep_minus_input_odd",
             "deep_minus_input_even", "upper_deep_minus_input_all", "lower_deep_minus_input_all",
             "shortfall_deep_minus_input", "early_deep_shortfall_minus_input"]].sort_values(
        "deep_minus_input_all").round(2).to_string())
    print(R.round(4).to_string(index=False))
    for s, g in P.groupby("subset", sort=False):
        print(f"\n{s} (bin = mm above the boundary):")
        print(g.set_index("bin")[["n", "shortfall", "onset_rel"]].round(3).T.to_string())


if __name__ == "__main__":
    main()
