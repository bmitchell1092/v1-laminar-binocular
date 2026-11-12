"""Check AD: laminar results on alignment v4 tiers, coarse and fine compartments (docs/09, registered 2026-11-12).

Usage: .venv\\Scripts\\python scripts\\ad_compartments_v4.py
Relabels check X's channels (x_channels.parquet, x_latency.parquet) with metadata/laminar_alignment_v4.csv: height
above the v3 anchor h = -rel_anchor_mm, scaled by the penetration's squish factor s (tier 1; 1 otherwise).
Coarse: supragranular h > 0.59 s, granular 0 < h <= 0.59 s, infragranular -0.49 s <= h < 0.
Fine: L2/3 > 0.59 s; L4A/B (0.28 s, 0.59 s]; L4C (0, 0.28 s]; L5 [-0.22 s, 0); L6 [-0.49 s, -0.22 s).
Tests AD1-AD5 on tier 1 (primary), tiers 1 + 2 and all 22. Writes docs/qc/ad_tests.csv and docs/qc/ad_groups.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from v1lb.config import metadata_dir, out_root
from x_depth_v3 import SPLITS, ai_groups, shortfall_groups, wil

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
COARSE = ("supragranular", "granular", "infragranular")
FINE = ("L2/3", "L4A/B", "L4C", "L5", "L6")


def label(h, s):
    coarse = np.select([h > 0.59 * s, h > 0, h >= -0.49 * s - 1e-9], ["supragranular", "granular", "infragranular"], "")
    fine = np.select([h > 0.59 * s, h > 0.28 * s, h > 0, h >= -0.22 * s - 1e-9, h >= -0.49 * s - 1e-9],
                     ["L2/3", "L4A/B", "L4C", "L5", "L6"], "")
    return coarse, fine


def relabel(df, V):
    df = df[df.grp != ""].merge(V[["penetration", "tier", "s"]], on="penetration")
    h = -df.rel_anchor_mm.to_numpy()
    df["coarse"], df["fine"] = label(h, df.s.to_numpy())
    df["h_scaled"] = h / df.s
    return df


def diff(tab, a, b):
    return tab[a] - tab[b] if a in tab and b in tab else pd.Series(dtype=float)


def tests(C, L, name):
    R = []
    tag = lambda t: f"{t} [{name}]"                                       # noqa: E731
    # AD1 coarse shortfall and index
    S = {s: shortfall_groups(C, grp="coarse", split=s) for s in SPLITS}
    d = {s: diff(S[s], "infragranular", "granular") for s in SPLITS}
    R.append(wil(d["all"], "greater", tag("AD1 shortfall infragranular - granular > 0"), (d["odd"], d["even"])))
    u = {s: diff(S[s], "supragranular", "granular") for s in SPLITS}
    R.append(wil(u["all"], "two-sided", tag("AD1 (reported) shortfall supragranular - granular"), (u["odd"], u["even"])))
    for g in COARSE:
        R.append(wil(S["all"][g], "less", tag(f"AD1 (reported) shortfall < 0, {g}"), (S["odd"][g], S["even"][g])))
    G = ai_groups(C, grp="coarse")
    gap = (G.AI - G.AI_sum).unstack()
    R.append(wil(diff(gap, "infragranular", "granular"), "greater", tag("AD1 (AI - AI_sum) infragranular - granular > 0")))
    R.append(wil(diff(gap, "supragranular", "granular"), "two-sided", tag("AD1 (reported) (AI - AI_sum) supra - granular")))
    # AD2 fine shortfall relative to L4C
    F = {s: shortfall_groups(C, grp="fine", split=s) for s in SPLITS}
    for g, alt in (("L2/3", "two-sided"), ("L4A/B", "two-sided"), ("L5", "greater"), ("L6", "greater")):
        v = {s: diff(F[s], g, "L4C") for s in SPLITS}
        R.append(wil(v["all"], alt, tag(f"AD2 shortfall {g} - L4C"), (v["odd"], v["even"])))
    # AD3 onsets
    for grp, pairs in (("coarse", (("supragranular", "granular", "greater"), ("infragranular", "granular", "two-sided"))),
                       ("fine", (("L2/3", "L4C", "greater"), ("L4A/B", "L4C", "greater"), ("L5", "L4C", "greater"),
                                 ("L6", "L4C", "two-sided")))):
        on = {s: L[(L.split == s) & (L[grp] != "")].groupby(["penetration", grp]).dom_90.median().unstack()
              for s in SPLITS}
        for a, b, alt in pairs:
            v = {s: diff(on[s], a, b) for s in SPLITS}
            R.append(wil(v["all"], alt, tag(f"AD3 onset dominant eye 90%, {a} - {b} (ms)"), (v["odd"], v["even"])))
    # AD4 binocular advance per coarse compartment
    adv = {s: L[(L.split == s) & (L.coarse != "")].groupby(["penetration", "coarse"]).adv_45.median().unstack()
           for s in SPLITS}
    for g in COARSE:
        R.append(wil(adv["all"][g], "less", tag(f"AD4 binocular advance at 45% < 0, {g} (ms)"),
                     (adv["odd"][g], adv["even"][g])))
    # AD5 continuous scaled depth
    X = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip") & (C.coarse != "")]
    X = X.assign(S=X.added - X.predicted)
    fit = smf.mixedlm("S ~ h_scaled", X, groups=X.penetration).fit(reml=True)
    est, p = float(fit.params["h_scaled"]), float(fit.pvalues["h_scaled"])
    R.append(dict(test=tag("AD5 shortfall ~ scaled height above the anchor (per mm), slope < 0"),
                  n=X.penetration.nunique(), n_channels=len(X), median=est, p=p / 2 if est < 0 else 1 - p / 2,
                  holds=bool(est < 0 and p / 2 < .05)))
    return R


def groups(C, L):
    out = []
    for grp in ("coarse", "fine"):
        sf = shortfall_groups(C, grp=grp).stack().rename("A_minus_P")
        ai = ai_groups(C, grp=grp)
        on = L[(L.split == "all") & (L[grp] != "")].groupby(["penetration", grp])[["dom_90", "adv_45"]].median()
        n = C[(C.split == "all") & np.isclose(C.c, 0.45) & (C.window == "dip") & (C[grp] != "")].groupby(
            ["penetration", grp]).size().rename("n_channels")
        g = pd.concat([n, sf, ai, on], axis=1).reset_index()
        out.append(g.rename(columns={g.columns[1]: "compartment"}).assign(scheme=grp))
    return pd.concat(out)


def main():
    V = pd.read_csv(metadata_dir() / "laminar_alignment_v4.csv")
    V = V[V.s.notna()]
    B = out_root() / "binocular"
    C = relabel(pd.read_parquet(B / "x_channels.parquet"), V)
    L = relabel(pd.read_parquet(B / "x_latency.parquet"), V)
    sets = {"tier 1": (1,), "tiers 1+2": (1, 2), "all 22": (0, 1, 2)}
    R = []
    for name, tiers in sets.items():
        R += tests(C[C.tier.isin(tiers)], L[L.tier.isin(tiers)], name)
    R = pd.DataFrame(R)
    R.to_csv(QC / "ad_tests.csv", index=False)
    G = groups(C, L).merge(V[["penetration", "tier", "s"]], on="penetration")
    G.to_csv(QC / "ad_groups.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(R.drop(columns=[c for c in ("n_channels",) if c in R]).round(4).to_string(index=False))
    for name, tiers in sets.items():
        g = G[G.tier.isin(tiers)]
        print(f"\n{name}: medians over penetrations")
        print(g.groupby(["scheme", "compartment"])[["n_channels", "A_minus_P", "AI", "AI_sum", "dom_90",
                                                     "adv_45"]].median().round(3).to_string())
        print(g.groupby(["scheme", "compartment"]).penetration.nunique().rename("n_pen").to_string())


if __name__ == "__main__":
    main()
