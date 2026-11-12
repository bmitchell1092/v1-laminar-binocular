"""Check AD: laminar results on alignment v4 tiers, coarse and fine compartments (docs/09, registered 2026-11-12).

Usage: .venv\\Scripts\\python scripts\\ad_compartments_v4.py
Relabels check X's channels (x_channels.parquet, x_latency.parquet) with metadata/laminar_alignment_{ALIGN}.csv (v5
since 2026-11-13: anchor at the bottom of the initial sink; v4 outputs kept as ad_*.csv, v5 as ad5_*.csv): height
above the v4 anchor from contact numbers and the probe's spacing, h = -sign (contact - anchor) * spacing, scaled by the
penetration's squish factor s (tier 1; 1 otherwise); first-sheet contacts only (at or below the top of cortex).
Probes with 50 um spacing (docs/05 C18) reach contacts X never measured, so their channels are recomputed with X's
per-channel code (x_depth_v3.penetration, spacing-aware) and cached in <out_root>/binocular/ad_spacing_*.parquet.
Coarse: supragranular h > 0.59 s, granular 0 < h <= 0.59 s, infragranular -0.49 s <= h < 0.
Fine: L2/3 > 0.59 s; L4A/B (0.28 s, 0.59 s]; L4C (0, 0.28 s]; L5 [-0.22 s, 0); L6 [-0.49 s, -0.22 s).
Tests AD1-AD5 on tier 1 (primary), tiers 1 + 2 and all 22. Writes docs/qc/ad_tests.csv and docs/qc/ad_groups.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from v1lb.config import metadata_dir, out_root
from x_depth_v3 import SPLITS, ai_groups, penetration, shortfall_groups, wil

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
COARSE = ("supragranular", "granular", "infragranular")
ALIGN = "v5"          # laminar alignment for AD, AE and the laminar figures (2026-11-13; v4 results kept as ad_*/ae_*)
AD_PFX = {"v4": "ad_", "v5": "ad5_"}
AE_PFX = {"v4": "ae_", "v5": "ae5_"}


def alignment(version=ALIGN):
    """The laminar alignment table (metadata/laminar_alignment_{version}.csv)."""
    return pd.read_csv(metadata_dir() / f"laminar_alignment_{version}.csv")
FINE = ("L2/3", "L4A/B", "L4C", "L5", "L6")


def label(h, s):
    coarse = np.select([h > 0.59 * s, h > 0, h >= -0.49 * s - 1e-9], ["supragranular", "granular", "infragranular"], "")
    fine = np.select([h > 0.59 * s, h > 0.28 * s, h > 0, h >= -0.22 * s - 1e-9, h >= -0.49 * s - 1e-9],
                     ["L2/3", "L4A/B", "L4C", "L5", "L6"], "")
    return coarse, fine


def relabel(df, V):
    df = df.merge(V[["penetration", "tier", "s", "orientation", "anchor_contact", "top_contact", "spacing_mm"]],
                  on="penetration")
    sign = np.where(df.orientation == "ascending", 1, -1)
    df = df[sign * (df.contact - df.top_contact) >= 0].copy()
    sign = np.where(df.orientation == "ascending", 1, -1)
    h = -sign * (df.contact - df.anchor_contact).to_numpy() * df.spacing_mm.to_numpy()
    df["coarse"], df["fine"] = label(h, df.s.to_numpy())
    df["h_scaled"] = h / df.s
    return df.drop(columns=["orientation", "anchor_contact", "top_contact"])


def spacing_cache(V):
    """Per-channel rows, onsets and space-time arrays of the non-100-um probes, recomputed with X's per-channel code
    (x_depth_v3.penetration, spacing-aware) and cached in <out_root>/binocular/ad_spacing_*."""
    B = out_root() / "binocular"
    odd = sorted(V[V.anchor_contact.notna() & ~np.isclose(V.spacing_mm, 0.1)].penetration)
    f = {k: B / f"ad_spacing_{k}" for k in ("channels.parquet", "latency.parquet", "meta.parquet", "spacetime.npz")}
    if not odd:
        return odd, None
    if not (all(x.exists() for x in f.values()) and sorted(set(pd.read_parquet(f["channels.parquet"]).penetration)) == odd):
        A3 = pd.read_csv(metadata_dir() / "laminar_alignment_v3.csv").set_index("penetration")
        done = [penetration((p, A3.loc[p].to_dict())) for p in odd]
        pd.DataFrame([r for rows, _, _ in done for r in rows]).to_parquet(f["channels.parquet"])
        pd.DataFrame([r for _, lat, _ in done for r in lat]).assign(
            adv_45=lambda x: x.bin_45 - x.dom_45).to_parquet(f["latency.parquet"])
        st = [x for _, _, x in done if x]
        pd.DataFrame([m for x in st for m in x["meta"]]).to_parquet(f["meta.parquet"])
        np.savez_compressed(f["spacetime.npz"], onset=np.concatenate([np.stack(x["onset"]) for x in st]).astype(np.float32),
                            stim=np.concatenate([np.stack(x["stim"]) for x in st]).astype(np.float32))
    return odd, f


def load_channels(V):
    """X's channel rows and onsets, with those of non-100-um probes replaced by the recomputed ones."""
    B = out_root() / "binocular"
    C, L = pd.read_parquet(B / "x_channels.parquet"), pd.read_parquet(B / "x_latency.parquet")
    odd, f = spacing_cache(V)
    if odd:
        C = pd.concat([C[~C.penetration.isin(odd)], pd.read_parquet(f["channels.parquet"])], ignore_index=True)
        L = pd.concat([L[~L.penetration.isin(odd)], pd.read_parquet(f["latency.parquet"])], ignore_index=True)
    return C, L


def load_spacetime(V):
    """X's space-time arrays (onset and stimulus clocks) and their channel table, with the non-100-um probes
    replaced. Returns onset [ch, 14, REL], stim [ch, 14, T_STIM], meta, rel, t."""
    d = out_root() / "figures" / "data"
    Z = np.load(d / "spacetime_v3.npz")
    M = pd.read_parquet(d / "spacetime_v3_meta.parquet")
    keep = ~M.penetration.isin(sorted(V[~np.isclose(V.spacing_mm.fillna(0.1), 0.1)].penetration)).to_numpy()
    on, stim, M = Z["onset"][keep], Z["stim"][keep], M[keep]
    odd, f = spacing_cache(V)
    if odd:
        S = np.load(f["spacetime.npz"])
        on, stim = np.concatenate([on, S["onset"]]), np.concatenate([stim, S["stim"]])
        M = pd.concat([M, pd.read_parquet(f["meta.parquet"])], ignore_index=True)
    return on, stim, M.reset_index(drop=True), Z["rel"], Z["t"]

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
    V = alignment()
    V = V[V.s.notna()]
    C0, L0 = load_channels(V)
    C, L = relabel(C0, V), relabel(L0, V)
    sets = {"tier 1": (1,), "tiers 1+2": (1, 2), "all anchored": (0, 1, 2)}
    R = []
    for name, tiers in sets.items():
        R += tests(C[C.tier.isin(tiers)], L[L.tier.isin(tiers)], name)
    R = pd.DataFrame(R)
    R.to_csv(QC / f"{AD_PFX[ALIGN]}tests.csv", index=False)
    G = groups(C, L).merge(V[["penetration", "tier", "s"]], on="penetration")
    G.to_csv(QC / f"{AD_PFX[ALIGN]}groups.csv", index=False)
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
