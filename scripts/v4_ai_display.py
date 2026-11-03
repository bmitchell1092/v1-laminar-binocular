"""V4 check V: the additivity index as the displayed measure (docs/09, registered 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\v4_ai_display.py
AI = mean(B) / mean(D + O) and AI_sum = mean(R_D(c_D + a c_O)) / mean(D + O) in windows after each channel's response
onset (0-20, 30-100, 130-200 ms), where D + O > 0.12 of the dominant-eye peak at 90%.
* per penetration and cell (both contrasts >= 11%, c_D + a c_O <= 90%) from population.npz (channel means, the
  penetration's median eye weight);
* per penetration and depth group (above / at, within 2 contacts / below the earliest sink) at 45% from
  spacetime.npz (channel means within the group; each channel's own eye weight for the summation prediction).
V1 matched 45%, 30-100 ms, AI < AI_sum; V2 (AI - AI_sum) below - at > 0 (45%, 30-100 ms); V3 130-200 ms, cells with
the other eye stronger (mean over cells per penetration), AI < AI_sum; V4 (reported) 0-20 ms, AI - AI_sum.
One-sided Wilcoxon over penetrations. Writes docs/qc/v4_V.csv, v4_V_cells.csv and v4_V_groups.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from v1lb.binocular import CONTRASTS, response_at
from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
WINDOWS = {"early": (0, 20), "dip": (30, 100), "late": (130, 200)}
CS = (0.11, 0.225, 0.45, 0.9)
FLOOR = 0.12


def cells():
    pop = dict(np.load(out_root() / "figures" / "data" / "population.npz"))
    meta = pd.read_parquet(out_root() / "figures" / "data" / "population_meta.parquet").set_index("penetration")
    O, rel = pop["onset"], pop["rel"]
    rows = []
    for k, p in enumerate(pop["pens"]):
        a = float(meta.loc[str(p), "a"])
        for cd in CS:
            for co in CS:
                if cd + a * co > CONTRASTS[-1]:
                    continue
                iD, iO = list(CONTRASTS).index(cd), list(CONTRASTS).index(co)
                B, S = O[k, iD, iO], O[k, iD, 0] + O[k, 0, iO]
                P = response_at(O[k, :, 0], cd + a * co)
                for w, (w0, w1) in WINDOWS.items():
                    sel = (rel >= w0) & (rel < w1)
                    s = float(np.nanmean(S[sel]))
                    rows.append(dict(penetration=str(p), monkey=str(p)[7], c_dom=cd, c_oth=co, window=w,
                                     AI=float(np.nanmean(B[sel])) / s if s > FLOOR else np.nan,
                                     AI_sum=float(np.nanmean(P[sel])) / s if s > FLOOR else np.nan))
    return pd.DataFrame(rows)


def groups(c=0.45):
    Z = np.load(out_root() / "figures" / "data" / "spacetime.npz")
    meta = pd.read_parquet(out_root() / "figures" / "data" / "spacetime_meta.parquet")
    X, rel = Z["onset"], Z["rel"]
    i, j = list(CONTRASTS).index(c), CS.index(c)
    D, Oe, B = X[:, i], X[:, 6 + j], X[:, 10 + j]
    P = np.stack([response_at(X[q, :6], c * (1 + meta.a.iloc[q])) for q in range(len(meta))])
    rs = meta.rel_sink.to_numpy()
    grp = np.where(np.abs(rs) <= 2, "at", np.where(rs < -2, "above", "below"))
    rows = []
    for (p, g), idx in pd.Series(range(len(meta))).groupby([meta.penetration, grp]):
        q = idx.to_numpy()
        for w, (w0, w1) in WINDOWS.items():
            sel = (rel >= w0) & (rel < w1)
            s = float(np.nanmean((D[q] + Oe[q])[:, sel]))
            rows.append(dict(penetration=p, monkey=p[7], grp=g, window=w, n_channels=len(q),
                             AI=float(np.nanmean(B[q][:, sel])) / s if s > FLOOR else np.nan,
                             AI_sum=float(np.nanmean(P[q][:, sel])) / s if s > FLOOR else np.nan))
    return pd.DataFrame(rows)


def wil(v, alt, name, extra=None):
    v = pd.Series(v).dropna()
    p = stats.wilcoxon(v, alternative=alt).pvalue
    want = -1 if alt == "less" else 1
    return dict(test=name, n=len(v), median=float(v.median()), n_in_direction=int((np.sign(v) == want).sum()), p=p,
                holds=bool(p < .05) if alt != "two-sided" else np.nan, **(extra or {}))


def main():
    C = cells()
    C.to_csv(QC / "v4_V_cells.csv", index=False)
    G = groups()
    G.to_csv(QC / "v4_V_groups.csv", index=False)
    C["gap"], G["gap"] = C.AI - C.AI_sum, G.AI - G.AI_sum
    rows = []
    m = C[np.isclose(C.c_dom, 0.45) & np.isclose(C.c_oth, 0.45)].set_index(["penetration", "window"])
    rows.append(wil(m.xs("dip", level="window").gap, "less", "V1 matched 45%, 30-100 ms: AI - AI_sum < 0",
                    dict(AI=float(m.xs("dip", level="window").AI.median()),
                         AI_sum=float(m.xs("dip", level="window").AI_sum.median()))))
    piv = G[G.window == "dip"].pivot(index="penetration", columns="grp", values="gap")
    rows.append(wil(piv["below"] - piv["at"], "greater", "V2 (AI - AI_sum) below - at the sink > 0, 45%, 30-100 ms"))
    hi = C[(C.window == "late") & (C.c_oth > C.c_dom)].groupby("penetration")[["gap", "AI", "AI_sum"]].mean()
    rows.append(wil(hi.gap, "less", "V3 other eye stronger, 130-200 ms: AI - AI_sum < 0",
                    dict(AI=float(hi.AI.median()), AI_sum=float(hi.AI_sum.median()))))
    e = m.xs("early", level="window")
    rows.append(wil(e.gap, "two-sided", "(reported) V4 matched 45%, 0-20 ms: AI - AI_sum",
                    dict(AI=float(e.AI.median()), AI_sum=float(e.AI_sum.median()))))
    for g in ("above", "at", "below"):
        x = G[(G.window == "dip") & (G.grp == g)]
        rows.append(wil(x.gap, "less", f"(reported) {g} sink, 45%, 30-100 ms: AI - AI_sum < 0",
                        dict(AI=float(x.AI.median()), AI_sum=float(x.AI_sum.median()))))
    T = pd.DataFrame(rows)
    T.to_csv(QC / "v4_V.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
