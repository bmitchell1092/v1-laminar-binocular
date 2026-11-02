"""Supplementary Table S1: the main statistics per monkey (docs/09, figure set for V3).

Usage: .venv\\Scripts\\python scripts\\figures\\tableS1_monkeys.py
Every row is recomputed from saved outputs with the statistic of its figure, separately for monkey B, monkey J and
both. Value = median over penetrations; sign = penetrations in the direction of the paper's claim / penetrations;
p = Wilcoxon signed-rank test over penetrations (one-sided in the claimed direction; two-sided for the plain early
and late values, which carry no test), or a one-sided sign test for the slopes, as in Fig 3b and Supplementary Fig S2b.
* linearity index early (0-20 ms after onset), late (50-200 ms) and early - late, both eyes at 11, 22.5 and 45%
  (Fig 2d: per penetration, mean of the windows' channel medians; k_windows.parquet);
* slope of the transition-window linearity index (window starts 20-35 ms) on log contrast (Fig 3b, S2b;
  normalization_timing.parquet, all presentations);
* iso-response contour exponent, late - early, at 11, 22.5 and 45% (Fig S3; v3_C_penetration.parquet);
* delayed vs static normalization held out on the interaction loss, % less error (Fig 4d; v3_B2.parquet);
* pool weight - drive weight of the other eye, conductance and fixed-delay fits (Fig 4e; fig5_cond.parquet,
  fig5_model.parquet);
* conductance vs fixed-delay held out, % less error (Fig S5f; fig5_cond.parquet). Held-out comparisons leave out
  211008, whose even half has no usable condition (docs/09).
Writes docs/qc/tableS1_monkeys.csv; prints the table.
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from v1lb.config import out_root

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
DIOPTIC = [0.11, 0.225, 0.45]
EARLY, LATE = [0, 10], [50, 75, 100, 150]             # window starts, as Fig 2d
NO_EVEN_HALF = "211008_B_eA"


def monkey_of(index) -> pd.Series:
    return pd.Series([p.split("_")[1] for p in index], index=index)


def rows_for(name: str, v: pd.Series, direction: str | None, test: str = "wilcoxon", pct: bool = False) -> dict:
    """One table row: median, sign count and p per monkey and overall. direction: 'pos', 'neg' or None (no test)."""
    v = v.dropna()
    mk = monkey_of(v.index)
    row = {"statistic": name}
    for grp, x in (("B", v[mk == "B"]), ("J", v[mk == "J"]), ("all", v)):
        val = 100 * (1 - np.exp(-x)) if pct else x
        row[f"{grp}_n"] = len(x)
        row[f"{grp}_value"] = float(np.median(val)) if len(x) else np.nan
        if direction is None or len(x) == 0:
            row[f"{grp}_sign"], row[f"{grp}_p"] = "", np.nan
            continue
        k = int((x > 0).sum()) if direction == "pos" else int((x < 0).sum())
        row[f"{grp}_sign"] = f"{k}/{len(x)}"
        alt = "greater" if direction == "pos" else "less"
        if test == "sign":
            row[f"{grp}_p"] = stats.binomtest(k, len(x), 0.5, alternative="greater").pvalue
        else:
            row[f"{grp}_p"] = stats.wilcoxon(x, alternative=alt).pvalue if len(x) > 1 else np.nan
    return row


def early_late_rows() -> list[dict]:
    W = pd.read_parquet(out_root() / "binocular" / "k_windows.parquet")
    out = []
    for c in DIOPTIC:
        per = W[W.c == c].groupby(["penetration", "w0"]).L.median().unstack()
        e, l_ = per[EARLY].mean(axis=1), per[LATE].mean(axis=1)
        ok = e.notna() & l_.notna()
        e, l_ = e[ok], l_[ok]
        lab = f"{100 * c:g}%"
        out.append(rows_for(f"linearity index, early (0-20 ms), {lab}", e, None))
        out.append(rows_for(f"linearity index, late (50-200 ms), {lab}", l_, None))
        out.append(rows_for(f"linearity index, early - late, {lab}", e - l_, "pos"))
    return out


def slope_row() -> dict:
    R = pd.read_parquet(out_root() / "binocular" / "normalization_timing.parquet")
    per = R[R.split == "all"].groupby(["name", "c", "s"]).L.median().unstack("s")
    tr = per.loc[:, (per.columns >= 20) & (per.columns <= 35)].median(axis=1).unstack("c").reindex(columns=DIOPTIC)
    sl = tr.dropna().apply(lambda r: np.polyfit(np.log(DIOPTIC), r.values, 1)[0], axis=1)
    return rows_for("slope of transition-window (20-50 ms) linearity index on log contrast", sl, "neg", test="sign")


def contour_rows() -> list[dict]:
    S = pd.read_parquet(out_root() / "binocular" / "v3_C_penetration.parquet")
    out = []
    for c in DIOPTIC:
        m = S[S.cstar == c].pivot(index="penetration", columns="window", values="m")
        out.append(rows_for(f"iso-response contour exponent, late - early, {100 * c:g}%", m["late"] - m["early"],
                            "pos"))
    return out


def model_rows() -> list[dict]:
    b = out_root() / "binocular"
    B2 = pd.read_parquet(b / "v3_B2.parquet")
    C = pd.read_parquet(b / "fig5_cond.parquet").set_index("penetration")
    M = pd.read_parquet(b / "fig5_model.parquet").set_index("penetration")
    held = C.drop(index=NO_EVEN_HALF, errors="ignore")
    return [
        rows_for("delayed vs static, held-out interaction error, % less", np.log(B2.cv_static / B2.cv_full), "pos",
                 pct=True),
        rows_for("pool weight - drive weight of the other eye, conductance fit", C.cond_p_w - C.cond_p_a, "pos"),
        rows_for("pool weight - drive weight of the other eye, fixed-delay fit", M.p_w - M.p_a, "pos"),
        rows_for("conductance vs fixed-delay, held-out error, % less", np.log(held.cv_full / held.cv_cond), "pos",
                 pct=True),
    ]


def main():
    rows = early_late_rows() + [slope_row()] + contour_rows() + model_rows()
    T = pd.DataFrame(rows)
    cols = ["statistic"] + [f"{g}_{k}" for g in ("B", "J", "all") for k in ("n", "value", "sign", "p")]
    T = T[cols]
    path = REPO / "docs" / "qc" / "tableS1_monkeys.csv"
    T.to_csv(path, index=False, float_format="%.4g")
    print(path)
    with pd.option_context("display.width", 250, "display.max_colwidth", 80):
        print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
