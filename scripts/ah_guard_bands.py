"""Check AH: guard bands at the coarse compartment boundaries (docs/09, registered 2026-11-13).

Usage: .venv\\Scripts\\python scripts\\ah_guard_bands.py
On laminar alignment v5's coarse compartments (scaled heights; ad_compartments_v4.relabel), channels within g mm of a
compartment boundary (0 and 0.59 mm) are dropped, g = 0.1 and 0.2 mm. AD1's tests on the rest: infragranular -
granular > 0 in the shortfall (A - P, 45%, 30-100 ms; odd and even halves) and in AI - AI_sum; supragranular -
granular two-sided. Sets: tiers 1 + 2 (main), tier 1, all anchored. Writes docs/qc/ah_tests.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from ad_compartments_v4 import alignment, diff, load_channels, relabel
from x_depth_v3 import SPLITS, ai_groups, shortfall_groups, wil

REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
BOUNDS = (0.0, 0.59)
SETS = {"tiers 1+2": (1, 2), "tier 1": (1,), "all anchored": (0, 1, 2)}


def main():
    V = alignment()
    V = V[V.s.notna()]
    C0, _ = load_channels(V)
    C = relabel(C0, V)
    R = []
    for g in (0.0, 0.1, 0.2):
        keep = np.ones(len(C), bool)
        for b in BOUNDS:
            keep &= np.abs(C.h_scaled.to_numpy() - b) >= g - 1e-9
        Cg = C[keep]
        for name, tiers in SETS.items():
            X = Cg[Cg.tier.isin(tiers)]
            tag = f"[{name}, guard {g:.1f} mm]"
            S = {s: shortfall_groups(X, grp="coarse", split=s) for s in SPLITS}
            d = {s: diff(S[s], "infragranular", "granular") for s in SPLITS}
            R.append(wil(d["all"], "greater", f"AH shortfall infragranular - granular > 0 {tag}", (d["odd"], d["even"])))
            u = diff(S["all"], "supragranular", "granular")
            R.append(wil(u, "two-sided", f"AH (reported) shortfall supragranular - granular {tag}"))
            G = ai_groups(X, grp="coarse")
            gap = (G.AI - G.AI_sum).unstack()
            R.append(wil(diff(gap, "infragranular", "granular"), "greater",
                         f"AH (AI - AI_sum) infragranular - granular > 0 {tag}"))
            R[-1]["n_channels"] = int(X[(X.split == "all") & np.isclose(X.c, 0.45) & (X.window == "dip")
                                        & (X.coarse != "")].shape[0])
    R = pd.DataFrame(R)
    R.to_csv(QC / "ah_tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(R[["test", "n", "median", "n_in_direction", "p", "holds", "n_channels"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
