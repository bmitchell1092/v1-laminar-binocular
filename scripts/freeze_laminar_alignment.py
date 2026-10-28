"""Freeze the laminar alignment table (metadata/laminar_alignment_v1.csv) from the laminar reports.

Decision (Blake, 2026-10-28): keep tiers A and B, set tier C aside. Claude reviewed the A/B
calls and applied the exceptions below, each with its reason, so the dataset stays solid.
Columns:
  include, tier, orientation ('ascending' = contact number increases with depth),
  top_of_cortex, l45, spiking_bottom (advisory: deepest spiking contact before a silent run of
  >= 3; underestimates on small-spike probes, so it is not a cortical boundary), deep_limit
  (hard limit, only where a second fold is confirmed), basis (rule / override), reason.
Depth of contact c relative to L4/5, in contacts (100 um), positive = superficial:
  ascending: l45 - c;  descending: c - l45.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from v1lb.config import out_root

REPO = Path(__file__).resolve().parents[1]

# L4/5 overrides: the rule's step was contradicted by both cross-checks and by anatomy.
OVERRIDE_L45 = {
    "211208_B_eA": (13.5, "spiking at 6-18, silent 19-28, spiking again 29-32: cortex, white matter, "
                          "second fold. First-column phase step 13.5 matches spectral crossing 13.9 and "
                          "earliest sink 13; the rule's 28.5 is the entry into the second fold."),
}
# Confirmed second folds: contacts deeper than this belong to another column.
DEEP_LIMIT = {"211208_B_eA": 18}
# A/B penetrations set aside, with reasons.
SET_ASIDE = {
    "221128_J_eA": "tier A, but SPI and spectral both call it descending against the legacy ascending sort; "
                   "the sink sits only ~5 contacts below the first spiking contact in either reading.",
    "211103_B_eA": "two candidate L4/5 steps (13.5, 29.5): cross-checks (15.3, 12) favour 13.5, which is "
                   "only 2.5 contacts below the top of cortex; 29.5 is 18.5 contacts below it.",
    "221209_J_eB": "spectral crossing (14.7) and earliest sink (14) agree with each other but with "
                   "neither phase step (18.5, 21.5).",
    "221214_J_eB": "L4/5 (27.5) only 3.5 contacts below the top of cortex (24): anatomically implausible.",
}
MIN_SILENT_RUN, RATE_FLOOR = 3, 2.5


def bottom_of_cortex(rate, orientation, toc):
    n = len(rate)
    order = np.arange(n) if orientation == "ascending" else np.arange(n)[::-1]
    start = int(np.flatnonzero(order == toc - 1)[0])
    last, silent = start, 0
    for i in range(start, n):
        if rate[order[i]] > RATE_FLOOR:
            last, silent = i, 0
        else:
            silent += 1
            if silent >= MIN_SILENT_RUN:
                break
    return int(order[last]) + 1


def main():
    src = out_root() / "laminar"
    rows = []
    for p in sorted(src.glob("*_laminar.json")):
        r = json.loads(p.read_text())
        pen = r["penetration"]
        rate = np.load(src / f"{pen}_laminar.npz")["rate"]
        l45, basis, reason = r["l45"], "rule", ""
        include = r["tier"] in ("A", "B")
        if pen in OVERRIDE_L45:
            l45, reason = OVERRIDE_L45[pen]
            basis = "override"
        if pen in SET_ASIDE:
            include, reason = False, SET_ASIDE[pen]
        elif not include:
            reason = "tier C: " + re.sub(r"np\.float64\(([^)]*)\)", r"\1", "; ".join(r["flags"][:3]))
        toc = r["top_of_cortex"]
        boc = bottom_of_cortex(rate, r["orientation"], toc) if toc > 0 else -1
        rows.append(dict(penetration=pen, include=include, tier=r["tier"], confidence=r["confidence"],
                         complexity=r["complexity"], orientation=r["orientation"], top_of_cortex=toc,
                         l45=l45, spiking_bottom=boc, deep_limit=DEEP_LIMIT.get(pen, np.nan),
                         spectral_crossing=r["spectral_crossing"],
                         earliest_sink=r["earliest_sink"], basis=basis, reason=reason))
    df = pd.DataFrame(rows).sort_values(["include", "tier", "confidence"], ascending=[False, True, False])
    out = REPO / "metadata" / "laminar_alignment_v1.csv"
    df.to_csv(out, index=False)
    print(df[["penetration", "include", "tier", "orientation", "top_of_cortex", "l45", "deep_limit", "basis"]]
          .to_string(index=False))
    print(f"\nincluded {int(df.include.sum())} of {len(df)} -> {out}")


if __name__ == "__main__":
    main()
