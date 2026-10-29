"""Laminar alignment v2: the flash-evoked CSD layer-4 sink is the primary landmark (Blake, 2026-10-28).

Usage: .venv\\Scripts\\python scripts\\freeze_laminar_alignment_v2.py
Inputs: output/flash_csd (run_flash_csd.py, register_flash_csd.py --all), output/laminar (reports),
metadata/penetrations.csv, metadata/laminar_alignment_v1.csv, docs/qc/laminar_calls_2026-10-28.csv.
Writes metadata/laminar_alignment_v2.csv.

Rules
- Orientation: majority of SPI (template fit), spectrolaminar peaks and the legacy by-eye sort direction.
- Landmark per evp: strongest sink of the mean flash CSD in 35-65 ms (`evp.dominant_early_sink`), moved into
  bminteroc contact numbers by the registration shift (shift used only if its margin >= 0.1, else 0).
  L4/5 = the boundary just below the sink's deepest contact.
- Usable evp: sink z >= Z_MIN (window-mean peak / SD of the baseline CSD), deepest contact not the probe's
  last contact in the deep direction, sink peak at or below the top of cortex.
- Several usable evps: median L4/5; flagged if their range exceeds 3 contacts.
- No usable evp: SPI L4/5 from v1, only for penetrations included in v1 (tier B, 'spi fallback').
- Tier A: unanimous orientation, usable flash landmark, and either >= 2 consistent evps or one with z >= Z_STRONG.
  Tier B: other usable flash landmarks (orientation at least 2 of 3) or SPI fallback. Tier C: everything else.
  Included: A and B.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from v1lb import laminar as lam
from v1lb.config import metadata_dir, out_root
from v1lb.evp import dominant_early_sink

REPO = Path(__file__).resolve().parents[1]
Z_MIN, Z_STRONG, MAX_RANGE = 5.0, 10.0, 3.0


def main():
    fl = out_root() / "flash_csd"
    reg = pd.read_csv(fl / "flash_csd_registered.csv")
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    v1 = pd.read_csv(metadata_dir() / "laminar_alignment_v1.csv").set_index("penetration")
    calls = pd.read_csv(REPO / "docs" / "qc" / "laminar_calls_2026-10-28.csv").set_index("penetration")
    rows, evp_rows = [], []
    for pen in pens.index:
        votes = [calls.loc[pen, "orientation"], calls.loc[pen, "spectral_orientation"], pens.loc[pen, "legacy_sort"]]
        orientation = max(set(votes), key=votes.count)
        unanimous = votes.count(orientation) == 3
        lr = np.load(out_root() / "laminar" / f"{pen}_laminar.npz")
        n = len(lr["rate"])
        toc, top_in, _ = lam.top_of_cortex(lr["rate"], lr["evoked_t"], orientation)
        sgn = 1 if orientation == "ascending" else -1
        usable = []
        for _, r in reg[reg.penetration == pen].iterrows():
            d = np.load(fl / f"{pen}_{r.evp}.npz")
            csd, t = d["csd"], d["t_ms"]
            peak, top, bottom, val = dominant_early_sink(csd, t)
            base = csd[(t >= -50) & (t < 0)]
            z = abs(val) / np.nanstd(base) if peak >= 0 else 0.0
            shift = int(r["shift"]) if r.reg_margin >= 0.1 else 0
            ok, why = True, ""
            if peak < 0 or z < Z_MIN:
                ok, why = False, f"weak sink (z {z:.1f})"
            else:
                c_peak = int(d["contacts"][peak]) + shift
                deep = (int(d["contacts"][bottom]) if orientation == "ascending" else int(d["contacts"][top])) + shift
                if (orientation == "ascending" and deep >= n) or (orientation == "descending" and deep <= 1):
                    ok, why = False, "sink reaches the deep end of the probe"
                elif not 1 <= deep <= n:
                    ok, why = False, "registered sink falls off the probe"
                elif toc > 0 and sgn * (c_peak - toc) < 0:
                    ok, why = False, "sink above the top of cortex"
            l45 = deep + 0.5 * sgn if ok else np.nan
            evp_rows.append(dict(penetration=pen, evp=r.evp, minutes=r.minutes_from_bminteroc, shift=shift,
                                 reg_margin=r.reg_margin, z=round(z, 1), usable=ok, reason=why, l45=l45))
            if ok:
                usable.append((l45, z, abs(r.minutes_from_bminteroc), r.evp))
        flags = []
        if usable:
            l45s = np.array([u[0] for u in usable])
            spread = float(l45s.max() - l45s.min())
            if spread > MAX_RANGE:
                closest = min(usable, key=lambda u: u[2])
                l45, basis = closest[0], f"flash {closest[3]} (closest; evps disagree by {spread:g})"
                flags.append(f"evps disagree by {spread:g} contacts")
                consistent = False
            else:
                l45, basis = float(np.median(l45s)), "flash " + "+".join(u[3] for u in usable)
                consistent = True
            strong = (len(usable) >= 2 and consistent) or max(u[1] for u in usable) >= Z_STRONG
            tier = "A" if unanimous and consistent and strong else "B"
        elif v1.loc[pen, "include"]:
            l45, basis, tier = v1.loc[pen, "l45"], "spi fallback (no usable flash landmark)", "B"
            flags.append("no usable flash landmark")
        else:
            l45, basis, tier = np.nan, "none", "C"
        if not unanimous:
            flags.append("orientation by 2-of-3 vote: " + ", ".join(votes))
        if top_in:
            flags.append("top contact already in cortex")
        if tier == "A" and flags:
            tier = "B"
        rows.append(dict(penetration=pen, include=tier in ("A", "B"), tier=tier, orientation=orientation,
                         orientation_votes="/".join(v[:4] for v in votes), top_of_cortex=toc, l45=l45,
                         l45_v1=v1.loc[pen, "l45"], basis=basis, n_usable_evps=len(usable),
                         deep_limit=np.nan, flags="; ".join(flags)))
    df = pd.DataFrame(rows).sort_values(["include", "tier", "penetration"], ascending=[False, True, True])
    df.to_csv(metadata_dir() / "laminar_alignment_v2.csv", index=False)
    pd.DataFrame(evp_rows).to_csv(REPO / "docs" / "qc" / "flash_landmarks_v2.csv", index=False)
    print(df.drop(columns=["flags"]).to_string(index=False))
    print(f"\nincluded {int(df.include.sum())} (A {int((df.tier == 'A').sum())}, B {int((df.tier == 'B').sum())}); "
          f"by monkey: {df[df.include].penetration.str[7].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
