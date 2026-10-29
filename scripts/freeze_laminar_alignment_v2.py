"""Laminar alignment v2: the flash-evoked CSD layer-4 sink is the primary landmark (Blake, 2026-10-28).

Usage: .venv\\Scripts\\python scripts\\freeze_laminar_alignment_v2.py
Inputs: output/flash_csd (run_flash_csd.py, register_flash_csd.py --all), output/laminar (reports),
metadata/penetrations.csv, metadata/laminar_alignment_v1.csv, docs/qc/laminar_calls_2026-10-28.csv.
Writes metadata/laminar_alignment_v2.csv and docs/qc/flash_landmarks_v2.csv.

Rules
- Orientation: SPI vs depth is the arbiter (Blake): the SPI template call when its fit margin >= SPI_DECISIVE;
  otherwise the majority of SPI, spectrolaminar peaks and the legacy by-eye sort. A decisive SPI call that both
  other markers contradict is kept as the orientation but the penetration is set aside (tier C, not depth-usable).
- Flash landmark per evp: outlier contacts interpolated and 3-point Hamming smoothing before the CSD (the raw CSD
  alternates sink/source between neighbouring contacts in many penetrations: contact gain mismatch, 05 C12); then
  the strongest sink of the mean CSD in 35-65 ms (`evp.dominant_early_sink`), moved into bminteroc contact numbers
  by the registration shift (used only if its margin >= 0.1). L4/5 = the boundary just below the sink's deepest
  contact. Usable: z >= Z_MIN (window-mean peak / SD of the baseline CSD), not at the probe's deep end, on the probe
  after the shift, peak not above the top of cortex. Several usable evps: median; if their range exceeds MAX_RANGE,
  the evp closest in time to bminteroc (flagged).
- Plausible depth: L4/5 between DIST_MIN and DIST_MAX contacts below the top of cortex (spontaneous MUA; ~L1/L2),
  or at most DIST_MAX when the top contact is already in cortex (the distance is then a lower bound). Macaque V1
  L4/5 lies ~1.1-1.3 mm below the pia; the window allows for top-of-cortex error and probe angle.
- L4/5 choice: flash if plausible; flash if implausible but the SPI step agrees within AGREE contacts (two
  independent markers outweigh the anatomical prior); else the SPI step if plausible; else none.
- Tier A: flash plausible, SPI agrees, unanimous orientation, consistent evps (>= 2, or one with z >= Z_STRONG),
  no flags. Tier B: any other L4/5. Tier C: no L4/5 or orientation conflict. Included (layer labels): A and B.
  depth_usable (depth below the top of cortex, the fallback axis): every penetration without an orientation conflict.
- deep_limit carried over from v1 (confirmed folds) unless L4/5 lies beyond it.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from v1lb import laminar as lam
from v1lb.config import metadata_dir, out_root
from v1lb.evp import dominant_early_sink

REPO = Path(__file__).resolve().parents[1]
Z_MIN, Z_STRONG, MAX_RANGE = 5.0, 10.0, 3.0
SPI_DECISIVE = 0.1
DIST_MIN, DIST_MAX = 5.0, 17.0
AGREE = 3.0


def flash_landmark(pen, r, orientation, n, toc):
    """L4/5 from one evp (bminteroc contact numbers) and the per-evp QC row."""
    d = np.load(out_root() / "flash_csd" / f"{pen}_{r.evp}.npz")
    lfp, t, con = d["lfp"], d["t_ms"], d["contacts"]
    bad = lam.bad_channels(lfp[(t >= 20) & (t <= 150)])
    csd = lam.csd(lfp, smooth=True, bad=bad)
    peak, top, bottom, val = dominant_early_sink(csd, t)
    z = abs(val) / np.nanstd(csd[(t >= -50) & (t < 0)]) if peak >= 0 else 0.0
    shift = int(r["shift"]) if r.reg_margin >= 0.1 else 0
    sgn = 1 if orientation == "ascending" else -1
    ok, why, l45 = True, "", np.nan
    if peak < 0 or z < Z_MIN:
        ok, why = False, f"weak sink (z {z:.1f})"
    else:
        ends = (int(con[top]) + shift, int(con[bottom]) + shift)
        deep = max(ends) if sgn == 1 else min(ends)
        c_peak = int(con[peak]) + shift
        if (sgn == 1 and deep >= n) or (sgn == -1 and deep <= 1):
            ok, why = False, "sink reaches the deep end of the probe"
        elif not 1 <= deep <= n:
            ok, why = False, "registered sink falls off the probe"
        elif toc > 0 and sgn * (c_peak - toc) < 0:
            ok, why = False, "sink above the top of cortex"
        else:
            l45 = deep + 0.5 * sgn
    return dict(penetration=pen, evp=r.evp, minutes=r.minutes_from_bminteroc, shift=shift, reg_margin=r.reg_margin,
                bad_contacts=",".join(str(int(c) + shift) for c in con[bad]), z=round(float(z), 1),
                usable=ok, reason=why, l45=l45)


def main():
    fl = out_root() / "flash_csd"
    reg = pd.read_csv(fl / "flash_csd_registered.csv")
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    v1 = pd.read_csv(metadata_dir() / "laminar_alignment_v1.csv").set_index("penetration")
    calls = pd.read_csv(REPO / "docs" / "qc" / "laminar_calls_2026-10-28.csv").set_index("penetration")
    rows, evp_rows = [], []
    for pen in pens.index:
        c = calls.loc[pen]
        votes = [c["orientation"], c["spectral_orientation"], pens.loc[pen, "legacy_sort"]]
        spi_decisive = c["spi_fit_margin"] >= SPI_DECISIVE
        orientation = votes[0] if spi_decisive else max(set(votes), key=votes.count)
        unanimous = votes.count(orientation) == 3
        conflict = votes.count(orientation) == 1
        sgn = 1 if orientation == "ascending" else -1

        lr = np.load(out_root() / "laminar" / f"{pen}_laminar.npz")
        n = len(lr["rate"])
        toc, top_in, _ = lam.top_of_cortex(lr["rate"], lr["evoked_t"], orientation)

        def dist(l45):
            return sgn * (l45 - toc)

        def plausible(l45):
            return np.isfinite(l45) and dist(l45) <= DIST_MAX and (top_in or dist(l45) >= DIST_MIN)

        ev = [flash_landmark(pen, r, orientation, n, toc) for _, r in reg[reg.penetration == pen].iterrows()]
        evp_rows += ev
        usable = [e for e in ev if e["usable"]]
        flags = []
        flash, spread, consistent, strong = np.nan, np.nan, False, False
        if usable:
            l45s = np.array([e["l45"] for e in usable])
            spread = float(l45s.max() - l45s.min())
            consistent = spread <= MAX_RANGE
            if consistent:
                flash = float(np.median(l45s))
            else:
                flash = min(usable, key=lambda e: abs(e["minutes"]))["l45"]
                flags.append(f"evps disagree by {spread:g} contacts (closest evp used)")
            strong = (len(usable) >= 2 and consistent) or max(e["z"] for e in usable) >= Z_STRONG
        spi = float(c["l45"])
        spi_agrees = np.isfinite(flash) and np.isfinite(spi) and abs(spi - flash) <= AGREE

        if plausible(flash):
            l45, basis = flash, "flash"
        elif spi_agrees:
            l45, basis = flash, "flash (SPI agrees)"
            flags.append(f"L4/5 {dist(flash):g} contacts below top of cortex (outside {DIST_MIN:g}-{DIST_MAX:g})")
        elif plausible(spi):
            l45, basis = spi, "spi (flash " + ("implausible depth)" if np.isfinite(flash) else "unusable)")
        else:
            l45, basis = np.nan, "none"
        if not unanimous:
            flags.append(("decisive SPI contradicted by spectral and legacy: " if conflict else
                          f"orientation by {'SPI (decisive)' if spi_decisive else '2-of-3 vote'}: ") + ", ".join(votes))
        if top_in:
            flags.append("top contact already in cortex")

        deep_limit = v1.loc[pen, "deep_limit"]
        if pd.notna(deep_limit) and np.isfinite(l45) and sgn * (l45 - deep_limit) > 0:
            flags.append(f"v1 deep limit {deep_limit:g} dropped (L4/5 below it)")
            deep_limit = np.nan

        if conflict or not np.isfinite(l45):
            tier = "C"
        elif basis == "flash" and spi_agrees and unanimous and consistent and strong and not flags:
            tier = "A"
        else:
            tier = "B"
        rows.append(dict(penetration=pen, include=tier in ("A", "B"), depth_usable=not conflict, tier=tier,
                         orientation=orientation, orientation_votes="/".join(v[:4] for v in votes),
                         spi_margin=round(float(c["spi_fit_margin"]), 3), top_of_cortex=toc,
                         top_contact_in_cortex=bool(top_in), l45=l45,
                         l45_minus_toc=dist(l45) if np.isfinite(l45) else np.nan,
                         l45_flash=flash, l45_spi=spi, l45_v1=v1.loc[pen, "l45"], basis=basis,
                         n_usable_evps=len(usable), evp_spread=spread, deep_limit=deep_limit, flags="; ".join(flags)))
    df = pd.DataFrame(rows).sort_values(["include", "tier", "penetration"], ascending=[False, True, True])
    df.to_csv(metadata_dir() / "laminar_alignment_v2.csv", index=False)
    pd.DataFrame(evp_rows).to_csv(REPO / "docs" / "qc" / "flash_landmarks_v2.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["flags", "orientation_votes", "l45_v1"]).to_string(index=False))
    inc = df[df.include]
    print(f"\nincluded {len(inc)} (A {int((df.tier == 'A').sum())}, B {int((df.tier == 'B').sum())}); "
          f"by monkey: {inc.penetration.str[7].value_counts().to_dict()}; basis: {inc.basis.value_counts().to_dict()}; "
          f"depth-usable {int(df.depth_usable.sum())}")


if __name__ == "__main__":
    main()
