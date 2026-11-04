"""Laminar alignment v3: depth from a cross-checked input/deep anchor (docs/09 W and X; Blake, 2026-11-03).

Usage: .venv\\Scripts\\python scripts\\freeze_laminar_alignment_v3.py
Inputs: docs/qc/laminar_audit.csv (scripts/laminar_audit.py), output/laminar/{pen}_laminar.npz (contact count),
metadata/laminar_alignment_v2.csv (orientation used to build the datasets). Writes metadata/laminar_alignment_v3.csv.
Laminar markers only; nothing binocular is read.

Rules
- Orientation: the audit's majority of three votes (phase template, spectrolaminar peaks, SPI vs contact).
- Three estimates of the input/deep boundary, in mm along the probe: the phase flip (the Davis template boundary,
  if on the probe and its step >= MIN_STEP rad); the earliest monocular sink + SINK_OFFSET; the alpha/beta-gamma
  crossing + CROSS_OFFSET (offsets = the medians over the 14 clear flips, docs/09 W).
- Anchor: the flip if the sink or crossing estimate lies within TOL of it; else the mean of the sink and crossing
  estimates if they lie within TOL of each other; else none (left out of the depth analyses). Rounded to the nearest
  boundary between contacts (ties to the deeper side).
- Top of cortex: the evoked-MUA top (the spontaneous-rate top is kept for reference).
- Deep cap: contacts more than DEEP_CAP mm below the anchor are excluded (no white-matter call; second sheets).
- Sensitivity anchor: the boundary 2.5 contacts below the earliest sink (input group = 5 contacts centred on it).
"""

from pathlib import Path

import numpy as np
import pandas as pd

from v1lb.config import metadata_dir, out_root

REPO = Path(__file__).resolve().parents[1]
MIN_STEP, TOL = 0.5, 0.4               # rad; mm
SINK_OFFSET, CROSS_OFFSET = 0.3, 0.1   # mm below the marker
DEEP_CAP = 0.7                         # mm below the anchor
# Blake delegated the four cards flagged in W; each decision keeps the rule's outcome (docs/09 X)
NOTES = {
    "221215_J_eA": "orientation as v2 (the audit's 'legacy' was the laminar-report call); flip + crossing + earliest "
                   "MUA agree; the sink (0.75 mm below) and reversed step lie beyond the deep cap",
    "221214_J_eA": "sink and crossing agree with an on-probe phase step the template missed; the alpha/beta "
                   "inversion comes from contacts beyond the deep cap",
    "221118_J_eA": "left out: flip, sink and crossing spread over 1.1 mm",
    "211217_B_eA": "no second sheet: the deep early zone is low-amplitude and patchy, no reversed phase step; "
                   "beyond the deep cap",
    "211025_B_eA": "flip 0.25 mm below the top of cortex (not L4/5); sink and crossing agree 1.1 mm deeper",
    "221109_J_eA": "flip in unresponsive tissue at the deep end; sink, crossing and earliest MUA agree 2.3 mm higher",
}


def to_boundary(x: float, sign: int) -> float:
    """Nearest boundary between contacts (k + 0.5); ties go to the deeper side."""
    lo = np.floor(x - 0.5) + 0.5
    if np.isclose(x - lo, 0.5):
        return lo + 1 if sign > 0 else lo
    return lo if x - lo < 0.5 else lo + 1


def anchor(r) -> tuple[float, str]:
    """Anchor in mm relative to the audit's template reversal (+ = deeper), and its basis."""
    est = {}
    if r.reversal_on_probe and abs(r.reversal_step_rad) >= MIN_STEP:
        est["flip"] = 0.0
    if np.isfinite(r.sink_depth_mm):
        est["sink"] = r.sink_depth_mm + SINK_OFFSET
    if np.isfinite(r.crossing_depth_mm):
        est["crossing"] = r.crossing_depth_mm + CROSS_OFFSET
    near = lambda a, b: a in est and b in est and abs(est[a] - est[b]) <= TOL + 1e-9   # noqa: E731
    agree = [k for k in ("sink", "crossing") if near("flip", k)]
    if agree:
        return 0.0, "flip (" + ", ".join(agree) + " agree)"
    if near("sink", "crossing"):
        return (est["sink"] + est["crossing"]) / 2, "sink and crossing (flip absent or disagrees)"
    return np.nan, "none: no two estimates within %.1f mm" % TOL


def main():
    W = pd.read_csv(REPO / "docs" / "qc" / "laminar_audit.csv")
    V2 = pd.read_csv(metadata_dir() / "laminar_alignment_v2.csv").set_index("penetration")
    rows = []
    for _, r in W.iterrows():
        pen = r.penetration
        n = len(np.load(out_root() / "laminar" / f"{pen}_laminar.npz")["rate"])
        sign = 1 if r.orientation == "ascending" else -1           # contact number increases with depth if ascending
        at = lambda mm: r.reversal + sign * 10 * mm                   # noqa: E731  (mm from the reversal -> contact)
        off, basis = anchor(r)
        a = to_boundary(at(off), sign) if np.isfinite(off) else np.nan
        flip_ok = bool(r.reversal_on_probe and abs(r.reversal_step_rad) >= MIN_STEP)
        rows.append(dict(
            penetration=pen, depth_usable=bool(np.isfinite(a)), orientation=r.orientation, orientation_votes=r.votes,
            orientation_v2=V2.orientation.get(pen, ""), n_contacts=n, anchor_contact=a, anchor_basis=basis,
            flip_contact=r.reversal if flip_ok else np.nan, flip_step_rad=r.reversal_step_rad if flip_ok else np.nan,
            sink_contact=r.sink, crossing_contact=round(at(r.crossing_depth_mm), 1) if np.isfinite(r.crossing_depth_mm)
            else np.nan, top_contact=r.top_mua, top_rate_contact=r.top_rate,
            top_to_anchor_mm=round(sign * (a - r.top_mua) / 10, 2) if np.isfinite(a) else np.nan,
            deep_cap_mm=DEEP_CAP,
            sink_anchor_contact=r.sink + sign * 2.5 if np.isfinite(r.sink) else np.nan,
            note=NOTES.get(pen, "")))
    A = pd.DataFrame(rows)
    A.to_csv(metadata_dir() / "laminar_alignment_v3.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(A.drop(columns=["note", "orientation_votes"]).to_string(index=False))
    u = A[A.depth_usable]
    print(f"\n{len(u)}/{len(A)} anchored: {u.anchor_basis.str.startswith('flip').sum()} on the flip, "
          f"{(~u.anchor_basis.str.startswith('flip')).sum()} on the sink and crossing; "
          f"orientation differs from v2 in {', '.join(A.penetration[A.orientation != A.orientation_v2]) or 'none'}")
    print("top of cortex to anchor (mm): median %.2f, IQR %.2f-%.2f" % tuple(u.top_to_anchor_mm.quantile([.5, .25, .75])))


if __name__ == "__main__":
    main()
