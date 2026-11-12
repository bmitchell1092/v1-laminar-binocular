"""Laminar alignment v4: confidence tiers and compression-scaled compartments (docs/09 AD; Blake, 2026-11-12).

Usage: .venv\\Scripts\\python scripts\\freeze_laminar_alignment_v4.py
Inputs: metadata/laminar_alignment_v3.csv (anchor and markers), the laminar audit's monocular cache via
laminar_audit.markers (evoked-MUA amplitude per contact). Laminar markers and monocular MUA only; nothing binocular is
read. Writes metadata/laminar_alignment_v4.csv.

Rationale (reports/Macaque V1 laminar markers anatomy.md): no histology fixes the sink-to-4C/5 distance, but the layer
budget does constrain the anchor. From the evoked-MUA top (top of layer 2) to the 4C/5 border is about 1.24 mm
(Zhu et al. 2026: L2/3 650, 4A/B 311, 4C 281 um; Cox et al. 2019: 1.2 mm). Cortex is compressed or stretched by
different amounts each day (probe advance, retraction, settling time; Blake), so the ratio of the observed span to
1.24 mm is taken as the penetration's squish factor s. A plausible s confirms the anchor only up to that tolerance, and
the compartments are scaled by s. The bottom of grey matter is not used: an automated call from MUA or CSD fired 0.2 mm
below the anchor in three probes (impossible) and MUA continues 1.3-1.7 mm below the anchor in several J probes.

Rules
- T1 top observed: the evoked-MUA top is not the probe's superficial-end contact, and the contact above it is undriven
  (audit amplitude < 0.2 of the probe's 90th percentile).
- T2 plausible span: s = (top to anchor) / 1.24 mm in [S_MIN, S_MAX] (daily squish or stretch up to ~25-35%, plus one
  contact of rounding).
- T3 cross-checked anchor: two independent boundary estimates (flip with step >= 0.5 rad, sink + 0.3 mm,
  crossing + 0.1 mm) agree within TOL mm, and the anchor is one of them or their mean (v3 anchors unchanged).
- T4 room below: at least 0.49 * s mm of probe below the anchor.
- Tier 1 (confident): T1-T4. Tier 2 (consistent, top off the probe): T3, T4, and the span to the probe's end, a lower
  bound, gives s >= S_MIN. Otherwise tier 0 (left out of the laminar analyses; kept for the all-penetration reference).
- s for labelling: the T2 value in tier 1; 1 otherwise.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import laminar_audit as LA
from v1lb.config import metadata_dir

REPO = Path(__file__).resolve().parents[1]
SPAN_MM = 1.24                      # evoked-MUA top (top of L2) to the 4C/5 border
S_MIN, S_MAX = 0.75, 1.35
TOL = 0.3                           # mm; literature tolerances 0.2-0.3 (Davis et al. 2023 flag >= 0.2)
MIN_STEP, SINK_OFFSET, CROSS_OFFSET = 0.5, 0.3, 0.1
L56_MM = 0.49                       # L5 + L6 (Zhu et al. 2026)


def estimates(r):
    """Boundary estimates in mm relative to the v3 anchor (+ = deeper)."""
    sign = 1 if r.orientation == "ascending" else -1
    rel = lambda c: sign * (c - r.anchor_contact) / 10              # noqa: E731
    est = {}
    if np.isfinite(r.flip_contact) and np.isfinite(r.flip_step_rad) and abs(r.flip_step_rad) >= MIN_STEP:
        est["flip"] = rel(r.flip_contact)
    if np.isfinite(r.sink_contact):
        est["sink"] = rel(r.sink_contact) + SINK_OFFSET
    if np.isfinite(r.crossing_contact):
        est["crossing"] = rel(r.crossing_contact) + CROSS_OFFSET
    return est


def cross_checked(est, basis):
    """T3: the estimates the anchor rests on agree within TOL."""
    if basis.startswith("flip"):
        return any(k in est and abs(est[k] - est["flip"]) <= TOL + 1e-9 for k in ("sink", "crossing"))
    if basis.startswith("sink and crossing"):
        return "sink" in est and "crossing" in est and abs(est["sink"] - est["crossing"]) <= TOL + 1e-9
    return False


def main():
    A = pd.read_csv(metadata_dir() / "laminar_alignment_v3.csv")
    rows = []
    for r in A.itertuples(index=False):
        row = dict(penetration=r.penetration, monkey=r.penetration.split("_")[1], anchor_contact=r.anchor_contact,
                   orientation=r.orientation, top_contact=r.top_contact)
        if not r.depth_usable:
            rows.append(row | dict(tier=0, s=np.nan, reason="no v3 anchor"))
            continue
        _, arr = LA.markers(r.penetration)
        n, amp = arr["n"], arr["amp_n"]
        asc = r.orientation == "ascending"
        pos = lambda c: (c - 1) if asc else (n - c)                     # noqa: E731  contact -> position, 0 = superficial end
        p_top, p_anc = pos(r.top_contact), pos(r.anchor_contact)
        top_obs = bool(p_top > 0 and amp[p_top - 1] < 0.2)
        span = (p_anc - p_top) / 10                                     # top contact centre to the anchor
        s_top = span / SPAN_MM
        below = (n - 1 - p_anc) / 10 + 0.05
        est = estimates(r)
        t3 = cross_checked(est, r.anchor_basis)
        s_lab = s_top if top_obs else 1.0
        t4 = below >= L56_MM * s_lab - 1e-9
        reasons = []
        if not t3:
            reasons.append(f"estimates disagree by > {TOL} mm")
        if not t4:
            reasons.append(f"only {below:.2f} mm of probe below the anchor")
        if top_obs:
            t2 = S_MIN <= s_top <= S_MAX
            if not t2:
                reasons.append(f"span {span:.2f} mm (s = {s_top:.2f})")
            tier = 1 if (t2 and t3 and t4) else 0
        else:
            reasons.append(f"top off the probe (span >= {span:.2f} mm)")
            tier = 2 if (t3 and t4 and s_top >= S_MIN) else 0
            if s_top < S_MIN:
                reasons.append("lower-bound span too short to check")
        rows.append(row | dict(tier=tier, s=round(s_lab, 3) if tier == 1 else 1.0, top_observed=top_obs,
                               span_mm=round(span, 2), s_top=round(s_top, 3), probe_below_anchor_mm=round(below, 2),
                               **{f"est_{k}_mm": round(v, 2) for k, v in est.items()},
                               anchor_basis=r.anchor_basis, reason="; ".join(reasons)))
    out = pd.DataFrame(rows)
    out.to_csv(metadata_dir() / "laminar_alignment_v4.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print(out[["penetration", "tier", "s", "top_observed", "span_mm", "est_flip_mm", "est_sink_mm", "est_crossing_mm",
               "reason"]].sort_values(["tier", "penetration"]).to_string(index=False))
    print(out.groupby(["tier", "monkey"]).size())


if __name__ == "__main__":
    main()
