"""Laminar alignment v5: the 4C/5 boundary at the bottom of the initial sink (docs/09 AF and v5; Blake, 2026-11-13).

Usage: .venv\\Scripts\\python scripts\\freeze_laminar_alignment_v5.py
Inputs: metadata/laminar_alignment_v3.csv (contact-numbered markers: earliest sink, phase flip, crossing, top of
cortex), the laminar audit's monocular cache (<out_root>/laminar_audit/{pen}_mono.npz) and laminar_audit.markers
(evoked-MUA amplitude). Laminar markers and monocular responses only. Writes metadata/laminar_alignment_v5.csv.

Rules
- Sink bottom: the mean monocular CSD (one eye at >= 45%, eyes pooled, baseline-subtracted) at the time of the initial
  sink's peak (20-70 ms) at the earliest-sink contact; walking deeper from that contact, the boundary between the last
  sink contact and the first source contact (Maier-lab convention: the 4C/5 border is the bottom of the initial sink).
- Cross-checks: the phase flip (|step| >= 0.5 rad) and the alpha/beta-gamma crossing, each shifted by its median
  offset from the sink bottom over the clear flips (|step| >= 1 rad, as in audit W).
- Anchor: the sink bottom if the flip or crossing estimate lies within TOL of it; else the flip estimate if the
  crossing estimate agrees within TOL; else none.
- Tiers, s, T1-T4: as v4 (freeze_laminar_alignment_v4), with the v5 anchor.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import laminar_audit as LA
from freeze_laminar_alignment_v4 import L56_MM, S_MAX, S_MIN, SPAN_MM, TOL
from v1lb.config import metadata_dir, out_root
from v1lb.dataset import contact_spacing_mm

REPO = Path(__file__).resolve().parents[1]
MIN_STEP, CLEAR_STEP = 0.5, 1.0


def sink_bottom(pen, orientation, sink_contact):
    """Contact number (k + 0.5) of the boundary below the initial sink, and the sink's peak time (ms)."""
    M = np.load(out_root() / "laminar_audit" / f"{pen}_mono.npz")
    t, c = M["t"], M["contact"]
    csd = (M["csd_contra"] + M["csd_ipsi"]) / 2
    csd = csd - np.nanmean(csd[t < 0], axis=0)
    sign = 1 if orientation == "ascending" else -1
    o = np.argsort(sign * c)                                   # superficial -> deep
    cs, cc = csd[:, o], c[o]
    k = int(np.flatnonzero(cc == sink_contact)[0])
    w = (t >= 20) & (t < 70)
    tpk = float(t[w][np.nanargmin(cs[w, k])])
    prof = cs[t == tpk][0]
    j = k
    while j + 1 < len(prof) and np.isfinite(prof[j + 1]) and prof[j + 1] < 0:
        j += 1
    if j + 1 >= len(prof):
        return np.nan, tpk
    return float(cc[j] + cc[j + 1]) / 2, tpk


def main():
    A = pd.read_csv(metadata_dir() / "laminar_alignment_v3.csv")
    rows = []
    for r in A.itertuples(index=False):
        h = contact_spacing_mm(r.penetration)
        sign = 1 if r.orientation == "ascending" else -1
        sb, tpk = sink_bottom(r.penetration, r.orientation, r.sink_contact) if np.isfinite(r.sink_contact) \
            else (np.nan, np.nan)
        rows.append(dict(r._asdict(), h=h, sign=sign, sink_bottom=sb, sink_peak_ms=tpk))
    A = pd.DataFrame(rows)
    flip_ok = A.flip_contact.notna() & (A.flip_step_rad.abs() >= MIN_STEP)
    clear = A.flip_contact.notna() & (A.flip_step_rad.abs() >= CLEAR_STEP)
    up = lambda a, b, r: (r.sign * (b - a)) * r.h               # noqa: E731  mm by which contact a lies above b
    off_flip = np.nanmedian([up(r.sink_bottom, r.flip_contact, r) for r in A[clear].itertuples()])
    off_cross = np.nanmedian([up(r.sink_bottom, r.crossing_contact, r) for r in A[clear].itertuples()])
    print(f"offsets over {int(clear.sum())} clear flips: sink bottom {off_flip:.3f} mm above the flip, "
          f"{off_cross:.3f} mm above the crossing")
    out = []
    for r, fok in zip(A.itertuples(), flip_ok):
        est = {}
        if np.isfinite(r.sink_bottom):
            est["sink"] = r.sink_bottom
        if fok:
            est["flip"] = r.flip_contact - r.sign * off_flip / r.h
        if np.isfinite(r.crossing_contact):
            est["crossing"] = r.crossing_contact - r.sign * off_cross / r.h
        near = lambda a, b: a in est and b in est and abs(est[a] - est[b]) * r.h <= TOL + 1e-9   # noqa: E731
        agree = [k for k in ("flip", "crossing") if near("sink", k)]
        if agree:
            anchor, basis = est["sink"], "sink bottom (" + ", ".join(agree) + " agree)"
        elif near("flip", "crossing"):
            x = est["flip"]
            lo = np.floor(x - 0.5) + 0.5
            anchor = lo if x - lo < 0.5 else lo + 1
            basis = "flip and crossing (sink disagrees)"
        else:
            anchor, basis = np.nan, "none: no two estimates within %.1f mm" % TOL
        row = dict(penetration=r.penetration, monkey=r.penetration.split("_")[1], orientation=r.orientation,
                   spacing_mm=r.h, top_contact=r.top_contact, anchor_contact=anchor, anchor_basis=basis,
                   anchor_v4_minus_v5_mm=np.nan, sink_bottom=r.sink_bottom, sink_peak_ms=r.sink_peak_ms,
                   **{f"est_{k}_mm": (round(r.sign * (v - anchor) * r.h, 2) if np.isfinite(anchor) else np.nan)
                      for k, v in est.items()})
        if not np.isfinite(anchor):
            out.append(row | dict(tier=0, s=np.nan, reason="no anchor"))
            continue
        _, arr = LA.markers(r.penetration)
        n, amp = arr["n"], arr["amp_n"]
        pos = (lambda c: c - 1) if r.orientation == "ascending" else (lambda c: n - c)
        p_top, p_anc = pos(r.top_contact), pos(anchor)
        top_obs = bool(p_top > 0 and amp[p_top - 1] < 0.2)
        span = (p_anc - p_top) * r.h
        s_top = span / SPAN_MM
        below = (n - 1 - p_anc) * r.h + r.h / 2
        s_lab = s_top if top_obs else 1.0
        t4 = below >= L56_MM * s_lab - 1e-9
        reasons = []
        if not t4:
            reasons.append(f"only {below:.2f} mm of probe below the anchor")
        if top_obs:
            t2 = S_MIN <= s_top <= S_MAX
            if not t2:
                reasons.append(f"span {span:.2f} mm (s = {s_top:.2f})")
            tier = 1 if (t2 and t4) else 0
        else:
            reasons.append(f"top off the probe (span >= {span:.2f} mm)")
            tier = 2 if (t4 and S_MIN <= s_top <= S_MAX) else 0       # the lower bound must lie in the plausible range
            if s_top < S_MIN:
                reasons.append("lower-bound span too short to check")
            elif s_top > S_MAX:
                reasons.append("lower-bound span already too long")
        out.append(row | dict(tier=tier, s=round(s_lab, 3) if tier == 1 else 1.0, top_observed=top_obs,
                              span_mm=round(span, 2), s_top=round(s_top, 3), probe_below_anchor_mm=round(below, 2),
                              reason="; ".join(reasons)))
    V = pd.DataFrame(out)
    V4 = pd.read_csv(metadata_dir() / "laminar_alignment_v4.csv").set_index("penetration")
    V["anchor_v4_minus_v5_mm"] = [round(r.sign * (V4.anchor_contact.get(r.penetration, np.nan) - r.anchor_contact)
                                        * r.h, 2) for r in A.assign(anchor_contact=V.anchor_contact).itertuples()]
    V["tier_v4"] = V.penetration.map(V4.tier)
    V.to_csv(metadata_dir() / "laminar_alignment_v5.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print(V[["penetration", "tier", "tier_v4", "s", "spacing_mm", "anchor_contact", "anchor_basis", "span_mm",
             "anchor_v4_minus_v5_mm", "reason"]].sort_values(["tier", "penetration"]).to_string(index=False))
    print(V.groupby(["tier", "monkey"]).size())


if __name__ == "__main__":
    main()
