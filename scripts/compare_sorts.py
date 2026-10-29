"""Kilosort 4 pilot vs the 2023 phy curation, same penetration, same metrics.

Usage: .venv\\Scripts\\python scripts\\compare_sorts.py <penetration> [...] [--tag name]
Needs <out_root>/kilosort4/<pen>[_<tag>]/ (run_ks4.py) and <out_root>/units/<pen>_units.parquet (run_units.py).
Writes units.parquet / spikes.parquet into the Kilosort 4 folder and appends to docs/qc/ks4_pilot.csv.
Matching: a Kilosort 4 unit matches the 2023 cluster that shares the most of its spikes (within 10 samples).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from v1lb import units as U
from v1lb.config import metadata_dir, out_root, recording_files
from v1lb.dataset import ALIGNMENT_FILE
from v1lb.io import NsxFile

REPO = Path(__file__).resolve().parents[1]


def shared_fraction(a: np.ndarray, b: np.ndarray, tol: int = 10) -> float:
    """Fraction of spikes in a with a spike of b within tol ticks."""
    if len(a) == 0 or len(b) == 0:
        return 0.0
    i = np.clip(np.searchsorted(b, a), 1, len(b) - 1)
    return float(np.mean(np.minimum(np.abs(b[i] - a), np.abs(b[i - 1] - a)) <= tol))


def summary(tab: pd.DataFrame, label: str) -> dict:
    wi = tab[tab.well_isolated & tab.in_cortex]
    return {f"{label}_clusters": len(tab), f"{label}_well_isolated": len(wi),
            f"{label}_wi_over_1hz": int((wi.rate_hz >= 1).sum()),
            f"{label}_wi_contacts": int(wi.peak_contact.nunique()),
            f"{label}_median_snr": round(float(tab.snr.median()), 2),
            f"{label}_median_isi": round(float(tab.isi_contamination.median()), 3),
            f"{label}_median_ampcut": round(float(tab.amplitude_cutoff.median()), 3)}


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    tag = opts.get("tag", "")
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    al = pd.read_csv(metadata_dir() / ALIGNMENT_FILE).set_index("penetration")
    qc_path = REPO / "docs" / "qc" / "ks4_pilot.csv"
    qc = pd.read_csv(qc_path) if qc_path.exists() else pd.DataFrame()
    for pen in args:
        p = pens.loc[pen]
        res = out_root() / "kilosort4" / (pen + (f"_{tag}" if tag else ""))
        ns = NsxFile(recording_files(p.session, p.bminteroc)[".ns6"])
        tab, spikes, rep = U.units(p.session, p.bminteroc, bank=p.bank, sort=U.load_ks4(res, ns.labels))
        tab = U.with_depth(tab, al.loc[pen])
        tab.insert(0, "penetration", pen)
        tab.to_parquet(res / "units.parquet")
        pd.DataFrame({"cluster_id": np.concatenate([np.full(len(v), k) for k, v in spikes.items()]),
                      "tick": np.concatenate(list(spikes.values()))}).to_parquet(res / "spikes.parquet")
        row = dict(penetration=pen, tag=tag, offset=rep["offset"], offset_measured=rep["measured"],
                   offset_consistent=rep["consistent"]) | summary(tab, "ks4")
        old_path = out_root() / "units" / f"{pen}_units.parquet"
        if old_path.exists():
            old = pd.read_parquet(old_path)
            row |= summary(old, "ks2023")
            S = pd.read_parquet(out_root() / "units" / f"{pen}_spikes.parquet")
            old_sp = {c: np.sort(g.tick.to_numpy()) for c, g in S.groupby("cluster_id")}
            match = []
            for _, u in tab[tab.well_isolated].iterrows():
                best = max(old_sp, key=lambda c: shared_fraction(spikes[u.cluster_id], old_sp[c]))
                o = old.set_index("cluster_id").loc[best]
                match.append(dict(ks4=u.cluster_id, ks2023=best, shared=shared_fraction(spikes[u.cluster_id], old_sp[best]),
                                  ks2023_group=o.group, ks2023_wi=bool(o.well_isolated)))
            m = pd.DataFrame(match)
            if len(m):
                row |= dict(ks4_wi_matched_2023=int((m.shared >= 0.5).sum()),
                            ks4_wi_matched_2023_wi=int(((m.shared >= 0.5) & m.ks2023_wi).sum()))
                m.to_csv(res / "match_2023.csv", index=False)
        print(pd.Series(row).to_string(), flush=True)
        qc = pd.concat([qc[(qc.get("penetration") != pen) | (qc.get("tag").fillna("") != tag)] if len(qc) else qc,
                        pd.DataFrame([row])], ignore_index=True)
    qc.to_csv(qc_path, index=False)
