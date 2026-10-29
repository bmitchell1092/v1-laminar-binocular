"""Import the 2023 Kilosort/phy curation for every depth-usable penetration that has a sort.

Usage: .venv\Scripts\python scripts\run_units.py [penetration ...]
Writes <out_root>/units/<pen>_units.parquet (one row per good/mua cluster of the penetration's bank, with
isolation metrics and depth coordinates) and <pen>_spikes.parquet (cluster_id, tick on the NSP clock), and
docs/qc/kilosort_offsets.csv (measured spike-time offset per sort).
"""

import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from v1lb.config import metadata_dir, out_root
from v1lb.dataset import ALIGNMENT_FILE, channel_table
from v1lb.units import sort_dir, units

REPO = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    al = pd.read_csv(metadata_dir() / ALIGNMENT_FILE).set_index("penetration")
    todo = [a for a in sys.argv[1:] if not a.startswith("--")] or al.index[al.depth_usable].tolist()
    out = out_root() / "units"
    out.mkdir(parents=True, exist_ok=True)
    qc_path = REPO / "docs" / "qc" / "kilosort_offsets.csv"
    qc = pd.read_csv(qc_path).set_index("penetration") if qc_path.exists() else pd.DataFrame()
    for pen in todo:
        p, a = pens.loc[pen], al.loc[pen]
        if sort_dir(p.session, p.bminteroc) is None:
            print(f"{pen}: no sort", flush=True)
            continue
        t = time.time()
        try:
            tab, spikes, rep = units(p.session, p.bminteroc, bank=p.bank)
        except Exception:
            print(f"{pen}: FAILED\n{traceback.format_exc()}", flush=True)
            continue
        ch = channel_table(32, a.orientation, float(a.l45), int(a.top_of_cortex),
                           float(a.deep_limit) if pd.notna(a.deep_limit) else np.nan).set_index("contact")
        for col in ("depth_l45", "depth_um", "depth_below_toc", "in_cortex", "beyond_deep_limit"):
            tab[col] = ch.loc[tab.peak_contact.fillna(tab.contact).astype(int), col].to_numpy()
        tab.insert(0, "penetration", pen)
        tab.to_parquet(out / f"{pen}_units.parquet")
        pd.DataFrame({"cluster_id": np.concatenate([np.full(len(v), k) for k, v in spikes.items()]),
                      "tick": np.concatenate(list(spikes.values()))}).to_parquet(out / f"{pen}_spikes.parquet")
        qc.loc[pen, ["sort", "offset_samples", "consistent", "unit_lags"]] = [
            f"{p.session}_{p.bminteroc}", rep["offset"], rep["consistent"], str(rep["lags"])]
        g = tab[tab.group == "good"]
        print(f"{pen}: {time.time() - t:.0f}s offset {rep['offset']}; good {len(g)}, well isolated "
              f"{int(tab.well_isolated.sum())}", flush=True)
    qc.index.name = "penetration"
    qc.sort_index().to_csv(qc_path)
