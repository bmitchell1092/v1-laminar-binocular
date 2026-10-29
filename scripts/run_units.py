"""Import single units for every depth-usable penetration: the 2023 Kilosort/phy curation where it is usable,
otherwise the Kilosort 4 sort in <out_root>/kilosort4/<pen>/ (run_ks4.py; pilot showed equal yield, 05 C14).

Usage: .venv\\Scripts\\python scripts\\run_units.py [penetration ...]
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

from v1lb.config import metadata_dir, out_root, recording_files
from v1lb.dataset import ALIGNMENT_FILE
from v1lb.io import NsxFile
from v1lb.units import BROKEN_SORTS, load_ks4, sort_dir, units, with_depth

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
        rec, ks4 = f"{p.session}_{p.bminteroc}", out_root() / "kilosort4" / pen
        if sort_dir(p.session, p.bminteroc) is not None and rec not in BROKEN_SORTS:
            sorter, sort = "kilosort2_phy2023", None
        elif (ks4 / "spike_times.npy").exists():
            sorter = "kilosort4"
            sort = load_ks4(ks4, NsxFile(recording_files(p.session, p.bminteroc)[".ns6"]).labels)
        else:
            print(f"{pen}: no usable sort", flush=True)
            continue
        t = time.time()
        try:
            tab, spikes, rep = units(p.session, p.bminteroc, bank=p.bank, sort=sort)
        except Exception:
            print(f"{pen}: FAILED\n{traceback.format_exc()}", flush=True)
            continue
        tab = with_depth(tab, a)
        tab.insert(0, "penetration", pen)
        tab.insert(1, "sorter", sorter)
        tab.to_parquet(out / f"{pen}_units.parquet")
        pd.DataFrame({"cluster_id": np.concatenate([np.full(len(v), k) for k, v in spikes.items()]),
                      "tick": np.concatenate(list(spikes.values()))}).to_parquet(out / f"{pen}_spikes.parquet")
        qc.loc[pen, ["sorter", "sort", "offset_samples", "measured_samples", "consistent", "unit_lags"]] = [
            sorter, rec, rep["offset"], rep["measured"], rep["consistent"], str(rep["lags"])]
        g = tab[tab.group == "good"]
        print(f"{pen}: {time.time() - t:.0f}s offset {rep['offset']}; good {len(g)}, well isolated "
              f"{int(tab.well_isolated.sum())}", flush=True)
    qc.index.name = "penetration"
    qc.sort_index().to_csv(qc_path)
