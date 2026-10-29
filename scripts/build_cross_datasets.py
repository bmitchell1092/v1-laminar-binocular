"""Datasets for the dual-probe cross-recordings: each probe during the bminteroc run whose grating sat on the
other probe's RF (held-out data for the candidate-1 confirmation, docs/10). Written to <out_root>/datasets_cross/
so no discovery script picks them up.

Usage: .venv\\Scripts\\python scripts\\build_cross_datasets.py
"""

import time
import traceback

import pandas as pd

from v1lb.config import metadata_dir, out_root
from v1lb.dataset import ALIGNMENT_FILE, build_penetration

if __name__ == "__main__":
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    al = pd.read_csv(metadata_dir() / ALIGNMENT_FILE).set_index("penetration")
    out = out_root() / "datasets_cross"
    for session, g in pens.groupby("session"):
        if g.bank.nunique() < 2:
            continue
        runs = set(g.bminteroc)
        for pen, row in g.iterrows():
            if not al.loc[pen, "depth_usable"]:
                continue
            for rec in sorted(runs - {row.bminteroc}):
                t = time.time()
                try:
                    ds = build_penetration(pen, out, recording=rec)
                    print(f"{pen} in {rec}: {time.time() - t:.0f}s, {ds.sizes['presentation']} presentations", flush=True)
                except Exception:
                    print(f"{pen} in {rec}: FAILED\n{traceback.format_exc()}", flush=True)
