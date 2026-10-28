"""Build the per-penetration data files for every included penetration in the laminar alignment.

Usage: .venv\\Scripts\\python scripts\\build_datasets.py [penetration ...] [--overwrite]
Writes <out_root>/datasets/<pen>.nc and <pen>_presentations.parquet.
"""

import sys
import time
import traceback

import pandas as pd

from v1lb.config import metadata_dir, out_root
from v1lb.dataset import ALIGNMENT_FILE, build_penetration

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    overwrite = "--overwrite" in sys.argv
    align = pd.read_csv(metadata_dir() / ALIGNMENT_FILE)
    pens = args or align.loc[align.include, "penetration"].tolist()
    out = out_root() / "datasets"
    for pen in pens:
        if (out / f"{pen}.nc").exists() and not overwrite:
            print(f"{pen}: exists, skipped", flush=True)
            continue
        t = time.time()
        try:
            ds = build_penetration(pen, out)
            size = (out / f"{pen}.nc").stat().st_size / 1e6
            print(f"{pen}: {time.time() - t:.0f}s, {ds.sizes['presentation']} presentations x "
                  f"{ds.sizes['channel']} channels, {size:.0f} MB", flush=True)
        except Exception:
            print(f"{pen}: FAILED\n{traceback.format_exc()}", flush=True)
