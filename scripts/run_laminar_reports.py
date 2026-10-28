"""Run the laminar report for every penetration in metadata/penetrations.csv.

Usage: .venv\\Scripts\\python scripts\\run_laminar_reports.py [penetration ...] [--overwrite]
Writes <out_root>/laminar/<pen>_laminar.{png,npz,json} and laminar_summary.csv.
"""

import json
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

from v1lb.config import out_root
from v1lb.laminar_report import laminar_report

REPO = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    overwrite = "--overwrite" in sys.argv
    pens = pd.read_csv(REPO / "metadata" / "penetrations.csv", dtype=str).fillna("")
    if args:
        pens = pens[pens.penetration.isin(args)]
    out = out_root() / "laminar"
    rows = []
    for _, r in pens.iterrows():
        js = out / f"{r.penetration}_laminar.json"
        if js.exists() and not overwrite:
            rows.append(json.loads(js.read_text()))
            continue
        t = time.time()
        try:
            s = laminar_report(r.to_dict(), out, legacy_table=REPO / "docs" / "legacy_TuneList_all_versions.csv")
            rows.append(s)
            print(f"{r.penetration}: {time.time() - t:.0f}s  {s['orientation']} ToC {s['top_of_cortex']} "
                  f"tier {s['tier']} ({s['confidence']}, {s['complexity']}) L4/5 {s['l45']} steps {s['phase_steps']}", flush=True)
        except Exception:
            print(f"{r.penetration}: FAILED\n{traceback.format_exc()}", flush=True)
    pd.DataFrame(rows).to_csv(out / "laminar_summary.csv", index=False)
