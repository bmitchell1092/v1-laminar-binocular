"""Register each evp recording to its penetration's bminteroc recording and move the flash
initial-sink bottom into bminteroc contact numbers.

Usage: .venv\\Scripts\\python scripts\\register_flash_csd.py [--all]
Needs <out_root>/flash_csd/flash_csd_summary.csv (run_flash_csd.py). Default: included penetrations only.
Writes <out_root>/flash_csd/flash_csd_registered.csv.
"""

import sys
import time

import numpy as np
import pandas as pd

from v1lb.config import metadata_dir, out_root
from v1lb.registration import depth_profiles, register

if __name__ == "__main__":
    src = out_root() / "flash_csd"
    summ = pd.read_csv(src / "flash_csd_summary.csv")
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v1.csv").set_index("penetration")
    if "--all" not in sys.argv:
        summ = summ[summ.penetration.isin(al.index[al.include])]
    rows = []
    for pen, grp in summ.groupby("penetration"):
        p, a = pens.loc[pen], al.loc[pen]
        t = time.time()
        ref = depth_profiles(p.session, p.bminteroc, p.bank)
        for _, s in grp.iterrows():
            reg = register(depth_profiles(p.session, s.evp, p.bank), ref)
            row = dict(s, shift=reg.shift, reg_score=round(reg.score, 3), reg_margin=round(reg.margin, 3),
                       orientation=a.orientation, l45=a.l45)
            if pd.notna(s.sink_contact):
                deepest = s.sink_bottom_contact if a.orientation == "ascending" else s.sink_top_contact
                c = deepest + reg.shift
                boundary = c + 0.5 if a.orientation == "ascending" else c - 0.5
                row["registered_sink_bottom_contact"] = c
                row["registered_sink_bottom_depth"] = (a.l45 - boundary) if a.orientation == "ascending" else (boundary - a.l45)
            rows.append(row)
            print(f"{pen} {s.evp}: shift {reg.shift:+d} (r {reg.score:.2f}, margin {reg.margin:.2f}) -> "
                  f"sink bottom depth {row.get('registered_sink_bottom_depth', np.nan):+.1f} vs L4/5", flush=True)
        print(f"  ({time.time() - t:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(src / "flash_csd_registered.csv", index=False)
