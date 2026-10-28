"""Flash-evoked CSD for every evp recording of every penetration, in contact-number order.

Usage: .venv\\Scripts\\python scripts\\run_flash_csd.py
Writes <out_root>/flash_csd/<pen>_<evp>.npz and flash_csd_summary.csv. Contacts are in
label order (1..n); depth relative to L4/5 is applied later from the alignment table.
"""

import time
import traceback

import numpy as np
import pandas as pd

from v1lb.config import metadata_dir, out_root, recording_files
from v1lb.evp import flash_csd
from v1lb.io import NevFile

if __name__ == "__main__":
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).fillna("")
    out = out_root() / "flash_csd"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in pens.iterrows():
        b0 = NevFile(recording_files(r.session, r.bminteroc)[".nev"]).origin_utc
        for evp in [e for e in r.evp.split(";") if e]:
            t = time.time()
            try:
                f = flash_csd(r.session, evp, r.bank, "ascending")
                dt_min = (NevFile(recording_files(r.session, evp)[".nev"]).origin_utc - b0).total_seconds() / 60
                np.savez_compressed(out / f"{r.penetration}_{evp}.npz", t_ms=f.t_ms, contacts=f.contacts,
                                    lfp=f.lfp, csd=f.csd, mua=f.mua)
                ok = f.sink_index >= 0
                rows.append(dict(penetration=r.penetration, evp=evp, minutes_from_bminteroc=round(dt_min),
                                 n_flashes=f.n_flashes,
                                 sink_contact=int(f.contacts[f.sink_index]) if ok else None,
                                 sink_onset_ms=f.sink_onset_ms,
                                 sink_top_contact=int(f.contacts[f.sink_top]) if ok else None,
                                 sink_bottom_contact=int(f.contacts[f.sink_bottom]) if ok else None))
                print(f"{r.penetration} {evp}: {time.time() - t:.0f}s {rows[-1]}", flush=True)
            except Exception:
                print(f"{r.penetration} {evp}: FAILED\n{traceback.format_exc()}", flush=True)
    pd.DataFrame(rows).to_csv(out / "flash_csd_summary.csv", index=False)
