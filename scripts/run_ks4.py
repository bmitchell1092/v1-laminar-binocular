"""Kilosort 4 on a penetration's bminteroc recording, read in place from the .ns6 (pilot started 2026-10-28).

Run in the Kilosort environment (CUDA PyTorch + kilosort 4.1, v1lb installed with --no-deps):
  .venv-ks4\\Scripts\\python scripts\\run_ks4.py <penetration> [...] [--th 9,8] [--tag name]
Writes <out_root>/kilosort4/<pen>[_<tag>]/ (phy-format results + run.json).

The .ns6 data packet is memory-mapped read-only as [samples, channels] int16 (all channels); the probe map picks
the penetration's bank: 32 contacts in one column, y = (contact - 1) x 100 um. Recordings with several data
packets (gaps of a few samples, 05 A4) are sorted as one stream (Kilosort's BinaryFileGroup); run.json keeps the
packet table so `v1lb.units` maps samples back to NSP ticks.
nblocks = 0 (no drift correction): Kilosort's docs advise skipping it for probes with <= 64 channels or
>= 50 um spacing, where drift estimates are unreliable. Other settings are Kilosort 4 defaults unless given.
"""

import json
import sys
import time

import numpy as np
import pandas as pd
import torch
from kilosort import DEFAULT_SETTINGS, run_kilosort
from kilosort.io import BinaryFileGroup

from v1lb.config import metadata_dir, out_root, recording_files
from v1lb.io.nsx import NsxFile


def probe_for(ns: NsxFile, bank: str) -> dict:
    rows = sorted((int(lab[2:]), i) for i, lab in enumerate(ns.labels) if lab.startswith(f"e{bank}") and lab[2:].isdigit())
    contacts = np.array([c for c, _ in rows])
    return dict(chanMap=np.array([i for _, i in rows]), xc=np.zeros(len(rows)), yc=(contacts - 1) * 100.0,
                kcoords=np.zeros(len(rows)), n_chan=len(rows))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    th = [float(x) for x in opts.get("th", "9,8").split(",")]
    tag = opts.get("tag", "")
    pens = pd.read_csv(metadata_dir() / "penetrations.csv", dtype=str).set_index("penetration")
    for pen in args:
        p = pens.loc[pen]
        f = recording_files(p.session, p.bminteroc)[".ns6"]
        ns = NsxFile(f)
        maps = [np.memmap(f, dtype="<i2", mode="r", offset=pk.offset, shape=(pk.n_samples, len(ns.labels)))
                for pk in ns.packets]
        data = maps[0] if len(maps) == 1 else BinaryFileGroup(file_objects=maps)
        probe = probe_for(ns, p.bank)
        settings = dict(DEFAULT_SETTINGS) | dict(n_chan_bin=len(ns.labels), fs=30000, nblocks=0,
                                                 Th_universal=th[0], Th_learned=th[1])
        res = out_root() / "kilosort4" / (pen + (f"_{tag}" if tag else ""))
        res.mkdir(parents=True, exist_ok=True)
        t = time.time()
        # `filename` is required but only recorded (data come from `file_object`; nothing is written next to it)
        run_kilosort(settings=settings, probe=probe, filename=f, file_object=data, results_dir=res,
                     data_dtype="int16", device=torch.device("cuda"))
        meta = dict(penetration=pen, recording=f"{p.session}_{p.bminteroc}", ns6=str(f),
                    packets=[[pk.timestamp, pk.n_samples] for pk in ns.packets], bank=p.bank, labels=[ns.labels[i] for i in probe["chanMap"]],
                    settings={k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
                              for k, v in settings.items()},
                    kilosort=__import__("kilosort").__version__, torch=torch.__version__,
                    minutes=round((time.time() - t) / 60, 1))
        (res / "run.json").write_text(json.dumps(meta, indent=1))
        n = len(pd.read_csv(res / "cluster_KSLabel.tsv", sep="\t"))
        print(f"{pen}: {meta['minutes']} min, {n} clusters -> {res}", flush=True)


if __name__ == "__main__":
    main()
