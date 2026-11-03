"""QC (docs/09 check J): is the 1 kHz LFP stream (.ns2) aligned with the 30 kHz stream (.ns6) when both are read by NSP tick?
Same channel, 20 s; the ns6 signal low-passed (a) with the amplifier's own LFP filter (causal 4th-order 500 Hz:
the lag should be 0 if the clocks agree) and (b) zero-phase (the lag is then the hardware filter's delay).

Usage: .venv\Scripts\python scripts\qc_nsx_clock.py   (reads raw files on D:, read-only; two sessions per monkey)
Result 2026-11-03: lag 0.0 ms with the causal filter (r = 1.000), 0.9 ms against zero-phase, in all four sessions.
"""
import glob, os
import numpy as np
from scipy import signal
from v1lb.io.nsx import NsxFile

files = []
for root in [r"D:/NEUROPHYS DATA/2021/Raw Data B", r"D:/NEUROPHYS DATA/2021/Raw Data J"]:
    fs = sorted(glob.glob(os.path.join(root, "**", "*bminteroc*.ns2"), recursive=True))
    files += [fs[0], fs[len(fs) // 2]]
for f2 in files:
    f6 = f2[:-1] + "6"
    if not os.path.exists(f6):
        continue
    n2, n6 = NsxFile(f2), NsxFile(f6)
    lab = [c.label for c in n2.channels if c.label.startswith("e")][10]
    start = max(n2.t_start, n6.t_start) + 30000 * 60
    stop = start + 30000 * 20
    x2, t2 = n2.read(start, stop, [lab])
    x6, t6 = n6.read(start - 30000, stop + 30000, [lab])
    x2, x6 = x2[:, 0].astype(float), x6[:, 0].astype(float)
    out = {}
    for name, filt in (("amplifier LFP filter, causal", "causal"), ("zero-phase 500 Hz", "zero")):
        sos = signal.butter(4, 500, "lowpass", fs=30000, output="sos")
        y = signal.sosfilt(sos, x6) if filt == "causal" else signal.sosfiltfilt(sos, x6)
        lags = np.arange(-150, 151)                        # ticks, +-5 ms
        r = []
        for L in lags:
            yi = np.interp(t2 - L, t6, y)                    # ns6 signal at (ns2 tick - L)
            r.append(np.corrcoef(yi - yi.mean(), x2 - x2.mean())[0, 1])
        r = np.array(r)
        out[name] = (lags[np.argmax(r)] / 30.0, r.max())
    print(os.path.basename(f2), lab, {k: (round(v[0], 3), round(v[1], 4)) for k, v in out.items()})
