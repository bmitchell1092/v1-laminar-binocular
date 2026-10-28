"""Compare the flash-evoked initial sink with laminar alignment v1.

Usage: .venv\\Scripts\\python scripts\\compare_flash_csd.py
Reads <out_root>/flash_csd (from run_flash_csd.py) and metadata/laminar_alignment_v1.csv.
Writes flash_csd_vs_alignment.csv (one row per penetration x evp) and flash_csd_grid.png
(closest evp per penetration, CSD on depth relative to L4/5).
Sink-bottom depth: depth (contacts, + = superficial) of the boundary just below the deepest
contact of the initial sink. 0 = it coincides with the L4/5 call; > 0 = sink bottom lies above it.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from v1lb.config import metadata_dir, out_root  # noqa: E402


def depth_of(c, l45, orientation):
    return (l45 - c) if orientation == "ascending" else (c - l45)


if __name__ == "__main__":
    src = out_root() / "flash_csd"
    summ = pd.read_csv(src / "flash_csd_summary.csv")
    al = pd.read_csv(metadata_dir() / "laminar_alignment_v1.csv").set_index("penetration")
    rows = []
    for _, s in summ.iterrows():
        a = al.loc[s.penetration]
        if pd.isna(s.sink_contact) or pd.isna(a.l45):
            rows.append(dict(s, include=bool(a.include), orientation=a.orientation, l45=a.l45, sink_bottom_depth=np.nan))
            continue
        deepest = s.sink_bottom_contact if a.orientation == "ascending" else s.sink_top_contact
        boundary = deepest + 0.5 if a.orientation == "ascending" else deepest - 0.5
        rows.append(dict(s, include=bool(a.include), orientation=a.orientation, l45=a.l45,
                         sink_depth=depth_of(s.sink_contact, a.l45, a.orientation),
                         sink_bottom_depth=depth_of(boundary, a.l45, a.orientation)))
    df = pd.DataFrame(rows)
    df["closest"] = df.groupby("penetration").minutes_from_bminteroc.transform(lambda m: m.abs() == m.abs().min())
    df.to_csv(src / "flash_csd_vs_alignment.csv", index=False)

    inc = df[df.include & df.closest].sort_values("penetration")
    n = len(inc)
    cols = 4
    fig, axes = plt.subplots(int(np.ceil(n / cols)), cols, figsize=(4.2 * cols, 4.6 * np.ceil(n / cols)), squeeze=False)
    for ax, (_, r) in zip(axes.ravel(), inc.iterrows()):
        d = np.load(src / f"{r.penetration}_{r.evp}.npz")
        a = al.loc[r.penetration]
        depth = depth_of(d["contacts"], a.l45, a.orientation)
        order = np.argsort(-depth)                      # superficial at top
        csd = d["csd"][:, order]
        lim = np.nanpercentile(np.abs(csd[(d["t_ms"] > 0) & (d["t_ms"] < 150)]), 99)
        ax.imshow(csd.T, aspect="auto", cmap="RdBu", vmin=-lim, vmax=lim,
                  extent=[d["t_ms"][0], d["t_ms"][-1] + 1, depth[order][-1] - 0.5, depth[order][0] + 0.5])
        ax.axhline(0, color="m", lw=2)
        if np.isfinite(r.sink_bottom_depth):
            ax.axhline(r.sink_bottom_depth, color="k", ls="--", lw=1.2)
            ax.plot(r.sink_onset_ms, r.sink_depth, "k*", ms=11)
        if a.top_of_cortex > 0:
            ax.axhline(depth_of(a.top_of_cortex, a.l45, a.orientation) + 0.5, color="g", lw=1.5)
        ax.set(xlim=(-20, 150), title=f"{r.penetration} {r.evp} ({r.minutes_from_bminteroc:+.0f} min)\n"
                                       f"sink bottom {r.sink_bottom_depth:+.1f} vs L4/5, onset {r.sink_onset_ms:.0f} ms",
               xlabel="ms from flash")
        ax.title.set_fontsize(9)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes[0, 0].set_ylabel("contacts from L4/5 call (+ superficial)")
    fig.suptitle("Flash-evoked CSD (red = sink) vs laminar alignment v1: magenta = L4/5 call, dashed = initial-sink "
                 "bottom, star = earliest sink, green = top of cortex", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(src / "flash_csd_grid.png", dpi=70)
    print(inc[["penetration", "evp", "minutes_from_bminteroc", "n_flashes", "sink_onset_ms", "sink_depth",
               "sink_bottom_depth"]].to_string(index=False))
    print("\nall evps of included penetrations (drift check):")
    print(df[df.include].pivot_table(index="penetration", columns="evp", values="sink_bottom_depth").round(1).to_string())
