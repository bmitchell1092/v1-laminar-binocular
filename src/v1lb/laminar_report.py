"""One-page laminar report per penetration: markers, rule-based calls, review flags.

Decision rules (Blake, 2026-10-28):
1. Orientation from SPI vs depth (template fit to Davis et al. 2023), cross-checked by the
   spectrolaminar peaks; disagreement is flagged.
2. Top of cortex = first active contact from the superficial end (spontaneous MUA above the
   noise floor, or evoked analog MUA).
3. L4/5 boundary = 2nd phase step below the top of cortex if there are two or more, else the
   only one. Extra deeper steps (white matter, a second fold of V1/V2) are expected on some days.
4. Spectrolaminar crossing and earliest CSD sink are reported as cross-checks.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import laminar as lam  # noqa: E402
from .bminteroc import presentation_table  # noqa: E402
from .config import out_root, recording_files  # noqa: E402
from .io import NsxFile  # noqa: E402
from .signals import analog_mua_epochs, probe_labels, read_lfp, threshold_spikes  # noqa: E402

TICKS = 30000


def legacy_limits(penetration: str, table: Path) -> dict[str, tuple[float, float, float]]:
    """Legacy [ToC, SinkBtm, BoC] by version, in raw contact numbers (flipped versions converted)."""
    if not table.exists():
        return {}
    df = pd.read_csv(table, dtype=str)
    df = df[df.Penetration == penetration]
    out = {}
    for _, r in df.iterrows():
        try:
            v = tuple(float(r[k]) for k in ("ToC", "SinkBtm", "BoC"))
        except (TypeError, ValueError):
            continue
        name = r.source.replace("ERR.mat in ", "").replace("git BM-ephys 1. preprocess/", "")
        if r.SortDirection == "descending":
            v = tuple(33 - x for x in v)
            name += " (flipped->raw)"
        out[f"{name} {r.source_date[:7]}"] = v
    return out


def laminar_report(row: dict, out_dir: Path | None = None, max_minutes: float = 30,
                   legacy_table: Path | None = None) -> dict:
    pen, session, bank, rec = row["penetration"], row["session"], row["bank"], row["bminteroc"]
    out_dir = Path(out_dir or out_root() / "laminar")
    out_dir.mkdir(parents=True, exist_ok=True)
    f = recording_files(session, rec)
    ns2, ns6 = NsxFile(f[".ns2"]), NsxFile(f[".ns6"])
    labels = probe_labels(ns2, bank)
    n = len(labels)
    start = ns6.t_start + 30 * TICKS
    stop = min(ns6.t_stop, start + int(max_minutes * 60 * TICKS))
    minutes = (stop - start) / TICKS / 60

    lfp, lt = read_lfp(ns2, labels, start, stop)
    phase = lam.generalized_phase(lam.bandpass(np.nan_to_num(lfp), 1000, (5, 50)), 1000, lp=0)
    spikes = threshold_spikes(ns6, labels, start, stop)
    rate = np.array([s.size for s in spikes]) / (minutes * 60)
    sp = lam.spike_phase(spikes, phase, lt)
    prof = lam.phase_profile(sp, min_spikes=500)
    T = lam.davis_templates()
    fit = lam.fit_template(prof, T["phase"], "phase", (T["n_phase"] >= 10).astype(float))
    freqs, relpow = lam.relative_power(lfp, 1000)
    ab, gm = lam.band_profile(freqs, relpow, (10, 19)), lam.band_profile(freqs, relpow, (75, 150))
    sc = lam.spectral_crossing(ab, gm)

    P, meta = presentation_table(session, rec)
    sel = P[P.completed & (P.pd_on_tick > 0) & (P.contrast_L + P.contrast_R >= 0.9)]
    ev = sel.pd_on_tick.to_numpy()[:300]
    env = analog_mua_epochs(ns6, labels, ev, pre_ms=100, post_ms=200)
    d = env[:, 140:220].mean(axis=1) - env[:, 20:100].mean(axis=1)        # 40-120 ms vs -80-0 ms
    evoked_t = d.mean(axis=0) / (d.std(axis=0, ddof=1) / np.sqrt(len(d)))
    ep = np.stack([ns2.read(t - 50 * 30, t + 150 * 30, labels)[0] for t in ev])
    ep = ep - ep[:, :50].mean(axis=1, keepdims=True)
    csd_avg = lam.csd(ep.mean(axis=0))
    t_ms = np.arange(-50, 150).astype(float)
    sink_onset, sink_first = lam.earliest_sink(csd_avg, t_ms)

    flags = []
    orientation = fit.orientation
    if sc.orientation != fit.orientation:
        flags.append(f"orientation disagreement: SPI {fit.orientation} vs spectrolaminar {sc.orientation}")
    if fit.margin < 0.1:
        flags.append(f"weak SPI orientation margin ({fit.margin:.2f})")
    if (sp.n >= 500).sum() < 8:
        flags.append(f"few spike channels with >=500 spikes ({int((sp.n >= 500).sum())})")
    toc, top_in, toc_by = lam.top_of_cortex(rate, evoked_t, orientation)
    if toc < 0:
        flags.append("top of cortex not found")
    if top_in:
        flags.append("top contact already in cortex: the supra->input step may lie above the probe")
    order = lam.superficial_order(n, orientation)
    active_idx = [int(i) for i in order if rate[i] > 2.5 or evoked_t[i] > 6.0]
    deepest_active = (active_idx[-1] + 1) if active_idx else None
    steps = lam.phase_steps(prof, orientation, start_contact=max(toc, 1), stop_contact=deepest_active)
    l45 = lam.l45_from_steps(steps)
    if not np.isfinite(l45):
        flags.append("no phase step found below the top of cortex")
    # evidence that the probe crossed a second fold of cortex (V1 or V2) below the first column
    reversed_steps = [s for s in steps if lam.step_kind(s[2]) == "reversed" and abs(s[2]) >= 0.8]
    gaps = [(float(a[0]), float(b[0])) for a, b in zip(steps[:-1], steps[1:]) if lam.mua_gap(rate, a[0], b[0])]
    second_fold = bool(reversed_steps or gaps)
    if reversed_steps:
        flags.append(f"reversed phase step(s) at {[float(s[0]) for s in reversed_steps]}: possible second fold")
    if gaps:
        flags.append(f"MUA gap between steps {gaps}: possible white matter / second fold")
    sink = sink_first + 1 if sink_first >= 0 else np.nan
    agree_spec = bool(np.isfinite(l45) and np.isfinite(sc.crossing) and abs(sc.crossing - l45) <= 2)
    agree_sink = bool(np.isfinite(l45) and np.isfinite(sink) and abs(sink - l45) <= 3)
    if np.isfinite(l45) and not agree_spec:
        flags.append(f"spectrolaminar crossing ({sc.crossing:.1f}) >2 contacts from L4/5 call ({l45:.1f})")
    if np.isfinite(l45) and not agree_sink:
        flags.append(f"earliest sink ({sink}) >3 contacts from L4/5 call ({l45:.1f})")

    # complexity and confidence, for ranking penetrations (Blake, 2026-10-28)
    checks = {
        "orientation markers agree": sc.orientation == fit.orientation,
        "SPI orientation margin >= 0.1": fit.margin >= 0.1,
        "top of cortex on the probe": toc > 0 and not top_in,
        "L4/5 step found": bool(np.isfinite(l45)),
        "spectral crossing within 2": agree_spec,
        "earliest sink within 3": agree_sink,
        ">= 16 spike channels with >= 500 spikes": int((sp.n >= 500).sum()) >= 16,
        "no second-fold evidence": not second_fold,
    }
    confidence = sum(checks.values()) / len(checks)
    essential = checks["orientation markers agree"] and checks["L4/5 step found"]
    tier = "A" if essential and confidence >= 0.85 else "B" if essential and confidence >= 0.6 else "C"
    complexity = ("possible second fold" if second_fold else
                  "top contact in cortex" if top_in else "single column")

    summary = dict(penetration=pen, recording=rec, minutes=round(minutes, 1), n_contacts=n,
                   tier=tier, confidence=round(confidence, 2), complexity=complexity,
                   orientation=orientation, spi_fit_margin=round(fit.margin, 3), spi_fit_score=round(fit.score, 3),
                   spectral_orientation=sc.orientation, spectral_crossing=round(sc.crossing, 2),
                   gamma_peak=sc.gamma_peak, alpha_peak=sc.alpha_peak,
                   top_of_cortex=toc, top_of_cortex_by=toc_by, top_contact_in_cortex=top_in,
                   phase_steps=[(round(s, 1), round(m, 2), round(g, 2)) for s, m, g in steps], l45=l45,
                   earliest_sink=None if np.isnan(sink) else int(sink),
                   earliest_sink_ms=float(sink_onset[sink_first]) if sink_first >= 0 else None,
                   n_csd_trials=len(ev), checks=checks, flags=flags)
    legacy = legacy_limits(pen, legacy_table) if legacy_table else {}
    _figure(pen, labels, sp, prof, rate, evoked_t, ab, gm, csd_avg, t_ms, summary, legacy, out_dir / f"{pen}_laminar.png")
    np.savez_compressed(out_dir / f"{pen}_laminar.npz", spi=sp.spi, theta=sp.theta, ppc=sp.ppc, n=sp.n,
                        rate=rate, evoked_t=evoked_t, profile=prof, alpha_beta=ab, gamma=gm, csd=csd_avg,
                        labels=np.array(labels))
    (out_dir / f"{pen}_laminar.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


def _figure(pen, labels, sp, prof, rate, evoked_t, ab, gm, csd_avg, t_ms, s, legacy, path):
    n = len(labels)
    ch = np.arange(1, n + 1)
    ext = [0.5, n + 0.5, n + 0.5, 0.5]
    fig, ax = plt.subplots(1, 7, figsize=(30, 7), gridspec_kw=dict(width_ratios=[1.2, 1.2, 0.8, 0.6, 0.6, 0.6, 1.2]))
    im = ax[0].imshow(sp.theta, cmap="twilight", vmin=-np.pi, vmax=np.pi, extent=ext)
    plt.colorbar(im, ax=ax[0], fraction=0.046)
    ax[0].set(title=f"{pen}\npreferred phase (rows: spike ch, cols: LFP ch)", xlabel="LFP contact", ylabel="spike contact")
    im = ax[1].imshow(sp.spi, cmap="hot", extent=ext)
    plt.colorbar(im, ax=ax[1], fraction=0.046)
    ax[1].set(title="SPI", xlabel="LFP contact")
    ax[2].plot(prof, ch, "o-")
    ax[2].set(xlim=(-np.pi, np.pi), title="preferred phase vs LFP contact\n(Davis: supra/input ±π, deep ≈ −1.2)")
    ax[2].axvline(-1.2, color="0.6", ls=":")
    ax[3].barh(ch, rate)
    ax[3].axvline(2.5, color="0.5", ls=":")
    ax[3].set(title="MUA rate (Hz)")
    ax[4].barh(ch, np.clip(evoked_t, -5, 40))
    ax[4].axvline(6, color="0.5", ls=":")
    ax[4].set(title="evoked MUA (t)")
    ax[5].plot(ab, ch, label="10-19 Hz")
    ax[5].plot(gm, ch, label="75-150 Hz")
    ax[5].legend(fontsize=7)
    ax[5].set(title="relative power")
    lim = np.nanmax(np.abs(csd_avg)) * 0.8
    ax[6].imshow(csd_avg.T, aspect="auto", cmap="RdBu", vmin=-lim, vmax=lim, extent=[t_ms[0], t_ms[-1] + 1, n + 0.5, 0.5])
    ax[6].set(title="CSD (red = sink), high-contrast onsets", xlabel="ms from photodiode")
    if s["earliest_sink"]:
        ax[6].plot(s["earliest_sink_ms"], s["earliest_sink"], "k*", ms=14)
    for a in ax[2:]:
        a.set_ylim(n + 0.5, 0.5)
        if s["top_of_cortex"] > 0:
            a.axhline(s["top_of_cortex"], color="g", lw=2)
        for b, _, g in s["phase_steps"]:
            a.axhline(b, color={"canonical": "m", "reversed": "c"}.get(lam.step_kind(g), "0.4"), lw=0.8, ls="--")
        if np.isfinite(s["l45"]):
            a.axhline(s["l45"], color="m", lw=2.5)
        if np.isfinite(s["spectral_crossing"]):
            a.axhline(s["spectral_crossing"], color="orange", lw=1.2, ls="-.")
        for k, (name, (toc, snk, boc)) in enumerate(legacy.items()):
            a.axhline(snk, color=f"C{k + 4}", lw=1, ls=":")
    txt = (f"TIER {s['tier']} (confidence {s['confidence']:.2f}; {s['complexity']})   "
           f"orientation: {s['orientation']} (SPI margin {s['spi_fit_margin']:+.2f}; spectral {s['spectral_orientation']})   "
           f"top of cortex (green): {s['top_of_cortex']} [{s['top_of_cortex_by']}]   "
           f"steps (pos, |rad|, signed; cyan = reversed): {s['phase_steps']}   L4/5 (magenta): {s['l45']}   "
           f"spectral crossing (orange): {s['spectral_crossing']:.1f}   earliest sink: {s['earliest_sink']} "
           f"({s['earliest_sink_ms']} ms)\nlegacy SinkBtm (dotted): "
           + "; ".join(f"C{k + 4} {n_}: {v[1]:g}" for k, (n_, v) in enumerate(legacy.items()))
           + ("\nFLAGS: " + " | ".join(s["flags"]) if s["flags"] else ""))
    fig.text(0.01, 0.005, txt, fontsize=9, va="bottom", color="darkred" if s["flags"] else "black")
    plt.tight_layout(rect=(0, 0.08, 1, 1))
    plt.savefig(path, dpi=60)
    plt.close(fig)
