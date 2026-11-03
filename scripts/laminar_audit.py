"""Laminar audit (docs/09 W): one report card per penetration from laminar markers only, blind to binocular outcomes.

Usage: .venv\\Scripts\\python scripts\\laminar_audit.py [cache|cards|all]
cache: monocular grating CSD and MUA per contact (clean presentations, one eye at 45 or 90%, the other eye 0; contra
       and ipsi kept apart), -50 to 150 ms from stimulus onset -> <out_root>/laminar_audit/{pen}_mono.npz.
cards: the markers from output/laminar/{pen}_laminar.npz (preferred phase, SPI and PPC, spontaneous rate,
       relative power) and the monocular cache; the automated calls (orientation, top of cortex, phase reversal,
       input layer, bottom of the first sheet, second sheet) -> one PNG card per penetration and
       docs/qc/laminar_audit.csv. The calls are drafts for Blake's review (docs/09 W).
"""

import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
import matplotlib.ticker  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from v1lb import laminar as lam  # noqa: E402
from v1lb.config import out_root  # noqa: E402
from v1lb.dataset import clean_presentations, open_penetration  # noqa: E402

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
QC = REPO / "docs" / "qc"
OUT = out_root() / "laminar_audit"
WIN = (-50, 150)
RATE_FLOOR = 2.5          # spontaneous threshold crossings (Hz) above which a contact is in grey matter
MIN_SPIKES = 500          # spike channels used for SPI and the phase profile


def pens():
    A = pd.read_csv(REPO / "metadata" / "laminar_alignment_v2.csv")
    return A[A.depth_usable.astype(bool)].penetration.tolist()


# ------------------------------------------------------------------------------------------------ cache
def cache_one(pen):
    d = out_root() / "datasets"
    ds = open_penetration(pen)
    P = pd.read_parquet(d / f"{pen}_presentations.parquet")
    E = pd.read_parquet(d / f"{pen}_presentation_eye.parquet")
    keep = clean_presentations(ds, P, E)
    t = ds.time.values
    w = (t >= WIN[0]) & (t < WIN[1])
    cc, ci = P.contrast_contra.to_numpy(), P.contrast_ipsi.to_numpy()
    sel = {"contra": keep & (cc >= 0.45) & (ci == 0), "ipsi": keep & (ci >= 0.45) & (cc == 0)}
    mua, csd = ds.mua.values[:, w].astype(float), ds.csd.values[:, w].astype(float)
    tt = t[w]
    base = tt < 0
    out = dict(t=tt, contact=ds.contact.values, depth_below_toc=ds.depth_below_toc.values,
               in_cortex=ds.in_cortex.values, bad_lfp=ds.bad_lfp.values)
    for eye, s in sel.items():
        m = mua[s] - mua[s][:, base].mean(axis=1, keepdims=True)
        c = csd[s] - csd[s][:, base].mean(axis=1, keepdims=True)
        ev = m[:, (tt >= 40) & (tt < 120)].mean(axis=1)
        out[f"mua_{eye}"] = m.mean(axis=0)
        out[f"csd_{eye}"] = c.mean(axis=0)
        out[f"evoked_t_{eye}"] = ev.mean(axis=0) / (ev.std(axis=0, ddof=1) / np.sqrt(len(ev)))
        out[f"n_{eye}"] = int(s.sum())
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{pen}_mono.npz", **out)
    return pen, out["n_contra"], out["n_ipsi"]


# ---------------------------------------------------------------------------------------------- markers
def half_max_onset(x, t, lo=10, hi=120):
    w = (t >= lo) & (t < hi)
    pk = np.nanmax(x[w])
    j = np.flatnonzero(w & (x >= 0.5 * pk))[0]
    if not w[j - 1]:
        return float(t[j])
    return float(t[j - 1] + (0.5 * pk - x[j - 1]) / (x[j] - x[j - 1]))


def run_median(x, k=3):
    return pd.Series(x).rolling(k, center=True, min_periods=2).median().to_numpy()


def first_run(mask, start, length):
    """First index >= start that begins `length` consecutive True values, or None."""
    for i in range(max(int(start), 0), len(mask) - length + 1):
        if mask[i:i + length].all():
            return i
    return None


def sink_onsets(csd, t, ok, window=(20, 100), min_rel=0.4, frac=0.25):
    """Per contact, onset of the evoked sink: from the sink's peak in `window`, back to the last sample above
    `frac` x peak (interpolated), so noise before the sink cannot set it. Contacts count when `ok` and their peak
    reaches `min_rel` of the strongest sink among `ok` contacts. csd is [time, contact], sinks negative."""
    w = np.flatnonzero((t >= window[0]) & (t <= window[1]))
    X = np.where(np.isfinite(csd[w]), csd[w], np.inf)
    ok = ok & np.isfinite(csd[w]).all(axis=0)
    pk_i = w[np.argmin(X, axis=0)]
    pk = csd[pk_i, np.arange(csd.shape[1])]
    on = np.full(csd.shape[1], np.nan)
    if not ok.any() or np.nanmin(pk[ok]) >= 0:
        return on
    strongest = np.nanmin(pk[ok])
    for c in np.flatnonzero(ok & (pk <= min_rel * strongest)):
        j = pk_i[c]
        while j > w[0] and csd[j - 1, c] <= frac * pk[c]:
            j -= 1
        x0, x1 = csd[j - 1, c], csd[j, c]
        on[c] = t[j - 1] + (frac * pk[c] - x0) / (x1 - x0) if j > w[0] else t[j]
    return on


def spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 6:
        return np.nan
    return float(stats.spearmanr(x[ok], y[ok]).statistic)


def template_curve(n, fo, orientation):
    """The fitted Davis phase template at each contact (NaN off the template), contact order."""
    T = lam.davis_templates()
    tmpl, wt = T["phase"], (T["n_phase"] >= 10).astype(float)
    good = np.isfinite(tmpl) & (wt > 0)
    lo, hi = np.flatnonzero(good)[[0, -1]]
    core = np.where(good[lo:hi + 1], tmpl[lo:hi + 1], np.nan)
    ext = np.r_[np.full(n, core[0]), core, np.full(n, core[-1])]
    c = np.arange(n)
    k = fo.shift + c if orientation == "ascending" else fo.shift - c
    return np.where((k >= 0) & (k < ext.size), ext[np.clip(k, 0, ext.size - 1)], np.nan)


def markers(pen):
    L = np.load(out_root() / "laminar" / f"{pen}_laminar.npz")
    J = json.loads((out_root() / "laminar" / f"{pen}_laminar.json").read_text())
    M = np.load(OUT / f"{pen}_mono.npz")
    n = len(L["rate"])
    contact = np.arange(1, n + 1)
    T = lam.davis_templates()
    wt = (T["n_phase"] >= 10).astype(float)
    good_sp = L["n"] >= MIN_SPIKES
    spi_c = np.where(good_sp, np.nanmean(L["spi"], axis=1), np.nan)
    ppc_c = np.where(good_sp, np.nanmean(L["ppc"], axis=1), np.nan)
    ab_c, gm_c = L["alpha_beta"].astype(float), L["gamma"].astype(float)

    # 1. orientation: phase template, spectrolaminar peaks, SPI falling with contact order
    fit = lam.fit_template(L["profile"], T["phase"], "phase", wt)
    sc = lam.spectral_crossing(ab_c, gm_c)
    act = L["rate"] > RATE_FLOOR
    rho_c = spearman(contact[act].astype(float), spi_c[act])
    votes = {"phase": fit.orientation, "spectral": sc.orientation,
             "spi": "ascending" if rho_c < 0 else "descending"}
    asc = sum(v == "ascending" for v in votes.values())
    orientation = "ascending" if asc >= 2 else "descending"
    order = lam.superficial_order(n, orientation)
    to_pos = (lambda b: b - 1) if orientation == "ascending" else (lambda b: n - b)   # contact number -> position
    to_contact = lambda p: int(order[int(p)]) + 1                                    # noqa: E731

    # everything below in probe position, 0 = superficial end, 0.1 mm per step
    pos = np.arange(n)
    rate, spi, ppc, ab, gm = (x[order] for x in (L["rate"], spi_c, ppc_c, ab_c, gm_c))
    prof, good_sp = L["profile"][order], good_sp[order]
    t = M["t"]
    raw = np.argsort(M["contact"])          # datasets store contacts superficial -> deep (v2 orientation)
    mua = (M["mua_contra"] + M["mua_ipsi"])[:, raw][:, order]
    csd = ((M["csd_contra"] + M["csd_ipsi"]) / 2)[:, raw][:, order]
    ev_t = np.maximum(M["evoked_t_contra"], M["evoked_t_ipsi"])[raw][order]
    amp = run_median(np.nanmax(mua[(t >= 30) & (t < 120)], axis=0))
    amp_n = amp / np.nanpercentile(amp, 90)

    # 2. top of cortex: spontaneous rate (existing rule) and evoked MUA (>= 20% of the probe's 90th percentile, t >= 6)
    toc, top_in, _ = lam.top_of_cortex(L["rate"], L["evoked_t"], orientation)
    p_top_rate = to_pos(toc) if toc > 0 else np.nan
    p_top = first_run((amp_n >= 0.2) & (ev_t >= 6), 0, 3)
    p_top = 0 if p_top is None else p_top

    # 3. phase reversal: Davis template boundary with the orientation fixed; its step across 3 contacts each side
    fo = lam.fit_template(L["profile"], T["phase"], "phase", wt, orientations=(orientation,))
    p_rev = to_pos(fo.boundary)
    on_probe = 1 <= p_rev <= n - 2
    step = np.nan
    if on_probe:
        i = int(np.floor(p_rev))
        a, b = np.exp(1j * prof[max(i - 2, 0):i + 1]), np.exp(1j * prof[i + 1:i + 4])
        step = float(np.angle(np.nanmean(b) * np.conj(np.nanmean(a))))
    clear = bool(on_probe and abs(step) >= 1.0)
    steps = [(to_pos(s[0]), s[1], s[2]) for s in lam.phase_steps(L["profile"], orientation)]

    # 4. bottom of the first sheet: below the reversal (or 0.8 mm below the top), first 2 contacts under 35%
    start = int(np.ceil(p_rev)) if clear else p_top + 8
    quiet = amp_n < 0.35
    g = first_run(quiet, start, 2)
    bottom, gap_end = (g - 1, None) if g is not None else (n - 1, None)
    if g is not None:
        e = g
        while e < n and quiet[e]:
            e += 1
        gap_end = e if e < n else None
    second_active = gap_end is not None and first_run(amp_n >= 0.5, gap_end, 2) is not None
    rev_below = [s for s in steps if lam.step_kind(s[2]) == "reversed" and abs(s[2]) >= 0.8 and s[0] > p_rev + 1]

    # 5. input layer: earliest MUA onset and earliest sink in the first sheet; a second early zone below it
    driven = (ev_t > 6) & (amp_n >= 0.2) & (pos >= p_top)
    on_mua = np.full(n, np.nan)
    for p in np.flatnonzero(driven):
        on_mua[p] = half_max_onset(mua[:, p], t)
    om = np.where(driven, run_median(on_mua), np.nan)
    first = (pos >= p_top) & (pos <= bottom)
    p_mua = int(np.nanargmin(np.where(first, om, np.nan))) if np.isfinite(om[first]).any() else np.nan
    deep = pos > bottom + 2
    p_mua2 = int(np.nanargmin(np.where(deep, om, np.nan))) if np.isfinite(om[deep]).any() else np.nan
    early2 = bool(np.isfinite(p_mua2) and np.isfinite(p_mua) and om[p_mua2] <= om[p_mua] + 3
                  and np.nanmax(om[p_mua:p_mua2]) >= om[p_mua] + 5)      # separated by a ridge of later onsets
    s_on = sink_onsets(csd, t, pos >= p_top)
    s_first = np.where(first, s_on, np.nan)
    p_sink = int(np.nanargmin(s_first)) if np.isfinite(s_first).any() else np.nan
    p_sink_all = int(np.nanargmin(s_on)) if np.isfinite(s_on).any() else np.nan

    # 6. Blake's checks within the first sheet
    rho_spi = spearman(pos[first & good_sp].astype(float), spi[first & good_sp])
    rho_ppc = spearman(pos[first & good_sp].astype(float), ppc[first & good_sp])
    rho_ab = spearman(pos[first].astype(float), ab[first] - gm[first])
    p_cross = to_pos(sc.crossing) if np.isfinite(sc.crossing) else np.nan

    rel = lambda p: round((p - p_rev) / 10, 2) if np.isfinite(p) else np.nan       # noqa: E731
    row = dict(
        penetration=pen, orientation=orientation, votes="/".join(v[:4] for v in votes.values()),
        n_votes_agree=max(asc, 3 - asc), legacy_orientation=J["orientation"],
        top_rate=toc, top_mua=to_contact(p_top), top_rate_minus_mua_mm=round((p_top_rate - p_top) / 10, 1),
        reversal=round(float(fo.boundary), 1), reversal_on_probe=on_probe, reversal_step_rad=round(step, 2),
        reversal_clear=clear, reversal_fit=round(fo.score, 3),
        top_depth_mm=rel(p_top), bottom=to_contact(bottom), bottom_depth_mm=rel(bottom + 0.5),
        first_sheet_mm=round((bottom - p_top + 1) / 10, 1), white_matter_dip=g is not None,
        second_active_zone=second_active, reversed_steps_below=len(rev_below),
        sink=to_contact(p_sink) if np.isfinite(p_sink) else np.nan, sink_depth_mm=rel(p_sink),
        sink_onset_ms=round(float(s_on[p_sink]), 1) if np.isfinite(p_sink) else np.nan,
        sink_whole_probe_depth_mm=rel(p_sink_all),
        mua=to_contact(p_mua) if np.isfinite(p_mua) else np.nan, mua_depth_mm=rel(p_mua),
        mua_onset_ms=round(float(om[p_mua]), 1) if np.isfinite(p_mua) else np.nan,
        sink_minus_mua_mm=round((p_sink - p_mua) / 10, 1) if np.isfinite(p_sink) and np.isfinite(p_mua) else np.nan,
        second_early_zone=early2, second_zone_depth_mm=rel(p_mua2) if early2 else np.nan,
        second_zone_onset_ms=round(float(om[p_mua2]), 1) if early2 else np.nan,
        crossing_depth_mm=rel(p_cross), spi_rho=round(rho_spi, 2), ppc_rho=round(rho_ppc, 2),
        alpha_beta_rho=round(rho_ab, 2), legacy_l45=J["l45"], legacy_sink=J["earliest_sink"])
    arrays = dict(n=n, order=order, prof=prof, tmpl=template_curve(n, fo, orientation)[order], steps=steps,
                  spi=spi, rate=rate, ab=ab, gm=gm, t=t, mua=mua, csd=csd, s_on=s_on, on_mua=on_mua, om=om,
                  amp_n=amp_n, p_top=p_top, p_top_rate=p_top_rate, p_rev=p_rev, bottom=bottom, p_sink=p_sink,
                  p_mua=p_mua, p_mua2=p_mua2 if early2 else np.nan, p_cross=p_cross)
    return row, arrays


# ------------------------------------------------------------------------------------------------- card
INK, GREY, MAG, BLUE, ORANGE = "#0b0b0b", "#8a8985", "#d6336c", "#2a78d6", "#eb6834"


def card(row, A, path):
    n = A["n"]
    y = np.arange(n)
    hl = [(A["p_top"] - 0.5, INK, "--"), (A["p_top_rate"] - 0.5, GREY, ":"), (A["p_rev"], MAG, "-"),
          (A["bottom"] + 0.5, INK, "--"), (A["p_sink"], BLUE, ":"), (A["p_mua"], ORANGE, ":"),
          (A["p_mua2"], ORANGE, "-."), (A["p_cross"], "#1baf7a", ":")]
    fig, axs = plt.subplots(1, 6, figsize=(16, 6.6), sharey=True)
    fig.subplots_adjust(left=0.045, right=0.965, top=0.76, bottom=0.08, wspace=0.14)
    ax = axs[0]
    ax.plot(A["tmpl"], y, color=GREY, lw=3, alpha=0.5, label="Davis template")
    ax.plot(A["prof"], y, "o-", ms=3, color=INK, lw=1)
    for p, _, s in A["steps"]:
        k = lam.step_kind(s)
        ax.plot([np.pi * 0.92], [p], marker="<", ms=5, color=MAG if k == "canonical" else BLUE if k == "reversed"
                else GREY)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks([-np.pi, 0, np.pi], ["-π", "0", "π"])
    ax.set_title("preferred phase (grey: Davis template)\n◄ steps: magenta canonical, blue reversed", fontsize=7)
    ax.set_ylabel("position from superficial end (contacts, 0.1 mm)")
    ax = axs[1]
    ax.plot(A["spi"], y, "o-", ms=3, color=INK, lw=1)
    a2 = ax.twiny()
    a2.semilogx(np.maximum(A["rate"], 0.3), y, color=ORANGE, lw=0.9)
    a2.axvline(RATE_FLOOR, color=ORANGE, lw=0.6, ls=":")
    a2.set_xticks([1, 3, 10, 30], ["1", "3", "10", "30"])
    a2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a2.tick_params(labelsize=6, colors=ORANGE)
    ax.set_title(f"SPI (black); spontaneous rate, Hz (orange, log)\nρ(SPI, depth) in first sheet {row['spi_rho']}",
                 fontsize=7, pad=16)
    ax = axs[2]
    ax.plot(A["ab"], y, color=BLUE)
    ax.plot(A["gm"], y, color=ORANGE)
    ax.set_title(f"power 10–19 Hz (blue), 75–150 Hz (orange)\nρ(α/β − γ, depth) in first sheet {row['alpha_beta_rho']}",
                 fontsize=7)
    t = A["t"]
    w = (t >= 0) & (t <= 60)
    ax = axs[3]
    C = A["csd"][w].T
    v = np.nanpercentile(np.abs(C), 98)
    ax.imshow(C, aspect="auto", cmap="RdBu", vmin=-v, vmax=v, extent=(0, 60, n - 0.5, -0.5), interpolation="bilinear")
    ax.plot(A["s_on"], y, "k.", ms=4)
    ax.set_title(f"monocular CSD (sink blue; · onsets)\nearliest sink {row['sink_onset_ms']} ms", fontsize=7)
    ax.set_xlabel("ms from stimulus onset", fontsize=7)
    ax = axs[4]
    X = A["mua"]
    pk = np.nanmax(X[(t >= 10) & (t < 120)], axis=0)
    ax.imshow((X / np.where(pk > 0, pk, np.nan))[w].T, aspect="auto", cmap="Greys", vmin=0, vmax=1,
              extent=(0, 60, n - 0.5, -0.5), interpolation="bilinear")
    ax.plot(A["on_mua"], y, "o", ms=2.5, color=ORANGE)
    ax.plot(A["om"], y, color=ORANGE, lw=1)
    ax.set_title(f"monocular MUA (own peak = 1; · onsets)\nearliest onset {row['mua_onset_ms']} ms", fontsize=7)
    ax.set_xlabel("ms from stimulus onset", fontsize=7)
    ax = axs[5]
    ax.plot(A["amp_n"], y, "o-", ms=3, color=INK, lw=1)
    for x in (0.2, 0.35, 0.5):
        ax.axvline(x, color=GREY, lw=0.5, ls=":")
    ax.set_xlim(0, 1.6)
    ax.set_title("evoked MUA amplitude\n(probe's 90th percentile = 1)", fontsize=7.5)
    ax.yaxis.set_label_position("right")
    rax = ax.secondary_yaxis("right")
    rax.set_yticks(y[::2], [str(int(A["order"][p]) + 1) for p in y[::2]], fontsize=6)
    rax.set_ylabel("contact", fontsize=7)
    for ax in axs:
        for p, col, ls in hl:
            if np.isfinite(p):
                ax.axhline(p, color=col, lw=0.9, ls=ls)
        ax.tick_params(labelsize=7)
    axs[0].set_ylim(n - 0.5, -0.5)
    flags = []
    if row["n_votes_agree"] < 3:
        flags.append(f"orientation 2-of-3 ({row['votes']})")
    if row["orientation"] != row["legacy_orientation"]:
        flags.append("orientation differs from legacy")
    if not row["reversal_clear"]:
        flags.append("no clear phase reversal on the probe")
    if abs(row["top_rate_minus_mua_mm"]) > 0.3:
        flags.append(f"tops disagree by {row['top_rate_minus_mua_mm']} mm")
    if row["first_sheet_mm"] > 2.0:
        flags.append(f"first sheet {row['first_sheet_mm']} mm (> 2.0)")
    if row["second_active_zone"] or row["reversed_steps_below"] or row["second_early_zone"]:
        flags.append("second sheet evidence: " + ", ".join(k for k, v in (
            ("active zone below a dip", row["second_active_zone"]),
            ("reversed phase step", row["reversed_steps_below"]),
            (f"early MUA zone {row['second_zone_onset_ms']} ms", row["second_early_zone"])) if v))
    if np.isfinite(row["sink_minus_mua_mm"]) and abs(row["sink_minus_mua_mm"]) > 0.2:
        flags.append(f"sink and MUA {row['sink_minus_mua_mm']:+.1f} mm apart")
    head = (f"{row['penetration']} — orientation {row['orientation']} (votes phase/spectral/SPI {row['votes']}; "
            f"legacy {row['legacy_orientation']});  depth from the phase reversal (contact {row['reversal']}, "
            f"step {row['reversal_step_rad']} rad):  top {row['top_depth_mm']} mm,  earliest sink "
            f"{row['sink_depth_mm']} mm,  earliest MUA {row['mua_depth_mm']} mm,  α/β–γ crossing "
            f"{row['crossing_depth_mm']} mm,  bottom of first sheet {row['bottom_depth_mm']} mm")
    key = ("lines: top (evoked MUA) black dashed · top (spontaneous rate) grey dotted · phase reversal magenta · "
           "bottom of first sheet black dashed · earliest sink blue dotted · earliest MUA orange dotted · "
           "second early zone orange dash-dot · α/β–γ crossing green dotted")
    fig.text(0.01, 0.975, head, fontsize=8, va="top", wrap=True)
    fig.text(0.01, 0.925, "flags: " + ("; ".join(flags) if flags else "none"), fontsize=8, va="top",
             color="#b3261e" if flags else INK, wrap=True)
    fig.text(0.01, 0.89, key, fontsize=6.5, va="top", color="#52514e")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return flags


def summary(R, path):
    """All penetrations on one axis: depth of each marker from the phase reversal."""
    R = R.sort_values(["reversal_clear", "first_sheet_mm"], ascending=[False, True]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 0.28 * len(R) + 1.2))
    for k, r in R.iterrows():
        ax.plot([r.top_depth_mm, r.bottom_depth_mm], [k, k], color="#d9d8d3", lw=4, solid_capstyle="butt")
        for col, c, m in (("sink_depth_mm", BLUE, "s"), ("mua_depth_mm", ORANGE, "o"),
                          ("crossing_depth_mm", "#1baf7a", "^"), ("second_zone_depth_mm", ORANGE, "x")):
            ax.plot([r[col]], [k], marker=m, ms=5, color=c, mfc=c if m != "^" else "none", lw=0)
        ax.text(-3.05, k, r.penetration + ("" if r.reversal_clear else "  (no clear reversal)"), fontsize=6.5,
                va="center", ha="left")
    ax.axvline(0, color=MAG, lw=1)
    ax.set_xlim(-3.1, 3.2)
    ax.set_ylim(len(R) - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlabel("depth from the phase reversal (mm; negative = above)")
    ax.set_title("grey: first sheet (top to bottom) · ■ earliest sink · ● earliest MUA · △ α/β–γ crossing · "
                 "× second early MUA zone", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    P = pens()
    if stage in ("cache", "all"):
        todo = [p for p in P if not (OUT / f"{p}_mono.npz").exists()]
        with ProcessPoolExecutor(4) as ex:
            for r in ex.map(cache_one, todo):
                print("cached", *r, flush=True)
    if stage in ("cards", "all"):
        rows = []
        for p in P:
            row, A = markers(p)
            row["flags"] = "; ".join(card(row, A, OUT / f"{p}_card.png"))
            rows.append(row)
        R = pd.DataFrame(rows)
        R.to_csv(QC / "laminar_audit.csv", index=False)
        summary(R, OUT / "laminar_audit_summary.png")
        pd.set_option("display.width", 250)
        cols = ["penetration", "orientation", "votes", "reversal", "reversal_step_rad", "reversal_clear",
                "top_depth_mm", "bottom_depth_mm", "first_sheet_mm", "sink_depth_mm", "sink_onset_ms",
                "mua_depth_mm", "mua_onset_ms", "crossing_depth_mm", "second_early_zone", "second_zone_depth_mm",
                "second_active_zone", "reversed_steps_below", "spi_rho", "alpha_beta_rho"]
        print(R[cols].to_string(index=False))


if __name__ == "__main__":
    main()
