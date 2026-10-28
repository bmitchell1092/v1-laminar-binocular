"""Build the laminar review page (HTML + report images) from <out_root>/laminar/*.json.

Usage: .venv\\Scripts\\python scripts\\build_laminar_review.py <dest_dir>
Writes <dest_dir>/index.html and <dest_dir>/img/<pen>.png. Verdicts are not part of the
page source: the published page keeps them in its `db` store under verdicts/<penetration>.
"""

import html
import json
import shutil
import sys
from pathlib import Path

from v1lb.config import out_root
from v1lb.laminar_report import legacy_limits

REPO = Path(__file__).resolve().parents[1]
TIER_ORDER = {"A": 0, "B": 1, "C": 2}


def fmt(v, nd=1):
    if v is None:
        return "–"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return html.escape(str(v))
    if f != f:
        return "–"
    return f"{f:.{nd}f}".rstrip("0").rstrip(".") if nd else str(int(f))


def legacy_text(pen):
    lim = legacy_limits(pen, REPO / "docs" / "legacy_TuneList_all_versions.csv")
    if not lim:
        return "–", ""
    items = sorted(lim.items(), key=lambda kv: kv[0][-7:])
    latest = items[-1]
    title = "; ".join(f"{k}: SinkBtm {v[1]:g} (ToC {v[0]:g}, BoC {v[2]:g})" for k, v in items)
    return f"{latest[1][1]:g}", html.escape(title)


def main(dest: Path):
    src = out_root() / "laminar"
    rows = [json.loads(p.read_text()) for p in sorted(src.glob("*_laminar.json"))]
    rows.sort(key=lambda r: (TIER_ORDER.get(r["tier"], 9), -r["confidence"], r["penetration"]))
    (dest / "img").mkdir(parents=True, exist_ok=True)
    for r in rows:
        shutil.copyfile(src / f"{r['penetration']}_laminar.png", dest / "img" / f"{r['penetration']}.png")
    counts = {t: sum(r["tier"] == t for r in rows) for t in "ABC"}
    items = []
    for i, r in enumerate(rows, 1):
        pen = r["penetration"]
        leg, leg_title = legacy_text(pen)
        agree = r["orientation"] == r["spectral_orientation"]
        steps = ", ".join(f"{s[0]:g} ({s[2]:+.2f})" for s in r["phase_steps"]) or "none"
        checks = "".join(f'<li class="{"ok" if v else "no"}">{html.escape(k)}</li>' for k, v in r["checks"].items())
        flags = "".join(f"<li>{html.escape(f)}</li>" for f in r["flags"]) or "<li>none</li>"
        items.append(f"""
<details class="pen" data-tier="{r['tier']}" data-pen="{pen}" id="{pen}">
  <summary>
    <span class="rank">{i}</span>
    <span class="name">{pen}</span>
    <span class="tier t{r['tier']}">{r['tier']}</span>
    <span class="num conf">{r['confidence']:.2f}</span>
    <span class="cx">{html.escape(r['complexity'])}</span>
    <span class="ori{'' if agree else ' warn'}">{r['orientation'][:4]}{'' if agree else ' ≠ ' + r['spectral_orientation'][:4]}</span>
    <span class="num">{fmt(r['top_of_cortex'], 0)}</span>
    <span class="num l45">{fmt(r['l45'])}</span>
    <span class="num">{fmt(r['spectral_crossing'])}</span>
    <span class="num">{fmt(r['earliest_sink'], 0)}</span>
    <span class="num legacy" title="{leg_title}">{leg}</span>
    <span class="verdict-chip" data-chip="{pen}">not reviewed</span>
  </summary>
  <div class="body">
    <div class="figwrap"><img src="img/{pen}.png" alt="Laminar report for {pen}" loading="lazy"></div>
    <div class="facts">
      <section>
        <h3>Calls</h3>
        <dl>
          <dt>Recording</dt><dd>{html.escape(r['recording'])}, {fmt(r['minutes'])} min, {r['n_contacts']} contacts</dd>
          <dt>Orientation</dt><dd>SPI {r['orientation']} (margin {r['spi_fit_margin']:+.2f}); spectrolaminar {r['spectral_orientation']}</dd>
          <dt>Top of cortex</dt><dd>contact {fmt(r['top_of_cortex'], 0)} ({html.escape(r['top_of_cortex_by'])})</dd>
          <dt>Phase steps</dt><dd>{steps} <span class="hint">position (signed rad)</span></dd>
          <dt>L4/5</dt><dd>{fmt(r['l45'])}</dd>
          <dt>Cross-checks</dt><dd>spectral crossing {fmt(r['spectral_crossing'])} (γ peak {r['gamma_peak']}, αβ peak {r['alpha_peak']}); earliest sink {fmt(r['earliest_sink'], 0)} at {fmt(r['earliest_sink_ms'], 0)} ms ({r['n_csd_trials']} trials)</dd>
          <dt>Legacy SinkBtm</dt><dd>{leg_title or '–'}</dd>
        </dl>
      </section>
      <section><h3>Checks</h3><ul class="checks">{checks}</ul></section>
      <section><h3>Flags</h3><ul class="flags">{flags}</ul></section>
      <section class="review">
        <h3>Your verdict</h3>
        <div class="choices" role="group" aria-label="Verdict for {pen}">
          <button type="button" data-v="keep" data-for="{pen}">Keep</button>
          <button type="button" data-v="unsure" data-for="{pen}">Unsure</button>
          <button type="button" data-v="exclude" data-for="{pen}">Exclude</button>
        </div>
        <label for="note-{pen}">Note (corrected boundary, reason to exclude, …)</label>
        <textarea id="note-{pen}" data-note="{pen}" rows="2"></textarea>
        <p class="status" data-status="{pen}"></p>
      </section>
    </div>
  </div>
</details>""")
    page = TEMPLATE.replace("{{ITEMS}}", "\n".join(items)).replace("{{N}}", str(len(rows))) \
        .replace("{{A}}", str(counts["A"])).replace("{{B}}", str(counts["B"])).replace("{{C}}", str(counts["C"]))
    (dest / "index.html").write_text(page, encoding="utf-8")
    print(f"wrote {dest / 'index.html'} with {len(rows)} penetrations")


TEMPLATE = r"""<title>V1 Laminar Calls</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: a ranked list of penetrations; each row opens to its laminar report and a verdict panel. */
:root {
  --bg: #f6f7f9; --surface: #ffffff; --ink: #1b1f27; --muted: #5d6573; --line: #dde1e8;
  --accent: #a3238f; --accent-ink: #ffffff;
  --tA: #1f7a4d; --tB: #a86400; --tC: #b3261e; --keep: #1f7a4d; --exclude: #b3261e; --unsure: #a86400;
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "Cascadia Mono", Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #12151b; --surface: #1a1e26; --ink: #e7eaf0; --muted: #9aa3b2; --line: #2c323d;
  --accent: #e46fd2; --accent-ink: #12151b; --tA: #5cc98f; --tB: #e0a54a; --tC: #f07a70;
  --keep: #5cc98f; --exclude: #f07a70; --unsure: #e0a54a; color-scheme: dark } }
:root[data-theme="dark"] {
  --bg: #12151b; --surface: #1a1e26; --ink: #e7eaf0; --muted: #9aa3b2; --line: #2c323d;
  --accent: #e46fd2; --accent-ink: #12151b; --tA: #5cc98f; --tB: #e0a54a; --tC: #f07a70;
  --keep: #5cc98f; --exclude: #f07a70; --unsure: #e0a54a; color-scheme: dark }
body { background: var(--bg); color: var(--ink); font-family: var(--body); line-height: 1.45;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 24px 48px; }
main { max-width: 1400px; margin: 0 auto; display: grid; gap: 20px; }
h1 { font-family: var(--display); font-weight: 600; font-size: 1.9rem; margin: 0; text-wrap: balance; }
h3 { font-family: var(--display); font-weight: 600; font-size: 0.8rem; text-transform: uppercase;
  letter-spacing: 0.06em; color: var(--muted); margin: 0 0 6px; }
.lede { color: var(--muted); max-width: 70ch; margin: 6px 0 0; }
.strip { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }
.stat { background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 8px 12px;
  display: grid; gap: 2px; min-width: 7rem; }
.stat b { font-family: var(--mono); font-size: 1.25rem; font-variant-numeric: tabular-nums; }
.stat span { font-size: 0.78rem; color: var(--muted); }
.filters { display: flex; flex-wrap: wrap; gap: 6px; margin-left: auto; }
.filters button { font: inherit; font-size: 0.85rem; border: 1px solid var(--line); background: var(--surface);
  color: var(--ink); border-radius: 999px; padding: 4px 12px; cursor: pointer; }
.filters button[aria-pressed="true"] { background: var(--ink); color: var(--bg); border-color: var(--ink); }
.legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 0.82rem; color: var(--muted); }
.legend i { display: inline-block; width: 22px; height: 0; border-top: 3px solid; vertical-align: middle; margin-right: 6px; }
.list { overflow-x: auto; }
.head, summary { display: grid; align-items: center; gap: 10px;
  grid-template-columns: 2.2rem 8.5rem 2.2rem 3.4rem minmax(9rem, 1fr) 6.5rem 3rem 3.4rem 3.6rem 3rem 3.6rem 7.5rem;
  min-width: 960px; }
.head { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted);
  padding: 0 12px 6px; border-bottom: 1px solid var(--line); }
.pen { background: var(--surface); border: 1px solid var(--line); border-radius: 6px; margin-top: 6px; min-width: 960px; }
.pen[hidden] { display: none !important; }
summary { list-style: none; cursor: pointer; padding: 10px 12px; }
summary::-webkit-details-marker { display: none; }
summary:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
.pen[open] summary { border-bottom: 1px solid var(--line); }
.rank { color: var(--muted); font-family: var(--mono); font-size: 0.85rem; }
.name { font-family: var(--mono); font-weight: 500; }
.num { font-family: var(--mono); font-variant-numeric: tabular-nums; text-align: right; }
.l45 { color: var(--accent); font-weight: 500; }
.cx { color: var(--muted); font-size: 0.88rem; }
.ori { font-family: var(--mono); font-size: 0.85rem; }
.ori.warn { color: var(--tC); font-weight: 500; }
.tier { font-family: var(--mono); font-weight: 600; text-align: center; border-radius: 4px; padding: 1px 0;
  color: var(--surface); }
.tA { background: var(--tA); } .tB { background: var(--tB); } .tC { background: var(--tC); }
.verdict-chip { font-size: 0.8rem; color: var(--muted); border: 1px dashed var(--line); border-radius: 999px;
  padding: 2px 8px; text-align: center; }
.verdict-chip.keep { color: var(--keep); border: 1px solid var(--keep); }
.verdict-chip.exclude { color: var(--exclude); border: 1px solid var(--exclude); }
.verdict-chip.unsure { color: var(--unsure); border: 1px solid var(--unsure); }
.body { display: grid; gap: 16px; padding: 14px 12px 16px; }
.figwrap { overflow-x: auto; border: 1px solid var(--line); border-radius: 4px; background: #fff; }
.figwrap img { max-width: none; width: 1800px; display: block; }
.facts { display: grid; gap: 18px; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); }
.facts section { min-width: 0; }
dl { display: grid; grid-template-columns: max-content 1fr; gap: 4px 12px; margin: 0; font-size: 0.9rem; }
dt { color: var(--muted); } dd { margin: 0; min-width: 0; overflow-wrap: anywhere; }
.hint { color: var(--muted); font-size: 0.8rem; }
ul { margin: 0; padding-left: 18px; font-size: 0.9rem; display: grid; gap: 3px; }
.checks { list-style: none; padding: 0; }
.checks li::before { font-family: var(--mono); margin-right: 8px; }
.checks li.ok::before { content: "✓"; color: var(--tA); }
.checks li.no::before { content: "✗"; color: var(--tC); }
.review { display: grid; gap: 8px; align-content: start; }
.choices { display: flex; gap: 6px; flex-wrap: wrap; }
.choices button { font: inherit; border: 1px solid var(--line); background: var(--bg); color: var(--ink);
  border-radius: 4px; padding: 6px 14px; cursor: pointer; }
.choices button:focus-visible, textarea:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
.choices button[aria-pressed="true"][data-v="keep"] { background: var(--keep); border-color: var(--keep); color: var(--surface); }
.choices button[aria-pressed="true"][data-v="exclude"] { background: var(--exclude); border-color: var(--exclude); color: var(--surface); }
.choices button[aria-pressed="true"][data-v="unsure"] { background: var(--unsure); border-color: var(--unsure); color: var(--surface); }
.choices button:disabled { opacity: 0.5; cursor: not-allowed; }
label { font-size: 0.82rem; color: var(--muted); }
textarea { font: inherit; font-size: 0.9rem; background: var(--bg); color: var(--ink); border: 1px solid var(--line);
  border-radius: 4px; padding: 6px 8px; resize: vertical; width: 100%; box-sizing: border-box; }
.status { font-size: 0.8rem; color: var(--muted); margin: 0; min-height: 1.2em; }
.notice { font-size: 0.85rem; color: var(--muted); margin: 0; }
@media (prefers-reduced-motion: reduce) { * { transition: none !important; } }
</style>
<main>
  <header>
    <h1>V1 Laminar Calls</h1>
    <p class="lede">Rule-based laminar alignment for {{N}} penetrations (monkeys B and J), ranked by confidence. Orientation comes from spike-phase coupling (SPI) checked against the spectrolaminar peaks. Top of cortex is the first spiking contact. L4/5 is the 2nd non-reversed phase step below it. Depth units are contacts (100 µm).</p>
  </header>
  <div class="strip">
    <div class="stat"><b>{{A}}</b><span>tier A, clear</span></div>
    <div class="stat"><b>{{B}}</b><span>tier B, caveats</span></div>
    <div class="stat"><b>{{C}}</b><span>tier C, unclear</span></div>
    <div class="stat"><b id="n-reviewed">0</b><span>reviewed of {{N}}</span></div>
    <div class="filters" role="group" aria-label="Filter by tier">
      <button type="button" data-filter="all" aria-pressed="true">All</button>
      <button type="button" data-filter="A" aria-pressed="false">A</button>
      <button type="button" data-filter="B" aria-pressed="false">B</button>
      <button type="button" data-filter="C" aria-pressed="false">C</button>
      <button type="button" data-filter="todo" aria-pressed="false">Not reviewed</button>
    </div>
  </div>
  <div class="legend" aria-label="Report figure legend">
    <span><i style="border-color:#2ca02c"></i>top of cortex</span>
    <span><i style="border-color:#c71fbf"></i>L4/5 call</span>
    <span><i style="border-color:#c71fbf;border-top-style:dashed;border-top-width:2px"></i>canonical step</span>
    <span><i style="border-color:#17becf;border-top-style:dashed;border-top-width:2px"></i>reversed step</span>
    <span><i style="border-color:#ff8c00;border-top-style:dashed;border-top-width:2px"></i>spectral crossing</span>
    <span><i style="border-color:#8c564b;border-top-style:dotted;border-top-width:2px"></i>legacy SinkBtm versions</span>
    <span>★ earliest CSD sink</span>
  </div>
  <p class="notice" id="store-notice">Verdicts are saved for everyone with access to this page.</p>
  <div class="list">
    <div class="head" aria-hidden="true">
      <span>#</span><span>Penetration</span><span>Tier</span><span>Conf.</span><span>Complexity</span><span>Orient.</span>
      <span>ToC</span><span>L4/5</span><span>Spec. ×</span><span>Sink</span><span>Legacy</span><span>Verdict</span>
    </div>
    {{ITEMS}}
  </div>
</main>
<script>
(() => {
  const verdicts = {};
  const LABEL = { keep: "keep", exclude: "exclude", unsure: "unsure" };
  let db = null, filter = "all";

  function render() {
    let reviewed = 0;
    document.querySelectorAll(".pen").forEach(el => {
      const pen = el.dataset.pen, v = verdicts[pen];
      if (v && v.verdict) reviewed++;
      const chip = el.querySelector("[data-chip]");
      chip.className = "verdict-chip" + (v && v.verdict ? " " + v.verdict : "");
      chip.textContent = v && v.verdict ? LABEL[v.verdict] : "not reviewed";
      el.querySelectorAll(".choices button").forEach(b => b.setAttribute("aria-pressed", String(!!v && v.verdict === b.dataset.v)));
      const ta = el.querySelector("[data-note]");
      if (document.activeElement !== ta) ta.value = (v && v.note) || "";
      const show = filter === "all" || (filter === "todo" ? !(v && v.verdict) : el.dataset.tier === filter);
      el.hidden = !show;
    });
    document.getElementById("n-reviewed").textContent = String(reviewed);
  }

  async function save(pen, patch, statusEl) {
    if (!db) { statusEl.textContent = "Saving is not available in this view."; return; }
    const body = Object.assign({ verdict: "", note: "" }, verdicts[pen] || {}, patch, { updated: new Date().toISOString() });
    statusEl.textContent = "Saving…";
    try {
      await db.doc("verdicts/" + pen).set(body);
      verdicts[pen] = body; render();
      statusEl.textContent = "Saved.";
    } catch (e) {
      statusEl.textContent = e && e.code === "invalid_argument"
        ? "You can view this page but not record verdicts."
        : "Could not save (" + ((e && e.code) || "error") + "). Try again.";
    }
  }

  document.querySelectorAll(".filters button").forEach(b => b.addEventListener("click", () => {
    filter = b.dataset.filter;
    document.querySelectorAll(".filters button").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
    render();
  }));
  document.querySelectorAll(".choices button").forEach(b => b.addEventListener("click", () => {
    const pen = b.dataset.for, current = verdicts[pen] && verdicts[pen].verdict;
    save(pen, { verdict: current === b.dataset.v ? "" : b.dataset.v }, document.querySelector(`[data-status="${pen}"]`));
  }));
  document.querySelectorAll("[data-note]").forEach(ta => {
    let timer = null;
    ta.addEventListener("input", () => {
      clearTimeout(timer);
      timer = setTimeout(() => {
        const pen = ta.dataset.note;
        if (((verdicts[pen] && verdicts[pen].note) || "") !== ta.value)
          save(pen, { note: ta.value }, document.querySelector(`[data-status="${pen}"]`));
      }, 900);
    });
  });
  render();

  (async () => {
    const notice = document.getElementById("store-notice");
    try { db = window.claude && window.claude.use ? await window.claude.use("db") : null; } catch (e) { db = null; }
    if (!db) {
      notice.textContent = "Verdicts can't be saved in this view; open the page on claude.ai to record them.";
      document.querySelectorAll(".choices button, [data-note]").forEach(el => el.disabled = true);
      return;
    }
    db.collection("verdicts").onSnapshot(snap => {
      snap.docs.forEach(d => { if (d.exists) verdicts[d.id] = d.data(); });
      snap.docChanges().forEach(c => { if (c.type === "removed") delete verdicts[c.doc.id]; });
      render();
    }, () => { notice.textContent = "Lost the connection to saved verdicts; reload the page."; });
  })();
})();
</script>
"""

if __name__ == "__main__":
    main(Path(sys.argv[1]))
