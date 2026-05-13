"""Build an offline-friendly HTML report comparing simulator results.

Reads ``*_results.json`` files written by ``run_biomodels`` and writes a small
``index.html`` **into the same folder** so the existing JSON files double as
the report's data layer (no copying). Plotly JS is embedded once into the
index; each model's series are fetched on click. Designed to scale to
thousands of models — bundle size grows linearly in data, but only one
model's data is loaded at a time.

Sharing: zip / tar the results folder. Browsers block ``fetch()`` from
``file://`` URLs, so the recipient runs ``python3 -m http.server`` in the
folder; the index shows clear in-page instructions if a fetch fails. Pass
``--zip`` to bundle the folder into a single archive automatically.

Exposes:
- ``build_report(results_dir, output_path=None)`` — function / CLI
- ``ResultsAnalysisStep`` — a process-bigraph Step
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from plotly.offline import get_plotlyjs

from process_bigraph import Step


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load_results(results_dir: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for f in sorted(Path(results_dir).glob("*_results.json")):
        try:
            out.append(json.loads(f.read_text()))
        except json.JSONDecodeError as e:
            print(f"⚠ Skipping {f.name}: {e}")
    return out


# ---------------------------------------------------------------------------
# Per-model summary + bucketing
# ---------------------------------------------------------------------------

# Normalized-RMSE thresholds (mean across shared species).
_BUCKET_THRESHOLDS = [
    ("good",       0.01,  "Good (≤1%)"),
    ("borderline", 0.10,  "Borderline (1–10%)"),
    ("large",      math.inf, "Large diff (>10%)"),
]
_NO_COMPARISON_BUCKET = ("none", "No comparison")


def _union_species(engines: Dict[str, Any]) -> List[str]:
    species: List[str] = []
    seen = set()
    for payload in engines.values():
        if not payload:
            continue
        for c in payload.get("columns", []):
            if c not in seen:
                seen.add(c)
                species.append(c)
    return species


def _compute_model_summary(model: Dict[str, Any]) -> Dict[str, Any]:
    """Compute per-species RMSE, normalized RMSE, and a quality bucket."""
    engines = model.get("engines") or {}
    species = _union_species(engines)

    rmse_by_species: Dict[str, float] = {}
    nrmse_by_species: Dict[str, float] = {}
    n_shared = 0

    for sp in species:
        ys_per_engine: Dict[str, List[float]] = {}
        for eng_name, payload in engines.items():
            if not payload:
                continue
            cols = payload.get("columns", [])
            if sp not in cols:
                continue
            j = cols.index(sp)
            ys_per_engine[eng_name] = [row[j] for row in payload.get("values", [])]
        if len(ys_per_engine) < 2:
            continue
        n_shared += 1
        keys = sorted(ys_per_engine.keys())
        y1, y2 = ys_per_engine[keys[0]], ys_per_engine[keys[1]]
        n = min(len(y1), len(y2))
        if n == 0:
            continue
        rmse = math.sqrt(sum((y1[k] - y2[k]) ** 2 for k in range(n)) / n)
        rmse_by_species[sp] = rmse
        denom = max(
            (abs(v) for v in y1[:n]), default=0.0
        )
        denom = max(denom, max((abs(v) for v in y2[:n]), default=0.0))
        if denom > 0:
            nrmse_by_species[sp] = rmse / denom

    mean_nrmse: Optional[float]
    if nrmse_by_species:
        mean_nrmse = sum(nrmse_by_species.values()) / len(nrmse_by_species)
    else:
        mean_nrmse = None

    if mean_nrmse is None:
        bucket_id, bucket_label = _NO_COMPARISON_BUCKET
    else:
        for bid, threshold, label in _BUCKET_THRESHOLDS:
            if mean_nrmse <= threshold:
                bucket_id, bucket_label = bid, label
                break

    return {
        "n_shared": n_shared,
        "rmse_by_species": rmse_by_species,
        "nrmse_by_species": nrmse_by_species,
        "mean_nrmse": mean_nrmse,
        "bucket": bucket_id,
        "bucket_label": bucket_label,
    }


# ---------------------------------------------------------------------------
# HTML scaffolding
# ---------------------------------------------------------------------------

def _metadata_card_html(model: Dict[str, Any], summary: Dict[str, Any]) -> str:
    bid = model.get("biomodel_id", "")
    md = model.get("metadata") or {}
    name = md.get("name") or bid
    n_species = md.get("n_species", "—")
    n_reactions = md.get("n_reactions", "—")
    n_parameters = md.get("n_parameters", "—")
    compartments = md.get("compartments") or []
    compartments_str = ", ".join(compartments) if compartments else "—"
    n_points = md.get("n_points")
    duration = model.get("duration")
    time_unit = model.get("time_unit") or ""
    if duration is None:
        sim_str = "—"
    else:
        n_pts = n_points if n_points is not None else "—"
        sim_str = f"{duration:g} {time_unit}".strip() + f", {n_pts} time points"
    url = md.get("biomodels_url")
    link_html = (
        f' <a href="{url}" target="_blank" rel="noreferrer noopener">BioModels page ↗</a>'
        if url
        else ""
    )
    mean_nrmse = summary.get("mean_nrmse")
    nrmse_str = f"{mean_nrmse:.4g}" if mean_nrmse is not None else "—"
    rows = [
        ("Name", name),
        ("BioModels ID", f"{bid}{link_html}"),
        ("Species", n_species),
        ("Reactions", n_reactions),
        ("Compartments", compartments_str),
        ("Parameters", n_parameters),
        ("Simulation", sim_str),
        ("Mean nRMSE (shared species)", f"{nrmse_str} — {summary.get('bucket_label','')}"),
    ]
    body = "\n".join(
        f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows
    )
    return f'<table class="meta">{body}</table>'


def _summary_table_html(
    models: List[Dict[str, Any]], summaries: List[Dict[str, Any]]
) -> str:
    headers = [
        ("Model",            "string"),
        ("Engines",          "string"),
        ("#species copasi",  "number"),
        ("#species tellurium", "number"),
        ("#shared",          "number"),
        ("mean nRMSE",       "number"),
        ("Bucket",           "string"),
    ]
    head = '<table id="summary"><thead><tr>' + "".join(
        f'<th class="sortable" data-sort-type="{t}">{label}<span class="sort-arrow"></span></th>'
        for label, t in headers
    ) + "</tr></thead><tbody>"

    body_rows = []
    for m, s in zip(models, summaries):
        engines = m.get("engines") or {}
        nc = len((engines.get("copasi") or {}).get("columns", [])) if engines.get("copasi") else 0
        nt = len((engines.get("tellurium") or {}).get("columns", [])) if engines.get("tellurium") else 0
        engines_present = ", ".join(k for k, v in engines.items() if v) or "—"
        bid = m.get("biomodel_id", "")
        mn = s.get("mean_nrmse")
        mn_str = f"{mn:.4g}" if mn is not None else "—"
        mn_sort = mn if mn is not None else float("inf")
        body_rows.append(
            "<tr>"
            f'<td data-sort-value="{bid}"><a href="#" class="nav-target" data-target="model-{bid}">{bid}</a></td>'
            f'<td data-sort-value="{engines_present}">{engines_present}</td>'
            f'<td data-sort-value="{nc}">{nc}</td>'
            f'<td data-sort-value="{nt}">{nt}</td>'
            f'<td data-sort-value="{s["n_shared"]}">{s["n_shared"]}</td>'
            f'<td data-sort-value="{mn_sort}">{mn_str}</td>'
            f'<td data-sort-value="{s["bucket"]}">{s.get("bucket_label","")}</td>'
            "</tr>"
        )
    return head + "\n".join(body_rows) + "</tbody></table>"


def _sidebar_html(
    models: List[Dict[str, Any]], summaries: List[Dict[str, Any]]
) -> str:
    """Render sidebar grouped by RMSE bucket with section counts."""
    bucket_order = ["large", "borderline", "good", "none"]
    bucket_labels = {
        "good":       "Good (≤1%)",
        "borderline": "Borderline (1–10%)",
        "large":      "Large diff (>10%)",
        "none":       "No comparison",
    }
    by_bucket: Dict[str, List[Tuple[Dict[str, Any], Dict[str, Any]]]] = {b: [] for b in bucket_order}
    for m, s in zip(models, summaries):
        by_bucket.setdefault(s["bucket"], []).append((m, s))

    sections: List[str] = []
    for b in bucket_order:
        items = by_bucket.get(b) or []
        if not items:
            continue
        # within a bucket, worst (highest nRMSE) first; None at the bottom
        items.sort(key=lambda x: -(x[1].get("mean_nrmse") or -1))
        links = "\n".join(
            f'<li><a href="#" class="nav-target" data-target="model-{m["biomodel_id"]}" '
            f'data-name="{(m.get("metadata") or {}).get("name", "")}">'
            f'{m["biomodel_id"]}</a></li>'
            for m, _ in items
        )
        sections.append(
            f'<div class="group" data-bucket="{b}">'
            f'<div class="group-header">{bucket_labels[b]} <span class="count">{len(items)}</span></div>'
            f'<ul>{links}</ul>'
            f'</div>'
        )
    return "\n".join(sections)


# ---------------------------------------------------------------------------
# CSS / JS
# ---------------------------------------------------------------------------

_CSS = """
:root { color-scheme: light; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 0; display: flex; min-height: 100vh; color: #222; }
nav { width: 280px; background: #f5f6f8; padding: 1rem; box-sizing: border-box; border-right: 1px solid #d8dbe0; position: sticky; top: 0; height: 100vh; overflow-y: auto; }
nav h1 { font-size: 1rem; margin: 0 0 0.6rem 0; color: #444; }
nav input[type="search"] { width: 100%; box-sizing: border-box; padding: 0.4rem 0.6rem; border: 1px solid #c8cdd5; border-radius: 4px; font-size: 0.9rem; margin-bottom: 0.8rem; }
nav .group { margin-bottom: 0.6rem; }
nav .group-header { font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.04em; color: #5d6573; padding: 0.25rem 0.5rem; }
nav .group-header .count { float: right; background: #d8dbe0; color: #333; border-radius: 9px; padding: 0 6px; font-weight: 600; }
nav .group[data-bucket="large"] .group-header { color: #b3261e; }
nav .group[data-bucket="borderline"] .group-header { color: #b8741a; }
nav .group[data-bucket="good"] .group-header { color: #1b6e3c; }
nav ul { list-style: none; padding: 0; margin: 0; }
nav a { display: block; padding: 0.35rem 0.7rem; color: #2a5dbf; text-decoration: none; border-radius: 4px; font-size: 0.9rem; margin-bottom: 1px; }
nav a:hover { background: #e7ebf5; }
nav a.active { background: #2a5dbf; color: white; }
nav .overview-link { font-weight: 600; margin-bottom: 0.6rem; }
main { flex: 1; padding: 1.5rem 2rem; overflow-x: auto; }
.panel { display: none; }
.panel.active { display: block; }
h2 { margin-top: 0; }
.placeholder { color: #888; font-style: italic; padding: 1rem; }
table { border-collapse: collapse; margin-bottom: 1rem; font-size: 0.92rem; }
th, td { border: 1px solid #d8dbe0; padding: 5px 10px; text-align: left; }
th { background: #f0f2f5; font-weight: 600; }
th.sortable { cursor: pointer; user-select: none; }
th.sortable:hover { background: #e6e9ee; }
.sort-arrow { display: inline-block; width: 0.8em; color: #999; }
th[data-sort-dir="asc"] .sort-arrow::after { content: "▲"; color: #333; }
th[data-sort-dir="desc"] .sort-arrow::after { content: "▼"; color: #333; }
td a { color: #2a5dbf; text-decoration: none; }
td a:hover { text-decoration: underline; }
table.meta { max-width: 600px; }
table.meta th { background: #fafbfd; width: 200px; }
.error-card { border: 1px solid #d8a; background: #fff5f7; padding: 0.8rem 1rem; border-radius: 6px; max-width: 720px; }
.error-card pre { background: #f3f4f6; padding: 0.5rem 0.7rem; border-radius: 4px; overflow-x: auto; font-size: 0.85rem; }
.error-card code { background: #f3f4f6; padding: 0 4px; border-radius: 3px; }
"""

# IMPORTANT: this is JS embedded in an f-string-built HTML page, but the JS
# block itself is NOT an f-string — keep raw braces for JS object literals.
_JS = r"""
(function () {
  var FETCH_MAP = {};
  try { FETCH_MAP = JSON.parse(document.getElementById('fetch-map').textContent); } catch (e) {}
  var RENDERED = new Set();
  var INFLIGHT = new Map();

  function speciesLabel(sp, units) {
    var u = units && units[sp];
    return u ? sp + ' (' + u + ')' : sp;
  }

  function renderFetchError(target, modelId, path, err) {
    var fname = (window.location.pathname.split('/').pop() || 'index.html');
    target.innerHTML =
      '<div class="error-card">' +
      '<p><strong>Could not load <code>' + path + '</code>:</strong> ' + (err && err.message || err) + '</p>' +
      '<p>Browsers block <code>fetch()</code> from <code>file://</code> URLs. Serve this folder over HTTP, then re-open in the browser:</p>' +
      '<pre>cd "' + decodeURIComponent(window.location.pathname.replace(/[^\/]*$/, '')) + '"\npython3 -m http.server 8000</pre>' +
      '<p>Then visit <a href="http://localhost:8000/' + fname + '">http://localhost:8000/' + fname + '</a></p>' +
      '</div>';
  }

  // Convert the raw {time, columns, values} per-engine shape (from the runner)
  // into the {species, time_unit, species_units, engines{time,values_by_species}}
  // shape Plotly construction expects.
  function normalizeResults(raw) {
    var engines = raw.engines || {};
    var species = [];
    var seen = {};
    Object.keys(engines).forEach(function (engName) {
      var p = engines[engName];
      if (!p) return;
      (p.columns || []).forEach(function (sp) {
        if (!seen[sp]) { seen[sp] = 1; species.push(sp); }
      });
    });
    var outEngines = {};
    Object.keys(engines).forEach(function (engName) {
      var p = engines[engName];
      if (!p) return;
      var cols = p.columns || [];
      var values = p.values || [];
      var vbs = {};
      cols.forEach(function (sp, j) {
        var col = new Array(values.length);
        for (var k = 0; k < values.length; k++) col[k] = values[k][j];
        vbs[sp] = col;
      });
      outEngines[engName] = { time: p.time || [], values_by_species: vbs };
    });
    return {
      species: species,
      species_units: raw.species_units || {},
      time_unit: raw.time_unit,
      engines: outEngines,
    };
  }

  function plotData(target, data) {
    var species = data.species || [];
    if (species.length === 0) {
      target.innerHTML = '<p class="placeholder">No data for this model.</p>';
      return;
    }
    var cols = 3;
    var rows = Math.ceil(species.length / cols);
    var traces = [];
    var layout = {
      grid: { rows: rows, columns: cols, pattern: 'independent' },
      height: 240 * rows + 80,
      legend: { orientation: 'h', y: 1.04, x: 0 },
      margin: { t: 60, b: 40, l: 60, r: 20 },
    };
    var seenLegend = new Set();
    species.forEach(function (sp, i) {
      var r = Math.floor(i / cols);
      var c = i % cols;
      var idx = i + 1;
      var xKey = 'xaxis' + (idx === 1 ? '' : idx);
      var yKey = 'yaxis' + (idx === 1 ? '' : idx);
      var xRef = idx === 1 ? 'x' : 'x' + idx;
      var yRef = idx === 1 ? 'y' : 'y' + idx;
      layout[yKey] = { title: { text: speciesLabel(sp, data.species_units) } };
      var lastRowForCol = Math.floor((species.length - 1 - c) / cols);
      if (r === lastRowForCol) {
        layout[xKey] = { title: { text: data.time_unit ? 'time (' + data.time_unit + ')' : 'time' } };
      } else {
        layout[xKey] = {};
      }
      Object.keys(data.engines || {}).forEach(function (engName) {
        var eng = data.engines[engName];
        var ys = eng && eng.values_by_species && eng.values_by_species[sp];
        if (!ys) return;
        var showInLegend = !seenLegend.has(engName);
        seenLegend.add(engName);
        traces.push({
          x: eng.time,
          y: ys,
          mode: 'lines',
          name: engName,
          legendgroup: engName,
          showlegend: showInLegend,
          xaxis: xRef,
          yaxis: yRef,
        });
      });
    });
    Plotly.newPlot(target, traces, layout, { responsive: true, displaylogo: false });
  }

  function buildPlot(modelId) {
    if (RENDERED.has(modelId)) return;
    if (INFLIGHT.has(modelId)) return;
    var target = document.getElementById('plot-' + modelId);
    if (!target) return;
    var path = FETCH_MAP[modelId];
    if (!path) {
      target.innerHTML = '<p class="placeholder">No data file mapped for ' + modelId + '.</p>';
      RENDERED.add(modelId);
      return;
    }
    target.innerHTML = '<p class="placeholder">Loading…</p>';
    var p = fetch(path)
      .then(function (r) {
        if (!r.ok) throw new Error('HTTP ' + r.status);
        return r.json();
      })
      .then(function (raw) {
        plotData(target, normalizeResults(raw));
        RENDERED.add(modelId);
      })
      .catch(function (err) { renderFetchError(target, modelId, path, err); })
      .finally(function () { INFLIGHT.delete(modelId); });
    INFLIGHT.set(modelId, p);
  }

  // Navigation
  var navLinks = document.querySelectorAll('nav a.nav-target');
  var allLinks = document.querySelectorAll('a.nav-target');
  var panels = document.querySelectorAll('section.panel');
  function show(id) {
    panels.forEach(function (p) { p.classList.toggle('active', p.id === id); });
    navLinks.forEach(function (a) { a.classList.toggle('active', a.dataset.target === id); });
    if (id && id.indexOf('model-') === 0) {
      buildPlot(id.substring('model-'.length));
    }
    window.scrollTo(0, 0);
  }
  allLinks.forEach(function (a) {
    a.addEventListener('click', function (e) {
      e.preventDefault();
      show(a.dataset.target);
    });
  });

  // Sidebar search
  var search = document.getElementById('sidebar-search');
  if (search) {
    search.addEventListener('input', function () {
      var q = search.value.trim().toLowerCase();
      document.querySelectorAll('nav .group').forEach(function (grp) {
        var visible = 0;
        grp.querySelectorAll('a.nav-target').forEach(function (a) {
          var hay = (a.textContent + ' ' + (a.dataset.name || '')).toLowerCase();
          var match = !q || hay.indexOf(q) !== -1;
          a.parentNode.style.display = match ? '' : 'none';
          if (match) visible++;
        });
        grp.style.display = (q && visible === 0) ? 'none' : '';
        var counter = grp.querySelector('.count');
        if (counter) counter.textContent = visible;
      });
    });
  }

  // Sortable summary table
  document.querySelectorAll('table#summary th.sortable').forEach(function (th) {
    th.addEventListener('click', function () {
      var table = th.closest('table');
      var tbody = table.querySelector('tbody');
      var rows = Array.from(tbody.querySelectorAll('tr'));
      var headerCells = Array.from(th.parentNode.children);
      var idx = headerCells.indexOf(th);
      var type = th.dataset.sortType || 'string';
      var dir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';
      headerCells.forEach(function (other) { other.removeAttribute('data-sort-dir'); });
      th.dataset.sortDir = dir;
      rows.sort(function (a, b) {
        var av = a.children[idx].dataset.sortValue;
        var bv = b.children[idx].dataset.sortValue;
        if (av === undefined) av = a.children[idx].textContent;
        if (bv === undefined) bv = b.children[idx].textContent;
        if (type === 'number') {
          av = parseFloat(av); bv = parseFloat(bv);
          if (isNaN(av)) av = Infinity;
          if (isNaN(bv)) bv = Infinity;
        }
        var cmp = av < bv ? -1 : av > bv ? 1 : 0;
        return dir === 'asc' ? cmp : -cmp;
      });
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });
})();
"""


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def _discover_models(results_dir: Path) -> List[Tuple[Dict[str, Any], Path]]:
    """Return ``[(model_dict, source_file_path), ...]`` from a results dir."""
    out: List[Tuple[Dict[str, Any], Path]] = []
    for f in sorted(results_dir.glob("*_results.json")):
        try:
            d = json.loads(f.read_text())
        except json.JSONDecodeError as e:
            print(f"⚠ Skipping {f.name}: {e}")
            continue
        if d.get("biomodel_id"):
            out.append((d, f))
    return out


def build_report(
    results_dir: str,
    output_path: Optional[str] = None,
    make_zip: bool = False,
) -> Path:
    """Write ``index.html`` into the results folder so existing ``*_results.json``
    files double as the data layer. Returns the path to the index file.

    - ``output_path=None`` → ``<results_dir>/index.html``
    - ``output_path`` ending in ``.html`` → write index there; data files are
      referenced relative to the index.
    - any other ``output_path`` → treated as a directory; the index goes in it
      and points back at ``results_dir`` via relative paths.

    If ``make_zip`` is true, also write ``<bundle_dir>.zip`` containing the
    index plus every referenced ``*_results.json`` for one-file sharing.
    """
    results_dir = Path(results_dir)
    discovered = _discover_models(results_dir)
    if not discovered:
        raise ValueError(f"No *_results.json files found in {results_dir}")

    if output_path is None:
        index_path = results_dir / "index.html"
    else:
        op = Path(output_path)
        if op.suffix.lower() == ".html":
            index_path = op
        else:
            index_path = op / "index.html"
    index_dir = index_path.parent
    index_dir.mkdir(parents=True, exist_ok=True)

    models = [m for m, _ in discovered]
    summaries = [_compute_model_summary(m) for m in models]

    fetch_map: Dict[str, str] = {}
    for m, src in discovered:
        bid = m["biomodel_id"]
        # filenames are stored relative to where index.html lives
        fetch_map[bid] = os.path.relpath(src.resolve(), index_dir.resolve())

    panels_html: List[str] = []
    for m, s in zip(models, summaries):
        bid = m["biomodel_id"]
        panels_html.append(
            f'<section class="panel" id="model-{bid}">'
            f'<h2>{bid}</h2>'
            f'{_metadata_card_html(m, s)}'
            f'<div class="plot-target" id="plot-{bid}"><p class="placeholder">Click to render plot.</p></div>'
            f'</section>'
        )

    sidebar = _sidebar_html(models, summaries)
    summary_table = _summary_table_html(models, summaries)
    plotly_js = get_plotlyjs()
    fetch_map_json = json.dumps(fetch_map).replace("</", "<\\/")

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Biomodels regression report</title>
<style>{_CSS}</style>
<script>{plotly_js}</script>
</head>
<body>
<nav>
  <h1>Biomodels regression</h1>
  <input id="sidebar-search" type="search" placeholder="Filter models…" autocomplete="off">
  <a href="#" class="nav-target overview-link active" data-target="overview">Overview</a>
  {sidebar}
</nav>
<main>
  <section class="panel active" id="overview">
    <h2>Overview</h2>
    <p>{len(models)} models. Sidebar groups models by mean normalized RMSE between simulators; click any model to lazy-fetch and render its plot.</p>
    {summary_table}
  </section>
  {''.join(panels_html)}
</main>
<script type="application/json" id="fetch-map">{fetch_map_json}</script>
<script>{_JS}</script>
</body>
</html>
"""

    index_path.write_text(html, encoding="utf-8")

    if make_zip:
        import zipfile
        bundle_dir = index_dir
        zip_path = bundle_dir.with_suffix(".zip") if bundle_dir.suffix else Path(str(bundle_dir) + ".zip")
        zip_root = bundle_dir.name or "report"
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            zf.write(index_path, arcname=f"{zip_root}/{index_path.name}")
            for m, src in discovered:
                arc = f"{zip_root}/{os.path.relpath(src.resolve(), index_dir.resolve())}"
                zf.write(src, arcname=arc)
        print(f"📦 Bundled to {zip_path}")

    return index_path


# ---------------------------------------------------------------------------
# Step
# ---------------------------------------------------------------------------

class ResultsAnalysisStep(Step):
    config_schema = {
        "results_dir": "string",
        # When empty, the index is written into ``results_dir`` itself.
        "output_path": {"_type": "string", "_default": ""},
        "make_zip": {"_type": "boolean", "_default": False},
    }

    def inputs(self) -> Dict[str, Any]:
        return {}

    def outputs(self) -> Dict[str, Any]:
        return {"report_path": "string"}

    def update(self, state) -> Dict[str, Any]:
        path = build_report(
            self.config["results_dir"],
            self.config.get("output_path") or None,
            make_zip=bool(self.config.get("make_zip", False)),
        )
        print(f"📊 Wrote report: {path}")
        return {"report_path": str(path)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build a sharable comparison report bundled with simulator results.",
    )
    parser.add_argument("--results-dir", default="out_biomodels")
    parser.add_argument(
        "--output",
        default=None,
        help="Where to write index.html. Default: <results-dir>/index.html.",
    )
    parser.add_argument(
        "--zip",
        action="store_true",
        help="Also produce <bundle>.zip containing the index and result files.",
    )
    args = parser.parse_args()
    p = build_report(args.results_dir, args.output, make_zip=args.zip)
    bundle_dir = p.parent
    print(f"Wrote {p}")
    print()
    print("To view (browsers block file:// fetches, so a tiny server is needed):")
    print(f"  cd {bundle_dir}")
    print( "  python3 -m http.server 8000")
    print( "  open http://localhost:8000/")
    if not args.zip:
        print()
        print("To package for sharing:")
        print(f"  python -m biomodels_regression.analysis --results-dir {args.results_dir} --zip")
