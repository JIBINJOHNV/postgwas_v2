"""Standalone reporting for integrated PoPS gene-level results."""

from __future__ import annotations

from html import escape
import json
from pathlib import Path
from typing import Mapping, Sequence

from postgwas.core.io.reports import write_html_report


def write_integrated_gene_report(
    records: Sequence[Mapping[str, object]],
    columns: Sequence[str],
    summary: Mapping[str, object],
    *,
    dataset_id: str,
    tsv_path: Path,
    report_path: Path,
    page_size: int,
    null_value: str,
) -> Path:
    """Write a searchable, sortable, paginated HTML view of the merged table."""
    safe_records = [
        [record.get(column) for column in columns]
        for record in records
    ]
    payload = json.dumps(
        {"columns": list(columns), "rows": safe_records},
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    ).replace("</", "<\\/")
    category_counts = summary.get("status_counts", {})
    cards = "".join(
        '<div class="card"><span>%s</span><strong>%s</strong></div>'
        % (escape(str(label)), f"{int(count):,}")
        for label, count in category_counts.items()
    )
    document = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%(dataset)s integrated PoPS gene results</title>
<style>
:root{font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif;color:#172033;background:#f5f7fb}
body{margin:0;padding:24px}.shell{max-width:1600px;margin:auto}.panel{background:#fff;border:1px solid #dce2eb;border-radius:12px;box-shadow:0 4px 18px #26334d12;padding:20px;margin-bottom:16px}
h1{font-size:1.55rem;margin:0 0 8px}p{color:#566174;line-height:1.45}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:16px}.card{border:1px solid #e2e7ef;border-radius:9px;padding:12px}.card span{display:block;color:#687489;font-size:.82rem}.card strong{font-size:1.35rem}
.tools{display:flex;gap:12px;align-items:center;flex-wrap:wrap}.tools input{min-width:280px;flex:1;padding:9px;border:1px solid #bbc4d2;border-radius:7px}.tools button,.tools a{padding:8px 12px;border:1px solid #9da8b8;border-radius:7px;background:#fff;color:#172033;text-decoration:none;cursor:pointer}.tools button:disabled{opacity:.45}.table-wrap{overflow:auto;max-height:68vh;border:1px solid #dce2eb;border-radius:8px;margin-top:14px}table{border-collapse:collapse;width:max-content;min-width:100%%;font-size:.79rem}th,td{padding:8px 10px;border-bottom:1px solid #e8ecf2;text-align:left;white-space:nowrap}th{position:sticky;top:0;background:#eef2f7;cursor:pointer;z-index:1}tr:nth-child(even){background:#fafbfd}.muted{color:#667085}.pager{display:flex;justify-content:space-between;align-items:center;margin-top:12px;gap:12px}
</style>
</head>
<body><main class="shell">
<section class="panel"><h1>Integrated PoPS gene results</h1>
<p>Dataset: <strong>%(dataset)s</strong>. This full gene union combines the official PoPS predictions with the input target-gene results, compatibility decisions, PoPS model-use flags, and available annotated MAGMA fields. Missing values are displayed as <code>%(null)s</code>.</p>
<div class="cards">%(cards)s</div></section>
<section class="panel postgwas-data-browser" data-source="data" data-page-size="%(page_size)d" data-null-value="%(null_attr)s"><div class="tools"><input id="search" type="search" placeholder="Search any gene or value" aria-label="Search table"><a href="%(tsv)s" download>Download TSV</a></div>
<div class="table-wrap"><table><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div>
<div class="pager"><button id="prev">Previous</button><span id="page" class="muted"></span><button id="next">Next</button></div></section>
</main>
<script id="data" type="application/json">%(payload)s</script>
</body></html>
""" % {
        "dataset": escape(dataset_id),
        "null": escape(null_value),
        "null_attr": escape(null_value, quote=True),
        "cards": cards,
        "tsv": escape(tsv_path.name, quote=True),
        "payload": payload,
        "page_size": page_size,
    }
    return write_html_report(document, report_path)


__all__ = ["write_integrated_gene_report"]
