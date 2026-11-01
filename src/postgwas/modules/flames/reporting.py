"""Presentation-only unified reporting for validated FLAMES results."""

from __future__ import annotations

from html import escape
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from postgwas.core.io.reports import collect_html_reports, write_html_report
from postgwas.modules.flames.errors import FlamesError
from postgwas.modules.flames.stages import FLAMES_STAGES


def _clean(value: Any) -> Any:
    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


def _read_report_table(
    path: Path,
    *,
    delimiter: str,
    required_columns: Sequence[str],
    allow_empty: bool,
) -> dict[str, Any]:
    try:
        frame = pd.read_csv(path, sep=delimiter, dtype=object)
    except (OSError, UnicodeError, ValueError, pd.errors.ParserError) as exc:
        raise FlamesError(
            "Cannot read validated FLAMES report input %s: %s" % (path, exc)
        ) from exc
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise FlamesError(
            "Validated FLAMES report input %s is missing configured columns: %s"
            % (path, ", ".join(missing))
        )
    if frame.empty and not allow_empty:
        raise FlamesError("Validated FLAMES report input contains no rows: %s" % path)
    columns = list(required_columns)
    return {
        "columns": columns,
        "rows": [
            {column: _clean(row[column]) for column in columns}
            for row in frame.to_dict(orient="records")
        ],
        "total_rows": int(len(frame)),
    }


def build_flames_report(
    module,
    *,
    dataset_id: str,
    raw_scores: Path,
    prioritized_genes: Path,
    report_path: Path,
    input_metrics: Mapping[str, Any],
    credible_set_metrics: Sequence[Mapping[str, Any]],
    annotation_metrics: Mapping[str, Any],
    score_metrics: Mapping[str, Any],
    resources: Mapping[str, Any],
    output_paths: Mapping[str, Path],
    upstream_results: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a report model without changing or recalculating FLAMES results."""
    schema = module.result_schema
    raw_columns = list(schema.model_dump().values())
    prediction_columns = [
        schema.locus_column,
        schema.filename_column,
        schema.symbol_column,
        schema.gene_column,
        schema.scaled_score_column,
        schema.raw_score_column,
        schema.precision_column,
    ]
    raw_table = _read_report_table(
        raw_scores,
        delimiter=module.input_schema.table_delimiter,
        required_columns=raw_columns,
        allow_empty=False,
    )
    prediction_table = _read_report_table(
        prioritized_genes,
        delimiter=module.input_schema.table_delimiter,
        required_columns=prediction_columns,
        allow_empty=True,
    )
    if raw_table["total_rows"] != int(score_metrics["scored_genes"]):
        raise FlamesError(
            "FLAMES report raw-score row count disagrees with validated score metrics"
        )
    if prediction_table["total_rows"] != int(score_metrics["prioritized_genes"]):
        raise FlamesError(
            "FLAMES report prioritized-gene row count disagrees with validated score metrics"
        )

    warnings = []
    for source in ("magma", "pops"):
        count = int(annotation_metrics["gene_rows_without_%s_evidence" % source])
        if count:
            warnings.append({
                "check": "Rows without %s evidence" % source.upper(),
                "count": count,
                "detail": (
                    "The published FLAMES implementation assigned zero for this "
                    "absent study-specific evidence; the affected gene IDs are "
                    "listed below."
                ),
            })
    zero_features = list(annotation_metrics["features_zero_in_every_locus"])
    if zero_features:
        warnings.append({
            "check": "Features zero in every locus",
            "count": len(zero_features),
            "detail": (
                "Review annotation coverage and resource compatibility. A zero-only "
                "feature is not by itself proof that the run is invalid."
            ),
        })
    if module.vep_mode == "api" or module.cadd_mode == "api":
        warnings.append({
            "check": "Live annotation service",
            "count": None,
            "detail": (
                "At least one annotation source used a live API. Requests are "
                "audited, but remote content and availability are not version-pinned."
            ),
        })

    stage_details = (
        "%s credible sets; %s variants; %s shared MAGMA/PoPS genes"
        % (
            int(input_metrics["credible_sets"]),
            int(input_metrics["credible_set_variants"]),
            int(input_metrics["shared_magma_pops_genes"]),
        ),
        "%s model features; %s annotation files inventoried; VEP=%s; CADD=%s"
        % (
            len(resources["feature_names"]),
            len(resources["annotation_inventory"]),
            module.vep_mode,
            module.cadd_mode,
        ),
        "%s locus files; %s annotated locus-gene rows"
        % (
            int(annotation_metrics["annotation_files"]),
            int(annotation_metrics["annotated_genes"]),
        ),
        "%s scored rows across %s loci; %s prioritized genes"
        % (
            int(score_metrics["scored_genes"]),
            int(score_metrics["loci"]),
            int(score_metrics["prioritized_genes"]),
        ),
        "%s published and provenance files listed, including this HTML report"
        % len(output_paths),
    )
    return {
        "dataset_id": dataset_id,
        "genome_build": module.genome_build.value,
        "locus_window_bp": int(module.locus_window_bp),
        "minimum_credible_set_coverage": float(
            module.minimum_credible_set_coverage
        ),
        "report_path": str(report_path),
        "report_config": module.html_report.model_dump(mode="json"),
        "summary": {
            "credible_sets": int(input_metrics["credible_sets"]),
            "credible_set_variants": int(input_metrics["credible_set_variants"]),
            "annotated_rows": int(annotation_metrics["annotated_genes"]),
            "scored_rows": int(score_metrics["scored_genes"]),
            "scored_loci": int(score_metrics["loci"]),
            "prioritized_genes": int(score_metrics["prioritized_genes"]),
        },
        "stages": [
            {"Stage": number, "Step": label, "Status": "Completed", "Evidence": detail}
            for number, (label, detail) in enumerate(
                zip(FLAMES_STAGES, stage_details), start=1,
            )
        ],
        "warnings": warnings,
        "credible_sets": [dict(row) for row in credible_set_metrics],
        "raw_table": raw_table,
        "prediction_table": prediction_table,
        "missing_magma_genes": list(annotation_metrics["genes_without_magma_evidence"]),
        "missing_pops_genes": list(annotation_metrics["genes_without_pops_evidence"]),
        "zero_features": zero_features,
        "resources": [
            {"Resource": "FLAMES implementation", "Value": str(resources["script"])},
            {"Resource": "Classifier model", "Value": str(resources["model"])},
            {"Resource": "Model feature manifest", "Value": str(resources["features"])},
            {"Resource": "Annotation bundle", "Value": module.upstream.annotation_bundle_id},
            {"Resource": "Annotation directory", "Value": str(resources["annotation_directory"])},
            {"Resource": "Python runtime", "Value": str(resources["python"])},
            {"Resource": "VEP mode", "Value": module.vep_mode},
            {"Resource": "CADD mode", "Value": module.cadd_mode},
        ],
        "inputs": [
            {"Input": "Fine-mapping index", "Path": str(resources["interchange"]["index_path"])},
            {"Input": "MAGMA gene results", "Path": str(resources["magma"])},
            {"Input": "MAGMA gene-property results", "Path": str(resources["magma_covariate"])},
            {"Input": "PoPS scores", "Path": str(resources["pops"])},
        ],
        "outputs": [
            {"Output": name, "Path": str(path)}
            for name, path in output_paths.items()
        ],
        "upstream_reports": collect_html_reports(upstream_results),
    }


def _title(value: str) -> str:
    labels = {
        "finemap": "Fine-mapping",
        "ld_clump": "LD clumping",
        "magma": "MAGMA",
        "pops": "PoPS",
    }
    if str(value).lower() in labels:
        return labels[str(value).lower()]
    abbreviations = {
        "cadd": "CADD",
        "ld": "LD",
        "magma": "MAGMA",
        "pip": "PIP",
        "pops": "PoPS",
        "xgb": "XGB",
    }
    return " ".join(
        abbreviations.get(part.lower(), part.capitalize())
        for part in str(value).replace("-", "_").split("_")
    )


def _value(value: Any, digits: int, null_value: str) -> str:
    value = _clean(value)
    if value is None:
        return null_value
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if not math.isfinite(value):
            return null_value
        return format(value, ".%dg" % digits)
    return str(value)


def _static_table(
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    *,
    report_path: Path,
    digits: int,
    null_value: str,
    path_columns: Sequence[str] = (),
) -> str:
    if not rows:
        return '<p class="empty">No records were produced for this section.</p>'
    header = "".join("<th>%s</th>" % escape(_title(column)) for column in columns)
    body = []
    for row in rows:
        cells = []
        for column in columns:
            raw = _clean(row.get(column))
            rendered = escape(_value(raw, digits, null_value))
            if raw is not None and column in path_columns:
                href = os.path.relpath(str(raw), report_path.parent)
                rendered = '<a href="%s">%s</a>' % (
                    escape(href, quote=True), rendered,
                )
            cells.append("<td>%s</td>" % rendered)
        body.append("<tr>%s</tr>" % "".join(cells))
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (header, "".join(body))
    )


def _interactive_table(
    table: Mapping[str, Any],
    *,
    table_id: str,
    page_size: int,
    digits: int,
    null_value: str,
) -> str:
    columns = list(table["columns"])
    rows = [
        [_value(row.get(column), digits, null_value) for column in columns]
        for row in table["rows"]
    ]
    if not rows:
        return '<p class="empty">No genes met the published FLAMES prioritization rule.</p>'
    payload = json.dumps(
        {"columns": columns, "rows": rows},
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    ).replace("</", "<\\/")
    header = "".join("<th>%s</th>" % escape(_title(column)) for column in columns)
    return """
<div class="postgwas-data-browser" data-source="%(id)s-data" data-page-size="%(page)d" data-null-value="%(null)s">
  <div class="table-tools"><input type="search" class="table-search" placeholder="Search any gene, locus, or value" aria-label="Search table"><button type="button" class="previous">Previous</button><span class="page-state"></span><button type="button" class="next-page">Next</button></div>
  <div class="table-wrap"><table><thead><tr>%(header)s</tr></thead><tbody></tbody></table></div>
</div><script type="application/json" id="%(id)s-data">%(payload)s</script>""" % {
        "id": escape(table_id, quote=True),
        "page": page_size,
        "null": escape(null_value, quote=True),
        "header": header,
        "payload": payload,
    }


def _gene_list(title: str, values: Sequence[str]) -> str:
    if not values:
        return "<h3>%s</h3><p>None.</p>" % escape(title)
    return "<h3>%s</h3><p><code>%s</code></p>" % (
        escape(title), escape(", ".join(values)),
    )


def render_flames_html(report: Mapping[str, Any]) -> str:
    """Render a self-contained unified FLAMES report."""
    config = report["report_config"]
    digits = int(config["significant_digits"])
    page_size = int(config["page_size"])
    null_value = str(config["null_value"])
    report_path = Path(report["report_path"])
    summary = report["summary"]
    warning_rows = report["warnings"]
    status_class = "warning" if warning_rows else "good"
    status_text = "Completed — review limitations" if warning_rows else "Completed"
    cards = "".join(
        '<div class="metric"><span>%s</span><strong>%s</strong></div>'
        % (escape(_title(key)), escape(_value(value, digits, null_value)))
        for key, value in summary.items()
    )
    upstream = []
    for item in report["upstream_reports"]:
        available = Path(item["path"]).is_file()
        upstream.append({
            "Module": _title(item["module"]),
            "Report": item["path"],
            "Availability": "Available" if available else "Recorded report file not found",
            "_available": available,
        })
    upstream_rows = []
    for row in upstream:
        path_columns = ("Report",) if row.pop("_available") else ()
        upstream_rows.append(_static_table(
            [row], ("Module", "Report", "Availability"),
            report_path=report_path, digits=digits, null_value=null_value,
            path_columns=path_columns,
        ))
    warning_table = _static_table(
        warning_rows, ("check", "count", "detail"), report_path=report_path,
        digits=digits, null_value=null_value,
    )
    credible_table = _static_table(
        report["credible_sets"], ("row_number", "chromosome", "variants", "pip_mass"),
        report_path=report_path, digits=digits, null_value=null_value,
    )
    styles = """
:root{font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif;color:#172033;background:#f5f7fb}*{box-sizing:border-box}body{margin:0}.hero{background:#14213d;color:#fff;padding:28px max(24px,calc((100% - 1500px)/2))}.hero h1{margin:0 0 7px}.hero p{margin:3px 0;color:#dbe6ff}.nav{background:#fff;border-bottom:1px solid #d9e0eb;padding:10px max(24px,calc((100% - 1500px)/2));position:sticky;top:0;z-index:3}.nav a{margin-right:16px;color:#234a91;text-decoration:none}main{max-width:1500px;margin:auto;padding:22px}.panel{background:#fff;border:1px solid #dce2eb;border-radius:12px;padding:20px;margin-bottom:16px;box-shadow:0 3px 14px #26334d10}h2{margin-top:0}.notice{border-left:5px solid #238636;background:#effaf2;padding:13px 16px;margin-bottom:16px}.notice.warning{border-color:#b7791f;background:#fff8e6}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:10px}.metric{border:1px solid #e1e6ee;border-radius:9px;padding:12px}.metric span{display:block;color:#687489;font-size:.79rem}.metric strong{font-size:1.3rem}.table-wrap{overflow:auto;border:1px solid #dce2eb;border-radius:8px;margin-top:10px}table{border-collapse:collapse;width:max-content;min-width:100%;font-size:.79rem}th,td{padding:8px 10px;border-bottom:1px solid #e8ecf2;text-align:left;vertical-align:top;white-space:nowrap}th{background:#eef2f7}tr:nth-child(even){background:#fafbfd}.table-tools{display:flex;gap:9px;align-items:center;flex-wrap:wrap}.table-tools input{min-width:280px;flex:1;padding:8px;border:1px solid #aeb9c8;border-radius:6px}.table-tools button{padding:7px 10px}.empty{color:#667085}.concepts{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}.concepts article{border:1px solid #e1e6ee;border-radius:9px;padding:14px}.concepts h3{margin-top:0}.footer{max-width:1500px;margin:auto;padding:0 22px 28px;color:#667085;font-size:.85rem}code{white-space:normal;overflow-wrap:anywhere}a{color:#2457a6}
"""
    return """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>%(dataset)s unified FLAMES report</title><style>%(styles)s</style></head><body>
<header class="hero"><h1>Unified FLAMES gene-prioritization report</h1><p>Dataset: <strong>%(dataset)s</strong> · genome build: <strong>%(build)s</strong> · published model window: <strong>%(window)s bp</strong></p></header>
<nav class="nav"><a href="#overview">Overview</a><a href="#workflow">Workflow</a><a href="#prioritized">Prioritized genes</a><a href="#all-scores">All scores</a><a href="#qc">QC</a><a href="#provenance">Provenance</a></nav><main>
<div class="notice %(status_class)s"><strong>%(status)s.</strong> Every displayed score row passed the configured schema, finite-value, range, within-locus normalization, and cross-file checks. Completion does not establish that a prioritized gene is causal.</div>
<section class="panel" id="overview"><h2>Results at a glance</h2><div class="metrics">%(cards)s</div></section>
<section class="panel" id="workflow"><h2>Complete FLAMES workflow</h2><p>These are the five validated stages executed by the PostGWAS FLAMES boundary. Pipeline prerequisite reports are linked below when the pipeline supplied them explicitly.</p>%(stages)s<h3>Recorded prerequisite reports</h3>%(upstream)s</section>
<section class="panel" id="prioritized"><h2>Prioritized genes</h2><p>These rows exactly reproduce the validated <code>.pred</code> output. The published FLAMES selection rule is represented by the output itself; this report does not apply another cutoff.</p>%(predictions)s</section>
<section class="panel" id="all-scores"><h2>All scored locus-gene rows</h2><p><code>FLAMES_scaled</code> is normalized to sum to one within each locus and is intended for within-locus ranking. It is not a posterior probability that a gene is causal.</p>%(raw)s</section>
<section class="panel" id="qc"><h2>Scientific QC and limitations</h2>%(warnings)s<h3>Credible-set input audit</h3><p>Each row is one indexed credible set. A sum of marginal PIPs can exceed one for multi-effect fine-mapping models and is not forced to one.</p>%(credible)s%(missing_magma)s%(missing_pops)s%(zero_features)s</section>
<section class="panel"><h2>How to interpret the scores</h2><div class="concepts"><article><h3>FLAMES model score</h3><p>The XGBoost classifier integrates configured variant-to-gene features. Its validity depends on the published model, feature manifest, annotation resources, gene universe, and genome-build compatibility.</p></article><article><h3>Combined ranking score</h3><p><code>FLAMES_raw</code> combines model and PoPS evidence using the published implementation. <code>FLAMES_scaled</code> divides that score by its locus total, so values are comparable for ranking genes within a locus, not as calibrated causal probabilities across loci.</p></article><article><h3>Estimated cumulative precision</h3><p>This is a calibration-derived estimate from the FLAMES benchmark. Transferability is conditional on similarity to the training and evaluation setting; it is not a study-specific posterior probability or proof of causality.</p></article><article><h3>Prioritized genes</h3><p>A prioritized gene is a model-supported candidate for follow-up. Fine-mapping uncertainty, LD quality, incomplete annotations, tissue relevance, and missing upstream evidence remain important limitations.</p></article></div></section>
<section class="panel" id="provenance"><h2>Inputs, resources, and outputs</h2><h3>Scientific inputs</h3>%(inputs)s<h3>Model and annotation resources</h3>%(resources)s<h3>Published outputs</h3>%(outputs)s</section>
</main><footer class="footer">PostGWAS unified FLAMES report · consult the completion manifest, resolved YAML, source tables, API audit, and module logs for the full reproducibility record.</footer></body></html>""" % {
        "dataset": escape(report["dataset_id"]),
        "build": escape(report["genome_build"]),
        "window": f"{report['locus_window_bp']:,}",
        "styles": styles,
        "status_class": status_class,
        "status": escape(status_text),
        "cards": cards,
        "stages": _static_table(
            report["stages"], ("Stage", "Step", "Status", "Evidence"),
            report_path=report_path, digits=digits, null_value=null_value,
        ),
        "upstream": "".join(upstream_rows) if upstream_rows else '<p class="empty">No prerequisite HTML reports were recorded for this direct FLAMES invocation.</p>',
        "predictions": _interactive_table(
            report["prediction_table"], table_id="prioritized-table",
            page_size=page_size, digits=digits, null_value=null_value,
        ),
        "raw": _interactive_table(
            report["raw_table"], table_id="raw-score-table",
            page_size=page_size, digits=digits, null_value=null_value,
        ),
        "warnings": warning_table,
        "credible": credible_table,
        "missing_magma": _gene_list("Genes without MAGMA evidence", report["missing_magma_genes"]),
        "missing_pops": _gene_list("Genes without PoPS evidence", report["missing_pops_genes"]),
        "zero_features": _gene_list("Features zero in every locus", report["zero_features"]),
        "inputs": _static_table(
            report["inputs"], ("Input", "Path"), report_path=report_path,
            digits=digits, null_value=null_value, path_columns=("Path",),
        ),
        "resources": _static_table(
            report["resources"], ("Resource", "Value"), report_path=report_path,
            digits=digits, null_value=null_value,
        ),
        "outputs": _static_table(
            report["outputs"], ("Output", "Path"), report_path=report_path,
            digits=digits, null_value=null_value, path_columns=("Path",),
        ),
    }


def write_flames_html_report(module, *, destination: Path | None = None, **kwargs) -> Path:
    """Build, render, and atomically write the unified FLAMES HTML report."""
    report = build_flames_report(module, **kwargs)
    output = Path(destination) if destination is not None else Path(kwargs["report_path"])
    return write_html_report(render_flames_html(report), output)


__all__ = [
    "build_flames_report",
    "render_flames_html",
    "write_flames_html_report",
]
