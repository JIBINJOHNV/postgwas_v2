"""Standalone HTML views of persisted harmonisation evidence.

The renderer is deliberately presentation-only.  It accepts the already-built
run-summary records and dataset manifests; it never opens a summary-statistics
file or VCF and it does not derive scientific metrics.
"""

from __future__ import annotations

from html import escape
import json
import math
import os
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence
from urllib.parse import quote, urlsplit

from postgwas.core.io.reports import write_html_report
from postgwas.core.pipeline_logging import _elapsed
from postgwas.modules.harmonisation.report_steps import (
    chromosome_steps, dataset_steps, post_merge_steps, summarise_chromosome_evidence,
)


_STATUS_STYLES = {
    "OK": "ok",
    "PARTIAL": "warning",
    "RUNNING": "running",
    "NOT_RUN": "muted",
    "COMPLETED": "ok",
    "PASS": "ok",
    "WARNING": "warning",
    "COMPLETED WITH WARNINGS": "warning",
    "NOT REQUESTED": "muted",
    "NOT NEEDED": "muted",
    "NOT RECORDED": "muted",
    "EVIDENCE RECORDED": "running",
    "NOT RUN — INCOMPLETE INPUT": "muted",
    "FAILED — CONTINUED BY POLICY": "error",
}


def _text(value: Any) -> str:
    if value is None or value == "":
        return "Not recorded"
    if value is True:
        return "Yes"
    if value is False:
        return "No"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.6g}" if math.isfinite(value) else "Non-finite value recorded"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _cell(value: Any, *, percent: bool = False) -> str:
    if percent and isinstance(value, (int, float)) and not isinstance(value, bool):
        return escape(_percentage(value))
    return escape(_text(value))


def _status(value: Any) -> str:
    return str(value or "NOT_RUN").upper()


def _status_badge(value: Any) -> str:
    status = _status(value)
    css_class = _STATUS_STYLES.get(status, "error")
    label = status if "_" in str(value) or status in {"OK", "PARTIAL", "RUNNING"} else status.lower().capitalize()
    return '<span class="badge %s">%s</span>' % (
        css_class, escape(label.replace("_", " ")),
    )


def _safe_id(value: Any) -> str:
    # Reversible encoding keeps distinct dataset names distinct; ':' cannot
    # occur in an unencoded ID, so encoded and unencoded names cannot collide.
    text = str(value)
    return text if re.fullmatch(r"[A-Za-z0-9_.-]+", text) else "u:" + text.encode().hex()


def _href(value: Any, report_path: str | Path | None = None) -> str | None:
    if value is None or value == "":
        return None
    text = str(value)
    parsed = urlsplit(text)
    if parsed.scheme in {"file", "http", "https"}:
        return text
    if parsed.scheme:
        return None
    try:
        path = Path(text).expanduser()
        if report_path is not None:
            parent = Path(report_path).expanduser().absolute().parent
            path = path if path.is_absolute() else parent / path
            return quote(os.path.relpath(path, parent), safe="/.-_")
        return path.absolute().as_uri()
    except (OSError, ValueError):
        return None


def _path_value(value: Any, report_path: str | Path | None = None, *, check_exists: bool = False) -> str:
    if value is None or value == "":
        return '<span class="muted-text">Not recorded</span>'
    raw = str(value)
    label = escape(Path(raw).name)
    full_path = '<details class="file-path"><summary>Full path</summary><span class="path">%s</span></details>' % escape(raw)
    href = _href(value, report_path)
    if href is None:
        return '<span class="path">%s (unsupported link)</span>' % escape(raw)
    if check_exists and not urlsplit(raw).scheme:
        path = Path(raw).expanduser()
        if not path.is_absolute() and report_path is not None:
            path = Path(report_path).parent / path
        if not path.is_file():
            return '<span class="muted-text">%s — not present at report generation</span>%s' % (label, full_path)
    return '<a class="path" href="%s">%s</a>%s' % (escape(href, quote=True), label, full_path)


def _percentage(value: Any) -> str:
    """Do not round a nonzero finding to zero, or a partial result to 100%."""
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
        return _text(value)
    displayed = f"{value:.4f}"
    if 0 < value < 100 and displayed == "0.0000":
        displayed = "<0.0001"
    elif 0 < value < 100 and displayed == "100.0000":
        displayed = ">99.9999"
    return displayed + "%"


def _fraction(value: Any) -> str:
    """Display a recorded fraction as a percentage; do not invent a denominator."""
    return _percentage(value * 100) if isinstance(value, (int, float)) and not isinstance(value, bool) else _text(value)


def _ratio(numerator: Any, denominator: Any, basis: str) -> str:
    """Presentation arithmetic only, using compatible recorded counts."""
    text = f"{_text(numerator)} / {_text(denominator)}"
    if isinstance(numerator, (int, float)) and isinstance(denominator, (int, float)) and not isinstance(numerator, bool) and not isinstance(denominator, bool) and math.isfinite(numerator) and math.isfinite(denominator) and denominator > 0:
        text += f" ({_percentage(100 * numerator / denominator)} of {basis})"
    return text


def _metric_cards(items: Iterable[tuple[str, Any, str]]) -> str:
    cards = []
    for label, value, kind in items:
        if kind == "status":
            rendered = _status_badge(value)
        elif kind == "percent":
            rendered = _cell(value, percent=True)
        else:
            rendered = _cell(value)
        cards.append(
            '<div class="metric"><div class="metric-label">%s</div>'
            '<div class="metric-value">%s</div></div>'
            % (escape(label), rendered)
        )
    return '<div class="metrics">%s</div>' % "".join(cards)


def _key_value_table(
    items: Iterable[tuple[str, Any]],
    *,
    path_labels: Iterable[str] = (),
    percent_labels: Iterable[str] = (),
) -> str:
    path_names = set(path_labels)
    percent_names = set(percent_labels)
    rows = []
    for label, value in items:
        if label in path_names:
            rendered = _path_value(value)
        else:
            rendered = _cell(value, percent=label in percent_names)
        rows.append(
            "<tr><th scope=\"row\">%s</th><td>%s</td></tr>"
            % (escape(label), rendered)
        )
    return '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody></table></div>' % "".join(rows)


def _section(title: str, body: str, *, subtitle: str | None = None) -> str:
    description = "" if not subtitle else '<p class="section-note">%s</p>' % escape(subtitle)
    return '<section><h4>%s</h4>%s%s</section>' % (escape(title), description, body)


def _mapping_details(title: str, mapping: Mapping[str, Any]) -> str:
    if not mapping:
        return ""
    rows = _key_value_table(
        ("Most similar reference population" if key == "closest_population" else str(key).replace("_", " ").capitalize(), "Not set (null)" if value is None else "Empty string" if value == "" else value)
        for key, value in mapping.items() if not isinstance(value, (Mapping, list, tuple))
    )
    children = []
    for key, value in mapping.items():
        label = str(key).replace("_", " ").capitalize()
        if isinstance(value, Mapping):
            children.append(_mapping_details(label, value))
        elif isinstance(value, (list, tuple)):
            if not value:
                children.append(_key_value_table([(label, "None recorded")]))
            elif all(not isinstance(item, Mapping) for item in value):
                children.append(_key_value_table([(label, "; ".join("Empty string" if item == "" else _text(item) for item in value))]))
            else:
                children.append(_mapping_details(label, {str(i): item for i, item in enumerate(value, 1)}))
    return '<details><summary>%s</summary>%s%s</details>' % (escape(title), rows, "".join(children))


def _stage_card(step: Mapping[str, Any], *, level: int = 5) -> str:
    """One shared view for dataset, chromosome and post-merge step evidence."""
    summary = "".join('<p>%s</p>' % escape(str(line)) for line in step.get("summary", []))
    rules = step.get("policies") or {}
    details = dict(step.get("details") or {})
    recorded = details.pop("Recorded stage", {})
    # The merged evidence already contains step extras. Avoid printing them a
    # second time; retain the call identity and failures for reproducibility.
    call = {key: recorded[key] for key in ("func", "error", "failure_hint") if key in recorded}
    if call:
        details["Recorded call"] = call
    return (
        '<div class="stage-card"><h{level}>{number}. {title} {badge}</h{level}>'
        '<p class="section-note">{purpose}</p>{summary}{metrics}{groups}{rules}{details}</div>'
    ).format(
        level=level, number=step["number"], title=escape(step["title"]),
        badge=_status_badge(step["status"]), purpose=escape(step["purpose"]),
        summary=summary or '<p class="muted-text">No step-specific outcome was recorded.</p>',
        metrics=_key_value_table(("Duration" if key == "Elapsed seconds" else key, _elapsed(value) if key == "Elapsed seconds" and value is not None else value)
                                 for key, value in (step.get("metrics") or {}).items()),
        rules=_mapping_details("Rule applied — recorded configuration", rules),
        details=_mapping_details("Supporting evidence and recorded decisions", details),
        groups="".join('<h{level}>{title}</h{level}><p>{note}</p>{metrics}'.format(
            level=min(level + 1, 6), title=escape(group["title"]),
            note=escape(group.get("note", "")), metrics=_key_value_table(group.get("metrics", {}).items()),
        ) for group in step.get("groups", [])),
    )


def _input_section(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    values = manifest.get("pre_vcf_summary") or {}
    rows = (
        ("Input variants", values.get("total_variant_infile", record.get("input_variants"))),
        ("Rows read", values.get("total_variant_read", record.get("parsed_variants"))),
        ("Missing required values", values.get("total_variant_removed_missing_values")),
        ("Invalid coordinates", values.get("total_variant_removed_null_coords")),
        ("Unsupported chromosomes", values.get("total_variant_removed_unsupported_chromosomes")),
        ("Non-standard alleles", values.get("total_variant_removed_non_standard_alleles")),
        ("Duplicate variants", values.get("total_variant_removed_duplicates")),
        ("Ready for harmonisation", values.get("total_variant_remaining_for_harmonisation", record.get("input_ready_variants"))),
        ("Ready SNPs", values.get("total_variant_ready_snps", record.get("input_ready_snps"))),
        ("Ready indels / other", values.get("total_variant_ready_indels_or_other", record.get("input_ready_indels_or_other_variants"))),
        ("Chromosome-stage removals", values.get("total_variant_removed_chromosome_harmonisation")),
        ("All rejected variants", values.get("total_variant_rejected_end_to_end", record.get("rejected_variants"))),
        ("Sent to GWAS-to-VCF", values.get("total_variant_in_vcf_input", record.get("harmonised_variants"))),
        ("End-to-end accounting balanced", values.get("row_accounting_balanced", record.get("row_accounting_balanced"))),
        ("End-to-end accounting complete", values.get("row_accounting_complete", record.get("row_accounting_complete"))),
    )
    return _section(
        "Input and harmonisation accounting",
        _key_value_table(rows),
        subtitle="Counts are copied from the persisted pre-VCF and reconciliation evidence.",
    )


def _decision_section(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    dataset = manifest.get("dataset") or {}
    decisions = dataset.get("study_decisions") or {}
    build = dataset.get("genome_build") or {}
    strand_consensus = decisions.get("strand_consensus") or {}
    rows = (
        ("Genome build", record.get("genome_build")),
        ("Genome-build decision source", record.get("genome_build_source")),
        ("Build matches / testable", "%s / %s" % (_text((build.get("matches") or {}).get(record.get("genome_build"))), _text(build.get("testable_variants")))),
        ("Effect type", record.get("effect_type")),
        ("Effect-type decision source", record.get("effect_type_source")),
        ("Effect scale", record.get("effect_scale")),
        ("SE scale", decisions.get("se_scale")),
        ("SE-scale decision source", decisions.get("se_scale_source")),
        ("P-value type", record.get("p_value_type")),
        ("P-value decision source", record.get("p_value_type_source")),
        ("Frequency type", record.get("frequency_type")),
        ("Frequency decision source", record.get("frequency_type_source")),
        ("Strand consensus", record.get("strand_consensus")),
        ("Strand consensus support", _fraction(strand_consensus.get("dominant_fraction"))),
        ("Initial MAF suspicion", decisions.get("eaf_is_maf_initial")),
        ("Final frequency interpretation", decisions.get("frequency_type")),
    )
    return _section(
        "Study-wide scientific decisions",
        _key_value_table(rows) + _mapping_details("Full recorded study-decision evidence", decisions),
        subtitle="Effect, SE, p-value and strand decisions are resolved study-wide and reused by chromosome workers. Initial MAF suspicion is provisional; the final frequency interpretation also incorporates subsequent chromosome-level reference validation.",
    )


def _strand_section(record: Mapping[str, Any]) -> str:
    rows = (
        ("Status", record.get("strand_status")),
        ("Metadata complete", record.get("strand_metadata_complete")),
        ("Reference panel", record.get("strand_reference_panel")),
        ("Reference population", record.get("strand_reference_population")),
        ("Chromosomes expected", record.get("strand_chromosomes_expected")),
        ("Chromosomes completed", record.get("strand_chromosomes_completed")),
        ("Variants evaluated", record.get("strand_variants_evaluated")),
        ("Matched", record.get("strand_matched")),
        ("Forward", record.get("strand_forward")),
        ("Forward swapped", record.get("strand_forward_swapped")),
        ("Reverse complement", record.get("strand_reverse_complement")),
        ("Reverse-complement swapped", record.get("strand_reverse_complement_swapped")),
        ("Reference unmatched detected", record.get("strand_reference_unmatched_detected")),
        ("Reference unmatched retained by opt-in policy", record.get("strand_reference_unmatched_retained")),
        ("Reference unmatched removed", record.get("strand_reference_unmatched")),
        ("Palindromic frequency conflict", record.get("strand_palindromic_frequency_conflict")),
        ("Palindromic orientation unavailable", record.get("strand_palindromic_orientation_unavailable")),
        ("Palindromic frequency discordant", record.get("strand_palindromic_frequency_discordant")),
        ("Palindromic ambiguous", record.get("strand_palindromic_ambiguous")),
        ("Reference ambiguous", record.get("strand_reference_ambiguous")),
        ("Total removed during orientation", record.get("strand_removed_total")),
        ("Orientation accounting balanced", record.get("strand_accounting_balanced")),
    )
    references = record.get("strand_reference_files")
    if isinstance(references, str):
        try:
            references = json.loads(references)
        except (ValueError, TypeError):
            references = [references]
    return _section(
        "Strand and allele orientation", _key_value_table(rows)
        + _mapping_details("All strand-reference files", {"Paths": references or []}),
        subtitle="Combined recorded chromosome evidence. Strand removals are included in the strand/EAF step total, not additional removals. Palindromic conflicts mean strand and frequency evidence supported opposite orientations.",
    )


def _vcf_section(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    assessment = manifest.get("qc_assessment") or {}
    rows = (
        ("Genome build assessed", manifest.get("qc_genome_build") or record.get("genome_build")),
        ("Final merged VCF records", record.get("final_vcf_variants")),
        ("SNPs", record.get("final_vcf_snps")),
        ("Indels / other variants", record.get("final_vcf_indels_or_other_variants")),
        ("Active QC rules", record.get("qc_active_rule_count")),
        ("All variants passing QC", record.get("qc_passed_variants")),
        ("SNPs passing QC", record.get("qc_passed_snps")),
        ("SNPs failing at least one rule", record.get("qc_failed_snps")),
        ("SNP pass percentage", record.get("qc_passed_snp_percent")),
        ("SNP fail percentage", record.get("qc_failed_snp_percent")),
        ("All variants failing at least one rule", assessment.get("excluded")),
        ("All-variant QC-pass percentage", _fraction(assessment.get("retained_fraction"))),
        ("QC definition", assessment.get("definition")),
    )
    return _section(
        "VCF content and virtual QC",
        _key_value_table(
            rows,
            percent_labels={"SNP pass percentage", "SNP fail percentage"},
        ),
        subtitle="Virtual QC only: the merged VCF remains unfiltered; no filtered VCF is created. Rule failures can overlap and must not be added together. Counts refer only to the build shown here.",
    )


def _population_comparison_table(population: Mapping[str, Any], *, compact: bool = False) -> str:
    """Share recorded population comparisons between the overview and detail."""
    rows = []
    fields = population.get("fields") or {}
    for label, values in (population.get("populations") or {}).items():
        missing = fields.get(label) or {}
        rows.append((
            label, values.get("comparable_variants"),
            values.get("pearson_correlation"), values.get("mean_absolute_difference"),
            _ratio(missing.get("missing"), missing.get("total_records"), "VCF records")
            if compact else values.get("inverted_mean_absolute_difference"),
        ))
    if not rows:
        return '<p class="section-note">Population-frequency comparisons not recorded. This is not a passing similarity result.</p>'
    return _data_table((
        "Population", "Comparable variants", "Pearson correlation",
        "Mean absolute AF difference",
        "Missing reference AF / VCF records" if compact else "Inverted mean absolute AF difference",
    ), rows)


def _frequency_section(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    population = manifest.get("population_frequency_qc") or {}
    rows = (
        ("Frequency type", record.get("frequency_type")),
        ("Most similar reference population", record.get("closest_population")),
        ("Comparison reference", "%s / %s" % (_text(record.get("comparison_af_panel")), _text(record.get("comparison_af_population")))),
        ("Comparable variants", record.get("af_comparable_variants")),
        ("Concordant variants", record.get("af_concordant_variants")),
        ("Concordance percentage", record.get("af_concordance_percent")),
        ("Mismatched variants", record.get("af_mismatched_variants")),
        ("Mismatch percentage", record.get("af_mismatch_percent")),
        ("Missing study AF", record.get("study_af_missing")),
        ("Missing reference AF", record.get("reference_af_missing")),
        ("Absolute-difference cutoff", record.get("af_difference_cutoff")),
        ("Population-frequency QC status", population.get("status")),
        ("Population-frequency decision", population.get("decision_reason")),
    )
    field_rows = []
    for label, values in (population.get("fields") or {}).items():
        field_rows.append((
            label,
            values.get("total_records"),
            values.get("missing"),
            _fraction(values.get("missing_fraction")),
            values.get("invalid"),
            _fraction(values.get("invalid_fraction")),
        ))
    field_table = ""
    if field_rows:
        headings = (
            "Field", "Total records", "Missing", "Missing percentage",
            "Invalid", "Invalid percentage",
        )
        field_table = _data_table(headings, field_rows)

    return _section(
        "Allele-frequency QC",
        _key_value_table(
            rows,
            percent_labels={"Concordance percentage", "Mismatch percentage"},
        ) + _population_comparison_table(population) + field_table,
        subtitle="Population similarity is measured on the unfiltered, input-build merged VCF. Most similar reference population is based on allele-frequency similarity among the reference populations compared; it is not confirmation of ancestry. Population correlations use the common eligible comparison set; the selected-reference AF check has its own denominator. Inverted difference is a diagnostic comparison with 1 − AF, not a correction applied by this report.",
    )


def _statistical_section(record: Mapping[str, Any]) -> str:
    rows = (
        ("Effect scale", record.get("effect_scale")),
        ("Pre-VCF invalid effect-statistic removals", record.get("pre_vcf_invalid_effect_statistic_removals")),
        ("QC-passed missing / invalid Neff", record.get("qc_passed_missing_or_invalid_neff")),
        ("QC-passed missing imputation score", record.get("qc_passed_missing_imputation_score")),
        ("Low-Neff reference percentile", _fraction(record.get("neff_reference_quantile"))),
        ("Low-Neff reference value", record.get("neff_reference_value")),
        ("Low-Neff minimum relative to reference", _fraction(record.get("neff_minimum_fraction_of_reference"))),
        ("Low-Neff threshold", record.get("neff_minimum_threshold")),
        ("Raw low-Neff variants", record.get("raw_low_neff_variants")),
        ("Raw low-Neff percentage", record.get("raw_low_neff_percent")),
        ("QC-passed low-Neff variants", record.get("qc_passed_low_neff_variants")),
        ("QC-passed low-Neff percentage", record.get("qc_passed_low_neff_percent")),
        ("QC-passed Neff upper-tail diagnostic", record.get("qc_passed_neff_upper_outliers")),
    )
    return _section(
        "Statistical quality",
        _key_value_table(
            rows,
            percent_labels={
                "Raw low-Neff percentage",
                "QC-passed low-Neff percentage",
            },
        ),
        subtitle="Neff means effective sample size. Low-Neff and upper-tail findings are diagnostics unless an active rule explicitly excludes them; passing other QC rules does not remove these concerns. Low-Neff percentages use variants with available effective sample size.",
    )


def _chromosome_sort_key(value: Any) -> tuple[int, int | str]:
    text = str(value).upper().removeprefix("CHR")
    if text.isdigit():
        return (0, int(text))
    return (1, text)


def _data_table(
    headings: Sequence[str], rows: Sequence[Sequence[Any]],
) -> str:
    """Render already-recorded tabular values without transforming them."""
    head = "".join(
        '<th scope="col">%s</th>' % escape(label) for label in headings
    )
    body = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % _cell(value) for value in row)
        for row in rows
    )
    return (
        '<div class="table-wrap subtable"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (head, body)
    )


def _chromosome_section(manifest: Mapping[str, Any]) -> str:
    dataset = manifest.get("dataset") or {}
    summaries = dataset.get("chromosome_summaries") or {}
    states = dataset.get("chromosomes") or {}
    headers = (
        "Chromosome", "Status", "Attempt", "Rows in", "Rows out", "Rejected",
        "Balanced", "Duration", "Warnings", "Errors",
    )
    rows = []
    drilldowns = []
    chromosomes = sorted(set(summaries) | set(states), key=_chromosome_sort_key)
    for chromosome in chromosomes:
        summary = summaries.get(chromosome) or {}
        state = states.get(chromosome) or {}
        values = (
            chromosome,
            summary.get("status", state.get("status")),
            summary.get("attempt", state.get("attempts")),
            summary.get("rows_in"),
            summary.get("rows_out"),
            summary.get("rejected"),
            summary.get("balanced"),
            _elapsed(summary["elapsed"]) if summary.get("elapsed") is not None else None,
            summary.get("warnings"),
            summary.get("errors"),
        )
        rows.append(values)
        cards = "".join(_stage_card(step) for step in chromosome_steps(summary, manifest.get("policies")))
        strand = ((summary.get("stage_qc") or {}).get("eaf_qc") or {}).get("strand_orientation") or {}
        orientation = _key_value_table([
            ("Unmatched retained", strand.get("reference_unmatched_retained")),
            ("Unmatched removed", strand.get("reference_unmatched")),
        ])
        drilldowns.append('<details class="chromosome-detail"><summary>Chromosome %s — %s · all 16 steps</summary>%s%s</details>' % (
            escape(str(chromosome)), escape(_text(summary.get("status", state.get("status")))), cards, orientation,
        ))
    if not chromosomes:
        drilldowns.append('<p class="empty">No chromosome outcomes were recorded. The workflow below describes the steps, not completed work.</p>')
        drilldowns.extend(_stage_card(step) for step in chromosome_steps({}, manifest.get("policies")))
    return _section(
        "Chromosome-wise processing",
        _data_table(headers, rows) + "".join(drilldowns),
        subtitle="Open a chromosome to follow its actual steps. Rows out here means harmonised rows before VCF creation, not target-build liftover output. Only recorded evidence is displayed; warnings do not necessarily mean values changed or variants were removed.",
    )


def _flatten_report_paths(prefix: str, value: Any) -> list[tuple[str, Any]]:
    if isinstance(value, Mapping):
        rows = []
        for key, nested in value.items():
            label = "%s — %s" % (prefix, str(key).replace("_", " "))
            rows.extend(_flatten_report_paths(label, nested))
        return rows
    if isinstance(value, (list, tuple)):
        return [("%s — %d" % (prefix, index), item) for index, item in enumerate(value, start=1)]
    return [(prefix, value)] if value not in (None, "") else []


def _output_paths(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> list[tuple[str, Any]]:
    """Share the recorded, retained output links between overview and detail."""
    dataset = manifest.get("dataset") or {}
    merged = manifest.get("merged_vcfs") or {}
    rows: list[tuple[str, Any]] = [
        ("Dataset manifest", record.get("manifest")),
        ("Screen report", record.get("screen_report")),
        ("Combined log", manifest.get("combined_log")),
        ("Rejected variants", dataset.get("rejected_variants_file")),
        ("Reject-reason table", dataset.get("reject_reason_table")),
        ("Pre-VCF column statistics", (manifest.get("gwas2vcf_summary_audit") or {}).get("path")),
        ("Population-frequency QC", (manifest.get("population_frequency_qc") or {}).get("report")),
    ]
    removed = manifest.get("gwas2vcf_intermediate") or {}
    for key, value in merged.items():
        if isinstance(value, str) and _href(value) and ("/" in value or "\\" in value):
            if key == "gwas2vcf" and removed.get("retained") is False and removed.get("removed_artifacts"):
                continue
            rows.append(("Merged output — " + key, value))
    rows.extend(_flatten_report_paths("VCF QC", (manifest.get("qc_assessment") or {}).get("reports") or {}))
    rows.extend(_flatten_report_paths("Concordance", (manifest.get("concordance_validation") or {}).get("reports") or {}))
    return rows


def _output_links(rows: Iterable[tuple[str, Any]], report_path: str | Path | None) -> str:
    return '<div class="table-wrap"><table><tbody>%s</tbody></table></div>' % "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>' % (
            escape(label), _path_value(value, report_path, check_exists=True),
        ) for label, value in rows
    )


def _outputs_section(record: Mapping[str, Any], manifest: Mapping[str, Any], report_path: str | Path | None) -> str:
    dataset = manifest.get("dataset") or {}
    merged = manifest.get("merged_vcfs") or {}
    removed = manifest.get("gwas2vcf_intermediate") or {}
    status_rows = (
        ("Dataset status", record.get("status") or manifest.get("status")),
        ("Completed chromosomes", dataset.get("completed")),
        ("Failed chromosomes", dataset.get("failed")),
        ("Merge status", merged.get("merge_status")),
        ("Required merge failures", merged.get("required_merge_failures")),
        ("Optional merge failures", merged.get("optional_merge_failures")),
        ("QC genome build", manifest.get("qc_genome_build")),
    )
    cleanup = ""
    if removed.get("retained") is False and removed.get("removed_artifacts"):
        cleanup = '<p>Raw GWAS-to-VCF intermediate: removed according to the recorded retention policy.</p>'
    return _section(
        "Outputs and provenance",
        _key_value_table(status_rows)
        + cleanup + _output_links(_output_paths(record, manifest), report_path)
        + _mapping_details("Software versions", manifest.get("versions") or {})
        + _mapping_details("Exact VCF provenance, input sources and transformations", manifest.get("vcf_provenance") or {})
        + _mapping_details("Effective scientific policies", manifest.get("policies") or {})
        + _mapping_details("Full resolved report configuration", (manifest.get("report_context") or {}).get("configuration") or {}),
        subtitle="Output links are relative to this report when possible. Missing files are labelled, not offered as downloads. Full reference paths are provenance; they may only be accessible on the original computer.",
    )


def _attention(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    """Surface observed findings without creating new scientific QC rules."""
    messages = []
    failure = record.get("failure_reason") or (manifest.get("failure") or {}).get("message")
    if failure:
        messages.append(("error", "Processing failure", str(failure)))
    status = _status(record.get("status") or manifest.get("status"))
    if not failure and status not in {"OK", "PARTIAL", "COMPLETED"}:
        messages.append(("warning", "Completion not established",
            f"Recorded processing status: {status}. Do not treat incomplete or unavailable results as validated outputs."))
    if record.get("effect_type_source") == "detector":
        messages.append(("warning", "Confirm the inferred effect type",
            f"Values were interpreted as {_text(record.get('effect_type'))}. Check the study documentation and explicitly declare the effect type when confirmed."))
    if record.get("p_value_type_source") == "detector":
        messages.append(("warning", "Confirm the inferred p-value type",
            f"P-values were interpreted as {_text(record.get('p_value_type'))}. Check the study documentation and set p_value_type explicitly when confirmed."))
    source = manifest.get("config") or {}
    if source.get("info_source") in {"external", "fixed_cli"}:
        messages.append(("warning", "INFO is not a study-measured score",
            f"Recorded source: {_text(source.get('info_source_detail') or source.get('info_source'))}. Passing the INFO rule does not establish measured variant-level imputation quality in this study."))
    if record.get("af_mismatched_variants"):
        messages.append(("warning", "Reference-frequency disagreements",
            f"{_ratio(record.get('af_mismatched_variants'), record.get('af_comparable_variants'), 'comparable variants')} disagreed under the recorded AF threshold. Review the reference population and AF QC below."))
    if record.get("reference_af_missing"):
        messages.append(("warning", "Reference frequencies missing",
            f"Variants lacking the selected reference AF: {_text(record['reference_af_missing'])}. Missing is different from frequency disagreement."))
    if record.get("qc_passed_low_neff_variants"):
        percent = record.get("qc_passed_low_neff_percent")
        proportion = f" ({_percentage(percent)} of those with available Neff)" if percent is not None else ""
        messages.append(("warning", "Low effective sample size remains",
            f"{_text(record['qc_passed_low_neff_variants'])} virtual-QC-passed variants{proportion} have Neff below {_text(record.get('neff_minimum_threshold'))}. This is a diagnostic, not an additional exclusion. Review suitability for the intended downstream method."))
    if str(record.get("status", "")).upper() == "PARTIAL":
        messages.append(("warning", "Incomplete dataset", "Some required work or outputs are incomplete. Interpret counts only for the recorded completed coverage."))
    for message in (manifest.get("report_context") or {}).get("sample_sheet_warnings") or []:
        messages.append(("warning", "Sample-sheet warning", str(message)))
    validation = manifest.get("concordance_validation") or {}
    if validation.get("status") and _status(validation["status"]) not in {"OK", "PASS", "PASSED"}:
        messages.append(("warning", "Review input–VCF concordance",
            f"Recorded outcome: {_text(validation['status'])}. Review the comparison below; input-only variants and statistical disagreements have different meanings."))
    if not messages:
        return '<p class="section-note">No headline alert is recorded here. Review stage evidence and chromosome warnings; this is not a statement that every scientific check passed.</p>'
    return "".join('<div class="callout %s"><strong>%s</strong><p>%s</p></div>' % (
        kind, escape(title), escape(message),
    ) for kind, title, message in messages)


def _inputs_section(manifest: Mapping[str, Any]) -> str:
    context = manifest.get("report_context") or {}
    inputs = context.get("input_mapping") or manifest.get("config") or {}
    rows = [(str(key).replace("_", " ").capitalize(), "Not provided" if value is None else value)
            for key, value in inputs.items() if not isinstance(value, (dict, list))]
    source = manifest.get("config") or {}
    return (
        _section("Study inputs and exact column mappings", _key_value_table(rows),
                 subtitle="Names are preserved exactly as supplied. An absent study column is not a failed check when a configured recovery path supplies its value later.")
        + _section("Resolved sources and initial checks", _key_value_table([
            ("Sample sheet", context.get("sample_sheet")),
            ("Resource directory", context.get("resource_directory") or source.get("resource_folder")),
            ("INFO source", source.get("info_source_detail")),
            ("User-assigned fixed INFO", source.get("fixed_info")),
        ]), subtitle="Study INFO, an external proxy and a user-assigned constant are distinct sources. Default strand/EAF validation references are also distinct from external study EAF and final VCF population annotations.")
        + _mapping_details("Exact resource checks and reference paths", (manifest.get("dataset") or {}).get("resource_preflight") or {})
    )


def _completeness_section(manifest: Mapping[str, Any]) -> str:
    evidence = (manifest.get("dataset") or {}).get("field_completeness") or {}
    conversion = evidence.get("numeric_conversion_by_column") or {}
    rows = []
    for field in evidence.get("fields") or []:
        column = field.get("input_column")
        numeric = conversion.get(column) or {}
        rows.append((field.get("field"), column, field.get("source_label"),
                     _ratio(field.get("post_parse_missing"), field.get("input_rows"), "input rows"),
                     numeric.get("invalid") if numeric else "Not assessed as numeric",
                     field.get("final_missing"), field.get("lifecycle_note")))
    return _section("Field completeness and recovery", _data_table(
        ("Field", "Study column", "Source", "Missing after parsing", "Numeric conversion failures", "Missing at final check", "Handling"), rows,
    ) + _mapping_details("Per-field recovery and final-check evidence", evidence),
        subtitle="Missing after parsing and failed numeric conversion are separate measurements. Final missing counts apply only to the assessed chromosomes and fields; absent evidence is never replaced by zero.")


def _qc_rules_section(manifest: Mapping[str, Any]) -> str:
    assessment = (manifest.get("report_context") or {}).get("qc_assessment") or {}
    blocks = []
    for rule in assessment.get("rules") or []:
        blocks.append('<div class="stage-card"><h5>%s</h5><p>%s</p>%s%s</div>' % (
            escape(str(rule.get("label") or rule.get("key"))), escape(str(rule.get("purpose") or "")),
            _key_value_table([
                ("Rule applied", rule.get("criterion")), ("Configured action", rule.get("decision")),
                ("Failing raw records", _ratio(rule.get("failed_raw"), (assessment.get("raw") or {}).get("num_records"), "raw VCF records")),
                ("Failing percentage of raw records", _fraction(rule.get("failed_fraction_raw"))),
                ("Failing this rule only", rule.get("unique_only_raw")),
                ("Also failing other rules", rule.get("overlap_raw")),
            ]), _mapping_details("Missing values and individual failure conditions", {"Conditions": rule.get("details") or []}),
        ))
    if not blocks:
        blocks.append('<p class="empty">Individual rule evidence is not recorded in this report view. Consult the linked QC report; no thresholds or passing counts are assumed.</p>')
    return _section("Active QC rules and configured actions", "".join(blocks)
        + _mapping_details("Inactive rules", {"Rules": assessment.get("inactive_rules") or []})
        + _mapping_details("Full recorded raw and QC-passed statistics", {
            "Raw merged VCF": assessment.get("raw") or {},
            "Virtual QC-passed subset": assessment.get("qc_passed") or {},
        }), subtitle="Failures can overlap. Rule-specific percentages use all raw VCF records, not the QC-passed subset. Missing-value actions are shown in each recorded criterion.")


def _concordance_section(manifest: Mapping[str, Any], *, compact: bool = False) -> str:
    context = manifest.get("report_context") or {}
    validation = manifest.get("concordance_validation") or {}
    requested = context.get("concordance_requested")
    status = validation.get("status") or ("Not requested" if requested is False else "Not recorded")
    if not validation:
        return _status_badge(status) + '<p>No input–VCF comparison result is recorded. This is not a passing concordance result.</p>'
    summary = validation.get("summary") or {}
    rows = []
    for name, values in (summary.get("variant_types") or {}).items():
        rows.append((name.replace("_", " "), values.get("input_unique_variants"), values.get("vcf_unique_variants"),
                     _ratio(values.get("exact_matched_variants"), values.get("variant_union"), "input/VCF variant union"),
                     values.get("input_only_variants"), values.get("vcf_only_variants")))
    body = _status_badge(status) + _data_table(
        ("Variant type", "Input unique", "VCF unique", "Allele-aware common", "Input only", "VCF only"), rows,
    )
    groups = [("All allele-aware matches", validation.get("metrics") or {})]
    if not compact:
        groups.extend((str(name).replace("_", " "), values) for name, values in (validation.get("metrics_by_variant_type") or {}).items())
    for title, metrics in groups:
        rows = [(str(name).replace("_", " "), _ratio(value.get("concordant"), value.get("checked"), "checked values"),
                 value.get("mismatches"), value.get("unavailable_in_input"), value.get("missing_in_vcf"), value.get("excluded_orientation"))
                for name, value in metrics.items() if isinstance(value, Mapping)]
        table = _data_table(
            ("Statistic", "Concordant / checked", "Disagreeing", "Unavailable in input", "Missing in VCF", "Orientation excluded"), rows)
        body += ('<h5>%s</h5>%s' % (escape(title), table)) if compact else _section(title, table)
    if compact:
        return body + '<p class="section-note">Values are compared only where the statistic and allele orientation are usable. Input-only or VCF-only records do not by themselves prove a statistical error. SNP and indel detail is available in section 7.</p>'
    body += '<p>Position diagnostics apply only to allele-unmatched variants. They do not prove allele equivalence: signed-effect magnitudes and folded frequencies may be compared when orientation is unknown.</p>'
    return body + _mapping_details("Variant matching, palindromic handling and position diagnostics", summary) \
        + _mapping_details("Position-only statistical diagnostics", validation.get("position_metrics_by_variant_type") or {}) \
        + _mapping_details("Configured concordance tolerances and policies", ((context.get("configuration") or {}).get("harmonisation") or {}).get("concordance_validation") or {})


def _workflow_sections(record: Mapping[str, Any], manifest: Mapping[str, Any], report_path: str | Path | None) -> list[tuple[str, str]]:
    dataset = manifest.get("dataset") or {}
    preparation = []
    for step in dataset_steps(manifest):
        preparation.append(_stage_card(step, level=4))
        if step["number"] == 3:
            preparation.append(_completeness_section(manifest))
        if step["number"] == 7:
            preparation.append(_decision_section(record, manifest))
    post = []
    for step in post_merge_steps(manifest):
        post.append(_stage_card(step, level=4))
        if step["number"] == 2:
            post.append(_frequency_section(record, manifest))
        if step["number"] == 4:
            post.extend((_vcf_section(record, manifest), _qc_rules_section(manifest), _statistical_section(record)))
    accounting = _input_section(record, manifest) + _mapping_details("Chromosome scheduling, completion and retries", {
        "Parallelism": dataset.get("parallelism") or {}, "Completed": dataset.get("completed"),
        "Failed": dataset.get("failed"), "Rounds": dataset.get("rounds"),
    }) + _mapping_details("Disjoint row accounting and recorded rejection reasons", {
        "Reconciliation": dataset.get("reconciliation") or {}, "Reject counts": dataset.get("reject_counts") or {},
    })
    return [
        ("Inputs and initial checks", _inputs_section(manifest)),
        ("Whole-study preparation", "".join(preparation)),
        ("Chromosome harmonisation", '<p>These 16 steps repeat for each chromosome. Effect-scale conversion precedes allele swapping; study-wide decisions are reused, not inferred again.</p>' + _strand_section(record) + _chromosome_section(manifest)),
        ("Chromosome completion and variant accounting", accounting),
        ("Merging and final dataset QC", "".join(post)),
        ("Input–VCF concordance", '<p>Optional comparison of the original input with the final VCF in the same genome build, after harmonisation and merged-VCF QC.</p>' + _concordance_section(manifest)),
        ("Final interpretation and downloads", '<p>Before downstream analysis, confirm the inferred study properties, chromosome coverage, reasons for removal, reference-frequency agreement and any remaining statistical warnings. A completed run is not a guarantee of suitability for every downstream method.</p>' + _outputs_section(record, manifest, report_path)),
    ]


def _decision_source(value: Any) -> str:
    """Translate recorded decision origins, never infer a missing origin."""
    return {
        "detector": "automatically inferred", "sample_sheet": "declared in sample sheet",
        "policy": "configured in YAML", "dataset_z_pvalue_cross_check": "study-wide Z/P cross-check",
        "study_level_statistic": "study-wide frequency distribution",
    }.get(value, _text(value).replace("_", " "))


def _evidence_count(group: Mapping[str, Any], key: str) -> str:
    item = (group.get("counts") or {}).get(key) or {}
    result = _text(item.get("value"))
    if item.get("value") is not None and not item.get("complete"):
        result += " — partial evidence from %s / %s chromosomes" % (
            _text(item.get("recorded_chromosomes")), _text(item.get("expected_chromosomes")),
        )
    return result


def _evidence_sources(group: Mapping[str, Any], key: str) -> str:
    values = (group.get("sources") or {}).get(key) or []
    descriptions = {
        "calculated_from_case_control": "Effective N calculated from case and control counts",
        "calculated_from_beta_se": "Calculated from harmonised BETA / SE",
        "standard_info": "Standard INFO scale",
        "warn": "warn — report disagreements and retain them at this check",
        "reject": "reject — remove affected variants at this check",
    }
    result = "; ".join(descriptions.get(value, _text(value)) for value in values) if values else "Not recorded"
    coverage = (group.get("source_coverage") or {}).get(key) or {}
    if values and coverage and not coverage.get("complete"):
        result += " — recorded for %s / %s chromosomes" % (
            _text(coverage.get("recorded_chromosomes")), _text(coverage.get("expected_chromosomes")),
        )
    return result


def _overview_inputs(record: Mapping, manifest: Mapping) -> list[tuple[str, str, int]]:
    dataset = manifest.get("dataset") or {}
    source = manifest.get("config") or {}
    pre = manifest.get("pre_vcf_summary") or {}
    chromosomes = sorted((dataset.get("chromosomes") or {}), key=_chromosome_sort_key)
    run = _key_value_table([
        ("Dataset", record.get("dataset_id") or manifest.get("sample_id")),
        ("Input file", record.get("input_file") or manifest.get("sumstat_file")),
        ("Trait type supplied", source.get("trait_type")),
        ("Recorded trait interpretation", (manifest.get("vcf_provenance") or {}).get("trait_type")),
        ("PostGWAS version", (manifest.get("vcf_provenance") or {}).get("postgwas_version")),
        ("Started", manifest.get("started")),
        ("Completed", manifest.get("completed_at")),
        ("Duration", _elapsed(manifest.get("elapsed_seconds")) if manifest.get("elapsed_seconds") is not None else None),
        ("Chromosomes included", ", ".join(chromosomes) if chromosomes else None),
        ("Chromosomes completed", _ratio(record.get("strand_chromosomes_completed"), record.get("strand_chromosomes_expected"), "expected chromosomes")),
        ("Failed chromosomes", "None" if dataset.get("failed") == [] else dataset.get("failed")),
        ("Resource validation", (dataset.get("resource_preflight") or {}).get("status")),
        ("Resource directory", source.get("resource_folder")),
        ("Dataset output directory", manifest.get("output_folder")),
    ], path_labels={"Input file", "Resource directory", "Dataset output directory"})
    input_rows = pre.get("total_variant_read", record.get("parsed_variants"))
    ready = pre.get("total_variant_remaining_for_harmonisation", record.get("input_ready_variants"))
    rows = [("Input variants", record.get("input_variants")), ("Rows read", input_rows)]
    rows.extend((label, _ratio(pre.get(key), input_rows, "rows read")) for label, key in (
        ("Missing mandatory values — removed", "total_variant_removed_missing_values"),
        ("Invalid coordinates — removed", "total_variant_removed_null_coords"),
        ("Unsupported chromosomes — removed", "total_variant_removed_unsupported_chromosomes"),
        ("Non-standard alleles — removed", "total_variant_removed_non_standard_alleles"),
        ("Duplicate variants — removed", "total_variant_removed_duplicates"),
    ))
    rows.extend([
        ("Ready for harmonisation", _ratio(ready, input_rows, "rows read")),
        ("Ready SNPs", _ratio(record.get("input_ready_snps"), ready, "ready variants")),
        ("Ready indels / other variants", _ratio(record.get("input_ready_indels_or_other_variants"), ready, "ready variants")),
    ])
    return [("Study and run", run, 2), ("Input validation", _key_value_table(rows), 3)]


def _overview_alignment(record: Mapping, manifest: Mapping, evidence: Mapping) -> list[tuple[str, str, int]]:
    dataset = manifest.get("dataset") or {}
    decisions = dataset.get("study_decisions") or {}
    build = dataset.get("genome_build") or {}
    consensus = decisions.get("strand_consensus") or {}
    source = manifest.get("config") or {}
    provenance = manifest.get("vcf_provenance") or {}
    ptype = record.get("p_value_type")
    reference_checks = decisions.get("eaf_reference_decisions")
    if (decisions.get("eaf_is_maf_initial") is False and isinstance(reference_checks, Mapping)
            and reference_checks and all(value == 0 for value in reference_checks.values())):
        reference_checks = "No second-level MAF check needed: the input was not flagged as MAF-like"
    rows = [
        ("Genome build", _text(record.get("genome_build")) + " · " + _decision_source(record.get("genome_build_source"))),
        ("Build support", _ratio((build.get("matches") or {}).get(record.get("genome_build")), build.get("testable_variants"), "build-testable variants")),
        ("Effect type", _text(record.get("effect_type")) + " · " + _decision_source(record.get("effect_type_source"))),
        ("Effect transformation", record.get("effect_scale")),
        ("SE scale", _text(decisions.get("se_scale")) + " · " + _decision_source(decisions.get("se_scale_source"))),
        ("P-value type", _text(ptype) + " · " + _decision_source(record.get("p_value_type_source"))),
        ("P-value harmonisation", {"raw": "Raw probability scale retained; range policies still apply", "neglog10": "−log10(P) converted to raw P; range policies still apply"}.get(ptype)),
        ("Final frequency interpretation", _text(record.get("frequency_type")) + " · " + _decision_source(record.get("frequency_type_source"))),
        ("Initially MAF-like", decisions.get("eaf_is_maf_initial")),
        ("Frequency reference-check outcomes", reference_checks),
        ("P-value VCF representation", provenance.get("pvalue_vcf_formula")),
    ]
    strand = [
        ("Strand reference", _text(record.get("strand_reference_panel")) + " / " + _text(record.get("strand_reference_population"))),
        ("Study strand consensus", record.get("strand_consensus")),
        ("Consensus support", _fraction(consensus.get("dominant_fraction")) + " among " + _text(consensus.get("informative")) + " informative variants"),
        ("Allele matches", _ratio(record.get("strand_matched"), record.get("strand_variants_evaluated"), "evaluated variants")),
    ]
    strand.extend((label, _ratio(record.get(key), record.get("strand_variants_evaluated"), "evaluated variants")) for label, key in (
        ("Forward / unchanged", "strand_forward"), ("Forward / allele-swapped", "strand_forward_swapped"),
        ("Reverse-complement", "strand_reverse_complement"), ("Reverse-complement + swap", "strand_reverse_complement_swapped"),
        ("Reference unmatched — retained", "strand_reference_unmatched_retained"),
        ("Reference unmatched — removed", "strand_reference_unmatched"),
        ("Palindromic consensus/AF conflict — removed", "strand_palindromic_frequency_conflict"),
        ("Palindromic ambiguous — removed", "strand_palindromic_ambiguous"),
        ("Reference ambiguous — removed", "strand_reference_ambiguous"),
        ("Total removed at strand alignment", "strand_removed_total"),
    ))
    frequency = evidence.get("frequency") or {}
    strand.extend([
        ("EAF source", _text(decisions.get("eaf_provenance")) + " · " + _text(source.get("eaf_col") or source.get("eafcolumn"))),
        ("External EAF file / template", source.get("eaffile") if source.get("eaffile") else "Not used" if decisions.get("eaf_provenance") == "study_supplied" else None),
        ("External EAF column", source.get("eafcolumn") if source.get("eaffile") else "Not used" if decisions.get("eaf_provenance") == "study_supplied" else None),
        ("Non-palindromic AF comparisons", _evidence_count(frequency, "strand_non_palindromic_af_comparable")),
        ("Non-palindromic AF disagreements", _evidence_count(frequency, "strand_non_palindromic_af_discordant")),
        ("Non-palindromic disagreement action", _evidence_sources(frequency, "strand_af_discordance_action")),
        ("Palindromic AF comparisons", _evidence_count(frequency, "strand_palindromic_af_comparable")),
        ("Palindromic AF disagreements", _evidence_count(frequency, "strand_palindromic_af_discordant")),
        ("Palindromic disagreement action", _evidence_sources(frequency, "strand_palindromic_af_discordance_action")),
    ])
    note = '<p class="section-note">Strand and AF checks overlap within one stage; do not add their removal totals. A warning does not mean AF was changed. External EAF is a reference estimate, not a measurement from the study.</p>'
    return [("Study-wide decisions", _key_value_table(rows), 3),
            ("Strand and allele-frequency harmonisation", _key_value_table(strand) + note, 4)]


def _overview_statistics(record: Mapping, manifest: Mapping, evidence: Mapping) -> str:
    source = manifest.get("config") or {}
    provenance = manifest.get("vcf_provenance") or {}
    samples = evidence.get("sample_size") or {}
    effects = evidence.get("effects") or {}
    pvalues = evidence.get("pvalues") or {}
    info = evidence.get("info") or {}
    rows = [
        ("Sample-size handling", _evidence_sources(samples, "Neff_status")),
        ("Sample-size formula", provenance.get("sample_size_formula")),
        ("Case-count source", _evidence_sources(samples, "ncase_source")),
        ("Control / total-count source", _evidence_sources(samples, "ncontrol_source")),
        ("Missing sample sizes detected", _evidence_count(samples, "variants_affected")),
        ("Missing sample-size policy", _evidence_sources(samples, "missing_action")),
        ("Missing sample sizes filled", _evidence_count(samples, "variants_imputed")),
        ("Missing sample sizes removed", _evidence_count(samples, "variants_removed")),
        ("BETA source", provenance.get("effect_source")),
        ("Effect input column", source.get("beta_or_col")),
        ("SE source", provenance.get("se_source")),
        ("SE input column", source.get("se_col")),
        ("Z source", provenance.get("z_source") or _evidence_sources(effects, "z_source")),
    ]
    rows.extend((label, _evidence_count(effects, key)) for label, key in (
        ("BETA reconstructed from Z", "beta_derived_from_z"),
        ("SE reconstructed from Z-only inputs", "se_reconstructed_from_z"),
        ("Missing SE recovered from Z", "se_recovered_from_z"),
        ("Missing SE recovered from P", "se_missing_recovered_from_pvalue"),
        ("SE using an explicit zero-P approximation", "se_from_zero_p_approximation"),
    ))
    rows.extend((label, _evidence_count(pvalues, key)) for label, key in (
        ("P-values clipped below the limit", "variants_with_pvalues_clipped_low"),
        ("P-values clipped above the limit", "variants_with_pvalues_clipped_high"),
    ))
    rows.extend([
        ("Invalid effect-statistic removals before VCF", record.get("pre_vcf_invalid_effect_statistic_removals")),
        ("INFO source", source.get("info_source_detail") or _evidence_sources(info, "source")),
        ("Fixed INFO assigned by user", "Not used" if source.get("info_source") in {"internal", "external"} else source.get("fixed_info")),
        ("INFO metric type", _evidence_sources(info, "info_score_type")),
        ("Missing INFO at stage end", _evidence_count(info, "missing_info_after")),
        ("INFO rounding corrections", _evidence_count(info, "rescaled_within_tolerance")),
        ("Out-of-range INFO rejected", _evidence_count(info, "rejected_out_of_range")),
        ("Missing INFO rejected", _evidence_count(info, "rejected_missing")),
    ])
    body = _key_value_table(rows)
    validation = evidence.get("validation") or {}
    checks = []
    for key, title in (("beta_se_z", "Z versus BETA / SE"), ("z_p", "Z versus P")):
        group = validation.get(key) or {}
        checks.append((title, *(_evidence_count(group, name) for name in ("checked", "discordant", "removed", "retained_discordant")), _evidence_sources(group, "action")))
    body += '<h5>Statistical consistency</h5>' + _data_table(
        ("Check", "Compared", "Disagreed", "Removed", "Retained by this check", "Policy"), checks,
    )
    ranges = []
    for column, bounds in (samples.get("ranges") or {}).items():
        coverage = bounds.get("coverage") or {}
        ranges.append((column, bounds.get("min"), bounds.get("max"),
                       _ratio(coverage.get("recorded_chromosomes"), coverage.get("expected_chromosomes"), "expected chromosomes")))
    body += '<h5>Recorded sample-size ranges</h5>' + (_data_table(("Column", "Minimum", "Maximum", "Coverage"), ranges)
        if ranges else '<p class="section-note">Sample-size ranges not recorded.</p>')
    body += _mapping_details("Missing sample-size plan", samples.get("missing_plan") or {})
    return body + '<p class="section-note">Counts refer to their own stages. Supplied and reconstructed statistics, external INFO proxies and fixed constants are different sources. A zero is displayed only when supported by recorded evidence; partial coverage is labelled.</p>'


def _overview_export(record: Mapping, manifest: Mapping, evidence: Mapping) -> str:
    dataset = manifest.get("dataset") or {}
    pre = manifest.get("pre_vcf_summary") or {}
    merged = manifest.get("merged_vcfs") or {}
    vcf = evidence.get("vcf") or {}
    configuration = ((manifest.get("report_context") or {}).get("configuration") or {}).get("harmonisation") or {}
    target_build = (configuration.get("vcf_processing") or {}).get("target_builds", {}).get(record.get("genome_build"))
    rows = [
        ("Exported to GWAS-to-VCF", pre.get("total_variant_in_vcf_input", record.get("harmonised_variants"))),
        ("Pre-VCF rejections", _ratio(record.get("rejected_variants"), record.get("parsed_variants", record.get("input_variants")), "rows read")),
        ("Unprocessed rows from failed chromosomes", record.get("unprocessed_failed_chromosome_variants")),
        ("Accounting balanced", pre.get("row_accounting_balanced", record.get("row_accounting_balanced"))),
        ("Accounting complete", pre.get("row_accounting_complete", record.get("row_accounting_complete"))),
        ("Assessed merged VCF — " + _text(manifest.get("qc_genome_build") or record.get("genome_build")), record.get("final_vcf_variants")),
        ("SNPs in assessed merged VCF", _ratio(record.get("final_vcf_snps"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("Indels / other in assessed merged VCF", _ratio(record.get("final_vcf_indels_or_other_variants"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("Target genome build", target_build),
    ]
    rows.extend((label, _evidence_count(vcf, key)) for label, key in (
        ("Adapter output — recorded chromosome total", "adapter_output"),
        ("After normalization — recorded chromosome total", "normalization_output"),
        ("Entering liftover", "liftover_input"), ("Liftover plugin rejections", "liftover_rejected"),
        ("Successfully lifted but excluded by swap policy", "liftover_swap_excluded"),
        ("Target output — recorded chromosome total", "target_output"),
    ))
    rows.extend([
        ("Liftover swap policy", _evidence_sources(vcf, "swap_policy")),
        ("Merge outcome", merged.get("merge_status")),
        ("Required merge failures", "None" if merged.get("required_merge_failures") == [] else merged.get("required_merge_failures")),
    ])
    reasons = dataset.get("reject_counts") or {}
    losses = _data_table(("Recorded pre-VCF rejection reason", "Rejected / rows read"), [
        (str(reason).replace("_", " "), _ratio(count, record.get("parsed_variants", record.get("input_variants")), "rows read"))
        for reason, count in reasons.items() if isinstance(count, (int, float)) and count > 0
    ]) if reasons else '<p class="section-note">Detailed rejection counts not recorded.</p>'
    return _key_value_table(rows) + losses + '<p class="section-note">Target counts above are sums of recorded chromosome liftover outcomes, not a new count of the merged target VCF. Successfully lifted variants excluded by policy are not liftover failures. Pre-VCF rejections, liftover exclusions and virtual-QC failures are separate; do not add them together.</p>'


def _overview_final_qc(record: Mapping, manifest: Mapping) -> list[tuple[str, str, int]]:
    population = manifest.get("population_frequency_qc") or {}
    selected = (population.get("populations") or {}).get(record.get("closest_population")) or {}
    assessment = manifest.get("qc_assessment") or {}
    detail = (manifest.get("report_context") or {}).get("qc_assessment") or {}
    raw = detail.get("raw") or {}
    passed = detail.get("qc_passed") or {}
    frequencies = _key_value_table([
        ("Assessment genome build", manifest.get("qc_genome_build") or record.get("genome_build")),
        ("Comparison panel / population", _text(record.get("comparison_af_panel")) + " / " + _text(record.get("comparison_af_population"))),
        ("AF agreement criterion", "Absolute AF difference ≤ " + _text(record.get("af_difference_cutoff"))),
        ("Concordant AF", _ratio(record.get("af_concordant_variants"), record.get("af_comparable_variants"), "comparable variants")),
        ("Discordant AF", _ratio(record.get("af_mismatched_variants"), record.get("af_comparable_variants"), "comparable variants")),
        ("Missing study AF", _ratio(record.get("study_af_missing"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("Missing reference AF", _ratio(record.get("reference_af_missing"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("Most similar reference population", record.get("closest_population")),
        ("Reference-population AF correlation", selected.get("pearson_correlation")),
        ("Population-correlation comparable variants", selected.get("comparable_variants", population.get("comparable_variants"))),
        ("Population comparison status", population.get("status")),
        ("Population comparison decision", population.get("decision_reason")),
        ("Population comparison VCF", population.get("raw_merged_vcf")),
    ], path_labels={"Population comparison VCF"})
    frequencies += '<h5>GWAS AF versus every reference population</h5>' + _population_comparison_table(population, compact=True)
    frequencies += '<p class="section-note">Population correlations use the common eligible records from the unfiltered input-build merged VCF; missingness uses all its records. Higher correlation and smaller mean AF difference indicate greater similarity, not confirmed ancestry. If study AF came from a reference panel, similarity to that panel is not independent evidence. The selected-reference AF agreement check above has its own denominator.</p>'
    qc = _key_value_table([
        ("Genome build assessed", manifest.get("qc_genome_build") or record.get("genome_build")),
        ("Active QC rules", record.get("qc_active_rule_count")),
        ("All variants passing", _ratio(record.get("qc_passed_variants"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("All variants failing ≥1 rule", _ratio(assessment.get("excluded"), record.get("final_vcf_variants"), "assessed VCF records")),
        ("SNPs passing", _ratio(record.get("qc_passed_snps"), record.get("final_vcf_snps"), "SNPs in the assessed VCF")),
        ("SNPs failing ≥1 rule", _ratio(record.get("qc_failed_snps"), record.get("final_vcf_snps"), "SNPs in the assessed VCF")),
        ("Indels / other passing", _ratio(passed.get("num_non_snps"), raw.get("num_non_snps"), "non-SNP records")),
        ("QC-passed missing / invalid Neff", record.get("qc_passed_missing_or_invalid_neff")),
        ("QC-passed missing INFO", record.get("qc_passed_missing_imputation_score")),
        ("QC-passed low Neff", _ratio(record.get("qc_passed_low_neff_variants"), passed.get("effective_sample_size_available"), "QC-passed variants with available Neff")),
        ("Low-Neff threshold", record.get("neff_minimum_threshold")),
        ("QC-passed Neff upper outliers", record.get("qc_passed_neff_upper_outliers")),
    ])
    rules = [(rule.get("label"), _ratio(rule.get("failed_raw"), record.get("final_vcf_variants"), "assessed VCF records"), rule.get("criterion")) for rule in detail.get("rules") or []]
    qc += _data_table(("Active rule", "Failing records", "Criterion"), rules) if rules else '<p class="section-note">Individual QC-rule outcomes not recorded.</p>'
    qc += '<p class="important-note"><strong>Virtual QC leaves the merged VCF unfiltered.</strong> It assesses a subset; it does not create a QC-filtered VCF. Failures of different rules may overlap. Low-Neff and upper-tail findings are diagnostics unless an active rule explicitly excludes them. No target-build QC pass is implied by a source-build assessment.</p>'
    return [("Reference-frequency QC", frequencies, 6), ("Final VCF QC", qc, 6)]


def _results_overview(record: Mapping, manifest: Mapping, report_path: str | Path | None) -> str:
    evidence = summarise_chromosome_evidence(manifest)
    blocks = _overview_inputs(record, manifest) + _overview_alignment(record, manifest, evidence)
    blocks.extend([
        ("Sample size and statistical handling", _overview_statistics(record, manifest, evidence), 4),
        ("Export, VCF creation, and liftover", _overview_export(record, manifest, evidence), 5),
    ])
    blocks.extend(_overview_final_qc(record, manifest))
    blocks.append(("Input–VCF concordance", _concordance_section(manifest, compact=True), 7))
    links = _output_paths(record, manifest)
    downloads = [(label, path) for label, path in links if label in {"Dataset manifest", "Rejected variants", "Reject-reason table", "VCF QC — html"} or label.startswith("Merged output — ")]
    blocks.append(("Final takeaway", _attention(record, manifest) + '<h5>Outputs and supporting evidence</h5>' + _output_links(downloads, report_path), 8))
    prefix = "dataset-" + _safe_id(record.get("dataset_id") or manifest.get("sample_id") or "Dataset")
    cards = _metric_cards((
        ("Processing status", record.get("status") or manifest.get("status"), "status"),
        ("Input variants", record.get("input_variants"), "number"),
        ("Pre-VCF rejected", record.get("rejected_variants"), "number"),
        ("VCF — " + _text(manifest.get("qc_genome_build") or record.get("genome_build")), record.get("final_vcf_variants"), "number"),
        ("Virtual QC passed", record.get("qc_passed_variants"), "number"),
        ("Virtual QC failed", (manifest.get("qc_assessment") or {}).get("excluded"), "number"),
    ))
    cards += _metric_cards((
        ("Inferred genome build", record.get("genome_build_detected"), "text"),
        ("Detected effect type", record.get("effect_type_detected"), "text"),
        ("Detected P-value type", record.get("p_value_type_detected"), "text"),
        ("Most similar reference population", record.get("closest_population"), "text"),
    ))
    navigation = '<nav class="overview-nav" aria-label="Results summary"><ol>%s</ol></nav>' % "".join(
        '<li><a href="#%s-overview-%d">%s</a></li>' % (prefix, number, escape(title))
        for number, (title, _, _) in enumerate(blocks, 1)
    )
    panels = "".join('<section class="overview-panel" id="%s-overview-%d"><header><h4>1.%d %s</h4><a href="#%s-section-%d">Detailed evidence →</a></header>%s</section>' % (
        prefix, number, number, escape(title), prefix, target, body,
    ) for number, (title, body, target) in enumerate(blocks, 1))
    return cards + '<p class="important-note">Processing completion and scientific QC are different outcomes. Detected values are recorded inference results, not proof of the study semantics; the settings actually used and their sources appear in 1.3. Counts describe the named stage and build. Missing evidence is not a zero or a passing result.</p>' + navigation + '<div class="overview-panels">' + panels + '</div>'


def _dataset_article(
    record: Mapping[str, Any],
    manifest: Mapping[str, Any],
    *,
    heading_level: int,
    include_report_link: bool,
    report_path: str | Path | None = None,
) -> str:
    dataset_id = record.get("dataset_id") or manifest.get("sample_id") or "Dataset"
    title_tag = "h%d" % heading_level
    report_link = ""
    if include_report_link and record.get("html_report"):
        report_link = '<a class="button" href="%s">Open dataset-only report</a>' % escape(
            _href(record.get("html_report"), report_path) or "#", quote=True,
        )
    contents = [("Results at a glance", _results_overview(record, manifest, report_path))]
    contents.extend(_workflow_sections(record, manifest, report_path))
    prefix = "dataset-" + _safe_id(dataset_id)
    navigation = '<nav class="section-nav" aria-label="Report sections"><ol>%s</ol></nav>' % "".join(
        '<li><a href="#%s-section-%d">%s</a></li>' % (prefix, index, escape(title))
        for index, (title, _body) in enumerate(contents, 1))
    sections = "".join('<section class="report-section" id="%s-section-%d"><h3>%d. %s</h3>%s</section>' % (
        prefix, index, index, escape(title), body,
    ) for index, (title, body) in enumerate(contents, 1))
    return (
        '<article id="dataset-%s" class="dataset"><header class="dataset-header">'
        '<div><%s>%s</%s><p>Harmonisation dataset report</p></div>%s</header>'
        '%s%s</article>'
        % (
            _safe_id(dataset_id), title_tag, escape(str(dataset_id)), title_tag,
            report_link, navigation, sections,
        )
    )


_STYLES = """
:root{--ink:#172b3a;--muted:#536777;--line:#d8e2e8;--panel:#fff;--canvas:#f3f6f8;--brand:#165c64;--ok:#166534;--warn:#854d0e;--error:#991b1b}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:20px}
body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.65 system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:#165b8c;text-underline-offset:3px}a:hover{text-decoration-thickness:2px}a:focus-visible,summary:focus-visible,input:focus-visible{outline:3px solid #126f86;outline-offset:4px}
.page-header{padding:32px max(24px,calc((100vw - 1280px)/2));color:#fff;background:#173e4c;border-bottom:5px solid #6fc6bd}
.eyebrow{margin:0 0 8px;font-size:12px;font-weight:700;letter-spacing:.16em;text-transform:uppercase}
.page-header h1{margin:0;font-size:clamp(24px,3.4vw,40px);line-height:1.2;overflow-wrap:anywhere}.page-header .lede{max-width:900px;margin:12px 0 0;color:#e0edf0}
.container{max-width:1280px;margin:auto;padding:24px 24px 56px}.provenance-note{padding:12px 16px;border:1px solid var(--line);border-left:4px solid var(--brand);background:#fff;border-radius:6px;font-size:13px}
.run-nav,.overview,.dataset{background:var(--panel);border:1px solid var(--line);border-radius:12px}
.run-nav,.overview{padding:20px;margin:20px 0}.run-nav ul{display:flex;flex-wrap:wrap;gap:10px;padding:0;list-style:none}.run-nav label{display:block;font-weight:700}
input[type=search]{width:min(100%,540px);padding:10px 12px;border:1px solid #819aaa;border-radius:6px;font:inherit}
.button,.run-nav a{display:inline-block;padding:7px 12px;background:#e8f4f3;border-radius:6px;color:#164d52;text-decoration:none}
.dataset{margin:20px 0;padding:24px}.dataset-header{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;border-bottom:1px solid var(--line);padding-bottom:18px}
.dataset-header h2{margin:0;font-size:24px;overflow-wrap:anywhere}.dataset-header p{margin:5px 0;color:var(--muted)}
.section-nav{padding:18px 0;border-bottom:1px solid var(--line)}.section-nav ol{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:7px 30px;margin:0;padding-left:22px}
.report-section{margin:32px 0 0;padding-top:8px}.report-section>h3{font-size:23px;border-bottom:2px solid #a6cfce;padding-bottom:10px;margin:0 0 18px}
.overview-nav{padding:16px 20px;background:#edf5f5;border:1px solid #c7dfdf;border-radius:8px;margin:20px 0}
.overview-nav ol{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:7px 32px;margin:0;padding-left:22px}
.overview-panel{border:1px solid var(--line);border-radius:9px;margin:20px 0;padding:20px;background:#fff}
.overview-panel>header{display:flex;flex-wrap:wrap;align-items:baseline;justify-content:space-between;gap:8px 20px;border-bottom:1px solid var(--line);padding-bottom:12px;margin-bottom:16px}
.overview-panel>header h4{margin:0}.overview-panel>header a{font-size:13px}.overview-panel h5{margin-top:20px}
.important-note{padding:14px 18px;background:#eef5f9;border-left:4px solid #317798;border-radius:5px;color:#254d63}
h4{font-size:18px;margin:18px 0 8px}h5{font-size:16px;margin:0 0 8px}.section-note,.muted-text,.empty{color:var(--muted)}
.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:20px 0}
.metric{padding:14px;border:1px solid var(--line);border-top:3px solid var(--brand);border-radius:7px;background:#fbfdfd;min-width:0}
.metric-label{font-size:12px;font-weight:700;color:var(--muted)}.metric-value{font-size:24px;font-weight:750;overflow-wrap:anywhere;font-variant-numeric:tabular-nums}
.badge{display:inline-block;padding:3px 8px;border-radius:5px;font-size:11px;font-weight:750;letter-spacing:.02em;vertical-align:middle}
.badge.ok{color:var(--ok);background:#e2f4e8}.badge.warning{color:var(--warn);background:#fef3c7}.badge.error{color:var(--error);background:#fee2e2}.badge.running{color:#075985;background:#e0f2fe}.badge.muted{color:#475569;background:#edf1f4}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:7px;max-width:100%;margin:10px 0}.subtable{margin-top:14px}
table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:9px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;overflow-wrap:anywhere;font-variant-numeric:tabular-nums}
thead th{background:#edf3f6;font-size:12px;font-weight:750;position:sticky;top:0}tbody tr:last-child th,tbody tr:last-child td{border-bottom:0}.key-value th{width:34%;color:#445c6c;background:#f8fafb;font-weight:600}
.stage-card{border:1px solid var(--line);border-left:3px solid #5a9f9e;border-radius:8px;margin:16px 0;padding:16px;background:#fff}.stage-card p{margin:8px 0}
.stage-card .badge{margin-left:8px}.path{overflow-wrap:anywhere;word-break:break-word}.file-path{font-size:12px}
.callout{padding:14px 16px;border-radius:7px;margin:12px 0;border:1px solid}.callout p{margin:5px 0 0}.callout.error{color:#7f1d1d;background:#fef2f2;border-color:#fecaca}.callout.warning{color:#713f12;background:#fffbeb;border-color:#f0d797}
details{margin:12px 0;border:1px solid var(--line);border-radius:7px;background:#fbfcfd}summary{cursor:pointer;padding:10px 12px;font-weight:650}details[open]>summary{border-bottom:1px solid var(--line)}
details>.table-wrap,details>.stage-card,details>details,details>p{margin:12px}details>.dataset{margin:0;border:0}.dataset-group>summary{font-size:18px;padding:16px}.chromosome-detail>summary{background:#edf5f5}.file-path{border:0;background:none;margin:3px 0}.file-path summary{padding:0;font-weight:400}
.footer{max-width:1280px;margin:auto;padding:0 24px 32px;color:var(--muted);font-size:13px}
[hidden]{display:none!important}
@media(max-width:700px){.overview-nav ol{grid-template-columns:1fr}.overview-panel{padding:12px}.overview-panel .key-value th{width:42%}}
@media(max-width:720px){.container{padding:12px 8px 32px}.dataset{padding:14px}.dataset-header{display:block}.section-nav ol{grid-template-columns:1fr}.page-header{padding:24px 18px}.key-value th{width:40%}.stage-card{padding:10px}.report-section>h3{font-size:20px}th,td{padding:8px}}
@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0;border:0}.page-header .lede{color:var(--muted)}.container{max-width:none;padding:0}.section-nav,.run-nav,.button,input{display:none}.dataset{border:0;padding:0}.table-wrap{overflow:visible}.stage-card{break-inside:avoid}details::details-content{display:block;content-visibility:visible}a{color:inherit}.report-section{break-before:auto}thead{display:table-header-group}}
"""


def _document(title: str, header: str, body: str) -> str:
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<style>%s</style>
</head>
<body>
%s
<main class="container">
<div class="provenance-note"><strong>Evidence reuse:</strong> this report is rendered from the existing harmonisation run summary, dataset manifest, and QC report references. It does not reread summary statistics or VCFs and does not calculate new scientific metrics.</div>
%s
</main>
<footer class="footer">PostGWAS harmonisation report · review the linked canonical manifest and QC tables for machine-readable provenance.</footer>
</body>
</html>
""" % (escape(title), _STYLES, header, body)


def render_dataset_report(
    record: Mapping[str, Any], manifest: Mapping[str, Any],
    *, report_path: str | Path | None = None,
) -> str:
    """Render one complete dataset report from already-persisted evidence."""
    dataset_id = record.get("dataset_id") or manifest.get("sample_id") or "Dataset"
    header = (
        '<header class="page-header"><p class="eyebrow">PostGWAS</p>'
        '<h1>%s</h1><p class="lede">Harmonisation input, scientific decisions, '
        'chromosome processing, VCF content, QC, and output provenance.</p></header>'
        % escape(str(dataset_id))
    )
    return _document(
        "%s — PostGWAS harmonisation report" % dataset_id,
        header,
        _dataset_article(
            record, manifest, heading_level=2, include_report_link=False,
            report_path=report_path or record.get("html_report"),
        ),
    )


def _run_overview(records: Sequence[Mapping[str, Any]], report_path: str | Path | None = None) -> str:
    headers = (
        "Dataset", "Status", "Input variants", "Build", "Effect / P-value",
        "Frequency", "Final VCF", "SNPs passing virtual QC",
        "Dataset report",
    )
    rows = []
    for record in records:
        report_link = record.get("html_report")
        if report_link:
            rendered_link = '<a href="%s">Open report</a>' % escape(
                _href(report_link, report_path) or "#", quote=True,
            )
        else:
            rendered_link = '<span class="muted-text">Not available</span>'
        values = (
            record.get("dataset_id"),
            None,
            record.get("input_variants"),
            record.get("genome_build"),
            _text(record.get("effect_type")) + " / " + _text(record.get("p_value_type")),
            record.get("frequency_type"),
            record.get("final_vcf_variants"),
            _ratio(record.get("qc_passed_snps"), record.get("final_vcf_snps"), "VCF SNPs"),
            None,
        )
        cells = []
        for index, value in enumerate(values):
            if index == 1:
                rendered = _status_badge(record.get("status"))
            elif index == len(values) - 1:
                rendered = rendered_link
            else:
                rendered = _cell(value)
            cells.append("<td>%s</td>" % rendered)
        rows.append('<tr data-search="%s">%s</tr>' % (
            escape(str(record.get("dataset_id")) + " " + _status(record.get("status")), quote=True), "".join(cells),
        ))
    head = "".join("<th scope=\"col\">%s</th>" % escape(label) for label in headers)
    return '<section class="overview"><h2>All selected datasets</h2><div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div></section>' % (head, "".join(rows))


def render_run_report(
    records: Sequence[Mapping[str, Any]],
    manifests: Mapping[str, Mapping[str, Any]],
    *, report_path: str | Path | None = None,
) -> str:
    """Render one multi-dataset report containing every selected dataset."""
    navigation = '<nav class="run-nav" aria-label="Dataset reports"><label for="dataset-search">Find a dataset by name or status</label><input type="search" id="dataset-search" placeholder="Dataset name or status"><p>All selected datasets are included below. Open any dataset for its complete report.</p><ul>%s</ul></nav>' % "".join(
        '<li><a href="#dataset-%s">%s</a></li>'
        % (_safe_id(record.get("dataset_id")), escape(str(record.get("dataset_id"))))
        for record in records
    )
    articles = "".join(
        '<details class="dataset-group" data-search="%s"><summary>%s · %s</summary>%s</details>' % (
            escape(str(record.get("dataset_id")) + " " + _status(record.get("status")), quote=True),
            escape(str(record.get("dataset_id"))), _status_badge(record.get("status")), _dataset_article(
            record,
            manifests.get(str(record.get("dataset_id"))) or {},
            heading_level=2,
            include_report_link=True,
            report_path=report_path,
        ))
        for record in records
    )
    header = (
        '<header class="page-header"><p class="eyebrow">PostGWAS</p>'
        '<h1>Harmonisation run report</h1><p class="lede">%d selected dataset%s '
        'with input statistics, study-wide decisions, chromosome evidence, VCF '
        'statistics, QC, and output provenance.</p></header>'
        % (len(records), "" if len(records) == 1 else "s")
    )
    return _document(
        "PostGWAS harmonisation run report",
        header,
        navigation + _run_overview(records, report_path) + articles + '''<script type="text/javascript">
const search = document.getElementById('dataset-search');
search.addEventListener('input', () => {
  const query = search.value.toLocaleLowerCase();
  document.querySelectorAll('[data-search]').forEach(item => {
    item.hidden = !item.dataset.search.toLocaleLowerCase().includes(query);
  });
});
document.querySelectorAll('.run-nav a').forEach(link => link.addEventListener('click', () => {
  const target = document.getElementById(link.hash.slice(1));
  if (target) { const group = target.closest('details'); group.open = true; group.hidden = false; }
}));
</script>''',
    )


def write_dataset_report(
    path: str | Path,
    record: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> Path:
    """Atomically write one dataset HTML report."""
    return write_html_report(render_dataset_report(record, manifest, report_path=path), path)


def write_run_report(
    path: str | Path,
    records: Sequence[Mapping[str, Any]],
    manifests: Mapping[str, Mapping[str, Any]],
) -> Path:
    """Atomically write the combined HTML report in sample-sheet order."""
    return write_html_report(render_run_report(records, manifests, report_path=path), path)


__all__ = [
    "render_dataset_report",
    "render_run_report",
    "write_dataset_report",
    "write_run_report",
]
