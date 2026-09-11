"""Machine-readable and human-readable reports for validated filtering runs.

The report builders are presentation-only. They consume the reconciled counts
produced by the filtering engine and never reopen or reinterpret a VCF.
"""

from __future__ import annotations

from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from rich.cells import cell_len

from postgwas.core.io.reports import write_delimited_report, write_html_report
from postgwas.core.validation_reporting import (
    consolidate_validation_fields, flush_file_validation_display,
)
from postgwas.core.ui.screen import (
    SYMBOLS,
    render_action_plan,
    screen_field,
    screen_line,
    screen_section_lines,
)


SUMMARY_COLUMNS = (
    "record_type",
    "dataset_id",
    "genome_build",
    "input_vcf",
    "input_validation_status",
    "declared_contigs",
    "required_declared_fields_status",
    "required_declared_fields_present",
    "required_declared_fields_total",
    "order",
    "filter_id",
    "category",
    "reason",
    "action",
    "variants_before",
    "variants_after",
    "variants_removed",
    "variants_retained_percent",
    "soft_filter_variants",
    "variants_matching_reason",
    "variants_removed_for_this_reason",
    "variants_failing_multiple_rules",
    "overlapping_rule_matches",
    "reconciliation_status",
    "requested_threads",
    "polars_thread_pool_size",
    "thread_budget_enforced",
    "bcftools_expression",
    "filtered_vcf",
    "soft_filtered_vcf",
    "reason_summary",
    "summary_csv",
    "html_report",
    "filter_log",
)


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        return value != value
    except Exception:
        return True


def _count(value: Any, missing: str = "unavailable") -> str:
    if _is_missing(value):
        return missing
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def _rule_removes(rule: Mapping[str, Any]) -> bool:
    """Return whether a rule participates in hard-removal attribution."""
    return bool(rule.get("removes", rule.get("action") == "remove"))


def _display_rule_action(rule: Mapping[str, Any]) -> str:
    """Return the user-facing action without changing the internal contract."""
    return "EXCLUDE" if _rule_removes(rule) else "KEEP"


def _display_rule_label(rule: Mapping[str, Any]) -> str:
    """Return the prepared rule label or derive the missing-field label."""
    label = rule.get("label")
    if label:
        return str(label)
    field = rule.get("field")
    return "Missing %s" % field if field else "Unnamed filtering rule"


def _display_rule_group(rule: Mapping[str, Any]) -> tuple[str | None, str]:
    """Return an optional rule-family heading and a validated screen kind."""
    group = rule.get("display_group")
    kind = str(rule.get("display_group_kind") or "analysis")
    return (
        str(group) if group else None,
        kind if kind in SYMBOLS else "analysis",
    )


def _earlier_rule_overlap(rule: Mapping[str, Any]) -> int | None:
    """Return failures already assigned to an earlier removal rule."""
    if (
        not _rule_removes(rule)
        or _is_missing(rule.get("count"))
        or _is_missing(rule.get("primary_count"))
    ):
        return None
    return int(rule["count"]) - int(rule["primary_count"])


def _retained_percent(before: Any, after: Any) -> float | None:
    if _is_missing(before) or _is_missing(after):
        return None
    before = int(before)
    return 100.0 if before == 0 and int(after) == 0 else (
        None if before == 0 else 100.0 * int(after) / before
    )


def _required_field_counts(
    validation: Mapping[str, Any],
) -> tuple[int, int]:
    fields = validation.get("required_fields") or ()
    return (
        sum(bool(field.get("declared")) for field in fields),
        len(fields),
    )


def _required_field_summary(validation: Mapping[str, Any]) -> str:
    if not validation.get("required_fields_evaluated", True):
        return "not evaluated"
    present, total = _required_field_counts(validation)
    if total == 0:
        return "none required by the active filters"
    if present == total:
        return "all %d present" % total
    missing = total - present
    return "%d of %d present; %d missing" % (present, total, missing)


def _required_field_status(validation: Mapping[str, Any]) -> str:
    if not validation.get("required_fields_evaluated", True):
        return "NOT_EVALUATED"
    present, total = _required_field_counts(validation)
    return "PASS" if present == total else "FAILED"


def _header_contract_summary(validation: Mapping[str, Any]) -> str:
    """Summarize the complete input-header decision without hiding field counts."""
    status = str(validation.get("status") or "UNKNOWN").upper()
    decision = "PASSED" if status == "PASS" else status
    if not validation.get("required_fields_evaluated", True):
        return "%s · required INFO/FORMAT fields not evaluated" % decision
    present, total = _required_field_counts(validation)
    if total == 0:
        return "%s · no INFO/FORMAT fields required by active rules" % decision
    return "%s · %d/%d required INFO/FORMAT fields declared" % (
        decision, present, total,
    )


def build_filtering_summary(
    *,
    dataset_id: str,
    genome_build: str,
    missing_checks: Sequence[Mapping[str, Any]],
    condition_checks: Sequence[Mapping[str, Any]],
    reason_statistics: Mapping[str, Any],
    variants_before: int | None,
    variants_after: int | None,
    soft_filter_variants: int | None,
    soft_filter_enabled: bool,
    input_validation: Mapping[str, Any] | None,
    data_flow: Mapping[str, Any] | None,
    outputs: Mapping[str, Any],
    resolved_configuration: Mapping[str, Any],
    runtime_seconds: float,
    display_missing_counts: bool,
) -> dict[str, Any]:
    """Build the single evidence object used by every filtering presentation."""
    removed = (
        None
        if variants_before is None or variants_after is None
        else int(variants_before) - int(variants_after)
    )
    attributed = reason_statistics.get("primary_removed_total")
    reconciled = bool(
        reason_statistics.get("reconciled")
        and removed is not None
        and attributed is not None
        and int(attributed) == removed
    )
    recorded_rules = reason_statistics.get("checks") or ()
    if recorded_rules:
        rules = [dict(rule) for rule in recorded_rules]
    else:
        rules = [
            dict(rule, category="missing_value") for rule in missing_checks
        ] + [
            dict(rule, category="filter_condition", action="remove")
            for rule in condition_checks
        ]
    rules.sort(key=lambda item: item.get("priority", 999))
    thread_runtime = dict(reason_statistics.get("aggregation") or {})
    validation = dict(input_validation or {})
    validation["required_fields"] = [
        dict(field) for field in validation.get("required_fields") or ()
    ]
    return {
        "dataset_id": str(dataset_id),
        "genome_build": str(genome_build),
        "attribution_method": (
            "A variant can fail more than one active rule. Each removed variant "
            "is assigned to the first failed rule in the recorded audit order so "
            "primary-removal counts remain mutually exclusive."
        ),
        "variants_before": variants_before,
        "variants_after": variants_after,
        "variants_removed": removed,
        "variants_retained_percent": _retained_percent(
            variants_before, variants_after,
        ),
        "soft_filter_enabled": bool(soft_filter_enabled),
        "soft_filter_variants": soft_filter_variants,
        "input_validation": validation,
        "primary_removed_total": attributed,
        "overlap_variants": reason_statistics.get("overlap_variants"),
        "extra_rule_matches": reason_statistics.get("extra_rule_matches"),
        "reconciled": reconciled,
        "reconciliation_status": "PASS" if reconciled else "FAILED",
        "thread_runtime": thread_runtime,
        "display_missing_counts": bool(display_missing_counts),
        "rules": rules,
        "data_flow": dict(data_flow or {}),
        "outputs": dict(outputs),
        "resolved_configuration": dict(resolved_configuration),
        "runtime_seconds": float(runtime_seconds),
    }


def filtering_summary_records(summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return one overall row plus one row per active rule for the CSV."""
    outputs = summary["outputs"]
    common = {
        "dataset_id": summary["dataset_id"],
        "genome_build": summary["genome_build"],
    }
    thread_runtime = summary.get("thread_runtime") or {}
    validation = summary.get("input_validation") or {}
    fields_present, fields_total = _required_field_counts(validation)
    records = [{
        **common,
        "record_type": "overall",
        "input_vcf": validation.get("input_vcf"),
        "input_validation_status": validation.get("status"),
        "declared_contigs": validation.get("declared_contigs"),
        "required_declared_fields_status": _required_field_status(validation),
        "required_declared_fields_present": fields_present,
        "required_declared_fields_total": fields_total,
        "variants_before": summary.get("variants_before"),
        "variants_after": summary.get("variants_after"),
        "variants_removed": summary.get("variants_removed"),
        "variants_retained_percent": (
            None
            if summary.get("variants_retained_percent") is None
            else "%.6f" % float(summary["variants_retained_percent"])
        ),
        "soft_filter_variants": summary.get("soft_filter_variants"),
        "variants_failing_multiple_rules": summary.get("overlap_variants"),
        "overlapping_rule_matches": summary.get("extra_rule_matches"),
        "reconciliation_status": summary.get("reconciliation_status"),
        "requested_threads": thread_runtime.get("requested_threads"),
        "polars_thread_pool_size": thread_runtime.get(
            "polars_thread_pool_size"
        ),
        "thread_budget_enforced": thread_runtime.get(
            "thread_budget_enforced"
        ),
        "filtered_vcf": outputs.get("filtered_vcf"),
        "soft_filtered_vcf": outputs.get("soft_filtered_vcf"),
        "reason_summary": outputs.get("reason_summary"),
        "summary_csv": outputs.get("summary_csv"),
        "html_report": outputs.get("html_report"),
        "filter_log": outputs.get("filter_log"),
    }]
    for order, rule in enumerate(summary.get("rules") or (), 1):
        removes = _rule_removes(rule)
        records.append({
            **common,
            "record_type": "rule",
            "order": order,
            "filter_id": rule.get("tag") if removes else None,
            "category": rule.get("category"),
            "reason": rule.get("label"),
            "action": rule.get("action", "remove"),
            "variants_matching_reason": rule.get("count"),
            "variants_removed_for_this_reason": rule.get("primary_count"),
            "bcftools_expression": rule.get("expr"),
        })
    return records


def write_filtering_summary_csv(
    summary: Mapping[str, Any], destination: str | Path,
) -> Path:
    """Write the reconciled overall and rule-level summary atomically."""
    return write_delimited_report(
        filtering_summary_records(summary),
        destination,
        fieldnames=SUMMARY_COLUMNS,
        delimiter=",",
        null_value="",
    )


def render_input_vcf_validation(
    validation: Mapping[str, Any], *, label_width: int,
) -> str:
    """Render input-contract evidence for successful and failed preflight."""
    if not validation:
        return ""
    status = str(validation.get("status") or "UNKNOWN").upper()
    status_kind = "success" if status == "PASS" else "error"
    input_vcf = validation.get("input_vcf")
    card_fields = [
        (
            "info", "Input VCF",
            os.path.basename(str(input_vcf)) if input_vcf else "unavailable",
        ),
        (
            "genetic", "Genome build",
            validation.get("genome_build") or "not validated",
        ),
        (
            "count", "Declared contigs",
            _count(validation.get("declared_contigs")),
        ),
        (
            "count", "Total variants",
            _count(validation.get("total_variants"), missing="not counted"),
        ),
        (
            status_kind, "Header contract", _header_contract_summary(validation),
        ),
    ]
    fields = validation.get("required_fields") or ()
    fields_ok = _required_field_status(validation) == "PASS"
    provenance = validation.get("field_provenance") or {}
    provenance_status = str(
        provenance.get("status") or "UNAVAILABLE"
    ).upper()
    provenance_detail = provenance.get("detail")
    provenance_value = provenance_status
    if provenance_detail:
        provenance_value += " · %s" % provenance_detail
    card_fields.append((
        "success" if provenance_status == "AVAILABLE" else "warning",
        "Field provenance",
        provenance_value,
    ))
    for field_summary in validation.get("field_summaries") or ():
        card_fields.append((
            field_summary.get("kind", "info"),
            field_summary["label"],
            field_summary["value"],
        ))
    if status != "PASS" or not fields_ok:
        for field in fields:
            declared = bool(field.get("declared"))
            required_by = "; ".join(
                str(reason) for reason in field.get("required_by") or ()
            )
            value = "present" if declared else "missing"
            if required_by:
                value += "; required by %s" % required_by
            card_fields.append((
                "success" if declared else "error",
                field.get("field", "Unknown field"),
                value,
            ))
    if validation.get("error"):
        card_fields.append(("error", "Validation error", validation["error"]))
    card_fields = consolidate_validation_fields(card_fields)
    flush_file_validation_display()
    if not card_fields:
        return ""
    return "\n".join((
        "", screen_line(status_kind, "Input VCF validation", indent=4), "",
        *(screen_field(kind, label, value, indent=8, label_width=label_width)
          for kind, label, value in card_fields), "",
    ))


def render_filtering_summary(
    summary: Mapping[str, Any], *, label_width: int,
) -> str:
    """Render one aligned terminal summary from reconciled evidence."""
    indent = 8

    def begin_section(lines: list[str], kind: str, title: str) -> None:
        """Separate summary sections and their contents consistently."""
        lines.extend(screen_section_lines(kind, title))

    lines = [
        "",
        screen_line("analysis", "Variant filtering summary", indent=4),
        "",
        screen_field(
            "info", "Dataset", summary["dataset_id"],
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "genetic", "Genome build", summary["genome_build"],
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "info", "Attribution method", summary["attribution_method"],
            indent=indent, label_width=label_width,
        ),
    ]
    flow_fields = (
        ("Input summary statistics", "total_variant_infile"),
        ("Read by harmonisation", "total_variant_read"),
        ("Invalid coordinates removed", "total_variant_removed_null_coords"),
        (
            "Non-standard alleles removed",
            "total_variant_removed_non_standard_alleles",
        ),
        (
            "Remaining for harmonisation",
            "total_variant_remaining_for_harmonisation",
        ),
        ("Missing effect frequency", "total_variant_with_missing_eaf"),
        ("Invalid effect statistics", "total_variant_with_invalid_beta_se"),
        ("Used for VCF creation", "total_variant_in_vcf_input"),
    )
    data_flow = summary.get("data_flow") or {}
    if any(key in data_flow for _, key in flow_fields):
        begin_section(lines, "count", "Upstream harmonisation flow")
        for label, key in flow_fields:
            if key in data_flow:
                lines.append(screen_field(
                    "loss"
                    if "removed" in label.lower() or "invalid" in label.lower()
                    else "count",
                    label,
                    _count(data_flow.get(key)),
                    indent=indent,
                    label_width=label_width,
                ))

    configured_rules = list(summary.get("rules") or ())
    show_missing = bool(summary.get("display_missing_counts"))
    result_rules = [
        rule for rule in configured_rules
        if rule.get("category") != "missing_value" or show_missing
    ]
    missing_results_hidden = any(
        rule.get("category") == "missing_value" for rule in configured_rules
    ) and not show_missing
    begin_section(lines, "analysis", "Filtering results by rule")
    if result_rules:
        current_group = None
        for index, rule in enumerate(result_rules):
            group, group_kind = _display_rule_group(rule)
            if group is not None and group != current_group:
                lines.extend((screen_line(
                    group_kind, group, indent=indent,
                ), ""))
            current_group = group
            rule_indent = indent + 4 if group is not None else indent
            breakdown_indent = " " * (rule_indent + 8)
            removes = _rule_removes(rule)
            lines.append(screen_line(
                "loss" if removes else "info",
                _display_rule_label(rule),
                indent=rule_indent,
            ))
            if removes:
                lines.extend((
                    screen_line(
                        "info",
                        "%s variants failed this rule and were removed:"
                        % _count(rule.get("count")),
                        indent=rule_indent + 4,
                    ),
                    "%s• %s failed no earlier EXCLUDE rule"
                    % (breakdown_indent, _count(rule.get("primary_count"))),
                    "%s• %s also failed one or more earlier EXCLUDE rules"
                    % (breakdown_indent, _count(_earlier_rule_overlap(rule))),
                ))
            else:
                lines.append(screen_line(
                    "info",
                    "%s variants matched this KEEP rule; this rule did not "
                    "exclude any variants"
                    % _count(rule.get("count")),
                    indent=rule_indent + 4,
                ))
            if index < len(result_rules) - 1:
                lines.append("")
    if missing_results_hidden:
        if result_rules:
            lines.append("")
        lines.append(screen_line(
            "info", "Missing-value rule results are hidden by configuration",
            indent=indent,
        ))
    elif not result_rules:
        lines.append(screen_line(
            "info", "No filtering rule is active", indent=indent,
        ))

    begin_section(lines, "count", "Exact removal accounting")
    lines.extend((
        screen_field(
            "analysis", "Variants failing multiple rules",
            _count(summary.get("overlap_variants")),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "analysis", "Rule failures beyond the first",
            _count(summary.get("extra_rule_matches")),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "loss", "Removed by all filters", _count(summary.get("variants_removed")),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "success" if summary.get("reconciled") else "error",
            "Removal count check",
            "%s assigned to reasons %s %s removed in total"
            % (
                _count(summary.get("primary_removed_total")),
                "=" if summary.get("reconciled") else "!=",
                _count(summary.get("variants_removed")),
            ),
            indent=indent,
            label_width=label_width,
        ),
    ))
    begin_section(lines, "analysis", "Saved reports")
    for kind, label, key in (
        ("info", "Detailed reason report", "reason_summary"),
        ("success", "Summary CSV", "summary_csv"),
        ("success", "Detailed HTML report", "html_report"),
    ):
        value = summary["outputs"].get(key)
        if value:
            lines.append(screen_field(
                kind, label, os.path.basename(str(value)),
                indent=indent, label_width=label_width,
            ))
    begin_section(lines, "count", "Final filtering outcome")
    if summary.get("soft_filter_enabled"):
        lines.append(screen_field(
            "info", "Soft-filter audit VCF",
            "%s variants retained in %s"
            % (
                _count(summary.get("soft_filter_variants")),
                os.path.basename(str(summary["outputs"].get("soft_filtered_vcf"))),
            ),
            indent=indent, label_width=label_width,
        ))
    retained_percent = summary.get("variants_retained_percent")
    retained_value = (
        "unavailable"
        if retained_percent is None
        else "%.2f%%" % float(retained_percent)
    )
    lines.extend((
        screen_field(
            "count", "VCF before filtering", _count(summary.get("variants_before")),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "count", "VCF after filtering", _count(summary.get("variants_after")),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "success", "Variants retained", retained_value,
            indent=indent, label_width=label_width,
        ),
    ))
    lines.append("")
    return "\n".join(lines)


def render_filtering_plan(rules: Sequence[Mapping[str, Any]]) -> str:
    """Render active policies once in their exact first-failure audit order."""
    return render_action_plan(
        "Filtering plan · applied in this order",
        [
            {
                "action": _display_rule_action(rule),
                "label": _display_rule_label(rule),
            }
            for rule in rules
        ],
        empty_message=(
            "No removal or missing-value rule is active; all records pass"
        ),
    )


def _html_value(value: Any) -> str:
    if value is None or value == "":
        return "Not available"
    if value is True:
        return "Yes"
    if value is False:
        return "No"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.6g}"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _cell(value: Any) -> str:
    return escape(_html_value(value))


def _path_cell(value: Any) -> str:
    if value is None or value == "":
        return '<span class="muted">Not generated</span>'
    path = Path(str(value)).expanduser()
    try:
        href = path.resolve().as_uri()
    except (OSError, ValueError):
        href = ""
    label = escape(str(value))
    if not href:
        return '<span class="path">%s</span>' % label
    return '<a class="path" href="%s">%s</a>' % (
        escape(href, quote=True), label,
    )


def _table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    head = "".join(
        '<th scope="col">%s</th>' % escape(header) for header in headers
    )
    body = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % _cell(value) for value in row)
        for row in rows
    )
    if not rows:
        body = '<tr><td colspan="%d" class="muted">No rows</td></tr>' % len(headers)
    return '<div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (
        head, body,
    )


def _key_value_table(
    rows: Sequence[tuple[str, Any]], *, path_labels: Sequence[str] = (),
) -> str:
    paths = set(path_labels)
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (escape(label), _path_cell(value) if label in paths else _cell(value))
        for label, value in rows
    )
    return '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody></table></div>' % body


def _rule_count_interpretation(summary: Mapping[str, Any]) -> str:
    """Explain first-failure attribution with one observed run example."""
    removal_rules = [
        rule for rule in summary.get("rules") or () if _rule_removes(rule)
    ]
    example_rule = next(
        (
            rule
            for rule in removal_rules
            if (
                _earlier_rule_overlap(rule) is not None
                and _earlier_rule_overlap(rule) > 0
            )
        ),
        None,
    )
    if example_rule is None:
        worked_example = (
            '<div class="worked-example"><h3>Worked example</h3>'
            '<p>No active rule in this run had an earlier-rule overlap greater '
            'than zero. For every removal rule, the independent failure count '
            'therefore equals its primary-removal count.</p></div>'
        )
    else:
        failed = example_rule.get("count")
        primary = example_rule.get("primary_count")
        overlap = _earlier_rule_overlap(example_rule)
        worked_example = (
            '<div class="worked-example"><h3>Worked example from this run</h3>'
            '<p><strong>%s</strong></p><dl class="example-counts">'
            '<div><dt>Variants failing this rule</dt><dd>%s</dd></div>'
            '<div><dt>Primary removals assigned here</dt><dd>%s</dd></div>'
            '<div><dt>Overlapping earlier removal rules</dt><dd>%s</dd></div>'
            '</dl><p>Every variant in the first count fails this rule and is '
            'excluded from the hard-filtered VCF. The second and third counts '
            'partition that failure count according to whether a primary '
            'removal reason had already been assigned.</p></div>'
            % (
                escape(str(example_rule.get("label", "Removal rule"))),
                _cell(failed),
                _cell(primary),
                _cell(overlap),
            )
        )
    return (
        '<section id="rule-count-interpretation"><h2>How to read rule counts</h2>'
        '<p class="section-note">Each active removal rule is evaluated '
        'independently, but final removal accounting assigns every excluded '
        'variant to exactly one primary reason.</p>'
        '<div class="definition-grid">'
        '<article class="definition"><h3>1. Variants failing this rule</h3>'
        '<p>The independent number of variants satisfying this rule\'s failure '
        'expression. A variant may appear in this count for several rules, so '
        'these values must not be summed to estimate total removals.</p></article>'
        '<article class="definition"><h3>2. Primary removals assigned here</h3>'
        '<p>The subset first encountered at this rule after applying the '
        'recorded audit order. These counts are mutually exclusive and their '
        'sum must equal the number removed from the hard-filtered VCF.</p></article>'
        '<article class="definition"><h3>3. Overlapping earlier removal rules</h3>'
        '<p>Variants failing this rule whose primary reason was already assigned '
        'to an earlier rule. They are still removed; they are simply not counted '
        'a second time under this rule.</p></article></div>'
        '<div class="formula"><strong>Per-rule identity</strong><code>'
        'overlapping earlier removal rules = variants failing this rule − '
        'primary removals assigned here</code></div>'
        '%s'
        '<div class="callout"><strong>Accounting rule:</strong> sum primary '
        'removals, not independent rule-failure counts. The reconciliation '
        'status confirms that primary removal assignments equal actual '
        'hard-filter removals.</div></section>'
        % worked_example
    )


def _input_validation_html(summary: Mapping[str, Any]) -> str:
    """Explain the validated input contract and each active field dependency."""
    validation = summary.get("input_validation") or {}
    if not validation:
        return ""
    status = str(validation.get("status") or "UNKNOWN").upper()
    status_class = "ok" if status == "PASS" else "fail"
    field_rows = [
        (
            field.get("field"),
            "; ".join(str(value) for value in field.get("required_by") or ()),
            "Present" if field.get("declared") else "Missing",
        )
        for field in validation.get("required_fields") or ()
    ]
    field_table = _table(
        ("Required VCF field", "Required by active filter", "Header status"),
        field_rows,
    )
    error = ""
    if validation.get("error"):
        error = (
            '<div class="notice validation-error"><strong>Validation error:</strong> '
            "%s</div>" % _cell(validation["error"])
        )
    provenance = validation.get("field_provenance") or {}
    provenance_value = str(
        provenance.get("status") or "UNAVAILABLE"
    ).upper()
    if provenance.get("detail"):
        provenance_value += " · %s" % provenance["detail"]
    validation_rows = [
        ("Input VCF", validation.get("input_vcf")),
        ("Genome build", validation.get("genome_build")),
        ("Declared contigs", validation.get("declared_contigs")),
        ("Total variants", validation.get("total_variants")),
        ("Header contract", _header_contract_summary(validation)),
        ("Field provenance", provenance_value),
    ]
    validation_rows.extend(
        (field["label"], field["value"])
        for field in validation.get("field_summaries") or ()
    )
    return (
        '<section id="input-vcf-validation"><h2>Input VCF validation</h2>'
        '<p class="section-note">PostGWAS checks the authoritative genome-build '
        'declaration, contig declarations, and every INFO/FORMAT header field '
        'needed by the selected filters before variant filtering starts.</p>'
        '<p><strong>Contract status:</strong> <span class="badge %s">%s</span></p>'
        "%s"
        '<h3>Required field declarations</h3><p class="section-note">%s. A '
        '<strong>Present</strong> result means the field is defined in the VCF '
        'header. It does not mean every record contains a value; record-level '
        'missing values are measured separately in the filtering-rule table.</p>%s'
        "%s</section>"
        % (
            status_class,
            escape(status),
            _key_value_table(
                validation_rows, path_labels=("Input VCF",),
            ),
            escape(_required_field_summary(validation)),
            field_table,
            error,
        )
    )


_STYLES = """
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ee;--canvas:#f4f7fb;--panel:#fff;--brand:#155e75;--brand2:#0891b2;--ok:#166534;--warn:#92400e}*{box-sizing:border-box}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page-header{padding:36px max(24px,calc((100vw - 1320px)/2));color:#fff;background:linear-gradient(135deg,#164e63,#0e7490)}.eyebrow{margin:0 0 5px;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase;opacity:.86}.page-header h1{margin:0;font-size:clamp(28px,4vw,44px);line-height:1.1}.page-header p{max-width:900px;margin:12px 0 0;color:#cffafe}.container{max-width:1320px;margin:0 auto;padding:26px 24px 58px}.notice,.callout,section{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:0 7px 24px rgba(15,23,42,.05)}.notice{margin-bottom:20px;padding:14px 16px;border-left:4px solid var(--brand2)}.notice.validation-error{margin:16px 0 0;border-left-color:#dc2626;background:#fef2f2}.callout{margin:16px 0 0;padding:14px 16px;border-left:4px solid var(--brand2)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:12px;margin:0 0 22px}.metric{padding:15px;border:1px solid var(--line);border-radius:11px;background:#fff}.metric-label{font-size:12px;font-weight:800;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.metric-value{margin-top:4px;font-size:22px;font-weight:780;overflow-wrap:anywhere}.badge{display:inline-flex;padding:5px 9px;border-radius:999px;font-size:12px;font-weight:850}.badge.ok{color:var(--ok);background:#dcfce7}.badge.fail{color:#991b1b;background:#fee2e2}section{margin:18px 0;padding:22px}section h2{margin:0 0 5px;font-size:21px}section h3{margin:18px 0 5px;font-size:17px}.section-note{margin:0 0 14px;color:var(--muted)}.definition-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.definition{padding:16px;border:1px solid var(--line);border-radius:10px;background:#f8fafc}.definition h3,.worked-example h3{margin:0 0 7px;font-size:16px}.definition p,.worked-example p{margin:0}.worked-example p+p{margin-top:8px}.formula,.worked-example{margin-top:14px;padding:15px 16px;border-radius:10px;background:#ecfeff}.formula strong{display:block;margin-bottom:5px;color:var(--brand)}.formula code{font:600 13px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;white-space:normal}.example-counts{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px;margin:12px 0}.example-counts div{padding:10px;border:1px solid #a5f3fc;border-radius:8px;background:#fff}.example-counts dt{font-size:12px;font-weight:800;color:var(--muted);text-transform:uppercase}.example-counts dd{margin:3px 0 0;font-size:20px;font-weight:800;color:var(--brand)}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px}table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}thead th{background:#ecfeff;color:#334155;font-size:12px;text-transform:uppercase;letter-spacing:.04em}.key-value th{width:290px;background:#f8fafc;color:#475569}.path{overflow-wrap:anywhere;word-break:break-word;color:#075985}.muted{color:var(--muted)}.footer{max-width:1320px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:900px){.definition-grid,.example-counts{grid-template-columns:1fr}}@media(max-width:700px){.container{padding:16px 10px 40px}section{padding:14px}.key-value th{width:42%}th,td{padding:8px}.page-header{padding:28px 18px}}@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0}.page-header p{color:var(--muted)}.container{max-width:none;padding:0}.notice,section{box-shadow:none;break-inside:avoid}.table-wrap{overflow:visible}a{color:inherit;text-decoration:none}}
"""


def render_filtering_html(summary: Mapping[str, Any]) -> str:
    """Render a self-contained filtering report from reconciled evidence."""
    retained = summary.get("variants_retained_percent")
    cards = (
        ("Input variants", _count(summary.get("variants_before"))),
        ("Retained variants", _count(summary.get("variants_after"))),
        ("Removed variants", _count(summary.get("variants_removed"))),
        (
            "Retained",
            "Not available" if retained is None else "%.2f%%" % float(retained),
        ),
        (
            "Variants failing multiple rules",
            _count(summary.get("overlap_variants")),
        ),
    )
    metrics = '<div class="metrics">%s</div>' % "".join(
        '<div class="metric"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div></div>'
        % (escape(label), escape(value))
        for label, value in cards
    )
    status_class = "ok" if summary.get("reconciled") else "fail"
    reconciliation = (
        '<div class="notice"><strong>Removal reconciliation:</strong> '
        '<span class="badge %s">%s</span> %s primary removal assignments and '
        '%s hard removals.</div>'
        % (
            status_class,
            escape(str(summary.get("reconciliation_status"))),
            _cell(summary.get("primary_removed_total")),
            _cell(summary.get("variants_removed")),
        )
    )

    rule_rows = []
    for order, rule in enumerate(summary.get("rules") or (), 1):
        removes = _rule_removes(rule)
        failure_count = (
            rule.get("count")
            if removes
            else "%s observed; not a removal rule" % _count(rule.get("count"))
        )
        rule_rows.append((
            order,
            rule.get("tag") if removes else "Not applied",
            str(rule.get("category", "")).replace("_", " ").title(),
            _display_rule_group(rule)[0] or "Not specified",
            _display_rule_label(rule),
            _display_rule_action(rule),
            failure_count,
            rule.get("primary_count") if removes else "Not applicable",
            _earlier_rule_overlap(rule) if removes else "Not applicable",
            rule.get("expr"),
        ))
    rule_table = _table(
        (
            "Order", "FILTER ID", "Category", "Rule family", "Rule", "Action",
            "Variants failing this rule", "Primary removals assigned here",
            "Overlapping earlier removal rules", "bcftools expression",
        ),
        rule_rows,
    )

    flow_rows = [
        (key.replace("_", " ").title(), value)
        for key, value in (summary.get("data_flow") or {}).items()
    ]
    flow_rows.extend((
        ("VCF before filtering", summary.get("variants_before")),
        ("VCF after filtering", summary.get("variants_after")),
        ("Removed by all filters", summary.get("variants_removed")),
        ("Variants failing multiple rules", summary.get("overlap_variants")),
        ("Rule failures beyond the first", summary.get("extra_rule_matches")),
        ("Soft-filter VCF variants", summary.get("soft_filter_variants")),
    ))

    module = (summary.get("resolved_configuration") or {}).get("filtering") or {}
    execution = (summary.get("resolved_configuration") or {}).get("execution") or {}
    settings = (
        ("Minimum MAF", module.get("maf_min")),
        ("Minimum INFO", module.get("info_min")),
        ("Maximum INFO", module.get("info_max")),
        ("Minimum -log10(P)", module.get("minimum_neglog10_p")),
        ("Maximum study/reference AF difference", module.get("frequency_difference_max")),
        ("Reference population tag", module.get("reference_population_tag")),
        ("Include indels and other non-SNPs", module.get("include_indels")),
        ("Remove ambiguous palindromic SNPs", module.get("remove_palindromic")),
        ("Remove MHC", module.get("remove_mhc")),
        ("Write soft-filter VCF", module.get("write_soft_filter_vcf")),
        ("Missing AF action", module.get("missing_af_action")),
        ("Missing INFO action", module.get("missing_info_action")),
        ("Missing P-value action", module.get("missing_pvalue_action")),
        ("Threads", execution.get("threads")),
        ("Memory (GB)", execution.get("memory_gb")),
    )
    outputs = summary["outputs"]
    output_rows = (
        ("Filtered VCF", outputs.get("filtered_vcf")),
        ("Filtered VCF index", outputs.get("filtered_vcf_index")),
        ("Soft-filter VCF", outputs.get("soft_filtered_vcf")),
        ("Soft-filter VCF index", outputs.get("soft_filtered_vcf_index")),
        ("Detailed reason TSV", outputs.get("reason_summary")),
        ("Summary CSV", outputs.get("summary_csv")),
        ("Detailed HTML report", outputs.get("html_report")),
        ("Filtering log", outputs.get("filter_log")),
        ("MHC exclusion BED", outputs.get("mhc_exclusion_bed")),
    )
    executable_rows = tuple(
        (name, value)
        for name, value in (
            (summary.get("resolved_configuration") or {}).get("executables") or {}
        ).items()
    )
    sections = "".join((
        _input_validation_html(summary),
        _rule_count_interpretation(summary),
        '<section><h2>Variant flow and exact accounting</h2><p class="section-note">Counts are copied from the validated filtering audit and hard-filter output.</p>%s</section>'
        % _key_value_table(tuple(flow_rows)),
        '<section><h2>Filtering plan and rule results</h2><p class="section-note">The table lists every applied EXCLUDE or KEEP policy in the recorded first-failure audit order, followed by its observed result. Independent rule-failure counts may overlap; primary-removal counts are mutually exclusive. Rows with action KEEP report an observed condition count and mark removal attribution as not applicable.</p>%s</section>'
        % rule_table,
        '<section><h2>Resolved scientific policy</h2><p class="section-note">These are the effective schema-validated settings used for this run.</p>%s</section>'
        % _key_value_table(settings),
        '<section><h2>Outputs</h2>%s</section>'
        % _key_value_table(output_rows, path_labels=tuple(label for label, _ in output_rows)),
        '<section><h2>Execution provenance</h2>%s</section>'
        % _key_value_table((
            ("Dataset", summary.get("dataset_id")),
            ("Genome build", summary.get("genome_build")),
            ("Filtering runtime seconds", summary.get("runtime_seconds")),
            ("Requested threads", (
                summary.get("thread_runtime") or {}
            ).get("requested_threads")),
            ("Effective Polars threads", (
                summary.get("thread_runtime") or {}
            ).get("polars_thread_pool_size")),
            ("Polars thread budget enforced", (
                summary.get("thread_runtime") or {}
            ).get("thread_budget_enforced")),
            *executable_rows,
        ), path_labels=tuple(name for name, _ in executable_rows)),
    ))
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s — PostGWAS filtering report</title>
<style>%s</style>
</head>
<body>
<header class="page-header"><p class="eyebrow">PostGWAS</p><h1>Variant filtering report</h1><p>%s · %s · hard-filter output, rule-level attribution, QC accounting, and run provenance.</p></header>
<main class="container">
<div class="notice"><strong>Evidence reuse:</strong> this report is rendered from the reconciled statistics already used for the terminal summary and CSV. It does not reread the VCF or calculate a second set of scientific metrics.</div>
%s
%s
%s
</main>
<footer class="footer">PostGWAS filtering report · consult the linked CSV, reason TSV, and canonical log for machine-readable evidence.</footer>
</body>
</html>
""" % (
        escape(str(summary["dataset_id"])),
        _STYLES,
        escape(str(summary["dataset_id"])),
        escape(str(summary["genome_build"])),
        metrics,
        reconciliation,
        sections,
    )


def write_filtering_html_report(
    summary: Mapping[str, Any], destination: str | Path,
) -> Path:
    """Write the self-contained detailed filtering report atomically."""
    return write_html_report(render_filtering_html(summary), destination)


__all__ = [
    "SUMMARY_COLUMNS",
    "build_filtering_summary",
    "filtering_summary_records",
    "render_filtering_html",
    "render_filtering_plan",
    "render_filtering_summary",
    "render_input_vcf_validation",
    "write_filtering_html_report",
    "write_filtering_summary_csv",
]
