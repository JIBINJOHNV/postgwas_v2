"""Human-readable formatter reports from already-validated result metadata.

The screen and HTML renderers are presentation-only. They never reopen the
GWAS-VCF, exported tables, or LDSC reference, and therefore cannot calculate a
second, potentially divergent set of scientific results.
"""

from __future__ import annotations

from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import quote

from postgwas.core.input_validation import current_validation_session
from postgwas.core.io.reports import write_html_report
from postgwas.core.ui.screen import screen_field, screen_line

from .contracts import (
    CUSTOM_OUTPUT_TARGET,
    FORMAT_CONTRACTS,
    formatter_target_display_name,
)


def _count(value: Any) -> str:
    if value is None:
        return "Not available"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def _number(value: Any) -> str:
    if value is None:
        return "Not available"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.10g}"
    return str(value)


def _quantity(value: Any, singular: str, plural: str | None = None) -> str:
    count = int(value)
    label = singular if count == 1 else (plural or singular + "s")
    return "%s %s" % (_count(count), label)


def _percentage(numerator: Any, denominator: Any) -> float | None:
    try:
        numerator_value = float(numerator)
        denominator_value = float(denominator)
    except (TypeError, ValueError):
        return None
    if denominator_value <= 0:
        return None
    return 100.0 * numerator_value / denominator_value


def _percentage_text(numerator: Any, denominator: Any) -> str:
    value = _percentage(numerator, denominator)
    if value is None:
        return "percentage unavailable"
    if value in (0.0, 100.0):
        return "%.0f%%" % value
    precision = 3 if value < 0.1 else 2
    return ("%%.%df%%%%" % precision) % value


def _target_name(
    target: str,
    variant_id_observations: Mapping[str, Mapping] | None = None,
) -> str:
    return formatter_target_display_name(target, variant_id_observations)


def _identifier_label(result: Mapping[str, Any]) -> str:
    return (
        "rsID"
        if result.get("variant_id_type") == "rsid"
        else "configured coordinate-and-allele ID"
    )


def _exclusion_rows(
    result: Mapping[str, Any],
) -> list[tuple[str, int]]:
    """Return disjoint formatter-stage exclusions in processing order."""
    identifier = (
        "Missing or invalid rsID"
        if result.get("variant_id_type") == "rsid"
        else "Missing or invalid configured identifier"
    )
    rows = [
        (
            identifier,
            int(result.get(
                "identifier_missing_rows_excluded",
                result.get("identifier_rows_excluded", 0),
            )),
        ),
        (
            "Rows removed by duplicate policy",
            int(result.get("identifier_duplicate_rows_excluded", 0)),
        ),
        (
            "No rsID match in supplied LDSC SNP/allele reference",
            int(result.get("rows_excluded_not_in_reference", 0)),
        ),
        (
            "Allele mismatch with supplied LDSC SNP/allele reference",
            int(result.get("rows_excluded_reference_allele_mismatch", 0)),
        ),
        (
            "Outside configured chromosomes",
            int(result.get("rows_excluded_unconfigured_chromosome", 0)),
        ),
        (
            "Missing or invalid required output values",
            int(result.get("schema_rows_excluded", 0)),
        ),
    ]
    categorized = sum(count for _, count in rows)
    other = int(result.get("rows_excluded", 0)) - categorized
    if other > 0:
        rows.append(("Other recorded exclusions", other))
    return [(label, count) for label, count in rows if count]


def _common_input_rows(results: Mapping[str, Mapping[str, Any]]) -> int | None:
    values = {
        int(result["rows_in"])
        for result in results.values()
        if result.get("rows_in") is not None
    }
    return next(iter(values)) if len(values) == 1 else None


def ldsc_sample_prevalence_screen_fields(
    results,
    case_count_column,
    control_count_column,
    label_width,
    *,
    indent=6,
):
    """Render the case fraction already calculated by the LDSC exporter."""
    result = results.get("ldsc")
    if result is None:
        return []
    if result.get("trait_type") != "binary":
        return [screen_field(
            "analysis",
            "GWAS-VCF case fraction",
            "not applicable to a quantitative trait",
            indent=indent,
            label_width=label_width,
        )]

    value = result.get("sample_prev")
    aggregation = result.get("sample_prevalence_aggregation")
    variants = result.get("sample_prevalence_variants")
    minimum = result.get("sample_prevalence_minimum")
    maximum = result.get("sample_prevalence_maximum")
    if any(
        item is None
        for item in (value, aggregation, variants, minimum, maximum)
    ):
        return [screen_field(
            "warning",
            "GWAS-VCF case fraction",
            "unavailable because the formatter result metadata is incomplete",
            indent=indent,
            label_width=label_width,
        )]

    fields = [
        screen_field(
            "analysis",
            "Calculated case fraction",
            "%.10g (%.2f%%)" % (float(value), float(value) * 100.0),
            indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "Formula",
            "%s / (%s + %s)"
            % (case_count_column, case_count_column, control_count_column),
            indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "Aggregation",
            "%s across retained LDSC-input variants" % aggregation,
            indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "count",
            "Variants used",
            _count(variants),
            indent=indent,
            label_width=label_width,
        ),
    ]
    if float(minimum) == float(maximum):
        fields.append(screen_field(
            "info",
            "Case-fraction variation",
            "none — %.10g for every retained variant" % float(minimum),
            indent=indent,
            label_width=label_width,
        ))
    else:
        fields.append(screen_field(
            "info",
            "Case-fraction range",
            "%.10g to %.10g" % (float(minimum), float(maximum)),
            indent=indent,
            label_width=label_width,
        ))

    count_ranges = (
        result.get("sample_prevalence_case_count_minimum"),
        result.get("sample_prevalence_case_count_maximum"),
        result.get("sample_prevalence_control_count_minimum"),
        result.get("sample_prevalence_control_count_maximum"),
    )
    if any(item is None for item in count_ranges):
        fields.append(screen_field(
            "warning",
            "%s/%s values" % (case_count_column, control_count_column),
            "unavailable in stored completion metadata; rerun with --overwrite",
            indent=indent,
            label_width=label_width,
        ))
        return fields

    case_minimum, case_maximum, control_minimum, control_maximum = map(
        float, count_ranges,
    )
    if case_minimum == case_maximum and control_minimum == control_maximum:
        total = case_minimum + control_minimum
        fields.extend([
            screen_field(
                "count",
                "Reported %s" % case_count_column,
                "%s (%s)"
                % (
                    format(case_minimum, ",.10g"),
                    _percentage_text(case_minimum, total),
                ),
                indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "count",
                "Reported %s" % control_count_column,
                "%s (%s)"
                % (
                    format(control_minimum, ",.10g"),
                    _percentage_text(control_minimum, total),
                ),
                indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "count",
                "Reported total sample count",
                format(total, ",.10g"),
                indent=indent,
                label_width=label_width,
            ),
        ])
    else:
        fields.extend([
            screen_field(
                "count",
                "%s range" % case_count_column,
                "%s to %s"
                % (
                    format(case_minimum, ",.10g"),
                    format(case_maximum, ",.10g"),
                ),
                indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "count",
                "%s range" % control_count_column,
                "%s to %s"
                % (
                    format(control_minimum, ",.10g"),
                    format(control_maximum, ",.10g"),
                ),
                indent=indent,
                label_width=label_width,
            ),
        ])
    return fields


def _duplicate_screen_lines(
    result,
    label_width,
    *,
    heading_indent=6,
    field_indent=6,
):
    groups = int(result.get("identifier_duplicate_groups", 0))
    rows = int(result.get("identifier_duplicate_rows", 0))
    input_rows = int(result.get("rows_in", 0))
    validated = result.get("identifier_uniqueness_validated") is True
    if groups == 0:
        return [
            "",
            screen_line(
                "analysis", "Identifier validation", indent=heading_indent,
            ),
            screen_field(
                "success" if validated else "warning",
                "Final identifier validation",
                (
                    "Passed — all %s retained identifiers are unique"
                    % _count(result.get("rows_out"))
                    if validated
                    else "uniqueness validation was not recorded"
                ),
                indent=field_indent,
                label_width=label_width,
            ),
        ]

    lines = [
        "",
        screen_line(
            "analysis",
            "Duplicate %s handling" % _identifier_label(result),
            indent=heading_indent,
        ),
        screen_field(
            "count",
            "Duplicate groups found",
            "%s containing %s (%s of GWAS-VCF rows)"
            % (
                _quantity(groups, "group"),
                _quantity(rows, "row"),
                _percentage_text(rows, input_rows),
            ),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "Configured duplicate policy",
            result.get("identifier_duplicate_policy", "not recorded"),
            indent=field_indent,
            label_width=label_width,
        ),
    ]
    exact_collapsed = int(result.get(
        "identifier_exact_duplicate_rows_collapsed", 0,
    ))
    if exact_collapsed:
        lines.append(screen_field(
            "info",
            "Exact repeated rows collapsed",
            "%s (%s of GWAS-VCF rows)"
            % (
                _count(exact_collapsed),
                _percentage_text(exact_collapsed, input_rows),
            ),
            indent=field_indent,
            label_width=label_width,
        ))
    if result.get("reference_selection_applied") is True:
        resolved = int(result.get(
            "identifier_duplicate_groups_resolved_by_reference", 0,
        ))
        ambiguous = int(result.get(
            "identifier_duplicate_groups_unresolved_after_reference", 0,
        ))
        not_retained = max(0, groups - resolved - ambiguous)
        lines.extend([
            screen_field(
                "success",
                "Resolved by supplied LDSC SNP/allele reference",
                "%s (%s of duplicate groups)"
                % (_quantity(resolved, "group"), _percentage_text(resolved, groups)),
                indent=field_indent,
                label_width=label_width,
            ),
            screen_field(
                "info",
                "Not retained after reference matching",
                "%s (%s)"
                % (_quantity(not_retained, "group"), _percentage_text(not_retained, groups)),
                indent=field_indent,
                label_width=label_width,
            ),
            screen_field(
                "info",
                "Remaining ambiguous groups",
                "%s (%s)"
                % (_quantity(ambiguous, "group"), _percentage_text(ambiguous, groups)),
                indent=field_indent,
                label_width=label_width,
            ),
        ])
    else:
        resolved = int(result.get(
            "identifier_duplicate_groups_resolved_by_policy", 0,
        ))
        unresolved = int(result.get(
            "identifier_duplicate_groups_unresolved", 0,
        ))
        lines.extend([
            screen_field(
                "info",
                "Resolved by duplicate policy",
                "%s (%s)"
                % (_quantity(resolved, "group"), _percentage_text(resolved, groups)),
                indent=field_indent,
                label_width=label_width,
            ),
            screen_field(
                "info",
                "Unresolved groups excluded",
                "%s (%s)"
                % (_quantity(unresolved, "group"), _percentage_text(unresolved, groups)),
                indent=field_indent,
                label_width=label_width,
            ),
        ])
    lines.extend([
        screen_field(
            "info",
            "Rows removed by duplicate policy",
            "%s (%s of GWAS-VCF rows)"
            % (
                _quantity(
                    result.get("identifier_duplicate_rows_excluded", 0),
                    "row",
                ),
                _percentage_text(
                    result.get("identifier_duplicate_rows_excluded", 0),
                    input_rows,
                ),
            ),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            "success" if validated else "warning",
            "Final identifier validation",
            (
                "Passed — all %s retained identifiers are unique"
                % _count(result.get("rows_out"))
                if validated
                else "uniqueness validation was not recorded"
            ),
            indent=field_indent,
            label_width=label_width,
        ),
    ])
    return lines


def render_formatter_screen_summary(
    dataset_id,
    vcf,
    study_design,
    results,
    log_path,
    html_report_path,
    case_count_column,
    control_count_column,
    label_width,
    variant_id_observations,
    output_destinations,
    input_evidence,
    minimum_p_value,
):
    """Render a concise but scientifically explicit completion summary."""
    field_indent = 10
    nested_field_indent = 14
    nested_label_width = label_width - (nested_field_indent - field_indent)
    prepared_for = [
        _target_name(target, variant_id_observations) for target in results
    ]
    completion_heading = (
        "%s input files prepared" % prepared_for[0]
        if len(prepared_for) == 1
        else "Downstream-analysis input files prepared"
    )
    lines = [
        "",
        screen_line("analysis", completion_heading, indent=2),
        "",
        screen_line("analysis", "Source GWAS-VCF", indent=6),
        screen_field(
            "info", "Dataset", dataset_id,
            indent=field_indent, label_width=label_width,
        ),
        screen_field(
            "genetic", "Input GWAS-VCF", vcf,
            indent=field_indent, label_width=label_width, path_value=True,
        ),
        screen_field(
            "count",
            "Total variants in input GWAS-VCF",
            _count(_common_input_rows(results)),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            "genetic",
            "Genome build",
            input_evidence.get("genome_build", "not recorded"),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "VCF embedded dataset/sample",
            input_evidence.get("postgwas_dataset_id", "not recorded"),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            (
                "success"
                if str(input_evidence.get("postgwas_dataset_id")) == dataset_id
                else "warning"
            ),
            "Run ID versus VCF identity",
            (
                "match"
                if str(input_evidence.get("postgwas_dataset_id")) == dataset_id
                else "different — verify that this is the intended GWAS-VCF"
            ),
            indent=field_indent,
            label_width=label_width,
        ),
        screen_field(
            "success",
            "VCF provenance status",
            "%s; PostGWAS %s"
            % (
                input_evidence.get("postgwas_status", "not recorded"),
                input_evidence.get("postgwas_version", "not recorded"),
            ),
            indent=field_indent,
            label_width=label_width,
        ),
    ]
    if study_design is None:
        lines.append(screen_field(
            "info",
            "Trait-type inference",
            "not needed to prepare the selected input files",
            indent=field_indent,
            label_width=label_width,
        ))
    else:
        trait_label = (
            "Binary case-control"
            if study_design.trait_type == "binary"
            else "Quantitative"
        )
        lines.extend([
            screen_field(
                "analysis",
                "Trait type",
                "%s (inferred from sample-count fields)" % trait_label,
                indent=field_indent,
                label_width=label_width,
            ),
            screen_field(
                (
                    "info"
                    if study_design.trait_type == "quantitative"
                    else (
                        "success"
                        if study_design.case_counts_present == study_design.rows
                        else "warning"
                    )
                ),
                "%s completeness" % case_count_column,
                "%s / %s variants (%s)"
                % (
                    _count(study_design.case_counts_present),
                    _count(study_design.rows),
                    _percentage_text(
                        study_design.case_counts_present, study_design.rows,
                    ),
                ),
                indent=field_indent,
                label_width=label_width,
            ),
            screen_field(
                (
                    "success"
                    if study_design.control_counts_present == study_design.rows
                    else "warning"
                ),
                "%s completeness" % control_count_column,
                "%s / %s variants (%s)"
                % (
                    _count(study_design.control_counts_present),
                    _count(study_design.rows),
                    _percentage_text(
                        study_design.control_counts_present, study_design.rows,
                    ),
                ),
                indent=field_indent,
                label_width=label_width,
            ),
        ])
    for target, observation in variant_id_observations.items():
        identifier = (
            "rsIDs"
            if observation["variant_id_type"] == "rsid"
            else "configured coordinate-and-allele unique IDs"
        )
        lines.append(screen_field(
            "genetic",
            "BIM identifier convention · %s"
            % _target_name(target, variant_id_observations),
            (
                identifier
                if current_validation_session() is not None
                else "%s detected after scanning %s BIM variants"
                % (identifier, _count(observation["variants"]))
            ),
            indent=field_indent,
            label_width=label_width,
        ))
        lines.extend([
            screen_field(
                "info",
                "BIM file · %s"
                % _target_name(target, variant_id_observations),
                observation.get("bim_file", "not recorded"),
                indent=field_indent,
                label_width=label_width,
                path_value=True,
            ),
            screen_field(
                "warning",
                "GWAS–BIM variant overlap · %s"
                % _target_name(target, variant_id_observations),
                "not assessed by formatter; the downstream module applies its "
                "resolved reference-matching policy",
                indent=field_indent,
                label_width=label_width,
            ),
        ])

    for target, result in results.items():
        target_name = _target_name(target, variant_id_observations)
        input_rows = int(result.get("rows_in", 0))
        retained = int(result.get("rows_out", 0))
        excluded = int(result.get("rows_excluded", input_rows - retained))
        heading = (
            "LDSC files created"
            if target == "ldsc"
            else "%s files created" % target_name
        )
        lines.extend([
            "",
            screen_line("analysis", heading, indent=6),
            "",
            screen_line(
                "analysis", "Variant accounting", indent=10,
            ),
            screen_field(
                "count",
                "Variants read from GWAS-VCF",
                _count(input_rows),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
        ])
        reference_variants = result.get("reference_variants")
        if reference_variants is not None:
            lines.append(screen_field(
                "genetic",
                "Supplied LDSC SNP/allele reference",
                "%s unique rsIDs" % _count(reference_variants),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ))
        retained_value = "%s (%s of GWAS-VCF variants)" % (
            _count(retained), _percentage_text(retained, input_rows),
        )
        if reference_variants is not None:
            retained_value += (
                "; %s of supplied LDSC SNP/allele reference variants" %
                _percentage_text(retained, reference_variants)
            )
        lines.extend([
            screen_field(
                "success",
                "Variants written to %s files"
                % ("LDSC" if target == "ldsc" else target_name),
                retained_value,
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
            screen_field(
                "count",
                "Variants excluded before writing",
                "%s (%s of GWAS-VCF variants)"
                % (_count(excluded), _percentage_text(excluded, input_rows)),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
            screen_field(
                "success" if int(result.get("schema_rows_excluded", 0)) == 0 else "warning",
                "Required-value validation",
                "%s retained; %s excluded for missing or invalid required values"
                % (
                    _count(result.get("rows_out", 0)),
                    _count(result.get("schema_rows_excluded", 0)),
                ),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
            screen_field(
                (
                    "success"
                    if int(result.get("p_values_bounded", 0)) == 0
                    else "warning"
                ),
                "P values bounded at %g" % minimum_p_value,
                _count(result.get("p_values_bounded", 0)),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
            screen_field(
                "info",
                "Interpretation",
                "formatter retention only; the downstream analysis can retain "
                "fewer variants after its own reference, scope, and mapping checks",
                indent=nested_field_indent,
                label_width=nested_label_width,
            ),
        ])
        for label, count in _exclusion_rows(result):
            lines.append(screen_field(
                "info",
                label,
                "%s (%s of GWAS-VCF variants)"
                % (_count(count), _percentage_text(count, input_rows)),
                indent=nested_field_indent,
                label_width=nested_label_width,
            ))
        if target == "ldsc":
            lines.append(screen_field(
                "info",
                "Important",
                "LDSC munging may apply additional INFO, MAF, sample-size, "
                "and statistical-quality filters",
                indent=nested_field_indent,
                label_width=nested_label_width,
            ))
        lines.extend(_duplicate_screen_lines(
            result,
            nested_label_width,
            heading_indent=10,
            field_indent=nested_field_indent,
        ))
        target_destinations = [
            (label, path)
            for label, path in output_destinations.items()
            if label == target or label.startswith(target + ".")
        ]
        if target_destinations:
            lines.extend([
                "",
                screen_line("analysis", "Prepared scientific files", indent=10),
            ])
            for label, path in target_destinations:
                output_name = (
                    label.split(".", 1)[1].replace("_", " ")
                    if "." in label
                    else "input"
                )
                output_name = output_name.replace("snp", "SNP").replace(
                    "p values", "P-value",
                )
                lines.append(screen_field(
                    "success",
                    "%s file" % output_name,
                    path,
                    indent=nested_field_indent,
                    label_width=nested_label_width,
                    path_value=True,
                ))
        if target == "ldsc":
            lines.extend([
                "",
                screen_line(
                    "analysis", "GWAS-VCF case fraction", indent=10,
                ),
            ])
            lines.extend(ldsc_sample_prevalence_screen_fields(
                results,
                case_count_column,
                control_count_column,
                nested_label_width,
                indent=nested_field_indent,
            ))
    lines.extend([
        "",
        screen_line("analysis", "Reports and logs", indent=6),
        screen_field(
            "success",
            "Detailed HTML report",
            html_report_path,
            indent=field_indent,
            label_width=label_width,
            path_value=True,
        ),
        screen_field(
            "info",
            "Full formatter log",
            log_path,
            indent=field_indent,
            label_width=label_width,
            path_value=True,
        ),
        "",
    ])
    return "\n".join(lines)


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
        return f"{value:,.10g}"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _cell(value: Any) -> str:
    return escape(_html_value(value))


def _path_cell(
    value: Any,
    *,
    output_directory: Path,
    report_path: Path,
) -> str:
    if value is None or value == "":
        return '<span class="muted">Not generated</span>'
    path = Path(str(value)).expanduser()
    try:
        resolved = path.resolve()
    except (OSError, ValueError):
        return '<span class="path">%s</span>' % escape(str(value))
    try:
        relative_to_root = resolved.relative_to(output_directory.resolve())
        display = (
            "Formatter output root"
            if not relative_to_root.parts
            else relative_to_root.as_posix()
        )
        relative = os.path.relpath(resolved, report_path.parent.resolve())
        href = quote(Path(relative).as_posix())
    except ValueError:
        display = str(resolved)
        try:
            href = resolved.as_uri()
        except ValueError:
            href = ""
    if not href:
        return '<span class="path">%s</span>' % escape(display)
    return '<a class="path" href="%s">%s</a>' % (
        escape(href, quote=True), escape(display),
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
        body = '<tr><td colspan="%d" class="muted">No rows</td></tr>' % len(
            headers
        )
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (head, body)
    )


def _key_value_table(rows: Sequence[tuple[str, Any]]) -> str:
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (escape(label), _cell(value))
        for label, value in rows
    )
    return (
        '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody>'
        '</table></div>' % body
    )


def _path_table(
    rows: Sequence[tuple[str, Any]],
    *,
    output_directory: Path,
    report_path: Path,
) -> str:
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (
            escape(label),
            _path_cell(
                value,
                output_directory=output_directory,
                report_path=report_path,
            ),
        )
        for label, value in rows
    )
    return (
        '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody>'
        '</table></div>' % body
    )


def _metric_cards(rows: Sequence[tuple[str, Any]]) -> str:
    return '<div class="metrics">%s</div>' % "".join(
        '<div class="metric"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div></div>'
        % (escape(label), _cell(value))
        for label, value in rows
    )


def _study_design_html(
    study_design,
    case_count_column,
    control_count_column,
) -> str:
    if study_design is None:
        return _key_value_table((
            ("Trait inference", "Not required for the selected formats"),
        ))
    trait = (
        "Binary case-control"
        if study_design.trait_type == "binary"
        else "Quantitative"
    )
    return _key_value_table((
        ("Trait type", "%s (inferred from sample-count fields)" % trait),
        (
            "%s completeness" % case_count_column,
            "%s / %s (%s)"
            % (
                _count(study_design.case_counts_present),
                _count(study_design.rows),
                _percentage_text(
                    study_design.case_counts_present, study_design.rows,
                ),
            ),
        ),
        (
            "%s completeness" % control_count_column,
            "%s / %s (%s)"
            % (
                _count(study_design.control_counts_present),
                _count(study_design.rows),
                _percentage_text(
                    study_design.control_counts_present, study_design.rows,
                ),
            ),
        ),
        (
            "Inference rule",
            "Quantitative only when the configured case-count field is absent "
            "for every variant; otherwise binary",
        ),
    ))


def _target_overview_rows(results, variant_id_observations=None):
    rows = []
    for target, result in results.items():
        input_rows = result.get("rows_in")
        retained = result.get("rows_out")
        excluded = result.get("rows_excluded")
        rows.append((
            _target_name(target, variant_id_observations),
            input_rows,
            retained,
            _percentage_text(retained, input_rows),
            excluded,
            _percentage_text(excluded, input_rows),
            _identifier_label(result),
            result.get("identifier_duplicate_policy"),
            "Passed" if result.get("identifier_uniqueness_validated") else "Not recorded",
        ))
    return rows


def _duplicate_html(result: Mapping[str, Any]) -> str:
    groups = int(result.get("identifier_duplicate_groups", 0))
    resolved_reference = int(result.get(
        "identifier_duplicate_groups_resolved_by_reference", 0,
    ))
    ambiguous_reference = int(result.get(
        "identifier_duplicate_groups_unresolved_after_reference", 0,
    ))
    not_retained_reference = max(
        0, groups - resolved_reference - ambiguous_reference,
    )
    rows = [
        ("Selected identifier type", _identifier_label(result)),
        ("Configured duplicate policy", result.get("identifier_duplicate_policy")),
        ("Duplicate groups detected", groups),
        ("Rows in duplicate groups", result.get("identifier_duplicate_rows", 0)),
        (
            "Rows in duplicate groups (% of input)",
            _percentage_text(
                result.get("identifier_duplicate_rows", 0), result.get("rows_in"),
            ),
        ),
        (
            "Exact repeated rows collapsed",
            result.get("identifier_exact_duplicate_rows_collapsed", 0),
        ),
    ]
    if result.get("reference_selection_applied") is True:
        rows.extend([
            (
                "Groups resolved by supplied LDSC SNP/allele reference",
                resolved_reference,
            ),
            (
                "Resolved by reference (% of duplicate groups)",
                _percentage_text(resolved_reference, groups),
            ),
            (
                "Groups not retained after reference matching",
                not_retained_reference,
            ),
            (
                "Ambiguous groups remaining after reference matching",
                ambiguous_reference,
            ),
        ])
    rows.extend([
        (
            "Conflicting groups reaching duplicate policy",
            result.get("identifier_conflicting_duplicate_groups", 0),
        ),
        (
            "Groups resolved by duplicate policy",
            result.get("identifier_duplicate_groups_resolved_by_policy", 0),
        ),
        (
            "Unresolved groups excluded by duplicate policy",
            result.get("identifier_duplicate_groups_unresolved", 0),
        ),
        (
            "Rows excluded by duplicate policy",
            result.get("identifier_duplicate_rows_excluded", 0),
        ),
        (
            "Final uniqueness invariant",
            (
                "Passed: %s retained rows and %s unique identifiers"
                % (
                    _count(result.get("identifier_validation_rows")),
                    _count(result.get("identifier_unique_values")),
                )
                if result.get("identifier_uniqueness_validated") is True
                else "Not recorded"
            ),
        ),
        ("Identifier selection reused", result.get("identifier_selection_reused")),
    ])
    return _key_value_table(tuple(rows))


def _schema_html(report: Mapping[str, Any]) -> str:
    rows = []
    for output in report.get("outputs") or ():
        for column in output.get("columns") or ():
            rows.append((
                str(output.get("output", "table")).replace("_", " ").title(),
                column.get("source"),
                column.get("saved_as"),
                column.get("transformation") or "Copied without transformation",
            ))
    return _table(
        ("Output", "Canonical source", "Saved column", "Transformation"),
        rows,
    )


def _ldsc_html(
    result: Mapping[str, Any],
    case_count_column: str,
    control_count_column: str,
) -> str:
    reference = ""
    if result.get("reference_selection_applied") is True:
        reference_rows = (
            ("Supplied LDSC SNP/allele reference variants", result.get("reference_variants")),
            ("GWAS-VCF rows entering reference matching", result.get("reference_rows_in")),
            ("Candidate rows retained by rsID and allele matching", result.get("reference_rows_out")),
            ("Unique variants in final LDSC formatter input", result.get("rows_out")),
            (
                "Final LDSC-input coverage of supplied SNP/allele reference",
                _percentage_text(
                    result.get("rows_out"), result.get("reference_variants"),
                ),
            ),
            (
                "No rsID match in supplied LDSC SNP/allele reference",
                result.get("rows_excluded_not_in_reference"),
            ),
            (
                "Allele mismatch with supplied LDSC SNP/allele reference",
                result.get("rows_excluded_reference_allele_mismatch"),
            ),
            (
                "Duplicate groups resolved uniquely by supplied SNP/allele reference",
                result.get("identifier_duplicate_groups_resolved_by_reference"),
            ),
            (
                "Duplicate groups still ambiguous after supplied SNP/allele reference",
                result.get(
                    "identifier_duplicate_groups_unresolved_after_reference"
                ),
            ),
        )
        reference = (
            '<h3>LDSC SNP/allele reference matching</h3>%s'
            '<div class="callout"><strong>Scope:</strong> this is the supplied '
            'SNP/allele list used for formatter and <code>munge_sumstats.py</code> '
            'selection. It is distinct from the chromosome-split LD-score files '
            'supplied later through <code>--ref-ld-chr</code> and '
            '<code>--w-ld-chr</code>.</div>'
            % _key_value_table(reference_rows)
        )

    if result.get("trait_type") != "binary":
        prevalence = _key_value_table((
            ("GWAS-VCF case fraction", "Not applicable to a quantitative trait"),
        ))
    else:
        case_minimum = result.get("sample_prevalence_case_count_minimum")
        case_maximum = result.get("sample_prevalence_case_count_maximum")
        control_minimum = result.get("sample_prevalence_control_count_minimum")
        control_maximum = result.get("sample_prevalence_control_count_maximum")
        prevalence = _key_value_table((
            (
                "Calculated case fraction",
                "%s (%s)"
                % (
                    _number(result.get("sample_prev")),
                    _percentage_text(result.get("sample_prev"), 1),
                ),
            ),
            (
                "Formula",
                "%s / (%s + %s)"
                % (case_count_column, case_count_column, control_count_column),
            ),
            ("Aggregation", result.get("sample_prevalence_aggregation")),
            ("Retained variants used", result.get("sample_prevalence_variants")),
            (
                "Observed case-fraction range",
                "%s to %s"
                % (
                    _number(result.get("sample_prevalence_minimum")),
                    _number(result.get("sample_prevalence_maximum")),
                ),
            ),
            (
                "%s range" % case_count_column,
                "%s to %s" % (_number(case_minimum), _number(case_maximum)),
            ),
            (
                "%s range" % control_count_column,
                "%s to %s"
                % (_number(control_minimum), _number(control_maximum)),
            ),
        ))
    return (
        reference
        + '<h3>GWAS-VCF case fraction</h3>%s'
        '<div class="callout"><strong>Interpretation:</strong> this data-derived '
        'case fraction is calculated only from variants retained in the formatter '
        'LDSC table. LDSC munging can apply additional INFO, MAF, sample-size, and '
        'statistical-quality filters before heritability estimation.</div>'
        % prevalence
    )


def _target_details_html(
    target: str,
    result: Mapping[str, Any],
    schema_report: Mapping[str, Any],
    output_destinations: Mapping[str, str],
    *,
    output_directory: Path,
    report_path: Path,
    case_count_column: str,
    control_count_column: str,
    variant_id_observations: Mapping[str, Mapping] | None = None,
) -> str:
    input_rows = int(result.get("rows_in", 0))
    exclusions = _exclusion_rows(result)
    exclusion_table = _table(
        ("Exclusion reason", "Variants", "% of GWAS-VCF input"),
        [
            (label, count, _percentage_text(count, input_rows))
            for label, count in exclusions
        ],
    )
    flow = _key_value_table((
        ("GWAS-VCF variants", result.get("rows_in")),
        ("Variants retained for this formatter input", result.get("rows_out")),
        (
            "Retained percentage",
            _percentage_text(result.get("rows_out"), result.get("rows_in")),
        ),
        ("Variants excluded", result.get("rows_excluded")),
        (
            "Excluded percentage",
            _percentage_text(result.get("rows_excluded"), result.get("rows_in")),
        ),
        ("P values bounded at configured minimum", result.get("p_values_bounded", 0)),
        ("Sample-size mode", result.get("sample_size_mode")),
        ("Written columns", result.get("columns")),
    ))
    target_paths = [
        (label, path)
        for label, path in output_destinations.items()
        if (
            label == target
            or label.startswith(target + ".")
            or label.startswith(target + "[")
        )
    ]
    source = (
        "Configured PostGWAS custom-column contract"
        if target == CUSTOM_OUTPUT_TARGET
        else FORMAT_CONTRACTS[target].source
    )
    special = (
        _ldsc_html(result, case_count_column, control_count_column)
        if target == "ldsc"
        else ""
    )
    return (
        '<section id="target-%s"><h2>%s</h2>'
        '<p class="section-note">Validated formatter evidence for this target. '
        'Percentages use all extracted GWAS-VCF variants as the denominator.</p>'
        '<h3>Variant flow</h3>%s<h3>Exclusion accounting</h3>%s'
        '<h3>Identifier and duplicate handling</h3>%s%s'
        '<h3>Configured output schema</h3>%s'
        '<h3>Scientific statistic semantics</h3>%s'
        '<h3>Generated target artifacts</h3>%s'
        '<p class="source"><strong>Scientific contract source:</strong> %s</p>'
        '</section>'
        % (
            escape(target, quote=True),
            escape(_target_name(target, variant_id_observations)),
            flow,
            exclusion_table,
            _duplicate_html(result),
            special,
            _schema_html(schema_report),
            _key_value_table((
                ("P-value representation", schema_report.get("p_value")),
                ("Frequency representation", schema_report.get("frequency")),
                ("Frequency interpretation", schema_report.get("frequency_interpretation")),
                ("Sample-size representation", schema_report.get("sample_size")),
            )),
            _path_table(
                target_paths,
                output_directory=output_directory,
                report_path=report_path,
            ),
            escape(source),
        )
    )


_STYLES = """
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ee;--canvas:#f4f7fb;--panel:#fff;--brand:#155e75;--brand2:#0891b2;--ok:#166534;--warn:#92400e}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page-header{padding:36px max(24px,calc((100vw - 1360px)/2));color:#fff;background:linear-gradient(135deg,#164e63,#0e7490)}.eyebrow{margin:0 0 5px;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase;opacity:.86}.page-header h1{margin:0;font-size:clamp(28px,4vw,44px);line-height:1.1}.page-header p{max-width:950px;margin:12px 0 0;color:#cffafe}.container{max-width:1360px;margin:0 auto;padding:26px 24px 58px}.notice,.callout,section{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:0 7px 24px rgba(15,23,42,.05)}.notice{margin-bottom:20px;padding:14px 16px;border-left:4px solid var(--brand2)}.callout{margin:16px 0;padding:14px 16px;border-left:4px solid var(--brand2)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:12px;margin:0 0 22px}.metric{padding:15px;border:1px solid var(--line);border-radius:11px;background:#fff}.metric-label{font-size:12px;font-weight:800;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.metric-value{margin-top:4px;font-size:22px;font-weight:780;overflow-wrap:anywhere}.navigation{margin:18px 0;padding:14px 18px;background:#ecfeff;border:1px solid #a5f3fc;border-radius:12px}.navigation strong{margin-right:12px}.navigation a{display:inline-block;margin:3px 10px 3px 0;color:#075985;font-weight:700}section{margin:18px 0;padding:22px}section h2{margin:0 0 5px;font-size:23px}section h3{margin:22px 0 8px;font-size:17px}.section-note{margin:0 0 14px;color:var(--muted)}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px}table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}thead th{background:#ecfeff;color:#334155;font-size:12px;text-transform:uppercase;letter-spacing:.04em}.key-value th{width:315px;background:#f8fafc;color:#475569}.path{overflow-wrap:anywhere;word-break:break-word;color:#075985}.muted,.source{color:var(--muted)}code{padding:1px 4px;border-radius:4px;background:#e2e8f0;font:13px ui-monospace,SFMono-Regular,Menlo,monospace}.footer{max-width:1360px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:700px){.container{padding:16px 10px 40px}section{padding:14px}.key-value th{width:42%}th,td{padding:8px}.page-header{padding:28px 18px}}@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0}.page-header p{color:var(--muted)}.container{max-width:none;padding:0}.notice,section{box-shadow:none;break-inside:avoid}.table-wrap{overflow:visible}a{color:inherit;text-decoration:none}}
"""


def render_formatter_html_report(
    *,
    dataset_id: str,
    input_vcf: str | Path,
    output_directory: str | Path,
    report_path: str | Path,
    study_design,
    results: Mapping[str, Mapping[str, Any]],
    schema_reports: Mapping[str, Mapping[str, Any]],
    output_destinations: Mapping[str, str],
    log_path: str | Path,
    resolved_config_path: str | Path,
    completion_manifest_path: str | Path,
    case_count_column: str,
    control_count_column: str,
    variant_id_observations: Mapping[str, Mapping[str, Any]],
    resolved_bcftools: str | Path,
    module,
) -> str:
    """Render one self-contained, detailed formatter evidence report."""
    output_root = Path(output_directory).expanduser().resolve()
    destination = Path(report_path).expanduser().resolve()
    input_rows = _common_input_rows(results)
    ldsc = results.get("ldsc")
    metrics = [
        ("GWAS-VCF variants", _count(input_rows)),
        ("Formatter targets", len(results)),
        (
            "Trait type",
            (
                "Not inferred"
                if study_design is None
                else study_design.trait_type
            ),
        ),
    ]
    if ldsc is not None:
        metrics.extend([
            ("LDSC-input variants", _count(ldsc.get("rows_out"))),
            (
                "LDSC retained",
                _percentage_text(ldsc.get("rows_out"), ldsc.get("rows_in")),
            ),
        ])

    navigation = '<nav class="navigation"><strong>Report sections</strong>%s</nav>' % "".join(
        '<a href="#target-%s">%s</a>'
        % (
            escape(target, quote=True),
            escape(_target_name(target, variant_id_observations)),
        )
        for target in results
    )
    overview = (
        '<section><h2>Target export overview</h2>'
        '<p class="section-note">Retained and excluded percentages use the '
        'extracted GWAS-VCF variant count for each target.</p>%s</section>'
        % _table(
            (
                "Target", "GWAS-VCF variants", "Retained", "Retained %",
                "Excluded", "Excluded %", "Identifier", "Duplicate policy",
                "Final uniqueness",
            ),
            _target_overview_rows(results, variant_id_observations),
        )
    )
    reference_observations = ""
    if variant_id_observations:
        reference_observations = (
            '<h3>External reference identifier observations</h3>%s'
            % _table(
                ("Target", "Identifier type", "Reference variants scanned"),
                [
                    (
                        _target_name(target, variant_id_observations),
                        observation.get("variant_id_type"),
                        observation.get("variants"),
                    )
                    for target, observation in variant_id_observations.items()
                ],
            )
        )
    input_section = (
        '<section><h2>Input and study-design evidence</h2>'
        '<p class="section-note">Study design is inferred only for formatter '
        'targets whose sample-size contract differs between binary and '
        'quantitative traits.</p>%s%s</section>'
        % (
            _path_table(
                (("Input GWAS-VCF", input_vcf),),
                output_directory=output_root,
                report_path=destination,
            )
            + _study_design_html(
                study_design, case_count_column, control_count_column,
            ),
            reference_observations,
        )
    )
    details = "".join(
        _target_details_html(
            target,
            result,
            schema_reports[target],
            output_destinations,
            output_directory=output_root,
            report_path=destination,
            case_count_column=case_count_column,
            control_count_column=control_count_column,
            variant_id_observations=variant_id_observations,
        )
        for target, result in results.items()
    )
    provenance_paths = [
        *sorted(output_destinations.items()),
        ("formatter.log", log_path),
        ("formatter.resolved_configuration", resolved_config_path),
        ("formatter.completion_manifest", completion_manifest_path),
    ]
    output_section = (
        '<section><h2>Complete artifact inventory</h2>'
        '<p class="section-note">Formatter artifacts use paths relative to the '
        'output root, so links remain valid when a completed pipeline stage is '
        'published or the run directory is copied and checksum-validated.</p>%s'
        '</section>'
        % _path_table(
            provenance_paths,
            output_directory=output_root,
            report_path=destination,
        )
    )
    settings = (
        ("Minimum representable raw P", module.minimum_p_value),
        ("VCF include expression", module.vcf_include_expression),
        ("Configured chromosomes", module.chromosomes),
        ("Default identifier type", module.variant_identifiers.default_type),
        (
            "Default duplicate policy",
            module.variant_identifiers.default_duplicate_policy,
        ),
        (
            "Per-target identifier types",
            module.variant_identifiers.target_types,
        ),
        (
            "Per-target duplicate policies",
            module.variant_identifiers.target_duplicate_policies,
        ),
        (
            "LDSC case-fraction aggregation",
            module.ldsc_sample_prevalence.aggregation,
        ),
        (
            "LDSC SNP/allele reference",
            module.ldsc_reference.merge_alleles_file,
        ),
        ("Canonical VCF input contract", module.input_contract.model_dump(mode="json")),
    )
    settings_section = (
        '<section><h2>Resolved formatter policy</h2>'
        '<p class="section-note">These schema-validated settings describe the '
        'effective formatter run. Scientific output schemas are shown within '
        'each target section.</p>%s</section>' % _key_value_table(settings)
    )
    provenance_section = (
        '<section><h2>Execution provenance</h2>%s%s</section>'
        % (
            _key_value_table((("Dataset", dataset_id),)),
            _path_table(
                (
                    ("Output root", output_root),
                    ("Resolved bcftools executable", resolved_bcftools),
                    ("Canonical formatter log", log_path),
                    ("Resolved configuration", resolved_config_path),
                    ("Checksum completion manifest", completion_manifest_path),
                    ("This HTML report", destination),
                ),
                output_directory=output_root,
                report_path=destination,
            ),
        )
    )
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s — PostGWAS formatter report</title>
<style>%s</style>
</head>
<body>
<header class="page-header"><p class="eyebrow">PostGWAS</p><h1>Formatter report</h1><p>%s · GWAS-VCF input evidence, target-specific variant accounting, reference matching, duplicate handling, output schemas, sample-size semantics, and reproducible artifact provenance.</p></header>
<main class="container">
<div class="notice"><strong>Evidence reuse:</strong> this report is rendered from the same validated result metadata used by the terminal summary, canonical log, and completion manifest. It does not reread scientific inputs or recalculate formatter metrics.</div>
%s
%s
%s
%s
%s
%s
%s
</main>
<footer class="footer">PostGWAS formatter report · review the canonical log and checksum completion manifest for machine-readable provenance.</footer>
</body>
</html>
""" % (
        escape(dataset_id),
        _STYLES,
        escape(dataset_id),
        _metric_cards(metrics),
        navigation,
        input_section,
        overview,
        details,
        settings_section,
        output_section + provenance_section,
    )


def write_formatter_html_report(
    destination: str | Path,
    **values,
) -> Path:
    """Atomically write one complete detailed formatter HTML report."""
    return write_html_report(
        render_formatter_html_report(report_path=destination, **values),
        destination,
    )


__all__ = [
    "ldsc_sample_prevalence_screen_fields",
    "render_formatter_html_report",
    "render_formatter_screen_summary",
    "write_formatter_html_report",
]
