"""Streaming QC and consistent reports for population LD annotation."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Any, Mapping, Sequence

from postgwas.core.io.reports import write_delimited_report, write_html_report
from postgwas.core.processes import run_checked_line_command
from postgwas.core.ui.screen import screen_field, screen_line


@dataclass(frozen=True)
class LDReferenceSummary:
    """Small BED inventory needed for annotation coverage reporting."""

    path: Path
    contigs: frozenset[str]
    labels: frozenset[str]
    block_count: int


SUMMARY_COLUMNS = (
    "record_type",
    "dataset_id",
    "genome_build",
    "population",
    "contig",
    "total_variants",
    "annotated_variants",
    "unassigned_variants",
    "annotation_percent",
    "any_population_annotated",
    "all_populations_annotated",
    "unassigned_on_bed_contigs",
    "unassigned_outside_bed_contigs",
    "bed_blocks",
    "unique_bed_labels",
    "blocks_used",
    "empty_blocks",
    "duplicate_bed_labels",
    "unexpected_annotation_labels",
    "bed_file",
)


def _percentage(count: int, total: int) -> float:
    return 0.0 if total == 0 else 100.0 * int(count) / int(total)


def calculate_ld_annotation_summary(
    vcf: Path,
    *,
    bcftools: str,
    dataset_id: str,
    genome_build: str,
    populations: Sequence[str],
    info_ids: Mapping[str, str],
    references: Mapping[str, LDReferenceSummary],
    expected_variant_count: int,
    vcf_contigs: Sequence[str],
    logger=None,
) -> dict:
    """Count annotation coverage in one checked, constant-memory VCF stream."""
    selected = tuple(str(population) for population in populations)
    annotated_counts = {population: 0 for population in selected}
    observed_labels = {population: set() for population in selected}
    unassigned_by_contig: dict[str, int] = {}
    total_variants = 0
    any_population_annotated = 0
    all_populations_annotated = 0
    expected_fields = len(selected) + 1

    def consume(line: str) -> None:
        nonlocal total_variants
        nonlocal any_population_annotated
        nonlocal all_populations_annotated
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != expected_fields or not fields[0]:
            raise ValueError(
                "LD-annotation QC query returned %d fields; expected %d"
                % (len(fields), expected_fields)
            )
        total_variants += 1
        assigned = []
        for population, value in zip(selected, fields[1:], strict=True):
            present = value not in {"", "."}
            assigned.append(present)
            if present:
                annotated_counts[population] += 1
                observed_labels[population].add(value)
        if any(assigned):
            any_population_annotated += 1
        else:
            contig = fields[0]
            unassigned_by_contig[contig] = unassigned_by_contig.get(contig, 0) + 1
        if all(assigned):
            all_populations_annotated += 1

    query_format = "%CHROM" + "".join(
        "\t%%INFO/%s" % info_ids[population] for population in selected
    ) + "\n"
    streamed_lines = run_checked_line_command(
        [bcftools, "query", "-f", query_format, str(vcf)],
        "Calculating LD-annotation QC summary",
        line_consumer=consume,
        logger=logger,
        error_type=ValueError,
    )
    if streamed_lines != total_variants:
        raise ValueError(
            "LD-annotation QC stream accounting is inconsistent: %d lines, %d variants"
            % (streamed_lines, total_variants)
        )
    if total_variants != int(expected_variant_count):
        raise ValueError(
            "LD-annotation QC counted %d variants; the validated TBI contains %d"
            % (total_variants, expected_variant_count)
        )

    bed_contigs = set().union(
        *(references[population].contigs for population in selected)
    )
    unassigned_on_bed_contigs = sum(
        count
        for contig, count in unassigned_by_contig.items()
        if contig in bed_contigs
    )
    ordered_unassigned = {
        contig: unassigned_by_contig[contig]
        for contig in vcf_contigs
        if contig in unassigned_by_contig
    }
    ordered_unassigned.update(
        (contig, count)
        for contig, count in sorted(unassigned_by_contig.items())
        if contig not in ordered_unassigned
    )

    population_summaries = {}
    for population in selected:
        reference = references[population]
        labels = observed_labels[population]
        annotated = annotated_counts[population]
        population_summaries[population] = {
            "annotated_variants": annotated,
            "unassigned_variants": total_variants - annotated,
            "annotation_percent": _percentage(annotated, total_variants),
            "bed_blocks": reference.block_count,
            "unique_bed_labels": len(reference.labels),
            "blocks_used": len(labels.intersection(reference.labels)),
            "empty_blocks": len(reference.labels.difference(labels)),
            "duplicate_bed_labels": reference.block_count - len(reference.labels),
            "unexpected_annotation_labels": len(labels.difference(reference.labels)),
            "bed_file": str(reference.path),
        }

    fully_unassigned = total_variants - any_population_annotated
    summary = {
        "dataset_id": dataset_id,
        "genome_build": genome_build,
        "total_variants": total_variants,
        "any_population_annotated": any_population_annotated,
        "all_populations_annotated": all_populations_annotated,
        "fully_unassigned_variants": fully_unassigned,
        "unassigned_on_bed_contigs": unassigned_on_bed_contigs,
        "unassigned_outside_bed_contigs": (
            fully_unassigned - unassigned_on_bed_contigs
        ),
        "unassigned_by_contig": ordered_unassigned,
        "bed_contigs": [
            contig for contig in vcf_contigs if contig in bed_contigs
        ],
        "populations": population_summaries,
    }
    if logger is not None and hasattr(logger, "record"):
        logger.record("OBSERVED", "ld_annotation_summary", **summary)
    return summary


def _summary_records(summary: Mapping) -> list[dict]:
    total = int(summary["total_variants"])
    any_annotated = int(summary["any_population_annotated"])
    records = [
        {
            "record_type": "overall",
            "dataset_id": summary["dataset_id"],
            "genome_build": summary["genome_build"],
            "population": "ALL",
            "contig": "ALL",
            "total_variants": total,
            "annotated_variants": any_annotated,
            "unassigned_variants": summary["fully_unassigned_variants"],
            "annotation_percent": "%.6f" % _percentage(any_annotated, total),
            "any_population_annotated": any_annotated,
            "all_populations_annotated": summary[
                "all_populations_annotated"
            ],
            "unassigned_on_bed_contigs": summary[
                "unassigned_on_bed_contigs"
            ],
            "unassigned_outside_bed_contigs": summary[
                "unassigned_outside_bed_contigs"
            ],
        }
    ]
    for population, values in summary["populations"].items():
        records.append(
            {
                "record_type": "population",
                "dataset_id": summary["dataset_id"],
                "genome_build": summary["genome_build"],
                "population": population,
                "contig": "ALL",
                "total_variants": total,
                "annotated_variants": values["annotated_variants"],
                "unassigned_variants": values["unassigned_variants"],
                "annotation_percent": "%.6f" % values["annotation_percent"],
                "bed_blocks": values["bed_blocks"],
                "unique_bed_labels": values["unique_bed_labels"],
                "blocks_used": values["blocks_used"],
                "empty_blocks": values["empty_blocks"],
                "duplicate_bed_labels": values["duplicate_bed_labels"],
                "unexpected_annotation_labels": values[
                    "unexpected_annotation_labels"
                ],
                "bed_file": values["bed_file"],
            }
        )
    for contig, count in summary["unassigned_by_contig"].items():
        records.append(
            {
                "record_type": "unassigned_contig",
                "dataset_id": summary["dataset_id"],
                "genome_build": summary["genome_build"],
                "population": "ALL",
                "contig": contig,
                "unassigned_variants": count,
            }
        )
    return records


def write_ld_annotation_summary(summary: Mapping, destination: Path) -> Path:
    """Write the validated summary as one atomic, machine-readable CSV."""
    return write_delimited_report(
        _summary_records(summary),
        destination,
        fieldnames=SUMMARY_COLUMNS,
        delimiter=",",
        null_value="",
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
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


def _html_cell(value: Any) -> str:
    return escape(_html_value(value))


def _html_path(value: Any) -> str:
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
        escape(href, quote=True),
        label,
    )


def _html_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    *,
    path_columns: Sequence[int] = (),
    empty_message: str = "No rows",
) -> str:
    paths = set(path_columns)
    heading = "".join(
        '<th scope="col">%s</th>' % escape(header) for header in headers
    )
    body_rows = []
    for row in rows:
        cells = "".join(
            "<td>%s</td>"
            % (_html_path(value) if column in paths else _html_cell(value))
            for column, value in enumerate(row)
        )
        body_rows.append("<tr>%s</tr>" % cells)
    body = "".join(body_rows)
    if not rows:
        body = (
            '<tr><td colspan="%d" class="muted">%s</td></tr>'
            % (len(headers), escape(empty_message))
        )
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (heading, body)
    )


def _html_key_value_table(
    rows: Sequence[tuple[str, Any]],
    *,
    path_labels: Sequence[str] = (),
) -> str:
    paths = set(path_labels)
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (
            escape(label),
            _html_path(value) if label in paths else _html_cell(value),
        )
        for label, value in rows
    )
    return (
        '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody>'
        "</table></div>" % body
    )


def _bed_contig_availability(present: int, requested: int) -> str:
    """Describe contig presence without implying interval-level coverage."""
    present = int(present)
    requested = int(requested)
    if requested < 1 or present < 0 or present > requested:
        raise ValueError("BED-contig availability counts are inconsistent")
    if requested == 1:
        if present:
            return "Contig present in the requested BED file"
        return "Contig absent from the requested BED file"
    if present == 0:
        return "Contig absent from all %d requested BED files" % requested
    if present == requested:
        return "Contig present in all %d requested BED files" % requested
    return "Contig present in %d of %d requested BED files" % (
        present,
        requested,
    )


_HTML_STYLES = """
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ee;--canvas:#f4f7fb;--panel:#fff;--brand:#155e75;--brand2:#0891b2;--ok:#166534;--warn:#92400e}*{box-sizing:border-box}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page-header{padding:36px max(24px,calc((100vw - 1320px)/2));color:#fff;background:linear-gradient(135deg,#164e63,#0e7490)}.eyebrow{margin:0 0 5px;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase;opacity:.86}.page-header h1{margin:0;font-size:clamp(28px,4vw,44px);line-height:1.1}.page-header p{max-width:900px;margin:12px 0 0;color:#cffafe}.container{max-width:1320px;margin:0 auto;padding:26px 24px 58px}.notice,.callout,section{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:0 7px 24px rgba(15,23,42,.05)}.notice{margin-bottom:20px;padding:14px 16px;border-left:4px solid var(--brand2)}.callout{margin:16px 0 0;padding:14px 16px;border-left:4px solid var(--warn)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:0 0 22px}.metric{padding:15px;border:1px solid var(--line);border-radius:11px;background:#fff}.metric-label{font-size:12px;font-weight:800;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.metric-value{margin-top:4px;font-size:22px;font-weight:780;overflow-wrap:anywhere}.metric-detail{margin-top:2px;color:var(--muted);font-size:13px}section{margin:18px 0;padding:22px}section h2{margin:0 0 5px;font-size:21px}.section-note{margin:0 0 14px;color:var(--muted)}.definition-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.definition{padding:16px;border:1px solid var(--line);border-radius:10px;background:#f8fafc}.definition h3{margin:0 0 7px;font-size:16px}.definition p{margin:0}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px}table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}thead th{background:#ecfeff;color:#334155;font-size:12px;text-transform:uppercase;letter-spacing:.04em}.key-value th{width:290px;background:#f8fafc;color:#475569}.path{overflow-wrap:anywhere;word-break:break-word;color:#075985}.muted{color:var(--muted)}.footer{max-width:1320px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:900px){.definition-grid{grid-template-columns:1fr}}@media(max-width:700px){.container{padding:16px 10px 40px}section{padding:14px}.key-value th{width:42%}th,td{padding:8px}.page-header{padding:28px 18px}}@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0}.page-header p{color:var(--muted)}.container{max-width:none;padding:0}.notice,section{box-shadow:none;break-inside:avoid}.table-wrap{overflow:visible}a{color:inherit;text-decoration:none}}
"""


def render_ld_annotation_html(
    summary: Mapping,
    *,
    input_vcf: Path,
    vcf_contigs: Sequence[str],
    info_ids: Mapping[str, str],
    references: Mapping[str, LDReferenceSummary],
    threads: int,
    bcftools: str,
    outputs: Mapping[str, Any],
) -> str:
    """Render one self-contained report from the validated coverage evidence."""
    total = int(summary["total_variants"])
    any_annotated = int(summary["any_population_annotated"])
    all_annotated = int(summary["all_populations_annotated"])
    fully_unassigned = int(summary["fully_unassigned_variants"])
    unassigned_on_bed_contigs = int(summary["unassigned_on_bed_contigs"])
    unassigned_outside_bed_contigs = int(
        summary["unassigned_outside_bed_contigs"]
    )
    cards = (
        (
            "Total VCF variants",
            f"{total:,}",
            "Records in the validated annotated VCF",
        ),
        (
            "Assigned in any population",
            f"{any_annotated:,}",
            "%.2f%% of all variants" % _percentage(any_annotated, total),
        ),
        (
            "Assigned in every population",
            f"{all_annotated:,}",
            "%.2f%% of all variants" % _percentage(all_annotated, total),
        ),
        (
            "No assignment in any population",
            f"{fully_unassigned:,}",
            "%.2f%% of all variants" % _percentage(fully_unassigned, total),
        ),
        (
            "No assignment on contigs present in BEDs",
            f"{unassigned_on_bed_contigs:,}",
            "%.2f%% of fully unassigned"
            % _percentage(unassigned_on_bed_contigs, fully_unassigned),
        ),
        (
            "No assignment on contigs absent from all BEDs",
            f"{unassigned_outside_bed_contigs:,}",
            "%.2f%% of fully unassigned"
            % _percentage(unassigned_outside_bed_contigs, fully_unassigned),
        ),
    )
    metrics = '<div class="metrics">%s</div>' % "".join(
        '<div class="metric"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div><div class="metric-detail">%s</div>'
        "</div>"
        % (escape(label), escape(value), escape(detail))
        for label, value, detail in cards
    )

    population_rows = []
    for population, values in summary["populations"].items():
        population_rows.append(
            (
                population,
                info_ids.get(population),
                values["annotated_variants"],
                "%.2f%%" % float(values["annotation_percent"]),
                values["unassigned_variants"],
                values["blocks_used"],
                values["empty_blocks"],
                values["unexpected_annotation_labels"],
            )
        )
    population_table = _html_table(
        (
            "Population",
            "VCF INFO field",
            "Variants with block assignment",
            "Annotation rate",
            "Variants without block assignment",
            "Reference blocks used",
            "Reference blocks with no variants",
            "Labels not found in BED",
        ),
        population_rows,
    )

    reference_rows = []
    for population, values in summary["populations"].items():
        reference = references[population]
        ordered_contigs = [
            contig for contig in vcf_contigs if contig in reference.contigs
        ]
        ordered_contigs.extend(
            sorted(reference.contigs.difference(ordered_contigs))
        )
        reference_rows.append(
            (
                population,
                values["bed_file"],
                values["bed_blocks"],
                values["unique_bed_labels"],
                values["duplicate_bed_labels"],
                len(reference.contigs),
                tuple(ordered_contigs),
            )
        )
    reference_table = _html_table(
        (
            "Population",
            "BED file",
            "BED interval rows",
            "Unique block labels",
            "Duplicate labels",
            "BED contig count",
            "BED contigs",
        ),
        reference_rows,
        path_columns=(1,),
    )

    requested_bed_count = len(references)
    unassigned_rows = [
        (
            contig,
            count,
            "%.2f%%" % _percentage(count, fully_unassigned),
            _bed_contig_availability(
                sum(
                    contig in reference.contigs
                    for reference in references.values()
                ),
                requested_bed_count,
            ),
        )
        for contig, count in summary["unassigned_by_contig"].items()
    ]
    unassigned_table = _html_table(
        (
            "Contig",
            "Variants with no LD-block assignment",
            "Share of all unassigned variants",
            "LD-block BED availability",
        ),
        unassigned_rows,
        empty_message=(
            "Every output-VCF variant is assigned in at least one requested "
            "population."
        ),
    )

    output_rows = (
        ("Annotated VCF", outputs.get("annotated_vcf")),
        ("VCF index", outputs.get("annotated_index")),
        ("Summary CSV", outputs.get("summary_file")),
        ("Detailed HTML report", outputs.get("html_report")),
        ("Canonical log", outputs.get("log_file")),
    )
    output_table = _html_key_value_table(
        output_rows,
        path_labels=tuple(label for label, _value in output_rows),
    )
    run_table = _html_key_value_table(
        (
            ("Dataset", summary["dataset_id"]),
            ("Genome build declared by VCF", summary["genome_build"]),
            ("Input VCF", input_vcf),
            ("Declared VCF contigs", len(vcf_contigs)),
            ("Contig labels", tuple(vcf_contigs)),
            ("Requested populations", tuple(summary["populations"])),
            (
                "Annotation behavior",
                "Non-filtering; unmatched variants remain in the output VCF",
            ),
            (
                "VCF/BED contig matching",
                "Exact compatible labels required; no automatic chr-prefix conversion",
            ),
            ("Compression threads", int(threads)),
            ("bcftools executable", bcftools),
        ),
        path_labels=("Input VCF", "bcftools executable"),
    )
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s — PostGWAS LD-block annotation report</title>
<style>%s</style>
</head>
<body>
<header class="page-header"><p class="eyebrow">PostGWAS</p><h1>LD-block annotation report</h1><p>%s · %s · population-specific coverage, block use, unassigned-variant distribution, and run provenance.</p></header>
<main class="container">
<div class="notice"><strong>Evidence reuse:</strong> this report uses the same validated, constant-memory coverage summary written to the terminal, CSV, and canonical log. Rendering the HTML does not rescan the VCF or calculate a second set of scientific metrics.</div>
%s
<section><h2>How to interpret annotation coverage</h2><p class="section-note">LD annotation does not remove variants: unmatched variants remain in the output VCF with the requested LD-block INFO fields empty. Population-specific counts overlap and must not be added together.</p><div class="definition-grid"><article class="definition"><h3>Assigned in any population</h3><p>The union: a variant has an LD-block label for at least one requested population.</p></article><article class="definition"><h3>Assigned in every population</h3><p>The intersection: a variant has an LD-block label for every requested population.</p></article><article class="definition"><h3>No assignment in any population</h3><p>Every requested population LD-block INFO field is empty for the variant. Its contig-level distribution is shown below.</p></article></div></section>
<section><h2>Population coverage and block use</h2><p class="section-note">Each row refers only to the named population. A reference block is used when it contains at least one output-VCF variant; a reference block with no variants is present in the BED but contains none of the output-VCF variants.</p>%s</section>
<section><h2>Population LD-block reference inventory</h2><p class="section-note">These files passed LD-annotation file and contig-name preflight and are listed for provenance. BED content alone cannot prove its genome build or population origin; verify each resource's provenance independently.</p>%s</section>
<section><h2>Variants with no LD-block assignment in any requested population</h2><p class="section-note">Each row counts output-VCF variants for which every requested population LD-block INFO field is empty. The rows sum to <strong>%s</strong>. LD-block BED availability describes whether the contig occurs in each reference file, not whether a listed variant overlaps an interval.</p>%s</section>
<section><h2>Run and input provenance</h2><p class="section-note">The genome build is the declaration read from the VCF header. BED files do not encode enough provenance for PostGWAS to independently prove their genome build.</p>%s</section>
<section><h2>Validated outputs</h2>%s</section>
<div class="callout"><strong>Scientific interpretation:</strong> LD blocks are population- and genome-build-specific approximations of LD structure. A block label is not a causal locus, and an empty or unassigned interval is not evidence that association is absent.</div>
</main>
<footer class="footer">PostGWAS LD-block annotation report · use the linked CSV and canonical log for machine-readable evidence and complete command provenance.</footer>
</body>
</html>
""" % (
        escape(str(summary["dataset_id"])),
        _HTML_STYLES,
        escape(str(summary["dataset_id"])),
        escape(str(summary["genome_build"])),
        metrics,
        population_table,
        reference_table,
        escape(f"{fully_unassigned:,}"),
        unassigned_table,
        run_table,
        output_table,
    )


def write_ld_annotation_html_report(
    summary: Mapping,
    destination: Path,
    *,
    input_vcf: Path,
    vcf_contigs: Sequence[str],
    info_ids: Mapping[str, str],
    references: Mapping[str, LDReferenceSummary],
    threads: int,
    bcftools: str,
    outputs: Mapping[str, Any],
) -> Path:
    """Write the detailed, self-contained report atomically."""
    document = render_ld_annotation_html(
        summary,
        input_vcf=input_vcf,
        vcf_contigs=vcf_contigs,
        info_ids=info_ids,
        references=references,
        threads=threads,
        bcftools=bcftools,
        outputs=outputs,
    )
    return write_html_report(document, destination)


def render_ld_annotation_summary(
    summary: Mapping,
    summary_file: Path,
    *,
    label_width: int,
    separator_column: int | None = None,
    annotated_vcf: Path | None = None,
    annotated_index: Path | None = None,
    html_report: Path | None = None,
    log_file: Path | None = None,
) -> str:
    """Render the completed QC findings with the shared MAGMA-style fields."""
    total = int(summary["total_variants"])
    fully_unassigned = int(summary["fully_unassigned_variants"])
    summary_indent = 6
    lines = [
        "",
        screen_line("analysis", "LD-block annotation summary", indent=2),
        screen_field(
            "info", "Dataset", summary["dataset_id"], indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "genetic", "Genome build", summary["genome_build"],
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "count", "Total variants", f"{total:,}", indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "genetic",
            "Assigned in any population",
            "%s · %.2f%%"
            % (
                f"{int(summary['any_population_annotated']):,}",
                _percentage(summary["any_population_annotated"], total),
            ),
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "genetic",
            "Assigned in every population",
            "%s · %.2f%%"
            % (
                f"{int(summary['all_populations_annotated']):,}",
                _percentage(summary["all_populations_annotated"], total),
            ),
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "warning" if fully_unassigned else "success",
            "Unassigned in every population",
            "%s · %.2f%%"
            % (f"{fully_unassigned:,}", _percentage(fully_unassigned, total)),
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
        screen_field(
            "count",
            "Unassigned on BED contigs",
            f"{int(summary['unassigned_on_bed_contigs']):,}",
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        ),
    ]
    bed_contigs = set(summary["bed_contigs"])
    for contig, count in summary["unassigned_by_contig"].items():
        if contig in bed_contigs:
            continue
        lines.append(
            screen_field(
                "warning",
                "Contig %s without LD blocks" % contig,
                f"{int(count):,} unassigned variants",
                indent=summary_indent,
                label_width=label_width,
                separator_column=separator_column,
            )
        )
    lines.extend(
        [
            "",
            screen_line(
                "analysis", "Population coverage", indent=summary_indent
            ),
        ]
    )
    for population_index, (population, values) in enumerate(
        summary["populations"].items()
    ):
        if population_index:
            lines.append("")
        lines.extend(
            [
                screen_field(
                    "genetic",
                    "%s annotated" % population,
                    "%s/%s (%.2f%%)"
                    % (
                        f"{int(values['annotated_variants']):,}",
                        f"{total:,}",
                        float(values["annotation_percent"]),
                    ),
                    indent=summary_indent,
                    label_width=label_width,
                    separator_column=separator_column,
                ),
                screen_field(
                    "warning" if values["unassigned_variants"] else "success",
                    "%s unassigned" % population,
                    "%s/%s (%.2f%%)"
                    % (
                        f"{int(values['unassigned_variants']):,}",
                        f"{total:,}",
                        _percentage(values["unassigned_variants"], total),
                    ),
                    indent=summary_indent,
                    label_width=label_width,
                    separator_column=separator_column,
                ),
                screen_field(
                    "count",
                    "%s LD blocks" % population,
                    "%s used / %s total; %s empty"
                    % (
                        f"{int(values['blocks_used']):,}",
                        f"{int(values['bed_blocks']):,}",
                        f"{int(values['empty_blocks']):,}",
                    ),
                    indent=summary_indent,
                    label_width=label_width,
                    separator_column=separator_column,
                ),
            ]
        )
        if values["duplicate_bed_labels"]:
            lines.append(
                screen_field(
                    "warning",
                    "%s duplicate BED labels" % population,
                    f"{int(values['duplicate_bed_labels']):,}",
                    indent=summary_indent,
                    label_width=label_width,
                    separator_column=separator_column,
                )
            )
        if values["unexpected_annotation_labels"]:
            lines.append(
                screen_field(
                    "warning",
                    "%s labels absent from BED" % population,
                    f"{int(values['unexpected_annotation_labels']):,}",
                    indent=summary_indent,
                    label_width=label_width,
                    separator_column=separator_column,
                )
            )
    lines.extend(
        [
            "",
            screen_line("analysis", "Saved outputs", indent=summary_indent),
        ]
    )
    if annotated_vcf is not None:
        lines.append(
            screen_field(
                "success", "Annotated VCF", annotated_vcf.name,
                indent=summary_indent,
                label_width=label_width,
                separator_column=separator_column,
            )
        )
    if annotated_index is not None:
        lines.append(
            screen_field(
                "success", "VCF index", annotated_index.name,
                indent=summary_indent,
                label_width=label_width,
                separator_column=separator_column,
            )
        )
    lines.append(
        screen_field(
            "success", "Summary CSV", summary_file.name,
            indent=summary_indent,
            label_width=label_width,
            separator_column=separator_column,
        )
    )
    if html_report is not None:
        lines.append(
            screen_field(
                "success", "Detailed HTML report", html_report.name,
                indent=summary_indent,
                label_width=label_width,
                separator_column=separator_column,
            )
        )
    if log_file is not None:
        lines.append(
            screen_field(
                "success", "Detailed log", log_file.name,
                indent=summary_indent,
                label_width=label_width,
                separator_column=separator_column,
            )
        )
    lines.extend(("", ""))
    return "\n".join(lines)


__all__ = [
    "LDReferenceSummary",
    "SUMMARY_COLUMNS",
    "calculate_ld_annotation_summary",
    "render_ld_annotation_html",
    "render_ld_annotation_summary",
    "write_ld_annotation_html_report",
    "write_ld_annotation_summary",
]
