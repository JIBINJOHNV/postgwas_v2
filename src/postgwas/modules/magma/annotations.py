"""Validation and post-processing for functional MAGMA annotations."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from postgwas.core.gene_coordinates import read_gene_coordinates
from postgwas.core.input_validation import record_file_validation
from postgwas.core.io.delimiters import open_text
from postgwas.core.io.reports import write_delimited_report, write_text_lines_report
from postgwas.core.paths import require_nonempty_file
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.reference import read_reference_bim_matches


def read_gene_locations(
    gene_location_file: str | Path,
    module_config,
    label: str = "MAGMA gene-location reference",
    *,
    has_header: bool | None = None,
) -> dict[str, tuple[str, int, int, str, str | None]]:
    """Apply MAGMA's chromosome policy to the shared immutable gene reader."""
    input_config = module_config.input
    rows = read_gene_coordinates(
        gene_location_file,
        column_roles=["gene" if role == "gene_id" else role for role in input_config.gene_location_columns],
        delimiter_pattern=input_config.table_delimiter_pattern,
        has_header=input_config.gene_location_has_header if has_header is None else has_header,
        comment_prefix=module_config.annotation_validation.comment_prefix,
        compressed=True,
        label=label,
        error_type=MagmaError,
    )
    invalid_chromosomes = set(input_config.invalid_chromosome_labels)
    for row in rows:
        if row.chromosome in invalid_chromosomes:
            message = "%s has an invalid gene ID, chromosome, interval, or strand for gene %s" % (label, row.gene)
            record_file_validation(gene_location_file, label, status="failed", message=message)
            raise MagmaError(message)
    return {
        row.gene: (row.chromosome, row.start, row.end, row.strand, row.alternate)
        for row in rows
    }


def write_pathway_compatible_gene_locations(
    gene_location_file: str | Path,
    output_file: str | Path,
    module_config,
    logger,
) -> tuple[Path, dict]:
    """Write a positional MAGMA location file keyed by alternate gene IDs."""
    locations = read_gene_locations(gene_location_file, module_config)
    roles = module_config.input.gene_location_columns
    alternate_column = (
        roles.index("alternate_gene_id") + 1
        if "alternate_gene_id" in roles else None
    )
    grouped: dict[str, list[tuple[str, str, int, int, str]]] = {}
    missing_alternate_ids = 0
    for primary, (chromosome, start, end, strand, alternate) in locations.items():
        if alternate is None:
            missing_alternate_ids += 1
            continue
        grouped.setdefault(alternate, []).append(
            (primary, chromosome, start, end, strand)
        )
    if not grouped:
        raise MagmaError(
            "The configured gene-location file contains no non-empty alternate "
            "gene identifiers%s"
            % (
                " in column %d" % alternate_column
                if alternate_column is not None else ""
            )
        )
    ambiguous = {
        alternate: candidates
        for alternate, candidates in grouped.items()
        if len(candidates) > 1
    }
    policy = module_config.gene_sets.alternate_id_duplicate_policy
    if ambiguous and policy == "error":
        raise MagmaError(
            "Gene-location column %d contains %d identifiers assigned to multiple "
            "intervals; examples: %s. Change "
            "gene_sets.alternate_id_duplicate_policy or provide unique identifiers."
            % (
                alternate_column,
                len(ambiguous),
                ", ".join(list(ambiguous)[:5]),
            )
        )
    selected = []
    for alternate, candidates in grouped.items():
        # The configured longest-interval policy is deterministic: primary ID
        # breaks equal-length ties without depending on input row order.
        primary, chromosome, start, end, strand = sorted(
            candidates,
            key=lambda value: (-(value[3] - value[2] + 1), value[0]),
        )[0]
        selected.append((alternate, chromosome, start, end, strand, primary))
    delimiter = module_config.input.output_table_delimiter
    destination = write_text_lines_report(
        [delimiter.join(map(str, record)) for record in selected],
        output_file,
    )
    metrics = {
        "input_rows": len(locations),
        "missing_alternate_id_rows": missing_alternate_ids,
        "unique_alternate_ids": len(grouped),
        "ambiguous_alternate_ids": len(ambiguous),
        "duplicate_rows_removed": sum(
            len(candidates) - 1 for candidates in ambiguous.values()
        ),
        "duplicate_policy": policy,
        "retained_rows": len(selected),
    }
    logger.record(
        "OUTPUT",
        "magma_pathway_compatible_gene_locations",
        path=str(destination),
        **metrics,
    )
    return destination, metrics


def validate_gene_annotation(
    annotation_file: str | Path,
    ld_reference_prefix: str | Path,
    module_config,
    logger,
) -> dict:
    """Validate MAGMA annotation structure and exact BIM identifier overlap."""
    source = require_nonempty_file(
        annotation_file, "MAGMA gene annotation", error_type=MagmaError,
    )
    policy = module_config.annotation_validation
    coordinate = re.compile(policy.coordinate_pattern)
    delimiter = re.compile(module_config.input.table_delimiter_pattern)
    genes: set[str] = set()
    variants: set[str] = set()
    assignments = 0
    excluded_placeholders = 0
    try:
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text or text.startswith(policy.comment_prefix):
                    continue
                fields = delimiter.split(text)
                if len(fields) < 3:
                    raise MagmaError(
                        "MAGMA annotation line %d has fewer than three fields: %s"
                        % (line_number, source)
                    )
                gene, location, *assigned = fields
                match = coordinate.fullmatch(location)
                if match is None:
                    raise MagmaError(
                        "MAGMA annotation line %d has an invalid configured coordinate: %s"
                        % (line_number, location)
                    )
                start = int(match.group("start"))
                end = int(match.group("end"))
                if start < 1 or end < start:
                    raise MagmaError(
                        "MAGMA annotation line %d has an invalid interval: %s"
                        % (line_number, location)
                    )
                if gene in genes:
                    raise MagmaError(
                        "MAGMA annotation repeats gene identifier %s at line %d"
                        % (gene, line_number)
                    )
                usable = [
                    value for value in assigned
                    if value not in policy.invalid_variant_identifiers
                ]
                excluded_placeholders += len(assigned) - len(usable)
                genes.add(gene)
                variants.update(usable)
                assignments += len(usable)
    except (OSError, UnicodeError, ValueError) as exc:
        raise MagmaError("Cannot validate MAGMA annotation %s: %s" % (source, exc)) from exc
    if not genes or not variants:
        raise MagmaError("MAGMA annotation contains no gene-to-variant assignments")

    matched, reference_variants = read_reference_bim_matches(
        ld_reference_prefix, module_config.input, variants,
    )
    matched_variants = len(matched)
    overlap = matched_variants / len(variants)
    minimum = policy.minimum_bim_variant_overlap_fraction
    if overlap < minimum:
        raise MagmaError(
            "Only %d/%d unique annotation variants (%.2f%%) occur in BIM field 2; "
            "the configured minimum is %.2f%%. The annotation and LD reference "
            "must use the same variant identifiers and release."
            % (matched_variants, len(variants), overlap * 100, minimum * 100)
        )
    result = {
        "annotation_file": str(source),
        "genes": len(genes),
        "unique_annotation_variants": len(variants),
        "variant_assignments": assignments,
        "excluded_placeholder_assignments": excluded_placeholders,
        "reference_variants": reference_variants,
        "matched_annotation_variants": matched_variants,
        "bim_overlap_fraction": overlap,
        "minimum_bim_overlap_fraction": minimum,
    }
    logger.record("RESULT", "magma_annotation_validation", **result)
    return result


def prepare_scoped_gene_annotation(
    annotation_file: str | Path,
    output_file: str | Path,
    excluded_genes_file: str | Path,
    summary_file: str | Path,
    module_config,
    analysis_scope: dict,
    logger,
) -> dict:
    """Create a run-owned MAGMA annotation restricted to the configured scope."""
    source = require_nonempty_file(
        annotation_file, "MAGMA gene annotation", error_type=MagmaError,
    )
    policy = module_config.annotation_validation
    coordinate = re.compile(policy.coordinate_pattern)
    delimiter = re.compile(module_config.input.table_delimiter_pattern)
    chromosome_exclusions = set(analysis_scope["exclude_chromosomes"])
    exclude_mhc = analysis_scope.get("exclude_mhc_genes", False)
    mhc = analysis_scope["mhc_region"]
    retained: list[str] = []
    exclusions: list[tuple[str, str, str, int, int]] = []
    total = 0
    report = module_config.exclusion_reporting
    try:
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text or text.startswith(policy.comment_prefix):
                    continue
                fields = delimiter.split(text)
                if len(fields) < 2:
                    raise MagmaError(
                        "MAGMA annotation line %d has fewer than two fields: %s"
                        % (line_number, source)
                    )
                match = coordinate.fullmatch(fields[1])
                if match is None:
                    raise MagmaError(
                        "MAGMA annotation line %d has an invalid configured "
                        "coordinate: %s" % (line_number, fields[1])
                    )
                chromosome = re.sub(
                    module_config.input.chromosome_prefix_pattern,
                    "",
                    match.group("chromosome"),
                ).upper()
                chromosome = module_config.input.chromosome_aliases.get(
                    chromosome, chromosome,
                )
                start = int(match.group("start"))
                end = int(match.group("end"))
                total += 1
                reason = None
                if chromosome in chromosome_exclusions:
                    reason = report.chromosome_reason
                elif (
                    exclude_mhc
                    and mhc is not None
                    and chromosome == mhc["chromosome"]
                    and start <= mhc["end"]
                    and end >= mhc["start"]
                ):
                    reason = report.mhc_reason
                if reason is None:
                    retained.append(text)
                else:
                    exclusions.append((fields[0], reason, chromosome, start, end))
    except (OSError, UnicodeError, ValueError) as exc:
        raise MagmaError(
            "Cannot apply MAGMA gene exclusions to %s: %s" % (source, exc)
        ) from exc
    if not retained:
        raise MagmaError(
            "MAGMA chromosome and MHC policies excluded every annotated unit"
        )
    destination = Path(output_file)
    excluded_destination = Path(excluded_genes_file)
    summary_destination = Path(summary_file)
    for path in (destination, excluded_destination, summary_destination):
        path.parent.mkdir(parents=True, exist_ok=True)
    try:
        destination.write_text("\n".join(retained) + "\n", encoding="utf-8")
        excluded_destination.write_text(
            "%s\t%s\t%s\t%s\t%s\n"
            % (
                module_config.result_schema.gene_id_column,
                report.reason_column,
                module_config.input.chromosome_column,
                report.start_column,
                report.end_column,
            )
            + "".join(
                "%s\t%s\t%s\t%d\t%d\n" % record for record in exclusions
            ),
            encoding="utf-8",
        )
    except OSError as exc:
        raise MagmaError(
            "Cannot write MAGMA scoped annotation or exclusion audit: %s" % exc
        ) from exc
    by_reason = {
        reason: sum(record[1] == reason for record in exclusions)
        for reason in (report.chromosome_reason, report.mhc_reason)
    }
    summary = [
        {
            report.scope_column: report.annotated_units_scope,
            report.reason_column: reason,
            report.excluded_count_column: count,
            report.input_count_column: total,
            report.retained_count_column: len(retained),
        }
        for reason, count in by_reason.items()
    ]
    write_delimited_report(
        summary,
        summary_destination,
        fieldnames=[
            report.scope_column,
            report.reason_column,
            report.excluded_count_column,
            report.input_count_column,
            report.retained_count_column,
        ],
        delimiter=module_config.result_schema.report_delimiter,
        null_value=module_config.result_schema.report_null_value,
    )
    result = {
        "source_annotation": str(source),
        "analysis_annotation": str(destination),
        "input_units": total,
        "retained_units": len(retained),
        "excluded_units": len(exclusions),
        "excluded_by_reason": by_reason,
        "excluded_gene_ids": [record[0] for record in exclusions],
        "excluded_genes_file": str(excluded_destination),
        "exclusion_summary": str(summary_destination),
    }
    logger.record("RESULT", "magma_gene_scope", **result)
    return result


def merge_gene_annotations(
    annotation_files: list[str | Path],
    gene_location_file: str | Path,
    output_file: str | Path,
    module_config,
    logger,
) -> dict:
    """Union nMAGMA SNP assignments using its canonical gene coordinates."""
    policy = module_config.annotation_validation
    input_config = module_config.input
    coordinate = re.compile(policy.coordinate_pattern)
    delimiter = re.compile(input_config.table_delimiter_pattern)
    source_locations = read_gene_locations(
        gene_location_file, module_config, "nMAGMA gene-location reference",
    )
    gene_locations = {
        gene: "%s:%d:%d" % values[:3]
        for gene, values in source_locations.items()
    }

    genes: dict[str, tuple[list[str], set[str]]] = {}
    source_assignments = []
    absent_genes: set[str] = set()
    coordinate_disagreements = 0
    excluded_noncanonical_assignments = 0
    invalid_component_coordinates = 0
    for annotation_file in annotation_files:
        source = require_nonempty_file(
            annotation_file, "MAGMA annotation component", error_type=MagmaError,
        )
        assignments = 0
        placeholders = 0
        source_coordinate_disagreements = 0
        source_absent_genes: set[str] = set()
        source_excluded_noncanonical_assignments = 0
        source_invalid_component_coordinates = 0
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text or text.startswith(policy.comment_prefix):
                    continue
                fields = delimiter.split(text)
                if len(fields) < 3:
                    raise MagmaError(
                        "MAGMA annotation component %s has an invalid record at "
                        "line %d" % (source, line_number)
                    )
                gene, location, *variants = fields
                usable = [
                    variant
                    for variant in variants
                    if variant not in policy.invalid_variant_identifiers
                ]
                placeholders += len(variants) - len(usable)
                match = coordinate.fullmatch(location)
                valid_coordinate = match is not None
                if match is not None:
                    start = int(match.group("start"))
                    end = int(match.group("end"))
                    valid_coordinate = start >= 1 and end >= start
                if not valid_coordinate:
                    source_invalid_component_coordinates += 1
                    invalid_component_coordinates += 1
                if gene not in gene_locations:
                    source_absent_genes.add(gene)
                    absent_genes.add(gene)
                    source_excluded_noncanonical_assignments += len(usable)
                    excluded_noncanonical_assignments += len(usable)
                    continue
                if location != gene_locations[gene]:
                    source_coordinate_disagreements += 1
                    coordinate_disagreements += 1
                if gene not in genes:
                    genes[gene] = ([], set())
                ordered, observed = genes[gene]
                for variant in usable:
                    if variant not in observed:
                        ordered.append(variant)
                        observed.add(variant)
                        assignments += 1
        source_assignments.append(
            {
                "path": str(source),
                "new_unique_assignments": assignments,
                "excluded_placeholder_assignments": placeholders,
                "coordinate_disagreement_records": (
                    source_coordinate_disagreements
                ),
                "genes_absent_from_location_reference": len(source_absent_genes),
                "assignments_excluded_outside_location_reference": (
                    source_excluded_noncanonical_assignments
                ),
                "invalid_component_coordinate_records": (
                    source_invalid_component_coordinates
                ),
            }
        )
    retained_genes = {
        gene: genes[gene]
        for gene in gene_locations
        if gene in genes and genes[gene][0]
    }
    if not retained_genes:
        raise MagmaError("nMAGMA annotation components contain no genes")
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with destination.open("w", encoding="utf-8", newline="\n") as handle:
            for gene, (variants, _) in retained_genes.items():
                handle.write(
                    input_config.output_table_delimiter.join(
                        [gene, gene_locations[gene], *variants]
                    )
                    + "\n"
                )
    except OSError as exc:
        raise MagmaError("Cannot write merged nMAGMA annotation: %s" % exc) from exc
    destination = require_nonempty_file(
        destination, "Merged nMAGMA annotation", error_type=MagmaError,
    )
    result = {
        "annotation_components": source_assignments,
        "gene_location_file": str(Path(gene_location_file).expanduser().resolve()),
        "canonical_gene_locations": len(gene_locations),
        "genes": len(retained_genes),
        "unique_gene_variant_assignments": sum(
            len(values[0]) for values in retained_genes.values()
        ),
        "excluded_placeholder_assignments": sum(
            record["excluded_placeholder_assignments"]
            for record in source_assignments
        ),
        "coordinate_disagreement_records": coordinate_disagreements,
        "genes_absent_from_location_reference": len(absent_genes),
        "assignments_excluded_outside_location_reference": (
            excluded_noncanonical_assignments
        ),
        "invalid_component_coordinate_records": invalid_component_coordinates,
        "output": str(destination),
    }
    logger.record("RESULT", "n_magma_annotation_merge", **result)
    return result


def _regulatory_element_table(
    source: Path,
    schema,
    input_config,
    *,
    include_gene: bool,
) -> pd.DataFrame:
    delimiter = re.compile(
        schema.mapping_delimiter_pattern
        if include_gene
        else schema.location_delimiter_pattern
    )
    has_header = schema.mapping_has_header if include_gene else schema.location_has_header
    required_indexes = [
        schema.element_id_column,
        schema.chromosome_column,
        schema.start_column,
        schema.end_column,
    ]
    if include_gene:
        required_indexes.append(schema.gene_id_column)
        if schema.score_column is not None:
            required_indexes.append(schema.score_column)
    rows = []
    maximum_index = max(required_indexes)
    with open_text(source) as handle:
        for line_number, raw in enumerate(handle, 1):
            if has_header and line_number == 1:
                continue
            text = raw.strip()
            if not text:
                continue
            fields = delimiter.split(text)
            if len(fields) <= maximum_index:
                raise MagmaError(
                    "Element-to-gene mapping line %d has %d fields; configured "
                    "column %d is required."
                    % (line_number, len(fields), maximum_index)
                )
            element = fields[schema.element_id_column].strip()
            chromosome = re.sub(
                input_config.chromosome_prefix_pattern,
                "",
                fields[schema.chromosome_column].strip(),
            ).upper()
            chromosome = input_config.chromosome_aliases.get(
                chromosome, chromosome,
            )
            try:
                start = int(fields[schema.start_column])
                end = int(fields[schema.end_column])
            except ValueError as exc:
                raise MagmaError(
                    "Regulatory-element reference contains a non-integer interval "
                    "at line %d: %s" % (line_number, source)
                ) from exc
            if (
                not element
                or chromosome in input_config.invalid_chromosome_labels
                or start < 1
                or end < start
            ):
                raise MagmaError(
                    "Regulatory-element reference contains an invalid identifier or "
                    "interval at line %d: %s" % (line_number, source)
                )
            gene = fields[schema.gene_id_column].strip() if include_gene else None
            if include_gene and not gene:
                raise MagmaError(
                    "Element-to-gene mapping contains an empty gene at line %d."
                    % line_number
                )
            score = (
                fields[schema.score_column].strip()
                if include_gene and schema.score_column is not None
                else None
            )
            rows.append((element, chromosome, start, end, gene, score))
    if not rows:
        raise MagmaError("Regulatory-element reference contains no usable records")
    frame = pd.DataFrame(
        rows,
        columns=["element", "chromosome", "start", "end", "gene", "source_score"],
    )
    coordinate_counts = frame.groupby("element")[["chromosome", "start", "end"]].nunique()
    conflicts = coordinate_counts.max(axis=1).gt(1)
    if conflicts.any():
        raise MagmaError(
            "Regulatory-element reference assigns conflicting coordinates to: %s"
            % conflicts[conflicts].index[:5].tolist()
        )
    subset = ["element", "gene"] if include_gene else ["element"]
    return frame.drop_duplicates(subset, keep="first")


def map_regulatory_elements_to_genes(
    element_results: pd.DataFrame,
    element_location_file: str | Path,
    mapping_file: str | Path,
    output_file: str | Path,
    module_config,
    logger,
) -> pl.DataFrame:
    """Apply the published chromMAGMA lowest-element-p assignment per gene."""
    schema = module_config.chrom_magma_mapping
    element_column = schema.result_element_column
    p_column = schema.result_p_value_column
    missing = sorted({element_column, p_column} - set(element_results.columns))
    if missing or element_results.empty:
        raise MagmaError(
            "chromMAGMA element results are empty or missing columns: %s"
            % ", ".join(missing or [element_column, p_column])
        )
    elements = element_results.copy()
    elements[element_column] = elements[element_column].astype(str).str.strip()
    elements[p_column] = pd.to_numeric(elements[p_column], errors="coerce")
    invalid = (
        elements[element_column].eq("")
        | elements[p_column].isna()
        | ~np.isfinite(elements[p_column])
        | elements[p_column].lt(0)
        | elements[p_column].gt(1)
    )
    if invalid.any():
        raise MagmaError(
            "chromMAGMA element results contain %d invalid identifiers or p-values"
            % int(invalid.sum())
        )
    mapping_path = require_nonempty_file(
        mapping_file, "chromMAGMA element-to-gene mapping", error_type=MagmaError,
    )
    location_path = require_nonempty_file(
        element_location_file,
        "chromMAGMA regulatory-element locations",
        error_type=MagmaError,
    )
    mapping = _regulatory_element_table(
        mapping_path, schema, module_config.input, include_gene=True,
    )
    locations = _regulatory_element_table(
        location_path, schema, module_config.input, include_gene=False,
    )
    unique_results = set(elements[element_column])
    missing_locations = unique_results - set(locations["element"])
    if missing_locations:
        raise MagmaError(
            "%d tested chromMAGMA elements are absent from the configured "
            "regulatory-element location file; examples: %s"
            % (len(missing_locations), sorted(missing_locations)[:5])
        )
    mapped_results = unique_results & set(mapping["element"])
    overlap = len(mapped_results) / len(unique_results)
    if overlap < schema.minimum_element_mapping_fraction:
        raise MagmaError(
            "Only %d/%d tested regulatory elements (%.2f%%) occur in the "
            "element-to-gene mapping; the configured minimum is %.2f%%."
            % (
                len(mapped_results), len(unique_results), overlap * 100,
                schema.minimum_element_mapping_fraction * 100,
            )
        )
    tested_mapping = mapping[mapping["element"].isin(mapped_results)].merge(
        locations[["element", "chromosome", "start", "end"]],
        on="element",
        how="left",
        suffixes=("_mapping", "_location"),
        validate="many_to_one",
    )
    coordinate_mismatch = (
        tested_mapping["chromosome_mapping"]
        != tested_mapping["chromosome_location"]
    ) | (tested_mapping["start_mapping"] != tested_mapping["start_location"]) | (
        tested_mapping["end_mapping"] != tested_mapping["end_location"]
    )
    if coordinate_mismatch.any():
        raise MagmaError(
            "The element-to-gene mapping and regulatory-element location file "
            "disagree for %d tested elements; examples: %s"
            % (
                int(coordinate_mismatch.sum()),
                tested_mapping.loc[coordinate_mismatch, "element"].head(5).tolist(),
            )
        )
    joined = elements.merge(
        tested_mapping[["element", "gene", "source_score"]],
        left_on=element_column,
        right_on="element",
        how="inner",
        validate="one_to_many",
    )
    counts = joined.groupby("gene")["element"].nunique()
    joined = joined.sort_values(
        ["gene", p_column, "element"], kind="mergesort",
    )
    selected = joined.drop_duplicates("gene", keep="first").copy()
    original_columns = [
        column for column in element_results.columns
        if column not in {element_column, p_column}
    ]
    selected.rename(
        columns={
            column: schema.report_original_column_prefix + column.lower()
            for column in original_columns
        },
        inplace=True,
    )
    output = pd.DataFrame({
        module_config.result_schema.gene_id_column: selected["gene"],
        p_column: selected[p_column],
        schema.report_element_column: selected["element"],
        schema.report_element_count_column: selected["gene"].map(counts),
        schema.report_source_score_column: selected["source_score"],
    })
    for original in original_columns:
        output[schema.report_original_column_prefix + original.lower()] = selected[
            schema.report_original_column_prefix + original.lower()
        ]
    destination = write_delimited_report(
        output.to_dict(orient="records"),
        output_file,
        fieldnames=output.columns.tolist(),
        delimiter=module_config.result_schema.report_delimiter,
        null_value=module_config.result_schema.report_null_value,
    )
    result = {
        "tested_regulatory_elements": len(unique_results),
        "mapped_regulatory_elements": len(mapped_results),
        "element_mapping_fraction": overlap,
        "element_gene_links": len(joined),
        "genes": len(output),
        "assignment_policy": schema.gene_assignment_policy,
        "tie_breaker": schema.tie_breaker,
        "output": str(destination),
    }
    logger.record("RESULT", "chrom_magma_element_to_gene", **result)
    return pl.DataFrame(output.to_dict(orient="list"))


__all__ = [
    "map_regulatory_elements_to_genes",
    "merge_gene_annotations",
    "prepare_scoped_gene_annotation",
    "validate_gene_annotation",
]
