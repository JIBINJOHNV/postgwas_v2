"""Validated MAGMA gene and competitive gene-set analysis.

MAGMA protocol invariants retained in code are limited to its documented output
suffixes, metadata marker, and command modifiers. Runtime policies, paths,
thresholds, columns, and resource limits come from the canonical MAGMA YAML
configuration.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from postgwas.core.io.delimiters import open_text
from postgwas.core.io.reports import write_delimited_report
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.processes import run_checked_command
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magma.annotations import (
    map_regulatory_elements_to_genes,
    merge_gene_annotations,
    read_gene_locations,
    validate_gene_annotation,
)
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.reference import read_reference_bim_matches


# These suffixes are fixed by MAGMA's public file contract, not PostGWAS policy.
_ANNOTATION_SUFFIX = ".genes.annot"
_GENE_RAW_SUFFIX = ".genes.raw"
_GENE_RESULT_SUFFIX = ".genes.out"
_GENE_SET_SUFFIX = ".gsa.out"
# MAGMA prefixes metadata rows in its human-readable result tables with '#'.
_RESULT_COMMENT_PREFIX = "#"


@dataclass(frozen=True)
class MagmaPreflight:
    """Validated values that must be available before staging is created."""

    snp_location_file: Path
    p_value_file: Path
    ld_reference_prefix: str
    executable: str
    version: str
    mapping_gene_sets: dict[str, Path | None]


def _tool_output(prefix: str | Path, suffix: str) -> Path:
    return Path("%s%s" % (prefix, suffix))


def resolve_magma_output_paths(
    output_directory: str | Path,
    dataset_id: str,
    module_config,
    analysis_name: str | None = None,
) -> dict[str, Path]:
    """Resolve shared and mapping-specific outputs from canonical patterns."""
    output = Path(output_directory).expanduser().resolve()
    layout = module_config.output_layout
    selected = analysis_name or module_config.mapping.primary
    validate_filename_component(selected, "mapping analysis", error_type=MagmaError)
    paths = {
        "harmonised_p_values": configured_output_path(
            output,
            layout.harmonised_p_values,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "harmonised_snp_locations": configured_output_path(
            output,
            layout.harmonised_snp_locations,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "mapping_comparison": configured_output_path(
            output,
            layout.mapping_comparison,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "annotation_prefix": configured_output_path(
            output,
            layout.annotation_prefix,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "component_annotation_prefix": configured_output_path(
            output,
            layout.component_annotation_prefix,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "gene_prefix": configured_output_path(
            output,
            layout.gene_result_prefix,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
            upstream=module_config.gene_window_upstream_kb,
            downstream=module_config.gene_window_downstream_kb,
        ),
        "gene_batch_prefix": configured_output_path(
            output,
            layout.gene_batch_prefix,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
            upstream=module_config.gene_window_upstream_kb,
            downstream=module_config.gene_window_downstream_kb,
        ),
        "gene_set_prefix": configured_output_path(
            output,
            layout.gene_set_prefix,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "corrected_genes": configured_output_path(
            output,
            layout.corrected_genes,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "corrected_gene_sets": configured_output_path(
            output,
            layout.corrected_gene_sets,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "annotated_gene_sets": configured_output_path(
            output,
            layout.annotated_gene_sets,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "prepared_gene_sets": configured_output_path(
            output,
            layout.prepared_gene_sets,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "chrom_magma_genes": configured_output_path(
            output,
            layout.chrom_magma_genes,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
    }
    paths.update(
        {
            "gene_annotation": _tool_output(
                paths["annotation_prefix"], _ANNOTATION_SUFFIX,
            ),
            "component_gene_annotation": _tool_output(
                paths["component_annotation_prefix"], _ANNOTATION_SUFFIX,
            ),
            "genes_raw": _tool_output(paths["gene_prefix"], _GENE_RAW_SUFFIX),
            "genes_out": _tool_output(paths["gene_prefix"], _GENE_RESULT_SUFFIX),
            "gene_sets_raw": _tool_output(paths["gene_set_prefix"], _GENE_SET_SUFFIX),
        }
    )
    return paths


def _normalize_chromosome(values: pd.Series, input_config) -> pd.Series:
    normalized = (
        values.astype(str)
        .str.strip()
        .str.replace(input_config.chromosome_prefix_pattern, "", regex=True)
        .str.upper()
    )
    return normalized.replace(input_config.chromosome_aliases)


def _validate_p_values(
    values: pd.Series,
    source: str | Path,
    *,
    allow_zero: bool,
) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    lower_invalid = numeric < 0 if allow_zero else numeric <= 0
    invalid = numeric.isna() | ~np.isfinite(numeric) | lower_invalid | (numeric > 1)
    if invalid.any():
        interval = "[0, 1]" if allow_zero else "(0, 1]"
        raise MagmaError(
            "P-value column in %s contains %d invalid values; examples: %s. "
            "Values must be finite and within %s."
            % (source, int(invalid.sum()), values[invalid].head(5).tolist(), interval)
        )
    return numeric


def _read_table(
    path: str | Path,
    label: str,
    delimiter_pattern: str,
    *,
    comment_prefix: str | None = None,
    column_types: dict[str, str] | None = None,
) -> pd.DataFrame:
    file_path = require_nonempty_file(path, label, error_type=MagmaError)
    try:
        return pd.read_csv(
            file_path,
            sep=delimiter_pattern,
            comment=comment_prefix,
            engine="c",
            dtype=column_types,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise MagmaError("Cannot parse %s %s: %s" % (label, file_path, exc)) from exc


def resolve_gene_set_identifiers(
    location_file: str | Path,
    gene_sets: pl.DataFrame,
    minimum_overlap: float,
    module_config,
) -> tuple[pl.DataFrame, set[str], dict]:
    """Validate six-column locations and use column-six IDs only when needed."""
    locations = read_gene_locations(location_file, module_config)
    primary_ids = set(locations)
    alternate_ids: dict[str, list[str]] = {}
    for primary, (*_, alternate) in locations.items():
        if alternate is not None:
            alternate_ids.setdefault(alternate, []).append(primary)
    input_column = module_config.result_schema.report_input_genes_column
    genes = set().union(
        *(set(value.split(",")) for value in gene_sets[input_column].to_list())
    )
    primary_matches = genes & primary_ids
    primary_comparison = min(len(primary_ids), len(genes))
    primary_fraction = (
        len(primary_matches) / primary_comparison if primary_comparison else 0.0
    )
    if primary_fraction >= minimum_overlap:
        return gene_sets, primary_ids, {
            "identifier_source": "primary_gene_id",
            "input_unique_ids": len(genes),
            "matched_unique_ids": len(primary_matches),
            "comparison_unique_ids": primary_comparison,
            "match_fraction": primary_fraction,
            "translated_unique_ids": 0,
            "one_to_many_identifiers": 0,
        }

    translated = (genes - primary_ids) & set(alternate_ids)
    records = []
    for record in gene_sets.iter_rows(named=True):
        resolved = []
        for gene in record[input_column].split(","):
            resolved.extend(
                [gene] if gene in primary_ids else alternate_ids.get(gene, [gene])
            )
        record[input_column] = ",".join(dict.fromkeys(resolved))
        records.append(record)
    resolved_gene_sets = pl.DataFrame(records, schema=gene_sets.schema)
    effective_ids = set().union(
        *(set(value.split(",")) for value in resolved_gene_sets[input_column].to_list())
    )
    matched = effective_ids & primary_ids
    comparison_size = min(len(primary_ids), len(effective_ids))
    match_fraction = len(matched) / comparison_size if comparison_size else 0.0
    if match_fraction < minimum_overlap:
        raise MagmaError(
            "Gene-set identifiers match neither primary gene-location column 1 "
            "(%d/%d; %.2f%%) nor the effective column-1 IDs after translating "
            "column 6 (%d/%d; %.2f%%); the configured minimum is %.2f%%."
            % (
                len(primary_matches), primary_comparison, primary_fraction * 100,
                len(matched), comparison_size, match_fraction * 100,
                minimum_overlap * 100,
            )
        )
    return resolved_gene_sets, primary_ids, {
        "identifier_source": "alternate_gene_id",
        "input_unique_ids": len(genes),
        "matched_unique_ids": len(matched),
        "comparison_unique_ids": comparison_size,
        "match_fraction": match_fraction,
        "translated_unique_ids": len(translated),
        "one_to_many_identifiers": sum(
            len(alternate_ids[gene]) > 1 for gene in translated
        ),
    }


def _gene_ids_from_table(
    table_file: str | Path,
    label: str,
    delimiter_pattern: str,
    identifier_column: int,
    *,
    has_header: bool,
    comment_prefix: str,
) -> set[str]:
    """Read normalized IDs from a configured reference-table column."""
    source = require_nonempty_file(table_file, label, error_type=MagmaError)
    delimiter = re.compile(delimiter_pattern)
    identifiers: set[str] = set()
    header_pending = has_header
    with open_text(source) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(comment_prefix):
                continue
            if header_pending:
                header_pending = False
                continue
            fields = delimiter.split(text)
            if identifier_column >= len(fields):
                raise MagmaError(
                    "%s line %d has %d columns; configured gene-ID column is %d"
                    % (label, line_number, len(fields), identifier_column + 1)
                )
            identifier = _normalize_gene_id(fields[identifier_column])
            if identifier is not None:
                identifiers.add(identifier)
    if not identifiers:
        raise MagmaError("%s contains no gene identifiers" % label)
    return identifiers


def validate_gene_set_identifier_overlap(
    reference_gene_ids: set[str],
    gene_sets: pl.DataFrame,
    minimum_overlap: float | None,
    input_genes_column: str,
) -> dict:
    """Measure gene-set overlap and optionally enforce a compatibility floor."""
    gene_set_values = [
        set(value.split(",")) for value in gene_sets[input_genes_column].to_list()
    ]
    set_ids = set().union(*gene_set_values)
    overlap = reference_gene_ids & set_ids
    comparison_size = min(len(reference_gene_ids), len(set_ids))
    fraction = len(overlap) / comparison_size if comparison_size else 0.0
    usable_sets = sum(
        not reference_gene_ids.isdisjoint(values) for values in gene_set_values
    )
    result = {
        "reference_unique_ids": len(reference_gene_ids),
        "geneset_unique_ids": len(set_ids),
        "overlapping_unique_ids": len(overlap),
        "comparison_unique_ids": comparison_size,
        "overlap_fraction_of_smaller_universe": fraction,
        "minimum_required_overlap": minimum_overlap,
        "total_gene_sets": gene_sets.height,
        "gene_sets_with_at_least_one_reference_gene": usable_sets,
    }
    if minimum_overlap is not None and fraction < minimum_overlap:
        raise MagmaError(
            "Mapping-reference and gene-set identifiers overlap for %d/%d "
            "unique IDs (%.2f%%); the configured minimum is %.2f%%. Use one "
            "gene identifier system."
            % (len(overlap), comparison_size, fraction * 100, minimum_overlap * 100)
        )
    return result


def prepare_magma_variant_inputs(
    pval_file: str | Path,
    snp_loc_file: str | Path,
    ld_ref: str | Path,
    pval_output: str | Path,
    location_output: str | Path,
    module_config,
    logger,
) -> dict:
    """Prepare MAGMA tables and optionally retain their exact BIM intersection."""
    columns = module_config.input
    policy = module_config.snp_harmonisation
    resolve_variants = policy.resolve_variants_to_reference
    pval = _read_table(
        pval_file,
        "MAGMA p-value file",
        columns.table_delimiter_pattern,
        column_types={
            columns.variant_id_column: "string",
            columns.p_value_column: "string",
            columns.sample_size_column: "string",
        },
    )
    locations = _read_table(
        snp_loc_file,
        "MAGMA SNP-location file",
        columns.table_delimiter_pattern,
        column_types={
            columns.variant_id_column: "string",
            columns.chromosome_column: "string",
            columns.position_column: "string",
            columns.reference_allele_column: "string",
            columns.alternate_allele_column: "string",
        },
    )
    required_pval = {
        columns.variant_id_column,
        columns.p_value_column,
        columns.sample_size_column,
    }
    required_location = {
        columns.variant_id_column,
        columns.chromosome_column,
        columns.position_column,
        columns.reference_allele_column,
        columns.alternate_allele_column,
    }
    missing_pval = sorted(required_pval - set(pval.columns))
    missing_location = sorted(required_location - set(locations.columns))
    if missing_pval or missing_location:
        detail = []
        if missing_pval:
            detail.append("p-value table: %s" % ", ".join(missing_pval))
        if missing_location:
            detail.append("location table: %s" % ", ".join(missing_location))
        raise MagmaError("Required MAGMA input columns are missing from " + "; ".join(detail))
    if pval.empty:
        raise MagmaError("MAGMA p-value file contains no variants")

    variant = columns.variant_id_column
    original_columns = pval.columns.tolist()
    missing_pval_ids = pval[variant].isna()
    pval[variant] = pval[variant].astype(str).str.strip()
    if missing_pval_ids.any() or pval[variant].eq("").any():
        raise MagmaError("MAGMA p-value input contains empty variant identifiers")
    pval["_P_NUMERIC"] = _validate_p_values(
        pval[columns.p_value_column], pval_file, allow_zero=False,
    )
    sample_sizes = pd.to_numeric(pval[columns.sample_size_column], errors="coerce")
    invalid_n = sample_sizes.isna() | ~np.isfinite(sample_sizes) | (sample_sizes <= 0)
    if invalid_n.any():
        raise MagmaError(
            "Sample-size column %s contains %d missing, non-finite or non-positive values"
            % (columns.sample_size_column, int(invalid_n.sum()))
        )
    pval[columns.sample_size_column] = sample_sizes
    pval["_ROW_ORDER"] = np.arange(len(pval))

    location_columns = [
        variant, columns.chromosome_column, columns.position_column,
    ]
    allele_columns = [columns.reference_allele_column, columns.alternate_allele_column]
    location_columns.extend(allele_columns)
    locations = locations[location_columns].copy()
    missing_location_ids = locations[variant].isna()
    locations[variant] = locations[variant].astype(str).str.strip()
    if missing_location_ids.any() or locations[variant].eq("").any():
        raise MagmaError("MAGMA SNP-location input contains empty variant identifiers")
    raw_chromosomes = locations[columns.chromosome_column].astype("string").str.strip()
    invalid_chromosomes = raw_chromosomes.isna() | raw_chromosomes.isin(
        columns.invalid_chromosome_labels
    )
    if invalid_chromosomes.any():
        raise MagmaError(
            "MAGMA SNP-location input contains %d invalid chromosome labels."
            % int(invalid_chromosomes.sum())
        )
    locations["CHR_NORM"] = _normalize_chromosome(
        locations[columns.chromosome_column], columns,
    )
    positions = pd.to_numeric(
        locations[columns.position_column], errors="coerce",
    )
    invalid_positions = positions.isna() | positions.le(0) | positions.mod(1).ne(0)
    if invalid_positions.any():
        raise MagmaError(
            "MAGMA SNP-location input contains missing, non-positive, or non-integer "
            "positions"
        )
    locations["BP_NORM"] = positions.astype("int64")
    missing_alleles = locations[allele_columns].isna().any(axis=1)
    first = locations[columns.reference_allele_column].astype(str).str.upper().str.strip()
    second = locations[columns.alternate_allele_column].astype(str).str.upper().str.strip()
    if (
        missing_alleles.any()
        or first.isin(columns.invalid_allele_labels).any()
        or second.isin(columns.invalid_allele_labels).any()
    ):
        raise MagmaError("MAGMA SNP-location input contains invalid REF or ALT alleles")
    locations["REF_NORM"] = first
    locations["ALT_NORM"] = second
    locations["ALLELE_KEY"] = ["|".join(sorted(pair)) for pair in zip(first, second)]
    conflicts = (
        locations.groupby(variant)[["CHR_NORM", "BP_NORM", "ALLELE_KEY"]]
        .nunique()
        .max(axis=1)
        .gt(1)
    )
    if conflicts.any():
        raise MagmaError(
            "SNP-location input assigns conflicting coordinates or allele pairs "
            "to the same ID; examples: %s"
            % conflicts[conflicts].index[:5].tolist()
        )
    locations = locations.drop_duplicates(variant, keep="first")

    working = pval.merge(
        locations[[variant, *[name for name in locations.columns if name != variant]]],
        on=variant,
        how="left",
        validate="many_to_one",
    )
    missing_location_rows = working["BP_NORM"].isna()
    if missing_location_rows.any():
        raise MagmaError(
            "%d p-value rows have no matching SNP-location record, so their "
            "chromosome, position, and alleles cannot be validated; examples: %s"
            % (
                int(missing_location_rows.sum()),
                working.loc[missing_location_rows, variant].head(5).tolist(),
            )
        )
    input_unique = int(working[variant].nunique())
    reference_variant_count = None
    reference_id_match_rows = None
    reference_unique_id_matches = None
    reference_compatible_rows = None
    reference_unique_compatible_variants = None
    reference_coordinate_mismatch_rows = None
    reference_allele_mismatch_rows = None
    unmatched_rows = None
    overlap_fraction = None
    if resolve_variants:
        reference_by_id, reference_variant_count = read_reference_bim_matches(
            ld_ref, columns, set(working[variant]),
        )
        reference_ids = set(reference_by_id.index)
        id_matched = working[variant].isin(reference_ids)
        reference_id_match_rows = int(id_matched.sum())
        unmatched_rows = int((~id_matched).sum())
        reference_unique_id_matches = int(working.loc[id_matched, variant].nunique())

        selected = reference_by_id.loc[working.loc[id_matched, variant]]
        coordinate_matches = (
            working.loc[id_matched, "CHR_NORM"].to_numpy()
            == selected["CHR_NORM"].to_numpy()
        ) & (
            working.loc[id_matched, "BP_NORM"].to_numpy()
            == selected["BP_NORM"].to_numpy()
        )
        allele_matches = (
            working.loc[id_matched, "ALLELE_KEY"].to_numpy()
            == selected["ALLELE_KEY"].to_numpy()
        )
        compatible_matches = coordinate_matches & allele_matches
        compatible = pd.Series(False, index=working.index)
        compatible.loc[id_matched] = compatible_matches
        reference_coordinate_mismatch_rows = int((~coordinate_matches).sum())
        reference_allele_mismatch_rows = int(
            (coordinate_matches & ~allele_matches).sum()
        )
        reference_compatible_rows = int(compatible.sum())
        resolved_unique = int(working.loc[compatible, variant].nunique())
        reference_unique_compatible_variants = resolved_unique
        overlap_fraction = resolved_unique / input_unique if input_unique else 0.0
        if (
            resolved_unique == 0
            or overlap_fraction < policy.minimum_overlap_fraction
        ):
            raise MagmaError(
                "Only %d/%d unique formatter variants (%.2f%%) have compatible "
                "BIM identifiers, coordinates, and allele pairs; the configured "
                "minimum is %.2f%%. Of %d ID-matched rows, %d have coordinate "
                "mismatches and %d have allele-pair mismatches. Confirm that the "
                "formatter and LD reference use the same genome build and alleles."
                % (
                    resolved_unique,
                    input_unique,
                    overlap_fraction * 100,
                    policy.minimum_overlap_fraction * 100,
                    reference_id_match_rows,
                    reference_coordinate_mismatch_rows,
                    reference_allele_mismatch_rows,
                )
            )
        working = working.loc[compatible].copy()
    else:
        working = working.copy()

    resolved = working
    resolved = resolved.sort_values(
        [variant, "_P_NUMERIC", "_ROW_ORDER"], kind="mergesort",
    )
    before_deduplication = len(resolved)
    resolved = resolved.drop_duplicates(variant, keep="first")
    duplicate_rows = before_deduplication - len(resolved)
    resolved = resolved.sort_values("_ROW_ORDER", kind="mergesort")

    harmonised_pval = resolved[original_columns].copy()
    harmonised_pval[columns.p_value_column] = resolved["_P_NUMERIC"].to_numpy()
    harmonised_locations = pd.DataFrame(
        {
            variant: harmonised_pval[variant].to_numpy(),
            columns.chromosome_column: resolved["CHR_NORM"].to_numpy(),
            columns.position_column: resolved["BP_NORM"].to_numpy(),
        }
    )
    pval_destination = Path(pval_output)
    location_destination = Path(location_output)
    pval_destination.parent.mkdir(parents=True, exist_ok=True)
    location_destination.parent.mkdir(parents=True, exist_ok=True)
    harmonised_pval.to_csv(
        pval_destination,
        sep=columns.output_table_delimiter,
        index=False,
    )
    # MAGMA --snp-loc requires the first three columns and forbids a header.
    harmonised_locations.to_csv(
        location_destination,
        sep=columns.output_table_delimiter,
        index=False,
        header=False,
    )
    qc = {
        "input_rows": len(pval),
        "input_unique_variants": input_unique,
        "reference_intersection_enabled": resolve_variants,
        "reference_variant_count": reference_variant_count,
        "reference_id_match_rows": reference_id_match_rows,
        "reference_unique_id_matches": reference_unique_id_matches,
        "reference_compatible_rows": reference_compatible_rows,
        "reference_unique_compatible_variants": (
            reference_unique_compatible_variants
        ),
        "reference_coordinate_mismatch_rows": reference_coordinate_mismatch_rows,
        "reference_allele_mismatch_rows": reference_allele_mismatch_rows,
        "not_in_reference_rows": unmatched_rows,
        "overlap_fraction": overlap_fraction,
        "duplicates_resolved_by_lowest_p": int(duplicate_rows),
        "retained_rows": len(harmonised_pval),
    }
    logger.record("RESULT", "magma_variant_preparation", **qc)
    return {
        "pval_file": str(pval_destination),
        "snp_loc_file": str(location_destination),
        "qc": qc,
    }


def _version_tuple(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("."))


def require_supported_magma(executable: str, module_config, logger) -> str:
    output = run_checked_command(
        [executable, *module_config.version_arguments],
        "Read MAGMA version",
        logger=logger,
        error_type=MagmaError,
        stderr_to_stdout=True,
    )
    match = re.search(module_config.version_pattern, output)
    if match is None:
        raise MagmaError(
            "MAGMA version output did not match modules.magma.version_pattern: %r"
            % output.strip()
        )
    version = match.group(1)
    observed = _version_tuple(version)
    minimum = _version_tuple(module_config.minimum_magma_version)
    width = max(len(observed), len(minimum))
    if observed + (0,) * (width - len(observed)) < minimum + (0,) * (
        width - len(minimum)
    ):
        raise MagmaError(
            "MAGMA %s is too old; configured minimum is %s."
            % (version, module_config.minimum_magma_version)
        )
    logger.record("OBSERVED", "magma_version", version=version, executable=executable)
    return version


def _validate_ld_reference(ld_ref: str | Path, input_config) -> None:
    missing = [
        "%s%s" % (ld_ref, suffix)
        for suffix in input_config.required_reference_extensions
        if not Path("%s%s" % (ld_ref, suffix)).is_file()
        or Path("%s%s" % (ld_ref, suffix)).stat().st_size <= 0
    ]
    if missing:
        raise MagmaError("LD reference is incomplete or empty: %s" % ", ".join(missing))


def _count_annotation_genes(annotation_file: Path) -> int:
    with open_text(annotation_file) as handle:
        return sum(
            1
            for line in handle
            if line.strip() and not line.startswith(_RESULT_COMMENT_PREFIX)
        )


def _batch_plan(annotation_file: Path, configuration) -> tuple[int, int, int]:
    module = configuration.modules.magma
    genes = _count_annotation_genes(annotation_file)
    if genes < 1:
        raise MagmaError("MAGMA annotation contains no genes")
    if not module.batching.enabled:
        return genes, 1, 1
    memory_workers = max(
        1, int(configuration.execution.memory_gb // module.batching.memory_per_process_gb),
    )
    workers = max(1, min(configuration.execution.threads, memory_workers))
    gene_batches = max(1, genes // module.batching.minimum_genes_per_batch)
    batches = max(1, min(workers, gene_batches))
    return genes, batches, min(workers, batches)


def correct_gene_set_p_values(
    gsa_out_file: str | Path,
    output_file: str | Path,
    module_config,
    logger,
) -> pl.DataFrame:
    """Add configured global and explicitly named family corrections."""
    schema = module_config.result_schema
    source = Path(gsa_out_file).expanduser().resolve()
    frame = _read_table(
        source,
        "MAGMA gene-set output",
        schema.table_delimiter_pattern,
        comment_prefix=_RESULT_COMMENT_PREFIX,
    )
    name_column = schema.gene_set_name_column
    p_column = schema.gene_set_p_value_column
    full_name_column = schema.gene_set_full_name_column
    required = {name_column, p_column}
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise MagmaError(
            "MAGMA gene-set output is empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )
    frame[p_column] = _validate_p_values(frame[p_column], source, allow_zero=True)
    if full_name_column not in frame.columns:
        frame[full_name_column] = frame[name_column]
    else:
        frame[full_name_column] = frame[full_name_column].fillna(frame[name_column])
    frame[full_name_column] = frame[full_name_column].astype(str).str.strip()
    if frame[full_name_column].duplicated().any():
        raise MagmaError(
            "MAGMA gene-set output contains duplicate %s values" % full_name_column
        )

    values = frame[p_column].to_numpy(dtype=float)
    for method in module_config.multiple_testing.global_methods:
        column = schema.global_correction_column_pattern.format(method=method)
        frame[column] = adjust_p_values(values, method)
        logger.record(
            "PARAM", "multiple_testing_family",
            family="all_gene_sets", method=method, tests=len(values),
        )
    for name, family in module_config.multiple_testing.families.items():
        selected = frame[full_name_column].str.contains(
            family.pattern, regex=True, na=False,
        )
        family_values = frame.loc[selected, p_column].to_numpy(dtype=float)
        for method in family.methods:
            column = schema.family_correction_column_pattern.format(
                family=name, method=method,
            )
            frame[column] = np.nan
            if len(family_values):
                frame.loc[selected, column] = adjust_p_values(family_values, method)
            logger.record(
                "PARAM", "multiple_testing_family",
                family=name, pattern=family.pattern, method=method,
                tests=len(family_values),
            )
    frame = frame.sort_values(p_column, kind="mergesort")
    destination = write_delimited_report(
        frame.to_dict(orient="records"),
        output_file,
        fieldnames=frame.columns.tolist(),
        delimiter=schema.report_delimiter,
        null_value=schema.report_null_value,
    )
    logger.record("OUTPUT", "corrected_gene_sets", path=str(destination), rows=len(frame))
    return pl.DataFrame(frame.to_dict(orient="list"))


def _normalize_gene_id(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    normalized = str(value).strip()
    return normalized or None


def correct_gene_p_values(
    genes_out_file: str | Path,
    output_file: str | Path,
    module_config,
    logger,
) -> pl.DataFrame:
    """Add configured adjusted p-values to the complete MAGMA gene result."""
    schema = module_config.result_schema
    source = Path(genes_out_file).expanduser().resolve()
    frame = _read_table(
        source,
        "MAGMA gene output",
        schema.table_delimiter_pattern,
        comment_prefix=_RESULT_COMMENT_PREFIX,
    )
    gene_column = schema.gene_id_column
    p_column = schema.gene_p_value_column
    required = {gene_column, p_column}
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise MagmaError(
            "MAGMA gene output is empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )

    frame[p_column] = _validate_p_values(frame[p_column], source, allow_zero=True)
    normalized_genes = frame[gene_column].map(_normalize_gene_id)
    if normalized_genes.isna().any():
        raise MagmaError(
            "MAGMA gene output contains %d missing or empty gene identifiers"
            % int(normalized_genes.isna().sum())
        )
    if normalized_genes.duplicated().any():
        raise MagmaError("MAGMA gene output contains duplicate normalized gene IDs")

    values = frame[p_column].to_numpy(dtype=float)
    for method in module_config.multiple_testing.gene_methods:
        column = schema.global_correction_column_pattern.format(method=method)
        if column in frame.columns:
            raise MagmaError(
                "Configured gene correction column already exists in MAGMA output: %s"
                % column
            )
        frame[column] = adjust_p_values(values, method)
        logger.record(
            "PARAM",
            "multiple_testing_family",
            family="all_genes",
            method=method,
            tests=len(values),
        )

    destination = write_delimited_report(
        frame.to_dict(orient="records"),
        output_file,
        fieldnames=frame.columns.tolist(),
        delimiter=schema.report_delimiter,
        null_value=schema.report_null_value,
    )
    logger.record(
        "OUTPUT",
        "corrected_genes",
        path=str(destination),
        rows=len(frame),
        methods=list(module_config.multiple_testing.gene_methods),
    )
    return pl.DataFrame(frame.to_dict(orient="list"))


def _significance_outcome(
    frame: pl.DataFrame,
    p_column: str,
    methods: Sequence[str],
    module_config,
    singular: str,
    plural: str,
) -> tuple[str, dict, list[tuple[str, str, object]]]:
    """Summarise configured nominal and adjusted result thresholds once."""
    threshold = module_config.multiple_testing.reporting_significance_threshold
    schema = module_config.result_schema
    required = [p_column] + [
        schema.global_correction_column_pattern.format(method=method)
        for method in methods
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise MagmaError(
            "Cannot summarise significance because corrected result columns are "
            "missing: %s" % ", ".join(missing)
        )

    tested = frame.height
    nominal = int((frame[p_column] <= threshold).sum())
    adjusted = {
        method: int(
            (
                frame[
                    schema.global_correction_column_pattern.format(method=method)
                ]
                <= threshold
            ).sum()
        )
        for method in methods
    }
    labels = module_config.multiple_testing.reporting_method_labels
    subject = singular if tested == 1 else plural
    findings = [f"{nominal:,} nominally significant"]
    findings.extend(
        f"{adjusted[method]:,} significant after global {labels[method]} correction"
        for method in methods
    )
    message = (
        f"{tested:,} {subject} tested at the configured p ≤ {threshold:g} "
        f"threshold: {'; '.join(findings)}."
    )
    metrics = {
        "tested": tested,
        "significance_threshold": threshold,
        "nominal_significant": nominal,
        "adjusted_significant": adjusted,
    }
    fields = [
        ("count", f"{plural.capitalize()} tested", tested),
        ("analysis", "Reporting threshold", f"p ≤ {threshold:g}"),
        ("info", "Nominally significant", nominal),
    ]
    fields.extend(
        (
            "analysis",
            f"Significant after global {labels[method]}",
            adjusted[method],
        )
        for method in methods
    )
    return message, metrics, fields


def parse_gene_set_file(
    gene_set_path: str | Path,
    module_config,
    logger,
    *,
    input_format: str | None = None,
) -> pl.DataFrame:
    """Parse configured GMT, native MAGMA, or two-column memberships once."""
    source = require_nonempty_file(
        gene_set_path, "gene-set file", error_type=MagmaError,
    )
    configured_format = input_format or module_config.gene_sets.input_format
    schema = module_config.result_schema
    records = []
    seen: set[str] = set()
    detected_format: str | None = None
    table_delimiter = re.compile(module_config.input.table_delimiter_pattern)
    if configured_format == "membership":
        membership = module_config.gene_sets
        delimiter = re.compile(membership.membership_delimiter_pattern)
        grouped: dict[str, list[str]] = {}
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                if membership.membership_has_header and line_number == 1:
                    continue
                text = raw.strip()
                if not text or text.startswith(_RESULT_COMMENT_PREFIX):
                    continue
                fields = delimiter.split(text)
                maximum = max(
                    membership.membership_set_column,
                    membership.membership_gene_column,
                )
                if len(fields) <= maximum:
                    raise MagmaError(
                        "Malformed two-column gene-set membership at line %d"
                        % line_number
                    )
                name = fields[membership.membership_set_column].strip()
                gene = _normalize_gene_id(fields[membership.membership_gene_column])
                if not name or gene is None:
                    raise MagmaError(
                        "Empty set or gene identifier in membership line %d"
                        % line_number
                    )
                grouped.setdefault(name, []).append(gene)
        for name, values in grouped.items():
            genes = list(dict.fromkeys(values))
            joined_genes = ",".join(genes)
            records.append(
                {
                    schema.gene_set_full_name_column: name,
                    schema.report_gene_set_description_column: None,
                    schema.report_source_input_genes_column: joined_genes,
                    schema.report_input_genes_column: joined_genes,
                }
            )
        detected_format = "membership"
    else:
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.lstrip().startswith(
                    _RESULT_COMMENT_PREFIX
                ):
                    continue
                text = raw.rstrip("\r\n")
                tab_parts = text.split("\t")
                line_format = "gmt" if len(tab_parts) >= 3 else "magma"
                if detected_format is None:
                    detected_format = (
                        line_format if configured_format == "auto" else configured_format
                    )
                if configured_format == "auto" and line_format != detected_format:
                    raise MagmaError(
                        "Gene-set file mixes GMT and native MAGMA records at line %d"
                        % line_number
                    )
                if detected_format == "gmt":
                    if len(tab_parts) < 3 or not tab_parts[0].strip():
                        raise MagmaError("Malformed GMT gene-set line %d" % line_number)
                    name = tab_parts[0].strip()
                    description = tab_parts[1].strip() or None
                    gene_values = tab_parts[2:]
                else:
                    parts = table_delimiter.split(text)
                    if len(parts) < 2 or not parts[0].strip():
                        raise MagmaError(
                            "Malformed native MAGMA gene-set line %d" % line_number
                        )
                    name = parts[0].strip()
                    description = None
                    gene_values = parts[1:]
                if name in seen:
                    raise MagmaError("Duplicate gene-set name: %s" % name)
                seen.add(name)
                genes = list(
                    dict.fromkeys(
                        gene
                        for gene in (_normalize_gene_id(value) for value in gene_values)
                        if gene is not None
                    )
                )
                if not genes:
                    raise MagmaError("Gene set %s contains no valid gene IDs" % name)
                joined_genes = ",".join(genes)
                records.append(
                    {
                        schema.gene_set_full_name_column: name,
                        schema.report_gene_set_description_column: description,
                        schema.report_source_input_genes_column: joined_genes,
                        schema.report_input_genes_column: joined_genes,
                    }
                )
    if not records:
        raise MagmaError("Gene-set file contains no gene sets: %s" % source)
    logger.record(
        "DECIDE",
        "gene_set_input_format",
        configured=configured_format,
        detected=detected_format,
        gene_sets=len(records),
    )
    return pl.DataFrame(records)


def write_native_gene_set_file(
    gene_sets: pl.DataFrame,
    output_file: str | Path,
    module_config,
) -> Path:
    """Write parsed memberships in MAGMA's native set-annotation format."""
    schema = module_config.result_schema
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with destination.open("w", encoding="utf-8", newline="\n") as handle:
            for record in gene_sets.iter_rows(named=True):
                handle.write(
                    "%s %s\n"
                    % (
                        record[schema.gene_set_full_name_column],
                        record[schema.report_input_genes_column].replace(",", " "),
                    )
                )
    except OSError as exc:
        raise MagmaError(
            "Cannot write native MAGMA gene-set file %s: %s"
            % (destination, exc)
        ) from exc
    return require_nonempty_file(
        destination, "Prepared MAGMA gene-set file", error_type=MagmaError,
    )


def annotate_gene_set_results(
    genes_out_file: str | Path,
    gene_sets: pl.DataFrame,
    corrected: pl.DataFrame,
    output_file: str | Path,
    dataset_id: str,
    module_config,
    logger,
) -> str:
    """Annotate corrected gene-set rows with tested genes and gene p-values."""
    schema = module_config.result_schema
    source = Path(genes_out_file).expanduser().resolve()
    genes = _read_table(
        source,
        "MAGMA gene output",
        schema.table_delimiter_pattern,
        comment_prefix=_RESULT_COMMENT_PREFIX,
    )
    gene_column = schema.gene_id_column
    gene_p_column = schema.gene_p_value_column
    missing = sorted({gene_column, gene_p_column} - set(genes.columns))
    if missing or genes.empty:
        raise MagmaError(
            "MAGMA gene output is empty or missing columns: %s"
            % ", ".join(missing or [gene_column, gene_p_column])
        )
    genes[gene_p_column] = _validate_p_values(
        genes[gene_p_column], source, allow_zero=True,
    )
    genes[gene_column] = genes[gene_column].map(_normalize_gene_id)
    genes = genes.dropna(subset=[gene_column])
    if genes[gene_column].duplicated().any():
        raise MagmaError("MAGMA gene output contains duplicate normalized gene IDs")
    gene_p = dict(zip(genes[gene_column], genes[gene_p_column]))
    details = []
    for gene_string in gene_sets[schema.report_input_genes_column].to_list():
        values = list(dict.fromkeys(gene_string.split(",")))
        common = sorted(
            ((gene, gene_p[gene]) for gene in values if gene in gene_p),
            key=lambda item: item[0],
        )
        details.append(
            {
                schema.report_common_genes_column: (
                    ",".join(gene for gene, _ in common) or None
                ),
                schema.report_common_gene_p_values_column: (
                    ",".join(str(value) for _, value in common) or None
                ),
                schema.report_total_genes_column: len(values),
                schema.report_common_gene_count_column: len(common),
            }
        )
    annotated = (
        corrected.join(
            pl.concat([gene_sets, pl.DataFrame(details)], how="horizontal_extend"),
            on=schema.gene_set_full_name_column,
            how="left",
            coalesce=True,
        )
        .sort(schema.gene_set_p_value_column)
        .with_columns(pl.lit(dataset_id).alias(schema.report_dataset_column))
    )
    destination = write_delimited_report(
        annotated.to_dicts(),
        output_file,
        fieldnames=annotated.columns,
        delimiter=schema.report_delimiter,
        null_value=schema.report_null_value,
    )
    logger.record("OUTPUT", "annotated_gene_sets", path=str(destination), rows=annotated.height)
    return str(destination)


def _run_command(
    command: Sequence[str],
    purpose: str,
    expected_outputs: Sequence[Path],
    configuration,
    logger,
) -> None:
    run_checked_command(
        command,
        purpose,
        logger=logger,
        error_type=MagmaError,
        timeout_seconds=configuration.execution.timeout_seconds,
        expected_outputs=expected_outputs,
    )


def _run_gene_associations(
    executable: str,
    ld_ref: str,
    annotation_file: Path,
    paths: dict[str, Path],
    definition,
    configuration,
    logger,
) -> tuple[int, int, int]:
    """Run one configured MAGMA gene analysis with bounded batching."""
    module = configuration.modules.magma
    inputs = module.input
    genes, batches, workers = _batch_plan(annotation_file, configuration)
    common = [
        executable,
        "--bfile", ld_ref,
        "--gene-annot", str(annotation_file),
        "--pval", str(paths["harmonised_p_values"]),
        "use=%s,%s" % (inputs.variant_id_column, inputs.p_value_column),
        "ncol=%s" % inputs.sample_size_column,
        "duplicate=error",
        "--gene-model", module.gene_model,
        "--seed", str(configuration.execution.random_seed),
    ]
    if definition.gene_settings:
        common.extend(["--gene-settings", *definition.gene_settings])
    if batches == 1:
        _run_command(
            [*common, "--out", str(paths["gene_prefix"])],
            "MAGMA %s gene association" % definition.display_name,
            [paths["genes_raw"], paths["genes_out"]],
            configuration,
            logger,
        )
    else:
        def run_batch(number: int) -> None:
            expected = [
                _tool_output(
                    paths["gene_batch_prefix"],
                    ".batch%d_%d%s" % (number, batches, _GENE_RAW_SUFFIX),
                ),
                _tool_output(
                    paths["gene_batch_prefix"],
                    ".batch%d_%d%s" % (number, batches, _GENE_RESULT_SUFFIX),
                ),
            ]
            _run_command(
                [
                    *common,
                    "--out", str(paths["gene_batch_prefix"]),
                    "--batch", str(number), str(batches),
                ],
                "MAGMA %s batch %d/%d"
                % (definition.display_name, number, batches),
                expected,
                configuration,
                logger,
            )

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(run_batch, number)
                for number in range(1, batches + 1)
            ]
            try:
                for future in as_completed(futures):
                    future.result()
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        _run_command(
            [
                executable,
                "--merge", str(paths["gene_batch_prefix"]),
                "--out", str(paths["gene_prefix"]),
            ],
            "MAGMA %s batch merge" % definition.display_name,
            [paths["genes_raw"], paths["genes_out"]],
            configuration,
            logger,
        )
    return genes, batches, workers


def _mapping_gene_set(definition, inputs) -> Path | None:
    if definition.method == "chrom_magma":
        return None
    configured = definition.gene_set_file or inputs.gene_set_file
    return (
        require_nonempty_file(configured, "gene-set file", error_type=MagmaError)
        if configured is not None
        else None
    )


def _write_mapping_comparison(
    analyses: dict[str, dict],
    output_file: Path,
    module,
    logger,
) -> str:
    """Write one provenance-rich long table without pooling hypothesis families."""
    schema = module.result_schema
    frames = []
    for name in module.mapping.selected:
        definition = module.mapping.definitions[name]
        frame = pd.DataFrame(analyses[name]["corrected_gene_results"].to_dicts())
        metadata = {
            schema.report_mapping_name_column: name,
            schema.report_mapping_method_column: definition.method,
            schema.report_mapping_context_column: definition.context,
            schema.report_gene_id_type_column: definition.gene_id_type,
            schema.report_source_name_column: definition.source_name,
            schema.report_source_version_column: definition.source_version,
            schema.report_statistic_type_column: definition.result_statistic_type,
            schema.report_statistic_interpretation_column: (
                definition.result_statistic_interpretation
            ),
        }
        for column, value in reversed(list(metadata.items())):
            frame.insert(0, column, value)
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True, sort=False)
    destination = write_delimited_report(
        combined.to_dict(orient="records"),
        output_file,
        fieldnames=combined.columns.tolist(),
        delimiter=schema.report_delimiter,
        null_value=schema.report_null_value,
    )
    logger.record(
        "OUTPUT",
        "magma_mapping_comparison",
        path=str(destination),
        rows=len(combined),
        mapping_analyses=list(module.mapping.selected),
        statistical_pooling=False,
    )
    return str(destination)


def preflight_magma_analysis(
    dataset_id: str,
    configuration,
    logger,
) -> MagmaPreflight:
    """Validate run-level MAGMA resources without creating analysis outputs."""
    validate_filename_component(dataset_id, "dataset_id", error_type=MagmaError)
    module = configuration.modules.magma
    inputs = module.input
    if configuration.execution.random_seed < 1:
        raise MagmaError("MAGMA requires a positive execution.random_seed")
    required = {
        "SNP-location file": inputs.snp_location_file,
        "p-value file": inputs.p_value_file,
        "LD-reference prefix": inputs.ld_reference_prefix,
    }
    missing = [label for label, value in required.items() if value is None]
    if missing:
        raise MagmaError(
            "Required MAGMA inputs are missing: %s. Provide CLI options or set "
            "modules.magma.input in --run-config." % ", ".join(missing)
        )
    snp_location_file = require_nonempty_file(
        inputs.snp_location_file, "SNP-location file", error_type=MagmaError,
    )
    p_value_file = require_nonempty_file(
        inputs.p_value_file, "p-value file", error_type=MagmaError,
    )
    ld_reference_prefix = str(
        Path(inputs.ld_reference_prefix).expanduser().resolve()
    )
    _validate_ld_reference(ld_reference_prefix, inputs)
    executable = resolve_executable(
        configuration.resources.executables.magma,
        "MAGMA executable",
        error_type=MagmaError,
    )
    version = require_supported_magma(executable, module, logger)
    mapping_gene_sets = {
        name: _mapping_gene_set(module.mapping.definitions[name], inputs)
        for name in module.mapping.selected
    }
    return MagmaPreflight(
        snp_location_file=snp_location_file,
        p_value_file=p_value_file,
        ld_reference_prefix=ld_reference_prefix,
        executable=executable,
        version=version,
        mapping_gene_sets=mapping_gene_sets,
    )


def run_magma_analysis(
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger,
    *,
    preflight: MagmaPreflight | None = None,
) -> dict:
    """Run positional and configured functional MAGMA mappings independently."""
    prepared = preflight or preflight_magma_analysis(
        dataset_id, configuration, logger,
    )
    module = configuration.modules.magma
    inputs = module.input
    snp_loc = prepared.snp_location_file
    pval = prepared.p_value_file
    ld_ref = prepared.ld_reference_prefix
    executable = prepared.executable
    version = prepared.version
    mapping_gene_sets = prepared.mapping_gene_sets
    total_steps = 1 + sum(
        3 + (4 if mapping_gene_sets[name] is not None else 0)
        for name in module.mapping.selected
    )
    shared_paths = resolve_magma_output_paths(
        output_directory, dataset_id, module, module.mapping.primary,
    )
    for path in shared_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    step_number = 1
    variant_stage_title = (
        "Retain variants present in LD reference"
        if module.snp_harmonisation.resolve_variants_to_reference
        else "Prepare MAGMA input tables"
    )
    with logger.step(
        step_number, total_steps, variant_stage_title, "prepare_magma_variant_inputs",
    ) as step:
        variant_preparation = prepare_magma_variant_inputs(
            pval,
            snp_loc,
            ld_ref,
            shared_paths["harmonised_p_values"],
            shared_paths["harmonised_snp_locations"],
            module,
            logger,
        )
        qc = variant_preparation["qc"]
        removed = qc["input_rows"] - qc["retained_rows"]
        step.set_rows(qc["retained_rows"], removed=removed)
        variant_fields = [
            ("count", "Formatter rows", qc["input_rows"]),
            ("success", "Prepared variants", qc["retained_rows"]),
            (
                "info", "BIM intersection",
                "applied" if qc["reference_intersection_enabled"] else "not requested",
            ),
        ]
        if qc["reference_intersection_enabled"]:
            variant_fields.extend(
                [
                    (
                        "count", "Unique BIM ID matches",
                        qc["reference_unique_id_matches"],
                    ),
                    (
                        "success", "Compatible unique variants",
                        qc["reference_unique_compatible_variants"],
                    ),
                    ("warning", "Rows absent from BIM", qc["not_in_reference_rows"]),
                    (
                        "warning", "Rows with coordinate mismatch",
                        qc["reference_coordinate_mismatch_rows"],
                    ),
                    (
                        "warning", "Rows with allele-pair mismatch",
                        qc["reference_allele_mismatch_rows"],
                    ),
                ]
            )
        variant_fields.append(
            (
                "info", "Duplicate rows consolidated",
                "%s (lowest-p-value rule)"
                % f"{qc['duplicates_resolved_by_lowest_p']:,}",
            )
        )
        variant_outcome = (
            "%s formatter rows; %s ID-, coordinate-, and allele-compatible "
            "variants retained; %s rows were not retained."
            % (
                f"{qc['input_rows']:,}",
                f"{qc['retained_rows']:,}",
                f"{removed:,}",
            )
            if qc["reference_intersection_enabled"]
            else "%s formatter rows produced %s unique MAGMA variants; %s "
            "rows were not retained."
            % (
                f"{qc['input_rows']:,}",
                f"{qc['retained_rows']:,}",
                f"{removed:,}",
            )
        )
        step.outcome(
            variant_outcome,
            fields=variant_fields,
            **qc,
        )
    step_number += 1

    analyses: dict[str, dict] = {}
    for name in module.mapping.selected:
        definition = module.mapping.definitions[name]
        gene_set_file = mapping_gene_sets[name]
        parsed_gene_sets = None
        prepared_gene_sets = None
        gene_id_validation = None
        tested_gene_coverage = None
        paths = resolve_magma_output_paths(
            output_directory, dataset_id, module, name,
        )
        for path in paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        paths["harmonised_p_values"] = shared_paths["harmonised_p_values"]
        paths["harmonised_snp_locations"] = shared_paths["harmonised_snp_locations"]

        with logger.step(
            step_number,
            total_steps,
            "%s · prepare SNP-to-unit annotation" % definition.display_name,
            "magma_mapping_annotation",
        ) as step:
            if definition.method == "positional":
                configured_location = (
                    definition.gene_location_file or inputs.gene_location_file
                )
                if configured_location is None:
                    raise MagmaError(
                        "Mapping %s requires a gene-location file" % name
                    )
                location_reference = require_nonempty_file(
                    configured_location, "%s gene-location file" % name,
                    error_type=MagmaError,
                )
                upstream = (
                    definition.annotation_window_upstream_kb
                    if definition.annotation_window_upstream_kb is not None
                    else module.gene_window_upstream_kb
                )
                downstream = (
                    definition.annotation_window_downstream_kb
                    if definition.annotation_window_downstream_kb is not None
                    else module.gene_window_downstream_kb
                )
                annotation_file = paths["gene_annotation"]
                _run_command(
                    [
                        executable,
                        "--annotate", "window=%s,%s" % (upstream, downstream),
                        "--snp-loc", str(paths["harmonised_snp_locations"]),
                        "--gene-loc", str(location_reference),
                        "--out", str(paths["annotation_prefix"]),
                    ],
                    "MAGMA %s annotation" % definition.display_name,
                    [annotation_file],
                    configuration,
                    logger,
                )
                annotation_validation = {
                    "source": "generated_by_magma",
                    "upstream_kb": upstream,
                    "downstream_kb": downstream,
                }
            elif definition.method == "n_magma":
                location_reference = require_nonempty_file(
                    definition.gene_location_file,
                    "%s positional gene-location file" % name,
                    error_type=MagmaError,
                )
                _run_command(
                    [
                        executable,
                        "--annotate",
                        "window=%s,%s"
                        % (
                            definition.annotation_window_upstream_kb,
                            definition.annotation_window_downstream_kb,
                        ),
                        "--snp-loc", str(paths["harmonised_snp_locations"]),
                        "--gene-loc", str(location_reference),
                        "--out", str(paths["component_annotation_prefix"]),
                    ],
                    "MAGMA %s positional annotation component"
                    % definition.display_name,
                    [paths["component_gene_annotation"]],
                    configuration,
                    logger,
                )
                annotation_file = paths["gene_annotation"]
                merge_result = merge_gene_annotations(
                    [
                        paths["component_gene_annotation"],
                        *definition.gene_annotation_files,
                    ],
                    location_reference,
                    annotation_file,
                    module,
                    logger,
                )
                overlap_validation = validate_gene_annotation(
                    annotation_file, ld_ref, module, logger,
                )
                annotation_validation = {
                    "source": "merged_n_magma_components",
                    "positional_window_upstream_kb": (
                        definition.annotation_window_upstream_kb
                    ),
                    "positional_window_downstream_kb": (
                        definition.annotation_window_downstream_kb
                    ),
                    "merge": merge_result,
                    "bim_overlap": overlap_validation,
                }
            elif definition.method == "chrom_magma":
                location_reference = require_nonempty_file(
                    definition.regulatory_element_location_file,
                    "%s regulatory-element location file" % name,
                    error_type=MagmaError,
                )
                annotation_file = paths["gene_annotation"]
                _run_command(
                    [
                        executable,
                        "--annotate",
                        "window=%s,%s"
                        % (
                            definition.annotation_window_upstream_kb,
                            definition.annotation_window_downstream_kb,
                        ),
                        "--snp-loc", str(paths["harmonised_snp_locations"]),
                        "--gene-loc", str(location_reference),
                        "--out", str(paths["annotation_prefix"]),
                    ],
                    "MAGMA %s regulatory-element annotation"
                    % definition.display_name,
                    [annotation_file],
                    configuration,
                    logger,
                )
                annotation_validation = {"source": "generated_by_magma"}
            else:
                annotation_file = require_nonempty_file(
                    definition.gene_annotation_file,
                    "%s gene annotation" % name,
                    error_type=MagmaError,
                )
                annotation_validation = validate_gene_annotation(
                    annotation_file, ld_ref, module, logger,
                )
            annotated_units = _count_annotation_genes(annotation_file)
            annotation_fields = [
                ("analysis", "Mapping method", definition.method),
                ("genetic", "Annotated units", annotated_units),
                ("info", "Biological context", definition.context),
            ]
            if definition.method == "n_magma":
                annotation_fields.extend(
                    [
                        (
                            "genetic",
                            "Canonical gene locations",
                            merge_result["canonical_gene_locations"],
                        ),
                        (
                            "info",
                            "Component coordinate differences",
                            merge_result["coordinate_disagreement_records"],
                        ),
                        (
                            "info",
                            "Invalid component coordinates replaced",
                            merge_result["invalid_component_coordinate_records"],
                        ),
                        (
                            "warning",
                            "Genes outside location reference",
                            merge_result["genes_absent_from_location_reference"],
                        ),
                    ]
                )
            step.set_rows(annotated_units)
            step.outcome(
                "%s annotation contains %s tested units."
                % (definition.display_name, f"{annotated_units:,}"),
                fields=annotation_fields,
                mapping=name,
                method=definition.method,
                annotated_units=annotated_units,
                validation=annotation_validation,
            )
            step.output("gene_annotation", path=str(annotation_file))
        step_number += 1

        if gene_set_file is not None:
            with logger.step(
                step_number,
                total_steps,
                "%s · validate competitive gene sets" % definition.display_name,
                "parse_gene_set_file",
            ) as step:
                parsed_gene_sets = parse_gene_set_file(
                    gene_set_file,
                    module,
                    logger,
                    input_format=definition.gene_set_format,
                )
                if definition.method in {"positional", "n_magma"}:
                    parsed_gene_sets, reference_gene_ids, identifier_resolution = (
                        resolve_gene_set_identifiers(
                            location_reference,
                            parsed_gene_sets,
                            (
                                definition.minimum_gene_id_overlap_fraction
                                or module.gene_sets.minimum_gene_id_overlap_fraction
                            ),
                            module,
                        )
                    )
                else:
                    reference_gene_ids = _gene_ids_from_table(
                        annotation_file,
                        "%s gene annotation" % definition.display_name,
                        module.input.table_delimiter_pattern,
                        0,
                        has_header=False,
                        comment_prefix=module.annotation_validation.comment_prefix,
                    )
                    identifier_resolution = {
                        "identifier_source": "annotation_gene_id",
                        "translated_unique_ids": 0,
                        "one_to_many_identifiers": 0,
                    }
                gene_id_validation = validate_gene_set_identifier_overlap(
                    reference_gene_ids,
                    parsed_gene_sets,
                    (
                        definition.minimum_gene_id_overlap_fraction
                        or module.gene_sets.minimum_gene_id_overlap_fraction
                    ),
                    module.result_schema.report_input_genes_column,
                )
                gene_id_validation["identifier_resolution"] = identifier_resolution
                prepared_gene_sets = write_native_gene_set_file(
                    parsed_gene_sets, paths["prepared_gene_sets"], module,
                )
                logger.record(
                    "RESULT",
                    "gene_id_compatibility",
                    mapping=name,
                    reference_universe="complete_mapping_reference",
                    **gene_id_validation,
                )
                step.set_rows(parsed_gene_sets.height)
                step.outcome(
                    "%s validated %s competitive gene sets against the complete "
                    "mapping reference."
                    % (definition.display_name, f"{parsed_gene_sets.height:,}"),
                    fields=[
                        (
                            "genetic", "Reference genes",
                            gene_id_validation["reference_unique_ids"],
                        ),
                        (
                            "count", "Unique gene-set genes",
                            gene_id_validation["geneset_unique_ids"],
                        ),
                        (
                            "success", "Overlapping identifiers",
                            "%s/%s (%.2f%%)"
                            % (
                                f"{gene_id_validation['overlapping_unique_ids']:,}",
                                f"{gene_id_validation['comparison_unique_ids']:,}",
                                gene_id_validation[
                                    "overlap_fraction_of_smaller_universe"
                                ] * 100,
                            ),
                        ),
                        (
                            "success", "Gene sets represented in reference",
                            gene_id_validation[
                                "gene_sets_with_at_least_one_reference_gene"
                            ],
                        ),
                        (
                            "info", "Gene-set identifier source",
                            identifier_resolution["identifier_source"],
                        ),
                        (
                            "count", "Alternate IDs translated",
                            identifier_resolution["translated_unique_ids"],
                        ),
                        (
                            "info", "One-to-many IDs expanded",
                            identifier_resolution["one_to_many_identifiers"],
                        ),
                    ],
                    mapping=name,
                    reference_universe="complete_mapping_reference",
                    **gene_id_validation,
                )
            step_number += 1

        with logger.step(
            step_number,
            total_steps,
            "%s · calculate associations" % definition.display_name,
            "magma_gene_analysis",
        ) as step:
            genes, batches, workers = _run_gene_associations(
                executable,
                ld_ref,
                annotation_file,
                paths,
                definition,
                configuration,
                logger,
            )
            step.set_rows(genes)
            step.outcome(
                "%s submitted %s annotated units in %s batch(es)."
                % (definition.display_name, f"{genes:,}", f"{batches:,}"),
                fields=[
                    ("genetic", "Annotated units submitted", genes),
                    ("analysis", "Association batches", batches),
                    ("analysis", "Parallel workers", workers),
                ],
                mapping=name,
                annotated_units=genes,
                batches=batches,
                workers=workers,
            )
        step_number += 1

        with logger.step(
            step_number,
            total_steps,
            (
                "%s · map regulatory-element results to genes"
                if definition.method == "chrom_magma"
                else "%s · correct p-values"
            )
            % definition.display_name,
            "correct_gene_p_values",
        ) as step:
            chrom_gene_report = None
            if definition.method == "chrom_magma":
                element_results = _read_table(
                    paths["genes_out"],
                    "chromMAGMA regulatory-element output",
                    module.result_schema.table_delimiter_pattern,
                    comment_prefix=_RESULT_COMMENT_PREFIX,
                )
                mapped = map_regulatory_elements_to_genes(
                    element_results,
                    definition.regulatory_element_location_file,
                    definition.element_to_gene_file,
                    paths["chrom_magma_genes"],
                    module,
                    logger,
                )
                chrom_gene_report = str(paths["chrom_magma_genes"])
                result_frame = mapped
                gene_result_path = paths["chrom_magma_genes"]
                corrected_gene_path = None
                significance = {
                    "tested": mapped.height,
                    "significance_not_reported": True,
                    "reason": definition.result_statistic_interpretation,
                }
                step.set_rows(mapped.height)
                step.outcome(
                    "%s mapped %s genes. The minimum linked element p-value is "
                    "reported for ranking; gene-level significance correction is "
                    "not applied."
                    % (definition.display_name, f"{mapped.height:,}"),
                    fields=[
                        ("genetic", "Genes with mapped elements", mapped.height),
                        (
                            "info", "Statistical interpretation",
                            definition.result_statistic_interpretation,
                        ),
                    ],
                    mapping=name,
                    **significance,
                )
                step.output("mapped_gene_ranking", path=str(gene_result_path))
            else:
                corrected_genes = correct_gene_p_values(
                    paths["genes_out"], paths["corrected_genes"], module, logger,
                )
                result_frame = corrected_genes
                gene_result_path = paths["corrected_genes"]
                corrected_gene_path = str(paths["corrected_genes"])
                outcome, significance, fields = _significance_outcome(
                    corrected_genes,
                    module.result_schema.gene_p_value_column,
                    module.multiple_testing.gene_methods,
                    module,
                    "gene",
                    "genes",
                )
                step.set_rows(corrected_genes.height)
                step.outcome(outcome, fields=fields, mapping=name, **significance)
                step.output(
                    "corrected_gene_associations", path=str(paths["corrected_genes"]),
                )
        step_number += 1

        if parsed_gene_sets is not None:
            tested_gene_ids = {
                identifier
                for identifier in (
                    _normalize_gene_id(value)
                    for value in result_frame[
                        module.result_schema.gene_id_column
                    ].to_list()
                )
                if identifier is not None
            }
            tested_gene_coverage = validate_gene_set_identifier_overlap(
                tested_gene_ids,
                parsed_gene_sets,
                None,
                module.result_schema.report_input_genes_column,
            )
            if not tested_gene_coverage[
                "gene_sets_with_at_least_one_reference_gene"
            ]:
                raise MagmaError(
                    "%s produced no tested genes represented in the configured "
                    "gene sets; competitive analysis cannot be interpreted."
                    % definition.display_name
                )
            logger.record(
                "RESULT",
                "gene_set_tested_gene_coverage",
                mapping=name,
                reference_universe="study_tested_genes",
                **tested_gene_coverage,
            )

        analysis_result = {
            "mapping_name": name,
            "mapping_method": definition.method,
            "display_name": definition.display_name,
            "biological_context": definition.context,
            "gene_id_type": definition.gene_id_type,
            "annotation_source": definition.source_name,
            "annotation_version": definition.source_version,
            "result_statistic_type": definition.result_statistic_type,
            "result_statistic_interpretation": (
                definition.result_statistic_interpretation
            ),
            "gene_annotation": str(annotation_file),
            "magma_genes_prefix": str(paths["gene_prefix"]),
            "magma_genes_raw": str(paths["genes_raw"]),
            "magma_genes_out": str(paths["genes_out"]),
            "magma_genes_corrected": corrected_gene_path,
            "magma_gene_results": str(gene_result_path),
            "chrom_magma_gene_report": chrom_gene_report,
            "corrected_gene_results": result_frame,
            "annotation_validation": annotation_validation,
            "batching": {
                "annotated_units": genes,
                "batches": batches,
                "workers": workers,
            },
            "gene_id_validation": gene_id_validation,
            "tested_gene_coverage": tested_gene_coverage,
        }

        if gene_set_file is not None:
            with logger.step(
                step_number,
                total_steps,
                "%s · test competitive gene sets" % definition.display_name,
                "magma_gene_set_analysis",
            ) as step:
                _run_command(
                    [
                        executable,
                        "--gene-results", str(paths["genes_raw"]),
                        "--set-annot", str(prepared_gene_sets),
                        "--out", str(paths["gene_set_prefix"]),
                        "--seed", str(configuration.execution.random_seed),
                    ],
                    "MAGMA %s competitive gene-set analysis"
                    % definition.display_name,
                    [paths["gene_sets_raw"]],
                    configuration,
                    logger,
                )
                step.set_rows(parsed_gene_sets.height)
                step.outcome(
                    "%s submitted %s validated gene sets using %s tested genes."
                    % (
                        definition.display_name,
                        f"{parsed_gene_sets.height:,}",
                        f"{tested_gene_coverage['reference_unique_ids']:,}",
                    ),
                    fields=[
                        (
                            "count", "Validated gene sets submitted",
                            parsed_gene_sets.height,
                        ),
                        (
                            "genetic", "Study-tested genes",
                            tested_gene_coverage["reference_unique_ids"],
                        ),
                        (
                            "genetic", "Gene-set genes tested",
                            "%s/%s"
                            % (
                                f"{tested_gene_coverage['overlapping_unique_ids']:,}",
                                f"{tested_gene_coverage['geneset_unique_ids']:,}",
                            ),
                        ),
                        (
                            "success", "Gene sets containing tested genes",
                            tested_gene_coverage[
                                "gene_sets_with_at_least_one_reference_gene"
                            ],
                        ),
                    ],
                    mapping=name,
                    tested_gene_coverage=tested_gene_coverage,
                )
            step_number += 1

            with logger.step(
                step_number,
                total_steps,
                "%s · correct gene-set p-values" % definition.display_name,
                "correct_gene_set_p_values",
            ) as step:
                corrected_sets = correct_gene_set_p_values(
                    paths["gene_sets_raw"], paths["corrected_gene_sets"], module, logger,
                )
                outcome, significance, fields = _significance_outcome(
                    corrected_sets,
                    module.result_schema.gene_set_p_value_column,
                    module.multiple_testing.global_methods,
                    module,
                    "gene set",
                    "gene sets",
                )
                step.set_rows(corrected_sets.height)
                step.outcome(outcome, fields=fields, mapping=name, **significance)
            step_number += 1

            with logger.step(
                step_number,
                total_steps,
                "%s · annotate gene-set results" % definition.display_name,
                "annotate_gene_set_results",
            ) as step:
                pathway = annotate_gene_set_results(
                    paths["genes_out"],
                    parsed_gene_sets,
                    corrected_sets,
                    paths["annotated_gene_sets"],
                    dataset_id,
                    module,
                    logger,
                )
                step.set_rows(corrected_sets.height)
                step.outcome(
                    "%s wrote %s annotated competitive gene-set results."
                    % (definition.display_name, f"{corrected_sets.height:,}"),
                    fields=[
                        ("success", "Tested gene-set results", corrected_sets.height),
                    ],
                    mapping=name,
                )
            step_number += 1
            analysis_result.update(
                {
                    "prepared_gene_sets": str(prepared_gene_sets),
                    "magma_gene_sets_raw": str(paths["gene_sets_raw"]),
                    "magma_gene_sets_corrected": str(paths["corrected_gene_sets"]),
                    "magma_pathway": pathway,
                }
            )
        analyses[name] = analysis_result

    comparison = _write_mapping_comparison(
        analyses, shared_paths["mapping_comparison"], module, logger,
    )
    primary = analyses[module.mapping.primary]
    result = {
        "magma_version": version,
        "magma_executable": executable,
        "variant_preparation": variant_preparation,
        "mapping_analyses": {
            name: {
                key: value
                for key, value in analysis.items()
                if key != "corrected_gene_results"
            }
            for name, analysis in analyses.items()
        },
        "primary_mapping": module.mapping.primary,
        "magma_mapping_comparison": comparison,
        "magma_genes_raw": primary["magma_genes_raw"],
        "magma_genes_out": primary["magma_genes_out"],
        "magma_gene_results": primary["magma_gene_results"],
        "magma_genes_prefix": primary["magma_genes_prefix"],
        "magma_gene_annotation": primary["gene_annotation"],
        "gene_id_validation": primary["gene_id_validation"],
        "tested_gene_coverage": primary["tested_gene_coverage"],
        "batching": primary["batching"],
    }
    if primary["magma_genes_corrected"] is not None:
        result["magma_genes_corrected"] = primary["magma_genes_corrected"]
    if primary["chrom_magma_gene_report"] is not None:
        result["chrom_magma_gene_report"] = primary["chrom_magma_gene_report"]
    for key in (
        "magma_gene_sets_raw", "magma_gene_sets_corrected", "magma_pathway",
    ):
        if key in primary:
            result[key] = primary[key]
    return result


__all__ = [
    "MagmaPreflight",
    "annotate_gene_set_results",
    "correct_gene_p_values",
    "correct_gene_set_p_values",
    "prepare_magma_variant_inputs",
    "parse_gene_set_file",
    "preflight_magma_analysis",
    "require_supported_magma",
    "resolve_gene_set_identifiers",
    "resolve_magma_output_paths",
    "run_magma_analysis",
]
