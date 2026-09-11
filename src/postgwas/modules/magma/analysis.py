"""Validated MAGMA gene and competitive gene-set analysis.

MAGMA protocol invariants retained in code are limited to its documented output
suffixes, metadata marker, and command modifiers. Runtime policies, paths,
thresholds, columns, and resource limits come from the canonical MAGMA YAML
configuration.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from postgwas.core.gene_sets import GeneSetFormat, validate_gene_set_source
from postgwas.core.io.delimiters import open_text, pandas_csv_engine
from postgwas.core.io.reports import write_delimited_report
from postgwas.core.genomic_scope import resolve_genomic_analysis_scope
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.processes import run_checked_command
from postgwas.core.plink import validate_plink_files
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magma.annotations import (
    map_regulatory_elements_to_genes,
    merge_gene_annotations,
    prepare_scoped_gene_annotation,
    read_gene_locations,
    validate_gene_annotation,
    write_pathway_compatible_gene_locations,
)
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.reference import read_reference_bim_matches
from postgwas.modules.magma.reporting import (
    MAGMA_PIPELINE_STAGE_KEYS,
    magma_variant_input_outcome_fields,
)


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
    gene_set_plans: dict[str, dict]
    analysis_scope: dict


@dataclass(frozen=True)
class MagmaReferencePreflight:
    """Validated external resources independent of formatter-created tables."""

    ld_reference_prefix: str
    executable: str
    version: str
    gene_set_plans: dict[str, dict]
    analysis_scope: dict


def _tool_output(prefix: str | Path, suffix: str) -> Path:
    return Path("%s%s" % (prefix, suffix))


def _annotation_window_argument(upstream_kb: int, downstream_kb: int) -> str:
    """Return MAGMA's documented upstream-then-downstream window modifier."""
    return "window=%s,%s" % (upstream_kb, downstream_kb)


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
        "excluded_variants": configured_output_path(
            output,
            layout.excluded_variants,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "mapping_comparison": configured_output_path(
            output,
            layout.mapping_comparison,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "pipeline_summary_csv": configured_output_path(
            output,
            layout.pipeline_summary_csv,
            error_type=MagmaError,
            dataset_id=dataset_id,
        ),
        "pipeline_summary_html": configured_output_path(
            output,
            layout.pipeline_summary_html,
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
        "scoped_annotation": configured_output_path(
            output,
            layout.scoped_annotation,
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
        "pathway_compatible_gene_locations": configured_output_path(
            output,
            layout.pathway_compatible_gene_locations,
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
        "excluded_genes": configured_output_path(
            output,
            layout.excluded_genes,
            error_type=MagmaError,
            dataset_id=dataset_id,
            analysis=selected,
        ),
        "exclusion_summary": configured_output_path(
            output,
            layout.exclusion_summary,
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


def _analysis_scope(configuration) -> dict:
    """Resolve build-specific MAGMA exclusions once from validated configuration."""
    module = configuration.modules.magma
    def normalize(value: object) -> str:
        return str(
            _normalize_chromosome(
                pd.Series([str(value)]), module.input,
            ).iloc[0]
        )
    return resolve_genomic_analysis_scope(
        genome_build=module.genome_build,
        genomes=configuration.resources.genomes,
        mhc=module.mhc,
        chromosomes=module.chromosomes,
        normalize_chromosome=normalize,
        analysis_name="MAGMA",
        error_type=MagmaError,
    )


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
    logger,
    *,
    comment_prefix: str | None = None,
    column_types: dict[str, str] | None = None,
) -> pd.DataFrame:
    file_path = require_nonempty_file(path, label, error_type=MagmaError)
    engine = pandas_csv_engine(delimiter_pattern)
    logger.record(
        "PARAM",
        "table_parser",
        label=label,
        path=str(file_path),
        delimiter_pattern=delimiter_pattern,
        pandas_engine=engine,
    )
    try:
        return pd.read_csv(
            file_path,
            sep=delimiter_pattern,
            comment=comment_prefix,
            engine=engine,
            dtype=column_types,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise MagmaError("Cannot parse %s %s: %s" % (label, file_path, exc)) from exc


def _require_available_result_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    """Reject configured output-column collisions before changing a result."""
    counts = Counter(columns)
    duplicates = sorted(column for column, count in counts.items() if count > 1)
    if duplicates:
        raise MagmaError(
            "Configured %s output columns are not unique: %s"
            % (label, ", ".join(duplicates))
        )
    existing = sorted(set(columns) & set(frame.columns))
    if existing:
        raise MagmaError(
            "Configured %s output columns already exist in MAGMA output: %s"
            % (label, ", ".join(existing))
        )


def resolve_gene_set_identifiers(
    location_file: str | Path,
    gene_sets: pl.DataFrame,
    minimum_overlap: float,
    module_config,
    *,
    enforce: bool = True,
    allow_alternate_reference_ids: bool = False,
) -> tuple[pl.DataFrame, set[str], dict]:
    """Select compatible gene-reference IDs without changing pathway members."""
    locations = read_gene_locations(location_file, module_config)
    location_roles = module_config.input.gene_location_columns
    alternate_column = (
        location_roles.index("alternate_gene_id") + 1
        if "alternate_gene_id" in location_roles else None
    )
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
    alternate_candidate_matches = genes & set(alternate_ids)
    alternate_matches = (genes - primary_ids) & set(alternate_ids)
    ambiguous_alternates = {
        identifier for identifier, targets in alternate_ids.items()
        if len(targets) > 1
    }
    unmatched_input_ids = genes - primary_ids - set(alternate_ids)
    common_summary = {
        "location_rows": len(locations),
        "location_primary_unique_ids": len(primary_ids),
        "location_alternate_unique_ids": len(alternate_ids),
        "location_ambiguous_alternate_ids": len(ambiguous_alternates),
        "alternate_identifier_column": alternate_column,
        "alternate_identifiers_available": bool(alternate_ids),
        "direct_primary_matches": len(primary_matches),
        "alternate_candidate_matches": len(alternate_candidate_matches),
        "alternate_ids_matched": len(alternate_matches),
        "direct_primary_input_fraction": (
            len(primary_matches) / len(genes) if genes else 0.0
        ),
        "alternate_input_fraction": (
            len(alternate_candidate_matches) / len(genes) if genes else 0.0
        ),
        "alternate_reference_fraction": (
            len(alternate_candidate_matches) / len(alternate_ids)
            if alternate_ids else 0.0
        ),
        "ambiguous_input_alternate_ids": len(
            alternate_matches & ambiguous_alternates
        ),
        "unmatched_input_ids": len(unmatched_input_ids),
        "input_resolved_ids": len(primary_matches) + len(alternate_matches),
        "input_resolution_fraction": (
            (len(primary_matches) + len(alternate_matches)) / len(genes)
            if genes else 0.0
        ),
    }
    primary_comparison = min(len(primary_ids), len(genes))
    primary_fraction = (
        len(primary_matches) / primary_comparison if primary_comparison else 0.0
    )
    common_summary.update({
        "primary_comparison_unique_ids": primary_comparison,
        "primary_match_fraction": primary_fraction,
    })
    if primary_fraction >= minimum_overlap:
        return gene_sets, primary_ids, {
            **common_summary,
            "alternate_ids_matched": 0,
            "ambiguous_input_alternate_ids": 0,
            "unmatched_input_ids": len(genes - primary_ids),
            "input_resolved_ids": len(primary_matches),
            "input_resolution_fraction": (
                len(primary_matches) / len(genes) if genes else 0.0
            ),
            "identifier_source": "primary_gene_id",
            "input_unique_ids": len(genes),
            "matched_unique_ids": len(primary_matches),
            "comparison_unique_ids": primary_comparison,
            "match_fraction": primary_fraction,
            "pathway_identifiers_modified": False,
        }

    if allow_alternate_reference_ids:
        selected_reference_ids = set(alternate_ids)
        matched = genes & selected_reference_ids
        comparison_size = min(len(selected_reference_ids), len(genes))
        match_fraction = len(matched) / comparison_size if comparison_size else 0.0
        if enforce and match_fraction < minimum_overlap:
            alternate_label = (
                "alternate column %d" % alternate_column
                if alternate_column is not None and alternate_ids else
                "an available alternate identifier column"
            )
            raise MagmaError(
                "Gene-set identifiers match neither primary gene-location column 1 "
                "(%d/%d; %.2f%%) nor %s (%d/%d; %.2f%%); the "
                "configured minimum is %.2f%%."
                % (
                    len(primary_matches), primary_comparison,
                    primary_fraction * 100, alternate_label,
                    len(matched), comparison_size,
                    match_fraction * 100, minimum_overlap * 100,
                )
            )
        return gene_sets, selected_reference_ids, {
            **common_summary,
            "identifier_source": "alternate_gene_id",
            "input_unique_ids": len(genes),
            "matched_unique_ids": len(matched),
            "comparison_unique_ids": comparison_size,
            "match_fraction": match_fraction,
            "pathway_identifiers_modified": False,
        }

    if enforce:
        raise MagmaError(
            "Gene-set identifiers match primary gene-location column 1 for only "
            "%d/%d unique IDs (%.2f%%), below the configured %.2f%% minimum. "
            "Pathway identifiers are never rewritten."
            % (
                len(primary_matches), primary_comparison, primary_fraction * 100,
                minimum_overlap * 100,
            )
        )
    return gene_sets, primary_ids, {
        **common_summary,
        "identifier_source": "primary_gene_id",
        "input_unique_ids": len(genes),
        "matched_unique_ids": len(primary_matches),
        "comparison_unique_ids": primary_comparison,
        "match_fraction": primary_fraction,
        "pathway_identifiers_modified": False,
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
    represented_gene_counts = [
        len(reference_gene_ids & values) for values in gene_set_values
    ]
    usable_sets = sum(count > 0 for count in represented_gene_counts)
    result = {
        "reference_unique_ids": len(reference_gene_ids),
        "geneset_unique_ids": len(set_ids),
        "overlapping_unique_ids": len(overlap),
        "comparison_unique_ids": comparison_size,
        "overlap_fraction_of_smaller_universe": fraction,
        "minimum_required_overlap": minimum_overlap,
        "total_gene_sets": gene_sets.height,
        "gene_sets_with_at_least_one_reference_gene": usable_sets,
        "gene_sets_with_exactly_one_reference_gene": sum(
            count == 1 for count in represented_gene_counts
        ),
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
    *,
    analysis_scope: dict | None = None,
    excluded_output: str | Path | None = None,
) -> dict:
    """Validate MAGMA variant inputs and apply configured reference/scope policies."""
    columns = module_config.input
    policy = module_config.snp_harmonisation
    resolve_variants = policy.resolve_variants_to_reference
    pval = _read_table(
        pval_file,
        "MAGMA p-value file",
        columns.table_delimiter_pattern,
        logger,
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
        logger,
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
    reference_by_id, reference_variant_count = read_reference_bim_matches(
        ld_ref, columns, set(working[variant]),
    )
    reference_ids = set(reference_by_id.index)
    id_matched = working[variant].isin(reference_ids)
    reference_id_match_rows = int(id_matched.sum())
    unmatched_rows = int((~id_matched).sum())
    reference_unique_id_matches = int(
        working.loc[id_matched, variant].nunique()
    )
    overlap_fraction = (
        reference_unique_id_matches / input_unique if input_unique else 0.0
    )
    if (
        reference_unique_id_matches == 0
        or overlap_fraction < policy.minimum_overlap_fraction
    ):
        bim_identifier_column = columns.bim_columns.index("variant_id") + 1
        raise MagmaError(
            "Only %d/%d unique GWAS variant identifiers (%.2f%%) occur in BIM "
            "column %d; the configured minimum is %.2f%%. Confirm that the "
            "formatter identifier convention and LD reference are compatible."
            % (
                reference_unique_id_matches,
                input_unique,
                overlap_fraction * 100,
                bim_identifier_column,
                policy.minimum_overlap_fraction * 100,
            )
        )
    working = working.loc[id_matched].copy() if resolve_variants else working.copy()

    excluded_chromosome_rows = 0
    excluded_mhc_rows = 0
    scope = analysis_scope or {
        "mhc_policy": "include",
        "exclude_mhc_snps": False,
        "mhc_region": None,
        "mhc_source": None,
        "exclude_chromosomes": [],
    }
    chromosome_mask = working["CHR_NORM"].isin(scope["exclude_chromosomes"])
    excluded_chromosome_counts = {
        chromosome: int(working["CHR_NORM"].eq(chromosome).sum())
        for chromosome in scope["exclude_chromosomes"]
    }
    excluded_chromosome_rows = int(chromosome_mask.sum())
    mhc_mask = pd.Series(False, index=working.index)
    region = scope["mhc_region"]
    if scope.get("exclude_mhc_snps", False):
        if region is None:
            raise MagmaError("Resolved MAGMA MHC region is missing")
        mhc_mask = (
            working["CHR_NORM"].eq(region["chromosome"])
            & working["BP_NORM"].between(
                region["start"], region["end"], inclusive="both",
            )
        )
    # Chromosome exclusions are attributed first; this keeps assigned reasons
    # mutually exclusive even when a custom MHC region lies on an excluded contig.
    assigned_mhc_mask = mhc_mask & ~chromosome_mask
    excluded_mhc_rows = int(assigned_mhc_mask.sum())
    excluded_mask = chromosome_mask | assigned_mhc_mask
    excluded = working.loc[excluded_mask].copy()
    if not excluded.empty:
        reason_column = module_config.exclusion_reporting.reason_column
        excluded[reason_column] = np.where(
            chromosome_mask.loc[excluded.index],
            module_config.exclusion_reporting.chromosome_reason,
            module_config.exclusion_reporting.mhc_reason,
        )
    if excluded_output is not None:
        destination = Path(excluded_output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        report = pd.DataFrame({
            variant: excluded.get(variant),
            columns.chromosome_column: excluded.get("CHR_NORM"),
            columns.position_column: excluded.get("BP_NORM"),
            columns.reference_allele_column: excluded.get("REF_NORM"),
            columns.alternate_allele_column: excluded.get("ALT_NORM"),
            module_config.exclusion_reporting.reason_column: excluded.get(
                module_config.exclusion_reporting.reason_column
            ),
        })
        report.to_csv(
            destination,
            sep=columns.output_table_delimiter,
            index=False,
        )
    working = working.loc[~excluded_mask].copy()
    if working.empty:
        raise MagmaError(
            "MAGMA chromosome and MHC policies excluded every eligible variant"
        )

    resolved = working
    duplicate_mask = resolved[variant].duplicated(keep=False)
    duplicate_rows_detected = int(duplicate_mask.sum())
    duplicate_ids = resolved.loc[duplicate_mask, variant].drop_duplicates()
    duplicate_groups = int(len(duplicate_ids))
    duplicate_policy = policy.duplicate_policy
    if duplicate_groups and duplicate_policy == "err":
        raise MagmaError(
            "MAGMA input contains %d duplicated SNP IDs across %d rows after "
            "optional LD-reference filtering; examples: %s. Use unique inputs, "
            "or explicitly select --duplicate-policy lowest_p or remove."
            % (
                duplicate_groups,
                duplicate_rows_detected,
                duplicate_ids.head(5).tolist(),
            )
        )
    if duplicate_policy == "lowest_p":
        resolved = resolved.sort_values(
            [variant, "_P_NUMERIC", "_ROW_ORDER"], kind="mergesort",
        ).drop_duplicates(variant, keep="first")
    elif duplicate_policy == "remove":
        resolved = resolved.loc[~duplicate_mask].copy()
    # Schema validation makes any other policy unreachable.
    elif duplicate_policy != "err":
        raise MagmaError("Unsupported MAGMA duplicate policy: %s" % duplicate_policy)
    duplicate_rows_removed = int(len(working) - len(resolved))
    resolved = resolved.sort_values("_ROW_ORDER", kind="mergesort")
    if resolved.empty:
        raise MagmaError(
            "MAGMA duplicate policy %s removed every eligible variant; no analysis "
            "input remains." % duplicate_policy
        )

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
        "not_in_reference_rows": unmatched_rows,
        "overlap_fraction": overlap_fraction,
        "duplicate_policy": duplicate_policy,
        "duplicate_groups_detected": duplicate_groups,
        "duplicate_rows_detected": duplicate_rows_detected,
        "duplicate_rows_removed": duplicate_rows_removed,
        "duplicate_groups_resolved_by_lowest_p": (
            duplicate_groups if duplicate_policy == "lowest_p" else 0
        ),
        "duplicate_groups_removed": (
            duplicate_groups if duplicate_policy == "remove" else 0
        ),
        "mhc_policy": scope["mhc_policy"],
        "mhc_region": scope["mhc_region"],
        "mhc_region_source": scope["mhc_source"],
        "excluded_chromosomes": scope["exclude_chromosomes"],
        "excluded_chromosome_rows": excluded_chromosome_rows,
        "excluded_chromosome_counts": excluded_chromosome_counts,
        "excluded_mhc_rows": excluded_mhc_rows,
        "excluded_scope_rows": excluded_chromosome_rows + excluded_mhc_rows,
        "retained_rows": len(harmonised_pval),
    }
    logger.record("RESULT", "magma_variant_preparation", **qc)
    return {
        "pval_file": str(pval_destination),
        "snp_loc_file": str(location_destination),
        "excluded_variants": str(excluded_output) if excluded_output is not None else None,
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
        logger,
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

    global_columns = {
        method: schema.global_correction_column_pattern.format(method=method)
        for method in module_config.multiple_testing.global_methods
    }
    family_columns = {
        (name, method): schema.family_correction_column_pattern.format(
            family=name, method=method,
        )
        for name, family in module_config.multiple_testing.families.items()
        for method in family.methods
    }
    primary_columns = [
        schema.primary_correction_method_column,
        schema.primary_adjusted_p_value_column,
        schema.primary_significant_column,
    ]
    _require_available_result_columns(
        frame,
        [*global_columns.values(), *family_columns.values(), *primary_columns],
        "gene-set correction",
    )

    values = frame[p_column].to_numpy(dtype=float)
    for method, column in global_columns.items():
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
            column = family_columns[(name, method)]
            frame[column] = np.nan
            if len(family_values):
                frame.loc[selected, column] = adjust_p_values(family_values, method)
            logger.record(
                "PARAM", "multiple_testing_family",
                family=name, pattern=family.pattern, method=method,
                tests=len(family_values),
            )
    primary_method = module_config.multiple_testing.primary_method
    primary_source_column = global_columns[primary_method]
    threshold = module_config.multiple_testing.reporting_significance_threshold
    frame[schema.primary_correction_method_column] = primary_method
    frame[schema.primary_adjusted_p_value_column] = frame[primary_source_column]
    frame[schema.primary_significant_column] = (
        frame[schema.primary_adjusted_p_value_column] <= threshold
    )
    primary_significant = int(frame[schema.primary_significant_column].sum())
    logger.record(
        "RESULT",
        "primary_gene_set_correction",
        family="all_gene_sets",
        method=primary_method,
        method_label=module_config.multiple_testing.reporting_method_labels[
            primary_method
        ],
        threshold=threshold,
        tests=len(values),
        significant=primary_significant,
        source_column=primary_source_column,
        adjusted_p_value_column=schema.primary_adjusted_p_value_column,
        significance_column=schema.primary_significant_column,
    )
    frame = frame.sort_values(p_column, kind="mergesort")
    destination = write_delimited_report(
        frame.to_dict(orient="records"),
        output_file,
        fieldnames=frame.columns.tolist(),
        delimiter=schema.report_delimiter,
        null_value=schema.report_null_value,
    )
    logger.record(
        "OUTPUT",
        "corrected_gene_sets",
        path=str(destination),
        rows=len(frame),
        primary_method=primary_method,
        primary_significant=primary_significant,
    )
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
    *,
    gene_location_file: str | Path | None = None,
    gene_location_has_header: bool | None = None,
) -> pl.DataFrame:
    """Join configured gene metadata and add adjusted p-values to MAGMA results."""
    schema = module_config.result_schema
    source = Path(genes_out_file).expanduser().resolve()
    frame = _read_table(
        source,
        "MAGMA gene output",
        schema.table_delimiter_pattern,
        logger,
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
    frame[gene_column] = normalized_genes

    if gene_location_file is not None:
        locations = read_gene_locations(
            gene_location_file,
            module_config,
            has_header=gene_location_has_header,
        )
        annotation_columns = {
            "chromosome": schema.gene_reference_chromosome_column,
            "start": schema.gene_reference_start_column,
            "end": schema.gene_reference_end_column,
            "strand": schema.gene_reference_strand_column,
            "alternate_gene_id": schema.gene_reference_alternate_id_column,
        }
        _require_available_result_columns(
            frame, list(annotation_columns.values()), "gene-location annotation",
        )
        location_frame = pd.DataFrame.from_records(
            [
                {
                    gene_column: gene,
                    annotation_columns["chromosome"]: chromosome,
                    annotation_columns["start"]: start,
                    annotation_columns["end"]: end,
                    annotation_columns["strand"]: strand,
                    annotation_columns["alternate_gene_id"]: alternate,
                }
                for gene, (
                    chromosome, start, end, strand, alternate,
                ) in locations.items()
            ]
        )
        frame = frame.merge(
            location_frame, on=gene_column, how="left", validate="one_to_one",
        )
        missing_annotations = int(
            frame[annotation_columns["chromosome"]].isna().sum()
        )
        if missing_annotations:
            raise MagmaError(
                "%d MAGMA gene results are absent from the configured gene-location "
                "reference; examples: %s"
                % (
                    missing_annotations,
                    frame.loc[
                        frame[annotation_columns["chromosome"]].isna(), gene_column
                    ].head(5).tolist(),
                )
            )
        logger.record(
            "RESULT",
            "gene_result_location_annotation",
            gene_location_file=str(
                Path(gene_location_file).expanduser().resolve()
            ),
            annotated_genes=len(frame),
            missing_annotations=missing_annotations,
            columns=annotation_columns,
        )

    correction_columns = {
        method: schema.global_correction_column_pattern.format(method=method)
        for method in module_config.multiple_testing.gene_methods
    }
    _require_available_result_columns(
        frame, list(correction_columns.values()), "gene correction",
    )
    values = frame[p_column].to_numpy(dtype=float)
    for method, column in correction_columns.items():
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
    *,
    primary_method: str | None = None,
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
    if primary_method is not None:
        if primary_method not in adjusted:
            raise MagmaError(
                "Configured primary correction %s was not calculated" % primary_method
            )
        primary_count = adjusted[primary_method]
        secondary_methods = [
            method for method in methods if method != primary_method
        ]
        message = (
            f"{tested:,} {subject} tested at the configured p ≤ {threshold:g} "
            f"threshold: {primary_count:,} significant by the primary global "
            f"{labels[primary_method]} correction; {nominal:,} nominally "
            "significant (descriptive only)."
        )
        metrics = {
            "tested": tested,
            "significance_threshold": threshold,
            "nominal_significant": nominal,
            "adjusted_significant": adjusted,
            "primary_correction": primary_method,
            "primary_significant": primary_count,
        }
        fields = [
            ("count", f"{plural.capitalize()} tested", tested),
            ("analysis", "Reporting threshold", f"p ≤ {threshold:g}"),
            (
                "analysis",
                "Primary correction",
                f"Global {labels[primary_method]}",
            ),
            (
                "analysis",
                "Significant by primary correction",
                primary_count,
            ),
            ("info", "Nominally significant (descriptive)", nominal),
        ]
        if secondary_methods:
            fields.append(
                (
                    "info",
                    "Other global corrections in results",
                    ", ".join(labels[method] for method in secondary_methods),
                )
            )
        return message, metrics, fields

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
    return_metadata: bool = False,
) -> pl.DataFrame | tuple[pl.DataFrame, dict]:
    """Build configured gene-set columns from the shared validated reader."""
    source = require_nonempty_file(
        gene_set_path, "gene-set file", error_type=MagmaError,
    )
    configured_format = input_format or module_config.gene_sets.input_format
    schema = module_config.result_schema
    records = []
    membership = module_config.gene_sets
    is_membership = configured_format == "membership"

    def append_record(row):
        joined_genes = ",".join(row.genes)
        records.append({
            schema.gene_set_full_name_column: row.name,
            schema.report_gene_set_description_column: row.description or None,
            schema.report_source_input_genes_column: joined_genes,
            schema.report_input_genes_column: joined_genes,
        })

    validated = validate_gene_set_source(
        source,
        GeneSetFormat(
            input_format=configured_format,
            delimiter_pattern=(membership.membership_delimiter_pattern if is_membership
                               else module_config.input.table_delimiter_pattern),
            comment_prefix=_RESULT_COMMENT_PREFIX,
            compressed=True,
            empty_gene_policy="reject" if is_membership else "drop",
            duplicate_gene_policy="deduplicate",
            reject_name_whitespace=False,
            membership_has_header=membership.membership_has_header if is_membership else False,
            membership_set_column=membership.membership_set_column if is_membership else None,
            membership_gene_column=membership.membership_gene_column if is_membership else None,
        ),
        error_type=MagmaError,
        on_row=append_record,
    )
    detected_format = validated.detected_format
    logger.record(
        "DECIDE",
        "gene_set_input_format",
        configured=configured_format,
        detected=detected_format,
        gene_sets=len(records),
    )
    frame = pl.DataFrame(
        records,
        schema={
            schema.gene_set_full_name_column: pl.String,
            schema.report_gene_set_description_column: pl.String,
            schema.report_source_input_genes_column: pl.String,
            schema.report_input_genes_column: pl.String,
        },
    )
    if return_metadata:
        return frame, {
            "configured_format": configured_format,
            "detected_format": detected_format,
            "gene_sets": len(records),
        }
    return frame


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
        logger,
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
                schema.report_untested_genes_column: (
                    ",".join(gene for gene in values if gene not in gene_p) or None
                ),
                schema.report_untested_gene_count_column: sum(
                    gene not in gene_p for gene in values
                ),
            }
        )
    annotated = (
        corrected.join(
            pl.concat([gene_sets, pl.DataFrame(details)], how="horizontal"),
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


def _prepare_gene_set_plans(mapping_gene_sets, module, logger) -> dict[str, dict]:
    """Validate optional pathway inputs before any MAGMA scientific command."""
    plans: dict[str, dict] = {}
    for name in module.mapping.selected:
        source = mapping_gene_sets[name]
        if source is None:
            plans[name] = {"status": "not_requested", "reason": None}
            continue
        definition = module.mapping.definitions[name]
        gene_set_format = (
            module.gene_sets.input_format
            if definition.gene_set_format == "auto"
            else definition.gene_set_format
        )
        parsed, input_metadata = parse_gene_set_file(
            source,
            module,
            logger,
            input_format=gene_set_format,
            return_metadata=True,
        )
        if (
            input_metadata["detected_format"] == "membership"
            and module.gene_sets.membership_has_header
        ):
            raise MagmaError(
                "MAGMA column-based --set-annot input cannot contain a header. "
                "Provide the original headerless membership file; PostGWAS does "
                "not rewrite pathway inputs."
            )
        minimum = (
            definition.minimum_gene_id_overlap_fraction
            or module.gene_sets.minimum_gene_id_overlap_fraction
        )
        if definition.method in {"positional", "n_magma"}:
            location = (
                definition.gene_location_file or module.input.gene_location_file
            )
            if location is None:
                raise MagmaError(
                    "Mapping %s requires a gene-location file" % name
                )
            reference = require_nonempty_file(
                location,
                "%s gene-location file" % name,
                error_type=MagmaError,
            )
            parsed, reference_ids, resolution = resolve_gene_set_identifiers(
                reference,
                parsed,
                minimum,
                module,
                enforce=False,
                allow_alternate_reference_ids=(definition.method == "positional"),
            )
        else:
            reference = require_nonempty_file(
                definition.gene_annotation_file,
                "%s source gene annotation" % definition.display_name,
                error_type=MagmaError,
            )
            reference_ids = _gene_ids_from_table(
                reference,
                "%s source gene annotation" % definition.display_name,
                module.input.table_delimiter_pattern,
                0,
                has_header=False,
                comment_prefix=module.annotation_validation.comment_prefix,
            )
            resolution = {
                "identifier_source": "annotation_gene_id",
                "pathway_identifiers_modified": False,
            }
        validation = validate_gene_set_identifier_overlap(
            reference_ids,
            parsed,
            None,
            module.result_schema.report_input_genes_column,
        )
        validation["identifier_resolution"] = resolution
        observed = validation["overlap_fraction_of_smaller_universe"]
        if observed < minimum:
            alternate_column = resolution.get("alternate_identifier_column")
            alternate_reason = (
                "The optional alternate gene-identifier column is absent or "
                "contains no usable identifiers. "
                if not resolution.get("alternate_identifiers_available", True)
                else ""
            )
            incompatibility = (
                "Gene-set identifiers for %s overlap the mapping reference for only "
                "%d/%d unique IDs (%.2f%%), below the configured %.2f%% minimum. "
                "%sUse a gene-set file "
                "whose identifiers match gene-location column 1 or its optional "
                "%s."
                % (
                    definition.display_name,
                    validation["overlapping_unique_ids"],
                    validation["comparison_unique_ids"],
                    observed * 100,
                    minimum * 100,
                    alternate_reason,
                    (
                        "column %d aliases" % alternate_column
                        if alternate_column is not None else
                        "alternate identifier column"
                    ),
                )
            )
            if module.gene_sets.identifier_mismatch_action == "error":
                raise MagmaError(
                    "%s The configured action is error, so the complete MAGMA "
                    "run was stopped." % incompatibility
                )
            reason = (
                "Competitive gene-set analysis was not performed. %s "
                "Gene-association analysis will continue." % incompatibility
            )
            logger.warn(reason)
            logger.record(
                "SKIP",
                "magma_gene_set_analysis",
                mapping=name,
                reason="incompatible_gene_identifiers",
                message=reason,
                **validation,
            )
            plans[name] = {
                "status": "skipped",
                "reason": reason,
                "pathway_file": str(source),
                "gene_reference_file": str(reference),
                "minimum_overlap": minimum,
                "input_metadata": input_metadata,
                "validation": validation,
            }
            continue
        plans[name] = {
            "status": "ready",
            "reason": None,
            "pathway_file": str(source),
            "gene_reference_file": str(reference),
            "minimum_overlap": minimum,
            "input_metadata": input_metadata,
            "parsed_gene_sets": parsed,
            "validation": validation,
        }
    return plans


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
            schema.report_mhc_policy_column: analyses[name]["analysis_scope"][
                "mhc_policy"
            ],
            schema.report_mhc_region_column: analyses[name]["analysis_scope"][
                "mhc_region_label"
            ],
            schema.report_excluded_chromosomes_column: ",".join(
                analyses[name]["analysis_scope"]["exclude_chromosomes"]
            ),
            schema.report_excluded_units_column: analyses[name]["gene_scope"][
                "excluded_units"
            ],
            schema.report_gene_set_status_column: analyses[name][
                "gene_set_analysis"
            ]["status"],
            schema.report_gene_set_reason_column: analyses[name][
                "gene_set_analysis"
            ]["reason"],
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
    *,
    references: MagmaReferencePreflight | None = None,
    include_gene_sets: bool = True,
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
    reference = references or preflight_magma_references(
        dataset_id, configuration, logger, include_gene_sets=include_gene_sets,
    )
    return MagmaPreflight(
        snp_location_file=snp_location_file,
        p_value_file=p_value_file,
        ld_reference_prefix=reference.ld_reference_prefix,
        executable=reference.executable,
        version=reference.version,
        gene_set_plans=reference.gene_set_plans,
        analysis_scope=reference.analysis_scope,
    )


def preflight_magma_references(
    dataset_id: str,
    configuration,
    logger,
    *,
    include_gene_sets: bool = True,
) -> MagmaReferencePreflight:
    """Validate all MAGMA resources that precede GWAS-VCF extraction."""
    validate_filename_component(dataset_id, "dataset_id", error_type=MagmaError)
    ld_reference_prefix, executable, version, analysis_scope = (
        validate_magma_ld_reference_and_runtime(configuration, logger)
    )
    gene_set_plans = validate_magma_gene_set_references(
        configuration, logger, include_gene_sets=include_gene_sets,
    )
    return MagmaReferencePreflight(
        ld_reference_prefix=ld_reference_prefix,
        executable=executable,
        version=version,
        gene_set_plans=gene_set_plans,
        analysis_scope=analysis_scope,
    )


def validate_magma_ld_reference_and_runtime(configuration, logger):
    """Validate the PLINK LD reference and configured MAGMA executable."""
    module = configuration.modules.magma
    inputs = module.input
    if configuration.execution.random_seed < 1:
        raise MagmaError("MAGMA requires a positive execution.random_seed")
    if inputs.ld_reference_prefix is None:
        raise MagmaError("Required MAGMA input is missing: LD-reference prefix")
    ld_reference_prefix = str(Path(inputs.ld_reference_prefix).expanduser().resolve())
    validate_plink_files(
        ld_reference_prefix, inputs.required_reference_extensions, error_type=MagmaError,
        missing_message="LD reference is incomplete or empty",
    )
    executable = resolve_executable(
        configuration.resources.executables.magma,
        "MAGMA executable",
        error_type=MagmaError,
    )
    version = require_supported_magma(executable, module, logger)
    return ld_reference_prefix, executable, version, _analysis_scope(configuration)


def validate_magma_gene_set_references(
    configuration, logger, *, include_gene_sets: bool = True,
):
    """Validate pathway files against each selected mapping's gene reference."""
    module = configuration.modules.magma
    inputs = module.input
    mapping_gene_sets = {
        name: (
            _mapping_gene_set(module.mapping.definitions[name], inputs)
            if include_gene_sets else None
        )
        for name in module.mapping.selected
    }
    return _prepare_gene_set_plans(mapping_gene_sets, module, logger)


def run_magma_analysis(
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger,
    *,
    preflight: MagmaPreflight | None = None,
    pipeline_progress=None,
    pipeline_stage_numbers: Mapping[str, int] | None = None,
    formatter_result=None,
    variant_id_observation=None,
) -> dict:
    """Run positional and configured functional MAGMA mappings independently."""
    include_gene_sets = (
        pipeline_stage_numbers is None
        or "pathway_analysis" in pipeline_stage_numbers
    )
    prepared = preflight or preflight_magma_analysis(
        dataset_id,
        configuration,
        logger,
        include_gene_sets=include_gene_sets,
    )
    module = configuration.modules.magma
    inputs = module.input
    snp_loc = prepared.snp_location_file
    pval = prepared.p_value_file
    ld_ref = prepared.ld_reference_prefix
    executable = prepared.executable
    version = prepared.version
    gene_set_plans = prepared.gene_set_plans
    analysis_scope = prepared.analysis_scope
    active_pipeline_stages = dict(pipeline_stage_numbers or {})
    if pipeline_progress is not None and not active_pipeline_stages:
        active_pipeline_stages = {
            key: number
            for number, key in enumerate(MAGMA_PIPELINE_STAGE_KEYS, 1)
        }

    def start_pipeline_stage(key: str) -> None:
        if pipeline_progress is not None and key in active_pipeline_stages:
            pipeline_progress.start(active_pipeline_stages[key])

    def complete_pipeline_stage(key: str, **values) -> None:
        if pipeline_progress is not None and key in active_pipeline_stages:
            pipeline_progress.complete(active_pipeline_stages[key], **values)

    total_steps = 1 + sum(
        3 + (2 if gene_set_plans[name]["status"] == "ready" else 0)
        for name in module.mapping.selected
    )
    shared_paths = resolve_magma_output_paths(
        output_directory, dataset_id, module, module.mapping.primary,
    )
    for path in shared_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    step_number = 1
    if pipeline_progress is not None and pipeline_progress.current == 0:
        start_pipeline_stage("variant_inputs")
    variant_stage_title = (
        "Assess GWAS-BIM identifier overlap and apply configured variant policies"
    )
    with logger.step(
        step_number, total_steps, variant_stage_title, "prepare_magma_variant_inputs",
    ) as step:
        variant_preparation = (
            formatter_result.get("variant_preparation")
            if isinstance(formatter_result, Mapping) else None
        )
        if variant_preparation is None:
            variant_preparation = prepare_magma_variant_inputs(
                pval,
                snp_loc,
                ld_ref,
                shared_paths["harmonised_p_values"],
                shared_paths["harmonised_snp_locations"],
                module,
                logger,
                analysis_scope=analysis_scope,
                excluded_output=shared_paths["excluded_variants"],
            )
        else:
            shared_paths["harmonised_p_values"] = Path(
                variant_preparation["pval_file"]
            )
            shared_paths["harmonised_snp_locations"] = Path(
                variant_preparation["snp_loc_file"]
            )
            logger.record(
                "SKIP",
                "prepare_magma_variant_inputs",
                reason="completed_by_formatter_in_memory",
                p_value_file=variant_preparation["pval_file"],
                snp_location_file=variant_preparation["snp_loc_file"],
            )
        qc = variant_preparation["qc"]
        removed = qc["input_rows"] - qc["retained_rows"]
        step.set_rows(qc["retained_rows"], removed=removed)
        variant_fields = [
            ("count", "Formatter rows", qc["input_rows"]),
            ("success", "Prepared variants", qc["retained_rows"]),
            ("count", "BIM reference variants", qc["reference_variant_count"]),
            ("success", "GWAS identifiers found in BIM", (
                "%s/%s (%.2f%%)"
                % (
                    f"{qc['reference_unique_id_matches']:,}",
                    f"{qc['input_unique_variants']:,}",
                    qc["overlap_fraction"] * 100,
                )
            )),
            (
                "analysis", "Minimum identifier overlap",
                "%.2f%%" % (
                    module.snp_harmonisation.minimum_overlap_fraction * 100
                ),
            ),
            (
                "info", "BIM intersection",
                "applied" if qc["reference_intersection_enabled"] else "not requested",
            ),
            (
                "warning" if qc["not_in_reference_rows"] else "success",
                "Rows absent from BIM", qc["not_in_reference_rows"],
            ),
        ]
        variant_fields.extend(
            [
                ("analysis", "MHC policy", qc["mhc_policy"]),
                (
                    "genetic", "MHC region",
                    (
                        "%s:%s-%s"
                        % (
                            qc["mhc_region"]["chromosome"],
                            qc["mhc_region"]["start"],
                            qc["mhc_region"]["end"],
                        )
                        if qc["mhc_region"] is not None else "not used"
                    ),
                ),
                (
                    "genetic", "Excluded chromosomes",
                    ", ".join(qc["excluded_chromosomes"]) or "none",
                ),
                (
                    (
                        "loss" if qc["excluded_chromosome_counts"]
                        else "success"
                    ),
                    "Variants excluded by chromosome",
                    ", ".join(
                        "%s=%s" % (chromosome, count)
                        for chromosome, count
                        in qc["excluded_chromosome_counts"].items()
                    ) or "none",
                ),
                (
                    "loss" if qc["excluded_mhc_rows"] else "success",
                    "Variants excluded from MHC",
                    qc["excluded_mhc_rows"],
                ),
            ]
        )
        variant_fields.extend(
            [
                ("info", "Duplicate policy", qc["duplicate_policy"]),
                (
                    "count", "Duplicate ID groups detected",
                    qc["duplicate_groups_detected"],
                ),
                (
                    "info",
                    (
                        "Rows removed after keeping lowest P"
                        if qc["duplicate_policy"] == "lowest_p"
                        else "Rows removed with duplicated IDs"
                    ),
                    qc["duplicate_rows_removed"],
                ),
            ]
        )
        variant_outcome = (
            "%s formatter rows; %s BIM-matched identifiers retained; %s rows "
            "were not retained."
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
    if pipeline_progress is not None:
        complete_pipeline_stage(
            "variant_inputs",
            outcome_fields=magma_variant_input_outcome_fields(
                configuration,
                formatter_result or {"rows_in": qc["input_rows"]},
                variant_id_observation or {},
                qc,
                {
                    "snp_loc_file": variant_preparation["snp_loc_file"],
                    "pval_file": variant_preparation["pval_file"],
                },
            ),
        )
    step_number += 1

    analyses: dict[str, dict] = {}
    for name in module.mapping.selected:
        definition = module.mapping.definitions[name]
        gene_set_plan = gene_set_plans[name]
        parsed_gene_sets = gene_set_plan.get("parsed_gene_sets")
        gene_id_validation = gene_set_plan.get("validation")
        tested_gene_coverage = None
        location_reference = None
        pathway_location_preparation = None
        gene_significance = None
        pathway_significance = None
        paths = resolve_magma_output_paths(
            output_directory, dataset_id, module, name,
        )
        for path in paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        paths["harmonised_p_values"] = shared_paths["harmonised_p_values"]
        paths["harmonised_snp_locations"] = shared_paths["harmonised_snp_locations"]

        start_pipeline_stage("gene_annotation")
        with logger.step(
            step_number,
            total_steps,
            "%s · create the MAGMA SNP-to-gene annotation"
            % definition.display_name,
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
                resolution = (
                    (gene_set_plan.get("validation") or {}).get(
                        "identifier_resolution"
                    ) or {}
                )
                if (
                    gene_set_plan["status"] == "ready"
                    and resolution.get("identifier_source") == "alternate_gene_id"
                ):
                    location_reference, pathway_location_preparation = (
                        write_pathway_compatible_gene_locations(
                            location_reference,
                            paths["pathway_compatible_gene_locations"],
                            module,
                            logger,
                        )
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
                        "--annotate",
                        _annotation_window_argument(upstream, downstream),
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
                    "pathway_compatible_gene_location": (
                        pathway_location_preparation
                    ),
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
                        _annotation_window_argument(
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
                        _annotation_window_argument(
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
            scope_result = prepare_scoped_gene_annotation(
                annotation_file,
                paths["scoped_annotation"],
                paths["excluded_genes"],
                paths["exclusion_summary"],
                module,
                analysis_scope,
                logger,
            )
            source_annotation_file = annotation_file
            annotation_file = paths["scoped_annotation"]
            annotated_units = scope_result["retained_units"]
            annotation_fields = [
                ("analysis", "Mapping method", definition.method),
                ("genetic", "Annotated units", annotated_units),
                (
                    (
                        "loss"
                        if scope_result["excluded_by_reason"][
                            module.exclusion_reporting.chromosome_reason
                        ] else "success"
                    ),
                    "Units excluded by chromosome",
                    scope_result["excluded_by_reason"][
                        module.exclusion_reporting.chromosome_reason
                    ],
                ),
                (
                    (
                        "loss"
                        if scope_result["excluded_by_reason"][
                            module.exclusion_reporting.mhc_reason
                        ] else "success"
                    ),
                    "Units excluded from MHC",
                    scope_result["excluded_by_reason"][
                        module.exclusion_reporting.mhc_reason
                    ],
                ),
                ("info", "Biological context", definition.context),
            ]
            if pathway_location_preparation is not None:
                annotation_fields.extend([
                    (
                        "success", "Pathway gene identifiers",
                        "unchanged",
                    ),
                    (
                        "genetic", "Gene-location identifier",
                        "alternate ID from configured column %d"
                        % (
                            inputs.gene_location_columns.index(
                                "alternate_gene_id"
                            ) + 1
                        ),
                    ),
                    (
                        "count", "Derived gene-location rows",
                        pathway_location_preparation["retained_rows"],
                    ),
                    (
                        (
                            "warning"
                            if pathway_location_preparation[
                                "ambiguous_alternate_ids"
                            ] else "success"
                        ),
                        "Alternate IDs with multiple intervals",
                        pathway_location_preparation[
                            "ambiguous_alternate_ids"
                        ],
                    ),
                ])
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
                scope=scope_result,
            )
            step.output("gene_annotation", path=str(annotation_file))
        if pipeline_progress is not None:
            complete_pipeline_stage(
                "gene_annotation",
                outcome_fields=[
                    ("genetic", "Annotated units", annotated_units),
                    ("count", "Excluded units", scope_result["excluded_units"]),
                ],
            )
        step_number += 1

        if gene_set_plan["status"] == "ready":
            logger.record(
                "RESULT",
                "gene_id_compatibility",
                mapping=name,
                reference_universe="complete_mapping_reference",
                pathway_file=str(gene_set_plan["pathway_file"]),
                pathway_file_modified=False,
                **gene_id_validation,
            )

        start_pipeline_stage("gene_analysis")
        with logger.step(
            step_number,
            total_steps,
            "%s · calculate gene-association statistics" % definition.display_name,
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
        if pipeline_progress is not None:
            complete_pipeline_stage(
                "gene_analysis",
                outcome_fields=[
                    ("genetic", "Genes submitted", genes),
                    ("analysis", "Parallel workers", workers),
                ],
            )
        step_number += 1

        start_pipeline_stage("gene_results")
        with logger.step(
            step_number,
            total_steps,
            (
                "%s · map regulatory-element results to genes"
                if definition.method == "chrom_magma"
                else "%s · annotate gene results and adjust p-values for multiple testing"
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
                    logger,
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
                gene_significance = significance
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
                    paths["genes_out"],
                    paths["corrected_genes"],
                    module,
                    logger,
                    gene_location_file=location_reference,
                    gene_location_has_header=(
                        False if pathway_location_preparation is not None else None
                    ),
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
                gene_significance = significance
                if corrected_genes.height > genes:
                    raise MagmaError(
                        "MAGMA returned more gene results (%d) than annotated "
                        "genes submitted (%d) for %s"
                        % (corrected_genes.height, genes, definition.display_name)
                    )
                fields = [
                    ("genetic", "Genes submitted to MAGMA", genes),
                    *fields,
                    (
                        "warning", "Genes without a valid MAGMA result",
                        genes - corrected_genes.height,
                    ),
                ]
                step.set_rows(corrected_genes.height)
                step.outcome(outcome, fields=fields, mapping=name, **significance)
                step.output(
                    "annotated_gene_associations",
                    path=str(paths["corrected_genes"]),
                )
        if pipeline_progress is not None:
            complete_pipeline_stage(
                "gene_results",
                outcome_fields=[
                    ("genetic", "Genes tested", gene_significance["tested"]),
                    (
                        "info", "Nominally significant",
                        gene_significance.get("nominal_significant", "not applicable"),
                    ),
                ],
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
            "gene_id_type": (
                inputs.alternate_gene_id_type
                if pathway_location_preparation is not None
                else definition.gene_id_type
            ),
            "effective_gene_id_source": (
                ((gene_set_plan.get("validation") or {}).get(
                    "identifier_resolution"
                ) or {}).get("identifier_source", "mapping_definition")
            ),
            "annotation_source": definition.source_name,
            "annotation_version": definition.source_version,
            "result_statistic_type": definition.result_statistic_type,
            "result_statistic_interpretation": (
                definition.result_statistic_interpretation
            ),
            "gene_annotation": str(annotation_file),
            "source_gene_annotation": str(source_annotation_file),
            "gene_location_reference": (
                str(location_reference) if location_reference is not None else None
            ),
            "pathway_compatible_gene_location": (
                str(paths["pathway_compatible_gene_locations"])
                if pathway_location_preparation is not None else None
            ),
            "gene_scope": scope_result,
            "analysis_scope": {
                **analysis_scope,
                "mhc_region_label": (
                    "%s:%s-%s"
                    % (
                        analysis_scope["mhc_region"]["chromosome"],
                        analysis_scope["mhc_region"]["start"],
                        analysis_scope["mhc_region"]["end"],
                    )
                    if analysis_scope["mhc_region"] is not None else None
                ),
            },
            "magma_genes_prefix": str(paths["gene_prefix"]),
            "magma_genes_raw": str(paths["genes_raw"]),
            "magma_genes_out": str(paths["genes_out"]),
            "magma_genes_annotated": corrected_gene_path,
            "magma_gene_results": str(gene_result_path),
            "chrom_magma_gene_report": chrom_gene_report,
            "corrected_gene_results": result_frame,
            "annotation_validation": annotation_validation,
            "batching": {
                "annotated_units": genes,
                "batches": batches,
                "workers": workers,
            },
            "gene_significance": gene_significance,
            "gene_id_validation": gene_id_validation,
            "tested_gene_coverage": tested_gene_coverage,
            "gene_set_analysis": {
                "status": gene_set_plan["status"],
                "reason": gene_set_plan["reason"],
            },
        }

        if gene_set_plan["status"] == "ready":
            start_pipeline_stage("pathway_analysis")
            with logger.step(
                step_number,
                total_steps,
                "%s · run competitive pathway analysis with the original pathway file"
                % definition.display_name,
                "magma_gene_set_analysis",
            ) as step:
                set_annotation_arguments = [
                    "--set-annot", str(gene_set_plan["pathway_file"]),
                ]
                if gene_set_plan["input_metadata"]["detected_format"] == "membership":
                    set_annotation_arguments.append(
                        "col=%d,%d"
                        % (
                            module.gene_sets.membership_gene_column + 1,
                            module.gene_sets.membership_set_column + 1,
                        )
                    )
                _run_command(
                    [
                        executable,
                        "--gene-results", str(paths["genes_raw"]),
                        *set_annotation_arguments,
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
                    "%s submitted the original pathway file containing %s gene "
                    "sets using %s tested genes."
                    % (
                        definition.display_name,
                        f"{parsed_gene_sets.height:,}",
                        f"{tested_gene_coverage['reference_unique_ids']:,}",
                    ),
                    fields=[
                        (
                            "info", "Pathway file",
                            Path(gene_set_plan["pathway_file"]).name,
                        ),
                        (
                            "success", "Pathway file handling",
                            "passed directly to MAGMA without modification",
                        ),
                        (
                            "count", "Pathway records submitted",
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
                    pathway_file=str(gene_set_plan["pathway_file"]),
                    pathway_file_modified=False,
                    tested_gene_coverage=tested_gene_coverage,
                )
            if pipeline_progress is not None:
                complete_pipeline_stage(
                    "pathway_analysis",
                    outcome_fields=[
                        ("count", "Pathways submitted", parsed_gene_sets.height),
                        (
                            "genetic", "Study-tested genes",
                            tested_gene_coverage["reference_unique_ids"],
                        ),
                    ],
                )
            step_number += 1

            start_pipeline_stage("pathway_results")
            with logger.step(
                step_number,
                total_steps,
                "%s · annotate pathway results and adjust pathway p-values"
                % definition.display_name,
                "annotate_and_correct_gene_set_results",
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
                    primary_method=module.multiple_testing.primary_method,
                )
                pathway_significance = significance
                if corrected_sets.height > parsed_gene_sets.height:
                    raise MagmaError(
                        "MAGMA returned more pathway results (%d) than pathways "
                        "submitted (%d) for %s"
                        % (
                            corrected_sets.height,
                            parsed_gene_sets.height,
                            definition.display_name,
                        )
                    )
                fields = [
                    (
                        "count", "Pathways submitted to MAGMA",
                        parsed_gene_sets.height,
                    ),
                    *fields,
                    (
                        "warning", "Pathways without a valid MAGMA result",
                        parsed_gene_sets.height - corrected_sets.height,
                    ),
                ]
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
                    outcome,
                    fields=[
                        *fields,
                        (
                            "success", "Annotated pathway results",
                            corrected_sets.height,
                        ),
                    ],
                    mapping=name,
                    **significance,
                )
            if pipeline_progress is not None:
                complete_pipeline_stage(
                    "pathway_results",
                    outcome_fields=[
                        ("count", "Pathways tested", pathway_significance["tested"]),
                        (
                            "success", "Primary-adjusted significant",
                            pathway_significance["primary_significant"],
                        ),
                    ],
                )
            step_number += 1
            analysis_result.update(
                {
                    "pathway_file": str(gene_set_plan["pathway_file"]),
                    "magma_gene_sets_raw": str(paths["gene_sets_raw"]),
                    "magma_gene_sets_corrected": str(paths["corrected_gene_sets"]),
                    "magma_pathway": pathway,
                    "pathway_significance": pathway_significance,
                    "gene_set_analysis": {
                        "status": "completed",
                        "reason": None,
                    },
                }
            )
        elif (
            pipeline_progress is not None
            and "pathway_analysis" in active_pipeline_stages
        ):
            start_pipeline_stage("pathway_analysis")
            complete_pipeline_stage(
                "pathway_analysis",
                outcome=(
                    "Pathway analysis was not requested or gene identifiers "
                    "were incompatible"
                ),
            )
            start_pipeline_stage("pathway_results")
            complete_pipeline_stage(
                "pathway_results",
                outcome=(
                    "No pathway result required annotation or p-value adjustment"
                ),
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
        "analysis_scope": analysis_scope,
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
        "gene_set_analysis": primary["gene_set_analysis"],
        "batching": primary["batching"],
    }
    if primary["magma_genes_annotated"] is not None:
        result["magma_genes_annotated"] = primary["magma_genes_annotated"]
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
    "MagmaReferencePreflight",
    "annotate_gene_set_results",
    "correct_gene_p_values",
    "correct_gene_set_p_values",
    "prepare_magma_variant_inputs",
    "parse_gene_set_file",
    "preflight_magma_analysis",
    "preflight_magma_references",
    "require_supported_magma",
    "resolve_gene_set_identifiers",
    "resolve_magma_output_paths",
    "run_magma_analysis",
    "validate_magma_gene_set_references",
    "validate_magma_ld_reference_and_runtime",
]
