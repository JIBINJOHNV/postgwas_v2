"""Final MAGMA p-value and SNP-location export."""

from pathlib import Path

import polars as pl

from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path
from postgwas.modules.formatting.contracts import NAMED_OUTPUT_RESULT_KEYS
from postgwas.modules.formatting.table import (
    complete_rows,
    count_bounded_negative_log10_values,
    mapped_table,
    transformation_source,
    validation_expression,
    FormattingError,
)


def export_magma(
    frame,
    output_directory,
    dataset_id,
    config,
    *,
    overwrite,
    magma_config=None,
    ld_reference_prefix=None,
    analysis_scope=None,
    identifier_qc=None,
):
    """Write final MAGMA inputs, including pipeline reference selection."""
    schema = config.exports["magma"]
    result_keys = NAMED_OUTPUT_RESULT_KEYS["magma"]
    usable, excluded = complete_rows(
        frame,
        schema.validation.required_columns,
        validation_expression(schema.validation),
    )
    location_schema = schema.outputs["snp_location"]
    p_value_schema = schema.outputs["p_values"]
    p_source = transformation_source(
        p_value_schema, "negative_log10_to_raw_p",
    )
    p_bounded = count_bounded_negative_log10_values(
        usable, p_source, config.minimum_p_value,
    )
    preparation_qc = None
    excluded_variants = None
    if magma_config is not None:
        if ld_reference_prefix is None:
            raise FormattingError(
                "MAGMA pipeline export requires the resolved LD-reference prefix"
            )
        if analysis_scope is None:
            raise FormattingError(
                "MAGMA pipeline export requires the resolved analysis scope"
            )
        columns = magma_config.input
        policy = magma_config.snp_harmonisation
        variant = config.canonical_columns.resolved_variant_id
        chromosome = config.canonical_columns.chromosome
        position = config.canonical_columns.position
        bim = Path("%s%s" % (
            ld_reference_prefix, columns.bim_extension,
        )).expanduser().resolve()
        roles = list(columns.bim_columns)
        reference = (
            pl.scan_csv(
                bim,
                separator=columns.bim_delimiter,
                has_header=False,
                new_columns=roles,
                schema_overrides=[pl.String] * len(roles),
                low_memory=True,
            )
            .select(pl.col("variant_id").str.strip_chars().alias(variant))
            .collect(streaming=True)
        )
        if reference.is_empty() or reference[variant].is_null().any() or (
            reference[variant] == ""
        ).any():
            raise FormattingError(
                "MAGMA LD-reference BIM contains no usable variant identifiers: %s"
                % bim
            )
        input_unique = usable[variant].n_unique()
        input_identifiers = usable.select(variant).unique()
        repeated_reference_ids = (
            reference.filter(pl.col(variant).is_duplicated())
            .select(variant)
            .unique()
            .join(input_identifiers, on=variant, how="semi")
        )
        if not repeated_reference_ids.is_empty():
            examples = repeated_reference_ids[variant].head(5).to_list()
            raise FormattingError(
                "MAGMA LD-reference BIM assigns duplicate records to %d GWAS "
                "variant identifiers; examples: %s"
                % (repeated_reference_ids.height, examples)
            )
        matched_ids = (
            input_identifiers
            .join(reference.unique(), on=variant, how="semi")
        )
        matched = matched_ids.height
        overlap = matched / input_unique if input_unique else 0.0
        if matched == 0 or overlap < policy.minimum_overlap_fraction:
            raise FormattingError(
                "Only %d/%d unique GWAS variant identifiers (%.2f%%) occur in "
                "the MAGMA BIM reference; the configured minimum is %.2f%%."
                % (
                    matched,
                    input_unique,
                    overlap * 100,
                    policy.minimum_overlap_fraction * 100,
                )
            )
        selected = (
            usable.join(matched_ids, on=variant, how="semi")
            if policy.resolve_variants_to_reference else usable
        )
        scope = analysis_scope
        chromosome_mask = pl.col(chromosome).is_in(scope["exclude_chromosomes"])
        excluded_chromosome_counts = {
            label: selected.filter(pl.col(chromosome) == label).height
            for label in scope["exclude_chromosomes"]
        }
        region = scope["mhc_region"]
        mhc_mask = pl.lit(False)
        if scope.get("exclude_mhc_snps", False):
            if region is None:
                raise FormattingError("Resolved MAGMA MHC region is missing")
            mhc_mask = (
                (pl.col(chromosome) == region["chromosome"])
                & pl.col(position).is_between(
                    region["start"], region["end"], closed="both",
                )
            )
        reason = magma_config.exclusion_reporting.reason_column
        scoped = selected.with_columns(
            pl.when(chromosome_mask)
            .then(pl.lit(magma_config.exclusion_reporting.chromosome_reason))
            .when(mhc_mask)
            .then(pl.lit(magma_config.exclusion_reporting.mhc_reason))
            .otherwise(None)
            .alias(reason)
        )
        removed = scoped.filter(pl.col(reason).is_not_null())
        excluded_chromosome_rows = removed.filter(
            pl.col(reason) == magma_config.exclusion_reporting.chromosome_reason
        ).height
        excluded_mhc_rows = removed.filter(
            pl.col(reason) == magma_config.exclusion_reporting.mhc_reason
        ).height
        usable = scoped.filter(pl.col(reason).is_null()).drop(reason)
        if usable.is_empty():
            raise FormattingError(
                "MAGMA chromosome and MHC policies excluded every eligible variant"
            )
        excluded_path = configured_output_path(
            output_directory,
            magma_config.output_layout.excluded_variants,
            dataset_id=dataset_id,
            error_type=FormattingError,
        )
        excluded_variants = write_dataframe_table(
            removed.select(
                variant,
                chromosome,
                position,
                config.canonical_columns.reference_allele,
                config.canonical_columns.alternate_allele,
                reason,
            ),
            excluded_path,
            overwrite=overwrite,
            runtime=config.runtime,
            error_type=FormattingError,
        )
        identifier_qc = identifier_qc or {}
        preparation_qc = {
            "input_rows": int(frame.height - excluded),
            "input_unique_variants": int(input_unique),
            "reference_intersection_enabled": (
                policy.resolve_variants_to_reference
            ),
            "reference_variant_count": int(reference.height),
            "reference_id_match_rows": int(matched),
            "reference_unique_id_matches": int(matched),
            "not_in_reference_rows": int(input_unique - matched),
            "overlap_fraction": float(overlap),
            "duplicate_policy": policy.duplicate_policy,
            "duplicate_groups_detected": int(
                identifier_qc.get("identifier_duplicate_groups", 0)
            ),
            "duplicate_rows_detected": int(
                identifier_qc.get("identifier_duplicate_rows", 0)
            ),
            "duplicate_rows_removed": int(
                identifier_qc.get("identifier_duplicate_rows_excluded", 0)
            ),
            "duplicate_groups_resolved_by_lowest_p": int(
                identifier_qc.get(
                    "identifier_duplicate_groups_resolved_by_policy", 0,
                )
                if policy.duplicate_policy == "lowest_p" else 0
            ),
            "duplicate_groups_removed": int(
                identifier_qc.get("identifier_duplicate_groups", 0)
                if policy.duplicate_policy == "remove" else 0
            ),
            "mhc_policy": scope["mhc_policy"],
            "mhc_region": scope["mhc_region"],
            "mhc_region_source": scope["mhc_source"],
            "excluded_chromosomes": scope["exclude_chromosomes"],
            "excluded_chromosome_rows": excluded_chromosome_rows,
            "excluded_chromosome_counts": excluded_chromosome_counts,
            "excluded_mhc_rows": excluded_mhc_rows,
            "excluded_scope_rows": (
                excluded_chromosome_rows + excluded_mhc_rows
            ),
            "retained_rows": usable.height,
        }
    location_table = mapped_table(
        usable, location_schema, config.minimum_p_value,
    )
    location_header = location_schema.include_header
    if magma_config is not None:
        location_table = location_table.select(
            magma_config.input.variant_id_column,
            magma_config.input.chromosome_column,
            magma_config.input.position_column,
        )
        # MAGMA --annotate requires a headerless SNP, chromosome, position file.
        location_header = False
    p_value_table = mapped_table(
        usable, p_value_schema, config.minimum_p_value,
    )
    location = write_dataframe_table(
        location_table,
        configured_output_path(
            output_directory, location_schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        include_header=location_header,
        error_type=FormattingError,
    )
    p_values = write_dataframe_table(
        p_value_table,
        configured_output_path(
            output_directory, p_value_schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        include_header=p_value_schema.include_header,
        error_type=FormattingError,
    )
    result = {
        result_keys["snp_location"]: location,
        result_keys["p_values"]: p_values,
        "rows_in": frame.height,
        "rows_out": usable.height,
        "rows_excluded": excluded,
        "p_values_bounded": p_bounded,
        "columns": {
            "snp_location": location_table.columns,
            "p_values": p_value_table.columns,
        },
    }
    if preparation_qc is not None:
        result["variant_preparation"] = {
            "pval_file": p_values,
            "snp_loc_file": location,
            "excluded_variants": excluded_variants,
            "qc": preparation_qc,
        }
    return result


__all__ = ["export_magma"]
