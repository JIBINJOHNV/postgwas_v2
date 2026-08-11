"""MiXeR univariate summary-statistics export."""

import polars as pl
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path

from postgwas.modules.formatting.table import (
    FormattingError,
    complete_rows,
    infer_study_design,
    mapped_table,
    validation_expression,
)


def export_mixer(
    frame, output_directory, dataset_id, config, *, overwrite, study_design=None,
):
    """Write the official MiXeR ``SNP CHR BP A1 A2 N Z`` schema."""
    settings = config.mixer
    schema = config.exports["mixer"]
    design = study_design
    if design is None:
        design = infer_study_design(
            frame,
            config.study_design.case_count_column,
            config.study_design.control_count_column,
        )
    required = list(schema.validation.required_columns)
    if settings.minimum_info is not None:
        required.append(settings.info_column)

    condition = (
        validation_expression(schema.validation)
        & pl.col(settings.chromosome_column).is_in(config.chromosomes)
        & (
            pl.col(settings.effect_allele_column)
            != pl.col(settings.other_allele_column)
        )
    )
    if settings.minimum_info is not None:
        condition &= pl.col(settings.info_column) >= settings.minimum_info
    if settings.snps_only:
        condition &= (
            (pl.col(settings.effect_allele_column).str.len_chars() == 1)
            & (pl.col(settings.other_allele_column).str.len_chars() == 1)
            & pl.col(settings.effect_allele_column).is_in(settings.allowed_alleles)
            & pl.col(settings.other_allele_column).is_in(settings.allowed_alleles)
        )

    # MiXeR's recommendation is relative to the study-wide median N, not the
    # median after INFO or variant-type exclusions have already been applied.
    sample_sizes, _ = complete_rows(
        frame,
        [settings.chromosome_column, settings.sample_size_column],
        (pl.col(settings.sample_size_column) > 0)
        & pl.col(settings.chromosome_column).is_in(config.chromosomes),
    )
    median_n = sample_sizes.select(pl.col(settings.sample_size_column).median()).item()
    if median_n is None or median_n <= 0:
        raise FormattingError("MiXeR effective sample size has no positive median.")
    minimum_n = float(median_n) * settings.minimum_sample_size_fraction

    usable, _ = complete_rows(frame, required, condition)
    usable = usable.filter(pl.col(settings.sample_size_column) >= minimum_n)
    if usable.is_empty():
        raise FormattingError(
            "No variants remain after the configured MiXeR sample-size check."
        )

    output = mapped_table(usable, schema, config.minimum_p_value)
    path = write_dataframe_table(
        output,
        configured_output_path(
            output_directory, schema.output_file, dataset_id=dataset_id,
        ),
        overwrite=overwrite,
        runtime=config.runtime,
        error_type=FormattingError,
    )
    return {
        "mixer_input": path,
        "trait_type": design.trait_type,
        "sample_size_mode": (
            "effective_n_from_case_control_counts"
            if design.trait_type == "binary"
            else "total_n_from_nco"
        ),
        "median_sample_size": float(median_n),
        "minimum_sample_size": minimum_n,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": frame.height - output.height,
        "columns": output.columns,
    }


__all__ = ["export_mixer"]
