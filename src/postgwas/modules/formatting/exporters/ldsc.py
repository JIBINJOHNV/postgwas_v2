"""LDSC munge-sumstats input export."""

import polars as pl
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.paths import configured_output_path

from postgwas.modules.formatting.table import (
    FormattingError,
    complete_rows,
    count_bounded_negative_log10_values,
    infer_study_design,
    mapped_table,
    transformation_source,
    validation_expression,
)


def export_ldsc(
    frame, output_directory, dataset_id, config, *, overwrite, study_design=None,
):
    """Write the sample-size schema understood by CBIIT ``munge_sumstats``."""
    schema = config.exports["ldsc"]
    design = study_design
    if design is None:
        design = infer_study_design(
            frame,
            config.study_design.case_count_column,
            config.study_design.control_count_column,
        )
    case_control = design.trait_type == "binary"
    sample_columns = schema.trait_required_columns[design.trait_type]
    required = [*schema.validation.required_columns, *sample_columns]
    usable, excluded = complete_rows(
        frame,
        required,
        validation_expression(schema.validation, positive_columns=sample_columns),
    )
    p_source = transformation_source(
        schema, "negative_log10_to_raw_p",
    )
    p_bounded = count_bounded_negative_log10_values(
        usable, p_source, config.minimum_p_value,
    )

    sample_prevalence = None
    sample_prevalence_aggregation = None
    sample_prevalence_variants = 0
    sample_prevalence_minimum = None
    sample_prevalence_maximum = None
    sample_prevalence_case_count_minimum = None
    sample_prevalence_case_count_maximum = None
    sample_prevalence_control_count_minimum = None
    sample_prevalence_control_count_maximum = None
    if case_control:
        case_column = config.study_design.case_count_column
        control_column = config.study_design.control_count_column
        prevalence = (
            pl.col(case_column)
            / (pl.col(case_column) + pl.col(control_column))
        )
        sample_prevalence_aggregation = (
            config.ldsc_sample_prevalence.aggregation
        )
        if sample_prevalence_aggregation == "median":
            aggregate = prevalence.median()
        elif sample_prevalence_aggregation == "mean":
            aggregate = prevalence.mean()
        else:  # Retain an actionable boundary error if validation is bypassed.
            raise FormattingError(
                "Unsupported LDSC sample-prevalence aggregation: %s"
                % sample_prevalence_aggregation
            )
        summary = usable.select(
            aggregate.alias("value"),
            prevalence.min().alias("minimum"),
            prevalence.max().alias("maximum"),
            pl.col(case_column).min().alias("case_count_minimum"),
            pl.col(case_column).max().alias("case_count_maximum"),
            pl.col(control_column).min().alias("control_count_minimum"),
            pl.col(control_column).max().alias("control_count_maximum"),
        ).to_dicts()[0]
        sample_prevalence = float(summary["value"])
        sample_prevalence_variants = usable.height
        sample_prevalence_minimum = float(summary["minimum"])
        sample_prevalence_maximum = float(summary["maximum"])
        sample_prevalence_case_count_minimum = float(
            summary["case_count_minimum"]
        )
        sample_prevalence_case_count_maximum = float(
            summary["case_count_maximum"]
        )
        sample_prevalence_control_count_minimum = float(
            summary["control_count_minimum"]
        )
        sample_prevalence_control_count_maximum = float(
            summary["control_count_maximum"]
        )

    mapping = dict(schema.columns)
    mapping.update(schema.trait_columns[design.trait_type])
    mapping.update(schema.trailing_columns)
    if not case_control:
        sample_size_mode = "quantitative_n_from_nco"
    else:
        sample_size_mode = "case_control_counts"
    output = mapped_table(
        usable, schema, config.minimum_p_value, columns=mapping,
    )
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
        "ldsc_file": path,
        "sample_prev": sample_prevalence,
        "sample_prevalence_aggregation": sample_prevalence_aggregation,
        "sample_prevalence_variants": sample_prevalence_variants,
        "sample_prevalence_minimum": sample_prevalence_minimum,
        "sample_prevalence_maximum": sample_prevalence_maximum,
        "sample_prevalence_case_count_minimum": (
            sample_prevalence_case_count_minimum
        ),
        "sample_prevalence_case_count_maximum": (
            sample_prevalence_case_count_maximum
        ),
        "sample_prevalence_control_count_minimum": (
            sample_prevalence_control_count_minimum
        ),
        "sample_prevalence_control_count_maximum": (
            sample_prevalence_control_count_maximum
        ),
        "trait_type": design.trait_type,
        "sample_size_mode": sample_size_mode,
        "columns": output.columns,
        "rows_in": frame.height,
        "rows_out": output.height,
        "rows_excluded": excluded,
        "p_values_bounded": p_bounded,
    }


__all__ = ["export_ldsc"]
