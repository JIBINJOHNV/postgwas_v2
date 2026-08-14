"""Typed configuration for the upstream PoPS v0.2 integration."""

from __future__ import annotations

import re
from pathlib import Path
from string import Formatter
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GeneRankingReportingConfig,
    GenomeBuild,
    ModuleConfig,
    StrictModel,
)


PopsMethod = Literal["ridge", "lasso", "linreg"]
PopsGeneUniversePolicy = Literal["strict", "intersect"]


def _optional_nonempty(value: str | None) -> str | None:
    if value is not None and not value.strip():
        raise ValueError("must be null or non-empty")
    return value


class PopsInputSchema(StrictModel):
    table_delimiter_pattern: str
    target_table_delimiter: str
    gene_annotation_id_column: str
    gene_annotation_name_column: str
    gene_annotation_chromosome_column: str
    gene_annotation_tss_column: str
    magma_gene_id_column: str
    magma_score_column: str
    target_gene_id_column: str
    target_score_column: str
    feature_rows_suffix: str
    feature_columns_pattern: str
    feature_matrix_pattern: str
    magma_genes_out_suffix: str
    magma_genes_raw_suffix: str
    prediction_score_column: str
    prediction_training_column: str
    marginal_selected_column: str

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("target_table_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator(
        "gene_annotation_id_column", "gene_annotation_name_column",
        "gene_annotation_chromosome_column",
        "gene_annotation_tss_column", "magma_gene_id_column",
        "magma_score_column", "target_gene_id_column", "target_score_column",
        "prediction_score_column", "prediction_training_column",
        "marginal_selected_column",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator(
        "feature_rows_suffix", "magma_genes_out_suffix", "magma_genes_raw_suffix",
    )
    @classmethod
    def safe_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value

    @field_validator("feature_columns_pattern", "feature_matrix_pattern")
    @classmethod
    def valid_chunk_pattern(cls, value: str) -> str:
        try:
            fields = [field for _, field, _, _ in Formatter().parse(value) if field]
            rendered = value.format(chunk=0)
        except (KeyError, ValueError) as exc:
            raise ValueError("must be a valid pattern containing {chunk}") from exc
        if fields != ["chunk"] or not rendered.startswith(".") or Path(rendered).name != rendered:
            raise ValueError("must be a filename suffix pattern containing {chunk} once")
        return value


class PopsOutputLayout(StrictModel):
    output_prefix: str
    staging_directory: str
    resolved_config_file: str
    completion_manifest: str
    service_log_file: str
    predictions_suffix: str
    coefficients_suffix: str
    marginals_suffix: str
    upstream_log_suffix: str
    training_data_suffix: str
    matrix_data_suffix: str
    compatible_genes_out_suffix: str
    compatible_genes_raw_suffix: str
    excluded_genes_out_suffix: str
    excluded_genes_raw_suffix: str
    gene_compatibility_table_suffix: str
    gene_compatibility_report_suffix: str

    @field_validator("output_prefix", "staging_directory", "service_log_file")
    @classmethod
    def dataset_pattern(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or "{dataset_id}" not in value:
            raise ValueError("must be a relative path containing {dataset_id}")
        return value

    @field_validator("resolved_config_file", "completion_manifest")
    @classmethod
    def safe_resolved_config_path(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value

    @field_validator(
        "predictions_suffix", "coefficients_suffix", "marginals_suffix",
        "upstream_log_suffix", "training_data_suffix", "matrix_data_suffix",
        "compatible_genes_out_suffix", "compatible_genes_raw_suffix",
        "excluded_genes_out_suffix", "excluded_genes_raw_suffix",
        "gene_compatibility_table_suffix", "gene_compatibility_report_suffix",
    )
    @classmethod
    def safe_output_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value

    @model_validator(mode="after")
    def unique_suffixes(self):
        suffixes = [
            self.predictions_suffix, self.coefficients_suffix,
            self.marginals_suffix, self.upstream_log_suffix,
            self.training_data_suffix, self.matrix_data_suffix,
            self.compatible_genes_out_suffix, self.compatible_genes_raw_suffix,
            self.excluded_genes_out_suffix, self.excluded_genes_raw_suffix,
            self.gene_compatibility_table_suffix,
            self.gene_compatibility_report_suffix,
        ]
        if len(suffixes) != len(set(suffixes)):
            raise ValueError("output suffixes must be unique")
        return self


PopsReporting = GeneRankingReportingConfig


class PopsConfig(ModuleConfig):
    genome_build: GenomeBuild | None = None
    minimum_gene_count: int = Field(ge=1)
    gene_universe_policy: PopsGeneUniversePolicy
    magma_association_prefix: str | None = None
    feature_matrix_prefix: str | None = None
    feature_matrix_chunks: int = Field(ge=1)
    gene_location_file: str | None = None
    control_features_file: str | None = None
    target_score_file: str | None = None
    target_covariates_file: str | None = None
    target_error_covariance_file: str | None = None
    use_magma_covariates: bool
    use_magma_error_covariance: bool
    covariate_projection_chromosomes: list[str] | None = None
    remove_hla_during_covariate_projection: bool
    feature_subset_file: str | None = None
    feature_selection_chromosomes: list[str] | None = None
    feature_selection_p_cutoff: float = Field(gt=0, le=1, allow_inf_nan=False)
    maximum_selected_features: int | None = Field(default=None, ge=1)
    forward_selected_features: int | None = Field(default=None, ge=1)
    remove_hla_during_feature_selection: bool
    training_chromosomes: list[str] | None = None
    remove_hla_during_training: bool
    method: PopsMethod
    save_matrix_files: bool
    verbose: bool
    reporting: GeneRankingReportingConfig
    input_schema: PopsInputSchema
    output_layout: PopsOutputLayout

    @field_validator(
        "magma_association_prefix", "feature_matrix_prefix", "gene_location_file",
        "control_features_file", "target_score_file", "target_covariates_file",
        "target_error_covariance_file", "feature_subset_file",
    )
    @classmethod
    def optional_nonempty_path(cls, value: str | None) -> str | None:
        return _optional_nonempty(value)

    @field_validator(
        "covariate_projection_chromosomes", "feature_selection_chromosomes",
        "training_chromosomes",
    )
    @classmethod
    def valid_chromosome_list(cls, values: list[str] | None) -> list[str] | None:
        if values is None:
            return None
        normalized = [str(value).strip() for value in values]
        if not normalized or any(not value for value in normalized):
            raise ValueError("must be null or contain one or more chromosome labels")
        if len(normalized) != len(set(normalized)):
            raise ValueError("must not contain duplicate chromosome labels")
        return normalized

    @model_validator(mode="after")
    def target_companions_require_target_scores(self):
        if self.magma_association_prefix is not None and self.target_score_file is not None:
            raise ValueError(
                "magma_association_prefix and target_score_file are mutually exclusive"
            )
        if self.target_score_file is None and (
            self.target_covariates_file is not None
            or self.target_error_covariance_file is not None
        ):
            raise ValueError(
                "target covariates and covariance require target_score_file"
            )
        if (
            self.gene_universe_policy == "intersect"
            and self.target_score_file is not None
        ):
            raise ValueError(
                "gene_universe_policy=intersect is supported only for MAGMA "
                "targets; custom targets must already be aligned"
            )
        for label, output_suffix, raw_suffix in (
            (
                "compatible",
                self.output_layout.compatible_genes_out_suffix,
                self.output_layout.compatible_genes_raw_suffix,
            ),
            (
                "excluded",
                self.output_layout.excluded_genes_out_suffix,
                self.output_layout.excluded_genes_raw_suffix,
            ),
        ):
            genes_out_suffix = self.input_schema.magma_genes_out_suffix
            genes_raw_suffix = self.input_schema.magma_genes_raw_suffix
            if not output_suffix.endswith(genes_out_suffix):
                raise ValueError(
                    "%s derived MAGMA output suffix must end with %s"
                    % (label, genes_out_suffix)
                )
            prefix = output_suffix[:-len(genes_out_suffix)]
            if raw_suffix != prefix + genes_raw_suffix:
                raise ValueError(
                    "%s derived MAGMA output suffixes must share one prefix"
                    % label
                )
        return self


__all__ = [
    "PopsConfig",
    "PopsGeneUniversePolicy",
    "PopsMethod",
    "PopsReporting",
]
