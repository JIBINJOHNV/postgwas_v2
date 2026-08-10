"""Typed configuration for the upstream K-POPS v1.0.0 integration."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GeneRankingReportingConfig,
    GenomeBuild,
    ModuleConfig,
    StrictModel,
)


class KPopsInputSchema(StrictModel):
    table_delimiter_pattern: str
    predictions_table_delimiter: str
    published_table_delimiter: str
    gene_id_column: str
    gene_name_column: str
    chromosome_column: str
    tss_column: str
    magma_gene_id_column: str
    magma_score_column: str
    predictions_gene_id_column: str
    predictions_score_column: str
    published_gene_id_column: str
    predictions_gene_id_from_index: bool
    kernel_genes_suffix: str
    kernel_matrix_suffix: str
    magma_genes_out_suffix: str
    magma_genes_raw_suffix: str
    kernel_float_bytes: int = Field(ge=1)

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("predictions_table_delimiter", "published_table_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator(
        "gene_id_column", "gene_name_column", "chromosome_column", "tss_column",
        "magma_gene_id_column", "magma_score_column",
        "predictions_gene_id_column", "predictions_score_column",
        "published_gene_id_column",
    )
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator(
        "kernel_genes_suffix", "kernel_matrix_suffix", "magma_genes_out_suffix",
        "magma_genes_raw_suffix",
    )
    @classmethod
    def safe_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value


class KPopsOutputLayout(StrictModel):
    output_prefix: str
    staging_directory: str
    resolved_config_file: str
    completion_manifest: str
    service_log_file: str
    predictions_suffix: str
    coefficients_suffix: str
    attribution_suffix: str
    attribution_rows_suffix: str
    attribution_columns_suffix: str

    @field_validator("output_prefix", "staging_directory", "service_log_file")
    @classmethod
    def dataset_pattern(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or "{dataset_id}" not in value:
            raise ValueError("must be a relative path containing {dataset_id}")
        return value

    @field_validator("resolved_config_file", "completion_manifest")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value

    @field_validator(
        "predictions_suffix", "coefficients_suffix", "attribution_suffix",
        "attribution_rows_suffix", "attribution_columns_suffix",
    )
    @classmethod
    def safe_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value


class KPopsConfig(ModuleConfig):
    genome_build: GenomeBuild | None = None
    script_path: str
    gene_annotation_file: str | None = None
    kernel_matrix_prefix: str | None = None
    magma_association_prefix: str | None = None
    minimum_gene_count: int = Field(ge=2)
    kernel_validation_chunk_rows: int = Field(ge=1)
    kernel_symmetry_relative_tolerance: float = Field(ge=0, allow_inf_nan=False)
    kernel_symmetry_absolute_tolerance: float = Field(ge=0, allow_inf_nan=False)
    use_magma_covariates: bool
    covariate_projection_chromosomes: list[str] | None = None
    remove_hla_during_covariate_projection: bool
    training_chromosomes: list[str]
    remove_hla_during_training: bool
    remove_hla_during_testing: bool
    device: Literal["cpu", "cuda", "mps"]
    top_contributor_gene_count: int = Field(ge=1)
    anchor_genes: list[str]
    anchor_gene_type: Literal["ENSGID", "NAME"]
    save_attribution_files: bool
    verbose: bool
    reporting: GeneRankingReportingConfig
    input_schema: KPopsInputSchema
    output_layout: KPopsOutputLayout

    @field_validator(
        "script_path", "gene_annotation_file", "kernel_matrix_prefix",
        "magma_association_prefix",
    )
    @classmethod
    def nonempty_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator("covariate_projection_chromosomes", "training_chromosomes")
    @classmethod
    def valid_chromosomes(cls, values: list[str] | None) -> list[str] | None:
        if values is None:
            return None
        normalized = [str(value).strip() for value in values]
        if not normalized or any(not value for value in normalized):
            raise ValueError("must contain one or more chromosome labels")
        if len(normalized) != len(set(normalized)):
            raise ValueError("must not contain duplicate chromosome labels")
        return normalized

    @field_validator("anchor_genes")
    @classmethod
    def unique_anchor_genes(cls, values: list[str]) -> list[str]:
        normalized = [value.strip() for value in values]
        if any(not value for value in normalized) or len(normalized) != len(set(normalized)):
            raise ValueError("must contain unique non-empty gene identifiers")
        return normalized

    @model_validator(mode="after")
    def valid_training_mode(self):
        special = {"loco", "all"}.intersection(self.training_chromosomes)
        if special and len(self.training_chromosomes) != 1:
            raise ValueError("training_chromosomes 'loco' or 'all' must be used alone")
        return self


__all__ = ["KPopsConfig"]
