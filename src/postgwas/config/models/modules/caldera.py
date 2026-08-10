"""Typed configuration for CALDERA causal-gene prioritisation."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

from pydantic import field_validator, model_validator

from postgwas.config.models.common import GenomeBuild, ModuleConfig, StrictModel


class CalderaInputSchema(StrictModel):
    table_delimiter: str
    pops_gene_id_column: str
    pops_score_column: str
    credible_set_locus_column: str
    credible_set_chromosome_column: str
    credible_set_position_column: str
    credible_set_probability_column: str
    credible_set_coding_gene_column: str
    credible_set_rsid_column: str

    @field_validator("table_delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator(
        "pops_gene_id_column", "pops_score_column", "credible_set_locus_column",
        "credible_set_chromosome_column", "credible_set_position_column",
        "credible_set_probability_column", "credible_set_coding_gene_column",
        "credible_set_rsid_column",
    )
    @classmethod
    def nonempty_names(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value


class CalderaPipelineInput(StrictModel):
    index_file_names: list[str]
    index_filename_column: str
    index_locus_column: str
    credible_variant_column: str
    credible_probability_column: str
    variant_identifier_pattern: str
    locus_separator: str

    @field_validator("index_file_names")
    @classmethod
    def safe_index_names(cls, values: list[str]) -> list[str]:
        if (
            not values
            or len(values) != len(set(values))
            or any(Path(value).name != value or not value for value in values)
        ):
            raise ValueError("must contain unique plain filenames")
        return values

    @field_validator(
        "index_filename_column", "index_locus_column", "credible_variant_column",
        "credible_probability_column", "locus_separator",
    )
    @classmethod
    def nonempty_value(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("variant_identifier_pattern")
    @classmethod
    def valid_variant_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        required = {"chromosome", "position"}
        if not required.issubset(compiled.groupindex):
            raise ValueError("must define named chromosome and position groups")
        return value


class CalderaResultSchema(StrictModel):
    locus_column: str
    locus_position_column: str
    gene_name_column: str
    normalized_probability_column: str
    raw_probability_column: str
    locus_gene_count_column: str
    distance_column: str
    pops_score_column: str
    coding_probability_column: str
    imputed_pops_column: str
    gene_id_column: str

    @field_validator("*")
    @classmethod
    def nonempty_column(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def unique_columns(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("result columns must be unique")
        return self


class CalderaOutputLayout(StrictModel):
    results_file: str
    pipeline_credible_sets_file: str
    staging_directory: str
    resolved_config_file: str
    completion_manifest: str
    service_log_file: str

    @field_validator("*")
    @classmethod
    def safe_relative_pattern(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value

    @model_validator(mode="after")
    def dataset_paths(self):
        for field in ("results_file", "pipeline_credible_sets_file", "staging_directory", "service_log_file"):
            if "{dataset_id}" not in getattr(self, field):
                raise ValueError("%s must contain {dataset_id}" % field)
        return self


class CalderaConfig(ModuleConfig):
    genome_build: GenomeBuild
    repository_path: str | None = None
    adapter_script_path: str | None = None
    installed_repository_relative_path: str
    pops_file: str | None = None
    credible_set_file: str | None = None
    # The pinned upstream implementation fixes this boundary internally and
    # exposes no function argument with which PostGWAS could change it.
    minimum_credible_set_coverage: Literal[0.95]
    assembly_by_genome_build: dict[GenomeBuild, int]
    upstream_script_relative_path: str
    coding_variants_relative_path: str
    gene_locations_relative_pattern: str
    model_relative_path: str
    input_schema: CalderaInputSchema
    pipeline_input: CalderaPipelineInput
    result_schema: CalderaResultSchema
    output_layout: CalderaOutputLayout

    @field_validator(
        "repository_path", "adapter_script_path", "pops_file", "credible_set_file",
    )
    @classmethod
    def nonempty_path(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator(
        "upstream_script_relative_path", "coding_variants_relative_path",
        "gene_locations_relative_pattern", "model_relative_path",
        "installed_repository_relative_path",
    )
    @classmethod
    def safe_repository_path(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the CALDERA repository")
        return value

    @field_validator("gene_locations_relative_pattern")
    @classmethod
    def assembly_pattern(cls, value: str) -> str:
        if value.count("{assembly}") != 1:
            raise ValueError("must contain {assembly} exactly once")
        return value

    @model_validator(mode="after")
    def assembly_map_is_complete(self):
        required = set(GenomeBuild)
        if set(self.assembly_by_genome_build) != required:
            raise ValueError("assembly_by_genome_build must define GRCh37 and GRCh38")
        return self


__all__ = ["CalderaConfig"]
