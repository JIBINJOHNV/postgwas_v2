"""Typed configuration for single-trait CBIIT LDSC heritability."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GenomeBuild,
    ModuleConfig,
    Population,
    StrictModel,
)


class LDSCReferenceLayout(StrictModel):
    """Filename suffixes required by CBIIT chromosome-split LD references."""

    ld_score_suffix: str
    m_suffix: str
    m_5_50_suffix: str

    @field_validator("*")
    @classmethod
    def safe_suffix(cls, value: str) -> str:
        if not value.startswith(".") or Path(value).name != value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value

    @model_validator(mode="after")
    def unique_suffixes(self):
        values = (self.ld_score_suffix, self.m_suffix, self.m_5_50_suffix)
        if len(values) != len(set(values)):
            raise ValueError("reference suffixes must be unique")
        return self


class LDSCOutputLayout(StrictModel):
    """All LDSC-owned output names and temporary locations."""

    output_prefix: str
    staging_directory: str
    resolved_config_file: str
    service_log_file: str
    munged_sumstats_suffix: str
    upstream_log_suffix: str
    observed_prefix_suffix: str
    liability_prefix_suffix: str
    covariance_suffix: str
    delete_values_suffix: str
    partitioned_delete_values_suffix: str

    @field_validator("output_prefix", "staging_directory", "service_log_file")
    @classmethod
    def dataset_pattern(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or "{dataset_id}" not in value:
            raise ValueError("must be a relative path containing {dataset_id}")
        return value

    @field_validator("resolved_config_file")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("must be a relative path below the output directory")
        return value

    @field_validator(
        "munged_sumstats_suffix", "upstream_log_suffix",
        "observed_prefix_suffix", "liability_prefix_suffix",
        "covariance_suffix", "delete_values_suffix",
        "partitioned_delete_values_suffix",
    )
    @classmethod
    def safe_suffix(cls, value: str) -> str:
        if not value or Path(value).name != value:
            raise ValueError("must be a non-empty filename suffix")
        return value

    @model_validator(mode="after")
    def unique_suffixes(self):
        analysis_outputs = [
            prefix + suffix
            for prefix in (
                self.observed_prefix_suffix, self.liability_prefix_suffix,
            )
            for suffix in (
                self.upstream_log_suffix, self.covariance_suffix,
                self.delete_values_suffix,
                self.partitioned_delete_values_suffix,
            )
        ]
        outputs = [
            self.munged_sumstats_suffix,
            self.upstream_log_suffix,
            *analysis_outputs,
        ]
        if len(outputs) != len(set(outputs)):
            raise ValueError("rendered LDSC output suffixes must be unique")
        return self


class LDSCConfig(ModuleConfig):
    """Supported upstream defaults for formatter-to-heritability execution.

    CBIIT LDSC reads chromosome-split references for the 22 autosomes. That
    chromosome count is an upstream protocol invariant, not a tunable PostGWAS
    default, so the schema deliberately accepts only 22.
    """

    population: Population | None
    genome_build: GenomeBuild | None
    chromosomes: Literal[22]
    intercept: float | None = Field(default=None, allow_inf_nan=False)
    minimum_info: float = Field(ge=0, le=1, allow_inf_nan=False)
    minimum_maf: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    minimum_n: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    chunksize: int = Field(ge=1)
    keep_maf: bool
    two_step: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    chisq_max: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    n_blocks: int = Field(gt=1)
    use_m_5_50: bool
    print_covariance: bool
    print_delete_values: bool
    sample_prevalence: float | None = Field(
        default=None, gt=0, lt=1, allow_inf_nan=False,
    )
    population_prevalence: float | None = Field(
        default=None, gt=0, lt=1, allow_inf_nan=False,
    )
    reference_layout: LDSCReferenceLayout
    output_layout: LDSCOutputLayout

    @model_validator(mode="after")
    def compatible_estimator_settings(self):
        if self.intercept is not None and self.two_step is not None:
            raise ValueError(
                "intercept and two_step cannot both be set for LDSC heritability"
            )
        return self
