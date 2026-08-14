"""Schema-validated configuration for genotype-free GWAS-VCF QC."""

import re
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import (
    GenomeBuild,
    GenomicRegion,
    ModuleConfig,
    StrictModel,
)
from postgwas.core.paths import validate_filename_component
from postgwas.core.vcf import VCF_TAG


class QCSummaryInputsConfig(StrictModel):
    vcf: Path | None = None
    dataset_id: str | None = None

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return validate_filename_component(
            value, "modules.qc_summary.inputs.dataset_id"
        )


class QCSummaryMHCRegionConfig(GenomicRegion):
    start: int = Field(ge=1)


class QCSummaryRulesConfig(StrictModel):
    minimum_neglog10_p: float | None = Field(
        default=None, ge=0, allow_inf_nan=False,
    )
    missing_pvalue_action: Literal["keep", "remove"]
    maf_min: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    missing_af_action: Literal["keep", "remove"]
    info_min: float = Field(ge=0, le=1, allow_inf_nan=False)
    info_max: float = Field(ge=0, le=100, allow_inf_nan=False)
    missing_info_action: Literal["keep", "remove"]
    maximum_af_difference: float = Field(ge=0, le=1, allow_inf_nan=False)
    include_indels: bool
    remove_palindromic: bool
    palindromic_lower: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    palindromic_upper: float = Field(ge=0.5, le=1, allow_inf_nan=False)
    remove_mhc: bool
    mhc_regions: dict[GenomeBuild, QCSummaryMHCRegionConfig]
    sample_size_outlier_standard_deviations: float = Field(
        gt=0, allow_inf_nan=False,
    )

    @model_validator(mode="after")
    def validate_intervals(self):
        if self.info_max < self.info_min:
            raise ValueError("info_max must be greater than or equal to info_min")
        if self.palindromic_upper <= self.palindromic_lower:
            raise ValueError(
                "palindromic_upper must be greater than palindromic_lower"
            )
        return self


class QCSummaryVcfFieldsConfig(StrictModel):
    chromosome: str
    position: str
    reference_allele: str
    alternate_allele: str
    study_info_af: str
    external_info_af: str
    study_format_af: str
    imputation_format: str
    log_pvalue_format: str
    effective_sample_size_format: str

    @field_validator(
        "chromosome", "position", "reference_allele", "alternate_allele",
    )
    @classmethod
    def fixed_field_query(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"%[A-Z][A-Z0-9_]*", value):
            raise ValueError("must be a bcftools fixed-field query")
        return value

    @field_validator("study_info_af")
    @classmethod
    def info_field_query(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", value):
            raise ValueError("must be a bcftools INFO-field query")
        return value

    @field_validator(
        "study_format_af", "imputation_format", "log_pvalue_format",
        "effective_sample_size_format",
    )
    @classmethod
    def format_field_query(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"\[%[A-Za-z][A-Za-z0-9_.-]*\]", value):
            raise ValueError("must be a bracketed bcftools FORMAT-field query")
        return value

    @field_validator("external_info_af")
    @classmethod
    def external_af_template(cls, value: str) -> str:
        value = value.strip()
        if value.count("{external_af}") != 1:
            raise ValueError(
                "must contain exactly one {external_af} placeholder"
            )
        try:
            rendered = value.format(external_af="TAG")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not re.fullmatch(r"%INFO/[A-Za-z][A-Za-z0-9_.-]*", rendered):
            raise ValueError("must render a valid bcftools INFO-field query")
        return value


class QCSummaryTableConfig(StrictModel):
    delimiter: str
    null_values: list[str]
    null_output: str
    temporary_suffix: str
    io_buffer_bytes: int = Field(gt=0)

    @field_validator("delimiter")
    @classmethod
    def one_character_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("must contain exactly one character")
        return value

    @field_validator("null_values")
    @classmethod
    def unique_null_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique values")
        return values

    @field_validator("temporary_suffix")
    @classmethod
    def valid_temporary_suffix(cls, value: str) -> str:
        if not value.startswith(".") or "/" in value:
            raise ValueError("must be a filename suffix beginning with '.'")
        return value


class QCSummaryOutputLayoutConfig(StrictModel):
    metric_report: str
    rule_report: str
    assessment_json: str
    temporary_table_prefix: str
    resolved_config_file: str
    service_log_file: str
    completion_manifest: str

    @field_validator("*")
    @classmethod
    def safe_output_pattern(cls, value: str) -> str:
        value = value.strip()
        if not value or "{dataset_id}" not in value:
            raise ValueError("must be non-empty and contain {dataset_id}")
        try:
            rendered = Path(value.format(dataset_id="dataset", build="GRCh37"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError("must stay inside the configured output directory")
        return value

    @model_validator(mode="after")
    def unique_build_specific_outputs(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("QC output path patterns must be unique")
        for field in type(self).model_fields:
            if "{build}" not in getattr(self, field):
                raise ValueError("%s must contain {build}" % field)
        return self


class QCSummaryConfig(ModuleConfig):
    inputs: QCSummaryInputsConfig
    output_directory: Path | None = None
    target_build: GenomeBuild
    reference_af_column: str
    rules: QCSummaryRulesConfig
    vcf_fields: QCSummaryVcfFieldsConfig
    table: QCSummaryTableConfig
    output_layout: QCSummaryOutputLayoutConfig

    @field_validator("reference_af_column")
    @classmethod
    def valid_reference_af_column(cls, value: str) -> str:
        value = value.strip()
        if not VCF_TAG.fullmatch(value):
            raise ValueError(
                "must be a valid VCF INFO tag beginning with a letter and then "
                "containing only letters, digits, dot, underscore, or hyphen"
            )
        return value

    def mhc_region(self, genome_build: str | GenomeBuild) -> QCSummaryMHCRegionConfig:
        try:
            return self.rules.mhc_regions[GenomeBuild(genome_build)]
        except (KeyError, ValueError) as exc:
            raise ValueError(
                "No MHC region is configured for genome build %s" % genome_build
            ) from exc
