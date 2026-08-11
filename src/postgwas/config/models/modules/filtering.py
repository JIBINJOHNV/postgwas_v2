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


class FilteringInputsConfig(StrictModel):
    vcf: Path | None = None
    dataset_id: str | None = None

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return validate_filename_component(
            value, "modules.filtering.inputs.dataset_id"
        )


class FilteringMHCConfig(GenomicRegion):
    genome_build: GenomeBuild
    start: int = Field(ge=1)


class FilteringVcfFieldsConfig(StrictModel):
    chromosome: str
    position: str
    reference_allele: str
    alternate_allele: str
    variant_type: str
    study_af_format: str
    imputation_quality_format: str
    log_pvalue_format: str
    study_af_info: str
    external_af_info: str

    @field_validator(
        "chromosome",
        "position",
        "reference_allele",
        "alternate_allele",
        "variant_type",
    )
    @classmethod
    def valid_fixed_field(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", value):
            raise ValueError("must be a bcftools fixed-field expression")
        return value

    @field_validator(
        "study_af_format",
        "imputation_quality_format",
        "log_pvalue_format",
    )
    @classmethod
    def valid_format_field(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"FORMAT/[A-Za-z][A-Za-z0-9_.-]*", value):
            raise ValueError("must be a valid bcftools FORMAT field expression")
        return value

    @field_validator("study_af_info")
    @classmethod
    def valid_info_field(cls, value: str) -> str:
        value = value.strip()
        if not re.fullmatch(r"INFO/[A-Za-z][A-Za-z0-9_.-]*", value):
            raise ValueError("must be a valid bcftools INFO field expression")
        return value

    @field_validator("external_af_info")
    @classmethod
    def valid_external_info_template(cls, value: str) -> str:
        value = value.strip()
        if value.count("{reference_population_tag}") != 1:
            raise ValueError(
                "must contain exactly one {reference_population_tag} placeholder"
            )
        try:
            rendered = value.format(reference_population_tag="TAG")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not re.fullmatch(r"INFO/[A-Za-z][A-Za-z0-9_.-]*", rendered):
            raise ValueError("must render a valid bcftools INFO field")
        return value


class FilteringOutputLayoutConfig(StrictModel):
    filtered_vcf: str
    log_file: str
    reason_summary: str
    mhc_exclusion_bed: str

    @field_validator("*")
    @classmethod
    def safe_output_pattern(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("output path patterns must not be empty")
        if "{dataset_id}" not in value or "{genome_build}" not in value:
            raise ValueError(
                "output path patterns must contain {dataset_id} and {genome_build}"
            )
        try:
            rendered = Path(
                value.format(dataset_id="dataset", genome_build="build")
            )
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if rendered.is_absolute() or ".." in rendered.parts:
            raise ValueError(
                "output path patterns must stay inside the output directory"
            )
        return value

    @model_validator(mode="after")
    def validate_destinations(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("filtering output path patterns must be unique")
        if not self.filtered_vcf.endswith(".vcf.gz"):
            raise ValueError("filtered_vcf must end with .vcf.gz for tabix indexing")
        return self


class FilteringConfig(ModuleConfig):
    inputs: FilteringInputsConfig
    output_directory: Path | None = None
    genome_build: GenomeBuild
    genome_build_header_tokens: dict[GenomeBuild, list[str]]
    maf_min: float | None = Field(default=None, ge=0, le=0.5, allow_inf_nan=False)
    missing_af_action: Literal["keep", "remove"]
    info_min: float | None = Field(default=None, ge=0, le=1, allow_inf_nan=False)
    info_max: float | None = Field(default=None, ge=0, le=1, allow_inf_nan=False)
    missing_info_action: Literal["keep", "remove"]
    minimum_neglog10_p: float | None = Field(
        default=None, ge=0, allow_inf_nan=False,
    )
    missing_pvalue_action: Literal["keep", "remove"]
    reference_population_tag: str
    frequency_difference_max: float | None = Field(
        default=None, ge=0, le=1, allow_inf_nan=False,
    )
    include_indels: bool
    remove_palindromic: bool
    palindromic_lower: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    palindromic_upper: float = Field(ge=0.5, le=1, allow_inf_nan=False)
    remove_mhc: bool
    mhc: FilteringMHCConfig
    empty_expression_action: Literal["match_all", "match_none"]
    sort_output: bool
    report_missing_counts: bool
    vcf_fields: FilteringVcfFieldsConfig
    output_layout: FilteringOutputLayoutConfig

    @field_validator("reference_population_tag")
    @classmethod
    def valid_reference_population_tag(cls, value: str) -> str:
        value = value.strip()
        if not VCF_TAG.fullmatch(value):
            raise ValueError(
                "must be a valid VCF tag: letters followed by letters, digits, "
                "dot, underscore, or hyphen"
            )
        return value

    @field_validator("genome_build_header_tokens")
    @classmethod
    def valid_build_header_tokens(
        cls, values: dict[GenomeBuild, list[str]],
    ) -> dict[GenomeBuild, list[str]]:
        if not values:
            raise ValueError("must define at least one genome-build header token")
        normalized: dict[GenomeBuild, list[str]] = {}
        for build, tokens in values.items():
            cleaned = [str(token).strip() for token in tokens]
            if not cleaned or any(not token for token in cleaned):
                raise ValueError(
                    "every genome build must define non-empty header tokens"
                )
            if len(cleaned) != len(set(token.lower() for token in cleaned)):
                raise ValueError("genome-build header tokens must be unique per build")
            normalized[build] = cleaned
        return normalized

    @model_validator(mode="after")
    def validate_intervals(self):
        if (
            self.info_min is not None
            and self.info_max is not None
            and self.info_max < self.info_min
        ):
            raise ValueError("info_max must be greater than or equal to info_min")
        if self.palindromic_upper <= self.palindromic_lower:
            raise ValueError("palindromic_upper must be greater than palindromic_lower")
        if self.mhc.genome_build != self.genome_build:
            raise ValueError("mhc.genome_build must match genome_build")
        if self.genome_build not in self.genome_build_header_tokens:
            raise ValueError(
                "genome_build_header_tokens must define the configured genome_build"
            )
        return self
