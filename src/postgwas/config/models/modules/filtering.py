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
from postgwas.core.variant_qc import MAXIMUM_SUPPORTED_INFO_SCORE
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
    start: int = Field(ge=1)


class FilteringMHCRegionOverrideConfig(StrictModel):
    chromosome: str | None = None
    start: int | None = Field(default=None, ge=1)
    end: int | None = Field(default=None, gt=0)

    @field_validator("chromosome")
    @classmethod
    def nonempty_chromosome(cls, value: str | None) -> str | None:
        if value is None:
            return None
        value = value.strip()
        if not value:
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def valid_selected_interval(self):
        if (
            self.start is not None
            and self.end is not None
            and self.end <= self.start
        ):
            raise ValueError("end must be greater than start")
        return self


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
    soft_filtered_vcf: str
    log_file: str
    preflight_log: str
    reason_summary: str
    summary_csv: str
    html_report: str
    mhc_exclusion_bed: str

    @field_validator("*")
    @classmethod
    def safe_output_pattern(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("output path patterns must not be empty")
        if "{dataset_id}" not in value:
            raise ValueError("output path patterns must contain {dataset_id}")
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
        build_specific = (
            self.filtered_vcf,
            self.soft_filtered_vcf,
            self.log_file,
            self.reason_summary,
            self.summary_csv,
            self.html_report,
            self.mhc_exclusion_bed,
        )
        if any("{genome_build}" not in pattern for pattern in build_specific):
            raise ValueError(
                "build-specific filtering outputs must contain {genome_build}"
            )
        if "{genome_build}" in self.preflight_log:
            raise ValueError(
                "preflight_log must not contain {genome_build} because it records "
                "failures that occur before the VCF build is inferred"
            )
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("filtering output path patterns must be unique")
        for name in ("filtered_vcf", "soft_filtered_vcf"):
            if not getattr(self, name).endswith(".vcf.gz"):
                raise ValueError("%s must end with .vcf.gz for tabix indexing" % name)
        if not self.reason_summary.endswith(".tsv"):
            raise ValueError("reason_summary must end with .tsv")
        if not self.summary_csv.endswith(".csv"):
            raise ValueError("summary_csv must end with .csv")
        if not self.html_report.endswith(".html"):
            raise ValueError("html_report must end with .html")
        return self


class FilteringReasonIdsConfig(StrictModel):
    empty_expression: str
    missing_pvalue: str
    pvalue_below_threshold: str
    missing_study_af: str
    maf_outside_range: str
    missing_imputation_quality: str
    imputation_quality_outside_range: str
    missing_study_info_af: str
    missing_external_af: str
    frequency_difference: str
    non_snp: str
    palindromic: str
    mhc: str

    @field_validator("*")
    @classmethod
    def valid_filter_id(cls, value: str) -> str:
        value = value.strip()
        if not VCF_TAG.fullmatch(value) or value == "PASS":
            raise ValueError("must be a valid non-PASS VCF FILTER ID")
        return value

    @model_validator(mode="after")
    def unique_filter_ids(self):
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("filter reason IDs must be unique")
        return self


class FilteringConfig(ModuleConfig):
    inputs: FilteringInputsConfig
    output_directory: Path | None = None
    maf_min: float | None = Field(default=None, ge=0, le=0.5, allow_inf_nan=False)
    missing_af_action: Literal["keep", "remove"]
    info_min: float | None = Field(default=None, ge=0, le=1, allow_inf_nan=False)
    info_max: float | None = Field(
        default=None,
        ge=0,
        le=MAXIMUM_SUPPORTED_INFO_SCORE,
        allow_inf_nan=False,
    )
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
    mhc_regions: dict[GenomeBuild, FilteringMHCConfig]
    mhc_region_override: FilteringMHCRegionOverrideConfig
    empty_expression_action: Literal["match_all", "match_none"]
    sort_output: bool
    write_soft_filter_vcf: bool
    report_missing_counts: bool
    filter_reason_ids: FilteringReasonIdsConfig
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

    @field_validator("mhc_regions")
    @classmethod
    def configured_mhc_regions(
        cls, values: dict[GenomeBuild, FilteringMHCConfig],
    ) -> dict[GenomeBuild, FilteringMHCConfig]:
        if not values:
            raise ValueError("must define at least one build-specific MHC region")
        return values

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
        override_values = self.mhc_region_override.model_dump().values()
        if any(value is not None for value in override_values) and not self.remove_mhc:
            raise ValueError(
                "mhc_region_override requires remove_mhc: true or --remove-mhc"
            )
        return self
