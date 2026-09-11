from typing import Literal

from pydantic import ConfigDict, Field, field_validator, model_validator

from postgwas.config.models.common import ModuleConfig


class ManhattanConfig(ModuleConfig):
    model_config = ConfigDict(allow_inf_nan=False)
    significance_threshold: float = Field(gt=0, le=1)
    suggestive_threshold: float = Field(gt=0, le=1)
    autosome_count: int = Field(gt=0)
    minimum_af: float = Field(ge=0, le=0.5)
    minimum_neglog10_p: float = Field(ge=0)
    loglog_pvalue: float = Field(gt=1)
    cytoband_ratio: float = Field(gt=0)
    maximum_height: float | None = Field(gt=0)
    chromosome_spacing: int = Field(ge=0)
    width: float = Field(gt=0)
    height: float | None = Field(gt=0)
    font_size: int = Field(gt=0)
    allelic_shift: bool
    flag_coding: bool
    file_format: Literal["pdf", "png"]
    output_file: str = Field(min_length=1)
    log_file: str = Field(min_length=1)
    plot_data_file: str = Field(min_length=1)
    single_chromosome_height: float = Field(gt=0)
    genomewide_height: float = Field(gt=0)
    png_dpi: int = Field(gt=0)
    pvalue_field: str
    frequency_field: str
    allelic_shift_field: str
    consequence_field: str
    consequence_info_field: str
    coding_terms: list[str] = Field(min_length=1)
    chromosome_label_policy: Literal["exact", "strip_prefix"]
    chromosome_prefix: str = Field(min_length=1)
    axis_label_rows: int = Field(ge=1)
    caption_font_size: float = Field(gt=0)
    caption_template: str = Field(min_length=1)

    @field_validator("pvalue_field", "frequency_field", "allelic_shift_field", "consequence_field", "consequence_info_field")
    @classmethod
    def valid_vcf_field(cls, value):
        from postgwas.core.vcf import VCF_TAG

        if not VCF_TAG.fullmatch(value):
            raise ValueError("must be a valid VCF field name")
        return value

    @field_validator("coding_terms")
    @classmethod
    def valid_coding_terms(cls, values):
        if len(values) != len(set(values)) or any(not value.strip() or "," in value for value in values):
            raise ValueError("coding_terms must contain unique nonempty terms without commas")
        return values

    @model_validator(mode="after")
    def valid_plot_ranges(self):
        if self.allelic_shift and self.minimum_af > 0:
            raise ValueError("allelic_shift cannot be combined with minimum_af > 0")
        if self.maximum_height is not None and self.maximum_height < self.minimum_neglog10_p:
            raise ValueError("maximum_height must be at least minimum_neglog10_p")
        return self
