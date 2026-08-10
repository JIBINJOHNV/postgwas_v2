from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from postgwas.config.models.common import GenomicRegion, ModuleConfig, StrictModel


class FilteringInputsConfig(StrictModel):
    vcf: Path | None = None
    dataset_id: str | None = None


class FilteringConfig(ModuleConfig):
    inputs: FilteringInputsConfig
    output_directory: Path | None = None
    maf_min: float = Field(ge=0, le=0.5)
    info_min: float = Field(ge=0, le=1)
    info_max: float | None = Field(default=None, ge=0, le=1)
    missing_info_action: Literal["keep", "remove"]
    minimum_neglog10_p: float | None = Field(default=None, ge=0)
    reference_population_tag: str
    frequency_difference_max: float = Field(ge=0, le=1)
    include_indels: bool
    remove_palindromic: bool
    palindromic_lower: float = Field(ge=0, le=0.5)
    palindromic_upper: float = Field(ge=0.5, le=1)
    remove_mhc: bool
    mhc: GenomicRegion

    @model_validator(mode="after")
    def validate_intervals(self):
        if self.info_max is not None and self.info_max < self.info_min:
            raise ValueError("info_max must be greater than or equal to info_min")
        if self.palindromic_upper <= self.palindromic_lower:
            raise ValueError("palindromic_upper must be greater than palindromic_lower")
        return self
