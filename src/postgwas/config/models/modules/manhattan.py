from typing import Literal

from pydantic import Field

from postgwas.config.models.common import ModuleConfig


class ManhattanConfig(ModuleConfig):
    significance_threshold: float = Field(gt=0, le=1)
    suggestive_threshold: float = Field(gt=0, le=1)
    autosome_count: int = Field(gt=0)
    minimum_af: float = Field(ge=0, le=0.5)
    minimum_neglog10_p: float = Field(ge=0)
    loglog_pvalue: float = Field(gt=0)
    cytoband_ratio: float = Field(gt=0)
    maximum_height: float | None = Field(gt=0)
    chromosome_spacing: int = Field(ge=0)
    width: float = Field(gt=0)
    height: float | None = Field(gt=0)
    font_size: int = Field(gt=0)
    allelic_shift: bool
    flag_coding: bool
    file_format: Literal["pdf", "png"]
