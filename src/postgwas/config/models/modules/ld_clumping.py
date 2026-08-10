from pydantic import Field

from postgwas.config.models.common import GenomicRegion, ModuleConfig, Population


class LDClumpingConfig(ModuleConfig):
    population: Population
    lead_pvalue: float = Field(gt=0, le=1)
    clump_r2: float = Field(ge=0, le=1)
    lead_r2: float = Field(ge=0, le=1)
    window_kb: int = Field(gt=0)
    merge_distance_bp: int = Field(ge=0)
    remove_mhc: bool
    mhc: GenomicRegion
