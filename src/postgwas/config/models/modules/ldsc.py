from pydantic import Field

from postgwas.config.models.common import ModuleConfig, Population


class LDSCConfig(ModuleConfig):
    population: Population
    intercept: float | None = None
    minimum_info: float = Field(ge=0, le=1)
    minimum_maf: float = Field(ge=0, le=0.5)
    chromosomes: int = Field(ge=1, le=25)
