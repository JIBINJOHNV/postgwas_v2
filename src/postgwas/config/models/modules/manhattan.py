from pydantic import Field

from postgwas.config.models.common import ModuleConfig


class ManhattanConfig(ModuleConfig):
    significance_threshold: float = Field(gt=0, le=1)
    suggestive_threshold: float = Field(gt=0, le=1)
    width: float = Field(gt=0)
    height: float = Field(gt=0)
    file_format: str
