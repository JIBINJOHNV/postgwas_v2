from typing import Literal

from pydantic import Field, model_validator

from postgwas.config.models.common import ModuleConfig, Population


class AlleleOrientationConfig(ModuleConfig):
    population: Population
    frequency_tolerance: float = Field(ge=0, le=1)
    palindromic_lower: float = Field(ge=0, le=0.5)
    palindromic_upper: float = Field(ge=0.5, le=1)
    unresolved_action: Literal["reject", "keep", "fail"]

    @model_validator(mode="after")
    def validate_palindromic_interval(self):
        if self.palindromic_upper <= self.palindromic_lower:
            raise ValueError("palindromic_upper must be greater than palindromic_lower")
        return self
