from pathlib import Path

from pydantic import field_validator

from postgwas.config.models.common import StrictModel
from postgwas.core.paths import validate_filename_component


class RunConfig(StrictModel):
    output_directory: Path
    dataset_id: str
    overwrite: bool
    resume: bool

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str) -> str:
        return validate_filename_component(value, "run.dataset_id")
