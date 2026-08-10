from pathlib import Path

from pydantic import Field

from postgwas.config.models.common import StrictModel


class ExecutionConfig(StrictModel):
    threads: int = Field(ge=1)
    memory_gb: float = Field(gt=0)
    reserve_threads: int = Field(ge=0)
    usable_memory_fraction: float = Field(gt=0, le=1)
    random_seed: int = Field(ge=0)
    temporary_directory: Path | None = None
    retries: int = Field(ge=0)
    timeout_seconds: int | None = Field(default=None, gt=0)
    fail_fast: bool = True
