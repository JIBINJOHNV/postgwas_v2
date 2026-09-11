from pathlib import Path

from pydantic import field_validator

from postgwas.config.models.common import StrictModel


class PipelineValidationConfig(StrictModel):
    """Presentation of mandatory pipeline validation; no scientific defaults."""

    report_file: Path

    @field_validator("report_file")
    @classmethod
    def relative_report_file(cls, value: Path) -> Path:
        if value.is_absolute() or value == Path(".") or ".." in value.parts:
            raise ValueError("pipeline.validation.report_file must be a safe relative file path")
        return value


class PipelineConfig(StrictModel):
    modules: list[str]
    fail_fast: bool
    stop_after: str | None = None
    validation: PipelineValidationConfig
