from pathlib import PurePath
from typing import Literal

from pydantic import Field, field_validator

from postgwas.config.models.common import StrictModel


LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]


class LoggingConfig(StrictModel):
    console_level: LogLevel
    file_level: LogLevel
    show_screen: bool
    show_progress: bool
    progress_refresh_seconds: float = Field(gt=0, le=60)
    terminal_label_width: int = Field(ge=20, le=80)
    filename: str
    screen_log_file: str
    capture_external_tools: bool
    include_timestamps: bool

    @field_validator("screen_log_file")
    @classmethod
    def safe_relative_log_path(cls, value: str) -> str:
        text = str(value).strip()
        path = PurePath(text)
        if (
            not text
            or path.is_absolute()
            or ".." in path.parts
            or path.name in ("", ".")
        ):
            raise ValueError("must be a non-empty relative file path without '..'")
        return text
