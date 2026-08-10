from typing import Literal

from pydantic import Field

from postgwas.config.models.common import StrictModel


LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]


class LoggingConfig(StrictModel):
    console_level: LogLevel
    file_level: LogLevel
    show_progress: bool
    terminal_label_width: int = Field(ge=20, le=80)
    filename: str
    capture_external_tools: bool
    include_timestamps: bool
