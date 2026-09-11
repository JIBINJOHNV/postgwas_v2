from pathlib import PurePath
from string import Formatter
from typing import Literal, get_args

from pydantic import Field, field_validator, model_validator

from postgwas.config.models.common import StrictModel


LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

# Presentation roles, not analysis decisions. Keep the role inventory identical
# to the shared screen renderer; the schema test checks the complete mapping.
TerminalRole = Literal[
    "genetic", "count", "analysis", "decision", "attention", "info", "loss",
    "success", "warning", "error", "run", "retry", "text",
]


class TerminalStyleConfig(StrictModel):
    """One run-wide terminal theme; defaults live only in application.yaml."""

    color: Literal["auto", "never"]
    styles: dict[TerminalRole, str]
    heading_indent: int = Field(ge=0, le=12)
    field_indent: int = Field(ge=0, le=20)
    indent_step: int = Field(ge=1, le=8)

    @field_validator("styles")
    @classmethod
    def complete_valid_styles(cls, value: dict[str, str]) -> dict[str, str]:
        from rich.style import Style
        from rich.errors import StyleSyntaxError

        missing = set(get_args(TerminalRole)) - set(value)
        if missing:
            raise ValueError("Missing terminal style roles: %s" % ", ".join(sorted(missing)))
        for role, style in value.items():
            try:
                Style.parse(style)
            except StyleSyntaxError as exc:
                raise ValueError("Invalid terminal style for %s: %s" % (role, style)) from exc
        return value


class FileValidationDisplayConfig(StrictModel):
    """Shared presentation only; these limits never control validation."""

    max_screen_files: int | None = Field(ge=1)
    max_checks_per_file: int = Field(ge=1)
    max_list_items: int = Field(ge=1)
    metric_labels: dict[str, str]
    role_labels: dict[str, str]
    outcome_metric_aliases: dict[str, str | list[str]]
    direct_report_file: str

    @model_validator(mode="after")
    def valid_metric_aliases(self):
        targets = [key for value in self.outcome_metric_aliases.values()
                   for key in ([value] if isinstance(value, str) else value)]
        unknown = set(targets) - set(self.metric_labels)
        if unknown or any(not label.strip() for label in self.outcome_metric_aliases):
            raise ValueError("outcome_metric_aliases must map non-empty labels to configured metric_labels keys")
        if any(not value for value in self.outcome_metric_aliases.values()):
            raise ValueError("outcome_metric_aliases must not contain empty target lists")
        return self

    @field_validator("direct_report_file")
    @classmethod
    def safe_report_pattern(cls, value: str) -> str:
        path = PurePath(value)
        if (not value.strip() or path.is_absolute() or ".." in path.parts
                or path.name in ("", ".")):
            raise ValueError("direct_report_file must be a safe relative file path")
        try:
            for _, name, specification, conversion in Formatter().parse(value):
                if name is not None and (name != "command" or specification or conversion):
                    raise ValueError("Only {command} is supported")
        except ValueError as exc:
            raise ValueError("direct_report_file permits only the {command} placeholder") from exc
        return value


class LoggingConfig(StrictModel):
    console_level: LogLevel
    file_level: LogLevel
    show_screen: bool
    show_progress: bool
    progress_refresh_seconds: float = Field(gt=0, le=60)
    terminal_label_width: int = Field(ge=20, le=80)
    terminal_style: TerminalStyleConfig
    filename: str
    screen_log_file: str
    capture_external_tools: bool
    include_timestamps: bool
    file_validation: FileValidationDisplayConfig

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
