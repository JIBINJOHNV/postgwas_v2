from pathlib import Path
from string import Formatter
from typing import Literal

from pydantic import field_validator

from postgwas.config.models.common import StrictModel
from postgwas.core.paths import validate_filename_component


class ResumePolicyConfig(StrictModel):
    """Global decisions for checksum-validated checkpoint recovery."""

    checkpoint_validation: Literal["sha256"]
    checkpoint_directory: Path
    direct_manifest: str
    pipeline_stage_manifest: str
    audit_log: str
    partial_results: Literal["resume_validated_stages"]
    changed_parameters: Literal["warn_and_restart", "error"]
    changed_inputs: Literal["warn_and_restart", "error"]
    unvalidated_outputs: Literal["warn_and_restart", "error"]

    @field_validator("checkpoint_directory")
    @classmethod
    def relative_checkpoint_directory(cls, value: Path) -> Path:
        if (
            value.is_absolute()
            or value in (Path("."), Path(".."))
            or ".." in value.parts
        ):
            raise ValueError(
                "run.resume_policy.checkpoint_directory must be a safe relative directory"
            )
        return value

    @staticmethod
    def _validate_pattern(
        value: str,
        *,
        field_name: str,
        required: set[str],
    ) -> str:
        pattern = str(value).strip()
        path = Path(pattern)
        if not pattern or path.is_absolute() or ".." in path.parts:
            raise ValueError(
                "run.resume_policy.%s must be a safe relative path pattern"
                % field_name
            )
        try:
            fields = {
                name for _, name, _, _ in Formatter().parse(pattern) if name
            }
        except ValueError as exc:
            raise ValueError(
                "run.resume_policy.%s is not a valid path pattern" % field_name
            ) from exc
        if fields != required:
            raise ValueError(
                "run.resume_policy.%s must contain exactly: %s"
                % (
                    field_name,
                    ", ".join("{%s}" % name for name in sorted(required)),
                )
            )
        return pattern

    @field_validator("direct_manifest")
    @classmethod
    def valid_direct_manifest(cls, value: str) -> str:
        return cls._validate_pattern(
            value, field_name="direct_manifest", required={"command"},
        )

    @field_validator("pipeline_stage_manifest")
    @classmethod
    def valid_pipeline_stage_manifest(cls, value: str) -> str:
        return cls._validate_pattern(
            value,
            field_name="pipeline_stage_manifest",
            required={"module", "stage_number"},
        )

    @field_validator("audit_log")
    @classmethod
    def relative_audit_log(cls, value: str) -> str:
        return cls._validate_pattern(
            value, field_name="audit_log", required=set(),
        )


class RunConfig(StrictModel):
    output_directory: Path
    dataset_id: str
    overwrite: bool
    resume: bool
    resume_policy: ResumePolicyConfig

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str) -> str:
        return validate_filename_component(value, "run.dataset_id")
