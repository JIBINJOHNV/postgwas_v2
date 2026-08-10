"""Configuration primitives shared by GCTA-backed scientific modules."""

from __future__ import annotations

import re
from pathlib import Path

from pydantic import field_validator, model_validator

from postgwas.config.models.common import ModuleConfig, StrictModel


class GctaReferenceConfig(StrictModel):
    prefix: str | None = None
    population: str | None = None
    required_extensions: list[str]
    bim_columns: list[str]
    table_delimiter_pattern: str

    @field_validator("prefix", "population")
    @classmethod
    def optional_nonempty_text(cls, value: str | None) -> str | None:
        if value is not None and not value.strip():
            raise ValueError("must be null or non-empty")
        return value

    @field_validator("required_extensions", "bim_columns")
    @classmethod
    def unique_nonempty_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique values")
        if any(not str(value).strip() for value in values):
            raise ValueError("must not contain empty values")
        return values

    @field_validator("table_delimiter_pattern")
    @classmethod
    def valid_delimiter_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @model_validator(mode="after")
    def required_bim_roles(self):
        # ``--bfile`` is a GCTA/PLINK protocol invariant: these three files
        # share one prefix and their standard suffixes are not user policies.
        if set(self.required_extensions) != {".bed", ".bim", ".fam"}:
            raise ValueError(
                "required_extensions must contain .bed, .bim, and .fam exactly once"
            )
        required = {
            "chromosome", "variant_id", "genetic_distance", "position",
            "allele1", "allele2",
        }
        if set(self.bim_columns) != required:
            raise ValueError(
                "bim_columns must map the six PLINK BIM roles: %s"
                % ", ".join(sorted(required))
            )
        return self


class GctaBackedModuleConfig(ModuleConfig):
    """Executable-version contract shared by every GCTA-backed module."""

    minimum_gcta_version: str
    version_arguments: list[str]
    version_pattern: str
    version_probe_output_argument: str
    version_probe_output_name: str
    version_probe_temporary_prefix: str

    @field_validator("minimum_gcta_version")
    @classmethod
    def semantic_version(cls, value: str) -> str:
        if not re.fullmatch(r"\d+(?:\.\d+){1,2}", value):
            raise ValueError("must be a dotted numeric version")
        return value

    @field_validator("version_arguments")
    @classmethod
    def safe_version_arguments(cls, values: list[str]) -> list[str]:
        if not values or any(not value.strip() for value in values):
            raise ValueError("must contain one or more non-empty arguments")
        return values

    @field_validator("version_pattern")
    @classmethod
    def compilable_version_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        if compiled.groups < 1:
            raise ValueError("must capture the numeric version in group 1")
        return value

    @field_validator("version_probe_output_argument")
    @classmethod
    def safe_output_argument(cls, value: str) -> str:
        if not value.strip() or any(character.isspace() for character in value):
            raise ValueError("must be one non-empty command argument")
        return value

    @field_validator("version_probe_output_name")
    @classmethod
    def safe_output_name(cls, value: str) -> str:
        if not value.strip() or Path(value).name != value:
            raise ValueError("must be a non-empty filename")
        return value

    @field_validator("version_probe_temporary_prefix")
    @classmethod
    def safe_temporary_prefix(cls, value: str) -> str:
        if not value.strip() or "/" in value or "\\" in value:
            raise ValueError("must be a non-empty filename prefix")
        return value


__all__ = ["GctaBackedModuleConfig", "GctaReferenceConfig"]
