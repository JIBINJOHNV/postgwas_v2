from pathlib import Path
from string import Formatter
from typing import Literal

from pydantic import Field, field_validator

from postgwas.config.models.common import (
    GenomeBuild,
    ModuleConfig,
    Population,
    StrictModel,
)


class PredLDConfig(StrictModel):
    minimum_r2: float = Field(ge=0, le=1)
    minimum_maf: float = Field(ge=0, le=0.5)
    mode: Literal["TOP_LD"]
    correlation_method: Literal["pearson", "spearman"]
    imputed_dataset_suffix: str
    memory_gb_per_worker: float = Field(gt=0, allow_inf_nan=False)
    free_memory_threshold_gb: float = Field(gt=0, allow_inf_nan=False)
    free_memory_threshold_fraction: float = Field(
        gt=0, le=1, allow_inf_nan=False,
    )
    memory_poll_seconds: float = Field(gt=0, allow_inf_nan=False)
    worker_poll_seconds: float = Field(gt=0, allow_inf_nan=False)
    large_chromosomes: list[str]
    preferred_chromosome_order: list[str] = Field(min_length=1)
    reference_subdirectory_template: str
    reference_file_template: str
    reference_file_kinds: list[str] = Field(min_length=1)

    @field_validator("imputed_dataset_suffix")
    @classmethod
    def safe_dataset_suffix(cls, value: str) -> str:
        permitted = (
            "_-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        )
        if not value or any(character not in permitted for character in value):
            raise ValueError(
                "must contain only letters, digits, underscores, or hyphens"
            )
        return value

    @field_validator("large_chromosomes", "preferred_chromosome_order")
    @classmethod
    def unique_normalised_chromosomes(cls, value: list[str]) -> list[str]:
        normalised = [
            str(chromosome).upper().removeprefix("CHR")
            for chromosome in value
        ]
        if any(not chromosome for chromosome in normalised):
            raise ValueError("must not contain an empty chromosome label")
        if len(set(normalised)) != len(normalised):
            raise ValueError("must contain unique chromosome labels")
        return normalised

    @field_validator("reference_subdirectory_template")
    @classmethod
    def safe_reference_subdirectory_template(cls, value: str) -> str:
        cls._validate_template(
            value,
            allowed={"mode", "population"},
            required={"mode", "population"},
            label="reference_subdirectory_template",
            allow_directories=True,
        )
        return value

    @field_validator("reference_file_template")
    @classmethod
    def safe_reference_file_template(cls, value: str) -> str:
        cls._validate_template(
            value,
            allowed={"population", "chromosome", "kind"},
            required={"population", "chromosome", "kind"},
            label="reference_file_template",
            allow_directories=False,
        )
        return value

    @field_validator("reference_file_kinds")
    @classmethod
    def safe_reference_file_kinds(cls, value: list[str]) -> list[str]:
        if len(set(value)) != len(value):
            raise ValueError("must contain unique reference file kinds")
        if any(
            not kind
            or Path(kind).name != kind
            or any(character.isspace() for character in kind)
            for kind in value
        ):
            raise ValueError("must contain safe non-empty filename components")
        return value

    @staticmethod
    def _validate_template(
        value: str,
        *,
        allowed: set[str],
        required: set[str],
        label: str,
        allow_directories: bool,
    ) -> None:
        fields = set()
        try:
            parsed = list(Formatter().parse(value))
        except ValueError as exc:
            raise ValueError("contains invalid format syntax") from exc
        for _literal, field, format_spec, conversion in parsed:
            if field is None:
                continue
            if field not in allowed or format_spec or conversion:
                raise ValueError(
                    "%s permits only plain placeholders: %s"
                    % (label, ", ".join(sorted(allowed)))
                )
            fields.add(field)
        if fields != required:
            raise ValueError(
                "%s must contain exactly these placeholders: %s"
                % (label, ", ".join(sorted(required)))
            )
        rendered = value.format(**{field: "VALUE" for field in allowed})
        path = Path(rendered)
        if (
            not value
            or path.is_absolute()
            or ".." in path.parts
            or (not allow_directories and path.name != rendered)
        ):
            raise ValueError("%s must be a safe relative path template" % label)


class ImputationEnginesConfig(StrictModel):
    pred_ld: PredLDConfig


class ImputationConfig(ModuleConfig):
    engine: Literal["pred_ld"]
    input_directory: Path | None = None
    ld_reference_directory: Path | None = None
    genome_build: GenomeBuild | None = None
    population: Population
    post_harmonisation_directory: str
    engines: ImputationEnginesConfig

    @field_validator("post_harmonisation_directory")
    @classmethod
    def safe_output_directory(cls, value: str) -> str:
        path = Path(value)
        if not value or path.is_absolute() or ".." in path.parts:
            raise ValueError(
                "must be a non-empty relative path below the pipeline output"
            )
        return value
