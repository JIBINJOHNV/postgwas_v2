import re
from pathlib import Path

from pydantic import ValidationInfo, field_validator

from postgwas.config.models.common import ModuleConfig, Population, StrictModel
from postgwas.core.paths import validate_filename_component


class LDAnnotationInputsConfig(StrictModel):
    """Run-specific inputs for population LD-block annotation."""

    vcf: Path | None = None
    ld_region_dir: Path | None = None
    dataset_id: str | None = None

    @field_validator("vcf", "ld_region_dir", mode="before")
    @classmethod
    def empty_paths_are_missing(cls, value):
        if isinstance(value, str) and not value.strip():
            return None
        return value

    @field_validator("dataset_id")
    @classmethod
    def safe_dataset_id(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return validate_filename_component(
            value, "modules.ld_annotation.inputs.dataset_id"
        )


class LDAnnotationConfig(ModuleConfig):
    inputs: LDAnnotationInputsConfig
    output_directory: Path | None = None
    populations: list[Population]
    include_unassigned: bool
    bed_filename_template: str
    info_field_template: str
    info_description_template: str
    output_filename_template: str
    summary_filename_template: str
    html_report_filename_template: str
    canonical_log_filename_template: str

    @field_validator("populations")
    @classmethod
    def unique_populations(cls, values: list[Population]) -> list[Population]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain one or more unique populations")
        return values

    @field_validator("bed_filename_template")
    @classmethod
    def valid_bed_filename_template(cls, value: str) -> str:
        value = value.strip()
        required = ("{genome_build}", "{population}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain {genome_build} and {population} exactly once"
            )
        try:
            rendered = Path(
                value.format(genome_build="BUILD", population="POPULATION")
            )
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if (
            rendered.is_absolute()
            or rendered.name != str(rendered)
            or "/" in value
            or "\\" in value
        ):
            raise ValueError("must render one filename inside --ld-region-dir")
        # The BED suffix is a protocol invariant: bcftools uses it to apply BED's
        # zero-based, half-open coordinate convention to annotation intervals.
        if not str(rendered).endswith(".bed.gz"):
            raise ValueError("must render a .bed.gz filename")
        return value

    @field_validator("info_field_template")
    @classmethod
    def valid_info_field_template(cls, value: str) -> str:
        value = value.strip()
        if value.count("{population}") != 1:
            raise ValueError("must contain {population} exactly once")
        try:
            rendered = value.format(population="EUR")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", rendered):
            raise ValueError("must render a valid VCF INFO identifier")
        return value

    @field_validator("info_description_template")
    @classmethod
    def valid_info_description_template(cls, value: str) -> str:
        value = value.strip()
        required = ("{population}", "{genome_build}")
        if any(value.count(token) != 1 for token in required):
            raise ValueError(
                "must contain {population} and {genome_build} exactly once"
            )
        try:
            rendered = value.format(population="EUR", genome_build="GRCh37")
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if not rendered or any(
            character in rendered for character in ('"', "\r", "\n")
        ):
            raise ValueError(
                "must render a nonempty single-line VCF description without quotes"
            )
        return value

    @field_validator(
        "output_filename_template",
        "summary_filename_template",
        "html_report_filename_template",
        "canonical_log_filename_template",
    )
    @classmethod
    def valid_output_filename_template(
        cls,
        value: str,
        info: ValidationInfo,
    ) -> str:
        value = value.strip()
        if value.count("{dataset_id}") != 1:
            raise ValueError("must contain {dataset_id} exactly once")
        try:
            rendered = Path(value.format(dataset_id="dataset"))
        except (KeyError, ValueError) as exc:
            raise ValueError("contains an unsupported placeholder: %s" % exc) from exc
        if (
            rendered.is_absolute()
            or rendered.name != str(rendered)
            or "/" in value
            or "\\" in value
        ):
            raise ValueError("must render one filename inside the output directory")
        required_suffix = {
            "output_filename_template": ".vcf.gz",
            "summary_filename_template": ".csv",
            "html_report_filename_template": ".html",
            "canonical_log_filename_template": ".log",
        }[info.field_name]
        if not str(rendered).endswith(required_suffix):
            raise ValueError("must render a %s filename" % required_suffix)
        return value
