import re
from typing import Literal, get_args

from pydantic import Field, RootModel, field_validator, model_validator

from postgwas.config.models.common import ModuleConfig, StrictModel
from postgwas.core.variant_identifiers import (
    IDENTIFIER_TEMPLATE_FIELDS,
    VariantIdentifierType,
    identifier_template_parts,
)


FormattingTarget = Literal[
    "magma", "susie", "finemap", "pred_ld", "ldsc", "mixer", "gcta_gene",
]
FormattingTransform = Literal[
    "negative_log10_to_raw_p", "effect_frequency_to_minor_frequency",
]


class FormattingVcfFields(RootModel[dict[str, str]]):
    """Canonical table column -> bcftools query expression."""

    @model_validator(mode="after")
    def validate_projection(self):
        if not self.root:
            raise ValueError("vcf_fields must contain at least one configured field")
        invalid = [
            key for key, value in self.root.items()
            if not str(key).strip() or not str(value).strip()
        ]
        if invalid:
            raise ValueError("vcf_fields contains an empty column name or query")
        return self


class FormattingValidationConfig(StrictModel):
    required_columns: list[str] = Field(default_factory=list)
    positive_columns: list[str] = Field(default_factory=list)
    nonnegative_columns: list[str] = Field(default_factory=list)
    open_unit_interval_columns: list[str] = Field(default_factory=list)
    closed_unit_interval_columns: list[str] = Field(default_factory=list)

    @field_validator("*")
    @classmethod
    def unique_columns(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)):
            raise ValueError("must not contain duplicate column names")
        return values


class FormattingTableConfig(StrictModel):
    output_file: str
    columns: dict[str, str]
    transformations: dict[str, FormattingTransform] = Field(default_factory=dict)
    integer_columns: list[str] = Field(default_factory=list)
    validation: FormattingValidationConfig = Field(
        default_factory=FormattingValidationConfig
    )

    @model_validator(mode="after")
    def validate_mapping(self):
        if not self.output_file.strip():
            raise ValueError("output_file must not be empty")
        if not self.columns:
            raise ValueError("columns must contain at least one source-to-output mapping")
        unknown = sorted(set(self.transformations) - set(self.columns))
        if unknown:
            raise ValueError(
                "transformations reference unmapped source columns: %s"
                % ", ".join(unknown)
            )
        unknown = sorted(set(self.integer_columns) - set(self.columns))
        if unknown:
            raise ValueError(
                "integer_columns reference unmapped source columns: %s"
                % ", ".join(unknown)
            )
        if len(set(self.columns.values())) != len(self.columns):
            raise ValueError("output column names must be unique")
        return self


class FormattingExportConfig(StrictModel):
    output_file: str | None = None
    output_directory: str | None = None
    partition_file: str | None = None
    partition_by: str | None = None
    columns: dict[str, str] = Field(default_factory=dict)
    trait_columns: dict[Literal["binary", "quantitative"], dict[str, str]] = Field(
        default_factory=dict
    )
    trailing_columns: dict[str, str] = Field(default_factory=dict)
    transformations: dict[str, FormattingTransform] = Field(default_factory=dict)
    integer_columns: list[str] = Field(default_factory=list)
    validation: FormattingValidationConfig = Field(
        default_factory=FormattingValidationConfig
    )
    trait_required_columns: dict[
        Literal["binary", "quantitative"], list[str]
    ] = Field(default_factory=dict)
    outputs: dict[str, FormattingTableConfig] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_export(self):
        if not self.columns and not self.outputs:
            raise ValueError("an export needs columns or named outputs")
        if self.columns and not self.output_file and not self.partition_file:
            raise ValueError("a mapped export needs output_file or partition_file")
        mapped = set(self.columns) | set(self.trailing_columns)
        for values in self.trait_columns.values():
            mapped.update(values)
        unknown = sorted(set(self.transformations) - mapped)
        if unknown:
            raise ValueError(
                "transformations reference unmapped source columns: %s"
                % ", ".join(unknown)
            )
        unknown = sorted(set(self.integer_columns) - mapped)
        if unknown:
            raise ValueError(
                "integer_columns reference unmapped source columns: %s"
                % ", ".join(unknown)
            )
        return self


class MixerFormattingConfig(StrictModel):
    minimum_info: float | None = Field(default=None, ge=0, le=1)
    minimum_sample_size_fraction: float = Field(gt=0, le=1)
    snps_only: bool
    chromosome_column: str
    position_column: str
    effect_allele_column: str
    other_allele_column: str
    sample_size_column: str
    info_column: str
    allowed_alleles: list[str]


class FormattingCanonicalColumns(StrictModel):
    chromosome: str
    position: str
    variant_id: str
    reference_allele: str
    alternate_allele: str
    resolved_variant_id: str


class FormattingChromosomeLabels(StrictModel):
    """Configured normalization used before coordinates enter any export."""

    prefix_pattern: str
    aliases: dict[str, str]

    @field_validator("prefix_pattern")
    @classmethod
    def valid_prefix_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("aliases")
    @classmethod
    def nonempty_aliases(cls, values: dict[str, str]) -> dict[str, str]:
        if any(
            not source.strip() or not target.strip()
            for source, target in values.items()
        ):
            raise ValueError("must not contain empty chromosome labels")
        return values


class FormattingVariantIdentifierConfig(StrictModel):
    """General formatter policy for selecting downstream variant identifiers."""

    default_type: VariantIdentifierType
    target_types: dict[FormattingTarget, VariantIdentifierType]
    rsid_pattern: str
    rsid_extraction_pattern: str
    unique_id_template: str

    @field_validator("rsid_pattern")
    @classmethod
    def valid_rsid_pattern(cls, value: str) -> str:
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        return value

    @field_validator("rsid_extraction_pattern")
    @classmethod
    def valid_rsid_extraction_pattern(cls, value: str) -> str:
        try:
            compiled = re.compile(value)
        except re.error as exc:
            raise ValueError("must be a valid regular expression") from exc
        if compiled.groups != 1:
            raise ValueError("must contain exactly one capturing group for the rsID")
        return value

    @field_validator("unique_id_template")
    @classmethod
    def valid_unique_id_template(cls, value: str) -> str:
        try:
            parsed = identifier_template_parts(value)
        except ValueError as exc:
            raise ValueError("must be a valid format template") from exc
        fields = [field for _, field, _, _ in parsed if field is not None]
        required = set(IDENTIFIER_TEMPLATE_FIELDS)
        if len(fields) != len(required) or set(fields) != required:
            raise ValueError(
                "must contain chromosome, position, reference_allele, and "
                "alternate_allele exactly once"
            )
        if any(
            format_spec or conversion
            for _, field, format_spec, conversion in parsed
            if field is not None
        ):
            raise ValueError("must not use conversions or format specifiers")
        return value


class FormattingStudyDesign(StrictModel):
    case_count_column: str
    control_count_column: str


class FormattingRuntimeConfig(StrictModel):
    log_file: str
    resolved_config_file: str
    completion_manifest_file: str
    temporary_table_prefix: str
    temporary_table_suffix: str
    atomic_output_suffix: str
    atomic_uncompressed_suffix: str
    table_delimiter: str
    input_null_values: list[str]
    output_null_value: str
    compressed_suffixes: list[str]
    io_buffer_bytes: int = Field(gt=0)

    @field_validator(
        "log_file", "resolved_config_file", "completion_manifest_file",
        "temporary_table_prefix",
        "temporary_table_suffix", "atomic_output_suffix",
        "atomic_uncompressed_suffix", "output_null_value",
    )
    @classmethod
    def nonempty_runtime_value(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("must not be empty")
        return value

    @field_validator("table_delimiter")
    @classmethod
    def single_delimiter(cls, value: str) -> str:
        if len(value) != 1:
            raise ValueError("table_delimiter must be exactly one character")
        return value

    @field_validator("input_null_values", "compressed_suffixes")
    @classmethod
    def nonempty_unique_values(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("must contain unique values and must not be empty")
        return values

    @model_validator(mode="after")
    def valid_temporary_prefix(self):
        if "{dataset_id}" not in self.temporary_table_prefix:
            raise ValueError("temporary_table_prefix must contain {dataset_id}")
        if "/" in self.temporary_table_prefix or "\\" in self.temporary_table_prefix:
            raise ValueError("temporary_table_prefix must be a filename prefix")
        for suffix in (self.atomic_output_suffix, self.atomic_uncompressed_suffix):
            if "/" in suffix or "\\" in suffix:
                raise ValueError("atomic output suffixes must not contain directories")
        if self.atomic_output_suffix == self.atomic_uncompressed_suffix:
            raise ValueError("atomic output suffixes must be different")
        return self


class FormattingResolvedConfig(StrictModel):
    """Configuration paths retained in target-scoped formatter metadata."""

    common_fields: list[str]
    format_fields: dict[FormattingTarget, list[str]]

    @field_validator("common_fields")
    @classmethod
    def unique_common_fields(cls, values: list[str]) -> list[str]:
        if not values or len(values) != len(set(values)):
            raise ValueError("common_fields must contain unique values")
        return values

    @field_validator("format_fields")
    @classmethod
    def complete_format_fields(
        cls, values: dict[FormattingTarget, list[str]],
    ) -> dict[FormattingTarget, list[str]]:
        expected = set(get_args(FormattingTarget))
        if set(values) != expected:
            raise ValueError(
                "format_fields must contain every formatter target exactly once"
            )
        for target, paths in values.items():
            if not paths or len(paths) != len(set(paths)):
                raise ValueError(
                    "format_fields.%s must contain unique values" % target
                )
            if "exports.%s" % target not in paths:
                raise ValueError(
                    "format_fields.%s must include exports.%s" % (target, target)
                )
        return values


class FormattingConfig(ModuleConfig):
    formats: list[FormattingTarget]
    format_order: list[FormattingTarget]
    module_formats: dict[str, list[FormattingTarget]]
    chromosomes: list[str]
    minimum_p_value: float = Field(gt=0, le=1)
    vcf_fields: FormattingVcfFields
    numeric_columns: list[str]
    canonical_columns: FormattingCanonicalColumns
    chromosome_labels: FormattingChromosomeLabels
    variant_identifiers: FormattingVariantIdentifierConfig
    study_design: FormattingStudyDesign
    runtime: FormattingRuntimeConfig
    vcf_include_expression: str | None = None
    exports: dict[FormattingTarget, FormattingExportConfig]
    mixer: MixerFormattingConfig
    resolved_config: FormattingResolvedConfig

    @field_validator("formats", "format_order", "chromosomes")
    @classmethod
    def unique_values(cls, values: list[str]) -> list[str]:
        if len(values) != len(set(values)):
            raise ValueError("must not contain duplicate values")
        return values

    @model_validator(mode="after")
    def all_targets_have_schemas(self):
        expected = set(get_args(FormattingTarget))
        missing = sorted(expected - set(self.exports))
        if missing:
            raise ValueError("exports is missing formatter targets: %s" % ", ".join(missing))
        if set(self.format_order) != expected:
            raise ValueError("format_order must contain every formatter target exactly once")
        available = self.model_dump(mode="python", exclude={"resolved_config"})
        for target, paths in self.resolved_config.format_fields.items():
            overlap = set(paths) & set(self.resolved_config.common_fields)
            if overlap:
                raise ValueError(
                    "resolved_config paths are repeated for %s: %s"
                    % (target, ", ".join(sorted(overlap)))
                )
        for path in [
            *self.resolved_config.common_fields,
            *(
                path
                for paths in self.resolved_config.format_fields.values()
                for path in paths
            ),
        ]:
            current = available
            for part in path.split("."):
                if not isinstance(current, dict) or part not in current:
                    raise ValueError(
                        "resolved_config references an unknown field: %s" % path
                    )
                current = current[part]
        for module, formats in self.module_formats.items():
            if not module.strip():
                raise ValueError("module_formats contains an empty module name")
            if len(formats) != len(set(formats)):
                raise ValueError(
                    "module_formats.%s must not contain duplicate formats" % module
                )
        known = set(self.vcf_fields.root) | {self.canonical_columns.resolved_variant_id}
        references = set(self.numeric_columns)
        references.update(self.canonical_columns.model_dump().values())
        references.update(self.study_design.model_dump().values())
        references.update({
            self.mixer.chromosome_column,
            self.mixer.position_column,
            self.mixer.effect_allele_column,
            self.mixer.other_allele_column,
            self.mixer.sample_size_column,
            self.mixer.info_column,
        })
        for schema in self.exports.values():
            references.update(schema.columns)
            references.update(schema.trailing_columns)
            validation = schema.validation
            references.update(validation.required_columns)
            references.update(validation.positive_columns)
            references.update(validation.nonnegative_columns)
            references.update(validation.open_unit_interval_columns)
            references.update(validation.closed_unit_interval_columns)
            for values in schema.trait_columns.values():
                references.update(values)
            for values in schema.trait_required_columns.values():
                references.update(values)
            for output in schema.outputs.values():
                references.update(output.columns)
                references.update(output.validation.required_columns)
                references.update(output.validation.positive_columns)
                references.update(output.validation.nonnegative_columns)
                references.update(output.validation.open_unit_interval_columns)
                references.update(output.validation.closed_unit_interval_columns)
        unknown = sorted(references - known)
        if unknown:
            raise ValueError(
                "formatter schemas reference unknown canonical columns: %s"
                % ", ".join(unknown)
            )
        return self
