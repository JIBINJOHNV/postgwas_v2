"""Versioned, normalized sample-sheet schema for harmonisation.

This preflight deliberately never opens a GWAS or external data file. It only
normalizes the sample sheet, validates row semantics, and checks file metadata.
Header and scientific inference checks belong to the later execution preflight.
"""

from __future__ import annotations

import csv
import re
import unicodedata
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, PrivateAttr, field_validator, model_validator

from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.core.values import MISSING_TEXT


INFERABLE_FIELDS = ("trait_type", "effect_type", "p_value_type", "delimiter")

ENUM_ALIASES = {
    "trait_type": {
        "case-control": "case_control",
        "case control": "case_control",
        "binary": "case_control",
        "continuous": "quantitative",
    },
    "effect_type": {
        "odds ratio": "odds_ratio",
        "odds-ratio": "odds_ratio",
        "or": "odds_ratio",
    },
    "p_value_type": {
        "p": "raw",
        "pvalue": "raw",
        "p-value": "raw",
        "-log10": "neglog10",
        "-log10p": "neglog10",
        "mlogp": "neglog10",
        "-ln": "negln",
    },
    "delimiter": {
        "\\t": "tab",
        "tab-separated": "tab",
        "tab separated": "tab",
        "csv": "comma",
        "comma-separated": "comma",
        "white space": "whitespace",
    },
}

ALLOWED_ENUMS = {
    "trait_type": {"auto", "quantitative", "case_control"},
    "effect_type": {"auto", "beta", "odds_ratio"},
    "p_value_type": {"auto", "raw", "neglog10", "negln"},
    "delimiter": {"auto", "tab", "comma", "semicolon", "space", "whitespace"},
}


def _unicode(value: Any) -> Any:
    return unicodedata.normalize("NFKC", value).strip() if isinstance(value, str) else value


def normalise_header(value: str) -> str:
    normalized = _unicode(value).lower().replace("-", " ")
    return re.sub(r"\s+", "_", normalized)


def _missing(value: Any) -> bool:
    return value is None or (
        isinstance(value, str) and _unicode(value).upper() in MISSING_TEXT
    )


def _normalise_inferable(field: str, value: Any, warnings: list[str], row_number: int) -> str:
    if _missing(value):
        return "auto"
    text = _unicode(value).lower()
    text = ENUM_ALIASES[field].get(text, text.replace("-", "_").replace(" ", "_"))
    if text not in ALLOWED_ENUMS[field]:
        warnings.append(
            "Row %d: invalid %s=%r was normalized to 'auto'." % (row_number, field, value)
        )
        return "auto"
    return text


def _normalise_count(value: Any, field: str, row_number: int) -> int | None:
    if _missing(value):
        return None
    text = _unicode(str(value)).replace(",", "")
    if not re.fullmatch(r"[0-9]+", text):
        raise ConfigurationError(
            "Row %d: %s must be a positive whole number; received %r."
            % (row_number, field, value)
        )
    result = int(text)
    if result <= 0:
        raise ConfigurationError("Row %d: %s must be greater than zero." % (row_number, field))
    return result


def _normalise_path(value: Any, sample_sheet: Path) -> str | None:
    if _missing(value):
        return None
    candidate = Path(_unicode(str(value))).expanduser()
    if not candidate.is_absolute():
        candidate = sample_sheet.parent / candidate
    return str(candidate.resolve(strict=False))


def _require_nonempty_file(path: Path, field: str, row_number: int) -> None:
    """Check metadata only; never open or read the target file."""
    try:
        stat = path.stat()
    except OSError as exc:
        raise ConfigurationError(
            "Row %d: %s does not exist or cannot be inspected: %s (%s)"
            % (row_number, field, path, exc)
        ) from exc
    if not path.is_file():
        raise ConfigurationError("Row %d: %s is not a regular file: %s" % (row_number, field, path))
    if stat.st_size <= 0:
        raise ConfigurationError("Row %d: %s is empty: %s" % (row_number, field, path))


class HarmonisationSampleSheetRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    config_version: int
    dataset_id: str
    input_file: Path
    trait_type: Literal["auto", "quantitative", "case_control"] = "auto"
    chromosome_column: str | None = None
    position_column: str | None = None
    chromosome_position_column: str | None = None
    variant_id_column: str | None = None
    effect_allele_column: str
    other_allele_column: str
    effect_allele_frequency_column: str | None = None
    effect_column: str | None = None
    effect_type: Literal["auto", "beta", "odds_ratio"] = "auto"
    standard_error_column: str | None = None
    z_score_column: str | None = None
    p_value_column: str
    p_value_type: Literal["auto", "raw", "neglog10", "negln"] = "auto"
    control_count_column: str | None = None
    case_count_column: str | None = None
    control_count: int | None = None
    case_count: int | None = None
    imputation_info_column: str | None = None
    external_info_file: Path | None = None
    external_info_column: str | None = None
    external_eaf_file: Path | None = None
    external_eaf_column: str | None = None
    delimiter: Literal["auto", "tab", "comma", "semicolon", "space", "whitespace"] = "auto"
    _normalisation_warnings: list[str] = PrivateAttr(default_factory=list)

    @field_validator("*", mode="before")
    @classmethod
    def normalise_missing_and_text(cls, value):
        if _missing(value):
            return None
        return _unicode(value)

    @field_validator("dataset_id")
    @classmethod
    def validate_dataset_id(cls, value: str) -> str:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
            raise ValueError("use only letters, numbers, '.', '_' and '-'")
        return value

    @model_validator(mode="after")
    def validate_alternatives(self):
        if self.config_version != 2:
            raise ValueError("config_version must be 2")
        if not ((self.chromosome_column and self.position_column) or self.chromosome_position_column):
            raise ValueError(
                "provide chromosome_column plus position_column, or chromosome_position_column"
            )
        if not self.effect_column and not self.z_score_column:
            raise ValueError("provide effect_column or z_score_column")

        internal_eaf = self.effect_allele_frequency_column is not None
        external_eaf = self.external_eaf_file is not None or self.external_eaf_column is not None
        if external_eaf and not (self.external_eaf_file and self.external_eaf_column):
            raise ValueError("external EAF requires both external_eaf_file and external_eaf_column")
        if internal_eaf == external_eaf:
            raise ValueError("provide exactly one EAF source: internal column XOR external file+column")

        internal_info = self.imputation_info_column is not None
        external_info = self.external_info_file is not None or self.external_info_column is not None
        if (
            not internal_info
            and external_info
            and not (self.external_info_file and self.external_info_column)
        ):
            raise ValueError("external INFO requires both external_info_file and external_info_column")

        has_controls = self.control_count_column is not None or self.control_count is not None
        has_cases = self.case_count_column is not None or self.case_count is not None
        if not has_controls:
            raise ValueError("provide control_count_column or control_count")
        if self.trait_type == "case_control" and not has_cases:
            raise ValueError("case-control traits require case_count_column or case_count")
        if self.trait_type == "quantitative" and has_cases:
            raise ValueError("quantitative traits must not provide case-count fields")
        return self

    @property
    def normalisation_warnings(self) -> tuple[str, ...]:
        return tuple(self._normalisation_warnings)


def _read_records(sample_sheet: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with sample_sheet.open("r", encoding="utf-8", newline="") as handle:
            sample = handle.read(8192)
            if not sample.strip():
                raise ConfigurationError("Sample sheet is empty: %s" % sample_sheet)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            except csv.Error:
                dialect = csv.excel
            handle.seek(0)
            reader = csv.reader(handle, dialect)
            raw_header = next(reader)
            records = []
            for line_number, row in enumerate(reader, start=2):
                if not any(cell.strip() for cell in row):
                    continue
                if len(row) != len(raw_header):
                    raise ConfigurationError(
                        "Sample-sheet row %d has %d values but the header has %d columns."
                        % (line_number, len(row), len(raw_header))
                    )
                records.append(dict(zip(raw_header, row)))
    except (OSError, csv.Error) as exc:
        raise ConfigurationError("Cannot parse sample sheet %s: %s" % (sample_sheet, exc)) from exc
    if not records:
        raise ConfigurationError("Sample sheet has no dataset rows: %s" % sample_sheet)
    return raw_header, records


def load_harmonisation_sample_sheet(path: str | Path) -> list[HarmonisationSampleSheetRow]:
    """Normalize and validate rows without opening any referenced data file."""
    sample_sheet = Path(path).expanduser().resolve()
    raw_header, raw_records = _read_records(sample_sheet)
    normalized_header = [normalise_header(name) for name in raw_header]
    duplicates = sorted({name for name in normalized_header if normalized_header.count(name) > 1})
    if duplicates:
        raise ConfigurationError("Duplicate normalized sample-sheet columns: " + ", ".join(duplicates))

    rows: list[HarmonisationSampleSheetRow] = []
    seen_ids: dict[str, int] = {}
    seen_inputs: dict[Path, int] = {}
    normalized_names = dict(zip(raw_header, normalized_header))
    for row_number, raw in enumerate(raw_records, start=2):
        record = {normalized_names[key]: value for key, value in raw.items()}
        warnings: list[str] = []
        for field in INFERABLE_FIELDS:
            record[field] = _normalise_inferable(field, record.get(field), warnings, row_number)
        for field in ("control_count", "case_count"):
            record[field] = _normalise_count(record.get(field), field, row_number)
        for field in ("input_file", "external_info_file", "external_eaf_file"):
            record[field] = _normalise_path(record.get(field), sample_sheet)
        try:
            row = HarmonisationSampleSheetRow.model_validate(record)
        except Exception as exc:
            raise ConfigurationError(
                "Invalid sample-sheet row %d in %s: %s" % (row_number, sample_sheet, exc)
            ) from exc
        row._normalisation_warnings.extend(warnings)
        if row.imputation_info_column and (
            row.external_info_file or row.external_info_column
        ):
            row._normalisation_warnings.append(
                "Row %d: internal INFO column %r has priority; the listed "
                "external INFO source will not be used."
                % (row_number, row.imputation_info_column)
            )

        _require_nonempty_file(row.input_file, "input_file", row_number)
        if row.external_eaf_file:
            _require_nonempty_file(row.external_eaf_file, "external_eaf_file", row_number)
        if row.external_info_file and row.imputation_info_column is None:
            _require_nonempty_file(row.external_info_file, "external_info_file", row_number)

        normalized_id = row.dataset_id.casefold()
        if normalized_id in seen_ids:
            raise ConfigurationError(
                "Duplicate dataset_id %r in rows %d and %d."
                % (row.dataset_id, seen_ids[normalized_id], row_number)
            )
        seen_ids[normalized_id] = row_number

        normalized_input = row.input_file.resolve(strict=False)
        if normalized_input in seen_inputs:
            raise ConfigurationError(
                "Duplicate input_file %s in rows %d and %d."
                % (normalized_input, seen_inputs[normalized_input], row_number)
            )
        seen_inputs[normalized_input] = row_number
        rows.append(row)
    return rows


def to_harmonisation_input(
    row: HarmonisationSampleSheetRow,
    *,
    resource_directory: str | Path,
    output_directory: str | Path,
    output_layout: dict[str, str],
    fixed_info: float | None = None,
    fixed_info_column: str | None = None,
) -> dict[str, Any]:
    """Translate a validated sample-sheet row into the engine input contract."""
    value = lambda item: None if item is None else str(item)
    internal_info = row.imputation_info_column is not None
    external_info = bool(row.external_info_file and row.external_info_column)
    if internal_info:
        info_source = "internal"
        info_source_detail = "internal column %s" % row.imputation_info_column
        info_file = None
        info_column = None
        resolved_fixed_info = None
    elif external_info:
        info_source = "external"
        info_source_detail = "external file %s, column %s" % (
            row.external_info_file,
            row.external_info_column,
        )
        info_file = row.external_info_file
        info_column = row.external_info_column
        resolved_fixed_info = None
    elif fixed_info is not None:
        if fixed_info_column is None or not str(fixed_info_column).strip():
            raise ConfigurationError(
                "The fixed INFO working column is not configured."
            )
        info_source = "fixed_cli"
        info_source_detail = "--fixed-info %g" % float(fixed_info)
        info_file = None
        info_column = None
        resolved_fixed_info = float(fixed_info)
    else:
        raise ConfigurationError(
            "Dataset %s has no imputation-quality source. Provide "
            "imputation_info_column, provide external_info_file together with "
            "external_info_column, or explicitly supply --fixed-info VALUE."
            % row.dataset_id
        )
    return {
        "sumstat_file": str(row.input_file),
        "gwas_outputname": row.dataset_id,
        "trait_type": row.trait_type,
        "declared_effect_type": row.effect_type,
        "declared_pvalue_type": row.p_value_type,
        "delimiter": row.delimiter,
        "chr_col": value(row.chromosome_column),
        "pos_col": value(row.position_column),
        "chr_pos_col": value(row.chromosome_position_column),
        "snp_id_col": value(row.variant_id_column),
        "ea_col": row.effect_allele_column,
        "oa_col": row.other_allele_column,
        "eaf_col": value(row.effect_allele_frequency_column),
        "beta_or_col": value(row.effect_column),
        "se_col": value(row.standard_error_column),
        "imp_z_col": value(row.z_score_column),
        "pval_col": row.p_value_column,
        "ncontrol_col": value(row.control_count_column),
        "ncase_col": value(row.case_count_column),
        "ncontrol": value(row.control_count),
        "ncase": value(row.case_count),
        "imp_info_col": value(row.imputation_info_column),
        "infofile": value(info_file),
        "infocolumn": value(info_column),
        "fixed_info": resolved_fixed_info,
        "fixed_info_column": (
            str(fixed_info_column) if resolved_fixed_info is not None else None
        ),
        "info_source": info_source,
        "info_source_detail": info_source_detail,
        "provided_external_info_file": value(row.external_info_file),
        "eaffile": value(row.external_eaf_file),
        "eafcolumn": value(row.external_eaf_column),
        "provided_external_eaf_file": value(row.external_eaf_file),
        "resource_folder": str(Path(resource_directory).expanduser().resolve()),
        "output_folder": str(configured_output_path(
            output_directory,
            output_layout["dataset_directory"],
            error_type=ConfigurationError,
            dataset_id=row.dataset_id,
        )),
    }
