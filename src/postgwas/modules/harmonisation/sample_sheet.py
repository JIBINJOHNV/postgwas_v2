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

from pydantic import (
    BaseModel,
    ConfigDict,
    PrivateAttr,
    ValidationInfo,
    ValidationError,
    field_validator,
    model_validator,
)

from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.core.values import MISSING_TEXT
from postgwas.modules.harmonisation.resource_paths import (
    validate_external_resource_spec,
)


INFERABLE_FIELDS = ("trait_type", "effect_type", "p_value_type", "delimiter")
SAMPLE_SHEET_VERSION = 2

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
    "p_value_type": {"auto", "raw", "neglog10"},
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


def _require_external_resource(path: Path, field: str, row_number: int) -> None:
    """Accept one non-empty file or an explicit chromosome/build template."""
    try:
        validate_external_resource_spec(path)
    except ValueError as exc:
        raise ConfigurationError(
            "Row %d: %s is invalid. %s" % (row_number, field, exc)
        ) from exc


class HarmonisationSampleSheetRow(BaseModel):
    # User metadata may accompany mappings without entering the analysis contract.
    model_config = ConfigDict(extra="ignore")

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
    p_value_type: Literal["auto", "raw", "neglog10"] = "auto"
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
    _row_number: int | None = PrivateAttr(default=None)

    @field_validator("*", mode="before")
    @classmethod
    def normalise_missing_and_text(cls, value):
        if _missing(value):
            return None
        return _unicode(value)

    @field_validator("dataset_id")
    @classmethod
    def validate_dataset_id(cls, value: str) -> str:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._+-]*", value):
            raise ValueError("use only letters, numbers, '.', '_', '+' and '-'")
        return value

    @model_validator(mode="after")
    def validate_alternatives(self, info: ValidationInfo):
        problems = []
        if self.config_version != SAMPLE_SHEET_VERSION:
            problems.append("config_version must be %d" % SAMPLE_SHEET_VERSION)
        draft = bool(info.context and info.context.get("allow_missing_mappings"))
        if not draft and not ((self.chromosome_column and self.position_column) or self.chromosome_position_column):
            problems.append(
                "provide chromosome_column plus position_column, or chromosome_position_column"
            )
        if not draft and not self.effect_column and not self.z_score_column:
            problems.append("provide effect_column or z_score_column")

        internal_eaf = self.effect_allele_frequency_column is not None
        external_eaf = self.external_eaf_file is not None or self.external_eaf_column is not None
        if external_eaf and not (self.external_eaf_file and self.external_eaf_column):
            missing = (
                "external_eaf_file"
                if self.external_eaf_file is None
                else "external_eaf_column"
            )
            problems.append(
                "Dataset %r has an incomplete external EAF source: %s is missing. "
                "Set both external_eaf_file and external_eaf_column, or set both "
                "to NA and provide effect_allele_frequency_column instead."
                % (self.dataset_id, missing)
            )
        draft_allows_missing_eaf = bool(
            info.context and info.context.get("allow_missing_eaf")
        )
        if internal_eaf and external_eaf:
            problems.append(
                "Dataset %r provides two allele-frequency sources: internal "
                "effect_allele_frequency_column=%r and external_eaf_file plus "
                "external_eaf_column. Keep exactly one source and set the other "
                "alternative to NA."
                % (self.dataset_id, self.effect_allele_frequency_column)
            )
        if not internal_eaf and not external_eaf and not draft_allows_missing_eaf:
            problems.append(
                "Dataset %r has no allele-frequency source. Complete the sample "
                "sheet by setting either effect_allele_frequency_column to the "
                "exact study EAF/MAF column name, or both external_eaf_file and "
                "external_eaf_column. All three fields are currently NA."
                % self.dataset_id
            )

        internal_info = self.imputation_info_column is not None
        external_info = self.external_info_file is not None or self.external_info_column is not None
        if (
            not internal_info
            and external_info
            and not (self.external_info_file and self.external_info_column)
        ):
            missing = (
                "external_info_file"
                if self.external_info_file is None
                else "external_info_column"
            )
            problems.append(
                "Dataset %r has an incomplete external INFO source: %s is "
                "missing. Set both external_info_file and external_info_column, "
                "or set both to NA and provide imputation_info_column or the "
                "--fixed-info option instead."
                % (self.dataset_id, missing)
            )

        has_controls = self.control_count_column is not None or self.control_count is not None
        has_cases = self.case_count_column is not None or self.case_count is not None
        draft_allows_missing_sample_size = bool(
            info.context and info.context.get("allow_missing_sample_size")
        )
        if not has_controls and not draft_allows_missing_sample_size:
            problems.append("provide control_count_column or control_count")
        if self.trait_type == "case_control" and not has_cases and not draft_allows_missing_sample_size:
            problems.append("case-control traits require case_count_column or case_count")
        if self.trait_type == "quantitative" and has_cases:
            problems.append("quantitative traits must not provide case-count fields")
        if problems:
            raise ValueError("; ".join(problems))
        return self

    @classmethod
    def validate_generated_draft(cls, value: Any) -> "HarmonisationSampleSheetRow":
        """Validate a generator draft while leaving user-completable fields unresolved."""
        return _HarmonisationSampleSheetDraftRow.model_validate(
            value,
            context={
                "allow_missing_mappings": True,
                "allow_missing_sample_size": True,
                "allow_missing_eaf": True,
            },
        )

    @property
    def normalisation_warnings(self) -> tuple[str, ...]:
        return tuple(self._normalisation_warnings)


class _HarmonisationSampleSheetDraftRow(HarmonisationSampleSheetRow):
    """Unresolved mapping cells are permitted only during draft publication.

    Normal sample-sheet loading continues to use the required fields above.
    """

    effect_allele_column: str | None = None
    other_allele_column: str | None = None
    p_value_column: str | None = None


def _read_records(sample_sheet: Path, issues: list[dict]) -> tuple[list[str], list[tuple[int, dict[str, str]]]]:
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
            for row in reader:
                line_number = reader.line_num
                if not any(cell.strip() for cell in row):
                    continue
                if len(row) != len(raw_header):
                    issues.append({
                        "row": line_number, "dataset_id": "", "input_file": "",
                        "stage": "sample sheet", "status": "failed",
                        "problem": "Sample-sheet row %d has %d values but the header has %d columns."
                        % (line_number, len(row), len(raw_header)),
                    })
                    continue
                records.append((line_number, dict(zip(raw_header, row))))
    except (OSError, csv.Error) as exc:
        raise ConfigurationError("Cannot parse sample sheet %s: %s" % (sample_sheet, exc)) from exc
    if not records and not issues:
        raise ConfigurationError("Sample sheet has no dataset rows: %s" % sample_sheet)
    return raw_header, records


def _validation_error_text(exc: ValidationError) -> str:
    """Render Pydantic row errors without implementation details or help URLs."""
    messages: list[str] = []
    for error in exc.errors(
        include_url=False,
        include_context=False,
        include_input=False,
    ):
        location = ".".join(str(part) for part in error.get("loc", ()))
        message = str(error.get("msg") or "invalid value")
        if message.startswith("Value error, "):
            message = message[len("Value error, "):]
        detail = "%s: %s" % (location, message) if location else message
        if detail not in messages:
            messages.append(detail)
    return "; ".join(messages) or str(exc)


def load_harmonisation_sample_sheet(
    path: str | Path, *, issues: list[dict] | None = None,
) -> list[HarmonisationSampleSheetRow]:
    """Check every row; optionally collect issues for the full-run preflight.

    Without a collector, preserve the public fail-on-invalid-sheet contract,
    but raise only after every independently readable row has been checked.
    """
    collected = [] if issues is None else issues
    sample_sheet = Path(path).expanduser().resolve()
    raw_header, raw_records = _read_records(sample_sheet, collected)
    normalized_header = [normalise_header(name) for name in raw_header]
    duplicates = sorted({name for name in normalized_header if normalized_header.count(name) > 1})
    if duplicates:
        raise ConfigurationError("Duplicate normalized sample-sheet columns: " + ", ".join(duplicates))

    rows: list[HarmonisationSampleSheetRow] = []
    seen_ids: dict[str, int] = {}
    seen_inputs: dict[str, int] = {}
    normalized_names = dict(zip(raw_header, normalized_header))
    for row_number, raw in raw_records:
        record = {normalized_names[key]: value for key, value in raw.items()}
        row_problems = []
        warnings: list[str] = []
        for field in INFERABLE_FIELDS:
            record[field] = _normalise_inferable(field, record.get(field), warnings, row_number)
        for field in ("control_count", "case_count"):
            try:
                record[field] = _normalise_count(record.get(field), field, row_number)
            except ConfigurationError as exc:
                row_problems.append(str(exc))
                record[field] = None
        for field in ("input_file", "external_info_file", "external_eaf_file"):
            try:
                record[field] = _normalise_path(record.get(field), sample_sheet)
            except (OSError, ValueError, RuntimeError) as exc:
                row_problems.append("Invalid %s: %s" % (field, exc))
                record[field] = None
        try:
            row = HarmonisationSampleSheetRow.model_validate(record)
        except ValidationError as exc:
            row = None
            row_problems.append(_validation_error_text(exc))
        # File availability and duplicate identities do not depend on a valid
        # scientific mapping. Check them even when another cell is invalid.
        for field in ("input_file", "external_eaf_file", "external_info_file"):
            value = record.get(field)
            if value is None or (field == "external_info_file" and not _missing(record.get("imputation_info_column"))):
                continue
            try:
                check = _require_nonempty_file if field == "input_file" else _require_external_resource
                check(Path(value), field, row_number)
            except (ConfigurationError, OSError, ValueError, RuntimeError) as exc:
                row_problems.append(str(exc))
        for field, seen, identity in (
            ("dataset_id", seen_ids, str(record.get("dataset_id") or "").casefold()),
            ("input_file", seen_inputs, record.get("input_file")),
        ):
            if not identity:
                continue
            if identity in seen:
                row_problems.append("Duplicate %s %r in rows %d and %d." % (
                    field, record[field], seen[identity], row_number,
                ))
            else:
                seen[identity] = row_number
        if row_problems:
            for problem in row_problems:
                collected.append({
                    "row": row_number, "dataset_id": record.get("dataset_id") or "",
                    "input_file": record.get("input_file") or "",
                    "stage": "sample sheet", "status": "failed",
                    "problem": "Invalid sample-sheet row %d for dataset %r: %s"
                    % (row_number, record.get("dataset_id"), problem),
                })
            continue
        row._row_number = row_number
        row._normalisation_warnings.extend(warnings)
        if row.imputation_info_column and (row.external_info_file or row.external_info_column):
            row._normalisation_warnings.append(
                "Row %d: internal INFO column %r has priority; the listed external INFO source will not be used."
                % (row_number, row.imputation_info_column)
            )
        rows.append(row)
    if issues is None and collected:
        raise ConfigurationError("Sample-sheet validation failed:\n" + "\n".join(
            item["problem"] for item in collected
        ))
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
        "output_root": str(Path(output_directory).expanduser().resolve()),
        "output_folder": str(configured_output_path(
            output_directory,
            output_layout["dataset_directory"],
            error_type=ConfigurationError,
            dataset_id=row.dataset_id,
        )),
    }
