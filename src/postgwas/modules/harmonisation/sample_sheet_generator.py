"""Generate a reviewed v2 harmonisation sample-sheet draft from GWAS headers."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import os
from pathlib import Path
import re
import tempfile
import unicodedata
from typing import Any, Sequence

from pydantic import ValidationError

from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.required_arguments import RequiredArgument, require_resolved_arguments
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples, mark_cli_required_help
from postgwas.core.io.delimiters import (
    delimiter_character,
    detect_delimiter,
    open_text,
)

from .resource_paths import validate_external_resource_spec
from .policies import load_policies
from .sample_sheet import (
    HarmonisationSampleSheetRow,
    SAMPLE_SHEET_VERSION,
    _validation_error_text,
    load_harmonisation_sample_sheet,
)


# Sample-sheet protocol fields: both sources share file/column pairs.
_EXTERNAL_SOURCES = {
    "eaf": "effect_allele_frequency_column",
    "info": "imputation_info_column",
}


@dataclass(frozen=True)
class SampleSheetDatasetStatus:
    dataset_id: str
    input_file: Path
    missing_eaf: bool
    missing_sample_size: bool
    missing_info: bool
    missing_mappings: tuple[str, ...]
    has_case_control_columns: bool
    has_coordinates: bool
    warnings: tuple[str, ...]
    eaf_source: str = "missing"
    info_source: str = "missing"

    @property
    def needs_attention(self) -> bool:
        return bool(
            self.missing_eaf
            or self.missing_sample_size
            or self.missing_info
            or self.missing_mappings
            or self.warnings
        )


@dataclass(frozen=True)
class SampleSheetFileRejection:
    input_file: Path
    reason: str


@dataclass(frozen=True)
class SampleSheetGenerationResult:
    output_file: Path
    rejection_report_file: Path
    review_report_file: Path
    candidate_count: int
    dataset_count: int
    trait_type: str
    effect_type: str
    p_value_type: str
    warnings: tuple[str, ...]
    requires_completion: bool
    datasets: tuple[SampleSheetDatasetStatus, ...]
    rejected_files: tuple[SampleSheetFileRejection, ...]


def _normalise_column_name(value: str) -> str:
    text = unicodedata.normalize("NFKC", str(value)).strip().upper()
    return re.sub(r"[^A-Z0-9]+", "_", text).strip("_")


def _candidate_files(
    input_directory: Path,
    excluded_files: Sequence[Path],
    suffixes: Sequence[str],
) -> list[Path]:
    if not input_directory.is_dir():
        raise ConfigurationError(
            "Summary-statistics input directory does not exist: %s"
            % input_directory
        )
    excluded = {
        path.resolve(strict=False)
        for path in excluded_files
    }
    candidates = [
        path.resolve()
        for path in input_directory.iterdir()
        if path.is_file()
        and not path.name.startswith(".")
        and path.resolve() not in excluded
        and any(path.name.lower().endswith(suffix) for suffix in suffixes)
    ]
    candidates.sort(key=lambda path: (path.name.casefold(), str(path)))
    if not candidates:
        raise ConfigurationError(
            "No supported summary-statistics files were found in %s. "
            "Configured suffixes: %s"
            % (input_directory, ", ".join(suffixes))
        )
    return candidates


def _dataset_id(path: Path, suffixes: Sequence[str]) -> str:
    lowered = path.name.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if lowered.endswith(suffix):
            return path.name[: -len(suffix)].rstrip(".")
    return path.stem


def _delimiter_name(path: Path, policies) -> tuple[str, str]:
    candidates = list(policies.get("input.delimiter_candidates"))
    result = detect_delimiter(
        path,
        candidates=candidates,
        minimum_columns=int(policies.get("input.delimiter_min_columns")),
        maximum_columns=int(policies.get("input.delimiter_max_columns")),
        sample_lines=int(policies.get("input.delimiter_sample_rows")),
        comment_prefix=(
            "##" if bool(policies.get("input.strip_double_hash_lines")) else None
        ),
    )
    for name in candidates:
        if delimiter_character(name) == result.value:
            return name, str(result.value)
    raise ConfigurationError(
        "Delimiter detection returned an unconfigured separator for %s" % path
    )


def _read_header(path: Path, delimiter: str, strip_double_hash: bool) -> list[str]:
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line.strip() or (strip_double_hash and line.startswith("##")):
                continue
            headers = [value.strip() for value in next(csv.reader([line], delimiter=delimiter))]
            if not headers or any(not value for value in headers):
                raise ConfigurationError(
                    "The detected header contains an empty column name: %s" % path
                )
            normalized = [_normalise_column_name(value) for value in headers]
            duplicates = sorted({
                value for value in normalized if normalized.count(value) > 1
            })
            if duplicates:
                raise ConfigurationError(
                    "The header contains ambiguous duplicate column names after "
                    "normalization (%s): %s"
                    % (", ".join(duplicates), path)
                )
            return headers
    raise ConfigurationError("No header line was found in %s" % path)


def _header_index(headers: Sequence[str]) -> dict[str, str]:
    return {_normalise_column_name(value): value for value in headers}


def _matches(index: dict[str, str], aliases: Sequence[str]) -> list[str]:
    return [
        index[normalized]
        for alias in aliases
        if (normalized := _normalise_column_name(alias)) in index
    ]


def _one_match(
    index: dict[str, str], aliases: Sequence[str], label: str,
    path: Path, warnings: list[str],
) -> str | None:
    """Select by canonical alias priority and record alternatives for review."""
    matches = list(dict.fromkeys(_matches(index, aliases)))
    if len(matches) > 1:
        warnings.append(
            "%s: selected %s for %s by configured alias priority; alternatives: %s."
            % (path.name, matches[0], label, ", ".join(matches[1:]))
        )
    return matches[0] if matches else None


def _coordinates(index, aliases, path, warnings):
    chromosome = _one_match(index, aliases.chromosome_column, "chromosome", path, warnings)
    position = _one_match(index, aliases.position_column, "position", path, warnings)
    combined = _one_match(
        index, aliases.chromosome_position_column, "combined-coordinate", path, warnings,
    )
    if chromosome and position:
        if combined:
            warnings.append(
                "%s: separate chromosome/position columns were selected; combined "
                "column %s was not used." % (path.name, combined)
            )
        return chromosome, position, None
    if combined:
        return None, None, combined
    return chromosome, position, None


def _frequency(index, aliases, prefixes, alternative_frequency, effect_allele, path, warnings):
    # Category order retains explicit EAF priority. Within each category the
    # canonical YAML alias order is authoritative; the selected header is literal.
    ordered = (
        aliases.effect_allele_frequency_column
        + alternative_frequency.frequency_columns
        + aliases.maf_frequency_column
        + aliases.ambiguous_frequency_column
        + aliases.reference_frequency_column
    )
    selected = _one_match(index, ordered, "frequency", path, warnings)
    if selected is None:
        prefix_aliases = [
            original for normalized, original in index.items()
            if any(normalized.startswith(prefix) for prefix in prefixes)
        ]
        selected = _one_match(index, prefix_aliases, "frequency", path, warnings)
    if selected is None:
        return None
    normalized = _normalise_column_name(selected)
    if normalized in {_normalise_column_name(v) for v in aliases.maf_frequency_column}:
        warnings.append(
            "%s: %s is MAF-like, not confirmed EAF; allele alignment requires "
            "harmonisation/reference checks." % (path.name, selected)
        )
    elif normalized in {_normalise_column_name(v) for v in aliases.reference_frequency_column}:
        warnings.append(
            "%s: %s is a reference-panel frequency candidate retained in the draft; "
            "its provenance and allele alignment are not established by the header."
            % (path.name, selected)
        )
    elif normalized in {_normalise_column_name(v) for v in aliases.ambiguous_frequency_column}:
        warnings.append(
            "%s: frequency column %s does not state its allele orientation."
            % (path.name, selected)
        )
    elif normalized in {_normalise_column_name(v) for v in alternative_frequency.frequency_columns}:
        if _normalise_column_name(effect_allele) not in {
            _normalise_column_name(v) for v in alternative_frequency.allele_columns
        }:
            warnings.append(
                "%s: ALT frequency %s was retained, but effect allele %s is not "
                "ALT-labelled; confirm allele alignment before harmonisation."
                % (path.name, selected, effect_allele)
            )
    return selected


def _sample_size(index, aliases, path, warnings, on_missing):
    controls = _one_match(index, aliases.control_count_column, "control-count", path, warnings)
    cases = _one_match(index, aliases.case_count_column, "case-count", path, warnings)
    total = _one_match(
        index, aliases.total_sample_size_column + aliases.effective_sample_size_column,
        "total/effective sample size", path, warnings,
    )
    if controls and total:
        warnings.append(
            "%s: explicit control column %s was selected; %s was not used as a control count."
            % (path.name, controls, total)
        )
    elif not controls and cases and total:
        warnings.append(
            "%s: %s cannot be assumed to contain controls when case count %s is present; "
            "that total/effective column was not used for controls." % (path.name, total, cases)
        )
    elif not controls and not cases:
        controls = total
        if total in _matches(index, aliases.effective_sample_size_column):
            warnings.append(
                "%s: %s is an effective sample-size candidate; confirm its meaning."
                % (path.name, total)
            )
    if not controls:
        controls = _one_match(
            index, aliases.fallback_sample_size_column, "fallback sample-size", path, warnings,
        )
        if controls:
            warnings.append(
                "%s: %s was mapped to control_count_column as a fallback; confirm "
                "whether it represents total N, effective N or actual controls before "
                "harmonisation. The header alone does not establish its meaning."
                % (path.name, controls)
            )
        if not controls and on_missing == "fail":
            raise ConfigurationError(
                "%s has no recognized control, total, or effective sample-size column." % path
            )
    return controls, cases


def _missing_mappings(row):
    """Describe the analysis contract without requiring a draft to satisfy it."""
    missing = [
        field for field in ("effect_allele_column", "other_allele_column", "p_value_column")
        if not row.get(field)
    ]
    if not row.get("chromosome_position_column"):
        missing.extend(
            field for field in ("chromosome_column", "position_column") if not row.get(field)
        )
    if not row.get("effect_column") and not row.get("z_score_column"):
        missing.append("effect_column or z_score_column")
    return tuple(missing)


def _row_for_file(path: Path, generator_config, policies) -> tuple[dict[str, Any], list[str]]:
    delimiter_name, delimiter = _delimiter_name(path, policies)
    headers = _read_header(path, delimiter, bool(policies.get("input.strip_double_hash_lines")))
    index = _header_index(headers)
    aliases = generator_config.column_aliases
    warnings: list[str] = []
    chromosome, position, combined = _coordinates(index, aliases, path, warnings)
    row = {
        "config_version": SAMPLE_SHEET_VERSION,
        "dataset_id": _dataset_id(path, generator_config.supported_suffixes),
        "input_file": str(path),
        "trait_type": generator_config.trait_type,
        "effect_type": generator_config.effect_type,
        "p_value_type": generator_config.p_value_type,
        "chromosome_column": chromosome,
        "position_column": position,
        "chromosome_position_column": combined,
        "delimiter": delimiter_name,
    }
    # Equivalent mapping fields share their schema names with configured aliases.
    for field in (
        "variant_id_column", "effect_allele_column", "other_allele_column",
        "standard_error_column", "z_score_column", "p_value_column", "imputation_info_column",
    ):
        row[field] = _one_match(index, getattr(aliases, field), field, path, warnings)
    row["effect_column"] = _one_match(
        index, aliases.beta_effect_column + aliases.odds_ratio_effect_column
        + aliases.ambiguous_effect_column, "effect", path, warnings,
    )
    row["effect_allele_frequency_column"] = _frequency(
        index, aliases, generator_config.effect_frequency_prefixes,
        generator_config.alternative_frequency, row["effect_allele_column"], path, warnings,
    )
    row["control_count_column"], row["case_count_column"] = _sample_size(
        index, aliases, path, warnings, generator_config.on_missing_sample_size,
    )
    # Combined coordinates and Z scores satisfy their respective screening roles.
    screening = dict(row)
    if combined:
        screening.update(chromosome_column=combined, position_column=combined)
    if row["z_score_column"] and not row["effect_column"]:
        screening["effect_column"] = row["z_score_column"]
    mapped = sum(bool(screening.get(field)) for field in generator_config.screening_columns)
    if mapped < generator_config.minimum_mapped_columns:
        raise ConfigurationError(
            "%s has only %d/%d screening mappings; at least %d are required for a draft."
            % (path, mapped, len(generator_config.screening_columns), generator_config.minimum_mapped_columns)
        )
    return row, warnings


def _has_sample_size(row: dict[str, Any]) -> bool:
    return bool(row.get("control_count_column") or row.get("control_count"))


def _source_kind(row: dict[str, Any], source: str) -> str:
    if row.get(_EXTERNAL_SOURCES[source]):
        return "internal"
    if row.get(f"external_{source}_file") and row.get(f"external_{source}_column"):
        return "external"
    return "missing"


def _external_defaults(config, *, config_file, overrides):
    """Validate shared references once; retain templates and exclude actual files."""
    values = config.model_dump()
    requirements = []
    for source in _EXTERNAL_SOURCES:
        fields = [f"external_{source}_{part}" for part in ("file", "column")]
        if any(values[field] for field in fields):
            requirements.extend(RequiredArgument(
                "--" + field.replace("_", "-"),
                "modules.harmonisation.sample_sheet_generator." + field,
                values[field],
            ) for field in fields)
    require_resolved_arguments(requirements)
    defaults, excluded = {}, set()
    for source in _EXTERNAL_SOURCES:
        file_field, column_field = f"external_{source}_file", f"external_{source}_column"
        if not values[file_field]:
            continue
        path = Path(values[file_field]).expanduser()
        # YAML paths are relative to that YAML; CLI paths are relative to cwd.
        if not path.is_absolute() and config_file and file_field not in overrides:
            path = Path(config_file).expanduser().resolve().parent / path
        path = path.resolve()
        try:
            matches = validate_external_resource_spec(path)
            for match in matches:
                with match.open("rb") as handle:
                    handle.read(1)
            excluded.update(matches)
        except (ValueError, OSError) as exc:
            raise ConfigurationError("Invalid argument --%s: %s" % (
                file_field.replace("_", "-"), exc,
            )) from exc
        defaults[file_field] = str(path)
        defaults[column_field] = values[column_field]
    return defaults, excluded


def _write_sample_sheet(
    rows: list[dict[str, Any]],
    output_file: Path,
    *,
    allow_incomplete: bool,
) -> None:
    for row in rows:
        if allow_incomplete:
            HarmonisationSampleSheetRow.validate_generated_draft(row)
        else:
            HarmonisationSampleSheetRow.model_validate(row)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    fields = list(HarmonisationSampleSheetRow.model_fields)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".%s." % output_file.name,
        suffix=".tmp",
        dir=output_file.parent,
        text=True,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    field: "NA" if row.get(field) is None else row.get(field)
                    for field in fields
                })
        if not allow_incomplete:
            load_harmonisation_sample_sheet(temporary)
        os.replace(temporary, output_file)
    finally:
        temporary.unlink(missing_ok=True)


def _write_table_report(
    rows: Sequence[dict[str, Any]],
    fields: Sequence[str],
    output_file: Path,
) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".%s." % output_file.name,
        suffix=".tmp",
        dir=output_file.parent,
        text=True,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=fields,
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output_file)
    finally:
        temporary.unlink(missing_ok=True)


def _write_generation_reports(dataset_statuses, rejected_files, suffixes, rejection_report, review_report):
    """Persist all findings, including when no sample-sheet draft can be written."""
    _write_table_report(
        [{"input_file": str(item.input_file), "reason": item.reason} for item in rejected_files],
        ("input_file", "reason"), rejection_report,
    )
    review_rows = [{
        "dataset_id": item.dataset_id,
        "input_file": str(item.input_file),
        "status": "needs_attention" if item.needs_attention else "ready",
        "missing_mappings": "; ".join(item.missing_mappings),
        "missing_info": item.missing_info,
        "missing_sample_size": item.missing_sample_size,
        "missing_frequency": item.missing_eaf,
        "eaf_source": item.eaf_source,
        "info_source": item.info_source,
        "review_notes": "; ".join(item.warnings),
    } for item in dataset_statuses]
    review_rows.extend({
        "dataset_id": _dataset_id(item.input_file, suffixes),
        "input_file": str(item.input_file), "status": "rejected", "review_notes": item.reason,
    } for item in rejected_files)
    review_rows.sort(key=lambda row: (Path(row["input_file"]).name.casefold(), row["input_file"]))
    _write_table_report(
        review_rows,
        ("dataset_id", "input_file", "status", "missing_mappings", "missing_info",
         "missing_sample_size", "missing_frequency", "eaf_source", "info_source", "review_notes"),
        review_report,
    )


def generate_sample_sheet(
    input_directory: str | Path,
    output_file: str | Path,
    *,
    config_file: str | Path | None = None,
    external_eaf_file: str | Path | None = None,
    external_eaf_column: str | None = None,
    external_info_file: str | Path | None = None,
    external_info_column: str | None = None,
) -> SampleSheetGenerationResult:
    """Validate detected mappings and atomically publish a v2 sample-sheet draft."""
    overrides = {key: str(value) for key, value in {
        "external_eaf_file": external_eaf_file, "external_eaf_column": external_eaf_column,
        "external_info_file": external_info_file, "external_info_column": external_info_column,
    }.items() if value is not None}
    module_config = load_module_configuration("harmonisation", config_file, cli_overrides={
        "sample_sheet_generator." + key: value for key, value in overrides.items()
    })
    generator_config = module_config.sample_sheet_generator
    defaults, reference_files = _external_defaults(
        generator_config, config_file=config_file, overrides=overrides,
    )
    policies = load_policies(module_config.policies)
    source_directory = Path(input_directory).expanduser().resolve()
    destination = Path(output_file).expanduser().resolve()
    rejection_report = destination.with_name(
        destination.name + generator_config.rejection_report_suffix
    )
    review_report = destination.with_name(destination.name + generator_config.review_report_suffix)
    if any(path.resolve() in reference_files for path in (destination, rejection_report, review_report)):
        raise ConfigurationError("Output paths must not overwrite supplied external reference files.")
    files = _candidate_files(
        source_directory,
        (destination, rejection_report, review_report, *reference_files),
        generator_config.supported_suffixes,
    )

    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    dataset_statuses: list[SampleSheetDatasetStatus] = []
    rejected_files: list[SampleSheetFileRejection] = []
    for path in files:
        try:
            row, row_warnings = _row_for_file(path, generator_config, policies)
            for source, internal_field in _EXTERNAL_SOURCES.items():
                if not row.get(internal_field) and defaults.get(f"external_{source}_file"):
                    row.update({f"external_{source}_{part}": defaults[f"external_{source}_{part}"]
                                for part in ("file", "column")})
                    row_warnings.append(
                        "%s: external %s assigned because no internal mapping was detected; "
                        "confirm build, population and source provenance before harmonisation."
                        % (path.name, source.upper())
                    )
            HarmonisationSampleSheetRow.validate_generated_draft(row)
            rows.append(row)
            warnings.extend(row_warnings)
            dataset_statuses.append(SampleSheetDatasetStatus(
                dataset_id=str(row["dataset_id"]),
                input_file=path,
                missing_eaf=_source_kind(row, "eaf") == "missing",
                missing_sample_size=not _has_sample_size(row),
                missing_info=_source_kind(row, "info") == "missing",
                eaf_source=_source_kind(row, "eaf"),
                info_source=_source_kind(row, "info"),
                missing_mappings=_missing_mappings(row),
                has_case_control_columns=bool(row.get("case_count_column") and row.get("control_count_column")),
                has_coordinates=bool(row.get("chromosome_position_column") or
                                     (row.get("chromosome_column") and row.get("position_column"))),
                warnings=tuple(row_warnings),
            ))
        except (ConfigurationError, OSError, ValueError, csv.Error) as exc:
            rejected_files.append(SampleSheetFileRejection(
                input_file=path,
                reason=(
                    _validation_error_text(exc)
                    if isinstance(exc, ValidationError) else str(exc)
                ),
            ))
    _write_generation_reports(
        dataset_statuses, rejected_files, generator_config.supported_suffixes,
        rejection_report, review_report,
    )
    if rejected_files and generator_config.failure_policy == "fail_all":
        raise ConfigurationError(
            "Sample-sheet generation rejected %d/%d files; fail_all requested, "
            "no sample sheet was written. See %s and %s."
            % (len(rejected_files), len(files), rejection_report, review_report)
        )
    if not rows:
        raise ConfigurationError(
            "None of the %d candidate files passed draft screening; no sample sheet "
            "was written. See %s and %s."
            % (len(files), rejection_report, review_report)
        )

    dataset_ids: dict[str, str] = {}
    for row in rows:
        dataset_id = str(row["dataset_id"])
        normalized = dataset_id.casefold()
        if normalized in dataset_ids:
            raise ConfigurationError(
                "Input filenames produce duplicate dataset_id %r: %s and %s"
                % (dataset_id, dataset_ids[normalized], row["input_file"])
            )
        dataset_ids[normalized] = str(row["input_file"])

    requires_completion = any(
        not _has_sample_size(row) or _source_kind(row, "eaf") == "missing" or bool(_missing_mappings(row))
        for row in rows
    )
    _write_sample_sheet(rows, destination, allow_incomplete=requires_completion)
    return SampleSheetGenerationResult(
        output_file=destination,
        rejection_report_file=rejection_report,
        review_report_file=review_report,
        candidate_count=len(files),
        dataset_count=len(rows),
        trait_type=generator_config.trait_type,
        effect_type=generator_config.effect_type,
        p_value_type=generator_config.p_value_type,
        warnings=tuple(warnings),
        requires_completion=requires_completion,
        datasets=tuple(dataset_statuses),
        rejected_files=tuple(rejected_files),
    )


def _generation_summary(result: SampleSheetGenerationResult) -> str:
    """Keep the terminal aggregate-only; individual findings live in TSV reports."""
    datasets = result.datasets
    attention = sum(dataset.needs_attention for dataset in datasets)
    return "\n".join((
        "Sample-sheet generation summary",
        "  Input files       : %d" % result.candidate_count,
        "  Rows written      : %d" % result.dataset_count,
        "  Rejected files    : %d" % len(result.rejected_files),
        "  Ready             : %d" % (result.dataset_count - attention),
        "  Need attention    : %d" % attention,
        "  Missing INFO      : %d" % sum(d.missing_info for d in datasets),
        "  External EAF      : %d" % sum(d.eaf_source == "external" for d in datasets),
        "  External INFO     : %d" % sum(d.info_source == "external" for d in datasets),
        "  No internal INFO  : %d" % sum(d.info_source != "internal" for d in datasets),
        "  Missing sample N  : %d" % sum(d.missing_sample_size for d in datasets),
        "  Missing frequency : %d" % sum(d.missing_eaf for d in datasets),
        "  Frequency mapped  : %d" % sum(not d.missing_eaf for d in datasets),
        "  Case/control cols : %d" % sum(d.has_case_control_columns for d in datasets),
        "  Coordinates mapped: %d" % sum(d.has_coordinates for d in datasets),
        "  Missing mappings  : %d" % sum(bool(d.missing_mappings) for d in datasets),
        "  Mapping review    : %d" % sum(bool(d.warnings) for d in datasets),
        "  Trait type        : %s" % result.trait_type,
        "  Effect type       : %s" % result.effect_type,
        "  P-value type      : %s" % result.p_value_type,
        "  Output            : %s" % result.output_file,
        "  Review report     : %s" % result.review_report_file,
        "  Rejection report  : %s" % result.rejection_report_file,
        "  Mapping counts describe retained rows, not validated study designs or frequency types.",
    ))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas sample_sheet_generator",
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            ("Draft a sample sheet and review reports:", "postgwas sample_sheet_generator",
             ("--input-directory /path/to/raw_data", "--output /path/to/studies.csv")),
            ("Fill missing EAF and INFO mappings:", "postgwas sample_sheet_generator",
             ("--input-directory /path/to/raw_data", "--output /path/to/studies.csv",
              "--external-eaf-file '/references/{build}/eaf_chr{chromosome}.tsv.gz'",
              "--external-eaf-column AF",
              "--external-info-file '/references/{build}/info_chr{chromosome}.tsv.gz'",
              "--external-info-column INFO")),
            ("Use custom mapping aliases:", "postgwas sample_sheet_generator",
             ("--input-directory /path/to/raw_data", "--output /path/to/studies.csv",
              "--run-config /path/to/harmonisation.yaml")),
            ("Export harmonisation configuration:", "postgwas config export",
             ("--module harmonisation", "--style full", "--output harmonisation.yaml")),
        ),
        description=(
            "Generate a version-2 PostGWAS harmonisation sample-sheet draft "
            "from summary-statistics headers."
        )
    )
    parser.add_argument(
        "--input-directory",
        metavar="PATH",
        help="Directory containing summary-statistics files.",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        help="CSV sample sheet to create atomically.",
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        help="Optional harmonisation or full-run YAML overriding header aliases.",
    )
    for source in _EXTERNAL_SOURCES:
        for part in ("file", "column"):
            parser.add_argument(
                f"--external-{source}-{part}", default=argparse.SUPPRESS,
                metavar="PATH" if part == "file" else "NAME",
                help=(f"External {source.upper()} {part} used only when its internal mapping is missing. "
                      f"Requires --external-{source}-{'column' if part == 'file' else 'file'}."
                      + (" Supports {build} and {chromosome}." if part == "file" else "")),
            )
    mark_cli_required_help(parser, ("input_directory", "output"))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        require_resolved_arguments((
            RequiredArgument("--input-directory", None, args.input_directory),
            RequiredArgument("--output", None, args.output),
        ))
        result = generate_sample_sheet(
            args.input_directory,
            args.output,
            config_file=args.run_config,
            **{key: value for key, value in vars(args).items() if key.startswith("external_")},
        )
    except ConfigurationError as exc:
        parser.error(str(exc))
    print(_generation_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
