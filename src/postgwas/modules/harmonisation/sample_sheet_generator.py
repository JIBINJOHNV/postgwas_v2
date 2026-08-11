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

from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.io.delimiters import (
    delimiter_character,
    detect_delimiter,
    open_text,
)

from .policies import load_policies
from .sample_sheet import (
    HarmonisationSampleSheetRow,
    SAMPLE_SHEET_VERSION,
    load_harmonisation_sample_sheet,
)


@dataclass(frozen=True)
class SampleSheetDatasetStatus:
    dataset_id: str
    input_file: Path
    missing_eaf: bool
    missing_sample_size: bool
    missing_info: bool
    warnings: tuple[str, ...]

    @property
    def needs_attention(self) -> bool:
        return bool(
            self.missing_eaf
            or self.missing_sample_size
            or self.missing_info
            or self.warnings
        )


@dataclass(frozen=True)
class SampleSheetGenerationResult:
    output_file: Path
    dataset_count: int
    trait_type: str
    effect_type: str
    p_value_type: str
    warnings: tuple[str, ...]
    requires_completion: bool
    datasets: tuple[SampleSheetDatasetStatus, ...]


def _normalise_column_name(value: str) -> str:
    text = unicodedata.normalize("NFKC", str(value)).strip().upper()
    return re.sub(r"[^A-Z0-9]+", "_", text).strip("_")


def _candidate_files(
    input_directory: Path,
    output_file: Path,
    suffixes: Sequence[str],
) -> list[Path]:
    if not input_directory.is_dir():
        raise ConfigurationError(
            "Summary-statistics input directory does not exist: %s"
            % input_directory
        )
    resolved_output = output_file.resolve(strict=False)
    candidates = [
        path.resolve()
        for path in input_directory.iterdir()
        if path.is_file()
        and not path.name.startswith(".")
        and path.resolve() != resolved_output
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
    index: dict[str, str],
    aliases: Sequence[str],
    label: str,
    path: Path,
) -> str | None:
    matches = list(dict.fromkeys(_matches(index, aliases)))
    if len(matches) > 1:
        raise ConfigurationError(
            "%s has multiple candidate %s columns (%s). Select the intended "
            "column manually in the sample sheet."
            % (path, label, ", ".join(matches))
        )
    return matches[0] if matches else None


def _coordinates(index, aliases, path, warnings):
    chromosome = _one_match(
        index, aliases.chromosome_column, "chromosome", path,
    )
    position = _one_match(index, aliases.position_column, "position", path)
    combined = _one_match(
        index, aliases.chromosome_position_column, "combined-coordinate", path,
    )
    if chromosome and position:
        if combined:
            warnings.append(
                "%s: separate chromosome/position columns were selected; combined "
                "column %s was not used." % (path.name, combined)
            )
        return chromosome, position, None
    if combined and not chromosome and not position:
        return None, None, combined
    raise ConfigurationError(
        "%s must contain chromosome plus position columns, or one combined-coordinate "
        "column; detected chromosome=%r, position=%r, combined=%r."
        % (path, chromosome, position, combined)
    )


def _frequency(
    index,
    aliases,
    prefixes,
    alternative_frequency,
    effect_allele,
    path,
    warnings,
):
    explicit = _one_match(
        index,
        aliases.effect_allele_frequency_column,
        "effect-allele-frequency",
        path,
    )
    prefix_matches = [
        original
        for normalized, original in index.items()
        if any(
            normalized == prefix or normalized.startswith(prefix + "_")
            for prefix in prefixes
        )
    ]
    alternative = _one_match(
        index,
        alternative_frequency.frequency_columns,
        "ALT-allele-frequency",
        path,
    )
    direct_matches = list(dict.fromkeys(
        ([explicit] if explicit else [])
        + prefix_matches
        + ([alternative] if alternative else [])
    ))
    if len(direct_matches) > 1:
        raise ConfigurationError(
            "%s has multiple candidate effect-frequency columns (%s)."
            % (path, ", ".join(direct_matches))
        )
    if direct_matches:
        if alternative:
            alternative_alleles = {
                _normalise_column_name(value)
                for value in alternative_frequency.allele_columns
            }
            if _normalise_column_name(effect_allele) not in alternative_alleles:
                raise ConfigurationError(
                    "%s contains ALT-frequency column %s, but the selected effect-"
                    "allele column is %s rather than an ALT-labelled column. Do not "
                    "use ALT frequency as EAF unless the reported effect is aligned "
                    "to ALT; select the correct columns manually."
                    % (path, alternative, effect_allele)
                )
        return direct_matches[0]

    maf = _one_match(index, aliases.maf_frequency_column, "MAF", path)
    ambiguous = _one_match(
        index, aliases.ambiguous_frequency_column, "ambiguous frequency", path,
    )
    reference = _one_match(
        index, aliases.reference_frequency_column, "reference frequency", path,
    )
    possible = [value for value in (maf, ambiguous) if value]
    if len(possible) > 1:
        raise ConfigurationError(
            "%s has both MAF-like and ambiguous frequency candidates (%s)."
            % (path, ", ".join(possible))
        )
    if maf:
        warnings.append(
            "%s: %s is MAF-like, not confirmed EAF. Harmonisation must confirm "
            "its allele alignment against the configured reference panel."
            % (path.name, maf)
        )
        return maf
    if ambiguous:
        warnings.append(
            "%s: frequency column %s does not state its allele orientation. Review "
            "the generated mapping before harmonisation." % (path.name, ambiguous)
        )
        return ambiguous
    if reference:
        warnings.append(
            "%s: reference-panel frequency column %s was not used as study EAF. "
            "Set effect_allele_frequency_column to a valid study-frequency column, "
            "or set both external_eaf_file and external_eaf_column before "
            "harmonisation."
            % (path.name, reference)
        )
    return None


def _effect_column(index, aliases, z_score, path):
    beta = _one_match(index, aliases.beta_effect_column, "beta", path)
    odds_ratio = _one_match(
        index, aliases.odds_ratio_effect_column, "odds-ratio", path,
    )
    ambiguous = _one_match(
        index, aliases.ambiguous_effect_column, "ambiguous effect", path,
    )
    categories = [value for value in (beta, odds_ratio, ambiguous) if value]
    if len(categories) > 1:
        raise ConfigurationError(
            "%s has multiple candidate effect columns (%s). Select the intended "
            "effect column manually." % (path, ", ".join(categories))
        )
    if categories:
        return categories[0]
    if z_score:
        return None
    raise ConfigurationError(
        "%s has neither a recognized effect column nor a Z-score column." % path
    )


def _sample_size(index, aliases, path, warnings, on_missing):
    controls = _one_match(
        index, aliases.control_count_column, "control-count", path,
    )
    cases = _one_match(index, aliases.case_count_column, "case-count", path)
    total = _one_match(
        index, aliases.total_sample_size_column, "total-sample-size", path,
    )
    effective = _one_match(
        index, aliases.effective_sample_size_column, "effective-sample-size", path,
    )
    invalid = _one_match(
        index, aliases.invalid_sample_size_column, "unsafe sample-size", path,
    )
    if controls:
        ignored = [value for value in (total, effective) if value]
        if ignored:
            warnings.append(
                "%s: explicit control column %s was selected; %s was not used as a "
                "control count." % (path.name, controls, ", ".join(ignored))
            )
    else:
        alternatives = [value for value in (total, effective) if value]
        if len(alternatives) > 1:
            raise ConfigurationError(
                "%s has both total and effective sample-size columns (%s); select "
                "the scientifically appropriate one manually."
                % (path, ", ".join(alternatives))
            )
        if cases and alternatives:
            raise ConfigurationError(
                "%s has case count %s but no explicit control count. Column %s cannot "
                "be assumed to contain controls." % (path, cases, alternatives[0])
            )
        controls = alternatives[0] if alternatives else None
        if effective:
            warnings.append(
                "%s: %s is treated as the lone effective/total sample-size input; "
                "review this mapping before harmonisation." % (path.name, effective)
            )
    if cases and not controls:
        raise ConfigurationError(
            "%s has a case-count column but no usable control-count column." % path
        )
    if not controls:
        if invalid:
            raise ConfigurationError(
                "%s contains %s, but a meta-analysis weight is not a sample size."
                % (path, invalid)
            )
        if on_missing == "fail":
            raise ConfigurationError(
                "%s has no recognized control, total, or effective sample-size "
                "column." % path
            )
        return None, None
    return controls, cases


def _row_for_file(path: Path, generator_config, policies) -> tuple[dict[str, Any], list[str]]:
    delimiter_name, delimiter = _delimiter_name(path, policies)
    headers = _read_header(
        path,
        delimiter,
        bool(policies.get("input.strip_double_hash_lines")),
    )
    index = _header_index(headers)
    aliases = generator_config.column_aliases
    warnings: list[str] = []

    chromosome, position, combined = _coordinates(index, aliases, path, warnings)
    effect_allele = _one_match(
        index, aliases.effect_allele_column, "effect-allele", path,
    )
    other_allele = _one_match(
        index, aliases.other_allele_column, "other-allele", path,
    )
    if not effect_allele or not other_allele:
        raise ConfigurationError(
            "%s must contain recognized effect and other allele columns; detected "
            "effect=%r, other=%r." % (path, effect_allele, other_allele)
        )

    variant_id = _one_match(
        index, aliases.variant_id_column, "variant-identifier", path,
    )
    frequency = _frequency(
        index,
        aliases,
        generator_config.effect_frequency_prefixes,
        generator_config.alternative_frequency,
        effect_allele,
        path,
        warnings,
    )
    standard_error = _one_match(
        index, aliases.standard_error_column, "standard-error", path,
    )
    z_score = _one_match(index, aliases.z_score_column, "Z-score", path)
    effect = _effect_column(index, aliases, z_score, path)
    p_value = _one_match(index, aliases.p_value_column, "p-value", path)
    if not p_value:
        raise ConfigurationError("%s has no recognized p-value column." % path)
    controls, cases = _sample_size(
        index,
        aliases,
        path,
        warnings,
        generator_config.on_missing_sample_size,
    )
    info = _one_match(
        index, aliases.imputation_info_column, "imputation-quality", path,
    )

    row = {
        "config_version": SAMPLE_SHEET_VERSION,
        "dataset_id": _dataset_id(path, generator_config.supported_suffixes),
        "input_file": str(path),
        # Headers select sample-size fields but cannot establish study design.
        # Full-dataset harmonisation makes that scientific decision.
        "trait_type": generator_config.trait_type,
        "chromosome_column": chromosome,
        "position_column": position,
        "chromosome_position_column": combined,
        "variant_id_column": variant_id,
        "effect_allele_column": effect_allele,
        "other_allele_column": other_allele,
        "effect_allele_frequency_column": frequency,
        "effect_column": effect,
        # Header aliases select columns only. Full-dataset harmonisation decides
        # both statistical scales unless the user edits the generated sheet.
        "effect_type": generator_config.effect_type,
        "standard_error_column": standard_error,
        "z_score_column": z_score,
        "p_value_column": p_value,
        "p_value_type": generator_config.p_value_type,
        "control_count_column": controls,
        "case_count_column": cases,
        "control_count": None,
        "case_count": None,
        "imputation_info_column": info,
        "external_info_file": None,
        "external_info_column": None,
        "external_eaf_file": None,
        "external_eaf_column": None,
        "delimiter": delimiter_name,
    }
    return row, warnings


def _has_sample_size(row: dict[str, Any]) -> bool:
    return bool(row.get("control_count_column") or row.get("control_count"))


def _has_eaf(row: dict[str, Any]) -> bool:
    internal = bool(row.get("effect_allele_frequency_column"))
    external = bool(row.get("external_eaf_file") and row.get("external_eaf_column"))
    return internal or external


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


def generate_sample_sheet(
    input_directory: str | Path,
    output_file: str | Path,
    *,
    config_file: str | Path | None = None,
) -> SampleSheetGenerationResult:
    """Validate detected mappings and atomically publish a v2 sample-sheet draft."""
    module_config = load_module_configuration("harmonisation", config_file)
    generator_config = module_config.sample_sheet_generator
    policies = load_policies(module_config.policies)
    source_directory = Path(input_directory).expanduser().resolve()
    destination = Path(output_file).expanduser().resolve()
    files = _candidate_files(
        source_directory,
        destination,
        generator_config.supported_suffixes,
    )

    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    dataset_statuses: list[SampleSheetDatasetStatus] = []
    problems: list[str] = []
    for path in files:
        try:
            row, row_warnings = _row_for_file(path, generator_config, policies)
            rows.append(row)
            warnings.extend(row_warnings)
            dataset_statuses.append(SampleSheetDatasetStatus(
                dataset_id=str(row["dataset_id"]),
                input_file=path,
                missing_eaf=not _has_eaf(row),
                missing_sample_size=not _has_sample_size(row),
                missing_info=row.get("imputation_info_column") is None,
                warnings=tuple(row_warnings),
            ))
        except (ConfigurationError, OSError, ValueError, csv.Error) as exc:
            problems.append("%s: %s" % (path.name, exc))
    if problems:
        raise ConfigurationError(
            "Sample-sheet generation failed for %d file(s); no output was written:\n- %s"
            % (len(problems), "\n- ".join(problems))
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
        not _has_sample_size(row) or not _has_eaf(row)
        for row in rows
    )
    _write_sample_sheet(
        rows,
        destination,
        allow_incomplete=requires_completion,
    )
    return SampleSheetGenerationResult(
        output_file=destination,
        dataset_count=len(rows),
        trait_type=generator_config.trait_type,
        effect_type=generator_config.effect_type,
        p_value_type=generator_config.p_value_type,
        warnings=tuple(warnings),
        requires_completion=requires_completion,
        datasets=tuple(dataset_statuses),
    )


def _generation_summary(result: SampleSheetGenerationResult) -> str:
    attention = [dataset for dataset in result.datasets if dataset.needs_attention]
    ready = result.dataset_count - len(attention)
    lines = [
        "Sample-sheet generation summary",
        "  Input files       : %d" % result.dataset_count,
        "  Rows written      : %d" % result.dataset_count,
        "  Ready             : %d" % ready,
        "  Need attention    : %d" % len(attention),
        "  Output            : %s" % result.output_file,
        "  Trait type        : %s" % result.trait_type,
        "  Effect type       : %s" % result.effect_type,
        "  P-value type      : %s" % result.p_value_type,
    ]
    if not attention:
        return "\n".join(lines)

    lines.extend(("", "Datasets needing attention"))
    for dataset in attention:
        issues: list[str] = []
        if dataset.missing_eaf:
            issues.append("allele frequency")
        if dataset.missing_sample_size:
            issues.append("sample size")
        if dataset.missing_info:
            issues.append("imputation INFO")
        if dataset.warnings:
            issues.append("mapping review")
        lines.append(
            "  %s (%s): %s"
            % (dataset.dataset_id, dataset.input_file.name, ", ".join(issues))
        )

    if any(dataset.missing_eaf for dataset in attention):
        lines.extend((
            "",
            "How to provide allele frequency",
            "  Internal column : set effect_allele_frequency_column to the exact "
            "study EAF or MAF header.",
            "  External file   : set external_eaf_file AND external_eaf_column.",
            "  Required        : complete one of these alternatives before running "
            "harmonisation; do not invent a fixed frequency.",
        ))
    if any(dataset.missing_sample_size for dataset in attention):
        lines.extend((
            "",
            "How to complete sample size",
            "  Quantitative study",
            "    Set control_count_column to the total-N column, or set "
            "control_count to the fixed total N.",
            "  Case-control study",
            "    Set control_count_column/control_count AND "
            "case_count_column/case_count.",
        ))
    if any(dataset.missing_info for dataset in attention):
        lines.extend((
            "",
            "How to provide imputation INFO",
            "  Internal column  : set imputation_info_column",
            "  External file    : set external_info_file and external_info_column",
            "  Fixed run value  : use harmonisation --fixed-info VALUE",
        ))

    review_notes = [
        (dataset, warning)
        for dataset in attention
        for warning in dataset.warnings
    ]
    if review_notes:
        lines.extend(("", "Other mapping review notes"))
        for dataset, warning in review_notes:
            prefix = "%s: " % dataset.input_file.name
            detail = warning.removeprefix(prefix)
            lines.append("  %s: %s" % (dataset.dataset_id, detail))
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a version-2 PostGWAS harmonisation sample-sheet draft "
            "from summary-statistics headers."
        )
    )
    parser.add_argument(
        "--input-directory",
        required=True,
        metavar="PATH",
        help="Directory containing summary-statistics files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="PATH",
        help="CSV sample sheet to create atomically.",
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        help="Optional harmonisation or full-run YAML overriding header aliases.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = generate_sample_sheet(
            args.input_directory,
            args.output,
            config_file=args.run_config,
        )
    except ConfigurationError as exc:
        parser.error(str(exc))
    print(_generation_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
