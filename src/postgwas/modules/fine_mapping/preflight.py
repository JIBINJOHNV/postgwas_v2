"""Mode-aware input, LD-reference, and runtime validation for fine-mapping."""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from postgwas.core.paths import require_nonempty_file, resolve_executable
from postgwas.core.plink import (
    count_plink_samples,
    validate_plink_bed_dimensions,
    validate_plink_files,
)
from postgwas.core.preflight import (
    PreflightFileIdentity,
    capture_preflight_file_identities,
    require_unchanged_preflight_files,
)
from postgwas.modules.fine_mapping.engines.susie.main import (
    SUSIE_REQUIRED_R_PACKAGES,
    resolve_susie_r_runtime,
)
from postgwas.modules.fine_mapping.engines.susie.summary_preparation import (
    REQUIRED_SUMSTAT_COLUMNS,
    normalize_chromosome,
)


# BED/BIM/FAM names and BIM field counts are fixed by
# the public PLINK binary-file contract rather than PostGWAS run policy.
_PLINK_REFERENCE_SUFFIXES = (".bed", ".bim", ".fam")
_PLINK_BIM_FIELDS = 6
_CHROMOSOME_ALIASES = {
    "X": "23",
    "Y": "24",
    "XY": "25",
    "M": "26",
    "MT": "26",
}
# These names are the fixed formatter-to-engine table contract. They are not
# user-selectable column mappings; changing them requires changing both ends of
# the internal interface and its regression tests.
_FINEMAP_SUMSTAT_COLUMNS = (
    "rsid",
    "chromosome",
    "position",
    "allele1",
    "allele2",
    "maf",
    "beta",
    "se",
    "NEF",
)


class FineMappingPreflightError(ValueError):
    """A run-level input or resource failed before locus processing."""


@dataclass(frozen=True)
class FineMappingInputEvidence:
    summary_statistics_file: Path
    summary_statistics_rows: int
    summary_statistics_columns: tuple[str, ...]
    summary_statistics_chromosomes: tuple[str, ...]
    locus_summary_variants: int
    locus_file: Path
    input_loci: int
    eligible_loci: int
    loci_below_lp_threshold: int
    loci_excluded_by_mhc: int
    locus_chromosomes: tuple[str, ...]


@dataclass(frozen=True)
class FineMappingReferenceEvidence:
    prefix: str
    bed_file: Path
    bim_file: Path
    fam_file: Path
    variants: int
    samples: int
    chromosomes: tuple[str, ...]
    bed_bytes: int
    locus_variant_id_matches: int
    locus_variant_ids_requested: int
    coordinate_concordant_variant_matches: int
    loci_with_reference_variants: int
    loci_requested: int


@dataclass(frozen=True)
class FineMappingToolEvidence:
    name: str
    path: str
    version: str


@dataclass(frozen=True)
class FineMappingResourcePreflight:
    """External LD panel and tool evidence available before pipeline inputs."""

    engine: str
    genome_build: str
    reference_prefix: str
    reference_files: tuple[tuple[str, Path], ...]
    reference_identities: tuple[PreflightFileIdentity, ...]
    tools: tuple[FineMappingToolEvidence, ...]
    r_runtime: dict | None


@dataclass(frozen=True)
class FineMappingPreflight:
    engine: str
    genome_build: str
    inputs: FineMappingInputEvidence
    reference: FineMappingReferenceEvidence
    tools: tuple[FineMappingToolEvidence, ...]
    r_runtime: dict | None

    def input_screen_fields(self) -> list[tuple[str, str, object]]:
        engine_label = "SuSiE-RSS" if self.engine == "susie" else "FINEMAP"
        return [
            ("analysis", "Fine-mapping input", engine_label),
            ("info", "Summary-statistics file", self.inputs.summary_statistics_file.name),
            ("count", "Summary-statistic variants", self.inputs.summary_statistics_rows),
            ("count", "Variants inside eligible loci", self.inputs.locus_summary_variants),
            (
                "success",
                "Required summary columns",
                "%d/%d present"
                % (
                    len(self.inputs.summary_statistics_columns),
                    len(self.inputs.summary_statistics_columns),
                ),
            ),
            ("info", "Locus-definition file", self.inputs.locus_file.name),
            ("count", "Input loci", self.inputs.input_loci),
            ("success", "Loci eligible for preparation", self.inputs.eligible_loci),
            (
                "info" if not self.inputs.loci_below_lp_threshold else "loss",
                "Loci below LP threshold",
                self.inputs.loci_below_lp_threshold,
            ),
            (
                "info" if not self.inputs.loci_excluded_by_mhc else "loss",
                "Loci excluded by MHC policy",
                self.inputs.loci_excluded_by_mhc,
            ),
            (
                "genetic",
                "Eligible chromosomes",
                ", ".join(self.inputs.locus_chromosomes),
            ),
            ("success", "File-structure and value validation", "passed"),
        ]

    def resource_screen_fields(self) -> list[tuple[str, str, object]]:
        fields: list[tuple[str, str, object]] = [
            ("analysis", "PLINK LD-reference panel", Path(self.reference.prefix).name),
            ("genetic", "Declared genome build", self.genome_build),
            ("success", "BED/BIM/FAM files", "complete and internally consistent"),
            ("count", "Reference variants", self.reference.variants),
            ("count", "Reference samples", self.reference.samples),
            (
                (
                    "success"
                    if self.reference.locus_variant_id_matches
                    == self.reference.locus_variant_ids_requested
                    else "warning"
                ),
                "Locus variant IDs found in BIM",
                "%s/%s (%.2f%%)"
                % (
                    f"{self.reference.locus_variant_id_matches:,}",
                    f"{self.reference.locus_variant_ids_requested:,}",
                    100.0
                    * self.reference.locus_variant_id_matches
                    / self.reference.locus_variant_ids_requested,
                ),
            ),
            (
                (
                    "success"
                    if self.reference.coordinate_concordant_variant_matches
                    == self.reference.locus_variant_id_matches
                    else "warning"
                ),
                "Matched IDs with concordant coordinates",
                "%s/%s"
                % (
                    f"{self.reference.coordinate_concordant_variant_matches:,}",
                    f"{self.reference.locus_variant_id_matches:,}",
                ),
            ),
            (
                (
                    "success"
                    if self.reference.loci_with_reference_variants
                    == self.reference.loci_requested
                    else "warning"
                ),
                "Eligible loci with reference variants",
                "%s/%s"
                % (
                    f"{self.reference.loci_with_reference_variants:,}",
                    f"{self.reference.loci_requested:,}",
                ),
            ),
            (
                "genetic",
                "Reference chromosomes",
                ", ".join(self.reference.chromosomes),
            ),
        ]
        for tool in self.tools:
            fields.append(("info", tool.name, tool.version))
        if self.r_runtime is not None:
            fields.append(
                (
                    "success",
                    "Required R packages",
                    "%d/%d available"
                    % (len(SUSIE_REQUIRED_R_PACKAGES), len(SUSIE_REQUIRED_R_PACKAGES)),
                )
            )
        fields.append(("success", "Reference and runtime validation", "passed"))
        return fields

    def audit_rows(self) -> list[dict[str, object]]:
        rows = [
            {
                "category": "input",
                "check": "summary_statistics",
                "status": "passed",
                "value": self.inputs.summary_statistics_rows,
                "path": str(self.inputs.summary_statistics_file),
                "detail": "required columns and row values validated",
            },
            {
                "category": "input",
                "check": "locus_definitions",
                "status": "passed",
                "value": self.inputs.eligible_loci,
                "path": str(self.inputs.locus_file),
                "detail": (
                    "input=%d;below_lp=%d;mhc_excluded=%d"
                    % (
                        self.inputs.input_loci,
                        self.inputs.loci_below_lp_threshold,
                        self.inputs.loci_excluded_by_mhc,
                    )
                ),
            },
            {
                "category": "reference",
                "check": "plink_bed_bim_fam",
                "status": "passed",
                "value": self.reference.variants,
                "path": self.reference.prefix,
                "detail": "samples=%d;bed_bytes=%d"
                % (
                    self.reference.samples,
                    self.reference.bed_bytes,
                ),
            },
            {
                "category": "reference",
                "check": "summary_reference_overlap",
                "status": (
                    "passed"
                    if self.reference.loci_with_reference_variants
                    == self.reference.loci_requested
                    else "warning"
                ),
                "value": self.reference.coordinate_concordant_variant_matches,
                "path": self.reference.prefix,
                "detail": "locus_variant_id_matches=%d/%d"
                ";coordinate_concordant_matches=%d"
                ";loci_with_reference_variants=%d/%d"
                % (
                    self.reference.locus_variant_id_matches,
                    self.reference.locus_variant_ids_requested,
                    self.reference.coordinate_concordant_variant_matches,
                    self.reference.loci_with_reference_variants,
                    self.reference.loci_requested,
                ),
            },
        ]
        rows.extend(
            {
                "category": "runtime",
                "check": tool.name,
                "status": "passed",
                "value": tool.version,
                "path": tool.path,
                "detail": "executable resolved before analysis",
            }
            for tool in self.tools
        )
        return rows


def _normalise_chromosome_value(value: object) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.endswith(".0"):
        text = text[:-2]
    text = text.upper()
    return _CHROMOSOME_ALIASES.get(text, text)


def _chromosome_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    try:
        return 0, int(text)
    except ValueError:
        return 1, text


def _merge_intervals(
    intervals: list[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return (
        np.asarray([interval[0] for interval in merged], dtype=np.int64),
        np.asarray([interval[1] for interval in merged], dtype=np.int64),
    )


def _require_executable(value: str | Path, label: str) -> str:
    resolved = resolve_executable(
        value, label, error_type=FineMappingPreflightError,
    )
    if not os.access(resolved, os.X_OK):
        raise FineMappingPreflightError(
            "%s is not executable: %s" % (label, resolved)
        )
    return resolved


def _tool_version(
    executable: str,
    label: str,
    timeout_seconds: float,
    *,
    arguments: tuple[str, ...],
    version_pattern: str,
) -> str:
    """Validate one native informational command and its tool-specific grammar."""
    try:
        result = subprocess.run(
            [executable, *arguments],
            check=False,
            capture_output=True,
            text=True,
            timeout=float(timeout_seconds),
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise FineMappingPreflightError(
            "Cannot execute %s for runtime validation: %s" % (label, exc)
        ) from exc
    output = "\n".join((result.stdout or "", result.stderr or ""))
    if result.returncode != 0 or re.search(r"(?im)^\s*(?:!!\s*)?Error\b", output):
        raise FineMappingPreflightError(
            "%s runtime version probe failed (exit %s): %s %s"
            % (label, result.returncode, executable, " ".join(arguments))
        )
    versions = re.findall(version_pattern, output)
    if len(versions) != 1:
        raise FineMappingPreflightError(
            "%s runtime version probe did not report exactly one recognised version: %s %s"
            % (label, executable, " ".join(arguments))
        )
    return versions[0]


def _summary_file(args) -> Path:
    value = (
        args.susie_input_file
        if args.finemap_method == "susie"
        else args.finemap_in_files
    )
    label = (
        "SuSiE summary-statistics file"
        if args.finemap_method == "susie"
        else "FINEMAP summary-statistics file"
    )
    return require_nonempty_file(
        value, label, error_type=FineMappingPreflightError,
    )


def _summary_invalid_rows(chunk: pd.DataFrame, engine: str) -> pd.Series:
    if engine == "susie":
        chromosome = normalize_chromosome(chunk["CHR"])
        position = pd.to_numeric(chunk["BP"], errors="coerce")
        z_score = pd.to_numeric(chunk["EZ"], errors="coerce")
        sample_size = pd.to_numeric(chunk["NEF"], errors="coerce")
        lp_value = pd.to_numeric(chunk["LP"], errors="coerce")
        snp = chunk["SNP"].astype("string").str.strip()
        reference = chunk["REF"].astype("string").str.strip()
        alternate = chunk["ALT"].astype("string").str.strip()
        return (
            chromosome.isna()
            | chromosome.eq("")
            | position.isna()
            | position.lt(1)
            | position.mod(1).ne(0)
            | snp.isna()
            | snp.eq("")
            | reference.isna()
            | reference.eq("")
            | alternate.isna()
            | alternate.eq("")
            | reference.str.upper().eq(alternate.str.upper())
            | ~np.isfinite(z_score)
            | ~np.isfinite(sample_size)
            | sample_size.le(0)
            | ~np.isfinite(lp_value)
            | lp_value.lt(0)
        )

    chromosome = normalize_chromosome(chunk["chromosome"])
    position = pd.to_numeric(chunk["position"], errors="coerce")
    maf = pd.to_numeric(chunk["maf"], errors="coerce")
    beta = pd.to_numeric(chunk["beta"], errors="coerce")
    standard_error = pd.to_numeric(chunk["se"], errors="coerce")
    sample_size = pd.to_numeric(chunk["NEF"], errors="coerce")
    snp = chunk["rsid"].astype("string").str.strip()
    first = chunk["allele1"].astype("string").str.strip()
    second = chunk["allele2"].astype("string").str.strip()
    return (
        chromosome.isna()
        | chromosome.eq("")
        | position.isna()
        | position.lt(1)
        | position.mod(1).ne(0)
        | snp.isna()
        | snp.eq("")
        | first.isna()
        | first.eq("")
        | second.isna()
        | second.eq("")
        | first.str.upper().eq(second.str.upper())
        | ~np.isfinite(maf)
        | maf.le(0)
        | maf.gt(0.5)
        | ~np.isfinite(beta)
        | ~np.isfinite(standard_error)
        | standard_error.le(0)
        | ~np.isfinite(sample_size)
        | sample_size.le(0)
    )


def _validate_summary_statistics(
    args,
    eligible_intervals: tuple[tuple[str, int, int], ...],
) -> tuple[
    Path,
    int,
    tuple[str, ...],
    tuple[str, ...],
    dict[str, tuple[str, int]],
]:
    source = _summary_file(args)
    engine = args.finemap_method
    required = (
        tuple(REQUIRED_SUMSTAT_COLUMNS)
        if engine == "susie"
        else _FINEMAP_SUMSTAT_COLUMNS
    )
    settings = args.summary_statistics_preparation
    compression = settings["compression"]
    if compression == "none":
        compression = None
    row_count = 0
    chromosomes: set[str] = set()
    relevant_variants: dict[str, tuple[str, int]] = {}
    raw_intervals_by_chromosome: dict[str, list[tuple[int, int]]] = {}
    for chromosome, start, end in eligible_intervals:
        raw_intervals_by_chromosome.setdefault(chromosome, []).append((start, end))
    intervals_by_chromosome = {
        chromosome: _merge_intervals(intervals)
        for chromosome, intervals in raw_intervals_by_chromosome.items()
    }
    observed_columns: tuple[str, ...] | None = None
    try:
        chunks = pd.read_csv(
            source,
            sep=settings["separator"],
            compression=compression,
            chunksize=int(settings["chunk_rows"]),
            keep_default_na=False,
        )
        for chunk_number, chunk in enumerate(chunks, 1):
            if observed_columns is None:
                observed_columns = tuple(str(column) for column in chunk.columns)
                missing = [column for column in required if column not in chunk.columns]
                if missing:
                    raise FineMappingPreflightError(
                        "%s summary statistics are missing required columns: %s"
                        % (engine, ", ".join(missing))
                    )
            invalid = _summary_invalid_rows(chunk, engine)
            if invalid.any():
                raise FineMappingPreflightError(
                    "%s summary statistics contain %d invalid row(s) in input "
                    "chunk %d; coordinates, alleles, statistics, frequencies, "
                    "and sample sizes must satisfy the documented engine contract"
                    % (engine, int(invalid.sum()), chunk_number)
                )
            chromosome_column = "CHR" if engine == "susie" else "chromosome"
            position_column = "BP" if engine == "susie" else "position"
            identifier_column = "SNP" if engine == "susie" else "rsid"
            normalized_chromosomes = normalize_chromosome(
                chunk[chromosome_column]
            ).astype(str)
            positions = pd.to_numeric(chunk[position_column], errors="coerce")
            chromosomes.update(normalized_chromosomes.unique())
            relevant = pd.Series(False, index=chunk.index)
            for chromosome, (starts, ends) in intervals_by_chromosome.items():
                same_chromosome = normalized_chromosomes.eq(chromosome)
                if not same_chromosome.any():
                    continue
                row_indexes = np.flatnonzero(same_chromosome.to_numpy())
                selected_positions = positions.iloc[row_indexes].to_numpy(dtype=np.int64)
                interval_indexes = np.searchsorted(
                    starts, selected_positions, side="right"
                ) - 1
                inside = interval_indexes >= 0
                inside[inside] &= selected_positions[inside] <= ends[
                    interval_indexes[inside]
                ]
                relevant.iloc[row_indexes[inside]] = True
            selected_ids = chunk.loc[relevant, identifier_column].astype(str)
            duplicated = set(selected_ids[selected_ids.duplicated(keep=False)])
            duplicated.update(set(selected_ids).intersection(relevant_variants))
            if duplicated:
                examples = ", ".join(sorted(duplicated)[:5])
                raise FineMappingPreflightError(
                    "Summary statistics repeat variant IDs inside eligible loci; "
                    "examples: %s" % examples
                )
            for index, variant_id in selected_ids.items():
                relevant_variants[variant_id] = (
                    str(normalized_chromosomes.loc[index]),
                    int(positions.loc[index]),
                )
            row_count += len(chunk)
    except FineMappingPreflightError:
        raise
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise FineMappingPreflightError(
            "Cannot parse %s summary statistics %s: %s" % (engine, source, exc)
        ) from exc
    if observed_columns is None or row_count < 1:
        raise FineMappingPreflightError(
            "%s summary-statistics file contains no variants: %s" % (engine, source)
        )
    if not relevant_variants:
        raise FineMappingPreflightError(
            "Summary statistics contain no variants inside any eligible locus boundary"
        )
    return (
        source,
        row_count,
        required,
        tuple(sorted(chromosomes, key=_chromosome_sort_key)),
        relevant_variants,
    )


def _validate_loci(args) -> tuple[
    FineMappingInputEvidence,
    dict[str, tuple[str, int]],
    tuple[tuple[str, int, int], ...],
]:
    source = require_nonempty_file(
        args.locus_file,
        "fine-mapping locus-definition file",
        error_type=FineMappingPreflightError,
    )
    try:
        loci = pd.read_csv(source, sep=None, engine="python")
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise FineMappingPreflightError(
            "Cannot parse fine-mapping locus-definition file %s: %s" % (source, exc)
        ) from exc
    loci.columns = [str(column).strip().upper() for column in loci.columns]
    if loci.columns.duplicated().any():
        raise FineMappingPreflightError(
            "Locus-definition file contains duplicate column names"
        )
    required = (
        ("CHR", "POS", "LP")
        if args.locus_type == "point"
        else ("CHR", "START", "END", "LP")
    )
    missing = [column for column in required if column not in loci.columns]
    if missing:
        raise FineMappingPreflightError(
            "Locus-definition file is missing required %s-mode columns: %s"
            % (args.locus_type, ", ".join(missing))
        )
    if loci.empty:
        raise FineMappingPreflightError(
            "Locus-definition file contains no loci: %s" % source
        )
    chromosomes = normalize_chromosome(loci["CHR"])
    lp_values = pd.to_numeric(loci["LP"], errors="coerce")
    invalid = chromosomes.isna() | chromosomes.eq("") | ~np.isfinite(lp_values)
    invalid |= lp_values.lt(0)
    flank = int(args.window_kb) * int(
        args.fine_mapping_runtime_defaults["bases_per_kilobase"]
    )
    minimum_position = int(args.fine_mapping_runtime_defaults["minimum_position"])
    if args.locus_type == "point":
        positions = pd.to_numeric(loci["POS"], errors="coerce")
        invalid |= positions.isna() | positions.lt(minimum_position) | positions.mod(1).ne(0)
        starts = (positions - flank).clip(lower=minimum_position)
        ends = positions + flank
        duplicate_columns = ["CHR", "POS"]
    else:
        starts_raw = pd.to_numeric(loci["START"], errors="coerce")
        ends_raw = pd.to_numeric(loci["END"], errors="coerce")
        invalid |= (
            starts_raw.isna()
            | ends_raw.isna()
            | starts_raw.lt(minimum_position)
            | ends_raw.lt(starts_raw)
            | starts_raw.mod(1).ne(0)
            | ends_raw.mod(1).ne(0)
        )
        starts = (starts_raw - flank).clip(lower=minimum_position)
        ends = ends_raw + flank
        duplicate_columns = ["CHR", "START", "END"]
    if invalid.any():
        raise FineMappingPreflightError(
            "Locus-definition file contains %d invalid row(s); chromosomes, "
            "coordinates, and LP values must satisfy the documented locus contract"
            % int(invalid.sum())
        )
    normalized_loci = loci.copy()
    normalized_loci["CHR"] = chromosomes
    if normalized_loci.duplicated(duplicate_columns).any():
        raise FineMappingPreflightError(
            "Locus-definition file contains duplicate locus coordinates"
        )
    if "GENOMICLOCUS" in normalized_loci.columns:
        locus_names = normalized_loci["GENOMICLOCUS"].astype("string").str.strip()
        if locus_names.isna().any() or locus_names.eq("").any() or locus_names.duplicated().any():
            raise FineMappingPreflightError(
                "GenomicLocus values must be non-empty and unique"
            )
    passes_lp = lp_values.ge(float(args.lp_threshold))
    excluded_by_mhc = pd.Series(False, index=loci.index)
    if args.finemap_skip_mhc:
        mhc_chromosome = _normalise_chromosome_value(args.finemap_mhc_chromosome)
        excluded_by_mhc = (
            passes_lp
            & chromosomes.eq(mhc_chromosome)
            & starts.lt(int(args.finemap_mhc_end))
            & ends.gt(int(args.finemap_mhc_start))
        )
    eligible = passes_lp & ~excluded_by_mhc
    if not eligible.any():
        raise FineMappingPreflightError(
            "No loci remain after applying LP >= %s and the configured MHC policy"
            % args.lp_threshold
        )
    locus_chromosomes = tuple(
        sorted(set(chromosomes[eligible].astype(str)), key=_chromosome_sort_key)
    )
    eligible_intervals = tuple(
        (
            str(chromosomes.loc[index]),
            int(starts.loc[index]),
            int(ends.loc[index]),
        )
        for index in loci.index[eligible]
    )
    (
        summary_file,
        summary_rows,
        summary_columns,
        summary_chromosomes,
        relevant_variants,
    ) = _validate_summary_statistics(args, eligible_intervals)
    missing_chromosomes = sorted(set(locus_chromosomes) - set(summary_chromosomes))
    if missing_chromosomes:
        raise FineMappingPreflightError(
            "Summary statistics contain no variants for eligible locus chromosome(s): %s"
            % ", ".join(missing_chromosomes)
        )
    return FineMappingInputEvidence(
        summary_statistics_file=summary_file,
        summary_statistics_rows=summary_rows,
        summary_statistics_columns=summary_columns,
        summary_statistics_chromosomes=summary_chromosomes,
        locus_summary_variants=len(relevant_variants),
        locus_file=source,
        input_loci=len(loci),
        eligible_loci=int(eligible.sum()),
        loci_below_lp_threshold=int((~passes_lp).sum()),
        loci_excluded_by_mhc=int(excluded_by_mhc.sum()),
        locus_chromosomes=locus_chromosomes,
    ), relevant_variants, eligible_intervals


def _reference_prefix(value: str | Path) -> str:
    text = str(Path(value).expanduser())
    for suffix in _PLINK_REFERENCE_SUFFIXES:
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    return str(Path(text).resolve())


def _reference_files(prefix_value: str | Path) -> dict[str, Path]:
    prefix = _reference_prefix(prefix_value)
    paths = validate_plink_files(
        prefix, _PLINK_REFERENCE_SUFFIXES, error_type=FineMappingPreflightError,
    )
    return {suffix: paths[suffix.lstrip(".")] for suffix in _PLINK_REFERENCE_SUFFIXES}


def _validated_resource_files(
    resource_preflight: FineMappingResourcePreflight,
) -> dict[str, Path]:
    files = _reference_files(resource_preflight.reference_prefix)
    if tuple(files.items()) != resource_preflight.reference_files:
        raise FineMappingPreflightError(
            "Resolved PLINK LD-reference paths changed after pipeline preflight"
        )
    require_unchanged_preflight_files(
        resource_preflight.reference_identities,
        error_type=FineMappingPreflightError,
        label="PLINK LD-reference files",
    )
    return files


def _validate_plink_reference(
    prefix_value: str | Path,
    required_chromosomes: tuple[str, ...],
    relevant_variants: dict[str, tuple[str, int]],
    eligible_intervals: tuple[tuple[str, int, int], ...],
    *,
    validated_files: dict[str, Path] | None = None,
) -> FineMappingReferenceEvidence:
    prefix = _reference_prefix(prefix_value)
    files = validated_files or _reference_files(prefix)
    samples = count_plink_samples(files[".fam"], error_type=FineMappingPreflightError)

    variants = 0
    chromosomes: set[str] = set()
    matched_variant_ids: set[str] = set()
    concordant_variant_ids: set[str] = set()
    try:
        with files[".bim"].open("r", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                fields = raw.split()
                if len(fields) != _PLINK_BIM_FIELDS:
                    raise FineMappingPreflightError(
                        "PLINK BIM line %d has an invalid field count" % line_number
                    )
                chromosome, variant_id, _, position, first, second = fields
                try:
                    coordinate = int(position)
                    float(fields[2])
                except ValueError as exc:
                    raise FineMappingPreflightError(
                        "PLINK BIM line %d has an invalid coordinate" % line_number
                    ) from exc
                if (
                    coordinate < 1
                    or not variant_id.strip()
                    or not first.strip()
                    or not second.strip()
                ):
                    raise FineMappingPreflightError(
                        "PLINK BIM line %d has an empty ID/allele or invalid position"
                        % line_number
                    )
                normalized_chromosome = _normalise_chromosome_value(chromosome)
                chromosomes.add(normalized_chromosome)
                if variant_id in relevant_variants:
                    if variant_id in matched_variant_ids:
                        raise FineMappingPreflightError(
                            "PLINK BIM repeats locus variant ID %s" % variant_id
                        )
                    matched_variant_ids.add(variant_id)
                    summary_chromosome, summary_position = relevant_variants[variant_id]
                    if (
                        normalized_chromosome == summary_chromosome
                        and coordinate == summary_position
                    ):
                        concordant_variant_ids.add(variant_id)
                variants += 1
    except UnicodeDecodeError as exc:
        raise FineMappingPreflightError(
            "PLINK BIM is not valid UTF-8 text: %s" % files[".bim"]
        ) from exc
    if variants < 1:
        raise FineMappingPreflightError("PLINK BIM contains no variants")
    absent = sorted(set(required_chromosomes) - chromosomes)
    if absent:
        raise FineMappingPreflightError(
            "PLINK LD reference does not contain eligible locus chromosome(s): %s"
            % ", ".join(absent)
        )
    if not matched_variant_ids:
        raise FineMappingPreflightError(
            "No summary-statistic variant IDs inside eligible loci were found in "
            "the PLINK BIM; check build, ancestry, and identifier convention"
        )
    if not concordant_variant_ids:
        raise FineMappingPreflightError(
            "Summary-statistic IDs found in the PLINK BIM have no concordant "
            "chromosome and position; check genome build and identifier convention"
        )

    concordant_positions: dict[str, list[int]] = {}
    for variant_id in concordant_variant_ids:
        chromosome, position = relevant_variants[variant_id]
        concordant_positions.setdefault(chromosome, []).append(position)
    sorted_positions = {
        chromosome: np.asarray(sorted(positions), dtype=np.int64)
        for chromosome, positions in concordant_positions.items()
    }
    loci_with_reference_variants = 0
    for chromosome, start, end in eligible_intervals:
        positions = sorted_positions.get(chromosome)
        if positions is None:
            continue
        first = int(np.searchsorted(positions, start, side="left"))
        if first < len(positions) and positions[first] <= end:
            loci_with_reference_variants += 1

    observed_bed_bytes = validate_plink_bed_dimensions(
        files[".bed"], variants=variants, samples=samples,
        error_type=FineMappingPreflightError,
    )
    return FineMappingReferenceEvidence(
        prefix=prefix,
        bed_file=files[".bed"],
        bim_file=files[".bim"],
        fam_file=files[".fam"],
        variants=variants,
        samples=samples,
        chromosomes=tuple(sorted(chromosomes, key=_chromosome_sort_key)),
        bed_bytes=observed_bed_bytes,
        locus_variant_id_matches=len(matched_variant_ids),
        locus_variant_ids_requested=len(relevant_variants),
        coordinate_concordant_variant_matches=len(concordant_variant_ids),
        loci_with_reference_variants=loci_with_reference_variants,
        loci_requested=len(eligible_intervals),
    )


def _validate_runtime_tools(
    args,
) -> tuple[tuple[FineMappingToolEvidence, ...], dict | None]:
    """Resolve and version-check tools selected by the fine-mapping engine."""
    timeout = float(args.fine_mapping_runtime_defaults["tool_version_timeout_seconds"])
    tools: list[FineMappingToolEvidence] = []
    plink = _require_executable(args.plink, "PLINK executable")
    plink_version = _tool_version(
        plink, "PLINK executable", timeout,
        arguments=("--version",), version_pattern=r"(?m)^(PLINK v[0-9][^\r\n]+)$",
    )
    if args.finemap_method == "finemap" and "PLINK v2" not in plink_version:
        raise FineMappingPreflightError(
            "FINEMAP BGEN generation requires PLINK 2; %s reported: %s"
            % (plink, plink_version)
        )
    args.plink = plink
    args.plink_version = plink_version
    tools.append(FineMappingToolEvidence("PLINK", plink, plink_version))

    r_runtime = None
    if args.finemap_method == "susie":
        try:
            r_runtime = resolve_susie_r_runtime(args.rscript, timeout)
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
            raise FineMappingPreflightError(str(exc)) from exc
        args.rscript = r_runtime["rscript"]
        tools.append(
            FineMappingToolEvidence("R / susieR", args.rscript, r_runtime["version"])
        )
    else:
        # These flags and banner grammars are native CLI protocol invariants,
        # not selectable analysis parameters or executable-path defaults.
        for attribute, label, arguments, version_pattern in (
            ("bgenix", "BGENIX", ("-help",),
             r"(?m)^Welcome to bgenix\s*\n\(version:\s*([0-9][^,\s)]*)[,)]"),
            ("ldstore", "LDstore", ("--help",),
             r"(?m)^\|\s*Welcome to (LDstore v[0-9][^\s|]*)\s*\|$"),
            ("finemap_executable", "FINEMAP", ("--help",),
             r"(?m)^\|\s*Welcome to (FINEMAP v[0-9][^\s|]*)\s*\|$"),
        ):
            executable = _require_executable(getattr(args, attribute), "%s executable" % label)
            setattr(args, attribute, executable)
            tools.append(
                FineMappingToolEvidence(
                    label,
                    executable,
                    _tool_version(
                        executable, "%s executable" % label, timeout,
                        arguments=arguments, version_pattern=version_pattern,
                    ),
                )
            )
    return tuple(tools), r_runtime


def _apply_runtime_evidence(
    args,
    tools: tuple[FineMappingToolEvidence, ...],
    r_runtime: dict | None,
) -> None:
    """Restore validated runtime paths without repeating version subprocesses."""
    by_name = {tool.name: tool for tool in tools}
    expected = (
        {"PLINK", "R / susieR"}
        if args.finemap_method == "susie"
        else {"PLINK", "BGENIX", "LDstore", "FINEMAP"}
    )
    if set(by_name) != expected or len(by_name) != len(tools):
        raise FineMappingPreflightError(
            "Fine-mapping runtime evidence is incomplete or contains duplicate tools"
        )
    args.plink = by_name["PLINK"].path
    args.plink_version = by_name["PLINK"].version
    if args.finemap_method == "susie":
        if r_runtime is None:
            raise FineMappingPreflightError(
                "SuSiE pipeline preflight did not retain its validated R runtime"
            )
        args.rscript = str(r_runtime["rscript"])
        # The inherited process environment may contain large non-path values,
        # credentials, and machine-specific state. Keep it private so generic
        # checkpoint state and input discovery never serialise or scan it.
        args._susie_r_environment = dict(r_runtime["environment"])
        return
    if r_runtime is not None:
        raise FineMappingPreflightError(
            "FINEMAP pipeline preflight contains an unexpected R runtime"
        )
    args.bgenix = by_name["BGENIX"].path
    args.ldstore = by_name["LDstore"].path
    args.finemap_executable = by_name["FINEMAP"].path


def run_fine_mapping_resource_preflight(args) -> FineMappingResourcePreflight:
    """Validate external resources before formatter and locus files exist."""
    files = _reference_files(args.finemap_ld_reference)
    tools, r_runtime = _validate_runtime_tools(args)
    prefix = _reference_prefix(args.finemap_ld_reference)
    args.finemap_ld_reference = prefix
    return FineMappingResourcePreflight(
        engine=args.finemap_method,
        genome_build=args.genome_build,
        reference_prefix=prefix,
        reference_files=tuple(files.items()),
        reference_identities=capture_preflight_file_identities(
            files.values(),
            error_type=FineMappingPreflightError,
            label="PLINK LD-reference files",
        ),
        tools=tools,
        r_runtime=r_runtime,
    )


def run_fine_mapping_preflight(
    args,
    *,
    resource_preflight: FineMappingResourcePreflight | None = None,
) -> FineMappingPreflight:
    """Validate every mode-specific input and resource before locus analysis."""
    inputs, relevant_variants, eligible_intervals = _validate_loci(args)
    if resource_preflight is None:
        reference = _validate_plink_reference(
            args.finemap_ld_reference,
            inputs.locus_chromosomes,
            relevant_variants,
            eligible_intervals,
        )
        tools, r_runtime = _validate_runtime_tools(args)
    else:
        if resource_preflight.engine != args.finemap_method:
            raise FineMappingPreflightError(
                "Fine-mapping engine changed after pipeline resource preflight"
            )
        if resource_preflight.genome_build != args.genome_build:
            raise FineMappingPreflightError(
                "Fine-mapping genome build changed after pipeline resource preflight"
            )
        if resource_preflight.reference_prefix != _reference_prefix(
            args.finemap_ld_reference
        ):
            raise FineMappingPreflightError(
                "Fine-mapping LD-reference prefix changed after pipeline preflight"
            )
        reference = _validate_plink_reference(
            args.finemap_ld_reference,
            inputs.locus_chromosomes,
            relevant_variants,
            eligible_intervals,
            validated_files=_validated_resource_files(resource_preflight),
        )
        tools = resource_preflight.tools
        r_runtime = resource_preflight.r_runtime
    _apply_runtime_evidence(args, tools, r_runtime)
    args.finemap_ld_reference = reference.prefix
    return FineMappingPreflight(
        engine=args.finemap_method,
        genome_build=args.genome_build,
        inputs=inputs,
        reference=reference,
        tools=tools,
        r_runtime=r_runtime,
    )


__all__ = [
    "FineMappingInputEvidence",
    "FineMappingPreflight",
    "FineMappingPreflightError",
    "FineMappingReferenceEvidence",
    "FineMappingResourcePreflight",
    "FineMappingToolEvidence",
    "run_fine_mapping_preflight",
    "run_fine_mapping_resource_preflight",
]
