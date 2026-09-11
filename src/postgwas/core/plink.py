"""Reusable PLINK file checks without study-specific matching policies.

The six FAM fields and SNP-major BED encoding are protocol invariants:
https://www.cog-genomics.org/plink/1.9/formats#bed
Companion suffixes are supplied by the caller's validated configuration.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from postgwas.core.input_validation import record_file_validation, validate_once


# PLINK BED stores four two-bit genotypes per byte after this fixed header.
_SNP_MAJOR_BED_HEADER = bytes((0x6C, 0x1B, 0x01))
_GENOTYPES_PER_BYTE = 4
_FAM_FIELD_COUNT = 6
# The order and whitespace separation are specified by the PLINK BIM protocol;
# callers with existing schema-validated mappings continue to supply those.
PLINK_BIM_COLUMN_ROLES = (
    "chromosome", "variant_id", "genetic_distance", "position", "allele1", "allele2",
)
PLINK_TABLE_DELIMITER_PATTERN = r"\s+"


def validate_plink_files(
    prefix: str | Path,
    required_extensions: Sequence[str],
    *,
    error_type: type[Exception] = ValueError,
    label: str = "PLINK LD reference",
    missing_message: str | None = None,
) -> dict[str, Path]:
    """Validate configured companions, not a genotype bundle when only BIM is used."""
    resolved_prefix = Path(prefix).expanduser().resolve()
    paths = {
        suffix.lstrip("."): Path(str(resolved_prefix) + suffix)
        for suffix in required_extensions
    }
    invalid = []
    for path in paths.values():
        if not path.is_file() or path.stat().st_size <= 0:
            invalid.append(str(path))
            record_file_validation(
                path, "PLINK companion", checks=("regular nonempty file",),
                status="failed", message="Required companion is missing or empty.",
            )
    if invalid:
        message = missing_message or "%s files are missing or empty" % label
        raise error_type("%s: %s" % (message, ", ".join(invalid)))

    def inspect():
        for path in paths.values():
            try:
                with path.open("rb") as handle:
                    handle.read(1)
            except OSError as exc:
                record_file_validation(
                    path, "PLINK companion", checks=("readability",),
                    status="failed", message=str(exc),
                )
                raise error_type("Cannot read %s file %s: %s" % (label, path, exc)) from exc
            record_file_validation(
                path, "PLINK companion", checks=("regular nonempty file", "readability"),
                message="Availability checked; content checks are reported separately.",
            )
        return paths

    return validate_once(
        paths.values(),
        {"validator": "plink_companion_availability", "version": 1},
        inspect,
        error_type=error_type,
    )


def count_plink_samples(
    fam_file: str | Path,
    *,
    allow_blank_rows: bool = False,
    error_type: type[Exception] = ValueError,
) -> int:
    """Count samples while validating the six-field FAM structure in one pass."""
    path = Path(fam_file).expanduser().resolve()

    def inspect():
        count = 0
        blank_rows = 0
        first_blank_row = None
        try:
            with path.open("r", encoding="utf-8", errors="strict") as handle:
                for line_number, raw in enumerate(handle, 1):
                    if not raw.strip():
                        blank_rows += 1
                        if first_blank_row is None:
                            first_blank_row = line_number
                        continue
                    fields = raw.split()
                    if len(fields) != _FAM_FIELD_COUNT:
                        raise ValueError(
                            "PLINK FAM line %d has %d fields; expected %d: %s"
                            % (line_number, len(fields), _FAM_FIELD_COUNT, path)
                        )
                    count += 1
            if not count:
                raise ValueError("PLINK FAM contains no reference samples: %s" % path)
        except (OSError, ValueError) as exc:
            record_file_validation(
                path, "PLINK FAM", checks=("all-row field counts",),
                status="failed", message=str(exc), metrics={"samples_read": count},
            )
            raise error_type(str(exc)) from exc
        record_file_validation(
            path, "PLINK FAM", checks=("nonblank-row field counts", "sample count", "blank-row inventory"),
            metrics={"samples": count, "blank_rows": blank_rows},
        )
        return count, blank_rows, first_blank_row

    count, blank_rows, first_blank_row = validate_once(
        (path,),
        {"validator": "plink_fam_structure", "version": 2},
        inspect,
        error_type=error_type,
    )
    # A caller's blank-row policy does not require another full FAM scan.
    if blank_rows and not allow_blank_rows:
        message = "PLINK FAM line %d has 0 fields; expected %d: %s" % (
            first_blank_row, _FAM_FIELD_COUNT, path,
        )
        record_file_validation(
            path, "PLINK FAM", checks=("caller blank-row policy",),
            status="failed", message=message,
        )
        raise error_type(message)
    record_file_validation(
        path, "PLINK FAM", checks=("caller blank-row policy",),
        metrics={"allow_blank_rows": allow_blank_rows, "blank_rows": blank_rows},
    )
    return count


def validate_plink_bed_dimensions(
    bed_file: str | Path,
    *,
    variants: int,
    samples: int,
    error_type: type[Exception] = ValueError,
) -> int:
    """Validate SNP-major encoding and exact BIM/FAM dimensions, not genotype QC."""
    path = Path(bed_file).expanduser().resolve()

    def inspect():
        try:
            if variants < 1 or samples < 1:
                raise ValueError("PLINK BED validation requires positive BIM/FAM counts.")
            with path.open("rb") as handle:
                header = handle.read(len(_SNP_MAJOR_BED_HEADER))
            if header != _SNP_MAJOR_BED_HEADER:
                raise ValueError(
                    "PLINK BED does not have the required SNP-major binary header: %s" % path
                )
            observed = path.stat().st_size
            bytes_per_variant = (samples + _GENOTYPES_PER_BYTE - 1) // _GENOTYPES_PER_BYTE
            expected = len(_SNP_MAJOR_BED_HEADER) + variants * bytes_per_variant
            if observed != expected:
                raise ValueError(
                    "PLINK BED size is inconsistent with BIM/FAM dimensions: observed=%d, "
                    "expected=%d, variants=%d, samples=%d"
                    % (observed, expected, variants, samples)
                )
        except (OSError, ValueError) as exc:
            record_file_validation(
                path, "PLINK BED", checks=("SNP-major header", "BIM/FAM dimensions"),
                status="failed", message=str(exc),
            )
            raise error_type(str(exc)) from exc
        record_file_validation(
            path, "PLINK BED", checks=("SNP-major header", "BIM/FAM dimensions"),
            metrics={"variants": variants, "samples": samples, "bed_bytes": observed},
            message="Binary layout validated; genotype quality is not inferred.",
        )
        return observed

    return validate_once(
        (path,),
        {"validator": "plink_bed_dimensions", "version": 1, "variants": variants, "samples": samples},
        inspect,
        error_type=error_type,
    )


def validate_plink_bundle_dimensions(
    files: Mapping[str, Path],
    *,
    variants: int,
    samples: int | None = None,
    error_type: type[Exception] = ValueError,
) -> dict[str, int]:
    """Check an actual genotype consumer's bundle using an already scanned BIM.

    ``samples`` may be the result of an existing full FAM structural check, so
    COJO's previously validated sample count does not require a second scan.
    """
    normalized = {suffix.lstrip("."): path for suffix, path in files.items()}
    missing = {"bed", "bim", "fam"} - set(normalized)
    if missing:
        raise error_type(
            "PLINK genotype reference requires configured companions: %s"
            % ", ".join(sorted(missing))
        )

    def inspect():
        sample_count = samples
        if sample_count is None:
            sample_count = count_plink_samples(normalized["fam"], error_type=error_type)
        bed_bytes = validate_plink_bed_dimensions(
            normalized["bed"], variants=variants, samples=sample_count, error_type=error_type,
        )
        return {"variants": variants, "samples": sample_count, "bed_bytes": bed_bytes}

    return validate_once(
        (normalized[role] for role in ("bed", "bim", "fam")),
        {"validator": "plink_bundle_dimensions", "version": 1, "variants": variants, "samples": samples},
        inspect,
        error_type=error_type,
    )
