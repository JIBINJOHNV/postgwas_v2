"""Fast, aggregated validation of chromosome-specific harmonisation resources.

The genome build and observed chromosomes are not known until the dataset-level
checks have run.  This module validates the resulting exact resource set once,
before chromosome partitions are written or worker processes are started.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Mapping, Sequence

from postgwas.core.io.delimiters import open_text, resolve_delimiter
from postgwas.core.processes import run_checked_command
from postgwas.core.vcf import annotation_vcf_fields, declared_vcf_tag_definitions


_RESOURCE_LABELS = {
    "default_eaf_file": "default allele-frequency table",
    "default_comparison_af_file": "comparison allele-frequency VCF",
    "user_eaf_file": "external allele-frequency table",
    "user_info_file": "external INFO table",
    "dbsnp_file": "dbSNP VCF",
    "genome_fasta_file": "source-build FASTA",
    "target_fasta": "target-build FASTA",
    "annot_path": "gene annotation GFF",
    "chain_file": "liftover chain",
}

_ALWAYS_REQUIRED_KEYS = (
    "default_comparison_af_file",
    "dbsnp_file",
    "genome_fasta_file",
    "target_fasta",
    "annot_path",
    "chain_file",
)


class ResourcePreflightError(RuntimeError):
    """One or more resolved chromosome resources cannot support the run."""

    def __init__(self, issues: Sequence[Mapping[str, Any]]):
        self.issues = [dict(issue) for issue in issues]
        super().__init__(format_resource_preflight_issues(self.issues))


def _chromosome_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def _add_issue(
    issues: dict[tuple[str, str, str], set[str]],
    chromosome: str,
    resource: str,
    path: str | Path | None,
    problem: str,
) -> None:
    key = (resource, "not resolved" if path is None else str(path), problem)
    issues.setdefault(key, set()).add(str(chromosome))


def _issue_records(
    issues: Mapping[tuple[str, str, str], set[str]],
) -> list[dict[str, Any]]:
    records = []
    for (resource, path, problem), chromosomes in issues.items():
        records.append({
            "resource": resource,
            "path": path,
            "problem": problem,
            "chromosomes": sorted(chromosomes, key=_chromosome_sort_key),
        })
    return sorted(
        records,
        key=lambda item: (
            _chromosome_sort_key(item["chromosomes"][0]),
            item["resource"],
            item["path"],
        ),
    )


def format_resource_preflight_issues(issues: Sequence[Mapping[str, Any]]) -> str:
    """Render every resource problem without hiding later chromosomes."""
    lines = [
        "Harmonisation resource preflight found %d problem%s:"
        % (len(issues), "" if len(issues) == 1 else "s")
    ]
    for issue in issues:
        chromosomes = ", ".join(str(value) for value in issue["chromosomes"])
        lines.extend([
            "- %s (chromosome%s %s)"
            % (
                issue["resource"],
                "" if len(issue["chromosomes"]) == 1 else "s",
                chromosomes,
            ),
            "  Path: %s" % issue["path"],
            "  Problem: %s" % issue["problem"],
        ])
    return "\n".join(lines)


def _basic_file_problem(path: str | Path) -> str | None:
    candidate = Path(path)
    try:
        if not candidate.is_file():
            return "file does not exist or is not a regular file"
        if candidate.stat().st_size <= 0:
            return "file is empty"
    except OSError as exc:
        return "file metadata cannot be read: %s" % exc
    return None


def _index_path(path: Path, suffixes: Sequence[str]) -> Path | None:
    for suffix in suffixes:
        candidate = Path(str(path) + suffix)
        if _basic_file_problem(candidate) is None:
            return candidate
    return None


def _run(arguments: Sequence[str], purpose: str, logger=None) -> str:
    command = [str(value) for value in arguments]
    if logger is not None and hasattr(logger, "record"):
        logger.record("INPUT", "external_command", purpose=purpose, command=command)
    return run_checked_command(
        command,
        purpose,
        logger=None,
        error_type=RuntimeError,
    )


def _vcf_index_state(
    path: Path, bcftools: str, logger=None,
) -> tuple[dict[str, int | None], str | None]:
    index = _index_path(path, (".tbi", ".csi"))
    if index is None:
        return {}, (
            "BGZF VCF index is missing or empty; create a non-empty %s.tbi or "
            "%s.csi index before rerunning" % (path, path)
        )
    try:
        output = _run(
            (bcftools, "index", "--stats", str(path)),
            "inspect harmonisation reference VCF index",
            logger=logger,
        )
    except RuntimeError as exc:
        return {}, "VCF or its index is unreadable: %s" % exc
    contigs = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3 or not fields[0].strip():
            return {}, "the VCF index produced a malformed contig-statistics line"
        try:
            records = None if fields[2] == "." else int(fields[2])
        except ValueError:
            return {}, "the VCF index produced a non-integer record count"
        contigs[fields[0].strip()] = records
    if not contigs:
        return {}, "the VCF index reports no contigs"
    return contigs, None


def _vcf_info_definitions(
    path: Path, bcftools: str, logger=None,
) -> tuple[dict[str, dict[str, str]], str | None]:
    try:
        header = _run(
            (bcftools, "view", "--header-only", str(path)),
            "inspect harmonisation reference VCF header",
            logger=logger,
        )
    except RuntimeError as exc:
        return {}, "VCF header is unreadable: %s" % exc
    try:
        return declared_vcf_tag_definitions(header, "INFO"), None
    except ValueError as exc:
        return {}, "VCF header is ambiguous: %s" % exc


def _fasta_index_state(path: Path) -> tuple[set[str], str | None]:
    index = Path(str(path) + ".fai")
    problem = _basic_file_problem(index)
    if problem is not None:
        return set(), (
            "FASTA index %s is unavailable: %s; create a valid .fai index "
            "before rerunning"
            % (index, problem)
        )
    contigs = set()
    try:
        with index.open("r", encoding="utf-8", errors="replace") as handle:
            for line_number, raw in enumerate(handle, start=1):
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 5 or not fields[0]:
                    return set(), (
                        "FASTA index has a malformed line at %d" % line_number
                    )
                try:
                    length = int(fields[1])
                    int(fields[2])
                    int(fields[3])
                    int(fields[4])
                except ValueError:
                    return set(), (
                        "FASTA index has non-integer metadata at line %d"
                        % line_number
                    )
                if length <= 0:
                    return set(), (
                        "FASTA index reports a non-positive sequence length at line %d"
                        % line_number
                    )
                contigs.add(fields[0])
    except OSError as exc:
        return set(), "FASTA index cannot be read: %s" % exc
    if not contigs:
        return set(), "FASTA index contains no sequences"
    return contigs, None


def _table_header(
    path: Path,
    configured_delimiter: str,
    policies,
) -> tuple[set[str], str | None]:
    try:
        detected = resolve_delimiter(
            path,
            configured_delimiter,
            candidates=list(policies.get("input.delimiter_candidates")),
            minimum_columns=int(policies.get("input.delimiter_min_columns")),
            maximum_columns=int(policies.get("input.delimiter_max_columns")),
            sample_lines=int(policies.get("input.delimiter_sample_rows")),
        )
        with open_text(path) as handle:
            records = []
            for line in handle:
                if not line.strip() or line.startswith("##"):
                    continue
                records.append(line.rstrip("\r\n"))
                if len(records) == 2:
                    break
        if not records:
            return set(), "table has no readable header"
        if len(records) == 1:
            return set(), "table contains a header but no data rows"
        header = records[0]
        if detected.value == " ":
            columns = header.split()
        else:
            columns = next(csv.reader([header], delimiter=detected.value))
    except (OSError, UnicodeError, ValueError, csv.Error) as exc:
        return set(), "table header cannot be read: %s" % exc
    if not columns or any(not str(column).strip() for column in columns):
        return set(), "table header is empty or contains an unnamed column"
    normalized = [str(column).strip() for column in columns]
    duplicates = sorted({column for column in normalized if normalized.count(column) > 1})
    if duplicates:
        return set(), "table header contains duplicate columns: %s" % ", ".join(
            duplicates
        )
    return set(normalized), None


def _gff_problem(path: Path) -> str | None:
    try:
        with open_text(path) as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line == "##FASTA":
                    break
                if line.startswith("#"):
                    continue
                if len(raw.rstrip("\r\n").split("\t")) != 9:
                    return "first GFF data row does not contain exactly 9 columns"
                return None
    except (OSError, UnicodeError, ValueError) as exc:
        return "GFF file cannot be read: %s" % exc
    return "GFF file contains no data records"


def _chain_mapping_state(
    path: Path,
) -> tuple[dict[str, set[str]], str | None]:
    """Read source-to-target contig names from UCSC chain headers once."""
    mappings: dict[str, set[str]] = {}
    try:
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith("#") or not line.startswith("chain"):
                    continue
                fields = line.split()
                if fields[0] != "chain" or len(fields) != 13:
                    return {}, (
                        "chain header at line %d is not a valid 13-field UCSC "
                        "chain record" % line_number
                    )
                mappings.setdefault(fields[2], set()).add(fields[7])
    except (OSError, UnicodeError, ValueError) as exc:
        return {}, "chain file cannot be read: %s" % exc
    if not mappings:
        return {}, "chain file contains no chain headers"
    return mappings, None


def _mapping_columns(
    mapping: Mapping[str, Any], value_column: Any,
) -> tuple[set[str], str | None]:
    keys = ("chr", "pos", "a1", "a2", "delimiter")
    missing_keys = [key for key in keys if not str(mapping.get(key, "")).strip()]
    if missing_keys:
        return set(), "reference mapping is missing: %s" % ", ".join(missing_keys)
    if value_column is None or not str(value_column).strip():
        return set(), "the configured value column is missing"
    structural = {
        str(mapping["chr"]),
        str(mapping["pos"]),
        str(mapping["a1"]),
        str(mapping["a2"]),
    }
    value = str(value_column)
    if value in structural:
        return set(), (
            "the configured value column %r is also a structural coordinate/allele "
            "column" % value
        )
    return structural | {value}, None


def validate_harmonisation_resource_maps(
    resource_maps: Mapping[str, Mapping[str, Any]],
    *,
    bcftools: str,
    vcf_config: Mapping[str, Any],
    default_eaf_colmap: Mapping[str, Any],
    external_eaf_colmap: Mapping[str, Any],
    external_info_colmap: Mapping[str, Any],
    policies,
    require_default_eaf: bool,
    logger=None,
) -> dict[str, Any]:
    """Validate every resolved resource and return serialisable provenance."""
    issues: dict[tuple[str, str, str], set[str]] = {}
    file_cache: dict[str, str | None] = {}
    vcf_cache: dict[str, tuple[dict[str, int | None], str | None]] = {}
    vcf_header_cache: dict[
        str, tuple[dict[str, dict[str, str]], str | None]
    ] = {}
    fasta_cache: dict[str, tuple[set[str], str | None]] = {}
    table_cache: dict[tuple[str, str, tuple[str, ...]], str | None] = {}
    gff_cache: dict[str, str | None] = {}
    chain_cache: dict[str, tuple[dict[str, set[str]], str | None]] = {}

    expected_info_tags = {
        field.split("/", 1)[1]
        for field in annotation_vcf_fields(
            vcf_config["external_frequency_columns"], "INFO",
        )
    }
    required_keys = list(_ALWAYS_REQUIRED_KEYS)
    if require_default_eaf:
        required_keys.insert(0, "default_eaf_file")

    for chromosome, resources in resource_maps.items():
        active_keys = list(required_keys)
        for optional_key in ("user_eaf_file", "user_info_file"):
            if resources.get(optional_key) is not None:
                active_keys.append(optional_key)

        for key in active_keys:
            label = _RESOURCE_LABELS[key]
            raw_path = resources.get(key)
            if raw_path is None:
                _add_issue(
                    issues, chromosome, label, None,
                    "resource path was not resolved",
                )
                continue
            path = Path(str(raw_path))
            if str(path) not in file_cache:
                file_cache[str(path)] = _basic_file_problem(path)
            problem = file_cache[str(path)]
            if problem is not None:
                _add_issue(issues, chromosome, label, path, problem)
                continue

            if key in ("default_comparison_af_file", "dbsnp_file"):
                if str(path) not in vcf_cache:
                    vcf_cache[str(path)] = _vcf_index_state(
                        path, bcftools, logger=logger,
                    )
                contigs, problem = vcf_cache[str(path)]
                if problem is not None:
                    _add_issue(issues, chromosome, label, path, problem)
                elif str(chromosome) not in contigs:
                    _add_issue(
                        issues, chromosome, label, path,
                        "indexed contigs do not include the exact chromosome label %r; "
                        "bcftools allele annotation requires matching chromosome names"
                        % str(chromosome),
                    )
                elif contigs[str(chromosome)] == 0:
                    _add_issue(
                        issues, chromosome, label, path,
                        "the indexed chromosome contains zero VCF records",
                    )
                if key == "default_comparison_af_file":
                    if str(path) not in vcf_header_cache:
                        vcf_header_cache[str(path)] = _vcf_info_definitions(
                            path, bcftools, logger=logger,
                        )
                    definitions, header_problem = vcf_header_cache[str(path)]
                    if header_problem is not None:
                        _add_issue(issues, chromosome, label, path, header_problem)
                    else:
                        missing = sorted(expected_info_tags - set(definitions))
                        if missing:
                            _add_issue(
                                issues, chromosome, label, path,
                                "VCF header is missing configured INFO tag%s: %s"
                                % ("" if len(missing) == 1 else "s", ", ".join(missing)),
                            )
                        invalid = sorted(
                            tag for tag in expected_info_tags.intersection(definitions)
                            if definitions[tag].get("type") != "Float"
                            or definitions[tag].get("number") != "A"
                        )
                        if invalid:
                            details = ", ".join(
                                "%s(Number=%s,Type=%s)" % (
                                    tag,
                                    definitions[tag].get("number", "missing"),
                                    definitions[tag].get("type", "missing"),
                                )
                                for tag in invalid
                            )
                            _add_issue(
                                issues, chromosome, label, path,
                                "configured frequency INFO tags must be alternate-"
                                "allele values (Number=A, Type=Float): %s"
                                % details,
                            )

            elif key in ("genome_fasta_file", "target_fasta"):
                if str(path) not in fasta_cache:
                    fasta_cache[str(path)] = _fasta_index_state(path)
                contigs, problem = fasta_cache[str(path)]
                if problem is not None:
                    _add_issue(issues, chromosome, label, path, problem)
                elif key == "genome_fasta_file" and str(chromosome) not in contigs:
                    _add_issue(
                        issues, chromosome, label, path,
                        "FASTA index does not include the exact chromosome label %r"
                        % str(chromosome),
                    )

            elif key in ("default_eaf_file", "user_eaf_file", "user_info_file"):
                if key == "default_eaf_file":
                    mapping = default_eaf_colmap
                    value_column = resources.get("default_eaf_reference_column")
                elif key == "user_eaf_file":
                    mapping = external_eaf_colmap
                    value_column = resources.get("user_eaf_column")
                else:
                    mapping = external_info_colmap
                    value_column = resources.get("user_info_column")
                required_columns, mapping_problem = _mapping_columns(
                    mapping, value_column,
                )
                if mapping_problem is not None:
                    _add_issue(issues, chromosome, label, path, mapping_problem)
                    continue
                cache_key = (
                    str(path), str(mapping["delimiter"]),
                    tuple(sorted(required_columns)),
                )
                if cache_key not in table_cache:
                    columns, header_problem = _table_header(
                        path, str(mapping["delimiter"]), policies,
                    )
                    if header_problem is not None:
                        table_cache[cache_key] = header_problem
                    else:
                        missing = sorted(required_columns - columns)
                        table_cache[cache_key] = (
                            None if not missing else
                            "table header is missing configured column%s: %s"
                            % ("" if len(missing) == 1 else "s", ", ".join(missing))
                        )
                problem = table_cache[cache_key]
                if problem is not None:
                    _add_issue(issues, chromosome, label, path, problem)

            elif key == "annot_path":
                if str(path) not in gff_cache:
                    gff_cache[str(path)] = _gff_problem(path)
                problem = gff_cache[str(path)]
                if problem is not None:
                    _add_issue(issues, chromosome, label, path, problem)
            elif key == "chain_file":
                if str(path) not in chain_cache:
                    chain_cache[str(path)] = _chain_mapping_state(path)
                _mappings, problem = chain_cache[str(path)]
                if problem is not None:
                    _add_issue(issues, chromosome, label, path, problem)

        chain_path = resources.get("chain_file")
        target_fasta_path = resources.get("target_fasta")
        if chain_path is not None and target_fasta_path is not None:
            chain_state = chain_cache.get(str(chain_path))
            target_fasta_state = fasta_cache.get(str(target_fasta_path))
            if (
                chain_state is not None
                and chain_state[1] is None
                and target_fasta_state is not None
                and target_fasta_state[1] is None
            ):
                target_labels = chain_state[0].get(str(chromosome), set())
                if not target_labels:
                    _add_issue(
                        issues, chromosome, _RESOURCE_LABELS["chain_file"],
                        chain_path,
                        "chain file has no mapping from source chromosome label %r"
                        % str(chromosome),
                    )
                elif not target_labels.intersection(target_fasta_state[0]):
                    _add_issue(
                        issues, chromosome, _RESOURCE_LABELS["target_fasta"],
                        target_fasta_path,
                        "FASTA index contains none of the target contig labels "
                        "mapped by the chain for source chromosome %r: %s"
                        % (str(chromosome), ", ".join(sorted(target_labels))),
                    )

    records = _issue_records(issues)
    if records:
        raise ResourcePreflightError(records)

    return {
        "status": "passed",
        "chromosomes": sorted(resource_maps, key=_chromosome_sort_key),
        "require_default_eaf": bool(require_default_eaf),
        "unique_resource_files_checked": len(file_cache),
        "indexed_vcfs_checked": len(vcf_cache),
        "frequency_vcf_headers_checked": len(vcf_header_cache),
        "fasta_indexes_checked": len(fasta_cache),
        "table_headers_checked": len(table_cache),
        "structured_text_files_checked": len(gff_cache) + len(chain_cache),
        "resources_by_chromosome": {
            str(chromosome): dict(resources)
            for chromosome, resources in resource_maps.items()
        },
    }


def recheck_preflighted_resource_map(
    chromosome: str,
    resources: Mapping[str, Any],
    *,
    require_default_eaf: bool,
) -> None:
    """Cheap worker-side guard against resources disappearing after preflight."""
    issues: dict[tuple[str, str, str], set[str]] = {}
    active_keys = list(_ALWAYS_REQUIRED_KEYS)
    if require_default_eaf:
        active_keys.insert(0, "default_eaf_file")
    active_keys.extend(
        key for key in ("user_eaf_file", "user_info_file")
        if resources.get(key) is not None
    )
    for key in active_keys:
        path = resources.get(key)
        if path is None:
            _add_issue(
                issues, chromosome, _RESOURCE_LABELS[key], None,
                "resource path was not resolved",
            )
            continue
        problem = _basic_file_problem(path)
        if problem is not None:
            _add_issue(issues, chromosome, _RESOURCE_LABELS[key], path, problem)
            continue
        if key in ("default_comparison_af_file", "dbsnp_file") and _index_path(
            Path(path), (".tbi", ".csi"),
        ) is None:
            _add_issue(
                issues, chromosome, _RESOURCE_LABELS[key], path,
                "VCF index disappeared after dataset preflight",
            )
        if key in ("genome_fasta_file", "target_fasta") and _basic_file_problem(
            Path(str(path) + ".fai")
        ) is not None:
            _add_issue(
                issues, chromosome, _RESOURCE_LABELS[key], path,
                "FASTA index disappeared after dataset preflight",
            )
    records = _issue_records(issues)
    if records:
        raise ResourcePreflightError(records)


__all__ = [
    "ResourcePreflightError",
    "format_resource_preflight_issues",
    "recheck_preflighted_resource_map",
    "validate_harmonisation_resource_maps",
]
