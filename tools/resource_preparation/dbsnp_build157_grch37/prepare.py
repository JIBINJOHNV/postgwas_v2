#!/usr/bin/env python3
"""Prepare current and historical dbSNP Build 157 GRCh37 resources."""

from __future__ import annotations

import argparse
import bz2
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator

from postgwas.core.resource_preparation import (
    JsonLinesRunLogger,
    ResourcePreparationError,
    md5,
    require_executable,
    resumable_download,
    run_command,
    sha256,
)


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ResourceConfig(StrictModel):
    name: str = Field(min_length=1)
    genome_build: Literal["GRCh37.p13"]
    assembly_accession: Literal["GCF_000001405.25"]
    dbsnp_build: PositiveInt
    output_directory: Path


class DownloadConfig(StrictModel):
    name: str = Field(min_length=1)
    filename: str = Field(min_length=1)
    url: str = Field(pattern=r"^https://")
    md5_filename: str | None = None
    md5_url: str | None = Field(default=None, pattern=r"^https://")

    @model_validator(mode="after")
    def validate_checksum_pair(self) -> "DownloadConfig":
        if (self.md5_filename is None) != (self.md5_url is None):
            raise ValueError("md5_filename and md5_url must be provided together")
        for value in (self.filename, self.md5_filename):
            if value is not None and Path(value).name != value:
                raise ValueError("download filenames must not contain directories")
        return self


class InputConfig(StrictModel):
    grch37_vcf: str = Field(min_length=1)
    grch37_vcf_index: str = Field(min_length=1)
    merged_json: str = Field(min_length=1)
    withdrawn_json: str = Field(min_length=1)
    assembly_report: str = Field(min_length=1)


class ChromosomePolicyConfig(StrictModel):
    assembly_report_sequence_column: str = Field(min_length=1)
    assembly_report_role_column: str = Field(min_length=1)
    assembly_report_accession_column: str = Field(min_length=1)
    assembly_report_length_column: str = Field(min_length=1)
    included_sequence_role: str = Field(min_length=1)
    expected_primary_labels: list[str] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_labels(self) -> "ChromosomePolicyConfig":
        if len(self.expected_primary_labels) != len(set(self.expected_primary_labels)):
            raise ValueError("expected primary chromosome labels must be unique")
        return self


class ProcessingConfig(StrictModel):
    curl_executable: str = Field(min_length=1)
    bcftools_executable: str = Field(min_length=1)
    download_workers: PositiveInt
    download_retries: int = Field(ge=0)
    bcftools_threads: PositiveInt
    partial_suffix: str = Field(pattern=r"^\.[^/]+$")
    output_index_type: Literal["csi", "tbi"]
    alias_compression_level: int = Field(ge=1, le=9)


class DirectoryConfig(StrictModel):
    downloads: str = Field(min_length=1)
    derived: str = Field(min_length=1)
    metadata: str = Field(min_length=1)
    logs: str = Field(min_length=1)


class OutputConfig(StrictModel):
    standardized_vcf: str = Field(min_length=1)
    chromosome_map: str = Field(min_length=1)
    primary_targets: str = Field(min_length=1)
    merged_alias_table: str = Field(min_length=1)
    withdrawn_status_table: str = Field(min_length=1)
    manifest: str = Field(min_length=1)
    log: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_names(self) -> "OutputConfig":
        values = list(self.model_dump().values())
        if len(values) != len(set(values)):
            raise ValueError("output filenames must be unique")
        if any(Path(value).name != value for value in values):
            raise ValueError("output filenames must not contain directories")
        return self


class WorkflowConfig(StrictModel):
    resource: ResourceConfig
    downloads: list[DownloadConfig] = Field(min_length=1)
    inputs: InputConfig
    chromosome_policy: ChromosomePolicyConfig
    processing: ProcessingConfig
    directories: DirectoryConfig
    outputs: OutputConfig

    @model_validator(mode="after")
    def validate_cross_references(self) -> "WorkflowConfig":
        names = [item.name for item in self.downloads]
        filenames = [item.filename for item in self.downloads]
        checksum_filenames = [
            item.md5_filename for item in self.downloads if item.md5_filename is not None
        ]
        if len(names) != len(set(names)):
            raise ValueError("download names must be unique")
        if len(filenames + checksum_filenames) != len(set(filenames + checksum_filenames)):
            raise ValueError("download and checksum filenames must be unique")
        unknown = set(self.inputs.model_dump().values()) - set(names)
        if unknown:
            raise ValueError(f"input references unknown downloads: {sorted(unknown)}")
        directories = list(self.directories.model_dump().values())
        if len(directories) != len(set(directories)):
            raise ValueError("directory names must be unique")
        return self


def load_config(path: Path) -> WorkflowConfig:
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ResourcePreparationError(f"configuration must be a YAML mapping: {path}")
    config = WorkflowConfig.model_validate(raw)
    config.resource.output_directory = config.resource.output_directory.expanduser().resolve()
    return config


def build_layout(config: WorkflowConfig) -> dict[str, Path]:
    root = config.resource.output_directory
    result = {"root": root}
    for name, value in config.directories.model_dump().items():
        candidate = Path(value)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ResourcePreparationError(
                f"directory {name!r} must be relative to output_directory"
            )
        result[name] = root / candidate
    return result


def downloads_by_name(config: WorkflowConfig) -> dict[str, DownloadConfig]:
    return {item.name: item for item in config.downloads}


def download_path(
    config: WorkflowConfig, layout: dict[str, Path], logical_name: str,
) -> Path:
    return layout["downloads"] / downloads_by_name(config)[logical_name].filename


def download_all(
    config: WorkflowConfig,
    layout: dict[str, Path],
    curl: str,
    logger: JsonLinesRunLogger,
) -> None:
    jobs: list[tuple[str, Path]] = []
    for item in config.downloads:
        jobs.append((item.url, layout["downloads"] / item.filename))
        if item.md5_url is not None and item.md5_filename is not None:
            jobs.append((item.md5_url, layout["downloads"] / item.md5_filename))
    with ThreadPoolExecutor(max_workers=config.processing.download_workers) as executor:
        futures = [
            executor.submit(
                resumable_download,
                curl=curl,
                url=url,
                destination=destination,
                retries=config.processing.download_retries,
                partial_suffix=config.processing.partial_suffix,
                logger=logger,
            )
            for url, destination in jobs
        ]
        for future in as_completed(futures):
            future.result()


def read_expected_md5(path: Path, expected_filename: str) -> str:
    fields = path.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2:
        raise ResourcePreparationError(f"unexpected upstream MD5 format: {path}")
    checksum, filename = fields
    if filename.lstrip("*") != expected_filename:
        raise ResourcePreparationError(
            f"upstream MD5 names {filename!r}, expected {expected_filename!r}: {path}"
        )
    if len(checksum) != 32 or any(char not in "0123456789abcdefABCDEF" for char in checksum):
        raise ResourcePreparationError(f"invalid MD5 checksum in {path}")
    return checksum.lower()


def validate_downloads(
    config: WorkflowConfig, layout: dict[str, Path], logger: JsonLinesRunLogger,
) -> dict[str, dict[str, Any]]:
    metrics: dict[str, dict[str, Any]] = {}
    for item in config.downloads:
        path = layout["downloads"] / item.filename
        if not path.is_file() or path.stat().st_size == 0:
            raise ResourcePreparationError(f"download is missing or empty: {path}")
        record: dict[str, Any] = {"path": str(path), "bytes": path.stat().st_size}
        if item.md5_filename is not None:
            checksum_path = layout["downloads"] / item.md5_filename
            expected = read_expected_md5(checksum_path, item.filename)
            observed = md5(path)
            if observed != expected:
                raise ResourcePreparationError(
                    f"MD5 mismatch for {path}: observed {observed}, expected {expected}"
                )
            record["upstream_md5"] = expected
        metrics[item.name] = record
        logger.write("download_validated", logical_name=item.name, **record)
    return metrics


def parse_assembly_report(
    path: Path,
    policy: ChromosomePolicyConfig,
) -> list[tuple[str, str, int]]:
    header: list[str] | None = None
    selected: list[tuple[str, str, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith("# Sequence-Name\t"):
                header = line[2:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise ResourcePreparationError(
                    f"assembly report data precede the column header at line {line_number}"
                )
            fields = line.split("\t")
            if len(fields) != len(header):
                raise ResourcePreparationError(
                    f"assembly report line {line_number} has {len(fields)} columns; "
                    f"expected {len(header)}"
                )
            record = dict(zip(header, fields, strict=True))
            if record.get(policy.assembly_report_role_column) != policy.included_sequence_role:
                continue
            label = record.get(policy.assembly_report_sequence_column, "")
            accession = record.get(policy.assembly_report_accession_column, "")
            try:
                length = int(record[policy.assembly_report_length_column])
            except (KeyError, ValueError) as exc:
                raise ResourcePreparationError(
                    f"assembly report line {line_number} has an invalid configured length"
                ) from exc
            if not label or not accession or accession == "na":
                raise ResourcePreparationError(
                    f"assembly report line {line_number} lacks a primary label or accession"
                )
            selected.append((accession, label, length))
    observed = [label for _, label, _ in selected]
    if observed != policy.expected_primary_labels:
        raise ResourcePreparationError(
            f"primary assembly labels differ: observed {observed}, expected "
            f"{policy.expected_primary_labels}"
        )
    if len({accession for accession, _, _ in selected}) != len(selected):
        raise ResourcePreparationError("primary RefSeq accessions are not unique")
    return selected


def write_chromosome_files(
    mapping: list[tuple[str, str, int]],
    layout: dict[str, Path],
    outputs: OutputConfig,
) -> tuple[Path, Path]:
    layout["metadata"].mkdir(parents=True, exist_ok=True)
    mapping_path = layout["metadata"] / outputs.chromosome_map
    targets_path = layout["metadata"] / outputs.primary_targets
    mapping_path.write_text(
        "".join(f"{accession}\t{label}\n" for accession, label, _ in mapping),
        encoding="utf-8",
    )
    targets_path.write_text(
        "".join(f"{label}\t1\t{length}\n" for _, label, length in mapping),
        encoding="utf-8",
    )
    return mapping_path, targets_path


def vcf_index_stats(bcftools: str, path: Path) -> list[dict[str, int | str]]:
    output = run_command([bcftools, "index", "--stats", str(path)], capture=True)
    result: list[dict[str, int | str]] = []
    for line_number, line in enumerate(output.splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 3:
            raise ResourcePreparationError(
                f"unexpected bcftools index --stats line {line_number} for {path}: "
                f"{line!r}"
            )
        try:
            records = int(fields[2])
        except ValueError as exc:
            raise ResourcePreparationError(
                f"invalid indexed record count for {path}: {fields[2]!r}"
            ) from exc
        result.append({"contig": fields[0], "records": records})
    if not result:
        raise ResourcePreparationError(f"VCF index contains no contigs: {path}")
    return result


def source_vcf_index_contigs(bcftools: str, path: Path) -> list[str]:
    """Read contigs from the upstream Tabix index without count metadata."""
    adjacent_tabix = Path(bcftools).with_name("tabix")
    tabix = (
        str(adjacent_tabix)
        if adjacent_tabix.is_file()
        else require_executable("tabix")
    )
    output = run_command([tabix, "--list-chroms", str(path)], capture=True)
    contigs = [line for line in output.splitlines() if line]
    if not contigs:
        raise ResourcePreparationError(f"source VCF index contains no contigs: {path}")
    return contigs


def validate_source_vcf(
    config: WorkflowConfig,
    source: Path,
    bcftools: str,
) -> dict[str, Any]:
    header = run_command([bcftools, "view", "--header-only", str(source)], capture=True)
    expected_build = f"##dbSNP_BUILD_ID={config.resource.dbsnp_build}"
    expected_reference = f"##reference={config.resource.genome_build}"
    if expected_build not in header.splitlines():
        raise ResourcePreparationError(f"source VCF lacks {expected_build}: {source}")
    if expected_reference not in header.splitlines():
        raise ResourcePreparationError(f"source VCF lacks {expected_reference}: {source}")
    column_headers = [line for line in header.splitlines() if line.startswith("#CHROM\t")]
    if len(column_headers) != 1 or len(column_headers[0].split("\t")) != 8:
        raise ResourcePreparationError("dbSNP source VCF must be an eight-column sites VCF")
    return {
        "records": None,
        "record_count_status": "unavailable_in_upstream_tbi",
        "contigs": source_vcf_index_contigs(bcftools, source),
    }


def validate_standardized_vcf(
    path: Path,
    bcftools: str,
    expected_contigs: list[str],
    index_type: str,
) -> dict[str, Any]:
    suffix = ".csi" if index_type == "csi" else ".tbi"
    index_path = Path(str(path) + suffix)
    if not path.is_file() or path.stat().st_size == 0 or not index_path.is_file():
        raise ResourcePreparationError(f"standardized VCF or index is missing: {path}")
    stats = vcf_index_stats(bcftools, path)
    contigs = [str(item["contig"]) for item in stats]
    if contigs != expected_contigs:
        raise ResourcePreparationError(
            f"standardized VCF contigs differ: observed {contigs}, expected "
            f"{expected_contigs}"
        )
    records = sum(int(item["records"]) for item in stats)
    if records < 1:
        raise ResourcePreparationError("standardized VCF contains no records")
    samples = run_command(
        [bcftools, "query", "--list-samples", str(path)], capture=True,
    ).splitlines()
    if samples:
        raise ResourcePreparationError("standardized dbSNP VCF unexpectedly has samples")
    header = run_command(
        [bcftools, "view", "--header-only", str(path)], capture=True,
    )
    if any(line.startswith("##INFO=<") for line in header.splitlines()):
        raise ResourcePreparationError(
            "standardized dbSNP VCF unexpectedly retains INFO definitions"
        )
    alt_query = subprocess.Popen(
        [bcftools, "query", "--format", "%ALT\n", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if alt_query.stdout is None:
        alt_query.kill()
        alt_query.wait()
        raise ResourcePreparationError("could not inspect ALT alleles")
    has_multiallelic = any(b"," in line for line in alt_query.stdout)
    alt_query.stdout.close()
    stderr = alt_query.stderr.read().decode("utf-8", errors="replace") if alt_query.stderr else ""
    return_code = alt_query.wait()
    if return_code != 0:
        raise ResourcePreparationError(f"could not inspect ALT alleles: {stderr}")
    if has_multiallelic:
        raise ResourcePreparationError(
            "standardized dbSNP VCF contains multiallelic records"
        )
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "records": records,
        "contigs": contigs,
        "multiallelic_records": 0,
    }


def convert_vcf(
    *,
    source: Path,
    destination: Path,
    chromosome_map: Path,
    primary_targets: Path,
    expected_contigs: list[str],
    bcftools: str,
    threads: int,
    index_type: str,
) -> dict[str, Any]:
    try:
        return validate_standardized_vcf(
            destination, bcftools, expected_contigs, index_type,
        )
    except (ResourcePreparationError, subprocess.CalledProcessError):
        pass
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.", dir=destination.parent,
    ) as temporary_name:
        temporary = Path(temporary_name)
        prepared = temporary / destination.name
        with (
            tempfile.TemporaryFile() as annotate_error,
            tempfile.TemporaryFile() as norm_error,
            tempfile.TemporaryFile() as view_error,
        ):
            annotate = subprocess.Popen(
                [
                    bcftools,
                    "annotate",
                    "--rename-chrs",
                    str(chromosome_map),
                    "--remove",
                    "INFO",
                    "--threads",
                    str(threads),
                    "--output-type",
                    "u",
                    "--no-version",
                    str(source),
                ],
                stdout=subprocess.PIPE,
                stderr=annotate_error,
            )
            if annotate.stdout is None:
                raise ResourcePreparationError("bcftools annotate stdout pipe is unavailable")
            norm = subprocess.Popen(
                [
                    bcftools,
                    "norm",
                    "--multiallelics",
                    "-any",
                    "--threads",
                    str(threads),
                    "--output-type",
                    "u",
                    "--no-version",
                    "-",
                ],
                stdin=annotate.stdout,
                stdout=subprocess.PIPE,
                stderr=norm_error,
            )
            annotate.stdout.close()
            if norm.stdout is None:
                norm.kill()
                annotate.wait()
                raise ResourcePreparationError("bcftools norm stdout pipe is unavailable")
            view = subprocess.Popen(
                [
                    bcftools,
                    "view",
                    "--targets-file",
                    str(primary_targets),
                    "--threads",
                    str(threads),
                    "--output-type",
                    "z",
                    "--output",
                    str(prepared),
                    "--no-version",
                    "-",
                ],
                stdin=norm.stdout,
                stderr=view_error,
            )
            norm.stdout.close()
            view_return_code = view.wait()
            norm_return_code = norm.wait()
            annotate_return_code = annotate.wait()
            if (
                annotate_return_code != 0
                or norm_return_code != 0
                or view_return_code != 0
            ):
                annotate_error.seek(0)
                norm_error.seek(0)
                view_error.seek(0)
                raise ResourcePreparationError(
                    "bcftools chromosome conversion failed; annotate stderr: "
                    f"{annotate_error.read().decode('utf-8', errors='replace')}; "
                    "norm stderr: "
                    f"{norm_error.read().decode('utf-8', errors='replace')}; "
                    "view stderr: "
                    f"{view_error.read().decode('utf-8', errors='replace')}"
                )
        run_command([
            bcftools,
            "index",
            "--force",
            f"--{index_type}",
            "--threads",
            str(threads),
            str(prepared),
        ])
        metrics = validate_standardized_vcf(
            prepared, bcftools, expected_contigs, index_type,
        )
        suffix = ".csi" if index_type == "csi" else ".tbi"
        os.replace(prepared, destination)
        os.replace(Path(str(prepared) + suffix), Path(str(destination) + suffix))
    metrics["path"] = str(destination)
    metrics["bytes"] = destination.stat().st_size
    return metrics


def build_merged_alias_table(
    source: Path, destination: Path, compression_level: int,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    records = 0
    records_without_current_target = 0
    rows_without_current_target = 0
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", dir=destination.parent, delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        with bz2.open(source, "rt", encoding="utf-8") as source_handle, gzip.open(
            temporary, "wt", encoding="utf-8", compresslevel=compression_level,
        ) as output:
            output.write(
                "PREVIOUS_RSID\tMERGED_INTO_RSID\tRELATION\tMERGE_BUILD\tMERGE_DATE\n"
            )
            for line_number, raw in enumerate(source_handle, 1):
                try:
                    record = json.loads(raw)
                    previous = str(record["refsnp_id"])
                    merged = record["merged_snapshot_data"]
                    targets = [str(value) for value in merged["merged_into"]]
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    raise ResourcePreparationError(
                        f"invalid merged RefSNP JSON at line {line_number}"
                    ) from exc
                if not previous.isdigit() or any(not value.isdigit() for value in targets):
                    raise ResourcePreparationError(
                        f"invalid merged RefSNP identifiers at line {line_number}"
                    )
                if targets:
                    for target in targets:
                        output.write(
                            f"rs{previous}\trs{target}\tmerged_record\t"
                            f"{merged.get('proxy_build_id', '')}\t"
                            f"{merged.get('proxy_time', '')}\n"
                        )
                        rows += 1
                else:
                    output.write(
                        f"rs{previous}\t\tmerged_without_current_target\t"
                        f"{merged.get('proxy_build_id', '')}\t"
                        f"{merged.get('proxy_time', '')}\n"
                    )
                    rows += 1
                    records_without_current_target += 1
                    rows_without_current_target += 1
                for historical in record.get("dbsnp1_merges", []):
                    alias = str(historical.get("merged_rsid", ""))
                    if not alias.isdigit():
                        raise ResourcePreparationError(
                            f"invalid historical merged rsID at line {line_number}"
                        )
                    historical_targets = targets or [""]
                    for target in historical_targets:
                        relation = (
                            "historical_merge_member"
                            if target
                            else "historical_member_without_current_target"
                        )
                        current = f"rs{target}" if target else ""
                        output.write(
                            f"rs{alias}\t{current}\t{relation}\t"
                            f"{historical.get('revision', '')}\t"
                            f"{historical.get('merge_date', '')}\n"
                        )
                        rows += 1
                        if not target:
                            rows_without_current_target += 1
                records += 1
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    if records < 1 or rows < records:
        raise ResourcePreparationError("merged alias table has invalid record counts")
    return {
        "path": str(destination),
        "source_records": records,
        "output_rows": rows,
        "records_without_current_target": records_without_current_target,
        "rows_without_current_target": rows_without_current_target,
    }


def build_withdrawn_status_table(
    source: Path, destination: Path, compression_level: int,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    records = 0
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", dir=destination.parent, delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        with bz2.open(source, "rt", encoding="utf-8") as source_handle, gzip.open(
            temporary, "wt", encoding="utf-8", compresslevel=compression_level,
        ) as output:
            output.write(
                "RSID\tWITHDRAWN_ROOT_RSID\tRELATION\tWITHDRAWN_DATE\t"
                "LAST_UPDATE_BUILD\n"
            )
            for line_number, raw in enumerate(source_handle, 1):
                try:
                    record = json.loads(raw)
                    root = str(record["refsnp_id"])
                    withdrawn = record["withdrawn_snapshot_data"]
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    raise ResourcePreparationError(
                        f"invalid withdrawn RefSNP JSON at line {line_number}"
                    ) from exc
                if not root.isdigit():
                    raise ResourcePreparationError(
                        f"invalid withdrawn RefSNP identifier at line {line_number}"
                    )
                date = withdrawn.get("withdrawn_time", "")
                build = record.get("last_update_build_id", "")
                output.write(f"rs{root}\trs{root}\twithdrawn_record\t{date}\t{build}\n")
                rows += 1
                for historical in record.get("dbsnp1_merges", []):
                    alias = str(historical.get("merged_rsid", ""))
                    if not alias.isdigit():
                        raise ResourcePreparationError(
                            f"invalid historical withdrawn rsID at line {line_number}"
                        )
                    output.write(
                        f"rs{alias}\trs{root}\thistorical_merge_member\t{date}\t{build}\n"
                    )
                    rows += 1
                records += 1
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    if records < 1 or rows < records:
        raise ResourcePreparationError("withdrawn status table has invalid record counts")
    return {"path": str(destination), "source_records": records, "output_rows": rows}


def write_manifest(
    *,
    config: WorkflowConfig,
    config_path: Path,
    layout: dict[str, Path],
    download_metrics: dict[str, dict[str, Any]],
    source_vcf_metrics: dict[str, Any],
    standardized_vcf_metrics: dict[str, Any],
    merged_metrics: dict[str, Any],
    withdrawn_metrics: dict[str, Any],
    bcftools: str,
) -> Path:
    derived_paths = [
        Path(standardized_vcf_metrics["path"]),
        Path(merged_metrics["path"]),
        Path(withdrawn_metrics["path"]),
        layout["metadata"] / config.outputs.chromosome_map,
        layout["metadata"] / config.outputs.primary_targets,
    ]
    output_records = int(standardized_vcf_metrics["records"])
    manifest = {
        "schema_version": "dbsnp_build157_grch37_resource.v1",
        "resource": config.resource.model_dump(mode="json"),
        "resolved_configuration": config.model_dump(mode="json"),
        "configuration_path": str(config_path),
        "configuration_sha256": sha256(config_path),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "bcftools_version": run_command(
            [bcftools, "--version"], capture=True,
        ).splitlines()[0],
        "downloads": download_metrics,
        "chromosome_conversion": {
            "source_contigs": source_vcf_metrics["contigs"],
            "output_contigs": standardized_vcf_metrics["contigs"],
            "source_records": source_vcf_metrics["records"],
            "source_record_count_status": source_vcf_metrics[
                "record_count_status"
            ],
            "output_records": output_records,
            "excluded_non_primary_records": None,
            "excluded_record_count_status": "unavailable_in_upstream_tbi",
            "info_fields": "removed_all",
            "policy": (
                "retain only configured assembly-report role "
                f"{config.chromosome_policy.included_sequence_role!r} and rename "
                "configured source accessions to declared canonical chromosome labels; "
                "split all multiallelic records into biallelic records, and remove "
                "all INFO annotations from the derived VCF"
            ),
            "multiallelic_records": "split_to_biallelic",
        },
        "alias_tables": {
            "merged": merged_metrics,
            "withdrawn": withdrawn_metrics,
        },
        "derived_outputs": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in derived_paths
        ],
        "scientific_sources": [item.url for item in config.downloads],
    }
    manifest_path = layout["metadata"] / config.outputs.manifest
    with manifest_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)
    return manifest_path


def transform(
    config: WorkflowConfig,
    config_path: Path,
    layout: dict[str, Path],
    bcftools: str,
    logger: JsonLinesRunLogger,
    download_metrics: dict[str, dict[str, Any]],
) -> None:
    report = download_path(config, layout, config.inputs.assembly_report)
    mapping = parse_assembly_report(report, config.chromosome_policy)
    chromosome_map, primary_targets = write_chromosome_files(
        mapping, layout, config.outputs,
    )
    source_vcf = download_path(config, layout, config.inputs.grch37_vcf)
    source_vcf_metrics = validate_source_vcf(config, source_vcf, bcftools)
    logger.write("source_vcf_validated", **source_vcf_metrics)
    standardized_vcf = layout["derived"] / config.outputs.standardized_vcf
    standardized_metrics = convert_vcf(
        source=source_vcf,
        destination=standardized_vcf,
        chromosome_map=chromosome_map,
        primary_targets=primary_targets,
        expected_contigs=config.chromosome_policy.expected_primary_labels,
        bcftools=bcftools,
        threads=config.processing.bcftools_threads,
        index_type=config.processing.output_index_type,
    )
    standardized_metrics["excluded_non_primary_records"] = None
    standardized_metrics["excluded_record_count_status"] = (
        "unavailable_in_upstream_tbi"
    )
    standardized_metrics["info_fields"] = "removed_all"
    logger.write("standardized_vcf_completed", **standardized_metrics)
    merged_metrics = build_merged_alias_table(
        download_path(config, layout, config.inputs.merged_json),
        layout["derived"] / config.outputs.merged_alias_table,
        config.processing.alias_compression_level,
    )
    logger.write("merged_alias_table_completed", **merged_metrics)
    withdrawn_metrics = build_withdrawn_status_table(
        download_path(config, layout, config.inputs.withdrawn_json),
        layout["derived"] / config.outputs.withdrawn_status_table,
        config.processing.alias_compression_level,
    )
    logger.write("withdrawn_status_table_completed", **withdrawn_metrics)
    manifest = write_manifest(
        config=config,
        config_path=config_path,
        layout=layout,
        download_metrics=download_metrics,
        source_vcf_metrics=source_vcf_metrics,
        standardized_vcf_metrics=standardized_metrics,
        merged_metrics=merged_metrics,
        withdrawn_metrics=withdrawn_metrics,
        bcftools=bcftools,
    )
    logger.write("manifest_completed", path=str(manifest))


def validate_final(
    config: WorkflowConfig, layout: dict[str, Path], bcftools: str,
) -> dict[str, Any]:
    vcf_metrics = validate_standardized_vcf(
        layout["derived"] / config.outputs.standardized_vcf,
        bcftools,
        config.chromosome_policy.expected_primary_labels,
        config.processing.output_index_type,
    )
    required = [
        layout["derived"] / config.outputs.merged_alias_table,
        layout["derived"] / config.outputs.withdrawn_status_table,
        layout["metadata"] / config.outputs.manifest,
        layout["metadata"] / config.outputs.chromosome_map,
        layout["metadata"] / config.outputs.primary_targets,
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ResourcePreparationError(f"final resource files are missing or empty: {missing}")
    return {"vcf": vcf_metrics, "required_files": [str(path) for path in required]}


def prepare(config_path: Path, stage: str) -> None:
    config = load_config(config_path)
    layout = build_layout(config)
    curl = require_executable(config.processing.curl_executable)
    bcftools = require_executable(config.processing.bcftools_executable)
    layout["logs"].mkdir(parents=True, exist_ok=True)
    logger = JsonLinesRunLogger(layout["logs"] / config.outputs.log)
    logger.record_unfinished_previous_run()
    logger.write(
        "run_started",
        stage=stage,
        config=str(config_path),
        resolved_configuration=config.model_dump(mode="json"),
    )
    try:
        if stage in {"download", "all"}:
            download_all(config, layout, curl, logger)
        if stage in {"transform", "validate", "all"}:
            download_metrics = validate_downloads(config, layout, logger)
        else:
            download_metrics = {}
        if stage in {"transform", "all"}:
            transform(
                config,
                config_path,
                layout,
                bcftools,
                logger,
                download_metrics,
            )
        if stage == "validate":
            logger.write("final_validation_completed", **validate_final(config, layout, bcftools))
        logger.write("run_completed", stage=stage, status="success")
    except BaseException as exc:
        logger.write(
            "run_completed",
            stage=stage,
            status="failure",
            error_type=type(exc).__name__,
            error=str(exc),
        )
        raise


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--stage", required=True, choices=("download", "transform", "validate", "all"),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        prepare(args.config.expanduser().resolve(), args.stage)
    except (
        OSError,
        ResourcePreparationError,
        subprocess.CalledProcessError,
        yaml.YAMLError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
