#!/usr/bin/env python3
"""Prepare population-specific 1000 Genomes Phase 3 genotype VCFs.

The workflow is resumable and writes only validated, indexed files to their
final locations. All release-specific policy is supplied by a validated YAML
configuration rather than inferred from filenames or genomic coordinates.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator

from postgwas.core.resource_preparation import (
    JsonLinesRunLogger,
    require_executable,
    resumable_download,
    run_command,
    sha256,
)


# Three concurrent jobs times four bcftools threads uses the 12 logical CPUs
# available on the preparation host without unbounded oversubscription.
CONTIG_PROCESSING_WORKERS = 3
CONCATENATION_WORKERS = 3


class PreparationError(RuntimeError):
    """Raised when a resource cannot be prepared without ambiguity or data loss."""


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ResourceConfig(StrictModel):
    name: str = Field(min_length=1)
    genome_build: Literal["GRCh37"]
    release: str = Field(min_length=1)
    source_url: str = Field(pattern=r"^https://")
    output_directory: Path


class PanelConfig(StrictModel):
    filename: str = Field(min_length=1)
    sample_column: str = Field(min_length=1)
    grouping_column: str = Field(min_length=1)
    groups: list[str] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_groups(self) -> "PanelConfig":
        if len(self.groups) != len(set(self.groups)):
            raise ValueError("panel groups must be unique")
        return self


class ContigConfig(StrictModel):
    label: str = Field(min_length=1)
    filename: str = Field(min_length=1)


class OutputSetConfig(StrictModel):
    name: str = Field(min_length=1)
    contigs: list[str] = Field(min_length=1)
    filename_template: str = Field(min_length=1)


class ProcessingConfig(StrictModel):
    bcftools_executable: str = Field(min_length=1)
    curl_executable: str = Field(min_length=1)
    bcftools_threads: PositiveInt
    download_workers: PositiveInt
    download_retries: int = Field(ge=0)
    index_suffix: Literal[".tbi", ".csi"]
    output_index_type: Literal["csi", "tbi"]
    recalculate_info_tags: list[Literal["AC", "AN", "AF", "NS"]]
    remove_per_contig_outputs_after_success: bool
    keep_downloads_after_success: bool


class DirectoryConfig(StrictModel):
    downloads: str = Field(min_length=1)
    metadata: str = Field(min_length=1)
    per_contig: str = Field(min_length=1)
    final: str = Field(min_length=1)
    logs: str = Field(min_length=1)


class NamingConfig(StrictModel):
    partial_suffix: str = Field(pattern=r"^\.[^/]+$")
    group_sample_list_template: str = Field(min_length=1)
    per_contig_vcf_template: str = Field(min_length=1)
    split_groups_filename: str = Field(min_length=1)
    split_prepared_vcf_template: str = Field(min_length=1)
    concat_input_list_filename: str = Field(min_length=1)
    log_filename: str = Field(min_length=1)
    manifest_filename: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_templates(self) -> "NamingConfig":
        template_fields = (
            "group_sample_list_template",
            "per_contig_vcf_template",
            "split_prepared_vcf_template",
        )
        for field in template_fields:
            if "{group}" not in getattr(self, field):
                raise ValueError(f"naming.{field} must contain '{{group}}'")
        for field, value in self.model_dump().items():
            if field != "partial_suffix" and Path(value).name != value:
                raise ValueError(f"naming.{field} must be a filename, not a path")
        return self


class WorkflowConfig(StrictModel):
    resource: ResourceConfig
    panel: PanelConfig
    contigs: list[ContigConfig] = Field(min_length=1)
    output_sets: list[OutputSetConfig] = Field(min_length=1)
    processing: ProcessingConfig
    directories: DirectoryConfig
    naming: NamingConfig

    @model_validator(mode="after")
    def validate_cross_references(self) -> "WorkflowConfig":
        labels = [item.label for item in self.contigs]
        if len(labels) != len(set(labels)):
            raise ValueError("contig labels must be unique")
        names = [item.name for item in self.output_sets]
        if len(names) != len(set(names)):
            raise ValueError("output-set names must be unique")
        known = set(labels)
        used: dict[str, str] = {}
        for output_set in self.output_sets:
            if "{group}" not in output_set.filename_template:
                raise ValueError(
                    f"output set {output_set.name!r} filename_template must contain "
                    "'{group}'"
                )
            unknown = set(output_set.contigs) - known
            if unknown:
                raise ValueError(
                    f"output set {output_set.name!r} has unknown contigs: "
                    f"{sorted(unknown)}"
                )
            for label in output_set.contigs:
                if label in used:
                    raise ValueError(
                        f"contig {label!r} occurs in output sets {used[label]!r} and "
                        f"{output_set.name!r}"
                    )
                used[label] = output_set.name
        missing = known - set(used)
        if missing:
            raise ValueError(f"contigs are absent from output sets: {sorted(missing)}")
        directories = self.directories.model_dump().values()
        if len(set(directories)) != len(tuple(directories)):
            raise ValueError("directory names must be unique")
        return self


def load_config(path: Path) -> WorkflowConfig:
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise PreparationError(f"configuration must be a YAML mapping: {path}")
    config = WorkflowConfig.model_validate(raw)
    config.resource.output_directory = config.resource.output_directory.expanduser().resolve()
    return config


def ensure_relative_filename(value: str, label: str) -> None:
    if Path(value).name != value:
        raise PreparationError(f"{label} must be a filename, not a path: {value!r}")


def layout(config: WorkflowConfig) -> dict[str, Path]:
    root = config.resource.output_directory
    result = {"root": root}
    for key, value in config.directories.model_dump().items():
        candidate = Path(value)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise PreparationError(f"directory {key!r} must be relative to output_directory")
        result[key] = root / candidate
    return result


def source_url(config: WorkflowConfig, filename: str) -> str:
    ensure_relative_filename(filename, "source filename")
    return f"{config.resource.source_url.rstrip('/')}/{filename}"


def download_sources(
    config: WorkflowConfig,
    paths: dict[str, Path],
    curl: str,
    logger: JsonLinesRunLogger,
) -> None:
    requested = [config.panel.filename]
    for contig in config.contigs:
        requested.extend((contig.filename, contig.filename + config.processing.index_suffix))
    jobs = [
        (source_url(config, filename), paths["downloads"] / filename)
        for filename in requested
    ]
    with ThreadPoolExecutor(max_workers=config.processing.download_workers) as executor:
        futures = {
            executor.submit(
                resumable_download,
                curl=curl,
                url=url,
                destination=destination,
                retries=config.processing.download_retries,
                partial_suffix=config.naming.partial_suffix,
                logger=logger,
            ): destination
            for url, destination in jobs
        }
        for future in as_completed(futures):
            future.result()


def read_panel(config: WorkflowConfig, panel_path: Path) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {group: [] for group in config.panel.groups}
    seen: set[str] = set()
    with panel_path.open("r", encoding="utf-8") as handle:
        header_line = handle.readline()
        if not header_line:
            raise PreparationError(f"sample panel is empty: {panel_path}")
        header = header_line.rstrip("\r\n").split("\t")
        try:
            sample_index = header.index(config.panel.sample_column)
            group_index = header.index(config.panel.grouping_column)
        except ValueError as exc:
            raise PreparationError(
                f"panel does not contain configured columns {config.panel.sample_column!r} "
                f"and {config.panel.grouping_column!r}: {header}"
            ) from exc
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if max(sample_index, group_index) >= len(fields):
                raise PreparationError(f"panel line {line_number} has too few columns")
            sample = fields[sample_index].strip()
            group = fields[group_index].strip()
            if not sample or sample in seen:
                raise PreparationError(
                    f"panel line {line_number} has an empty or duplicate sample: {sample!r}"
                )
            seen.add(sample)
            if group in groups:
                groups[group].append(sample)
    empty = [group for group, samples in groups.items() if not samples]
    if empty:
        raise PreparationError(f"configured groups have no panel samples: {empty}")
    return groups


def vcf_samples(bcftools: str, path: Path) -> list[str]:
    output = run_command(
        [bcftools, "query", "--list-samples", str(path)], capture=True,
    )
    return [line for line in output.splitlines() if line]


def source_vcf_index_contigs(bcftools: str, path: Path) -> list[str]:
    """Read contigs from an upstream Tabix index lacking count metadata."""
    adjacent_tabix = Path(bcftools).with_name("tabix")
    tabix = (
        str(adjacent_tabix)
        if adjacent_tabix.is_file()
        else require_executable("tabix")
    )
    output = run_command([tabix, "--list-chroms", str(path)], capture=True)
    contigs = [line for line in output.splitlines() if line]
    if not contigs:
        raise PreparationError(f"source VCF index contains no contigs: {path}")
    return contigs


def vcf_index_stats(bcftools: str, path: Path) -> list[dict[str, int | str]]:
    output = run_command([bcftools, "index", "--stats", str(path)], capture=True)
    result: list[dict[str, int | str]] = []
    for line_number, line in enumerate(output.splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 3:
            raise PreparationError(
                f"unexpected bcftools index --stats line {line_number} for {path}: "
                f"{line!r}"
            )
        try:
            records = int(fields[2])
        except ValueError as exc:
            raise PreparationError(
                f"invalid indexed record count for {path}: {fields[2]!r}"
            ) from exc
        result.append({"contig": fields[0], "records": records})
    if not result:
        raise PreparationError(f"VCF index contains no contigs: {path}")
    return result


def validate_indexed_vcf(
    bcftools: str,
    path: Path,
    index_type: str,
    expected_contigs: list[str] | None = None,
) -> dict[str, Any]:
    suffix = ".csi" if index_type == "csi" else ".tbi"
    index = Path(str(path) + suffix)
    if not path.is_file() or path.stat().st_size == 0 or not index.is_file():
        raise PreparationError(f"VCF or {index_type.upper()} index is missing: {path}")
    index_stats = vcf_index_stats(bcftools, path)
    observed_contigs = [str(item["contig"]) for item in index_stats]
    if expected_contigs is not None and observed_contigs != expected_contigs:
        raise PreparationError(
            f"indexed contigs differ for {path}: observed {observed_contigs}, "
            f"expected {expected_contigs}"
        )
    records = sum(int(item["records"]) for item in index_stats)
    if records < 1:
        raise PreparationError(f"VCF contains no indexed records: {path}")
    return {
        "samples": len(vcf_samples(bcftools, path)),
        "records": records,
        "contigs": observed_contigs,
        "bytes": path.stat().st_size,
    }


def write_group_files(
    groups: dict[str, list[str]],
    metadata_directory: Path,
    filename_template: str,
) -> None:
    metadata_directory.mkdir(parents=True, exist_ok=True)
    for group, samples in groups.items():
        path = metadata_directory / filename_template.format(group=group)
        path.write_text("".join(f"{sample}\n" for sample in samples), encoding="utf-8")


def split_contig(
    config: WorkflowConfig,
    contig: ContigConfig,
    groups: dict[str, list[str]],
    paths: dict[str, Path],
    bcftools: str,
    logger: JsonLinesRunLogger,
) -> dict[str, dict[str, Any]]:
    source = paths["downloads"] / contig.filename
    source_samples = vcf_samples(bcftools, source)
    source_sample_set = set(source_samples)
    group_lookup = {
        sample: group
        for group, samples in groups.items()
        for sample in samples
        if sample in source_sample_set
    }
    missing = {
        group: len(set(samples) - source_sample_set)
        for group, samples in groups.items()
    }
    expected = {group: sum(sample in source_sample_set for sample in samples) for group, samples in groups.items()}
    if any(count == 0 for count in expected.values()):
        raise PreparationError(
            f"contig {contig.label} has no samples for configured groups: "
            f"{[group for group, count in expected.items() if count == 0]}"
        )
    destination_directory = paths["per_contig"] / contig.label
    existing: dict[str, dict[str, Any]] = {}
    for group in config.panel.groups:
        destination = destination_directory / config.naming.per_contig_vcf_template.format(
            group=group,
        )
        try:
            metrics = validate_indexed_vcf(
                bcftools,
                destination,
                config.processing.output_index_type,
                [contig.label],
            )
            if metrics["samples"] != expected[group]:
                raise PreparationError(
                    f"existing {destination} has {metrics['samples']} samples; "
                    f"expected {expected[group]}"
                )
            existing[group] = metrics
        except (PreparationError, subprocess.CalledProcessError):
            continue
    if len(existing) == len(config.panel.groups):
        logger.write("split_skipped", contig=contig.label, reason="validated_outputs_exist")
        return existing
    remaining_groups = [
        group for group in config.panel.groups if group not in existing
    ]

    destination_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{contig.label}.", dir=destination_directory,
    ) as temporary_name:
        temporary = Path(temporary_name)
        group_file = temporary / config.naming.split_groups_filename
        with group_file.open("w", encoding="utf-8") as handle:
            for sample in source_samples:
                group = group_lookup.get(sample)
                if group in remaining_groups:
                    handle.write(f"{sample}\t-\t{group}\n")
        logger.write(
            "split_started",
            contig=contig.label,
            source=str(source),
            source_samples=len(source_samples),
            selected_samples=expected,
            panel_samples_absent_from_source=missing,
            reused_groups=sorted(existing),
            remaining_groups=remaining_groups,
        )
        run_command([
            bcftools,
            "+split",
            str(source),
            "--groups-file",
            str(group_file),
            "--output",
            str(temporary),
            "--output-type",
            "z",
        ])
        metrics_by_group = dict(existing)
        for group in remaining_groups:
            split_path = temporary / config.naming.per_contig_vcf_template.format(
                group=group,
            )
            if not split_path.is_file():
                raise PreparationError(f"bcftools +split did not create {split_path}")
            prepared = temporary / config.naming.split_prepared_vcf_template.format(
                group=group,
            )
            if config.processing.recalculate_info_tags:
                run_command([
                    bcftools,
                    "+fill-tags",
                    str(split_path),
                    "--output-type",
                    "z",
                    "--output",
                    str(prepared),
                    "--threads",
                    str(config.processing.bcftools_threads),
                    "--",
                    "--tags",
                    ",".join(config.processing.recalculate_info_tags),
                ])
            else:
                prepared = split_path
            run_command([
                bcftools,
                "index",
                "--force",
                f"--{config.processing.output_index_type}",
                "--threads",
                str(config.processing.bcftools_threads),
                str(prepared),
            ])
            metrics = validate_indexed_vcf(
                bcftools,
                prepared,
                config.processing.output_index_type,
                [contig.label],
            )
            if metrics["samples"] != expected[group]:
                raise PreparationError(
                    f"split {contig.label}/{group} has {metrics['samples']} samples; "
                    f"expected {expected[group]}"
                )
            destination = destination_directory / config.naming.per_contig_vcf_template.format(
                group=group,
            )
            index_suffix = ".csi" if config.processing.output_index_type == "csi" else ".tbi"
            os.replace(prepared, destination)
            os.replace(Path(str(prepared) + index_suffix), Path(str(destination) + index_suffix))
            metrics_by_group[group] = metrics
        logger.write("split_completed", contig=contig.label, groups=metrics_by_group)
        return metrics_by_group


def concat_outputs(
    config: WorkflowConfig,
    output_set: OutputSetConfig,
    group: str,
    paths: dict[str, Path],
    bcftools: str,
    logger: JsonLinesRunLogger,
) -> Path:
    per_contig_name = config.naming.per_contig_vcf_template.format(group=group)
    inputs = [paths["per_contig"] / label / per_contig_name for label in output_set.contigs]
    sample_lists = [vcf_samples(bcftools, path) for path in inputs]
    if any(samples != sample_lists[0] for samples in sample_lists[1:]):
        raise PreparationError(
            f"cannot concatenate output set {output_set.name!r} for {group}: "
            "the source contigs have different ordered sample lists"
        )
    filename = output_set.filename_template.format(group=group)
    ensure_relative_filename(filename, f"output set {output_set.name!r} filename")
    destination = paths["final"] / filename
    try:
        metrics = validate_indexed_vcf(
            bcftools,
            destination,
            config.processing.output_index_type,
            output_set.contigs,
        )
        if metrics["samples"] != len(sample_lists[0]):
            raise PreparationError("sample count differs")
        logger.write(
            "concat_skipped",
            output_set=output_set.name,
            group=group,
            path=str(destination),
            reason="validated_output_exists",
        )
        return destination
    except (PreparationError, subprocess.CalledProcessError):
        pass

    paths["final"].mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{group}.{output_set.name}.", dir=paths["final"]) as name:
        temporary = Path(name)
        file_list = temporary / config.naming.concat_input_list_filename
        file_list.write_text("".join(f"{path}\n" for path in inputs), encoding="utf-8")
        prepared = temporary / filename
        logger.write(
            "concat_started",
            output_set=output_set.name,
            group=group,
            inputs=[str(path) for path in inputs],
        )
        run_command([
            bcftools,
            "concat",
            "--file-list",
            str(file_list),
            "--output-type",
            "z",
            "--output",
            str(prepared),
            "--threads",
            str(config.processing.bcftools_threads),
        ])
        run_command([
            bcftools,
            "index",
            "--force",
            f"--{config.processing.output_index_type}",
            "--threads",
            str(config.processing.bcftools_threads),
            str(prepared),
        ])
        metrics = validate_indexed_vcf(
            bcftools,
            prepared,
            config.processing.output_index_type,
            output_set.contigs,
        )
        if metrics["samples"] != len(sample_lists[0]):
            raise PreparationError(
                f"final output has {metrics['samples']} samples; expected {len(sample_lists[0])}"
            )
        index_suffix = ".csi" if config.processing.output_index_type == "csi" else ".tbi"
        os.replace(prepared, destination)
        os.replace(Path(str(prepared) + index_suffix), Path(str(destination) + index_suffix))
    logger.write(
        "concat_completed",
        output_set=output_set.name,
        group=group,
        path=str(destination),
        metrics=metrics,
    )
    return destination


def validate_sources(
    config: WorkflowConfig,
    paths: dict[str, Path],
    bcftools: str,
    logger: JsonLinesRunLogger,
) -> None:
    for contig in config.contigs:
        source = paths["downloads"] / contig.filename
        index = Path(str(source) + config.processing.index_suffix)
        if not source.is_file() or not index.is_file():
            raise PreparationError(f"source VCF or index is missing for contig {contig.label}")
        samples = vcf_samples(bcftools, source)
        observed_contigs = source_vcf_index_contigs(bcftools, source)
        if observed_contigs != [contig.label]:
            raise PreparationError(
                f"source VCF contigs differ for {source}: observed "
                f"{observed_contigs}, expected {[contig.label]}"
            )
        if not samples:
            raise PreparationError(f"source VCF has no samples: {source}")
        logger.write(
            "source_validated",
            contig=contig.label,
            path=str(source),
            samples=len(samples),
            records=None,
            record_count_status="unavailable_in_upstream_tbi",
            observed_contigs=observed_contigs,
            bytes=source.stat().st_size,
        )


def write_manifest(
    config: WorkflowConfig,
    config_path: Path,
    outputs: list[Path],
    paths: dict[str, Path],
    bcftools: str,
) -> None:
    manifest = {
        "resource": config.resource.model_dump(mode="json"),
        "panel": config.panel.model_dump(mode="json"),
        "resolved_configuration": config.model_dump(mode="json"),
        "configuration_path": str(config_path),
        "configuration_sha256": sha256(config_path),
        "bcftools_version": run_command(
            [bcftools, "--version"], capture=True,
        ).splitlines()[0],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        ],
    }
    paths["metadata"].mkdir(parents=True, exist_ok=True)
    manifest_path = paths["metadata"] / config.naming.manifest_filename
    with manifest_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False)


def prepare(config_path: Path, stage: str) -> None:
    config = load_config(config_path)
    paths = layout(config)
    bcftools = require_executable(config.processing.bcftools_executable)
    curl = require_executable(config.processing.curl_executable)
    paths["logs"].mkdir(parents=True, exist_ok=True)
    log_path = paths["logs"] / config.naming.log_filename
    logger = JsonLinesRunLogger(log_path)
    logger.record_unfinished_previous_run()
    logger.write(
        "run_started",
        stage=stage,
        config=str(config_path),
        resolved_configuration=config.model_dump(mode="json"),
    )
    try:
        if stage in {"download", "all"}:
            download_sources(config, paths, curl, logger)
        if stage in {"prepare", "validate", "all"}:
            panel_path = paths["downloads"] / config.panel.filename
            if not panel_path.is_file():
                raise PreparationError(
                    f"panel is missing; run the download stage first: {panel_path}"
                )
            groups = read_panel(config, panel_path)
            write_group_files(
                groups,
                paths["metadata"],
                config.naming.group_sample_list_template,
            )
            validate_sources(config, paths, bcftools, logger)
        else:
            groups = {}
        outputs: list[Path] = []
        if stage in {"prepare", "all"}:
            logger.write(
                "parallel_processing_started",
                contig_workers=CONTIG_PROCESSING_WORKERS,
                concatenation_workers=CONCATENATION_WORKERS,
                bcftools_threads_per_process=config.processing.bcftools_threads,
            )
            with ThreadPoolExecutor(
                max_workers=CONTIG_PROCESSING_WORKERS,
            ) as executor:
                futures = {
                    executor.submit(
                        split_contig,
                        config,
                        contig,
                        groups,
                        paths,
                        bcftools,
                        logger,
                    ): contig.label
                    for contig in config.contigs
                }
                for future in as_completed(futures):
                    future.result()
            with ThreadPoolExecutor(
                max_workers=CONCATENATION_WORKERS,
            ) as executor:
                futures = [
                    executor.submit(
                        concat_outputs,
                        config,
                        output_set,
                        group,
                        paths,
                        bcftools,
                        logger,
                    )
                    for output_set in config.output_sets
                    for group in config.panel.groups
                ]
                outputs = [future.result() for future in as_completed(futures)]
            write_manifest(config, config_path, outputs, paths, bcftools)
            if config.processing.remove_per_contig_outputs_after_success:
                shutil.rmtree(paths["per_contig"])
            if not config.processing.keep_downloads_after_success:
                shutil.rmtree(paths["downloads"])
        elif stage == "validate":
            for output_set in config.output_sets:
                for group in config.panel.groups:
                    filename = output_set.filename_template.format(group=group)
                    output = paths["final"] / filename
                    metrics = validate_indexed_vcf(
                        bcftools,
                        output,
                        config.processing.output_index_type,
                        output_set.contigs,
                    )
                    logger.write("final_output_validated", path=str(output), metrics=metrics)
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
        "--stage",
        required=True,
        choices=("download", "prepare", "validate", "all"),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        prepare(args.config.expanduser().resolve(), args.stage)
    except (PreparationError, subprocess.CalledProcessError, OSError, yaml.YAMLError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
