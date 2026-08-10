#!/usr/bin/env python3
"""Install pinned MAGMA functional-mapping resources and emit valid YAML."""

from __future__ import annotations

import csv
import fnmatch
import os
from pathlib import Path
import re
import shutil
import tarfile
import zipfile
from collections import Counter, defaultdict

import yaml

from postgwas.config import load_configuration
from postgwas.core.resource_preparation import (
    JsonLinesRunLogger,
    ResourcePreparationError,
    create_staging_directory,
    require_executable,
    resumable_download,
    run_command,
    sha256,
)


def _safe_members(names: list[str], destination: Path) -> None:
    root = destination.resolve()
    for name in names:
        target = (destination / name).resolve()
        if target != root and root not in target.parents:
            raise ResourcePreparationError(
                "archive contains a path outside its destination: %s" % name
            )


def _extract(
    archive: Path,
    archive_type: str,
    destination: Path,
    excluded_members: list[str],
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    if archive_type == "zip":
        with zipfile.ZipFile(archive) as handle:
            members = [
                member for member in handle.infolist()
                if not _excluded(Path(member.filename), excluded_members)
            ]
            _safe_members([member.filename for member in members], destination)
            handle.extractall(destination, members=members)
        return
    if archive_type == "tar.gz":
        with tarfile.open(archive, "r:gz") as handle:
            members = [
                member for member in handle.getmembers()
                if not _excluded(Path(member.name), excluded_members)
            ]
            _safe_members([member.name for member in members], destination)
            handle.extractall(destination, members=members, filter="data")
        return
    raise ResourcePreparationError("unsupported archive type: %s" % archive_type)


def _relative_source(extracted: Path, archive_root: str) -> Path:
    source = extracted / archive_root if archive_root else extracted
    if not source.is_dir():
        raise ResourcePreparationError(
            "configured archive root does not exist: %s" % source
        )
    return source


def _excluded(relative: Path, patterns: list[str]) -> bool:
    name = relative.as_posix()
    return any(fnmatch.fnmatch(name, pattern) for pattern in patterns)


def _copy_selected(source: Path, destination: Path, specification: dict) -> int:
    selected: dict[str, Path] = {}
    for pattern in specification["include"]:
        for path in source.glob(pattern):
            if path.is_file():
                relative = path.relative_to(source)
                if not _excluded(relative, specification.get("exclude", [])):
                    selected[relative.as_posix()] = path
    if not selected:
        raise ResourcePreparationError(
            "download %s selected no reference files" % specification["archive_name"]
        )
    renames = specification.get("renames", {})
    for relative_name, path in sorted(selected.items()):
        target_name = renames.get(relative_name, relative_name)
        target = destination / target_name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    return len(selected)


def _slug(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    if not normalized or not re.match(r"[a-z]", normalized):
        raise ResourcePreparationError("cannot create a safe mapping name from %r" % value)
    return normalized


def _definition_common(group: dict, context: str) -> dict:
    definition = {
        "method": group["method"],
        "display_name": "%s · %s" % (
            group["display_prefix"], context.replace("_", " "),
        ),
        "genome_build": group["genome_build"],
        "population": group["population"],
        "gene_id_type": group["gene_id_type"],
        "context": context,
        "source_name": group["source_name"],
        "source_version": group["source_version"],
        "source_url": group["source_url"],
        "result_statistic_type": group["result_statistic_type"],
        "result_statistic_interpretation": group[
            "result_statistic_interpretation"
        ],
        "gene_settings": group["gene_settings"],
    }
    if "minimum_gene_id_overlap_fraction" in group:
        definition["minimum_gene_id_overlap_fraction"] = group[
            "minimum_gene_id_overlap_fraction"
        ]
    return definition


class _ReadableYamlDumper(yaml.SafeDumper):
    def ignore_aliases(self, data):
        return True


def _write_yaml(path: Path, value: dict) -> None:
    path.write_text(
        yaml.dump(
            value,
            Dumper=_ReadableYamlDumper,
            sort_keys=False,
            width=120,
            allow_unicode=True,
        ),
        encoding="utf-8",
    )


def _split_record(text: str, delimiter: re.Pattern[str], line_number: int, label: str) -> list[str]:
    fields = delimiter.split(text.strip())
    if not fields or any(not value.strip() for value in fields):
        raise ResourcePreparationError(
            "%s contains an empty field at line %d" % (label, line_number)
        )
    return fields


def _read_gene_identifier_bridge(
    exported_file: Path,
    gene_location_file: Path,
    specification: dict,
) -> list[tuple[str, str, str]]:
    source_to_bridge: dict[str, set[str]] = defaultdict(set)
    with exported_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            specification["rdata_source_gene_column"],
            specification["rdata_bridge_column"],
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ResourcePreparationError(
                "RData identifier export is missing columns: %s"
                % ", ".join(sorted(required))
            )
        for record in reader:
            source = re.sub(
                specification["source_gene_version_pattern"],
                "",
                record[specification["rdata_source_gene_column"]].strip(),
            )
            bridge = record[specification["rdata_bridge_column"]].strip()
            if source and bridge:
                source_to_bridge[source].add(bridge)

    bridge_to_target: dict[str, set[str]] = defaultdict(set)
    delimiter = re.compile(specification["gene_location_delimiter_pattern"])
    required_index = max(
        specification["gene_location_target_column"],
        specification["gene_location_bridge_column"],
    )
    with gene_location_file.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(specification["comment_prefix"]):
                continue
            fields = _split_record(
                text, delimiter, line_number, "gene-location reference",
            )
            if len(fields) <= required_index:
                raise ResourcePreparationError(
                    "gene-location reference line %d does not contain configured "
                    "column %d" % (line_number, required_index)
                )
            bridge_to_target[
                fields[specification["gene_location_bridge_column"]]
            ].add(fields[specification["gene_location_target_column"]])

    crosswalk = []
    for source, bridges in source_to_bridge.items():
        if len(bridges) != 1:
            continue
        bridge = next(iter(bridges))
        targets = bridge_to_target.get(bridge, set())
        if len(targets) == 1:
            crosswalk.append((source, bridge, next(iter(targets))))
    if not crosswalk:
        raise ResourcePreparationError(
            "identifier references produced no unambiguous Ensembl-to-Entrez mapping"
        )
    return sorted(crosswalk)


def _annotation_gene_ids(path: Path, specification: dict) -> set[str]:
    identifiers = set()
    required_index = specification["annotation_gene_column"]
    delimiter = re.compile(specification["annotation_delimiter_pattern"])
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(specification["comment_prefix"]):
                continue
            fields = delimiter.split(text)
            if len(fields) <= required_index:
                raise ResourcePreparationError(
                    "eMAGMA annotation %s line %d lacks configured gene column %d"
                    % (path, line_number, required_index)
                )
            identifiers.add(fields[required_index])
    if not identifiers:
        raise ResourcePreparationError(
            "eMAGMA annotation contains no gene identifiers: %s" % path
        )
    return identifiers


def _harmonise_emagma_network_identifiers(
    staging: Path,
    work: Path,
    configuration_directory: Path,
    specification: dict,
    logger,
) -> list[dict]:
    rscript = require_executable(specification["rscript_executable"])
    export_script = configuration_directory / specification["r_export_script"]
    if not export_script.is_file() or export_script.stat().st_size == 0:
        raise ResourcePreparationError(
            "configured RData export script does not exist or is empty: %s"
            % export_script
        )
    rdata = staging / specification["rdata_file"]
    gene_location = staging / specification["gene_location_file"]
    for path, label in (
        (rdata, "identifier RData"),
        (gene_location, "gene-location identifier reference"),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ResourcePreparationError("%s is missing or empty: %s" % (label, path))
    exported = work / specification["r_export_file"]
    exported.parent.mkdir(parents=True, exist_ok=True)
    run_command([
        rscript,
        str(export_script),
        str(rdata),
        str(exported),
        specification["rdata_object"],
        specification["rdata_source_gene_column"],
        specification["rdata_bridge_column"],
    ])
    crosswalk = _read_gene_identifier_bridge(
        exported, gene_location, specification,
    )
    source_to_target = {source: target for source, _, target in crosswalk}
    mapped_sources = set(source_to_target)
    crosswalk_destination = staging / specification["crosswalk_file"]
    crosswalk_destination.parent.mkdir(parents=True, exist_ok=True)
    crosswalk_columns = specification["crosswalk_columns"]
    with crosswalk_destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=specification["output_delimiter"])
        writer.writerow([
            crosswalk_columns["source_gene"],
            crosswalk_columns["bridge_gene"],
            crosswalk_columns["target_gene"],
        ])
        writer.writerows(crosswalk)

    source_directory = staging / specification["source_directory"]
    networks = sorted(source_directory.glob(specification["source_pattern"]))
    if len(networks) != specification["expected_network_files"]:
        raise ResourcePreparationError(
            "expected %d eMAGMA network files but found %d in %s"
            % (
                specification["expected_network_files"],
                len(networks),
                source_directory,
            )
        )
    destination_directory = staging / specification["destination_directory"]
    destination_directory.mkdir(parents=True, exist_ok=True)
    membership_delimiter = re.compile(
        specification["membership_delimiter_pattern"]
    )
    membership_index = max(
        specification["membership_set_column"],
        specification["membership_gene_column"],
    )
    summaries = []
    for network in networks:
        source_memberships = 0
        source_genes = set()
        retained_pairs: list[tuple[str, str]] = []
        observed_pairs = set()
        with network.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text or text.startswith(specification["comment_prefix"]):
                    continue
                fields = _split_record(
                    text, membership_delimiter, line_number, "eMAGMA network",
                )
                if len(fields) <= membership_index:
                    raise ResourcePreparationError(
                        "eMAGMA network %s line %d lacks configured membership "
                        "column %d" % (network, line_number, membership_index)
                    )
                source_memberships += 1
                source_gene = re.sub(
                    specification["source_gene_version_pattern"],
                    "",
                    fields[specification["membership_gene_column"]],
                )
                source_genes.add(source_gene)
                target_gene = source_to_target.get(source_gene)
                pair = (
                    fields[specification["membership_set_column"]], target_gene,
                )
                if target_gene is not None and pair not in observed_pairs:
                    retained_pairs.append(pair)
                    observed_pairs.add(pair)
        if not source_genes:
            raise ResourcePreparationError(
                "eMAGMA network contains no gene memberships: %s" % network
            )
        mapped_source_genes = source_genes & mapped_sources
        mapping_fraction = len(mapped_source_genes) / len(source_genes)
        if mapping_fraction < specification["minimum_unique_gene_mapping_fraction"]:
            raise ResourcePreparationError(
                "eMAGMA network %s maps only %d/%d unique Ensembl genes (%.2f%%); "
                "the configured minimum is %.2f%%"
                % (
                    network,
                    len(mapped_source_genes),
                    len(source_genes),
                    mapping_fraction * 100,
                    specification["minimum_unique_gene_mapping_fraction"] * 100,
                )
            )
        destination = destination_directory / network.name
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=specification["output_delimiter"])
            writer.writerows(retained_pairs)

        annotation = (
            staging
            / specification["annotation_directory"]
            / (network.stem + specification["annotation_suffix"])
        )
        if not annotation.is_file():
            raise ResourcePreparationError(
                "eMAGMA network has no tissue-matched annotation: %s" % network
            )
        annotation_genes = _annotation_gene_ids(annotation, specification)
        target_genes = {target for _, target in retained_pairs}
        overlap = target_genes & annotation_genes
        comparison_size = min(len(target_genes), len(annotation_genes))
        overlap_fraction = len(overlap) / comparison_size if comparison_size else 0.0
        if overlap_fraction < specification["minimum_annotation_overlap_fraction"]:
            raise ResourcePreparationError(
                "eMAGMA network %s has only %d/%d unique genes in its annotation "
                "(%.2f%%); the configured minimum is %.2f%%"
                % (
                    network,
                    len(overlap),
                    comparison_size,
                    overlap_fraction * 100,
                    specification["minimum_annotation_overlap_fraction"] * 100,
                )
            )
        summaries.append(
            {
                "network": network.name,
                "source_memberships": source_memberships,
                "retained_memberships": len(retained_pairs),
                "source_unique_genes": len(source_genes),
                "mapped_unique_genes": len(mapped_source_genes),
                "unique_gene_mapping_fraction": mapping_fraction,
                "annotation_unique_genes": len(annotation_genes),
                "overlapping_annotation_genes": len(overlap),
                "annotation_overlap_fraction": overlap_fraction,
            }
        )

    summary_destination = staging / specification["summary_file"]
    summary_destination.parent.mkdir(parents=True, exist_ok=True)
    summary_columns = specification["summary_columns"]
    with summary_destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[summary_columns[role] for role in summary_columns],
            delimiter=specification["output_delimiter"],
        )
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {summary_columns[role]: summary[role] for role in summary_columns}
            )
    logger.write(
        "emagma_network_identifiers_harmonised",
        source_gene_id_type=specification["source_gene_id_type"],
        target_gene_id_type=specification["target_gene_id_type"],
        crosswalk_records=len(crosswalk),
        networks=len(summaries),
        summary=specification["summary_file"],
    )
    return summaries


def _generated_module_config(staging: Path, final_root: Path, configuration: dict) -> dict:
    build = configuration["genome_build"]
    population = configuration["population"]
    positional = configuration["positional"]
    definitions = {
        "positional": {
            "method": "positional",
            "display_name": positional["display_name"],
            "genome_build": build,
            "population": population,
            "gene_id_type": positional["gene_id_type"],
            "context": positional["context"],
            "source_name": positional["source_name"],
            "source_version": positional["source_version"],
            "source_url": positional["source_url"],
            "result_statistic_type": positional["result_statistic_type"],
            "result_statistic_interpretation": positional[
                "result_statistic_interpretation"
            ],
            "gene_location_file": str(final_root / configuration["gene_location_file"]),
        }
    }
    for configured_group in configuration["annotation_groups"]:
        group = {
            **configured_group,
            "genome_build": build,
            "population": population,
        }
        excluded = group.get("exclude_globs", [])
        for annotation in sorted(staging.glob(group["glob"])):
            if not annotation.is_file() or any(
                fnmatch.fnmatch(annotation.name, pattern) for pattern in excluded
            ):
                continue
            source_context = (
                re.sub(r"(?:\.genes\.annot(?:\.gz)?|\.annot)$", "", annotation.name)
                if group["context_from_filename"]
                else group["display_prefix"]
            )
            context = group.get("context_aliases", {}).get(
                source_context, source_context,
            )
            name = "%s_%s" % (group["name_prefix"], _slug(source_context))
            definition = _definition_common(group, context)
            definition["gene_annotation_file"] = str(
                final_root / annotation.relative_to(staging)
            )
            network_directory = group.get("gene_set_directory")
            if network_directory:
                network = staging / network_directory / (
                    source_context + group["gene_set_suffix"]
                )
                if network.is_file():
                    definition["gene_set_file"] = str(
                        final_root / network.relative_to(staging)
                    )
                    definition["gene_set_format"] = group["gene_set_format"]
            definitions[name] = definition

    n_config = configuration["n_magma"]
    for tissue, values in n_config["tissues"].items():
        components = [staging / value for value in values["components"]]
        missing = [str(path) for path in components if not path.is_file()]
        if missing:
            raise ResourcePreparationError(
                "nMAGMA components are missing for %s: %s"
                % (tissue, ", ".join(missing))
            )
        definitions["n_magma_%s" % _slug(tissue)] = {
            "method": "n_magma",
            "display_name": "%s · %s" % (
                n_config["display_prefix"], tissue.replace("_", " "),
            ),
            "genome_build": build,
            "population": population,
            "gene_id_type": n_config["gene_id_type"],
            "context": tissue,
            "source_name": n_config["source_name"],
            "source_version": n_config["source_version"],
            "source_url": n_config["source_url"],
            "result_statistic_type": n_config["result_statistic_type"],
            "result_statistic_interpretation": n_config[
                "result_statistic_interpretation"
            ],
            "gene_location_file": str(final_root / n_config["gene_location_file"]),
            "annotation_window_upstream_kb": n_config[
                "annotation_window_upstream_kb"
            ],
            "annotation_window_downstream_kb": n_config[
                "annotation_window_downstream_kb"
            ],
            "gene_annotation_files": [
                str(final_root / path.relative_to(staging)) for path in components
            ],
        }

    chrom = configuration["chrom_magma"]
    definitions[chrom["name"]] = {
        "method": "chrom_magma",
        "display_name": chrom["display_name"],
        "genome_build": build,
        "population": population,
        "gene_id_type": chrom["gene_id_type"],
        "context": chrom["context"],
        "source_name": chrom["source_name"],
        "source_version": chrom["source_version"],
        "source_url": chrom["source_url"],
        "result_statistic_type": chrom["result_statistic_type"],
        "result_statistic_interpretation": chrom[
            "result_statistic_interpretation"
        ],
        "regulatory_element_location_file": str(
            final_root / chrom["regulatory_element_location_file"]
        ),
        "element_to_gene_file": str(final_root / chrom["element_to_gene_file"]),
        "annotation_window_upstream_kb": chrom["annotation_window_upstream_kb"],
        "annotation_window_downstream_kb": chrom["annotation_window_downstream_kb"],
    }
    observed_counts = Counter(
        definition["method"] for definition in definitions.values()
    )
    expected_counts = configuration["expected_mapping_definition_counts"]
    if dict(observed_counts) != expected_counts:
        raise ResourcePreparationError(
            "generated mapping-definition counts do not match configuration: "
            "observed %s; expected %s" % (dict(observed_counts), expected_counts)
        )
    return {
        "genome_build": build,
        "population": population,
        "input": {
            "ld_reference_prefix": str(
                final_root / configuration["ld_reference_prefix"]
            ),
            "gene_location_file": str(
                final_root / configuration["gene_location_file"]
            ),
        },
        "mapping": {
            "selected": configuration["selected"],
            "primary": configuration["primary"],
            "definitions": definitions,
        },
    }


def _generated_run_config(
    staging: Path,
    final_root: Path,
    configuration: dict,
) -> dict:
    """Nest generated MAGMA values in the reloadable pipeline namespace."""
    return {
        "modules": {
            "magma": _generated_module_config(
                staging, final_root, configuration,
            ),
        },
    }


def _preparation_digests(config_path: Path, config: dict) -> dict[str, str]:
    configuration_directory = config_path.parent.resolve()
    required = {
        "config.yaml": config_path.resolve(),
        "preparation.py": Path(__file__).resolve(),
        config["identifier_harmonisation"]["emagma_networks"]["r_export_script"]: (
            configuration_directory
            / config["identifier_harmonisation"]["emagma_networks"]["r_export_script"]
        ),
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise ResourcePreparationError(
            "resource-preparation source files are missing: %s" % ", ".join(missing)
        )
    return {name: sha256(path) for name, path in required.items()}


def _manifest(staging: Path, config: dict, preparation_digests: dict[str, str]) -> dict:
    excluded = {config["manifest_file"], config["run_log_file"]}
    files = []
    for path in sorted(item for item in staging.rglob("*") if item.is_file()):
        relative = path.relative_to(staging).as_posix()
        if relative not in excluded:
            files.append(
                {"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)}
            )
    return {
        "schema_version": config["schema_version"],
        "genome_build": config["configuration"]["genome_build"],
        "population": config["configuration"]["population"],
        "preparation_file_sha256": preparation_digests,
        "files": files,
        "sources": [
            {
                "name": value["source_name"],
                "version": value["source_version"],
                "url": value["url"],
                "archive_sha256": value["sha256"],
                "license": value["license"],
            }
            for value in config["downloads"].values()
        ],
    }


def _validate_existing(
    output: Path,
    config: dict,
    preparation_digests: dict[str, str],
) -> bool:
    manifest_path = output / config["manifest_file"]
    generated_config = output / config["generated_config_file"]
    if not output.exists():
        return False
    if not manifest_path.is_file() or not generated_config.is_file():
        raise ResourcePreparationError(
            "existing output is incomplete and will not be replaced: %s" % output
        )
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("preparation_file_sha256") != preparation_digests:
        raise ResourcePreparationError(
            "existing output was generated with a different resource-preparation "
            "configuration or script and will not be reused: %s" % output
        )
    for record in manifest.get("files", []):
        path = output / record["path"]
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ResourcePreparationError(
                "existing installed resource failed checksum validation: %s" % path
            )
    load_configuration(generated_config)
    return True


def refresh_metadata(
    config_path: Path,
    output: Path,
    *,
    configuration: dict | None = None,
) -> Path:
    """Regenerate validated configuration metadata without changing resources."""
    config_path = config_path.expanduser().resolve()
    output = output.expanduser().resolve()
    config = configuration or _load_preparation_configuration(config_path)
    manifest_path = output / config["manifest_file"]
    generated_path = output / config["generated_config_file"]
    if not manifest_path.is_file():
        raise ResourcePreparationError(
            "cannot refresh an installation without its manifest: %s" % output
        )
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    records = manifest.get("files", [])
    generated_relative = config["generated_config_file"]
    generated_records = [
        record for record in records if record.get("path") == generated_relative
    ]
    if len(generated_records) != 1:
        raise ResourcePreparationError(
            "resource manifest must contain exactly one generated configuration"
        )
    for record in records:
        if record.get("path") == generated_relative:
            continue
        path = output / record["path"]
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ResourcePreparationError(
                "installed scientific resource failed checksum validation: %s" % path
            )

    temporary_config = generated_path.with_name(generated_path.name + ".refreshing")
    temporary_manifest = manifest_path.with_name(manifest_path.name + ".refreshing")
    logger = JsonLinesRunLogger(output / config["run_log_file"])
    logger.write("metadata_refresh_started", output=str(output))
    try:
        generated = _generated_run_config(
            output, output, config["configuration"],
        )
        _write_yaml(temporary_config, generated)
        load_configuration(temporary_config)
        generated_records[0].update({
            "bytes": temporary_config.stat().st_size,
            "sha256": sha256(temporary_config),
        })
        manifest["preparation_file_sha256"] = _preparation_digests(
            config_path, config,
        )
        _write_yaml(temporary_manifest, manifest)
        os.replace(temporary_config, generated_path)
        os.replace(temporary_manifest, manifest_path)
        logger.write(
            "metadata_refresh_completed",
            generated_config=str(generated_path),
        )
        return output
    except BaseException as exc:
        logger.write(
            "metadata_refresh_failed",
            error_type=type(exc).__name__,
            error=str(exc),
        )
        raise
    finally:
        temporary_config.unlink(missing_ok=True)
        temporary_manifest.unlink(missing_ok=True)


def prepare(
    config_path: Path,
    output: Path,
    cache: Path,
    *,
    configuration: dict | None = None,
) -> Path:
    config_path = config_path.expanduser().resolve()
    config = configuration or _load_preparation_configuration(config_path)
    output = output.expanduser().resolve()
    cache = cache.expanduser().resolve()
    preparation_digests = _preparation_digests(config_path, config)
    if _validate_existing(output, config, preparation_digests):
        return output
    staging = create_staging_directory(output)
    logger = JsonLinesRunLogger(staging / config["run_log_file"])
    logger.write("run_started", output=str(output), config=str(config_path.resolve()))
    work = staging / ".work"
    try:
        curl = require_executable(config["download_executable"])
        for name, specification in config["downloads"].items():
            archive = resumable_download(
                curl=curl,
                url=specification["url"],
                destination=cache / specification["archive_name"],
                retries=config["download_retries"],
                partial_suffix=config["partial_suffix"],
                logger=logger,
                expected_sha256=specification["sha256"],
            )
            extracted = work / name
            _extract(
                archive,
                specification["archive_type"],
                extracted,
                config["archive_member_exclude"],
            )
            source = _relative_source(extracted, specification["archive_root"])
            installed = _copy_selected(
                source, staging / specification["destination"], specification,
            )
            for nested in specification.get("nested", []):
                matches = sorted(source.glob(nested["pattern"]))
                if not matches:
                    raise ResourcePreparationError(
                        "nested archive pattern selected nothing: %s"
                        % nested["pattern"]
                    )
                for nested_archive in matches:
                    _extract(
                        nested_archive,
                        nested["archive_type"],
                        staging / nested["destination"],
                        config["archive_member_exclude"],
                    )
            logger.write("resource_installed", name=name, files=installed)

        _harmonise_emagma_network_identifiers(
            staging,
            work,
            config_path.parent,
            config["identifier_harmonisation"]["emagma_networks"],
            logger,
        )
        generated = _generated_run_config(
            staging, output, config["configuration"],
        )
        generated_path = staging / config["generated_config_file"]
        generated_path.parent.mkdir(parents=True, exist_ok=True)
        _write_yaml(generated_path, generated)
        load_configuration(generated_path)
        shutil.rmtree(work)
        manifest = _manifest(staging, config, preparation_digests)
        manifest_path = staging / config["manifest_file"]
        _write_yaml(manifest_path, manifest)
        logger.write(
            "run_completed",
            resources=len(manifest["files"]),
            generated_config=str(output / config["generated_config_file"]),
        )
        os.replace(staging, output)
        return output
    except BaseException as exc:
        logger.write("run_failed", error_type=type(exc).__name__, error=str(exc))
        raise


def _load_preparation_configuration(config_path: Path) -> dict:
    path = config_path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ResourcePreparationError(
            "resource configuration does not exist or is empty: %s" % path
        )
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ResourcePreparationError("unsupported resource configuration schema")
    return config


def packaged_configuration() -> tuple[Path, dict]:
    """Return the installed canonical YAML path and loaded mapping."""
    path = Path(__file__).with_name("config.yaml").resolve()
    return path, _load_preparation_configuration(path)


def run(
    action: str,
    *,
    config_path: Path | None = None,
    configuration: dict | None = None,
    output_directory: Path | None = None,
    cache_directory: Path | None = None,
) -> tuple[Path, Path]:
    """Resolve portable defaults once and execute one preparation operation."""
    if action not in {"prepare", "refresh"}:
        raise ResourcePreparationError(
            "resource action must be either 'prepare' or 'refresh'"
        )
    if config_path is None:
        selected_config, packaged = packaged_configuration()
        config = configuration or packaged
    else:
        selected_config = config_path.expanduser().resolve()
        config = configuration or _load_preparation_configuration(selected_config)
    output = (
        output_directory or Path(config["default_output_directory"])
    ).expanduser().resolve()
    cache = (
        cache_directory or Path(config["default_cache_directory"])
    ).expanduser().resolve()
    installed = (
        refresh_metadata(
            selected_config,
            output,
            configuration=config,
        )
        if action == "refresh"
        else prepare(
            selected_config,
            output,
            cache,
            configuration=config,
        )
    )
    return installed, installed / config["generated_config_file"]


__all__ = ["packaged_configuration", "prepare", "refresh_metadata", "run"]
