#!/usr/bin/env python3
"""Map GMT gene pathways to exact PLINK BIM IDs for GCTA fastBAT sets."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import io
from pathlib import Path
import shlex
import shutil
import sqlite3
import sys
import tempfile

import yaml

from postgwas.core.resource_preparation import (
    ResourcePreparationError,
    create_staging_directory,
    required_file,
    sha256,
    unique_nonempty_values,
    validate_output_filename,
    write_text,
)


OFFICIAL_FASTBAT_DOCUMENTATION = (
    "https://yanglab.westlake.edu.cn/software/gcta/#fastBAT"
)


@dataclass(frozen=True)
class Pathway:
    index: int
    name: str
    description: str
    genes: tuple[str, ...]
    duplicate_gene_entries: int


@dataclass(frozen=True)
class GeneInterval:
    chromosome: str
    start: int
    end: int
    gene: str


def _normalize_chromosome(value: str, policy: str) -> str:
    if policy == "strip_chr_prefix" and value.lower().startswith("chr"):
        return value[3:]
    return value


def _parse_gmt(path: Path, duplicate_gene_policy: str) -> list[Pathway]:
    pathways: list[Pathway] = []
    names: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ResourcePreparationError(
                    f"GMT line {line_number} must contain pathway, description, "
                    "and at least one gene"
                )
            name = fields[0].strip()
            description = fields[1].strip()
            genes = [value.strip() for value in fields[2:]]
            if not name or any(character.isspace() for character in name):
                raise ResourcePreparationError(
                    f"GMT line {line_number} has an empty or whitespace-containing "
                    "pathway ID"
                )
            if name == "END":
                raise ResourcePreparationError("GMT pathway ID END is reserved by GCTA")
            if name in names:
                raise ResourcePreparationError(f"GMT repeats pathway ID {name!r}")
            if any(not gene for gene in genes):
                raise ResourcePreparationError(
                    f"GMT pathway {name!r} contains an empty gene identifier"
                )
            unique_genes = tuple(dict.fromkeys(genes))
            duplicate_count = len(genes) - len(unique_genes)
            if duplicate_count and duplicate_gene_policy == "error":
                raise ResourcePreparationError(
                    f"GMT pathway {name!r} contains {duplicate_count} duplicate "
                    "gene entries"
                )
            names.add(name)
            pathways.append(
                Pathway(
                    len(pathways), name, description, unique_genes, duplicate_count,
                )
            )
    if not pathways:
        raise ResourcePreparationError(f"GMT file contains no pathways: {path}")
    return pathways


def _read_gene_intervals(
    path: Path,
    requested_genes: set[str],
    allowed_chromosomes: set[str],
    chromosome_policy: str,
    window_bp: int,
) -> tuple[dict[str, GeneInterval], dict[str, list[GeneInterval]], set[str]]:
    genes: set[str] = set()
    requested: dict[str, GeneInterval] = {}
    requested_chromosomes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) != 4:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has {len(fields)} fields; expected 4"
                )
            raw_chromosome, start_text, end_text, gene = fields
            chromosome = _normalize_chromosome(raw_chromosome, chromosome_policy)
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has non-integer coordinates"
                ) from exc
            if chromosome not in allowed_chromosomes:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has disallowed chromosome "
                    f"{raw_chromosome!r} after configured normalization"
                )
            if start < 1 or end < start or not gene:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has invalid coordinates or gene ID"
                )
            if gene in genes:
                raise ResourcePreparationError(
                    f"gene list repeats gene identifier {gene!r}"
                )
            genes.add(gene)
            if gene in requested_genes:
                requested[gene] = GeneInterval(
                    chromosome,
                    max(1, start - window_bp),
                    end + window_bp,
                    gene,
                )
                requested_chromosomes.add(chromosome)
    if not genes:
        raise ResourcePreparationError(f"gene list contains no genes: {path}")
    intervals: dict[str, list[GeneInterval]] = {}
    for interval in requested.values():
        intervals.setdefault(interval.chromosome, []).append(interval)
    for values in intervals.values():
        values.sort(key=lambda value: (value.start, value.end, value.gene))
    return requested, intervals, requested_chromosomes


def _pathway_memberships(pathways: list[Pathway]) -> dict[str, set[int]]:
    memberships: dict[str, set[int]] = {}
    for pathway in pathways:
        for gene in pathway.genes:
            memberships.setdefault(gene, set()).add(pathway.index)
    return memberships


def _create_mapping_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute(
        "CREATE TABLE variants (variant_id TEXT PRIMARY KEY) WITHOUT ROWID"
    )
    connection.execute(
        "CREATE TABLE memberships ("
        "pathway_index INTEGER NOT NULL, variant_id TEXT NOT NULL, "
        "variant_order INTEGER NOT NULL, "
        "PRIMARY KEY (pathway_index, variant_id)) WITHOUT ROWID"
    )
    connection.execute(
        "CREATE TABLE gene_memberships ("
        "pathway_index INTEGER NOT NULL, gene TEXT NOT NULL, "
        "variant_id TEXT NOT NULL, variant_order INTEGER NOT NULL, "
        "PRIMARY KEY (pathway_index, gene, variant_id)) WITHOUT ROWID"
    )
    return connection


def _map_bim_variants(
    bim: Path,
    intervals_by_chromosome: dict[str, list[GeneInterval]],
    gene_memberships: dict[str, set[int]],
    allowed_chromosomes: set[str],
    chromosome_policy: str,
    connection: sqlite3.Connection,
) -> dict:
    current_chromosome: str | None = None
    completed_chromosomes: set[str] = set()
    active: list[GeneInterval] = []
    next_interval = 0
    last_position = 0
    variant_count = 0
    mapped_variant_count = 0
    overlapping_memberships_collapsed = 0
    mapped_genes: set[str] = set()
    observed_chromosomes: set[str] = set()
    with bim.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                raise ResourcePreparationError(
                    f"PLINK BIM contains a blank line at {line_number}"
                )
            fields = raw.split()
            if len(fields) != 6:
                raise ResourcePreparationError(
                    f"PLINK BIM line {line_number} has {len(fields)} fields; expected 6"
                )
            raw_chromosome, variant_id, _, position_text, _, _ = fields
            chromosome = _normalize_chromosome(raw_chromosome, chromosome_policy)
            if chromosome not in allowed_chromosomes:
                raise ResourcePreparationError(
                    f"PLINK BIM line {line_number} has disallowed chromosome "
                    f"{raw_chromosome!r} after configured normalization"
                )
            if not variant_id or variant_id in {".", "END"}:
                raise ResourcePreparationError(
                    f"PLINK BIM line {line_number} has an empty, '.', or reserved "
                    "END variant ID"
                )
            try:
                position = int(position_text)
            except ValueError as exc:
                raise ResourcePreparationError(
                    f"PLINK BIM line {line_number} has a non-integer position"
                ) from exc
            if position < 1:
                raise ResourcePreparationError(
                    f"PLINK BIM line {line_number} has a non-positive position"
                )
            if chromosome != current_chromosome:
                if chromosome in completed_chromosomes:
                    raise ResourcePreparationError(
                        "PLINK BIM chromosome blocks are not contiguous; sort the BIM "
                        f"before mapping (repeated chromosome {raw_chromosome!r})"
                    )
                if current_chromosome is not None:
                    completed_chromosomes.add(current_chromosome)
                current_chromosome = chromosome
                active = []
                next_interval = 0
                last_position = 0
            if position < last_position:
                raise ResourcePreparationError(
                    f"PLINK BIM positions decrease within chromosome {raw_chromosome!r} "
                    f"at line {line_number}; sort the BIM before mapping"
                )
            last_position = position
            observed_chromosomes.add(chromosome)
            try:
                connection.execute(
                    "INSERT INTO variants (variant_id) VALUES (?)", (variant_id,),
                )
            except sqlite3.IntegrityError as exc:
                raise ResourcePreparationError(
                    f"PLINK BIM repeats variant ID {variant_id!r} at line {line_number}"
                ) from exc
            variant_count += 1

            chromosome_intervals = intervals_by_chromosome.get(chromosome, [])
            while (
                next_interval < len(chromosome_intervals)
                and chromosome_intervals[next_interval].start <= position
            ):
                active.append(chromosome_intervals[next_interval])
                next_interval += 1
            active = [interval for interval in active if interval.end >= position]
            gene_pathway_pairs: set[tuple[int, str]] = set()
            for interval in active:
                mapped_genes.add(interval.gene)
                gene_pathway_pairs.update(
                    (pathway_index, interval.gene)
                    for pathway_index in gene_memberships[interval.gene]
                )
            pathway_indices = {
                pathway_index for pathway_index, _ in gene_pathway_pairs
            }
            if pathway_indices:
                mapped_variant_count += 1
                overlapping_memberships_collapsed += (
                    len(gene_pathway_pairs) - len(pathway_indices)
                )
                connection.executemany(
                    "INSERT OR IGNORE INTO memberships "
                    "(pathway_index, variant_id, variant_order) VALUES (?, ?, ?)",
                    (
                        (pathway_index, variant_id, variant_count)
                        for pathway_index in pathway_indices
                    ),
                )
                connection.executemany(
                    "INSERT INTO gene_memberships "
                    "(pathway_index, gene, variant_id, variant_order) "
                    "VALUES (?, ?, ?, ?)",
                    (
                        (pathway_index, gene, variant_id, variant_count)
                        for pathway_index, gene in gene_pathway_pairs
                    ),
                )
    if variant_count == 0:
        raise ResourcePreparationError(f"PLINK BIM contains no variants: {bim}")
    connection.commit()
    return {
        "bim_variants": variant_count,
        "bim_variants_mapped_to_requested_genes": mapped_variant_count,
        "overlapping_gene_variant_memberships_collapsed": (
            overlapping_memberships_collapsed
        ),
        "genes_with_at_least_one_bim_variant": mapped_genes,
        "bim_chromosomes": observed_chromosomes,
    }


def _table_text(fieldnames: list[str], rows: list[dict]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _readme_text(manifest: dict) -> str:
    resource = manifest["resource"]
    validation = manifest["validation"]
    outputs = manifest["outputs"]
    return f"""# {resource['name']}

This directory contains a GCTA `--fastBAT-set-list` generated from a GMT
gene-set file, a build-matched gene-coordinate list, and a PLINK BIM LD
reference. Variant identifiers are copied verbatim from BIM column 2; they are
never inferred or converted to rsIDs.

Genome build: `{resource['genome_build']}`. The GMT mapping coordinates, BIM
reference, harmonised GWAS, and analysis configuration must use this same
build. A {resource['gene_window_kb']} kb window was applied to both sides of
each gene during gene-to-variant mapping.

Pathways written: {validation['pathways_written']}  
Pathways omitted: {validation['pathways_omitted']}  
GMT genes not found in the coordinate list: {validation['unmapped_gene_entries']}  
Unique BIM variants written across retained pathways: {validation['unique_variants_written']}

`{outputs['pathway_mapping']}` reports gene and variant counts for every input
pathway; `{outputs['gene_variant_mapping']}` records the exact gene-to-BIM-ID
assignments; and `{outputs['unmapped_genes']}` records each unresolved GMT gene.
Full provenance and checksums are in `{outputs['manifest']}` and
`{outputs['checksums']}`.

This resource is for fastBAT custom-set analysis. GCTA does not document an
mBAT-combo option equivalent to `--fastBAT-set-list`.

Scientific format source:

- {OFFICIAL_FASTBAT_DOCUMENTATION}
"""


def prepare_resource(args: argparse.Namespace) -> dict:
    gmt = required_file(args.gmt, "GMT pathway file")
    gene_list = required_file(args.gene_list, "GCTA gene-coordinate file")
    bim = required_file(args.bim, "PLINK BIM file")
    output_directory = Path(args.output_directory).expanduser().resolve()
    metadata_names = {
        args.pathway_mapping_name,
        args.gene_variant_mapping_name,
        args.unmapped_genes_name,
        args.manifest_name,
        args.readme_name,
        args.checksums_name,
    }
    validate_output_filename(args.output_name, metadata_names)
    for name in metadata_names:
        validate_output_filename(name, {args.output_name})
    if len(metadata_names) != 6:
        raise ResourcePreparationError("resource output filenames must be unique")
    if args.gene_window_kb < 0:
        raise ResourcePreparationError("--gene-window-kb must be nonnegative")
    allowed_chromosomes = unique_nonempty_values(
        args.allowed_chromosomes, "--allowed-chromosomes",
    )

    pathways = _parse_gmt(gmt, args.duplicate_gene_policy)
    gene_memberships = _pathway_memberships(pathways)
    requested_genes = set(gene_memberships)
    window_bp = args.gene_window_kb * 1000
    gene_coordinates, intervals, gene_chromosomes = _read_gene_intervals(
        gene_list,
        requested_genes,
        allowed_chromosomes,
        args.chromosome_label_policy,
        window_bp,
    )
    missing_genes = requested_genes - set(gene_coordinates)
    if missing_genes and args.unmapped_gene_policy == "error":
        examples = ", ".join(sorted(missing_genes)[:10])
        raise ResourcePreparationError(
            f"{len(missing_genes)} GMT genes are absent from the coordinate list; "
            f"examples: {examples}"
        )

    staging = create_staging_directory(output_directory)
    database_handle = tempfile.NamedTemporaryFile(dir=staging, delete=False)
    database = Path(database_handle.name)
    database_handle.close()
    database.unlink()
    connection: sqlite3.Connection | None = None
    try:
        connection = _create_mapping_database(database)
        bim_metrics = _map_bim_variants(
            bim,
            intervals,
            gene_memberships,
            allowed_chromosomes,
            args.chromosome_label_policy,
            connection,
        )
        shared_chromosomes = gene_chromosomes & bim_metrics["bim_chromosomes"]
        if not shared_chromosomes:
            raise ResourcePreparationError(
                "gene coordinates and PLINK BIM do not share chromosome labels"
            )

        mapping_rows: list[dict] = []
        unmapped_rows: list[dict] = []
        retained: list[Pathway] = []
        empty: list[str] = []
        total_memberships = 0
        genes_with_variants = bim_metrics["genes_with_at_least_one_bim_variant"]
        for pathway in pathways:
            mapped_genes = [gene for gene in pathway.genes if gene in gene_coordinates]
            missing = [gene for gene in pathway.genes if gene not in gene_coordinates]
            hit_genes = [gene for gene in mapped_genes if gene in genes_with_variants]
            variant_rows = connection.execute(
                "SELECT variant_id FROM memberships WHERE pathway_index = ? "
                "ORDER BY variant_order",
                (pathway.index,),
            ).fetchall()
            variant_ids = [row[0] for row in variant_rows]
            status = "retained" if variant_ids else "omitted_empty"
            if variant_ids:
                retained.append(pathway)
                total_memberships += len(variant_ids)
            else:
                empty.append(pathway.name)
            mapping_rows.append({
                "pathway": pathway.name,
                "description": pathway.description,
                "input_genes": len(pathway.genes),
                "duplicate_gene_entries_removed": pathway.duplicate_gene_entries,
                "genes_with_coordinates": len(mapped_genes),
                "genes_with_bim_variants": len(hit_genes),
                "unmapped_genes": len(missing),
                "unique_bim_variants": len(variant_ids),
                "status": status,
            })
            unmapped_rows.extend(
                {"pathway": pathway.name, "gene": gene} for gene in missing
            )
        if empty and args.empty_pathway_policy == "error":
            examples = ", ".join(empty[:10])
            raise ResourcePreparationError(
                f"{len(empty)} pathways have no mapped BIM variants; examples: "
                f"{examples}"
            )
        if not retained:
            raise ResourcePreparationError("no pathways retain mapped BIM variants")

        output_file = staging / args.output_name
        with output_file.open("w", encoding="utf-8") as handle:
            for pathway in retained:
                handle.write(pathway.name + "\n")
                rows = connection.execute(
                    "SELECT variant_id FROM memberships WHERE pathway_index = ? "
                    "ORDER BY variant_order",
                    (pathway.index,),
                )
                for (variant_id,) in rows:
                    handle.write(variant_id + "\n")
                handle.write("END\n\n")
        mapping_file = staging / args.pathway_mapping_name
        gene_variant_file = staging / args.gene_variant_mapping_name
        unmapped_file = staging / args.unmapped_genes_name
        write_text(
            mapping_file,
            _table_text(
                [
                    "pathway", "description", "input_genes",
                    "duplicate_gene_entries_removed", "genes_with_coordinates",
                    "genes_with_bim_variants", "unmapped_genes",
                    "unique_bim_variants", "status",
                ],
                mapping_rows,
            ),
        )
        write_text(
            unmapped_file,
            _table_text(["pathway", "gene"], unmapped_rows),
        )
        with gene_variant_file.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["pathway", "gene", "variant_id"])
            rows = connection.execute(
                "SELECT pathway_index, gene, variant_id FROM gene_memberships "
                "ORDER BY pathway_index, variant_order, gene"
            )
            for pathway_index, gene, variant_id in rows:
                writer.writerow([pathways[pathway_index].name, gene, variant_id])
        unique_variants_written = connection.execute(
            "SELECT COUNT(DISTINCT variant_id) FROM memberships"
        ).fetchone()[0]
        manifest = {
            "schema_version": "gcta_fastbat_pathway_resource.v1",
            "resource": {
                "name": args.resource_name,
                "genome_build": args.genome_build,
                "method": "fastbat_set",
                "gcta_option": "--fastBAT-set-list",
                "file": output_file.name,
                "format": "GCTA set blocks: set ID, exact BIM IDs, END",
                "gene_window_kb": args.gene_window_kb,
                "sha256": sha256(output_file),
            },
            "sources": {
                "gmt": {"path": str(gmt), "sha256": sha256(gmt)},
                "gene_coordinates": {
                    "path": str(gene_list), "sha256": sha256(gene_list),
                },
                "plink_bim": {"path": str(bim), "sha256": sha256(bim)},
            },
            "policies": {
                "allowed_chromosomes": sorted(allowed_chromosomes),
                "chromosome_label_policy": args.chromosome_label_policy,
                "duplicate_gene_policy": args.duplicate_gene_policy,
                "unmapped_gene_policy": args.unmapped_gene_policy,
                "empty_pathway_policy": args.empty_pathway_policy,
                "variant_identifier_policy": (
                    "copy PLINK BIM column 2 verbatim; never infer rsIDs"
                ),
            },
            "validation": {
                "input_pathways": len(pathways),
                "pathways_written": len(retained),
                "pathways_omitted": len(empty),
                "input_unique_genes": len(requested_genes),
                "genes_with_coordinates": len(gene_coordinates),
                "unmapped_unique_genes": len(missing_genes),
                "unmapped_gene_entries": len(unmapped_rows),
                "bim_variants": bim_metrics["bim_variants"],
                "bim_variants_mapped_to_requested_genes": (
                    bim_metrics["bim_variants_mapped_to_requested_genes"]
                ),
                "overlapping_gene_variant_memberships_collapsed": (
                    bim_metrics["overlapping_gene_variant_memberships_collapsed"]
                ),
                "total_pathway_variant_memberships": total_memberships,
                "unique_variants_written": unique_variants_written,
                "shared_chromosomes": sorted(shared_chromosomes),
            },
            "outputs": {
                "pathway_mapping": mapping_file.name,
                "gene_variant_mapping": gene_variant_file.name,
                "unmapped_genes": unmapped_file.name,
                "manifest": args.manifest_name,
                "readme": args.readme_name,
                "checksums": args.checksums_name,
            },
            "generation": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "command": getattr(args, "generation_command", None) or shlex.join(sys.argv),
                "utility": str(Path(__file__).resolve()),
            },
            "scientific_sources": [OFFICIAL_FASTBAT_DOCUMENTATION],
        }
        manifest_file = staging / args.manifest_name
        readme_file = staging / args.readme_name
        checksum_file = staging / args.checksums_name
        write_text(
            manifest_file,
            yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True),
        )
        write_text(readme_file, _readme_text(manifest))
        checksummed = [
            output_file, mapping_file, gene_variant_file, unmapped_file,
            manifest_file, readme_file,
        ]
        write_text(
            checksum_file,
            "".join(f"{sha256(path)}  {path.name}\n" for path in checksummed),
        )
        connection.close()
        connection = None
        database.unlink()
        staging.replace(output_directory)
    except BaseException:
        if connection is not None:
            connection.close()
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert GMT gene pathways to GCTA fastBAT SNP-set blocks using "
            "build-matched gene coordinates and exact PLINK BIM variant IDs."
        )
    )
    parser.add_argument("--gmt", required=True, type=Path, metavar="PATH")
    parser.add_argument("--gene-list", required=True, type=Path, metavar="PATH")
    parser.add_argument("--bim", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-directory", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-name", required=True, metavar="NAME")
    parser.add_argument("--resource-name", required=True, metavar="NAME")
    parser.add_argument("--genome-build", required=True, metavar="BUILD")
    parser.add_argument("--gene-window-kb", required=True, type=int, metavar="KB")
    parser.add_argument(
        "--allowed-chromosomes", required=True, nargs="+", metavar="CHR",
    )
    parser.add_argument(
        "--chromosome-label-policy",
        required=True,
        choices=("exact", "strip_chr_prefix"),
    )
    parser.add_argument(
        "--duplicate-gene-policy", required=True, choices=("error", "deduplicate"),
    )
    parser.add_argument(
        "--unmapped-gene-policy", required=True, choices=("error", "report"),
    )
    parser.add_argument(
        "--empty-pathway-policy", required=True, choices=("error", "omit"),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    from postgwas.config import load_configuration

    names = load_configuration().modules.gcta_gene.set_annotation.conversion.output_names
    args.pathway_mapping_name = names.pathway_mapping
    args.gene_variant_mapping_name = names.gene_variant_mapping
    args.unmapped_genes_name = names.unmapped_genes
    args.manifest_name = names.manifest
    args.readme_name = names.readme
    args.checksums_name = names.checksums
    for value in (args.output_name, args.resource_name, args.genome_build):
        if not str(value).strip():
            parser.error("output name, resource name, and genome build must be non-empty")
    try:
        manifest = prepare_resource(args)
    except (OSError, ResourcePreparationError, sqlite3.Error, yaml.YAMLError) as exc:
        parser.error(str(exc))
    print(
        "Created %s (%d pathways, %d unique BIM variants, %s)"
        % (
            Path(args.output_directory).expanduser().resolve() / args.output_name,
            manifest["validation"]["pathways_written"],
            manifest["validation"]["unique_variants_written"],
            manifest["resource"]["genome_build"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
