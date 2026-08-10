#!/usr/bin/env python3
"""Create a validated GCTA fastBAT/mBAT-combo gene list from a table.

This is a resource-preparation utility, not part of the runtime scientific
module. Column positions and resource metadata are explicit command-line
inputs so the utility does not infer a genome build or a source schema.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import shlex
import shutil
import sys
from typing import Iterable

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
OFFICIAL_MBAT_COMBO_DOCUMENTATION = (
    "https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo"
)


def _validate_column_positions(
    positions: Iterable[int], expected_source_columns: int,
) -> tuple[int, ...]:
    values = tuple(positions)
    if expected_source_columns < 1:
        raise ResourcePreparationError("--expected-source-columns must be positive")
    if len(set(values)) != len(values):
        raise ResourcePreparationError("source column positions must be unique")
    if any(value < 1 or value > expected_source_columns for value in values):
        raise ResourcePreparationError(
            "source column positions must be between 1 and "
            f"{expected_source_columns}"
        )
    return tuple(value - 1 for value in values)


def convert_gene_list(
    source: Path,
    *,
    expected_source_columns: int,
    chromosome_column: int,
    start_column: int,
    end_column: int,
    gene_column: int,
    allowed_chromosomes: set[str],
) -> tuple[str, dict]:
    """Return a headerless GCTA table and validation metrics.

    Coordinates are copied exactly. GCTA applies any requested flanking window
    later through ``--fastBAT-wind`` or ``--mBAT-wind``.
    """
    chromosome_index, start_index, end_index, gene_index = _validate_column_positions(
        (chromosome_column, start_column, end_column, gene_column),
        expected_source_columns,
    )
    if not allowed_chromosomes or any(not value.strip() for value in allowed_chromosomes):
        raise ResourcePreparationError(
            "--allowed-chromosomes must contain non-empty chromosome labels"
        )

    records: list[tuple[str, int, int, str]] = []
    genes: set[str] = set()
    record_keys: set[tuple[str, int, int, str]] = set()
    chromosomes: Counter[str] = Counter()
    blank_lines = 0
    with source.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                blank_lines += 1
                continue
            fields = raw.split()
            if len(fields) != expected_source_columns:
                raise ResourcePreparationError(
                    f"source line {line_number} has {len(fields)} fields; "
                    f"expected {expected_source_columns}"
                )
            chromosome = fields[chromosome_index]
            gene = fields[gene_index]
            try:
                start = int(fields[start_index])
                end = int(fields[end_index])
            except ValueError as exc:
                raise ResourcePreparationError(
                    f"source line {line_number} has non-integer coordinates"
                ) from exc
            if chromosome not in allowed_chromosomes:
                raise ResourcePreparationError(
                    f"source line {line_number} has disallowed chromosome "
                    f"{chromosome!r}"
                )
            if start < 1 or end < start:
                raise ResourcePreparationError(
                    f"source line {line_number} has invalid coordinates "
                    f"{start}-{end}"
                )
            if not gene:
                raise ResourcePreparationError(
                    f"source line {line_number} has an empty gene identifier"
                )
            if gene in genes:
                raise ResourcePreparationError(
                    f"source line {line_number} repeats gene identifier {gene!r}"
                )
            record = (chromosome, start, end, gene)
            if record in record_keys:
                raise ResourcePreparationError(
                    f"source line {line_number} repeats output record {record!r}"
                )
            genes.add(gene)
            record_keys.add(record)
            records.append(record)
            chromosomes[chromosome] += 1

    if not records:
        raise ResourcePreparationError(f"source contains no gene records: {source}")
    output = "".join(
        f"{chromosome}\t{start}\t{end}\t{gene}\n"
        for chromosome, start, end, gene in records
    )
    return output, {
        "source_rows": len(records),
        "output_rows": len(records),
        "excluded_rows": 0,
        "blank_source_lines": blank_lines,
        "unique_gene_identifiers": len(genes),
        "duplicate_gene_identifiers": 0,
        "duplicate_output_records": 0,
        "invalid_coordinate_rows": 0,
        "chromosome_counts": dict(sorted(chromosomes.items())),
    }


def _metadata_file(path: Path | None) -> dict | None:
    if path is None:
        return None
    source = required_file(path, "source metadata file")
    return {"path": str(source), "sha256": sha256(source)}


def _readme_text(manifest: dict) -> str:
    resource = manifest["resource"]
    source = manifest["source"]
    validation = manifest["validation"]
    outputs = manifest["outputs"]
    return f"""# {resource['name']}

This directory contains a headerless, tab-delimited GCTA gene list for
`fastbat_gene` and `mbat_combo` with columns:

```text
chromosome  start  end  gene
```

Genome build: `{resource['genome_build']}`. This build was declared explicitly
when the resource was generated; it was not inferred from coordinates. The
harmonised GWAS, PLINK LD reference, and this gene list must all use the same
genome build.

The output is an exact projection of source columns
{source['columns']['chromosome']}, {source['columns']['start']},
{source['columns']['end']}, and {source['columns']['gene']}. Coordinates were
not shifted, lifted over, expanded, or changed according to strand. Gene
windows are applied only by GCTA through `--fastBAT-wind` or `--mBAT-wind`.

Source: `{source['path']}`  
Source SHA-256: `{source['sha256']}`  
Rows written: {validation['output_rows']}  
Rows excluded: {validation['excluded_rows']}

Full provenance and validation metrics are recorded in
`{outputs['manifest']}`; file checksums are recorded in
`{outputs['checksums']}`.

Scientific format sources:

- {OFFICIAL_FASTBAT_DOCUMENTATION}
- {OFFICIAL_MBAT_COMBO_DOCUMENTATION}
"""


def prepare_resource(args: argparse.Namespace) -> dict:
    source = required_file(args.input, "source gene-location file")
    output_directory = Path(args.output_directory).expanduser().resolve()
    metadata_names = {args.manifest_name, args.readme_name, args.checksums_name}
    validate_output_filename(args.output_name, metadata_names)
    for name in metadata_names:
        validate_output_filename(name, {args.output_name})
    if len(metadata_names) != 3:
        raise ResourcePreparationError("resource output filenames must be unique")
    source_readme = _metadata_file(args.source_readme)
    source_report = _metadata_file(args.source_report)
    allowed_chromosomes = unique_nonempty_values(
        args.allowed_chromosomes, "--allowed-chromosomes",
    )

    output_text, validation = convert_gene_list(
        source,
        expected_source_columns=args.expected_source_columns,
        chromosome_column=args.chromosome_column,
        start_column=args.start_column,
        end_column=args.end_column,
        gene_column=args.gene_column,
        allowed_chromosomes=allowed_chromosomes,
    )
    staging = create_staging_directory(output_directory)
    output_file = staging / args.output_name
    manifest_file = staging / args.manifest_name
    readme_file = staging / args.readme_name
    checksum_file = staging / args.checksums_name
    try:
        write_text(output_file, output_text)
        output_sha256 = sha256(output_file)
        manifest = {
            "schema_version": "gcta_gene_list_resource.v1",
            "resource": {
                "name": args.resource_name,
                "genome_build": args.genome_build,
                "methods": ["fastbat_gene", "mbat_combo"],
                "file": output_file.name,
                "format": "headerless tab-delimited text",
                "columns": ["chromosome", "start", "end", "gene"],
                "sha256": output_sha256,
            },
            "source": {
                "label": args.source_label,
                "path": str(source),
                "sha256": sha256(source),
                "expected_columns": args.expected_source_columns,
                "columns": {
                    "chromosome": args.chromosome_column,
                    "start": args.start_column,
                    "end": args.end_column,
                    "gene": args.gene_column,
                },
                "readme": source_readme,
                "report": source_report,
            },
            "conversion": {
                "coordinate_policy": (
                    "preserved exactly; no liftover or window expansion"
                ),
                "gene_identifier_policy": (
                    "copied exactly from the configured source column"
                ),
                "genome_build_policy": (
                    "declared by caller; not inferred from coordinates"
                ),
            },
            "validation": validation,
            "outputs": {
                "manifest": manifest_file.name,
                "readme": readme_file.name,
                "checksums": checksum_file.name,
            },
            "generation": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "command": shlex.join(sys.argv),
                "utility": str(Path(__file__).resolve()),
            },
            "scientific_sources": [
                OFFICIAL_FASTBAT_DOCUMENTATION,
                OFFICIAL_MBAT_COMBO_DOCUMENTATION,
            ],
        }
        write_text(
            manifest_file,
            yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True),
        )
        write_text(readme_file, _readme_text(manifest))
        checksums = {
            output_file.name: sha256(output_file),
            manifest_file.name: sha256(manifest_file),
            readme_file.name: sha256(readme_file),
        }
        write_text(
            checksum_file,
            "".join(f"{digest}  {name}\n" for name, digest in checksums.items()),
        )
        staging.replace(output_directory)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Project a validated source gene-location table into the headerless "
            "four-column format shared by GCTA fastBAT and mBAT-combo."
        )
    )
    parser.add_argument("--input", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-directory", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-name", required=True, metavar="NAME")
    parser.add_argument("--resource-name", required=True, metavar="NAME")
    parser.add_argument("--source-label", required=True, metavar="TEXT")
    parser.add_argument("--genome-build", required=True, metavar="BUILD")
    parser.add_argument("--expected-source-columns", required=True, type=int, metavar="N")
    parser.add_argument("--chromosome-column", required=True, type=int, metavar="N")
    parser.add_argument("--start-column", required=True, type=int, metavar="N")
    parser.add_argument("--end-column", required=True, type=int, metavar="N")
    parser.add_argument("--gene-column", required=True, type=int, metavar="N")
    parser.add_argument(
        "--allowed-chromosomes", required=True, nargs="+", metavar="CHR",
    )
    parser.add_argument("--source-readme", type=Path, metavar="PATH")
    parser.add_argument("--source-report", type=Path, metavar="PATH")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    from postgwas.config import load_configuration

    names = load_configuration().modules.gcta_gene.gene_annotation.resource_output_names
    args.manifest_name = names.manifest
    args.readme_name = names.readme
    args.checksums_name = names.checksums
    for name in (args.output_name, args.resource_name, args.source_label, args.genome_build):
        if not str(name).strip():
            parser.error("resource names, labels, and genome build must be non-empty")
    try:
        manifest = prepare_resource(args)
    except (OSError, ResourcePreparationError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    print(
        "Created %s (%d genes, %s)"
        % (
            Path(args.output_directory).expanduser().resolve() / args.output_name,
            manifest["validation"]["output_rows"],
            manifest["resource"]["genome_build"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
