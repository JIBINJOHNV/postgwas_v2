#!/usr/bin/env python3
"""
annot_ldblock.py — Annotate GWAS VCFs with population-specific LD blocks.
"""

from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

from postgwas.core.input_validation import record_file_validation
from postgwas.core.interval_validation import validate_bed4_file, validate_vcf_bed_contigs
from postgwas.core.io.artifacts import publish_artifact_set
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
)
from postgwas.core.processes import run_checked_command, run_checked_pipeline
from postgwas.core.ui import StageProgress
from postgwas.core.vcf import (
    count_indexed_vcf_records,
    declared_vcf_tag_definitions,
    read_vcf_header,
    select_vcf_sample,
    validate_postgwas_vcf_provenance,
    validate_vcf_header_contract,
)
from postgwas.modules.ld_annotation.reporting import (
    LDReferenceSummary,
    calculate_ld_annotation_summary,
    write_ld_annotation_html_report,
    write_ld_annotation_summary,
)




def validate_ld_block_references(
    directory: str | Path,
    *,
    genome_build: str,
    populations: Sequence[str],
    bed_filename_template: str,
    vcf_contigs: Sequence[str],
) -> dict[str, LDReferenceSummary]:
    """Share exact resource selection and validation between preflight and use."""
    references = {}
    for population in populations:
        filename = bed_filename_template.format(
            genome_build=genome_build, population=population,
        )
        bed = require_nonempty_file(
            Path(directory).expanduser().resolve() / filename,
            "LD block BED file missing or empty",
            error_type=FileNotFoundError,
        )
        reference = validate_bed4_file(bed)
        try:
            validate_vcf_bed_contigs(vcf_contigs, reference.contigs, bed)
        except ValueError as exc:
            record_file_validation(
                bed, "LD-block chromosome compatibility", status="failed",
                message=str(exc),
            )
            raise
        record_file_validation(
            bed, "LD-block chromosome compatibility",
            checks=("exact chromosome naming against input VCF declarations",),
            metrics={"declared_genome_build": genome_build, "population": population},
            message="Build and population select the configured resource; BED does not independently prove either declaration.",
        )
        references[population] = LDReferenceSummary(
            path=reference.path,
            contigs=reference.contigs,
            labels=reference.labels,
            block_count=reference.block_count,
        )
    return references


def _validate_annotated_output(
    vcf: Path,
    *,
    bcftools: str,
    genome_build: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    expected_contigs: Sequence[str],
    info_ids: Sequence[str],
    logger=None,
) -> int:
    """Validate final VCF metadata and prove its TBI can count all records."""
    header = read_vcf_header(vcf, bcftools, logger=logger)
    observed_build, observed_contigs = validate_vcf_header_contract(
        header=header,
        genome_build_header=genome_build_header,
        supported_genome_builds=supported_genome_builds,
        required_fields=["INFO/%s" % info_id for info_id in info_ids],
    )
    if observed_build != genome_build:
        raise ValueError(
            "Annotated VCF declares genome build %s; expected %s"
            % (observed_build, genome_build)
        )
    if list(observed_contigs) != list(expected_contigs):
        raise ValueError(
            "Annotated VCF contig declarations differ from the input VCF"
        )
    definitions = declared_vcf_tag_definitions(header, "INFO")
    invalid = [
        info_id
        for info_id in info_ids
        if definitions.get(info_id) != {"number": "1", "type": "String"}
    ]
    if invalid:
        raise ValueError(
            "Annotated VCF has invalid LD-block INFO definitions: %s"
            % ", ".join(invalid)
        )
    count_text = run_checked_command(
        [bcftools, "index", "-n", str(vcf)],
        "Validating annotated VCF index",
        logger=logger,
    ).strip()
    try:
        record_count = int(count_text)
    except ValueError as exc:
        raise ValueError(
            "bcftools returned an invalid annotated-VCF record count: %r"
            % count_text
        ) from exc
    if record_count < 0:
        raise ValueError("Annotated VCF record count cannot be negative")
    return record_count


@contextmanager
def _annotation_stage(
    progress: StageProgress,
    logger,
    number: int,
    total: int,
    title: str,
    function_name: str,
):
    """Use the canonical logger when present and retain progress-only support."""
    if logger is not None:
        with logger.step(number, total, title, function_name) as step:
            yield step
        return
    with progress.step(number, total, title) as step:
        yield step


# ------------------------------------------------------------
# Main Function
# ------------------------------------------------------------
def annotate_ldblocks(
    vcf_path: str,
    output_directory: str,
    ld_dir: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    bed_filename_template: str,
    info_field_template: str,
    info_description_template: str,
    output_filename_template: str,
    summary_filename_template: str,
    html_report_filename_template: str,
    bcftools_bin: str,
    *,
    provenance_headers: Mapping[str, str],
    populations: Sequence[str],
    threads: int,
    dataset_id: str,
    logger=None,
    stage_progress: StageProgress | None = None,
) -> dict:
    """Annotate population LD blocks through one checked bcftools pipeline."""
    progress = stage_progress or StageProgress(
        "LD-block annotation progress",
        enabled=False,
    )
    total_stages = 6
    with _annotation_stage(
        progress,
        logger,
        1,
        total_stages,
        "Validate input VCF metadata",
        "validate_vcf_header_contract",
    ) as step:
        # Dependency, argument, and VCF-contract checks fail before any output
        # directory or scientific artifact is created.
        bcftools = resolve_executable(bcftools_bin, "bcftools")
        if int(threads) < 1:
            raise ValueError("LD annotation threads must be at least 1")
        selected_populations = tuple(
            getattr(population, "value", str(population))
            for population in populations
        )
        if not selected_populations or len(selected_populations) != len(
            set(selected_populations)
        ):
            raise ValueError(
                "LD annotation requires one or more unique populations"
            )

        current_vcf = Path(vcf_path).expanduser().resolve()
        if not current_vcf.is_file() or current_vcf.stat().st_size <= 0:
            raise FileNotFoundError(
                "Input VCF does not exist or is empty: %s" % current_vcf
            )

        # The VCF header is the authoritative genome-build declaration for this
        # operation. BED coordinates do not contain enough information to infer a
        # build, so the declared VCF build selects the expected reference filename.
        header = read_vcf_header(current_vcf, bcftools, logger=logger)
        validate_postgwas_vcf_provenance(
            header, provenance_headers, vcf_path=current_vcf, logger=logger,
        )
        genome_build, vcf_contigs = validate_vcf_header_contract(
            header=header,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            required_fields=(),
        )
        if not vcf_contigs:
            raise ValueError(
                "Input VCF header contains no ##contig declarations; LD-block "
                "BED chromosome compatibility cannot be validated."
            )
        input_variant_count = count_indexed_vcf_records(
            current_vcf,
            bcftools,
            logger=logger,
        )
        select_vcf_sample(
            current_vcf,
            dataset_id,
            bcftools,
            logger=logger,
        )
        step.outcome(
            "Validated the input VCF metadata, single-sample contract, and "
            "indexed record count.",
            fields=(
                ("genetic", "Genome build", genome_build),
                ("count", "Declared contigs", len(vcf_contigs)),
                ("count", "Total variants", input_variant_count),
                ("count", "VCF sample columns", 1),
                ("info", "Input VCF", current_vcf.name),
            ),
            input_vcf=str(current_vcf),
            genome_build=genome_build,
            declared_contigs=len(vcf_contigs),
            total_variants=input_variant_count,
            sample_columns=1,
        )

    # Preflight every selected population before creating output. BED files are
    # intentionally read directly because bcftools supports them without a
    # tabix index and LDetect resources contain only a small number of intervals.
    beds: dict[str, Path] = {}
    references: dict[str, LDReferenceSummary] = {}
    info_ids: dict[str, str] = {}
    with _annotation_stage(
        progress,
        logger,
        2,
        total_stages,
        "Validate population LD-block references",
        "validate_ld_block_references",
    ) as step:
        references = validate_ld_block_references(
            ld_dir,
            genome_build=genome_build,
            populations=selected_populations,
            bed_filename_template=bed_filename_template,
            vcf_contigs=vcf_contigs,
        )
        for population in selected_populations:
            beds[population] = references[population].path
            info_ids[population] = info_field_template.format(
                population=population
            )
        reference_fields: list[tuple] = [
            (
                "analysis",
                "Populations requested",
                ", ".join(selected_populations),
            ),
            (
                "success",
                "Required references",
                "%d/%d available"
                % (len(references), len(selected_populations)),
            ),
            ("success", "Chromosome naming", "exact match"),
        ]
        for population in selected_populations:
            reference = references[population]
            reference_fields.extend(
                (
                    ("genetic", "%s reference" % population),
                    ("count", "LD blocks", reference.block_count),
                    ("count", "Reference contigs", len(reference.contigs)),
                    ("info", "BED file", reference.path.name),
                )
            )
        step.outcome(
            "Validated every requested population LD-block reference.",
            fields=reference_fields,
            references={
                population: {
                    "bed_file": str(references[population].path),
                    "blocks": references[population].block_count,
                    "unique_labels": len(references[population].labels),
                    "contigs": len(references[population].contigs),
                }
                for population in selected_populations
            },
            populations=selected_populations,
            chromosome_naming="exact_match",
        )

    temporary_vcf: Path | None = None
    temporary_index: Path | None = None
    temporary_summary: Path | None = None
    temporary_html: Path | None = None
    try:
        with _annotation_stage(
            progress,
            logger,
            3,
            total_stages,
            "Annotate population LD blocks",
            "run_checked_pipeline",
        ) as step:
            step_output_directory = Path(output_directory).expanduser().resolve()
            step_output_directory.mkdir(parents=True, exist_ok=True)
            annotated_vcf = configured_output_path(
                step_output_directory,
                output_filename_template,
                dataset_id=dataset_id,
            )
            annotated_index = Path(str(annotated_vcf) + ".tbi")
            summary_file = configured_output_path(
                step_output_directory,
                summary_filename_template,
                dataset_id=dataset_id,
            )
            html_report = configured_output_path(
                step_output_directory,
                html_report_filename_template,
                dataset_id=dataset_id,
            )
            temporary_token = uuid4().hex
            temporary_vcf = annotated_vcf.with_name(
                ".%s.%s.tmp.vcf.gz" % (annotated_vcf.name, temporary_token)
            )
            temporary_index = Path(str(temporary_vcf) + ".tbi")
            temporary_summary = summary_file.with_name(
                ".%s.%s.tmp.csv" % (summary_file.name, temporary_token)
            )
            temporary_html = html_report.with_name(
                ".%s.%s.tmp.html" % (html_report.name, temporary_token)
            )

            commands: list[list[str]] = []
            for index, population in enumerate(selected_populations):
                info_id = info_ids[population]
                info_description = info_description_template.format(
                    population=population,
                    genome_build=genome_build,
                )
                # Each BED interval supplies one scalar label. Number=1/Type=String
                # and CHROM,FROM,TO are bcftools protocol invariants protected by
                # BED-contract and final-header validation tests.
                header_line = (
                    '##INFO=<ID=%s,Number=1,Type=String,Description="%s">'
                    % (info_id, info_description)
                )
                command = [
                    bcftools,
                    "annotate",
                    "-a",
                    str(beds[population]),
                    "-c",
                    "CHROM,FROM,TO,%s" % info_id,
                    "-H",
                    header_line,
                ]
                source = str(current_vcf) if index == 0 else "-"
                if index == len(selected_populations) - 1:
                    # TBI is the published index contract consumed downstream;
                    # only the final compressed stream needs an index.
                    command.extend(
                        [
                            "--threads",
                            str(threads),
                            "-Oz",
                            "--write-index=tbi",
                            "-o",
                            str(temporary_vcf),
                            source,
                        ]
                    )
                else:
                    command.extend(["-Ou", source])
                commands.append(command)

            run_checked_pipeline(
                commands,
                "Annotating population LD blocks",
                logger=logger,
                expected_outputs=(temporary_vcf, temporary_index),
            )
            step.outcome(
                "Added one independent INFO annotation for each population.",
                fields=(
                    (
                        "analysis",
                        "Populations annotated",
                        ", ".join(selected_populations),
                    ),
                    (
                        "success",
                        "INFO fields created",
                        ", ".join(info_ids.values()),
                    ),
                    ("info", "Population stages", len(commands)),
                    ("info", "Compression threads", int(threads)),
                ),
                populations=selected_populations,
                info_fields=tuple(info_ids.values()),
                population_stages=len(commands),
                compression_threads=int(threads),
                output_vcf=str(temporary_vcf),
            )

        with _annotation_stage(
            progress,
            logger,
            4,
            total_stages,
            "Validate annotated VCF and index",
            "validate_annotated_output",
        ) as step:
            variant_count = _validate_annotated_output(
                temporary_vcf,
                bcftools=bcftools,
                genome_build=genome_build,
                genome_build_header=genome_build_header,
                supported_genome_builds=supported_genome_builds,
                expected_contigs=vcf_contigs,
                info_ids=tuple(info_ids.values()),
                logger=logger,
            )
            if variant_count != input_variant_count:
                raise ValueError(
                    "LD-block annotation changed the VCF record count: input "
                    "VCF has %d variants but the annotated VCF has %d. "
                    "Annotation must retain every input variant."
                    % (input_variant_count, variant_count)
                )
            step.outcome(
                "Validated the annotated VCF metadata, INFO fields, and index.",
                fields=(
                    ("count", "VCF variants validated", variant_count),
                    ("success", "Variant count preserved", "yes"),
                    (
                        "success",
                        "INFO fields validated",
                        ", ".join(info_ids.values()),
                    ),
                    ("success", "Genome build preserved", genome_build),
                    ("success", "TBI index", "valid"),
                ),
                variants=variant_count,
                input_variants=input_variant_count,
                variant_count_preserved=True,
                info_fields=tuple(info_ids.values()),
                genome_build=genome_build,
                tbi_index=str(temporary_index),
            )

        with _annotation_stage(
            progress,
            logger,
            5,
            total_stages,
            "Calculate and save annotation reports",
            "calculate_ld_annotation_summary",
        ) as step:
            summary = calculate_ld_annotation_summary(
                temporary_vcf,
                bcftools=bcftools,
                dataset_id=dataset_id,
                genome_build=genome_build,
                populations=selected_populations,
                info_ids=info_ids,
                references=references,
                expected_variant_count=variant_count,
                vcf_contigs=vcf_contigs,
                logger=logger,
            )
            write_ld_annotation_summary(summary, temporary_summary)
            report_outputs = {
                "annotated_vcf": annotated_vcf,
                "annotated_index": annotated_index,
                "summary_file": summary_file,
                "html_report": html_report,
                "log_file": (
                    None if logger is None else Path(logger.log_path)
                ),
            }
            write_ld_annotation_html_report(
                summary,
                temporary_html,
                input_vcf=current_vcf,
                vcf_contigs=vcf_contigs,
                info_ids=info_ids,
                references=references,
                threads=int(threads),
                bcftools=bcftools,
                outputs=report_outputs,
            )
            step.outcome(
                "Calculated annotation coverage and rendered the CSV and HTML reports.",
                fields=(
                    (
                        "success",
                        "Coverage summary",
                        "calculated for %d population%s"
                        % (
                            len(selected_populations),
                            "" if len(selected_populations) == 1 else "s",
                        ),
                    ),
                    ("success", "Reports prepared", "CSV and HTML"),
                ),
                total_variants=summary["total_variants"],
                any_population_annotated=summary[
                    "any_population_annotated"
                ],
                all_populations_annotated=summary[
                    "all_populations_annotated"
                ],
                fully_unassigned_variants=summary[
                    "fully_unassigned_variants"
                ],
                unassigned_on_bed_contigs=summary[
                    "unassigned_on_bed_contigs"
                ],
                population_coverage=summary["populations"],
                summary_file=str(temporary_summary),
                html_report=str(temporary_html),
            )

        with _annotation_stage(
            progress,
            logger,
            6,
            total_stages,
            "Publish validated VCF, index, and reports",
            "publish_artifact_set",
        ) as step:
            artifact_set = (
                (temporary_vcf, annotated_vcf),
                (temporary_index, annotated_index),
                (temporary_summary, summary_file),
                (temporary_html, html_report),
            )
            publish_artifact_set(artifact_set)
            published_count = len(artifact_set)
            step.outcome(
                "Published the validated LD-annotation artifact set.",
                fields=(
                    (
                        "success",
                        "Analysis artifacts published",
                        "%d/%d" % (published_count, published_count),
                    ),
                ),
                annotated_vcf=str(annotated_vcf),
                annotated_index=str(annotated_index),
                summary_file=str(summary_file),
                html_report=str(html_report),
                log_file=(
                    None if logger is None else str(Path(logger.log_path))
                ),
            )
    finally:
        if temporary_vcf is not None:
            temporary_vcf.unlink(missing_ok=True)
        if temporary_index is not None:
            temporary_index.unlink(missing_ok=True)
        if temporary_summary is not None:
            temporary_summary.unlink(missing_ok=True)
        if temporary_html is not None:
            temporary_html.unlink(missing_ok=True)
    return {
        "annotated_vcf": annotated_vcf,
        "annotated_index": annotated_index,
        "summary_file": summary_file,
        "html_report": html_report,
        "summary": summary,
        "variant_count": variant_count,
    }
