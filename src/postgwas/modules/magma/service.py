"""Configuration, logging and publication boundary for MAGMA."""

from __future__ import annotations

import shutil
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.checkpointing import (
    ExecutionCheckpoint,
    decode_checkpoint_value,
    discover_input_files,
    software_identity,
)
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    PreflightFileIdentity,
    PreflightLogBuffer,
    PreflightLogEvent,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    replay_preflight_log,
    require_pipeline_input_vcf,
    require_unchanged_preflight_files,
)
from postgwas.core.paths import (
    configured_output_path,
    remove_owned_directory,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.plink import validate_plink_bundle_dimensions, validate_plink_files
from postgwas.core.ui import StageProgress, print_screen_block
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.annotations import read_gene_locations
from postgwas.modules.magma.analysis import (
    MagmaReferencePreflight,
    preflight_magma_analysis,
    preflight_magma_references,
    resolve_magma_output_paths,
    run_magma_analysis,
)
from postgwas.modules.magma.reporting import (
    build_magma_pipeline_summary,
    magma_pipeline_stage_numbers,
    magma_variant_input_outcome_fields,
    render_magma_results_screen,
    write_magma_pipeline_csv,
    write_magma_pipeline_html,
)


@dataclass(frozen=True)
class MagmaPipelineResources:
    """Validated MAGMA resources retained outside checkpointed run state."""

    configuration: Any
    reference: MagmaReferencePreflight
    include_gene_sets: bool
    primary_gene_location_metrics: Mapping[str, Any] | None
    file_identities: tuple[PreflightFileIdentity, ...]
    log_events: tuple[PreflightLogEvent, ...]


def _resolved_configuration(args):
    module_overrides = explicit_overrides(
        args,
        {
            "genome_build": "genome_build",
            "snp_location_file": "input.snp_location_file",
            "p_value_file": "input.p_value_file",
            "magma_ld_reference": "input.ld_reference_prefix",
            "gene_location_file": "input.gene_location_file",
            "gene_set_file": "input.gene_set_file",
            "magma_positional_gene_id_type": (
                "mapping.definitions.positional.gene_id_type"
            ),
            "magma_positional_source_name": (
                "mapping.definitions.positional.source_name"
            ),
            "magma_positional_source_version": (
                "mapping.definitions.positional.source_version"
            ),
            "magma_positional_source_url": (
                "mapping.definitions.positional.source_url"
            ),
            "magma_positional_context": (
                "mapping.definitions.positional.context"
            ),
            "gene_location_alternate_id_type": "input.alternate_gene_id_type",
            "sample_size_column": "input.sample_size_column",
            "resolve_variants_to_reference": (
                "snp_harmonisation.resolve_variants_to_reference"
            ),
            "minimum_snp_overlap": "snp_harmonisation.minimum_overlap_fraction",
            "duplicate_policy": "snp_harmonisation.duplicate_policy",
            "minimum_gene_id_overlap": "gene_sets.minimum_gene_id_overlap_fraction",
            "gene_set_identifier_mismatch": "gene_sets.identifier_mismatch_action",
            "alternate_gene_id_duplicate_policy": (
                "gene_sets.alternate_id_duplicate_policy"
            ),
            "window_upstream": "gene_window_upstream_kb",
            "window_downstream": "gene_window_downstream_kb",
            "gene_model": "gene_model",
            "magma_mapping": "mapping.selected",
            "primary_magma_mapping": "mapping.primary",
            "mhc_policy": "mhc.policy",
            "mhc_chrom": "mhc.region_override.chromosome",
            "mhc_start": "mhc.region_override.start",
            "mhc_end": "mhc.region_override.end",
            "exclude_chromosomes": "chromosomes.exclude",
            "magma_memory_per_worker_gb": "batching.memory_per_process_gb",
        },
    )
    global_overrides = explicit_overrides(
        args,
        {
            "dataset_id": "run.dataset_id",
            "output_directory": "run.output_directory",
            "threads": "execution.threads",
            "memory_gb": "execution.memory_gb",
            "seed": "execution.random_seed",
            "magma": "resources.executables.magma",
            "resume": "run.resume",
            "overwrite": "run.overwrite",
        },
    )
    return load_run_configuration_for_module(
        "magma",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _resolved_magma_metadata_paths(module) -> tuple[str, ...]:
    """Select all used MAGMA settings while omitting unselected definitions."""
    paths = [name for name in type(module).model_fields if name != "mapping"]
    paths.extend(("mapping.selected", "mapping.primary"))
    paths.extend(
        "mapping.definitions.%s" % name for name in module.mapping.selected
    )
    return tuple(paths)


def _expected_magma_artifacts(
    output: Path,
    dataset_id: str,
    module,
    *,
    include_gene_sets: bool = True,
) -> dict[str, Path]:
    """Return run outputs plus immutable annotations used by selected mappings."""
    paths = resolve_magma_output_paths(
        output, dataset_id, module, module.mapping.primary,
    )
    primary_definition = module.mapping.definitions[module.mapping.primary]
    expected = {
        "magma_genes_prefix": paths["gene_prefix"],
        "magma_genes_raw": paths["genes_raw"],
        "magma_genes_out": paths["genes_out"],
        "magma_gene_annotation": paths["scoped_annotation"],
        "magma_excluded_variants": paths["excluded_variants"],
        "magma_excluded_genes": paths["excluded_genes"],
        "magma_exclusion_summary": paths["exclusion_summary"],
        "magma_harmonised_p_values": paths["harmonised_p_values"],
        "magma_harmonised_snp_locations": paths["harmonised_snp_locations"],
        "magma_mapping_comparison": paths["mapping_comparison"],
        "magma_pipeline_summary_csv": paths["pipeline_summary_csv"],
        "magma_pipeline_summary_html": paths["pipeline_summary_html"],
    }
    if include_gene_sets:
        expected["magma_pathway_compatible_gene_locations"] = paths[
            "pathway_compatible_gene_locations"
        ]
    if primary_definition.method == "chrom_magma":
        expected["chrom_magma_gene_report"] = paths["chrom_magma_genes"]
    else:
        expected["magma_genes_annotated"] = paths["corrected_genes"]
    primary_gene_set = (
        primary_definition.gene_set_file or module.input.gene_set_file
    )
    if (
        include_gene_sets
        and primary_definition.method != "chrom_magma"
        and primary_gene_set is not None
    ):
        expected.update(
            {
                "magma_gene_sets_raw": paths["gene_sets_raw"],
                "magma_gene_sets_corrected": paths["corrected_gene_sets"],
                "magma_pathway": paths["annotated_gene_sets"],
            }
        )
    for name in module.mapping.selected:
        if name == module.mapping.primary:
            continue
        mapping_paths = resolve_magma_output_paths(output, dataset_id, module, name)
        definition = module.mapping.definitions[name]
        expected.update({
            "%s_magma_genes_raw" % name: mapping_paths["genes_raw"],
            "%s_magma_genes_out" % name: mapping_paths["genes_out"],
            "%s_magma_gene_annotation" % name: mapping_paths["scoped_annotation"],
            "%s_magma_excluded_genes" % name: mapping_paths["excluded_genes"],
            "%s_magma_exclusion_summary" % name: mapping_paths[
                "exclusion_summary"
            ],
        })
        if include_gene_sets:
            expected[
                "%s_magma_pathway_compatible_gene_locations" % name
            ] = mapping_paths["pathway_compatible_gene_locations"]
        if definition.method == "chrom_magma":
            expected["%s_chrom_magma_genes" % name] = mapping_paths[
                "chrom_magma_genes"
            ]
        else:
            expected["%s_magma_genes_annotated" % name] = mapping_paths[
                "corrected_genes"
            ]
        if include_gene_sets and definition.method != "chrom_magma" and (
            definition.gene_set_file or module.input.gene_set_file
        ):
            expected.update(
                {
                    "%s_magma_gene_sets_raw" % name: mapping_paths["gene_sets_raw"],
                    "%s_magma_gene_sets_corrected" % name: mapping_paths[
                        "corrected_gene_sets"
                    ],
                    "%s_magma_pathway" % name: mapping_paths["annotated_gene_sets"],
                }
            )
    return expected


def _relative_to_directory(path: Path, directory: Path) -> Path | None:
    try:
        return path.expanduser().resolve().relative_to(directory.expanduser().resolve())
    except (OSError, ValueError):
        return None


def _existing_primary_outputs(
    output: Path,
    dataset_id: str,
    module,
    *,
    include_gene_sets: bool = True,
) -> list[Path]:
    """Return only existing files owned by this run's output directory."""
    return [
        path
        for name, path in _expected_magma_artifacts(
            output,
            dataset_id,
            module,
            include_gene_sets=include_gene_sets,
        ).items()
        if (
            name != "magma_genes_prefix"
            and _relative_to_directory(path, output) is not None
            and path.exists()
        )
    ]


def _publish_staged_files(staging: Path, output: Path, overwrite: bool) -> list[Path]:
    sources = sorted(path for path in staging.rglob("*") if path.is_file())
    destinations = [output / source.relative_to(staging) for source in sources]
    conflicts = [path for path in destinations if path.exists()]
    if conflicts and not overwrite:
        raise MagmaError(
            "MAGMA output already exists: %s. Use --resume or --overwrite."
            % ", ".join(str(path) for path in conflicts[:5])
        )
    for source, destination in zip(sources, destinations):
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.replace(destination)
    shutil.rmtree(staging)
    return destinations


def _published_value(value: Any, staging: Path, output: Path) -> Any:
    if isinstance(value, dict):
        return {name: _published_value(item, staging, output) for name, item in value.items()}
    if isinstance(value, list):
        return [_published_value(item, staging, output) for item in value]
    if isinstance(value, str):
        path = Path(value)
        relative = _relative_to_directory(path, staging)
        if relative is None:
            return value
        return str((output / relative).resolve())
    return value


def _staging_has_material_entries(staging: Path) -> bool:
    """Return whether staging contains anything beyond empty directories."""
    if staging.is_symlink() or not staging.is_dir():
        return True
    return any(
        path.is_symlink() or not path.is_dir()
        for path in staging.rglob("*")
    )


def _pipeline_includes_magma_gene_sets(args, module) -> bool:
    """Mirror the validated single-target dependency contracts used at runtime."""
    requested = tuple(
        getattr(
            args,
            "_pipeline_requested_modules",
            getattr(args, "modules", ()),
        )
        or ()
    )
    gene_only_progress_available = (
        len(module.mapping.selected) == 1
        and module.mapping.definitions[module.mapping.selected[0]].method
        == "positional"
    )
    return not (
        requested in {("magmacovar",), ("pops",)}
        and gene_only_progress_available
    )


def preflight_magma_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate external MAGMA resources before VCF-derived tables are made."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = _resolved_configuration(args)
    module = configuration.modules.magma
    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    if module.genome_build.value != observed_build:
        raise MagmaError(
            "The harmonised GWAS-VCF declares genome build %s, but MAGMA "
            "resolves to %s. The GWAS, PLINK LD reference, and mapping "
            "resources must use the same build."
            % (observed_build, module.genome_build.value)
        )
    include_gene_sets = _pipeline_includes_magma_gene_sets(args, module)
    buffered_log = PreflightLogBuffer()
    reference = preflight_magma_references(
        configuration.run.dataset_id,
        configuration,
        buffered_log,
        include_gene_sets=include_gene_sets,
    )
    reference_files = validate_plink_files(
        reference.ld_reference_prefix,
        module.input.required_reference_extensions,
        error_type=MagmaError,
    )
    configure_reference_variant_identifiers(
        args,
        configuration.modules.formatting,
        [BimIdentifierRequirement(
            consumer="MAGMA",
            formatter_target="magma",
            bim_file=reference_files["bim"],
            column_roles=module.input.bim_columns,
            delimiter_pattern=module.input.table_delimiter_pattern,
        )],
    )
    validate_plink_bundle_dimensions(
        reference_files, variants=args.variant_id_observations["magma"]["variants"],
        error_type=MagmaError,
    )
    primary_name = module.mapping.primary
    primary_definition = module.mapping.definitions[primary_name]
    primary_gene_location_metrics = None
    if primary_definition.method == "positional":
        location_file = (
            primary_definition.gene_location_file
            or module.input.gene_location_file
        )
        locations = read_gene_locations(
            location_file,
            module,
            "%s gene-location file" % primary_name,
        )
        alternate_ids = [
            values[4] for values in locations.values() if values[4] is not None
        ]
        alternate_counts = Counter(alternate_ids)
        primary_gene_location_metrics = {
            "path": str(Path(location_file).expanduser().resolve()),
            "genes": len(locations),
            "alternate_unique_ids": len(alternate_counts),
            "ambiguous_alternate_ids": sum(
                count > 1 for count in alternate_counts.values()
            ),
        }
    selected_definitions = [
        module.mapping.definitions[name] for name in module.mapping.selected
    ]
    resource_files = [*discover_input_files(
        module.input,
        selected_definitions,
        reference.executable,
    ).values(), *reference_files.values()]
    resources = MagmaPipelineResources(
        configuration=configuration,
        reference=reference,
        include_gene_sets=include_gene_sets,
        primary_gene_location_metrics=primary_gene_location_metrics,
        file_identities=capture_preflight_file_identities(
            resource_files,
            error_type=MagmaError,
            label="MAGMA resource",
        ),
        log_events=buffered_log.events,
    )
    return pipeline_preflight_evidence(
        "magma",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate formatter-created SNP-location and association tables.",
            "Validate exact GWAS-variant overlap with the PLINK BIM.",
        ),
    )


def _magma_pipeline_execution_configuration(args, resources):
    """Add pipeline-generated inputs to the validated MAGMA configuration."""
    configuration = resources.configuration
    run = configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    magma = configuration.modules.magma
    inputs = magma.input.model_copy(update={
        "snp_location_file": str(args.snp_location_file),
        "p_value_file": str(args.p_value_file),
    })
    modules = configuration.modules.model_copy(update={
        "magma": magma.model_copy(update={"input": inputs}),
    })
    return configuration.model_copy(
        update={"run": run, "modules": modules}, deep=True,
    )


def run_magma_direct(
    args,
    ctx=None,
    *,
    pipeline_resources: MagmaPipelineResources | None = None,
):
    """Resolve configuration once, run MAGMA, and always finalize its log."""
    try:
        configuration = (
            _magma_pipeline_execution_configuration(args, pipeline_resources)
            if pipeline_resources is not None
            else _resolved_configuration(args)
        )
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
        try:
            dataset = validate_filename_component(raw_dataset, "dataset_id")
        except ValueError:
            dataset = fallback.run.dataset_id
        log_path = configured_output_path(
            output,
            fallback.modules.magma.output_layout.log_file,
            error_type=MagmaError,
            dataset_id=dataset,
        )
        write_log_record(
            log_path,
            "ERROR",
            "MAGMA configuration failed: %s: %s" % (type(exc).__name__, exc),
            sample_id=dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise

    module = configuration.modules.magma
    pipeline_stage_numbers = magma_pipeline_stage_numbers(args)
    include_gene_sets = (
        pipeline_stage_numbers is None
        or "pathway_analysis" in pipeline_stage_numbers
    )
    if (
        pipeline_resources is not None
        and pipeline_resources.include_gene_sets != include_gene_sets
    ):
        raise MagmaError(
            "MAGMA pathway-analysis scope changed after pipeline preflight"
        )
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=MagmaError,
            label="MAGMA resource",
        )
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = str(configuration.run.dataset_id).strip()
    log_path = configured_output_path(
        output,
        module.output_layout.log_file,
        error_type=MagmaError,
        dataset_id=dataset,
    )
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
        stage_progress=StageProgress(
            "MAGMA analysis progress",
            enabled=(
                configuration.logging.show_progress
                and getattr(args, "_pipeline_stage_progress", None) is None
            ),
            outcome_label_width=configuration.logging.terminal_label_width,
        ),
    )
    checkpoint = None
    analysis_started = False
    staging = None
    try:
        output.mkdir(parents=True, exist_ok=True)
        resolved_config_path = configured_output_path(
            output,
            module.output_layout.resolved_config_file,
            error_type=MagmaError,
            dataset_id=dataset,
        )
        logger.record(
            "PARAM",
            "magma_run",
            dataset_id=dataset,
            genome_build=module.genome_build.value,
            population=module.population.value,
            gene_model=module.gene_model,
            window_upstream_kb=module.gene_window_upstream_kb,
            window_downstream_kb=module.gene_window_downstream_kb,
            sample_size_column=module.input.sample_size_column,
            chromosome_prefix_pattern=module.input.chromosome_prefix_pattern,
            chromosome_aliases=module.input.chromosome_aliases,
            invalid_chromosome_labels=module.input.invalid_chromosome_labels,
            invalid_allele_labels=module.input.invalid_allele_labels,
            minimum_snp_overlap=module.snp_harmonisation.minimum_overlap_fraction,
            duplicate_policy=module.snp_harmonisation.duplicate_policy,
            resolve_variants_to_reference=(
                module.snp_harmonisation.resolve_variants_to_reference
            ),
            mhc_policy=module.mhc.policy,
            mhc_region=(
                module.mhc.region_override.model_dump()
                if module.mhc.region_override is not None
                else configuration.resources.genomes[
                    module.genome_build.value
                ].regions["mhc"].model_dump()
            ),
            excluded_chromosomes=module.chromosomes.exclude,
            analysis_contract=(
                "gene_and_pathway_association"
                if include_gene_sets else
                "gene_association_only_for_magmacovar"
            ),
            mapping_analyses=module.mapping.selected,
            primary_mapping=module.mapping.primary,
            threads=configuration.execution.threads,
            memory_gb=configuration.execution.memory_gb,
            random_seed=configuration.execution.random_seed,
            overwrite=configuration.run.overwrite,
            resume=configuration.run.resume,
        )

        completion_manifest = configured_output_path(
            output,
            module.output_layout.completion_manifest,
            error_type=MagmaError,
            dataset_id=dataset,
        )
        staging = configured_output_path(
            output,
            module.output_layout.staging_directory,
            error_type=MagmaError,
            dataset_id=dataset,
        )
        if (
            (staging.exists() or staging.is_symlink())
            and _staging_has_material_entries(staging)
            and not completion_manifest.is_file()
            and not configuration.run.overwrite
        ):
            raise MagmaError(
                "An isolated incomplete MAGMA run exists at %s without a valid "
                "completion checkpoint. Review it and use --overwrite to replace "
                "it." % staging
            )

        preflight = preflight_magma_analysis(
            dataset,
            configuration,
            logger,
            references=(
                pipeline_resources.reference
                if pipeline_resources is not None else None
            ),
            include_gene_sets=include_gene_sets,
        )
        if pipeline_resources is not None:
            replay_preflight_log(pipeline_resources.log_events, logger)
        formatter_result = (
            ctx.get("formatter", {}).get("magma", {})
            if isinstance(ctx, Mapping) else {}
        )
        variant_observation = (
            getattr(args, "variant_id_observations", None) or {}
        ).get("magma")
        if include_gene_sets:
            for name, plan in preflight.gene_set_plans.items():
                logger.record(
                    "VALIDATE",
                    "magma_gene_identifier_preflight",
                    mapping=name,
                    status=plan["status"],
                    pathway_file=plan.get("pathway_file"),
                    gene_reference_file=plan.get("gene_reference_file"),
                    input_metadata=plan.get("input_metadata"),
                    minimum_overlap=plan.get("minimum_overlap"),
                    validation=plan.get("validation"),
                    reason=plan.get("reason"),
                )
        logger.record(
            "RESULT",
            "magma_preflight_summary",
            formatter_variants=formatter_result.get("rows_out"),
            ld_reference_prefix=preflight.ld_reference_prefix,
            analysis_contract=(
                "gene_and_pathway_association"
                if include_gene_sets else "gene_association_only"
            ),
            **(
                {
                    "mapping_gene_set_status": {
                        name: plan["status"]
                        for name, plan in preflight.gene_set_plans.items()
                    }
                }
                if include_gene_sets else {}
            ),
        )
        selected_definitions = [
            module.mapping.definitions[name]
            for name in module.mapping.selected
        ]
        checkpoint = ExecutionCheckpoint(
            manifest_path=completion_manifest,
            output_root=output,
            artifact_root=output,
            identity={
                "scope": "direct_module",
                "stage": "magma",
                "dataset_id": dataset,
            },
            configuration={
                "module": module.model_dump(mode="json"),
                "include_gene_sets": include_gene_sets,
                "execution": configuration.execution.model_dump(mode="json"),
                "genome_resources": configuration.resources.genomes[
                    module.genome_build.value
                ].model_dump(mode="json"),
                "magma_executable": preflight.executable,
            },
            inputs=discover_input_files(
                module.input,
                selected_definitions,
                preflight.executable,
                excluded_roots=(output,),
            ),
            policy=configuration.run.resume_policy,
            resume=configuration.run.resume,
            overwrite=configuration.run.overwrite,
            logger=logger,
            software=software_identity(preflight.executable),
            error_type=MagmaError,
        )
        checkpoint_decision = checkpoint.prepare()
        if checkpoint_decision.action == "resume":
            resumed = decode_checkpoint_value(
                checkpoint_decision.document.get("state")
            )
            if not isinstance(resumed, dict):
                raise MagmaError(
                    "MAGMA completion checkpoint contains invalid result state"
                )
            resumed["resumed"] = True
            resumed["resume_mode"] = "validated_checkpoint"
            args.magma = resumed["magma_executable"]
            if ctx is not None:
                ctx["magma"] = resumed
            logger.record(
                "SKIP",
                "magma_analysis",
                reason="validated_checkpoint",
                outputs=resumed,
            )
            logger.record(
                "STATUS",
                "magma_run",
                status="COMPLETED",
                resumed=True,
                resume_mode="validated_checkpoint",
            )
            pipeline_progress = getattr(args, "_pipeline_stage_progress", None)
            if pipeline_progress is not None:
                if pipeline_progress.current == 0:
                    pipeline_progress.start(
                        pipeline_stage_numbers["variant_inputs"]
                    )
                variant_preparation = resumed["variant_preparation"]
                pipeline_progress.complete(
                    pipeline_stage_numbers["variant_inputs"],
                    outcome_fields=magma_variant_input_outcome_fields(
                        configuration,
                        formatter_result or {
                            "rows_in": variant_preparation["qc"]["input_rows"]
                        },
                        variant_observation or {},
                        variant_preparation["qc"],
                        {
                            "snp_loc_file": variant_preparation["snp_loc_file"],
                            "pval_file": variant_preparation["pval_file"],
                        },
                    ),
                )
                for key in (
                    "gene_annotation",
                    "gene_analysis",
                    "gene_results",
                    "pathway_analysis",
                    "pathway_results",
                ):
                    if key not in pipeline_stage_numbers:
                        continue
                    stage = pipeline_stage_numbers[key]
                    pipeline_progress.start(stage)
                    pipeline_progress.complete(
                        stage,
                        outcome="Reused checksum-validated MAGMA outputs",
                    )
                pipeline_progress.start(pipeline_stage_numbers["publish"])
                args._pipeline_stage_completion = {
                    "outcome_fields": [
                        ("success", "Published outputs", "validated and reused"),
                        ("success", "Completion checkpoint", "validated"),
                    ],
                }
            print_screen_block(render_magma_results_screen(
                resumed,
                configuration,
                dataset_id=dataset,
                output_directory=output,
                log_file=log_path,
                label_width=configuration.logging.terminal_label_width,
                resumed=True,
                include_pathway_results=include_gene_sets,
            ))
            return resumed

        existing = _existing_primary_outputs(
            output,
            dataset,
            module,
            include_gene_sets=include_gene_sets,
        )
        if existing and not configuration.run.overwrite:
            raise MagmaError(
                "Existing or incomplete MAGMA results were found without a "
                "valid completion checkpoint: %s. Use --overwrite to replace "
                "them or choose another output directory."
                % ", ".join(str(path) for path in existing)
            )

        write_resolved_configuration(
            configuration,
            resolved_config_path,
            modules=("magma",),
            resource_paths=(
                "executables.magma",
                "genomes.%s.regions.mhc" % module.genome_build.value,
            ),
            module_paths={"magma": _resolved_magma_metadata_paths(module)},
        )

        if staging.exists() or staging.is_symlink():
            if (
                not configuration.run.overwrite
                and _staging_has_material_entries(staging)
            ):
                raise MagmaError(
                    "An isolated incomplete MAGMA run exists at %s without a "
                    "valid completion checkpoint. Review it and use --overwrite "
                    "to replace it." % staging
                )
            remove_owned_directory(
                staging,
                output,
                "MAGMA staging directory",
                error_type=MagmaError,
            )
            if not configuration.run.overwrite:
                logger.record(
                    "ACTION",
                    "magma_empty_staging_removed",
                    path=str(staging),
                    reason="no_partial_files",
                )
        staging.mkdir(parents=True)
        analysis_started = True

        staged_result = run_magma_analysis(
            staging,
            dataset,
            configuration,
            logger,
            preflight=preflight,
            pipeline_progress=getattr(args, "_pipeline_stage_progress", None),
            pipeline_stage_numbers=pipeline_stage_numbers,
            formatter_result=formatter_result,
            variant_id_observation=variant_observation,
        )
        staged_paths = resolve_magma_output_paths(
            staging, dataset, module, module.mapping.primary,
        )
        summary_result = _published_value(staged_result, staging, output)
        summary_records = build_magma_pipeline_summary(
            dataset,
            configuration,
            preflight,
            summary_result,
            formatter_result=formatter_result,
            variant_id_observation=variant_observation,
            input_vcf=getattr(args, "vcf", None),
            include_pathway_stages=include_gene_sets,
        )
        write_magma_pipeline_csv(
            summary_records,
            staged_paths["pipeline_summary_csv"],
            module.pipeline_summary_schema,
        )
        write_magma_pipeline_html(
            summary_records,
            staged_paths["pipeline_summary_html"],
            dataset_id=dataset,
            schema=module.pipeline_summary_schema,
            result=staged_result,
            configuration=configuration,
            include_pathway_results=include_gene_sets,
        )
        staged_result["magma_pipeline_summary_csv"] = str(
            staged_paths["pipeline_summary_csv"]
        )
        staged_result["magma_pipeline_summary_html"] = str(
            staged_paths["pipeline_summary_html"]
        )
        logger.record(
            "OUTPUT",
            "magma_pipeline_execution_summary",
            stages=len(summary_records),
            csv=str(staged_paths["pipeline_summary_csv"]),
            html=str(staged_paths["pipeline_summary_html"]),
        )
        for summary_record in summary_records:
            logger.record(
                "RESULT",
                "magma_pipeline_summary_stage",
                **summary_record,
            )
        pipeline_progress = getattr(args, "_pipeline_stage_progress", None)
        if pipeline_progress is not None:
            pipeline_progress.start(pipeline_stage_numbers["publish"])
        published = _publish_staged_files(
            staging, output, configuration.run.overwrite,
        )
        result = _published_value(staged_result, staging, output)
        result["published_files"] = [str(path) for path in published]
        logger.record("OUTPUT", "magma_outputs", files=result["published_files"])
        logger.record("STATUS", "magma_run", status="COMPLETED")
        checkpoint.write(
            status="COMPLETED",
            declared_values={
                "result": result,
                "resolved_configuration": resolved_config_path,
            },
            state=result,
            metrics={
                "mapping_count": len(result["mapping_analyses"]),
                "retained_variants": result["variant_preparation"]["qc"][
                    "retained_rows"
                ],
            },
        )
        if pipeline_progress is not None:
            # The pipeline executor completes this final stage only after its
            # own checkpoint is safely written. This prevents a failed outer
            # checkpoint from being displayed as 100% complete.
            args._pipeline_stage_completion = {
                "outcome_fields": [
                    ("success", "Published files", len(published)),
                    ("success", "Completion checkpoint", "validated and written"),
                ],
            }
        # Downstream pipeline steps reuse the same validated executable.
        args.magma = result["magma_executable"]
        if ctx is not None:
            ctx["magma"] = result

        print_screen_block(render_magma_results_screen(
            result,
            configuration,
            dataset_id=dataset,
            output_directory=output,
            log_file=log_path,
            label_width=configuration.logging.terminal_label_width,
            include_pathway_results=include_gene_sets,
        ))
        return result
    except BaseException as exc:
        if checkpoint is not None and analysis_started:
            try:
                checkpoint.write(
                    status="PARTIAL",
                    declared_values=staging,
                    state={"error_type": type(exc).__name__, "message": str(exc)},
                )
            except BaseException as checkpoint_exc:
                logger.error(
                    "MAGMA partial checkpoint failed: %s: %s"
                    % (type(checkpoint_exc).__name__, checkpoint_exc)
                )
        if not logger.summary()["failed"]:
            logger.error("MAGMA analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = ["run_magma_direct"]
