"""Configuration, logging and publication boundary for MAGMA."""

from __future__ import annotations

import shutil
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
from postgwas.core.paths import (
    configured_output_path,
    remove_owned_directory,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.magma.errors import MagmaError
from postgwas.modules.magma.analysis import (
    preflight_magma_analysis,
    resolve_magma_output_paths,
    run_magma_analysis,
)


def _resolved_configuration(args):
    module_overrides = explicit_overrides(
        args,
        {
            "snp_location_file": "input.snp_location_file",
            "p_value_file": "input.p_value_file",
            "magma_ld_reference": "input.ld_reference_prefix",
            "gene_location_file": "input.gene_location_file",
            "gene_set_file": "input.gene_set_file",
            "sample_size_column": "input.sample_size_column",
            "resolve_variants_to_reference": (
                "snp_harmonisation.resolve_variants_to_reference"
            ),
            "minimum_snp_overlap": "snp_harmonisation.minimum_overlap_fraction",
            "minimum_gene_id_overlap": "gene_sets.minimum_gene_id_overlap_fraction",
            "window_upstream": "gene_window_upstream_kb",
            "window_downstream": "gene_window_downstream_kb",
            "gene_model": "gene_model",
            "magma_mapping": "mapping.selected",
            "primary_magma_mapping": "mapping.primary",
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


def _expected_magma_artifacts(output: Path, dataset_id: str, module) -> dict[str, Path]:
    """Return run outputs plus immutable annotations used by selected mappings."""
    paths = resolve_magma_output_paths(
        output, dataset_id, module, module.mapping.primary,
    )
    primary_definition = module.mapping.definitions[module.mapping.primary]
    primary_annotation = (
        Path(primary_definition.gene_annotation_file).expanduser().resolve()
        if primary_definition.gene_annotation_file is not None
        else paths["gene_annotation"]
    )
    expected = {
        "magma_genes_prefix": paths["gene_prefix"],
        "magma_genes_raw": paths["genes_raw"],
        "magma_genes_out": paths["genes_out"],
        "magma_gene_annotation": primary_annotation,
        "magma_harmonised_p_values": paths["harmonised_p_values"],
        "magma_harmonised_snp_locations": paths["harmonised_snp_locations"],
        "magma_mapping_comparison": paths["mapping_comparison"],
    }
    if primary_definition.method == "chrom_magma":
        expected["chrom_magma_gene_report"] = paths["chrom_magma_genes"]
    else:
        expected["magma_genes_corrected"] = paths["corrected_genes"]
    primary_gene_set = (
        primary_definition.gene_set_file or module.input.gene_set_file
    )
    if primary_definition.method != "chrom_magma" and primary_gene_set is not None:
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
        annotation = (
            Path(definition.gene_annotation_file).expanduser().resolve()
            if definition.gene_annotation_file is not None
            else mapping_paths["gene_annotation"]
        )
        expected.update({
            "%s_magma_genes_raw" % name: mapping_paths["genes_raw"],
            "%s_magma_genes_out" % name: mapping_paths["genes_out"],
            "%s_magma_gene_annotation" % name: annotation,
        })
        if definition.method == "chrom_magma":
            expected["%s_chrom_magma_genes" % name] = mapping_paths[
                "chrom_magma_genes"
            ]
        else:
            expected["%s_magma_genes_corrected" % name] = mapping_paths[
                "corrected_genes"
            ]
        if definition.method != "chrom_magma" and (
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


def _existing_primary_outputs(output: Path, dataset_id: str, module) -> list[Path]:
    """Return only existing files owned by this run's output directory."""
    return [
        path
        for name, path in _expected_magma_artifacts(output, dataset_id, module).items()
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


def run_magma_direct(args, ctx=None):
    """Resolve configuration once, run MAGMA, and always finalize its log."""
    try:
        configuration = _resolved_configuration(args)
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
            enabled=configuration.logging.show_progress,
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
            mapping_analyses=module.mapping.selected,
            primary_mapping=module.mapping.primary,
            threads=configuration.execution.threads,
            memory_gb=configuration.execution.memory_gb,
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
            dataset, configuration, logger,
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
                "execution": configuration.execution.model_dump(mode="json"),
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
            print(
                "\n%s\n"
                % screen_line(
                    "success",
                    "MAGMA: validated completed outputs reused",
                    indent=2,
                )
            )
            return resumed

        existing = _existing_primary_outputs(output, dataset, module)
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
            resource_paths=("executables.magma",),
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
        )
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
        # Downstream pipeline steps reuse the same validated executable.
        args.magma = result["magma_executable"]
        if ctx is not None:
            ctx["magma"] = result

        lines = [
            "",
            screen_line("analysis", "MAGMA analysis completed", indent=2),
            screen_field(
                "info", "Dataset", dataset, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            screen_field(
                "success",
                (
                    "Corrected gene results"
                    if "magma_genes_corrected" in result
                    else "chromMAGMA gene ranking"
                ),
                result["magma_gene_results"],
                indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            screen_field(
                "success",
                "Mapping result catalogue",
                result["magma_mapping_comparison"],
                indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
        ]
        if "magma_pathway" in result:
            lines.append(
                screen_field(
                    "success",
                    "Gene-set results",
                    result["magma_pathway"],
                    indent=6,
                    label_width=configuration.logging.terminal_label_width,
                )
            )
        lines.extend(
            [
                screen_field(
                    "info", "Full log", log_path, indent=6,
                    label_width=configuration.logging.terminal_label_width,
                ),
                "",
            ]
        )
        print("\n".join(lines))
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
