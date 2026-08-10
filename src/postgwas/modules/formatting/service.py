"""Application service for one-pass GWAS-VCF format export."""

from __future__ import annotations

from pathlib import Path
import tempfile

from rich import box
from rich.console import Console
from rich.table import Table

from postgwas.config import (
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.ui.screen import screen_field, screen_line

from .exporters.finemap import export_finemap
from .exporters.gcta_gene import export_gcta_gene
from .exporters.ldsc import export_ldsc
from .exporters.magma import export_magma
from .exporters.mixer import export_mixer
from .exporters.pred_ld import export_pred_ld
from .exporters.susie import export_susie
from .contracts import FORMAT_CONTRACTS
from .table import (
    FormattingError,
    infer_study_design,
    load_harmonised_vcf,
    select_variant_identifiers,
)
from .resume import (
    formatter_resolved_paths,
    resume_formatter_outputs,
    write_formatter_completion_manifest,
)


EXPORTERS = {
    "magma": export_magma,
    "gcta_gene": export_gcta_gene,
    "susie": export_susie,
    "finemap": export_finemap,
    "pred_ld": export_pred_ld,
    "ldsc": export_ldsc,
    "mixer": export_mixer,
}


def _resolved_configuration(args):
    module_overrides = explicit_overrides(
        args,
        {
            "format": "formats",
            "variant_id_type": "variant_identifiers.default_type",
            "variant_id_types": "variant_identifiers.target_types",
        },
    )
    if "variant_identifiers.default_type" in module_overrides:
        module_overrides["variant_identifiers.target_types"] = {}
    global_overrides = explicit_overrides(
        args,
        {
            "dataset_id": "run.dataset_id",
            "output_directory": "run.output_directory",
            "threads": "execution.threads",
            "memory_gb": "execution.memory_gb",
            "seed": "execution.random_seed",
            "bcftools": "resources.executables.bcftools",
            "resume": "run.resume",
            "overwrite": "run.overwrite",
        },
    )

    return load_run_configuration_for_module(
        "formatting",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _screen_summary(
    dataset_id,
    vcf,
    study_design,
    results,
    log_path,
    case_count_column,
    control_count_column,
    label_width,
    variant_id_observations,
):
    lines = [
        "",
        screen_line("analysis", "Formatter completed", indent=2),
        screen_field(
            "info", "Dataset", dataset_id, indent=6, label_width=label_width,
        ),
        screen_field(
            "genetic", "Input GWAS-VCF", vcf, indent=6, label_width=label_width,
        ),
        screen_field(
            "analysis",
            "Inferred trait type",
            "%s (%s present for %s/%s variants; %s present for %s/%s)"
            % (
                study_design.trait_type,
                case_count_column,
                f"{study_design.case_counts_present:,}",
                f"{study_design.rows:,}",
                control_count_column,
                f"{study_design.control_counts_present:,}",
                f"{study_design.rows:,}",
            ),
            indent=6,
            label_width=label_width,
        ),
    ]
    for target, observation in variant_id_observations.items():
        identifier_label = (
            "rsIDs"
            if observation["variant_id_type"] == "rsid"
            else "configured coordinate-and-allele unique IDs"
        )
        lines.append(screen_field(
            "genetic",
            "%s reference IDs" % target,
            "%s (%s BIM variants scanned)"
            % (identifier_label, f"{observation['variants']:,}"),
            indent=6,
            label_width=label_width,
        ))
    for target, result in results.items():
        lines.append(screen_field(
            "success",
            target,
            "%s variants exported using %s; %s excluded because the "
            "identifier or required values were missing or invalid"
            % (
                f"{result['rows_out']:,}",
                "rsIDs" if result["variant_id_type"] == "rsid" else "unique IDs",
                f"{result['rows_excluded']:,}",
            ),
            indent=6,
            label_width=label_width,
        ))
    lines.append(screen_field(
        "info", "Full log", log_path, indent=6, label_width=label_width,
    ))
    lines.append("")
    return "\n".join(lines)


def _configured_output_columns(schema, trait_type):
    """Render output columns directly from the resolved YAML schema."""
    if schema.outputs:
        return "; ".join(
            "%s: %s" % (name.replace("_", " "), ", ".join(output.columns.values()))
            for name, output in schema.outputs.items()
        )
    columns = list(schema.columns.values())
    columns.extend(schema.trait_columns.get(trait_type, {}).values())
    columns.extend(schema.trailing_columns.values())
    return ", ".join(columns)


def _print_contract_table(selected, module, trait_type):
    """Show the reviewed tool contracts without exposing implementation details."""
    table = Table(
        title="🔬  Validated downstream input schemas",
        box=box.SIMPLE_HEAVY,
        show_lines=True,
        header_style="bold cyan",
        padding=(0, 1),
    )
    table.add_column("Tool", style="bold", no_wrap=True)
    table.add_column("Columns written", ratio=4)
    table.add_column("Allele frequency", ratio=3)
    table.add_column("Sample size", ratio=3)
    for target in selected:
        contract = FORMAT_CONTRACTS[target]
        table.add_row(
            contract.name,
            _configured_output_columns(module.exports[target], trait_type),
            contract.frequency,
            contract.sample_size,
        )
    console = Console(width=128)
    console.print()
    console.print(table)


def run_formatter_direct(args, ctx=None, *, configuration=None):
    """Create all requested tool inputs from one checked VCF extraction."""
    try:
        configuration = configuration or _resolved_configuration(args)
    except BaseException as exc:
        fallback = load_run_configuration_for_module("formatting")
        failure_output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        failure_dataset = str(
            getattr(args, "dataset_id", None) or fallback.run.dataset_id
        ).strip()
        failure_log = configured_output_path(
            failure_output,
            fallback.modules.formatting.runtime.log_file,
            error_type=FormattingError,
            dataset_id=failure_dataset,
        )
        write_log_record(
            failure_log,
            "ERROR",
            "Formatter configuration failed: %s: %s" % (type(exc).__name__, exc),
            sample_id=failure_dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise
    module = configuration.modules.formatting
    output_directory = Path(configuration.run.output_directory).expanduser().resolve()
    dataset_id = str(configuration.run.dataset_id).strip()
    output_directory.mkdir(parents=True, exist_ok=True)
    log_path = configured_output_path(
        output_directory,
        module.runtime.log_file,
        error_type=FormattingError,
        dataset_id=dataset_id,
    )
    logger = PipelineLogger(
        dataset_id,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )
    work_table = None
    variant_id_observations = (
        getattr(args, "variant_id_observations", None) or {}
    )
    try:
        vcf_value = getattr(args, "vcf", None)
        if not vcf_value:
            raise FormattingError(
                "A harmonised GWAS-VCF is required; provide --vcf PATH."
            )
        vcf = Path(vcf_value).expanduser().resolve()
        selected = list(module.formats)
        if not selected:
            raise FormattingError(
                "No output format was selected. Provide --format followed by one or "
                "more of: %s, or set modules.formatting.formats in --run-config."
                % ", ".join(EXPORTERS)
            )
        if len(selected) != len(set(selected)):
            raise FormattingError("Each requested output format may appear only once.")
        unknown = [target for target in selected if target not in EXPORTERS]
        if unknown:
            raise FormattingError("Unknown output format: %s" % ", ".join(unknown))
        identifier_policy = module.variant_identifiers.model_copy(update={
            "target_types": {
                target: identifier_type
                for target, identifier_type
                in module.variant_identifiers.target_types.items()
                if target in selected
            },
        })
        module = module.model_copy(update={
            "variant_identifiers": identifier_policy,
        })
        configuration.modules.formatting = module

        bcftools = configuration.resources.executables.bcftools
        try:
            resolved_bcftools = resolve_executable(
                str(bcftools), "bcftools executable", error_type=FormattingError,
            )
        except FormattingError as exc:
            raise FormattingError(
                "%s. Set resources.executables.bcftools in --run-config or provide "
                "--bcftools PATH." % exc
            ) from exc

        logger.record("INPUT", "formatter_run", dataset=dataset_id, vcf=str(vcf))
        logger.record("PARAM", "formats", values=selected)
        logger.record("PARAM", "minimum_p_value", value=module.minimum_p_value)
        logger.record("PARAM", "vcf_fields", values=module.vcf_fields.model_dump())
        logger.record(
            "PARAM",
            "chromosome_labels",
            values=module.chromosome_labels.model_dump(mode="json"),
        )
        logger.record(
            "PARAM", "variant_identifiers",
            values=module.variant_identifiers.model_dump(mode="json"),
        )
        for target, observation in variant_id_observations.items():
            logger.record(
                "OBSERVED", "reference_variant_identifiers",
                target=target, **observation,
            )
        for target in selected:
            contract = FORMAT_CONTRACTS[target]
            schema = module.exports[target]
            logger.record(
                "PARAM", "formatter_schema", target=target,
                columns=schema.columns,
                trait_columns=schema.trait_columns,
                trailing_columns=schema.trailing_columns,
                named_outputs={
                    name: output.columns for name, output in schema.outputs.items()
                },
                frequency=contract.frequency,
                sample_size=contract.sample_size, source=contract.source,
            )
        resolved_config_path = configured_output_path(
            output_directory,
            module.runtime.resolved_config_file,
            error_type=FormattingError,
            dataset_id=dataset_id,
        )
        resolved_module_paths = {
            "formatting": formatter_resolved_paths(configuration, selected),
        }
        if configuration.run.resume and not configuration.run.overwrite:
            resumed = resume_formatter_outputs(
                output_directory=output_directory,
                dataset_id=dataset_id,
                vcf=vcf,
                selected=selected,
                configuration=configuration,
                log_path=log_path,
            )
            if resumed is not None:
                write_resolved_configuration(
                    configuration,
                    resolved_config_path,
                    modules=("formatting",),
                    resource_paths=("executables.bcftools",),
                    module_paths=resolved_module_paths,
                )
                if ctx is not None:
                    ctx["formatter"] = resumed
                logger.record(
                    "SKIP", "formatter_run", reason="validated_resume",
                    formats=selected,
                )
                logger.record(
                    "STATUS", "formatter_run", status="COMPLETED",
                    formats=selected, resumed=True,
                )
                print(screen_line(
                    "success",
                    "Formatter outputs validated; continuing from completed step",
                    indent=2,
                ))
                return resumed
        write_resolved_configuration(
            configuration,
            resolved_config_path,
            modules=("formatting",),
            resource_paths=("executables.bcftools",),
            module_paths=resolved_module_paths,
        )

        handle = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=module.runtime.temporary_table_suffix,
            prefix=module.runtime.temporary_table_prefix.format(dataset_id=dataset_id),
            dir=output_directory,
            delete=False,
        )
        work_table = Path(handle.name)
        handle.close()
        with logger.step(
            1, len(selected) + 1, "Read harmonised GWAS-VCF", "load_harmonised_vcf",
        ) as step:
            frame = load_harmonised_vcf(
                vcf,
                work_table,
                dataset_id,
                resolved_bcftools,
                module,
                logger=logger,
            )
            step.set_rows(frame.height, removed=0)
            step.output("canonical_variants", rows=frame.height)

        study_design = infer_study_design(
            frame,
            module.study_design.case_count_column,
            module.study_design.control_count_column,
        )
        logger.record(
            "DECIDE",
            "study_design",
            trait_type=study_design.trait_type,
            rule="quantitative_if_case_count_entirely_missing_otherwise_binary",
            case_count_column=module.study_design.case_count_column,
            control_count_column=module.study_design.control_count_column,
            case_count_present=study_design.case_counts_present,
            control_count_present=study_design.control_counts_present,
            variants=study_design.rows,
            metadata_header_used=False,
        )
        if (
            study_design.trait_type == "binary"
            and (
                study_design.case_counts_present < study_design.rows
                or study_design.control_counts_present < study_design.rows
            )
        ):
            logger.warning(
                "The study is binary because the configured case-count column contains "
                "at least one value, but case/control counts are incomplete for some "
                "variants. Count-dependent outputs will exclude those variants."
            )

        _print_contract_table(selected, module, study_design.trait_type)

        results = {}
        for number, target in enumerate(selected, 2):
            identifier_type = module.variant_identifiers.target_types.get(
                target, module.variant_identifiers.default_type,
            )
            target_frame, identifier_qc = select_variant_identifiers(
                frame, module, identifier_type,
            )
            logger.record(
                "DECIDE", "formatter_variant_identifiers",
                target=target, **identifier_qc,
            )
            with logger.step(
                number,
                len(selected) + 1,
                "Create %s input" % target,
                EXPORTERS[target].__name__,
                rows_in=frame.height,
            ) as step:
                result = EXPORTERS[target](
                    target_frame,
                    output_directory,
                    dataset_id,
                    module,
                    overwrite=configuration.run.overwrite,
                )
                excluded = int(result["rows_excluded"])
                identifier_excluded = int(identifier_qc["identifier_rows_excluded"])
                if identifier_excluded:
                    step.qc(
                        "%s variant identifiers" % identifier_type,
                        "Exclude records that cannot represent the selected "
                        "identifier convention.",
                        frame.height,
                        target_frame.height,
                        reason="missing_or_invalid_variant_identifier",
                        warn=True,
                    )
                step.qc(
                    "required formatter fields",
                    "Exclude records that cannot be represented correctly in this tool's input.",
                    target_frame.height,
                    int(result["rows_out"]),
                    reason="missing_or_invalid_required_value",
                    warn=excluded > 0,
                )
                bounded = int(result.get("p_values_bounded", 0))
                if bounded:
                    step.qc(
                        "raw p-value numeric bound",
                        "Bound raw p-values at the configured minimum to avoid floating-point underflow.",
                        int(result["rows_out"]),
                        int(result["rows_out"]),
                        changed=bounded,
                        warn=True,
                    )
                total_excluded = frame.height - int(result["rows_out"])
                step.set_rows(int(result["rows_out"]), removed=total_excluded)
                for name, value in result.items():
                    if name.endswith(("_file", "_input", "_folder")):
                        step.output(name, path=value)
                step.output(
                    "written_schema", target=target, columns=result.get("columns"),
                    sample_size_mode=result.get("sample_size_mode"),
                )
                result["log_file"] = str(log_path)
                result.update(identifier_qc)
                result["rows_in"] = frame.height
                result["schema_rows_excluded"] = excluded
                result["rows_excluded"] = total_excluded
                results[target] = result

        if ctx is not None:
            ctx["formatter"] = results
        write_formatter_completion_manifest(
            output_directory=output_directory,
            dataset_id=dataset_id,
            vcf=vcf,
            selected=selected,
            configuration=configuration,
            results=results,
        )
        logger.record(
            "STATUS", "formatter_run", status="COMPLETED",
            formats=selected, input_variants=frame.height,
        )
        print(_screen_summary(
            dataset_id,
            str(vcf),
            study_design,
            results,
            str(log_path),
            module.study_design.case_count_column,
            module.study_design.control_count_column,
            configuration.logging.terminal_label_width,
            variant_id_observations,
        ))
        return results
    except BaseException as exc:
        if not logger.summary()["failed"]:
            logger.error("Formatter failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        if work_table is not None:
            work_table.unlink(missing_ok=True)
        logger.close()


__all__ = ["EXPORTERS", "run_formatter_direct"]
