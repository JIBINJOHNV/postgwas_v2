"""Application service for one-pass GWAS-VCF format export."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import tempfile
from typing import get_args

from postgwas.config import (
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.formatting import FormattingSampleSizeRole
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.preflight import require_unchanged_preflight_files
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.vcf import vcf_query_field_label

from .contracts import (
    CUSTOM_OUTPUT_TARGET,
    DUPLICATE_POLICY_DESCRIPTIONS,
    FORMAT_CONTRACTS,
    formatter_result_targets,
    formatter_target_display_name,
)
from .exporters.custom import export_custom
from .exporters.finemap import export_finemap
from .exporters.gcta_gene import export_gcta_gene
from .exporters.ldsc import export_ldsc
from .exporters.magma import export_magma
from .exporters.mixer import export_mixer
from .exporters.pred_ld import export_pred_ld
from .exporters.susie import export_susie
from .ldsc_reference import (
    resolve_ldsc_merge_alleles_file,
    select_ldsc_reference_variants,
)
from .resume import (
    formatter_output_paths,
    formatter_resolved_paths,
    magma_formatter_content_paths,
    resume_formatter_outputs,
    write_formatter_completion_manifest,
)
from .reporting import (
    ldsc_sample_prevalence_screen_fields as _ldsc_sample_prevalence_screen_fields,
    render_formatter_screen_summary as _screen_summary,
    write_formatter_html_report,
)
from .table import (
    FormattingError,
    infer_study_design,
    load_harmonised_vcf,
    resolve_duplicate_identifiers,
    select_variant_identifiers,
    validate_unique_identifiers,
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


def _magma_variant_input_pipeline_progress(args, *, pipeline_mode: bool):
    """Resolve the shared formatter-to-MAGMA stage from the attached plan."""
    progress = getattr(args, "_pipeline_stage_progress", None)
    plan = getattr(args, "_pipeline_progress_plan", {}) or {}
    stage = (plan.get("magma_stage_numbers") or {}).get("variant_inputs")
    if not pipeline_mode or progress is None or stage is None:
        return None
    return progress, stage


def _identifier_selection_settings(module, target: str) -> tuple[str, str]:
    """Resolve one target's identifier type and duplicate policy."""
    return (
        module.variant_identifiers.target_types.get(
            target, module.variant_identifiers.default_type,
        ),
        module.variant_identifiers.target_duplicate_policies.get(
            target,
            module.variant_identifiers.default_duplicate_policy,
        ),
    )


def _validate_output_destinations(
    output_directory: Path,
    dataset_id: str,
    selected: list[str],
    module,
    *,
    configuration=None,
) -> dict[str, str]:
    """Resolve every selected output and reject collisions before extraction."""
    output_paths = formatter_output_paths(
        output_directory,
        dataset_id,
        selected,
        module,
        configuration=configuration,
    )
    destinations: dict[Path, list[str]] = {}
    for label, path in output_paths.items():
        destinations.setdefault(path, []).append(label)
    for label, pattern in {
        "formatter log": module.runtime.log_file,
        "resolved configuration": module.runtime.resolved_config_file,
        "completion manifest": module.runtime.completion_manifest_file,
    }.items():
        path = configured_output_path(
            output_directory,
            pattern,
            error_type=FormattingError,
            dataset_id=dataset_id,
        )
        destinations.setdefault(path, []).append(label)

    collisions = [
        (path, labels)
        for path, labels in destinations.items()
        if len(labels) > 1
    ]
    if collisions:
        details = "; ".join(
            "%s resolve to %s" % (" and ".join(labels), path)
            for path, labels in collisions
        )
        raise FormattingError(
            "Formatter output collision: %s. Configure a unique output_file or "
            "partition_file for each selected output." % details
        )
    return {label: str(path) for label, path in output_paths.items()}


def _resolved_configuration(args):
    module_overrides = explicit_overrides(
        args,
        {
            "format": "formats",
            "variant_id_type": "variant_identifiers.default_type",
            "variant_id_types": "variant_identifiers.target_types",
            "variant_id_duplicate_policies": (
                "variant_identifiers.target_duplicate_policies"
            ),
            "duplicate_id_policy": (
                "variant_identifiers.default_duplicate_policy"
            ),
            "merge_alleles": "ldsc_reference.merge_alleles_file",
        },
    )
    if "variant_identifiers.default_type" in module_overrides:
        module_overrides["variant_identifiers.target_types"] = {}
    if "variant_identifiers.default_duplicate_policy" in module_overrides:
        module_overrides["variant_identifiers.target_duplicate_policies"] = {}
    if hasattr(args, "custom_output_file"):
        module_overrides["custom_output.output_file"] = args.custom_output_file
    if hasattr(args, "custom_columns"):
        module_overrides["custom_output.columns"] = dict(args.custom_columns)
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


def _column_records(mapping, transformations):
    """Return ordered source-to-saved mappings with configured transformations."""
    return [
        {
            "source": source,
            "saved_as": saved_as,
            "transformation": transformations.get(source),
        }
        for source, saved_as in mapping.items()
    ]


def _configured_output_schemas(module, target, trait_type):
    """Resolve the exact output mappings used by one selected formatter target."""
    if target == CUSTOM_OUTPUT_TARGET:
        columns = []
        for role, saved_as in module.custom_output.columns.items():
            contract = module.custom_output.field_contracts[role]
            columns.append({
                "source": contract.source,
                "saved_as": saved_as,
                "transformation": contract.transformation,
            })
        return [{"output": "custom", "columns": columns}]

    schema = module.exports[target]
    if schema.outputs:
        return [
            {
                "output": name,
                "columns": _column_records(
                    output.columns, output.transformations,
                ),
            }
            for name, output in schema.outputs.items()
        ]
    mapping = dict(schema.columns)
    if schema.trait_columns:
        if trait_type not in schema.trait_columns:
            raise FormattingError(
                "Cannot report the %s output schema without its binary or "
                "quantitative study design." % target
            )
        mapping.update(schema.trait_columns[trait_type])
    mapping.update(schema.trailing_columns)
    return [{
        "output": "table",
        "columns": _column_records(mapping, schema.transformations),
    }]


def _render_column_mappings(outputs):
    multiple_outputs = len(outputs) > 1
    rendered = []
    for output in outputs:
        mappings = ", ".join(
            "%s → %s" % (column["source"], column["saved_as"])
            for column in output["columns"]
        )
        if multiple_outputs:
            mappings = "%s: %s" % (
                output["output"].replace("_", " "), mappings,
            )
        rendered.append(mappings)
    return "; ".join(rendered)


def _render_saved_statistic(
    outputs,
    source,
    *,
    unchanged,
    transformed_by,
    transformed,
):
    multiple_outputs = len(outputs) > 1
    rendered = []
    for output in outputs:
        for column in output["columns"]:
            if (
                column["source"] != source
                and column["transformation"] != transformed_by
            ):
                continue
            transformation = column["transformation"]
            if transformation is None:
                representation = unchanged
            elif transformation == transformed_by:
                representation = transformed
            else:
                representation = "configured transformation: %s" % transformation
            label = "%s → %s (%s)" % (
                column["source"], column["saved_as"], representation,
            )
            if multiple_outputs:
                label = "%s: %s" % (
                    output["output"].replace("_", " "), label,
                )
            rendered.append(label)
    return "; ".join(rendered) or "not saved"


def _render_saved_sample_sizes(outputs, module, target, trait_type):
    """Report resolved source, numerical meaning, destination, and tool use."""
    canonical = module.canonical_columns
    reporting = module.sample_size_reporting
    source_roles = {
        getattr(canonical, role): role
        for role in get_args(FormattingSampleSizeRole)
    }
    multiple_outputs = len(outputs) > 1
    rendered = []
    for output in outputs:
        for column in output["columns"]:
            source = column["source"]
            role = source_roles.get(source)
            if role is None:
                continue
            meaning = reporting.source_semantics[role].resolve(trait_type)
            vcf_source = vcf_query_field_label(module.vcf_fields.root[source])
            transformation = column["transformation"]
            representation = (
                "copied"
                if transformation is None
                else "configured transformation: %s" % transformation
            )
            label = "%s: %s → %s (%s)" % (
                meaning, vcf_source, column["saved_as"], representation,
            )
            if multiple_outputs:
                label = "%s · %s" % (
                    output["output"].replace("_", " "), label,
                )
            rendered.append(label)
    if not rendered:
        return "not saved"
    if target != CUSTOM_OUTPUT_TARGET:
        rendered.append(reporting.target_notes[target].resolve(trait_type))
    return "; ".join(rendered)


def _configured_schema_reports(
    targets,
    module,
    trait_types=None,
    variant_id_observations=None,
):
    """Build screen/log reports solely from resolved formatter configuration."""
    trait_types = trait_types or {}
    reports = {}
    canonical = module.canonical_columns
    unique_id_label = "coordinate-and-allele ID (%s)" % (
        module.variant_identifiers.unique_id_template.format(
            chromosome=canonical.chromosome,
            position=canonical.position,
            reference_allele=canonical.reference_allele,
            alternate_allele=canonical.alternate_allele,
        )
    )
    for target in targets:
        outputs = _configured_output_schemas(
            module, target, trait_types.get(target),
        )
        variant_id_type, duplicate_id_policy = _identifier_selection_settings(
            module, target,
        )
        name = formatter_target_display_name(
            target, variant_id_observations,
        )
        if target == CUSTOM_OUTPUT_TARGET:
            frequency_interpretation = "custom field semantics"
        else:
            contract = FORMAT_CONTRACTS[target]
            frequency_interpretation = contract.frequency
        reports[target] = {
            "name": name,
            "variant_id_type": variant_id_type,
            "duplicate_id_policy": duplicate_id_policy,
            "duplicate_id_policy_description": (
                DUPLICATE_POLICY_DESCRIPTIONS[duplicate_id_policy]
            ),
            "variant_id_type_label": (
                "rsID"
                if variant_id_type == "rsid"
                else unique_id_label
            ),
            "outputs": outputs,
            "column_mappings": _render_column_mappings(outputs),
            "p_value": _render_saved_statistic(
                outputs,
                canonical.negative_log10_p_value,
                unchanged="-log10(P)",
                transformed_by="negative_log10_to_raw_p",
                transformed="raw P (10^(-%s))"
                % canonical.negative_log10_p_value,
            ),
            "frequency": _render_saved_statistic(
                outputs,
                canonical.effect_allele_frequency,
                unchanged="effect-allele frequency (EAF)",
                transformed_by="effect_frequency_to_minor_frequency",
                transformed=(
                    "minor-allele frequency (MAF = min(%s, 1-%s))"
                    % (
                        canonical.effect_allele_frequency,
                        canonical.effect_allele_frequency,
                    )
                ),
            ),
            "frequency_interpretation": frequency_interpretation,
            "sample_size": _render_saved_sample_sizes(
                outputs, module, target, trait_types.get(target),
            ),
        }
    return reports


def _log_schema_reports(logger, reports):
    for target, report in reports.items():
        logger.record(
            "DECIDE",
            "formatter_saved_schema",
            target=target,
            downstream_consumer=report["name"],
            variant_id_type=report["variant_id_type"],
            duplicate_id_policy=report["duplicate_id_policy"],
            outputs=report["outputs"],
            p_value_saved=report["p_value"],
            allele_frequency_saved=report["frequency"],
            frequency_interpretation=report["frequency_interpretation"],
            sample_size_saved=report["sample_size"],
        )


def _print_contract_summary(reports, label_width):
    """Explain the configured file preparation without implying analysis ran."""
    lines = [
        "",
        screen_line(
            "analysis",
            "Formatter plan · prepare files required by downstream analyses",
            indent=2,
        ),
    ]
    for report in reports.values():
        lines.extend([
            "",
            screen_line("analysis", report["name"], indent=6),
            screen_field(
                "info",
                "What this step does",
                "Create validated PostGWAS intermediate files for this tool; "
                "the downstream analysis is not run in this step",
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "genetic",
                "Variant identifiers written",
                report["variant_id_type_label"],
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "info",
                "Duplicate-ID handling",
                "%s — %s"
                % (
                    report["duplicate_id_policy"],
                    report["duplicate_id_policy_description"],
                ),
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "info",
                "Columns and transformations",
                report["column_mappings"],
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "analysis",
                "P value written",
                report["p_value"],
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "genetic",
                "Allele frequency written",
                report["frequency"],
                indent=10,
                label_width=label_width,
            ),
            screen_field(
                "count",
                "Sample size written",
                report["sample_size"],
                indent=10,
                label_width=label_width,
            ),
        ])
    print("\n".join(lines))


def _log_ldsc_sample_prevalence(logger, result, module):
    """Record the configured binary-trait prevalence reduction or its skip."""
    aggregation = module.ldsc_sample_prevalence.aggregation
    trait_type = result.get("trait_type")
    if trait_type != "binary":
        logger.record(
            "SKIP",
            "ldsc_sample_prevalence",
            reason="not_applicable_to_quantitative_trait",
            trait_type=trait_type,
            configured_aggregation=aggregation,
        )
        return
    logger.record(
        "DECIDE",
        "ldsc_sample_prevalence",
        aggregation=result.get("sample_prevalence_aggregation"),
        formula="%s/(%s+%s)"
        % (
            module.study_design.case_count_column,
            module.study_design.case_count_column,
            module.study_design.control_count_column,
        ),
        variants_used=result.get("sample_prevalence_variants"),
        minimum=result.get("sample_prevalence_minimum"),
        maximum=result.get("sample_prevalence_maximum"),
        case_count_minimum=result.get(
            "sample_prevalence_case_count_minimum"
        ),
        case_count_maximum=result.get(
            "sample_prevalence_case_count_maximum"
        ),
        control_count_minimum=result.get(
            "sample_prevalence_control_count_minimum"
        ),
        control_count_maximum=result.get(
            "sample_prevalence_control_count_maximum"
        ),
        value=result.get("sample_prev"),
    )


def run_formatter_direct(
    args,
    ctx=None,
    *,
    configuration=None,
    emit_terminal_summary: bool = True,
):
    """Create all requested tool inputs from one checked VCF extraction."""
    pipeline_mode = ctx is not None
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
    magma_pipeline_evidence = (
        ctx.validation("magma")
        if pipeline_mode and hasattr(ctx, "validation") else None
    )
    magma_pipeline_resources = (
        magma_pipeline_evidence.resources
        if magma_pipeline_evidence is not None else None
    )
    magma_formatter_configuration = (
        magma_pipeline_resources.configuration
        if magma_pipeline_resources is not None else None
    )
    if pipeline_mode and hasattr(ctx, "validation"):
        for module_name in ctx.validation_modules():
            evidence = ctx.validation(module_name)
            resources = getattr(evidence, "resources", None)
            identities = getattr(resources, "file_identities", None)
            if identities is not None:
                require_unchanged_preflight_files(
                    identities,
                    error_type=FormattingError,
                    label="%s resource" % module_name.replace("_", " "),
                )
    if magma_formatter_configuration is not None:
        configuration.modules.magma = (
            magma_formatter_configuration.modules.magma.model_copy(
                update={"enabled": True},
            )
        )
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
    progress = None
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
        targets = formatter_result_targets(module, selected)
        if not targets:
            raise FormattingError(
                "No output was selected. Provide --format followed by one or more "
                "of: %s, or provide --custom-output with --id and the requested "
                "column-name options."
                % ", ".join(EXPORTERS)
            )
        if len(selected) != len(set(selected)):
            raise FormattingError("Each requested output format may appear only once.")
        unknown = [target for target in selected if target not in EXPORTERS]
        if unknown:
            raise FormattingError("Unknown output format: %s" % ", ".join(unknown))
        merge_alleles_file = resolve_ldsc_merge_alleles_file(module, selected)
        if merge_alleles_file is not None:
            ldsc_reference = module.ldsc_reference.model_copy(update={
                "merge_alleles_file": merge_alleles_file,
            })
            module = module.model_copy(update={
                "ldsc_reference": ldsc_reference,
            })
        identifier_policy = module.variant_identifiers.model_copy(update={
            "target_types": {
                target: identifier_type
                for target, identifier_type
                in module.variant_identifiers.target_types.items()
                if target in selected
            },
            "target_duplicate_policies": {
                target: duplicate_policy
                for target, duplicate_policy
                in module.variant_identifiers.target_duplicate_policies.items()
                if target in selected
            },
        })
        module = module.model_copy(update={
            "variant_identifiers": identifier_policy,
        })
        configuration.modules.formatting = module
        output_destinations = _validate_output_destinations(
            output_directory,
            dataset_id,
            selected,
            module,
            configuration=configuration,
        )
        html_report_path = Path(
            output_destinations["formatter.html_report"]
        )
        completion_manifest_path = configured_output_path(
            output_directory,
            module.runtime.completion_manifest_file,
            error_type=FormattingError,
            dataset_id=dataset_id,
        )

        current_vcf_validation = (
            ctx.validation("current_vcf", {})
            if pipeline_mode and hasattr(ctx, "validation")
            else {}
        )
        indexed_vcf_validation = (
            current_vcf_validation.get("indexed")
            if isinstance(current_vcf_validation, dict)
            else None
        )
        executable_identity = (
            current_vcf_validation.get("bcftools_identity")
            if isinstance(current_vcf_validation, dict)
            else None
        )
        if executable_identity is not None:
            require_unchanged_preflight_files(
                executable_identity,
                error_type=FormattingError,
                label="bcftools executable",
            )
        if indexed_vcf_validation is not None:
            resolved_bcftools = indexed_vcf_validation.bcftools
        else:
            bcftools = configuration.resources.executables.bcftools
            try:
                resolved_bcftools = resolve_executable(
                    str(bcftools),
                    "bcftools executable",
                    error_type=FormattingError,
                )
            except FormattingError as exc:
                raise FormattingError(
                    "%s. Install bcftools on PATH or set "
                    "resources.executables.bcftools in --run-config." % exc
                ) from exc

        logger.record("INPUT", "formatter_run", dataset=dataset_id, vcf=str(vcf))
        logger.record("PARAM", "formats", values=targets)
        logger.record(
            "VALIDATE",
            "formatter_output_destinations",
            status="PASSED",
            values=output_destinations,
        )
        logger.record("PARAM", "minimum_p_value", value=module.minimum_p_value)
        logger.record(
            "PARAM",
            "formatter_input_contract",
            values=module.input_contract.model_dump(mode="json"),
        )
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
        if merge_alleles_file is not None:
            logger.record(
                "INPUT",
                "ldsc_merge_alleles",
                path=str(merge_alleles_file),
                matching="rsid_and_strand_unambiguous_allele_pair",
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
                transformations=schema.transformations,
                named_outputs={
                    name: {
                        "columns": output.columns,
                        "transformations": output.transformations,
                    }
                    for name, output in schema.outputs.items()
                },
                frequency=contract.frequency,
                source=contract.source,
            )
        if module.custom_output.active:
            logger.record(
                "PARAM",
                "formatter_schema",
                target=CUSTOM_OUTPUT_TARGET,
                output_file=module.custom_output.output_file,
                columns=module.custom_output.columns,
                field_contracts={
                    role: module.custom_output.field_contracts[role].model_dump(
                        mode="json"
                    )
                    for role in module.custom_output.columns
                },
                required_field_policy="exclude_rows_missing_any_requested_field",
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
        resolved_modules = ["formatting"]
        if "magma" in selected and magma_formatter_configuration is not None:
            resolved_modules.append("magma")
            resolved_module_paths["magma"] = magma_formatter_content_paths()
        if configuration.run.resume and not configuration.run.overwrite:
            resumed = resume_formatter_outputs(
                output_directory=output_directory,
                dataset_id=dataset_id,
                vcf=vcf,
                selected=selected,
                configuration=configuration,
                log_path=log_path,
                logger=logger,
            )
            if resumed is not None:
                write_resolved_configuration(
                    configuration,
                    resolved_config_path,
                    modules=tuple(resolved_modules),
                    resource_paths=("executables.bcftools",),
                    module_paths=resolved_module_paths,
                )
                if ctx is not None:
                    ctx["formatter"] = resumed
                resumed_trait_types = {
                    target: result.get("trait_type")
                    for target, result in resumed.items()
                    if result.get("trait_type") is not None
                }
                schema_reports = _configured_schema_reports(
                    targets,
                    module,
                    resumed_trait_types,
                    variant_id_observations,
                )
                _log_schema_reports(logger, schema_reports)
                if emit_terminal_summary and not pipeline_mode:
                    _print_contract_summary(
                        schema_reports,
                        configuration.logging.terminal_label_width,
                    )
                if "ldsc" in resumed:
                    _log_ldsc_sample_prevalence(
                        logger, resumed["ldsc"], module,
                    )
                logger.record(
                    "SKIP", "formatter_run", reason="validated_resume",
                    formats=targets,
                )
                logger.record(
                    "STATUS", "formatter_run", status="COMPLETED",
                    formats=targets, resumed=True,
                )
                resume_lines = [screen_line(
                    "success",
                    "Formatter outputs validated; continuing from completed step",
                    indent=2,
                )]
                if "ldsc" in resumed:
                    resume_lines.extend([
                        "",
                        screen_line(
                            "analysis", "GWAS-VCF case fraction", indent=6,
                        ),
                    ])
                    resume_lines.extend(_ldsc_sample_prevalence_screen_fields(
                        resumed,
                        module.study_design.case_count_column,
                        module.study_design.control_count_column,
                        configuration.logging.terminal_label_width,
                        indent=10,
                    ))
                resume_lines.extend([
                    "",
                    screen_line("analysis", "Reports and logs", indent=6),
                    screen_field(
                        "success",
                        "Detailed HTML report",
                        str(html_report_path),
                        indent=10,
                        label_width=(
                            configuration.logging.terminal_label_width
                        ),
                        path_value=True,
                    ),
                    screen_field(
                        "info",
                        "Full formatter log",
                        str(log_path),
                        indent=10,
                        label_width=(
                            configuration.logging.terminal_label_width
                        ),
                        path_value=True,
                    ),
                ])
                if emit_terminal_summary:
                    print(
                        resume_lines[0]
                        if pipeline_mode else "\n".join(resume_lines)
                    )
                magma_progress = _magma_variant_input_pipeline_progress(
                    args, pipeline_mode=pipeline_mode,
                )
                if magma_progress is not None:
                    pipeline_progress, variant_input_stage = magma_progress
                    pipeline_progress.start(variant_input_stage)
                return resumed
        write_resolved_configuration(
            configuration,
            resolved_config_path,
            modules=tuple(resolved_modules),
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
        stage_total = len(targets) + 1
        pipeline_progress = getattr(args, "_pipeline_stage_progress", None)
        detailed_pipeline_progress = (
            pipeline_mode and pipeline_progress is not None
        )
        magma_progress = _magma_variant_input_pipeline_progress(
            args, pipeline_mode=pipeline_mode,
        )
        progress = StageProgress(
            "Formatter preparation progress",
            enabled=(
                emit_terminal_summary and configuration.logging.show_progress
                and not detailed_pipeline_progress
            ),
            outcome_label_width=configuration.logging.terminal_label_width,
        )
        load_title = "Validate and read the harmonised GWAS-VCF"
        if magma_progress is not None:
            pipeline_progress, variant_input_stage = magma_progress
            pipeline_progress.start(variant_input_stage)
        progress.start_step(1, stage_total, load_title)
        try:
            validated_header_evidence = (
                current_vcf_validation.get("harmonised")
                if isinstance(current_vcf_validation, dict)
                else None
            )
            with logger.step(
                1, stage_total, load_title, "load_harmonised_vcf",
            ) as step:
                loaded = load_harmonised_vcf(
                    vcf,
                    work_table,
                    resolved_bcftools,
                    module,
                    logger=logger,
                    return_header_evidence=True,
                    validated_header_evidence=validated_header_evidence,
                    validated_sample=(
                        indexed_vcf_validation.sample
                        if indexed_vcf_validation is not None
                        else None
                    ),
                )
                if isinstance(loaded, tuple):
                    frame, input_evidence = loaded
                else:
                    frame = loaded
                    input_evidence = {
                        "genome_build": "not recorded",
                        "postgwas_dataset_id": dataset_id,
                        "postgwas_version": "not recorded",
                        "postgwas_status": "not recorded",
                    }
                step.set_rows(frame.height, removed=0)
                step.output("canonical_variants", rows=frame.height)
        except BaseException:
            progress.fail_step(1, stage_total, load_title)
            raise
        embedded_dataset = str(input_evidence["postgwas_dataset_id"])
        if embedded_dataset != dataset_id:
            logger.warning(
                "Run dataset ID %s differs from the GWAS-VCF embedded dataset "
                "and sample ID %s. The run ID controls output naming; verify "
                "that this is the intended GWAS-VCF."
                % (dataset_id, embedded_dataset)
            )
        progress.complete_step(
            1,
            stage_total,
            load_title,
            outcome_fields=[
                ("count", "Total variants in input GWAS-VCF", frame.height),
                ("genetic", "Genome build", input_evidence["genome_build"]),
                ("info", "VCF embedded dataset/sample", embedded_dataset),
                ("success", "VCF structural validation", "passed"),
            ],
        )
        study_design_required_by = [
            target for target in selected
            if target in module.study_design.required_formats
        ]
        study_design = None
        if study_design_required_by:
            study_design = infer_study_design(
                frame,
                module.study_design.case_count_column,
                module.study_design.control_count_column,
            )
            logger.record(
                "DECIDE",
                "study_design",
                trait_type=study_design.trait_type,
                required_by=study_design_required_by,
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
                    "The study is binary because the configured case-count column "
                    "contains at least one value, but case/control counts are "
                    "incomplete for some variants. Count-dependent outputs will "
                    "exclude those variants."
                )
        else:
            logger.record(
                "SKIP",
                "study_design",
                reason="not_required_for_selected_formats",
                formats=selected,
            )

        trait_types = {
            target: study_design.trait_type
            for target in study_design_required_by
        } if study_design is not None else {}
        schema_reports = _configured_schema_reports(
            targets,
            module,
            trait_types,
            variant_id_observations,
        )
        _log_schema_reports(logger, schema_reports)
        if emit_terminal_summary and not pipeline_mode:
            _print_contract_summary(
                schema_reports,
                configuration.logging.terminal_label_width,
            )

        results = {}
        exporters = {**EXPORTERS, CUSTOM_OUTPUT_TARGET: export_custom}
        identifier_cache_uses = Counter(
            _identifier_selection_settings(module, target)
            for target in targets
            if not (
                target == "ldsc" and merge_alleles_file is not None
            )
        )
        identifier_cache = {}
        for number, target in enumerate(targets, 2):
            identifier_type, duplicate_policy = _identifier_selection_settings(
                module, target,
            )
            reference_selection = (
                target == "ldsc" and merge_alleles_file is not None
            )
            if reference_selection and identifier_type != "rsid":
                raise FormattingError(
                    "LDSC --merge-alleles matching requires rsid identifiers; "
                    "remove --variant-id-type unique or select rsid."
                )
            if reference_selection:
                target_frame, identifier_qc = select_variant_identifiers(
                    frame,
                    module,
                    identifier_type,
                    duplicate_policy="allow",
                )
                target_frame, reference_qc = select_ldsc_reference_variants(
                    target_frame,
                    module,
                    merge_alleles_file,
                )
                identifier_qc.update(reference_qc)
                reference_frame_rows = target_frame.height
                target_frame, duplicate_qc = resolve_duplicate_identifiers(
                    target_frame,
                    module,
                    identifier_type,
                    duplicate_policy=duplicate_policy,
                )
                identifier_qc.update(duplicate_qc)
                identifier_qc["identifier_rows_out"] = target_frame.height
                identifier_qc["identifier_rows_excluded"] = (
                    int(identifier_qc["identifier_missing_rows_excluded"])
                    + int(identifier_qc["identifier_duplicate_rows_excluded"])
                )
                identifier_qc.update(validate_unique_identifiers(
                    target_frame,
                    module,
                    identifier_type,
                    duplicate_policy,
                ))
                identifier_selection_reused = False
                logger.record(
                    "DECIDE",
                    "ldsc_reference_selection",
                    target=target,
                    **reference_qc,
                )
            else:
                cache_key = (identifier_type, duplicate_policy)
                cached = identifier_cache.get(cache_key)
                if cached is None:
                    target_frame, identifier_qc = select_variant_identifiers(
                        frame,
                        module,
                        identifier_type,
                        duplicate_policy=duplicate_policy,
                    )
                    identifier_qc.update(validate_unique_identifiers(
                        target_frame,
                        module,
                        identifier_type,
                        duplicate_policy,
                    ))
                    if identifier_cache_uses[cache_key] > 1:
                        identifier_cache[cache_key] = (
                            target_frame,
                            dict(identifier_qc),
                        )
                    identifier_selection_reused = False
                else:
                    target_frame, cached_qc = cached
                    identifier_qc = dict(cached_qc)
                    identifier_selection_reused = True
                identifier_cache_uses[cache_key] -= 1
                if identifier_cache_uses[cache_key] == 0:
                    identifier_cache.pop(cache_key, None)
                reference_frame_rows = target_frame.height
            identifier_qc["identifier_selection_reused"] = (
                identifier_selection_reused
            )
            logger.record(
                "DECIDE", "formatter_variant_identifiers",
                target=target, **identifier_qc,
            )
            logger.record(
                "VALIDATE",
                "formatter_unique_identifiers",
                status="PASSED",
                target=target,
                variant_id_type=identifier_type,
                duplicate_policy=duplicate_policy,
                variants=target_frame.height,
                selection_reused=identifier_selection_reused,
            )
            target_name = formatter_target_display_name(
                target, variant_id_observations,
            )
            stage_title = "Create and validate %s input files" % target_name
            progress.start_step(number, stage_total, stage_title)
            step_context = logger.step(
                number,
                stage_total,
                stage_title,
                exporters[target].__name__,
                rows_in=frame.height,
            )
            step = None
            try:
                step = step_context.__enter__()
                exporter_kwargs = {
                    "overwrite": configuration.run.overwrite,
                }
                if target == "magma" and magma_pipeline_resources is not None:
                    magma_configuration = magma_pipeline_resources.configuration
                    exporter_kwargs.update({
                        "magma_config": magma_configuration.modules.magma,
                        "ld_reference_prefix": (
                            magma_pipeline_resources.reference.ld_reference_prefix
                        ),
                        "analysis_scope": (
                            magma_pipeline_resources.reference.analysis_scope
                        ),
                        "identifier_qc": identifier_qc,
                    })
                if target in study_design_required_by:
                    exporter_kwargs["study_design"] = study_design
                result = exporters[target](
                    target_frame,
                    output_directory,
                    dataset_id,
                    module,
                    **exporter_kwargs,
                )
                if result.get("variant_preparation") is not None:
                    logger.record(
                        "RESULT",
                        "magma_variant_preparation",
                        **result["variant_preparation"]["qc"],
                    )
                    step.output(
                        "excluded_variants",
                        path=result["variant_preparation"][
                            "excluded_variants"
                        ],
                    )
                if target == "ldsc":
                    _log_ldsc_sample_prevalence(logger, result, module)
                excluded = int(result["rows_excluded"])
                chromosome_excluded = int(
                    result.get("rows_excluded_unconfigured_chromosome", 0)
                )
                rows_on_configured_chromosomes = (
                    target_frame.height - chromosome_excluded
                )
                missing_identifier_excluded = int(identifier_qc.get(
                    "identifier_missing_rows_excluded",
                    identifier_qc["identifier_rows_excluded"],
                ))
                rows_after_missing_identifiers = (
                    frame.height - missing_identifier_excluded
                )
                duplicate_policy_rows_in = (
                    reference_frame_rows
                    if reference_selection
                    else rows_after_missing_identifiers
                )
                if missing_identifier_excluded:
                    step.qc(
                        "%s variant identifiers" % identifier_type,
                        "Exclude records that cannot represent the selected "
                        "identifier convention.",
                        frame.height,
                        rows_after_missing_identifiers,
                        reason="missing_or_invalid_variant_identifier",
                        warn=True,
                    )
                reference_excluded = (
                    int(identifier_qc.get("rows_excluded_not_in_reference", 0))
                    + int(identifier_qc.get(
                        "rows_excluded_reference_allele_mismatch", 0,
                    ))
                )
                if reference_excluded:
                    step.qc(
                        "LDSC merge-alleles reference",
                        "Retain records matching both a reference rsID and its "
                        "strand-unambiguous allele pair.",
                        rows_after_missing_identifiers,
                        reference_frame_rows,
                        reason="not_in_reference_or_allele_mismatch",
                        warn=True,
                    )
                duplicate_identifier_excluded = int(identifier_qc.get(
                    "identifier_duplicate_rows_excluded", 0,
                ))
                if duplicate_identifier_excluded:
                    step.qc(
                        "duplicated %s identifiers" % identifier_type,
                        DUPLICATE_POLICY_DESCRIPTIONS[duplicate_policy],
                        duplicate_policy_rows_in,
                        target_frame.height,
                        reason="duplicate_variant_identifier",
                        warn=True,
                    )
                if chromosome_excluded:
                    step.observed(
                        "excluded_unconfigured_chromosomes",
                        target=target,
                        chromosomes=result.get("excluded_chromosomes", {}),
                        variants=chromosome_excluded,
                    )
                    step.qc(
                        "configured formatter chromosomes",
                        "Exclude records whose canonical chromosome is not "
                        "configured for this output.",
                        target_frame.height,
                        rows_on_configured_chromosomes,
                        reason="unconfigured_chromosome",
                        warn=True,
                    )
                step.qc(
                    "required formatter fields",
                    "Exclude records that cannot be represented correctly in this tool's input.",
                    rows_on_configured_chromosomes,
                    (
                        target_frame.height - excluded
                        if result.get("variant_preparation") is not None
                        else int(result["rows_out"])
                    ),
                    reason="missing_or_invalid_required_value",
                    warn=excluded > 0,
                )
                bounded = int(result.get("p_values_bounded", 0))
                if bounded:
                    step.qc(
                        "raw p-value numeric bound",
                        "Bound raw p-values at the configured minimum to avoid "
                        "floating-point underflow.",
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
                result["input_vcf_metadata"] = dict(input_evidence)
                result["rows_in"] = frame.height
                result["schema_rows_excluded"] = excluded
                result["rows_excluded"] = total_excluded
                results[target] = result
            except BaseException as exc:
                try:
                    if step is not None:
                        step_context.__exit__(type(exc), exc, exc.__traceback__)
                finally:
                    progress.fail_step(number, stage_total, stage_title)
                raise
            else:
                step_context.__exit__(None, None, None)
            progress.complete_step(
                number,
                stage_total,
                stage_title,
                outcome_fields=[
                    ("success", "Variants written", int(result["rows_out"])),
                    ("count", "Variants excluded", total_excluded),
                ],
            )

        for result in results.values():
            result["html_report"] = str(html_report_path)
        write_formatter_html_report(
            html_report_path,
            dataset_id=dataset_id,
            input_vcf=vcf,
            output_directory=output_directory,
            study_design=study_design,
            results=results,
            schema_reports=schema_reports,
            output_destinations=output_destinations,
            log_path=log_path,
            resolved_config_path=resolved_config_path,
            completion_manifest_path=completion_manifest_path,
            case_count_column=module.study_design.case_count_column,
            control_count_column=module.study_design.control_count_column,
            variant_id_observations=variant_id_observations,
            resolved_bcftools=resolved_bcftools,
            module=module,
        )
        logger.record(
            "OUTPUT",
            "formatter_html_report",
            path=str(html_report_path),
            evidence_source="validated_formatter_result_metadata",
        )
        write_formatter_completion_manifest(
            output_directory=output_directory,
            dataset_id=dataset_id,
            vcf=vcf,
            selected=selected,
            configuration=configuration,
            results=results,
        )
        if ctx is not None:
            ctx["formatter"] = results
        logger.record(
            "STATUS", "formatter_run", status="COMPLETED",
            formats=targets, input_variants=frame.height,
        )
        if emit_terminal_summary and not pipeline_mode:
            print(_screen_summary(
                dataset_id,
                str(vcf),
                study_design,
                results,
                str(log_path),
                str(html_report_path),
                module.study_design.case_count_column,
                module.study_design.control_count_column,
                configuration.logging.terminal_label_width,
                variant_id_observations,
                output_destinations,
                input_evidence,
                module.minimum_p_value,
            ))
        return results
    except BaseException as exc:
        if not logger.summary()["failed"]:
            logger.error("Formatter failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        if progress is not None:
            progress.close()
        if work_table is not None:
            work_table.unlink(missing_ok=True)
        logger.close()


__all__ = ["EXPORTERS", "run_formatter_direct"]
