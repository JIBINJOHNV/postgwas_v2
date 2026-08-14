"""Application service for one-pass GWAS-VCF format export."""

from __future__ import annotations

from pathlib import Path
import tempfile
from typing import get_args

from rich import box
from rich.console import Console
from rich.table import Table
from rich.text import Text

from postgwas.config import (
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.formatting import FormattingSampleSizeRole
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.ui.screen import screen_field, screen_line

from .contracts import (
    CUSTOM_OUTPUT_TARGET,
    FORMAT_CONTRACTS,
    formatter_result_targets,
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
    resume_formatter_outputs,
    write_formatter_completion_manifest,
)
from .table import (
    FormattingError,
    infer_study_design,
    load_harmonised_vcf,
    resolve_duplicate_identifiers,
    select_variant_identifiers,
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


def _ldsc_sample_prevalence_screen_fields(
    results,
    case_count_column,
    control_count_column,
    label_width,
):
    """Render the sample prevalence already calculated by the LDSC exporter."""
    result = results.get("ldsc")
    if result is None:
        return []
    if result.get("trait_type") != "binary":
        return [screen_field(
            "analysis",
            "LDSC sample prevalence",
            "not applicable to a quantitative trait",
            indent=6,
            label_width=label_width,
        )]

    value = result.get("sample_prev")
    aggregation = result.get("sample_prevalence_aggregation")
    variants = result.get("sample_prevalence_variants")
    minimum = result.get("sample_prevalence_minimum")
    maximum = result.get("sample_prevalence_maximum")
    if any(
        item is None
        for item in (value, aggregation, variants, minimum, maximum)
    ):
        return [screen_field(
            "warning",
            "LDSC sample prevalence",
            "unavailable because the formatter result metadata is incomplete",
            indent=6,
            label_width=label_width,
        )]

    fields = [
        screen_field(
            "analysis",
            "LDSC sample prevalence",
            "%.10g (%.2f%%)" % (float(value), float(value) * 100.0),
            indent=6,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "Prevalence calculation",
            "%s of %s / (%s + %s)"
            % (
                aggregation,
                case_count_column,
                case_count_column,
                control_count_column,
            ),
            indent=6,
            label_width=label_width,
        ),
        screen_field(
            "count",
            "Prevalence variants",
            "%s exported variants" % f"{int(variants):,}",
            indent=6,
            label_width=label_width,
        ),
        screen_field(
            "info",
            "Prevalence range",
            "%.10g to %.10g" % (float(minimum), float(maximum)),
            indent=6,
            label_width=label_width,
        ),
    ]
    count_ranges = (
        result.get("sample_prevalence_case_count_minimum"),
        result.get("sample_prevalence_case_count_maximum"),
        result.get("sample_prevalence_control_count_minimum"),
        result.get("sample_prevalence_control_count_maximum"),
    )
    if any(value is None for value in count_ranges):
        fields.append(screen_field(
            "warning",
            "%s/%s ranges" % (case_count_column, control_count_column),
            "unavailable in stored completion metadata; rerun with --overwrite",
            indent=6,
            label_width=label_width,
        ))
        return fields

    case_minimum, case_maximum, control_minimum, control_maximum = count_ranges
    fields.extend([
        screen_field(
            "count",
            "%s range" % case_count_column,
            "%s to %s"
            % (
                format(float(case_minimum), ",.10g"),
                format(float(case_maximum), ",.10g"),
            ),
            indent=6,
            label_width=label_width,
        ),
        screen_field(
            "count",
            "%s range" % control_count_column,
            "%s to %s"
            % (
                format(float(control_minimum), ",.10g"),
                format(float(control_maximum), ",.10g"),
            ),
            indent=6,
            label_width=label_width,
        ),
    ])
    return fields


def _validate_output_destinations(
    output_directory: Path,
    dataset_id: str,
    selected: list[str],
    module,
) -> dict[str, str]:
    """Resolve every selected output and reject collisions before extraction."""
    output_paths = formatter_output_paths(
        output_directory, dataset_id, selected, module,
    )
    destinations: dict[Path, list[str]] = {}
    for label, path in output_paths.items():
        destinations.setdefault(path, []).append(label)
    if module.custom_output.active:
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
    ]
    if study_design is None:
        lines.append(screen_field(
            "info",
            "Study-design inference",
            "not required for the selected formats",
            indent=6,
            label_width=label_width,
        ))
    else:
        lines.append(screen_field(
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
        ))
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
        exclusion_counts = [
            (
                int(result.get(
                    "identifier_missing_rows_excluded",
                    result.get("identifier_rows_excluded", 0),
                )),
                "missing or invalid identifiers",
            ),
            (
                int(result.get("identifier_duplicate_rows_excluded", 0)),
                "rows with duplicated identifiers",
            ),
            (
                int(result.get("rows_excluded_not_in_reference", 0)),
                "not in the LDSC merge-alleles reference",
            ),
            (
                int(result.get(
                    "rows_excluded_reference_allele_mismatch", 0,
                )),
                "incompatible with LDSC reference alleles",
            ),
            (
                int(result.get(
                    "rows_excluded_unconfigured_chromosome", 0,
                )),
                "outside configured chromosomes",
            ),
            (
                int(result.get("schema_rows_excluded", 0)),
                "missing or invalid required values",
            ),
        ]
        exclusion_details = [
            "%s %s" % (f"{count:,}", reason)
            for count, reason in exclusion_counts
            if count
        ]
        categorized = sum(count for count, _ in exclusion_counts)
        uncategorized = int(result["rows_excluded"]) - categorized
        if uncategorized > 0:
            exclusion_details.append(
                "%s other exclusions" % f"{uncategorized:,}"
            )
        exclusion_summary = "%s excluded" % f"{result['rows_excluded']:,}"
        if exclusion_details:
            exclusion_summary += " (%s)" % "; ".join(exclusion_details)
        lines.append(screen_field(
            "success",
            target,
            "%s variants exported using %s; duplicate policy %s; %s"
            % (
                f"{result['rows_out']:,}",
                "rsIDs" if result["variant_id_type"] == "rsid" else "unique IDs",
                result.get("identifier_duplicate_policy", "not recorded"),
                exclusion_summary,
            ),
            indent=6,
            label_width=label_width,
        ))
    lines.extend(_ldsc_sample_prevalence_screen_fields(
        results,
        case_count_column,
        control_count_column,
        label_width,
    ))
    lines.append(screen_field(
        "info", "Full log", log_path, indent=6, label_width=label_width,
    ))
    lines.append("")
    return "\n".join(lines)


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


def _vcf_query_label(query):
    """Render a configured simple FORMAT query without hiding custom queries."""
    query = str(query).strip()
    if query.startswith("[%") and query.endswith("]"):
        tag = query[2:-1]
        if tag and all(
            character.isalnum() or character == "_" for character in tag
        ):
            return "FORMAT/%s" % tag
    return query


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
            vcf_source = _vcf_query_label(module.vcf_fields.root[source])
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


def _configured_schema_reports(targets, module, trait_types=None):
    """Build screen/log reports solely from resolved formatter configuration."""
    trait_types = trait_types or {}
    reports = {}
    canonical = module.canonical_columns
    for target in targets:
        outputs = _configured_output_schemas(
            module, target, trait_types.get(target),
        )
        variant_id_type = module.variant_identifiers.target_types.get(
            target,
            module.variant_identifiers.default_type,
        )
        duplicate_id_policy = (
            module.variant_identifiers.target_duplicate_policies.get(
                target,
                module.variant_identifiers.default_duplicate_policy,
            )
        )
        if target == CUSTOM_OUTPUT_TARGET:
            name = "Custom CLI table"
            frequency_interpretation = "custom field semantics"
        else:
            contract = FORMAT_CONTRACTS[target]
            name = contract.name
            frequency_interpretation = contract.frequency
        reports[target] = {
            "name": name,
            "variant_id_type": variant_id_type,
            "duplicate_id_policy": duplicate_id_policy,
            "variant_id_type_label": (
                "rsID"
                if variant_id_type == "rsid"
                else "unique"
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
            variant_id_type=report["variant_id_type"],
            duplicate_id_policy=report["duplicate_id_policy"],
            outputs=report["outputs"],
            p_value_saved=report["p_value"],
            allele_frequency_saved=report["frequency"],
            frequency_interpretation=report["frequency_interpretation"],
            sample_size_saved=report["sample_size"],
        )


def _print_contract_table(reports):
    """Show exact configured column mappings and statistic representations."""
    table = Table(
        title="🔬  Validated downstream input schemas",
        box=box.SIMPLE_HEAVY,
        show_lines=True,
        header_style="bold cyan",
        padding=(0, 1),
    )
    table.add_column("Tool", style="bold", no_wrap=True)
    table.add_column("Variant ID type / duplicate policy", width=18)
    table.add_column("Source → saved column", ratio=6)
    table.add_column("P value saved", ratio=2)
    table.add_column("Frequency saved", ratio=2)
    table.add_column("Sample size saved", ratio=3)
    for report in reports.values():
        table.add_row(
            report["name"],
            Text(
                "%s / %s"
                % (
                    report["variant_id_type_label"],
                    report["duplicate_id_policy"],
                )
            ),
            Text(report["column_mappings"]),
            Text(report["p_value"]),
            Text(report["frequency"]),
            Text(report["sample_size"]),
        )
    console = Console(width=144)
    console.print()
    console.print(table)


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
        )

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
        logger.record("PARAM", "formats", values=targets)
        logger.record(
            "VALIDATE",
            "formatter_output_destinations",
            status="PASSED",
            values=output_destinations,
        )
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
                    modules=("formatting",),
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
                    targets, module, resumed_trait_types,
                )
                _log_schema_reports(logger, schema_reports)
                _print_contract_table(schema_reports)
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
                resume_lines.extend(_ldsc_sample_prevalence_screen_fields(
                    resumed,
                    module.study_design.case_count_column,
                    module.study_design.control_count_column,
                    configuration.logging.terminal_label_width,
                ))
                print("\n".join(resume_lines))
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
            1, len(targets) + 1, "Read harmonised GWAS-VCF", "load_harmonised_vcf",
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
            targets, module, trait_types,
        )
        _log_schema_reports(logger, schema_reports)
        _print_contract_table(schema_reports)

        results = {}
        exporters = {**EXPORTERS, CUSTOM_OUTPUT_TARGET: export_custom}
        for number, target in enumerate(targets, 2):
            identifier_type = module.variant_identifiers.target_types.get(
                target, module.variant_identifiers.default_type,
            )
            reference_selection = (
                target == "ldsc" and merge_alleles_file is not None
            )
            if reference_selection and identifier_type != "rsid":
                raise FormattingError(
                    "LDSC --merge-alleles matching requires rsid identifiers; "
                    "remove --variant-id-type unique or select rsid."
                )
            duplicate_policy = (
                module.variant_identifiers.target_duplicate_policies.get(
                    target,
                    module.variant_identifiers.default_duplicate_policy,
                )
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
                logger.record(
                    "DECIDE",
                    "ldsc_reference_selection",
                    target=target,
                    **reference_qc,
                )
            else:
                target_frame, identifier_qc = select_variant_identifiers(
                    frame,
                    module,
                    identifier_type,
                    duplicate_policy=duplicate_policy,
                )
                reference_frame_rows = target_frame.height
            logger.record(
                "DECIDE", "formatter_variant_identifiers",
                target=target, **identifier_qc,
            )
            with logger.step(
                number,
                len(targets) + 1,
                "Create %s input" % target,
                exporters[target].__name__,
                rows_in=frame.height,
            ) as step:
                exporter_kwargs = {
                    "overwrite": configuration.run.overwrite,
                }
                if target in study_design_required_by:
                    exporter_kwargs["study_design"] = study_design
                result = exporters[target](
                    target_frame,
                    output_directory,
                    dataset_id,
                    module,
                    **exporter_kwargs,
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
                    duplicate_policy_descriptions = {
                        "exclude_all": (
                            "Exclude every conflicting record in a duplicated "
                            "identifier group."
                        ),
                        "error": (
                            "Collapse exact repeated records; stop if any "
                            "conflicting duplicated-identifier group remains."
                        ),
                        "most_significant": (
                            "Retain only a unique largest configured -log10(P); "
                            "exclude competing rows and whole tied or missing-rank groups."
                        ),
                        "highest_maf": (
                            "Retain only a unique largest MAF derived from EAF; "
                            "exclude competing rows and whole tied or missing-rank groups."
                        ),
                        "highest_info": (
                            "Retain only a unique largest imputation INFO value; "
                            "exclude competing rows and whole tied or missing-rank groups."
                        ),
                    }
                    step.qc(
                        "duplicated %s identifiers" % identifier_type,
                        duplicate_policy_descriptions[duplicate_policy],
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
                        "Exclude records whose normalized chromosome is not "
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
                    int(result["rows_out"]),
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
            formats=targets, input_variants=frame.height,
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
