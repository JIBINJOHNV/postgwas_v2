"""
PostGWAS Pipeline Runners (v56, Dynamic Directories & Logic Fixes)

Changes (NO LOGIC CHANGE):
• setup_subdir: robust step prefix normalization (02 vs 2), always creates folders
• runners: never reference `outputs` in finally if step failed
• restore args.output_directory safely even on error
"""

import os
import shutil
from pathlib import Path

from postgwas.config import load_configuration, load_module_configuration
from postgwas.config.models.modules.single_cell import (
    SINGLE_CELL_TOOL_PIPELINE_DEPENDENCIES,
    single_cell_formatter_targets,
)
from postgwas.modules.formatting.contracts import required_formats
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
    configure_required_variant_identifier_type,
)
from postgwas.core.input_validation import current_validation_session
from postgwas.core.paths import configured_output_path
from postgwas.core.paths import resolve_executable
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.values import format_count, format_percentage


def _startup_resource_outcome_fields(fields):
    """Avoid repeating the common startup inventory; retain its cautions.

    This applies only to static resource presentation. Callers must continue
    their validation and progress accounting, and must not use it for generated
    inputs, scientific compatibility decisions, or final result provenance.
    """
    if current_validation_session() is None:
        return fields
    return [field for field in fields if field[0] in {"warning", "error", "loss"}]


# ============================================================
# HELPER: Directory Manager (DYNAMIC)
# ============================================================


def _normalize_step_prefix(step):
    """
    Convert step to a 2-digit string prefix.
    Accepts: int(2), "2", "02", None
    """
    if step is None:
        return "00"
    # if already something like "02"
    try:
        # step may be "02" or "2" or 2
        i = int(str(step))
        return f"{i:02d}"
    except Exception:
        # fallback: if user injected something weird, keep as string
        s = str(step).strip()
        if s == "":
            return "00"
        return s


def setup_subdir(args, base_name):
    """
    Dynamically creates subdirectories based on execution order.
    Example: If base_name is 'formatter' and step is 2, creates '02_formatter'.
    """
    prefix = _normalize_step_prefix(getattr(args, "_step_num", None))
    folder_name = f"{prefix}_{base_name}"

    original_output_directory = args.output_directory
    if original_output_directory is None:
        raise ValueError("args.output_directory is None")

    # Ensure root outdir exists (important in docker-mounted paths)
    os.makedirs(original_output_directory, exist_ok=True)

    step_output_directory = os.path.join(original_output_directory, folder_name)
    os.makedirs(step_output_directory, exist_ok=True)

    # Update the canonical output directory temporarily for this pipeline step.
    args.output_directory = step_output_directory
    return original_output_directory


def _subdir_path(args, base_name, *, step=None):
    """Resolve a pipeline stage directory without changing the CLI namespace."""
    prefix = _normalize_step_prefix(
        getattr(args, "_step_num", None) if step is None else step
    )
    return Path(args.output_directory).expanduser().resolve() / (
        "%s_%s" % (prefix, base_name)
    )


# ============================================================
# HELPER: Dependency Debugger
# ============================================================

def check_and_resolve_binaries(args, required_tools):
    def inject_to_path(binary_path):
        directory = os.path.dirname(binary_path)
        if directory and directory not in os.environ.get("PATH", ""):
            os.environ["PATH"] = directory + os.pathsep + os.environ.get("PATH", "")

    if "bcftools" in required_tools:
        configuration = load_configuration(getattr(args, "run_config", None))
        cmd = configuration.resources.executables.bcftools
        resolved = shutil.which(cmd)
        if resolved:
            # Internal execution state only; bcftools is not a public CLI option.
            args.bcftools = resolved
            inject_to_path(resolved)
        else:
            raise FileNotFoundError(f"Missing executable: {cmd}")

    if "plink" in required_tools:
        cmd = getattr(args, "plink", "plink")
        resolved = shutil.which(cmd)
        if not resolved and cmd == "plink":
            resolved = shutil.which("plink2")
        if resolved:
            args.plink = resolved
            inject_to_path(resolved)
        else:
            raise FileNotFoundError(f"Missing executable: {cmd}")

    for tool in ["tabix", "bgzip"]:
        if tool in required_tools:
            path = shutil.which(tool)
            if path:
                inject_to_path(path)
            else:
                raise FileNotFoundError(f"Missing tool: {tool}")


def _validate_current_pipeline_vcf(args, ctx):
    """Validate and publish the current VCF after every data transition."""
    from postgwas.core.vcf import (
        FormattingError,
        validate_harmonised_indexed_vcf,
    )
    from postgwas.core.preflight import require_unchanged_preflight_files

    if not (
        hasattr(ctx, "validation") and hasattr(ctx, "publish_validation")
    ):
        return None
    configuration = load_configuration(getattr(args, "run_config", None))
    module = configuration.modules.formatting
    current = ctx.validation("current_vcf", {})
    cached = current.get("indexed") if isinstance(current, dict) else None
    entry = ctx.validation("input_vcf", {})
    entry_indexed = (
        entry.get("indexed") if isinstance(entry, dict) else None
    )
    executable_identity = next(
        (
            evidence.get("bcftools_identity")
            for evidence in (current, entry)
            if isinstance(evidence, dict)
            and evidence.get("bcftools_identity") is not None
        ),
        None,
    )
    if executable_identity is not None:
        require_unchanged_preflight_files(
            executable_identity,
            error_type=FormattingError,
            label="bcftools executable",
        )
    validated_executable = next(
        (
            evidence.bcftools
            for evidence in (cached, entry_indexed)
            if evidence is not None and str(evidence.bcftools).strip()
        ),
        None,
    )
    bcftools = validated_executable or resolve_executable(
        configuration.resources.executables.bcftools,
        "bcftools executable",
        error_type=FormattingError,
    )
    evidence = validate_harmonised_indexed_vcf(
        args.vcf,
        getattr(args, "dataset_id", None) or configuration.run.dataset_id,
        bcftools,
        module,
        cached=cached,
    )
    if executable_identity is not None:
        evidence["bcftools_identity"] = executable_identity
    ctx.publish_validation("current_vcf", evidence)
    return evidence


def _validated_magma_gene_result(ctx, consumer):
    """Return a calibrated primary MAGMA result or fail before invalid reuse."""
    magma = ctx.get("magma")
    if not isinstance(magma, dict):
        raise ValueError("%s requires a completed MAGMA analysis" % consumer)
    primary = magma.get("primary_mapping")
    analyses = magma.get("mapping_analyses")
    definition = analyses.get(primary) if isinstance(analyses, dict) else None
    statistic = (
        definition.get("result_statistic_type")
        if isinstance(definition, dict)
        else None
    )
    if statistic != "calibrated_gene_p_value":
        raise ValueError(
            "%s requires calibrated MAGMA gene results, but primary mapping %r "
            "provides %r. Choose a positional, eMAGMA, H-MAGMA, or nMAGMA "
            "mapping with --primary-magma-mapping."
            % (consumer, primary, statistic)
        )
    return magma


def _formatter_consumers_for_current_vcf(args, ctx):
    """Select only consumers of the VCF state handled by this formatter stage.

    Imputation is a real data transition.  Its first formatter stage must create
    only PRED-LD input from the original VCF; a later formatter stage must create
    downstream inputs from the newly harmonised imputed VCF.  Formatting every
    active consumer at both boundaries duplicates full-VCF extraction and can
    leave scientifically obsolete pre-imputation tables beside final inputs.
    """
    active_modules = list(dict.fromkeys(
        getattr(args, "modules", None) or (),
    ))
    if (
        getattr(args, "apply_imputation", False)
        and "imputation" not in active_modules
    ):
        active_modules.append("imputation")
    if "imputation" not in active_modules:
        return active_modules, False
    if "imputation" not in ctx:
        return ["imputation"], True
    return [
        module for module in active_modules
        if module != "imputation"
    ], False


# ============================================================
# RUNNERS
# ============================================================

def run_annot_ldblock_runner(args, ctx):
    from postgwas.modules.ld_annotation.service import run_annot_ldblock

    root = setup_subdir(args, "annot_ldblock")
    outputs = None
    try:
        outputs = run_annot_ldblock(args)
        ctx["annot_ldblock"] = outputs
        # keep behavior: next steps can rely on args.vcf being updated
        if isinstance(outputs, dict) and "annotated_vcf" in outputs:
            args.vcf = outputs["annotated_vcf"]
            _validate_current_pipeline_vcf(args, ctx)
        return outputs
    finally:
        args.output_directory = root


def run_sumstat_filter_runner(args, ctx):
    from postgwas.modules.filtering.service import run_sumstat_filter_direct

    # Determine base name based on context
    base_name = "filter_post_imp" if "imputation" in ctx else "filter_pre_imp"

    root = setup_subdir(args, base_name)
    outputs = None
    try:
        outputs = run_sumstat_filter_direct(args)
        if "imputation" in ctx:
            ctx["post_imputation_filter"] = outputs
        else:
            ctx["sumstat_filter"] = outputs

        # keep behavior: args.vcf updated
        if isinstance(outputs, dict) and "filtered_vcf" in outputs:
            args.vcf = outputs["filtered_vcf"]
            _validate_current_pipeline_vcf(args, ctx)
        return outputs
    finally:
        args.output_directory = root


def run_formatter_runner(args, ctx):
    from postgwas.modules.formatting.service import (
        run_formatter_direct,
    )
    from postgwas.modules.magma.service import MagmaPipelineResources
    from postgwas.modules.magma.reporting import magma_pipeline_stage_numbers
    from postgwas.modules.gcta_gene.service import GctaGenePipelineResources

    active_modules, pre_imputation_formatter = (
        _formatter_consumers_for_current_vcf(args, ctx)
    )
    current_vcf_evidence = _validate_current_pipeline_vcf(args, ctx)
    magma_validation = (
        ctx.validation("magma") if hasattr(ctx, "validation") else None
    )
    magma_resources = (
        magma_validation.resources if magma_validation is not None else None
    )
    if magma_resources is not None and not isinstance(
        magma_resources, MagmaPipelineResources,
    ):
        raise RuntimeError("Pipeline MAGMA preflight evidence has an invalid type")
    gcta_validation = (
        ctx.validation("gcta_gene") if hasattr(ctx, "validation") else None
    )
    gcta_resources = (
        gcta_validation.resources if gcta_validation is not None else None
    )
    if "gcta_gene" in active_modules and not isinstance(
        gcta_resources, GctaGenePipelineResources,
    ):
        raise RuntimeError(
            "GCTA gene formatter preparation requires completed pipeline "
            "resource preflight evidence."
        )
    if "ld_clump" in active_modules:
        ld_clumping = load_module_configuration(
            "ld_clumping",
            getattr(args, "run_config", None),
            cli_overrides=(
                {"methods": args.clumping_methods}
                if hasattr(args, "clumping_methods")
                else None
            ),
        )
        if (
            "cojo-slct" in ld_clumping.methods
            and "gcta_cojo" not in active_modules
        ):
            # The formatter needs the existing GCTA consumer contract, but the
            # standalone GCTA pipeline module is not added to the execution plan.
            active_modules.append("gcta_cojo")
    formatting_config = load_module_configuration(
        "formatting", getattr(args, "run_config", None),
    )
    requested_formats = (
        []
        if pre_imputation_formatter
        else list(getattr(args, "format", None) or formatting_config.formats)
    )
    if "single_cell" in active_modules:
        single_cell = load_module_configuration(
            "single_cell",
            getattr(args, "run_config", None),
            cli_overrides=(
                {"tools": args.tools} if hasattr(args, "tools") else None
            ),
        )
        requested_formats.extend(
            single_cell_formatter_targets(single_cell.tools)
        )
    args.format = required_formats(
        formatting_config,
        active_modules,
        requested_formats,
    )
    detailed_progress = getattr(args, "_pipeline_stage_progress", None)
    detailed_plan = getattr(args, "_pipeline_progress_plan", {}) or {}
    detailed_kind = detailed_plan.get("kind")
    if detailed_kind is None and detailed_progress is not None:
        # Preserve programmatic callers that attach the historical MAGMA
        # controller directly rather than going through the executor.
        detailed_kind = "magma"
    magma_only_progress = (
        not pre_imputation_formatter
        and detailed_kind in {"magma", "magmacovar", "pops"}
    )
    magma_stages = magma_pipeline_stage_numbers(args) if magma_only_progress else None
    gcta_only_progress = (
        not pre_imputation_formatter and detailed_kind == "gcta_gene"
    )
    gcta_configuration = (
        gcta_resources.configuration if gcta_resources is not None else None
    )
    gcta_module = (
        gcta_configuration.modules.gcta_gene
        if gcta_configuration is not None else None
    )
    gcta_log_path = None
    if magma_only_progress or gcta_only_progress:
        if magma_only_progress:
            detailed_progress.start(magma_stages["vcf"])
        else:
            from postgwas.modules.gcta_gene.reporting import (
                gcta_gene_overlap_outcome_fields,
                gcta_gene_reference_outcome_fields,
                gcta_gmt_outcome_fields,
                gcta_ld_reference_outcome_fields,
                gcta_ma_input_outcome_fields,
                gcta_set_source_outcome_fields,
            )
            from postgwas.modules.gcta_gene.stages import (
                configure_pipeline_stage_callbacks,
                start_pipeline_stage,
            )

            if gcta_resources is None:
                raise RuntimeError(
                    "Detailed GCTA progress requires completed pipeline "
                    "resource preflight evidence."
                )
            next_step = int(getattr(args, "_step_num", 1)) + 1
            gcta_output = _subdir_path(args, "gcta_gene", step=next_step)
            gcta_log_path = configured_output_path(
                gcta_output,
                gcta_module.output_layout.log_file,
                error_type=RuntimeError,
                dataset_id=gcta_configuration.run.dataset_id,
                method=gcta_module.method,
            )

            def start_gcta_preflight(key):
                with PipelineLogger(
                    gcta_configuration.run.dataset_id,
                    "run",
                    str(gcta_log_path.parent),
                    level=gcta_configuration.logging.file_level,
                    screen_level=gcta_configuration.logging.console_level,
                    log_path=str(gcta_log_path),
                ) as preflight_logger:
                    start_pipeline_stage(args, key, preflight_logger)

            def complete_gcta_preflight(key, fields, *, startup_resource=False):
                from postgwas.modules.gcta_gene.stages import (
                    complete_pipeline_stage,
                )

                with PipelineLogger(
                    gcta_configuration.run.dataset_id,
                    "run",
                    str(gcta_log_path.parent),
                    level=gcta_configuration.logging.file_level,
                    screen_level=gcta_configuration.logging.console_level,
                    log_path=str(gcta_log_path),
                ) as preflight_logger:
                    complete_pipeline_stage(
                        args,
                        key,
                        outcome_fields=fields,
                        screen_outcome_fields=(
                            _startup_resource_outcome_fields(fields)
                            if startup_resource else None
                        ),
                        logger=preflight_logger,
                    )

            configure_pipeline_stage_callbacks(
                args, gcta_configuration, gcta_log_path,
            )
            start_gcta_preflight("vcf")
    if not args.vcf or not os.path.exists(args.vcf):
        raise FileNotFoundError("Input VCF missing: %s" % args.vcf)
    if magma_only_progress or gcta_only_progress:
        if current_vcf_evidence is None:
            raise RuntimeError(
                "Detailed pipeline progress requires reusable VCF validation evidence"
            )
        evidence = current_vcf_evidence["harmonised"]
        if magma_only_progress:
            # Pipeline MAGMA derives its analysis build from the validated VCF.
            # The MAGMA schema subsequently verifies that selected mapping
            # resources declare this same build before scientific analysis.
            args.genome_build = evidence["genome_build"]
        elif str(gcta_module.genome_build) != str(evidence["genome_build"]):
            raise ValueError(
                "The harmonised GWAS-VCF declares genome build %s, but "
                "--genome-build resolves to %s. The GWAS, PLINK LD reference, "
                "and gene-coordinate file must use the same build."
                % (evidence["genome_build"], gcta_module.genome_build)
            )
        variant_count = current_vcf_evidence["indexed"].variant_count
        if magma_only_progress:
            args._magma_vcf_variant_count = variant_count
        else:
            args._gcta_vcf_variant_count = variant_count
        vcf_fields = [
                ("analysis", "Input summary-statistics VCF"),
                ("info", "Input file", Path(args.vcf).name),
                ("count", "Total variants", variant_count),
                ("genetic", "Genome build inferred from VCF", evidence["genome_build"]),
                (
                    "info", "VCF embedded dataset/sample",
                    evidence["postgwas_dataset_id"],
                ),
                ("success", "VCF structural validation", "passed"),
        ]
        if magma_only_progress:
            detailed_progress.complete(
                magma_stages["vcf"], outcome_fields=vcf_fields,
            )
        else:
            complete_gcta_preflight("vcf", vcf_fields)
    requirements = []
    if any(
        "magma" in formatting_config.module_formats.get(module, ())
        for module in active_modules
    ):
        if magma_resources is None:
            raise RuntimeError(
                "MAGMA formatter preparation requires completed pipeline "
                "resource preflight evidence."
            )
        magma_config = magma_resources.configuration.modules.magma
        reference_prefix = magma_resources.reference.ld_reference_prefix
        requirements.append(BimIdentifierRequirement(
            consumer="MAGMA",
            formatter_target="magma",
            bim_file="%s%s" % (
                reference_prefix, magma_config.input.bim_extension,
            ),
            column_roles=magma_config.input.bim_columns,
            delimiter_pattern=magma_config.input.table_delimiter_pattern,
        ))

    for module_name, label in (("gcta_gene", "GCTA gene"), ("gcta_cojo", "GCTA COJO")):
        if module_name not in active_modules:
            continue
        if module_name == "gcta_cojo":
            from postgwas.modules.gcta_cojo.service import (
                GctaCojoPipelineResources,
                gcta_cojo_analysis_label,
            )

            cojo_validation = (
                ctx.validation("gcta_cojo")
                if hasattr(ctx, "validation") else None
            )
            cojo_resources = (
                cojo_validation.resources
                if cojo_validation is not None else None
            )
            if isinstance(cojo_resources, GctaCojoPipelineResources):
                module_config = cojo_resources.configuration.modules.gcta_cojo
                cojo_reference = cojo_resources.reference
            elif cojo_resources is None:
                ld_validation = (
                    ctx.validation("ld_clump")
                    if hasattr(ctx, "validation") else None
                )
                ld_resources = (
                    ld_validation.resources
                    if ld_validation is not None else None
                )
                if getattr(ld_resources, "cojo_reference_validation", None):
                    module_config = ld_resources.configuration.modules.gcta_cojo
                    cojo_reference = ld_resources.cojo_reference_validation
                else:
                    module_config = None
                    cojo_reference = None
            else:
                raise RuntimeError(
                    "Pipeline GCTA-COJO preflight evidence has an invalid type"
                )
            if module_config is None or cojo_reference is None:
                raise RuntimeError(
                    "GCTA-COJO formatter preparation requires completed "
                    "pipeline resource preflight evidence."
                )
            label = gcta_cojo_analysis_label(module_config)
        else:
            if gcta_resources is None:
                raise RuntimeError(
                    "GCTA gene formatter preparation requires completed "
                    "pipeline resource preflight evidence."
                )
            module_config = gcta_resources.configuration.modules.gcta_gene
        argument_name = (
            "gcta_reference_prefix"
            if module_name == "gcta_gene"
            else "cojo_reference_prefix"
        )
        reference_prefix = (
            cojo_reference.prefix
            if module_name == "gcta_cojo"
            else gcta_resources.reference_prefix
        )
        if not reference_prefix:
            raise ValueError(
                "%s requires --%s before formatting so the BIM identifier "
                "convention can be selected."
                % (label, argument_name.replace("_", "-"))
            )
        bim_extensions = [
            suffix for suffix in module_config.reference.required_extensions
            if suffix.lower() == ".bim"
        ]
        if len(bim_extensions) != 1:
            raise ValueError(
                "%s reference.required_extensions must contain exactly one .bim suffix."
                % label
            )
        requirements.append(BimIdentifierRequirement(
            consumer=label,
            formatter_target="gcta_gene",
            bim_file="%s%s" % (reference_prefix, bim_extensions[0]),
            column_roles=module_config.reference.bim_columns,
            delimiter_pattern=module_config.reference.table_delimiter_pattern,
        ))

    if requirements:
        if gcta_only_progress:
            start_gcta_preflight("reference")
            gcta_reference_prefix = gcta_resources.reference_prefix
            gcta_reference_paths = gcta_resources.reference_paths
        if magma_only_progress:
            detailed_progress.start(magma_stages["ld_reference"])
            if magma_resources is None:
                raise RuntimeError(
                    "Detailed MAGMA progress requires completed pipeline "
                    "resource preflight evidence."
                )
            magma_configuration = magma_resources.configuration
            module = magma_configuration.modules.magma
            references = magma_resources.reference
            ld_reference_prefix = references.ld_reference_prefix
        configure_reference_variant_identifiers(
            args, formatting_config, requirements,
        )
        if gcta_only_progress:
            observation = gcta_resources.identifier_observation
            fields = gcta_ld_reference_outcome_fields(
                gcta_module,
                gcta_reference_prefix,
                gcta_reference_paths,
                reference_metrics={
                    "reference_variants": observation["variants"],
                },
                variant_id_type=observation["variant_id_type"],
            )
            complete_gcta_preflight(
                "reference", fields, startup_resource=True,
            )
        if magma_only_progress:
            observation = args.variant_id_observations["magma"]
            identifier_label = (
                "rsID"
                if observation["variant_id_type"] == "rsid"
                else "chromosome-position-allele ID"
            )
            detailed_progress.complete(
                magma_stages["ld_reference"],
                outcome_fields=_startup_resource_outcome_fields([
                    ("analysis", "PLINK LD reference"),
                    (
                        "info", "Reference prefix",
                        Path(ld_reference_prefix).name,
                    ),
                    (
                        "success", "Required companion files",
                        "%s present and non-empty"
                        % ", ".join(
                            extension.removeprefix(".").upper()
                            for extension in module.input.required_reference_extensions
                        ),
                    ),
                    ("count", "Total BIM variants", observation["variants"]),
                    ("genetic", "Variant identifiers", identifier_label),
                    ("genetic", "Declared genome build", module.genome_build.value),
                    ("genetic", "Declared population", module.population.value),
                    ("success", "BIM structural validation", "passed"),
                    (
                        "warning", "Build and population provenance",
                        "taken from configuration; genome build and population "
                        "are not independently verifiable from PLINK file contents",
                    ),
                ]),
            )
            detailed_progress.start(magma_stages["gene_location"])
            mapping_name = module.mapping.primary
            definition = module.mapping.definitions[mapping_name]
            location_metrics = magma_resources.primary_gene_location_metrics
            if location_metrics is None:
                raise RuntimeError(
                    "Detailed positional MAGMA progress requires validated "
                    "gene-location metrics."
                )
            location_file = location_metrics["path"]
            location_roles = module.input.gene_location_columns
            alternate_column = (
                location_roles.index("alternate_gene_id") + 1
                if "alternate_gene_id" in location_roles else None
            )
            detailed_progress.complete(
                magma_stages["gene_location"],
                outcome_fields=[
                    ("analysis", "Gene-location reference"),
                    ("info", "Reference file", Path(location_file).name),
                    (
                        "info", "Biological context",
                        definition.context or "not provided",
                    ),
                    ("count", "Reference genes", location_metrics["genes"]),
                    (
                        "genetic", "Primary IDs · column 1",
                        "%s unique" % format_count(location_metrics["genes"]),
                    ),
                    (
                        "genetic",
                        (
                            "Alternate IDs · column %d" % alternate_column
                            if alternate_column is not None else
                            "Alternate IDs"
                        ),
                        (
                            "%s unique"
                            % format_count(
                                location_metrics["alternate_unique_ids"]
                            )
                            if alternate_column is not None else "not configured"
                        ),
                    ),
                    (
                        (
                            "warning"
                            if location_metrics["ambiguous_alternate_ids"]
                            else "success"
                        ),
                        "Duplicated alternate IDs",
                        location_metrics["ambiguous_alternate_ids"],
                    ),
                    ("success", "File-structure validation", "passed"),
                ],
            )
            include_gene_sets = "pathway_file" in magma_stages
            if include_gene_sets:
                detailed_progress.start(magma_stages["pathway_file"])
            if magma_resources.include_gene_sets != include_gene_sets:
                raise RuntimeError(
                    "MAGMA pathway-analysis scope changed after pipeline preflight"
                )
            duplicate_policy = {
                "err": "error",
                "lowest_p": "lowest_p",
                "remove": "exclude_all",
            }[module.snp_harmonisation.duplicate_policy]
            duplicate_policies = dict(
                getattr(args, "variant_id_duplicate_policies", None) or {}
            )
            duplicate_policies["magma"] = duplicate_policy
            args.variant_id_duplicate_policies = duplicate_policies
            plan = references.gene_set_plans[mapping_name]
            pathway_metadata = plan.get("input_metadata") or {}
            resolution = (plan.get("validation") or {}).get(
                "identifier_resolution", {},
            )
            pathway_requested = plan["status"] != "not_requested"
            if include_gene_sets:
                detailed_progress.complete(
                    magma_stages["pathway_file"],
                    outcome_fields=(
                        [
                            ("analysis", "Pathway file"),
                            (
                                "info", "Pathway file",
                                Path(plan["pathway_file"]).name,
                            ),
                            (
                                "info", "File format",
                                str(pathway_metadata["detected_format"]).upper(),
                            ),
                            (
                                "count", "Total pathway records",
                                pathway_metadata["gene_sets"],
                            ),
                            (
                                "genetic", "Unique gene identifiers",
                                resolution["input_unique_ids"],
                            ),
                            ("success", "File-structure validation", "passed"),
                        ]
                        if pathway_requested else
                        [
                            ("analysis", "Pathway file"),
                            ("info", "Pathway analysis", "not requested"),
                        ]
                    ),
                )
                detailed_progress.start(magma_stages["pathway_identifiers"])
            if not pathway_requested:
                identifier_fields = [
                    ("analysis", "Gene-identifier compatibility"),
                    ("info", "Assessment", "not required; no pathway file supplied"),
                ]
            else:
                primary_matches = resolution["direct_primary_matches"]
                input_ids = resolution["input_unique_ids"]
                alternate_matches = resolution["alternate_candidate_matches"]
                alternate_reference_ids = resolution[
                    "location_alternate_unique_ids"
                ]
                resolved_alternate_column = resolution[
                    "alternate_identifier_column"
                ]
                alternate_label = (
                    "Gene-location alternate IDs · column %d"
                    % resolved_alternate_column
                    if resolved_alternate_column is not None else
                    "Gene-location alternate IDs"
                )
                alternate_pathway_label = (
                    "Pathway-file IDs found in gene-location column %d"
                    % resolved_alternate_column
                    if resolved_alternate_column is not None else
                    "Pathway-file IDs found among gene-location alternate IDs"
                )
                alternate_reference_label = (
                    "Gene-location column %d IDs found in pathway file"
                    % resolved_alternate_column
                    if resolved_alternate_column is not None else
                    "Gene-location alternate IDs found in pathway file"
                )
                selected_source = resolution["identifier_source"]
                if plan["status"] == "ready" and selected_source == "alternate_gene_id":
                    decision = (
                        "use gene-location alternate IDs from column %d as the "
                        "MAGMA primary gene IDs" % resolved_alternate_column
                    )
                    reason = (
                        "primary-ID overlap is below the configured minimum; "
                        "alternate-ID overlap passes"
                    )
                    derived_reference = (
                        "column %d IDs will become primary IDs; original column 1 "
                        "IDs will be retained as alternate metadata"
                        % resolved_alternate_column
                    )
                elif plan["status"] == "ready":
                    decision = "use gene-location primary IDs from column 1"
                    reason = "direct primary-ID overlap passes the configured minimum"
                    derived_reference = "not required"
                else:
                    decision = "no compatible gene-location identifier column found"
                    reason = (
                        "primary-ID overlap is below the configured minimum and "
                        "no usable alternate-ID column is available"
                        if not resolution.get("alternate_identifiers_available") else
                        "both primary-ID and alternate-ID overlap are below the "
                        "configured minimum"
                    )
                    derived_reference = "not created"
                primary_match_kind = (
                    "success"
                    if plan["status"] == "ready"
                    and selected_source == "primary_gene_id" else "loss"
                )
                alternate_match_kind = (
                    "success"
                    if plan["status"] == "ready"
                    and selected_source == "alternate_gene_id" else "info"
                )
                duplicate_count = resolution[
                    "location_ambiguous_alternate_ids"
                ]
                unmatched_count = resolution["unmatched_input_ids"]
                identifier_fields = [
                    ("analysis", "Gene-identifier compatibility"),
                    (
                        "count", "Pathway identifier universe",
                        "%s unique identifiers" % format_count(input_ids),
                    ),
                    (
                        "genetic", "Gene-location primary IDs · column 1",
                        "%s unique identifiers" % format_count(
                            resolution["location_primary_unique_ids"]
                        ),
                    ),
                    (
                        "genetic",
                        alternate_label,
                        "%s unique identifiers" % format_count(
                            alternate_reference_ids
                        ),
                    ),
                    ("analysis", "Direct compatibility"),
                    (
                        primary_match_kind,
                        "Pathway-file IDs found in gene-location column 1",
                        "%s/%s (%s)"
                        % (
                            format_count(primary_matches), format_count(input_ids),
                            format_percentage(primary_matches, input_ids),
                        ),
                    ),
                    (
                        alternate_match_kind,
                        alternate_pathway_label,
                        "%s/%s (%s)"
                        % (
                            format_count(alternate_matches), format_count(input_ids),
                            format_percentage(alternate_matches, input_ids),
                        ),
                    ),
                    (
                        alternate_match_kind,
                        alternate_reference_label,
                        "%s/%s (%s)"
                        % (
                            format_count(alternate_matches),
                            format_count(alternate_reference_ids),
                            format_percentage(
                                alternate_matches, alternate_reference_ids,
                            ),
                        ),
                    ),
                    ("decision", "Identifier decision"),
                    (
                        "success" if plan["status"] == "ready" else "warning",
                        "Selected action",
                        decision,
                    ),
                    ("info", "Reason", reason),
                    ("info", "Derived gene-location reference", derived_reference),
                    (
                        "warning" if duplicate_count else "success",
                        "Duplicated alternate IDs",
                        duplicate_count,
                    ),
                    (
                        "info", "Duplicate handling",
                        module.gene_sets.alternate_id_duplicate_policy,
                    ),
                    (
                        "warning" if unmatched_count else "success",
                        "Unmatched pathway identifiers",
                        unmatched_count,
                    ),
                    (
                        "info", "Required compatibility",
                        format_percentage(plan["minimum_overlap"], 1),
                    ),
                    ("success", "Original pathway identifiers", "unchanged"),
                    (
                        "success" if plan["status"] == "ready" else "warning",
                        "Pathway analysis",
                        (
                            "available"
                            if plan["status"] == "ready" else
                            "not performed; gene-association analysis will continue"
                        ),
                    ),
                ]
            if include_gene_sets:
                detailed_progress.complete(
                    magma_stages["pathway_identifiers"],
                    outcome_fields=identifier_fields,
                )
    if gcta_only_progress:
        from postgwas.modules.gcta_gene.stages import pipeline_stage_number

        if pipeline_stage_number(args, "gene_list") is not None:
            gene_list = gcta_resources.gene_list
            gene_metrics = gcta_resources.gene_annotation_metrics
            if gene_list is None or gene_metrics is None:
                raise RuntimeError(
                    "GCTA gene-coordinate preflight evidence is incomplete."
                )
            start_gcta_preflight("gene_list")
            complete_gcta_preflight(
                "gene_list",
                gcta_gene_reference_outcome_fields(
                    gcta_module,
                    gene_list,
                    gene_metrics,
                ),
            )

        if pipeline_stage_number(args, "gmt") is not None:
            gmt = gcta_resources.gmt
            pathway_preflight = gcta_resources.pathway_gene_preflight
            if gmt is None or pathway_preflight is None:
                raise RuntimeError("GCTA GMT preflight evidence is incomplete.")
            start_gcta_preflight("gmt")
            complete_gcta_preflight(
                "gmt",
                gcta_gmt_outcome_fields(
                    gcta_module,
                    gmt,
                    pathway_preflight.pathway_count,
                    len(pathway_preflight.requested_genes),
                ),
            )
            start_gcta_preflight("gene_overlap")
            complete_gcta_preflight(
                "gene_overlap",
                gcta_gene_overlap_outcome_fields(
                    gcta_module,
                    pathway_preflight,
                ),
            )

        if pipeline_stage_number(args, "set_list") is not None:
            set_list = gcta_resources.set_list
            set_metrics = gcta_resources.set_source_metrics
            if set_list is None or set_metrics is None:
                raise RuntimeError(
                    "GCTA fastBAT set-list preflight evidence is incomplete."
                )
            start_gcta_preflight("set_list")
            complete_gcta_preflight(
                "set_list",
                gcta_set_source_outcome_fields(set_list, set_metrics),
            )
    if "ldsc" in args.format:
        configure_required_variant_identifier_type(
            args,
            formatting_config,
            consumer="LDSC HapMap3 summary-statistic munging",
            formatter_target="ldsc",
            required_type="rsid",
        )

    root = setup_subdir(args, "formatter")
    try:
        if gcta_only_progress:
            start_gcta_preflight("formatted_input")
        formatter_result = run_formatter_direct(args, ctx)
        if not gcta_only_progress:
            return formatter_result

        from postgwas.core.ui import print_screen_block
        from postgwas.modules.gcta_gene.service import validate_pipeline_gcta_input

        formatted = formatter_result.get("gcta_gene", {})
        input_file = formatted.get("summary_statistics_input_file")
        if not input_file:
            raise ValueError(
                "Formatter did not return the required GCTA .ma input artifact."
            )
        gcta_input_preflight = validate_pipeline_gcta_input(
            input_file,
            gcta_configuration,
            pipeline_resources=gcta_resources,
        )
        ctx.publish_validation("gcta_gene_input", gcta_input_preflight)
        input_metrics = gcta_input_preflight.input_metrics
        reference_metrics = gcta_input_preflight.reference_metrics
        absent = reference_metrics["input_variants_absent_from_reference"]
        fields = gcta_ma_input_outcome_fields(
            gcta_configuration,
            input_file,
            input_metrics,
            formatter_result=formatted,
            variant_id_type=formatted.get("variant_id_type"),
            reference_metrics=reference_metrics,
        )
        complete_gcta_preflight("formatted_input", fields)
        with PipelineLogger(
            gcta_configuration.run.dataset_id,
            "run",
            str(gcta_log_path.parent),
            level=gcta_configuration.logging.file_level,
            screen_level=gcta_configuration.logging.console_level,
            log_path=str(gcta_log_path),
        ) as preflight_logger:
            preflight_logger.record(
                "RESULT",
                "gcta_input_validation_summary",
                summary_statistic_unique_ids=len(gcta_input_preflight.identifiers),
                plink_bim_unique_ids=reference_metrics["reference_variants"],
                exact_ids_shared=reference_metrics["overlapping_variants"],
                summary_ids_absent_from_bim=absent,
                bim_ids_absent_from_summary=reference_metrics[
                    "reference_variants_absent_from_input"
                ],
                compatible_allele_pairs=reference_metrics[
                    "compatible_allele_pairs"
                ],
                input_rewritten=False,
            )
        label_width = gcta_configuration.logging.terminal_label_width
        print_screen_block("\n".join([
            "",
            screen_line("analysis", "GCTA input-validation summary", indent=2),
            screen_field(
                "count", "Summary-statistic unique IDs",
                len(gcta_input_preflight.identifiers), indent=6,
                label_width=label_width,
            ),
            screen_field(
                "count", "PLINK BIM unique IDs",
                reference_metrics["reference_variants"], indent=6,
                label_width=label_width,
            ),
            screen_field(
                "success", "Exact IDs shared",
                reference_metrics["overlapping_variants"], indent=6,
                label_width=label_width,
            ),
            screen_field(
                "warning" if absent else "success",
                "Summary IDs absent from BIM",
                "%s; will not be used by GCTA" % format(absent, ",")
                if absent else "0",
                indent=6, label_width=label_width,
            ),
            screen_field(
                "count", "BIM IDs absent from summary",
                reference_metrics["reference_variants_absent_from_input"],
                indent=6, label_width=label_width,
            ),
            screen_field(
                "success", "BIM-compatible .ma",
                "validated without rewriting", indent=6,
                label_width=label_width,
            ),
            "",
        ]))
        return formatter_result
    finally:
        args.output_directory = root


def run_imputation_runner(args, ctx):
    from postgwas.modules.imputation.service import (
        PredLDPipelineResources,
        run_sumstat_imputation_direct,
    )

    root = setup_subdir(args, "imputation")

    if "formatter" not in ctx or "pred_ld" not in ctx["formatter"] or not ctx["formatter"]["pred_ld"]:
        raise ValueError("Imputation requires the 'pred_ld' format from formatter.")

    pred_ld_folder = ctx["formatter"]["pred_ld"].get("pred_ld_folder")
    if not pred_ld_folder:
        raise ValueError("Formatter did not return a valid 'pred_ld_folder' path.")

    args.pred_ld_input_directory = pred_ld_folder
    preflight = (
        ctx.validation("imputation") if hasattr(ctx, "validation") else None
    )
    resources = preflight.resources if preflight is not None else None
    if not isinstance(resources, PredLDPipelineResources):
        raise RuntimeError(
            "Imputation execution requires completed pipeline resource "
            "preflight evidence."
        )

    try:
        current_validation = (
            ctx.validation("current_vcf", {})
            if hasattr(ctx, "validation")
            else {}
        )
        indexed_input = (
            current_validation.get("indexed")
            if isinstance(current_validation, dict)
            else None
        )
        outputs = run_sumstat_imputation_direct(
            args,
            pipeline_resources=resources,
        )
        output_build = (
            indexed_input.genome_build
            if indexed_input is not None
            else getattr(args, "genome_build", None)
        )
        output_vcf = (
            outputs.get(output_build)
            if isinstance(outputs, dict) and output_build is not None
            else None
        )
        if not output_vcf:
            raise ValueError(
                "Post-imputation harmonisation did not produce a VCF for the "
                "pipeline input genome build %r." % output_build
            )
        args.vcf = output_vcf
        _validate_current_pipeline_vcf(args, ctx)
        ctx["imputation"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_ld_clump_runner(args, ctx):
    from postgwas.modules.ld_clumping.service import run_ld_clump_direct

    root = setup_subdir(args, "ld_clump")
    try:
        current = _validate_current_pipeline_vcf(args, ctx)
        outputs = run_ld_clump_direct(
            args,
            ctx=ctx,
            pipeline=True,
            cached_vcf=(current["indexed"] if current is not None else None),
        )
        ctx["ld_clump"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_finemap_runner(args, ctx):
    from postgwas.modules.fine_mapping.service import run_fine_mapping

    root = setup_subdir(args, "finemap")
    standard = ctx.get("ld_clump", {}).get("ld_clump_standard")
    locus_file = standard.get("ldpruned_sig_file") if standard else None
    if not locus_file:
        raise ValueError(
            "Fine-mapping cannot start because standard LD clumping produced no "
            "genomic risk loci. Check the lead P threshold, LD reference overlap, "
            "and standard clumping log."
        )
    args.locus_file = locus_file
    preflight = (
        ctx.validation("finemap")
        if hasattr(ctx, "validation")
        else None
    )
    resources = preflight.resources if preflight is not None else None
    method = getattr(args, "finemap_method", None)
    if method is None:
        from postgwas.config import load_module_configuration

        method = load_module_configuration(
            "fine_mapping", getattr(args, "run_config", None)
        ).engine
    upstream_results = {name: value for name, value in ctx.items() if name != "finemap"}
    try:
        if method == "susie":
            args.susie_input_file = ctx["formatter"]["susie"]["susie_input"]
            outputs = run_fine_mapping(
                args, resource_preflight=resources,
                upstream_results=upstream_results,
            )
            ctx["finemap"] = outputs
            return outputs
        elif method == "finemap":
            if "finemap" not in ctx["formatter"]:
                raise KeyError(
                    "Formatter output for 'finemap' is missing. Did the "
                    "formatter run correctly?"
                )

            args.finemap_in_files = ctx["formatter"]["finemap"]["finemap_input"]
            outputs = run_fine_mapping(
                args, resource_preflight=resources,
                upstream_results=upstream_results,
            )
            ctx["finemap"] = outputs
            return outputs
        else:
            raise ValueError(
                f"Unknown fine-mapping method '{method}'. Expected 'susie' "
                "or 'finemap'."
            )
    finally:
        args.output_directory = root


def run_magma_runner(args, ctx):
    from postgwas.modules.magma.service import run_magma_direct

    root = setup_subdir(args, "magma")

    magma_inputs = ctx["formatter"]["magma"]
    args.snp_location_file = magma_inputs["snp_loc_file"]
    args.p_value_file = magma_inputs["pval_file"]
    preflight = (
        ctx.validation("magma") if hasattr(ctx, "validation") else None
    )
    resources = preflight.resources if preflight is not None else None

    try:
        outputs = run_magma_direct(
            args,
            ctx,
            pipeline_resources=resources,
        )
        ctx["magma"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_magmacovar_runner(args, ctx):
    from postgwas.modules.magmacovar.service import run_magma_covar_direct

    root = setup_subdir(args, "magma_covar")

    magma = _validated_magma_gene_result(ctx, "MAGMA gene-property analysis")
    args.magma_gene_results_file = magma["magma_genes_raw"]

    try:
        outputs = run_magma_covar_direct(args, ctx)
        ctx["magma_covar"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_single_cell_runner(args, ctx):
    from postgwas.modules.single_cell.service import (
        resolve_single_cell_configuration,
        run_single_cell_direct,
    )

    root = setup_subdir(args, "single_cell")
    try:
        module = resolve_single_cell_configuration(args).modules.single_cell
        magma_tools = {
            tool for tool in module.tools
            if "magma" in SINGLE_CELL_TOOL_PIPELINE_DEPENDENCIES[tool]
        }
        magma = None
        if magma_tools:
            magma = _validated_magma_gene_result(ctx, "single-cell analysis")
        if "magma_celltype" in module.tools and magma is not None:
            args.magma_gene_results_file = magma["magma_genes_raw"]
        if "scdrs" in module.tools and magma is not None:
            primary = magma["primary_mapping"]
            actual_gene_id_type = magma["mapping_analyses"][primary].get(
                "gene_id_type"
            )
            declared_gene_id_type = (
                module.scdrs.magma_gene_set.source_identifier_type
            )
            if actual_gene_id_type != declared_gene_id_type:
                raise ValueError(
                    "scDRS declares MAGMA source identifiers as %r, but pipeline "
                    "mapping %r produced %r; update the scDRS identifier mapping "
                    "configuration before running"
                    % (declared_gene_id_type, primary, actual_gene_id_type)
                )
            args.scdrs_gene_set_source = "magma"
            args.scdrs_magma_gene_results_file = magma["magma_genes_out"]
        if "ldsc_celltype" in module.tools:
            formatter = ctx.get("formatter")
            ldsc = (
                formatter.get("ldsc") if isinstance(formatter, dict) else None
            )
            sumstats_file = (
                ldsc.get("ldsc_file") if isinstance(ldsc, dict) else None
            )
            if not sumstats_file:
                raise ValueError(
                    "LDSC cell-type analysis requires the formatter's validated "
                    "LDSC output"
                )
            args.ldsc_celltype_sumstats_source = "formatter"
            args.ldsc_celltype_sumstats_file = sumstats_file
        return run_single_cell_direct(args, ctx)
    finally:
        args.output_directory = root


def run_pops_runner(args, ctx):
    from postgwas.modules.pops.service import run_pops_direct

    root = setup_subdir(args, "pops")

    magma = _validated_magma_gene_result(ctx, "PoPS")
    if not magma.get("magma_genes_prefix"):
        raise ValueError(
            "PoPS requires the validated MAGMA gene-result prefix from the "
            "preceding pipeline step."
        )
    args.magma_association_prefix = magma["magma_genes_prefix"]
    args.magma_annotated_results_file = magma.get("magma_genes_annotated")

    try:
        outputs = run_pops_direct(args, ctx)
        return outputs
    finally:
        args.output_directory = root


def run_kpops_runner(args, ctx):
    from postgwas.modules.kpops.service import run_kpops_direct

    root = setup_subdir(args, "kpops")
    magma = _validated_magma_gene_result(ctx, "K-POPS")
    prefix = magma.get("magma_genes_prefix")
    if not prefix:
        raise ValueError("K-POPS requires the validated MAGMA gene-result prefix")
    args.magma_association_prefix = prefix
    preflight = ctx.validation("kpops") if hasattr(ctx, "validation") else None
    resources = preflight.resources if preflight is not None else None
    try:
        return run_kpops_direct(args, ctx, pipeline_resources=resources)
    finally:
        args.output_directory = root


def run_caldera_runner(args, ctx):
    from postgwas.modules.caldera.service import run_caldera_direct

    root = setup_subdir(args, "caldera")
    pops_file = ctx.get("pops_output")
    finemap = ctx.get("finemap", {})
    credible_sets_directory = (
        finemap.get("flames_input") if isinstance(finemap, dict) else None
    )
    if not pops_file:
        raise ValueError("CALDERA requires validated PoPS predictions")
    if not credible_sets_directory:
        raise ValueError("CALDERA requires validated fine-mapping credible sets")
    args.pops_file = pops_file
    args.finemap_credible_sets_directory = credible_sets_directory
    preflight = ctx.validation("caldera") if hasattr(ctx, "validation") else None
    resources = preflight.resources if preflight is not None else None
    try:
        return run_caldera_direct(args, ctx, pipeline_resources=resources)
    finally:
        args.output_directory = root


def run_flames_runner(args, ctx):
    finemap = ctx.get("finemap", {})
    credible_sets_directory = (
        finemap.get("flames_input") if isinstance(finemap, dict) else None
    )
    if not credible_sets_directory:
        status = (
            finemap.get("status", "not available")
            if isinstance(finemap, dict)
            else "not available"
        )
        raise ValueError(
            "FLAMES requires at least one retained fine-mapping credible set; "
            "fine-mapping status is %s." % status
        )

    from postgwas.modules.flames.service import run_flames_direct

    root = setup_subdir(args, "flames")

    args.credible_sets_directory = credible_sets_directory
    magma = _validated_magma_gene_result(ctx, "FLAMES")
    args.magma_gene_results_file = magma["magma_genes_out"]
    magmacovar = ctx.get("magma_covar")
    if not isinstance(magmacovar, dict) or not magmacovar.get("raw_results"):
        raise ValueError(
            "FLAMES requires validated raw MAGMAcovar results from the "
            "magmacovar pipeline step"
        )
    args.magma_covariate_results_file = magmacovar["raw_results"]
    args.pops_scores_file = ctx["pops_output"]

    try:
        outputs = run_flames_direct(args)
        ctx["flames"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_heritability_runner(args, ctx):
    from postgwas.modules.ldsc.service import (
        LDSCPipelineResources,
        run_ldsc_direct,
    )

    root = setup_subdir(args, "heritability")

    args.ldsc_input = ctx["formatter"]["ldsc"]["ldsc_file"]
    preflight = (
        ctx.validation("heritability") if hasattr(ctx, "validation") else None
    )
    resources = preflight.resources if preflight is not None else None
    if not isinstance(resources, LDSCPipelineResources):
        raise RuntimeError(
            "LDSC execution requires completed pipeline resource preflight evidence."
        )

    try:
        outputs = run_ldsc_direct(
            args,
            ctx,
            pipeline_resources=resources,
        )
        ctx["heritability"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_manhattan_runner(args, ctx):
    from postgwas.modules.manhattan.service import run_assoc_plot_direct

    root = setup_subdir(args, "manhattan")
    try:
        outputs = run_assoc_plot_direct(args)
        ctx["manhattan"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_qc_summary_runner(args, ctx):
    from postgwas.modules.qc_summary.service import run_qc_summary_direct

    root = setup_subdir(args, "qc_summary")
    try:
        outputs = run_qc_summary_direct(args)
        ctx["qc_summary"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_mixer_runner(args, ctx):
    from postgwas.modules.mixer.service import (
        MixerPipelineResources,
        run_mixer_direct,
    )

    root = setup_subdir(args, "mixer")
    try:
        formatter = ctx.get("formatter", {})
        mixer_input = formatter.get("mixer", {}).get("mixer_input")
        if not mixer_input:
            raise ValueError("MiXeR requires the 'mixer' artifact from formatter.")
        args.mixer_input_file = mixer_input
        preflight = (
            ctx.validation("mixer") if hasattr(ctx, "validation") else None
        )
        resources = preflight.resources if preflight is not None else None
        if not isinstance(resources, MixerPipelineResources):
            raise RuntimeError(
                "MiXeR execution requires completed pipeline resource "
                "preflight evidence."
            )
        outputs = run_mixer_direct(
            args,
            ctx,
            pipeline_resources=resources,
        )
        ctx["mixer"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_gcta_gene_runner(args, ctx):
    from postgwas.modules.gcta_gene.service import (
        GctaGenePipelineResources,
        run_gcta_gene_direct,
    )

    root = setup_subdir(args, "gcta_gene")
    preflight = (
        ctx.validation("gcta_gene") if hasattr(ctx, "validation") else None
    )
    resources = preflight.resources if preflight is not None else None
    if not isinstance(resources, GctaGenePipelineResources):
        raise RuntimeError(
            "GCTA gene execution requires completed pipeline resource "
            "preflight evidence."
        )
    try:
        return run_gcta_gene_direct(
            args,
            ctx,
            pipeline_resources=resources,
        )
    finally:
        args.output_directory = root


def run_gcta_cojo_runner(args, ctx):
    from postgwas.modules.gcta_cojo.service import (
        GctaCojoPipelineResources,
        run_gcta_cojo_direct,
    )

    root = setup_subdir(args, "gcta_cojo")
    formatter = ctx.get("formatter", {})
    gcta_input = formatter.get("gcta_gene", {}).get(
        "summary_statistics_input_file"
    )
    if not gcta_input:
        raise ValueError(
            "GCTA COJO requires the shared GCTA .ma artifact from formatter."
        )
    args.gcta_cojo_input_file = gcta_input
    preflight = (
        ctx.validation("gcta_cojo") if hasattr(ctx, "validation") else None
    )
    resources = preflight.resources if preflight is not None else None
    if not isinstance(resources, GctaCojoPipelineResources):
        raise RuntimeError(
            "GCTA-COJO execution requires completed pipeline resource preflight "
            "evidence."
        )
    try:
        return run_gcta_cojo_direct(
            args,
            ctx,
            pipeline_resources=resources,
        )
    finally:
        args.output_directory = root


# Runner registration lives exclusively in ``pipeline.registry``.  This module
# now contains implementations only; keeping a second name-to-runner mapping
# here was the source of CLI/planner/executor drift.
