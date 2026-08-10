"""
PostGWAS Pipeline Runners (v56, Dynamic Directories & Logic Fixes)

Changes (NO LOGIC CHANGE):
• setup_subdir: robust step prefix normalization (02 vs 2), always creates folders
• runners: never reference `outputs` in finally if step failed
• restore args.output_directory safely even on error
"""

import os
import shutil

from postgwas.config import load_module_configuration
from postgwas.modules.formatting.contracts import required_formats
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)

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


# ============================================================
# HELPER: Dependency Debugger
# ============================================================

def check_and_resolve_binaries(args, required_tools):
    def inject_to_path(binary_path):
        directory = os.path.dirname(binary_path)
        if directory and directory not in os.environ.get("PATH", ""):
            os.environ["PATH"] = directory + os.pathsep + os.environ.get("PATH", "")

    if "bcftools" in required_tools or hasattr(args, "bcftools"):
        cmd = getattr(args, "bcftools", "bcftools")
        resolved = shutil.which(cmd)
        if resolved:
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


# ============================================================
# HELPER: VCF Extractor
# ============================================================

def get_vcf_from_context(ctx):
    if "post_imputation_filter" in ctx:
        return ctx["post_imputation_filter"]["filtered_vcf"]

    if "imputation" in ctx:
        data = ctx["imputation"]
        if isinstance(data, dict):
            return data.get("GRCh38", data.get("GRCh37"))
        return data

    if "sumstat_filter" in ctx:
        return ctx["sumstat_filter"]["filtered_vcf"]

    if "annot_ldblock" in ctx:
        return ctx["annot_ldblock"]["annotated_vcf"]

    return None


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
        return outputs
    finally:
        args.output_directory = root


def run_formatter_runner(args, ctx):
    from postgwas.modules.formatting.service import run_formatter_direct

    active_modules = list(getattr(args, "modules", None) or ())
    if getattr(args, "apply_imputation", False):
        active_modules.append("imputation")
    formatting_config = load_module_configuration(
        "formatting", getattr(args, "run_config", None),
    )
    args.format = required_formats(
        formatting_config,
        active_modules,
        getattr(args, "format", None) or (),
    )
    requirements = []
    if any(
        "magma" in formatting_config.module_formats.get(module, ())
        for module in active_modules
    ):
        magma_config = load_module_configuration(
            "magma", getattr(args, "run_config", None),
        )
        reference_prefix = (
            getattr(args, "magma_ld_reference", None)
            or magma_config.input.ld_reference_prefix
        )
        if not reference_prefix:
            raise ValueError(
                "MAGMA requires --magma-ld-reference before formatting so the "
                "BIM identifier convention can be selected."
            )
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
        module_config = load_module_configuration(
            module_name, getattr(args, "run_config", None),
        )
        argument_name = (
            "gcta_reference_prefix"
            if module_name == "gcta_gene"
            else "cojo_reference_prefix"
        )
        reference_prefix = (
            getattr(args, argument_name, None) or module_config.reference.prefix
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
        configure_reference_variant_identifiers(
            args, formatting_config, requirements,
        )

    # 2. Resolve Binaries
    check_and_resolve_binaries(args, ["bcftools"])

    # 3. Validate Input VCF exists
    if args.vcf and not os.path.exists(args.vcf):
        raise FileNotFoundError(f"Input VCF missing: {args.vcf}")

    root = setup_subdir(args, "formatter")
    try:
        return run_formatter_direct(args, ctx)
    finally:
        args.output_directory = root


def run_imputation_runner(args, ctx):
    from postgwas.modules.imputation.service import run_sumstat_imputation_direct

    root = setup_subdir(args, "imputation")

    if "formatter" not in ctx or "pred_ld" not in ctx["formatter"] or not ctx["formatter"]["pred_ld"]:
        raise ValueError("Imputation requires the 'pred_ld' format from formatter.")

    pred_ld_folder = ctx["formatter"]["pred_ld"].get("pred_ld_folder")
    if not pred_ld_folder:
        raise ValueError("Formatter did not return a valid 'pred_ld_folder' path.")

    args.pred_ld_input_directory = pred_ld_folder

    try:
        outputs = run_sumstat_imputation_direct(args)
        # keep behavior: set vcf to GRCh37 output (as you had)
        if isinstance(outputs, dict) and "GRCh37" in outputs:
            args.vcf = outputs["GRCh37"]
        ctx["imputation"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_ld_clump_runner(args, ctx):
    from postgwas.modules.ld_clumping.service import run_ld_clump_direct

    root = setup_subdir(args, "ld_clump")
    args.ld_mode = "by_regions"
    try:
        outputs = run_ld_clump_direct(args)
        ctx["ld_clump"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_finemap_runner(args, ctx):
    from postgwas.modules.fine_mapping.service import run_fine_mapping

    root = setup_subdir(args, "finemap")
    args.locus_file = ctx["ld_clump"]["ld_clump_standard"]["ldpruned_sig_file"]
    method = getattr(args, "finemap_method", None)
    if method is None:
        from postgwas.config import load_module_configuration

        method = load_module_configuration(
            "fine_mapping", getattr(args, "run_config", None)
        ).engine
    try:
        if method == "susie":
            args.susie_input_file = ctx["formatter"]["susie"]["susie_input"]
            outputs = run_fine_mapping(args)
            ctx["finemap"] = outputs
            return outputs
        elif method == "finemap":
            if "finemap" not in ctx["formatter"]:
                raise KeyError(
                    "Formatter output for 'finemap' is missing. Did the "
                    "formatter run correctly?"
                )

            args.finemap_in_files = ctx["formatter"]["finemap"]["finemap_input"]
            outputs = run_fine_mapping(args)
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

    try:
        outputs = run_magma_direct(args, ctx)
        ctx["magma"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_magmacovar_runner(args, ctx):
    from postgwas.modules.magmacovar.service import run_magma_covar_direct

    root = setup_subdir(args, "magma_covar")

    # FIX: key is "magma"
    magma = _validated_magma_gene_result(ctx, "MAGMA gene-property analysis")
    args.magma_gene_results_file = magma["magma_genes_raw"]

    try:
        outputs = run_magma_covar_direct(args)
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
        magma = _validated_magma_gene_result(ctx, "single-cell analysis")
        module = resolve_single_cell_configuration(args).modules.single_cell
        if "magma_celltype" in module.tools:
            args.magma_gene_results_file = magma["magma_genes_raw"]
        if "scdrs" in module.tools:
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
    try:
        return run_kpops_direct(args, ctx)
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
    try:
        return run_caldera_direct(args, ctx)
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
    args.magma_covariate_results_file = ctx["magma_covar"]
    args.pops_scores_file = ctx["pops_output"]

    try:
        outputs = run_flames_direct(args)
        ctx["flames"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_heritability_runner(args, ctx):
    from postgwas.modules.ldsc.service import run_ldsc_direct

    root = setup_subdir(args, "heritability")

    args.ldsc_input = ctx["formatter"]["ldsc"]["ldsc_file"]

    if args.samp_prev is None and args.pop_prev is not None:
        print("Since population prevalence is provided and sample prevalence is not provided, using inferred sample prevalence from summary statistics")

        print("Sample Prevalence = median(N_CAS) / (median(N_CAS) + median(N_CON))")

        sprev = ctx["formatter"]["ldsc"]["sample_prev"]

        print(f"Sample Prevalence = {sprev:.4f}")

        args.samp_prev = sprev

    try:
        outputs = run_ldsc_direct(args)
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
    from postgwas.modules.mixer.service import run_mixer_direct

    root = setup_subdir(args, "mixer")
    try:
        formatter = ctx.get("formatter", {})
        mixer_input = formatter.get("mixer", {}).get("mixer_input")
        if not mixer_input:
            raise ValueError("MiXeR requires the 'mixer' artifact from formatter.")
        args.mixer_input_file = mixer_input
        outputs = run_mixer_direct(args, ctx)
        ctx["mixer"] = outputs
        return outputs
    finally:
        args.output_directory = root


def run_gcta_gene_runner(args, ctx):
    from postgwas.modules.gcta_gene.service import run_gcta_gene_direct

    root = setup_subdir(args, "gcta_gene")
    try:
        return run_gcta_gene_direct(args, ctx)
    finally:
        args.output_directory = root


def run_gcta_cojo_runner(args, ctx):
    from postgwas.modules.gcta_cojo.service import run_gcta_cojo_direct

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
    try:
        return run_gcta_cojo_direct(args, ctx)
    finally:
        args.output_directory = root


# Runner registration lives exclusively in ``pipeline.registry``.  This module
# now contains implementations only; keeping a second name-to-runner mapping
# here was the source of CLI/planner/executor drift.
