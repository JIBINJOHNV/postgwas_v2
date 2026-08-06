"""
PostGWAS Pipeline Runners (v56, Dynamic Directories & Logic Fixes)

Changes (NO LOGIC CHANGE):
• setup_subdir: robust step prefix normalization (02 vs 2), always creates folders
• runners: never reference `outputs` in finally if step failed
• restore args.outdir safely even on error
"""

import os
import shutil
import sys
import traceback
import multiprocessing
from pathlib import Path

import pandas as pd  # Required for MAGMA batch calculation

# Workflow Imports
from postgwas.annot_ldblock.workflows import run_annot_ldblock
from postgwas.formatter.workflows import run_formatter_direct
from postgwas.sumstat_filter.workflows import run_sumstat_filter_direct
from postgwas.imputation.workflows import run_sumstat_imputation_direct
from postgwas.finemap.susie.workflows import run_parallel_susie
from postgwas.ld_clump.workflows import run_ld_clump_direct
from postgwas.gene_assoc.workflows import run_magma_direct
from postgwas.magmacovar.workflows import run_magma_covar_direct
from postgwas.pops.workflows import run_pops_direct
from postgwas.flames.workflows import run_flames_direct
from postgwas.h2_rg.workflows import run_ldsc_direct
from postgwas.manhattan.workflows import run_assoc_plot_direct
from postgwas.qc_summary.workflows import run_qc_summary_direct
from postgwas.finemap.finemap.workflows import run_finemap_pipeline

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

    original_outdir = args.outdir
    if original_outdir is None:
        raise ValueError("args.outdir is None")

    # Ensure root outdir exists (important in docker-mounted paths)
    os.makedirs(original_outdir, exist_ok=True)

    new_outdir = os.path.join(original_outdir, folder_name)
    os.makedirs(new_outdir, exist_ok=True)

    # Update args.outdir temporarily for the runner
    args.outdir = new_outdir
    return original_outdir


# ============================================================
# HELPER: Batch Calculator (MAGMA)
# ============================================================

def get_optimal_magma_batches(annot_file, user_requested_batches=None):
    """
    Calculates the optimal number of MAGMA batches based on Data and CPU limits.
    Prevents 'batch X is empty' errors on small datasets.
    """
    sys_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))

    try:
        df = pd.read_csv(annot_file, delim_whitespace=True, header=None, usecols=[1])
        data_chroms = df[1].nunique()
    except Exception:
        data_chroms = 22

    user_limit = int(user_requested_batches) if user_requested_batches else 999
    optimal = min(data_chroms, sys_cpus, user_limit)
    optimal = max(1, optimal)

    print(f"   ℹ️  [MAGMA] Auto-adjusted batches to: {optimal} (Chromosomes: {data_chroms})")
    return optimal


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


# ============================================================
# RUNNERS
# ============================================================

def run_annot_ldblock_runner(args, ctx):
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
        # always restore outdir
        args.outdir = root


def run_sumstat_filter_runner(args, ctx):
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
        args.outdir = root


def run_formatter_runner(args, ctx):
    root = setup_subdir(args, "formatter")

    # 1. Determine Formats
    if not hasattr(args, "format") or args.format is None:
        active_formats = set()
    else:
        active_formats = set(args.format)

    active_modules = set(args.modules) if args.modules else set()

    if getattr(args, "apply_imputation", False) or "imputation" in active_modules:
        active_formats.add("ldpred")
    if "heritability" in active_modules:
        active_formats.add("ldsc")
    if any(m in active_modules for m in ["magma", "magmacovar", "pops", "flames"]):
        active_formats.add("magma")
    if "finemap" in active_modules or "flames" in active_modules:
        active_formats.add("finemap_susie")
        active_formats.add("finemap_finemap")

    args.format = list(active_formats)

    # 2. Resolve Binaries
    required = ["bcftools", "tabix", "bgzip"]
    if any(fmt in args.format for fmt in ["ldpred", "magma", "finemap"]):
        required.append("plink")
    check_and_resolve_binaries(args, required)

    # 3. Validate Input VCF exists
    if args.vcf and not os.path.exists(args.vcf):
        print(f"\n❌ File Not Found: input VCF does not exist:\n   {args.vcf}")
        raise FileNotFoundError(f"Input VCF missing: {args.vcf}")

    try:
        outputs = run_formatter_direct(args)
        ctx["formatter"] = outputs
        return outputs
    except Exception:
        print("\n\n🔥 [CRITICAL FORMATTER FAILURE] 🔥")
        traceback.print_exc()
        print("---------------------------------------------------")
        raise
    finally:
        args.outdir = root


def run_imputation_runner(args, ctx):
    root = setup_subdir(args, "imputation")

    if "formatter" not in ctx or "ldpred" not in ctx["formatter"] or not ctx["formatter"]["ldpred"]:
        raise ValueError("Imputation requires 'ldpred' format from formatter.")

    ldpred_folder = ctx["formatter"]["ldpred"].get("ldpred_folder")
    if not ldpred_folder:
        raise ValueError("Formatter did not return a valid 'ldpred_folder' path.")

    args.predld_input_dir = ldpred_folder

    try:
        outputs = run_sumstat_imputation_direct(args)
        # keep behavior: set vcf to GRCh37 output (as you had)
        if isinstance(outputs, dict) and "GRCh37" in outputs:
            args.vcf = outputs["GRCh37"]
        ctx["imputation"] = outputs
        return outputs
    finally:
        args.outdir = root


def run_ld_clump_runner(args, ctx):
    root = setup_subdir(args, "ld_clump")
    args.ld_mode = "by_regions"
    try:
        outputs = run_ld_clump_direct(args)
        ctx["ld_clump"] = outputs
        print("ld clumping completed")
        print(ctx["ld_clump"])
        return outputs
    finally:
        args.outdir = root


def run_finemap_runner(args, ctx):
    root = setup_subdir(args, "finemap")
    print("Finemap analysis started")
    #args.locus_file = ctx["ld_clump"]["ldpruned_sig_file"]
    args.locus_file = ctx["ld_clump"]["ld_clump_standard"]["ldpruned_sig_file"]
    print(args.locus_file )
    # --- FIX: Normalize the argument (Handle List vs String) ---
    method = args.finemap_method
    if isinstance(method, list):
        # Create a clean string like "susie, finemap"
        methods_str = ", ".join(method)
        print(f"\t\t\t ℹ️  Multiple methods provided ({methods_str}). Defaulting to the first one: '{method[0]}'")
        method = method[0]
    try:
        # Check against the normalized 'method' variable
        if method == "susie":
            args.susie_input_file = ctx["formatter"]["finemap_susie"]["susie_input"]
            outputs = run_parallel_susie(args)
            ctx["finemap"] = outputs
            return outputs
        elif method == "finemap":
            print("\t\t\t🔹 Using FINEMAP method for finemapping")
            # Ensure formatter output exists
            if "finemap_finemap" not in ctx["formatter"]:
                 raise KeyError("Formatter output for 'finemap_finemap' is missing. Did the formatter run correctly?")

            args.finemap_in_files = ctx["formatter"]["finemap_finemap"]["finemap_input"]
            
            # Use the correct pipeline function you identified
            # Ensure 'run_finemap_pipeline' is imported at the top of runners.py!
            outputs = run_finemap_pipeline(args)
            
            ctx["finemap"] = outputs
            print(f"   ✅ Output: {outputs}")
            return outputs

        # --- SAFETY NET ---
        else:
            raise ValueError(f"CRITICAL ERROR: Unknown finemap method '{method}'. Options are 'susie' or 'finemap'.")
    finally:
        args.outdir = root


def run_magma_runner(args, ctx):
    root = setup_subdir(args, "magma")

    magma_inputs = ctx["formatter"]["magma"]
    args.snp_loc_file = magma_inputs["snp_loc_file"]
    args.pval_file = magma_inputs["pval_file"]

    try:
        outputs = run_magma_direct(args, ctx)
        ctx["magma"] = outputs
        return outputs
    finally:
        args.outdir = root


def run_magmacovar_runner(args, ctx):
    root = setup_subdir(args, "magma_covar")

    # FIX: key is "magma"
    args.magama_gene_assoc_raw = ctx["magma"]["magma_genes_raw"]

    try:
        outputs = run_magma_covar_direct(args)
        ctx["magma_covar"] = outputs
        return outputs
    finally:
        args.outdir = root


def run_pops_runner(args, ctx):
    root = setup_subdir(args, "pops")

    if "magma" in ctx and isinstance(ctx["magma"], dict) and "magma_genes_prefix" in ctx["magma"]:
        args.magma_assoc_prefix = ctx["magma"]["magma_genes_prefix"]
    else:
        print("⚠️ Warning: MAGMA prefix not found in context. PoPS might fail.")

    args.pops_verbose = True

    try:
        outputs = run_pops_direct(args, ctx)
        # keep your behavior: store pops_file into ctx["pops_output"]
        if "pops_output" not in ctx and isinstance(outputs, dict):
            ctx["pops_output"] = outputs.get("pops_file")
        else:
            # if already set elsewhere, leave it
            pass
        return outputs
    finally:
        args.outdir = root


def run_flames_runner(args, ctx):
    root = setup_subdir(args, "flames")

    args.finemap_cred_dir = ctx["finemap"]["flames_input"]
    args.magma_genes_out = ctx["magma"]["magma_genes_out"]
    args.magma_tissue_covar_results = ctx["magma_covar"]
    args.pops_score_file = ctx["pops_output"]

    try:
        outputs = run_flames_direct(args)
        ctx["flames"] = outputs
        return outputs
    finally:
        args.outdir = root

 
# def run_heritability_runner(args, ctx):
#     root = setup_subdir(args, "heritability")

#     # (typo preserved from your code: ldsc_inut)
#     args.ldsc_inut = ctx["formatter"]["ldsc"]["ldsc_file"] 
#     if args.sample_prev is None and args.population_prev is not None:
#         print("Since Population prevalence provided and sample prevalence is not provided, using inferred sample prevalence from summary statistics")
#         print(f"Sample Prevalence = median(N_CAS) / (median(N_CAS) + median(N_CON))")
#         print(f"Sample Prevalence = {ctx['formatter']['ldsc']['sample_prev']}")
#         args.sample_prev = ctx["formatter"]["ldsc"]["sample_prev"]
#     try:
#         outputs = run_ldsc_direct(args)
#         ctx["heritability"] = outputs
#         return outputs
#     finally:
#         args.outdir = root

def run_heritability_runner(args, ctx):
    root = setup_subdir(args, "heritability")

    args.ldsc_inut = ctx["formatter"]["ldsc"]["ldsc_file"]

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
        args.outdir = root
        

def run_manhattan_runner(args, ctx):
    root = setup_subdir(args, "manhattan")
    try:
        outputs = run_assoc_plot_direct(args)
        ctx["manhattan"] = outputs
        return outputs
    finally:
        args.outdir = root


def run_qc_summary_runner(args, ctx):
    root = setup_subdir(args, "qc_summary")
    try:
        outputs = run_qc_summary_direct(args)
        ctx["qc_summary"] = outputs
        return outputs
    finally:
        args.outdir = root


# ============================================================
# PIPELINE REGISTRY
# ============================================================

RUNNERS = {
    "annot_ldblock": run_annot_ldblock_runner,
    "sumstat_filter": run_sumstat_filter_runner,
    "post_imputation_filter": run_sumstat_filter_runner,
    "formatter": run_formatter_runner,
    "imputation": run_imputation_runner,
    "finemap": run_finemap_runner,
    "ld_clump": run_ld_clump_runner,
    "magma": run_magma_runner,
    "magmacovar": run_magmacovar_runner,
    "pops": run_pops_runner,
    "flames": run_flames_runner,
    "heritability": run_heritability_runner,
    "manhattan": run_manhattan_runner,
    "qc_summary": run_qc_summary_runner,
}
