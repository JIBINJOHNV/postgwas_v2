#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v52, Error Deduplication Fix)

Changes:
• FIX: Groups missing arguments to prevent duplicate error messages.
  (e.g., reports '--vcf' once, listing all modules that need it).
"""

import argparse
import sys
from rich_argparse import RichHelpFormatter
from rich.console import Console
from rich.table import Table

# Initialize Rich Console Globally
console = Console()

# =====================================================================
# IMPORT ORCHESTRATOR
# =====================================================================
from postgwas.pipeline.execution import execute_pipeline 

# =====================================================================
# CUSTOM PARSER
# =====================================================================

class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)

# =====================================================================
# PARSER BUILDERS
# =====================================================================

from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_genomeversion_parser,
    get_common_out_parser,
    get_common_sumstat_filter_parser,
    get_formatter_parser,
    get_plink_binary_parser,
    get_common_magma_assoc_parser,
    get_common_magma_covar_parser,
    get_common_pops_parser,
    get_flames_common_parser,
    get_common_imputation_parser,
    get_ldsc_common_parser,
    get_annot_ldblock_parser,
    get_assoc_plot_parser,
    get_magma_binary_parser,
    get_bcftools_binary_parser,
    get_plink_binary_parser,
    sumstat_summary_arg_parser,
    get_ld_clump_parser
)
from postgwas.clis.finemap_cli import (
    get_finemap_common_parser,
    get_common_susie_arguments,
    get_common_finemap_finemap_arguments
)

def get_qc_parser():
    parser = argparse.ArgumentParser(add_help=False)
    g = parser.add_argument_group("QC Summary Arguments")
    g.add_argument("--qc_output_prefix")
    return parser

# =====================================================================
# MODULE CONFIGURATION
# =====================================================================

MODULE_DESCRIPTIONS = {
    "harmonisation": ["Standardises summary statistics and converts them to VCF format."],
    "annot_ldblock": ["Annotates variants with LD blocks."],
    "formatter": ["Converts harmonised VCF to tool-specific formats."],
    "sumstat_filter": ["Filters summary statistics: INFO, MAF, QC."],
    "imputation": ["Performs summary statistics imputation using PRED-LD."],
    "finemap": ["Runs SuSiE fine-mapping."],
    "ld_clump": ["Runs PLINK LD-clumping."],
    "magma": ["Runs MAGMA gene-level association analysis."],
    "magmacovar": ["Runs MAGMA gene-property analysis."],
    "pops": ["Computes PoPS gene prioritisation scores."],
    "flames": ["Runs FLAMES integrative effector-gene scoring."],
    "heritability": ["Runs LDSC heritability."],
    "manhattan": ["Generates Manhattan & QQ plots."],
    "qc_summary": ["Generates QC summaries."],
}

#"post_imputation_filter": ["SECONDARY FILTER: Runs automatically after Harmonisation."],

MODULE_DEPENDENCIES = {
    "annot_ldblock": [],
    "harmonisation": [],
    "sumstat_filter": [],
    "formatter":     [],
    "magma":         [],
    "heritability":  [],
    "manhattan":     [],
    "qc_summary":    [],
    "imputation": ["formatter"],
    "post_imputation_filter": [],
    "ld_clump": ["annot_ldblock"], 
    "finemap": ["ld_clump"],
    "magmacovar": ["magma"],
    "pops": ["magma"],
    "flames": ["annot_ldblock", "finemap", "ld_clump", "magma", "magmacovar", "pops"],
}

MODULE_PARSERS = {
    "harmonisation": [],
    "annot_ldblock": [get_defaultresourse_parser, get_inputvcf_parser, get_genomeversion_parser, get_common_out_parser, get_annot_ldblock_parser],
    "formatter": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_formatter_parser, get_bcftools_binary_parser],
    "finemap": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_finemap_common_parser, get_common_susie_arguments,get_common_finemap_finemap_arguments, get_bcftools_binary_parser],
    "ld_clump": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_plink_binary_parser, get_bcftools_binary_parser,get_ld_clump_parser],
    "magma": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_common_magma_assoc_parser, get_magma_binary_parser, get_bcftools_binary_parser],
    "magmacovar": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_common_magma_assoc_parser, get_common_magma_covar_parser, get_magma_binary_parser, get_bcftools_binary_parser],
    "pops": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_common_pops_parser, get_common_magma_assoc_parser, get_common_magma_covar_parser, get_bcftools_binary_parser, get_magma_binary_parser],
    "flames": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_common_magma_assoc_parser, get_common_magma_covar_parser, get_common_pops_parser, get_flames_common_parser, get_bcftools_binary_parser, get_magma_binary_parser],
    "sumstat_filter": [get_defaultresourse_parser, get_inputvcf_parser, get_genomeversion_parser, get_common_out_parser, get_common_sumstat_filter_parser, get_bcftools_binary_parser],
    "post_imputation_filter": [get_defaultresourse_parser, get_inputvcf_parser, get_genomeversion_parser, get_common_out_parser, get_common_sumstat_filter_parser, get_bcftools_binary_parser],
    "imputation": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_common_imputation_parser, get_bcftools_binary_parser],
    "heritability": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_ldsc_common_parser],
    "manhattan": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_assoc_plot_parser],
    "qc_summary": [get_defaultresourse_parser, get_inputvcf_parser, get_common_out_parser, get_qc_parser,sumstat_summary_arg_parser, get_bcftools_binary_parser],
}

# =====================================================================
# HELPERS
# =====================================================================

def resolve_dependencies_unique(mods):
    """Returns unique list of modules for Argument Parsing"""
    resolved, visited = [], set()
    def dfs(x):
        if x in visited: return
        visited.add(x)
        for dep in MODULE_DEPENDENCIES.get(x, []):
            dfs(dep)
        resolved.append(x)
    for m in mods: dfs(m)
    return resolved

def print_full_pipeline_help(modules, parser):
    console.print("\n❗ Required arguments, resources, or dependencies are missing.")
    console.print("   Displaying full pipeline help:\n")
    console.print("📘 Detailed Pipeline Execution Plan:\n")
    for idx, m in enumerate(modules, 1):
        if m in MODULE_DESCRIPTIONS:
            console.print(f" {idx}) [cyan]{m}[/cyan]")
            for line in MODULE_DESCRIPTIONS[m]:
                console.print(f"      • {line}")
        print("")
    console.print("👇 Below are the arguments relevant to this pipeline:\n")
    parser.print_help()

# =====================================================================
# MAIN
# =====================================================================

def main():

    # 1. Minimal Parse
    mini = argparse.ArgumentParser(add_help=False)
    mini.add_argument("--modules", nargs="*")
    mini.add_argument("--apply-filter", action="store_true")
    mini.add_argument("--apply-imputation", action="store_true")
    mini.add_argument("--apply-manhattan", action="store_true")
    mini.add_argument("--heritability", action="store_true")
    mini.add_argument("-h", "--help", action="store_true")

    a1, _ = mini.parse_known_args()

    # 2. Show Table if Empty
    if not a1.modules and not (a1.apply_filter or a1.apply_imputation or a1.heritability):
        console.print("\nEach module operates directly on your harmonised GWAS VCF.\n")
        console.print("[bold cyan]Usage Example[/bold cyan]\n  postgwas pipeline --modules finemap --apply-filter\n")
        t = Table(title="PostGWAS Pipeline — Available Modules", header_style="bold cyan", show_lines=True)
        t.add_column("Module", style="magenta", no_wrap=True)
        t.add_column("Description", style="white")
        for k in MODULE_DESCRIPTIONS:
            t.add_row(k, MODULE_DESCRIPTIONS[k][0])
        console.print(t)
        sys.exit(0)

    # 3. Resolve Modules (Stage 1: Identify Active Components)
    raw_modules = list(a1.modules) if a1.modules else []
    
    # 3a. Add Flags to list
    if a1.apply_filter and "sumstat_filter" not in raw_modules:
        raw_modules.insert(0, "sumstat_filter")
    if a1.apply_imputation and "imputation" not in raw_modules:
        raw_modules.append("imputation")
    if a1.heritability and "heritability" not in raw_modules:
        raw_modules.append("heritability")
    if a1.apply_manhattan and "manhattan" not in raw_modules:
        raw_modules.append("manhattan")

    # 3b. Calculate ALL dependencies immediately
    all_active_modules = set(resolve_dependencies_unique(raw_modules))
    is_imputing = "imputation" in all_active_modules
    
    # 3c. Construct Explicit Execution Order
    execution_modules = []
    
    # --- PHASE 1: PRE-IMPUTATION FILTER ---
    if "sumstat_filter" in all_active_modules:
        execution_modules.append("sumstat_filter")
        
    # --- PHASE 2: IMPUTATION BLOCK ---
    if is_imputing:
        # Pre-imp formatting
        if "formatter" not in execution_modules: 
            execution_modules.append("formatter")
        
        execution_modules.append("imputation")
        
        # Post-imp filtering
        if "sumstat_filter" in all_active_modules or "post_imputation_filter" in all_active_modules:
            execution_modules.append("post_imputation_filter")

    # --- PHASE 3: ANALYSIS PREP (Annotation & Formatting) ---
    if "annot_ldblock" in all_active_modules:
        execution_modules.append("annot_ldblock")

    downstream_tools = ["ld_clump", "finemap", "magma", "magmacovar", "pops", "flames", "heritability"]
    
    # LOGIC FIX: Only run 2nd formatter if downstream analysis exists OR explicit user request
    has_downstream = any(m in all_active_modules for m in downstream_tools)
    user_explicitly_asked_format = "formatter" in raw_modules

    if has_downstream or (is_imputing and user_explicitly_asked_format):
        execution_modules.append("formatter")

    # --- PHASE 4: DOWNSTREAM ANALYSIS ---
    for tool in downstream_tools:
        if tool in all_active_modules:
            execution_modules.append(tool)
            
    if "manhattan" in all_active_modules:
        execution_modules.append("manhattan")
    if "qc_summary" in all_active_modules:
        execution_modules.append("qc_summary")

    # 3d. Build Unique List for Argument Parser
    parser_modules = resolve_dependencies_unique(execution_modules)

    # 4. Build Full Parser
    parent_parsers = []
    seen = set()
    for m in parser_modules:
        if m in MODULE_PARSERS:
            for fn in MODULE_PARSERS[m]:
                if fn not in seen:
                    parent_parsers.append(fn())
                    seen.add(fn)

    custom_usage = "postgwas pipeline [options] --modules " + " ".join(parser_modules)
    
    parser = HelpOnErrorArgumentParser(
        prog="postgwas pipeline",
        formatter_class=RichHelpFormatter,
        parents=parent_parsers,
        conflict_handler='resolve',
        usage=custom_usage
    )

    # -------------------------------------------------------------
    # Context-Sensitive Arguments
    # -------------------------------------------------------------
    
    is_formatter_explicit = "formatter" in raw_modules
    
    for action in parser._actions:
        if action.dest == "format":
            if not is_formatter_explicit:
                action.required = False
                action.help = argparse.SUPPRESS 

    # Re-add CLI Flags
    parser.add_argument("--modules", nargs="*")
    parser.add_argument("--apply-filter", action="store_true")
    parser.add_argument("--apply-imputation", action="store_true")
    parser.add_argument("--apply-manhattan", action="store_true")
    parser.add_argument("--heritability", action="store_true")

    if a1.help:
        print_full_pipeline_help(execution_modules, parser)
        sys.exit(0)

    # 5. Parse & Validate
    try:
        args = parser.parse_args()
    except ValueError as e:
        console.print(f"\n❌ [bold red]Argument Parsing Error:[/bold red] {e}")
        print_full_pipeline_help(execution_modules, parser)
        sys.exit(2)

    REQUIRED_ARGS = {
        "annot_ldblock": ["vcf", "genome-version","sample_id","outdir","ld-region-dir"],
        "sumstat_filter": ["vcf","sample_id","outdir"],
        "post_imputation_filter": ["vcf","sample_id","outdir"], 
        "imputation": ["vcf","sample_id","outdir","ref_ld","gwas2vcf_resource"],
        "formatter": ["vcf","sample_id","outdir"],
        "finemap": ["vcf","sample_id","outdir","ld-region-dir","finemap_ld_ref"],
        "ld_clump": ["vcf", "genome-version","sample_id","outdir","ld-region-dir",'ld_clump_population'],
        "magma": ["vcf","sample_id","outdir","ld_ref","gene_loc_file"],
        "magmacovar": ["vcf","sample_id","outdir","ld_ref","gene_loc_file","covariates"],
        "pops": ["vcf","sample_id","outdir","ld_ref","feature_mat_prefix","pops_gene_loc_file"],
        "flames": ["vcf","sample_id","outdir","ld_ref","covariates","feature_mat_prefix","pops_gene_loc_file","flames_annot_dir"],
        "heritability": ["vcf","sample_id","outdir","merge-alleles","ref-ld-chr","w-ld-chr"],
        "manhattan": ["vcf","sample_id","outdir"],
        "qc_summary": ["vcf","sample_id","outdir"],
    }

    # -----------------------------------------------------------
    # NEW LOGIC: Deduplicate Error Messages
    # -----------------------------------------------------------
    missing_args_map = {} # arg_name -> list of modules

    check_list = parser_modules.copy()
    if not is_formatter_explicit and "formatter" in check_list:
        check_list.remove("formatter")

    for m in check_list:
        if m in REQUIRED_ARGS:
            for arg_name in REQUIRED_ARGS[m]:
                attr_name = arg_name.replace("-", "_")
                val = getattr(args, attr_name, None)
                if not val:
                    if arg_name not in missing_args_map:
                        missing_args_map[arg_name] = []
                    missing_args_map[arg_name].append(m)

    if missing_args_map:
        console.print("\n❌ [bold red]Missing Required Arguments:[/bold red]")
        for arg, modules in missing_args_map.items():
            mod_list = ", ".join([f"[cyan]{m}[/cyan]" for m in modules])
            console.print(f"   • Argument [bold red]--{arg}[/bold red] is required by module(s): {mod_list}")

        print_full_pipeline_help(execution_modules, parser)
        sys.exit(2)

    # 6. Execute
    try:
        args.modules = parser_modules 
        execute_pipeline(args, execution_modules)
        
    except (ValueError, OSError, AttributeError, TypeError, SystemExit) as e:
        if isinstance(e, SystemExit):
            if e.code != 0:
                print_full_pipeline_help(execution_modules, parser)
                sys.exit(e.code)
            sys.exit(0)

        err_msg = str(e).lower()
        if "required" in err_msg or "executable not found" in err_msg:
             console.print(f"\n❌ [bold red]Execution Error:[/bold red] {e}") 
             print_full_pipeline_help(execution_modules, parser)
             sys.exit(2)
        raise

    console.print("\n🎉 Pipeline complete.\n")

if __name__ == "__main__":
    main()