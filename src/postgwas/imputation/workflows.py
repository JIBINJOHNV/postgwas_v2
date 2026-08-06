import argparse
import sys
import os
from pathlib import Path
from rich_argparse import RichHelpFormatter



from postgwas.harmonisation.cli import run_harmonisation


# =========================================================
# DISPATCHERS
# =========================================================
def run_sumstat_imputation_direct(args,ctx=None):
    """
    DIRECT mode:
      - Run ONLY the chosen imputation tool.
      - No mandatory post-processing.
    """
    tool = args.imputation_tool
    if tool == "pred_ld":
        # Lazy import: only bring in PRED-LD when needed
        from postgwas.imputation.pred_ld.pred_ld_runner import (
            run_pred_ld_parallel,process_pred_ld_results_all_parallel
        )
        run_pred_ld_parallel(
            predld_input_dir=args.predld_input_dir,
            pred_ld_ref=args.ref_ld,
            output_folder=args.outdir,
            output_prefix=args.sample_id,
            r2threshold=args.r2threshold,
            maf=args.maf,
            population=args.population,
            ref=args.ref,
            threads=args.nthreads
        )
        print("         ✅ Sumstat Imputation using PRED-LD software completed.")
        print(" ")
    else:
        raise ValueError(f"Imputation tool not yet implemented for DIRECT mode: {tool}")
    combined_df, corr_df,config_path=process_pred_ld_results_all_parallel(
        folder_path=args.outdir,
        output_path=args.outdir,
        gwas2vcf_resource_folder=args.gwas2vcf_resource,
        output_prefix=args.sample_id,
        corr_method= "pearson",
        threads=args.nthreads,
    )
    args.config=config_path
    args.defaults=args.gwas2vcf_default_config
    ## perform harmonisation
    outputs = run_harmonisation(args)
    if ctx is not None:
        ctx["imputation"] = outputs
    return outputs
