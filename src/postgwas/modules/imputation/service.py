import argparse
import sys
import os
from pathlib import Path
from rich_argparse import RichHelpFormatter



from postgwas.modules.harmonisation.cli import run_harmonisation
from postgwas.cli.compute import resolve_compute_args


# =========================================================
# DISPATCHERS
# =========================================================
def run_sumstat_imputation_direct(args,ctx=None):
    """
    DIRECT mode:
      - Run ONLY the chosen imputation tool.
      - No mandatory post-processing.
    """
    resolve_compute_args(args)
    tool = args.imputation_engine
    if tool == "pred_ld":
        # Lazy import: only bring in PRED-LD when needed
        from postgwas.modules.imputation.engines.pred_ld.pred_ld_runner import (
            run_pred_ld_parallel,process_pred_ld_results_all_parallel
        )
        run_pred_ld_parallel(
            predld_input_dir=args.pred_ld_input_directory,
            pred_ld_ref=args.imputation_ld_reference,
            output_folder=args.output_directory,
            output_prefix=args.dataset_id,
            r2threshold=args.imputation_r2_threshold,
            maf=args.imputation_minimum_maf,
            population=args.population,
            ref=args.ref,
            threads=args.threads
        )
        print("         ✅ Sumstat Imputation using PRED-LD software completed.")
        print(" ")
    else:
        raise ValueError(f"Imputation tool not yet implemented for DIRECT mode: {tool}")
    combined_df, corr_df, sample_sheet_path = process_pred_ld_results_all_parallel(
        folder_path=args.output_directory,
        output_path=args.output_directory,
        output_prefix=args.dataset_id,
        corr_method=args.corr_method,
        threads=args.threads,
    )
    args.sample_sheet = sample_sheet_path
    args.output_directory = str(
        Path(args.output_directory).parent / "imputed_harmonised"
    )
    outputs = run_harmonisation(args)
    if ctx is not None:
        ctx["imputation"] = outputs
    return outputs
