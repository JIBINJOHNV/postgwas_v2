
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# Internal imports
from postgwas.pops.pops import pops_main, get_pops_args



def run_pops_direct(args: argparse.Namespace, ctx=None):
    """
    Run PoPS directly using provided arguments.
    This bypasses the full PostGWAS pipeline and runs PoPS only.
    """
    # --------------------------------------------
    # Create output folder
    # --------------------------------------------
    pops_folder = Path(args.outdir)
    pops_folder.mkdir(parents=True, exist_ok=True)

    pops_out_prefix = pops_folder / f"{args.sample_id}_pops"

    # --------------------------------------------
    # Build argument list for PoPS
    # --------------------------------------------
    pops_argv = [
        "--gene_annot_path", args.pops_gene_loc_file,
        "--feature_mat_prefix", args.feature_mat_prefix,
        "--num_feature_chunks", str(args.num_feature_chunks),
        "--magma_prefix", args.magma_assoc_prefix,
        "--out_prefix", str(pops_out_prefix),
    ]

    # Optional: Control features
    if getattr(args, "control_features_path", None):
        pops_argv += ["--control_features_path", args.control_features_path]

    # Optional flags: use or ignore MAGMA covariates
    pops_argv.append(
        "--use_magma_covariates" if args.use_magma_covariates else "--ignore_magma_covariates"
    )
    pops_argv.append(
        "--use_magma_error_cov" if args.use_magma_error_cov else "--ignore_magma_error_cov"
    )

    # Optional Y data
    if getattr(args, "y_path", None):
        pops_argv += ["--y_path", args.y_path]
    if getattr(args, "y_covariates_path", None):
        pops_argv += ["--y_covariates_path", args.y_covariates_path]
    if getattr(args, "y_error_cov_path", None):
        pops_argv += ["--y_error_cov_path", args.y_error_cov_path]

    # Project-out chromosomes
    if getattr(args, "project_out_covariates_chromosomes", None):
        pops_argv += [
            "--project_out_covariates_chromosomes",
            *args.project_out_covariates_chromosomes
        ]

    pops_argv.append(
        "--project_out_covariates_remove_hla"
        if args.project_out_covariates_remove_hla
        else "--project_out_covariates_keep_hla"
    )

    # Subset features
    if getattr(args, "subset_features_path", None):
        pops_argv += ["--subset_features_path", args.subset_features_path]

    # Feature selection chromosomes
    if getattr(args, "feature_selection_chromosomes", None):
        pops_argv += [
            "--feature_selection_chromosomes",
            *args.feature_selection_chromosomes,
        ]

    pops_argv += [
        "--feature_selection_p_cutoff",
        str(args.feature_selection_p_cutoff)
    ]

    if args.feature_selection_max_num is not None:
        pops_argv += [
            "--feature_selection_max_num",
            str(args.feature_selection_max_num),
        ]

    if args.feature_selection_fss_num_features is not None:
        pops_argv += [
            "--feature_selection_fss_num_features",
            str(args.feature_selection_fss_num_features),
        ]

    pops_argv.append(
        "--feature_selection_remove_hla"
        if args.feature_selection_remove_hla
        else "--feature_selection_keep_hla"
    )

    if getattr(args, "training_chromosomes", None):
        pops_argv += [
            "--training_chromosomes",
            *args.training_chromosomes,
        ]

    pops_argv.append(
        "--training_remove_hla"
        if args.training_remove_hla
        else "--training_keep_hla"
    )

    pops_argv += [
        "--method", args.method,
        "--random_seed", str(args.seed),
    ]

    pops_argv.append(
        "--save_matrix_files"
        if args.save_matrix_files
        else "--no_save_matrix_files"
    )

    pops_argv.append(
        "--verbose"
        if args.pops_verbose
        else "--no_verbose"
    )
    # --------------------------------------------
    # Run PoPS
    # --------------------------------------------
    try:
        pops_args = get_pops_args(pops_argv)
        pops_main(vars(pops_args))
    except SystemExit as e:
        raise RuntimeError("PoPS exited prematurely") from e

    print("\n           🎉 PoPS Analysis Completed Successfully!")
    print()

    pops_out_file = f"{pops_out_prefix}.preds"

    result = {
        "status": "success",
        "pops_file": pops_out_file,
        "output_dir": str(pops_folder),
        "out_prefix": str(pops_out_prefix),
    }
    ctx["pops_output"] = pops_out_file
    return result
