#!/usr/bin/env python3
import argparse
import sys

# ---------------------------------------------------------
# IMPORTS
# Ensure postgwas.pops.pipeline exists and contains run_pipeline_logic
# ---------------------------------------------------------
from postgwas.pops.main import run_pipeline_logic
from postgwas.pops.pops import pops_main, get_pops_args


def run_direct_mode(remaining):
    """Direct passthrough to PoPS."""
    if not remaining or remaining == ["--help"] or remaining == ["-h"]:
        get_pops_args(["--help"])
        return
    pops_args = get_pops_args(remaining)
    return pops_main(vars(pops_args))


def main():
    if len(sys.argv) == 1:
        print("\n🔧 postgwas pops — PoPS Analysis Tool\n")
        print("   postgwas pops direct --help")
        print("   postgwas pops pipeline --help\n")
        sys.exit(0)

    parser = argparse.ArgumentParser(
        prog="postgwas pops",
        description="Run PoPS directly or in pipeline mode (VCF -> Magma -> PoPS).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="cmd")

    # ---------------- DIRECT MODE ----------------
    p_direct = subparsers.add_parser("direct", help="Run PoPS directly (pass through args).", add_help=False)
    p_direct.add_argument("remaining", nargs=argparse.REMAINDER)

    # ---------------- PIPELINE MODE ----------------
    p_pipe = subparsers.add_parser(
        "pipeline",
        help="Run full VCF → MAGMA → PoPS workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 1. Required pipeline inputs
    grp_req = p_pipe.add_argument_group("1. General Pipeline Inputs")
    grp_req.add_argument("--sumstat_vcf", required=True, help="GWAS VCF file.")
    grp_req.add_argument("--output_folder", required=True, help="Directory where results will be saved.")
    grp_req.add_argument("--sample_name", required=True, help="Sample identifier used for naming output files.")

    # 2. MAGMA args
    grp_magma = p_pipe.add_argument_group("2. MAGMA Analysis Parameters")
    grp_magma.add_argument("--ld_ref", required=True, help="Path to LD reference data (bed/bim/fam prefix).")
    grp_magma.add_argument(
        "--magma_gene_loc_file",
        required=True,
        help="User-supplied MAGMA gene location file (.loc)."
    )
    grp_magma.add_argument("--geneset_file", default=None, help="Optional geneset file for MAGMA analysis.")
    grp_magma.add_argument("--window_upstream", type=int, default=35, help="Window upstream (kb).")
    grp_magma.add_argument("--window_downstream", type=int, default=10, help="Window downstream (kb).")
    grp_magma.add_argument("--gene_model", default="snp-wise=mean", help="MAGMA gene model definition.")
    grp_magma.add_argument("--n_sample_col", default="N_COL", help="Column name for sample size in VCF.")
    grp_magma.add_argument("--num_batches", type=int, default=6, help="Number of batches for processing.")
    grp_magma.add_argument("--num_cores", type=int, default=6, help="Number of cores to use.")
    grp_magma.add_argument("--seed", type=int, default=10, help="Random seed for MAGMA.")
    grp_magma.add_argument("--magma", default="magma", help="Path to magma executable.")

    # 3. PoPS Required Files
    grp_feats = p_pipe.add_argument_group("3. PoPS Required Files")
    grp_feats.add_argument(
        "--feature_mat_prefix", 
        required=True,
        help="Prefix to the split feature matrix files. There must be .mat.*.npy files, .cols.*.txt files, and a .rows.txt file."
    )
    grp_feats.add_argument(
        "--num_feature_chunks", 
        type=int, 
        required=True,
        help="The number of feature matrix chunks."
    )
    grp_feats.add_argument(
        "--pops_gene_loc_file",
        required=True,
        help="Path to tab-separated gene annotation file for PoPS (must contain ENSGID, CHR, and TSS columns)."
    )
    
    # 4. PoPS Covariates
    grp_cov = p_pipe.add_argument_group("4. PoPS Covariates")
    grp_cov.add_argument('--use_magma_covariates', dest='use_magma_covariates', action='store_true', help="Project out MAGMA covariates before fitting.")
    grp_cov.add_argument('--ignore_magma_covariates', dest='use_magma_covariates', action='store_false', help="Ignore MAGMA covariates.")
    p_pipe.set_defaults(use_magma_covariates=True)

    grp_cov.add_argument('--use_magma_error_cov', dest='use_magma_error_cov', action='store_true', help="Use the MAGMA error covariance when fitting.")
    grp_cov.add_argument('--ignore_magma_error_cov', dest='use_magma_error_cov', action='store_false', help="Ignore the MAGMA error covariance when fitting.")
    p_pipe.set_defaults(use_magma_error_cov=True)

    grp_cov.add_argument("--y_path", help="Path to a custom target score. Must contain ENSGID and Score columns.")
    grp_cov.add_argument("--y_covariates_path", help="Optional path to covariates for custom target score provided in --y_path.")
    grp_cov.add_argument("--y_error_cov_path", help="Optional path to error covariance for custom target score provided in --y_path (npz or npy).")

    grp_cov.add_argument("--project_out_covariates_chromosomes", nargs="*", help="List chromosomes to consider when projecting out covariates.")
    grp_cov.add_argument('--project_out_covariates_remove_hla', dest='project_out_covariates_remove_hla', action='store_true', help="Remove HLA genes before projecting out covariates.")
    grp_cov.add_argument('--project_out_covariates_keep_hla', dest='project_out_covariates_remove_hla', action='store_false', help="Keep HLA genes when projecting out covariates.")
    p_pipe.set_defaults(project_out_covariates_remove_hla=True)

    # 5. PoPS Feature Selection
    grp_sel = p_pipe.add_argument_group("5. PoPS Feature Selection")
    grp_sel.add_argument("--subset_features_path", help="Optional path to list of features (one per line) to subset to.")
    grp_sel.add_argument("--control_features_path", help="Optional path to list of features (one per line) to always include.")
    grp_sel.add_argument("--feature_selection_chromosomes", nargs="*", help="List chromosomes to consider when performing feature selection.")
    grp_sel.add_argument("--feature_selection_p_cutoff", type=float, default=0.05, help="P-value cutoff to use when performing feature selection.")
    grp_sel.add_argument("--feature_selection_max_num", type=int, help="Maximum number of features to select, excluding control features.")
    grp_sel.add_argument("--feature_selection_fss_num_features", type=int, help="Number of features to select using forward stepwise selection (overrides other selection args).")
    
    grp_sel.add_argument('--feature_selection_remove_hla', dest='feature_selection_remove_hla', action='store_true', help="Remove HLA genes when performing feature selection.")
    grp_sel.add_argument('--feature_selection_keep_hla', dest='feature_selection_remove_hla', action='store_false', help="Keep HLA genes when performing feature selection.")
    p_pipe.set_defaults(feature_selection_remove_hla=True)

    # 6. PoPS Training
    grp_train = p_pipe.add_argument_group("6. PoPS Training")
    grp_train.add_argument("--training_chromosomes", nargs="*", help="List chromosomes to consider when computing model coefficients.")
    grp_train.add_argument('--training_remove_hla', dest='training_remove_hla', action='store_true', help="Remove HLA genes when computing model coefficients.")
    grp_train.add_argument('--training_keep_hla', dest='training_remove_hla', action='store_false', help="Keep HLA genes when computing model coefficients.")
    p_pipe.set_defaults(training_remove_hla=True)

    # 7. Runtime Settings
    grp_misc = p_pipe.add_argument_group("7. Runtime Settings")
    grp_misc.add_argument("--method", default="ridge", help="Regularization used (ridge, lasso, linreg).")
    
    grp_misc.add_argument('--save_matrix_files', dest='save_matrix_files', action='store_true', help="Save matrices used to compute model coefficients and PoP score.")
    grp_misc.add_argument('--no_save_matrix_files', dest='save_matrix_files', action='store_false', help="Do not save matrices.")
    p_pipe.set_defaults(save_matrix_files=False)
    
    grp_misc.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility.")
    
    grp_misc.add_argument('--verbose', dest='verbose', action='store_true', help="Get verbose output.")
    grp_misc.add_argument('--no_verbose', dest='verbose', action='store_false', help="Silence output.")
    p_pipe.set_defaults(verbose=False)

    # Dispatch
    args, unknown = parser.parse_known_args()

    if args.cmd == "direct":
        return run_direct_mode(unknown)
    elif args.cmd == "pipeline":
        if unknown:
            print(f"Warning: Unknown arguments provided to pipeline mode: {unknown}")
        
        # Call the logic function from pipeline.py
        return run_pipeline_logic(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()