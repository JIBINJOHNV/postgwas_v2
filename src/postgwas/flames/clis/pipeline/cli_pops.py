# postgwas/flames/clis/pipeline/cli_pops.py

from postgwas.utils.main import validate_path,validate_prefix_files


def add_pops_arguments(parser):
    """
    Add PoPS-specific settings, covariates, feature-selection,
    training parameters, and runtime behavior.
    """

    # -----------------------------------------------------
    # 4. PoPS SETTINGS
    # -----------------------------------------------------
    pops = parser.add_argument_group("4. PoPS Settings")

    pops.add_argument(
        "--pops_gene_annot",
        type=validate_path(must_exist=True, must_be_file=True),
        help="[bold bright_red]Required[/bold bright_red]: PoPS gene annotation file (must contain these columns ENSGID, CHR, TSS).",
        metavar=""
    )
    
    pops.add_argument(
        "--feature_prefix",
        type=lambda p: validate_prefix_files(p, ["*.npy", "*.txt"]),
        help="[bold bright_red]Required[/bold bright_red]: Prefix for PoPS feature matrix chunk files.",
        metavar=""
    )

    pops.add_argument(
        "--pops_cov",
        type=validate_path(must_exist=True, must_be_file=True),
        help="Optional PoPS covariance matrix (.npz).",
        metavar=""
    )

    # -----------------------------------------------------
    # 5. PoPS COVARIATES
    # -----------------------------------------------------
    cov = parser.add_argument_group("5. PoPS Covariates")

    # MAGMA Covariates
    cov.add_argument(
        "--use_magma_covariates",
        dest="use_magma_covariates",
        action="store_true",
        help="Optional Project out MAGMA covariates before fitting."
    )
    cov.add_argument(
        "--ignore_magma_covariates",
        dest="use_magma_covariates",
        action="store_false",
        help="Optional Ignore MAGMA covariates."
    )
    parser.set_defaults(use_magma_covariates=True)

    # MAGMA Error Covariance
    cov.add_argument(
        "--use_magma_error_cov",
        dest="use_magma_error_cov",
        action="store_true",
        help="Optional Use the MAGMA error covariance when fitting."
    )
    cov.add_argument(
        "--ignore_magma_error_cov",
        dest="use_magma_error_cov",
        action="store_false",
        help="Optional Ignore the MAGMA error covariance."
    )
    parser.set_defaults(use_magma_error_cov=True)

    # Custom target score
    cov.add_argument(
        "--y_path",
        type=str,
        metavar="",
        help="Optional custom target score. Must contain ENSGID and Score columns."
    )
    cov.add_argument(
        "--y_covariates_path",
        type=str,
        metavar="",
        help="Optional covariate matrix for custom target score."
    )
    cov.add_argument(
        "--y_error_cov_path",
        type=str,
        metavar="",
        help="Optional error covariance for custom target score (.npz or .npy)."
    )

    # Chromosome control
    cov.add_argument(
        "--project_out_covariates_chromosomes",
        nargs="*",
        metavar="",
        help="Optional Chromosomes to consider when projecting out covariates."
    )

    cov.add_argument(
        "--project_out_covariates_remove_hla",
        dest="project_out_covariates_remove_hla",
        action="store_true",
        help="Optional Remove HLA genes before projecting out covariates."
    )
    cov.add_argument(
        "--project_out_covariates_keep_hla",
        dest="project_out_covariates_remove_hla",
        action="store_false",
        help="Optional Keep HLA genes when projecting out covariates."
    )
    parser.set_defaults(project_out_covariates_remove_hla=True)

    # -----------------------------------------------------
    # 6. PoPS FEATURE SELECTION
    # -----------------------------------------------------
    sel = parser.add_argument_group("6. PoPS Feature Selection")

    sel.add_argument(
        "--subset_features_path",
        type=str,
        metavar="",
        help="Optional list of features (one per line) to subset."
    )
    sel.add_argument(
        "--control_features_path",
        type=str,
        metavar="",
        help="Optional list of always-include features."
    )
    sel.add_argument(
        "--feature_selection_chromosomes",
        nargs="*",
        metavar="",
        help="Optional Chromosomes to consider during feature selection."
    )
    sel.add_argument(
        "--feature_selection_p_cutoff",
        type=float,
        default=0.05,
        metavar="",
        help="Optional P-value cutoff for feature selection."
    )
    sel.add_argument(
        "--feature_selection_max_num",
        type=int,
        metavar="",
        help="Optional Maximum number of features to select (excluding controls)."
    )
    sel.add_argument(
        "--feature_selection_fss_num_features",
        type=int,
        metavar="",
        help="Optional Number of features for forward stepwise selection (FSS)."
    )

    sel.add_argument(
        "--feature_selection_remove_hla",
        dest="feature_selection_remove_hla",
        action="store_true",
        help="Optional Remove HLA genes during feature selection."
    )
    sel.add_argument(
        "--feature_selection_keep_hla",
        dest="feature_selection_remove_hla",
        action="store_false",
        help="Optional Keep HLA genes during feature selection."
    )
    parser.set_defaults(feature_selection_remove_hla=True)

    # -----------------------------------------------------
    # 7. PoPS TRAINING
    # -----------------------------------------------------
    train = parser.add_argument_group("7. PoPS Training")

    train.add_argument(
        "--training_chromosomes",
        nargs="*",
        metavar="",
        help="Optional Chromosomes used to train model coefficients."
    )

    train.add_argument(
        "--training_remove_hla",
        dest="training_remove_hla",
        action="store_true",
        help="Optional Remove HLA genes when computing model coefficients."
    )
    train.add_argument(
        "--training_keep_hla",
        dest="training_remove_hla",
        action="store_false",
        help="Optional Keep HLA genes when computing model coefficients."
    )
    parser.set_defaults(training_remove_hla=True)

    # -----------------------------------------------------
    # 8. PoPS RUNTIME SETTINGS
    # -----------------------------------------------------
    run = parser.add_argument_group("8. PoPS Runtime Settings")

    run.add_argument(
        "--method",
        default="ridge",
        metavar="",
        help="Optional Regularization method (ridge, lasso, linreg)."
    )

    run.add_argument(
        "--save_matrix_files",
        dest="save_matrix_files",
        action="store_true",
        help="Optional Save PoPS matrices used for coefficient and score calculation."
    )
    run.add_argument(
        "--no_save_matrix_files",
        dest="save_matrix_files",
        action="store_false",
        help="Optional Do not save PoPS matrices."
    )
    parser.set_defaults(save_matrix_files=False)

    run.add_argument(
        "--random_seed",
        type=int,
        default=42,
        metavar="",
        help="Optional Random seed for reproducibility."
    )

    run.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Optional Enable verbose output."
    )
    run.add_argument(
        "--no_verbose",
        dest="verbose",
        action="store_false",
        help="Optional Disable verbose output."
    )
    parser.set_defaults(verbose=False)

    return parser
