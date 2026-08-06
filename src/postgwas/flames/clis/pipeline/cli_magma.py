# postgwas/flames/clis/pipeline/cli_magma.py

from postgwas.utils.main import validate_path


def add_magma_arguments(parser):
    magma = parser.add_argument_group("2. MAGMA Settings")

    magma.add_argument(
        "--gene_loc",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".loc"]
        ),
        help=(
            "[bold red]Required[/bold red]: MAGMA gene location file (.loc).\n"
        ),
        metavar=""
    )

    magma.add_argument(
        "--pops_gene_loc",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".loc"]
        ),
        help=(
            "[bold red]Required[/bold red]: MAGMA gene location file (.loc).\n"
            "Genes should be same in PoPS gene annotation file and this file. must contain these columns ENSGID, CHR, START, END, NAME"
        ),
        metavar=""
    )
    
    magma.add_argument(
        "--window_upstream",
        type=int,
        default=35,
        help="Upstream locus extension (kb).",
        metavar=""
    )

    magma.add_argument(
        "--window_downstream",
        type=int,
        default=10,
        help="Downstream locus extension (kb).",
        metavar=""
    )

    magma.add_argument(
        "--magma_exe",
        type=str,
        help="[bold red]Required[/bold red]: Path to the MAGMA executable.",
        metavar=""
    )

    magma.add_argument(
        "--gene_model",
        choices=["snp-wise=mean", "snp-wise=top", "snp-wise=all"],
        default="snp-wise=mean",
        metavar="",
        help="MAGMA gene model for aggregating SNP statistics."
    )

    magma.add_argument(
        "--covar_file",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=["log2TPM.txt", "log2TPM.txt.gz"]
        ),
        metavar="",
        help=(
            "MAGMA covariate file (e.g., GTEx v8 expression matrix).\n"
            "Available at:\n"
            "  https://github.com/Marijn-Schipper/FLAMES\n"
            "[bold red]Required[/bold red]"
        )
    )

    return parser
