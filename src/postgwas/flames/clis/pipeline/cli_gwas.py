# postgwas/flames/clis/pipeline/cli_gwas.py

from postgwas.utils.main import validate_path, validate_plink_prefix


def add_gwas_arguments(parser):
    gwas = parser.add_argument_group("1. GWAS Input")

    gwas.add_argument(
        "--vcf",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".vcf", ".vcf.gz"]
        ),
        help=(
            "QC-passed GWAS summary statistics VCF.\n"
            "[bold yellow]Must come from postgwas harmonisation.[/bold yellow]\n"
            "[bold red]Required[/bold red]"
        ),
        metavar=""
    )

    gwas.add_argument(
        "--ld-ref",
        type=validate_plink_prefix,
        dest="ld_ref_panel",
        metavar="",
        help=(
            "PLINK reference panel prefix (without .bed/.bim/.fam).\n"
            "[bold cyan]Example:[/bold cyan] /ref/1000G_EUR\n"
            "[bold red]Required[/bold red]"
        )
    )

    return parser
