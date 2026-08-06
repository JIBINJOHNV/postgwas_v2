import argparse
from postgwas.utils.main import validate_path
from rich_argparse import RichHelpFormatter

def get_geneset_common_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for gene set over representation analysis.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("pathway enrichment Inputs")

    grp.add_argument(
        "--gene_inputfile",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help="Gene file (first column = gene symbols)",
    )

    grp.add_argument(
        "--biogrid-key",
        required=True,
        help=(
            "BioGRID API key (required). "
            "Used to fetch protein–protein interactions from BioGRID. "
            "Register and obtain a key at: https://webservice.thebiogrid.org/"
        ),
    )

    grp.add_argument(
        "--david-email",
        required=True,
        help=(
            "Registered email for DAVID enrichment (required). "
            "Must be the same email used to register with DAVID. "
            "Register at: https://david.ncifcrf.gov/webservice/register.html"
        ),
    )
    grp.add_argument("--dsigdb-gmt", help="DSigDB GMT file; https://dsigdb.tanlab.org/DSigDBv1.0/download.html")
    grp.add_argument("--reference-set", default="genome_protein-coding")
    return parser