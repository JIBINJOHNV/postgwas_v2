import argparse
from postgwas.core.execution.runtime import validate_path

def get_geneset_common_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for gene set over representation analysis.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Pathway-enrichment inputs")

    grp.add_argument(
        "--gene-input-file",
        metavar="PATH",
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
    grp.add_argument(
        "--dsigdb-gmt", metavar="PATH",
        help="Optional DSigDB GMT file downloaded from dsigdb.tanlab.org.",
    )
    grp.add_argument("--reference-set", default="genome_protein-coding")
    return parser
