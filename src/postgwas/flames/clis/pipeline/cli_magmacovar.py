# postgwas/flames/clis/pipeline/cli_magma_covar.py

from __future__ import annotations
from argparse import ArgumentTypeError
from postgwas.utils.main import validate_path


# -----------------------------
# VALIDATE MODEL STRING
# -----------------------------
def validate_model_string(model_str: str) -> str:
    """
    Validate MAGMA model string formats.

    Accepts:
      - analyse
      - joint
      - interaction
      - condition=<v1,v2,...>
      - condition-hide=<v1,v2,...>
      - condition-residualize=<v1,v2,...>
      - condition-interaction=<v1,v2,...>
      - interaction-each=<v1,v2,...>
      - interaction-all=<v1,v2,...>
    """
    if model_str is None:
        return None

    allowed = [
        "condition",
        "condition-hide",
        "condition-residualize",
        "condition-interaction",
        "interaction-each",
        "interaction-all",
        "analyse",
        "joint",
        "interaction",
    ]

    # Case: keyword-only models
    if model_str in ("analyse", "joint", "interaction"):
        return model_str

    # Must contain "="
    if "=" not in model_str:
        raise ArgumentTypeError(
            f"❌ Invalid model format '{model_str}'. Expected condition-hide=Average"
        )

    prefix, values = model_str.split("=", 1)

    if prefix not in allowed:
        raise ArgumentTypeError(
            f"❌ Invalid model type '{prefix}'. Allowed: {', '.join(allowed)}"
        )

    vals = values.split(",")
    if any(v.strip() == "" for v in vals):
        raise ArgumentTypeError(
            f"❌ Invalid model list '{model_str}': contains empty entries"
        )

    return model_str


# -----------------------------
# ADD ARGUMENTS
# -----------------------------
def add_magma_covar_arguments(parser):
    """
    Add ONLY:
      • --covariates
      • --model
      • --direction

    This file is meant to be imported into the main pipeline CLI.
    """

    magma_covar = parser.add_argument_group("MAGMA Covariates")

    magma_covar.add_argument(
        "--covariates",
        type=lambda p: validate_path(
            must_exist=True,
            must_be_file=True,
            allowed_suffixes=[".tsv", ".txt"],
            must_not_be_empty=True,
        ),
        help="Covariate matrix file (TSV/TXT). Must not be empty.",
        metavar=""
    )

    # ----------------------- MODEL OPTIONS ------------------------
    #model_group = parser.add_argument_group("MAGMA Model Options")

    magma_covar.add_argument(
        "--model",
        metavar="",
        type=validate_model_string,
        help=(
            "MAGMA gene model.\n"
            "Formats:\n"
            "  condition=<v1,v2,...>\n"
            "  condition-hide=<v1>\n"
            "  condition-residualize=<v1,v2>\n"
            "  condition-interaction=<v1>\n"
            "  interaction-each=<v1,v2>\n"
            "  interaction-all=<v1>\n"
            "  analyse | joint | interaction\n"
        )
    )

    magma_covar.add_argument(
        "--direction",
        choices=["two-sided", "greater", "less"],
        default="two-sided",
        help=(
            "MAGMA test direction. Passed to MAGMA as a free token.\n"
            "Example:\n"
            "  --model condition-hide=Average --direction greater\n"
            "MAGMA sees:\n"
            "   condition-hide=Average direction=greater"
        ),
        metavar=""
    )

    return parser
