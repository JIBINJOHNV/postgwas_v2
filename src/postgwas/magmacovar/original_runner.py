# ==========================================================
# Model syntax validator
# ==========================================================
def validate_model_string(model_str: str) -> str:
    """
    Validate MAGMA model syntax for gene-property analysis.

    Allowed prefixes:
      condition, condition-hide, condition-residualize,
      condition-interaction, interaction-each, interaction-all,
      analyse, joint, interaction.

    Examples:
      condition-hide=Average
      condition=Brain_Cortex,Brain_Cerebellum
      analyse
      joint
      interaction
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

    # Keyword-only cases
    if model_str in ("analyse", "joint", "interaction"):
        return model_str

    # Must have '=' for other cases
    if "=" not in model_str:
        raise ArgumentTypeError(
            f"❌ Invalid model format '{model_str}'. "
            f"Expected something like condition-hide=Average"
        )

    prefix, values = model_str.split("=", 1)

    if prefix not in allowed:
        raise ArgumentTypeError(
            f"❌ Invalid model type '{prefix}'. Allowed: {', '.join(allowed)}"
        )

    # Check empty pieces in comma-separated values
    vals = values.split(",")
    if any(v.strip() == "" for v in vals):
        raise ArgumentTypeError(
            f"❌ Invalid model list '{model_str}': contains empty entries"
        )

    return model_str




# ==========================================================
# Helper: attach model arguments to a parser
# ==========================================================
def add_covar_model_arguments(parser: argparse.ArgumentParser) -> None:
    """
    Add MAGMA covariate model arguments to a parser.

    These are shared between DIRECT and PIPELINE modes.
    """
    model_group = parser.add_argument_group("MAGMA Covariate Model Options")

    model_group.add_argument(
        "--covar_model",
        metavar=" ",
        type=validate_model_string,
        help=(
            "Specify a single MAGMA covariate model.\n"
            "Supported formats:\n"
            "  condition=<v1,v2,...>\n"
            "  condition-hide=<v1,v2,...>\n"
            "  condition-residualize=<v1,v2,...>\n"
            "  condition-interaction=<v1,v2,...>\n"
            "  interaction-each=<v1,v2,...>\n"
            "  interaction-all=<v1,v2,...>\n"
            "  analyse\n"
            "  joint\n"
            "  interaction\n\n"
            "Examples:\n"
            "  --covar_model condition-hide=Average\n"
            "  --covar_model condition=Brain_Cortex,Brain_Cerebellum\n\n"
            "NOTE for FLAMES users:\n"
            "  FLAMES typically DOES NOT recommend using a custom --covar_model\n"
            "  unless you know exactly what you're doing."
        ),
    )
    
    

# ==========================================================
# PIPELINE MODE ENGINE
# ==========================================================
def _infer_gene_results_path(args: argparse.Namespace, magma_output) -> str:
    """
    Infer MAGMA gene results (.genes.raw) path from magma_analysis_pipeline output
    and/or conventions.

    Adjust this if your magma_analysis_pipeline uses a different naming scheme.
    """
    # First try to read from returned dict, if any
    gene_results_file = None
    if magma_output is not None:
        for key in ("gene_results_file", "gene_results", "genes_raw_file"):
            if isinstance(magma_output, dict) and key in magma_output:
                gene_results_file = magma_output[key]
                break

    # Fallback: standard PostGWAS-style naming inside outdir
    if gene_results_file is None:
        # You may want to change this to your actual MAGMA output template.
        candidate = Path(args.outdir) / f"{args.sample_id}.genes.raw"
        gene_results_file = str(candidate)

    return str(gene_results_file)

