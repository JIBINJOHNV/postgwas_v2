import argparse
from pathlib import Path


# ==========================================================
# Model syntax validator
# ==========================================================
def validate_model_string(model_str: str | None) -> str | None:
    """Validate one MAGMA --model modifier."""
    if model_str is None:
        return None

    model_str = model_str.strip()
    if not model_str:
        raise argparse.ArgumentTypeError("MAGMA model modifier cannot be empty")

    # Valid modifiers that do not accept a value.
    flag_modifiers = {
        "joint-pairs",
        "interaction-pairs",
    }

    # Valid modifiers accepting variable names or file mode.
    list_modifiers = {
        "condition",
        "condition-hide",
        "condition-residualize",
        "condition-interaction",
        "interaction-each",
        "interaction-all",
    }

    file_modifiers = {
        "joint",
        "interaction",
    }

    if model_str in flag_modifiers:
        return model_str

    if "=" not in model_str:
        if model_str in {"analyse", "joint", "interaction"}:
            raise argparse.ArgumentTypeError(
                f"'{model_str}' requires a value. Examples: "
                "analyse=sets, joint=models.txt, "
                "interaction=interactions.txt"
            )

        raise argparse.ArgumentTypeError(
            f"Unknown or incomplete MAGMA model modifier: '{model_str}'"
        )

    prefix, raw_values = model_str.split("=", 1)
    prefix = prefix.strip()
    values = [value.strip() for value in raw_values.split(",")]

    if prefix in flag_modifiers:
        raise argparse.ArgumentTypeError(
            f"'{prefix}' does not accept a value; use '{prefix}'"
        )

    valid_prefixes = list_modifiers | file_modifiers | {"analyse"}
    if prefix not in valid_prefixes:
        raise argparse.ArgumentTypeError(
            f"Invalid MAGMA model modifier '{prefix}'. "
            f"Allowed: {', '.join(sorted(valid_prefixes | flag_modifiers))}"
        )

    if not values or any(not value for value in values):
        raise argparse.ArgumentTypeError(
            f"Invalid MAGMA model modifier '{model_str}': "
            "contains an empty value"
        )

    # analyse=sets
    # analyse=covar
    # analyse=list,Variable1,Variable2
    # analyse=file,variables.txt
    if prefix == "analyse":
        single_modes = {"sets", "cov", "covar", "covariates"}

        if len(values) == 1 and values[0] in single_modes:
            return f"{prefix}={values[0]}"

        if values[0] == "list" and len(values) >= 2:
            return f"{prefix}={','.join(values)}"

        if values[0] == "file" and len(values) == 2:
            return f"{prefix}={','.join(values)}"

        raise argparse.ArgumentTypeError(
            f"Invalid analyse specification '{model_str}'. Valid forms:\n"
            "  analyse=sets\n"
            "  analyse=covar\n"
            "  analyse=list,Variable1,Variable2\n"
            "  analyse=file,variables.txt"
        )

    # joint and interaction require one model-definition file.
    if prefix in file_modifiers:
        if len(values) != 1:
            raise argparse.ArgumentTypeError(
                f"'{prefix}' requires exactly one model-definition file. "
                f"Example: {prefix}=models.txt"
            )
        return f"{prefix}={values[0]}"

    # File mode for condition and interaction-each/all modifiers:
    # condition=file,variables.txt
    if values[0] == "file":
        if len(values) != 2:
            raise argparse.ArgumentTypeError(
                f"Invalid file mode '{model_str}'. Expected: "
                f"{prefix}=file,filename"
            )
        return f"{prefix}={','.join(values)}"

    # Direct condition-interaction variables must occur in pairs.
    if prefix == "condition-interaction":
        if len(values) < 2 or len(values) % 2 != 0:
            raise argparse.ArgumentTypeError(
                "condition-interaction requires an even number of variables "
                "forming interaction pairs. Example: "
                "condition-interaction=Variable1,Variable2"
            )

    return f"{prefix}={','.join(values)}"


# ==========================================================
# Add model arguments to a parser
# ==========================================================
def add_covar_model_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    """Add repeatable MAGMA model modifiers to an argument parser."""
    model_group = parser.add_argument_group(
        "MAGMA Covariate Model Options"
    )

    model_group.add_argument(
        "--covar_model",
        action="extend",
        nargs="+",
        type=validate_model_string,
        metavar="MODEL",
        help=(
            "Specify one or more MAGMA --model modifiers.\n\n"
            "Examples:\n"
            "  --covar_model condition-hide=Average\n"
            "  --covar_model "
            "condition=Brain_Cortex,Brain_Cerebellum\n"
            "  --covar_model analyse=sets\n"
            "  --covar_model joint=models.txt\n"
            "  --covar_model interaction=interactions.txt\n"
            "  --covar_model joint-pairs\n"
            "  --covar_model interaction-pairs\n\n"
            "Multiple modifiers can be combined:\n"
            "  --covar_model "
            "analyse=sets condition-hide=Average\n\n"
            "NOTE for FLAMES users:\n"
            "  FLAMES typically does not recommend a custom model "
            "unless you understand its statistical consequences."
        ),
    )


# ==========================================================
# Infer MAGMA gene-results path
# ==========================================================
def _infer_gene_results_path(
    args: argparse.Namespace,
    magma_output,
) -> str:
    """Return and validate the MAGMA .genes.raw output path."""
    gene_results_file = None

    if isinstance(magma_output, dict):
        # Exact file paths returned by the pipeline functions.
        file_keys = (
            "magma_genes_raw",
            "genes_raw",
            "gene_results_file",
            "gene_results",
            "genes_raw_file",
        )

        for key in file_keys:
            value = magma_output.get(key)
            if value:
                gene_results_file = value
                break

        # Try a returned prefix if an exact filepath was not supplied.
        if gene_results_file is None:
            for key in ("magma_genes_prefix", "merged_prefix"):
                value = magma_output.get(key)
                if value:
                    gene_results_file = f"{value}.genes.raw"
                    break

    elif isinstance(magma_output, (str, Path)):
        gene_results_file = magma_output

    # Backward-compatible fallback using the current naming convention.
    if gene_results_file is None:
        outdir = getattr(args, "outdir", None)
        if outdir is None:
            outdir = getattr(args, "output_dir", None)

        sample_id = getattr(args, "sample_id", None)

        if not outdir:
            raise ValueError(
                "Cannot infer MAGMA gene-results path: "
                "args.outdir or args.output_dir is missing"
            )

        if not sample_id:
            raise ValueError(
                "Cannot infer MAGMA gene-results path: "
                "args.sample_id is missing"
            )

        upstream = getattr(args, "window_upstream", 35)
        downstream = getattr(args, "window_downstream", 10)

        gene_results_file = (
            Path(outdir)
            / (
                f"{sample_id}_magma_"
                f"{upstream}up_{downstream}down.genes.raw"
            )
        )

    gene_results_path = Path(gene_results_file).expanduser()

    if not gene_results_path.is_file():
        raise FileNotFoundError(
            "MAGMA gene-results file was not found.\n"
            f"Expected path: {gene_results_path}"
        )

    if gene_results_path.stat().st_size == 0:
        raise ValueError(
            "MAGMA gene-results file exists but is empty.\n"
            f"File: {gene_results_path}"
        )

    return str(gene_results_path.resolve())