"""Validated MAGMA gene-property execution."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

from postgwas.core.paths import require_nonempty_file
from postgwas.core.processes import run_checked_command
from postgwas.modules.magmacovar.errors import MagmaCovarError


def _magma_gene_ids(path: Path) -> set[str]:
    """Read the gene rows from MAGMA's documented text ``.genes.raw`` format."""
    genes: set[str] = set()
    version_found = False
    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    if line[1:].strip().upper().startswith("VERSION"):
                        version_found = True
                    continue
                fields = line.split()
                if len(fields) < 7:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw row %d has fewer than seven required fields"
                        % line_number
                    )
                chromosome = fields[1].strip()
                try:
                    start = int(fields[2])
                    stop = int(fields[3])
                except ValueError as exc:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw row %d has invalid gene coordinates"
                        % line_number
                    ) from exc
                if not chromosome:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw row %d has an empty chromosome" % line_number
                    )
                if start < 0 or stop < start:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw row %d has an invalid start/stop interval"
                        % line_number
                    )
                gene = fields[0].strip()
                if not gene:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw contains an empty gene ID at line %d"
                        % line_number
                    )
                if gene in genes:
                    raise MagmaCovarError(
                        "MAGMA .genes.raw contains duplicate gene ID %r" % gene
                    )
                genes.add(gene)
    except UnicodeDecodeError as exc:
        raise MagmaCovarError(
            "MAGMA gene results are not valid UTF-8 text: %s" % path
        ) from exc
    except OSError as exc:
        raise MagmaCovarError(
            "Cannot read MAGMA gene results %s: %s" % (path, exc)
        ) from exc
    if not version_found:
        raise MagmaCovarError(
            "MAGMA gene results do not contain the required '# VERSION' metadata: %s"
            % path
        )
    if not genes:
        raise MagmaCovarError(
            "MAGMA gene results contain no readable gene rows: %s" % path
        )
    return genes


def validate_magma_covariate_table(
    covariates_file: str | Path,
    *,
    eligible_gene_ids: set[str] | None = None,
    minimum_genes: int,
    maximum_missing_fraction: float,
    missing_genes: str,
) -> dict:
    """Validate the covariate table, optionally against MAGMA's gene universe."""
    covariates = require_nonempty_file(
        covariates_file, "MAGMA gene-covariate file", error_type=MagmaCovarError,
    )
    header: list[str] | None = None
    covariate_genes: set[str] = set()
    missing_by_property: list[int] = []
    observed_values: list[set[float]] = []
    overlap = 0
    try:
        with covariates.open("r", encoding="utf-8", errors="strict") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split()
                if header is None:
                    header = fields
                    if len(header) < 2:
                        raise MagmaCovarError(
                            "MAGMA gene-covariate header must contain a gene-ID column "
                            "and at least one numeric property"
                        )
                    if len(set(header)) != len(header):
                        raise MagmaCovarError(
                            "MAGMA gene-covariate header contains duplicate column names"
                        )
                    missing_by_property = [0] * (len(header) - 1)
                    observed_values = [set() for _ in header[1:]]
                    continue
                if len(fields) != len(header):
                    raise MagmaCovarError(
                        "MAGMA gene-covariate row %d has %d fields; the header has %d"
                        % (line_number, len(fields), len(header))
                    )
                gene = fields[0].strip()
                if not gene:
                    raise MagmaCovarError(
                        "MAGMA gene-covariate row %d has an empty gene ID" % line_number
                    )
                if gene in covariate_genes:
                    raise MagmaCovarError(
                        "MAGMA gene-covariate file contains duplicate gene ID %r" % gene
                    )
                covariate_genes.add(gene)
                shared = eligible_gene_ids is None or gene in eligible_gene_ids
                if shared:
                    overlap += 1
                for index, raw_value in enumerate(fields[1:]):
                    if raw_value == "NA":
                        if shared:
                            missing_by_property[index] += 1
                        continue
                    try:
                        value = float(raw_value)
                    except ValueError as exc:
                        raise MagmaCovarError(
                            "MAGMA gene-covariate value at line %d, column %r is "
                            "neither numeric nor the required missing token NA"
                            % (line_number, header[index + 1])
                        ) from exc
                    if not math.isfinite(value):
                        raise MagmaCovarError(
                            "MAGMA gene-covariate value at line %d, column %r is not finite"
                            % (line_number, header[index + 1])
                        )
                    if shared and len(observed_values[index]) < 2:
                        observed_values[index].add(value)
    except UnicodeDecodeError as exc:
        raise MagmaCovarError(
            "MAGMA gene-covariate file is not valid UTF-8 text: %s" % covariates
        ) from exc
    except OSError as exc:
        raise MagmaCovarError(
            "Cannot read MAGMA gene covariates %s: %s" % (covariates, exc)
        ) from exc

    if header is None or not covariate_genes:
        raise MagmaCovarError(
            "MAGMA gene-covariate file has no data rows: %s" % covariates
        )
    if overlap < minimum_genes:
        if eligible_gene_ids is None:
            raise MagmaCovarError(
                "MAGMA gene-covariate file contains %d genes; at least %d are required"
                % (overlap, minimum_genes)
            )
        raise MagmaCovarError(
            "Only %d gene IDs overlap between MAGMA gene results and covariates; "
            "at least %d are required" % (overlap, minimum_genes)
        )

    absent_genes = (
        len(eligible_gene_ids - covariate_genes)
        if eligible_gene_ids is not None else 0
    )
    denominator = (
        len(eligible_gene_ids)
        if eligible_gene_ids is not None and missing_genes == "fill"
        else overlap
    )
    property_missing = []
    for index, name in enumerate(header[1:]):
        missing = missing_by_property[index]
        if missing_genes == "fill":
            missing += absent_genes
        missing_fraction = missing / denominator
        if missing_fraction > maximum_missing_fraction:
            raise MagmaCovarError(
                "Gene property %r has %.3f missing values among eligible genes, "
                "above the configured MAGMA maximum of %.3f"
                % (name, missing_fraction, maximum_missing_fraction)
            )
        if len(observed_values[index]) < 2:
            raise MagmaCovarError(
                "Gene property %r is constant or has fewer than two finite values "
                "among overlapping genes" % name
            )
        property_missing.append({
            "property": name,
            "missing_genes": missing,
            "missing_fraction": missing_fraction,
        })

    return {
        "covariates_file": str(covariates),
        "covariate_genes": len(covariate_genes),
        "overlapping_genes": overlap,
        "gene_results_without_covariates": absent_genes,
        "covariate_genes_without_gene_results": (
            len(covariate_genes - eligible_gene_ids)
            if eligible_gene_ids is not None else 0
        ),
        "properties": len(header) - 1,
        "property_names": header[1:],
        "property_missingness": property_missing,
    }


def validate_magma_covariate_inputs(
    gene_results_file: str | Path,
    covariates_file: str | Path,
    *,
    minimum_genes: int,
    maximum_missing_fraction: float,
    missing_genes: str,
) -> dict:
    """Validate MAGMA inputs and their gene-ID overlap without a dataframe."""
    gene_results = require_nonempty_file(
        gene_results_file, "MAGMA .genes.raw file", error_type=MagmaCovarError,
    )
    genes = _magma_gene_ids(gene_results)
    if len(genes) < minimum_genes:
        raise MagmaCovarError(
            "MAGMA gene results contain %d genes; at least %d are required"
            % (len(genes), minimum_genes)
        )
    summary = validate_magma_covariate_table(
        covariates_file,
        eligible_gene_ids=genes,
        minimum_genes=minimum_genes,
        maximum_missing_fraction=maximum_missing_fraction,
        missing_genes=missing_genes,
    )
    return {
        "gene_results_file": str(gene_results),
        "gene_results_genes": len(genes),
        **summary,
    }


def build_magma_covariate_command(
    magma_bin: str,
    gene_results_file: str | Path,
    covariates_file: str | Path,
    output_prefix: str | Path,
    *,
    model: Sequence[str],
    direction: str,
    missing_values: str,
    maximum_missing_fraction: float,
    missing_genes: str,
) -> list[str]:
    """Build the documented MAGMA command as an argument vector."""
    model_options = []
    for raw_option in model:
        option = str(raw_option).strip()
        if not option or option.startswith("-"):
            raise MagmaCovarError(
                "Each modules.magmacovar.model entry must be a MAGMA --model "
                "modifier without a leading dash"
            )
        model_options.append(option)
    gene_covar = [
        "--gene-covar",
        str(covariates_file),
        "missing-values=%s" % missing_values,
        "max-miss=%.12g" % maximum_missing_fraction,
    ]
    if missing_genes == "fill":
        gene_covar.append("missing-genes=fill")
    return [
        str(magma_bin),
        "--gene-results", str(gene_results_file),
        *gene_covar,
        "--model", *model_options, "direction-covar=%s" % direction,
        "--out", str(output_prefix),
    ]


def read_magma_covariate_results(
    results_file: str | Path, *, minimum_genes: int,
) -> list[dict]:
    """Return validated documented COVAR rows from MAGMA ``.gsa.out`` output."""
    path = require_nonempty_file(
        results_file, "MAGMA gene-property output", error_type=MagmaCovarError,
    )
    header: list[str] | None = None
    rows = []
    required = {"VARIABLE", "TYPE", "NGENES", "BETA", "BETA_STD", "SE", "P"}
    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for raw_line in handle:
                fields = raw_line.split()
                if not fields:
                    continue
                if header is None and required.issubset(fields):
                    header = fields
                    continue
                if header is None or len(fields) < len(header):
                    continue
                values = dict(zip(header, fields))
                if values.get("TYPE") != "COVAR":
                    continue
                try:
                    genes = int(values["NGENES"])
                    beta = float(values["BETA"])
                    beta_standardized = float(values["BETA_STD"])
                    standard_error = float(values["SE"])
                    p_value = float(values["P"])
                except (KeyError, ValueError) as exc:
                    raise MagmaCovarError(
                        "MAGMA gene-property output contains an invalid COVAR row"
                    ) from exc
                if genes < minimum_genes:
                    raise MagmaCovarError(
                        "MAGMA result %r used %d genes; at least %d are required"
                        % (values.get("VARIABLE", "unknown"), genes, minimum_genes)
                    )
                if (
                    not all(math.isfinite(value) for value in (
                        beta, beta_standardized, standard_error, p_value,
                    ))
                    or standard_error < 0
                    or not 0 <= p_value <= 1
                ):
                    raise MagmaCovarError(
                        "MAGMA result %r contains a non-finite or out-of-range statistic"
                        % values.get("VARIABLE", "unknown")
                    )
                rows.append(
                    {
                        "variable": values["VARIABLE"],
                        "genes": genes,
                        "beta": beta,
                        "beta_standardized": beta_standardized,
                        "standard_error": standard_error,
                        "p_value": p_value,
                    }
                )
    except UnicodeDecodeError as exc:
        raise MagmaCovarError(
            "MAGMA gene-property output is not valid UTF-8 text: %s" % path
        ) from exc
    except OSError as exc:
        raise MagmaCovarError(
            "Cannot read MAGMA gene-property output %s: %s" % (path, exc)
        ) from exc
    if header is None:
        raise MagmaCovarError(
            "MAGMA gene-property output is missing the documented result header"
        )
    if not rows:
        raise MagmaCovarError(
            "MAGMA gene-property output contains no COVAR results"
        )
    variables = [row["variable"] for row in rows]
    if len(variables) != len(set(variables)):
        raise MagmaCovarError(
            "MAGMA gene-property output contains duplicate COVAR variables"
        )
    return rows


def validate_magma_covariate_output(
    results_file: str | Path, *, minimum_genes: int,
) -> dict:
    """Validate and summarize COVAR rows in MAGMA's ``.gsa.out`` output."""
    rows = read_magma_covariate_results(
        results_file, minimum_genes=minimum_genes,
    )
    return {
        "tested_properties": len(rows),
        "minimum_result_genes": min(row["genes"] for row in rows),
        "maximum_result_genes": max(row["genes"] for row in rows),
    }


def run_magma_covariates(
    *,
    magma_bin: str,
    gene_results_file: str | Path,
    covariates_file: str | Path,
    output_prefix: str | Path,
    results_file: str | Path,
    native_log_file: str | Path,
    module,
    logger,
    timeout_seconds: float | None,
) -> dict:
    """Validate inputs, run MAGMA once, and validate its scientific output."""
    inputs = validate_magma_covariate_inputs(
        gene_results_file,
        covariates_file,
        minimum_genes=module.minimum_genes,
        maximum_missing_fraction=module.input.maximum_missing_fraction,
        missing_genes=module.input.missing_genes,
    )
    logger.record("OBSERVED", "magmacovar_inputs", **{
        key: value for key, value in inputs.items() if key != "property_missingness"
    })
    for missingness in inputs["property_missingness"]:
        logger.record("OBSERVED", "magmacovar_property_missingness", **missingness)

    command = build_magma_covariate_command(
        magma_bin,
        inputs["gene_results_file"],
        inputs["covariates_file"],
        output_prefix,
        model=module.model,
        direction=module.direction,
        missing_values=module.input.missing_values,
        maximum_missing_fraction=module.input.maximum_missing_fraction,
        missing_genes=module.input.missing_genes,
    )
    run_checked_command(
        command,
        "MAGMA gene-property analysis",
        logger=logger,
        error_type=MagmaCovarError,
        timeout_seconds=timeout_seconds,
        expected_outputs=[results_file, native_log_file],
    )
    output = validate_magma_covariate_output(
        results_file, minimum_genes=module.minimum_genes,
    )
    return {"inputs": inputs, "output": output, "command": command}


__all__ = [
    "build_magma_covariate_command",
    "read_magma_covariate_results",
    "run_magma_covariates",
    "validate_magma_covariate_inputs",
    "validate_magma_covariate_output",
    "validate_magma_covariate_table",
]
