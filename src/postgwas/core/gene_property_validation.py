"""Read-only MAGMA text-file contracts shared by gene-property consumers.

These checks preserve the native format and the caller's configured policies.
No rows, missing values, gene identifiers or hypothesis families are changed.
"""

from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.paths import require_nonempty_file


def _magma_gene_ids(path: Path, *, error_type) -> set[str]:
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
                    raise error_type(
                        "MAGMA .genes.raw row %d has fewer than seven required fields"
                        % line_number
                    )
                chromosome = fields[1].strip()
                try:
                    start = int(fields[2])
                    stop = int(fields[3])
                except ValueError as exc:
                    raise error_type(
                        "MAGMA .genes.raw row %d has invalid gene coordinates"
                        % line_number
                    ) from exc
                if not chromosome:
                    raise error_type(
                        "MAGMA .genes.raw row %d has an empty chromosome" % line_number
                    )
                if start < 0 or stop < start:
                    raise error_type(
                        "MAGMA .genes.raw row %d has an invalid start/stop interval"
                        % line_number
                    )
                gene = fields[0].strip()
                if not gene:
                    raise error_type(
                        "MAGMA .genes.raw contains an empty gene ID at line %d"
                        % line_number
                    )
                if gene in genes:
                    raise error_type(
                        "MAGMA .genes.raw contains duplicate gene ID %r" % gene
                    )
                genes.add(gene)
    except UnicodeDecodeError as exc:
        raise error_type(
            "MAGMA gene results are not valid UTF-8 text: %s" % path
        ) from exc
    except OSError as exc:
        raise error_type(
            "Cannot read MAGMA gene results %s: %s" % (path, exc)
        ) from exc
    if not version_found:
        raise error_type(
            "MAGMA gene results do not contain the required '# VERSION' metadata: %s"
            % path
        )
    if not genes:
        raise error_type(
            "MAGMA gene results contain no readable gene rows: %s" % path
        )
    return genes


def validate_magma_gene_results(
    gene_results_file: str | Path, *, minimum_genes: int, error_type=ValueError,
) -> tuple[Path, set[str]]:
    """Validate one native MAGMA ``.genes.raw`` file and its gene universe."""
    gene_results = require_nonempty_file(
        gene_results_file, "MAGMA .genes.raw file", error_type=error_type,
    )
    genes = validate_once(
        (gene_results,), {"validator": "magma_genes_raw"},
        lambda: _magma_gene_ids(gene_results, error_type=error_type),
        error_type=error_type,
    )
    if len(genes) < minimum_genes:
        raise error_type(
            "MAGMA gene results contain %d genes; at least %d are required"
            % (len(genes), minimum_genes)
        )
    record_file_validation(
        gene_results, "MAGMA gene results",
        checks=("native text format", "unique genes", "configured minimum gene count"),
        metrics={"genes": len(genes)},
    )
    return gene_results, set(genes)


def _inspect_magma_covariate_table(
    covariates_file: str | Path,
    *,
    eligible_gene_ids: set[str] | None = None,
    minimum_genes: int,
    maximum_missing_fraction: float,
    missing_genes: str,
    error_type,
) -> dict:
    """Validate the covariate table, optionally against MAGMA's gene universe."""
    covariates = require_nonempty_file(
        covariates_file, "MAGMA gene-covariate file", error_type=error_type,
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
                        raise error_type(
                            "MAGMA gene-covariate header must contain a gene-ID column "
                            "and at least one numeric property"
                        )
                    if len(set(header)) != len(header):
                        raise error_type(
                            "MAGMA gene-covariate header contains duplicate column names"
                        )
                    missing_by_property = [0] * (len(header) - 1)
                    observed_values = [set() for _ in header[1:]]
                    continue
                if len(fields) != len(header):
                    raise error_type(
                        "MAGMA gene-covariate row %d has %d fields; the header has %d"
                        % (line_number, len(fields), len(header))
                    )
                gene = fields[0].strip()
                if not gene:
                    raise error_type(
                        "MAGMA gene-covariate row %d has an empty gene ID" % line_number
                    )
                if gene in covariate_genes:
                    raise error_type(
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
                        raise error_type(
                            "MAGMA gene-covariate value at line %d, column %r is "
                            "neither numeric nor the required missing token NA"
                            % (line_number, header[index + 1])
                        ) from exc
                    if not math.isfinite(value):
                        raise error_type(
                            "MAGMA gene-covariate value at line %d, column %r is not finite"
                            % (line_number, header[index + 1])
                        )
                    if shared and len(observed_values[index]) < 2:
                        observed_values[index].add(value)
    except UnicodeDecodeError as exc:
        raise error_type(
            "MAGMA gene-covariate file is not valid UTF-8 text: %s" % covariates
        ) from exc
    except OSError as exc:
        raise error_type(
            "Cannot read MAGMA gene covariates %s: %s" % (covariates, exc)
        ) from exc

    if header is None or not covariate_genes:
        raise error_type(
            "MAGMA gene-covariate file has no data rows: %s" % covariates
        )
    if overlap < minimum_genes:
        if eligible_gene_ids is None:
            raise error_type(
                "MAGMA gene-covariate file contains %d genes; at least %d are required"
                % (overlap, minimum_genes)
            )
        raise error_type(
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
            denominator_description = (
                "all eligible .genes.raw genes"
                if missing_genes == "fill"
                else "gene IDs overlapping both input files"
            )
            policy_advice = (
                " If absent genes should be excluded instead, set "
                "modules.magmacovar.input.missing_genes to drop; for a "
                "magmacovar direct or pipeline command, use "
                "--covariate-missing-genes drop."
                if missing_genes == "fill"
                else ""
            )
            raise error_type(
                "Gene property %r has %d/%d missing values (%.3f) using "
                "missing-genes=%s with %s as the denominator, above the "
                "configured MAGMA maximum max-miss=%.3f. Correct or remove "
                "this property, or deliberately set "
                "modules.magmacovar.input.maximum_missing_fraction within 0 "
                "to 0.2; for a magmacovar direct or pipeline command, use "
                "--covariate-max-miss FRACTION. Changing the missing-values "
                "policy does not bypass max-miss.%s"
                % (
                    name,
                    missing,
                    denominator,
                    missing_fraction,
                    missing_genes,
                    denominator_description,
                    maximum_missing_fraction,
                    policy_advice,
                )
            )
        if len(observed_values[index]) < 2:
            raise error_type(
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


def validate_magma_covariate_table(
    covariates_file: str | Path, *, eligible_gene_ids: set[str] | None = None,
    minimum_genes: int, maximum_missing_fraction: float, missing_genes: str,
    error_type=ValueError,
) -> dict:
    """Reuse an unchanged table only for the same gene universe and policies.

    A later generated gene universe is a new compatibility check, not evidence
    covered by the startup structure-only validation. Store compact summaries,
    not a full gene-by-property matrix.
    """
    covariates = require_nonempty_file(
        covariates_file, "MAGMA gene-covariate file", error_type=error_type,
    )

    def inspect():
        result = _inspect_magma_covariate_table(
            covariates, eligible_gene_ids=eligible_gene_ids,
            minimum_genes=minimum_genes,
            maximum_missing_fraction=maximum_missing_fraction,
            missing_genes=missing_genes, error_type=error_type,
        )
        record_file_validation(
            covariates, "MAGMA gene covariates",
            checks=("all-row field counts", "unique genes", "finite numeric or native NA values",
                    "configured missingness and nonconstant properties"),
            metrics={key: value for key, value in result.items()
                     if key not in {"covariates_file", "property_names", "property_missingness"}},
            message=("Generated gene-universe compatibility checked."
                     if eligible_gene_ids is not None
                     else "Compatibility with generated gene results remains deferred."),
        )
        return result

    return deepcopy(validate_once((covariates,), {
        "validator": "magma_covariate_table",
        "eligible_gene_ids": sorted(eligible_gene_ids) if eligible_gene_ids is not None else None,
        "minimum_genes": minimum_genes,
        "maximum_missing_fraction": maximum_missing_fraction,
        "missing_genes": missing_genes,
    }, inspect, error_type=error_type))
