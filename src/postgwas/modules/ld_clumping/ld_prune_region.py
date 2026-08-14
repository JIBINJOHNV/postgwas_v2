"""Prune one strongest association per configured population LD block."""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.core.vcf import extract_vcf_table


class LDRegionClumpingError(RuntimeError):
    """Annotated-region clumping cannot produce scientifically valid output."""


def _normalised_chromosome_expression(column: str) -> pl.Expr:
    return (
        pl.col(column)
        .cast(pl.Utf8)
        .str.to_uppercase()
        .str.replace(r"^CHR", "")
    )


def _region_output_paths(output_directory, dataset_id, population, configuration):
    values = {"dataset_id": dataset_id, "population": population}
    return {
        name: configured_output_path(
            output_directory,
            pattern,
            error_type=LDRegionClumpingError,
            **values,
        )
        for name, pattern in configuration.output_layout.model_dump().items()
        if name.startswith("region_")
    }


def ld_clump_by_regions(
    sumstat_vcf: str,
    output_directory: str,
    dataset_id: str,
    population: str | None = None,
    bcftools: str | None = None,
    threads: int | None = None,
    *,
    configuration=None,
    logger=None,
):
    """Extract configured fields and select the lowest-P SNP in each LD block."""
    del threads  # bcftools query is a single streaming operation.
    application = None
    if configuration is None or population is None or bcftools is None:
        application = load_configuration()
    configuration = configuration or application.modules.ld_clumping
    population = population or configuration.population.value
    bcftools = bcftools or application.resources.executables.bcftools
    population = str(population).upper()
    if population not in {value.value for value in type(configuration.population)}:
        raise LDRegionClumpingError(
            "Population %s is not a supported canonical population" % population
        )
    output_root = Path(output_directory).expanduser().resolve()
    paths = _region_output_paths(
        output_root, dataset_id, population, configuration,
    )
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    fields = configuration.vcf_fields
    ld_field = "%s_LDblock" % population
    columns = {
        "CHR": fields.chromosome,
        "BP": fields.position,
        "REF": fields.reference_allele,
        "ALT": fields.alternate_allele,
        "BETA": fields.effect,
        "SE": fields.standard_error,
        "AF": fields.allele_frequency,
        "LP": fields.log_pvalue,
        ld_field: fields.ld_block.format(population=population),
    }
    try:
        extract_vcf_table(
            sumstat_vcf,
            paths["region_working_table"],
            dataset_id,
            columns,
            bcftools,
            delimiter=configuration.table.delimiter,
            io_buffer_bytes=configuration.table.io_buffer_bytes,
            include_expression=configuration.table.biallelic_include_expression,
            logger=logger,
            error_type=LDRegionClumpingError,
            purpose="Extracting annotated LD-block variants",
        )
        table = pl.read_csv(
            paths["region_working_table"],
            separator=configuration.table.delimiter,
            null_values=configuration.table.null_values,
            infer_schema_length=configuration.table.infer_schema_length,
        ).with_columns(
            pl.col("BP").cast(pl.Int64, strict=False),
            pl.col("LP").cast(pl.Float64, strict=False),
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise LDRegionClumpingError(
            "Cannot read extracted region-clumping table: %s" % exc
        ) from exc
    invalid = table.filter(
        pl.col("BP").is_null()
        | pl.col("LP").is_null()
        | pl.col("LP").is_nan()
        | (pl.col("LP") < 0)
    )
    if not invalid.is_empty():
        raise LDRegionClumpingError(
            "Region-clumping input contains %d invalid position or LP values"
            % invalid.height
        )

    before_mhc = table.height
    if configuration.remove_mhc:
        mhc = configuration.mhc_regions[configuration.genome_build]
        mhc_chromosome = str(mhc.chromosome).upper().removeprefix("CHR")
        normalised_chromosome = _normalised_chromosome_expression("CHR")
        table = table.filter(
            ~(
                (normalised_chromosome == mhc_chromosome)
                & pl.col("BP").is_between(mhc.start, mhc.end, closed="both")
            )
        )
    if logger is not None:
        logger.record(
            "ACTION",
            "region_mhc_exclusion",
            enabled=configuration.remove_mhc,
            rows_in=before_mhc,
            rows_out=table.height,
            removed=before_mhc - table.height,
            genome_build=configuration.genome_build.value,
        )

    annotated = table.filter(pl.col(ld_field).is_not_null())
    if table.height and annotated.is_empty():
        raise LDRegionClumpingError(
            "No variants contain %s annotations. Check that LD annotation used "
            "population %s and genome build %s."
            % (ld_field, population, configuration.genome_build.value)
        )
    prepared = annotated.with_columns(
        pl.concat_str(
            [
                _normalised_chromosome_expression("CHR"),
                pl.col("BP").cast(pl.Utf8),
                pl.col("REF").cast(pl.Utf8).str.to_uppercase(),
                pl.col("ALT").cast(pl.Utf8).str.to_uppercase(),
            ],
            separator="_",
        ).alias("SNP"),
        pl.col(ld_field).alias("LDblock"),
        (10.0 ** (-pl.col("LP"))).alias("P_value"),
        pl.col(ld_field).str.extract(r"_(\d+)_(\d+)$", 1).alias("START"),
        pl.col(ld_field).str.extract(r"_(\d+)_(\d+)$", 2).alias("END"),
    ).with_columns(
        pl.col("START").cast(pl.Int64, strict=False),
        pl.col("END").cast(pl.Int64, strict=False),
    )
    invalid_blocks = prepared.filter(
        pl.col("START").is_null()
        | pl.col("END").is_null()
        | (pl.col("END") <= pl.col("START"))
    )
    if not invalid_blocks.is_empty():
        examples = invalid_blocks["LDblock"].unique().head(5).to_list()
        raise LDRegionClumpingError(
            "Found %d variants with malformed LD-block coordinates; examples=%s"
            % (invalid_blocks.height, examples)
        )
    pruned = (
        prepared.sort(["P_value", "SNP"])
        .group_by("LDblock", maintain_order=True)
        .head(1)
    )
    significant = pruned.filter(pl.col("P_value") <= configuration.lead_pvalue)

    try:
        with gzip.open(paths["region_raw_table"], "wt", encoding="utf-8") as handle:
            prepared.write_csv(handle, separator=configuration.table.delimiter)
        pruned.write_csv(
            paths["region_pruned_table"], separator=configuration.table.delimiter,
        )
        significant.write_csv(
            paths["region_significant_table"],
            separator=configuration.table.delimiter,
        )
        with paths["region_log"].open("w", encoding="utf-8") as handle:
            handle.write("Annotated-region LD clumping completed\n")
            handle.write("input_vcf=%s\n" % sumstat_vcf)
            handle.write("population=%s\n" % population)
            handle.write("genome_build=%s\n" % configuration.genome_build.value)
            handle.write("lead_pvalue=%s\n" % configuration.lead_pvalue)
            handle.write("annotated_variants=%d\n" % prepared.height)
            handle.write("ld_blocks=%d\n" % pruned.height)
            handle.write("significant_blocks=%d\n" % significant.height)
            handle.write("allele_orientation=REF and ALT retained unchanged\n")
    except OSError as exc:
        raise LDRegionClumpingError("Cannot write region-clumping output: %s" % exc) from exc

    paths["region_working_table"].unlink(missing_ok=True)
    try:
        paths["region_working_table"].parent.rmdir()
    except OSError:
        pass
    if logger is not None:
        logger.record(
            "OUTPUT",
            "region_clumping",
            annotated_variants=prepared.height,
            ld_blocks=pruned.height,
            significant_blocks=significant.height,
            lead_pvalue=configuration.lead_pvalue,
            allele_orientation="REF and ALT retained exactly as queried",
            pruned_file=str(paths["region_pruned_table"]),
            significant_file=str(paths["region_significant_table"]),
        )
    return {
        "status": "completed",
        "ldpruned_sig_file": str(paths["region_significant_table"]),
        "ldpruned_file": str(paths["region_pruned_table"]),
        "raw_file": str(paths["region_raw_table"]),
        "log_file": str(paths["region_log"]),
        "annotated_variants": prepared.height,
        "ld_blocks": pruned.height,
        "significant_blocks": significant.height,
    }


__all__ = ["LDRegionClumpingError", "ld_clump_by_regions"]
