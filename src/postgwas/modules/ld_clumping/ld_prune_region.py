"""Prune one strongest association per configured population LD block."""

from __future__ import annotations

import gzip
import math
from pathlib import Path
import shutil
import subprocess

import polars as pl

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.core.vcf import extract_vcf_table
from postgwas.modules.ld_clumping.common import (
    LDRegionClumpingError,
    canonical_variant_expression,
    chromosome_expression,
    exclude_mhc,
)


# Stable audit schema. The source VCF expressions that populate these internal
# fields remain configuration-driven through ``vcf_fields``.
REGION_OUTSIDE_LD_COLUMNS = [
    "CHR",
    "BP",
    "SNP",
    "uniq_id",
    "ID",
    "REF",
    "ALT",
    "BETA",
    "SE",
    "AF",
    "LP",
    "P_value",
    "LDblock",
    "exclusion_reason",
]


def _write_compressed_table(frame, destination, *, delimiter, compressor, threads):
    """Write one gzip table, preferring the configured parallel compressor.

    Python's :mod:`gzip` is single-threaded and defaults to maximum compression,
    which dominates the runtime of this step on genome-wide tables.
    """
    executable = shutil.which(compressor) if compressor else None
    if executable is not None:
        try:
            with destination.open("wb") as handle:
                process = subprocess.Popen(
                    [executable, "-c", "-p", str(max(1, int(threads or 1)))],
                    stdin=subprocess.PIPE,
                    stdout=handle,
                )
                frame.write_csv(process.stdin, separator=delimiter)
                process.stdin.close()
                if process.wait() == 0:
                    return "%s -p %s" % (
                        Path(executable).name,
                        max(1, int(threads or 1)),
                    )
        except (OSError, ValueError):
            pass
    # Binary mode, not text. polars writes bytes through its own writer, and
    # handing it a text-mode gzip wrapper truncates the deflate stream and
    # appends the CSV uncompressed, producing a .gz that cannot be read back.
    with gzip.open(destination, "wb", compresslevel=6) as handle:
        frame.write_csv(handle, separator=delimiter)
    return "gzip"


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
    include_report_tables: bool = False,
    configuration=None,
    logger=None,
):
    """Extract configured fields and select the strongest SNP in each LD block."""
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
        "ID": fields.variant_id,
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
            # VCF projections contain identifiers such as X/Y/MT. Preserve all
            # queried tokens as text, then convert only validated numeric fields.
            schema_overrides={column: pl.String for column in columns},
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

    table, mhc_statistics = exclude_mhc(
        table, configuration, chromosome_column="CHR", position_column="BP",
    )
    if logger is not None:
        logger.record("ACTION", "region_mhc_exclusion", **mhc_statistics)

    reference_allele = pl.col("REF").cast(pl.Utf8).str.to_uppercase()
    alternate_allele = pl.col("ALT").cast(pl.Utf8).str.to_uppercase()
    enriched = table.with_columns(
        pl.concat_str(
            [
                chromosome_expression("CHR"),
                pl.col("BP").cast(pl.Utf8),
                reference_allele,
                alternate_allele,
            ],
            separator="_",
        ).alias("SNP"),
        # Allele-order-independent key, identical to the standard method's
        # uniq_id, so region and standard output can be joined directly.
        canonical_variant_expression("CHR", "BP", "REF", "ALT").alias(
            "uniq_id"
        ),
        pl.col(ld_field).alias("LDblock"),
        # P is presentation-only. LP remains the scientific threshold and
        # ranking value because very small P values underflow to zero.
        (10.0 ** (-pl.col("LP"))).alias("P_value"),
    )
    annotated = enriched.filter(pl.col("LDblock").is_not_null())
    if table.height and annotated.is_empty():
        raise LDRegionClumpingError(
            "No variants contain %s annotations. Check that LD annotation used "
            "population %s and genome build %s."
            % (ld_field, population, configuration.genome_build.value)
        )
    lead_lp_threshold = -math.log10(configuration.lead_pvalue)
    significant_outside_ld_regions = (
        enriched.filter(
            pl.col("LDblock").is_null()
            & (pl.col("LP") >= lead_lp_threshold)
        )
        .with_columns(
            pl.lit("missing_%s_annotation" % ld_field).alias(
                "exclusion_reason"
            )
        )
        .select(REGION_OUTSIDE_LD_COLUMNS)
        .sort(["CHR", "BP"])
    )
    prepared = annotated.with_columns(
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
    # Rank on LP, not on P. Sorting the underflowed P column would tie every
    # variant stronger than ~1e-308 at 0.0 and break the tie on the SNP string,
    # which selects a block lead by allele spelling rather than by association.
    pruned = (
        prepared.sort(["LP", "SNP"], descending=[True, False])
        .group_by("LDblock", maintain_order=True)
        .head(1)
    )
    significant = pruned.filter(pl.col("LP") >= lead_lp_threshold)

    try:
        compressor = _write_compressed_table(
            prepared,
            paths["region_raw_table"],
            delimiter=configuration.table.delimiter,
            compressor=configuration.table.compressor,
            threads=threads,
        )
        pruned.write_csv(
            paths["region_pruned_table"], separator=configuration.table.delimiter,
        )
        significant.write_csv(
            paths["region_significant_table"],
            separator=configuration.table.delimiter,
        )
        significant_outside_ld_regions.write_csv(
            paths["region_significant_outside_ld_regions"],
            separator=configuration.table.delimiter,
        )
        with paths["region_log"].open("w", encoding="utf-8") as handle:
            handle.write("Annotated-region LD clumping completed\n")
            handle.write("input_vcf=%s\n" % sumstat_vcf)
            handle.write("population=%s\n" % population)
            handle.write("genome_build=%s\n" % configuration.genome_build.value)
            handle.write("lead_pvalue=%s\n" % configuration.lead_pvalue)
            handle.write(
                "mhc_exclusion=%s rows_in=%d rows_out=%d removed=%d\n"
                % (
                    mhc_statistics["enabled"],
                    mhc_statistics["rows_in"],
                    mhc_statistics["rows_out"],
                    mhc_statistics["removed"],
                )
            )
            handle.write("annotated_variants=%d\n" % prepared.height)
            handle.write("ld_blocks=%d\n" % pruned.height)
            handle.write("significant_blocks=%d\n" % significant.height)
            handle.write(
                "genome_wide_significant_outside_ld_regions=%d\n"
                % significant_outside_ld_regions.height
            )
            handle.write(
                "outside_ld_region_definition=LP >= %s and %s is missing, "
                "after configured MHC exclusion\n"
                % (lead_lp_threshold, ld_field)
            )
            for row in significant_outside_ld_regions.iter_rows(named=True):
                handle.write(
                    "outside_ld_region_variant"
                    " | chromosome=%s | position=%s | canonical_id=%s"
                    " | input_variant_id=%s | ref=%s | alt=%s"
                    " | beta=%s | se=%s | af=%s | lp=%s | p_value=%s"
                    " | reason=%s\n"
                    % (
                        row["CHR"],
                        row["BP"],
                        row["uniq_id"],
                        row["ID"] or "none",
                        row["REF"],
                        row["ALT"],
                        row["BETA"],
                        row["SE"],
                        row["AF"],
                        row["LP"],
                        row["P_value"],
                        row["exclusion_reason"],
                    )
                )
            handle.write("block_lead_selected_by=LP descending\n")
            handle.write("raw_table_compressor=%s\n" % compressor)
            handle.write("allele_orientation=REF and ALT retained unchanged\n")
    except OSError as exc:
        raise LDRegionClumpingError(
            "Cannot write region-clumping output: %s" % exc
        ) from exc

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
            genome_wide_significant_outside_ld_regions=(
                significant_outside_ld_regions.height
            ),
            lead_pvalue=configuration.lead_pvalue,
            significant_outside_ld_regions_file=str(
                paths["region_significant_outside_ld_regions"]
            ),
            allele_orientation="REF and ALT retained exactly as queried",
            pruned_file=str(paths["region_pruned_table"]),
            significant_file=str(paths["region_significant_table"]),
        )
    result = {
        "status": "completed",
        "ldpruned_sig_file": str(paths["region_significant_table"]),
        "ldpruned_file": str(paths["region_pruned_table"]),
        "raw_file": str(paths["region_raw_table"]),
        "log_file": str(paths["region_log"]),
        "annotated_variants": prepared.height,
        "ld_blocks": pruned.height,
        "significant_blocks": significant.height,
        "genome_wide_significant_outside_ld_regions": (
            significant_outside_ld_regions.height
        ),
        "significant_outside_ld_regions_file": str(
            paths["region_significant_outside_ld_regions"]
        ),
        "output_files": {
            "raw_table": str(paths["region_raw_table"]),
            "pruned_table": str(paths["region_pruned_table"]),
            "significant_table": str(paths["region_significant_table"]),
            "significant_outside_ld_regions": str(
                paths["region_significant_outside_ld_regions"]
            ),
            "detailed_log": str(paths["region_log"]),
        },
    }
    if include_report_tables:
        # Presentation rows are transient. The service consumes them while
        # writing HTML, then removes them before return/checkpoint logging.
        result["_significant_outside_ld_region_details"] = {
            "columns": list(REGION_OUTSIDE_LD_COLUMNS),
            "rows": significant_outside_ld_regions.rows(),
            "total_rows": significant_outside_ld_regions.height,
            "output_file": str(paths["region_significant_outside_ld_regions"]),
        }
    return result


__all__ = ["LDRegionClumpingError", "ld_clump_by_regions"]
