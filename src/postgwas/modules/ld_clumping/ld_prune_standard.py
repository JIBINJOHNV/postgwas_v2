
import stat
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.core.ui import screen_field, screen_line
from postgwas.core.vcf import extract_vcf_table
from postgwas.modules.ld_clumping.common import (  # noqa: F401 - re-exported API
    PipelineStageError,
    canonical_variant_expression,
    canonical_variant_id,
    chromosome_expression,
    chromosome_sort_key,
    exclude_mhc,
    format_compact_reference_exclusion_counts,
    function_error,
    normalise_chromosome,
    reference_exclusion_analysis_limitation,
    reference_exclusion_reason_counts,
    reference_exclusion_reason_items,
    reference_exclusion_reason_label,
)


_STANDARD_FORMATTED_SCHEMA = {
    # Variant identity is textual by VCF contract; numeric values have already
    # passed explicit conversion validation before this table is published.
    "chrcol": pl.String,
    "poscol": pl.Int64,
    "neacol": pl.String,
    "eacol": pl.String,
    "rsIDcol": pl.String,
    "pcol": pl.Float64,
    # -log10(P) exactly as published by the VCF. P underflows to 0.0 below
    # ~1e-308, so this column, not pcol, orders variants by strength.
    "lpcol": pl.Float64,
    "becol": pl.Float64,
    "secol": pl.Float64,
    "eafcol": pl.Float64,
}


# Scientific settings decide the result and must come from the same module
# configuration the caller resolved. Environment settings locate tools and size
# the run, so filling them from the global configuration cannot change science.
_SCIENTIFIC_OPTION_SOURCES = {
    "candidate_p_threshold": ("module", "candidate_pvalue"),
    "candidate_p": ("module", "candidate_pvalue"),
    "window_kb": ("module", "window_kb"),
    "missing_index_action": ("module", "missing_index_action"),
    "summary_pvalue_thresholds": ("module", "summary_pvalue_thresholds"),
    "minimum_reference_maf": ("module", "minimum_reference_maf"),
    "reference_file_pattern": ("reference", "file_pattern"),
    "reference_reverse_file_pattern": ("reference", "reverse_file_pattern"),
    "reference_inventory_pattern": ("reference", "variant_inventory_pattern"),
    "reference_orientation": ("reference", "orientation"),
    "reference_index_suffix": ("reference", "index_suffix"),
}

_ENVIRONMENT_OPTION_SOURCES = {
    "tabix_bin": ("executables", "tabix"),
    "bcftools_bin": ("executables", "bcftools"),
    "threads": ("execution", "threads"),
    "memory_gb": ("execution", "memory_gb"),
    "screen_label_width": ("logging", "terminal_label_width"),
}

_RUNTIME_OPTION_SOURCES = {
    **_SCIENTIFIC_OPTION_SOURCES,
    **_ENVIRONMENT_OPTION_SOURCES,
}


def _resolve_standard_options(*, application=None, configuration=None, **provided):
    """Fill omitted runtime options from one configuration source.

    Option names are validated on every call, not only when something is
    missing. A caller that already resolved a module configuration may not have
    a *scientific* value quietly supplied from the packaged defaults, because
    that would mix two configuration sources within one analysis.
    """
    unknown = sorted(set(provided).difference(_RUNTIME_OPTION_SOURCES))
    if unknown:
        raise ConfigurationError(
            "Unknown LD-clumping runtime option(s): %s" % ", ".join(unknown)
        )
    missing = sorted(name for name, value in provided.items() if value is None)
    if not missing:
        return provided
    conflicting = [
        name
        for name in missing
        if application is None
        and configuration is not None
        and name in _SCIENTIFIC_OPTION_SOURCES
    ]
    if conflicting:
        raise ConfigurationError(
            "LD-clumping received a resolved module configuration but no value "
            "for %s. Pass the value explicitly rather than letting a scientific "
            "setting fall back to the packaged defaults." % ", ".join(conflicting)
        )
    if configuration is not None:
        # Environment-only fallback: keep the caller's module configuration.
        sources = {
            "module": configuration,
            "reference": configuration.reference,
            "compute": configuration.compute,
        }
        application = load_configuration()
        sources.update(
            executables=application.resources.executables,
            execution=application.execution,
            logging=application.logging,
        )
        return {
            **provided,
            **{
                name: getattr(sources[_RUNTIME_OPTION_SOURCES[name][0]],
                              _RUNTIME_OPTION_SOURCES[name][1])
                for name in missing
            },
        }
    application = application or load_configuration()
    module = application.modules.ld_clumping
    sources = {
        "module": module,
        "reference": module.reference,
        "compute": module.compute,
        "executables": application.resources.executables,
        "execution": application.execution,
        "logging": application.logging,
    }
    resolved = dict(provided)
    for name in missing:
        source, attribute = _RUNTIME_OPTION_SOURCES[name]
        resolved[name] = getattr(sources[source], attribute)
    return resolved


def add_canonical_ids(gwas, log=None):
    """Add canonical IDs and retain the strongest row for true duplicates."""
    required = {"chrcol", "poscol", "neacol", "eacol", "rsIDcol", "pcol"}
    missing = sorted(required.difference(gwas.columns))
    if missing:
        raise PipelineStageError(
            "03 canonical ID preparation",
            "add_canonical_ids",
            f"Missing required columns: {', '.join(missing)}",
        )
    emit = log or (lambda message: None)
    # Rank duplicates on LP when the lossless column is present, because P
    # underflows to 0.0 below ~1e-308 and would tie the strongest variants.
    has_lp = "lpcol" in gwas.columns
    strength = (
        pl.col("lpcol").cast(pl.Float64, strict=False)
        if has_lp
        else -pl.col("pcol").cast(pl.Float64, strict=False)
    )

    chrom = chromosome_expression("chrcol")
    ea = pl.col("eacol").cast(pl.Utf8).str.to_uppercase()
    nea = pl.col("neacol").cast(pl.Utf8).str.to_uppercase()
    gwas = gwas.with_row_index("_input_order").with_columns(
        chrom.alias("chrcol"),
        pl.col("rsIDcol").cast(pl.Utf8).alias("input_id"),
        pl.col("pcol").cast(pl.Float64, strict=False).alias("pcol"),
        ea.alias("_effect_allele_orientation"),
        nea.alias("_non_effect_allele_orientation"),
        strength.alias("_strength"),
        canonical_variant_expression(
            "chrcol", "poscol", "neacol", "eacol"
        ).alias("uniq_id"),
    )
    invalid_p = gwas.filter(pl.col("pcol").is_null() | pl.col("pcol").is_nan())
    if not invalid_p.is_empty():
        examples = invalid_p.select("chrcol", "poscol", "rsIDcol").head(3)
        raise PipelineStageError(
            "03 canonical ID preparation",
            "add_canonical_ids",
            f"Found {invalid_p.height} rows with missing or non-numeric P values; "
            f"examples={examples.to_dicts()}",
        )
    invalid = gwas.filter(
        pl.col("uniq_id").is_null()
        | pl.col("eacol").cast(pl.Utf8).is_in(["", ".", "NA"])
        | pl.col("neacol").cast(pl.Utf8).is_in(["", ".", "NA"])
    )
    if not invalid.is_empty():
        examples = invalid.select("chrcol", "poscol", "eacol", "neacol").head(3)
        raise PipelineStageError(
            "03 canonical ID preparation",
            "add_canonical_ids",
            f"Cannot create canonical IDs for {invalid.height} rows; "
            f"examples={examples.to_dicts()}",
        )

    duplicates = gwas.group_by("uniq_id").len().filter(pl.col("len") > 1)
    if not duplicates.is_empty():
        orientation_conflicts = (
            gwas.group_by("uniq_id")
            .agg(
                pl.struct(
                    "_effect_allele_orientation",
                    "_non_effect_allele_orientation",
                ).n_unique().alias("orientation_count")
            )
            .filter(pl.col("orientation_count") > 1)
        )
        if not orientation_conflicts.is_empty():
            examples = orientation_conflicts["uniq_id"].head(5).to_list()
            raise PipelineStageError(
                "03 canonical ID preparation",
                "add_canonical_ids",
                "Duplicate canonical variants have conflicting effect-allele "
                "orientation; automatic allele flipping is prohibited "
                f"| conflicting_variants={orientation_conflicts.height} "
                f"| examples={examples}",
            )
        examples = duplicates["uniq_id"].head(5).to_list()
        rows_removed = duplicates.select((pl.col("len") - 1).sum()).item()
        emit(
            "[WARNING] [STAGE: 03 canonical ID preparation] "
            "[FUNCTION: add_canonical_ids] FUMA-style duplicate resolution "
            "retained the strongest row for each canonical ID "
            f"| ranked_by={'LP' if has_lp else 'P'} "
            f"| duplicate_canonical_ids={duplicates.height} "
            f"| rows_removed={rows_removed} | examples={examples}"
        )
    return (
        gwas.sort(
            ["chrcol", "poscol", "_strength", "_input_order"],
            descending=[False, False, True, False],
        )
        .unique("uniq_id", keep="first", maintain_order=True)
        .drop(
            "_input_order",
            "_effect_allele_orientation",
            "_non_effect_allele_orientation",
            "_strength",
        )
    )


def _standard_conversion_paths(output_folder, sample_name, configuration=None):
    configuration = configuration or load_configuration().modules.ld_clumping
    output_dir = Path(output_folder).expanduser().resolve()
    values = {
        "dataset_id": sample_name,
        "population": configuration.population.value,
    }
    return (
        configured_output_path(
            output_dir,
            configuration.output_layout.standard_formatted_table,
            error_type=RuntimeError,
            **values,
        ),
        configured_output_path(
            output_dir,
            configuration.output_layout.standard_log,
            error_type=RuntimeError,
            **values,
        ),
    )


def vcf_to_standard_ldclump(
    sumstat_vcf: str,
    output_folder: str,
    sample_name: str,
    bcftools_path=None,
    threads=None,
    *,
    configuration=None,
    logger=None,
):
    """Extract and validate the configured GWAS-VCF projection without a shell."""
    del threads
    application = None
    if configuration is None or bcftools_path is None:
        application = load_configuration()
    configuration = configuration or application.modules.ld_clumping
    bcftools_path = bcftools_path or application.resources.executables.bcftools
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder).expanduser().resolve()
    output_file, log_file = _standard_conversion_paths(
        output_dir, sample_name, configuration,
    )
    working_file = configured_output_path(
        output_dir,
        configuration.output_layout.standard_working_table,
        error_type=RuntimeError,
        dataset_id=sample_name,
        population=configuration.population.value,
    )
    try:
        for path in (output_file, log_file, working_file):
            path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise function_error(
            "01 VCF conversion",
            "vcf_to_standard_ldclump",
            error,
            sample=sample_name,
            input_vcf=vcf_path,
            output_folder=output_dir,
        ) from error
    fields = configuration.vcf_fields
    columns = {
        "chrcol": fields.chromosome,
        "poscol": fields.position,
        "neacol": fields.reference_allele,
        "eacol": fields.alternate_allele,
        "rsIDcol": fields.variant_id,
        "lpcol": fields.log_pvalue,
        "becol": fields.effect,
        "secol": fields.standard_error,
        "eafcol": fields.allele_frequency,
    }
    try:
        extract_vcf_table(
            vcf_path,
            working_file,
            sample_name,
            columns,
            bcftools_path,
            delimiter=configuration.table.delimiter,
            io_buffer_bytes=configuration.table.io_buffer_bytes,
            include_expression=configuration.table.biallelic_include_expression,
            logger=logger,
            error_type=RuntimeError,
            purpose="Extracting the standard LD-clumping table",
        )
        extracted = pl.read_csv(
            working_file,
            separator=configuration.table.delimiter,
            null_values=configuration.table.null_values,
            infer_schema_length=configuration.table.infer_schema_length,
            # Preserve the lossless VCF projection first; numeric fields are
            # converted explicitly so invalid values remain actionable errors.
            schema_overrides={column: pl.String for column in columns},
        ).with_columns(
            pl.col("poscol").cast(pl.Int64, strict=False),
            pl.col("lpcol").cast(pl.Float64, strict=False),
            pl.col("becol").cast(pl.Float64, strict=False),
            pl.col("secol").cast(pl.Float64, strict=False),
            pl.col("eafcol").cast(pl.Float64, strict=False),
        )
    except (OSError, RuntimeError, pl.exceptions.PolarsError) as error:
        raise function_error(
            "01 VCF conversion",
            "vcf_to_standard_ldclump",
            error,
            sample=sample_name,
            input_vcf=vcf_path,
            output_file=output_file,
            log_file=log_file,
        ) from error
    invalid = extracted.filter(
        pl.col("poscol").is_null()
        | pl.col("lpcol").is_null()
        | pl.col("lpcol").is_nan()
        | (pl.col("lpcol") < 0)
    )
    if not invalid.is_empty():
        raise PipelineStageError(
            "01 VCF conversion",
            "vcf_to_standard_ldclump",
            "Extracted VCF contains %d invalid position or LP values"
            % invalid.height,
            sample=sample_name,
            input_vcf=vcf_path,
        )
    # lpcol is retained alongside pcol: P underflows to 0.0 below ~1e-308 and
    # cannot rank the strongest variants, which is exactly where lead SNPs are.
    formatted = extracted.with_columns(
        (10.0 ** (-pl.col("lpcol"))).alias("pcol")
    ).select(*_STANDARD_FORMATTED_SCHEMA)
    try:
        formatted.write_csv(output_file, separator=configuration.table.delimiter)
        with log_file.open("w", encoding="utf-8") as handle:
            handle.write("Standard LD-clumping VCF extraction completed\n")
            handle.write("input_vcf=%s\n" % vcf_path)
            handle.write("output_table=%s\n" % output_file)
            handle.write("variants=%d\n" % formatted.height)
            handle.write("allele_orientation=REF to neacol; ALT to eacol; unchanged\n")
    except OSError as error:
        raise function_error(
            "01 VCF conversion", "vcf_to_standard_ldclump", error,
            sample=sample_name, output_file=output_file, log_file=log_file,
        ) from error
    working_file.unlink(missing_ok=True)
    try:
        working_file.parent.rmdir()
    except OSError:
        pass
    # The frame is returned as well as written so the caller does not have to
    # serialise and re-parse a genome-wide table it already holds in memory.
    return str(output_file), formatted


LD_REFERENCE_COLUMNS = [
    "chr_a", "pos_a", "snp_a", "chr_b", "pos_b", "snp_b", "r2",
]

LD_VARIANT_INVENTORY_COLUMNS = [
    "chromosome",
    "position",
    "reference_id",
    "canonical_id",
    "allele_1",
    "allele_2",
    "minor_allele_frequency",
]

LD_REFERENCE_EXCLUSION_INPUT_SCHEMA = {
    "chromosome": pl.String,
    "position": pl.Int64,
    "canonical_id": pl.String,
    "input_variant_id": pl.String,
    "reason": pl.String,
    "observed_reference_ids": pl.String,
    "action": pl.String,
}

LD_REFERENCE_EXCLUSION_LOCUS_SCHEMA = {
    "reported_locus_boundary_status": pl.String,
    "overlapping_genomic_loci": pl.String,
    "overlapping_locus_chromosome": pl.String,
    "overlapping_locus_start": pl.Int64,
    "overlapping_locus_end": pl.Int64,
    "locus_membership_interpretation": pl.String,
}

LD_REFERENCE_EXCLUSION_SCHEMA = {
    **LD_REFERENCE_EXCLUSION_INPUT_SCHEMA,
    **LD_REFERENCE_EXCLUSION_LOCUS_SCHEMA,
}
LD_REFERENCE_EXCLUSION_INPUT_COLUMNS = list(
    LD_REFERENCE_EXCLUSION_INPUT_SCHEMA
)
LD_REFERENCE_EXCLUSION_COLUMNS = list(LD_REFERENCE_EXCLUSION_SCHEMA)

EXCLUSION_INSIDE_LOCUS_BOUNDARY = "inside_reported_locus_boundary"
EXCLUSION_OUTSIDE_LOCUS_BOUNDARIES = (
    "outside_all_reported_locus_boundaries"
)


@dataclass(frozen=True)
class ChromosomeReferenceData:
    """Validated, bidirectional LD adjacency for one chromosome.

    Index membership is established from the variant inventory, never inferred
    from the presence of a pair row. This preserves FUMA's distinction between
    an unsupported index and a supported index with no qualifying partners.
    """

    partners: dict[str, pl.DataFrame]
    verified_indexes: frozenset[str]
    exclusions: tuple[dict[str, object], ...]
    ld_rows: int
    ld_rows_before_maf: int = 0
    low_maf_partner_rows_excluded: int = 0
    low_maf_partner_variants_excluded: int = 0


def _empty_ld_partners() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "uniq_id": pl.String,
            "ld_id": pl.String,
            "ref_chr": pl.String,
            "ref_pos": pl.Int64,
            "r2": pl.Float64,
        }
    )


def reference_contigs(tabix_bin, ld_file_path):
    """List the contig labels the LD reference index actually contains."""
    try:
        result = subprocess.run(
            [tabix_bin, "--list-chroms", str(ld_file_path)],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError as error:
        raise PipelineStageError(
            "04 LD reference validation",
            "reference_contigs",
            "tabix executable was not found",
            ld_file=str(ld_file_path),
        ) from error
    except subprocess.CalledProcessError as error:
        stderr = (error.stderr or "").strip()
        raise PipelineStageError(
            "04 LD reference validation",
            "reference_contigs",
            "tabix could not list the reference contigs: %s"
            % (stderr or "exit code %s" % error.returncode),
            ld_file=str(ld_file_path),
        ) from error
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def resolve_reference_contig(chrom, ld_file_path, *, tabix_bin, contig_cache=None):
    """Return the reference's own label for ``chrom``.

    The summary statistics are normalised to bare labels while a reference may
    be built with ``chr``-prefixed contigs. Querying the wrong label returns no
    rows for every variant. Resolving the label from the index makes both
    conventions work and turns a genuine naming mismatch into one explicit
    reference-contract error.
    """
    key = str(ld_file_path)
    cache = contig_cache if contig_cache is not None else {}
    if key not in cache:
        contigs = reference_contigs(tabix_bin, ld_file_path)
        cache[key] = {normalise_chromosome(contig): contig for contig in contigs}
    available = cache[key]
    wanted = normalise_chromosome(chrom)
    if wanted in available:
        return available[wanted]
    raise PipelineStageError(
        "04 LD reference validation",
        "resolve_reference_contig",
        "The LD reference index contains no contig matching chromosome "
        f"{wanted} | reference_contigs={sorted(available.values())[:10]} "
        "| cause=the reference was built with a different contig naming "
        "convention or a different chromosome",
        chromosome=chrom,
        ld_file=key,
    )


def _batch_tabix_query(
    tabix_bin,
    data_path,
    chrom,
    positions,
    *,
    contig_cache,
    output_path=None,
):
    """Fetch distinct exact positions into a file or, for compatibility, text."""
    positions = sorted({int(position) for position in positions})
    if not positions:
        if output_path is not None:
            output_path = Path(output_path)
            output_path.write_text("", encoding="utf-8")
            return output_path
        return ""
    contig = resolve_reference_contig(
        chrom,
        data_path,
        tabix_bin=tabix_bin,
        contig_cache=contig_cache,
    )
    try:
        with TemporaryDirectory(prefix="postgwas_ld_regions_") as directory:
            regions_path = Path(directory) / "positions.tsv"
            with regions_path.open("w", encoding="utf-8") as handle:
                for position in positions:
                    handle.write(f"{contig}\t{position}\t{position}\n")
            command = [tabix_bin, "--regions", str(regions_path), str(data_path)]
            if output_path is None:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                )
            else:
                output_path = Path(output_path)
                with output_path.open("w", encoding="utf-8") as output_handle:
                    result = subprocess.run(
                        command,
                        stdout=output_handle,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True,
                    )
    except FileNotFoundError as error:
        raise PipelineStageError(
            "04 LD retrieval",
            "_batch_tabix_query",
            "tabix executable was not found",
            chromosome=chrom,
            reference_file=data_path,
        ) from error
    except subprocess.CalledProcessError as error:
        stderr = (error.stderr or "").strip()
        raise PipelineStageError(
            "04 LD retrieval",
            "_batch_tabix_query",
            "tabix returned exit code %s%s"
            % (
                error.returncode,
                ": %s" % stderr if stderr else "",
            ),
            chromosome=chrom,
            reference_file=data_path,
        ) from error
    return result.stdout if output_path is None else output_path


def _batch_tabix_table(
    tabix_bin,
    data_path,
    chrom,
    positions,
    *,
    contig_cache,
    parser,
):
    """Parse one batched tabix query without retaining its stdout in memory."""
    with TemporaryDirectory(prefix="postgwas_ld_payload_") as directory:
        source = _batch_tabix_query(
            tabix_bin,
            data_path,
            chrom,
            positions,
            contig_cache=contig_cache,
            output_path=Path(directory) / "payload.tsv",
        )
        return parser(source, reference_file=data_path, chromosome=chrom)


def _csv_payload_source(payload):
    """Return ``(is_empty, source)`` for a text payload or temporary file."""
    if isinstance(payload, Path):
        try:
            return payload.stat().st_size == 0, payload
        except OSError as error:
            raise RuntimeError(
                "Cannot inspect tabix payload %s: %s" % (payload, error)
            )
    return not payload.strip(), StringIO(payload)


def _parse_inventory_payload(payload, *, reference_file, chromosome):
    empty, source = _csv_payload_source(payload)
    if empty:
        return pl.DataFrame(
            schema={
                "chromosome": pl.String,
                "position": pl.Int64,
                "reference_id": pl.String,
                "canonical_id": pl.String,
                "allele_1": pl.String,
                "allele_2": pl.String,
                "minor_allele_frequency": pl.Float64,
            }
        )
    try:
        frame = pl.read_csv(
            source,
            separator="\t",
            has_header=False,
            new_columns=LD_VARIANT_INVENTORY_COLUMNS,
            schema_overrides={
                "chromosome": pl.String,
                "position": pl.Int64,
                "reference_id": pl.String,
                "canonical_id": pl.String,
                "allele_1": pl.String,
                "allele_2": pl.String,
                "minor_allele_frequency": pl.Float64,
            },
        )
    except Exception as error:
        raise function_error(
            "04 LD reference validation",
            "_parse_inventory_payload",
            error,
            chromosome=chromosome,
            reference_file=reference_file,
        ) from error
    return frame


def _parse_ld_payload(payload, *, reference_file, chromosome):
    empty, source = _csv_payload_source(payload)
    if empty:
        return pl.DataFrame(
            schema={
                "chr_a": pl.String,
                "pos_a": pl.Int64,
                "snp_a": pl.String,
                "chr_b": pl.String,
                "pos_b": pl.Int64,
                "snp_b": pl.String,
                "r2": pl.Float64,
            }
        )
    try:
        frame = pl.read_csv(
            source,
            separator="\t",
            has_header=False,
            new_columns=LD_REFERENCE_COLUMNS,
            schema_overrides={
                "chr_a": pl.String,
                "pos_a": pl.Int64,
                "snp_a": pl.String,
                "chr_b": pl.String,
                "pos_b": pl.Int64,
                "snp_b": pl.String,
                "r2": pl.Float64,
            },
        )
    except Exception as error:
        raise function_error(
            "04 LD retrieval",
            "_parse_ld_payload",
            error,
            chromosome=chromosome,
            reference_file=reference_file,
        ) from error
    invalid = frame.filter(
        pl.any_horizontal(
            [pl.col(column).is_null() for column in LD_REFERENCE_COLUMNS]
        )
        | (pl.col("snp_a") == "")
        | (pl.col("snp_b") == "")
        | (chromosome_expression("chr_a") != normalise_chromosome(chromosome))
        | (chromosome_expression("chr_b") != normalise_chromosome(chromosome))
        | (pl.col("pos_a") <= 0)
        | (pl.col("pos_b") <= 0)
        | (~pl.col("r2").is_finite())
        | (~pl.col("r2").is_between(0, 1, closed="both"))
    )
    if not invalid.is_empty():
        raise PipelineStageError(
            "04 LD reference validation",
            "_parse_ld_payload",
            "LD rows violate the chromosome, coordinate, or r2 contract",
            chromosome=chromosome,
            reference_file=reference_file,
            invalid_rows=invalid.height,
        )
    return frame


def _validate_inventory(inventory, *, chromosome, reference_file):
    expected = canonical_variant_expression(
        "chromosome", "position", "allele_1", "allele_2"
    )
    invalid = inventory.filter(
        pl.col("reference_id").is_null()
        | (pl.col("reference_id") == "")
        | pl.col("canonical_id").is_null()
        | (pl.col("canonical_id") == "")
        | (
            chromosome_expression("chromosome")
            != normalise_chromosome(chromosome)
        )
        | pl.col("position").is_null()
        | (pl.col("position") <= 0)
        | pl.col("minor_allele_frequency").is_null()
        | (~pl.col("minor_allele_frequency").is_finite())
        | (~pl.col("minor_allele_frequency").is_between(0, 0.5, closed="both"))
        | expected.is_null()
        | (pl.col("canonical_id") != expected)
    )
    if not invalid.is_empty():
        raise PipelineStageError(
            "04 LD reference validation",
            "_validate_inventory",
            "Variant inventory rows violate the chromosome, coordinate, MAF, "
            "or canonical allele-ID contract",
            chromosome=chromosome,
            reference_file=reference_file,
            examples=invalid["reference_id"].head(5).to_list(),
        )
    duplicate_ids = (
        inventory.group_by("reference_id")
        .agg(
            pl.col("canonical_id").n_unique().alias("canonical_count"),
            pl.col("minor_allele_frequency").n_unique().alias("maf_count"),
        )
        .filter(
            (pl.col("canonical_count") != 1)
            | (pl.col("maf_count") != 1)
        )
    )
    if not duplicate_ids.is_empty():
        raise PipelineStageError(
            "04 LD reference validation",
            "_validate_inventory",
            "A reference ID maps to multiple canonical variants or MAF values",
            chromosome=chromosome,
            reference_file=reference_file,
            examples=duplicate_ids["reference_id"].head(5).to_list(),
        )
    duplicate_variants = (
        inventory.group_by("canonical_id")
        .agg(pl.col("reference_id").n_unique().alias("reference_id_count"))
        .filter(pl.col("reference_id_count") != 1)
    )
    if not duplicate_variants.is_empty():
        raise PipelineStageError(
            "04 LD reference validation",
            "_validate_inventory",
            "A canonical allele-aware variant maps to multiple reference IDs",
            chromosome=chromosome,
            reference_file=reference_file,
            examples=duplicate_variants["canonical_id"].head(5).to_list(),
        )


def prepare_chromosome_reference(
    significant,
    chrom,
    ld_path,
    inventory_path,
    *,
    reverse_ld_path,
    orientation,
    tabix_bin,
    minimum_reference_maf,
    missing_index_action,
    log=None,
):
    """Validate all index variants and batch-load their LD adjacency."""
    emit = log or print
    contig_cache = {}
    positions = significant["poscol"].to_list()
    index_inventory = _batch_tabix_table(
        tabix_bin,
        inventory_path,
        chrom,
        positions,
        contig_cache=contig_cache,
        parser=_parse_inventory_payload,
    )
    _validate_inventory(
        index_inventory,
        chromosome=chrom,
        reference_file=inventory_path,
    )
    by_position = {
        int(position): frame
        for (position,), frame in index_inventory.partition_by(
            ["position"], as_dict=True
        ).items()
    }
    verified_rows = []
    exclusions = []
    for row in significant.iter_rows(named=True):
        position = int(row["poscol"])
        canonical_id = row["uniq_id"]
        observed = by_position.get(position)
        reason = None
        observed_ids = []
        if observed is None or observed.is_empty():
            reason = "missing_from_reference"
        else:
            observed_ids = observed["canonical_id"].unique().to_list()
            exact = observed.filter(pl.col("canonical_id") == canonical_id)
            if exact.is_empty():
                reason = "allele_mismatch"
            elif exact["minor_allele_frequency"].max() < minimum_reference_maf:
                reason = "below_reference_maf"
            else:
                verified_rows.append(exact)
        if reason is None:
            continue
        record = {
            "chromosome": normalise_chromosome(chrom),
            "position": position,
            "canonical_id": canonical_id,
            "input_variant_id": row.get("rsIDcol"),
            "reason": reason,
            "observed_reference_ids": ";".join(map(str, observed_ids)),
            "action": "error" if missing_index_action == "error" else "warning_skip",
        }
        exclusions.append(record)
        severity = "ERROR" if missing_index_action == "error" else "WARNING"
        emit(
            f"[{severity}] [STAGE: 04 LD reference validation] "
            "[FUNCTION: prepare_chromosome_reference] Index variant excluded "
            "| chromosome=%s | position=%s | expected_id=%s | reason=%s "
            "| observed_reference_ids=%s | action=%s"
            % (
                record["chromosome"],
                position,
                canonical_id,
                reason,
                record["observed_reference_ids"] or "none",
                record["action"],
            )
        )
    if exclusions and missing_index_action == "error":
        raise PipelineStageError(
            "04 LD reference validation",
            "prepare_chromosome_reference",
            "Significant index variants are absent from, allele-incompatible "
            "with, or below the MAF threshold of the LD reference",
            chromosome=chrom,
            excluded_indexes=len(exclusions),
            examples=[record["canonical_id"] for record in exclusions[:5]],
            action="error",
        )
    if not verified_rows:
        return ChromosomeReferenceData(
            partners={},
            verified_indexes=frozenset(),
            exclusions=tuple(exclusions),
            ld_rows=0,
        )
    verified_inventory = pl.concat(verified_rows).unique(
        ["reference_id", "canonical_id"]
    )
    verified_indexes = frozenset(verified_inventory["canonical_id"].to_list())
    verified_positions = verified_inventory["position"].to_list()
    ld_rows = _batch_tabix_table(
        tabix_bin,
        ld_path,
        chrom,
        verified_positions,
        contig_cache=contig_cache,
        parser=_parse_ld_payload,
    )
    if orientation == "upper_triangle_dual_index":
        reverse_rows = _batch_tabix_table(
            tabix_bin,
            reverse_ld_path,
            chrom,
            verified_positions,
            contig_cache=contig_cache,
            parser=_parse_ld_payload,
        )
        ld_rows = pl.concat([ld_rows, reverse_rows], how="vertical_relaxed")
    if ld_rows.is_empty():
        emit(
            "[INFO] [STAGE: 04 LD retrieval] "
            "[FUNCTION: prepare_chromosome_reference] All verified indexes have "
            "no stored LD pairs; each will be retained as self-only "
            "| chromosome=%s | verified_indexes=%s" % (chrom, len(verified_indexes))
        )
        return ChromosomeReferenceData(
            partners={},
            verified_indexes=verified_indexes,
            exclusions=tuple(exclusions),
            ld_rows=0,
        )

    ld_endpoints = pl.concat(
        [
            ld_rows.select(
                pl.col("pos_a").alias("position"),
                pl.col("snp_a").alias("reference_id"),
            ),
            ld_rows.select(
                pl.col("pos_b").alias("position"),
                pl.col("snp_b").alias("reference_id"),
            ),
        ]
    ).unique()
    endpoint_positions = ld_endpoints["position"].unique().to_list()
    endpoint_inventory = _batch_tabix_table(
        tabix_bin,
        inventory_path,
        chrom,
        endpoint_positions,
        contig_cache=contig_cache,
        parser=_parse_inventory_payload,
    )
    _validate_inventory(
        endpoint_inventory,
        chromosome=chrom,
        reference_file=inventory_path,
    )
    missing_endpoint_ids = (
        ld_endpoints.select("reference_id")
        .unique()
        .join(
            endpoint_inventory.select("reference_id").unique(),
            on="reference_id",
            how="anti",
        )
        .sort("reference_id")
    )
    if not missing_endpoint_ids.is_empty():
        raise PipelineStageError(
            "04 LD reference validation",
            "prepare_chromosome_reference",
            "LD pair endpoints are absent from the reference variant inventory",
            chromosome=chrom,
            missing_endpoint_ids=missing_endpoint_ids.height,
            examples=missing_endpoint_ids["reference_id"].head(5).to_list(),
            inventory_file=inventory_path,
        )
    index_mapping = verified_inventory.select(
        pl.col("reference_id").alias("snp_a"),
        pl.col("canonical_id").alias("_index_id"),
    ).unique("snp_a")
    endpoint_mapping = endpoint_inventory.select(
        pl.col("reference_id").alias("snp_b"),
        pl.col("canonical_id").alias("uniq_id"),
        pl.col("minor_allele_frequency").alias("_partner_maf"),
    ).unique("snp_b")
    resolved_before_maf = (
        ld_rows.join(index_mapping, on="snp_a", how="left", coalesce=True)
        .join(endpoint_mapping, on="snp_b", how="left", coalesce=True)
        .filter(pl.col("_index_id").is_not_null())
    )
    low_maf_rows = resolved_before_maf.filter(
        pl.col("_partner_maf") < minimum_reference_maf
    )
    low_maf_partner_rows_excluded = low_maf_rows.height
    low_maf_partner_variants_excluded = low_maf_rows["uniq_id"].n_unique()
    resolved_rows = (
        resolved_before_maf.filter(
            pl.col("_partner_maf") >= minimum_reference_maf
        )
        .select(
            "_index_id",
            "uniq_id",
            pl.col("snp_b").alias("ld_id"),
            pl.col("chr_b").alias("ref_chr"),
            pl.col("pos_b").alias("ref_pos"),
            "r2",
        )
    )
    partners = {
        index_id: frame.drop("_index_id")
        for (index_id,), frame in resolved_rows.partition_by(
            ["_index_id"], as_dict=True
        ).items()
    }
    emit(
        "[INFO] [STAGE: 04 LD retrieval] "
        "[FUNCTION: prepare_chromosome_reference] Batched reference loaded "
        "| chromosome=%s | orientation=%s | verified_indexes=%s "
        "| excluded_indexes=%s | ld_rows_before_maf=%s | ld_rows=%s "
        "| minimum_reference_maf=%s | low_maf_partner_rows_excluded=%s "
        "| low_maf_partner_variants_excluded=%s | tabix_data_queries=%s"
        % (
            chrom,
            orientation,
            len(verified_indexes),
            len(exclusions),
            resolved_before_maf.height,
            resolved_rows.height,
            minimum_reference_maf,
            low_maf_partner_rows_excluded,
            low_maf_partner_variants_excluded,
            4 if orientation == "upper_triangle_dual_index" else 3,
        )
    )
    return ChromosomeReferenceData(
        partners=partners,
        verified_indexes=verified_indexes,
        exclusions=tuple(exclusions),
        ld_rows=resolved_rows.height,
        ld_rows_before_maf=resolved_before_maf.height,
        low_maf_partner_rows_excluded=low_maf_partner_rows_excluded,
        low_maf_partner_variants_excluded=low_maf_partner_variants_excluded,
    )


def get_ld_partners(
    chrom,
    pos,
    canonical_id,
    ld_file_path,
    r2_threshold,
    log=None,
    *,
    window_kb=None,
    reference_data=None,
):
    """Return partners for one inventory-verified index variant."""
    resolved = _resolve_standard_options(
        window_kb=window_kb,
    )
    window_kb = resolved["window_kb"]
    emit = log or print
    if reference_data is None:
        raise PipelineStageError(
            "04 LD reference validation",
            "get_ld_partners",
            "A validated chromosome inventory and batched LD lookup are required",
            chromosome=chrom,
            variant=canonical_id,
            ld_file=ld_file_path,
        )
    if canonical_id not in reference_data.verified_indexes:
        raise PipelineStageError(
            "04 LD reference validation",
            "get_ld_partners",
            "An excluded index variant reached LD clumping",
            chromosome=chrom,
            variant=canonical_id,
            ld_file=ld_file_path,
        )
    self_snp = pl.DataFrame(
        {
            "uniq_id": [canonical_id],
            "ld_id": [canonical_id],
            "ref_chr": [str(chrom)],
            "ref_pos": [int(pos)],
            "r2": [1.0],
        }
    )
    available = reference_data.partners.get(canonical_id, _empty_ld_partners())
    window_bp = int(window_kb) * 1000
    partners = available.filter(
        (pl.col("r2") >= r2_threshold)
        & ((pl.col("ref_pos") - int(pos)).abs() <= window_bp)
    )
    if partners.is_empty():
        emit(
            "[INFO] [STAGE: 04 LD retrieval] "
            "[FUNCTION: get_ld_partners] The index is present in the variant "
            "inventory, but no partners passed the LD threshold; the index is "
            "retained as self-only | chromosome=%s | position=%s "
            "| expected_id=%s | r2_threshold=%s | action=self-only"
            % (chrom, pos, canonical_id, r2_threshold)
        )
    return (
        pl.concat([self_snp, partners], how="vertical_relaxed")
        .sort("r2", descending=True)
        .unique("uniq_id", keep="first", maintain_order=True)
    )


def _window_reachable_gwas_annotations(gwas, verified_significant, window_kb):
    """Keep GWAS annotations reachable by a verified FUMA index window.

    The later membership join has LD-reference partners on its left and uses
    this frame only to determine whether a partner is GWAS-tagged and, if so,
    to attach its association fields. Reference-only partners therefore remain
    available. Keeping every GWAS row in the union of verified index windows
    also ensures that GWAS-tagged rows above ``candidate_pvalue`` are still
    recognised and excluded after the join rather than misclassified as
    reference-only variants.
    """
    index_positions = (
        verified_significant.select(
            pl.col("poscol").alias("_nearest_index_position")
        )
        .unique()
        .sort("_nearest_index_position")
    )
    window_bp = int(window_kb) * 1000
    return (
        gwas.sort("poscol")
        .join_asof(
            index_positions,
            left_on="poscol",
            right_on="_nearest_index_position",
            strategy="nearest",
        )
        .filter(
            pl.col("_nearest_index_position").is_not_null()
            & (
                (pl.col("poscol") - pl.col("_nearest_index_position")).abs()
                <= window_bp
            )
        )
        .drop("_nearest_index_position")
    )


def find_ind_sig_snps(
    chr_df,
    ld_path,
    lead_p_threshold,
    r2_clump_threshold,
    threads=None,
    log=None,
    *,
    candidate_p_threshold=None,
    window_kb=None,
    reference_data=None,
):
    del threads
    resolved = _resolve_standard_options(
        candidate_p_threshold=candidate_p_threshold,
        window_kb=window_kb,
    )
    candidate_p_threshold = resolved["candidate_p_threshold"]
    window_kb = resolved["window_kb"]
    if reference_data is None:
        raise PipelineStageError(
            "04 LD reference validation",
            "find_ind_sig_snps",
            "Validated chromosome reference data were not provided",
        )
    emit = log or print
    # Rank on LP where available: P underflows to 0.0 below ~1e-308, which is
    # precisely the range lead SNPs occupy, and would order them arbitrarily.
    order = (
        [pl.col("lpcol").is_null(), pl.col("lpcol")]
        if "lpcol" in chr_df.columns
        else [pl.lit(False), -pl.col("pcol")]
    )
    candidates = (
        chr_df.filter(
            (pl.col("pcol") <= lead_p_threshold)
            & pl.col("uniq_id").is_in(
                list(reference_data.verified_indexes)
                if reference_data is not None
                else []
            )
        )
        .sort(order, descending=[False, True])
    )
    ind_sig_clumps = []
    chrom = chr_df.select(pl.col("chrcol").first()).item()
    emit(
        "[INFO] [STAGE: 05 independent SNP clumping] "
        "[FUNCTION: find_ind_sig_snps] Started "
        f"| chromosome={chrom} | significance_threshold={lead_p_threshold} "
        f"| significant_variants={candidates.height}"
    )
    # Walk the ranked candidates once and track depleted IDs in a set. Filtering
    # the frame per index SNP rescans every remaining row and is quadratic.
    ranked = candidates.iter_rows(named=True)
    candidate_ids = set(candidates["uniq_id"].to_list())
    depleted = set()
    remaining_count = candidates.height
    for top_snp in ranked:
        sid = top_snp["uniq_id"]
        if sid in depleted:
            continue
        variants_before = remaining_count
        partners = get_ld_partners(
            top_snp["chrcol"],
            top_snp["poscol"],
            sid,
            ld_path,
            r2_clump_threshold,
            log=emit,
            window_kb=window_kb,
            reference_data=reference_data,
        )
        members = (
            partners.join(chr_df, on="uniq_id", how="left", coalesce=True)
            .with_columns(
                pl.col("pcol").is_not_null().alias("is_gwas_tagged"),
                pl.coalesce("chrcol", "ref_chr").alias("chrcol"),
                pl.coalesce("poscol", "ref_pos").alias("poscol"),
            )
            .filter(
                (~pl.col("is_gwas_tagged"))
                | (pl.col("pcol") <= candidate_p_threshold)
                | (pl.col("uniq_id") == sid)
            )
            .with_columns(
                pl.lit(sid).alias("ind_sig_SNP_id"),
                pl.col("r2").alias("r2_with_IndSig"),
            )
            .drop("ref_chr", "ref_pos", "r2")
        )
        ind_sig_clumps.append(members)

        # Only deplete the significant-SNP pool. FUMA candidates may be linked
        # to more than one independent significant SNP, and partners include
        # reference-only variants that were never in the pool.
        newly_depleted = (
            set(partners["uniq_id"].to_list())
            .intersection(candidate_ids)
            .difference(depleted)
        )
        depleted.update(newly_depleted)
        remaining_count -= len(newly_depleted)
        variants_removed = variants_before - remaining_count
        emit(
            "[INFO] [STAGE: 05 independent SNP clumping] "
            "[FUNCTION: find_ind_sig_snps] Index SNP processed "
            f"| chromosome={chrom} | selected_variant={sid} "
            f"| selected_p={top_snp['pcol']} | selected_lp={top_snp.get('lpcol')} "
            f"| significant_variants_removed={variants_removed} "
            f"| significant_variants_remaining={remaining_count}"
        )
    emit(
        "[INFO] [STAGE: 05 independent SNP clumping] "
        "[FUNCTION: find_ind_sig_snps] Completed "
        f"| chromosome={chrom} "
        f"| independent_significant_snps={len(ind_sig_clumps)}"
    )
    return pl.concat(ind_sig_clumps) if ind_sig_clumps else pl.DataFrame()


def find_lead_snps(
    ind_sig_df,
    ld_path,
    r2_lead_threshold,
    threads=None,
    log=None,
    *,
    window_kb=None,
    reference_data=None,
):
    del threads
    if ind_sig_df.is_empty():
        return pl.DataFrame()
    resolved = _resolve_standard_options(
        window_kb=window_kb,
    )
    window_kb = resolved["window_kb"]
    if reference_data is None:
        raise PipelineStageError(
            "04 LD reference validation",
            "find_lead_snps",
            "Validated chromosome reference data were not provided",
        )
    emit = log or print
    # Same reason as the independent pass: order by LP, which does not underflow.
    order = (
        [pl.col("lpcol").is_null(), pl.col("lpcol")]
        if "lpcol" in ind_sig_df.columns
        else [pl.lit(False), -pl.col("pcol")]
    )
    ind_sig_heads = ind_sig_df.filter(
        pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
    ).sort(order, descending=[False, True])
    chrom = ind_sig_heads.select(pl.col("chrcol").first()).item()
    lead_clusters = []
    emit(
        "[INFO] [STAGE: 06 lead SNP clumping] "
        "[FUNCTION: find_lead_snps] Started "
        f"| chromosome={chrom} | r2_threshold={r2_lead_threshold} "
        f"| independent_significant_snps={ind_sig_heads.height}"
    )
    head_rows = {
        row["uniq_id"]: row for row in ind_sig_heads.iter_rows(named=True)
    }
    head_rank = {
        variant_id: rank for rank, variant_id in enumerate(head_rows)
    }
    unassigned_ids = set(head_rows)
    output_schema = {
        **ind_sig_heads.schema,
        "r2_with_Lead": pl.Float64,
        "lead_SNP_id": pl.String,
    }
    for top_row in ind_sig_heads.iter_rows(named=True):
        lid = top_row["uniq_id"]
        if lid not in unassigned_ids:
            continue
        partners = get_ld_partners(
            top_row["chrcol"],
            top_row["poscol"],
            lid,
            ld_path,
            r2_lead_threshold,
            log=emit,
            window_kb=window_kb,
            reference_data=reference_data,
        ).select("uniq_id", pl.col("r2").alias("r2_with_Lead"))
        partner_r2 = dict(partners.iter_rows())
        member_ids = sorted(
            unassigned_ids.intersection(partner_r2),
            key=head_rank.__getitem__,
        )
        members = pl.DataFrame(
            [
                {
                    **head_rows[variant_id],
                    "r2_with_Lead": partner_r2[variant_id],
                    "lead_SNP_id": lid,
                }
                for variant_id in member_ids
            ],
            schema=output_schema,
        )
        lead_clusters.append(members)
        unassigned_ids.difference_update(member_ids)
        remaining_count = len(unassigned_ids)
        emit(
            "[INFO] [STAGE: 06 lead SNP clumping] "
            "[FUNCTION: find_lead_snps] Lead SNP processed "
            f"| chromosome={chrom} | selected_variant={lid} "
            f"| independent_snps_assigned={members.height} "
            f"| independent_snps_remaining={remaining_count}"
        )
    emit(
        "[INFO] [STAGE: 06 lead SNP clumping] "
        "[FUNCTION: find_lead_snps] Completed "
        f"| chromosome={chrom} | lead_snps={len(lead_clusters)}"
    )
    return pl.concat(lead_clusters) if lead_clusters else pl.DataFrame()


def _lp_of_index_snp(columns):
    """Return the index SNP's LP, falling back to -log10(P) when LP is absent."""
    index_row = pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
    from_pvalue = (-pl.col("pcol").log10()).filter(index_row).first()
    if "lpcol" not in columns:
        return from_pvalue
    return pl.coalesce(pl.col("lpcol").filter(index_row).first(), from_pvalue)


def _stronger(candidate, incumbent):
    """True when ``candidate`` is the more significant of two loci."""
    left, right = candidate.get("lp"), incumbent.get("lp")
    if left is not None and right is not None:
        return left > right
    return candidate["p"] < incumbent["p"]


def define_genomic_risk_loci(
    ind_sig_df,
    leads_df,
    merge_dist,
    global_locus_start=1,
    *,
    summary_pvalue_thresholds=None,
):
    if leads_df.is_empty():
        return (
            pl.DataFrame(),
            pl.DataFrame(),
            pl.DataFrame(),
            pl.DataFrame(),
            global_locus_start,
        )

    summary_pvalue_thresholds = summary_pvalue_thresholds or (
        load_configuration().modules.ld_clumping.summary_pvalue_thresholds
    )
    candidate_p = {
        row["uniq_id"]: row["pcol"]
        for row in ind_sig_df.filter(pl.col("is_gwas_tagged"))
        .select("uniq_id", "pcol")
        .unique("uniq_id")
        .iter_rows(named=True)
    }

    is_clusters = (
        ind_sig_df.group_by("ind_sig_SNP_id")
        .agg(
            pl.col("chrcol").first().alias("chr"),
            pl.col("poscol").min().alias("start"),
            pl.col("poscol").max().alias("end"),
            pl.col("poscol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("pos"),
            pl.col("uniq_id").unique().alias("candidate_list"),
            pl.col("uniq_id")
            .filter(pl.col("is_gwas_tagged"))
            .unique()
            .alias("gwas_candidate_list"),
            pl.col("pcol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("p"),
            _lp_of_index_snp(ind_sig_df.columns).alias("lp"),
            pl.col("becol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("beta"),
            pl.col("secol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("se"),
            pl.col("rsIDcol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("rsID"),
            pl.col("eacol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("ea"),
            pl.col("neacol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("nea"),
            pl.col("eafcol")
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .first()
            .alias("eaf"),
            pl.concat_str(
                [
                    pl.col("uniq_id"),
                    pl.lit("("),
                    pl.col("r2_with_IndSig").round(3).cast(pl.Utf8),
                    pl.lit(")"),
                ]
            )
            .filter(pl.col("uniq_id") != pl.col("ind_sig_SNP_id"))
            .alias("ld_list"),
        )
        .with_columns(
            pl.col("candidate_list").list.len().alias("n_refsnps"),
            pl.col("gwas_candidate_list").list.len().alias("n_members"),
            pl.concat_str(
                [
                    pl.lit("{"),
                    pl.col("ind_sig_SNP_id"),
                    pl.lit(": "),
                    pl.col("ld_list").list.join("; "),
                    pl.lit("}"),
                ]
            ).alias("IndSig_Group"),
        )
        .sort(["chr", "start"])
    )

    lead_map = leads_df.select(
        pl.col("uniq_id").alias("ind_sig_SNP_id"),
        "lead_SNP_id",
        "r2_with_Lead",
    ).unique("ind_sig_SNP_id")

    initial = is_clusters.join(
        lead_map, on="ind_sig_SNP_id", how="left", coalesce=True
    )
    initial = (
        initial.group_by("lead_SNP_id")
        .agg(
            pl.col("chr").first(),
            pl.col("start").min(),
            pl.col("end").max(),
            pl.col("p").min(),
            pl.col("lp").max().alias("lp"),
            pl.col("pos")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_pos"),
            pl.col("beta")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_beta"),
            pl.col("se")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_se"),
            pl.col("rsID")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_rsid"),
            pl.col("ea")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_ea"),
            pl.col("nea")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_nea"),
            pl.col("eaf")
            .filter(pl.col("ind_sig_SNP_id") == pl.col("lead_SNP_id"))
            .first()
            .alias("l_eaf"),
            pl.col("ind_sig_SNP_id").alias("is_list"),
            pl.col("IndSig_Group").alias("is_groups"),
            pl.col("candidate_list").list.explode().unique().alias("candidate_list"),
            pl.col("gwas_candidate_list")
            .list.explode()
            .unique()
            .alias("gwas_candidate_list"),
            pl.concat_str(
                [
                    pl.col("ind_sig_SNP_id"),
                    pl.lit("("),
                    pl.col("r2_with_Lead").round(3).cast(pl.Utf8),
                    pl.lit(")"),
                ]
            )
            .filter(pl.col("ind_sig_SNP_id") != pl.col("lead_SNP_id"))
            .alias("l_ld_list"),
        )
        .with_columns(
            pl.col("candidate_list").list.len().alias("sum_refsnps"),
            pl.col("gwas_candidate_list").list.len().alias("sum_members"),
            pl.concat_str(
                [
                    pl.lit("{"),
                    pl.col("lead_SNP_id"),
                    pl.lit(": "),
                    pl.col("l_ld_list").list.join("; "),
                    pl.lit("}"),
                ]
            ).alias("Lead_Group"),
        )
        .sort(["chr", "start"])
    )

    merged = []
    idx = global_locus_start
    for row in initial.iter_rows(named=True):
        locus = {
            **row,
            "all_is": set(row["is_list"]),
            "all_is_g": set(row["is_groups"]),
            "all_l": {row["lead_SNP_id"]},
            "all_l_g": {row["Lead_Group"]},
            "all_candidates": set(row["candidate_list"]),
            "all_gwas_candidates": set(row["gwas_candidate_list"]),
            "m_dist": False,
        }
        if (
            merged
            and locus["chr"] == merged[-1]["chr"]
            and locus["start"] <= merged[-1]["end"] + merge_dist
        ):
            curr = merged[-1]
            curr["start"] = min(curr["start"], locus["start"])
            curr["end"] = max(curr["end"], locus["end"])
            curr["all_is"].update(locus["all_is"])
            curr["all_is_g"].update(locus["all_is_g"])
            curr["all_l"].update(locus["all_l"])
            curr["all_l_g"].update(locus["all_l_g"])
            curr["all_candidates"].update(locus["all_candidates"])
            curr["all_gwas_candidates"].update(locus["all_gwas_candidates"])
            curr["m_dist"] = True
            # Compare on LP, not P: two loci both stronger than ~1e-308 both
            # have p == 0.0, so a P comparison would keep whichever came first.
            if _stronger(locus, curr):
                for column in [
                    "lead_SNP_id",
                    "p",
                    "lp",
                    "l_pos",
                    "l_beta",
                    "l_se",
                    "l_rsid",
                    "l_ea",
                    "l_nea",
                    "l_eaf",
                ]:
                    curr[column] = locus[column]
        else:
            locus["Genomic_locus"] = idx
            idx += 1
            merged.append(locus)

    summary_rows = []
    l_map = []
    is_map = []
    for m in merged:
        p_values = [
            candidate_p[snp]
            for snp in m["all_gwas_candidates"]
            if snp in candidate_p and candidate_p[snp] is not None
        ]
        summary_rows.append(
            {
                "Genomic_locus": m["Genomic_locus"],
                "Locus_Length_kb": round((m["end"] - m["start"]) / 1000, 2),
                "Merged_by_Distance": int(m["m_dist"]),
                "Lead_uniqID": m["lead_SNP_id"],
                "Lead_rsID": m["l_rsid"],
                "CHR": m["chr"],
                "POS": m["l_pos"],
                "START": m["start"],
                "END": m["end"],
                "ea": m["l_ea"],
                "nea": m["l_nea"],
                "eaf": m["l_eaf"],
                "LP": m["lp"],
                "P_value": m["p"],
                "Beta": m["l_beta"],
                "SE": m["l_se"],
                "n_refsnps": len(m["all_candidates"]),
                "n_members": len(m["all_gwas_candidates"]),
                **{
                    name: sum(p <= threshold for p in p_values)
                    for name, threshold in summary_pvalue_thresholds.items()
                },
                "nIndSigSNPs": len(m["all_is"]),
                "nLeadSNPs": len(m["all_l"]),
                "IndSig_LD_Groups": " | ".join(sorted(m["all_is_g"])),
                "Lead_LD_Groups": " | ".join(sorted(m["all_l_g"])),
            }
        )
        for lid in m["all_l"]:
            l_map.append({"lead_SNP_id": lid, "Genomic_locus": m["Genomic_locus"]})
        for ind_sig_id in m["all_is"]:
            is_map.append(
                {
                    "ind_sig_SNP_id": ind_sig_id,
                    "Genomic_locus": m["Genomic_locus"],
                }
            )

    summary_df = pl.DataFrame(summary_rows) if summary_rows else pl.DataFrame()
    lead_locus_map = pl.DataFrame(l_map)
    ind_locus_map = pl.DataFrame(is_map)

    hierarchy_base = ind_sig_df.join(
        lead_map, on="ind_sig_SNP_id", how="left", coalesce=True
    ).join(ind_locus_map, on="ind_sig_SNP_id", how="left", coalesce=True)
    linked_ind_sigs = hierarchy_base.group_by(["Genomic_locus", "uniq_id"]).agg(
        pl.col("ind_sig_SNP_id").unique().alias("Linked_IndSigSNPs")
    )
    hierarchy = (
        hierarchy_base.sort("r2_with_IndSig", descending=True)
        .unique(["Genomic_locus", "uniq_id"], keep="first", maintain_order=True)
        .join(
            linked_ind_sigs,
            on=["Genomic_locus", "uniq_id"],
            how="left",
            coalesce=True,
        )
        .sort(["Genomic_locus", "poscol"])
    )

    is_clusters = (
        is_clusters.join(
            ind_locus_map, on="ind_sig_SNP_id", how="left", coalesce=True
        )
        .join(lead_map, on="ind_sig_SNP_id", how="left", coalesce=True)
        .sort(["Genomic_locus", "start"])
    )
    initial = initial.join(
        lead_locus_map, on="lead_SNP_id", how="left", coalesce=True
    ).sort(["Genomic_locus", "start"])
    return summary_df, hierarchy, is_clusters, initial, idx


def process_chromosome(
    chrom,
    gwas,
    ld_folder,
    pop,
    lead_p,
    r2_clump,
    r2_lead,
    merge_dist,
    threads=None,
    *,
    candidate_p=None,
    tabix_bin=None,
    window_kb=None,
    missing_index_action=None,
    minimum_reference_maf=None,
    reference_file_pattern=None,
    reference_reverse_file_pattern=None,
    reference_inventory_pattern=None,
    reference_orientation=None,
    reference_index_suffix=None,
    summary_pvalue_thresholds=None,
):
    """Clump one chromosome and return its loci plus ordered diagnostics.

    ``gwas`` may be the whole study: the chromosome subset is taken here rather
    than by the caller, so only the workers actually running hold a subset.
    """
    resolved = _resolve_standard_options(
        candidate_p=candidate_p,
        tabix_bin=tabix_bin,
        window_kb=window_kb,
        missing_index_action=missing_index_action,
        minimum_reference_maf=minimum_reference_maf,
        reference_file_pattern=reference_file_pattern,
        reference_reverse_file_pattern=reference_reverse_file_pattern,
        reference_inventory_pattern=reference_inventory_pattern,
        reference_orientation=reference_orientation,
        reference_index_suffix=reference_index_suffix,
        summary_pvalue_thresholds=summary_pvalue_thresholds,
    )
    candidate_p = resolved["candidate_p"]
    tabix_bin = resolved["tabix_bin"]
    window_kb = resolved["window_kb"]
    missing_index_action = resolved["missing_index_action"]
    minimum_reference_maf = resolved["minimum_reference_maf"]
    reference_file_pattern = resolved["reference_file_pattern"]
    reference_reverse_file_pattern = resolved["reference_reverse_file_pattern"]
    reference_inventory_pattern = resolved["reference_inventory_pattern"]
    reference_orientation = resolved["reference_orientation"]
    reference_index_suffix = resolved["reference_index_suffix"]
    summary_pvalue_thresholds = resolved["summary_pvalue_thresholds"]
    messages = []
    emit = messages.append
    # add_canonical_ids normalises the entire chromosome column once before
    # workers start. Repeating the string/regex expression over the full study
    # in every worker is redundant; only the scalar worker label needs normalising.
    canonical_chromosome = normalise_chromosome(chrom)
    gwas_subset = gwas.filter(pl.col("chrcol") == canonical_chromosome)
    significant = gwas_subset.filter(pl.col("pcol") <= lead_p)
    significant_variants = significant.height
    if significant_variants == 0:
        emit(
            "[INFO] [STAGE: 04 LD reference validation] "
            "[FUNCTION: process_chromosome] Chromosome has no variants at the "
            "configured lead threshold; no LD resource was required "
            f"| chromosome={chrom} | significance_threshold={lead_p}"
        )
        return {
            "chrom": chrom,
            "logs": messages,
            "progress": {
                "status": "no_significant_variants",
                "significant": 0,
                "input_significant": 0,
                "independent": 0,
                "lead": 0,
                "loci": 0,
            },
        }
    ld_path = Path(ld_folder) / reference_file_pattern.format(
        population=pop,
        chromosome=chrom,
    )
    if not ld_path.is_file() or ld_path.stat().st_size <= 0:
        error = PipelineStageError(
            "04 LD reference validation",
            "process_chromosome",
            "LD reference file was not found",
            chromosome=chrom,
            ld_file=str(ld_path),
        )
        error.chromosome_logs = messages
        raise error
    index_path = Path(str(ld_path) + reference_index_suffix)
    if not index_path.is_file() or index_path.stat().st_size <= 0:
        error = PipelineStageError(
            "04 LD reference validation",
            "process_chromosome",
            "LD reference tabix index was not found or is empty",
            chromosome=chrom,
            ld_file=str(ld_path),
            index_file=str(index_path),
        )
        error.chromosome_logs = messages
        raise error
    inventory_path = Path(ld_folder) / reference_inventory_pattern.format(
        population=pop,
        chromosome=chrom,
    )
    reverse_ld_path = (
        Path(ld_folder) / reference_reverse_file_pattern.format(
            population=pop,
            chromosome=chrom,
        )
        if reference_orientation == "upper_triangle_dual_index"
        else None
    )
    try:
        reference_data = prepare_chromosome_reference(
            significant,
            chrom,
            ld_path,
            inventory_path,
            reverse_ld_path=reverse_ld_path,
            orientation=reference_orientation,
            tabix_bin=tabix_bin,
            minimum_reference_maf=minimum_reference_maf,
            missing_index_action=missing_index_action,
            log=emit,
        )
    except PipelineStageError as error:
        error.chromosome_logs = messages
        raise
    reference_metrics = {
        "reference_ld_rows_before_maf": reference_data.ld_rows_before_maf,
        "reference_ld_rows_retained": reference_data.ld_rows,
        "low_maf_partner_rows_excluded": (
            reference_data.low_maf_partner_rows_excluded
        ),
        "low_maf_partner_variants_excluded": (
            reference_data.low_maf_partner_variants_excluded
        ),
    }
    analysed_significant = len(reference_data.verified_indexes)
    if analysed_significant == 0:
        emit(
            "[WARNING] [STAGE: 04 LD reference validation] "
            "[FUNCTION: process_chromosome] Every significant index variant "
            "was excluded by the validated LD-reference contract "
            f"| chromosome={chrom} | input_significant_variants={significant_variants} "
            f"| excluded_indexes={len(reference_data.exclusions)}"
        )
        return {
            "chrom": chrom,
            "logs": messages,
            "exclusions": list(reference_data.exclusions),
            "progress": {
                "status": "all_indexes_excluded",
                "significant": 0,
                "input_significant": significant_variants,
                "excluded_indexes": len(reference_data.exclusions),
                "independent": 0,
                "lead": 0,
                "loci": 0,
                **reference_metrics,
            },
        }
    verified_significant = significant.filter(
        pl.col("uniq_id").is_in(list(reference_data.verified_indexes))
    )
    if verified_significant.height != analysed_significant:
        error = PipelineStageError(
            "04 LD reference validation",
            "process_chromosome",
            "Verified LD-reference indexes do not map one-to-one to the "
            "canonical significant GWAS variants",
            chromosome=chrom,
            verified_reference_indexes=analysed_significant,
            matched_gwas_indexes=verified_significant.height,
        )
        error.chromosome_logs = messages
        raise error
    join_annotations = _window_reachable_gwas_annotations(
        gwas_subset,
        verified_significant,
        window_kb,
    )
    emit(
        "[INFO] [STAGE: 05 independent SNP clumping] "
        "[FUNCTION: _window_reachable_gwas_annotations] Prepared the GWAS "
        "annotation subset used by FUMA-style LD membership joins "
        f"| chromosome={chrom} | window_kb={window_kb} "
        f"| verified_indexes={analysed_significant} "
        f"| chromosome_gwas_variants={gwas_subset.height} "
        f"| join_annotation_variants={join_annotations.height} "
        "| scope=internal_annotation_join_only "
        "| reference_only_partners=preserved"
    )
    # Step A: Find independent significant SNPs.
    try:
        ind_sig = find_ind_sig_snps(
            join_annotations,
            str(ld_path),
            lead_p,
            r2_clump,
            log=emit,
            candidate_p_threshold=candidate_p,
            window_kb=window_kb,
            reference_data=reference_data,
        )
    except PipelineStageError as error:
        error.chromosome_logs = messages
        raise
    except Exception as error:
        wrapped_error = function_error(
            "05 independent SNP clumping",
            "find_ind_sig_snps",
            error,
            chromosome=chrom,
            ld_file=ld_path,
        )
        wrapped_error.chromosome_logs = messages
        raise wrapped_error from error
    if ind_sig.is_empty():
        emit(
            "[INFO] [STAGE: 05 independent SNP clumping] "
            "[FUNCTION: process_chromosome] No independent significant variants "
            "were produced after reference validation "
            f"| chromosome={chrom} | verified_indexes={analysed_significant}"
        )
        return {
            "chrom": chrom,
            "logs": messages,
            "exclusions": list(reference_data.exclusions),
            "progress": {
                "status": "no_independent_variants",
                "significant": analysed_significant,
                "input_significant": significant_variants,
                "excluded_indexes": len(reference_data.exclusions),
                "independent": 0,
                "lead": 0,
                "loci": 0,
                **reference_metrics,
            },
        }
    # Step B: Find lead SNPs.
    try:
        leads = find_lead_snps(
            ind_sig,
            str(ld_path),
            r2_lead,
            threads,
            log=emit,
            window_kb=window_kb,
            reference_data=reference_data,
        )
    except PipelineStageError as error:
        error.chromosome_logs = messages
        raise
    except Exception as error:
        wrapped_error = function_error(
            "06 lead SNP clumping",
            "find_lead_snps",
            error,
            chromosome=chrom,
            ld_file=ld_path,
        )
        wrapped_error.chromosome_logs = messages
        raise wrapped_error from error
    if leads.is_empty():
        emit(
            "[INFO] [STAGE: 06 lead SNP clumping] "
            "[FUNCTION: process_chromosome] Chromosome skipped because lead-SNP "
            "clumping produced no variants "
            f"| chromosome={chrom} | independent_significant_snps={ind_sig.height}"
        )
        return {
            "chrom": chrom,
            "logs": messages,
            "exclusions": list(reference_data.exclusions),
            "progress": {
                "status": "no_lead_variants",
                "significant": analysed_significant,
                "input_significant": significant_variants,
                "excluded_indexes": len(reference_data.exclusions),
                "independent": ind_sig.filter(
                    pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
                ).height,
                "lead": 0,
                "loci": 0,
                **reference_metrics,
            },
        }
    # Step C: Define boundaries
    emit(
        "[INFO] [STAGE: 07 locus definition] "
        "[FUNCTION: define_genomic_risk_loci] Started "
        f"| chromosome={chrom} | merge_distance={merge_dist}"
    )
    try:
        summ, hier, is_c, l_un, _ = define_genomic_risk_loci(
            ind_sig,
            leads,
            merge_dist,
            1,
            summary_pvalue_thresholds=summary_pvalue_thresholds,
        )
    except PipelineStageError as error:
        error.chromosome_logs = messages
        raise
    except Exception as error:
        wrapped_error = function_error(
            "07 locus definition",
            "define_genomic_risk_loci",
            error,
            chromosome=chrom,
            ld_file=ld_path,
        )
        wrapped_error.chromosome_logs = messages
        raise wrapped_error from error
    emit(
        "[INFO] [STAGE: 07 locus definition] "
        "[FUNCTION: define_genomic_risk_loci] Completed "
        f"| chromosome={chrom} | genomic_risk_loci={summ.height}"
    )
    return {
        "summ": summ,
        "hier": hier,
        "is_c": is_c,
        "l_un": l_un,
        "chrom": chrom,
        "logs": messages,
        "exclusions": list(reference_data.exclusions),
        "progress": {
            "status": "completed",
            "significant": analysed_significant,
            "input_significant": significant_variants,
            "excluded_indexes": len(reference_data.exclusions),
            "independent": ind_sig.filter(
                pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
            ).height,
            "lead": leads["lead_SNP_id"].n_unique(),
            "loci": summ.height,
            **reference_metrics,
        },
    }


def _append_chromosome_diagnostics(log_file, chrom, messages, error=None):
    """Append ordered worker diagnostics without echoing them to the terminal."""
    with Path(log_file).open("a", encoding="utf-8") as handle:
        handle.write(f"\n--- Chromosome {chrom} clumping diagnostics ---\n")
        for message in messages:
            handle.write(f"{message}\n")
        if error is not None:
            handle.write(f"[ERROR] {error}\n")
        handle.flush()


def _missing_chromosome_references(
    significant_counts,
    ld_folder,
    population,
    configuration,
):
    """Describe significant chromosomes lacking a usable LD file or index."""
    missing = []
    for row in significant_counts.iter_rows(named=True):
        chromosome = normalise_chromosome(row["chrcol"])
        ld_path = Path(ld_folder) / configuration.reference.file_pattern.format(
            population=population,
            chromosome=chromosome,
        )
        inventory_path = Path(
            ld_folder
        ) / configuration.reference.variant_inventory_pattern.format(
            population=population,
            chromosome=chromosome,
        )
        data_paths = [ld_path, inventory_path]
        if configuration.reference.orientation == "upper_triangle_dual_index":
            data_paths.append(
                Path(ld_folder)
                / configuration.reference.reverse_file_pattern.format(
                    population=population,
                    chromosome=chromosome,
                )
            )
        required_paths = [
            path
            for data_path in data_paths
            for path in (
                data_path,
                Path(str(data_path) + configuration.reference.index_suffix),
            )
        ]
        absent = []
        for path in required_paths:
            try:
                metadata = path.stat()
            except FileNotFoundError:
                absent.append(str(path))
                continue
            except OSError as error:
                raise PipelineStageError(
                    "04 LD reference validation",
                    "_missing_chromosome_references",
                    "Cannot inspect an LD reference resource: %s" % error,
                    chromosome=chromosome,
                    resource=str(path),
                ) from error
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_size <= 0:
                absent.append(str(path))
        if absent:
            missing.append(
                {
                    "chromosome": chromosome,
                    "significant_variants": int(row["significant"]),
                    "ld_file": str(ld_path),
                    "missing_resources": absent,
                }
            )
    return sorted(
        missing,
        key=lambda record: chromosome_sort_key(record["chromosome"]),
    )


def _warn_skipped_chromosome(record, detailed_log, logger=None):
    """Publish one explicit partial-reference warning to every audit surface."""
    chromosome = record["chromosome"]
    significant = record["significant_variants"]
    resource_count = len(record["missing_resources"])
    message = (
        "[WARNING] [STAGE: 04 LD reference validation] "
        "[FUNCTION: ld_clump_standard] Chromosome skipped because a required LD "
        "table, reverse sidecar, inventory, or tabix index is missing or empty "
        f"| chromosome={chromosome} | significant_variants_skipped={significant} "
        f"| missing_resources={record['missing_resources']} "
        "| action=skip_chromosome | consequence=results are not genome-wide"
    )
    _append_chromosome_diagnostics(detailed_log, chromosome, [message])
    print(
        screen_field(
            "warning",
            f"chr{chromosome} skipped",
            f"{significant:,} index "
            f"{'variant' if significant == 1 else 'variants'} · "
            f"{resource_count:,} LD "
            f"{'resource' if resource_count == 1 else 'resources'} "
            f"missing/empty · no chr{chromosome} loci",
            indent=6,
            label_width=24,
        ),
        flush=True,
    )
    if logger is not None:
        logger.record(
            "WARNING",
            "standard_ld_reference_chromosome_skipped",
            chromosome=chromosome,
            significant_variants_skipped=significant,
            missing_resources=record["missing_resources"],
            action="skip_chromosome",
            consequence="results_not_genome_wide",
        )


def _reference_coverage_metadata(chromosome_progress, missing_references):
    """Return machine-readable counts for analysed and skipped associations."""
    analysed_significant = sum(item["significant"] for item in chromosome_progress)
    input_significant = sum(
        item.get("input_significant", item["significant"])
        for item in chromosome_progress
    ) + sum(
        record["significant_variants"] for record in missing_references
    )
    skipped_significant = input_significant - analysed_significant
    metadata = {
        "reference_coverage_status": (
            "partial" if skipped_significant else "complete"
        ),
        "input_significant_variants": input_significant,
        "significant_variants": analysed_significant,
        "skipped_significant_variants": skipped_significant,
        "skipped_chromosomes": [
            record["chromosome"] for record in missing_references
        ],
        "missing_reference_resources": [
            resource
            for record in missing_references
            for resource in record["missing_resources"]
        ],
        "warnings": (
            sum(item["warnings"] for item in chromosome_progress)
            + len(missing_references)
        ),
        "other_warnings": sum(
            item.get("other_warnings", 0) for item in chromosome_progress
        ),
        "chromosomes_analysed": len(chromosome_progress),
        "chromosomes_with_risk_loci": [
            item["chrom"] for item in chromosome_progress if item["loci"] > 0
        ],
        "chromosomes_without_significant_variants": [
            item["chrom"]
            for item in chromosome_progress
            if item["status"] == "no_significant_variants"
        ],
        "chromosomes_significant_without_risk_locus": [
            item["chrom"]
            for item in chromosome_progress
            if item["status"] != "no_significant_variants"
            and item["loci"] == 0
        ],
        "independent_significant_snps": sum(
            item["independent"] for item in chromosome_progress
        ),
        "lead_snps": sum(item["lead"] for item in chromosome_progress),
        "genomic_risk_loci": sum(
            item["loci"] for item in chromosome_progress
        ),
        "chromosome_results": [dict(item) for item in chromosome_progress],
        "missing_reference_details": [
            dict(record) for record in missing_references
        ],
    }
    for field in (
        "reference_ld_rows_before_maf",
        "reference_ld_rows_retained",
        "low_maf_partner_rows_excluded",
        "low_maf_partner_variants_excluded",
    ):
        metadata[field] = sum(item.get(field, 0) for item in chromosome_progress)
    return metadata


def _missing_chromosome_exclusions(gwas, missing_references, lead_p):
    missing = {record["chromosome"] for record in missing_references}
    if not missing:
        return []
    rows = gwas.filter(
        (pl.col("pcol") <= lead_p)
        & chromosome_expression("chrcol").is_in(list(missing))
    )
    return [
        {
            "chromosome": normalise_chromosome(row["chrcol"]),
            "position": int(row["poscol"]),
            "canonical_id": row["uniq_id"],
            "input_variant_id": row.get("rsIDcol"),
            "reason": "missing_chromosome_reference",
            "observed_reference_ids": "",
            "action": "warning_skip",
        }
        for row in rows.iter_rows(named=True)
    ]


def _write_reference_exclusions(
    exclusions,
    output_directory,
    dataset_id,
    population,
    configuration,
):
    output_path = configured_output_path(
        output_directory,
        configuration.output_layout.standard_reference_exclusions,
        error_type=RuntimeError,
        dataset_id=dataset_id,
        population=population,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if exclusions:
        frame = pl.DataFrame(
            exclusions,
            schema=LD_REFERENCE_EXCLUSION_SCHEMA,
        )
    else:
        frame = pl.DataFrame(schema=LD_REFERENCE_EXCLUSION_SCHEMA)
    frame.write_csv(output_path, separator=configuration.table.delimiter)
    return output_path


def _classify_reference_exclusions_by_locus(exclusions, locus_summary=None):
    """Classify excluded indexes by final locus-boundary containment.

    Coordinate containment is descriptive only. An excluded index was never
    admitted to the validated LD graph, so containment must not be interpreted
    as LD membership or as evidence that the index contributed to the locus.
    Final merged loci are non-overlapping within a chromosome by construction;
    that invariant is checked before the sorted interval lookup.
    """
    if not exclusions:
        return []
    base = (
        pl.DataFrame(exclusions)
        .select(LD_REFERENCE_EXCLUSION_INPUT_COLUMNS)
        .with_row_index("_exclusion_order")
    )
    if locus_summary is None or locus_summary.is_empty():
        return (
            base.with_columns(
                pl.lit(EXCLUSION_OUTSIDE_LOCUS_BOUNDARIES).alias(
                    "reported_locus_boundary_status"
                ),
                pl.lit("").alias("overlapping_genomic_loci"),
                pl.lit("").alias("overlapping_locus_chromosome"),
                pl.lit(None, dtype=pl.Int64).alias(
                    "overlapping_locus_start"
                ),
                pl.lit(None, dtype=pl.Int64).alias(
                    "overlapping_locus_end"
                ),
                pl.lit("outside all reported locus boundaries").alias(
                    "locus_membership_interpretation"
                ),
            )
            .select(LD_REFERENCE_EXCLUSION_COLUMNS)
            .to_dicts()
        )
    required = {"Genomic_locus", "CHR", "START", "END"}
    missing = sorted(required.difference(locus_summary.columns))
    if missing:
        raise PipelineStageError(
            "09 result aggregation",
            "_classify_reference_exclusions_by_locus",
            "Final locus summary lacks columns required for excluded-index "
            "boundary classification",
            missing_columns=missing,
            available_columns=locus_summary.columns,
        )
    boundaries = (
        locus_summary.select(
            chromosome_expression("CHR").alias("_locus_chromosome"),
            pl.col("START").cast(pl.Int64, strict=False).alias("_locus_start"),
            pl.col("END").cast(pl.Int64, strict=False).alias("_locus_end"),
            pl.col("Genomic_locus").cast(pl.Int64, strict=False).alias(
                "_locus_id"
            ),
        )
        .sort(["_locus_chromosome", "_locus_start"])
        .with_columns(
            pl.col("_locus_end")
            .shift(1)
            .over("_locus_chromosome")
            .alias("_previous_locus_end")
        )
    )
    invalid = boundaries.filter(
        pl.any_horizontal(
            pl.col(column).is_null()
            for column in (
                "_locus_chromosome", "_locus_start", "_locus_end", "_locus_id"
            )
        )
        | (pl.col("_locus_start") <= 0)
        | (pl.col("_locus_end") < pl.col("_locus_start"))
        | (
            pl.col("_previous_locus_end").is_not_null()
            & (pl.col("_locus_start") <= pl.col("_previous_locus_end"))
        )
    )
    if not invalid.is_empty():
        raise PipelineStageError(
            "09 result aggregation",
            "_classify_reference_exclusions_by_locus",
            "Final reported locus boundaries are invalid or overlap within a "
            "chromosome; excluded-index containment would be ambiguous",
            invalid_loci=invalid["_locus_id"].head(10).to_list(),
        )
    classified = (
        base.sort(["chromosome", "position"])
        .join_asof(
            boundaries.drop("_previous_locus_end"),
            left_on="position",
            right_on="_locus_start",
            by_left="chromosome",
            by_right="_locus_chromosome",
            strategy="backward",
        )
        .with_columns(
            (
                pl.col("_locus_id").is_not_null()
                & (pl.col("position") <= pl.col("_locus_end"))
            ).alias("_inside_locus_boundary")
        )
        .with_columns(
            pl.when(pl.col("_inside_locus_boundary"))
            .then(pl.lit(EXCLUSION_INSIDE_LOCUS_BOUNDARY))
            .otherwise(pl.lit(EXCLUSION_OUTSIDE_LOCUS_BOUNDARIES))
            .alias("reported_locus_boundary_status"),
            pl.when(pl.col("_inside_locus_boundary"))
            .then(pl.col("_locus_id").cast(pl.String))
            .otherwise(pl.lit(""))
            .alias("overlapping_genomic_loci"),
            pl.when(pl.col("_inside_locus_boundary"))
            # Polars coalesces the right-side by-key after the as-of join;
            # the canonical left chromosome is the validated equivalent.
            .then(pl.col("chromosome"))
            .otherwise(pl.lit(""))
            .alias("overlapping_locus_chromosome"),
            pl.when(pl.col("_inside_locus_boundary"))
            .then(pl.col("_locus_start"))
            .otherwise(pl.lit(None, dtype=pl.Int64))
            .alias("overlapping_locus_start"),
            pl.when(pl.col("_inside_locus_boundary"))
            .then(pl.col("_locus_end"))
            .otherwise(pl.lit(None, dtype=pl.Int64))
            .alias("overlapping_locus_end"),
            pl.when(pl.col("_inside_locus_boundary"))
            .then(
                pl.lit(
                    "coordinate overlap only; excluded index was not used as "
                    "an LD member or locus seed"
                )
            )
            .otherwise(pl.lit("outside all reported locus boundaries"))
            .alias("locus_membership_interpretation"),
        )
        .sort("_exclusion_order")
        .select(LD_REFERENCE_EXCLUSION_COLUMNS)
    )
    return classified.to_dicts()


def _reference_exclusion_locus_counts(exclusions):
    counts = {
        "reference_exclusions_within_reported_locus_boundaries": 0,
        "reference_exclusions_outside_reported_locus_boundaries": 0,
    }
    for exclusion in exclusions:
        status = exclusion.get("reported_locus_boundary_status")
        if status == EXCLUSION_INSIDE_LOCUS_BOUNDARY:
            counts["reference_exclusions_within_reported_locus_boundaries"] += 1
        elif status == EXCLUSION_OUTSIDE_LOCUS_BOUNDARIES:
            counts["reference_exclusions_outside_reported_locus_boundaries"] += 1
        else:
            raise PipelineStageError(
                "09 result aggregation",
                "_reference_exclusion_locus_counts",
                "Excluded index lacks a valid reported-locus boundary status",
                canonical_id=exclusion.get("canonical_id"),
                boundary_status=status,
            )
    if sum(counts.values()) != len(exclusions):
        raise PipelineStageError(
            "09 result aggregation",
            "_reference_exclusion_locus_counts",
            "Excluded-index locus classification counts do not reconcile",
            excluded_indexes=len(exclusions),
            classified_indexes=sum(counts.values()),
        )
    return counts


def _append_reference_exclusion_locus_audit(
    detailed_log,
    exclusions,
    counts,
):
    """Write the summary and every classified excluded variant to the log."""
    with Path(detailed_log).open("a", encoding="utf-8") as handle:
        handle.write("\n--- Excluded-index relationship to reported loci ---\n")
        handle.write(
            "[INFO] [STAGE: 09 result aggregation] "
            "[FUNCTION: _classify_reference_exclusions_by_locus] Completed "
            "coordinate-only classification"
            " | excluded_indexes=%s | within_reported_locus_boundaries=%s"
            " | outside_reported_locus_boundaries=%s"
            " | interpretation=boundary overlap is not validated LD membership\n"
            % (
                len(exclusions),
                counts[
                    "reference_exclusions_within_reported_locus_boundaries"
                ],
                counts[
                    "reference_exclusions_outside_reported_locus_boundaries"
                ],
            )
        )
        for exclusion in exclusions:
            handle.write(
                "[INFO] [STAGE: 09 result aggregation] "
                "[FUNCTION: _classify_reference_exclusions_by_locus] "
                "Excluded index classified"
                " | chromosome=%s | position=%s | canonical_id=%s"
                " | input_variant_id=%s | exclusion_reason=%s"
                " | reported_locus_boundary_status=%s"
                " | overlapping_genomic_loci=%s"
                " | overlapping_locus_chromosome=%s"
                " | overlapping_locus_start=%s"
                " | overlapping_locus_end=%s"
                " | interpretation=%s\n"
                % (
                    exclusion["chromosome"],
                    exclusion["position"],
                    exclusion["canonical_id"],
                    exclusion.get("input_variant_id") or "none",
                    exclusion["reason"],
                    exclusion["reported_locus_boundary_status"],
                    exclusion.get("overlapping_genomic_loci") or "none",
                    exclusion.get("overlapping_locus_chromosome") or "none",
                    exclusion.get("overlapping_locus_start") or "none",
                    exclusion.get("overlapping_locus_end") or "none",
                    exclusion["locus_membership_interpretation"],
                )
            )
        handle.flush()


def _finalise_reference_exclusion_audit(
    exclusions,
    locus_summary,
    output_directory,
    dataset_id,
    population,
    configuration,
    detailed_log,
    logger=None,
):
    classified = _classify_reference_exclusions_by_locus(
        exclusions,
        locus_summary,
    )
    counts = _reference_exclusion_locus_counts(classified)
    output_path = _write_reference_exclusions(
        classified,
        output_directory,
        dataset_id,
        population,
        configuration,
    )
    _append_reference_exclusion_locus_audit(
        detailed_log,
        classified,
        counts,
    )
    if logger is not None:
        logger.record(
            "OUTPUT",
            "standard_reference_exclusions",
            output_file=str(output_path),
            excluded_indexes=len(classified),
            reasons=reference_exclusion_reason_counts(classified),
            **counts,
            interpretation=(
                "coordinate_boundary_overlap_only_not_validated_ld_membership"
            ),
        )
    return classified, output_path, counts


def _compact_chromosome_labels(chromosomes):
    """Format observed chromosome labels as compact, ordered ranges."""
    labels = sorted(
        {normalise_chromosome(chrom) for chrom in chromosomes},
        key=chromosome_sort_key,
    )
    numeric = [int(label) for label in labels if label.isdigit()]
    other = [label for label in labels if not label.isdigit()]
    compact = []
    if numeric:
        start = previous = numeric[0]
        for chromosome in numeric[1:] + [None]:
            if chromosome is not None and chromosome == previous + 1:
                previous = chromosome
                continue
            compact.append(
                f"chr{start}" if start == previous else f"chr{start}–{previous}"
            )
            if chromosome is not None:
                start = previous = chromosome
    compact.extend(f"chr{label}" for label in other)
    return ", ".join(compact) if compact else "none"


def _print_chromosome_progress(
    chrom,
    progress,
    warning_count,
    reference_exclusion_reasons=None,
):
    """Print only chromosome outcomes that require live user attention."""
    reference_exclusion_reasons = reference_exclusion_reasons or {}
    reference_exclusions = sum(reference_exclusion_reasons.values())
    other_warnings = max(0, warning_count - reference_exclusions)
    status = progress["status"]
    if status == "no_significant_variants" and warning_count == 0:
        return
    retained_significant = progress["significant"]
    input_significant = progress.get(
        "input_significant",
        retained_significant + reference_exclusions,
    )
    if input_significant != retained_significant + reference_exclusions:
        raise PipelineStageError(
            "08 chromosome execution",
            "_print_chromosome_progress",
            "Chromosome GWS input, retained, and exclusion counts do not "
            "reconcile",
            chromosome=chrom,
            input_significant=input_significant,
            retained_significant=retained_significant,
            excluded_significant=reference_exclusions,
        )
    clumping_outcome = (
        f"independent {progress['independent']:,} · "
        f"lead {progress['lead']:,} · loci {progress['loci']:,}"
    )
    lines = [
        screen_field(
            "warning" if warning_count or progress["loci"] == 0 else "success",
            f"chr{normalise_chromosome(chrom)} GWS variants",
            f"{input_significant:,} input · {retained_significant:,} retained · "
            f"{reference_exclusions:,} excluded",
            indent=6,
            label_width=24,
        ),
        screen_field(
            "genetic",
            "Clumping results",
            clumping_outcome,
            indent=8,
            label_width=22,
        ),
    ]
    if reference_exclusions:
        lines.append(
            screen_field(
                "warning",
                "Excluded GWS variants",
                format_compact_reference_exclusion_counts(
                    reference_exclusion_reasons
                ),
                indent=8,
                label_width=22,
            )
        )
    if other_warnings:
        lines.append(
            screen_field(
                "warning",
                "Other warnings",
                f"{other_warnings:,} · see detailed log",
                indent=8,
                label_width=22,
            )
        )
    print("\n".join(lines), end="\n\n", flush=True)


def _print_clumping_summary(
    chromosome_progress,
    detailed_log,
    skipped_references=(),
    *,
    label_width,
):
    """Print one concise end-of-run summary while the log retains full detail."""
    metrics = _reference_coverage_metadata(
        chromosome_progress, skipped_references,
    )
    no_significant = metrics["chromosomes_without_significant_variants"]
    with_loci = metrics["chromosomes_with_risk_loci"]
    significant_without_loci = metrics[
        "chromosomes_significant_without_risk_locus"
    ]
    other_warning_chromosomes = [
        progress["chrom"]
        for progress in chromosome_progress
        if progress.get("other_warnings", 0) > 0
    ]
    skipped_chromosomes = [
        record["chromosome"] for record in skipped_references
    ]
    totals = {
        "input_significant": metrics["input_significant_variants"],
        "significant": metrics["significant_variants"],
        "independent": metrics["independent_significant_snps"],
        "lead": metrics["lead_snps"],
        "loci": metrics["genomic_risk_loci"],
    }
    exclusion_reasons = {}
    for progress in chromosome_progress:
        for reason, count in (
            progress.get("reference_exclusion_reasons") or {}
        ).items():
            exclusion_reasons[reason] = exclusion_reasons.get(reason, 0) + count
    missing_chromosome_indexes = sum(
        record["significant_variants"] for record in skipped_references
    )
    if missing_chromosome_indexes:
        exclusion_reasons["missing_chromosome_reference"] = (
            exclusion_reasons.get("missing_chromosome_reference", 0)
            + missing_chromosome_indexes
        )
    reference_exclusions = sum(exclusion_reasons.values())
    other_warnings = metrics["other_warnings"]
    print()
    heading = (
        "Standard LD clumping finished with reference exclusions"
        if reference_exclusions
        else "Standard LD clumping completed"
    )
    print(screen_line("analysis", heading, indent=2))
    print(
        screen_field(
            "count",
            "Chromosomes analysed",
            f"{len(chromosome_progress):,}",
            indent=6,
            label_width=label_width,
        )
    )
    print(
        screen_field(
            "success" if with_loci else "info",
            "Chromosomes with risk loci",
            f"{_compact_chromosome_labels(with_loci)} ({len(with_loci):,})",
            indent=6,
            label_width=label_width,
        )
    )
    print(
        screen_field(
            "info",
            "No significant variants",
            f"{_compact_chromosome_labels(no_significant)} "
            f"({len(no_significant):,})",
            indent=6,
            label_width=label_width,
        )
    )
    if significant_without_loci:
        print(
            screen_field(
                "warning",
                "Significant but no risk locus",
                f"{_compact_chromosome_labels(significant_without_loci)} "
                f"({len(significant_without_loci):,})",
                indent=6,
                label_width=label_width,
            )
        )
    for kind, label, field in (
        ("count", "GWS variants supplied", "input_significant"),
        (
            "success",
            "GWS variants passing LD-reference checks",
            "significant",
        ),
    ):
        print(
            screen_field(
                kind,
                label,
                f"{totals[field]:,}",
                indent=6,
                label_width=label_width,
            )
        )
    if reference_exclusions:
        print(
            screen_field(
                "warning",
                "GWS variants excluded",
                f"{reference_exclusions:,}",
                indent=6,
                label_width=label_width,
            )
        )
        for reason, count in reference_exclusion_reason_items(exclusion_reasons):
            value = f"{count:,}"
            if reason == "missing_chromosome_reference" and skipped_chromosomes:
                value += f" · {_compact_chromosome_labels(skipped_chromosomes)}"
            print(
                screen_field(
                    "warning",
                    reference_exclusion_reason_label(reason),
                    value,
                    indent=8,
                    label_width=label_width - 2,
                )
            )
    print()
    for label, field in (
        ("Independent significant SNPs", "independent"),
        ("Lead SNPs", "lead"),
        ("Genomic risk loci", "loci"),
    ):
        print(
            screen_field(
                "genetic",
                label,
                f"{totals[field]:,}",
                indent=6,
                label_width=label_width,
            )
        )
    if reference_exclusions:
        print()
        print(
            screen_field(
                "warning",
                "Analysis limitation",
                reference_exclusion_analysis_limitation(reference_exclusions),
                indent=6,
                label_width=label_width,
            )
        )
    if other_warnings:
        warning_summary = f"{other_warnings:,}"
        if other_warning_chromosomes:
            warning_summary += (
                " across "
                f"{_compact_chromosome_labels(other_warning_chromosomes)}"
            )
        warning_summary += "; see detailed log"
        print(
            screen_field(
                "warning",
                "Other diagnostic warnings",
                warning_summary,
                indent=6,
                label_width=label_width,
            )
        )
    print(
        screen_field(
            "info",
            "Detailed log",
            Path(detailed_log).name,
            indent=6,
            label_width=label_width,
        )
    )


def ld_clump_standard(
    vcf_path,
    ld_folder,
    dataset_id,
    output_directory,
    pop=None,
    lead_p=None,
    r2_clump=None,
    r2_lead=None,
    merge_dist=None,
    bcftools_bin=None,
    threads=None,
    *,
    tabix_bin=None,
    memory_gb=None,
    screen_label_width=None,
    emit_terminal_summary=True,
    include_report_tables=False,
    configuration=None,
    reference_manifest=None,
    logger=None,
):
    application = None
    if configuration is None:
        application = load_configuration()
        configuration = application.modules.ld_clumping
    pop = pop or configuration.population.value
    if reference_manifest is not None:
        if reference_manifest.orientation != configuration.reference.orientation:
            raise PipelineStageError(
                "04 LD reference validation",
                "ld_clump_standard",
                "The validated manifest orientation changed before execution",
                manifest_orientation=reference_manifest.orientation,
                configured_orientation=configuration.reference.orientation,
            )
    resolution_source = None if application is None else application
    lead_p = configuration.lead_pvalue if lead_p is None else lead_p
    r2_clump = configuration.clump_r2 if r2_clump is None else r2_clump
    r2_lead = configuration.lead_r2 if r2_lead is None else r2_lead
    merge_dist = (
        configuration.merge_distance_bp if merge_dist is None else merge_dist
    )
    resolved = _resolve_standard_options(
        application=resolution_source,
        configuration=configuration if application is None else None,
        bcftools_bin=bcftools_bin,
        tabix_bin=tabix_bin,
        threads=threads,
        memory_gb=memory_gb,
        screen_label_width=screen_label_width,
    )
    bcftools_bin = resolved["bcftools_bin"]
    tabix_bin = resolved["tabix_bin"]
    threads = resolved["threads"]
    memory_gb = resolved["memory_gb"]
    screen_label_width = resolved["screen_label_width"]

    # 1. Conversion. The frame is returned with the published table so a
    #    genome-wide projection is not serialised and immediately re-parsed.
    try:
        tsv_path, gwas = vcf_to_standard_ldclump(
            vcf_path,
            output_directory,
            dataset_id,
            bcftools_bin,
            threads=threads,
            configuration=configuration,
            logger=logger,
        )
        _, detailed_log = _standard_conversion_paths(
            output_directory, dataset_id, configuration,
        )
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "01 VCF conversion",
            "vcf_to_standard_ldclump",
            error,
            sample=dataset_id,
            input_vcf=vcf_path,
        ) from error
    Path(detailed_log).parent.mkdir(parents=True, exist_ok=True)
    startup_messages = []
    try:
        gwas = add_canonical_ids(gwas, log=startup_messages.append)
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "03 canonical ID preparation",
            "add_canonical_ids",
            error,
            input_file=tsv_path,
        ) from error
    gwas, mhc_statistics = exclude_mhc(
        gwas, configuration, chromosome_column="chrcol", position_column="poscol",
    )
    startup_messages.append(
        "[INFO] [STAGE: 03 canonical ID preparation] [FUNCTION: exclude_mhc] "
        "MHC exclusion | %s"
        % " | ".join("%s=%s" % item for item in mhc_statistics.items())
    )
    _append_chromosome_diagnostics(detailed_log, "startup", startup_messages)
    if logger is not None:
        logger.record("ACTION", "standard_mhc_exclusion", **mhc_statistics)
        for message in startup_messages:
            if message.startswith("[WARNING]"):
                logger.record("WARNING", "standard_clumping_input", detail=message)
    # Calculate the lead-threshold counts once. Reference validation and worker
    # scheduling consume the same chromosome-level evidence on full GWAS data.
    significant_rows = (
        gwas.filter(pl.col("pcol") <= lead_p)
        .group_by("chrcol")
        .len()
        .rename({"len": "significant"})
    )
    missing_references = _missing_chromosome_references(
        significant_rows,
        ld_folder,
        pop,
        configuration,
    )
    if (
        missing_references
        and configuration.missing_chromosome_action == "error"
    ):
        missing_resources = [
            resource
            for record in missing_references
            for resource in record["missing_resources"]
        ]
        raise PipelineStageError(
            "04 LD reference validation",
            "ld_clump_standard",
            "Missing or empty LD resources for chromosomes containing significant "
            "variants: %s" % ", ".join(missing_resources),
            sample=dataset_id,
        )
    for record in missing_references:
        _warn_skipped_chromosome(record, detailed_log, logger=logger)
    skipped_chromosomes = {
        record["chromosome"] for record in missing_references
    }
    # --- BALANCED SCHEDULING ---
    try:
        # Cost is driven by the number of variants at the lead threshold, not by
        # total rows: each significant variant costs one clumping round and, at
        # worst, one tabix query. A dense chromosome must start first.
        chromosome_rows = gwas.group_by("chrcol").len().rename({"len": "rows"})
        chrom_stats = (
            chromosome_rows.join(
                significant_rows, on="chrcol", how="left", coalesce=True
            )
            .with_columns(pl.col("significant").fill_null(0))
            .sort(["significant", "rows"], descending=True)
        )
        # ThreadPoolExecutor dynamically assigns the next queued task to each
        # free worker. Submitting every costly chromosome first minimises the
        # long tail; alternating in the smallest chromosomes delays dense work.
        scheduled_chromosomes = [
            chromosome
            for chromosome in chrom_stats["chrcol"].to_list()
            if normalise_chromosome(chromosome) not in skipped_chromosomes
        ]
        # A worker holds one chromosome, not the whole study. Charging each
        # worker for the entire frame under-counts the fitting workers by
        # roughly the chromosome count and can serialise the run outright.
        total_rows = max(1, gwas.height)
        largest_share = int(chrom_stats["rows"].max() or 0) / total_rows
        worker_frame_gb = (gwas.estimated_size() / (1024**3)) * largest_share
        memory_per_worker_gb = max(
            configuration.compute.minimum_worker_memory_gb,
            worker_frame_gb * configuration.compute.input_memory_multiplier,
        )
        memory_adjusted_workers = max(1, int(memory_gb // memory_per_worker_gb))
        n_workers = max(
            1,
            min(threads, memory_adjusted_workers),
        )
    except Exception as error:
        raise function_error(
            "08 chromosome scheduling",
            "ld_clump_standard",
            error,
            input_file=tsv_path,
        ) from error
    # ---------------------------
    res_list = []
    reference_exclusions = _missing_chromosome_exclusions(
        gwas,
        missing_references,
        lead_p,
    )
    effective_workers = max(1, min(n_workers, len(scheduled_chromosomes)))
    print()
    print(screen_line("analysis", "Standard LD clumping", indent=2))
    print(
        screen_field(
            "count",
            "Chromosomes scheduled",
            f"{len(scheduled_chromosomes):,}",
            indent=6,
            label_width=24,
        )
    )
    print(
        screen_field(
            "analysis",
            "Parallel workers",
            f"{effective_workers:,}",
            indent=6,
            label_width=24,
        )
    )

    # 3. Concurrent Processing using Threads
    chromosome_progress = []
    try:
        with ThreadPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                chrom: executor.submit(
                    process_chromosome,
                    chrom,
                    # The full frame is passed, not a precomputed subset: this
                    # comprehension runs in the submitting thread, so slicing
                    # here would materialise every chromosome before any worker
                    # starts and roughly double peak memory.
                    gwas,
                    ld_folder,
                    pop,
                    lead_p,
                    r2_clump,
                    r2_lead,
                    merge_dist,
                    candidate_p=configuration.candidate_pvalue,
                    tabix_bin=tabix_bin,
                    window_kb=configuration.window_kb,
                    missing_index_action=configuration.missing_index_action,
                    minimum_reference_maf=configuration.minimum_reference_maf,
                    reference_file_pattern=configuration.reference.file_pattern,
                    reference_reverse_file_pattern=(
                        configuration.reference.reverse_file_pattern
                    ),
                    reference_inventory_pattern=(
                        configuration.reference.variant_inventory_pattern
                    ),
                    reference_orientation=configuration.reference.orientation,
                    reference_index_suffix=configuration.reference.index_suffix,
                    summary_pvalue_thresholds=(
                        configuration.summary_pvalue_thresholds
                    ),
                )
                for chrom in scheduled_chromosomes
            }
            ordered_chromosomes = sorted(futures, key=chromosome_sort_key)
            for chrom in ordered_chromosomes:
                future = futures[chrom]
                try:
                    result = future.result()
                    if result:
                        messages = result.pop("logs", [])
                        chromosome_exclusions = result.pop("exclusions", [])
                        reference_exclusions.extend(chromosome_exclusions)
                        exclusion_reasons = reference_exclusion_reason_counts(
                            chromosome_exclusions
                        )
                        _append_chromosome_diagnostics(
                            detailed_log, chrom, messages
                        )
                        warning_count = sum(
                            message.startswith("[WARNING]") for message in messages
                        )
                        progress = result.pop("progress")
                        other_warnings = max(
                            0, warning_count - len(chromosome_exclusions)
                        )
                        chromosome_progress.append(
                            {
                                "chrom": chrom,
                                "warnings": warning_count,
                                "other_warnings": other_warnings,
                                "reference_exclusions": len(
                                    chromosome_exclusions
                                ),
                                "reference_exclusion_reasons": exclusion_reasons,
                                **progress,
                            }
                        )
                        _print_chromosome_progress(
                            chrom,
                            progress,
                            warning_count,
                            exclusion_reasons,
                        )
                        if logger is not None:
                            # The canonical log otherwise records nothing about
                            # the work itself; full detail stays in detailed_log.
                            logger.record(
                                "RESULT",
                                "chromosome_clumping",
                                chromosome=chrom,
                                warnings=warning_count,
                                other_warnings=other_warnings,
                                reference_exclusions=len(
                                    chromosome_exclusions
                                ),
                                reference_exclusion_reasons=exclusion_reasons,
                                **progress,
                            )
                        if "summ" in result:
                            res_list.append(result)
                except PipelineStageError as error:
                    _append_chromosome_diagnostics(
                        detailed_log,
                        chrom,
                        getattr(error, "chromosome_logs", []),
                        error=error,
                    )
                    print(
                        screen_field(
                            "error",
                            f"chr{normalise_chromosome(chrom)} failed",
                            error,
                            indent=6,
                            label_width=24,
                        ),
                        flush=True,
                    )
                    print(
                        screen_field(
                            "info",
                            "Detailed log",
                            Path(detailed_log).name,
                            indent=6,
                            label_width=24,
                        ),
                        flush=True,
                    )
                    raise
                except Exception as error:
                    raise function_error(
                        "08 chromosome execution",
                        "process_chromosome",
                        error,
                        chromosome=chrom,
                    ) from error
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "08 chromosome execution",
            "ld_clump_standard",
            error,
            dataset=dataset_id,
        ) from error
    coverage_metadata = _reference_coverage_metadata(
        chromosome_progress,
        missing_references,
    )
    # 4. Final Aggregation
    aggregated_results = {}
    if res_list:
        try:
            res_list.sort(key=lambda result: chromosome_sort_key(result["chrom"]))
            final_res = {"summ": [], "hier": [], "is_c": [], "l_un": []}
            locus_counter = 1
            for r in res_list:
                unique_count = r["summ"]["Genomic_locus"].unique().len()
                for key, frames in final_res.items():
                    if "Genomic_locus" in r[key].columns:
                        r[key] = r[key].with_columns(
                            pl.col("Genomic_locus") + (locus_counter - 1)
                        )
                    frames.append(r[key])
                locus_counter += unique_count
            for key, frames in final_res.items():
                df = pl.concat(frames)
                # Ensure nested values match the complete TSV representation
                # before deriving the compact in-memory HTML views.
                aggregated_results[key] = df.select(
                    pl.col(column).list.join("; ")
                    if df[column].dtype.is_nested()
                    else pl.col(column)
                    for column in df.columns
                )
        except Exception as error:
            raise function_error(
                "09 result aggregation",
                "ld_clump_standard",
                error,
                dataset=dataset_id,
            ) from error
    try:
        (
            reference_exclusions,
            exclusion_path,
            exclusion_locus_counts,
        ) = _finalise_reference_exclusion_audit(
            reference_exclusions,
            aggregated_results.get("summ"),
            output_directory,
            dataset_id,
            pop,
            configuration,
            detailed_log,
            logger=logger,
        )
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "10 output writing",
            "_finalise_reference_exclusion_audit",
            error,
            dataset=dataset_id,
        ) from error
    coverage_metadata.update(
        reference_exclusions=len(reference_exclusions),
        reference_exclusions_file=str(exclusion_path),
        reference_exclusion_reasons=reference_exclusion_reason_counts(
            reference_exclusions
        ),
        **exclusion_locus_counts,
    )
    output_files = {
        "formatted_table": str(tsv_path),
        "detailed_log": str(detailed_log),
        "reference_exclusions": str(exclusion_path),
    }
    if (
        coverage_metadata["input_significant_variants"] > 0
        and coverage_metadata["significant_variants"] == 0
    ):
        error = PipelineStageError(
            "04 LD reference validation",
            "ld_clump_standard",
            "All significant variants were excluded because chromosome "
            "resources or exact allele-aware index variants were absent from "
            "the LD reference; no scientifically valid clumping result can be "
            "reported",
            sample=dataset_id,
            excluded_significant_variants=(
                coverage_metadata["skipped_significant_variants"]
            ),
            exclusions_file=exclusion_path,
        )
        error.validated_failure_outputs = (str(exclusion_path),)
        raise error
    if res_list:
        # These names connect the stable algorithm outputs to their report
        # presentation. Displayed columns remain configuration-controlled; the
        # full hierarchy and every omitted audit column stay in the TSV files.
        result_outputs = {
            "summ": (
                configuration.output_layout.standard_summary,
                "genomic_risk_loci",
            ),
            "l_un": (
                configuration.output_layout.standard_lead_clusters,
                "lead_snps",
            ),
            "is_c": (
                configuration.output_layout.standard_independent_clusters,
                "independent_significant_snps",
            ),
            "hier": (
                configuration.output_layout.standard_hierarchy,
                None,
            ),
        }
        report_tables = {}
        for key, (pattern, report_name) in result_outputs.items():
            output_file = configured_output_path(
                output_directory,
                pattern,
                error_type=RuntimeError,
                dataset_id=dataset_id,
                population=pop,
            )
            output_file.parent.mkdir(parents=True, exist_ok=True)
            df = aggregated_results[key]
            if include_report_tables and report_name is not None:
                report_columns = (
                    configuration.reporting.result_table_columns[report_name]
                )
                missing_columns = [
                    column for column in report_columns
                    if column not in df.columns
                ]
                if missing_columns:
                    raise PipelineStageError(
                        "09 result aggregation",
                        "ld_clump_standard",
                        "Configured HTML result columns are absent from the "
                        "validated clumping table",
                        result_table=report_name,
                        missing_columns=missing_columns,
                        available_columns=df.columns,
                    )
                report_tables[report_name] = {
                    "columns": list(report_columns),
                    "rows": df.select(report_columns).rows(),
                    "total_rows": df.height,
                    "output_file": str(output_file),
                }
            try:
                df.write_csv(
                    output_file, separator=configuration.table.delimiter,
                )
            except Exception as error:
                raise function_error(
                    "10 output writing",
                    "ld_clump_standard",
                    error,
                    result=key,
                    output_file=output_file,
                ) from error
            output_files[key] = str(output_file)
        if emit_terminal_summary:
            _print_clumping_summary(
                chromosome_progress,
                detailed_log,
                skipped_references=missing_references,
                label_width=screen_label_width,
            )
        summary_path = configured_output_path(
            output_directory,
            configuration.output_layout.standard_summary,
            error_type=RuntimeError,
            dataset_id=dataset_id,
            population=pop,
        )
        result = {
            "status": (
                "partial_reference"
                if coverage_metadata["reference_coverage_status"] == "partial"
                else "completed"
            ),
            "ldpruned_sig_file": str(summary_path),
            "log_file": str(detailed_log),
            "output_files": output_files,
            **coverage_metadata,
        }
        if logger is not None:
            logger.record(
                "OUTPUT", "standard_clumping", **result,
                allele_orientation=(
                    "REF/ALT retained; canonical matching ID is allele-order-independent"
                ),
            )
        if include_report_tables:
            # Presentation rows are intentionally transient: service.py consumes
            # them while writing HTML and removes them before return/checkpoint
            # logging, preventing large result tables from entering metadata.
            result["_report_tables"] = report_tables
            result["_reference_exclusion_details"] = reference_exclusions
        return result
    if emit_terminal_summary:
        _print_clumping_summary(
            chromosome_progress,
            detailed_log,
            skipped_references=missing_references,
            label_width=screen_label_width,
        )
    result = {
        "status": (
            "partial_reference"
            if coverage_metadata["reference_coverage_status"] == "partial"
            else "no_loci"
        ),
        "ldpruned_sig_file": None,
        "log_file": str(detailed_log),
        "output_files": output_files,
        **coverage_metadata,
    }
    if logger is not None:
        logger.record("RESULT", "standard_clumping_no_loci", **result)
    if include_report_tables:
        result["_reference_exclusion_details"] = reference_exclusions
    return result
