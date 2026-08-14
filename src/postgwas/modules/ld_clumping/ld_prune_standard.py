
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from pathlib import Path

import polars as pl

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.core.ui import screen_field, screen_line
from postgwas.core.vcf import extract_vcf_table


class PipelineStageError(RuntimeError):
    """Pipeline error with a stable stage, function, and execution context."""

    def __init__(self, stage, function, message, **context):
        details = " | ".join(
            f"{key}={value}" for key, value in context.items() if value is not None
        )
        text = f"[STAGE: {stage}] [FUNCTION: {function}] {message}"
        if details:
            text = f"{text} | {details}"
        super().__init__(text)
        self.stage = stage
        self.function = function
        self.context = context


def function_error(stage, function, error, **context):
    """Create a contextual error while retaining the original exception."""
    return PipelineStageError(
        stage,
        function,
        f"{type(error).__name__}: {error}",
        **context,
    )


def _resolve_standard_options(*, application=None, **provided):
    """Resolve omitted helper options once without reloading configured calls."""
    if all(value is not None for value in provided.values()):
        return provided
    application = application or load_configuration()
    module = application.modules.ld_clumping
    defaults = {
        "candidate_p_threshold": module.candidate_pvalue,
        "candidate_p": module.candidate_pvalue,
        "tabix_bin": application.resources.executables.tabix,
        "window_kb": module.window_kb,
        "missing_index_action": module.missing_index_action,
        "reference_file_pattern": module.reference.file_pattern,
        "reference_index_suffix": module.reference.index_suffix,
        "summary_pvalue_thresholds": module.summary_pvalue_thresholds,
        "bcftools_bin": application.resources.executables.bcftools,
        "threads": application.execution.threads,
        "memory_gb": application.execution.memory_gb,
    }
    unknown = sorted(set(provided).difference(defaults))
    if unknown:
        raise RuntimeError(
            "Unknown LD-clumping runtime option(s): %s" % ", ".join(unknown)
        )
    return {
        name: defaults[name] if value is None else value
        for name, value in provided.items()
    }


def canonical_variant_id(chrom, pos, allele_1, allele_2):
    """Create one allele-order-independent ID for summary and LD variants."""
    if any(value is None for value in (chrom, pos, allele_1, allele_2)):
        return None
    chrom = str(chrom).upper()
    if chrom.startswith("CHR"):  # noqa: FURB188 - support Python before 3.9
        chrom = chrom[3:]
    if chrom.isdigit():
        chrom = str(int(chrom))
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        return None
    alleles = sorted((str(allele_1).upper(), str(allele_2).upper()))
    if any(allele in {"", ".", "NA"} for allele in alleles):
        return None
    return f"{chrom}_{pos}_{alleles[0]}_{alleles[1]}"


def normalize_ld_id(ld_id):
    """Normalize chr:pos:a1:a2 or chr_pos_a1_a2 LD identifiers."""
    if ld_id is None:
        return None
    parts = str(ld_id).replace(":", "_").split("_")
    if len(parts) != 4:
        return None
    return canonical_variant_id(parts[0], parts[1], parts[2], parts[3])


def add_canonical_ids(gwas):
    """Add canonical IDs and retain the lowest-P row for true duplicates."""
    required = {"chrcol", "poscol", "neacol", "eacol", "rsIDcol", "pcol"}
    missing = sorted(required.difference(gwas.columns))
    if missing:
        raise PipelineStageError(
            "03 canonical ID preparation",
            "add_canonical_ids",
            f"Missing required columns: {', '.join(missing)}",
        )

    chromosome_text = (
        pl.col("chrcol").cast(pl.Utf8).str.to_uppercase().str.replace(r"^CHR", "")
    )
    chrom = (
        pl.when(chromosome_text.str.contains(r"^\d+$"))
        .then(chromosome_text.cast(pl.Int64, strict=False).cast(pl.Utf8))
        .otherwise(chromosome_text)
    )
    ea = pl.col("eacol").cast(pl.Utf8).str.to_uppercase()
    nea = pl.col("neacol").cast(pl.Utf8).str.to_uppercase()
    allele_1 = pl.when(ea <= nea).then(ea).otherwise(nea)
    allele_2 = pl.when(ea <= nea).then(nea).otherwise(ea)
    gwas = gwas.with_row_index("_input_order").with_columns(
        chrom.alias("chrcol"),
        pl.col("rsIDcol").cast(pl.Utf8).alias("input_id"),
        pl.col("pcol").cast(pl.Float64, strict=False).alias("pcol"),
        ea.alias("_effect_allele_orientation"),
        nea.alias("_non_effect_allele_orientation"),
        pl.concat_str(
            [
                chrom,
                pl.col("poscol").cast(pl.Int64).cast(pl.Utf8),
                allele_1,
                allele_2,
            ],
            separator="_",
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
        print(
            "[WARNING] [STAGE: 03 canonical ID preparation] "
            "[FUNCTION: add_canonical_ids] FUMA-style duplicate resolution "
            "retained the lowest-P row for each canonical ID "
            f"| duplicate_canonical_ids={duplicates.height} "
            f"| rows_removed={rows_removed} | examples={examples}"
        )
    return (
        gwas.sort(["chrcol", "poscol", "pcol", "_input_order"])
        .unique("uniq_id", keep="first", maintain_order=True)
        .drop(
            "_input_order",
            "_effect_allele_orientation",
            "_non_effect_allele_orientation",
        )
    )


def build_id_aliases(gwas):
    """Map canonical IDs and unambiguous input rsIDs to canonical IDs."""
    aliases = {}
    ambiguous = set()
    for row in gwas.select("uniq_id", "rsIDcol").unique().iter_rows(named=True):
        canonical_id = row["uniq_id"]
        for alias in (canonical_id, row["rsIDcol"]):
            if alias is None or str(alias) in {"", ".", "NA"}:
                continue
            alias = str(alias)
            if alias in aliases and aliases[alias] != canonical_id:
                ambiguous.add(alias)
            else:
                aliases[alias] = canonical_id
    for alias in ambiguous:
        aliases.pop(alias, None)
    return aliases


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
    formatted = (
        extracted.with_columns((10.0 ** (-pl.col("lpcol"))).alias("pcol"))
        .drop("lpcol")
        .select(
            "chrcol",
            "poscol",
            "neacol",
            "eacol",
            "rsIDcol",
            "pcol",
            "becol",
            "secol",
            "eafcol",
        )
    )
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
    return str(output_file)


def get_ld_partners(
    chrom,
    pos,
    canonical_id,
    ld_file_path,
    r2_threshold,
    id_aliases=None,
    log=None,
    *,
    tabix_bin=None,
    window_kb=None,
    missing_index_action=None,
):
    """Return LD partners from a symmetric first-endpoint tabix reference."""
    resolved = _resolve_standard_options(
        tabix_bin=tabix_bin,
        window_kb=window_kb,
        missing_index_action=missing_index_action,
    )
    tabix_bin = resolved["tabix_bin"]
    window_kb = resolved["window_kb"]
    missing_index_action = resolved["missing_index_action"]
    id_aliases = id_aliases or {}
    emit = log or print
    region = None
    try:
        self_snp = pl.DataFrame(
            {
                "uniq_id": [canonical_id],
                "ld_id": [canonical_id],
                "ref_chr": [chrom],
                "ref_pos": [pos],
                "r2": [1.0],
            }
        )
        region = f"{chrom}:{pos}-{pos}"
        result = subprocess.run(
            [tabix_bin, ld_file_path, region],
            capture_output=True,
            text=True,
            check=True,
        )
        if not result.stdout.strip():
            message = (
                "No LD rows were returned for the exact index-SNP position "
                f"| chromosome={chrom} | position={pos} "
                f"| expected_id={canonical_id} | ld_rows_returned=0 "
                f"| ld_file={ld_file_path} | region={region}"
            )
            if missing_index_action == "self_only":
                emit(
                    "[WARNING] [STAGE: 04 LD retrieval] "
                    "[FUNCTION: get_ld_partners] %s | action=self-only "
                    "| consequence=independence assumed explicitly" % message
                )
                return self_snp
            raise PipelineStageError(
                "04 LD retrieval",
                "get_ld_partners",
                message + " | action=error",
                chromosome=chrom,
                variant=canonical_id,
                ld_file=ld_file_path,
                region=region,
            )
        ld_df = pl.read_csv(
            StringIO(result.stdout),
            separator="\t",
            has_header=False,
            new_columns=["chr_a", "pos_a", "snp_a", "chr_b", "pos_b", "snp_b", "r2"],
        ).with_columns(
            pl.col("pos_a").cast(pl.Int64),
            pl.col("pos_b").cast(pl.Int64),
            pl.col("snp_a").cast(pl.Utf8),
            pl.col("snp_b").cast(pl.Utf8),
            pl.col("r2").cast(pl.Float64),
        )

        def resolve_id(raw_id):
            raw_id = str(raw_id)
            return id_aliases.get(raw_id) or normalize_ld_id(raw_id) or raw_id

        ld_df = ld_df.with_columns(
            pl.Series("match_a", [resolve_id(x) for x in ld_df["snp_a"]]),
            pl.Series("match_b", [resolve_id(x) for x in ld_df["snp_b"]]),
        )
        matched = ld_df.filter(pl.col("match_a") == canonical_id)
        if matched.is_empty():
            observed_index_ids = (
                ld_df["snp_a"].unique(maintain_order=True).head(5).to_list()
            )
            message = (
                "Reference ID/allele mismatch at the queried index-SNP position "
                f"| chromosome={chrom} | position={pos} "
                f"| expected_id={canonical_id} | ld_rows_returned={ld_df.height} "
                f"| observed_index_ids={observed_index_ids} "
                f"| ld_file={ld_file_path} | region={region}"
            )
            if missing_index_action == "self_only":
                emit(
                    "[WARNING] [STAGE: 04 LD retrieval] "
                    "[FUNCTION: get_ld_partners] %s | action=self-only "
                    "| consequence=independence assumed explicitly" % message
                )
                return self_snp
            raise PipelineStageError(
                "04 LD retrieval",
                "get_ld_partners",
                message + " | action=error",
                chromosome=chrom,
                variant=canonical_id,
                ld_file=ld_file_path,
                region=region,
            )
        window_bp = int(window_kb) * 1000
        filtered = matched.filter(
            (pl.col("r2") >= r2_threshold)
            & ((pl.col("pos_b") - int(pos)).abs() <= window_bp)
        )
        if not matched.is_empty() and filtered.is_empty():
            emit(
                "[INFO] [STAGE: 04 LD retrieval] "
                "[FUNCTION: get_ld_partners] The index SNP matched the LD "
                "reference, but no partners passed the r2 threshold "
                f"| chromosome={chrom} | position={pos} "
                f"| expected_id={canonical_id} | matched_ld_rows={matched.height} "
                f"| r2_threshold={r2_threshold} | action=self-only "
                "| consequence=no LD partners assigned "
                f"| ld_file={ld_file_path} | region={region}"
            )
        partners = filtered.select(
            [
                pl.col("match_b").alias("uniq_id"),
                pl.col("snp_b").alias("ld_id"),
                pl.lit(chrom).alias("ref_chr"),
                pl.col("pos_b").alias("ref_pos"),
                pl.col("r2"),
            ]
        )
        return (
            pl.concat([self_snp, partners], how="vertical_relaxed")
            .sort("r2", descending=True)
            .unique("uniq_id", keep="first", maintain_order=True)
        )
    except FileNotFoundError as error:
        raise PipelineStageError(
            "04 LD retrieval",
            "get_ld_partners",
            "tabix executable was not found",
            chromosome=chrom,
            variant=canonical_id,
            ld_file=ld_file_path,
            region=region,
        ) from error
    except subprocess.CalledProcessError as error:
        stderr = (error.stderr or "").strip()
        detail = f"tabix returned exit code {error.returncode}"
        if stderr:
            detail = f"{detail}: {stderr}"
        raise PipelineStageError(
            "04 LD retrieval",
            "get_ld_partners",
            detail,
            chromosome=chrom,
            variant=canonical_id,
            ld_file=ld_file_path,
            region=region,
        ) from error
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "04 LD retrieval",
            "get_ld_partners",
            error,
            chromosome=chrom,
            variant=canonical_id,
            ld_file=ld_file_path,
            region=region,
        ) from error


def find_ind_sig_snps(
    chr_df,
    ld_path,
    lead_p_threshold,
    r2_clump_threshold,
    threads=None,
    log=None,
    *,
    candidate_p_threshold=None,
    tabix_bin=None,
    window_kb=None,
    missing_index_action=None,
):
    del threads
    resolved = _resolve_standard_options(
        candidate_p_threshold=candidate_p_threshold,
        tabix_bin=tabix_bin,
        window_kb=window_kb,
        missing_index_action=missing_index_action,
    )
    candidate_p_threshold = resolved["candidate_p_threshold"]
    tabix_bin = resolved["tabix_bin"]
    window_kb = resolved["window_kb"]
    missing_index_action = resolved["missing_index_action"]
    emit = log or print
    remaining_leads = chr_df.filter(pl.col("pcol") <= lead_p_threshold).sort("pcol")
    ind_sig_clumps = []
    id_aliases = build_id_aliases(chr_df)
    chrom = chr_df.select(pl.col("chrcol").first()).item()
    emit(
        "[INFO] [STAGE: 05 independent SNP clumping] "
        "[FUNCTION: find_ind_sig_snps] Started "
        f"| chromosome={chrom} | significance_threshold={lead_p_threshold} "
        f"| significant_variants={remaining_leads.height}"
    )
    while not remaining_leads.is_empty():
        variants_before = remaining_leads.height
        top_snp = remaining_leads.row(0, named=True)
        sid = top_snp["uniq_id"]
        partners = get_ld_partners(
            top_snp["chrcol"],
            top_snp["poscol"],
            sid,
            ld_path,
            r2_clump_threshold,
            id_aliases,
            log=emit,
            tabix_bin=tabix_bin,
            window_kb=window_kb,
            missing_index_action=missing_index_action,
        )
        members = (
            partners.join(chr_df, on="uniq_id", how="left")
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
        members = members.with_columns(
            pl.lit(members["uniq_id"].n_unique()).alias("n_refsnps")
        )
        ind_sig_clumps.append(members)

        # Only deplete the significant-SNP pool. FUMA candidates may be linked
        # to more than one independent significant SNP.
        partner_ids = partners["uniq_id"].to_list()
        remaining_leads = remaining_leads.filter(~pl.col("uniq_id").is_in(partner_ids))
        variants_removed = variants_before - remaining_leads.height
        emit(
            "[INFO] [STAGE: 05 independent SNP clumping] "
            "[FUNCTION: find_ind_sig_snps] Index SNP processed "
            f"| chromosome={chrom} | selected_variant={sid} "
            f"| selected_p={top_snp['pcol']} "
            f"| significant_variants_removed={variants_removed} "
            f"| significant_variants_remaining={remaining_leads.height}"
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
    tabix_bin=None,
    window_kb=None,
    missing_index_action=None,
):
    del threads
    if ind_sig_df.is_empty():
        return pl.DataFrame()
    resolved = _resolve_standard_options(
        tabix_bin=tabix_bin,
        window_kb=window_kb,
        missing_index_action=missing_index_action,
    )
    tabix_bin = resolved["tabix_bin"]
    window_kb = resolved["window_kb"]
    missing_index_action = resolved["missing_index_action"]
    emit = log or print
    ind_sig_heads = ind_sig_df.filter(
        pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
    ).sort("pcol")
    chrom = ind_sig_heads.select(pl.col("chrcol").first()).item()
    lead_clusters = []
    remaining = ind_sig_heads
    id_aliases = build_id_aliases(ind_sig_heads)
    emit(
        "[INFO] [STAGE: 06 lead SNP clumping] "
        "[FUNCTION: find_lead_snps] Started "
        f"| chromosome={chrom} | r2_threshold={r2_lead_threshold} "
        f"| independent_significant_snps={remaining.height}"
    )
    while not remaining.is_empty():
        top_row = remaining.row(0, named=True)
        lid = top_row["uniq_id"]
        partners = get_ld_partners(
            top_row["chrcol"],
            top_row["poscol"],
            lid,
            ld_path,
            r2_lead_threshold,
            id_aliases,
            log=emit,
            tabix_bin=tabix_bin,
            window_kb=window_kb,
            missing_index_action=missing_index_action,
        ).select("uniq_id", pl.col("r2").alias("r2_with_Lead"))
        members = remaining.join(partners, on="uniq_id", how="inner").with_columns(
            pl.lit(lid).alias("lead_SNP_id")
        )
        lead_clusters.append(members)
        assigned_ids = members["uniq_id"].to_list()
        remaining = remaining.filter(~pl.col("uniq_id").is_in(assigned_ids))
        emit(
            "[INFO] [STAGE: 06 lead SNP clumping] "
            "[FUNCTION: find_lead_snps] Lead SNP processed "
            f"| chromosome={chrom} | selected_variant={lid} "
            f"| independent_snps_assigned={members.height} "
            f"| independent_snps_remaining={remaining.height}"
        )
    emit(
        "[INFO] [STAGE: 06 lead SNP clumping] "
        "[FUNCTION: find_lead_snps] Completed "
        f"| chromosome={chrom} | lead_snps={len(lead_clusters)}"
    )
    return pl.concat(lead_clusters) if lead_clusters else pl.DataFrame()


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

    initial = is_clusters.join(lead_map, on="ind_sig_SNP_id", how="left")
    initial = (
        initial.group_by("lead_SNP_id")
        .agg(
            pl.col("chr").first(),
            pl.col("start").min(),
            pl.col("end").max(),
            pl.col("p").min(),
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
            if locus["p"] < curr["p"]:
                for column in [
                    "lead_SNP_id",
                    "p",
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
                "LP": -math.log10(m["p"]) if m["p"] > 0 else 300.0,
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

    hierarchy_base = ind_sig_df.join(lead_map, on="ind_sig_SNP_id", how="left").join(
        ind_locus_map, on="ind_sig_SNP_id", how="left"
    )
    linked_ind_sigs = hierarchy_base.group_by(["Genomic_locus", "uniq_id"]).agg(
        pl.col("ind_sig_SNP_id").unique().alias("Linked_IndSigSNPs")
    )
    hierarchy = (
        hierarchy_base.sort("r2_with_IndSig", descending=True)
        .unique(["Genomic_locus", "uniq_id"], keep="first", maintain_order=True)
        .join(linked_ind_sigs, on=["Genomic_locus", "uniq_id"], how="left")
        .sort(["Genomic_locus", "poscol"])
    )

    is_clusters = (
        is_clusters.join(ind_locus_map, on="ind_sig_SNP_id", how="left")
        .join(lead_map, on="ind_sig_SNP_id", how="left")
        .sort(["Genomic_locus", "start"])
    )
    initial = initial.join(lead_locus_map, on="lead_SNP_id", how="left").sort(
        ["Genomic_locus", "start"]
    )
    return summary_df, hierarchy, is_clusters, initial, idx


def process_chromosome(
    chrom,
    gwas_subset,
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
    reference_file_pattern=None,
    reference_index_suffix=None,
    summary_pvalue_thresholds=None,
):
    """
    Worker function to process a single chromosome.
    Accepts the canonical thread budget for sub-functions.
    """
    resolved = _resolve_standard_options(
        candidate_p=candidate_p,
        tabix_bin=tabix_bin,
        window_kb=window_kb,
        missing_index_action=missing_index_action,
        reference_file_pattern=reference_file_pattern,
        reference_index_suffix=reference_index_suffix,
        summary_pvalue_thresholds=summary_pvalue_thresholds,
    )
    candidate_p = resolved["candidate_p"]
    tabix_bin = resolved["tabix_bin"]
    window_kb = resolved["window_kb"]
    missing_index_action = resolved["missing_index_action"]
    reference_file_pattern = resolved["reference_file_pattern"]
    reference_index_suffix = resolved["reference_index_suffix"]
    summary_pvalue_thresholds = resolved["summary_pvalue_thresholds"]
    messages = []
    emit = messages.append
    significant_variants = gwas_subset.filter(pl.col("pcol") <= lead_p).height
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
    # Step A: Find independent significant SNPs.
    try:
        ind_sig = find_ind_sig_snps(
            gwas_subset,
            str(ld_path),
            lead_p,
            r2_clump,
            log=emit,
            candidate_p_threshold=candidate_p,
            tabix_bin=tabix_bin,
            window_kb=window_kb,
            missing_index_action=missing_index_action,
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
            "[FUNCTION: process_chromosome] Chromosome skipped because no "
            "variants passed the significance threshold "
            f"| chromosome={chrom} | significance_threshold={lead_p} "
            f"| significant_variants=0 | status=normal"
        )
        return {
            "chrom": chrom,
            "logs": messages,
            "progress": {
                "status": "no_significant_variants",
                "significant": significant_variants,
                "independent": 0,
                "lead": 0,
                "loci": 0,
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
            tabix_bin=tabix_bin,
            window_kb=window_kb,
            missing_index_action=missing_index_action,
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
            "progress": {
                "status": "no_lead_variants",
                "significant": significant_variants,
                "independent": ind_sig.filter(
                    pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
                ).height,
                "lead": 0,
                "loci": 0,
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
        "progress": {
            "status": "completed",
            "significant": significant_variants,
            "independent": ind_sig.filter(
                pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
            ).height,
            "lead": leads["lead_SNP_id"].n_unique(),
            "loci": summ.height,
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


def _normalise_chromosome_label(chrom):
    chromosome = str(chrom).upper()
    if chromosome.startswith("CHR"):
        chromosome = chromosome[3:]
    return str(int(chromosome)) if chromosome.isdigit() else chromosome


def _chromosome_sort_key(chrom):
    chromosome = _normalise_chromosome_label(chrom)
    return (0, int(chromosome)) if chromosome.isdigit() else (1, chromosome)


def _compact_chromosome_labels(chromosomes):
    """Format observed chromosome labels as compact, ordered ranges."""
    labels = sorted(
        {_normalise_chromosome_label(chrom) for chrom in chromosomes},
        key=_chromosome_sort_key,
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


def _print_chromosome_progress(chrom, progress, warning_count):
    """Print only chromosome outcomes that require live user attention."""
    status = progress["status"]
    if status == "no_significant_variants" and warning_count == 0:
        return
    outcome = (
        f"{progress['significant']:,} significant · "
        f"{progress['independent']:,} independent · "
        f"{progress['lead']:,} lead · {progress['loci']:,} "
        f"{'risk locus' if progress['loci'] == 1 else 'risk loci'}"
    )
    if warning_count:
        outcome += (
            f" · {warning_count:,} "
            f"{'warning' if warning_count == 1 else 'warnings'}"
        )
    kind = "warning" if warning_count or progress["loci"] == 0 else "success"
    print(
        screen_field(
            kind,
            f"chr{_normalise_chromosome_label(chrom)}",
            outcome,
            indent=6,
            label_width=8,
        ),
        flush=True,
    )


def _print_clumping_summary(chromosome_progress, detailed_log):
    """Print one concise end-of-run summary while the log retains full detail."""
    no_significant = [
        progress["chrom"]
        for progress in chromosome_progress
        if progress["status"] == "no_significant_variants"
    ]
    with_loci = [
        progress["chrom"]
        for progress in chromosome_progress
        if progress["loci"] > 0
    ]
    significant_without_loci = [
        progress["chrom"]
        for progress in chromosome_progress
        if progress["status"] != "no_significant_variants"
        and progress["loci"] == 0
    ]
    warning_chromosomes = [
        progress["chrom"]
        for progress in chromosome_progress
        if progress["warnings"] > 0
    ]
    totals = {
        field: sum(progress[field] for progress in chromosome_progress)
        for field in ("significant", "independent", "lead", "loci", "warnings")
    }
    label_width = 29
    print()
    print(screen_line("analysis", "Standard LD clumping completed", indent=2))
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
        ("count", "Significant variants", "significant"),
        ("genetic", "Independent significant SNPs", "independent"),
        ("genetic", "Lead SNPs", "lead"),
        ("genetic", "Genomic risk loci", "loci"),
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
    warning_summary = f"{totals['warnings']:,}"
    if warning_chromosomes:
        warning_summary += (
            f" across {_compact_chromosome_labels(warning_chromosomes)}"
        )
    print(
        screen_field(
            "warning" if totals["warnings"] else "info",
            "Warnings",
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
    configuration=None,
    logger=None,
):
    application = None
    if configuration is None:
        application = load_configuration()
        configuration = application.modules.ld_clumping
    pop = pop or configuration.population.value
    lead_p = configuration.lead_pvalue if lead_p is None else lead_p
    r2_clump = configuration.clump_r2 if r2_clump is None else r2_clump
    r2_lead = configuration.lead_r2 if r2_lead is None else r2_lead
    merge_dist = (
        configuration.merge_distance_bp if merge_dist is None else merge_dist
    )
    resolved = _resolve_standard_options(
        application=application,
        bcftools_bin=bcftools_bin,
        tabix_bin=tabix_bin,
        threads=threads,
        memory_gb=memory_gb,
    )
    bcftools_bin = resolved["bcftools_bin"]
    tabix_bin = resolved["tabix_bin"]
    threads = resolved["threads"]
    memory_gb = resolved["memory_gb"]

    # 1. Conversion
    try:
        tsv_path = vcf_to_standard_ldclump(
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
    # 2. Load Data
    try:
        gwas = pl.read_csv(
            tsv_path, separator=configuration.table.delimiter,
        )
    except Exception as error:
        raise function_error(
            "02 summary loading",
            "ld_clump_standard",
            error,
            input_file=tsv_path,
        ) from error
    try:
        gwas = add_canonical_ids(gwas)
    except PipelineStageError:
        raise
    except Exception as error:
        raise function_error(
            "03 canonical ID preparation",
            "add_canonical_ids",
            error,
            input_file=tsv_path,
        ) from error
    before_mhc = gwas.height
    if configuration.remove_mhc:
        mhc = configuration.mhc_regions[configuration.genome_build]
        normalized_chromosome = (
            pl.col("chrcol")
            .cast(pl.Utf8)
            .str.to_uppercase()
            .str.replace(r"^CHR", "")
        )
        mhc_chromosome = str(mhc.chromosome).upper().removeprefix("CHR")
        gwas = gwas.filter(
            ~(
                (normalized_chromosome == mhc_chromosome)
                & pl.col("poscol").is_between(mhc.start, mhc.end, closed="both")
            )
        )
    Path(detailed_log).parent.mkdir(parents=True, exist_ok=True)
    with Path(detailed_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "MHC exclusion: enabled=%s genome_build=%s rows_in=%d rows_out=%d "
            "removed=%d\n"
            % (
                configuration.remove_mhc,
                configuration.genome_build.value,
                before_mhc,
                gwas.height,
                before_mhc - gwas.height,
            )
        )
    if logger is not None:
        logger.record(
            "ACTION",
            "standard_mhc_exclusion",
            enabled=configuration.remove_mhc,
            rows_in=before_mhc,
            rows_out=gwas.height,
            removed=before_mhc - gwas.height,
        )
    significant_chromosomes = (
        gwas.filter(pl.col("pcol") <= lead_p)["chrcol"].unique().to_list()
    )
    missing_resources = []
    for chromosome in significant_chromosomes:
        reference_path = Path(ld_folder) / configuration.reference.file_pattern.format(
            population=pop,
            chromosome=chromosome,
        )
        index_path = Path(str(reference_path) + configuration.reference.index_suffix)
        if not reference_path.is_file() or reference_path.stat().st_size <= 0:
            missing_resources.append(str(reference_path))
        if not index_path.is_file() or index_path.stat().st_size <= 0:
            missing_resources.append(str(index_path))
    if missing_resources:
        raise PipelineStageError(
            "04 LD reference validation",
            "ld_clump_standard",
            "Missing or empty LD resources for chromosomes containing significant "
            "variants: %s" % ", ".join(missing_resources),
            sample=dataset_id,
        )
    # --- BALANCED SCHEDULING ---
    try:
        chrom_stats = gwas.group_by("chrcol").len().sort("len", descending=True)
        balanced_chroms = []
        big_to_small = chrom_stats["chrcol"].to_list()
        while big_to_small:
            balanced_chroms.append(big_to_small.pop(0))
            if big_to_small:
                balanced_chroms.append(big_to_small.pop(-1))
        gwas_mem_gb = gwas.estimated_size() / (1024**3)
        memory_per_worker_gb = max(
            configuration.compute.minimum_worker_memory_gb,
            gwas_mem_gb * configuration.compute.input_memory_multiplier,
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
    effective_workers = max(1, min(n_workers, len(balanced_chroms)))
    print()
    print(screen_line("analysis", "Standard LD clumping", indent=2))
    print(
        screen_field(
            "count",
            "Chromosomes scheduled",
            f"{len(balanced_chroms):,}",
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
                    gwas.filter(pl.col("chrcol") == chrom),
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
                    reference_file_pattern=configuration.reference.file_pattern,
                    reference_index_suffix=configuration.reference.index_suffix,
                    summary_pvalue_thresholds=(
                        configuration.summary_pvalue_thresholds
                    ),
                )
                for chrom in balanced_chroms
            }
            ordered_chromosomes = sorted(futures, key=_chromosome_sort_key)
            for chrom in ordered_chromosomes:
                future = futures[chrom]
                try:
                    result = future.result()
                    if result:
                        messages = result.pop("logs", [])
                        _append_chromosome_diagnostics(
                            detailed_log, chrom, messages
                        )
                        warning_count = sum(
                            message.startswith("[WARNING]") for message in messages
                        )
                        progress = result.pop("progress")
                        chromosome_progress.append(
                            {
                                "chrom": chrom,
                                "warnings": warning_count,
                                **progress,
                            }
                        )
                        _print_chromosome_progress(
                            chrom,
                            progress,
                            warning_count,
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
                            f"chr{_normalise_chromosome_label(chrom)} failed",
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
    # 4. Final Aggregation
    if res_list:
        def chr_key(x):
            c = str(x["chrom"]).lower().replace("chr", "")
            return (0, int(c)) if c.isdigit() else (1, c)

        try:
            res_list.sort(key=chr_key)
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
        except Exception as error:
            raise function_error(
                "09 result aggregation",
                "ld_clump_standard",
                error,
                dataset=dataset_id,
            ) from error
        for key, pattern in {
            "summ": configuration.output_layout.standard_summary,
            "hier": configuration.output_layout.standard_hierarchy,
            "is_c": configuration.output_layout.standard_independent_clusters,
            "l_un": configuration.output_layout.standard_lead_clusters,
        }.items():
            output_file = configured_output_path(
                output_directory,
                pattern,
                error_type=RuntimeError,
                dataset_id=dataset_id,
                population=pop,
            )
            output_file.parent.mkdir(parents=True, exist_ok=True)
            try:
                df = pl.concat(final_res[key])
                # Ensure list types are converted to strings before CSV write
                df = df.select(
                    pl.col(c).list.join("; ") if df[c].dtype.is_nested() else pl.col(c)
                    for c in df.columns
                )
            except Exception as error:
                raise function_error(
                    "09 result aggregation",
                    "ld_clump_standard",
                    error,
                    result=key,
                    output_file=output_file,
                ) from error
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
        _print_clumping_summary(chromosome_progress, detailed_log)
        summary_path = configured_output_path(
            output_directory,
            configuration.output_layout.standard_summary,
            error_type=RuntimeError,
            dataset_id=dataset_id,
            population=pop,
        )
        result = {
            "status": "completed",
            "ldpruned_sig_file": str(summary_path),
            "log_file": str(detailed_log),
            "significant_variants": sum(
                item["significant"] for item in chromosome_progress
            ),
            "genomic_risk_loci": sum(item["loci"] for item in chromosome_progress),
        }
        if logger is not None:
            logger.record(
                "OUTPUT", "standard_clumping", **result,
                allele_orientation=(
                    "REF/ALT retained; canonical matching ID is allele-order-independent"
                ),
            )
        return result
    _print_clumping_summary(chromosome_progress, detailed_log)
    result = {
        "status": "no_loci",
        "ldpruned_sig_file": None,
        "log_file": str(detailed_log),
        "significant_variants": sum(
            item["significant"] for item in chromosome_progress
        ),
        "genomic_risk_loci": 0,
    }
    if logger is not None:
        logger.record("RESULT", "standard_clumping_no_loci", **result)
    return result
