"""FINEMAP workflow with strict allele, LD, and FLAMES validation."""

import logging
import math
import multiprocessing as mp
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import time

import numpy as np
import pandas as pd
import polars as pl

from postgwas.cli.compute import resolve_compute_args
from postgwas.core.execution.runtime import safe_thread_count
from postgwas.modules.fine_mapping.defaults import (
    DEFAULT_BASES_PER_KILOBASE,
    DEFAULT_BGEN_BITS,
    DEFAULT_COVERAGE_TOLERANCE,
    DEFAULT_CREDIBLE_SET_COVERAGE,
    DEFAULT_EXTERNAL_TOOL_THREADS,
    DEFAULT_FALLBACK_MEMORY_GB,
    DEFAULT_FINEMAP_RAM_PER_WORKER_GB,
    DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL,
    DEFAULT_GENOME_BUILD,
    DEFAULT_LD_CORRELATION_TOLERANCE,
    DEFAULT_LD_DIAGONAL_TOLERANCE,
    DEFAULT_LD_EIGENVALUE_TOLERANCE,
    DEFAULT_LD_MAX_CORRELATION,
    DEFAULT_LD_SYMMETRY_TOLERANCE,
    DEFAULT_MAX_MAF,
    DEFAULT_MIN_POSITION,
    DEFAULT_PLINK_MEMORY_MB,
    DEFAULT_SCHEMA_INFERENCE_LENGTH,
    SUPPORTED_GENOME_BUILDS,
    get_finemap_defaults,
)
from postgwas.modules.fine_mapping.engines.finemap.input_gen import (
    create_finemap_master,
    create_ldstore_master,
    write_finemap_z_file,
    write_snp_file,
)
from postgwas.modules.fine_mapping.engines.finemap.merge_results import (
    INVALID_MODEL_PROBABILITY_REASON,
    InvalidModelProbabilityError,
    finemap_credible_files,
    parse_cred_header,
    process_finemap_output,
)
from postgwas.modules.fine_mapping.engines.finemap.runner import (
    FinemapExternalToolTimeout,
    run_bgen_indexing,
    run_finemap_binary,
    run_ldstore,
    run_plink_extraction,
    run_plink_to_bgen,
)
from postgwas.modules.fine_mapping.progress import ProgressRecorder, write_run_configuration
from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths
from postgwas.modules.fine_mapping.preflight import _validate_runtime_tools
from postgwas.modules.fine_mapping.resource_guard import (
    VARIANT_LIMIT_FAILURE_REASON,
    combine_reasons,
    dense_ld_resource_qc,
    variant_limit_failure_detail,
)


logger = logging.getLogger("postgwas.modules.fine_mapping")


def _normalise_chromosome(values):
    normalised = (
        values.astype(str)
        .str.strip()
        .str.replace(r"^chr", "", regex=True, case=False)
        .str.replace(r"\.0$", "", regex=True)
        .str.upper()
    )
    return normalised.replace({"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"})


def _complement(allele):
    allele = str(allele).upper()
    if not allele or any(base not in "ACGT" for base in allele):
        return None
    return allele.translate(str.maketrans("ACGT", "TGCA"))


def _read_bim(bim_file):
    columns = ["ref_chromosome", "rsid", "cm", "ref_position", "ref_allele1", "ref_allele2"]
    ref_df = pd.read_csv(
        bim_file,
        sep=r"\s+",
        header=None,
        names=columns,
        dtype={"ref_chromosome": str, "rsid": str, "ref_allele1": str, "ref_allele2": str},
    )
    if ref_df.empty:
        raise ValueError(f"Reference BIM is empty: {bim_file}")
    ref_df["ref_chromosome"] = _normalise_chromosome(ref_df["ref_chromosome"])
    ref_df["ref_position"] = pd.to_numeric(ref_df["ref_position"], errors="coerce")
    ref_df["ref_allele1"] = ref_df["ref_allele1"].str.upper()
    ref_df["ref_allele2"] = ref_df["ref_allele2"].str.upper()
    if ref_df["rsid"].duplicated().any():
        duplicated = int(ref_df["rsid"].duplicated(keep=False).sum())
        raise ValueError(f"Reference BIM contains {duplicated} rows with duplicate variant IDs")
    return ref_df


def _harmonise_sumstats(sumstats, ref_df, max_maf=DEFAULT_MAX_MAF):
    """Align effect alleles and beta signs to BIM A1, retaining auditable QC counts."""
    required = [
        "rsid", "chromosome", "position", "allele1", "allele2",
        "maf", "beta", "se", "NEF",
    ]
    missing = [column for column in required if column not in sumstats.columns]
    if missing:
        raise ValueError(f"Summary statistics are missing required columns: {missing}")

    # Avoid an implicit pyarrow dependency: Polars' to_pandas() requires it
    # for many dtypes, while this conversion is sufficient for validated
    # scalar summary-statistic columns.
    frame = pd.DataFrame(sumstats.to_dicts()) if isinstance(sumstats, pl.DataFrame) else sumstats.copy()
    if frame.empty:
        raise ValueError("Summary-statistics file contains no variants")
    frame["rsid"] = frame["rsid"].astype(str).str.strip()
    if frame["rsid"].duplicated().any():
        duplicated = int(frame["rsid"].duplicated(keep=False).sum())
        raise ValueError(f"Summary statistics contain {duplicated} rows with duplicate variant IDs")

    frame["chromosome"] = _normalise_chromosome(frame["chromosome"])
    frame["position"] = pd.to_numeric(frame["position"], errors="coerce")
    for column in ["maf", "beta", "se", "NEF"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["allele1"] = frame["allele1"].astype(str).str.upper().str.strip()
    frame["allele2"] = frame["allele2"].astype(str).str.upper().str.strip()
    frame["maf"] = np.minimum(frame["maf"], 1.0 - frame["maf"])

    valid_numeric = (
        frame["position"].notna()
        & np.isfinite(frame["beta"])
        & np.isfinite(frame["se"])
        & (frame["se"] > 0)
        & np.isfinite(frame["NEF"])
        & (frame["NEF"] > 0)
        & np.isfinite(frame["maf"])
        & (frame["maf"] > 0)
        & (frame["maf"] <= float(max_maf))
        & frame["allele1"].ne("")
        & frame["allele2"].ne("")
    )
    n_invalid = int((~valid_numeric).sum())
    frame = frame.loc[valid_numeric].copy()

    joined = frame.merge(ref_df, on="rsid", how="inner", validate="one_to_one")
    n_id_match = len(joined)
    position_match = (
        joined["chromosome"].eq(joined["ref_chromosome"])
        & joined["position"].eq(joined["ref_position"])
    )

    effect = joined["allele1"]
    other = joined["allele2"]
    ref_effect = joined["ref_allele1"]
    ref_other = joined["ref_allele2"]
    comp_effect = effect.map(_complement)
    comp_other = other.map(_complement)

    direct = effect.eq(ref_effect) & other.eq(ref_other)
    swapped = effect.eq(ref_other) & other.eq(ref_effect)
    complement_direct = comp_effect.eq(ref_effect) & comp_other.eq(ref_other)
    complement_swapped = comp_effect.eq(ref_other) & comp_other.eq(ref_effect)
    palindromic = (effect + other).isin({"AT", "TA", "CG", "GC"})
    allele_match = (direct | swapped | complement_direct | complement_swapped) & ~palindromic
    keep = position_match & allele_match

    aligned = joined.loc[keep].copy()
    flip = (swapped | complement_swapped).loc[keep]
    aligned.loc[flip, "beta"] = -aligned.loc[flip, "beta"]
    aligned["chromosome"] = aligned["ref_chromosome"]
    aligned["position"] = aligned["ref_position"].astype(np.int64)
    aligned["allele1"] = aligned["ref_allele1"]
    aligned["allele2"] = aligned["ref_allele2"]
    aligned = aligned[required].sort_values(
        ["chromosome", "position", "allele1", "allele2", "rsid"], kind="stable"
    )
    if aligned.empty:
        raise ValueError(
            "No variants remained after reference ID, coordinate, and allele harmonisation"
        )
    variant_key = ["chromosome", "position", "allele1", "allele2"]
    if aligned.duplicated(variant_key).any():
        raise ValueError(
            "Multiple variant IDs map to the same chromosome, position, and allele pair"
        )

    qc = {
        "input_variants": int(len(sumstats)),
        "invalid_numeric_or_allele": n_invalid,
        "reference_id_matches": int(n_id_match),
        "coordinate_mismatches": int((~position_match).sum()),
        "ambiguous_palindromic": int((position_match & palindromic).sum()),
        "allele_mismatches": int((position_match & ~allele_match).sum()),
        "beta_sign_flips": int(flip.sum()),
        "retained_variants": int(len(aligned)),
    }
    return pl.DataFrame(aligned.to_dict(orient="list")), qc


def _prepare_loci(
    args,
    minimum_position=DEFAULT_MIN_POSITION,
    bases_per_kilobase=DEFAULT_BASES_PER_KILOBASE,
):
    loci_pd = pd.read_csv(args.locus_file, sep=None, engine="python")
    loci_df = pl.DataFrame(loci_pd.to_dict(orient="list"))
    loci_df = loci_df.rename({column: column.strip() for column in loci_df.columns})
    if loci_df.is_empty():
        raise ValueError(f"The input locus file is empty: {args.locus_file}")

    locus_type = str(args.locus_type).lower()
    required = ["CHR", "POS"] if locus_type == "point" else ["CHR", "START", "END"]
    required.append("LP")
    missing = [column for column in required if column not in loci_df.columns]
    if missing:
        raise ValueError(f"Locus file is missing required columns for {locus_type} mode: {missing}")

    expressions = [
        pl.col("CHR")
        .cast(pl.Utf8)
        .str.replace(r"(?i)^chr", "")
        .str.replace(r"\.0$", "")
        .str.to_uppercase()
        .replace({"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"})
        .alias("CHR"),
        pl.col("LP").cast(pl.Float64, strict=False).alias("LP"),
    ]
    if "P_value" in loci_df.columns:
        expressions.append(pl.col("P_value").cast(pl.Float64, strict=False).alias("P_value"))
    if locus_type == "point":
        expressions.append(pl.col("POS").cast(pl.Int64, strict=False).alias("POS"))
    else:
        expressions.extend([
            pl.col("START").cast(pl.Int64, strict=False).alias("START"),
            pl.col("END").cast(pl.Int64, strict=False).alias("END"),
        ])
    loci_df = loci_df.with_columns(expressions)

    flank_bp = int(float(args.window_kb) * int(bases_per_kilobase))
    if flank_bp < 0:
        raise ValueError("window_kb cannot be negative")
    if locus_type == "point":
        loci_df = loci_df.drop_nulls(subset=["CHR", "POS", "LP"]).with_columns([
            (pl.col("POS") - flank_bp).clip(lower_bound=int(minimum_position)).alias("START"),
            (pl.col("POS") + flank_bp).alias("END"),
        ])
    elif locus_type == "range":
        loci_df = loci_df.drop_nulls(subset=["CHR", "START", "END", "LP"])
        loci_df = loci_df.with_columns([
            (pl.col("START") - flank_bp).clip(lower_bound=int(minimum_position)).alias("START"),
            (pl.col("END") + flank_bp).alias("END"),
        ])
    else:
        raise ValueError("locus_type must be 'point' or 'range'")

    if loci_df.filter(
        (pl.col("START") < int(minimum_position)) | (pl.col("END") < pl.col("START"))
    ).height:
        raise ValueError("Locus file contains invalid coordinates")
    loci_df = loci_df.filter(pl.col("LP") >= float(args.lp_threshold))
    if loci_df.is_empty():
        raise ValueError(f"No loci passed LP >= {args.lp_threshold}")

    if getattr(args, "finemap_skip_mhc", False) and not getattr(args, "finemap_include_mhc", False):
        mhc_chrom = _normalise_chromosome(
            pd.Series([args.finemap_mhc_chromosome])
        ).iloc[0]
        loci_df = loci_df.filter(
            ~(
                (pl.col("CHR") == mhc_chrom)
                & (pl.col("START") < int(args.finemap_mhc_end))
                & (pl.col("END") > int(args.finemap_mhc_start))
            )
        )
        if loci_df.is_empty():
            raise ValueError("All loci were removed by MHC filtering")

    if "GenomicLocus" not in loci_df.columns:
        loci_df = loci_df.with_columns(
            pl.concat_str([
                pl.lit("chr"), pl.col("CHR"), pl.lit(":"),
                pl.col("START"), pl.lit("-"), pl.col("END"),
            ]).alias("GenomicLocus")
        )
    loci_df = loci_df.with_columns(pl.col("GenomicLocus").cast(pl.Utf8))
    if loci_df.select(pl.struct(["CHR", "START", "END"]).is_duplicated().any()).item():
        raise ValueError("Locus file contains duplicate CHR/START/END rows")
    if loci_df.select(pl.col("GenomicLocus").is_duplicated().any()).item():
        raise ValueError("GenomicLocus identifiers must be unique")
    return loci_df


def setup_directories(outdir_path, output_layout, dataset_id):
    """Create and return all pipeline-owned directories."""
    outdir = Path(outdir_path).resolve()
    output_paths = resolve_output_paths(outdir, output_layout, dataset_id)
    dirs = {
        **output_paths,
        "root": outdir,
        "temp": output_paths["finemap_inputs_directory"],
        "loci": output_paths["finemap_locus_work_directory"],
        "flames": output_paths["primary_flames_directory"],
        "inter": output_paths["finemap_selected_models_directory"],
    }
    for directory in [
        outdir,
        dirs["temp"],
        dirs["loci"],
        output_paths["quality_control_directory"],
        output_paths["run_metadata_directory"],
    ]:
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


def load_and_prep_inputs(
    args,
    qc_path,
    schema_inference_length=DEFAULT_SCHEMA_INFERENCE_LENGTH,
    max_maf=DEFAULT_MAX_MAF,
):
    """Load loci and strictly harmonise summary statistics to the LD reference."""
    logger.info("Loading and validating locus file: %s", args.locus_file)
    loci_df = _prepare_loci(args)

    logger.info("Loading summary statistics: %s", args.finemap_in_files)
    sumstats = pl.read_csv(
        args.finemap_in_files,
        separator="\t",
        null_values=["NA", ".", "nan", "NaN", "inf", "-inf"],
        infer_schema_length=int(schema_inference_length),
    )

    raw_ref = str(args.finemap_ld_reference)
    ld_ref_prefix = str(Path(raw_ref).with_suffix("")) if raw_ref.endswith((".bed", ".bim", ".fam")) else raw_ref
    for suffix in [".bed", ".bim", ".fam"]:
        if not Path(f"{ld_ref_prefix}{suffix}").is_file():
            raise FileNotFoundError(f"LD reference file is missing: {ld_ref_prefix}{suffix}")

    ref_df = _read_bim(f"{ld_ref_prefix}.bim")
    sumstats_filt, harmonisation_qc = _harmonise_sumstats(
        sumstats, ref_df, max_maf=max_maf
    )
    qc_path = Path(qc_path)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([harmonisation_qc]).to_csv(qc_path, sep="\t", index=False)
    logger.info(
        "Input ready: %d loci and %d harmonised variants",
        loci_df.height,
        sumstats_filt.height,
    )
    return loci_df, sumstats_filt, ld_ref_prefix


def _resolved_finemap_config(args, target_coverage):
    """Collect the already-resolved FINEMAP controls for tasks and audit."""
    config = {
        key: getattr(args, key)
        for key in [
            "algorithm", "n_causal_snps", "n_iter", "n_conv_sss",
            "prob_conv_sss_tol", "n_configs_top", "corr_config",
            "pvalue_snps", "cond_pvalue", "prior_std", "prior_k",
            "force_n_samples", "std_effects",
        ]
    }
    config["prob_cred_set"] = float(target_coverage)
    return config


def generate_tasks(
    loci_df,
    sumstats_df,
    ld_ref_prefix,
    dirs,
    args,
    task_skips_file,
    target_coverage=DEFAULT_CREDIBLE_SET_COVERAGE,
    runtime_defaults=None,
):
    """Create one validated FINEMAP task per locus."""
    finemap_config = _resolved_finemap_config(args, target_coverage)
    runtime_defaults = dict(runtime_defaults or get_finemap_defaults())

    sample_size_policy = args.sample_size_policy
    sample_size_summary_statistic = args.sample_size_summary_statistic
    warning_threshold = float(args.sample_size_relative_range_warning_threshold)
    maximum_variants = int(args.maximum_variants_per_locus)
    peak_matrix_multiplier = float(args.finemap_ld_peak_matrix_multiplier)
    reserved_worker_memory_gb = float(args.minimum_memory_per_worker_gb)
    if sample_size_policy != "warn" or sample_size_summary_statistic != "median":
        raise ValueError("Unsupported fine-mapping sample-size policy")
    if not math.isfinite(warning_threshold) or warning_threshold < 0:
        raise ValueError("Sample-size relative-range warning threshold must be non-negative")

    tasks = []
    skipped = []
    for row in loci_df.iter_rows(named=True):
        chrom = str(row["CHR"]).replace("chr", "")
        start, end = int(row["START"]), int(row["END"])
        locus_id = f"chr{chrom}_{start}_{end}"
        locus_coordinates = {
            "chromosome": chrom,
            "start": start,
            "end": end,
        }
        locus_ss = sumstats_df.filter(
            (pl.col("chromosome") == chrom)
            & pl.col("position").is_between(start, end, closed="both")
        )
        if locus_ss.is_empty():
            skipped.append({
                "status": "skipped", "locus_id": locus_id,
                "genomic_locus": str(row["GenomicLocus"]),
                "failure_reason": "no_harmonised_variants",
                **locus_coordinates,
            })
            continue

        resource_qc = dense_ld_resource_qc(
            locus_ss.height,
            maximum_variants,
            peak_matrix_multiplier,
            reserved_worker_memory_gb,
        )
        if resource_qc["resource_failure_reason"]:
            failure_detail = variant_limit_failure_detail(resource_qc)
            logger.error(
                "[STAGE] locus=%s stage=ld_resource_guard status=failed "
                "reason=%s variants=%d limit=%d estimated_peak_memory_gb=%.3f "
                "detail=%s",
                locus_id,
                resource_qc["resource_failure_reason"],
                resource_qc["input_variant_count"],
                resource_qc["maximum_variants_per_locus"],
                resource_qc["estimated_peak_ld_memory_gb"],
                failure_detail,
            )
            skipped.append({
                "status": "resource_limit_exceeded",
                "locus_id": locus_id,
                "genomic_locus": str(row["GenomicLocus"]),
                "failure_reason": resource_qc["resource_failure_reason"],
                "failure_detail": failure_detail,
                **locus_coordinates,
                **resource_qc,
            })
            continue

        nef_values = locus_ss.get_column("NEF").to_numpy().astype(float)
        if not len(nef_values) or not np.isfinite(nef_values).all() or (nef_values <= 0).any():
            skipped.append({
                "status": "skipped", "locus_id": locus_id,
                "genomic_locus": str(row["GenomicLocus"]),
                "failure_reason": "invalid_sample_size",
                **locus_coordinates,
                **resource_qc,
            })
            continue
        median_nef = float(np.median(nef_values))
        nef_min = float(np.min(nef_values))
        nef_max = float(np.max(nef_values))
        relative_range = (nef_max - nef_min) / median_nef
        sample_size_qc = {
            "nef_min": nef_min,
            "nef_max": nef_max,
            "nef_median": median_nef,
            "nef_relative_range": relative_range,
            "selected_nef": median_nef,
            "nef_policy": sample_size_policy,
            "nef_summary_statistic": sample_size_summary_statistic,
            "nef_warning_threshold": warning_threshold,
            "warning_reason": (
                combine_reasons(
                    "nef_relative_range_exceeds_threshold"
                    if relative_range > warning_threshold else None,
                    resource_qc["resource_warning_reason"],
                )
            ),
        }
        n_samples = int(round(float(median_nef)))
        if n_samples < 1:
            skipped.append({
                "status": "skipped", "locus_id": locus_id,
                "genomic_locus": str(row["GenomicLocus"]),
                "failure_reason": "invalid_sample_size",
                **locus_coordinates,
                **sample_size_qc,
                **resource_qc,
            })
            continue
        if sample_size_qc["warning_reason"]:
            logger.warning(
                "[STAGE] locus=%s stage=sample_size_qc status=warning reason=%s "
                "nef_min=%.12g nef_median=%.12g nef_max=%.12g relative_range=%.12g",
                locus_id, sample_size_qc["warning_reason"], nef_min,
                median_nef, nef_max, relative_range,
            )

        z_out = dirs["temp"] / f"{locus_id}.z"
        snp_out = dirs["temp"] / f"{locus_id}.snp"
        write_finemap_z_file(locus_ss, z_out)
        write_snp_file(locus_ss, snp_out)
        tasks.append({
            "locus_id": locus_id,
            "z_file": z_out,
            "snp_file": snp_out,
            "out_dir": dirs["loci"],
            "ld_ref": ld_ref_prefix,
            "n_samples": n_samples,
            "finemap_config": finemap_config,
            "ldstore_timeout_seconds": args.ldstore_timeout_seconds,
            "finemap_timeout_seconds": args.finemap_timeout_seconds,
            "termination_grace_seconds": args.termination_grace_seconds,
            "plink": args.plink,
            "bgenix": args.bgenix,
            "ldstore": args.ldstore,
            "finemap_executable": args.finemap_executable,
            "genomic_locus": str(row["GenomicLocus"]),
            **locus_coordinates,
            "runtime_defaults": runtime_defaults,
            **(
                {
                    "file_log_level": args.fine_mapping_logging[
                        "file_level"
                    ]
                }
                if hasattr(args, "fine_mapping_logging")
                else {}
            ),
            "sample_size_qc": sample_size_qc,
            "resource_qc": resource_qc,
        })

    skip_columns = [
        "status", "locus_id", "genomic_locus", "failure_reason",
        "failure_detail", "chromosome", "start", "end",
        "nef_min", "nef_max", "nef_median", "nef_relative_range",
        "selected_nef", "nef_policy", "nef_summary_statistic",
        "nef_warning_threshold", "warning_reason",
        "input_variant_count", "maximum_variants_per_locus",
        "dense_ld_matrix_gb", "ld_peak_matrix_multiplier",
        "estimated_peak_ld_memory_gb", "reserved_worker_memory_gb",
        "resource_warning_reason", "resource_failure_reason",
    ]
    pd.DataFrame(skipped).reindex(columns=skip_columns).to_csv(
        task_skips_file, sep="\t", index=False
    )
    return tasks, skipped


def _reconcile_extracted_variants(z_file, snp_file, extracted_bim):
    """Reorder/filter Z records to extracted BIM order and verify allele coding."""
    z_df = pd.read_csv(z_file, sep=r"\s+", dtype={"rsid": str, "chromosome": str})
    bim_df = _read_bim(extracted_bim)
    z_lookup = z_df.set_index("rsid", drop=False)
    records = []
    for row in bim_df.itertuples(index=False):
        if row.rsid not in z_lookup.index:
            raise ValueError(f"Extracted reference variant {row.rsid} is absent from the FINEMAP Z file")
        record = z_lookup.loc[row.rsid].copy()
        if isinstance(record, pd.DataFrame):
            raise ValueError(f"Duplicate Z-file identifier after extraction: {row.rsid}")
        if str(record["chromosome"]) != str(row.ref_chromosome) or int(record["position"]) != int(row.ref_position):
            raise ValueError(f"Coordinate mismatch after PLINK extraction for {row.rsid}")

        effect, other = str(record["allele1"]).upper(), str(record["allele2"]).upper()
        ref_effect, ref_other = row.ref_allele1, row.ref_allele2
        direct = (effect, other) == (ref_effect, ref_other)
        swapped = (effect, other) == (ref_other, ref_effect)
        comp = (_complement(effect), _complement(other))
        comp_direct = comp == (ref_effect, ref_other)
        comp_swapped = comp == (ref_other, ref_effect)
        if not (direct or swapped or comp_direct or comp_swapped):
            raise ValueError(f"Allele mismatch after PLINK extraction for {row.rsid}")
        if swapped or comp_swapped:
            record["beta"] = -float(record["beta"])
        record["chromosome"] = row.ref_chromosome
        record["position"] = int(row.ref_position)
        record["allele1"] = ref_effect
        record["allele2"] = ref_other
        records.append(record)

    if not records:
        raise ValueError("PLINK extracted zero variants")
    reconciled = pd.DataFrame(records)[z_df.columns]
    reconciled.to_csv(z_file, sep=" ", index=False, float_format="%.12g")
    reconciled[["rsid"]].to_csv(snp_file, index=False, header=False)
    return len(reconciled)


def _validate_ld_matrix(
    ld_matrix,
    n_variants,
    max_correlation=DEFAULT_LD_MAX_CORRELATION,
    correlation_tolerance=DEFAULT_LD_CORRELATION_TOLERANCE,
    symmetry_tolerance=DEFAULT_LD_SYMMETRY_TOLERANCE,
    diagonal_tolerance=DEFAULT_LD_DIAGONAL_TOLERANCE,
    eigenvalue_tolerance=DEFAULT_LD_EIGENVALUE_TOLERANCE,
):
    matrix = np.loadtxt(ld_matrix, dtype=float, ndmin=2)
    if matrix.shape != (n_variants, n_variants):
        raise ValueError(
            f"LD matrix shape {matrix.shape} does not match {n_variants} harmonised variants"
        )
    if not np.isfinite(matrix).all():
        raise ValueError("LD matrix contains non-finite values")
    maximum_correlation = float(np.max(np.abs(matrix)))
    if maximum_correlation > float(max_correlation) + float(correlation_tolerance):
        raise ValueError(
            f"LD matrix contains an invalid correlation magnitude ({maximum_correlation:.6g})"
        )
    asymmetry = float(np.max(np.abs(matrix - matrix.T)))
    diagonal_error = float(np.max(np.abs(np.diag(matrix) - 1.0)))
    if asymmetry > float(symmetry_tolerance):
        raise ValueError(f"LD matrix is not symmetric (maximum difference {asymmetry:.3g})")
    if diagonal_error > float(diagonal_tolerance):
        raise ValueError(f"LD matrix diagonal is not one (maximum error {diagonal_error:.3g})")
    min_eigenvalue = float(np.linalg.eigvalsh((matrix + matrix.T) / 2.0).min())
    if min_eigenvalue < -float(eigenvalue_tolerance):
        raise ValueError(f"LD matrix is not positive semidefinite (minimum eigenvalue {min_eigenvalue:.3g})")
    return {"ld_asymmetry": asymmetry, "ld_diagonal_error": diagonal_error, "ld_min_eigenvalue": min_eigenvalue}


def _effective_max_causal_snps(configured: int, available_variants: int) -> int:
    """Return a valid per-locus upper bound without changing its interpretation."""
    configured = int(configured)
    available_variants = int(available_variants)
    if configured < 1:
        raise ValueError("Configured FINEMAP n_causal_snps must be at least 1")
    if available_variants < 1:
        raise ValueError("A FINEMAP locus must contain at least one variant")
    return min(configured, available_variants)


def _reconcile_model_selection_failures(results, failures):
    """Replace only model-selection-failed locus outcomes with failure records."""
    failures_by_locus = {failure["locus_id"]: failure for failure in failures}
    reconciled = []
    for result in results:
        row = dict(result)
        failure = failures_by_locus.get(row.get("locus_id"))
        if failure is not None:
            row.update({
                "status": "failed",
                "analysis_status": "failed",
                "outcome_reason": failure["failure_reason"],
                "failure_reason": failure["failure_reason"],
                "failure_detail": failure.get("failure_detail"),
            })
        reconciled.append(row)
    return reconciled


def _process_single_locus(task):
    """Run PLINK, LDstore, and FINEMAP for one locus."""
    locus_id = task["locus_id"]
    runtime_defaults = dict(get_finemap_defaults())
    runtime_defaults.update(task.get("runtime_defaults", {}))
    target_coverage = float(runtime_defaults["credible_set_coverage"])
    locus_dir = Path(task["out_dir"]) / locus_id
    locus_dir.mkdir(parents=True, exist_ok=True)
    logger.info(
        "[STAGE] locus=%s stage=initialization status=completed output=%s",
        locus_id,
        locus_dir,
    )

    z_file = locus_dir / f"{locus_id}.z"
    snp_file = locus_dir / f"{locus_id}.snp"
    plink_out = locus_dir / f"{locus_id}_plink"
    bgen_prefix = locus_dir / f"{locus_id}_genotypes"
    bgen_file = locus_dir / f"{locus_id}_genotypes.bgen"
    bgi_file = locus_dir / f"{locus_id}_genotypes.bgen.bgi"
    ldstore_master = locus_dir / f"{locus_id}.ldstore.master"
    bcor_file = locus_dir / f"{locus_id}.bcor"
    ld_matrix = locus_dir / f"{locus_id}.ld"
    master_file = locus_dir / f"{locus_id}.master"
    config_file = locus_dir / f"{locus_id}.config"
    cred_file = locus_dir / f"{locus_id}.cred"
    log_file = locus_dir / f"{locus_id}.log"

    try:
        logger.info("[STAGE] locus=%s stage=stale_output_cleanup status=started", locus_id)
        for stale in locus_dir.glob(f"{locus_id}.cred*"):
            if stale.is_file():
                stale.unlink()
                logger.info(
                    "[STAGE] locus=%s stage=stale_output_cleanup status=removed file=%s",
                    locus_id,
                    stale,
                )
        logger.info("[STAGE] locus=%s stage=input_copy status=started", locus_id)
        shutil.copy2(task["z_file"], z_file)
        shutil.copy2(task["snp_file"], snp_file)
        logger.info("[STAGE] locus=%s stage=input_copy status=completed", locus_id)

        logger.info("[STAGE] locus=%s stage=plink_extraction status=started", locus_id)
        run_plink_extraction(
            task["ld_ref"],
            snp_file,
            plink_out,
            plink_binary=task["plink"],
            plink_memory_mb=runtime_defaults["plink_memory_mb"],
        )
        logger.info("[STAGE] locus=%s stage=plink_extraction status=completed", locus_id)
        logger.info("[STAGE] locus=%s stage=allele_reconciliation status=started", locus_id)
        n_snps_available = _reconcile_extracted_variants(
            z_file, snp_file, plink_out.with_suffix(".bim")
        )
        logger.info(
            "[STAGE] locus=%s stage=allele_reconciliation status=completed variants=%d",
            locus_id,
            n_snps_available,
        )
        logger.info("[STAGE] locus=%s stage=bgen_conversion status=started", locus_id)
        run_plink_to_bgen(
            plink_out,
            bgen_prefix,
            plink_binary=task["plink"],
            bgen_bits=runtime_defaults["bgen_bits"],
        )
        logger.info("[STAGE] locus=%s stage=bgen_conversion status=completed", locus_id)
        logger.info("[STAGE] locus=%s stage=bgen_indexing status=started", locus_id)
        run_bgen_indexing(bgen_file, bgenix_binary=task["bgenix"])
        logger.info("[STAGE] locus=%s stage=bgen_indexing status=completed", locus_id)

        extracted_fam = plink_out.with_suffix(".fam")
        with extracted_fam.open("rb") as handle:
            n_ld_samples = sum(1 for _ in handle)
        if n_ld_samples < 1:
            raise ValueError("Extracted LD reference contains no samples")

        logger.info(
            "[STAGE] locus=%s stage=ldstore status=started ld_samples=%d",
            locus_id,
            n_ld_samples,
        )
        create_ldstore_master(
            ldstore_master, z_file, bgen_file, bgi_file,
            bcor_file, ld_matrix, n_ld_samples,
        )
        run_ldstore(
            ldstore_master,
            threads=int(runtime_defaults["external_tool_threads"]),
            ldstore_binary=task["ldstore"],
            timeout_seconds=task["ldstore_timeout_seconds"],
            termination_grace_seconds=task["termination_grace_seconds"],
        )
        logger.info("[STAGE] locus=%s stage=ldstore status=completed", locus_id)
        logger.info("[STAGE] locus=%s stage=ld_validation status=started", locus_id)
        ld_qc = _validate_ld_matrix(
            ld_matrix,
            n_snps_available,
            max_correlation=runtime_defaults["ld_max_correlation"],
            correlation_tolerance=runtime_defaults["ld_correlation_tolerance"],
            symmetry_tolerance=runtime_defaults["ld_symmetry_tolerance"],
            diagonal_tolerance=runtime_defaults["ld_diagonal_tolerance"],
            eigenvalue_tolerance=runtime_defaults["ld_eigenvalue_tolerance"],
        )
        logger.info(
            "[STAGE] locus=%s stage=ld_validation status=completed min_eigenvalue=%.12g",
            locus_id,
            ld_qc["ld_min_eigenvalue"],
        )

        local_config = task["finemap_config"].copy()
        local_config["prob_cred_set"] = target_coverage
        configured_max_causal_snps = int(local_config["n_causal_snps"])
        effective_max_causal_snps = _effective_max_causal_snps(
            configured_max_causal_snps, n_snps_available
        )
        local_config["n_causal_snps"] = effective_max_causal_snps
        logger.info(
            "[STAGE] locus=%s stage=finemap_configuration status=completed "
            "configured_max_causal_snps=%d effective_max_causal_snps=%d "
            "available_variants=%d",
            locus_id,
            configured_max_causal_snps,
            effective_max_causal_snps,
            n_snps_available,
        )
        create_finemap_master(
            master_file, z_file, ld_matrix, snp_file,
            config_file, cred_file, log_file, task["n_samples"],
        )
        logger.info(
            "[STAGE] locus=%s stage=finemap status=started target_coverage=%.2f",
            locus_id,
            target_coverage,
        )
        success, msg = run_finemap_binary(
            master_file,
            local_config,
            threads=int(runtime_defaults["external_tool_threads"]),
            finemap_binary=task["finemap_executable"],
            timeout_seconds=task["finemap_timeout_seconds"],
            termination_grace_seconds=task["termination_grace_seconds"],
        )
        cred_files = finemap_credible_files(locus_dir, locus_name=locus_id)
        if not success or not cred_files:
            raise RuntimeError(
                msg if not success else "FINEMAP produced no credible-set file"
            )
        if not config_file.is_file() or config_file.stat().st_size == 0:
            raise RuntimeError("FINEMAP produced no non-empty configuration file")
        logger.info(
            "[STAGE] locus=%s stage=model_probability_validation "
            "status=started credible_files=%d",
            locus_id,
            len(cred_files),
        )
        model_probabilities = {
            path.name: parse_cred_header(path) for path in cred_files
        }
        logger.info(
            "[STAGE] locus=%s stage=model_probability_validation "
            "status=completed probabilities=%s",
            locus_id,
            ",".join(
                f"{name}:{probability:.12g}"
                for name, probability in model_probabilities.items()
            ),
        )
        logger.info(
            "[STAGE] locus=%s stage=finemap status=completed credible_files=%d",
            locus_id, len(cred_files),
        )

        logger.info("[STAGE] locus=%s stage=qc_write status=started", locus_id)
        pd.DataFrame([{
            "locus_id": locus_id,
            "n_gwas_samples": task["n_samples"],
            "n_ld_samples": n_ld_samples,
            "n_variants": n_snps_available,
            "configured_max_causal_snps": configured_max_causal_snps,
            "effective_max_causal_snps": effective_max_causal_snps,
            "target_coverage": target_coverage,
            "ldstore_timeout_seconds": task["ldstore_timeout_seconds"],
            "finemap_timeout_seconds": task["finemap_timeout_seconds"],
            "termination_grace_seconds": task["termination_grace_seconds"],
            "model_probability_validation": "valid",
            "n_model_probabilities": len(model_probabilities),
            **task["sample_size_qc"],
            **task["resource_qc"],
            **ld_qc,
        }]).to_csv(locus_dir / f"{locus_id}_QC.tsv", sep="\t", index=False)
        logger.info("[STAGE] locus=%s stage=qc_write status=completed", locus_id)

        logger.info("[STAGE] locus=%s stage=temporary_cleanup status=started", locus_id)
        for path in [
            plink_out.with_suffix(".bed"), plink_out.with_suffix(".bim"),
            plink_out.with_suffix(".fam"), bgen_file, bgi_file, bcor_file,
        ]:
            if path.exists():
                path.unlink()
                logger.info(
                    "[STAGE] locus=%s stage=temporary_cleanup status=removed file=%s",
                    locus_id,
                    path,
                )
        logger.info("[STAGE] locus=%s stage=completed status=success", locus_id)
        return {
            "status": "success",
            "locus_id": locus_id,
            "cred_file": str(cred_file),
            "config_file": str(config_file),
            "genomic_locus": task["genomic_locus"],
            "analysis_status": "success",
            "outcome_reason": "finemap_completed",
            "failure_reason": None,
            "failure_detail": None,
            "ldstore_timeout_seconds": task["ldstore_timeout_seconds"],
            "finemap_timeout_seconds": task["finemap_timeout_seconds"],
            "termination_grace_seconds": task["termination_grace_seconds"],
            "chromosome": task["chromosome"],
            "start": task["start"],
            "end": task["end"],
            **task["sample_size_qc"],
            **task["resource_qc"],
        }
    except Exception as exc:
        failure_reason = (
            exc.reason
            if isinstance(
                exc,
                (InvalidModelProbabilityError, FinemapExternalToolTimeout),
            )
            else str(exc)
        )
        incomplete_outputs_removed = []
        if isinstance(exc, FinemapExternalToolTimeout):
            timeout_outputs = (
                [bcor_file, ld_matrix]
                if exc.reason == "ldstore_timeout"
                else [
                    config_file,
                    log_file,
                    *locus_dir.glob(f"{locus_id}.cred*"),
                ]
            )
            for partial_output in timeout_outputs:
                if partial_output.is_file():
                    partial_output.unlink()
                    incomplete_outputs_removed.append(str(partial_output))
                    logger.warning(
                        "[STAGE] locus=%s stage=incomplete_output_cleanup "
                        "status=removed file=%s",
                        locus_id,
                        partial_output,
                    )
        logger.exception(
            "[STAGE] locus=%s stage=failed status=failed reason=%s detail=%s",
            locus_id,
            failure_reason,
            exc,
        )
        failure_record = {
            "status": "failed", "locus_id": locus_id,
            "genomic_locus": task["genomic_locus"],
            "analysis_status": "failed",
            "outcome_reason": failure_reason,
            "failure_reason": failure_reason,
            "failure_detail": str(exc),
            "timeout_stage": getattr(exc, "stage", None),
            "timeout_seconds": getattr(exc, "timeout_seconds", None),
            "timeout_elapsed_seconds": getattr(exc, "elapsed_seconds", None),
            "timeout_forced_termination": getattr(
                exc, "forced_termination", None
            ),
            "incomplete_outputs_removed": ";".join(
                incomplete_outputs_removed
            ) or None,
            "ldstore_timeout_seconds": task["ldstore_timeout_seconds"],
            "finemap_timeout_seconds": task["finemap_timeout_seconds"],
            "termination_grace_seconds": task["termination_grace_seconds"],
            "chromosome": task["chromosome"],
            "start": task["start"],
            "end": task["end"],
            **task["sample_size_qc"],
            **task["resource_qc"],
        }
        pd.DataFrame([failure_record]).to_csv(
            locus_dir / f"{locus_id}_QC.tsv", sep="\t", index=False
        )
        return failure_record


def process_single_locus(task):
    """Run one locus with a detailed file log and no worker console handler."""
    if "file_log_level" not in task:
        raise ValueError(
            "FINEMAP worker task is missing the resolved YAML file log level"
        )
    locus_id = task["locus_id"]
    locus_dir = Path(task["out_dir"]) / locus_id
    with detailed_file_logging(
        "postgwas.modules.fine_mapping",
        locus_dir / f"{locus_id}_debug.log",
        task["file_log_level"],
    ):
        return _process_single_locus(task)


def _finemap_run_configuration(
    args, dirs, ld_ref_prefix, genome_build, plink_version,
    runtime_defaults, resource_parameters,
):
    """Build the canonical FINEMAP run record for success or early failure."""
    return {
        "engine": "FINEMAP",
        "inputs": {
            "locus_file": args.locus_file,
            "summary_statistics": args.finemap_in_files,
            "ld_reference": ld_ref_prefix,
            "output_directory": dirs["root"],
        },
        "analysis_parameters": {
            "genome_build": genome_build,
            "locus_type": args.locus_type,
            "window_kb": args.window_kb,
            "lp_threshold": args.lp_threshold,
            "credible_set_coverage": runtime_defaults[
                "credible_set_coverage"
            ],
            "finemap_config": _resolved_finemap_config(
                args, runtime_defaults["credible_set_coverage"]
            ),
            "external_process_timeouts": {
                "ldstore_timeout_seconds": args.ldstore_timeout_seconds,
                "finemap_timeout_seconds": args.finemap_timeout_seconds,
                "termination_grace_seconds": (
                    args.termination_grace_seconds
                ),
            },
            "sample_size_policy": args.sample_size_policy,
            "sample_size_summary_statistic": (
                args.sample_size_summary_statistic
            ),
            "sample_size_relative_range_warning_threshold": (
                args.sample_size_relative_range_warning_threshold
            ),
            "maximum_variants_per_locus": args.maximum_variants_per_locus,
            "finemap_ld_peak_matrix_multiplier": (
                args.finemap_ld_peak_matrix_multiplier
            ),
        },
        "resource_parameters": resource_parameters,
        "output_layout": args.fine_mapping_output_layout,
        "software": {
            "plink_path": args.plink,
            "plink_version": plink_version,
            "bgenix_path": args.bgenix,
            "ldstore_path": args.ldstore,
            "finemap_path": args.finemap_executable,
        },
        "defaults": runtime_defaults,
    }


def _run_finemap_pipeline(args, progress, screen=None):
    """Implement the validated FINEMAP workflow with structured progress."""
    resolve_compute_args(args)
    outdir = Path(args.output_directory).resolve()
    start_time = time()
    pipeline_total = DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL
    runtime_defaults = dict(args.fine_mapping_runtime_defaults)
    runtime_defaults["credible_set_coverage"] = float(args.prob_cred_set)
    progress.record(
        "pipeline", "initialization", 1, pipeline_total, "completed",
        f"output_dir={outdir}",
    )
    preflight = getattr(args, "_fine_mapping_preflight", None)
    if screen and preflight is None:
        screen.complete(1, [
            ("analysis", "Engine", "FINEMAP"),
            ("info", "Dataset", args.dataset_id),
        ])

    genome_build = getattr(args, "genome_build", DEFAULT_GENOME_BUILD)
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError("genome_build must be GRCh37 or GRCh38")

    logger.info("[STAGE] stage=dependency_validation status=started")
    if preflight is not None:
        runtime_tools = preflight.tools
        logger.info(
            "[STAGE] stage=dependency_validation status=reused_preflight"
        )
    else:
        runtime_tools, _ = _validate_runtime_tools(args)
    tool_versions = {tool.name: tool.version for tool in runtime_tools}
    plink_version = tool_versions["PLINK"]
    logger.info(
        "[STAGE] stage=dependency_validation status=completed plink=%s version=%s",
        args.plink,
        plink_version,
    )
    progress.record(
        "pipeline", "dependency_validation", 2, pipeline_total, "completed",
        f"PLINK={plink_version}",
    )
    if screen and preflight is None:
        screen.complete(2, [
            ("success", "External tools", "validated"),
            ("info", "PLINK", plink_version),
            ("analysis", "Genome build", genome_build),
        ])
    dirs = setup_directories(
        args.output_directory,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    logger.info("[STAGE] stage=input_preparation status=started")
    loci_df, sumstats_filt, ld_ref_prefix = load_and_prep_inputs(
        args,
        dirs["finemap_input_qc_file"],
        schema_inference_length=runtime_defaults["schema_inference_length"],
        max_maf=runtime_defaults["max_maf"],
    )
    progress.record(
        "pipeline", "input_preparation", 3, pipeline_total, "completed",
        f"loci={loci_df.height}; harmonised_variants={sumstats_filt.height}",
    )
    if screen:
        screen.complete(3, [
            ("count", "Input loci", loci_df.height),
            ("count", "Harmonised variants", sumstats_filt.height),
            ("analysis", "Boundary flank", f"{args.window_kb:,} kb"),
        ])
    logger.info("[STAGE] stage=task_generation status=started")
    tasks, generation_skips = generate_tasks(
        loci_df,
        sumstats_filt,
        ld_ref_prefix,
        dirs,
        args,
        dirs["finemap_task_skips_file"],
        target_coverage=runtime_defaults["credible_set_coverage"],
        runtime_defaults=runtime_defaults,
    )
    if not tasks:
        generation_warning_count = sum(
            bool(row.get("warning_reason")) for row in generation_skips
        )
        pd.DataFrame(generation_skips).sort_values(
            "locus_id", kind="stable"
        ).to_csv(
            dirs["finemap_locus_status_file"], sep="\t", index=False
        )
        write_run_configuration(
            dirs["run_configuration_file"],
            _finemap_run_configuration(
                args, dirs, ld_ref_prefix, genome_build, plink_version,
                runtime_defaults,
                {
                    "requested_threads": args.threads,
                    "selected_workers": 0,
                    "minimum_memory_per_worker_gb": (
                        args.minimum_memory_per_worker_gb
                    ),
                },
            ),
        )
        logger.error(
            "FINEMAP final summary: successful=0 attempted=%d failed_or_skipped=%d "
            "sample_size_warnings=%d outcomes=%s",
            len(generation_skips), len(generation_skips),
            generation_warning_count,
            dirs["finemap_locus_status_file"],
        )
        for reason, count in pd.Series(
            [row["failure_reason"] for row in generation_skips]
        ).value_counts().sort_index().items():
            logger.error("  failure reason [%s]: %d", reason, int(count))
        n_variant_limit_failures = sum(
            row.get("failure_reason") == VARIANT_LIMIT_FAILURE_REASON
            for row in generation_skips
        )
        if n_variant_limit_failures:
            logger.error(
                "  Variant-limit failures: %d; override with "
                "--maximum-variants-per-locus COUNT or edit "
                "modules.fine_mapping.ld_resource_guard."
                "maximum_variants_per_locus in YAML",
                n_variant_limit_failures,
            )
        raise RuntimeError(
            "No FINEMAP loci were eligible for execution; inspect %s for "
            "variant limits, input eligibility, and sample-size failures."
            % dirs["finemap_locus_status_file"]
        )

    requested_threads = args.threads
    ram_per_worker_gb = max(
        float(runtime_defaults["finemap_ram_per_worker_gb"]),
        float(
            getattr(
                args,
                "minimum_memory_per_worker_gb",
                runtime_defaults["finemap_ram_per_worker_gb"],
            ) or runtime_defaults["finemap_ram_per_worker_gb"]
        ),
        max(
            float(task["resource_qc"]["estimated_peak_ld_memory_gb"])
            for task in tasks
        ),
    )
    if ram_per_worker_gb <= 0:
        raise ValueError("minimum_memory_per_worker_gb must be greater than zero")
    max_workers = safe_thread_count(
        min(int(requested_threads), len(tasks)),
        ram_per_worker_gb,
        available_ram_gb=args.memory_gb,
        enforce_memory_budget=True,
        reporter=None,
    )
    logger.info(
        "[STAGE] stage=task_generation status=completed tasks=%d workers=%d "
        "memory_budget_gb=%.3f ram_per_worker_gb=%.3f",
        len(tasks),
        max_workers,
        args.memory_gb,
        ram_per_worker_gb,
    )
    write_run_configuration(
        dirs["run_configuration_file"],
        _finemap_run_configuration(
            args, dirs, ld_ref_prefix, genome_build, plink_version,
            runtime_defaults,
            {
                "requested_threads": requested_threads,
                "selected_workers": max_workers,
                "memory_budget_gb": args.memory_gb,
                "ram_per_worker_gb": ram_per_worker_gb,
            },
        ),
    )
    progress.record(
        "pipeline", "task_generation", 4, pipeline_total, "completed",
        f"tasks={len(tasks)}; workers={max_workers}",
    )
    if screen:
        screen.complete(4, [
            ("count", "Eligible loci", len(tasks)),
            (
                "warning" if generation_skips else "info",
                "Skipped loci",
                len(generation_skips),
            ),
            ("analysis", "Parallel workers", max_workers),
        ])

    results = []
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=mp.get_context("spawn")) as executor:
        futures = {executor.submit(process_single_locus, task): task for task in tasks}
        for completed_count, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                logger.exception(
                    "[STAGE] locus=%s stage=worker_collection status=failed",
                    task["locus_id"],
                )
                result = {
                    "status": "failed",
                    "locus_id": task["locus_id"],
                    "genomic_locus": task["genomic_locus"],
                    "failure_reason": str(exc),
                    "chromosome": task["chromosome"],
                    "start": task["start"],
                    "end": task["end"],
                    **task["sample_size_qc"],
                    **task["resource_qc"],
                }
            results.append(result)
            progress.record(
                "loci",
                "finemap_execution",
                completed_count,
                len(tasks),
                result.get("status", "unknown"),
                f"locus={result.get('locus_id', task['locus_id'])}",
            )

    all_locus_outcomes = [*results, *generation_skips]
    results_df = pd.DataFrame(all_locus_outcomes).sort_values("locus_id", kind="stable")
    results_df.to_csv(dirs["finemap_locus_status_file"], sep="\t", index=False)
    successful_results = [result for result in results if result.get("status") == "success"]
    invalid_model_results = [
        result for result in results
        if result.get("failure_reason") == INVALID_MODEL_PROBABILITY_REASON
    ]
    if not successful_results and not invalid_model_results:
        logger.error(
            "FINEMAP final summary: successful=0 attempted=%d failed_or_skipped=%d "
            "outcomes=%s",
            len(all_locus_outcomes), len(all_locus_outcomes),
            dirs["finemap_locus_status_file"],
        )
        for reason, count in results_df["failure_reason"].dropna().astype(
            str
        ).value_counts().sort_index().items():
            logger.error("  failure reason [%s]: %d", reason, int(count))
        n_variant_limit_failures = int(
            results_df["failure_reason"].eq(VARIANT_LIMIT_FAILURE_REASON).sum()
        )
        if n_variant_limit_failures:
            logger.error(
                "  Variant-limit failures: %d; override with "
                "--maximum-variants-per-locus COUNT or edit "
                "modules.fine_mapping.ld_resource_guard."
                "maximum_variants_per_locus in YAML",
                n_variant_limit_failures,
            )
        raise RuntimeError(
            "All FINEMAP loci failed; inspect %s and per-locus logs in %s"
            % (
                dirs["finemap_locus_status_file"],
                dirs["finemap_locus_work_directory"],
            )
        )
    progress.record(
        "pipeline", "locus_execution", 5, pipeline_total, "completed",
        f"successful={len(successful_results)}; "
        f"failed_or_skipped={len(all_locus_outcomes) - len(successful_results)}",
    )
    if screen:
        screen.complete(5, [
            ("success", "Successful fits", len(successful_results)),
            (
                "warning"
                if len(all_locus_outcomes) > len(successful_results)
                else "info",
                "Failed or skipped",
                len(all_locus_outcomes) - len(successful_results),
            ),
        ])

    allowed_loci = {result["locus_id"] for result in successful_results}
    model_selection_loci = allowed_loci | {
        result["locus_id"] for result in invalid_model_results
    }
    locus_metadata = {
        task["locus_id"]: {
            "GenomicLocus": task["genomic_locus"],
            "chr": task["chromosome"],
            "start": task["start"],
            "end": task["end"],
            **task["sample_size_qc"],
            **task["resource_qc"],
            "analysis_status": "success",
            "outcome_reason": "finemap_completed",
            "failure_reason": None,
        }
        for task in tasks if task["locus_id"] in allowed_loci
    }
    output_result = process_finemap_output(
        raw_dir=dirs["loci"],
        inter_dir=str(dirs["inter"]),
        final_dir=str(dirs["flames"]),
        model_summary_file=dirs["finemap_model_summary_file"],
        annotation_directory=dirs[
            "primary_flames_annotations_directory"
        ],
        index_filename=args.overlap_resolution["index_filename"],
        manifest_filename=args.flames_manifest_filename,
        annotation_prefix=args.overlap_resolution["annotation_prefix"],
        allowed_loci=model_selection_loci,
        locus_metadata=locus_metadata,
        genome_build=genome_build,
        target_coverage=runtime_defaults["credible_set_coverage"],
        coverage_tolerance=runtime_defaults["coverage_tolerance"],
    )
    results = _reconcile_model_selection_failures(
        results, output_result["failures"]
    )
    all_locus_outcomes = [*results, *generation_skips]
    results_df = pd.DataFrame(all_locus_outcomes).sort_values(
        "locus_id", kind="stable"
    )
    results_df.to_csv(
        dirs["finemap_locus_status_file"], sep="\t", index=False
    )
    successful_results = [
        result for result in results
        if result.get("status") == "success"
        and result["locus_id"] in output_result["selected_loci"]
    ]
    allowed_loci = {result["locus_id"] for result in successful_results}
    locus_metadata = {
        locus: metadata
        for locus, metadata in locus_metadata.items()
        if locus in allowed_loci
    }
    if not successful_results:
        logger.error(
            "FINEMAP final summary: successful=0 attempted=%d "
            "failed_or_skipped=%d outcomes=%s",
            len(all_locus_outcomes),
            len(all_locus_outcomes),
            dirs["finemap_locus_status_file"],
        )
        for reason, count in results_df["failure_reason"].dropna().astype(
            str
        ).value_counts().sort_index().items():
            logger.error("  failure reason [%s]: %d", reason, int(count))
        raise RuntimeError(
            "All FINEMAP loci failed model-output validation; inspect %s and %s"
            % (
                dirs["finemap_locus_status_file"],
                dirs["finemap_model_summary_file"],
            )
        )
    genomic_loci_file = (
        dirs["flames"]
        / args.overlap_resolution["genomic_loci_filename"]
    )
    pd.DataFrame(list(locus_metadata.values())).to_csv(
        genomic_loci_file, sep="\t", index=False
    )
    progress.record(
        "pipeline", "output_formatting", 6, pipeline_total, "completed",
        f"flames_input={dirs['flames']}",
    )
    primary_index_file = (
        dirs["flames"] / args.overlap_resolution["index_filename"]
    )
    n_credible_sets = len(pd.read_csv(primary_index_file, sep="\t"))
    if screen:
        screen.complete(6, [
            ("success", "Validated loci", len(successful_results)),
            ("success", "Primary credible sets", n_credible_sets),
            ("success", "FLAMES handoff", "ready"),
        ])
    pd.DataFrame([
        {"software": "PLINK", "version": plink_version, "path": args.plink},
        {
            "software": "bgenix",
            "version": tool_versions["BGENIX"],
            "path": args.bgenix,
        },
        {
            "software": "LDstore",
            "version": tool_versions["LDstore"],
            "path": args.ldstore,
        },
        {
            "software": "FINEMAP",
            "version": tool_versions["FINEMAP"],
            "path": args.finemap_executable,
        },
        {"software": "genome_build", "version": genome_build, "path": ""},
        {
            "software": "credible_set_coverage",
            "version": str(runtime_defaults["credible_set_coverage"]),
            "path": "",
        },
    ]).to_csv(
        dirs["finemap_software_versions_file"], sep="\t", index=False
    )

    n_successful = len(successful_results)
    n_failed = len(all_locus_outcomes) - n_successful
    n_warned = int(results_df["warning_reason"].notna().sum()) if "warning_reason" in results_df else 0
    warning_reason_counts = {}
    if "warning_reason" in results_df:
        warnings = results_df["warning_reason"].dropna().astype(str)
        warnings = warnings[warnings.ne("")]
        warning_reason_counts = {
            str(reason): int(count)
            for reason, count in warnings.value_counts().sort_index().items()
        }
    failure_reasons = results_df["failure_reason"].dropna().astype(str)
    failure_reasons = failure_reasons[failure_reasons.ne("")]
    failure_reason_counts = {
        str(reason): int(count)
        for reason, count in failure_reasons.value_counts().sort_index().items()
    }
    n_variant_limit_failures = int(
        results_df["failure_reason"].eq(VARIANT_LIMIT_FAILURE_REASON).sum()
    )
    progress.record(
        "pipeline", "complete", pipeline_total, pipeline_total, "success",
        f"successful={n_successful}; failed={n_failed}; warnings={n_warned}; "
        f"elapsed_seconds={time() - start_time:.3f}",
    )
    logger.info(
        "FINEMAP final summary: elapsed=%.2fs successful=%d attempted=%d "
        "failed_or_skipped=%d sample_size_warnings=%d outcomes=%s",
        time() - start_time,
        n_successful,
        len(all_locus_outcomes),
        n_failed,
        n_warned,
        dirs["finemap_locus_status_file"],
    )
    if n_variant_limit_failures:
        logger.info(
            "  Variant-limit failures: %d; override with "
            "--maximum-variants-per-locus COUNT or edit "
            "modules.fine_mapping.ld_resource_guard.maximum_variants_per_locus "
            "in YAML",
            n_variant_limit_failures,
        )
    for column, label in (
        ("warning_reason", "warning reason"),
        ("failure_reason", "failure reason"),
    ):
        if column not in results_df:
            continue
        reasons = results_df[column].dropna().astype(str)
        reasons = reasons[reasons.ne("")]
        for reason, count in reasons.value_counts().sort_index().items():
            logger.info("  %s [%s]: %d", label, reason, int(count))
    if screen:
        screen.complete(7, [
            ("success", "Primary loci successful", n_successful),
            ("warning" if n_failed else "info", "Failed or skipped", n_failed),
            ("warning" if n_warned else "info", "Loci with warnings", n_warned),
            ("success", "Primary credible sets", n_credible_sets),
        ])
    return {
        "status": "success",
        "output_dir": str(dirs["root"]),
        "flames_input": str(dirs["flames"]),
        "credible_set_manifests": [
            str(dirs["flames"] / args.flames_manifest_filename)
        ],
        "flames_index": str(
            dirs["flames"] / args.overlap_resolution["index_filename"]
        ),
        "locus_status": str(dirs["finemap_locus_status_file"]),
        "genomic_loci": str(genomic_loci_file),
        "n_attempted": len(all_locus_outcomes),
        "n_successful": n_successful,
        "n_failed": n_failed,
        "n_warnings": n_warned,
        "warning_reason_counts": warning_reason_counts,
        "failure_reason_counts": failure_reason_counts,
        "n_credible_sets": n_credible_sets,
        "n_variant_limit_failures": n_variant_limit_failures,
    }


def run_finemap_pipeline(args, screen=None):
    """Coordinate FINEMAP and record success, failure, percentage, and remaining work."""
    outdir = Path(args.output_directory).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    output_paths = resolve_output_paths(
        outdir,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    output_paths["run_metadata_directory"].mkdir(parents=True, exist_ok=True)
    with detailed_file_logging(
        "postgwas.modules.fine_mapping",
        output_paths["pipeline_log_file"],
        args.fine_mapping_logging["file_level"],
        mode="a",
    ):
        progress = ProgressRecorder(
            output_paths["pipeline_progress_file"], logger=logger
        )
        try:
            result = _run_finemap_pipeline(args, progress, screen=screen)
            return {
                "status": result["status"],
                "output_dir": result["output_dir"],
                "flames_input": result["flames_input"],
                "credible_set_manifests": result[
                    "credible_set_manifests"
                ],
                "flames_index": result["flames_index"],
                "locus_status": result["locus_status"],
                "genomic_loci": result["genomic_loci"],
                "n_attempted": result["n_attempted"],
                "n_successful": result["n_successful"],
                "n_failed": result["n_failed"],
                "n_warnings": result["n_warnings"],
                "warning_reason_counts": result["warning_reason_counts"],
                "failure_reason_counts": result["failure_reason_counts"],
                "n_credible_sets": result["n_credible_sets"],
                "n_variant_limit_failures": result[
                    "n_variant_limit_failures"
                ],
            }
        except Exception as exc:
            latest = progress.latest.get("pipeline", {})
            completed = int(latest.get("completed", 0))
            progress.record(
                "pipeline", "failed", completed,
                DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL, "failed", str(exc)
            )
            logger.exception("FINEMAP pipeline failed: %s", exc)
            raise
