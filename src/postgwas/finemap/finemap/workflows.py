"""FINEMAP workflow with strict allele, LD, and FLAMES validation."""

import logging
import math
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import time

import numpy as np
import pandas as pd
import polars as pl

from postgwas.finemap.defaults import (
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
    DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS,
    DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS,
    SUPPORTED_GENOME_BUILDS,
    get_finemap_defaults,
)
from postgwas.finemap.finemap.input_gen import (
    create_finemap_master,
    create_ldstore_master,
    write_finemap_z_file,
    write_snp_file,
)
from postgwas.finemap.finemap.merge_results import process_finemap_output
from postgwas.finemap.finemap.runner import (
    run_bgen_indexing,
    run_finemap_binary,
    run_ldstore,
    run_plink_extraction,
    run_plink_to_bgen,
)
from postgwas.finemap.progress import ProgressRecorder, write_run_configuration


logger = logging.getLogger("postgwas.finemap")
TARGET_CREDIBLE_SET_COVERAGE = DEFAULT_CREDIBLE_SET_COVERAGE


def setup_logging(log_file=None, verbose=True):
    """Configure process-local console logging and an optional file log."""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w"))
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="\t\t\t\t%H:%M:%S",
        handlers=handlers,
        force=True,
    )


def _resolve_plink2(args, version_timeout_seconds=DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS):
    """Resolve PLINK 2 and document the exact version used."""
    supplied = getattr(args, "plink", None)
    candidate = str(Path(supplied).resolve()) if supplied and Path(supplied).is_file() else None
    candidate = candidate or shutil.which("plink2")
    if not candidate:
        raise EnvironmentError(
            "PLINK 2 was not found. Provide a PLINK 2 executable with --plink "
            "or place plink2 on PATH."
        )
    probe = subprocess.run(
        [candidate, "--version"], capture_output=True, text=True,
        timeout=float(version_timeout_seconds), check=False,
    )
    version_text = (probe.stdout or probe.stderr).strip()
    if probe.returncode != 0 or "PLINK v2" not in version_text:
        raise EnvironmentError(
            f"FINEMAP BGEN generation requires PLINK 2; '{candidate}' did not report PLINK v2."
        )
    return candidate, version_text.splitlines()[0]


def check_dependencies(plink_binary=None):
    """Ensure required external binaries are available."""
    required_tools = ["bgenix", "ldstore", "finemap"]
    missing = [tool for tool in required_tools if shutil.which(tool) is None]
    if plink_binary is None and shutil.which("plink2") is None:
        missing.append("plink2")
    elif plink_binary is not None and not Path(plink_binary).is_file():
        missing.append(str(plink_binary))
    if missing:
        raise EnvironmentError(
            f"Missing required tools: {', '.join(missing)}. Install them and ensure they are on PATH."
        )


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
            pd.Series([args.finemap_mhc_chrom])
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


def setup_directories(outdir_path):
    """Create and return all pipeline-owned directories."""
    outdir = Path(outdir_path).resolve()
    dirs = {
        "root": outdir,
        "temp": outdir / "temp_inputs",
        "loci": outdir / "loci_data",
        "final": outdir / "final_results",
        "flames": outdir / "flames_input",
        "inter": outdir / "finemap_best_results",
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


def load_and_prep_inputs(
    args,
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

    raw_ref = str(args.finemap_ld_ref)
    ld_ref_prefix = str(Path(raw_ref).with_suffix("")) if raw_ref.endswith((".bed", ".bim", ".fam")) else raw_ref
    for suffix in [".bed", ".bim", ".fam"]:
        if not Path(f"{ld_ref_prefix}{suffix}").is_file():
            raise FileNotFoundError(f"LD reference file is missing: {ld_ref_prefix}{suffix}")

    ref_df = _read_bim(f"{ld_ref_prefix}.bim")
    sumstats_filt, harmonisation_qc = _harmonise_sumstats(
        sumstats, ref_df, max_maf=max_maf
    )
    qc_path = Path(args.outdir).resolve() / "input_harmonization_qc.tsv"
    pd.DataFrame([harmonisation_qc]).to_csv(qc_path, sep="\t", index=False)
    logger.info(
        "Input ready: %d loci and %d harmonised variants",
        loci_df.height,
        sumstats_filt.height,
    )
    return loci_df, sumstats_filt, ld_ref_prefix


def generate_tasks(
    loci_df,
    sumstats_df,
    ld_ref_prefix,
    dirs,
    args,
    target_coverage=DEFAULT_CREDIBLE_SET_COVERAGE,
    runtime_defaults=None,
):
    """Create one validated FINEMAP task per locus."""
    finemap_config = {
        key: getattr(args, key)
        for key in [
            "n_causal_snps", "n_iter", "n_conv_sss",
            "prob_conv_sss_tol", "prior_k", "std_effects",
        ]
    }
    # Scientific and user requirement: every credible set is 95%.
    if float(target_coverage) != DEFAULT_CREDIBLE_SET_COVERAGE:
        raise ValueError("FINEMAP credible-set coverage must be exactly 0.95")
    finemap_config["prob_cred_set"] = float(target_coverage)
    runtime_defaults = dict(runtime_defaults or get_finemap_defaults())

    tasks = []
    skipped = []
    for row in loci_df.iter_rows(named=True):
        chrom = str(row["CHR"]).replace("chr", "")
        start, end = int(row["START"]), int(row["END"])
        locus_id = f"chr{chrom}_{start}_{end}"
        locus_ss = sumstats_df.filter(
            (pl.col("chromosome") == chrom)
            & pl.col("position").is_between(start, end, closed="both")
        )
        if locus_ss.is_empty():
            skipped.append({"locus_id": locus_id, "reason": "no_harmonised_variants"})
            continue

        median_nef = locus_ss.select(pl.col("NEF").median()).item()
        if median_nef is None or not math.isfinite(float(median_nef)) or float(median_nef) <= 0:
            skipped.append({"locus_id": locus_id, "reason": "invalid_sample_size"})
            continue
        n_samples = int(round(float(median_nef)))
        if n_samples < 1:
            skipped.append({"locus_id": locus_id, "reason": "invalid_sample_size"})
            continue

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
            "plink": args.plink,
            "genomic_locus": str(row["GenomicLocus"]),
            "chromosome": chrom,
            "start": start,
            "end": end,
            "runtime_defaults": runtime_defaults,
        })

    pd.DataFrame(skipped, columns=["locus_id", "reason"]).to_csv(
        dirs["root"] / "task_generation_skips.tsv", sep="\t", index=False
    )
    return tasks


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


def process_single_locus(task):
    """Run PLINK, LDstore, and FINEMAP for one locus."""
    locus_id = task["locus_id"]
    runtime_defaults = dict(get_finemap_defaults())
    runtime_defaults.update(task.get("runtime_defaults", {}))
    target_coverage = float(runtime_defaults["credible_set_coverage"])
    if target_coverage != DEFAULT_CREDIBLE_SET_COVERAGE:
        raise ValueError("FINEMAP credible-set coverage must be exactly 0.95")
    locus_dir = Path(task["out_dir"]) / locus_id
    locus_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(log_file=locus_dir / f"{locus_id}_debug.log")
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
        run_bgen_indexing(bgen_file)
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
            n_threads=int(runtime_defaults["external_tool_threads"]),
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
        local_config["n_causal_snps"] = max(
            1, min(int(local_config["n_causal_snps"]), n_snps_available)
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
            n_threads=int(runtime_defaults["external_tool_threads"]),
        )
        cred_files = sorted(
            path for path in locus_dir.glob(f"{locus_id}.cred*")
            if path.is_file() and path.stat().st_size > 0
        )
        if not success or not cred_files:
            raise RuntimeError(msg if not success else "FINEMAP produced no non-empty credible-set file")
        if not config_file.is_file() or config_file.stat().st_size == 0:
            raise RuntimeError("FINEMAP produced no non-empty configuration file")
        logger.info(
            "[STAGE] locus=%s stage=finemap status=completed credible_files=%d",
            locus_id,
            len(cred_files),
        )

        logger.info("[STAGE] locus=%s stage=qc_write status=started", locus_id)
        pd.DataFrame([{
            "locus_id": locus_id,
            "n_gwas_samples": task["n_samples"],
            "n_ld_samples": n_ld_samples,
            "n_variants": n_snps_available,
            "target_coverage": target_coverage,
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
        }
    except Exception as exc:
        logger.exception(
            "[STAGE] locus=%s stage=failed status=failed reason=%s",
            locus_id,
            exc,
        )
        return {"status": "failed", "locus_id": locus_id, "reason": str(exc)}


def _software_version(
    binary,
    args=("--version",),
    timeout_seconds=DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS,
):
    try:
        result = subprocess.run(
            [binary, *args], capture_output=True, text=True,
            timeout=float(timeout_seconds), check=False,
        )
        text = (result.stdout or result.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except Exception:
        return "unknown"


def _run_finemap_pipeline(args, progress):
    """Implement the validated FINEMAP workflow with structured progress."""
    outdir = Path(args.outdir).resolve()
    start_time = time()
    pipeline_total = DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL
    runtime_defaults = get_finemap_defaults()
    progress.record(
        "pipeline", "initialization", 1, pipeline_total, "completed",
        f"output_dir={outdir}",
    )

    genome_build = getattr(args, "genome_build", DEFAULT_GENOME_BUILD)
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError("genome_build must be GRCh37 or GRCh38")

    logger.info("[STAGE] stage=dependency_validation status=started")
    args.plink, plink_version = _resolve_plink2(
        args,
        version_timeout_seconds=runtime_defaults["tool_version_timeout_seconds"],
    )
    check_dependencies(args.plink)
    logger.info(
        "[STAGE] stage=dependency_validation status=completed plink=%s version=%s",
        args.plink,
        plink_version,
    )
    progress.record(
        "pipeline", "dependency_validation", 2, pipeline_total, "completed",
        f"PLINK={plink_version}",
    )
    dirs = setup_directories(args.outdir)
    logger.info("[STAGE] stage=input_preparation status=started")
    loci_df, sumstats_filt, ld_ref_prefix = load_and_prep_inputs(
        args,
        schema_inference_length=runtime_defaults["schema_inference_length"],
        max_maf=runtime_defaults["max_maf"],
    )
    progress.record(
        "pipeline", "input_preparation", 3, pipeline_total, "completed",
        f"loci={loci_df.height}; harmonised_variants={sumstats_filt.height}",
    )
    logger.info("[STAGE] stage=task_generation status=started")
    tasks = generate_tasks(
        loci_df,
        sumstats_filt,
        ld_ref_prefix,
        dirs,
        args,
        target_coverage=runtime_defaults["credible_set_coverage"],
        runtime_defaults=runtime_defaults,
    )
    if not tasks:
        raise RuntimeError("No loci contained usable harmonised variants and valid sample sizes")

    try:
        import psutil
        mem_gb = psutil.virtual_memory().total / (1024.0 ** 3)
    except ImportError:
        try:
            mem_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024.0 ** 3)
        except (AttributeError, ValueError):
            mem_gb = float(runtime_defaults["fallback_memory_gb"])
    requested_threads = getattr(
        args, "nthreads", getattr(args, "threads", os.cpu_count() or 1)
    ) or 1
    ram_per_worker_gb = max(
        float(runtime_defaults["finemap_ram_per_worker_gb"]),
        float(
            getattr(
                args,
                "min_ram_per_worker_gb",
                runtime_defaults["finemap_ram_per_worker_gb"],
            ) or runtime_defaults["finemap_ram_per_worker_gb"]
        ),
    )
    if ram_per_worker_gb <= 0:
        raise ValueError("min_ram_per_worker_gb must be greater than zero")
    max_workers = max(
        1,
        min(
            int(requested_threads),
            len(tasks),
            max(1, int(mem_gb // ram_per_worker_gb)),
        ),
    )
    logger.info(
        "[STAGE] stage=task_generation status=completed tasks=%d workers=%d "
        "detected_memory_gb=%.3f ram_per_worker_gb=%.3f",
        len(tasks),
        max_workers,
        mem_gb,
        ram_per_worker_gb,
    )
    write_run_configuration(
        dirs["root"] / "run_configuration.json",
        {
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
                "credible_set_coverage": runtime_defaults["credible_set_coverage"],
                "finemap_config": tasks[0]["finemap_config"],
            },
            "resource_parameters": {
                "requested_threads": requested_threads,
                "selected_workers": max_workers,
                "detected_memory_gb": mem_gb,
                "ram_per_worker_gb": ram_per_worker_gb,
            },
            "software": {"plink_path": args.plink, "plink_version": plink_version},
            "defaults": runtime_defaults,
        },
    )
    progress.record(
        "pipeline", "task_generation", 4, pipeline_total, "completed",
        f"tasks={len(tasks)}; workers={max_workers}",
    )

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
                    "reason": str(exc),
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

    results_df = pd.DataFrame(results).sort_values("locus_id", kind="stable")
    results_df.to_csv(dirs["root"] / "finemap_locus_status.tsv", sep="\t", index=False)
    successful_results = [result for result in results if result.get("status") == "success"]
    if not successful_results:
        raise RuntimeError("All FINEMAP loci failed; inspect finemap_locus_status.tsv and per-locus logs")
    progress.record(
        "pipeline", "locus_execution", 5, pipeline_total, "completed",
        f"successful={len(successful_results)}; failed={len(results) - len(successful_results)}",
    )

    allowed_loci = {result["locus_id"] for result in successful_results}
    locus_metadata = {
        task["locus_id"]: {
            "GenomicLocus": task["genomic_locus"],
            "chr": task["chromosome"],
            "start": task["start"],
            "end": task["end"],
        }
        for task in tasks if task["locus_id"] in allowed_loci
    }
    process_finemap_output(
        raw_dir=dirs["loci"],
        inter_dir=str(dirs["inter"]),
        final_dir=str(dirs["flames"]),
        allowed_loci=allowed_loci,
        locus_metadata=locus_metadata,
        genome_build=genome_build,
        target_coverage=runtime_defaults["credible_set_coverage"],
        coverage_tolerance=runtime_defaults["coverage_tolerance"],
    )
    pd.DataFrame(list(locus_metadata.values())).to_csv(
        dirs["flames"] / "genomic_loci.tsv", sep="\t", index=False
    )
    progress.record(
        "pipeline", "output_formatting", 6, pipeline_total, "completed",
        f"flames_input={dirs['flames']}",
    )
    pd.DataFrame([
        {"software": "PLINK", "version": plink_version, "path": args.plink},
        {"software": "bgenix", "version": _software_version("bgenix"), "path": shutil.which("bgenix")},
        {"software": "LDstore", "version": _software_version("ldstore"), "path": shutil.which("ldstore")},
        {"software": "FINEMAP", "version": _software_version("finemap"), "path": shutil.which("finemap")},
        {"software": "genome_build", "version": genome_build, "path": ""},
        {
            "software": "credible_set_coverage",
            "version": str(runtime_defaults["credible_set_coverage"]),
            "path": "",
        },
    ]).to_csv(dirs["root"] / "software_versions.tsv", sep="\t", index=False)

    n_successful = len(successful_results)
    n_failed = len(results) - n_successful
    progress.record(
        "pipeline", "complete", pipeline_total, pipeline_total, "success",
        f"successful={n_successful}; failed={n_failed}; elapsed_seconds={time() - start_time:.3f}",
    )
    logger.info(
        "Pipeline finished in %.2f seconds (%d/%d successful)",
        time() - start_time,
        n_successful,
        len(results),
    )
    return {
        "status": "success",
        "output_dir": str(dirs["root"]),
        "flames_input": str(dirs["flames"]),
        "n_attempted": len(results),
        "n_successful": n_successful,
        "n_failed": n_failed,
    }


def run_finemap_pipeline(args):
    """Coordinate FINEMAP and record success, failure, percentage, and remaining work."""
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(log_file=outdir / "pipeline_summary.log")
    progress = ProgressRecorder(outdir / "pipeline_progress.tsv", logger=logger)
    try:
        result = _run_finemap_pipeline(args, progress)
        return {
            "status": result["status"],
            "output_dir": result["output_dir"],
            "flames_input": result["flames_input"],
            "n_attempted": result["n_attempted"],
            "n_successful": result["n_successful"],
            "n_failed": result["n_failed"],
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
