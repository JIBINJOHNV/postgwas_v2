import gzip
import hashlib
import logging
import math
import os
import shutil
import subprocess
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.core.execution.runtime import validate_path, safe_thread_count
import pandas as pd
from postgwas.cli.compute import resolve_compute_args
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

from postgwas.modules.fine_mapping.defaults import (
    DEFAULT_BASES_PER_KILOBASE,
    DEFAULT_CREDIBLE_SET_COVERAGE,
    DEFAULT_GENOME_BUILD,
    DEFAULT_MHC_CHROMOSOME,
    DEFAULT_MHC_END,
    DEFAULT_MHC_START,
    DEFAULT_MIN_POSITION,
    DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS,
    DEFAULT_SUSIE_RAM_PER_WORKER_GB,
    DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL,
    DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS,
    SUPPORTED_GENOME_BUILDS,
    get_finemap_defaults,
)
from postgwas.modules.fine_mapping.progress import ProgressRecorder, write_run_configuration
from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.output_layout import (
    cleanup_successful_workers,
    resolve_output_paths,
)
from postgwas.modules.fine_mapping.resource_guard import (
    VARIANT_LIMIT_FAILURE_REASON,
    dense_ld_resource_qc,
)

# SuSiE backend (Python) ---------------------
from postgwas.modules.fine_mapping.engines.susie.main import (
    validate_locus_file,
    resolve_susie_r_runtime,
    run_susie,
)
from postgwas.modules.fine_mapping.engines.susie.summary_preparation import (
    prepare_locus_summary_statistics,
)


logger = logging.getLogger("postgwas.modules.fine_mapping.engines.susie")

SUSIE_STATUS_SUCCESS = "success"
SUSIE_STATUS_NO_CREDIBLE_SETS = "completed_no_credible_sets"


def _merge_text_tables(input_files, output_file, compressed=False):
    """Stream-merge tables while retaining exactly one validated header."""
    input_files = [Path(p) for p in input_files if Path(p).is_file() and Path(p).stat().st_size > 0]
    if not input_files:
        return False

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    open_out = gzip.open if compressed else open
    expected_header = None

    with open_out(output_file, "wt") as dst:
        for input_file in input_files:
            open_in = gzip.open if input_file.suffix == ".gz" else open
            with open_in(input_file, "rt") as src:
                header = src.readline()
                if not header:
                    continue
                normalized = header.rstrip("\r\n")
                if expected_header is None:
                    expected_header = normalized
                    dst.write(normalized + "\n")
                elif normalized != expected_header:
                    raise ValueError(
                        f"Worker output schema mismatch: {input_file}; "
                        f"expected '{expected_header}', found '{normalized}'"
                    )
                for line in src:
                    if line.rstrip("\r\n") == expected_header:
                        continue
                    dst.write(line)
    return expected_header is not None


def _write_flames_index(
    worker_dirs,
    flames_dir,
    annotation_dir,
    index_filename,
    annotation_prefix,
):
    """Create one FLAMES index row per exported configured credible set."""
    rows = []
    for worker_dir in worker_dirs:
        row_file = Path(worker_dir) / "flames_input" / "indexfile_rows.tsv"
        if not row_file.is_file() or row_file.stat().st_size == 0:
            continue
        try:
            rows.append(pd.read_csv(row_file, sep="\t", dtype=str))
        except pd.errors.EmptyDataError:
            continue

    flames_dir = Path(flames_dir)
    flames_dir.mkdir(parents=True, exist_ok=True)
    index_path = flames_dir / index_filename
    if index_path.exists():
        index_path.unlink()
    if not rows:
        return None

    index_df = pd.concat(rows, ignore_index=True).drop_duplicates()
    required_columns = ("Filename", "GenomicLocus")
    missing_columns = set(required_columns) - set(index_df.columns)
    if missing_columns:
        raise ValueError(
            "SuSiE FLAMES index fragments are missing required columns: %s"
            % ", ".join(sorted(missing_columns))
        )
    if index_df[list(required_columns)].isna().any().any():
        raise ValueError("SuSiE FLAMES index fragments contain empty required values")
    index_df["Filename"] = index_df["Filename"].map(
        lambda value: str((flames_dir / Path(value).name).resolve())
    )
    index_df = index_df[list(required_columns)]

    annotation_dir = Path(annotation_dir)
    index_df["Annotfiles"] = index_df["Filename"].map(
        lambda value: str(
            (
                annotation_dir
                / f"{annotation_prefix}{Path(value).stem}.txt"
            ).resolve()
        )
    )

    index_df.to_csv(index_path, sep="\t", index=False)
    return index_path


def _summarize_credible_set_handoff(index_file, successful_loci):
    """Distinguish converged SuSiE fits from retained credible sets."""
    successful_loci = {str(locus) for locus in successful_loci}
    if index_file is None:
        return {
            "status": SUSIE_STATUS_NO_CREDIBLE_SETS,
            "flames_input": None,
            "credible_set_loci": set(),
            "n_credible_sets": 0,
            "n_loci_with_credible_sets": 0,
            "n_converged_without_credible_sets": len(successful_loci),
        }

    index_file = Path(index_file)
    try:
        index = pd.read_csv(index_file, sep="\t", dtype=str)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(
            "Cannot validate generated SuSiE FLAMES index %s: %s"
            % (index_file, exc)
        ) from exc
    required_columns = ("Filename", "GenomicLocus", "Annotfiles")
    missing_columns = set(required_columns) - set(index.columns)
    if missing_columns or index.empty:
        detail = (
            "no credible-set rows"
            if index.empty
            else "missing columns: %s" % ", ".join(sorted(missing_columns))
        )
        raise ValueError(
            "Generated SuSiE FLAMES index is invalid (%s): %s"
            % (detail, index_file)
        )
    if index[list(required_columns)].isna().any().any():
        raise ValueError(
            "Generated SuSiE FLAMES index contains empty required values: %s"
            % index_file
        )
    if index["Filename"].duplicated().any():
        raise ValueError(
            "Generated SuSiE FLAMES index contains duplicate credible-set "
            "files: %s" % index_file
        )
    for filename in index["Filename"]:
        credible_set_file = Path(filename)
        if not credible_set_file.is_absolute():
            credible_set_file = index_file.parent / credible_set_file
        if (
            not credible_set_file.is_file()
            or credible_set_file.stat().st_size == 0
        ):
            raise FileNotFoundError(
                "Generated SuSiE FLAMES index references a missing or empty "
                "credible-set file: %s" % credible_set_file
            )

    credible_set_loci = set(index["GenomicLocus"].astype(str))
    unexpected_loci = credible_set_loci - successful_loci
    if unexpected_loci:
        raise ValueError(
            "SuSiE FLAMES index contains loci not recorded as converged: %s"
            % ", ".join(sorted(unexpected_loci))
        )
    return {
        "status": SUSIE_STATUS_SUCCESS,
        "flames_input": str(index_file.parent),
        "credible_set_loci": credible_set_loci,
        "n_credible_sets": int(len(index)),
        "n_loci_with_credible_sets": len(credible_set_loci),
        "n_converged_without_credible_sets": len(
            successful_loci - credible_set_loci
        ),
    }


def _summarize_recovery_audit(audit_file):
    """Validate and summarize the merged per-attempt SuSiE recovery audit."""
    audit_file = Path(audit_file)
    if not audit_file.is_file():
        raise FileNotFoundError(f"SuSiE recovery audit TSV is missing: {audit_file}")
    try:
        audit = pd.read_csv(audit_file, sep="\t", dtype=str)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(f"Cannot read SuSiE recovery audit TSV {audit_file}: {exc}") from exc
    required = {"stage", "status", "reason", "repairable", "genomic_locus"}
    missing = required - set(audit.columns)
    if missing:
        raise ValueError(
            "SuSiE recovery audit is missing required columns: "
            + ", ".join(sorted(missing))
        )
    recovery_fit = audit["stage"].eq("recovery_fit")
    repairable = audit["repairable"].str.lower().eq("true")
    ld_validation = audit["stage"].isin(["ld_validation", "ld_revalidation"])
    failed = audit["status"].eq("failed")
    return {
        "table": audit,
        "n_recovery_successes": int((recovery_fit & audit["status"].eq("success")).sum()),
        "n_recovery_failures": int((recovery_fit & failed).sum()),
        "n_fatal_ld_failures": int((ld_validation & failed & ~repairable).sum()),
    }


def _log_locus_reason_counts(qc_df, log_method):
    """Print compact warning and failure reason counts from canonical locus QC."""
    for column, label in (
        ("warning_reason", "warning reason"),
        ("failure_reason", "failure reason"),
    ):
        if column not in qc_df:
            continue
        reasons = qc_df[column].dropna().astype(str)
        reasons = reasons[reasons.ne("")]
        for reason, count in reasons.value_counts().sort_index().items():
            log_method("  %s [%s]: %d", label, reason, int(count))


def merge_workers(
    worker_dirs,
    output_paths,
    sample_id,
    recovery_audit_filename,
    index_filename,
    annotation_prefix,
    log=None,
):
    """Merge current worker outputs while logging every copied artifact."""
    log = log or logger
    artifact_destinations = {
        "plots": Path(output_paths["diagnostic_plots_directory"]),
        "flames_input": Path(output_paths["primary_flames_directory"]),
        "locus_files": Path(output_paths["primary_credible_sets_directory"]),
        "rds_files": Path(output_paths["fitted_models_directory"]),
        "logs": Path(output_paths["locus_logs_directory"]),
        "ld_matrix_related": Path(output_paths["ld_diagnostics_directory"]),
    }
    for w in worker_dirs:
        w = Path(w)
        log.info("[STAGE] stage=worker_merge status=started worker=%s", w)
        for sub, dst in artifact_destinations.items():
            src = w / sub

            if not src.exists():
                continue
            for f in src.iterdir():
                if f.name == "indexfile_rows.tsv":
                    continue
                out = dst / f.name
                if f.is_file():
                    dst.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, out)
                    log.info(
                        "[STAGE] stage=worker_merge status=copied source=%s destination=%s",
                        f,
                        out,
                    )
        log.info("[STAGE] stage=worker_merge status=completed worker=%s", w)
    log.info("[STAGE] stage=table_merge status=started")
    table_specs = [
        (
            f"{sample_id}_SuSiE_QC_summary.tsv",
            output_paths["susie_qc_file"],
            False,
        ),
        (
            f"{sample_id}_SuSiE_failed_loci.tsv",
            output_paths["susie_failed_loci_file"],
            False,
        ),
        (
            f"{sample_id}_SUSIE_combined_credibleset.csv",
            output_paths["susie_combined_credible_sets_file"],
            False,
        ),
        (
            f"{sample_id}_SUSIE_combined_results.csv.gz",
            output_paths["susie_combined_results_file"],
            True,
        ),
        (
            recovery_audit_filename,
            Path(output_paths["quality_control_directory"])
            / recovery_audit_filename,
            False,
        ),
    ]
    for filename, destination, compressed in table_specs:
        inputs = [Path(w) / filename for w in worker_dirs]
        merged = _merge_text_tables(inputs, destination, compressed=compressed)
        log.info(
            "[STAGE] stage=table_merge status=%s output=%s",
            "completed" if merged else "skipped",
            destination,
        )

    log_inputs = [Path(w) / f"{sample_id}_run_susie.log" for w in worker_dirs]
    log_inputs = [p for p in log_inputs if p.is_file()]
    combined_log = Path(output_paths["susie_combined_log_file"])
    combined_log.parent.mkdir(parents=True, exist_ok=True)
    with open(combined_log, "w") as dst:
        for input_file in log_inputs:
            with open(input_file, "r") as src:
                shutil.copyfileobj(src, dst)

    index_file = _write_flames_index(
        worker_dirs,
        output_paths["primary_flames_directory"],
        output_paths["primary_flames_annotations_directory"],
        index_filename,
        annotation_prefix,
    )
    log.info(
        "[STAGE] stage=flames_index status=%s output=%s",
        "completed" if index_file else "skipped",
        index_file or "none",
    )

    worker_progress = [Path(w) / "susie_locus_progress.tsv" for w in worker_dirs]
    _merge_text_tables(
        worker_progress,
        output_paths["susie_locus_progress_file"],
        compressed=False,
    )
    for source_name, destination in [
        ("run_configuration_r.tsv", output_paths["susie_r_configuration_file"]),
        (
            "software_versions_r.tsv",
            output_paths["susie_r_software_versions_file"],
        ),
    ]:
        source = next(
            (Path(w) / source_name for w in worker_dirs if (Path(w) / source_name).is_file()),
            None,
        )
        if source is not None:
            destination = Path(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            log.info(
                "[STAGE] stage=worker_metadata status=copied source=%s destination=%s",
                source,
                destination,
            )
    return index_file

def split_locus_file(
    locus_file: str,
    outdir: Path,
    n_chunks: int,
    sep: str = "\t",
    finemap_skip_mhc: bool = False,
    mhc_chr: str = DEFAULT_MHC_CHROMOSOME,
    mhc_start: int = DEFAULT_MHC_START,
    mhc_end: int = DEFAULT_MHC_END,
    log=None,
):
    """
    Split locus file into n_chunks row-wise.
    If finemap_skip_mhc=True, remove MHC loci BEFORE splitting.
    Ensures every chunk contains at least one locus.
    Returns list of chunk locus file paths.
    """
    log = log or logger
    outdir.mkdir(parents=True, exist_ok=True)

    # Auto-detect delimiter safely
    df = pd.read_csv(locus_file, sep=None, engine="python")

    # --------------------------------------------------
    # Optional MHC filtering (BEFORE splitting)
    # --------------------------------------------------
    if finemap_skip_mhc:
        before = len(df)

        df = df[
            ~(
                (df["CHR"].astype(str) == str(mhc_chr)) &
                (df["START"].astype(int) < int(mhc_end)) &
                (df["END"].astype(int) > int(mhc_start))
            )
        ].reset_index(drop=True)

        after = len(df)

        if after == 0:
            raise ValueError(
                "All loci removed after MHC filtering — nothing left to finemap."
            )

        log.info(
            "[STAGE] stage=mhc_filter status=completed before=%d after=%d",
            before,
            after,
        )

    # --------------------------------------------------
    # Adjust number of chunks to avoid empty files
    # --------------------------------------------------
    n = len(df)
    if n == 0:
        raise ValueError("The locus file contains no loci — nothing to finemap.")
    if int(n_chunks) < 1:
        raise ValueError("n_chunks must be at least 1.")
    n_chunks = min(n_chunks, n)   # 🔑 at least one locus per chunk
    chunk_size = math.ceil(n / n_chunks)

    chunk_files = []

    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, n)

        if start >= end:
            break  # safety

        chunk_df = df.iloc[start:end]

        # Extra safety: never write empty chunks
        if chunk_df.empty:
            continue

        chunk_path = outdir / f"locus_chunk_{i+1}.tsv"
        chunk_df.to_csv(chunk_path, sep=sep, index=False)
        chunk_files.append(chunk_path)
        log.info(
            "[STAGE] stage=locus_split status=written chunk=%s rows=%d",
            chunk_path,
            len(chunk_df),
        )

    log.info(
        "[STAGE] stage=locus_split status=completed chunks=%d loci=%d",
        len(chunk_files),
        n,
    )

    return chunk_files


def run_susie_worker(
    locus_chunk: str,
    summary_statistics_manifest: str,
    worker_id: int,
    args,
    base_output_dir: Path,
):
    worker_outdir = base_output_dir / f"worker_{worker_id}"
    worker_outdir.mkdir(parents=True, exist_ok=True)

    # Worker directories are pipeline-owned.  Clear only known generated
    # artifacts so a rerun cannot merge a stale locus or FLAMES file.
    for subdir in [
        "plots", "flames_input", "locus_files", "rds_files", "logs",
        "ld_matrix_related", "output",
    ]:
        generated = worker_outdir / subdir
        if generated.exists():
            shutil.rmtree(generated)
    for generated in worker_outdir.glob(f"{args.dataset_id}_*"):
        if generated.is_file():
            generated.unlink()

    return run_susie(
        locus_file=locus_chunk,                     # 🔑 CHUNKED LOCUS FILE
        sumstat_manifest=summary_statistics_manifest,
        sample_id=args.dataset_id,
        ld_ref=args.finemap_ld_reference,
        plink=args.plink,
        output_folder=str(worker_outdir),
        resolved_configuration_file=args.resolved_susie_configuration_file,
        rscript=args.rscript,
        r_environment=args.susie_r_environment,
    )


def _run_parallel_susie(args, progress, screen=None):
    """Implement parallel SuSiE with documented stages and resource controls."""
    resolve_compute_args(args)
    runtime_defaults = dict(args.fine_mapping_runtime_defaults)
    pipeline_total = DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL
    output_folder = Path(args.output_directory).resolve()
    output_paths = resolve_output_paths(
        output_folder,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    resolved_configuration_file = write_run_configuration(
        output_paths["resolved_susie_configuration_file"],
        args.resolved_fine_mapping_configuration,
    )
    args.resolved_susie_configuration_file = str(
        resolved_configuration_file.resolve()
    )
    resolved_configuration_sha256 = hashlib.sha256(
        resolved_configuration_file.read_bytes()
    ).hexdigest()
    progress.record(
        "pipeline", "initialization", 1, pipeline_total, "completed",
        f"output_dir={output_folder}",
    )
    if screen:
        screen.complete(1, [
            ("analysis", "Engine", "SuSiE-RSS"),
            ("info", "Dataset", args.dataset_id),
        ])
    final_plink_path = None

    logger.info("[STAGE] stage=dependency_validation status=started")
    r_runtime = resolve_susie_r_runtime(
        args.rscript,
        runtime_defaults["tool_version_timeout_seconds"],
    )
    args.rscript = r_runtime["rscript"]
    args.susie_r_environment = r_runtime["environment"]
    if args.plink:
        if Path(args.plink).exists():
            final_plink_path = args.plink
        else:
            logger.warning(
                "[STAGE] stage=dependency_validation status=warning "
                "reason=provided_plink_path_missing path=%s",
                args.plink,
            )
            final_plink_path = None

    if not final_plink_path:
        detected_plink = shutil.which("plink") or shutil.which("plink2")
        if detected_plink:
            final_plink_path = detected_plink

    if not final_plink_path:
        raise EnvironmentError(
            "PLINK executable not found; add plink/plink2 to PATH or provide --plink"
        )

    args.plink = final_plink_path
    try:
        plink_version = subprocess.run(
            [str(args.plink), "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=float(runtime_defaults["tool_version_timeout_seconds"]),
        )
        version_output = (plink_version.stdout or plink_version.stderr).strip()
        args.plink_version = version_output.splitlines()[0]
    except Exception:
        args.plink_version = "unknown"
    logger.info(
        "[STAGE] stage=dependency_validation status=completed "
        "plink=%s plink_version=%s rscript=%s r_version=%s r_libraries=%s",
        args.plink,
        args.plink_version,
        args.rscript,
        r_runtime["version"],
        os.pathsep.join(r_runtime["library_paths"]),
    )
    progress.record(
        "pipeline", "dependency_validation", 2, pipeline_total, "completed",
        f"PLINK={args.plink_version}",
    )
    if screen:
        screen.complete(2, [
            ("success", "External tools", "validated"),
            ("info", "PLINK", args.plink_version),
            ("analysis", "Genome build", args.genome_build),
        ])

    args.threads = safe_thread_count(
        requested_threads=args.threads,
        gb_per_thread=args.minimum_memory_per_worker_gb,
    )

    if args.finemap_include_mhc:
        args.finemap_skip_mhc = False

    genome_build = args.genome_build
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError("genome_build must be GRCh37 or GRCh38")
    validate_locus_file(args.locus_file, locus_type=args.locus_type)
    logger.info("[STAGE] stage=locus_preparation status=started")
    try:
        loci_df = pd.read_csv(args.locus_file, sep=None, engine="python")
    except Exception as e:
        raise ValueError(f"Could not read locus file: {e}") from e

    loci_df.columns = [x.strip().upper() for x in loci_df.columns]
    flank_bp = int(args.window_kb * runtime_defaults["bases_per_kilobase"])
    if flank_bp < 0:
        raise ValueError("window_kb cannot be negative")

    if args.locus_type == "point":
        if 'POS' in loci_df.columns and "CHR" in loci_df.columns:
            loci_df['POS'] = pd.to_numeric(loci_df['POS'])
            loci_df['START'] = (
                loci_df['POS'] - flank_bp
            ).clip(lower=int(runtime_defaults["minimum_position"])).astype(int)
            loci_df['END'] = (loci_df['POS'] + flank_bp).astype(int)
            modified_name = output_paths["point_locus_windows_file"]
            modified_name.parent.mkdir(parents=True, exist_ok=True)
            loci_df.to_csv(modified_name, sep="\t", index=False)
            args.locus_file = str(modified_name)
            logger.info(
                "[STAGE] stage=locus_preparation status=windowed window_kb=%s output=%s",
                args.window_kb,
                modified_name,
            )
        else:
            raise ValueError(
                "locus_type is 'point', but POS or CHR is missing; "
                f"found {list(loci_df.columns)}"
            )
    elif args.locus_type == "range":
        required = ["CHR", "START", "END"]
        if all(col in loci_df.columns for col in required):
            loci_df['START'] = pd.to_numeric(loci_df['START']).astype(int)
            loci_df['END'] = pd.to_numeric(loci_df['END']).astype(int)
            if flank_bp > 0:
                loci_df["START"] = (
                    loci_df["START"] - flank_bp
                ).clip(lower=int(runtime_defaults["minimum_position"])).astype(int)
                loci_df["END"] = (loci_df["END"] + flank_bp).astype(int)
                modified_name = output_paths["range_locus_windows_file"]
                modified_name.parent.mkdir(parents=True, exist_ok=True)
                loci_df.to_csv(modified_name, sep="\t", index=False)
                args.locus_file = str(modified_name)
                logger.info(
                    "[STAGE] stage=locus_preparation status=extended window_kb=%s output=%s",
                    args.window_kb,
                    modified_name,
                )
            else:
                logger.info(
                    "[STAGE] stage=locus_preparation status=exact_ranges source=%s",
                    args.locus_file,
                )
        else:
            raise ValueError(
                "locus_type is 'range', but CHR, START, or END is missing; "
                f"found {list(loci_df.columns)}"
            )
    logger.info(
        "[STAGE] stage=locus_preparation status=completed loci=%d file=%s",
        len(loci_df),
        args.locus_file,
    )
    progress.record(
        "pipeline", "locus_preparation", 3, pipeline_total, "completed",
        f"loci={len(loci_df)}; file={args.locus_file}",
    )
    if screen:
        screen.complete(3, [
            ("count", "Input loci", len(loci_df)),
            ("analysis", "Boundary flank", f"{args.window_kb:,} kb"),
            (
                "analysis",
                "MHC policy",
                "excluded" if args.finemap_skip_mhc else "retained",
            ),
        ])

    locus_chunk_dir = output_paths["locus_chunks_directory"]
    locus_chunks = split_locus_file(
        locus_file=args.locus_file,
        outdir=locus_chunk_dir,
        n_chunks=args.threads,
        sep="\t",
        finemap_skip_mhc=args.finemap_skip_mhc,
        mhc_chr=args.finemap_mhc_chromosome,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
        log=logger,
    )
    if not locus_chunks:
        raise RuntimeError("No locus chunks were generated — SuSiE was not run.")
    logger.info(
        "[STAGE] stage=summary_statistics_preparation status=started "
        "policy=%s",
        args.summary_statistics_preparation["policy"],
    )
    prepared_summary_statistics = prepare_locus_summary_statistics(
        source_file=args.susie_input_file,
        locus_files=locus_chunks,
        chromosome_cache_directory=output_paths[
            "chromosome_cache_directory"
        ],
        locus_input_directory=output_paths[
            "locus_summary_statistics_directory"
        ],
        manifest_directory=output_paths["summary_statistics_directory"],
        quality_control_directory=output_paths["quality_control_directory"],
        settings=args.summary_statistics_preparation,
        existing_chromosome_manifest=getattr(
            args, "summary_statistics_chromosome_manifest", None
        ),
    )
    logger.info(
        "[STAGE] stage=summary_statistics_preparation status=completed "
        "loci=%d variant_memberships=%d cache_reused=%s manifest=%s",
        prepared_summary_statistics["n_loci"],
        prepared_summary_statistics["n_variant_memberships"],
        prepared_summary_statistics["cache_reused"],
        prepared_summary_statistics["locus_manifest"],
    )
    prepared_manifest = pd.read_csv(
        prepared_summary_statistics["locus_manifest"], sep="\t"
    )
    largest_locus_variants = int(prepared_manifest["n_variants"].max())
    largest_locus_resource_qc = dense_ld_resource_qc(
        largest_locus_variants,
        args.maximum_variants_per_locus,
        args.susie_ld_peak_matrix_multiplier,
        args.minimum_memory_per_worker_gb,
    )
    memory_adjusted_threads = safe_thread_count(
        requested_threads=len(locus_chunks),
        gb_per_thread=max(
            args.minimum_memory_per_worker_gb,
            largest_locus_resource_qc["estimated_peak_ld_memory_gb"],
        ),
    )
    if memory_adjusted_threads < len(locus_chunks):
        shutil.rmtree(locus_chunk_dir)
        locus_chunks = split_locus_file(
            locus_file=args.locus_file,
            outdir=locus_chunk_dir,
            n_chunks=memory_adjusted_threads,
            sep="\t",
            finemap_skip_mhc=args.finemap_skip_mhc,
            mhc_chr=args.finemap_mhc_chromosome,
            mhc_start=args.finemap_mhc_start,
            mhc_end=args.finemap_mhc_end,
            log=logger,
        )
    write_run_configuration(
        output_paths["run_configuration_file"],
        {
            "engine": "SuSiE",
            "resolved_configuration": args.resolved_fine_mapping_configuration,
            "resolved_configuration_file": args.resolved_susie_configuration_file,
            "resolved_configuration_sha256": resolved_configuration_sha256,
            "inputs": {
                "locus_file": args.locus_file,
                "summary_statistics": args.susie_input_file,
                "summary_statistics_locus_manifest": (
                    prepared_summary_statistics["locus_manifest"]
                ),
                "summary_statistics_chromosome_manifest": (
                    prepared_summary_statistics["chromosome_manifest"]
                ),
                "ld_reference": args.finemap_ld_reference,
                "output_directory": output_folder,
                "sample_id": args.dataset_id,
            },
            "analysis_parameters": {
                "genome_build": genome_build,
                "locus_type": args.locus_type,
                "window_kb": args.window_kb,
                "lp_threshold": args.lp_threshold,
                "L": args.L,
                "minimum_purity": args.minimum_purity,
                "credible_set_coverage": args.resolved_fine_mapping_configuration[
                    "credible_set_coverage"
                ],
                "skip_mhc": args.finemap_skip_mhc,
                "mhc_chromosome": args.finemap_mhc_chromosome,
                "mhc_start": args.finemap_mhc_start,
                "mhc_end": args.finemap_mhc_end,
                "maximum_variants_per_locus": args.maximum_variants_per_locus,
            },
            "resource_parameters": {
                "selected_workers": len(locus_chunks),
                "minimum_memory_per_worker_gb": args.minimum_memory_per_worker_gb,
                "ld_timeout_seconds": args.ld_timeout_seconds,
                "susie_timeout_seconds": args.susie_timeout_seconds,
                "recovery_audit_filename": args.recovery_audit_filename,
                "sample_size_policy": args.sample_size_policy,
                "sample_size_summary_statistic": args.sample_size_summary_statistic,
                "sample_size_relative_range_warning_threshold": (
                    args.sample_size_relative_range_warning_threshold
                ),
                "susie_ld_peak_matrix_multiplier": (
                    args.susie_ld_peak_matrix_multiplier
                ),
                "largest_locus_variants": largest_locus_variants,
                "largest_locus_estimated_peak_ld_memory_gb": (
                    largest_locus_resource_qc["estimated_peak_ld_memory_gb"]
                ),
                "summary_statistics_preparation": (
                    args.summary_statistics_preparation
                ),
                "output_layout": args.fine_mapping_output_layout,
            },
            "software": {
                "plink_path": args.plink,
                "plink_version": args.plink_version,
                "rscript_path": args.rscript,
                "r_version": r_runtime["version"],
                "r_library_paths": r_runtime["library_paths"],
            },
            "defaults": runtime_defaults,
        },
    )
    progress.record(
        "pipeline", "locus_split", 4, pipeline_total, "completed",
        f"chunks={len(locus_chunks)}",
    )
    if screen:
        screen.complete(4, [
            ("count", "Prepared loci", prepared_summary_statistics["n_loci"]),
            (
                "count",
                "Variant memberships",
                prepared_summary_statistics["n_variant_memberships"],
            ),
            ("analysis", "Parallel workers", len(locus_chunks)),
        ])

    ctx = mp.get_context("spawn")

    results = []
    with ProcessPoolExecutor(
        max_workers=len(locus_chunks),
        mp_context=ctx
    ) as executor:

        futures = {
            executor.submit(
                run_susie_worker,
                str(chunk),
                prepared_summary_statistics["locus_manifest"],
                i + 1,
                args,
                output_paths["workers_directory"],
            ): chunk
            for i, chunk in enumerate(locus_chunks)
        }

        for completed_count, fut in enumerate(as_completed(futures), start=1):
            chunk = futures[fut]
            try:
                res = fut.result()
                results.append(res)
            except Exception as e:
                progress.record(
                    "workers", "susie_execution", completed_count,
                    len(locus_chunks), "failed", f"chunk={chunk}; reason={e}",
                )
                logger.exception(
                    "[STAGE] stage=susie_execution status=failed chunk=%s reason=%s",
                    chunk,
                    e,
                )
                raise
            progress.record(
                "workers", "susie_execution", completed_count,
                len(locus_chunks), res.get("status", "unknown"), f"chunk={chunk}",
            )

    progress.record(
        "pipeline", "worker_execution", 5, pipeline_total, "completed",
        f"workers={len(results)}",
    )
    if screen:
        screen.complete(5, [
            ("success", "Workers completed", len(results)),
            ("analysis", "Locus batches", len(locus_chunks)),
        ])

    current_worker_dirs = [
        output_paths["workers_directory"] / f"worker_{index + 1}"
        for index in range(len(locus_chunks))
    ]
    index_file = merge_workers(
        worker_dirs=current_worker_dirs,
        output_paths=output_paths,
        sample_id=args.dataset_id,
        recovery_audit_filename=args.recovery_audit_filename,
        index_filename=args.overlap_resolution["index_filename"],
        annotation_prefix=args.overlap_resolution["annotation_prefix"],
        log=logger,
    )

    qc_file = output_paths["susie_qc_file"]
    combined_file = output_paths["susie_combined_results_file"]
    recovery_audit_file = (
        output_paths["quality_control_directory"] / args.recovery_audit_filename
    )
    if not qc_file.is_file() or not recovery_audit_file.is_file():
        raise RuntimeError(
            "SuSiE workers completed but required reporting outputs are missing. "
            f"Expected {qc_file.name} and {recovery_audit_file.name}."
        )

    qc_df = pd.read_csv(qc_file, sep="\t")
    recovery_summary = _summarize_recovery_audit(recovery_audit_file)
    n_success = int(qc_df.get("converged", pd.Series(dtype=bool)).astype(str).str.lower().eq("true").sum())
    n_failed = int(len(qc_df) - n_success)
    n_warning_loci = int(
        qc_df.get("warning_reason", pd.Series(dtype=str)).fillna("").ne("").sum()
    )
    n_variant_limit_failures = int(
        qc_df.get("failure_reason", pd.Series(dtype=str))
        .fillna("")
        .eq(VARIANT_LIMIT_FAILURE_REASON)
        .sum()
    )
    recovery_successes = recovery_summary["n_recovery_successes"]
    recovery_failures = recovery_summary["n_recovery_failures"]
    fatal_ld_failures = recovery_summary["n_fatal_ld_failures"]
    if n_success == 0:
        logger.error("SuSiE final summary")
        logger.error("  loci attempted: %d", len(qc_df))
        logger.error("  successful loci: 0")
        logger.error("  failed or skipped loci: %d", n_failed)
        logger.error("  loci with warnings: %d", n_warning_loci)
        logger.error("  failed recovery fits: %d", recovery_failures)
        logger.error("  fatal LD validation failures: %d", fatal_ld_failures)
        logger.error("  variant-limit failures: %d", n_variant_limit_failures)
        if n_variant_limit_failures:
            logger.error(
                "  Override with --maximum-variants-per-locus COUNT or edit "
                "modules.fine_mapping.ld_resource_guard."
                "maximum_variants_per_locus in YAML"
            )
        logger.error("  QC summary TSV: %s", qc_file)
        logger.error("  recovery audit TSV: %s", recovery_audit_file)
        _log_locus_reason_counts(qc_df, logger.error)
        raise RuntimeError(
            "All SuSiE loci failed or were skipped; QC and recovery audit TSVs "
            "were finalized before failure."
        )
    if not combined_file.is_file():
        raise RuntimeError(
            "SuSiE recorded successful loci but the combined result is missing: "
            f"{combined_file}"
        )

    # FLAMES can consume these definitions with its -l option. Keep only loci
    # that both converged and yielded at least one retained credible set.
    used_loci = pd.read_csv(args.locus_file, sep=None, engine="python")
    used_loci.columns = [column.strip() for column in used_loci.columns]
    used_loci["CHR"] = used_loci["CHR"].astype(str).str.replace(
        r"^chr", "", regex=True, case=False
    )
    if "GenomicLocus" not in used_loci.columns:
        used_loci["GenomicLocus"] = (
            "chr" + used_loci["CHR"] + ":" + used_loci["START"].astype(str)
            + "-" + used_loci["END"].astype(str)
        )
    successful_loci = set(
        qc_df.loc[
            qc_df["converged"].astype(str).str.lower().eq("true"),
            "genomic_locus",
        ].astype(str)
    )
    handoff = _summarize_credible_set_handoff(index_file, successful_loci)
    used_loci = used_loci[
        used_loci["GenomicLocus"].astype(str).isin(handoff["credible_set_loci"])
    ]
    used_loci = used_loci.rename(columns={"CHR": "chr", "START": "start", "END": "end"})
    primary_flames_directory = output_paths["primary_flames_directory"]
    used_loci[["GenomicLocus", "chr", "start", "end"]].drop_duplicates().to_csv(
        primary_flames_directory
        / args.overlap_resolution["genomic_loci_filename"],
        sep="\t",
        index=False,
    )

    version_file = output_paths["susie_software_versions_file"]
    version_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"software": "PLINK", "version": args.plink_version, "path": str(args.plink)},
        {"software": "genome_build", "version": genome_build, "path": ""},
        {
            "software": "credible_set_coverage",
            "version": str(runtime_defaults["credible_set_coverage"]),
            "path": "",
        },
    ]).to_csv(version_file, sep="\t", index=False)

    progress.record(
        "pipeline", "merge_and_validation", 6, pipeline_total, "completed",
        f"successful={n_success}; failed={n_failed}; "
        f"credible_sets={handoff['n_credible_sets']}; "
        f"converged_without_credible_sets={handoff['n_converged_without_credible_sets']}",
    )
    if screen:
        screen.complete(6, [
            ("success", "Successful loci", n_success),
            ("warning" if n_failed else "info", "Failed or skipped", n_failed),
            ("success", "Primary credible sets", handoff["n_credible_sets"]),
            (
                "warning"
                if handoff["n_converged_without_credible_sets"]
                else "info",
                "Converged without sets",
                handoff["n_converged_without_credible_sets"],
            ),
        ])
    progress.record(
        "pipeline", "complete", pipeline_total, pipeline_total, handoff["status"],
        f"output_dir={output_folder}; successful={n_success}; failed={n_failed}; "
        f"credible_sets={handoff['n_credible_sets']}; "
        f"converged_without_credible_sets={handoff['n_converged_without_credible_sets']}",
    )
    logger.info(
        "SuSiE fine-mapping completed: output=%s status=%s successful=%d "
        "failed=%d credible_sets=%d converged_without_credible_sets=%d",
        output_folder,
        handoff["status"],
        n_success,
        n_failed,
        handoff["n_credible_sets"],
        handoff["n_converged_without_credible_sets"],
    )
    logger.info("SuSiE final summary")
    logger.info("  loci attempted: %d", len(qc_df))
    logger.info("  successful loci: %d", n_success)
    logger.info("  failed or skipped loci: %d", n_failed)
    logger.info("  loci with warnings: %d", n_warning_loci)
    logger.info("  successful recovery fits: %d", recovery_successes)
    logger.info("  failed recovery fits: %d", recovery_failures)
    logger.info("  fatal LD validation failures: %d", fatal_ld_failures)
    logger.info("  variant-limit failures: %d", n_variant_limit_failures)
    if n_variant_limit_failures:
        logger.info(
            "  Override with --maximum-variants-per-locus COUNT or edit "
            "modules.fine_mapping.ld_resource_guard.maximum_variants_per_locus "
            "in YAML"
        )
    logger.info("  recovery audit TSV: %s", recovery_audit_file)
    logger.info(
        "  exact locus summary-statistics files: %d",
        prepared_summary_statistics["n_loci"],
    )
    logger.info(
        "  summary-statistics chromosome cache reused: %s",
        prepared_summary_statistics["cache_reused"],
    )
    logger.info(
        "  summary-statistics preparation TSV: %s",
        prepared_summary_statistics["preparation_summary"],
    )
    _log_locus_reason_counts(qc_df, logger.info)
    warning_reasons = qc_df.get(
        "warning_reason", pd.Series(dtype=str)
    ).dropna().astype(str)
    warning_reasons = warning_reasons[warning_reasons.ne("")]
    warning_reason_counts = {
        str(reason): int(count)
        for reason, count in warning_reasons.value_counts().sort_index().items()
    }
    failure_reasons = qc_df.get(
        "failure_reason", pd.Series(dtype=str)
    ).dropna().astype(str)
    failure_reasons = failure_reasons[failure_reasons.ne("")]
    failure_reason_counts = {
        str(reason): int(count)
        for reason, count in failure_reasons.value_counts().sort_index().items()
    }
    workers_cleaned = cleanup_successful_workers(
        output_paths["workers_directory"],
        enabled=args.fine_mapping_output_layout["cleanup_successful_workers"],
    )
    if workers_cleaned:
        logger.info(
            "[STAGE] stage=worker_cleanup status=completed directory=%s",
            output_paths["workers_directory"],
        )
    if screen:
        screen.complete(7, [
            (
                "success"
                if handoff["status"] == SUSIE_STATUS_SUCCESS
                else "warning",
                "Primary analysis status",
                handoff["status"].replace("_", " "),
            ),
            (
                "warning" if n_warning_loci else "info",
                "Loci with warnings",
                n_warning_loci,
            ),
            ("success", "Primary credible sets", handoff["n_credible_sets"]),
            ("info", "Worker intermediates cleaned", workers_cleaned),
        ])
    return {
        "status": handoff["status"],
        "output_dir": str(output_folder),
        "flames_input": handoff["flames_input"],
        "credible_set_manifests": sorted(
            str(path) for path in primary_flames_directory.glob(
                "*_FLAMES_manifest.tsv"
            )
        ),
        "flames_index": str(index_file) if index_file else None,
        "locus_status": str(qc_file),
        "genomic_loci": str(
            primary_flames_directory
            / args.overlap_resolution["genomic_loci_filename"]
        ),
        "n_attempted": int(len(qc_df)),
        "n_successful": n_success,
        "n_failed": n_failed,
        "n_warnings": n_warning_loci,
        "warning_reason_counts": warning_reason_counts,
        "failure_reason_counts": failure_reason_counts,
        "n_credible_sets": handoff["n_credible_sets"],
        "n_loci_with_credible_sets": handoff["n_loci_with_credible_sets"],
        "n_converged_without_credible_sets": handoff[
            "n_converged_without_credible_sets"
        ],
        "recovery_audit": str(recovery_audit_file),
        "n_recovery_successes": recovery_successes,
        "n_recovery_failures": recovery_failures,
        "n_fatal_ld_failures": fatal_ld_failures,
        "n_variant_limit_failures": n_variant_limit_failures,
        "summary_statistics_chromosome_manifest": (
            prepared_summary_statistics["chromosome_manifest"]
        ),
        "summary_statistics_locus_manifest": (
            prepared_summary_statistics["locus_manifest"]
        ),
        "summary_statistics_preparation_summary": (
            prepared_summary_statistics["preparation_summary"]
        ),
    }


def run_parallel_susie(args, screen=None):
    """Coordinate SuSiE and retain existing return keys for downstream callers."""
    output_folder = Path(args.output_directory).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)
    output_paths = resolve_output_paths(
        output_folder,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    output_paths["run_metadata_directory"].mkdir(parents=True, exist_ok=True)
    with detailed_file_logging(
        "postgwas.modules.fine_mapping",
        output_paths["pipeline_log_file"],
        args.fine_mapping_logging["file_level"],
    ):
        progress = ProgressRecorder(
            output_paths["pipeline_progress_file"], logger=logger
        )
        try:
            result = _run_parallel_susie(args, progress, screen=screen)
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
                "n_loci_with_credible_sets": result[
                    "n_loci_with_credible_sets"
                ],
                "n_converged_without_credible_sets": result[
                    "n_converged_without_credible_sets"
                ],
                "recovery_audit": result["recovery_audit"],
                "n_recovery_successes": result["n_recovery_successes"],
                "n_recovery_failures": result["n_recovery_failures"],
                "n_fatal_ld_failures": result["n_fatal_ld_failures"],
                "n_variant_limit_failures": result[
                    "n_variant_limit_failures"
                ],
                "summary_statistics_chromosome_manifest": result[
                    "summary_statistics_chromosome_manifest"
                ],
                "summary_statistics_locus_manifest": result[
                    "summary_statistics_locus_manifest"
                ],
                "summary_statistics_preparation_summary": result[
                    "summary_statistics_preparation_summary"
                ],
            }
        except Exception as exc:
            latest = progress.latest.get("pipeline", {})
            completed = int(latest.get("completed", 0))
            progress.record(
                "pipeline", "failed", completed,
                DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL, "failed", str(exc)
            )
            logger.exception("SuSiE pipeline failed: %s", exc)
            raise
