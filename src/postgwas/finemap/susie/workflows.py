import gzip
import logging
import math
import shutil
import subprocess
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path, safe_thread_count
import pandas as pd
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

from postgwas.finemap.defaults import (
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
from postgwas.finemap.progress import ProgressRecorder, write_run_configuration

# SuSiE backend (Python) ---------------------
from postgwas.finemap.susie.main import (
    validate_locus_file,
    run_susie,
)


logger = logging.getLogger("postgwas.finemap.susie")


def setup_logging(log_file=None, verbose=True):
    """Configure SuSiE orchestration logs for console and persistent output."""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w"))
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


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


def _write_flames_index(worker_dirs, final_dir):
    """Create one FLAMES index row per exported 95% credible-set file."""
    rows = []
    for worker_dir in worker_dirs:
        row_file = Path(worker_dir) / "flames_input" / "indexfile_rows.tsv"
        if not row_file.is_file() or row_file.stat().st_size == 0:
            continue
        try:
            rows.append(pd.read_csv(row_file, sep="\t", dtype=str))
        except pd.errors.EmptyDataError:
            continue

    flames_dir = Path(final_dir) / "flames_input"
    annotation_dir = flames_dir / "annots"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    if not rows:
        credible_files = sorted(flames_dir.glob("*_CS_*.txt"))
        if not credible_files:
            return None
        index_df = pd.DataFrame({"Filename": [str(p.resolve()) for p in credible_files]})
    else:
        index_df = pd.concat(rows, ignore_index=True).drop_duplicates()
        index_df["Filename"] = index_df["Filename"].map(
            lambda value: str((flames_dir / Path(value).name).resolve())
        )
        keep_cols = [c for c in ["Filename", "GenomicLocus"] if c in index_df.columns]
        index_df = index_df[keep_cols]

    index_df["Annotfiles"] = index_df["Filename"].map(
        lambda value: str(
            (annotation_dir / f"annotated_{Path(value).stem}.txt").resolve()
        )
    )

    index_path = flames_dir / "indexfile.txt"
    index_df.to_csv(index_path, sep="\t", index=False)
    return index_path


def merge_workers(worker_dirs, final_dir, sample_id, log=None):
    """Merge current worker outputs while logging every copied artifact."""
    log = log or logger
    final_dir = Path(final_dir)
    for sub in [
        "plots", "flames_input", "locus_files",
        "rds_files", "logs", "ld_matrix_related"
    ]:
        (final_dir / sub).mkdir(parents=True, exist_ok=True)
    for w in worker_dirs:
        w = Path(w)
        log.info("[STAGE] stage=worker_merge status=started worker=%s", w)
        for sub in [
            "plots", "flames_input", "locus_files",
            "rds_files", "logs", "ld_matrix_related"
        ]:
            src = w / sub
            dst = final_dir / sub

            if not src.exists():
                continue
            for f in src.iterdir():
                if f.name == "indexfile_rows.tsv":
                    continue
                out = dst / f.name
                if f.is_file():
                    shutil.copy2(f, out)
                    log.info(
                        "[STAGE] stage=worker_merge status=copied source=%s destination=%s",
                        f,
                        out,
                    )
        log.info("[STAGE] stage=worker_merge status=completed worker=%s", w)
    log.info("[STAGE] stage=table_merge status=started")
    table_specs = [
        (f"{sample_id}_SuSiE_QC_summary.tsv", False),
        (f"{sample_id}_SuSiE_failed_loci.tsv", False),
        (f"{sample_id}_SUSIE_combined_credibleset.csv", False),
        (f"{sample_id}_SUSIE_combined_results.csv.gz", True),
    ]
    for filename, compressed in table_specs:
        inputs = [Path(w) / filename for w in worker_dirs]
        merged = _merge_text_tables(inputs, final_dir / filename, compressed=compressed)
        log.info(
            "[STAGE] stage=table_merge status=%s output=%s",
            "completed" if merged else "skipped",
            final_dir / filename,
        )

    log_inputs = [Path(w) / f"{sample_id}_run_susie.log" for w in worker_dirs]
    log_inputs = [p for p in log_inputs if p.is_file()]
    with open(final_dir / f"{sample_id}_run_susie.ALL.log", "w") as dst:
        for input_file in log_inputs:
            with open(input_file, "r") as src:
                shutil.copyfileobj(src, dst)

    index_file = _write_flames_index(worker_dirs, final_dir)
    log.info(
        "[STAGE] stage=flames_index status=%s output=%s",
        "completed" if index_file else "skipped",
        index_file or "none",
    )

    worker_progress = [Path(w) / "susie_locus_progress.tsv" for w in worker_dirs]
    _merge_text_tables(
        worker_progress,
        final_dir / f"{sample_id}_SuSiE_locus_progress.tsv",
        compressed=False,
    )
    for source_name, destination_name in [
        ("run_configuration_r.tsv", f"{sample_id}_run_configuration_r.tsv"),
        ("software_versions_r.tsv", f"{sample_id}_software_versions_r.tsv"),
    ]:
        source = next(
            (Path(w) / source_name for w in worker_dirs if (Path(w) / source_name).is_file()),
            None,
        )
        if source is not None:
            shutil.copy2(source, final_dir / destination_name)
            log.info(
                "[STAGE] stage=worker_metadata status=copied source=%s destination=%s",
                source,
                final_dir / destination_name,
            )

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
    for generated in worker_outdir.glob(f"{args.sample_id}_*"):
        if generated.is_file():
            generated.unlink()

    return run_susie(
        locus_file=locus_chunk,                     # 🔑 CHUNKED LOCUS FILE
        sumstat_file=args.susie_input_file,
        sample_id=args.sample_id,
        ld_ref=args.finemap_ld_ref,
        plink=args.plink,
        output_folder=str(worker_outdir),
        lp_threshold=args.lp_threshold,
        L=args.L,
        workers=1,                                  # 🔑 IMPORTANT
        min_ram_per_worker_gb=args.min_ram_per_worker_gb,
        timeout_ld_seconds=args.timeout_ld_seconds,
        timeout_susie_seconds=args.timeout_susie_seconds,
        skip_mhc=args.finemap_skip_mhc,
        finemap_mhc_chrom=args.finemap_mhc_chrom,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
        genome_build=getattr(args, "genome_build", DEFAULT_GENOME_BUILD),
    )



def _run_parallel_susie(args, progress):
    """Implement parallel SuSiE with documented stages and resource controls."""
    runtime_defaults = get_finemap_defaults()
    pipeline_total = DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL
    output_folder = Path(args.outdir).resolve()
    progress.record(
        "pipeline", "initialization", 1, pipeline_total, "completed",
        f"output_dir={output_folder}",
    )
    final_plink_path = None

    logger.info("[STAGE] stage=dependency_validation status=started")
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
        args.plink_version = (plink_version.stdout or plink_version.stderr).strip().splitlines()[0]
    except Exception:
        args.plink_version = "unknown"
    logger.info(
        "[STAGE] stage=dependency_validation status=completed plink=%s version=%s",
        args.plink,
        args.plink_version,
    )
    progress.record(
        "pipeline", "dependency_validation", 2, pipeline_total, "completed",
        f"PLINK={args.plink_version}",
    )

    args.min_ram_per_worker_gb = max(
        float(runtime_defaults["susie_ram_per_worker_gb"]),
        float(getattr(args, "min_ram_per_worker_gb", 0) or 0),
    )
    args.nthreads = safe_thread_count(
        requested_threads=args.nthreads,
        gb_per_thread=args.min_ram_per_worker_gb
    )

    if args.finemap_include_mhc:
        args.finemap_skip_mhc = False

    genome_build = getattr(args, "genome_build", DEFAULT_GENOME_BUILD)
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
            modified_name = output_folder / "locus_windows.tsv"
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
                modified_name = output_folder / "locus_ranges_extended.tsv"
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

    locus_chunk_dir = output_folder / "locus_chunks"
    locus_chunks = split_locus_file(
        locus_file=args.locus_file,
        outdir=locus_chunk_dir,
        n_chunks=args.nthreads,
        sep="\t",
        finemap_skip_mhc=args.finemap_skip_mhc,
        mhc_chr=args.finemap_mhc_chrom,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
        log=logger,
    )
    if not locus_chunks:
        raise RuntimeError("No locus chunks were generated — SuSiE was not run.")
    write_run_configuration(
        output_folder / "run_configuration.json",
        {
            "engine": "SuSiE",
            "inputs": {
                "locus_file": args.locus_file,
                "summary_statistics": args.susie_input_file,
                "ld_reference": args.finemap_ld_ref,
                "output_directory": output_folder,
                "sample_id": args.sample_id,
            },
            "analysis_parameters": {
                "genome_build": genome_build,
                "locus_type": args.locus_type,
                "window_kb": args.window_kb,
                "lp_threshold": args.lp_threshold,
                "L": args.L,
                "credible_set_coverage": runtime_defaults["credible_set_coverage"],
                "skip_mhc": args.finemap_skip_mhc,
                "mhc_chromosome": args.finemap_mhc_chrom,
                "mhc_start": args.finemap_mhc_start,
                "mhc_end": args.finemap_mhc_end,
            },
            "resource_parameters": {
                "selected_workers": len(locus_chunks),
                "min_ram_per_worker_gb": args.min_ram_per_worker_gb,
                "timeout_ld_seconds": args.timeout_ld_seconds,
                "timeout_susie_seconds": args.timeout_susie_seconds,
            },
            "software": {"plink_path": args.plink, "plink_version": args.plink_version},
            "defaults": runtime_defaults,
        },
    )
    progress.record(
        "pipeline", "locus_split", 4, pipeline_total, "completed",
        f"chunks={len(locus_chunks)}",
    )

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
                i + 1,
                args,
                output_folder,
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

    current_worker_dirs = [
        output_folder / f"worker_{index + 1}" for index in range(len(locus_chunks))
    ]
    merge_workers(
        worker_dirs=current_worker_dirs,
        final_dir=output_folder,
        sample_id=args.sample_id,
        log=logger,
    )

    qc_file = output_folder / f"{args.sample_id}_SuSiE_QC_summary.tsv"
    combined_file = output_folder / f"{args.sample_id}_SUSIE_combined_results.csv.gz"
    if not qc_file.is_file() or not combined_file.is_file():
        raise RuntimeError(
            "SuSiE workers completed but required combined outputs are missing. "
            f"Expected {qc_file.name} and {combined_file.name}."
        )

    qc_df = pd.read_csv(qc_file, sep="\t")
    n_success = int(qc_df.get("converged", pd.Series(dtype=bool)).astype(str).str.lower().eq("true").sum())
    n_failed = int(len(qc_df) - n_success)
    if n_success == 0:
        raise RuntimeError("All SuSiE loci failed or were skipped; no usable result was produced.")

    # FLAMES can consume these definitions with its -l option.  Keep only
    # loci that yielded a converged SuSiE result.
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
    used_loci = used_loci[used_loci["GenomicLocus"].astype(str).isin(successful_loci)]
    used_loci = used_loci.rename(columns={"CHR": "chr", "START": "start", "END": "end"})
    used_loci[["GenomicLocus", "chr", "start", "end"]].drop_duplicates().to_csv(
        output_folder / "flames_input" / "genomic_loci.tsv", sep="\t", index=False
    )

    version_file = output_folder / f"{args.sample_id}_software_versions.tsv"
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
        f"successful={n_success}; failed={n_failed}",
    )
    progress.record(
        "pipeline", "complete", pipeline_total, pipeline_total, "success",
        f"output_dir={output_folder}; successful={n_success}; failed={n_failed}",
    )
    logger.info(
        "SuSiE fine-mapping completed: output=%s successful=%d failed=%d",
        output_folder,
        n_success,
        n_failed,
    )
    flames_input=f'{output_folder}/flames_input/'
    return {
        "status": "success",
        "output_dir": str(output_folder),
        "flames_input": str(flames_input),
        "n_attempted": int(len(qc_df)),
        "n_successful": n_success,
        "n_failed": n_failed,
    }


def run_parallel_susie(args):
    """Coordinate SuSiE and retain existing return keys for downstream callers."""
    output_folder = Path(args.outdir).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)
    setup_logging(output_folder / "pipeline_summary.log")
    progress = ProgressRecorder(output_folder / "pipeline_progress.tsv", logger=logger)
    try:
        result = _run_parallel_susie(args, progress)
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
            DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL, "failed", str(exc)
        )
        logger.exception("SuSiE pipeline failed: %s", exc)
        raise
