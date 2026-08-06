import sys,math,shutil,argparse,os
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path, safe_thread_count
import pandas as pd
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed


# SuSiE backend (Python) ---------------------
from postgwas.finemap.susie.main import (
    validate_locus_file,
    run_susie,
)


def merge_workers(worker_dirs, final_dir, sample_id):
    final_dir = Path(final_dir)
    for sub in [
        "plots", "flames_input", "locus_files",
        "rds_files", "logs", "ld_matrix_related"
    ]:
        (final_dir / sub).mkdir(parents=True, exist_ok=True)
    for w in worker_dirs:
        w = Path(w)
        #print(f"\t\t\t\tMerging {w}")
        for sub in [
            "plots", "flames_input", "locus_files",
            "rds_files", "logs", "ld_matrix_related"
        ]:
            src = w / sub
            dst = final_dir / sub

            if not src.exists():
                continue
            for f in src.iterdir():
                out = dst / f.name
                if not out.exists():
                    shutil.copy2(f, out)
    # --------------------------------------------------
    # 2. Linux streaming merge for per-worker tables
    # --------------------------------------------------
    final_dir_str = str(final_dir)
    base = final_dir_str
    cmds = [
        f"mkdir -p {base}",
        f"awk 'NR==1 || FNR>1' {base}/worker_*/{sample_id}_SuSiE_QC_summary.tsv > {base}/{sample_id}_SuSiE_QC_summary.tsv",
        f"awk 'NR==1 || FNR>1' {base}/worker_*/{sample_id}_SuSiE_failed_loci.tsv > {base}/{sample_id}_SuSiE_failed_loci.tsv",
        f"awk 'NR==1 || FNR>1' {base}/worker_*/{sample_id}_SUSIE_combined_credibleset.csv > {base}/{sample_id}_SUSIE_combined_credibleset.csv",
        f"zcat {base}/worker_*/{sample_id}_SUSIE_combined_results.csv.gz | awk 'NR==1 || FNR>1' | gzip > {base}/{sample_id}_SUSIE_combined_results.csv.gz",
        f"cat {base}/worker_*/{sample_id}_run_susie.log > {base}/{sample_id}_run_susie.ALL.log",
    ]

    print(f"\t\t\t\tConcatining all batch analysis results")
    for c in cmds:
        try:
            os.system(f"{c} 2>/dev/null")
        except Exception as e:
            print(f"\t\t\t\t❌ Merging command failed: {c}\nError: {e}")
            #raise

def split_locus_file(
    locus_file: str,
    outdir: Path,
    n_chunks: int,
    sep: str = "\t",
    finemap_skip_mhc: bool = False,
    mhc_chr: str = "6",
    mhc_start: int = 25000000,
    mhc_end: int = 35000000,   # ✅ FIXED (35 Mb)
):
    """
    Split locus file into n_chunks row-wise.
    If finemap_skip_mhc=True, remove MHC loci BEFORE splitting.
    Ensures every chunk contains at least one locus.
    Returns list of chunk locus file paths.
    """
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

        print(f"\t\t\t [INFO] MHC filtering enabled: {before} → {after} loci")

    # --------------------------------------------------
    # Adjust number of chunks to avoid empty files
    # --------------------------------------------------
    n = len(df)
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

    print(f"\t\t\tSplit locus file into {len(chunk_files)} chunks")

    return chunk_files






def run_susie_worker(
    locus_chunk: str,
    worker_id: int,
    args,
    base_output_dir: Path,
):
    worker_outdir = base_output_dir / f"worker_{worker_id}"
    worker_outdir.mkdir(parents=True, exist_ok=True)

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
    )



def run_parallel_susie(args):

    final_plink_path = None

    # 1. Check if user provided a specific path
    if args.plink:
        if Path(args.plink).exists():
            final_plink_path = args.plink
        else:
            # Path provided but invalid → ignore and fall back to PATH search
            final_plink_path = None

    # 2. If we don't have a valid path yet, search PATH
    if not final_plink_path:
        detected_plink = shutil.which("plink") or shutil.which("plink2")
        if detected_plink:
            final_plink_path = detected_plink

    # 3. Final Validation: If still missing, CRASH.
    if not final_plink_path:
        print("\n❌ [ERROR] PLINK executable not found!")
        print("   The pipeline failed to find PLINK in the provided argument or the system $PATH.")
        print("   Please either:")
        print("   1. Add 'plink' to your system environment, OR")
        print("   2. Provide a valid full path using --plink /path/to/plink")
        sys.exit(1)

    # FIX: Update args.plink
    args.plink = final_plink_path
    # threads based on RAM
    args.min_ram_per_worker_gb = max(14, int(getattr(args, "min_ram_per_worker_gb", 0) or 0))
    args.nthreads = safe_thread_count(
        requested_threads=args.nthreads,
        gb_per_thread=args.min_ram_per_worker_gb
    )

    if args.finemap_include_mhc:
        args.finemap_skip_mhc = False

    validate_locus_file(args.locus_file)

    output_folder = Path(args.outdir)
    output_folder.mkdir(parents=True, exist_ok=True)


    # 2. Data Prep
    try:
        # sep=None with engine='python' is great for auto-detecting delimiters
        loci_df = pd.read_csv(args.locus_file, sep=None, engine="python")
    except Exception as e:
        sys.exit(f"ERROR: Could not read locus file: {e}")

    # Clean column names (strip whitespace and uppercase)
    loci_df.columns = [x.strip().upper() for x in loci_df.columns]
    flank_bp = int(args.window_kb * 1000)

    if args.locus_type == "point":
        if 'POS' in loci_df.columns and "CHR" in loci_df.columns:
            # 1. Ensure POS is numeric to avoid calculation errors
            loci_df['POS'] = pd.to_numeric(loci_df['POS'])
            
            # 2. Use .clip(lower=0) to ensure coordinates aren't negative
            # 3. Cast to int to ensure clean TSV output (no .0 decimals in coordinates)
            loci_df['START'] = (loci_df['POS'] - flank_bp).clip(lower=0).astype(int)
            loci_df['END'] = (loci_df['POS'] + flank_bp).astype(int)
            
            print(f"[*] Created ±{args.window_kb}kb windows around POS column.")
            
            # 4. Fix string formatting: use f-string with quotes correctly
            # We use os.path.splitext to avoid "file.csv.tsv" names
            base_path = os.path.splitext(args.locus_file)[0]
            modified_name = f"{base_path}_windows.tsv"
            
            loci_df.to_csv(modified_name, sep="\t", index=False)
            args.locus_file = modified_name
        else:
            sys.exit(f"ERROR: locus_type is 'point', but 'POS' or 'CHR' columns are missing. Found: {list(loci_df.columns)}")
    elif args.locus_type == "range":
        required = ["CHR", "START", "END"]
        if all(col in loci_df.columns for col in required):
            # Ensure START/END are integers even in range mode
            loci_df['START'] = pd.to_numeric(loci_df['START']).astype(int)
            loci_df['END'] = pd.to_numeric(loci_df['END']).astype(int)
            print("[*] Using exact CHR, START, END coordinates from file.")
        else:
            sys.exit(f"ERROR: locus_type is 'range', but 'CHR', 'START', or 'END' columns are missing. Found: {list(loci_df.columns)}")
    print(f"The locus file used in the analysis is {args.locus_file}")


    # ---- split locus file
    locus_chunk_dir = output_folder / "locus_chunks"
    locus_chunks = split_locus_file(
        locus_file=args.locus_file,
        outdir=locus_chunk_dir,
        n_chunks=args.nthreads,
        sep="\t",
        finemap_skip_mhc=args.finemap_skip_mhc,
        mhc_chr =args.finemap_mhc_chrom,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
    )
    # ---- multiprocessing (SAFE)
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

        for fut in as_completed(futures):
            chunk = futures[fut]
            try:
                res = fut.result()
                results.append(res)
            except Exception as e:
                print(f"❌ Worker failed for {chunk}: {e}")
                raise

    merge_workers(
        worker_dirs=sorted(output_folder.glob("worker_*")),
        final_dir=output_folder,
        sample_id=args.sample_id
    )

    print(f"                👉 Results are saved in: {output_folder.resolve()}\n")
    print(f"                🎉 Finemapping using SuSiE is Completed\n")
    flames_input=f'{output_folder}/flames_input/'
    return {
        "status": "success",
        "output_dir": str(output_folder),
        "flames_input": str(flames_input),
    }
