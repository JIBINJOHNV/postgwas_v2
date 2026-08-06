"""Reliable MAGMA gene and gene-set analysis helpers.

The module runs MAGMA, harmonizes summary-statistic SNP identifiers to the LD
reference, performs multiple-testing correction, and annotates gene-set results
with overlapping genes.
"""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from postgwas.utils.main import decide_magma_batches_from_annot
from statsmodels.stats.multitest import multipletests


class MagmaPipelineError(RuntimeError):
    """Error containing the exact MAGMA pipeline stage and sample."""

    def __init__(
        self,
        stage: str,
        sample_id: str,
        message: str,
        *,
        command: Sequence[str] | None = None,
        exit_code: int | None = None,
        stderr: str | None = None,
    ) -> None:
        details = [
            "MAGMA pipeline failed",
            f"Sample: {sample_id}",
            f"Stage: {stage}",
            f"Reason: {message}",
        ]
        if command:
            details.append(f"Command: {shlex.join([str(value) for value in command])}")
        if exit_code is not None:
            details.append(f"Exit code: {exit_code}")
        if stderr:
            details.append(f"MAGMA stderr: {stderr.strip()}")
        super().__init__("\n".join(details))
        self.stage = stage
        self.sample_id = sample_id


def setup_logger(log_file: str, name: str = "magma_logger") -> logging.Logger:
    """Create an isolated file logger for one resolved log path."""
    log_path = Path(log_file).expanduser()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"{name}:{log_path.resolve()}")
    logger.setLevel(logging.INFO)

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(file_handler)
    logger.propagate = False
    return logger


def _validate_sample_id(sample_id: str) -> None:
    if not sample_id or Path(sample_id).name != sample_id or "\x00" in sample_id:
        raise ValueError(
            "sample_id must be a non-empty filename component without path separators"
        )


def _require_nonempty_file(path: str | Path, label: str) -> Path:
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {file_path}")
    if file_path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {file_path}")
    return file_path


def validate_gene_identifier_overlap(
    gene_loc_file: str | Path,
    geneset_file: str | Path,
    minimum_overlap: float = 0.50,
) -> dict:
    """Confirm that gene-location and GMT files use compatible gene IDs.

    The first column of the MAGMA gene-location file is compared with the gene
    columns (third column onward) of the GMT file.  Complete equality is not
    expected because the files can cover different gene universes, so the
    overlap is calculated relative to the smaller unique-ID universe.
    """
    if not 0 <= minimum_overlap <= 1:
        raise ValueError("minimum gene-ID overlap must be between 0 and 1")

    gene_loc_path = _require_nonempty_file(gene_loc_file, "Gene-location file")
    gmt_path = _require_nonempty_file(geneset_file, "Gene-set file")

    gene_loc_ids: set[str] = set()
    with gene_loc_path.open("r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            gene_id = stripped.split()[0]
            if gene_id.upper() in {"GENE", "GENE_ID"}:
                continue
            gene_loc_ids.add(gene_id)

    if not gene_loc_ids:
        raise ValueError(
            f"No gene IDs were found in the first column of {gene_loc_path}"
        )

    geneset_ids: set[str] = set()
    geneset_count = 0
    usable_geneset_count = 0
    with gmt_path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"Malformed GMT line {line_number}: expected gene-set name, "
                    "description, and at least one gene ID"
                )
            set_gene_ids = {value.strip() for value in parts[2:] if value.strip()}
            if not set_gene_ids:
                raise ValueError(f"GMT line {line_number} contains no gene IDs")
            geneset_count += 1
            geneset_ids.update(set_gene_ids)
            if not gene_loc_ids.isdisjoint(set_gene_ids):
                usable_geneset_count += 1

    if not geneset_ids:
        raise ValueError(f"No gene IDs were found in GMT columns 3 onward: {gmt_path}")

    overlapping_ids = gene_loc_ids.intersection(geneset_ids)
    comparison_size = min(len(gene_loc_ids), len(geneset_ids))
    overlap_fraction = len(overlapping_ids) / comparison_size

    result = {
        "gene_location_unique_ids": len(gene_loc_ids),
        "geneset_unique_ids": len(geneset_ids),
        "comparison_unique_ids": comparison_size,
        "overlapping_unique_ids": len(overlapping_ids),
        "overlap_fraction_of_smaller_universe": overlap_fraction,
        "minimum_required_overlap": minimum_overlap,
        "total_gene_sets": geneset_count,
        "gene_sets_with_at_least_one_location_gene": usable_geneset_count,
    }

    if not overlapping_ids or overlap_fraction < minimum_overlap:
        gene_loc_examples = sorted(gene_loc_ids)[:5]
        geneset_examples = sorted(geneset_ids)[:5]
        raise ValueError(
            "Gene identifiers in the gene-location and GMT files are "
            "incompatible: "
            f"{len(overlapping_ids):,}/{comparison_size:,} unique IDs overlap "
            f"({overlap_fraction:.2%}); at least {minimum_overlap:.2%} is required. "
            "The gene-location first column and GMT columns 3 onward must use "
            "the same identifier system (for example, Entrez Gene IDs). "
            f"Gene-location examples: {gene_loc_examples}. "
            f"GMT examples: {geneset_examples}."
        )

    return result


def _resolve_magma_executable(magma: str) -> str:
    resolved = shutil.which(magma)
    if resolved is None:
        raise FileNotFoundError(
            f"MAGMA executable was not found or is not executable: {magma}"
        )
    return resolved


def _validate_ld_reference(ld_ref: str) -> None:
    missing = [
        f"{ld_ref}{suffix}"
        for suffix in (".bed", ".bim", ".fam")
        if not Path(f"{ld_ref}{suffix}").is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "LD reference is incomplete. Missing: " + ", ".join(missing)
        )


def _verify_outputs(outputs: Sequence[str | Path], stage: str, sample_id: str) -> None:
    for output in outputs:
        try:
            _require_nonempty_file(output, f"Expected output from {stage}")
        except (FileNotFoundError, ValueError) as error:
            raise MagmaPipelineError(stage, sample_id, str(error)) from error


def run_subprocess_with_logging(
    cmd: Sequence[str],
    logger: logging.Logger,
    *,
    stage: str,
    sample_id: str,
    expected_outputs: Sequence[str | Path] = (),
) -> subprocess.CompletedProcess[str]:
    """Run one command and preserve its exact stage, stdout, and stderr."""
    command = [str(value) for value in cmd]
    logger.info("[%s] Command: %s", stage, shlex.join(command))
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        logger.exception("[%s] Executable or input was not found", stage)
        raise MagmaPipelineError(
            stage,
            sample_id,
            str(error),
            command=command,
        ) from error
    except subprocess.CalledProcessError as error:
        logger.error("[%s] Command failed with exit code %s", stage, error.returncode)
        if error.stdout:
            logger.error("[%s] stdout:\n%s", stage, error.stdout.rstrip())
        if error.stderr:
            logger.error("[%s] stderr:\n%s", stage, error.stderr.rstrip())
        raise MagmaPipelineError(
            stage,
            sample_id,
            "MAGMA returned a non-zero exit code",
            command=command,
            exit_code=error.returncode,
            stderr=error.stderr,
        ) from error

    if result.stdout:
        logger.info("[%s] stdout:\n%s", stage, result.stdout.rstrip())
    if result.stderr:
        logger.warning("[%s] stderr:\n%s", stage, result.stderr.rstrip())
    _verify_outputs(expected_outputs, stage, sample_id)
    logger.info("[%s] Completed successfully", stage)
    return result


def _read_whitespace_table(path: str | Path, label: str) -> pd.DataFrame:
    file_path = _require_nonempty_file(path, label)
    try:
        return pd.read_csv(file_path, sep=r"\s+", engine="c")
    except Exception as error:
        raise ValueError(f"Failed to parse {label} {file_path}: {error}") from error


def _normalize_chromosome(values: pd.Series) -> pd.Series:
    normalized = (
        values.astype(str)
        .str.strip()
        .str.replace(r"^chr", "", case=False, regex=True)
        .str.upper()
    )
    return normalized.replace({"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"})


def _validate_p_values(values: pd.Series, source: str | Path) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    invalid = numeric.isna() | ~np.isfinite(numeric) | (numeric < 0) | (numeric > 1)
    if invalid.any():
        examples = values[invalid].head(5).tolist()
        raise ValueError(
            f"P column in {source} contains {int(invalid.sum())} invalid values; "
            f"examples: {examples}. Values must be finite numbers in [0, 1]."
        )
    return numeric


def harmonize_magma_inputs(
    pval_file: str,
    snp_loc_file: str,
    ld_ref: str,
    output_dir: str | Path,
    sample_id: str,
    n_sample_col: str,
    logger: logging.Logger,
    minimum_snp_overlap: float = 0.05,
) -> dict:
    """Map input SNP IDs to LD-reference IDs and retain the lowest-P duplicate.

    Direct SNP-ID matches are preferred. Remaining variants are mapped only when
    chromosome and position identify exactly one variant in the LD-reference BIM.
    Ambiguous positions are never guessed.
    """
    stage = "input_harmonization"
    if not 0 <= minimum_snp_overlap <= 1:
        raise ValueError("minimum_snp_overlap must be between 0 and 1")

    try:
        pval_df = _read_whitespace_table(pval_file, "MAGMA p-value file")
        snp_loc_df = _read_whitespace_table(snp_loc_file, "MAGMA SNP-location file")
        required_pval = {"SNP", "P", n_sample_col}
        required_loc = {"SNP", "CHR", "BP"}
        missing_pval = sorted(required_pval.difference(pval_df.columns))
        missing_loc = sorted(required_loc.difference(snp_loc_df.columns))
        if missing_pval:
            raise ValueError(
                f"P-value file is missing required columns: {', '.join(missing_pval)}"
            )
        if missing_loc:
            raise ValueError(
                f"SNP-location file is missing required columns: {', '.join(missing_loc)}"
            )
        if pval_df.empty:
            raise ValueError("P-value file contains no variants")

        original_columns = pval_df.columns.tolist()
        pval_df["SNP"] = pval_df["SNP"].astype(str).str.strip()
        pval_df["_P_NUMERIC"] = _validate_p_values(pval_df["P"], pval_file)
        sample_sizes = pd.to_numeric(pval_df[n_sample_col], errors="coerce")
        invalid_n = (
            sample_sizes.isna() | ~np.isfinite(sample_sizes) | (sample_sizes <= 0)
        )
        if invalid_n.any():
            raise ValueError(
                f"Sample-size column {n_sample_col!r} contains "
                f"{int(invalid_n.sum())} missing, non-numeric, or non-positive values"
            )
        pval_df["_ROW_ORDER"] = np.arange(len(pval_df))

        snp_loc_df = snp_loc_df[["SNP", "CHR", "BP"]].copy()
        snp_loc_df["SNP"] = snp_loc_df["SNP"].astype(str).str.strip()
        snp_loc_df["CHR_NORM"] = _normalize_chromosome(snp_loc_df["CHR"])
        snp_loc_df["BP_NORM"] = pd.to_numeric(snp_loc_df["BP"], errors="coerce")
        if snp_loc_df["BP_NORM"].isna().any():
            raise ValueError("SNP-location BP column contains non-numeric values")

        conflicting_locations = (
            snp_loc_df.groupby("SNP")[["CHR_NORM", "BP_NORM"]]
            .nunique()
            .max(axis=1)
            .gt(1)
        )
        if conflicting_locations.any():
            examples = conflicting_locations[conflicting_locations].index[:5].tolist()
            raise ValueError(
                "SNP-location file assigns conflicting positions to the same SNP; "
                f"examples: {examples}"
            )
        snp_loc_df = snp_loc_df.drop_duplicates("SNP", keep="first")

        bim_file = _require_nonempty_file(f"{ld_ref}.bim", "LD-reference BIM file")
        reference_df = pd.read_csv(
            bim_file,
            sep=r"\s+",
            engine="c",
            header=None,
            usecols=[0, 1, 3],
            dtype={0: str, 1: str},
        )
        reference_df.columns = ["REF_CHR", "REF_SNP", "REF_BP"]
        reference_df["REF_SNP"] = reference_df["REF_SNP"].astype(str).str.strip()
        reference_df["CHR_NORM"] = _normalize_chromosome(reference_df["REF_CHR"])
        reference_df["BP_NORM"] = pd.to_numeric(reference_df["REF_BP"], errors="coerce")
        if reference_df["BP_NORM"].isna().any():
            raise ValueError(
                "LD-reference BIM contains non-numeric base-pair positions"
            )
        duplicated_reference_ids = reference_df["REF_SNP"].duplicated(keep=False)
        if duplicated_reference_ids.any():
            examples = (
                reference_df.loc[duplicated_reference_ids, "REF_SNP"].head(5).tolist()
            )
            raise ValueError(
                f"LD-reference BIM contains duplicate SNP IDs; examples: {examples}"
            )

        working = pval_df.merge(
            snp_loc_df[["SNP", "CHR_NORM", "BP_NORM"]],
            on="SNP",
            how="left",
            validate="many_to_one",
        )
        reference_by_id = reference_df.rename(
            columns={
                "REF_SNP": "SNP",
                "CHR_NORM": "REF_CHR_NORM",
                "BP_NORM": "REF_BP_NORM",
            }
        )[["SNP", "REF_CHR_NORM", "REF_BP_NORM"]]
        working = working.merge(
            reference_by_id, on="SNP", how="left", validate="many_to_one"
        )
        direct_match = working["REF_BP_NORM"].notna()

        position_mismatch = (
            direct_match
            & working["BP_NORM"].notna()
            & (
                (working["CHR_NORM"] != working["REF_CHR_NORM"])
                | (working["BP_NORM"] != working["REF_BP_NORM"])
            )
        )
        if position_mismatch.any():
            examples = working.loc[position_mismatch, "SNP"].head(5).tolist()
            raise ValueError(
                "Directly matched SNP IDs have different coordinates in the input and "
                f"LD reference, suggesting a genome-build mismatch; examples: {examples}"
            )

        unique_coordinate_reference = reference_df[
            ~reference_df.duplicated(["CHR_NORM", "BP_NORM"], keep=False)
        ][["CHR_NORM", "BP_NORM", "REF_SNP"]].rename(
            columns={"REF_SNP": "COORDINATE_REF_SNP"}
        )
        working = working.merge(
            unique_coordinate_reference,
            on=["CHR_NORM", "BP_NORM"],
            how="left",
            validate="many_to_one",
        )
        working["CANONICAL_SNP"] = working["SNP"].where(
            direct_match, working["COORDINATE_REF_SNP"]
        )

        input_unique = int(working["SNP"].nunique())
        resolved_unique = int(
            working.loc[working["CANONICAL_SNP"].notna(), "SNP"].nunique()
        )
        overlap_fraction = resolved_unique / input_unique if input_unique else 0.0
        if resolved_unique == 0 or overlap_fraction < minimum_snp_overlap:
            raise ValueError(
                f"Only {resolved_unique:,}/{input_unique:,} unique input SNP IDs "
                f"({overlap_fraction:.2%}) could be matched to the LD reference; "
                f"required minimum is {minimum_snp_overlap:.2%}"
            )

        resolved = working[working["CANONICAL_SNP"].notna()].copy()
        resolved["SNP"] = resolved["CANONICAL_SNP"]
        resolved = resolved.sort_values(
            ["SNP", "_P_NUMERIC", "_ROW_ORDER"], kind="mergesort"
        )
        before_deduplication = len(resolved)
        resolved = resolved.drop_duplicates("SNP", keep="first")
        duplicates_removed = before_deduplication - len(resolved)
        resolved = resolved.sort_values("_ROW_ORDER", kind="mergesort")

        harmonized_pval = resolved[original_columns].copy()
        harmonized_pval["P"] = resolved["_P_NUMERIC"].to_numpy()
        reference_lookup = reference_df.set_index("REF_SNP")
        selected_reference = reference_lookup.loc[harmonized_pval["SNP"]]
        harmonized_snp_loc = pd.DataFrame(
            {
                "SNP": harmonized_pval["SNP"].to_numpy(),
                "CHR": selected_reference["REF_CHR"].to_numpy(),
                "BP": selected_reference["REF_BP"].to_numpy(),
            }
        )

        output_folder = Path(output_dir)
        output_folder.mkdir(parents=True, exist_ok=True)
        harmonized_pval_path = output_folder / f"{sample_id}_magma_pval_harmonized.tsv"
        harmonized_snp_loc_path = (
            output_folder / f"{sample_id}_magma_snp_loc_harmonized.tsv"
        )
        harmonized_pval.to_csv(harmonized_pval_path, sep="\t", index=False)
        harmonized_snp_loc.to_csv(harmonized_snp_loc_path, sep="\t", index=False)

        qc = {
            "input_rows": len(pval_df),
            "input_unique_snps": input_unique,
            "direct_id_matches": int(direct_match.sum()),
            "coordinate_matches": int(
                (~direct_match & working["CANONICAL_SNP"].notna()).sum()
            ),
            "unresolved_rows": int(working["CANONICAL_SNP"].isna().sum()),
            "overlap_fraction": overlap_fraction,
            "canonical_duplicates_removed": int(duplicates_removed),
            "retained_rows": len(harmonized_pval),
        }
        logger.info("[%s] QC: %s", stage, qc)
        logger.info(
            "[%s] Retained the lowest P value whenever multiple input records "
            "mapped to the same LD-reference SNP ID",
            stage,
        )
        return {
            "pval_file": str(harmonized_pval_path),
            "snp_loc_file": str(harmonized_snp_loc_path),
            "qc": qc,
        }
    except MagmaPipelineError:
        raise
    except Exception as error:
        logger.exception("[%s] Input harmonization failed", stage)
        raise MagmaPipelineError(stage, sample_id, str(error)) from error


def run_magma_analysis(
    magma_analysis_folder: str,
    sample_id: str,
    ld_ref: str,
    gene_loc_file: str,
    snp_loc_file: str,
    pval_file: str,
    log_file: str,
    window_upstream: int = 35,
    window_downstream: int = 10,
    gene_model: str = "snp-wise=mean",
    n_sample_col: str = "N_COL",
    num_cores: int = 6,
    seed: int = 10,
    magma: str = "magma",
    geneset_file: str | None = None,
    harmonize_snp_ids: bool = True,
    minimum_snp_overlap: float = 0.05,
    minimum_gene_id_overlap: float = 0.50,
) -> dict:
    """Run annotation, gene analysis, batch merge, and optional gene-set analysis."""
    _validate_sample_id(sample_id)
    if num_cores < 1:
        raise ValueError("num_cores must be at least 1")
    if window_upstream < 0 or window_downstream < 0:
        raise ValueError("MAGMA annotation windows cannot be negative")

    folder = Path(magma_analysis_folder)
    folder.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(log_file)

    gene_id_validation = None
    try:
        magma_executable = _resolve_magma_executable(magma)
        _validate_ld_reference(ld_ref)
        _require_nonempty_file(gene_loc_file, "Gene-location file")
        _require_nonempty_file(snp_loc_file, "SNP-location file")
        _require_nonempty_file(pval_file, "P-value file")
        if geneset_file is not None:
            _require_nonempty_file(geneset_file, "Gene-set file")
            gene_id_validation = validate_gene_identifier_overlap(
                gene_loc_file=gene_loc_file,
                geneset_file=geneset_file,
                minimum_overlap=minimum_gene_id_overlap,
            )
            logger.info(
                "[input_validation] Gene-ID compatibility: %s",
                gene_id_validation,
            )
            print(
                "                 [Preflight] Gene-ID overlap: "
                f"{gene_id_validation['overlapping_unique_ids']:,}/"
                f"{gene_id_validation['comparison_unique_ids']:,} "
                f"({gene_id_validation['overlap_fraction_of_smaller_universe']:.2%}); "
                "usable gene sets: "
                f"{gene_id_validation['gene_sets_with_at_least_one_location_gene']:,}/"
                f"{gene_id_validation['total_gene_sets']:,}"
            )
    except Exception as error:
        logger.exception("[input_validation] Validation failed")
        raise MagmaPipelineError("input_validation", sample_id, str(error)) from error

    harmonization_result = None
    if harmonize_snp_ids:
        harmonization_result = harmonize_magma_inputs(
            pval_file=pval_file,
            snp_loc_file=snp_loc_file,
            ld_ref=ld_ref,
            output_dir=folder,
            sample_id=sample_id,
            n_sample_col=n_sample_col,
            logger=logger,
            minimum_snp_overlap=minimum_snp_overlap,
        )
        pval_file = harmonization_result["pval_file"]
        snp_loc_file = harmonization_result["snp_loc_file"]

    gene_prefix = folder / f"{sample_id}_magma"
    merged_prefix = folder / (
        f"{sample_id}_magma_{window_upstream}up_{window_downstream}down"
    )
    annot_path = Path(f"{gene_prefix}.genes.annot")
    merged_raw = Path(f"{merged_prefix}.genes.raw")
    merged_out = Path(f"{merged_prefix}.genes.out")

    logger.info("[annotation] Running MAGMA annotation")
    print("                 [Step 1/4] Running MAGMA annotation...")
    annotate_cmd = [
        magma_executable,
        "--annotate",
        f"window={window_upstream},{window_downstream}",
        "--snp-loc",
        str(snp_loc_file),
        "--gene-loc",
        str(gene_loc_file),
        "--out",
        str(gene_prefix),
    ]
    run_subprocess_with_logging(
        annotate_cmd,
        logger,
        stage="annotation",
        sample_id=sample_id,
        expected_outputs=[annot_path],
    )

    try:
        with annot_path.open("r") as handle:
            annotation_lines = sum(1 for _ in handle)
        num_batches = (
            decide_magma_batches_from_annot(annot_file=str(annot_path))
            if annotation_lines > 1000
            else 1
        )
        num_batches = int(num_batches)
        if num_batches < 1:
            raise ValueError(f"Invalid number of MAGMA batches: {num_batches}")
    except Exception as error:
        logger.exception("[batch_planning] Failed to determine batch count")
        raise MagmaPipelineError("batch_planning", sample_id, str(error)) from error

    mode = "single run" if num_batches == 1 else f"{num_batches} batches"
    logger.info("[gene_analysis] Running MAGMA gene analysis (%s)", mode)
    print(f"                 [Step 2/4] Running MAGMA gene analysis ({mode})...")

    def run_batch(batch_number: int | None = None) -> None:
        stage = (
            "gene_analysis"
            if batch_number is None
            else f"gene_analysis_batch_{batch_number}_of_{num_batches}"
        )
        batch_cmd = [
            magma_executable,
            "--bfile",
            str(ld_ref),
            "--gene-annot",
            str(annot_path),
            "--pval",
            str(pval_file),
            f"ncol={n_sample_col}",
            "duplicate=error",
            "--gene-model",
            gene_model,
            "--out",
            str(gene_prefix),
        ]
        if num_batches > 1:
            batch_cmd.extend(["--batch", str(batch_number), str(num_batches)])
        expected = (
            [Path(f"{gene_prefix}.genes.raw"), Path(f"{gene_prefix}.genes.out")]
            if num_batches == 1
            else []
        )
        run_subprocess_with_logging(
            batch_cmd,
            logger,
            stage=stage,
            sample_id=sample_id,
            expected_outputs=expected,
        )

    if num_batches == 1:
        run_batch()
        try:
            Path(f"{gene_prefix}.genes.raw").replace(merged_raw)
            Path(f"{gene_prefix}.genes.out").replace(merged_out)
            _verify_outputs([merged_raw, merged_out], "gene_output_move", sample_id)
        except MagmaPipelineError:
            raise
        except Exception as error:
            logger.exception("[gene_output_move] Failed to move gene outputs")
            raise MagmaPipelineError(
                "gene_output_move", sample_id, str(error)
            ) from error
    else:
        with ThreadPoolExecutor(max_workers=min(num_cores, num_batches)) as executor:
            futures = {
                executor.submit(run_batch, batch): batch
                for batch in range(1, num_batches + 1)
            }
            for future in as_completed(futures):
                future.result()

        logger.info("[batch_merge] Merging MAGMA batch results")
        print("                 [Step 3/4] Merging MAGMA batch results...")
        merge_cmd = [
            magma_executable,
            "--merge",
            str(gene_prefix),
            "--out",
            str(merged_prefix),
        ]
        run_subprocess_with_logging(
            merge_cmd,
            logger,
            stage="batch_merge",
            sample_id=sample_id,
            expected_outputs=[merged_raw, merged_out],
        )
        for intermediate in folder.glob(f"{sample_id}_magma.batch*"):
            try:
                intermediate.unlink()
            except OSError as error:
                logger.warning(
                    "Could not delete intermediate %s: %s", intermediate, error
                )

    if num_batches == 1:
        logger.info("[batch_merge] Single run; merge was not required")
        print("                 [Step 3/4] Single run; merge step skipped.")

    gsa_out: Path | None = None
    if geneset_file is not None:
        logger.info("[gene_set_analysis] Running MAGMA gene-set analysis")
        print("                 [Step 4/4] Running MAGMA gene-set analysis...")
        gene_set_prefix = folder / sample_id
        gsa_out = Path(f"{gene_set_prefix}.gsa.out")
        geneset_cmd = [
            magma_executable,
            "--gene-results",
            str(merged_raw),
            "--set-annot",
            str(geneset_file),
            "--out",
            str(gene_set_prefix),
            "--seed",
            str(seed),
        ]
        run_subprocess_with_logging(
            geneset_cmd,
            logger,
            stage="gene_set_analysis",
            sample_id=sample_id,
            expected_outputs=[gsa_out],
        )
    else:
        logger.info("[gene_set_analysis] No gene-set file supplied; stage skipped")

    logger.info("MAGMA analysis completed successfully")
    return {
        # Preserve the original prefix-valued key for existing callers.
        "gene_annot": str(gene_prefix),
        "gene_annot_file": str(annot_path),
        "merged_prefix": str(merged_prefix),
        "genes_raw": str(merged_raw),
        "genes_out": str(merged_out),
        "gsa_out": str(gsa_out) if gsa_out else None,
        "snp_harmonization": harmonization_result,
        "gene_id_validation": gene_id_validation,
    }


def correct_p_values(
    magma_analysis_folder: str,
    sample_id: str,
    log_file: str,
    output_file: str | None = None,
    gsa_out_file: str | None = None,
) -> pl.DataFrame:
    """Validate MAGMA gene-set p-values and add global/category corrections."""
    stage = "p_value_correction"
    logger = setup_logger(log_file)
    folder = Path(magma_analysis_folder)
    output_path = (
        Path(output_file)
        if output_file is not None
        else folder / f"{sample_id}.gsa_corrected.tsv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gsa_out = (
        Path(gsa_out_file)
        if gsa_out_file is not None
        else folder / f"{sample_id}.gsa.out"
    )

    try:
        _require_nonempty_file(gsa_out, "MAGMA gene-set output")
        # MAGMA writes a whitespace-delimited table.  FULL_NAME is the final,
        # variable-width column and can be much longer than the names present
        # in the first few hundred rows.  read_fwf(infer_nrows=...) therefore
        # truncates later names to the inferred width and can create false
        # duplicates.  Gene-set identifiers cannot contain whitespace, so a
        # whitespace-delimited parser preserves every FULL_NAME exactly.
        pdf = pd.read_csv(gsa_out, sep=r"\s+", comment="#", engine="c")
        pdf = pdf.loc[:, [column for column in pdf.columns if str(column).strip()]]
        required = {"VARIABLE", "P"}
        missing = sorted(required.difference(pdf.columns))
        if missing:
            raise ValueError(
                f"MAGMA gene-set output is missing columns {missing}; found {pdf.columns.tolist()}"
            )
        if pdf.empty:
            raise ValueError("MAGMA gene-set output contains no result rows")
        pdf["P"] = _validate_p_values(pdf["P"], gsa_out)
        if "FULL_NAME" not in pdf.columns:
            pdf = pdf.rename(columns={"VARIABLE": "FULL_NAME"})
        else:
            pdf["FULL_NAME"] = pdf["FULL_NAME"].fillna(pdf["VARIABLE"])
        pdf["FULL_NAME"] = pdf["FULL_NAME"].astype(str).str.strip()
        duplicates = pdf["FULL_NAME"].duplicated(keep=False)
        if duplicates.any():
            examples = pdf.loc[duplicates, "FULL_NAME"].head(5).tolist()
            raise ValueError(f"Duplicate gene-set names found; examples: {examples}")

        # Construct from Python lists so pyarrow is not an undeclared requirement.
        df = pl.DataFrame(pdf.to_dict(orient="list"))
        number_of_tests = df.height
        p_values = df["P"].to_numpy()
        df = df.with_columns(
            [
                (pl.col("P") * number_of_tests)
                .clip(upper_bound=1.0)
                .alias("P_bonferroni_corr"),
                (1 - (1 - pl.col("P")).pow(number_of_tests)).alias("P_sidak_corr"),
                pl.Series("P_holm_corr", multipletests(p_values, method="holm")[1]),
                pl.Series("P_fdr_bh_corr", multipletests(p_values, method="fdr_bh")[1]),
            ]
        )

        categories = {
            "go_kegg_reactome": r"^(GOBP_|GOCC_|GOMF_|KEGG_|REACTOME_)",
            "go": r"^(GOBP_|GOCC_|GOMF_)",
            "reactome": r"^REACTOME_",
            "kegg": r"^KEGG_",
        }
        for prefix, pattern in categories.items():
            subset = df.filter(pl.col("FULL_NAME").str.contains(pattern))
            bonferroni_column = f"{prefix}_P_bonferroni_corr"
            fdr_column = f"{prefix}_P_fdr_bh_corr"
            if subset.height:
                subset_p = subset["P"].to_numpy()
                corrections = subset.select("FULL_NAME").with_columns(
                    [
                        pl.Series(
                            bonferroni_column,
                            np.minimum(subset_p * len(subset_p), 1.0),
                        ),
                        pl.Series(
                            fdr_column,
                            multipletests(subset_p, method="fdr_bh")[1],
                        ),
                    ]
                )
                df = df.join(corrections, on="FULL_NAME", how="left")
            else:
                df = df.with_columns(
                    [
                        pl.lit(None).alias(bonferroni_column),
                        pl.lit(None).alias(fdr_column),
                    ]
                )

        df = df.sort("P")
        df.write_csv(output_path, separator="\t")
        logger.info("[%s] Results saved to %s", stage, output_path)
        return df
    except MagmaPipelineError:
        raise
    except Exception as error:
        logger.exception("[%s] Failed", stage)
        raise MagmaPipelineError(stage, sample_id, str(error)) from error


def _normalize_gene_id(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    normalized = str(value).strip()
    return normalized or None


def parse_gmt_file(gmt_path: Path) -> pl.DataFrame:
    """Read GMT content, validate unique set names, and de-duplicate genes."""
    _require_nonempty_file(gmt_path, "GMT gene-set file")
    data = []
    seen_names: set[str] = set()
    with gmt_path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3 or not parts[0].strip():
                raise ValueError(
                    f"Malformed GMT line {line_number}: expected name, description, and genes"
                )
            name = parts[0].strip()
            if name in seen_names:
                raise ValueError(f"Duplicate gene-set name in GMT: {name}")
            seen_names.add(name)
            genes = list(
                dict.fromkeys(
                    gene
                    for gene in (_normalize_gene_id(value) for value in parts[2:])
                    if gene is not None
                )
            )
            if not genes:
                raise ValueError(f"GMT gene set {name!r} has no valid genes")
            data.append({"V1": name, "V2": parts[1].strip(), "V3": ",".join(genes)})
    if not data:
        raise ValueError(f"GMT file contains no gene sets: {gmt_path}")
    return pl.DataFrame(data)


def add_gene_details(gene_string: str, gene_p_values: dict[str, float]) -> dict:
    """Match one normalized GMT gene list to normalized MAGMA gene IDs."""
    genes = list(
        dict.fromkeys(
            gene
            for gene in (_normalize_gene_id(value) for value in gene_string.split(","))
            if gene is not None
        )
    )
    common = sorted(
        ((gene, gene_p_values[gene]) for gene in genes if gene in gene_p_values),
        key=lambda item: item[0],
    )
    return {
        "CommonGenes": ",".join(gene for gene, _ in common) or None,
        "CommonGenes_Pvalues": ",".join(str(p_value) for _, p_value in common) or None,
        "TotalGenes": len(genes),
        "CommonGeneCount": len(common),
    }


def process_gene_sets(
    magma_analysis_folder: str,
    sample_id: str,
    geneset_file: str,
    corrected_df: pl.DataFrame,
    log_file: str,
    genes_out_file: str | None = None,
    window_upstream: int = 35,
    window_downstream: int = 10,
    return_metadata: bool = False,
) -> pl.DataFrame | dict:
    """Add MAGMA gene overlaps to corrected gene-set results."""
    stage = "gene_set_annotation"
    logger = setup_logger(log_file)
    folder = Path(magma_analysis_folder)
    genes_out_path = (
        Path(genes_out_file)
        if genes_out_file is not None
        else folder
        / f"{sample_id}_magma_{window_upstream}up_{window_downstream}down.genes.out"
    )
    output_file = folder / f"{sample_id}.gsa_corrected_annotated.tsv"

    try:
        _require_nonempty_file(genes_out_path, "MAGMA gene output")
        genes_pdf = pd.read_fwf(genes_out_path, comment="#", infer_nrows=500)
        missing = sorted({"GENE", "P"}.difference(genes_pdf.columns))
        if missing:
            raise ValueError(
                f"MAGMA gene output is missing columns {missing}; found {genes_pdf.columns.tolist()}"
            )
        if genes_pdf.empty:
            raise ValueError("MAGMA gene output contains no genes")
        genes_pdf["P"] = _validate_p_values(genes_pdf["P"], genes_out_path)
        genes_pdf["GENE"] = genes_pdf["GENE"].map(_normalize_gene_id)
        genes_pdf = genes_pdf.dropna(subset=["GENE"])
        if genes_pdf["GENE"].duplicated().any():
            raise ValueError("MAGMA gene output contains duplicate normalized gene IDs")
        gene_p_values = dict(zip(genes_pdf["GENE"], genes_pdf["P"]))

        gene_sets = parse_gmt_file(Path(geneset_file))
        results = [
            add_gene_details(gene_string, gene_p_values)
            for gene_string in gene_sets["V3"].to_list()
        ]
        gene_sets_with_overlap = pl.concat(
            [gene_sets, pl.DataFrame(results)], how="horizontal"
        )
        merged_df = (
            corrected_df.join(
                gene_sets_with_overlap,
                left_on="FULL_NAME",
                right_on="V1",
                how="left",
            )
            .sort("P")
            .with_columns(pl.lit(sample_id).alias("sample_id"))
        )
        merged_df.write_csv(output_file, separator="\t")
        logger.info("[%s] Results saved to %s", stage, output_file)
        if return_metadata:
            return {"data": merged_df, "output_file": str(output_file)}
        return merged_df
    except MagmaPipelineError:
        raise
    except Exception as error:
        logger.exception("[%s] Failed", stage)
        raise MagmaPipelineError(stage, sample_id, str(error)) from error


def magma_analysis_pipeline(
    output_dir: str,
    sample_id: str,
    ld_ref: str,
    gene_loc_file: str,
    snp_loc_file: str,
    pval_file: str,
    geneset_file: str | None,
    log_file: str,
    threads: int = 6,
    window_upstream: int = 35,
    window_downstream: int = 10,
    gene_model: str = "snp-wise=mean",
    n_sample_col: str = "N_COL",
    seed: int = 10,
    magma: str = "magma",
    harmonize_snp_ids: bool = True,
    minimum_snp_overlap: float = 0.05,
    minimum_gene_id_overlap: float = 0.50,
) -> dict:
    """Run the complete MAGMA pipeline and return actual output paths."""
    try:
        magma_result = run_magma_analysis(
            magma_analysis_folder=output_dir,
            sample_id=sample_id,
            ld_ref=ld_ref,
            gene_loc_file=gene_loc_file,
            snp_loc_file=snp_loc_file,
            pval_file=pval_file,
            geneset_file=geneset_file,
            log_file=log_file,
            n_sample_col=n_sample_col,
            num_cores=threads,
            window_upstream=window_upstream,
            window_downstream=window_downstream,
            gene_model=gene_model,
            seed=seed,
            magma=magma,
            harmonize_snp_ids=harmonize_snp_ids,
            minimum_snp_overlap=minimum_snp_overlap,
            minimum_gene_id_overlap=minimum_gene_id_overlap,
        )

        result = {
            "magma_genes_raw": magma_result["genes_raw"],
            "magma_genes_out": magma_result["genes_out"],
            "magma_genes_prefix": magma_result["merged_prefix"],
            "snp_harmonization": magma_result["snp_harmonization"],
            "gene_id_validation": magma_result["gene_id_validation"],
        }
        if geneset_file is None:
            return result

        corrected_df = correct_p_values(
            magma_analysis_folder=output_dir,
            sample_id=sample_id,
            log_file=log_file,
            gsa_out_file=magma_result["gsa_out"],
        )
        gene_set_result = process_gene_sets(
            magma_analysis_folder=output_dir,
            sample_id=sample_id,
            geneset_file=geneset_file,
            corrected_df=corrected_df,
            log_file=log_file,
            genes_out_file=magma_result["genes_out"],
            window_upstream=window_upstream,
            window_downstream=window_downstream,
            return_metadata=True,
        )
        result["magma_pathway"] = gene_set_result["output_file"]
        return result
    except MagmaPipelineError:
        raise
    except Exception as error:
        logger = setup_logger(log_file)
        logger.exception("[pipeline] Unexpected failure")
        raise MagmaPipelineError("pipeline", sample_id, str(error)) from error
