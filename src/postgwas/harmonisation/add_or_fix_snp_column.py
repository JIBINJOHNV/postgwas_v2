
import polars as pl
from pathlib import Path
import io
from typing import Tuple, Dict

import polars as pl
import io
from pathlib import Path
from typing import Tuple, Dict


def add_or_fix_snp_column(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict
) -> Tuple[pl.DataFrame, dict]:
    """
    Ensure the dataframe has a valid SNP column:
      - If missing or set to 'NA', create it using chr_pos_ea_oa format.
      - Replace ':' with '_'.
      - Fill missing SNP IDs using chr_pos_ea_oa.
      - Multiprocessing-safe buffered logging (no print to screen).
    """

    # -------------------------------------------------------
    # Setup per-chromosome log file
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_snp.log"

    # Buffered logger (silent)
    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")   # <-- no printing to screen

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    log_print(f"📊 [Chr {chromosome}] Starting SNP column harmonization...")

    # Extract mappings
    chr_col = sample_column_dict.get("chr_col")
    pos_col = sample_column_dict.get("pos_col")
    ea_col  = sample_column_dict.get("ea_col")
    oa_col  = sample_column_dict.get("oa_col")
    snp_col = sample_column_dict.get("snp_id_col", "NA")

    # -------------------------------------------------------
    # Create SNP column if missing
    # -------------------------------------------------------
    if snp_col == "NA" or snp_col not in df.columns:
        log_print("🧬 SNP column not found — creating from chr, pos, ea, oa...")

        snp_col = "snp_col"

        df = df.with_columns(
            (
                pl.col(chr_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(pos_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(ea_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(oa_col).cast(pl.Utf8)
            ).alias(snp_col)
        )

        sample_column_dict["snp_id_col"] = snp_col

    else:
        log_print(f"✅ SNP column '{snp_col}' found — cleaning and checking missing values.")

    # -------------------------------------------------------
    # Clean SNP: replace ":" with "_"
    # -------------------------------------------------------
    df = df.with_columns(pl.col(snp_col).str.replace_all(":", "_"))

    # -------------------------------------------------------
    # Detect missing SNP IDs
    # -------------------------------------------------------
    missing_mask = df.select(pl.col(snp_col).is_null() | (pl.col(snp_col) == "")).to_series()
    n_missing = int(missing_mask.sum())

    log_print(f"🔍 Found {n_missing} missing SNP IDs.")

    # -------------------------------------------------------
    # Fill missing SNP IDs
    # -------------------------------------------------------
    if n_missing > 0:
        log_print("🧩 Filling missing SNP IDs from chr, pos, ea, oa columns...")

        filled_col = (
            pl.when(pl.col(snp_col).is_null() | (pl.col(snp_col) == ""))
            .then(
                pl.col(chr_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(pos_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(ea_col).cast(pl.Utf8)
                + pl.lit("_")
                + pl.col(oa_col).cast(pl.Utf8)
            )
            .otherwise(pl.col(snp_col))
        )

        df = df.with_columns(filled_col.alias(snp_col))
        log_print("✅ Missing SNP IDs filled successfully.")

    log_print(f"🎯 SNP column harmonization complete — total rows: {df.height:,}\n")

    # -------------------------------------------------------
    # WRITE LOGFILE
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, sample_column_dict
