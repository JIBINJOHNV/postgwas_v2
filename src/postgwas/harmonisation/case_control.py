import polars as pl
from typing import Tuple, Dict
from pathlib import Path
import io


def harmonize_sample_sizes(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Harmonize sample size columns (ncase, ncontrol) and compute Neff.
    All console output is suppressed — logs written ONLY to:

        logs/{gwas_outputname}_{chromosome}_harmonize_sample_sizes.log
    """

    # -------------------------------------------------------
    # Create log file
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_{chromosome}_harmonize_sample_sizes.log"

    # -------------------------------------------------------
    # Logging buffer (NO stdout printing)
    # -------------------------------------------------------
    log_buffer = io.StringIO()

    def log_print(*args, **kwargs):
        """Write ONLY to logfile, never to stdout."""
        log_buffer.write(" ".join(str(a) for a in args) + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC — uses log_print()
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print(f"\n📊 [Chr {chromosome}] Starting sample size harmonization...")

    # Extract mapping
    ncontrol_col = sample_column_dict.get("ncontrol_col", "NA")
    ncase_col    = sample_column_dict.get("ncase_col", "NA")
    ncontrol_val = sample_column_dict.get("ncontrol", "NA")
    ncase_val    = sample_column_dict.get("ncase", "NA")

    # ============================
    # CONTROLS
    # ============================
    if ncontrol_col != "NA" and ncontrol_col in df.columns:
        log_print(f"🔍 Using dataframe column for controls: '{ncontrol_col}' → Int64")
        df = df.with_columns(pl.col(ncontrol_col).cast(pl.Int64, strict=False))
        qc_info["ncontrol_source"] = f"column:{ncontrol_col}"

    elif ncontrol_val != "NA":
        log_print(f"ℹ️ Using provided constant for controls: {ncontrol_val} → df['ncontrol']")
        df = df.with_columns(pl.lit(int(ncontrol_val)).alias("ncontrol"))
        ncontrol_col = "ncontrol"
        sample_column_dict["ncontrol_col"] = "ncontrol"
        qc_info["ncontrol_source"] = f"fixed_value:{ncontrol_val}"

    else:
        log_print("⚠️ No control column/value found — skipping.")
        qc_info["ncontrol_source"] = "missing"

    # ============================
    # CASES
    # ============================
    if ncase_col != "NA" and ncase_col in df.columns:
        log_print(f"🔍 Using dataframe column for cases: '{ncase_col}' → Int64")
        df = df.with_columns(pl.col(ncase_col).cast(pl.Int64, strict=False))
        qc_info["ncase_source"] = f"column:{ncase_col}"

    elif ncase_val != "NA":
        log_print(f"ℹ️ Using provided constant for cases: {ncase_val} → df['ncase']")
        df = df.with_columns(pl.lit(int(ncase_val)).alias("ncase"))
        ncase_col = "ncase"
        sample_column_dict["ncase_col"] = "ncase"
        qc_info["ncase_source"] = f"fixed_value:{ncase_val}"

    else:
        log_print("⚠️ No case column/value found — skipping.")
        qc_info["ncase_source"] = "missing"

    # ============================
    # COMPUTE Neff
    # ============================
    has_controls = (ncontrol_col != "NA" and ncontrol_col in df.columns)
    has_cases    = (ncase_col != "NA" and ncase_col in df.columns)

    if has_controls and has_cases:
        try:
            log_print("🧮 Calculating Neff = 4 / (1/N_CASE + 1/N_CONTROL)")
            df = df.with_columns(
                (
                    4 / (
                        1 / pl.col(ncase_col).cast(pl.Float64) +
                        1 / pl.col(ncontrol_col).cast(pl.Float64)
                    )
                )
                .round(0)
                .cast(pl.Int64, strict=False)
                .alias("Neff")
            )
            qc_info["Neff_status"] = "calculated_from_case_control"

        except Exception as e:
            log_print(f"⚠️ Neff calculation failed: {e}")
            qc_info["Neff_status"] = "skipped_error"

    elif has_controls:
        log_print("ℹ️ Only controls available — Neff := ncontrol.")
        df = df.with_columns(pl.col(ncontrol_col).cast(pl.Int64).alias("Neff"))
        qc_info["Neff_status"] = f"fallback_ncontrol_only:{ncontrol_col}"

    elif has_cases:
        log_print("ℹ️ Only cases available — Neff := ncase.")
        df = df.with_columns(pl.col(ncase_col).cast(pl.Int64).alias("Neff"))
        qc_info["Neff_status"] = f"fallback_ncase_only:{ncase_col}"

    else:
        log_print("⚠️ Neither cases nor controls available — Neff not computed.")
        qc_info["Neff_status"] = "not_computed_no_inputs"

    sample_column_dict["neff_col"] = "Neff"

    # ============================
    # SUMMARY STATS
    # ============================
    cols_to_check = ["ncase_col", "ncontrol_col", "neff_col"]
    sample_stats = {}

    for key in cols_to_check:
        real_col = sample_column_dict.get(key)
        if real_col is None or real_col == "NA" or real_col not in df.columns:
            sample_stats[key] = "missing"
            continue

        stats_dict = (
            df.select([
                pl.col(real_col).min().alias("min"),
                pl.col(real_col).max().alias("max"),
                pl.col(real_col).mean().alias("mean"),
                pl.col(real_col).median().alias("median"),
                pl.len().alias("total")
            ])
            .to_dicts()[0]
        )

        sample_stats[real_col] = stats_dict

    qc_info["sample_size_stats"] = sample_stats
    qc_info["final_variants"] = df.height

    log_print("✅ Sample size harmonization completed.\n")

    # -------------------------------------------------------
    # Write log file
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, qc_info, sample_column_dict
