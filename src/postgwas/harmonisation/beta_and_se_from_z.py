from postgwas.harmonisation import print_qc_summary
from postgwas.harmonisation.detect_pvalue_type_process import log_print
from postgwas.harmonisation import export_gwas_sumstat_forgwas2vcf
import polars as pl
from typing import Tuple, Dict
from pathlib import Path
import io


def calculate_beta_and_se_from_z(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Compute BETA and SE from Z-score using robust GWAS-safe implementation.

    Features:
        • Multiprocessing-safe logging (no stdout)
        • Partitioned QC (like Beta/SE filtering block)
        • Handles EAF edge cases (0/1, NaN)
        • Tracks removed variants by reason
    """

    # -------------------------------------------------------
    # Setup
    # -------------------------------------------------------
    log_file = None
    sample_column_dict = sample_column_dict.copy()  # SAFE

    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")

    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_beta_se_from_z.log"

    invalid_dir = Path(output_dir) / "invalid_zscore_to_beta_se"
    invalid_dir.mkdir(parents=True, exist_ok=True)

    invalid_file = invalid_dir / f"{gwas_outputname}_chr{chromosome}.tsv"

    log_buffer = io.StringIO()
    screen_buffer = []

    def log_print(*args, indent: int = 0, screen: bool = True):
        # Use 2 spaces per indent level for cleaner terminal viewing
        prefix = "  " * indent 
        msg = f"{prefix}{' '.join(str(a) for a in args)}"
        log_buffer.write(msg + "\n")
        if screen:
            screen_buffer.append(msg)

    qc_info = {"initial_variants": df.height}

    try:
        # -------------------------------------------------------
        # Column mappings
        # -------------------------------------------------------
        imp_z_col = sample_column_dict.get("imp_z_col", "NA")
        beta_col = sample_column_dict.get("beta_or_col", "NA")
        se_col = sample_column_dict.get("se_col", "NA")
        eaf_col = sample_column_dict.get("eaf_col", "NA")

        log_print(
            f"📊 [Chr {chromosome}] Starting calculation of BETA/SE from Z-score column: {imp_z_col}...",
            indent=0,
            screen=True
        ) 

        # -------------------------------------------------------
        # Validate required columns
        # -------------------------------------------------------
        if beta_col == "NA" and imp_z_col == "NA":
            log_print(
                f"Beta column is not present, but Z-score column also missing, so not able to do harmonisation",
                indent=4,
                screen=True
            )
            qc_info["status"] = "skipped_no_zscore_or_effectsize_column"
            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            raise KeyError("Please provide sumstat file with Z-score or BETA/OR columns. Both columns are missing.") 
            
        # -------------------------------------------------------
        # Validate required columns
        # -------------------------------------------------------
        if beta_col != "NA" or se_col != "NA":
            if imp_z_col != "NA":
                log_print(
                    f"Beta and SE columns are present, so skipping the calculation of BETA and SE from Z-score column {imp_z_col}.",
                    indent=4,
                    screen=True
                )
            else:
                log_print(
                    f"Beta and SE columns are present, so skipping the calculation of BETA and SE from Z-score.",
                    indent=4,
                    screen=True
                )
                log_print(
                    f"In this file, the Z-score column is not present, it will be calculated in the next step from BETA and SE columns.\n",
                    indent=4,
                    screen=True
                )
            qc_info["status"] = "Beta and SE columns are present"

            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            return df, qc_info, sample_column_dict

        if beta_col != "NA" and se_col == "NA":
            log_print(
                f"Beta column is present, but SE column is missing, so skipping the calculation of BETA from Z-score column {imp_z_col}. SE will be calculated in the next step using `calculate_se_from_beta_pvalue` function.",
                indent=4,
                screen=True
            )
            qc_info["status"] = "Beta column is present, but SE column is missing"

            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            return df, qc_info, sample_column_dict

        if imp_z_col == "NA" or imp_z_col not in df.columns:
            log_print(
                "⚠️ Beta OR SE columns are not present, but Z-score column is missing, so not able to do harmonisation",
                indent=4,
                screen=True
            )
            qc_info["status"] = "skipped_no_z"

            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            raise KeyError("Please provide sumstat file with Z-score or BETA/OR and SE columns.")

        if eaf_col == "NA" or eaf_col not in df.columns:
            msg = f"❌ Missing EAF column: {eaf_col}"
            log_print(msg, indent=4, screen=True)

            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            raise KeyError(msg)

        if "Neff" not in df.columns:
            msg = "❌ Missing Neff column"
            log_print(msg, indent=4, screen=True)

            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())

            if screen_buffer:
                print("\n".join(screen_buffer), flush=True)

            raise KeyError(msg)

        # -------------------------------------------------------
        # Cast numeric
        # -------------------------------------------------------
        df = df.with_columns([
            pl.col(imp_z_col).cast(pl.Float64, strict=False),
            pl.col(eaf_col).cast(pl.Float64, strict=False),
            pl.col("Neff").cast(pl.Float64, strict=False),
        ])

        # -------------------------------------------------------
        # FILTERING (PARTITIONED QC — LIKE YOUR STYLE)
        # -------------------------------------------------------
        cond_z_null = pl.col(imp_z_col).is_null() | ~pl.col(imp_z_col).is_finite()
        cond_eaf_null = pl.col(eaf_col).is_null() | ~pl.col(eaf_col).is_finite()

        cond_eaf_lt_zero = pl.col(eaf_col) < 0
        cond_eaf_gt_one = pl.col(eaf_col) > 1

        cond_eaf_invalid = (
            cond_eaf_lt_zero |
            cond_eaf_gt_one
        )

        cond_neff_null = pl.col("Neff").is_null() | ~pl.col("Neff").is_finite()
        cond_neff_invalid = pl.col("Neff") <= 0

        # Partition
        df_z_null = df.filter(cond_z_null)

        df_eaf_null = df.filter(
            ~cond_z_null &
            cond_eaf_null
        )

        df_eaf_lt_zero = df.filter(
            ~cond_z_null &
            ~cond_eaf_null &
            cond_eaf_lt_zero
        )

        df_eaf_gt_one = df.filter(
            ~cond_z_null &
            ~cond_eaf_null &
            ~cond_eaf_lt_zero &
            cond_eaf_gt_one
        )

        df_eaf_invalid = pl.concat([
            df_eaf_lt_zero,
            df_eaf_gt_one
        ]).unique()

        df_neff_null = df.filter(
            ~cond_z_null &
            ~cond_eaf_null &
            ~cond_eaf_invalid &
            cond_neff_null
        )

        df_neff_invalid = df.filter(
            ~cond_z_null &
            ~cond_eaf_null &
            ~cond_eaf_invalid &
            ~cond_neff_null &
            cond_neff_invalid
        )

        # Add removal reason
        df_z_null = df_z_null.with_columns(pl.lit("null_z").alias("remove_reason"))
        df_eaf_null = df_eaf_null.with_columns(pl.lit("null_eaf").alias("remove_reason"))
        df_eaf_lt_zero = df_eaf_lt_zero.with_columns(pl.lit("eaf_lt_zero").alias("remove_reason"))
        df_eaf_gt_one = df_eaf_gt_one.with_columns(pl.lit("eaf_gt_one").alias("remove_reason"))
        df_neff_null = df_neff_null.with_columns(pl.lit("null_neff").alias("remove_reason"))
        df_neff_invalid = df_neff_invalid.with_columns(pl.lit("invalid_neff").alias("remove_reason"))

        # Combine removed
        df_removed = pl.concat([
            df_z_null,
            df_eaf_null,
            df_eaf_lt_zero,
            df_eaf_gt_one,
            df_neff_null,
            df_neff_invalid
        ], how="diagonal").unique()

        # Save invalid variants
        if df_removed.height > 0:
            df_removed.write_csv(invalid_file, separator="\t")

        # Valid condition
        condition_valid = (
            (~cond_z_null) &
            (~cond_eaf_null) &
            (~cond_eaf_invalid) &
            (~cond_neff_null) &
            (~cond_neff_invalid)
        )

        df = df.filter(condition_valid)

        # -------------------------------------------------------
        # COUNTS
        # -------------------------------------------------------
        n_initial = qc_info["initial_variants"]

        n_z_null = df_z_null.height
        n_eaf_null = df_eaf_null.height

        n_eaf_lt_zero = df_eaf_lt_zero.height
        n_eaf_gt_one = df_eaf_gt_one.height

        n_eaf_invalid = (
            n_eaf_lt_zero +
            n_eaf_gt_one
        )

        n_neff_null = df_neff_null.height
        n_neff_invalid = df_neff_invalid.height

        n_removed = n_z_null + n_eaf_null + n_eaf_invalid + n_neff_null + n_neff_invalid
        n_final = df.height

        pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

        qc_info.update({
            "removed_z_null": n_z_null,
            "removed_eaf_null": n_eaf_null,
            "removed_eaf_invalid": n_eaf_invalid,
            "removed_eaf_lt_zero": n_eaf_lt_zero,
            "removed_eaf_gt_one": n_eaf_gt_one,
            "removed_neff_null": n_neff_null,
            "removed_neff_invalid": n_neff_invalid,
            "variants_removed_total": n_removed,
            "variants_remaining": n_final,
            "invalid_variants_file": str(invalid_file) if df_removed.height > 0 else "NA"
        })

        # -------------------------------------------------------
        # COMPUTE DENOMINATOR (OPTIMIZED)
        # -------------------------------------------------------
        df = df.with_columns(
            (
                2 * pl.col(eaf_col) * (1 - pl.col(eaf_col)) *
                (pl.col("Neff") + pl.col(imp_z_col) ** 2)
            ).alias("_den")
        )

        # -------------------------------------------------------
        # COMPUTE BETA
        # -------------------------------------------------------
        if beta_col == "NA" or beta_col not in df.columns:
            df = df.with_columns(
                (pl.col(imp_z_col) / pl.col("_den").sqrt()).alias("BETA")
            )
            log_print(
                f"✅ BETA column not present in the file provided, so calculating BETA using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
                indent=4,
                screen=True
            )
            log_print(
                "⚠ If the effect allele frequency in the provided file is not from the ZScore calculated samples then beta calculation may not be accurate",
                indent=4,
                screen=True
            )
            sample_column_dict["beta_or_col"] = "BETA"
            beta_col = "BETA"
            qc_info["beta_computed"] = True
        else:
            log_print(
                f"✅ BETA column is present in the file provided, so not calculating BETA using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
                indent=4,
                screen=True
            )
            qc_info["beta_computed"] = False

        # -------------------------------------------------------
        # COMPUTE SE
        # -------------------------------------------------------
        if se_col == "NA" or se_col not in df.columns:
            df = df.with_columns(
                (1 / pl.col("_den").sqrt()).alias("SE")
            )
            log_print(
                f"✅ SE column not present in the file provided, so calculating SE using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
                indent=4,
                screen=True
            )
            sample_column_dict["se_col"] = "SE"
            se_col = "SE"
            qc_info["se_computed"] = True
        else:
            log_print(
                f"✅ SE column is present in the file provided, so not calculating SE using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
                indent=4,
                screen=True
            )
            qc_info["se_computed"] = False

        # -------------------------------------------------------
        # LOG SUMMARY
        # -------------------------------------------------------
        log_print(f"\n📊 [Chr {chromosome}] Z → Beta/SE QC Summary:", indent=4, screen=True)
        log_print(f"Initial: {n_initial}", indent=4, screen=True)
        log_print(f"Removed: {n_removed} ({pct_removed:.2f}%)", indent=4, screen=True)
        log_print(f"  - Null Z        : {n_z_null}", indent=4, screen=True)
        log_print(f"  - Null EAF      : {n_eaf_null}", indent=4, screen=True)
        log_print(f"  - Invalid EAF   : {n_eaf_invalid}", indent=4, screen=True)
        log_print(f"      - EAF < 0   : {n_eaf_lt_zero}", indent=4, screen=True)
        log_print(f"      - EAF > 1   : {n_eaf_gt_one}", indent=4, screen=True)
        log_print(f"  - Null Neff     : {n_neff_null}", indent=4, screen=True)
        log_print(f"  - Invalid Neff  : {n_neff_invalid}", indent=4, screen=True)
        log_print(f"Remaining: {n_final}", indent=4, screen=True)

        if df_removed.height > 0:
            log_print(f"Invalid variants saved: {invalid_file}", indent=4, screen=True)

        # Cleanup
        df = df.drop("_den")

        # -------------------------------------------------------
        # FINAL SUMMARY
        # -------------------------------------------------------
        summary = df.select([
            pl.col(beta_col).min().alias("beta_min"),
            pl.col(beta_col).max().alias("beta_max"),
            pl.col(beta_col).mean().alias("beta_mean"),
            pl.col(beta_col).std().alias("beta_std"),
            pl.col(se_col).min().alias("se_min"),
            pl.col(se_col).max().alias("se_max"),
            pl.col(se_col).mean().alias("se_mean"),
            pl.col(se_col).std().alias("se_std"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update(summary) 

        log_print("\n📈 BETA/SE Summary:", indent=4, screen=True) # Level 1 = 4 spaces
        
        # Use fixed-width formatting (:>10.6f) so decimals align perfectly
        beta_row = (f"BETA -> Min: {summary['beta_min']:>10.6f} | "
                    f"Max: {summary['beta_max']:>10.6f} | "
                    f"Mean: {summary['beta_mean']:>10.6f} | "
                    f"Std: {summary['beta_std']:>10.6f}")
                    
        se_row   = (f"SE   -> Min: {summary['se_min']:>10.6f} | "
                    f"Max: {summary['se_max']:>10.6f} | "
                    f"Mean: {summary['se_mean']:>10.6f} | "
                    f"Std: {summary['se_std']:>10.6f} | "
                    f"Total: {summary['total']}")

        log_print(beta_row, indent=4, screen=True) # Level 2 = 8 spaces
        log_print(se_row, indent=4, screen=True)

        # -------------------------------------------------------
        # SAVE LOG
        # -------------------------------------------------------
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        # print once at end
        if screen_buffer:
            print("\n".join(screen_buffer), flush=True)

        return df, qc_info, sample_column_dict

    except Exception:
        # write log + print buffered messages once on failure
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        if screen_buffer:
            print("\n".join(screen_buffer), flush=True)

        raise 
    
# def calculate_beta_and_se_from_z(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     """
#     Compute BETA and SE from Z-score using robust GWAS-safe implementation.

#     Features:
#         • Multiprocessing-safe logging (no stdout)
#         • Partitioned QC (like Beta/SE filtering block)
#         • Handles EAF edge cases (0/1, NaN)
#         • Tracks removed variants by reason
#     """

#     # -------------------------------------------------------
#     # Setup
#     # -------------------------------------------------------
#     sample_column_dict = sample_column_dict.copy()  # SAFE

#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir = sample_column_dict.get("output_folder", ".")

#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_beta_se_from_z.log"

#     log_buffer = io.StringIO()
#     screen_buffer = []   # <-- added (buffer screen messages)

#     def log_print(*args, indent: int = 0, screen: bool = True):
#         msg = ("    " * indent) + " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n")
#         if screen:
#             screen_buffer.append(msg)   # <-- changed (buffer, don't print now)

#     qc_info = {"initial_variants": df.height}

#     try:
#         # -------------------------------------------------------
#         # Column mappings
#         # -------------------------------------------------------
#         imp_z_col = sample_column_dict.get("imp_z_col", "NA")
#         beta_col = sample_column_dict.get("beta_or_col", "NA")
#         se_col = sample_column_dict.get("se_col", "NA")
#         eaf_col = sample_column_dict.get("eaf_col", "NA")

#         log_print(f"📊 [Chr {chromosome}] Starting BETA/SE from Z-score calculation...", indent=4, screen=True)
#         # -------------------------------------------------------
#         # Validate required columns
#         # -------------------------------------------------------
#         if beta_col != "NA" and se_col != "NA":
#             log_print(
#                 "Beta and SE columns are present, so skipping the calculation of BETA and SE from Z-score.",
#                 indent=4,
#                 screen=True
#             )
#             qc_info["status"] = "Beta and SE columns are present"

#             with open(log_file, "w") as f:
#                 f.write(log_buffer.getvalue())

#             if screen_buffer:
#                 print("\n".join(screen_buffer), flush=True)

#             return df, qc_info, sample_column_dict

#         if imp_z_col == "NA" or imp_z_col not in df.columns:
#             log_print(
#                 "⚠️ Beta OR SE columns are not present, but Z-score column is missing, so not able to do harmonisation",
#                 indent=4,
#                 screen=True
#             )
#             qc_info["status"] = "skipped_no_z"

#             with open(log_file, "w") as f:
#                 f.write(log_buffer.getvalue())

#             if screen_buffer:
#                 print("\n".join(screen_buffer), flush=True)

#             raise KeyError("Please provide sumstat file with Z-score or BETA/OR and SE columns.")

#         if eaf_col == "NA" or eaf_col not in df.columns:
#             msg = f"❌ Missing EAF column: {eaf_col}"
#             log_print(msg, indent=4, screen=True)

#             with open(log_file, "w") as f:
#                 f.write(log_buffer.getvalue())

#             if screen_buffer:
#                 print("\n".join(screen_buffer), flush=True)

#             raise KeyError(msg)

#         if "Neff" not in df.columns:
#             msg = "❌ Missing Neff column"
#             log_print(msg, indent=4, screen=True)

#             with open(log_file, "w") as f:
#                 f.write(log_buffer.getvalue())

#             if screen_buffer:
#                 print("\n".join(screen_buffer), flush=True)

#             raise KeyError(msg)

#         # -------------------------------------------------------
#         # Cast numeric
#         # -------------------------------------------------------
#         df = df.with_columns([
#             pl.col(imp_z_col).cast(pl.Float64, strict=False),
#             pl.col(eaf_col).cast(pl.Float64, strict=False),
#             pl.col("Neff").cast(pl.Float64, strict=False),
#         ])

#         # -------------------------------------------------------
#         # FILTERING (PARTITIONED QC — LIKE YOUR STYLE)
#         # -------------------------------------------------------
#         cond_z_null = pl.col(imp_z_col).is_null() | ~pl.col(imp_z_col).is_finite()
#         cond_eaf_null = pl.col(eaf_col).is_null() | ~pl.col(eaf_col).is_finite()
#         cond_eaf_invalid = (pl.col(eaf_col) <= 0) | (pl.col(eaf_col) >= 1)
#         cond_neff_null = pl.col("Neff").is_null() | ~pl.col("Neff").is_finite()
#         cond_neff_invalid = pl.col("Neff") <= 0

#         # Partition
#         df_z_null = df.filter(cond_z_null)
#         df_eaf_null = df.filter(~cond_z_null & cond_eaf_null)
#         df_eaf_invalid = df.filter(~cond_z_null & ~cond_eaf_null & cond_eaf_invalid)
#         df_neff_null = df.filter(~cond_z_null & ~cond_eaf_null & ~cond_eaf_invalid & cond_neff_null)
#         df_neff_invalid = df.filter(
#             ~cond_z_null &
#             ~cond_eaf_null &
#             ~cond_eaf_invalid &
#             ~cond_neff_null &
#             cond_neff_invalid
#         )

#         # Combine removed
#         df_removed = pl.concat([
#             df_z_null,
#             df_eaf_null,
#             df_eaf_invalid,
#             df_neff_null,
#             df_neff_invalid
#         ]).unique()

#         # Valid condition
#         condition_valid = (
#             (~cond_z_null) &
#             (~cond_eaf_null) &
#             (~cond_eaf_invalid) &
#             (~cond_neff_null) &
#             (~cond_neff_invalid)
#         )

#         df = df.filter(condition_valid)

#         # -------------------------------------------------------
#         # COUNTS
#         # -------------------------------------------------------
#         n_initial = qc_info["initial_variants"]

#         n_z_null = df_z_null.height
#         n_eaf_null = df_eaf_null.height
#         n_eaf_invalid = df_eaf_invalid.height
#         n_neff_null = df_neff_null.height
#         n_neff_invalid = df_neff_invalid.height

#         n_removed = n_z_null + n_eaf_null + n_eaf_invalid + n_neff_null + n_neff_invalid
#         n_final = df.height

#         pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#         # -------------------------------------------------------
#         # LOG SUMMARY
#         # -------------------------------------------------------
#         log_print(f"\n📊 [Chr {chromosome}] Z → Beta/SE QC Summary:", indent=6, screen=True)
#         log_print(f"Initial: {n_initial}", indent=6, screen=True)
#         log_print(f"Removed: {n_removed} ({pct_removed:.2f}%)", indent=6, screen=True)
#         log_print(f"  - Null Z        : {n_z_null}", indent=6, screen=True)
#         log_print(f"  - Null EAF      : {n_eaf_null}", indent=6, screen=True)
#         log_print(f"  - Invalid EAF   : {n_eaf_invalid}", indent=6, screen=True)
#         log_print(f"  - Null Neff     : {n_neff_null}", indent=6, screen=True)
#         log_print(f"  - Invalid Neff  : {n_neff_invalid}", indent=6, screen=True)
#         log_print(f"Remaining: {n_final}", indent=6, screen=True)

#         qc_info.update({
#             "removed_z_null": n_z_null,
#             "removed_eaf_null": n_eaf_null,
#             "removed_eaf_invalid": n_eaf_invalid,
#             "removed_neff_null": n_neff_null,
#             "removed_neff_invalid": n_neff_invalid,
#             "variants_removed_total": n_removed,
#             "variants_remaining": n_final
#         })

#         # -------------------------------------------------------
#         # COMPUTE DENOMINATOR (OPTIMIZED)
#         # -------------------------------------------------------
#         df = df.with_columns(
#             (
#                 2 * pl.col(eaf_col) * (1 - pl.col(eaf_col)) *
#                 (pl.col("Neff") + pl.col(imp_z_col) ** 2)
#             ).alias("_den")
#         )

#         # -------------------------------------------------------
#         # COMPUTE BETA
#         # -------------------------------------------------------
#         if beta_col == "NA" or beta_col not in df.columns:
#             df = df.with_columns(
#                 (pl.col(imp_z_col) / pl.col("_den").sqrt()).alias("BETA")
#             )
#             log_print(
#                 f"✅ BETA column not present in the file provided, so calculating BETA using Z-score column: {imp_z_col} and EAF column: {eaf_col}", 
#                 f"If the effect allele frequency in the provided file is not from the ZScore calcualted samples then beta calcualtion may not be accurate",
#                 indent=4,
#                 screen=True
#             )
#             sample_column_dict["beta_or_col"] = "BETA"
#             beta_col = "BETA"
#             qc_info["beta_computed"] = True
#         else:
#             log_print(
#                 f"✅ BETA column is present in the file provided, so not calculating BETA using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
#                 indent=4,
#                 screen=True
#             )
#             qc_info["beta_computed"] = False

#         # -------------------------------------------------------
#         # COMPUTE SE
#         # -------------------------------------------------------
#         if se_col == "NA" or se_col not in df.columns:
#             df = df.with_columns(
#                 (1 / pl.col("_den").sqrt()).alias("SE")
#             )
#             log_print(
#                 f"✅ SE column not present in the file provided, so calculating SE using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
#                 indent=4,
#                 screen=True
#             )
#             sample_column_dict["se_col"] = "SE"
#             se_col = "SE"
#             qc_info["se_computed"] = True
#         else:
#             log_print(
#                 f"✅ SE column is present in the file provided, so not calculating SE using Z-score column: {imp_z_col} and EAF column: {eaf_col}",
#                 indent=4,
#                 screen=True
#             )
#             qc_info["se_computed"] = False

#         # Cleanup
#         df = df.drop("_den")

#         # -------------------------------------------------------
#         # FINAL SUMMARY
#         # -------------------------------------------------------
#         summary = df.select([
#             pl.col(beta_col).min().alias("beta_min"),
#             pl.col(beta_col).max().alias("beta_max"),
#             pl.col(beta_col).mean().alias("beta_mean"),
#             pl.col(beta_col).std().alias("beta_std"),
#             pl.col(se_col).min().alias("se_min"),
#             pl.col(se_col).max().alias("se_max"),
#             pl.col(se_col).mean().alias("se_mean"),
#             pl.col(se_col).std().alias("se_std"),
#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         qc_info.update(summary)

#         log_print("\n📈 BETA/SE Summary:")
#         log_print(summary)

#         # -------------------------------------------------------
#         # SAVE LOG
#         # -------------------------------------------------------
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         # print once at end
#         if screen_buffer:
#             print("\n".join(screen_buffer), flush=True)

#         return df, qc_info, sample_column_dict

#     except Exception:
#         # write log + print buffered messages once on failure
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         if screen_buffer:
#             print("\n".join(screen_buffer), flush=True)

#         raise




# import polars as pl
# import numpy as np
# from typing import Tuple, Dict
# from pathlib import Path
# import io


# def calculate_beta_and_se_from_z(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     """
#     Compute BETA and SE from Z-score (imp_z_col) using zmaf & Neff formulas.
#     SAFE FOR MULTIPROCESSING:
#         • No console printing
#         • All messages logged to per-chromosome logfile
#     """

#     # -------------------------------------------------------
#     # Create logfile
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     eaf_col_name = sample_column_dict.get("eaf_col", "NA")

#     if eaf_col_name == "NA" or eaf_col_name not in df.columns:
#         msg = f"❌ CRITICAL ERROR: EAF column '{eaf_col_name}' is missing from the DataFrame. Check your input mappings."
#         log_print(msg)
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())
#         raise KeyError(msg)  # Force the process to stop 
 
#     if "Neff" not in df.columns:
#         raise KeyError("❌ CRITICAL ERROR: 'Neff' column is required for BETA/SE computation.")

#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_beta_se_from_z.log"

#     # -------------------------------------------------------
#     # Logging buffer (NO stdout output)
#     # -------------------------------------------------------
#     log_buffer = io.StringIO()

#     def log_print(*args, **kwargs):
#         """Write only to file buffer, never to stdout."""
#         log_buffer.write(" ".join(str(a) for a in args) + "\n")

#     # -------------------------------------------------------
#     # MAIN LOGIC
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     log_print("\n🧩 Starting computation of BETA and SE from Z-score...")

#     imp_z_col = sample_column_dict.get("imp_z_col", "NA")
#     beta_col  = sample_column_dict.get("beta_or_col", "NA")
#     se_col    = sample_column_dict.get("se_col", "NA")

#     # -------------------------------------------------------
#     # 1️⃣ Already present → skip
#     # -------------------------------------------------------
#     if beta_col in df.columns and se_col in df.columns:
#         log_print("✅ BETA+SE columns already present; skipping computation.")

#         summary = df.select([
#             pl.col(beta_col).min().alias("beta_min"),
#             pl.col(beta_col).max().alias("beta_max"),
#             pl.col(beta_col).mean().alias("beta_mean"),
#             pl.col(beta_col).std().alias("beta_std"),

#             pl.col(se_col).min().alias("se_min"),
#             pl.col(se_col).max().alias("se_max"),
#             pl.col(se_col).mean().alias("se_mean"),
#             pl.col(se_col).std().alias("se_std"),

#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         qc_info.update(summary)
#         qc_info["status"] = "existing_beta_se"

#         # write log
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         return df, qc_info, sample_column_dict

#     # -------------------------------------------------------
#     # 2️⃣ Check required inputs
#     # -------------------------------------------------------
#     if imp_z_col == "NA" or imp_z_col not in df.columns:
#         log_print("⚠️ No imputed Z-score column found — skipping.")
#         qc_info["status"] = "skipped_no_zscore"

#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         return df, qc_info, sample_column_dict

#     if eaf_col not in df.columns or "Neff" not in df.columns:
#         raise ValueError("❌ Missing required columns: 'eaf_col' and 'Neff' must be present.")

#     # -------------------------------------------------------
#     # 3️⃣ Ensure numeric types
#     # -------------------------------------------------------
#     df = df.with_columns([
#         pl.col(imp_z_col).cast(pl.Float64, strict=False),
#         pl.col(eaf_col).cast(pl.Float64, strict=False),
#         pl.col("Neff").cast(pl.Float64, strict=False),
#     ])

#     # -------------------------------------------------------
#     # 4️⃣ Z-score summary
#     # -------------------------------------------------------
#     z_summary = df.select([
#         pl.col(imp_z_col).min().alias("z_min"),
#         pl.col(imp_z_col).max().alias("z_max"),
#         pl.col(imp_z_col).mean().alias("z_mean"),
#         pl.col(imp_z_col).std().alias("z_std"),
#         pl.len().alias("total")
#     ]).to_dicts()[0]

#     log_print(
#         f"📊 Z summary: min={z_summary['z_min']:.6f}, "
#         f"max={z_summary['z_max']:.6f}, "
#         f"mean={z_summary['z_mean']:.6f}, std={z_summary['z_std']:.6f}"
#     )

#     qc_info.update(z_summary)

#     # -------------------------------------------------------
#     # 5️⃣ Compute BETA
#     # -------------------------------------------------------
#     if beta_col == "NA" or beta_col not in df.columns:
#         log_print("🧮 Computing BETA from Z-score... using formula:")
#         log_print("Beta = z / sqrt(2p(1− p)(n + z^2)) : PMID: 29763751")

#         df = df.with_columns(
#             (
#                 pl.col(imp_z_col) /
#                 (2 * pl.col(eaf_col) * (1 - pl.col(eaf_col)) *
#                  (pl.col("Neff") + pl.col(imp_z_col)**2)).sqrt()
#             ).alias("BETA")
#         )

#         beta_col = "BETA"
#         sample_column_dict["beta_or_col"] = "BETA"
#         qc_info["beta_computed"] = True

#     else:
#         log_print("ℹ️ BETA already present — skipping.")
#         qc_info["beta_computed"] = False

#     # -------------------------------------------------------
#     # 6️⃣ Compute SE
#     # -------------------------------------------------------
#     if se_col == "NA" or se_col not in df.columns:
#         log_print("🧮 Computing SE from Z-score... using formula:")
#         log_print("SE = 1 / sqrt(2p(1− p)(n + z^2)) :PMID: 29763751")

#         df = df.with_columns(
#             (
#                 1 /
#                 (2 * pl.col(eaf_col) * (1 - pl.col(eaf_col)) *
#                  (pl.col("Neff") + pl.col(imp_z_col)**2)).sqrt()
#             ).alias("SE")
#         )

#         se_col = "SE"
#         sample_column_dict["se_col"] = "SE"
#         qc_info["se_computed"] = True

#     else:
#         log_print("ℹ️ SE already present — skipping.")
#         qc_info["se_computed"] = False

#     # -------------------------------------------------------
#     # 7️⃣ Final BETA/SE summary
#     # -------------------------------------------------------
#     summary = df.select([
#         pl.col(beta_col).min().alias("beta_min"),
#         pl.col(beta_col).max().alias("beta_max"),
#         pl.col(beta_col).mean().alias("beta_mean"),
#         pl.col(beta_col).std().alias("beta_std"),

#         pl.col(se_col).min().alias("se_min"),
#         pl.col(se_col).max().alias("se_max"),
#         pl.col(se_col).mean().alias("se_mean"),
#         pl.col(se_col).std().alias("se_std"),

#         pl.len().alias("total")
#     ]).to_dicts()[0]

#     qc_info.update(summary)

#     log_print(
#         f"📈 BETA summary: min={summary['beta_min']:.6e}, "
#         f"max={summary['beta_max']:.6e}, mean={summary['beta_mean']:.6e}, std={summary['beta_std']:.6e}"
#     )
#     log_print(
#         f"📉 SE summary:   min={summary['se_min']:.6e}, "
#         f"max={summary['se_max']:.6e}, mean={summary['se_mean']:.6e}, std={summary['se_std']:.6e}"
#     )

#     log_print("🎯 Completed BETA/SE computation from Z-score.\n")

#     # -------------------------------------------------------
#     # Save logs
#     # -------------------------------------------------------
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return df, qc_info, sample_column_dict
