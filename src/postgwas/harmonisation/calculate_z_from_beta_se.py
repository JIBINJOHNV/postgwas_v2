import polars as pl
from pathlib import Path
import io
from typing import Tuple, Dict
import sys

def calculate_z_from_beta_se(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    remove_beta_zero: bool = True
) -> Tuple[pl.DataFrame, Dict, dict]:

    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # subfolder for invalid rows (beta, se, z)
    invalid_info_dir = Path(output_dir) / "invalid_beta_se_zscore"
    invalid_info_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

    log_buffer = io.StringIO()

    def log_print(*args, indent: int = 0):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(("\t" * indent) + msg + "\n")

    qc_info = {"initial_variants": df.height}

    try:
        log_print(f"📊 [Chr {chromosome}] Starting Z-score computation (BETA / SE)...")

        beta_col = sample_column_dict.get("beta_col", "NA")
        se_col   = sample_column_dict.get("se_col", "NA")

        # ==========================================================
        # CHECK EXISTING Z COLUMN
        # ==========================================================
        existing_z_col = sample_column_dict.get("imp_z_col", "NA")

        has_existing_z = (
            existing_z_col != "NA"
            and existing_z_col in df.columns
        )

        if has_existing_z:
            log_print(
                f"✅ Z-score already present in the file provided in column '{existing_z_col}' — BETA/SE QC will still be applied.",
                indent=4
            )

            df = df.with_columns(
                pl.col(existing_z_col)
                .cast(pl.Float64, strict=False)
                .alias(existing_z_col)
            )

            sample_column_dict["imp_z_col"] = existing_z_col

            log_print(
                f"🎯 Using existing Z-score column: {existing_z_col}",
                indent=4
            )

        # ==========================================================
        # ALWAYS APPLY BETA/SE QC
        # ==========================================================
        if (
            beta_col == "NA" or
            se_col == "NA" or
            beta_col not in df.columns or
            se_col not in df.columns
        ):
            raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

        # Cast
        df = df.with_columns([
            pl.col(beta_col).cast(pl.Float64, strict=False),
            pl.col(se_col).cast(pl.Float64, strict=False)
        ])

        if not has_existing_z:
            log_print(
                f"✅ Z-score not present in the file provided — calculating Z-score from BETA column {beta_col} and SE column {se_col}",
                indent=4
            )

        # Conditions
        cond_beta_null  = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
        cond_se_null    = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
        cond_se_invalid = (pl.col(se_col) <= 0)

        # OPTIONAL beta=0 removal
        cond_beta_zero = (pl.col(beta_col) == 0) if remove_beta_zero else pl.lit(False)

        # Partition removed
        df_beta_null  = df.filter(cond_beta_null)
        df_se_null    = df.filter(~cond_beta_null & cond_se_null)
        df_beta_zero  = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
        df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

        df_removed = pl.concat([
            df_beta_null,
            df_se_null,
            df_beta_zero,
            df_se_invalid
        ]).unique()

        # Keep valid
        condition_valid = (
            (~cond_beta_null) &
            (~cond_se_null) &
            (~cond_beta_zero) &
            (~cond_se_invalid)
        )

        df = df.filter(condition_valid)

        # Counts
        n_initial     = qc_info["initial_variants"]
        n_beta_null   = df_beta_null.height
        n_se_null     = df_se_null.height
        n_beta_zero   = df_beta_zero.height
        n_se_invalid  = df_se_invalid.height

        # Compute Z only if existing Z is not present
        if not has_existing_z:
            df = df.with_columns(
                pl.when(pl.col(se_col) > 1e-12)
                .then(pl.col(beta_col) / pl.col(se_col))
                .otherwise(None)
                .alias("imp_z_col")
            )

            sample_column_dict["imp_z_col"] = "imp_z_col"
        else:
            log_print(
                f"✅ Existing Z-score retained from column: {existing_z_col}",
                indent=4
            )

        # Remove invalid Z
        z_col = sample_column_dict["imp_z_col"]

        df_z_invalid = df.filter(
            pl.col(z_col).is_null() |
            ~pl.col(z_col).is_finite()
        )
        n_z_invalid = df_z_invalid.height

        if n_z_invalid > 0:
            df_removed = pl.concat([df_removed, df_z_invalid]).unique()

        df = df.filter(
            pl.col(z_col).is_not_null() &
            pl.col(z_col).is_finite()
        )

        n_removed    = df_removed.height
        n_final      = df.height
        pct_removed  = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

        # Save removed
        invalid_file = invalid_info_dir / f"{gwas_outputname}_chr{chromosome}_invalid_beta_se.tsv"
        if df_removed.height > 0:
            df_removed.write_csv(invalid_file, separator="\t")
            log_print(f"🧾 Invalid variants saved → {invalid_file}", indent=4)

        # Logging
        log_print(f"📊 [Chr {chromosome}] Removed {n_removed:,} variants", indent=4)
        log_print(f"Initial variants : {n_initial:,}", indent=4)
        log_print(f"Removed total    : {n_removed:,} ({pct_removed:.2f}%)", indent=4)
        log_print(f"Null BETA        : {n_beta_null:,}", indent=4)
        log_print(f"Null SE          : {n_se_null:,}", indent=4)
        log_print(f"BETA = 0         : {n_beta_zero:,}", indent=4)
        log_print(f"SE <= 0          : {n_se_invalid:,}", indent=4)
        log_print(f"Invalid Z        : {n_z_invalid:,}", indent=4)
        log_print(f"Remaining        : {n_final:,}", indent=4)

        qc_info.update({
            "variants_removed_due_to_null_beta": n_beta_null,
            "variants_removed_due_to_null_se": n_se_null,
            "variants_removed_due_to_zero_beta": n_beta_zero,
            "variants_removed_due_to_invalid_se": n_se_invalid,
            "variants_removed_due_to_invalid_z": n_z_invalid,
            "variants_removed_invalid_beta_se": n_removed,
            "variants_after_filter_invalid_beta_se": n_final,
            "beta_zero_removed_flag": remove_beta_zero
        })

        if df.height == 0:
            raise RuntimeError("❌ No valid variants left after filtering.")

        # ==========================================================
        # COMMON Z SUMMARY
        # ==========================================================
        z_col = sample_column_dict["imp_z_col"]

        z_summary = df.select([
            pl.col(z_col).min().alias("z_min"),
            pl.col(z_col).max().alias("z_max"),
            pl.col(z_col).mean().alias("z_mean"),
            pl.col(z_col).std().alias("z_std"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        log_print(
            f"✅ Z summary: min={z_summary['z_min']:.6f}, "
            f"max={z_summary['z_max']:.6f}, "
            f"mean={z_summary['z_mean']:.6f}, "
            f"std={z_summary['z_std']:.6f} "
            f"(n={z_summary['total']:,})",
            indent=4
        )

        qc_info.update({
            "z_min": z_summary["z_min"],
            "z_max": z_summary["z_max"],
            "z_mean": z_summary["z_mean"],
            "z_std": z_summary["z_std"],
            "total_variants": z_summary["total"],
        })

        log_print(f"🎯  Z-score computation completed for chromosome chr{chromosome}", indent=4)

    except Exception as e:
        log_print(f"❌ ERROR: {str(e)}")
        raise

    finally:
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        print(log_buffer.getvalue(), flush=True)

    return df, qc_info, sample_column_dict









# def calculate_z_from_beta_se(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     remove_beta_zero: bool = True
# ) -> Tuple[pl.DataFrame, Dict, dict]:

#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     # subfolder for invalid rows (beta, se, z)
#     invalid_info_dir = Path(output_dir) / "invalid_beta_se_zscore"
#     invalid_info_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

#     log_buffer = io.StringIO()

#     def log_print(*args, indent: int = 0):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(("\t" * indent) + msg + "\n")

#     qc_info = {"initial_variants": df.height}

#     try:
#         log_print(f"📊 [Chr {chromosome}] Starting Z-score computation (BETA / SE)...")

#         beta_col = sample_column_dict.get("beta_col", "NA")
#         se_col   = sample_column_dict.get("se_col", "NA")

#         # ==========================================================
#         # USE EXISTING Z COLUMN IF ALREADY PRESENT
#         # ==========================================================
#         existing_z_col = sample_column_dict.get("imp_z_col", "NA")
#         used_existing_z = False

#         if (
#             existing_z_col != "NA"
#             and existing_z_col in df.columns
#         ):
#             log_print(
#                 f"✅ Z-score already present in the file provided in column '{existing_z_col}' — skipping calculation.",
#                 indent=4
#             )

#             df = df.with_columns(
#                 pl.col(existing_z_col)
#                 .cast(pl.Float64, strict=False)
#                 .alias(existing_z_col)
#             )

#             sample_column_dict["imp_z_col"] = existing_z_col
#             used_existing_z = True

#             log_print(
#                 f"🎯 Using existing Z-score column: {existing_z_col}",
#                 indent=4
#             )

#         # ==========================================================
#         # CALCULATE Z ONLY IF NOT ALREADY PRESENT
#         # ==========================================================
#         if not used_existing_z:

#             if (
#                 beta_col == "NA" or
#                 se_col == "NA" or
#                 beta_col not in df.columns or
#                 se_col not in df.columns
#             ):
#                 raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

#             # Cast
#             df = df.with_columns([
#                 pl.col(beta_col).cast(pl.Float64, strict=False),
#                 pl.col(se_col).cast(pl.Float64, strict=False)
#             ])
#             log_print(
#                 f"✅ Z-score not present in the file provided — calculating Z-score from BETA column {beta_col} and SE column {se_col}",
#                 indent=4
#             )
#             # Conditions
#             cond_beta_null  = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
#             cond_se_null    = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
#             cond_se_invalid = (pl.col(se_col) <= 0)

#             # OPTIONAL beta=0 removal
#             cond_beta_zero = (pl.col(beta_col) == 0) if remove_beta_zero else pl.lit(False)

#             # Partition removed
#             df_beta_null  = df.filter(cond_beta_null)
#             df_se_null    = df.filter(~cond_beta_null & cond_se_null)
#             df_beta_zero  = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
#             df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

#             df_removed = pl.concat([
#                 df_beta_null,
#                 df_se_null,
#                 df_beta_zero,
#                 df_se_invalid
#             ]).unique()

#             # Keep valid
#             condition_valid = (
#                 (~cond_beta_null) &
#                 (~cond_se_null) &
#                 (~cond_beta_zero) &
#                 (~cond_se_invalid)
#             )

#             df = df.filter(condition_valid)

#             # Counts
#             n_initial     = qc_info["initial_variants"]
#             n_beta_null   = df_beta_null.height
#             n_se_null     = df_se_null.height
#             n_beta_zero   = df_beta_zero.height
#             n_se_invalid  = df_se_invalid.height

#             # Compute Z safely
#             df = df.with_columns(
#                 pl.when(pl.col(se_col) > 1e-12)
#                 .then(pl.col(beta_col) / pl.col(se_col))
#                 .otherwise(None)
#                 .alias("imp_z_col")
#             )

#             sample_column_dict["imp_z_col"] = "imp_z_col"

#             # Remove invalid Z
#             df_z_invalid = df.filter(
#                 pl.col("imp_z_col").is_null() |
#                 ~pl.col("imp_z_col").is_finite()
#             )
#             n_z_invalid = df_z_invalid.height

#             if n_z_invalid > 0:
#                 df_removed = pl.concat([df_removed, df_z_invalid]).unique()

#             df = df.filter(
#                 pl.col("imp_z_col").is_not_null() &
#                 pl.col("imp_z_col").is_finite()
#             )

#             n_removed    = df_removed.height
#             n_final      = df.height
#             pct_removed  = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#             # Save removed
#             invalid_file = invalid_info_dir / f"{gwas_outputname}_chr{chromosome}_invalid_beta_se.tsv"
#             if df_removed.height > 0:
#                 df_removed.write_csv(invalid_file, separator="\t")
#                 log_print(f"🧾 Invalid variants saved → {invalid_file}", indent=4)

#             # Logging
#             log_print(f"📊 [Chr {chromosome}] Removed {n_removed:,} variants", indent=4)
#             log_print(f"Initial variants : {n_initial:,}", indent=4)
#             log_print(f"Removed total    : {n_removed:,} ({pct_removed:.2f}%)", indent=4)
#             log_print(f"Null BETA        : {n_beta_null:,}", indent=4)
#             log_print(f"Null SE          : {n_se_null:,}", indent=4)
#             log_print(f"BETA = 0         : {n_beta_zero:,}", indent=4)
#             log_print(f"SE <= 0          : {n_se_invalid:,}", indent=4)
#             log_print(f"Invalid Z        : {n_z_invalid:,}", indent=4)
#             log_print(f"Remaining        : {n_final:,}", indent=4)

#             qc_info.update({
#                 "variants_removed_due_to_null_beta": n_beta_null,
#                 "variants_removed_due_to_null_se": n_se_null,
#                 "variants_removed_due_to_zero_beta": n_beta_zero,
#                 "variants_removed_due_to_invalid_se": n_se_invalid,
#                 "variants_removed_due_to_invalid_z": n_z_invalid,
#                 "variants_removed_invalid_beta_se": n_removed,
#                 "variants_after_filter_invalid_beta_se": n_final,
#                 "beta_zero_removed_flag": remove_beta_zero
#             })

#             if df.height == 0:
#                 raise RuntimeError("❌ No valid variants left after filtering.")

#         # ==========================================================
#         # COMMON Z SUMMARY (FOR BOTH EXISTING + COMPUTED Z)
#         # ==========================================================
#         z_col = sample_column_dict["imp_z_col"]

#         z_summary = df.select([
#             pl.col(z_col).min().alias("z_min"),
#             pl.col(z_col).max().alias("z_max"),
#             pl.col(z_col).mean().alias("z_mean"),
#             pl.col(z_col).std().alias("z_std"),
#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         log_print(
#             f"✅ Z summary: min={z_summary['z_min']:.6f}, "
#             f"max={z_summary['z_max']:.6f}, "
#             f"mean={z_summary['z_mean']:.6f}, "
#             f"std={z_summary['z_std']:.6f} "
#             f"(n={z_summary['total']:,})",
#             indent=4
#         )

#         qc_info.update({
#             "z_min": z_summary["z_min"],
#             "z_max": z_summary["z_max"],
#             "z_mean": z_summary["z_mean"],
#             "z_std": z_summary["z_std"],
#             "total_variants": z_summary["total"],
#         })

#         log_print(f"🎯  Z-score computation completed for chromosome chr{chromosome}", indent=4)

#     except Exception as e:
#         log_print(f"❌ ERROR: {str(e)}")
#         raise

#     finally:
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         print(log_buffer.getvalue(), flush=True)

#     return df, qc_info, sample_column_dict

# def calculate_z_from_beta_se(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     remove_beta_zero: bool = True   # 🔴 NEW OPTION
# ) -> Tuple[pl.DataFrame, Dict, dict]:

#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     # 🔴 subfolder for invalid rows (beta, se, z)
#     invalid_info_dir = Path(output_dir) / "invalid_beta_se_zscore"
#     invalid_info_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

#     log_buffer = io.StringIO()

#     def log_print(*args, indent: int = 0):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(("\t" * indent) + msg + "\n")

#     qc_info = {"initial_variants": df.height}

#     try:
#         log_print(f"📊 [Chr {chromosome}] Starting Z-score computation (BETA / SE)...")

#         beta_col = sample_column_dict.get("beta_col", "NA")
#         se_col   = sample_column_dict.get("se_col", "NA")

#         if (
#             beta_col == "NA" or
#             se_col == "NA" or
#             beta_col not in df.columns or
#             se_col not in df.columns
#         ):
#             raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

#         # Cast
#         df = df.with_columns([
#             pl.col(beta_col).cast(pl.Float64, strict=False),
#             pl.col(se_col).cast(pl.Float64, strict=False)
#         ])

#         # Conditions
#         cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
#         cond_se_null   = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
#         cond_se_invalid = (pl.col(se_col) <= 0)

#         # 🔴 OPTIONAL beta=0 removal
#         cond_beta_zero = (pl.col(beta_col) == 0) if remove_beta_zero else pl.lit(False)

#         # Partition removed
#         df_beta_null = df.filter(cond_beta_null)
#         df_se_null   = df.filter(~cond_beta_null & cond_se_null)
#         df_beta_zero = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
#         df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

#         df_removed = pl.concat([
#             df_beta_null,
#             df_se_null,
#             df_beta_zero,
#             df_se_invalid
#         ]).unique()

#         # Keep valid
#         condition_valid = (
#             (~cond_beta_null) &
#             (~cond_se_null) &
#             (~cond_beta_zero) &
#             (~cond_se_invalid)
#         )

#         df = df.filter(condition_valid)

#         # Counts
#         n_initial = qc_info["initial_variants"]

#         n_beta_null = df_beta_null.height
#         n_se_null   = df_se_null.height
#         n_beta_zero = df_beta_zero.height
#         n_se_invalid = df_se_invalid.height

#         # Compute Z safely
#         df = df.with_columns(
#             pl.when(pl.col(se_col) > 1e-12)
#             .then(pl.col(beta_col) / pl.col(se_col))
#             .otherwise(None)
#             .alias("imp_z_col")
#         )

#         sample_column_dict["imp_z_col"] = "imp_z_col"

#         # Remove invalid Z
#         df_z_invalid = df.filter(pl.col("imp_z_col").is_null() | ~pl.col("imp_z_col").is_finite())
#         n_z_invalid = df_z_invalid.height

#         if n_z_invalid > 0:
#             df_removed = pl.concat([df_removed, df_z_invalid]).unique()

#         df = df.filter(pl.col("imp_z_col").is_not_null() & pl.col("imp_z_col").is_finite())

#         n_removed = df_removed.height
#         n_final = df.height
#         pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#         # Save removed
#         invalid_file = invalid_info_dir / f"{gwas_outputname}_chr{chromosome}_invalid_beta_se.tsv"
#         if df_removed.height > 0:
#             df_removed.write_csv(invalid_file, separator="\t")
#             log_print(f"🧾 Invalid variants saved → {invalid_file}", indent=4)

#         # Logging
#         log_print(f"📊 [Chr {chromosome}] Removed {n_removed:,} variants", indent=4)
#         log_print(f"Initial variants : {n_initial:,}", indent=4)
#         log_print(f"Removed total    : {n_removed:,} ({pct_removed:.2f}%)", indent=4)
#         log_print(f"Null BETA        : {n_beta_null:,}", indent=4)
#         log_print(f"Null SE          : {n_se_null:,}", indent=4)
#         log_print(f"BETA = 0         : {n_beta_zero:,}", indent=4)
#         log_print(f"SE <= 0          : {n_se_invalid:,}", indent=4)
#         log_print(f"Invalid Z        : {n_z_invalid:,}", indent=4)
#         log_print(f"Remaining        : {n_final:,}", indent=4)

#         qc_info.update({
#             "variants_removed_due_to_null_beta": n_beta_null,
#             "variants_removed_due_to_null_se": n_se_null,
#             "variants_removed_due_to_zero_beta": n_beta_zero,
#             "variants_removed_due_to_invalid_se": n_se_invalid,
#             "variants_removed_due_to_invalid_z": n_z_invalid,
#             "variants_removed_invalid_beta_se": n_removed,
#             "variants_after_filter_invalid_beta_se": n_final,
#             "beta_zero_removed_flag": remove_beta_zero   # 🔴 NEW QC FLAG
#         })

#         if df.height == 0:
#             raise RuntimeError("❌ No valid variants left after filtering.")

#         # Z summary
#         z_summary = df.select([
#             pl.col("imp_z_col").min().alias("z_min"),
#             pl.col("imp_z_col").max().alias("z_max"),
#             pl.col("imp_z_col").mean().alias("z_mean"),
#             pl.col("imp_z_col").std().alias("z_std"),
#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         log_print(
#             f"✅ Z summary: min={z_summary['z_min']:.6f}, "
#             f"max={z_summary['z_max']:.6f}, "
#             f"mean={z_summary['z_mean']:.6f}, "
#             f"std={z_summary['z_std']:.6f} "
#             f"(n={z_summary['total']:,})",
#             indent=4
#         )

#         qc_info.update({
#             "z_min": z_summary["z_min"],
#             "z_max": z_summary["z_max"],
#             "z_mean": z_summary["z_mean"],
#             "z_std": z_summary["z_std"],
#             "total_variants": z_summary["total"],
#         })

#         log_print(f"🎯  Z-score computation completed for chromosome chr{chromosome}", indent=4)

#     except Exception as e:
#         log_print(f"❌ ERROR: {str(e)}")
#         raise

#     finally:
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         print(log_buffer.getvalue(), flush=True)

#     return df, qc_info, sample_column_dict



# def calculate_z_from_beta_se(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, Dict, dict]:

#     import io
#     import sys
#     from pathlib import Path

#     # -------------------------------------------------------
#     # Setup logfile
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

#     # -------------------------------------------------------
#     # Internal log buffer
#     # -------------------------------------------------------
#     log_buffer = io.StringIO()

#     def log_print(*args, indent: int = 0):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(("\t" * indent) + msg + "\n")

#     # -------------------------------------------------------
#     # MAIN LOGIC
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     try:
#         log_print(f"📊 [Chr {chromosome}] Starting Z-score computation (BETA / SE)...")

#         beta_col = sample_column_dict.get("beta_col", "NA")
#         se_col   = sample_column_dict.get("se_col", "NA")

#         # -------------------------------------------------------
#         # Validate input columns
#         # -------------------------------------------------------
#         if (
#             beta_col == "NA" or
#             se_col == "NA" or
#             beta_col not in df.columns or
#             se_col not in df.columns
#         ):
#             raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

#         # -------------------------------------------------------
#         # Cast to numeric
#         # -------------------------------------------------------
#         df = df.with_columns([
#             pl.col(beta_col).cast(pl.Float64, strict=False),
#             pl.col(se_col).cast(pl.Float64, strict=False)
#         ])

#         # -------------------------------------------------------
#         # Filtering conditions
#         # -------------------------------------------------------
#         cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
#         cond_se_null   = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
#         cond_beta_zero = (pl.col(beta_col) == 0)
#         cond_se_invalid = (pl.col(se_col) <= 0)

#         # -------------------------------------------------------
#         # Partition removed variants
#         # -------------------------------------------------------
#         df_beta_null = df.filter(cond_beta_null)
#         df_se_null   = df.filter(~cond_beta_null & cond_se_null)
#         df_beta_zero = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
#         df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

#         df_removed = pl.concat([
#             df_beta_null,
#             df_se_null,
#             df_beta_zero,
#             df_se_invalid
#         ]).unique()

#         # -------------------------------------------------------
#         # Keep valid variants
#         # -------------------------------------------------------
#         condition_valid = (
#             (~cond_beta_null) &
#             (~cond_se_null) &
#             (~cond_beta_zero) &
#             (~cond_se_invalid)
#         )

#         df = df.filter(condition_valid)

#         # -------------------------------------------------------
#         # Counts
#         # -------------------------------------------------------
#         n_initial = qc_info["initial_variants"]

#         n_beta_null = df_beta_null.height
#         n_se_null   = df_se_null.height
#         n_beta_zero = df_beta_zero.height
#         n_se_invalid = df_se_invalid.height

#         n_removed = n_beta_null + n_se_null + n_beta_zero + n_se_invalid
#         n_final = df.height

#         pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#         # -------------------------------------------------------
#         # Log summary (buffered)
#         # -------------------------------------------------------
#         log_print(f"📊 [Chr {chromosome}] Removed {n_removed:,} variants", indent=4)
#         log_print(f"Initial variants : {n_initial:,}", indent=4)
#         log_print(f"Removed total    : {n_removed:,} ({pct_removed:.2f}%)", indent=4)
#         log_print(f"Null BETA        : {n_beta_null:,}", indent=4)
#         log_print(f"Null SE          : {n_se_null:,}", indent=4)
#         log_print(f"BETA = 0         : {n_beta_zero:,}", indent=4)
#         log_print(f"SE <= 0          : {n_se_invalid:,}", indent=4)
#         log_print(f"Remaining        : {n_final:,}", indent=4)

#         qc_info.update({
#             "variants_removed_due_to_null_beta": n_beta_null,
#             "variants_removed_due_to_null_se": n_se_null,
#             "variants_removed_due_to_zero_beta": n_beta_zero,
#             "variants_removed_due_to_invalid_se": n_se_invalid,
#             "variants_removed_invalid_beta_se": n_removed,
#             "variants_after_filter_invalid_beta_se": n_final
#         })

#         # -------------------------------------------------------
#         # Compute Z
#         # -------------------------------------------------------
#         df = df.with_columns(
#             pl.when(pl.col(se_col) > 0)
#             .then(pl.col(beta_col) / pl.col(se_col))
#             .otherwise(None)
#             .alias("imp_z_col")
#         )

#         sample_column_dict["imp_z_col"] = "imp_z_col"

#         # -------------------------------------------------------
#         # Z summary
#         # -------------------------------------------------------
#         z_summary = df.select([
#             pl.col("imp_z_col").min().alias("z_min"),
#             pl.col("imp_z_col").max().alias("z_max"),
#             pl.col("imp_z_col").mean().alias("z_mean"),
#             pl.col("imp_z_col").std().alias("z_std"),
#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         log_print(
#             f"✅ Z summary: min={z_summary['z_min']:.6f}, "
#             f"max={z_summary['z_max']:.6f}, "
#             f"mean={z_summary['z_mean']:.6f}, "
#             f"std={z_summary['z_std']:.6f} "
#             f"(n={z_summary['total']:,})",
#             indent=4
#         )

#         qc_info.update({
#             "z_min": z_summary["z_min"],
#             "z_max": z_summary["z_max"],
#             "z_mean": z_summary["z_mean"],
#             "z_std": z_summary["z_std"],
#             "total_variants": z_summary["total"], 
#         })

#         log_print(f"🎯  Z-scores Computation completed for chromosome chr{chromosome}", indent=4)

#     except Exception as e:
#         log_print(f"❌ ERROR: {str(e)}")
#         raise

#     finally:
#         # -------------------------------------------------------
#         # ALWAYS write + print logs
#         # -------------------------------------------------------
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         print(log_buffer.getvalue(), flush=True)

#     return df, qc_info, sample_column_dict



# def calculate_z_from_beta_se(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     """
#     Compute imputed Z-score ('imp_z_col') from BETA and SE columns.

#     Formula:
#         imp_z_col = BETA / SE

#     This version:
#       • Silent (no screen printing except summary)
#       • Writes full logs to: logs/{gwas_outputname}_chr{chromosome}_calculate_z.log
#       • Safe for multiprocessing
#     """

#     # -------------------------------------------------------
#     # Setup logfile
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

#     # -------------------------------------------------------
#     # Internal log buffer (multiprocessing-safe)
#     # -------------------------------------------------------
#     log_buffer = io.StringIO()

#     def log_print(*args):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n")

#     # -------------------------------------------------------
#     # MAIN LOGIC
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     log_print("\n🧩 Starting computation of imputed Z-scores (imp_z_col)...")

#     beta_col = sample_column_dict.get("beta_col", "NA")
#     se_col   = sample_column_dict.get("se_col", "NA")

#     # -------------------------------------------------------
#     # Validate input columns
#     # -------------------------------------------------------
#     if (
#         beta_col == "NA" or
#         se_col == "NA" or
#         beta_col not in df.columns or
#         se_col not in df.columns
#     ):
#         raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

#     # -------------------------------------------------------
#     # Cast to numeric
#     # -------------------------------------------------------
#     df = df.with_columns([
#         pl.col(beta_col).cast(pl.Float64, strict=False),
#         pl.col(se_col).cast(pl.Float64, strict=False)
#     ])

#     # -------------------------------------------------------
#     # Filtering conditions (with finite checks)
#     # -------------------------------------------------------
#     cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
#     cond_se_null   = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
#     cond_beta_zero = (pl.col(beta_col) == 0)
#     cond_se_invalid = (pl.col(se_col) <= 0)

#     # -------------------------------------------------------
#     # Partition removed variants
#     # -------------------------------------------------------
#     df_beta_null = df.filter(cond_beta_null)
#     df_se_null   = df.filter(~cond_beta_null & cond_se_null)
#     df_beta_zero = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
#     df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

#     # -------------------------------------------------------
#     # Combine removed variants (safe)
#     # -------------------------------------------------------
#     df_removed = pl.concat([
#         df_beta_null,
#         df_se_null,
#         df_beta_zero,
#         df_se_invalid
#     ]).unique()

#     # -------------------------------------------------------
#     # Keep valid variants
#     # -------------------------------------------------------
#     condition_valid = (
#         (~cond_beta_null) &
#         (~cond_se_null) &
#         (~cond_beta_zero) &
#         (~cond_se_invalid)
#     )

#     df = df.filter(condition_valid)

#     # -------------------------------------------------------
#     # Counts
#     # -------------------------------------------------------
#     n_initial = qc_info["initial_variants"]

#     n_beta_null = df_beta_null.height
#     n_se_null   = df_se_null.height
#     n_beta_zero = df_beta_zero.height
#     n_se_invalid = df_se_invalid.height

#     n_removed = n_beta_null + n_se_null + n_beta_zero + n_se_invalid
#     n_final = df.height

#     pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#     # -------------------------------------------------------
#     # Print summary
#     # -------------------------------------------------------
#     indent = "\t" * 4  # 4-level indentation

#     print(f"\n{indent}📊 [Chr {chromosome}] Removed {n_removed:,} variants with invalid Beta/SE.")
#     print(f"{indent}" + "-"*50)

#     print(f"{indent}• Initial variants : {n_initial:,}")
#     print(f"{indent}• Removed total    : {n_removed:,} ({pct_removed:.2f}%)")

#     print(f"{indent}    ├─ Null BETA   : {n_beta_null:,}")
#     print(f"{indent}    ├─ Null SE     : {n_se_null:,}")
#     print(f"{indent}    ├─ BETA = 0    : {n_beta_zero:,}")
#     print(f"{indent}    └─ SE <= 0     : {n_se_invalid:,}")

#     print(f"{indent}• Remaining        : {n_final:,}")
#     print(f"{indent}" + "-"*50)

#     # -------------------------------------------------------
#     # Log summary
#     # -------------------------------------------------------
#     log_print(f"[Chr {chromosome}] QC Summary:")
#     log_print(f"Initial: {n_initial}")
#     log_print(f"Removed: {n_removed} ({pct_removed:.2f}%)")
#     log_print(f"Null BETA: {n_beta_null}")
#     log_print(f"Null SE: {n_se_null}")
#     log_print(f"BETA=0: {n_beta_zero}")
#     log_print(f"SE<=0: {n_se_invalid}")
#     log_print(f"Remaining: {n_final}")

#     # -------------------------------------------------------
#     # QC info update
#     # -------------------------------------------------------
#     qc_info.update({
#         "variants_removed_due_to_null_beta": n_beta_null,
#         "variants_removed_due_to_null_se": n_se_null,
#         "variants_removed_due_to_zero_beta": n_beta_zero,
#         "variants_removed_due_to_invalid_se": n_se_invalid,
#         "variants_removed_invalid_beta_se": n_removed,
#         "variants_after_filter_invalid_beta_se": n_final
#     })

#     # -------------------------------------------------------
#     # Compute Z safely
#     # -------------------------------------------------------
#     log_print("📈 Computing imp_z_col = BETA / SE …")

#     df = df.with_columns(
#         pl.when(pl.col(se_col) > 0)
#         .then(pl.col(beta_col) / pl.col(se_col))
#         .otherwise(None)
#         .alias("imp_z_col")
#     )

#     sample_column_dict["imp_z_col"] = "imp_z_col"

#     # -------------------------------------------------------
#     # Z summary
#     # -------------------------------------------------------
#     z_summary = df.select([
#         pl.col("imp_z_col").min().alias("z_min"),
#         pl.col("imp_z_col").max().alias("z_max"),
#         pl.col("imp_z_col").mean().alias("z_mean"),
#         pl.col("imp_z_col").std().alias("z_std"),
#         pl.len().alias("total")
#     ]).to_dicts()[0]

#     log_print(
#         f"✅ Imputed Z-score summary: "
#         f"min={z_summary['z_min']:.6f}, "
#         f"max={z_summary['z_max']:.6f}, "
#         f"mean={z_summary['z_mean']:.6f}, "
#         f"std={z_summary['z_std']:.6f} "
#         f"(n={z_summary['total']:,})"
#     )

#     qc_info.update({
#         "z_min": z_summary["z_min"],
#         "z_max": z_summary["z_max"],
#         "z_mean": z_summary["z_mean"],
#         "z_std": z_summary["z_std"],
#         "total_variants": z_summary["total"],
#     })

#     log_print("🎯 Computation of imp_z_col completed.\n")

#     # -------------------------------------------------------
#     # Write logfile
#     # -------------------------------------------------------
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return df, qc_info, sample_column_dict
