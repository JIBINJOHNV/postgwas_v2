import polars as pl
from scipy.stats import norm
from typing import Tuple, Dict
import numpy as np
from pathlib import Path
import io


def calculate_se_from_beta_pvalue(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    tail: int = 2
) -> Tuple[pl.DataFrame, Dict, dict]:

    # -------------------------------------------------------
    # INTERNAL LOG BUFFER
    # -------------------------------------------------------
    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # -------------------------------------------------------
    # LOG FILE SETUP
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_se.log"

    qc_info = {"initial_variants": df.height}

    try:
        # -------------------------------------------------------
        # MAIN
        # -------------------------------------------------------
        log_print(f"\n📊 [Chr {chromosome}] Starting SE calculation...")

        beta_col = sample_column_dict.get("beta_col", "NA")
        se_col   = sample_column_dict.get("se_col", "NA")
        pval_col = sample_column_dict.get("pval_col", "NA")

        # -------------------------------------------------------
        # VALIDATION
        # -------------------------------------------------------
        if (
            beta_col == "NA" or
            pval_col == "NA" or
            beta_col not in df.columns or
            pval_col not in df.columns
        ):
            raise ValueError("❌ Missing BETA or PVAL column.")

        # -------------------------------------------------------
        # SKIP IF SE EXISTS
        # -------------------------------------------------------
        if se_col != "NA" and se_col in df.columns:
            log_print(f"              ℹ️ SE column '{se_col}' exists → skipping.")

            stats = df.select([
                pl.col(se_col).min().alias("min"),
                pl.col(se_col).max().alias("max"),
                pl.col(se_col).mean().alias("mean"),
                pl.col(se_col).median().alias("median"),
                pl.len().alias("total")
            ]).to_dicts()[0]

            qc_info.update({
                "status": "SE already present",
                "se_min": stats["min"],
                "se_max": stats["max"],
                "se_mean": stats["mean"],
                "se_median": stats["median"],
                "total_variants": stats["total"],
            })

            return df, qc_info, sample_column_dict

        # -------------------------------------------------------
        # CAST
        # -------------------------------------------------------
        df = df.with_columns([
            pl.col(beta_col).cast(pl.Float64, strict=False),
            pl.col(pval_col).cast(pl.Float64, strict=False)
        ])

        # -------------------------------------------------------
        # QC
        # -------------------------------------------------------
        zero_pval = df.filter(pl.col(pval_col) == 0).height
        invalid_pval = df.filter(
            (pl.col(pval_col) <= 0) | (pl.col(pval_col) >= 1)
        ).height

        if zero_pval > 0 or invalid_pval > 0:
            log_print(f"              ⚠️ PVAL == 0           : {zero_pval:,}, chr{chromosome}")
            log_print(f"              ⚠️ Invalid PVAL        : {invalid_pval:,}, chr{chromosome}")

        # -------------------------------------------------------
        # SAFE PVAL
        # -------------------------------------------------------
        df = df.with_columns(
            pl.when(pl.col(pval_col) <= 0)
              .then(1e-300)
              .when(pl.col(pval_col) >= 1)
              .then(0.999999999)
              .otherwise(pl.col(pval_col))
              .alias("_pval_safe")
        )

        log_print(f"              📈 Computing SE using norm.ppf... chr{chromosome} using tail = {tail}")
        log_print(f"              ℹ️ Beta column: '{beta_col}', P-value column: '{pval_col}'")

        # ------------------------------------------------------- 
        # Compute Z
        # -------------------------------------------------------
        z = norm.ppf(df["_pval_safe"].to_numpy() / tail)
        z[z == 0] = np.nan
        se = np.abs(df[beta_col].to_numpy() / z)

        # -------------------------------------------------------
        # Attach back
        # -------------------------------------------------------
        df = df.with_columns(pl.Series("SE", se))
        sample_column_dict["se_col"] = "SE"

        # -------------------------------------------------------
        # QC SUMMARY
        # -------------------------------------------------------
        stats = df.select([
            pl.col("SE").min().alias("min"),
            pl.col("SE").max().alias("max"),
            pl.col("SE").mean().alias("mean"),
            pl.col("SE").median().alias("median"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update({
            "status": "SE calculated",
            "se_min": stats["min"],
            "se_max": stats["max"],
            "se_mean": stats["mean"],
            "se_median": stats["median"],
            "total_variants": stats["total"],
        })

        log_print(
            f"              ✅ SE stats → min={stats['min']:.6f}, max={stats['max']:.6f}, "
            f"mean={stats['mean']:.6f}, median={stats['median']:.6f}"
        )

        log_print("              🎯 Done.\n")

        return df, qc_info, sample_column_dict

    except Exception as e:
        log_print(f"❌ ERROR in chr{chromosome}: {str(e)}")
        raise

    finally:
        # -------------------------------------------------------
        # ALWAYS WRITE LOG FILE
        # -------------------------------------------------------
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        # -------------------------------------------------------
        # ALWAYS PRINT (NO MIXING)
        # -------------------------------------------------------
        print(log_buffer.getvalue(), flush=True)
        

# def calculate_se_from_beta_pvalue(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     tail: int = 2
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     """
#     Compute Standard Error (SE) from BETA and processed P-values when SE is missing.
#     This version includes dedicated multiprocessing-safe per-chromosome logging.
#     """

#     # -------------------------------------------------------
#     # INTERNAL LOG BUFFER (safe in multiprocessing)
#     # -------------------------------------------------------
#     log_buffer = io.StringIO()

#     def log_print(*args, verbose=False):
#         """Write to buffer; optionally print to console."""
#         msg = " ".join(str(a) for a in args)
#         # always write to buffer
#         log_buffer.write(msg + "\n")
#         # print only if verbose=True
#         if verbose:
#             print(msg, flush=True)
            

#     # -------------------------------------------------------
#     # Create logfile path
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_se_from_beta_and_pvalue.log"


#     # -------------------------------------------------------
#     # MAIN LOGIC
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     log_print("\n🧩 Starting SE calculation from BETA and processed p-values...")

#     beta_col = sample_column_dict.get("beta_or_col", "NA")
#     se_col   = sample_column_dict.get("se_col", "NA")
#     pval_col = sample_column_dict.get("pval_col", "NA")

#     # --- Validate required columns ---
#     if (
#         beta_col == "NA" or
#         pval_col == "NA" or
#         beta_col not in df.columns or
#         pval_col not in df.columns
#     ):
#         raise ValueError("❌ Missing required columns: BETA and/or p-value for SE calculation.")

#     # === Case: SE is already provided ===
#     if se_col != "NA" and se_col in df.columns:
#         log_print(f"ℹ️ SE column '{se_col}' already exists — skipping recalculation.")
#         qc_info["status"] = "SE already present"

#         stats = df.select([
#             pl.col(se_col).min().alias("min"),
#             pl.col(se_col).max().alias("max"),
#             pl.col(se_col).mean().alias("mean"),
#             pl.col(se_col).median().alias("median"),
#             pl.len().alias("total")
#         ]).to_dicts()[0]

#         qc_info.update({
#             "se_min": stats["min"],
#             "se_max": stats["max"],
#             "se_mean": stats["mean"],
#             "se_median": stats["median"],
#             "total_variants": stats["total"],
#         })

#         # WRITE LOG FILE
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         return df, qc_info, sample_column_dict

#     # --- Cast required columns safely ---
#     df = df.with_columns([
#         pl.col(beta_col).cast(pl.Float64, strict=False),
#         pl.col(pval_col).cast(pl.Float64, strict=False)
#     ])

#     log_print("📈 Computing SE using scipy.stats.norm.ppf()...")
#     log_print(f"🧮 Formula: SE = abs({beta_col} / norm.ppf({pval_col} / {tail}))")
#     log_print("ℹ️  Note: Using Inverse Normal Cumulative Distribution Function (Probit).")
#     # Convert to pandas for norm.ppf
#     pdf = df.select([beta_col, pval_col]).to_pandas()

#     # Protect norm.ppf against infinite results
#     pdf[pval_col] = pdf[pval_col].clip(lower=1e-300, upper=0.999999999)

#     # Main SE formula
#     pdf["SE"] = np.abs(pdf[beta_col] / norm.ppf(pdf[pval_col] / tail))

#     # Merge back into Polars
#     df = df.with_columns(pl.Series("SE", pdf["SE"].astype(float)))
#     sample_column_dict["se_col"] = "SE"

#     # --- QC summary ---
#     stats = df.select([
#         pl.col("SE").min().alias("min"),
#         pl.col("SE").max().alias("max"),
#         pl.col("SE").mean().alias("mean"),
#         pl.col("SE").median().alias("median"),
#         pl.len().alias("total")
#     ]).to_dicts()[0]

#     qc_info.update({
#         "se_min": stats["min"],
#         "se_max": stats["max"],
#         "se_mean": stats["mean"],
#         "se_median": stats["median"],
#         "total_variants": stats["total"],
#         "status": "SE calculated from BETA and p-value"
#     })

#     log_print(
#         f"✅ SE calculated: min={stats['min']:.6f}, max={stats['max']:.6f}, "
#         f"mean={stats['mean']:.6f}, median={stats['median']:.6f} "
#         f"over {stats['total']:,} variants."
#     )

#     log_print("🎯 SE calculation complete.\n")

#     # -------------------------------------------------------
#     # WRITE LOG FILE
#     # -------------------------------------------------------
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return df, qc_info, sample_column_dict
