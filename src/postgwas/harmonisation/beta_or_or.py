import polars as pl
from pathlib import Path
import io
from typing import Dict, Tuple

def is_beta_or_or(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict
) -> Tuple[pl.DataFrame, dict, dict]:
    """
    Determine whether the effect-size column contains Beta or Odds Ratio.
    If OR, convert to Beta using log(OR). Generate QC summary and full logfile.

    Logged version:
        • Dedicated logfile per chromosome
        • Multiprocessing-safe (no stdout)
        • All logs written only to logs/*.log
    """
    # -------------------------------------------------------
    # 1. Setup log buffer + file path
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_is_beta_or_or.log"

    buffer = io.StringIO()
    completed_successfully = False

    def log_print(*args, indent: int = 16):
        msg = " ".join(map(str, args))
        indent_str = " " * indent
        final = f"{indent_str}{msg}"
        buffer.write(final + "\n")

    try:
        log_print(f"📊 [Chr {chromosome}] Starting effect size detection [BETA / ODDS RATIO]...", indent=0)
        # -------------------------------------------------------
        # 2. Main logic 
        # -------------------------------------------------------
        qc_info = {"initial_variants": df.height}
        effect_col = sample_column_dict.get("beta_or_col", None)
        if not effect_col or effect_col not in df.columns:
            raise ValueError("❌ Effect size column ('beta_or_col') missing or not found in dataframe.")

        # Cast to float safely 
        # df = df.with_columns(pl.col(effect_col).cast(pl.Float64, strict=False)) 
        df = df.with_columns(pl.col(effect_col).cast(pl.Float64, strict=True))
        # -------------------------------------------------------
        # LOWEST 10 VALUES (sorted)
        # -------------------------------------------------------
        total = df.height
        negative_count = df.filter(pl.col(effect_col) < 0.0).height
        negative_frac = negative_count / total if total else 0.0

        qc_info.update({
            "effect_col": effect_col,
            "negative_value_count": negative_count,
            "negative_fraction": round(negative_frac, 6),
        })

        # -------------------------------------------------------
        # 3A. BETA detected
        # -------------------------------------------------------
        if negative_frac > 0.10:
            log_print(f"✅ Effect size Detected as Beta for chromosome {chromosome}: {negative_frac*100:.2f}% values are negative.")
            sample_column_dict["effect_type"] = "beta"
            sample_column_dict["beta_col"] = effect_col
            qc_info["effect_type"] = "beta"

            stats = (
                df.select([
                    pl.col(effect_col).min().alias("min_beta"),
                    pl.col(effect_col).max().alias("max_beta"),
                    pl.col(effect_col).mean().alias("mean_beta"),
                    pl.col(effect_col).std().alias("std_beta"),
                    pl.len().alias("total"),
                ])
                .to_dicts()[0]
            )

            qc_info.update({
                "min_beta": stats["min_beta"],
                "max_beta": stats["max_beta"],
                "mean_beta": stats["mean_beta"],
                "std_beta":  stats["std_beta"],
            })

            log_print(
                f"📈 Beta summary for chromosome {chromosome}: min={stats['min_beta']:.6f}, "
                f"max={stats['max_beta']:.6f}, mean={stats['mean_beta']:.6f}"
            )

        # -------------------------------------------------------
        # 3B. OR detected → convert
        # -------------------------------------------------------
        else:
            log_print(f"              ⚠️ Effect size  Detected as Odds Ratio: only {negative_frac*100:.2f}% are negative. for chromosome {chromosome}")
            # print(f"                ⚠️ Effect size Detected as Odds Ratio for chromosome {chromosome}: only {negative_frac*100:.2f}% are negative.")

            sample_column_dict["effect_type"] = "odds_ratio"
            qc_info["effect_type"] = "odds_ratio"

            # Pre-conversion stats
            pre_stats = (
                df.select([
                    pl.col(effect_col).min().alias("min_pre"),
                    pl.col(effect_col).max().alias("max_pre"),
                    pl.col(effect_col).mean().alias("mean_pre"),
                    pl.col(effect_col).std().alias("std_pre"),
                ])
                .to_dicts()[0]
            )

            qc_info.update({
                "pre_conversion_min": pre_stats["min_pre"],
                "pre_conversion_max": pre_stats["max_pre"],
                "pre_conversion_mean": pre_stats["mean_pre"],
                "pre_conversion_std": pre_stats["std_pre"],
            })

            log_print(
                f"For chromosome {chromosome }"
                f"📊 Pre-conversion OR: min={pre_stats['min_pre']:.6f}, "
                f"max={pre_stats['max_pre']:.6f}, mean={pre_stats['mean_pre']:.6f}"
            )

            # Convert OR → Beta
            log_print(f"🔄 Converting OR → log(OR) ... for chromosome {chromosome}")

            # df = df.with_columns(
            #     pl.col(effect_col).clip(lower_bound=1e-4).log().alias("beta")
            # )
            zero_count = df.filter(pl.col(effect_col) == 0.0).height
            total = df.height
            pct = (zero_count / total * 100) if total else 0

            log_print(
                f"⚠️ OR == 0 count: {zero_count:,} / {total:,} "
                f"({pct:.2f}%) in {effect_col} [chr{chromosome}] → set to None in beta"
            )
            df = df.with_columns(
                pl.when(pl.col(effect_col) > 0)
                .then(pl.col(effect_col).log())
                .otherwise(None)
                .alias("beta")
            ) 

            sample_column_dict["beta_col"] = "beta"
            qc_info["conversion"] = "OR_to_Beta_log_transform_applied"

            # Post-conversion stats
            post_stats = (
                df.select([
                    pl.col("beta").min().alias("min_post"),
                    pl.col("beta").max().alias("max_post"),
                    pl.col("beta").mean().alias("mean_post"),
                    pl.col("beta").std().alias("std_post"),
                ])
                .to_dicts()[0]
            )

            qc_info.update({
                "post_conversion_min": post_stats["min_post"],
                "post_conversion_max": post_stats["max_post"],
                "post_conversion_mean": post_stats["mean_post"],
                "post_conversion_std": post_stats["std_post"],
            })
            
            stats = (
                df.select([
                    pl.col(effect_col).min().alias("min_beta"),
                    pl.col(effect_col).max().alias("max_beta"),
                    pl.col(effect_col).mean().alias("mean_beta"),
                    pl.col(effect_col).std().alias("std_beta"),
                    pl.len().alias("total"),
                ])
                .to_dicts()[0]
            )

            qc_info.update({
                "min_beta": stats["min_beta"],
                "max_beta": stats["max_beta"],
                "mean_beta": stats["mean_beta"],
                "std_beta":  stats["std_beta"],
            })

            log_print(
                f"For chromosome {chromosome }"
                f"📈 Post-conversion Beta: min={post_stats['min_post']:.6f}, "
                f"max={post_stats['max_post']:.6f}, mean={post_stats['mean_post']:.6f}"
            )

        # -------------------------------------------------------
        # Finalize
        # -------------------------------------------------------
        qc_info["final_total"] = total
        log_print(f"🎯 Effect-size harmonization complete for chromosome {chromosome}. \n")

        # -------------------------------------------------------
        # Save logfile 
        # -------------------------------------------------------
        completed_successfully = True
        return df, qc_info, sample_column_dict
    except Exception as e:
        log_print(f"❌ ERROR chromosome {chromosome}: {str(e)}", indent=0)
        raise
    finally:
        with open(log_file, "w") as f:
            f.write(buffer.getvalue())

        if completed_successfully:
            print(buffer.getvalue(), end="")

# def is_beta_or_or(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, dict, dict]:
#     """
#     Determine whether the effect-size column contains Beta or Odds Ratio.
#     If OR, convert to Beta using log(OR). Generate QC summary and full logfile.

#     Logged version:
#         • Dedicated logfile per chromosome
#         • Multiprocessing-safe (no stdout)
#         • All logs written only to logs/*.log
#     """

#     # -------------------------------------------------------
#     # 1. Setup log buffer + file path
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_is_beta_or_or.log"
#     log_fh = open(log_file, "a")   # ✅ REQUIRED

#     def log_print(*args):
#         msg = " ".join(map(str, args))  
#         indent = " " * 16
#         final = f"{indent}{msg}"
#         print(final)                 # console
#         log_fh.write(final + "\n")   # file
#         log_fh.flush()
    
#     # -------------------------------------------------------
#     # 2. Main logic
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     log_print("\n🧩 Checking if effect column represents Beta or Odds Ratio (OR)...")

#     effect_col = sample_column_dict.get("beta_or_col", None)
#     if not effect_col or effect_col not in df.columns:
#         raise ValueError("❌ Effect size column ('beta_or_col') missing or not found in dataframe.")

#     # Cast to float safely
#     df = df.with_columns(pl.col(effect_col).cast(pl.Float64, strict=False))

#     total = df.height
#     negative_count = df.filter(pl.col(effect_col) < 0.0).height
#     negative_frac = negative_count / total if total else 0.0

#     qc_info.update({
#         "effect_col": effect_col,
#         "negative_value_count": negative_count,
#         "negative_fraction": round(negative_frac, 6),
#     })

#     # -------------------------------------------------------
#     # 3A. BETA detected
#     # -------------------------------------------------------
#     if negative_frac > 0.10:
#         log_print(f"✅ Detected as Beta: {negative_frac*100:.2f}% values are negative.")
#         print(f"                ✅ Effect size Detected as Beta for chromosome {chromosome}: {negative_frac*100:.2f}% values are negative.")
#         print(f" ")
#         sample_column_dict["effect_type"] = "beta"
#         sample_column_dict["beta_col"] = effect_col
#         qc_info["effect_type"] = "beta"

#         stats = (
#             df.select([
#                 pl.col(effect_col).min().alias("min_beta"),
#                 pl.col(effect_col).max().alias("max_beta"),
#                 pl.col(effect_col).mean().alias("mean_beta"),
#                 pl.col(effect_col).std().alias("std_beta"),
#                 pl.len().alias("total"),
#             ])
#             .to_dicts()[0]
#         )

#         qc_info.update({
#             "min_beta": stats["min_beta"],
#             "max_beta": stats["max_beta"],
#             "mean_beta": stats["mean_beta"],
#             "std_beta":  stats["std_beta"],
#         })

#         log_print(
#             f"📈 Beta summary: min={stats['min_beta']:.6f}, "
#             f"max={stats['max_beta']:.6f}, mean={stats['mean_beta']:.6f}"
#         )

#     # -------------------------------------------------------
#     # 3B. OR detected → convert
#     # -------------------------------------------------------
#     else:
#         log_print(f"⚠️ Detected as Odds Ratio: only {negative_frac*100:.2f}% are negative.")
#         print(f"                ⚠️ Effect size Detected as Odds Ratio for chromosome {chromosome}: only {negative_frac*100:.2f}% are negative.")

#         sample_column_dict["effect_type"] = "odds_ratio"
#         qc_info["effect_type"] = "odds_ratio"

#         # Pre-conversion stats
#         pre_stats = (
#             df.select([
#                 pl.col(effect_col).min().alias("min_pre"),
#                 pl.col(effect_col).max().alias("max_pre"),
#                 pl.col(effect_col).mean().alias("mean_pre"),
#                 pl.col(effect_col).std().alias("std_pre"),
#             ])
#             .to_dicts()[0]
#         )

#         qc_info.update({
#             "pre_conversion_min": pre_stats["min_pre"],
#             "pre_conversion_max": pre_stats["max_pre"],
#             "pre_conversion_mean": pre_stats["mean_pre"],
#             "pre_conversion_std": pre_stats["std_pre"],
#         })

#         log_print(
#             f"📊 Pre-conversion OR: min={pre_stats['min_pre']:.6f}, "
#             f"max={pre_stats['max_pre']:.6f}, mean={pre_stats['mean_pre']:.6f}"
#         )

#         # Convert OR → Beta
#         log_print("🔄 Converting OR → log(OR) ...")

#         df = df.with_columns(
#             pl.col(effect_col).clip(lower_bound=1e-4).log().alias("beta")
#         )

#         sample_column_dict["beta_col"] = "beta"
#         qc_info["conversion"] = "OR_to_Beta_log_transform_applied"

#         # Post-conversion stats
#         post_stats = (
#             df.select([
#                 pl.col("beta").min().alias("min_post"),
#                 pl.col("beta").max().alias("max_post"),
#                 pl.col("beta").mean().alias("mean_post"),
#                 pl.col("beta").std().alias("std_post"),
#             ])
#             .to_dicts()[0]
#         )

#         qc_info.update({
#             "post_conversion_min": post_stats["min_post"],
#             "post_conversion_max": post_stats["max_post"],
#             "post_conversion_mean": post_stats["mean_post"],
#             "post_conversion_std": post_stats["std_post"],
#         })
        
#         stats = (
#             df.select([
#                 pl.col(effect_col).min().alias("min_beta"),
#                 pl.col(effect_col).max().alias("max_beta"),
#                 pl.col(effect_col).mean().alias("mean_beta"),
#                 pl.col(effect_col).std().alias("std_beta"),
#                 pl.len().alias("total"),
#             ])
#             .to_dicts()[0]
#         )

#         qc_info.update({
#             "min_beta": stats["min_beta"],
#             "max_beta": stats["max_beta"],
#             "mean_beta": stats["mean_beta"],
#             "std_beta":  stats["std_beta"],
#         })

#         log_print(
#             f"📈 Post-conversion Beta: min={post_stats['min_post']:.6f}, "
#             f"max={post_stats['max_post']:.6f}, mean={post_stats['mean_post']:.6f}"
#         )

#     # -------------------------------------------------------
#     # Finalize
#     # -------------------------------------------------------
#     qc_info["final_total"] = total
#     log_print("🎯 Effect-size harmonization complete.\n")

#     # -------------------------------------------------------
#     # Save logfile
#     # -------------------------------------------------------
#     log_fh.close()
#     return df, qc_info, sample_column_dict
