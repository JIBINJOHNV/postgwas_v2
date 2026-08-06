import polars as pl
import math
from typing import Dict, Tuple, Any, Union
import io, sys
from pathlib import Path

# -------------------------------------------------------
# GLOBAL BUFFER
# -------------------------------------------------------
buffer = io.StringIO()

def log_print(buffer, *args):
    msg = " ".join(str(x) for x in args)
    buffer.write(msg + "\n")

# =======================================================
# DETECT PVAL TYPE (UNCHANGED)
# =======================================================
def detect_pval_type(
    df: pl.DataFrame,
    sample_column_dict,
    proportion_threshold: float = 0.001):

    pval_col = sample_column_dict.get("pval_col")
    if not pval_col or pval_col not in df.columns:
        raise ValueError("Missing or invalid p-value column in sample_column_dict.")

    pvals = (
        df.select(pl.col(pval_col).cast(pl.Float64, strict=False))
        .drop_nulls()
        .rename({pval_col: "p"})
    )

    total = pvals.height
    if total == 0:
        raise ValueError("No numeric P-values found for detection.")

    over_1p2 = pvals.filter(pl.col("p") > 1.005).height / total
    within_0_1 = pvals.filter((pl.col("p") >= 0) & (pl.col("p") <= 1.005)).height / total

    max_val = float(pvals.select(pl.col("p").max()).item())
    median_val = float(pvals.select(pl.col("p").median()).item())

    if over_1p2 >= proportion_threshold:
        detected_type = "mlogp"
    elif max_val > 2:
        detected_type = "mlogp"
    elif within_0_1 >= 0.98:
        detected_type = "pvalue"
    else:
        detected_type = "pvalue"

    return {
        "detected_type": detected_type,
        "n_values": total,
        "over_1.2_fraction": over_1p2,
        "within_[0,1.05]_fraction": within_0_1,
        "median": median_val,
        "max": max_val,
    }


# =======================================================
# PVAL → MLOGP (NO PRINT)
# =======================================================
def convert_pval_to_mlogp(
    df: pl.DataFrame,
    sample_column_dict,
    buffer,
    output_col: str = "LP",
    min_p: float = 1e-300,
    max_p: float = 0.9999):

    try:
        pval_col = sample_column_dict["pval_col"]
        if pval_col not in df.columns:
            raise ValueError(f"Column '{pval_col}' not found in DataFrame.")

        df = df.filter(
            (pl.col(pval_col).is_not_null()) &
            (pl.col(pval_col) > 0) &
            (pl.col(pval_col) < 1)
        )

        df = df.with_columns(
            pl.col(pval_col).clip(min_p, max_p).alias("__p_clipped")
        )

        df = df.with_columns(
            (-pl.col("__p_clipped").log10()).clip(0.0, None).alias(output_col)
        ).drop("__p_clipped")

        stats = df.select(
            pl.col(output_col).min().alias("min_LP"),
            pl.col(output_col).max().alias("max_LP"),
            pl.col(output_col).mean().alias("mean_LP"),
            pl.col(output_col).median().alias("median_LP"),
        ).to_dicts()[0]

        log_print(buffer, f"\t\t\t\t✅ Converted raw p-values to -log10(p) → '{output_col}'")

        sample_column_dict["pval_col"] = output_col

        return df, stats, sample_column_dict

    except Exception as e:
        log_print(buffer, f"\t\t\t\t❌ ERROR convert_pval_to_mlogp: {str(e)}")
        raise


# =======================================================
# MLOGP → PVAL (NO PRINT)
# =======================================================
# =======================================================
# MLOGP → PVAL (NO PRINT)
# =======================================================
def convert_mlogp_to_pval(
    df: pl.DataFrame,
    sample_column_dict,
    buffer,
    output_col: str = "PVAL",
    min_lp: float = 0.0,
    max_lp: float = 300.0
):

    try:
        pval_col = sample_column_dict["pval_col"]

        if pval_col not in df.columns:
            raise ValueError(f"Column '{pval_col}' not found in DataFrame.")

        n_clipped = df.filter(
            (pl.col(pval_col) < min_lp) | (pl.col(pval_col) > max_lp)
        ).height

        df = df.with_columns(
            pl.col(pval_col).clip(min_lp, max_lp).alias("__lp_clipped")
        )

        if n_clipped > 0:
            log_print(
                buffer,
                f"\t\t\t\t⚠️ Clipped {n_clipped} LP values to [{min_lp}, {max_lp}]"
            )

        df = df.with_columns(
            (10 ** (-pl.col("__lp_clipped"))).alias("__raw_p")
        )

        df = df.with_columns(
            pl.when(pl.col("__raw_p") == 0)
            .then(1e-300)
            .otherwise(pl.col("__raw_p"))
            .alias(output_col)
        ).drop(["__lp_clipped", "__raw_p"])

        # Compute summary stats after conversion
        mlogp_stats = df.select(
            pl.col(output_col).min().alias("min_PVAL"),
            pl.col(output_col).max().alias("max_PVAL"),
            pl.col(output_col).mean().alias("mean_PVAL"),
            pl.col(output_col).median().alias("median_PVAL"),
        ).to_dicts()[0]

        # Add clipping count without overwriting summary stats
        mlogp_stats["mlogp_values_total_clipped"] = n_clipped

        log_print(
            buffer,
            f"\t\t\t\t✅ Converted -log10(p) → raw p-value → '{output_col}'"
        )

        log_print(
            buffer,
            f"\t\t\t\tPVAL summary: "
            f"min={mlogp_stats['min_PVAL']:.3e}, "
            f"max={mlogp_stats['max_PVAL']:.3e}, "
            f"mean={mlogp_stats['mean_PVAL']:.3e}, "
            f"median={mlogp_stats['median_PVAL']:.3e}"
        )

        sample_column_dict["pval_col"] = output_col

        return df, mlogp_stats, sample_column_dict

    except Exception as e:
        log_print(buffer, f"\t\t\t\t❌ ERROR convert_mlogp_to_pval: {str(e)}")
        raise


# =======================================================
# MAIN FUNCTION
# =======================================================
def detect_and_convert_pval(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "LP",
    proportion_threshold: float = 0.001,
    max_allowed_lp: float = 300.0,
):

    buffer = io.StringIO()   # ✅ LOCAL buffer (multiprocessing-safe)

    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    logfile = log_dir / f"{gwas_outputname}_chr{chromosome}_detect_pval.log"

    status_msg = ""

    try:
        log_print(buffer, f"📊 [Chr {chromosome}] Starting P-value detection...")

        pval_col = sample_column_dict["pval_col"]

        df = df.with_columns(pl.col(pval_col).cast(pl.Float64, strict=False))

        # ==================================================
        # INITIAL QC (suffix retained exactly)
        # ==================================================
        initial_stats = df.select([
            pl.len().alias("initial_total_variants_with_pvalues"),
            pl.col(pval_col).is_null().sum().alias("initial_variants_with_null_pvalues"),
            (pl.col(pval_col) == 0).sum().alias("initial_variants_with_zero_pvalues"),
            (pl.col(pval_col) < 0).sum().alias("initial_variants_with_lt0_pvalues"),
            (pl.col(pval_col) > 1).sum().alias("initial_variants_with_gt1_pvalues"),
            (pl.col(pval_col) > 1.05).sum().alias("initial_variants_with_gt1_05_pvalues"),
            pl.col(pval_col).min().alias("initial_min_pvalues"),
            pl.col(pval_col).max().alias("initial_max_pvalues"),
            pl.col(pval_col).mean().alias("initial_mean_pvalues"),
            pl.col(pval_col).median().alias("initial_median_pvalues"),
        ]).to_dicts()[0]

        qc_info = initial_stats.copy()
        log_print(buffer, f"\t\t\t\t📊 Initial P-value stats:")
        log_print(buffer, f"\t\t\t\t\tTotal: {initial_stats['initial_total_variants_with_pvalues']:,}")
        log_print(buffer, f"\t\t\t\t\tNull: {initial_stats['initial_variants_with_null_pvalues']:,}")
        log_print(buffer, f"\t\t\t\t\tZero: {initial_stats['initial_variants_with_zero_pvalues']:,}")
        log_print(buffer, f"\t\t\t\t\t<0: {initial_stats['initial_variants_with_lt0_pvalues']:,}")
        log_print(buffer, f"\t\t\t\t\t>1: {initial_stats['initial_variants_with_gt1_pvalues']:,}")
        log_print(buffer, f"\t\t\t\t\t>1.05: {initial_stats['initial_variants_with_gt1_05_pvalues']:,}")

        # --------------------------------------------------
        # Detection
        # --------------------------------------------------
        detection_info = detect_pval_type(df, sample_column_dict, proportion_threshold)
        qc_info.update(detection_info)
        log_print(buffer, f"\t\t\t\t✅ P-value type detected: {detection_info['detected_type']} for chromosome {chromosome}")
        # ============================================================
        # RAW PVAL QC (fixed order)
        # ============================================================
        if detection_info["detected_type"] != "mlogp":

            total_before = df.height

            n_null = df.select(pl.col(pval_col).is_null().sum()).item()
            n_lt0  = df.select((pl.col(pval_col) < 0).sum()).item()
            n_gt1_05 = df.select((pl.col(pval_col) > 1.05).sum()).item()

            df = df.filter(
                (pl.col(pval_col).is_not_null()) &
                (pl.col(pval_col) >= 0) &
                (pl.col(pval_col) <= 1.05)
            )

            total_after_filter = df.height

            n_zero_after = df.select((pl.col(pval_col) == 0).sum()).item()

            df = df.with_columns(
                pl.when(pl.col(pval_col) == 0)
                .then(1e-300)
                .otherwise(pl.col(pval_col))
                .alias(pval_col)
            )

            df = df.with_columns(
                pl.col(pval_col).clip(1e-300, 0.999999).alias(pval_col)
            )

            df = df.with_columns(pl.col(pval_col).alias(output_col))

            sample_column_dict["pval_col"] = output_col
            pval_col = sample_column_dict["pval_col"]

            # ---------------- LOG ----------------
            if n_null > 0:
                log_print(buffer, f"\t\t\t\t⚠️ Removed {n_null} variants_with_null_pvalues")

            if n_lt0 > 0:
                log_print(buffer, f"\t\t\t\t⚠️ Removed {n_lt0} variants_with_lt0_pvalues")

            if n_gt1_05 > 0:
                log_print(buffer, f"\t\t\t\t⚠️ Removed {n_gt1_05} variants_with_gt1_05_pvalues")

            if n_zero_after > 0:
                log_print(buffer, f"\t\t\t\t⚠️ Replaced {n_zero_after} variants_with_zero_pvalues")

            if total_before != total_after_filter:
                log_print(buffer, f"\t\t\t\t📊 Raw_pvalues_QC: before filtering {total_before} → after filtering {total_after_filter}")

            qc_info.update({
                "raw_total_variants_with_pvalues": total_before,
                "after_filter_variants_with_pvalues": total_after_filter,
                "variants_with_null_pvalues_removed": n_null,
                "variants_with_lt0_pvalues_removed": n_lt0,
                "variants_with_gt1_05_pvalues_removed": n_gt1_05,
                "variants_with_zero_pvalues_replaced": n_zero_after,
            })

        else:
            total_before = df.height

            df, mlogp_stats, sample_column_dict = convert_mlogp_to_pval(
                df=df,
                sample_column_dict=sample_column_dict,
                buffer=buffer,
                output_col="PVAL",
                min_lp=0.0,
                max_lp=max_allowed_lp,
            )

            pval_col = sample_column_dict["pval_col"]

            total_after_filter = df.height

            qc_info.update({
                "mlogp_total_variants_with_pvalues": total_before,
                "after_mlogp_to_pval_variants_with_pvalues": total_after_filter,
                **mlogp_stats,
            })

        pval_col = sample_column_dict["pval_col"]

        final_stats = df.select([
            pl.len().alias("final_total_variants_with_pvalues"),
            pl.col(pval_col).is_null().sum().alias("final_variants_with_null_pvalues"),
            (pl.col(pval_col) == 0).sum().alias("final_variants_with_zero_pvalues"),
            (pl.col(pval_col) < 0).sum().alias("final_variants_with_lt0_pvalues"),
            (pl.col(pval_col) > 1).sum().alias("final_variants_with_gt1_pvalues"),
            (pl.col(pval_col) > 1.05).sum().alias("final_variants_with_gt1_05_pvalues"),
            pl.col(pval_col).min().alias("final_min_pvalues"),
            pl.col(pval_col).max().alias("final_max_pvalues"),
            pl.col(pval_col).mean().alias("final_mean_pvalues"),
            pl.col(pval_col).median().alias("final_median_pvalues"),
        ]).to_dicts()[0]

        qc_info.update(final_stats)
        status_msg = f"✅ Completed chromosome {chromosome}"

        return df, qc_info, sample_column_dict

    except Exception as e:
        log_print(buffer, f"\t\t\t\t❌ ERROR chromosome {chromosome}: {str(e)}")
        status_msg = f"❌ FAILED chromosome {chromosome}: {str(e)}"
        raise

    finally:
        with open(logfile, "w") as f:
            f.write(buffer.getvalue())

        print(buffer.getvalue(), end="")