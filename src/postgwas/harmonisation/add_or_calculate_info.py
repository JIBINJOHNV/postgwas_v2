import polars as pl
import os
from pathlib import Path
from typing import Tuple, Dict


def add_or_calculate_info(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    info_file: str = "NA",
    info_column: str = "info",
    default_info_file: str = "NA",
    default_info_column: str = "info",
    log_dir: str = ".",
    output_folder: str = ".",
    output_prefix: str = "gwas",
    enable_stdout=True,  # kept but ignored for controlled logging
) -> Tuple[pl.DataFrame, Dict, dict]:

    import io

    # ============================================================
    # SETUP
    # ============================================================
    output_folder = Path(output_folder).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    log_dir = Path(log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{output_prefix}_chr{chromosome}_info.log"

    # ------------------------------------------------------------
    # BUFFERED LOGGER
    # ------------------------------------------------------------
    log_buffer = io.StringIO()
    PREFIX = " " * 16

    def log_print(*args, prefix=True):
        msg = " ".join(map(str, args))
        if prefix:
            msg = PREFIX + msg
        log_buffer.write(msg + "\n")  # ✅ buffer only (no print)

    try:
        log_print(" " * 60, prefix=False)
        log_print(f"📊 [Chr {chromosome}] Starting INFO harmonisation...", prefix=False) 

        qc_info = {"initial_variants": df.height}

        chr_col = sample_column_dict["chr_col"]
        pos_col = sample_column_dict["pos_col"]
        ea_col = sample_column_dict["ea_col"]
        oa_col = sample_column_dict["oa_col"]
        imp_info_col = sample_column_dict.get("imp_info_col", "NA")

        # ============================================================
        # NORMALIZE INPUT
        # ============================================================
        df = df.with_columns([
            # Chr as String (cleaned)
            pl.col(chr_col).cast(pl.String).str.replace("^chr", "", literal=False).str.strip_chars(),
            # Pos as Int64 (strip whitespace first if it's currently a string)
            pl.col(pos_col).cast(pl.String).str.strip_chars().cast(pl.Int64, strict=False),
            # Alleles as String
            pl.col(ea_col).cast(pl.String).str.to_uppercase().str.strip_chars(),
            pl.col(oa_col).cast(pl.String).str.to_uppercase().str.strip_chars()
        ])

        # ============================================================
        # STEP 1: INTERNAL INFO
        # ============================================================
        if imp_info_col != "NA" and imp_info_col in df.columns:

            log_print(f"🔍 Using internal INFO column: {imp_info_col} ; for chr{chromosome}")

            df = df.with_columns(
                pl.col(imp_info_col).cast(pl.Float64, strict=False)
            )

            info_col = imp_info_col

        else:
            potential_chr_file = f"{info_file}_chr{chromosome}.tsv.gz"

            if info_file != "NA" and os.path.exists(potential_chr_file):
                info_file_to_use = potential_chr_file
                info_col_to_use = info_column
            elif info_file != "NA" and os.path.exists(info_file):
                info_file_to_use = info_file
                info_col_to_use = info_column
            elif default_info_file != "NA" and os.path.exists(default_info_file):
                info_file_to_use = default_info_file
                info_col_to_use = default_info_column
            else:
                log_print(f"❌ No INFO file found for chromosome {chromosome}")
                raise FileNotFoundError(f"No INFO file found for chromosome {chromosome}")

            log_print(f"📂 Using INFO file: {info_file_to_use} ; for chr{chromosome}")

            # ============================================================
            # STEP 2: LOAD INFO (External File)
            # ============================================================
            info_df = (
                pl.read_csv(
                    info_file_to_use,
                    separator="\t",
                    null_values=["NA", ".", "-"],
                    infer_schema_length=10000 # Help Polars guess types correctly
                )
                .select(["CHROM", "POS", "ALT", "REF", info_col_to_use])
                .with_columns([
                    # CHROM as String
                    pl.col("CHROM").cast(pl.String).str.replace("^chr", "", literal=False).str.strip_chars(),
                    # POS as Int64
                    pl.col("POS").cast(pl.String).str.strip_chars().cast(pl.Int64, strict=False),
                    # Alleles as String
                    pl.col("ALT").cast(pl.String).str.to_uppercase().str.strip_chars(),
                    pl.col("REF").cast(pl.String).str.to_uppercase().str.strip_chars(),
                    # Info Score as Float
                    pl.col(info_col_to_use).cast(pl.Float64, strict=False)
                ])
                .rename({
                    "CHROM": chr_col,
                    "POS": pos_col,
                    "ALT": ea_col,
                    "REF": oa_col
                })
            )

            log_print(f"📊 INFO variants loaded: {info_df.height:,}")

            df_direct = df.join(
                info_df,
                on=[chr_col, pos_col, ea_col, oa_col],
                how="left"
            )

            info_df_flip = info_df.rename({ea_col: oa_col, oa_col: ea_col})

            df_flip = df.join(
                info_df_flip,
                on=[chr_col, pos_col, ea_col, oa_col],
                how="left"
            )

            flip_col = f"{info_col_to_use}_flip"
            df_flip = df_flip.rename({info_col_to_use: flip_col})

            df_merged = df_direct.join(
                df_flip.select([chr_col, pos_col, ea_col, oa_col, flip_col]),
                on=[chr_col, pos_col, ea_col, oa_col],
                how="left"
            )

            df = df_merged.with_columns(
                pl.when(pl.col(info_col_to_use).is_null())
                .then(pl.col(flip_col))
                .otherwise(pl.col(info_col_to_use))
                .alias(info_col_to_use)
            ).drop(flip_col)

            info_col = info_col_to_use
            sample_column_dict["imp_info_col"] = info_col

        # ============================================================
        # QC METRICS
        # ============================================================
        total = df.height
        missing = df.filter(pl.col(info_col).is_null()).height
        lt_zero = df.filter(pl.col(info_col) < 0).height
        gt_one = df.filter(pl.col(info_col) > 1).height
        gt_1_05 = df.filter(pl.col(info_col) > 1.05).height

        invalid = df.filter(
            (pl.col(info_col) < 0) |
            (pl.col(info_col) > 1)
        ).height

        low_info = df.filter(
            (pl.col(info_col).is_not_null()) &
            (pl.col(info_col) <= 0.3)
        ).height

        good_info = df.filter(
            (pl.col(info_col) > 0.3) &
            (pl.col(info_col) <= 1.0)
        ).height

        stats = df.select([
            pl.col(info_col).min().alias("min"),
            pl.col(info_col).max().alias("max"),
            pl.col(info_col).mean().alias("mean"),
            pl.col(info_col).median().alias("median"),
        ]).to_dicts()[0]

        invalid_df = df.with_columns([
            pl.when(pl.col(info_col).is_null()).then(pl.lit("MISSING"))
            .when(pl.col(info_col) < 0).then(pl.lit("INFO<0"))
            .when(pl.col(info_col) > 1.05).then(pl.lit("INFO>1.05"))
            .when((pl.col(info_col) > 1) & (pl.col(info_col) <= 1.05)).then(pl.lit("1<INFO<=1.05"))
            .otherwise(pl.lit(None))
            .alias("INFO_QC_REASON")
        ]).filter(pl.col("INFO_QC_REASON").is_not_null())

        invalid_path = output_folder / f"{output_prefix}_chr{chromosome}_invalid_info.tsv"

        # ============================================================
        # LOG SUMMARY
        # ============================================================
        log_print("=" * 60)
        log_print(f"📊 INFO QC SUMMARY FOR Chr{chromosome}")
        log_print(f"Total variants        : {total:,}")
        log_print(f"Missing INFO          : {missing:,}")
        log_print(f"INFO < 0              : {lt_zero:,}")
        log_print(f"INFO > 1              : {gt_one:,}")
        log_print(f"INFO > 1.05           : {gt_1_05:,}")
        log_print(f"Invalid (total)       : {invalid:,}")
        log_print(f"Low INFO [0,0.3]      : {low_info:,}")
        log_print(f"Good INFO (0.3,1]     : {good_info:,}")
        log_print(f"Range                 : {stats['min']} → {stats['max']}")

        qc_info["total_variant_infile"] = total
        qc_info["missing_info"] = missing
        qc_info["lt_zero"] = lt_zero
        qc_info["gt_one"] = gt_one
        qc_info["gt_1_05"] = gt_1_05
        qc_info["invalid"] = invalid
        qc_info["low_info"] = low_info
        qc_info["good_info"] = good_info
        qc_info["min_info"] = stats["min"]
        qc_info["max_info"] = stats["max"]
        qc_info["mean_info"] = stats["mean"]
        qc_info["median_info"] = stats["median"]

        if invalid_df.height > 0:
            invalid_df.write_csv(invalid_path, separator="\t")
            log_print(f"⚠️ Invalid INFO saved → {invalid_path}.gz ({invalid_df.height:,})")
            os.system(f"gzip --force  {invalid_path} &")
        else:
            log_print("✅ No invalid INFO variants found")

        log_print("=" * 60)
        log_print("")

    except Exception as e:
        log_print(f"❌ ERROR: {str(e)}")
        raise

    finally:
        # ============================================================
        # ALWAYS write + print logs
        # ============================================================
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        print(log_buffer.getvalue(), flush=True)

    return df, qc_info, sample_column_dict