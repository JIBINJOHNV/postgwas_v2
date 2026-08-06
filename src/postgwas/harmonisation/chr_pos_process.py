from pathlib import Path
import io
import os
import polars as pl
from typing import Dict, Tuple


def fix_chr_pos_allele_column(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    drop_mt: bool = True
):
    import io
    from pathlib import Path
    import polars as pl
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_fix_chr_pos_column.log"
    log_buffer = io.StringIO()
    
    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")
    try:
        def clean_col(x):
            if x is None:
                return None
            x = str(x).strip()
            if x.upper() in {"", "NA", "NAN", "NONE"}:
                return None
            return x
        chr_col = clean_col(sample_column_dict.get("chr_col"))
        pos_col = clean_col(sample_column_dict.get("pos_col"))
        chr_pos_col = clean_col(sample_column_dict.get("chr_pos_col"))
        log_print("\n🧬 Starting chromosome/position harmonization...")
        # ----------------------------
        # SPLIT LOGIC
        # ----------------------------
        if chr_pos_col and chr_pos_col in df.columns and (not chr_col or not pos_col):
            log_print(f"🔹 Chromosome {chr_col} and/or position columns {pos_col} not found. Attempting to split combined chromosome-position column {chr_pos_col}...")
            non_null_before = df.select(
                pl.col(chr_pos_col).is_not_null().sum()
            ).item()
            dup_before = df.select(
                pl.col(chr_pos_col).is_duplicated().sum()
            ).item()
            split_success = False
            for sep in [":", "-", "_"]:
                try:
                    temp = df.with_columns([
                        pl.col(chr_pos_col)
                          .str.split_exact(sep, 1)
                          .struct.field("field_0")
                          .alias("Chr"),
                        pl.col(chr_pos_col)
                          .str.split_exact(sep, 1)
                          .struct.field("field_1")
                          .alias("Pos")
                    ])
                    non_null_after = temp.select(
                        pl.col("Pos").is_not_null().sum()
                    ).item()
                    if non_null_after == non_null_before:
                        df = temp
                        split_success = True
                        log_print(f"✅ Split using '{sep}'")
                        break
                except Exception as e:
                    log_print(f"Delimiter '{sep}' failed: {e}")
            if not split_success:
                raise RuntimeError(
                    f"❌ Failed to split '{chr_pos_col}' using ':', '-', '_' "
                    f"(non-null before={non_null_before})"
                )
            # Clean Pos
            df = df.with_columns([
                pl.col("Pos").str.replace_all(r"\D", "").alias("Pos")
            ])
            chr_col, pos_col = "Chr", "Pos"
            sample_column_dict.update({
                "chr_col": "Chr",
                "pos_col": "Pos",
                "chr_pos_col": None
            })
            dup_after = df.select(
                pl.struct(["Chr", "Pos"]).is_duplicated().sum()
            ).item()
            if dup_before != dup_after:
                raise RuntimeError(
                    f"❌ Duplicate mismatch after split: before={dup_before}, after={dup_after}"
                )
        # ----------------------------
        # NORMALIZATION
        # ----------------------------
        df = df.with_columns([
            pl.col(chr_col).cast(pl.Utf8).str.strip_chars()
              .str.replace(r"(?i)^chr", "", literal=False)
              .str.replace(r"\.0$", "", literal=False)
              .str.to_uppercase()
              .alias(chr_col),
            pl.col(pos_col).cast(pl.Float64, strict=False)
                           .cast(pl.Int64, strict=False)
                           .alias(pos_col)
        ])
        # ----------------------------
        # Chromosome standardization
        # ----------------------------
        df = df.with_columns(
            pl.col(chr_col).replace({
                "23": "X", "24": "Y", "25": "MT", "M": "MT"
            }, default=pl.col(chr_col)).alias(chr_col)
        )
        # ----------------------------
        # Filtering
        # ----------------------------
        allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]
        if not drop_mt:
            allowed.append("MT")
        df = df.filter(pl.col(chr_col).is_not_null())
        if pos_col in df.columns:
            df = df.filter(pl.col(pos_col).is_not_null())
        df = df.filter(pl.col(chr_col).is_in(allowed))
        log_print(f"✅ Completed: {df.height:,} rows")
        return df, sample_column_dict
    except Exception as e:
        log_print(f"❌ ERROR: {e}")
        raise
    finally:
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())


# def fix_chr_pos_allele_column(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     drop_mt: bool = True
# ) -> pl.DataFrame:
#     """
#     Process chromosome and position columns in a GWAS DataFrame.
#     """
#     import io  # Ensure io is imported inside or at top
#     print("\n🧬 Starting chromosome/position harmonization...")

#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir = sample_column_dict.get("output_folder", ".")
#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)
#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_fix_chr_pos_column.log"

#     log_buffer = io.StringIO()
#     def log_print(*args):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n")

#     log_print("\n🧬 Starting chromosome/position harmonization...") 

#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     chr_pos_col = sample_column_dict.get("chr_pos_col")

#     # 1. Handle combined Chr:Pos column
#     if chr_pos_col and chr_pos_col in df.columns and (not chr_col or not pos_col):
#         log_print("🔹 Splitting combined chromosome-position column...")
#         print("🔹 Splitting combined chromosome-position column...")
#         df = df.with_columns([
#             pl.col(chr_pos_col).str.split_exact(":", 1).struct.field("field_0").alias("Chr"),
#             pl.col(chr_pos_col).str.split_exact(":", 1).struct.field("field_1")
#                   .str.replace_all(r"\D", "").alias("Pos") # Remove non-digits from Pos
#         ])
#         chr_col, pos_col = "Chr", "Pos"
#         sample_column_dict.update({"chr_col": "Chr", "pos_col": "Pos", "chr_pos_col": None})

#     # 2. Initial cleanup and normalization
#     # We combine the .0 strip and the string cleanup into one step
#     df = df.with_columns([
#         pl.col(chr_col).cast(pl.Utf8).str.strip_chars()
#           .str.replace(r"(?i)^chr", "", literal=False)
#           .str.replace(r"\.0$", "", literal=False)
#           .str.to_uppercase()
#           .alias(chr_col),
#         pl.col(pos_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).alias(pos_col)
#     ])

#     # 3. Standardize Chromosome naming (X, Y, MT)
#     # Using 'replace' is safer and faster than multiple 'when/then'
#     df = df.with_columns(
#         pl.col(chr_col).replace({
#             "23": "X", "24": "Y", "25": "MT", "M": "MT"
#         }, default=pl.col(chr_col)).alias(chr_col)
#     )

#     # 4. Filter allowed chromosomes
#     allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]
#     if not drop_mt:
#         allowed.append("MT")

#     # CRITICAL: Filter out nulls BEFORE the 'is_in' check to avoid logic gaps
#     df = df.filter(pl.col(chr_col).is_not_null())
#     if pos_col in df.columns:
#         df = df.filter(pl.col(pos_col).is_not_null())

#     df = df.filter(pl.col(chr_col).is_in(allowed))

#     log_print(f"🧩 [chr{chromosome}] CHR,POS, EA,OA Normalization completed — {df.height:,} rows remain.\n")

#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return df, sample_column_dict
