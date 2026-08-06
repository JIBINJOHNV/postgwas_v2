import os
import polars as pl
from pathlib import Path
from typing import Tuple, Dict, Optional

def add_or_calculate_eaf(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    eaffile: str = "NA",
    default_eaf_file: str = "NA",
    default_eaf_eafcolumn: str = "EAF",
    maf_eaf_decision_cutoff: float = 0.95,
    external_eaf_colmap: Optional[Dict] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Harmonise EAF values using internal or external reference.
    STRONG FAIL-FAST VERSION.
    Writes logs IMMEDIATELY to disk to prevent data loss on crash.
    """

    # -------------------------------------------------------
    # 1. Setup Logging (Immediate Write Mode)
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_eaf.log"
    
    # Initialize file (overwrite old log)
    with open(log_file, "w") as f:
        f.write(f"--- Log started for {gwas_outputname} Chr {chromosome} ---\n")

    # Helper: Writes to disk immediately
    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        try:
            with open(log_file, "a") as f:
                f.write(msg + "\n")
        except Exception as e:
            print(f"CRITICAL LOGGING ERROR: {e}")

    qc_info = {"initial_variants": df.height}
    log_print("\n🧩 Starting allele frequency harmonization...")

    # -------------------------------------------------------
    # 2. Get Column Names
    # -------------------------------------------------------
    try:
        chr_col = sample_column_dict["chr_col"]
        pos_col = sample_column_dict["pos_col"]
        ea_col = sample_column_dict["ea_col"]
        oa_col = sample_column_dict["oa_col"]
        eaf_col = sample_column_dict.get("eaf_col", "NA")
        eafcolumn = sample_column_dict.get("eafcolumn", "EAF")
    except KeyError as e:
        log_print(f"❌ CRITICAL: Missing required column mapping in config: {e}")
        raise SystemExit(1)

    # -------------------------------------------------------
    # Helper — validate EAF column quality
    # -------------------------------------------------------
    def validate_eaf_column(
        df_check: pl.DataFrame,
        colname: str,
        label: str,
        maf_eaf_decision_cutoff=maf_eaf_decision_cutoff,
    ):
        if colname not in df_check.columns:
            msg = f"❌ {label} '{colname}' not found in dataframe for chromosome {chromosome}."
            log_print(msg)
            raise ValueError(msg)

        total = df_check.height
        out_of_range = df_check.filter(
            (pl.col(colname) < 0.0) | (pl.col(colname) > 1.05)
        ).height
        out_of_range_pct = out_of_range / total if total else 0

        low_freq = df_check.filter(pl.col(colname) <= 0.5).height
        low_freq_pct = low_freq / total if total else 0

        if out_of_range_pct > 0.05:
            msg = (f"❌ {label} '{colname}' contains {out_of_range_pct*100:.2f}% invalid values "
                   f"on chromosome {chromosome} → Pipeline stopped.")
            log_print(msg)
            raise ValueError(msg)

        if low_freq_pct > maf_eaf_decision_cutoff:
            msg = (f"🛑 CRITICAL ERROR: {label} '{colname}' behaves like MAF on chromosome {chromosome}.\n"
                   f"→ {low_freq_pct*100:.2f}% of values ≤ 0.5.\n"
                   f"→ This indicates the file provides MAF instead of EAF.\n"
                   f"Pipeline cannot continue.")
            log_print(msg)
            raise ValueError(msg)

        df_checked = df_check.with_columns(
            pl.col(colname).clip(0.000001, 1.0)
        )

        return df_checked, {
            "out_of_range_pct": out_of_range_pct,
            "low_freq_pct": low_freq_pct,
        }

    # -------------------------------------------------------
    # Helper — Load EXTERNAL EAF file
    # -------------------------------------------------------
    def load_eaf_file(
        path: str,
        eafcolumn: str,
        external_eaf_colmap: Optional[dict] = None,
        maf_eaf_decision_cutoff: float = maf_eaf_decision_cutoff,
    ) -> pl.DataFrame:

        if external_eaf_colmap is None:
            external_eaf_colmap = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}

        if path == "NA" or not os.path.exists(path):
            msg = f"❌ External EAF file '{path}' not found for chromosome {chromosome}."
            log_print(msg)
            raise FileNotFoundError(msg)

        log_print(f"📂 Loading external EAF: {os.path.basename(path)}")

        try:
            df_ext = pl.read_csv(path, separator="\t", has_header=True)
            if eafcolumn not in df_ext.columns: raise ValueError
            log_print("✔ Detected TAB-delimited format")
        except:
            try:
                df_ext = pl.read_csv(path, separator=" ", has_header=True)
                if eafcolumn not in df_ext.columns: raise ValueError
                log_print("✔ Detected SPACE-delimited format")
            except:
                log_print("⚠ Falling back to whitespace normalization via pandas")
                import pandas as pd
                df_ext = pl.from_pandas(pd.read_csv(path, delim_whitespace=True))

        for key in ["chr", "pos", "a1", "a2"]:
            if key not in external_eaf_colmap:
                msg = f"❌ Missing '{key}' in external_eaf_colmap."
                log_print(msg)
                raise ValueError(msg)
            if external_eaf_colmap[key] not in df_ext.columns:
                msg = f"❌ Expected column '{external_eaf_colmap[key]}' not found in external file."
                log_print(msg)
                raise ValueError(msg)

        chr_in = external_eaf_colmap["chr"]
        pos_in = external_eaf_colmap["pos"]
        a1_in = external_eaf_colmap["a1"]
        a2_in = external_eaf_colmap["a2"]

        log_print(f"✔ Using external EAF column map: {external_eaf_colmap}")

        df_ext = (
            df_ext.select([chr_in, pos_in, a1_in, a2_in, eafcolumn])
            .rename({chr_in: chr_col, pos_in: pos_col, a1_in: ea_col, a2_in: oa_col})
            .with_columns([
                pl.col(chr_col).str.replace("0X", "X"),
                pl.col(pos_col).cast(pl.Int64),
                pl.col(ea_col).str.to_uppercase(),
                pl.col(oa_col).str.to_uppercase(),
            ])
        )

        df_ext, _ = validate_eaf_column(df_ext, eafcolumn, "External EAF")
        return df_ext

    # -------------------------------------------------------
    # STEP 1 — Try internal EAF
    # -------------------------------------------------------
    eaf_valid = False

    if eaf_col != "NA" and eaf_col in df.columns:
        try:
            df, result = validate_eaf_column(df, eaf_col, "Internal EAF")
            qc_info.update(result)
            eaf_valid = True
            log_print("✔ Valid Internal EAF column found.")
        except ValueError as e:
            # We catch it here to log it, then exit
            log_print(f"🛑 Error with internal EAF: {e}")
            raise SystemExit(1)

    # -------------------------------------------------------
    # STEP 2 — Internal OK → compute zMAF
    # -------------------------------------------------------
    if eaf_valid:
        df_checked = df.with_columns(
            pl.when(pl.col(eaf_col) <= 0.5)
            .then(pl.col(eaf_col))
            .otherwise(1 - pl.col(eaf_col))
            .alias("zmaf")
        )
        sample_column_dict["eaf_col"] = eaf_col

    else:
        # ---------------------------------------------------
        # STEP 3 — Use external EAF (strict fail)
        # ---------------------------------------------------
        eaf_file_to_use = (
            eaffile if (eaffile != "NA" and os.path.exists(eaffile))
            else default_eaf_file
        )

        if eaffile != "NA":
            if os.path.exists(eaffile):
                eaf_file_to_use = eaffile
            elif os.path.exists(f"{eaffile}_chr{chromosome}.tsv.gz"):
                eaf_file_to_use = f"{eaffile}_chr{chromosome}.tsv.gz"
            else:
                msg = (f"\n🛑 PIPELINE STOPPED: Effect allele frequency (EAF) info missing.\n"
                       f"• No external EAF file was found for chromosome {chromosome}.\n"
                       f"• Expected: {eaffile}")
                log_print(msg)
                raise SystemExit(1)

        eafcolumn = (
            default_eaf_eafcolumn
            if eaf_file_to_use == default_eaf_file
            else sample_column_dict.get("eafcolumn", "EAF")
        )

        log_print(f"ℹ Using external EAF file: {eaf_file_to_use}")

        try:
            eaf_df = load_eaf_file(
                eaf_file_to_use,
                eafcolumn=eafcolumn,
                external_eaf_colmap=external_eaf_colmap,
                maf_eaf_decision_cutoff=maf_eaf_decision_cutoff,
            )
        except Exception as e:
            log_print(f"\n🛑 PIPELINE STOPPED: External EAF error on chromosome {chromosome}.\n{e}")
            raise SystemExit(1)

        df = df.with_columns([
            pl.col(ea_col).str.to_uppercase(),
            pl.col(oa_col).str.to_uppercase(),
        ])

        df1 = df.join(eaf_df, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

        eaf_flip = (
            eaf_df.rename({ea_col: oa_col, oa_col: ea_col})
            .with_columns((1 - pl.col(eafcolumn)).alias(eafcolumn))
        )
        df2 = df.join(eaf_flip, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

        df_checked = (
            pl.concat([df1, df2], how="vertical")
            .unique([chr_col, pos_col, ea_col, oa_col])
            .with_columns(
                pl.when(pl.col(eafcolumn) <= 0.5)
                .then(pl.col(eafcolumn))
                .otherwise(1 - pl.col(eafcolumn))
                .alias("zmaf")
            )
        )

        sample_column_dict["eaf_col"] = eafcolumn
        qc_info["external_eaf_file"] = eaf_file_to_use

    # -------------------------------------------------------
    # Final sanity clipping
    # -------------------------------------------------------
    eaf_used = sample_column_dict["eaf_col"]
    df_checked = df_checked.with_columns(
        pl.col(eaf_used).clip(0.000001, 1.0)
    )

    log_print(f"✅ EAF Harmonization Complete for Chr {chromosome}")
    
    return df_checked, qc_info, sample_column_dict