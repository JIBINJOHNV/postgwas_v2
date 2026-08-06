import polars as pl


def genome_build(t_df: pl.DataFrame, grch37_file: str, grch38_file: str, sample_column_dict: dict):
    """
    Infer genome build (GRCh37 or GRCh38) by matching variant positions and alleles.
    Parameters
    ----------
    t_df : pl.DataFrame
        Subset of GWAS summary statistics.
    grch37_file : str
        Path to GRCh37 reference variant file.
    grch38_file : str
        Path to GRCh38 reference variant file.
    sample_column_dict : dict
        Contains mapping for 'chr_col', 'pos_col', 'ea_col', 'oa_col'.
    Returns
    -------
    dict
        genome_build_info with counts, percentages, and inferred build.
    """
    # --- Load reference datasets ---
    grch37_df = (
        pl.read_csv(grch37_file, separator="\t", schema_overrides={"CHROM": pl.Utf8})
        .select(["CHROM", "POS", "REF", "ALT"])
        .with_columns([
            pl.col("CHROM")
            .cast(pl.Utf8)
            .str.replace_all("^0", ""),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
        ])
        .filter(
            ~pl.col("CHROM").str.to_lowercase().is_in(
                ["x", "y", "mt", "23", "24"]
            )
        )
    )

    grch38_df = (
        pl.read_csv(grch38_file, separator="\t", schema_overrides={"CHROM": pl.Utf8})
        .select(["CHROM", "POS", "REF", "ALT"])
        .with_columns([
            pl.col("CHROM").cast(pl.Utf8).str.replace_all("^0", ""),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
        ])
        .filter(
            ~pl.col("CHROM").str.to_lowercase().is_in(
                ["x", "y", "mt", "23", "24"]
            )
        )   
    )
    # --- Prepare GWAS dataframe ---
    chr_col = sample_column_dict["chr_col"]
    pos_col = sample_column_dict["pos_col"]
    ea_col  = sample_column_dict["ea_col"]
    oa_col  = sample_column_dict["oa_col"]
    
    t_df = (
        t_df.with_columns([
            pl.col(chr_col)
            .cast(pl.Utf8)
            .str.replace_all("(?i)^chr", "")  # remove chr / CHR
            .str.replace_all("^0", ""),       # remove leading 0
            pl.col(pos_col).cast(pl.Int64),
            pl.col(ea_col).cast(pl.Utf8).str.to_uppercase(),
            pl.col(oa_col).cast(pl.Utf8).str.to_uppercase(),
        ])
    )
    # --- Helper to count matches ---
    def count_matches(ref_df: pl.DataFrame) -> int:
        same = t_df.join(
            ref_df,
            left_on=[chr_col, pos_col, ea_col, oa_col],
            right_on=["CHROM", "POS", "REF", "ALT"],
            how="inner"
        ).height
        
        flip = t_df.join(
            ref_df,
            left_on=[chr_col, pos_col, oa_col, ea_col],
            right_on=["CHROM", "POS", "REF", "ALT"],
            how="inner"
        ).height
        return same + flip
    # --- Compute match statistics ---
    grch37_match = count_matches(grch37_df)
    grch38_match = count_matches(grch38_df)
    total = max(t_df.height, 1)
    perc37 = grch37_match / total
    perc38 = grch38_match / total

    w38 = perc38 / (perc38 + perc37)
    w37 = perc37 / (perc38 + perc37)

    # Normalized dominance weights
    denom = perc37 + perc38
    # Guard: no informative matches at all
    if denom == 0:
        inferred_build = "Ambiguous"
    else:
        w37 = perc37 / denom
        w38 = perc38 / denom
        if w38 >= 0.9:
            inferred_build = "GRCh38"
        elif w37 >= 0.9:
            inferred_build = "GRCh37"
        else:
            inferred_build = "Ambiguous"
    # --- Summary output ---
    genome_build_info = {
        "inferred_build": inferred_build,
        "grch37_matches": grch37_match,
        "grch38_matches": grch38_match,
        "percent_grch37": round(perc37 * 100, 2),
        "percent_grch38": round(perc38 * 100, 2),
        "total_variants": total
    }
    #print(f"📊 GRCh37: {grch37_match:,} ({perc37*100:.2f}%) | GRCh38: {grch38_match:,} ({perc38*100:.2f}%)")
    return genome_build_info
