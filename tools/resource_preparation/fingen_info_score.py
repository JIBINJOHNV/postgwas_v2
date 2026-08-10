import polars as pl
import os

input_file = "finnge_R12_annotated_variants_v1.gz"

# 1. Get schema with better inference
full_schema = pl.scan_csv(
    input_file,
    separator="\t",
    null_values="NA",
    infer_schema_length=10000,
).collect_schema()
all_columns = full_schema.names()

# 2. Identify needed columns
info_cols = [col for col in all_columns if col.startswith("INFO_")]
coord_cols = ["chr", "pos", "ref", "alt"]
required_cols = coord_cols + info_cols

print(f"Plan: Reading {len(required_cols)} columns, ignoring {len(all_columns) - len(required_cols)} columns.")

# 3. Force all INFO columns to Float64
required_schema = {}
for col in coord_cols:
    required_schema[col] = full_schema[col]

for col in info_cols:
    required_schema[col] = pl.Float64

# 4. Read only the required columns
print("Reading data... This may take a few minutes.")
df = pl.read_csv(
    input_file,
    separator="\t",
    columns=required_cols,
    schema_overrides=required_schema,
    null_values="NA",
    infer_schema_length=10000,
    low_memory=True,
)
print(f"Loaded {df.height} rows x {df.width} columns.")

# Replace chr 23 with "X" and cast entire column to string
df = df.with_columns(
    pl.when(pl.col("chr") == 23)
      .then(pl.lit("X"))
      .otherwise(pl.col("chr").cast(pl.Utf8))
      .alias("chr")
)

# 5. Compute horizontal stats
print("Computing Mean and Median across INFO columns...")
df = df.with_columns([
    pl.mean_horizontal(info_cols).alias("INFO_MEAN"),
    pl.concat_list(info_cols).list.eval(pl.element().median()).list.first().alias("INFO_MEDIAN"),
])

# 6. Process per chromosome
chroms = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]



# 6. Process per chromosome — now all filtering is string-based
chroms = [str(i) for i in range(1, 23)] + ["X"]


# 6. Process per chromosome
chroms = [str(i) for i in range(1, 23)] + ["X"]

for chrom in chroms:
    print(f"Processing Chromosome {chrom}...")
    try:
        chrom_df = df.filter(pl.col("chr") == chrom)
        if chrom_df.is_empty():
            print(f"  Chromosome {chrom}: no data found, skipping.")
            continue
        out = chrom_df.select([
            pl.col("chr").alias("CHROM"),
            pl.col("pos").alias("POS"),
            pl.col("ref").alias("REF"),
            pl.col("alt").alias("ALT"),
        ])
        # Write Mean
        out.with_columns(
            chrom_df["INFO_MEAN"].alias("INFO")
        ).write_csv(
            f"GRCh38_fingenR12_Mean_infoscore_chr{chrom}.tsv",
            separator="\t"
        )
        # Write Median
        out.with_columns(
            chrom_df["INFO_MEDIAN"].alias("INFO")
        ).write_csv(
            f"GRCh38_fingenR12_Median_infoscore_chr{chrom}.tsv",
            separator="\t"
        )
        os.system(f"bgzip -f GRCh38_fingenR12_Median_infoscore_chr{chrom}.tsv")
        os.system(f"tabix -s1 -b2 -e2 -S1 -f GRCh38_fingenR12_Median_infoscore_chr{chrom}.tsv.gz")
        os.system(f"bgzip -f GRCh38_fingenR12_Mean_infoscore_chr{chrom}.tsv")
        os.system(f"tabix -s1 -b2 -e2 -S1 -f GRCh38_fingenR12_Mean_infoscore_chr{chrom}.tsv.gz")
        print(f"  Chromosome {chrom}: {chrom_df.height} variants written.")
    except Exception as e:
        print(f"  Chromosome {chrom} error: {e}")

print("Successfully processed all specified chromosomes.")