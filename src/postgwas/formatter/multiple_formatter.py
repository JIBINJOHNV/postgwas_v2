import polars as pl
import os
import gc

def run_export_pipeline(input_df, configurations):
    """
    input_df: Your 20 million row Polars DataFrame
    configurations: A dictionary mapping filenames to column rename maps
    """
    # 1. Convert to LazyFrame
    # This ensures Polars only processes the columns needed for each specific file
    lazy_df = input_df.lazy()
    
    print(f"Total rows to process: {input_df.height:,}")
    print(f"Device detected: {'Apple Silicon' if 'arm' in os.uname().version.lower() else 'Intel Mac'}")
    print("-" * 30)

    for filename, rename_map in configurations.items():
        try:
            print(f"🚀 Starting export: {filename}...")
            
            # 2. Select and Rename
            # We use .select() with aliases for the fastest transformation
            output_lf = lazy_df.select([
                pl.col(old_name).alias(new_name) 
                for old_name, new_name in rename_map.items()
            ])
            
            # 3. Execute and Write
            # .collect() runs the actual computation
            # .write_csv() sends it to disk
            if filename.endswith(".csv"):
                output_lf.collect().write_csv(filename)
            elif filename.endswith((".txt", ".tsv")):
                output_lf.collect().write_csv(filename, separator="\t")
            elif filename.endswith(".parquet"):
                output_lf.collect().write_parquet(filename)
                
            print(f"✅ Successfully saved {filename}")
            
            # 4. Memory Management (Crucial for Mac)
            # Explicitly clear intermediate memory to keep 'Memory Pressure' low
            gc.collect()
            
        except Exception as e:
            print(f"❌ Error exporting {filename}: {e}")

# --- CONFIGURATION AREA ---

# Define your files here
# Format: "FileName": {"Original_Column": "New_Column"}
export_configs = {
    "GWAS_Summary_Portal.csv": {
        "chrcol": "CHR",
        "poscol": "BP",
        "pcol": "P-VALUE",
        "rsIDcol": "SNP"
    },
    "Locus_Analysis_Full.txt": {
        "uniq_id": "ID",
        "becol": "BETA",
        "secol": "STANDARD_ERROR",
        "eacol": "ALLELE_1",
        "neacol": "ALLELE_2"
    },
    "Internal_Archive.parquet": {
        "chrcol": "CHR",
        "poscol": "POS",
        "pcol": "P",
        "eafcol": "FREQ"
    }
}

# --- EXECUTION ---

if __name__ == "__main__":
    # Example: Loading your data (replace with your actual loading step)
    # df = pl.read_csv("your_massive_file.csv") 
    
    # Run the pipeline
    # run_export_pipeline(df, export_configs)
    
    print("\nAll tasks completed.")