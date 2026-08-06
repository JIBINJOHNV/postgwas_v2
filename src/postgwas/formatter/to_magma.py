import subprocess
import polars as pl
from pathlib import Path

from pathlib import Path
import subprocess
import textwrap
import polars as pl


def vcf_to_magma(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert a summary-statistics VCF file into MAGMA-compatible inputs.

    Outputs
    -------
    <output_folder>/<sample_name>_magma_snp_loc.tsv
    <output_folder>/<sample_name>_magma_P_val.tsv
    <output_folder>/<sample_name>_vcf_to_magma.log
    """

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # Files
    # ---------------------------------------------------------
    raw_output_file = output_dir / f"{sample_name}_magma_input.tsv"
    snp_loc_file = output_dir / f"{sample_name}_magma_snp_loc.tsv"
    pval_file = output_dir / f"{sample_name}_magma_P_val.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_magma.log"

    # ---------------------------------------------------------
    # Build bcftools pipeline
    # ---------------------------------------------------------
    cmd = textwrap.dedent(f"""
        set -o pipefail
        (
            printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tBETA\\tSE\\tN_COL\\tAF\\tLP\\n"
            bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\t[%LP]\\n' | \
            sed 's|:|_|g'
        ) > "{raw_output_file}"
    """)

    # ---------------------------------------------------------
    # Run command + logging
    # ---------------------------------------------------------
    try:
        with open(log_file, "w") as lf:
            lf.write("### vcf_to_magma command\n\n")
            lf.write(cmd + "\n")
            lf.write("\n### Command output (stdout / stderr)\n\n")

            subprocess.run(
                cmd,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=lf,
                stderr=lf,
            )

    except subprocess.CalledProcessError as e:
        with open(log_file, "a") as lf:
            lf.write("\n❌ bcftools pipeline failed\n")
            lf.write(str(e) + "\n")

        raise RuntimeError(
            f"vcf_to_magma failed for {vcf_path}. "
            f"See log file: {log_file}"
        )

    # ---------------------------------------------------------
    # Read intermediate file (Polars)
    # ---------------------------------------------------------
    try:
        df = pl.read_csv(
            raw_output_file,
            separator="\t",
            null_values=["."],
            schema_overrides={
                "CHR": pl.String,
                "BP": pl.Int64,
                "LP": pl.Float64,
                "N_COL": pl.Float64,
            },
        )
    except Exception as e:
        raise ValueError(
            f"Failed to read intermediate MAGMA file: {raw_output_file}\n{e}"
        )

    # ---------------------------------------------------------
    # Validate required columns
    # ---------------------------------------------------------
    required_cols = [
        "SNP", "CHR", "BP", "REF", "ALT",
        "LP", "BETA", "SE", "AF", "N_COL",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {raw_output_file}: {', '.join(missing)}"
        )

    # ---------------------------------------------------------
    # Convert -log10(P) → P (numerically safe)
    # ---------------------------------------------------------
    df = df.with_columns(
        pl.when(pl.col("LP") > 320)
        .then(1e-320)
        .otherwise(10 ** (-pl.col("LP")))
        .alias("P")
    )

    # ---------------------------------------------------------
    # Select and reorder columns
    # ---------------------------------------------------------
    df_clean = df.select([
        "SNP",
        "CHR",
        "BP",
        "REF",
        "ALT",
        "P",
        "BETA",
        "SE",
        "AF",
        "N_COL",
    ])

    # ---------------------------------------------------------
    # Write MAGMA inputs
    # ---------------------------------------------------------
    df_clean.select(["SNP", "CHR", "BP"]).write_csv(
        snp_loc_file, separator="\t"
    )

    df_clean.select(["SNP", "P", "N_COL"]).write_csv(
        pval_file, separator="\t"
    )

    # ---------------------------------------------------------
    # Final log summary
    # ---------------------------------------------------------
    with open(log_file, "a") as lf:
        lf.write("\n### MAGMA outputs generated successfully\n")
        lf.write(f"SNP location file: {snp_loc_file.name}\n")
        lf.write(f"P-value file: {pval_file.name}\n")

    return {
        "snp_loc_file": str(snp_loc_file),
        "pval_file": str(pval_file),
        "log_file": str(log_file),
    }
