import polars as pl
from pathlib import Path

from postgwas.finemap.defaults import DEFAULT_MAX_MAF

def write_finemap_z_file(
    locus_ss: pl.DataFrame,
    outfile: Path,
    max_maf: float = DEFAULT_MAX_MAF,
):
    """
    Writes the .z file required by FINEMAP.
    Format: Space-delimited, columns: rsid chromosome position allele1 allele2 maf beta se
    """
    cols = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    
    # Validation
    missing = [c for c in cols if c not in locus_ss.columns]
    if missing:
        raise ValueError(f"Input sumstats missing columns for FINEMAP .z generation: {missing}")

    if locus_ss.is_empty():
        raise ValueError("Cannot write an empty FINEMAP .z file")
    if locus_ss.select(pl.any_horizontal(pl.col(cols).is_null()).any()).item():
        raise ValueError("FINEMAP .z input contains null values")
    if locus_ss.select(pl.col("rsid").is_duplicated().any()).item():
        raise ValueError("FINEMAP .z input contains duplicate rsIDs")
    invalid_numeric = locus_ss.filter(
        ~pl.col("maf").is_finite()
        | (pl.col("maf") <= 0)
        | (pl.col("maf") > float(max_maf))
        | ~pl.col("beta").is_finite()
        | ~pl.col("se").is_finite()
        | (pl.col("se") <= 0)
    )
    if invalid_numeric.height:
        raise ValueError("FINEMAP .z input has invalid MAF, beta, or standard error")

    locus_ss.select(cols).write_csv(outfile, separator=" ")


def write_snp_file(locus_ss: pl.DataFrame, outfile: Path):
    """Writes the list of variants (RSIDs) for PLINK extraction."""
    if "rsid" not in locus_ss.columns or locus_ss.is_empty():
        raise ValueError("Cannot write a SNP file without rsIDs")
    if locus_ss.select(pl.col("rsid").is_null().any()).item():
        raise ValueError("SNP list contains null rsIDs")
    locus_ss.select("rsid").write_csv(outfile, has_header=False)


def create_ldstore_master(
    master_path: Path, 
    z_file: Path, 
    bgen_file: Path, 
    bgi_file: Path, 
    bcor_file: Path, 
    ld_matrix: Path, 
    n_samples: int
):
    """Creates the master file required by LDstore."""
    with master_path.open("w") as f:
        f.write("z;bgen;bgi;bcor;ld;n_samples\n")
        f.write(
            f"{z_file};{bgen_file};{bgi_file};"
            f"{bcor_file};{ld_matrix};{n_samples}\n"
        )


def create_finemap_master(
    master_path: Path,
    z_file: Path,
    ld_file: Path,
    snp_file: Path,
    config_file: Path,
    cred_file: Path,
    log_file: Path,
    n_samples: int
):
    """Creates the master file required by FINEMAP."""
    with master_path.open("w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(
            f"{z_file};{ld_file};{snp_file};"
            f"{config_file};{cred_file};{log_file};{n_samples}\n"
        )
