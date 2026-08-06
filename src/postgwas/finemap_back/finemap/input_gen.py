import polars as pl
from pathlib import Path

def write_finemap_z_file(locus_ss: pl.DataFrame, outfile: Path):
    """
    Writes the .z file required by FINEMAP.
    Format: Space-delimited, columns: rsid chromosome position allele1 allele2 maf beta se
    """
    cols = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    
    # Validation
    missing = [c for c in cols if c not in locus_ss.columns]
    if missing:
        raise ValueError(f"Input sumstats missing columns for FINEMAP .z generation: {missing}")

    locus_ss.select(cols).write_csv(outfile, separator=" ")


def write_snp_file(locus_ss: pl.DataFrame, outfile: Path):
    """Writes the list of variants (RSIDs) for PLINK extraction."""
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