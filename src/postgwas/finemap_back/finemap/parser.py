import polars as pl
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def parse_finemap_cred_file(cred_file: Path, locus_id: str) -> pl.DataFrame:
    """
    Parses the .cred file (SNP-level statistics).
    Returns a dataframe with columns: locus_id, rsid, prob, log10bf, etc.
    """
    if not cred_file.exists():
        logger.warning(f"{locus_id}: .cred file not found.")
        return pl.DataFrame()

    try:
        # FINEMAP .cred is space/tab separated with a header
        # Standard columns: index, rsid, chromosome, position, allele1, allele2, maf, beta, se, prob, log10bf, group, ...
        df = pl.read_csv(cred_file, separator=" ", null_values=["NA"])
        
        # Add locus identifier
        df = df.with_columns(pl.lit(locus_id).alias("locus_id"))
        
        # Rename 'prob' to 'PIP' (Posterior Inclusion Probability) for standard nomenclature
        if "prob" in df.columns:
            df = df.rename({"prob": "PIP"})
            
        return df
    except Exception as e:
        logger.error(f"Error parsing .cred file for {locus_id}: {e}")
        return pl.DataFrame()


def parse_finemap_config_file(config_file: Path, locus_id: str) -> pl.DataFrame:
    """
    Parses the .config file (Causal Set configurations).
    Returns a dataframe summarising the credible sets.
    """
    if not config_file.exists():
        return pl.DataFrame()
        
    try:
        df = pl.read_csv(config_file, separator=" ", null_values=["NA"])
        df = df.with_columns(pl.lit(locus_id).alias("locus_id"))
        return df
    except Exception as e:
        logger.error(f"Error parsing .config file for {locus_id}: {e}")
        return pl.DataFrame()


def parse_finemap_log_for_logBF(log_file: Path) -> float:
    """
    Optional: Scrapes the .log file to find the final LogBF of the model.
    """
    # Implementation depends on exact log format need
    return 0.0