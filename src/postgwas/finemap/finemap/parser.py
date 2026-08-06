import polars as pl
from pathlib import Path
import logging
import math
import re

import pandas as pd

logger = logging.getLogger(__name__)

def parse_finemap_cred_file(cred_file: Path, locus_id: str) -> pl.DataFrame:
    """
    Parse a FINEMAP .credK table while skipping its model-level comment lines.
    Multi-causal files retain their cred1/prob1, cred2/prob2, ... columns.
    """
    if not cred_file.exists():
        logger.warning(f"{locus_id}: .cred file not found.")
        return pl.DataFrame()

    try:
        pandas_df = pd.read_csv(
            cred_file, sep=r"\s+", comment="#", na_values=["NA", "."]
        )
        if pandas_df.empty:
            return pl.DataFrame()
        df = pl.DataFrame(pandas_df.to_dict(orient="list")).with_columns(
            pl.lit(locus_id).alias("locus_id")
        )
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
        pandas_df = pd.read_csv(
            config_file, sep=r"\s+", comment="#", na_values=["NA", "."]
        )
        if pandas_df.empty:
            return pl.DataFrame()
        df = pl.DataFrame(pandas_df.to_dict(orient="list")).with_columns(
            pl.lit(locus_id).alias("locus_id")
        )
        return df
    except Exception as e:
        logger.error(f"Error parsing .config file for {locus_id}: {e}")
        return pl.DataFrame()


def parse_finemap_log_for_logBF(log_file: Path) -> float:
    """Extract a reported log10 Bayes factor, returning NaN when unavailable."""
    if not Path(log_file).is_file():
        return math.nan
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    patterns = [
        re.compile(rf"log10\s*BF\s*[:=]\s*({number})", re.IGNORECASE),
        re.compile(rf"log10bf\s*[:=]?\s*({number})", re.IGNORECASE),
    ]
    try:
        value = math.nan
        for line in Path(log_file).read_text(errors="replace").splitlines():
            for pattern in patterns:
                match = pattern.search(line)
                if match:
                    value = float(match.group(1))
        return value
    except (OSError, ValueError):
        return math.nan
