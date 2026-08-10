from pathlib import Path
import pandas as pd


def load_gene_list(gene_file):
    """
    Load gene symbols from a file.
    Uses FIRST column only.
    """
    gene_file = Path(gene_file)

    if not gene_file.exists():
        raise FileNotFoundError(f"Gene file not found: {gene_file}")

    df = pd.read_csv(
        gene_file,
        sep=None,            # auto-detect delimiter
        engine="python",
        header=0
    )

    genes = (
        df.iloc[:, 0]
        .astype(str)
        .str.strip()
        .tolist()
    )

    genes = sorted({g for g in genes if g and g.lower() != "nan"})

    if not genes:
        raise ValueError("No valid gene symbols found in first column")

    return genes
