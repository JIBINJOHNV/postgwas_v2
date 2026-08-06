import numpy as np

def check_eaf_or_maf(df):
    if "maf" in df.columns:
        raise ValueError("MAF column provided — please use EAF instead.")
    if "eaf" not in df.columns:
        raise ValueError("EAF column missing and no external file provided.")
    print("EAF check passed.")
    return df

def compute_effective_n(df):
    if "ncases" in df.columns and "ncontrol" in df.columns:
        df["neff"] = 4 / (1/df["ncases"] + 1/df["ncontrol"])
    else:
        df["neff"] = np.nan
    print("Calculated effective sample size (neff).")
    return df
