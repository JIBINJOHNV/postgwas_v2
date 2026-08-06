import polars as pl
import os
import sys
from typing import Tuple, List

def validate_gwas_config(config: dict, df: pl.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validates GWAS configuration against DataFrame columns.
    
    Returns:
        Tuple[bool, List[str]]: 
          - bool: True if valid, False if failed.
          - List[str]: A list of specific error messages describing what is missing.
    """
    errors = []
    warnings = []
    
    def is_valid(value):
        return value not in [None, "", "NA", "na", "NaN"]
    
    # 1. Setup Rules
    or_pairs = [
        ("chr_col", "chr_pos_col"),
        ("pos_col", "chr_pos_col"),
        ("beta_or_col", "imp_z_col"),
        ("eaf_col", "eaffile"),
        ("imp_info_col", "infofile"),
        ("ncontrol_col", "ncontrol"),
    ]
    and_pairs = [("ea_col", "oa_col")]
    required = ["pval_col", "gwas_outputname", "sumstat_file"]
    
    df_columns = set(df.columns)
    
    # 2. Check OR pairs
    for a, b in or_pairs:
        a_val, b_val = config.get(a), config.get(b)
        
        # Rule 1: Config must have at least one defined
        if not (is_valid(a_val) or is_valid(b_val)):
            errors.append(f"Missing pair: You must provide either '{a}' or '{b}'. Both are NA.")
            continue

        # Rule 2: The defined one must exist in DataFrame (or be a valid file/scalar)
        valid_in_df = False
        
        # Check A
        if is_valid(a_val) and a_val in df_columns: valid_in_df = True
        
        # Check B
        if is_valid(b_val):
            # If B is a column, check existence
            if b_val in df_columns: 
                valid_in_df = True
            # If B is a file/number override (not a column), it's valid if simply defined
            elif "file" in b or b == "ncontrol":
                valid_in_df = True
        
        if not valid_in_df:
            errors.append(f"Missing column: Config defined '{a}'='{a_val}' or '{b}'='{b_val}', but neither were found in DataFrame.")

    # 3. Check AND pairs
    for a, b in and_pairs:
        a_val, b_val = config.get(a), config.get(b)
        
        if not (is_valid(a_val) and is_valid(b_val)):
            errors.append(f"Missing pair: Both '{a}' and '{b}' are required in config.")
        else:
            if a_val not in df_columns:
                errors.append(f"Missing column: '{a_val}' (for {a}) not found in DataFrame.")
            if b_val not in df_columns:
                errors.append(f"Missing column: '{b_val}' (for {b}) not found in DataFrame.")

    # 4. Check Required Singles
    for key in required:
        val = config.get(key)
        if not is_valid(val):
            errors.append(f"Missing key: Mandatory config '{key}' is empty/NA.")
        elif key == "pval_col" and val not in df_columns:
            errors.append(f"Missing column: P-value column '{val}' not found in DataFrame.")

    # 5. Check File Existence
    if is_valid(config.get("sumstat_file")) and not os.path.exists(config["sumstat_file"]):
        errors.append(f"File not found: {config['sumstat_file']}")

    # 6. Optional N_Case Check
    ncase_col = config.get("ncase_col")
    ncase_val = config.get("ncase")
    
    has_ncase = False
    if is_valid(ncase_col) and ncase_col in df_columns: has_ncase = True
    elif is_valid(ncase_val): has_ncase = True
    
    if not has_ncase:
        warnings.append("\t\t\t ℹ️  Note: Only Total N was provided. Assuming Quantitative Phenotype.")
        warnings.append("\t\t\t     ℹ️  Note: The VCF 'Total Sample Size' (SS) will be set to N.")
        warnings.append("\t\t\t     ℹ️  Note: The VCF 'Effective Sample Size' (NEF) will also be set to N.")
    # --- FINAL REPORTING ---
    # Print Warnings
    for w in warnings:
        print(w)

    # Print Errors
    if errors:
        print("\n" + "="*60)
        print(f"🛑 VALIDATION FAILED: Found {len(errors)} error(s)")
        print("="*60)
        for e in errors:
            print(f" • {e}")
        print("="*60 + "\n")
        sys.stdout.flush() 
        return False, errors

    print("\t\t\t\t ✅Configuration and Data validated successfully.")
    return True, []