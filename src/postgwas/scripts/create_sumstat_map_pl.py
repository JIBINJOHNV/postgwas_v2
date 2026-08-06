import os
import gzip
import bz2
import argparse
import polars as pl
from pathlib import Path
from typing import Dict, List, Tuple

# ---------------------------------------------------------
# 1. DEFINE THE COLUMN MAPPING CONSTANTS
# ---------------------------------------------------------
COLUMN_MAP_SOURCE = """
#CHR	chr_col
#CHROM	chr_col
CHR	chr_col
CHROM	chr_col
chromosome	chr_col
CHROMOSOME	chr_col
CHROMOSOMES	chr_col
CHROMS	chr_col
HG19CHRC	chr_col
SEQ	chr_col
SEQNAME	chr_col
SEQNAMES	chr_col
SEQS	chr_col
BASE	pos_col
BASE GRCH37	pos_col
BASE GRCH38	pos_col
BASE PAIR POSITION	pos_col
BASE_GRCH37	pos_col
BASE_GRCH38	pos_col
BASE_PAIR_LOCATION	pos_col
base_pair_location	pos_col
BASE-GRCH37	pos_col
BASE-GRCH38	pos_col
BASE.GRCH37	pos_col
BASE.GRCH38	pos_col
BASE.PAIR.LOCATION	pos_col
BASE.PAIR.POSITION	pos_col
BASEGRCH37	pos_col
BASEGRCH38	pos_col
BP	pos_col
BP_HG19	pos_col
BP_HG38	pos_col
BP-HG19	pos_col
BP-HG38	pos_col
BP.HG19	pos_col
BP.HG38	pos_col
BPHG19	pos_col
BPHG38	pos_col
BPOS	pos_col
GEN_POS	pos_col
GEN-POS	pos_col
GEN.POS	pos_col
GENOMEPOS	pos_col
GENPOS	pos_col
POS	pos_col
POS GRCH37	pos_col
POS GRCH38	pos_col
POS_B37	pos_col
POS_B38	pos_col
POS_GRCH37	pos_col
POS_GRCH38	pos_col
POS-GRCH37	pos_col
POS-GRCH38	pos_col
POS.GRCH37	pos_col
POS.GRCH38	pos_col
POS37	pos_col
POS38	pos_col
POSGEN	pos_col
POSGRCH37	pos_col
POSGRCH38	pos_col
POSITION	pos_col
POSITION_BP	pos_col
POSITION37	pos_col
POSITION38	pos_col
START	pos_col
#ID	snp_id_col
ID	snp_id_col
ID_DBSNP49	snp_id_col
MARKER	snp_id_col
MARKER NAME	snp_id_col
MARKER_NAME	snp_id_col
MARKER.NAME	snp_id_col
MARKERNAME	snp_id_col
PREDICTOR	snp_id_col
RS	snp_id_col
RS ID	snp_id_col
RS_ID	snp_id_col
RS_NUMBER	snp_id_col
RS_NUMBERS	snp_id_col
RS.ID	snp_id_col
RSID	snp_id_col
RSIDS	snp_id_col
rsid	snp_id_col
RSNUMBERS	snp_id_col
SNP	snp_id_col
SNP ID	snp_id_col
SNP_ID	snp_id_col
SNP.ID	snp_id_col
SNPID	snp_id_col
SNPNAME	snp_id_col
VARIANT	snp_id_col
VARIANT_ID	snp_id_col
A_1	ea_col
A-1	ea_col
A.1	ea_col
A1	ea_col
AL1	ea_col
ALLELE 1	ea_col
ALLELE_1	ea_col
ALLELE-1	ea_col
ALLELE.1	ea_col
ALLELE1	ea_col
ALT	ea_col
ALT.ALLELE	ea_col
ALTALLELE	ea_col
ALTERNATIVE	ea_col
ALTERNATIVE ALLELE	ea_col
ALTERNATIVE_ALLELE	ea_col
ALTERNATIVE-ALLELE	ea_col
ALTERNATIVE.ALLELE	ea_col
ALTERNATIVEALLELE	ea_col
EA	ea_col
EFF ALLELE	ea_col
EFF_ALLELE	ea_col
EFF-ALLELE	ea_col
EFF.ALLELE	ea_col
EFFALLELE	ea_col
EFFECT ALLELE	ea_col
EFFECT_ALLELE	ea_col
effect_allele	ea_col
EFFECT-ALLELE	ea_col
EFFECT.ALLELE	ea_col
EFFECTALLELE	ea_col
TESTED ALLELE	ea_col
TESTED_ALLELE	ea_col
TESTED-ALLELE	ea_col
TESTED.ALLELE	ea_col
TESTEDALLELE	ea_col
A_2	oa_col
A-2	oa_col
A.2	oa_col
A2	oa_col
AL2	oa_col
OA	oa_col
ALLELE0	oa_col
ALLELE 2	oa_col
ALLELE_2	oa_col
ALLELE-2	oa_col
ALLELE.2	oa_col
ALLELE2	oa_col
NEA	oa_col
NON EFFECT ALLELE	oa_col
NON_EFFECT_ALLELE	oa_col
NON-EFFECT-ALLELE	oa_col
NON.EFFECT.ALLELE	oa_col
NONEA	oa_col
NONEFFECT_ALLELE	oa_col
neffect_allele	oa_col
NONEFFECTALLELE	oa_col
OTHER ALLELE	oa_col
OTHER_ALLELE	oa_col
OTHER-ALLELE	oa_col
OTHER.ALLELE	oa_col
OTHERALLELE	oa_col
REF	oa_col
REFALLELE	oa_col
REFERENCE	oa_col
REFERENCE ALLELE	oa_col
REFERENCE_ALLELE	oa_col
REFERENCE-ALLELE	oa_col
REFERENCE.ALLELE    oa_col
REFERENCEALLELE oa_col
AlternateAllele oa_col
AF	eaf_col
A1FREQ	eaf_col
AF_ALT	eaf_col
AF_EFF	eaf_col
AF-ALT	eaf_col
AF.ALT	eaf_col
AF.EFF	eaf_col
ALL_AF	eaf_col
ALL_META_AF	eaf_col
ALLELE_FREQ	eaf_col
ALLELE_FREQUENCY	eaf_col
ALLELE_FRQ	eaf_col
ALLELEFREQ	eaf_col
ALT_AF	eaf_col
ALT_FREQ	eaf_col
ALT-AF	eaf_col
ALT.AF	eaf_col
ALTERNATIVE_AF	eaf_col
EAF	eaf_col
EAF_HRC	eaf_col
EAF_MAX	eaf_col
EFF_AF	eaf_col
EFFECT_AF	eaf_col
EFFECT_ALLELE_FREQ	eaf_col
EFFECT_ALLELE_FREQUENCY	eaf_col
EFFECT_ALLELE_FRQ	eaf_col
EFFECTALLELEFREQ	eaf_col
EFFECTALLELEMAXFREQ	eaf_col
EFFECTALLELEMINFREQ	eaf_col
effect_allele_frequency	eaf_col
EST.FRQ	eaf_col
F_U	eaf_col
FCON	eaf_col
FREQ	eaf_col
FREQ_EFFECT	eaf_col
FREQ_EFFECT_ALLELE	eaf_col
FREQ_EUROPEAN_1000GENOMES	eaf_col
FREQ_HAPMAP	eaf_col
FREQ_TESTED_ALLELE	eaf_col
FREQ_TESTED_ALLELE_IN_HRS	eaf_col
FREQUENCY	eaf_col
FRQ	eaf_col
FRQ_EFFECT_ALLELE	eaf_col
FRQ_TESTED_ALLELE	eaf_col
HRC_FRQ_A1	eaf_col
Freq1	eaf_col
FRQ_U*	eaf_col
FRQMAX	eaf_col
INC_AF	eaf_col
MAF	eaf_col
MINOR_AF	eaf_col
POOLED_ALT_AF	eaf_col
TESTED_AF	eaf_col
EffectAlleleFrequency	eaf_col
ALL_INV_VAR_META_BETA	beta_or_col
ALT_EFFSIZE	beta_or_col
B	beta_or_col
beta	beta_or_col
BETA	beta_or_col
BETA_SLEEPDURATION	beta_or_col
EFFECT	beta_or_col
EFFECT_BETA	beta_or_col
EFFECT_SIZE	beta_or_col
EFFECT_WEIGHT	beta_or_col
EFFECTS	beta_or_col
ES	beta_or_col
EST	beta_or_col
INV_VAR_META_BETA	beta_or_col
LOG ODDS	beta_or_col
LOG_ODDS	beta_or_col
LOG-ODDS	beta_or_col
LOG.ODDS	beta_or_col
LOGOR	beta_or_col
MTAG_BETA	beta_or_col
ODDS RATIO	beta_or_col
ODDS_RAT	beta_or_col
ODDS_RATIO	beta_or_col
ODDS-RAT	beta_or_col
ODDS-RATIO	beta_or_col
ODDS.RAT	beta_or_col
ODDS.RATIO	beta_or_col
ODDSRATIO	beta_or_col
OR	beta_or_col
STDBETA	beta_or_col
EffectSize.Beta	beta_or_col
EffectSize.OR	beta_or_col
ALL_INV_VAR_META_SEBETA	se_col
INV_VAR_META_SEBETA	se_col
MTAG_SE	se_col
SE	se_col
SE_DGC	se_col
SE_SLEEPDURATION	se_col
SEBETA	se_col
STANDARD ERROR	se_col
STANDARD_ERROR	se_col
STANDARD-ERROR	se_col
STANDARD.ERROR	se_col
STANDARDERROR	se_col
standard_error	se_col
STDERR	se_col
STDERRLOGOR	se_col
ALL_INV_VAR_META_P	pval_col
INV_VAR_META_P	pval_col
LOG10_P	pval_col
LOG10P	pval_col
MLOG10P	pval_col
MTAG_PVAL	pval_col
p_wald	pval_col
P	pval_col
P_BOLT_LMM	pval_col
P_DGC	pval_col
P_LINREG	pval_col
P_SLEEPDURATION	pval_col
P_VAL	pval_col
p_value	pval_col
P_VALUE	pval_col
P-VAL	pval_col
P-VALUE	pval_col
P.SE	pval_col
P.VAL	pval_col
P.VALUE	pval_col
PVAL	pval_col
PVALUE	pval_col
neg_log_10_p_value	pval_col
Pval_Estimate	pval_col
EZ	imp_z_col
z_score	imp_z_col
MTAG_Z	imp_z_col
Z	imp_z_col
Z-SCORE	imp_z_col
ZSCORE	imp_z_col
Z_Estimate	imp_z_col
CONTROLS_N	ncontrol_col
n	ncontrol_col
N	ncontrol_col
N_CON	ncontrol_col
N_CONTROL	ncontrol_col
N_CONTROLS	ncontrol_col
N_CTRL	ncontrol_col
NCO	ncontrol_col
NCON	ncontrol_col
NCONTROL	ncontrol_col
Neff	ncontrol_col
N_eff	ncontrol_col
Weight	ncontrol_col
CASES_N	ncase_col
N_CAS	ncase_col
N_CASE	ncase_col
N_CASES	ncase_col
N_EVENT	ncase_col
NC	ncase_col
NCA	ncase_col
NCAS	ncase_col
NCASE	ncase_col
NCASES	ncase_col
INFO	imp_info_col
IMPUTATIONACCURACY	imp_info_col
MEDIAN_INFO	imp_info_col
IMPINFO	imp_info_col
IMPUTATION	imp_info_col
R2HAT	imp_info_col
RSQ	imp_info_col
minINFO	imp_info_col
chr:pos	chr_pos_col
"""

# Parse column mapping
COLUMN_MAPPING = {}
for line in COLUMN_MAP_SOURCE.strip().split('\n'): 
    if line.strip(): 
        parts = line.split('\t')
        if len(parts) == 2:
            key, val = parts
            COLUMN_MAPPING[key.upper()] = val

COLUMN_MAPPING_ORDERED = []
for line in COLUMN_MAP_SOURCE.strip().split('\n'): 
    if line.strip(): 
        parts = line.split('\t')
        if len(parts) == 2:
            key, val = parts
            # Store as (Header_to_look_for, Target_Category)
            COLUMN_MAPPING_ORDERED.append((key.upper(), val))

# --- CORRECTED OUTPUT_COLUMNS (Fixed Spelling resource_folder) ---
OUTPUT_COLUMNS = [
    "sumstat_file", "gwas_outputname", "chr_col", "pos_col", "snp_id_col",
    "ea_col", "oa_col", "eaf_col", "beta_or_col", "se_col", "imp_z_col",
    "pval_col", "ncontrol_col", "ncase_col", "ncontrol", "ncase",
    "imp_info_col", "infofile", "infocolumn", "eaffile",
    "eafcolumn", "liftover", "chr_pos_col", "resourse_folder", "output_folder"
]

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

def get_file_handle(file_path: Path):
    if file_path.suffix == '.gz' or file_path.suffix == '.bgz':
        return gzip.open(file_path, 'rt', encoding='utf-8', errors='replace')
    elif file_path.suffix == '.bz2':
        return bz2.open(file_path, 'rt', encoding='utf-8', errors='replace')
    else:
        return open(file_path, 'r', encoding='utf-8', errors='replace')

def detect_delimiter(header_line: str) -> str:
    if '\t' in header_line:
        return '\t'
    elif ',' in header_line:
        return ','
    else:
        return r'\s+'

def get_header_and_delimiter(file_path: Path) -> Tuple[List[str], str]:
    try:
        with get_file_handle(file_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("##"):
                    continue
                delimiter = detect_delimiter(line)
                if delimiter == r'\s+':
                    headers = line.split()
                else:
                    headers = line.split(delimiter)
                headers = [h.strip('"').strip("'") for h in headers]
                return headers, delimiter
        return [], ''
    except Exception as e:
        print(f"Error reading header for {file_path}: {e}")
        return [], ''

# def map_header_to_standard(headers: List[str]) -> Dict[str, str]:
#     mapped_data = {target: "NA" for target in set(COLUMN_MAPPING.values())}
#     for h in headers:
#         h_upper = h.upper()
#         if h_upper in COLUMN_MAPPING:
#             target_col = COLUMN_MAPPING[h_upper]
#             if mapped_data[target_col] == "NA":
#                 mapped_data[target_col] = h
#         elif h_upper.startswith("FRQ_U"):
#             if mapped_data["eaf_col"] == "NA":
#                 mapped_data["eaf_col"] = h
#         elif h_upper.startswith("FREQ_U"):
#             if mapped_data["eaf_col"] == "NA":
#                 mapped_data["eaf_col"] = h
#     return mapped_data

def map_header_to_standard(headers: List[str]) -> Dict[str, str]:
    # Normalize file headers for comparison
    file_headers_upper = {h.upper(): h for h in headers}
    # Initialize result dictionary
    mapped_data = {target: "NA" for target in set(COLUMN_MAPPING.values())}
    # Track which categories are already 'locked' (first map found)
    assigned_categories = set()
    # Iterate through the PRIORITY list (COLUMN_MAP_SOURCE order)
    for source_variant, target_col in COLUMN_MAPPING_ORDERED:
        # If this category hasn't been filled yet...
        if target_col not in assigned_categories:
            # ...and if the current priority variant exists in the file
            if source_variant in file_headers_upper:
                mapped_data[target_col] = file_headers_upper[source_variant]
                # STOP further mapping for this category
                assigned_categories.add(target_col)
    # Handle special prefix-based logic for EAF if still NA
    if mapped_data["eaf_col"] == "NA":
        for h in headers:
            h_up = h.upper()
            if h_up.startswith("FRQ_U") or h_up.startswith("FREQ_U"):
                mapped_data["eaf_col"] = h
                break
    return mapped_data

def generate_sumstat_mapping_table(
    input_dir: str, 
    output_full_path: str,
    resource_folder: str,
    harmonisation_output_path: str
) -> None:
    input_path = Path(input_dir)
    harmonisation_base = Path(harmonisation_output_path).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    results = []
    # Extension whitelist: only process files that likely contain GWAS data
    valid_extensions = {'.gz', '.txt', '.tsv', '.bgz', '.sumstats', '.meta','.bz2','.csv'}
    
    for file_path in input_path.iterdir():
        # 1. Skip directories, hidden files, and non-GWAS extensions
        if (file_path.is_file() and 
            not file_path.name.startswith('.') and 
            file_path.suffix.lower() in valid_extensions):
            
            file_name = file_path.name
            gwas_name = file_name.split('.')[0] 
            
            headers, delimiter = get_header_and_delimiter(file_path)
            if not headers:
                continue
            
            mapped_cols = map_header_to_standard(headers)
            
            # NOTE: We do NOT create the directory here anymore!
            row = {
                "sumstat_file": str(file_path.absolute()),
                "gwas_outputname": gwas_name,                
                "chr_col": mapped_cols.get("chr_col", "NA"),
                "pos_col": mapped_cols.get("pos_col", "NA"),
                "snp_id_col": mapped_cols.get("snp_id_col", "NA"),
                "ea_col": mapped_cols.get("ea_col", "NA"),
                "oa_col": mapped_cols.get("oa_col", "NA"),
                "eaf_col": mapped_cols.get("eaf_col", "NA"),
                "beta_or_col": mapped_cols.get("beta_or_col", "NA"),
                "se_col": mapped_cols.get("se_col", "NA"),
                "imp_z_col": mapped_cols.get("imp_z_col", "NA"),
                "pval_col": mapped_cols.get("pval_col", "NA"),
                "ncontrol_col": mapped_cols.get("ncontrol_col", "NA"),
                "ncase_col": mapped_cols.get("ncase_col", "NA"),
                "imp_info_col": mapped_cols.get("imp_info_col", "NA"),
                "ncontrol": "NA",
                "ncase": "NA",
                "infofile": "NA",
                "infocolumn": "NA",
                "eaffile": "NA",
                "eafcolumn": "NA",
                "liftover": "Yes",
                "chr_pos_col": mapped_cols.get("chr_pos_col", "NA"),
                "resourse_folder": resource_folder,
            }
            results.append(row)

    if not results:
        print("⚠️ No valid summary statistics files found.")
        return

    df = pl.DataFrame(results)
    
    # Ensure all columns exist
    for col in OUTPUT_COLUMNS:
        if col not in df.columns:
            df = df.with_columns(pl.lit("NA").alias(col))

    # Filter out rows where the majority of content columns are "NA"
    content_cols = ["chr_col", "pos_col", "ea_col", "oa_col", "beta_or_col", "pval_col"]
    df = df.with_columns(
        na_count = pl.sum_horizontal([pl.col(c) == "NA" for c in content_cols])
    )

    original_count = df.height
    df = df.filter(pl.col("na_count") < (len(content_cols) // 2 + 1))
    
    removed_count = original_count - df.height
    if removed_count > 0:
        print(f"🧹 Removed {removed_count} non-GWAS or poorly mapped files.")

    if df.height == 0:
        print("❌ All files were filtered out.")
        return

    # --- FINAL STEP: Create folders and set 'output_folder' ONLY for valid files ---
    def finalize_row(row_dict):
        name = row_dict["gwas_outputname"]
        out_path = harmonisation_base / name
        out_path.mkdir(parents=True, exist_ok=True)
        row_dict["output_folder"] = str(out_path) + "/"
        return row_dict

    final_data = [finalize_row(r) for r in df.to_dicts()]
    df = pl.DataFrame(final_data)

    # Final selection and write
    df = df.select(OUTPUT_COLUMNS)
    df.write_csv(output_full_path)
    
    print(f"Successfully processed {df.height} files.")
    print(f"Configuration saved to: {output_full_path}")
    
# ---------------------------------------------------------
# 3. CLI SETUP
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Auto-detect sumstat columns.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("-r", "--resource-folder", 
                        default="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/")
    parser.add_argument("-ho", "--harmonisation-output-path", required=True)

    args = parser.parse_args()
    generate_sumstat_mapping_table(
        input_dir=args.input,
        output_full_path=args.output_path,
        resource_folder=args.resource_folder,
        harmonisation_output_path=args.harmonisation_output_path
    )

if __name__ == "__main__":
    main()