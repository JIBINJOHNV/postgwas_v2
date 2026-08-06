import time
import os
import argparse
import pandas as pd
import pred_ld_functions
import builtins

# Default global ref directory (used by pred_ld_functions)
if not hasattr(builtins, 'REF_DIR'):
    builtins.REF_DIR = 'ref'

start = time.time()

def main():
    global imp_snp_list
    imp_snp_list = []
    version = '1.0.0'
    print("---------------------------------------------------------------------------------")
    print("PRED-LD : GWAS Summary Statistics Imputation")
    print("Version " + version + "; January 2025")
    print("Copyright (C) 2025 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("---------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(description="Process data in chunks.")
    parser.add_argument('--file-path', type=str, required=True, help='Input file path')
    parser.add_argument('--r2threshold', type=float, required=True, help='R2 threshold')
    parser.add_argument('--pop', type=str, required=True, help='Population')
    parser.add_argument('--maf', type=float, required=True, help='MAF input value')
    parser.add_argument('--ref', type=str, required=False,
                        help='LD Reference files (Pheno_Scanner, TOP_LD, Hap_Map or all_panels)',
                        default='all_panels')
    parser.add_argument('--imp_list', type=str, required=False,
                        help='A filename to define SNPs to impute (each SNP has a new line, no header)')
    parser.add_argument('--ref_dir', type=str, required=False, default='ref',
                        help='Base directory for reference files (default: ref)')

    parser.add_argument('--out_dir', type=str, required=True,
                        help='Base directory for reference files (default: ref)')
    
    args = parser.parse_args()

    # Assign arguments
    file_path = args.file_path
    r2threshold = args.r2threshold
    population = args.pop
    maf_input = args.maf
    ref_file = args.ref
    imp_snp_list_path = args.imp_list
    ref_dir = args.ref_dir
    out_dir = args.out_dir

    # Update the global REF_DIR used inside pred_ld_functions
    builtins.REF_DIR = ref_dir
    builtins.OUT_DIR = out_dir

    if imp_snp_list_path is not None:
        imp_snp_list = list(pd.read_csv(imp_snp_list_path, header=None)[0])

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return

    pred_ld_functions.process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list)

if __name__ == "__main__":
    main()

end = time.time()
print(f"Total Time: {end - start} seconds")
