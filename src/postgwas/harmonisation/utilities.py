import os
import glob
from pathlib import Path

def combine_logs_per_chromosome(outdir):
    log_dir = Path(outdir)
    # ORDER YOU PROVIDED
    log_order = [
        'add_or_calculate_info.log', 'beta_se_from_z.log',
        'calculate_se_from_beta_and_pvalue.log', 'calculate_zscore_from_beta_and_se.log',
        'detect_pval.log', 'eaf.log','fix_chr_pos_column.log',
        'gwastovcf.log', 'is_beta_or_or.log', 'snp.log',
        'ADHD2022_iPSYCH_deCODE_PGC_chr11.log', 'bcftools_annot.log', 'export.log'
    ]
    # CHROMOSOMES
    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    for chrom in chromosomes:
        combined_log = log_dir / f"chr{chrom}_combined.log"
        with combined_log.open("w") as outfile:
            outfile.write(f"=== Combined Log for Chromosome {chrom} ===\n\n")
            for step in log_order:
                # -------------------------
                # CORRECTED PATTERNS
                # -------------------------
                pattern_1 = f"{log_dir}/*_chr{chrom}_{step}"
                pattern_2 = f"{log_dir}/*_{chrom}_harmonize_sample_sizes.log"
                patterns = [pattern_1, pattern_2]
                # -------------------------
                # COLLECT MATCHED FILES
                # -------------------------
                matched_files = []
                for p in patterns:
                    matched_files.extend(glob.glob(p))
                # -------------------------
                # WRITE LOG CONTENT
                # -------------------------
                for log_file in matched_files:
                    base = os.path.basename(log_file)
                    outfile.write(f"\n--- BEGIN {base} ---\n")
                    try:
                        with open(log_file, "r") as f:
                            outfile.write(f.read())
                    except Exception as e:
                        outfile.write(f"[ERROR READING FILE: {e}]")
                    outfile.write(f"\n--- END {base} ---\n")
                    # delete original log
                    Path(log_file).unlink(missing_ok=True)
