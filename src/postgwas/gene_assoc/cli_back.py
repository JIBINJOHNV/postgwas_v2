#!/usr/bin/env python3
"""
CLI for MAGMA full analysis pipeline.
"""

import argparse
from pathlib import Path
from postgwas.gene_tests.magma_main import magma_analysis_pipeline


def main():

    parser = argparse.ArgumentParser(
        prog="postgwas magma",
        description="Run full MAGMA gene & geneset analysis pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # ---------------------------------------------------------
    # Required input files
    # ---------------------------------------------------------
    parser.add_argument("--output_dir", required=True,
                        help="Folder where MAGMA outputs will be stored")

    parser.add_argument("--sample_id", required=True,
                        help="Sample / GWAS prefix (used for output filenames)")

    parser.add_argument("--ld_ref", required=True,
                        help="Prefix of PLINK reference panel for MAGMA (e.g., 1000G EUR)")

    parser.add_argument("--gene_loc_file", required=True,
                        help="MAGMA gene location file (.loc)")

    parser.add_argument("--snp_loc_file", required=True,
                        help="SNP location file generated from GWAS or VCF")

    parser.add_argument("--pval_file", required=True,
                        help="Input file containing SNP p-values for MAGMA")

    parser.add_argument("--geneset_file", required=True,
                        help="GMT file (MSigDB or custom pathways)")

    parser.add_argument("--log_file", required=True,
                        help="Log file prefix for saving pipeline logs")

    # ---------------------------------------------------------
    # Analysis parameters
    # ---------------------------------------------------------
    parser.add_argument("--threads", type=int, default=6,
                        help="Number of parallel MAGMA workers")

    parser.add_argument("--num_batches", type=int, default=6,
                        help="Number of MAGMA batches")

    parser.add_argument("--window_upstream", type=int, default=35,
                        help="MAGMA upstream annotation window (kb)")

    parser.add_argument("--window_downstream", type=int, default=10,
                        help="MAGMA downstream annotation window (kb)")

    parser.add_argument("--gene_model", default="snp-wise=mean",
                        help="MAGMA gene model (e.g., snp-wise=mean, top, all)")

    parser.add_argument("--n_sample_col", default="N_COL",
                        help="Column name with sample sizes for MAGMA pval input")

    parser.add_argument("--seed", type=int, default=10,
                        help="Random seed for MAGMA")

    parser.add_argument("--magma", default="magma",
                        help="Path to MAGMA binary")

    args = parser.parse_args()

    # ---------------------------------------------------------
    # Run pipeline
    # ---------------------------------------------------------
    print("\n🚀 Starting MAGMA full analysis pipeline...\n")

    result_df = magma_analysis_pipeline(
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        ld_ref=args.ld_ref,
        gene_loc_file=args.gene_loc_file,
        snp_loc_file=args.snp_loc_file,
        pval_file=args.pval_file,
        geneset_file=args.geneset_file,
        log_file=args.log_file,
        threads=args.threads,
        num_batches=args.num_batches,
        window_upstream=args.window_upstream,
        window_downstream=args.window_downstream,
        gene_model=args.gene_model,
        n_sample_col=args.n_sample_col,
        seed=args.seed,
        magma=args.magma,
    )

    print("\n🎉 MAGMA pipeline completed successfully.")
    print(f"Output saved inside: {args.output_dir}\n")

    return 0


if __name__ == "__main__":
    main()
