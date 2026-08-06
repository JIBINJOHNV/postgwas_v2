#!/usr/bin/env python3
import argparse
import sys
import os

from postgwas.qc_summary.main import run_full_vcf_qc_and_filter


def main():
    parser = argparse.ArgumentParser(
        description="Run full QC + filtering on a GWAS summary VCF."
    )

    # Required arguments
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to the input GWAS summary VCF (.vcf.gz)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output folder where filtered VCF and logs will be stored"
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files (e.g., myGWAS)"
    )

    # Optional filtering parameters
    parser.add_argument("--pval-cutoff", type=float, default=None,
                        help="Minimum -log10(P-value) required (FORMAT/LP). Default: None")
    parser.add_argument("--maf-cutoff", type=float, default=None,
                        help="MAF filter: keep AF between [maf, 1-maf]. Default: None")
    parser.add_argument("--allelefreq-diff-cutoff", type=float, default=0.2,
                        help="Max allowed |AF - external_AF| difference. Default: 0.2")
    parser.add_argument("--info-cutoff", type=float, default=0.7,
                        help="Min imputation INFO score (FORMAT/SI). Default: 0.7")
    parser.add_argument("--external-af-name", type=str, default="EUR",
                        help="External AF used for concordance filtering. Default: EUR")

    # Flags / toggles
    parser.add_argument("--no-indels", action="store_true",
                        help="Exclude indels (SNPs only). Default: include all")
    parser.add_argument("--no-palindromic", action="store_true",
                        help="Exclude palindromic SNPs with AF in [0.4, 0.6]")

    # Palindromic AF range
    parser.add_argument("--pal-lower", type=float, default=0.4,
                        help="Lower bound for palindromic AF exclusion. Default: 0.4")
    parser.add_argument("--pal-upper", type=float, default=0.6,
                        help="Upper bound for palindromic AF exclusion. Default: 0.6")

    # MHC region
    parser.add_argument("--remove-mhc", action="store_true",
                        help="Remove MHC region (default: off)")
    parser.add_argument("--mhc-chrom", default="6",
                        help="Chromosome for MHC region. Default: 6")
    parser.add_argument("--mhc-start", type=int, default=25000000,
                        help="Start of MHC region. Default: 25 Mb")
    parser.add_argument("--mhc-end", type=int, default=34000000,
                        help="End of MHC region. Default: 34 Mb")

    # Threads + memory
    parser.add_argument("--threads", type=int, default=5,
                        help="Number of threads for bcftools. Default: 5")
    parser.add_argument("--max-mem", type=str, default="5G",
                        help="Memory for bcftools sort --max-mem. Default: 5G")

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Input checks
    # -------------------------------------------------------------------------
    if not os.path.exists(args.vcf):
        print(f"❌ ERROR: Input VCF not found: {args.vcf}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # -------------------------------------------------------------------------
    # Run full QC
    # -------------------------------------------------------------------------
    try:
        output_vcf = run_full_vcf_qc_and_filter(
            vcf_path=args.vcf,
            output_folder=args.outdir,
            output_prefix=args.prefix,
            pval_cutoff=args.pval_cutoff,
            maf_cutoff=args.maf_cutoff,
            allelefreq_diff_cutoff=args.allelefreq_diff_cutoff,
            info_cutoff=args.info_cutoff,
            external_af_name=args.external_af_name,
            include_indels=not args.no_indels,
            include_palindromic=not args.no_palindromic,
            palindromic_af_lower=args.pal_lower,
            palindromic_af_upper=args.pal_upper,
            remove_mhc=args.remove_mhc,
            mhc_chrom=args.mhc_chrom,
            mhc_start=args.mhc_start,
            mhc_end=args.mhc_end,
            threads=args.threads,
            max_mem=args.max_mem,
        )

    except Exception as e:
        print(f"❌ QC/Filtering failed: {e}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # CLI prints ONLY the final file path
    # -------------------------------------------------------------------------
    print(output_vcf)


if __name__ == "__main__":
    main()
