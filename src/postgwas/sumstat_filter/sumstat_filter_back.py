import os
import subprocess
from typing import Optional, Dict, Any
import os, subprocess, io
from pathlib import Path
from postgwas.utils.main import run_cmd



def filter_gwas_vcf_bcftools(
    vcf_path: str,
    output_folder: str,
    output_prefix: str,
    pval_cutoff: Optional[float] = None,
    maf_cutoff: Optional[float] = None,
    allelefreq_diff_cutoff: Optional[float] = None,
    info_cutoff: Optional[float] = None,
    external_af_name: str = "EUR",
    include_indels: bool = True,
    exclude_palindromic: bool = False,
    palindromic_af_lower: float = 0.4,
    palindromic_af_upper: float = 0.6,
    remove_mhc: bool = False,
    mhc_chrom: str = "6",
    mhc_start: int = 25000000,
    mhc_end: int = 34000000,
    threads: int = 5,
    max_mem: str = "5G",
    extra_options: Optional[Dict[str, Any]] = None,
) -> Dict[str, str]:
    BOLD = "\033[1m"
    RESET = "\033[0m"
    
    vcf_name = os.path.basename(vcf_path)
    genomeversion = None
    if "GRCh38" in vcf_name:
        genomeversion = "GRCh38"
    elif "GRCh37" in vcf_name:
        genomeversion = "GRCh37"
    output_prefix = (
        f"{output_prefix}_{genomeversion}"
        if genomeversion
        else output_prefix)
    
    # ============================================================
    # Create log directory + log buffer
    # ============================================================
    step_dir = Path(output_folder)
    step_dir.mkdir(parents=True, exist_ok=True)

    log_dir = step_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{output_prefix}_filter_gwas_vcf_bcftools.log"

    log_buffer = io.StringIO()
    def log_print(*args):
        log_buffer.write(" ".join(str(a) for a in args) + "\n")

    # ============================================================
    # Header: Start
    # ============================================================
    if not os.path.exists(vcf_path):
        msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
        log_print(msg)
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())
        raise FileNotFoundError(msg)

    os.makedirs(output_folder, exist_ok=True)
    output_vcf = os.path.join(output_folder, f"{output_prefix}_filtered.vcf.gz")

    # ============================================================
    # Count variants before filtering
    # ============================================================
    pre_cmd = f"bcftools view --threads {threads} -H {vcf_path} | wc -l"
    pre = run_cmd(pre_cmd)
    try:
        pre_variants = int(pre.stdout.strip())
    except Exception:
        pre_variants = None
    
    if extra_options is not None:
        log_print(f"\n\t\t\t📊 Variants in the input file                        : {extra_options.get('total_variant_infile', 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants successfully read by harmonisation module: {extra_options.get('total_variant_read', 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants USED for VCF creation                    : {extra_options.get('total_variant_in_vcf_input', 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants BEFORE filtering                         : {pre_variants}")
        log_print("")

    # ============================================================
    # Build inclusion expression
    # ============================================================
    include_parts = []
    if pval_cutoff is not None:
        include_parts.append(f"(FORMAT/LP >= {pval_cutoff})")
    if maf_cutoff is not None:
        include_parts.append(f"(FORMAT/AF >= {maf_cutoff} & FORMAT/AF <= {1 - maf_cutoff})")
    if info_cutoff is not None:
        include_parts.append(f"(FORMAT/SI >= {info_cutoff})")
    if allelefreq_diff_cutoff is not None:
        include_parts.append(f"(abs(INFO/AF - INFO/{external_af_name}) <= {allelefreq_diff_cutoff})")

    include_expr = " & ".join(include_parts) if include_parts else "1"

    log_print("🔧 Include expression (bcftools -i):")
    log_print(include_expr)
    log_print("")

    # ============================================================
    # Build exclusion expression
    # ============================================================
    palindromic_logic = (
        '((REF=="A" & ALT=="T") | (REF=="T" & ALT=="A") | '
        '(REF=="C" & ALT=="G") | (REF=="G" & ALT=="C"))'
    )

    exclude_expr = None
    if exclude_palindromic:
        exclude_expr = (
            f"({palindromic_logic} & "
            f"(FORMAT/AF >= {palindromic_af_lower} & FORMAT/AF <= {palindromic_af_upper}))"
        )

    log_print("🔧 Exclude expression (bcftools -e):")
    log_print(str(exclude_expr))
    log_print("")

    # ============================================================
    # Build bcftools pipeline
    # ============================================================
    cmd_parts = []
    cmd_parts.append(f"bcftools view --threads {threads} -i '{include_expr}' {vcf_path}")

    if not include_indels:
        cmd_parts.append(f"bcftools view --threads {threads} --types snps")

    if exclude_expr:
        cmd_parts.append(f"bcftools view --threads {threads} -e '{exclude_expr}'")

    if remove_mhc:
        mhc_bed = os.path.join(output_folder, f"{output_prefix}_mhc_exclude.bed")
        with open(mhc_bed, "w") as f:
            f.write(f"{mhc_chrom}\t{mhc_start}\t{mhc_end}\n")
        cmd_parts.append(f"bcftools view --threads {threads} -T ^{mhc_bed}")

    cmd_parts.append(f"bcftools sort --temp-dir {output_folder} --max-mem {max_mem}")
    cmd_parts.append("bgzip -c")

    pipeline = " | ".join(cmd_parts) + f" > {output_vcf} && tabix -f -p vcf {output_vcf}"

    log_print("🚀 Full bcftools pipeline:")
    log_print(pipeline)
    log_print("")

    result = run_cmd(pipeline)
    if result.returncode != 0:
        log_print("❌ Error during bcftools pipeline execution.")
        log_print("STDERR:")
        log_print(result.stderr)
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())
        raise RuntimeError("bcftools pipeline failed")

    # ============================================================
    # Count variants after filtering
    # ============================================================
    post_cmd = f"bcftools view --threads {threads} -H {output_vcf} | wc -l"
    post = run_cmd(post_cmd)
    try:
        post_variants = int(post.stdout.strip())
    except Exception:
        post_variants = None

    log_print("📊 Variants filtering command used:")
    log_print(post_cmd)
    log_print(f"✅ Variants AFTER filtering                          : {post_variants}")
    log_print(f"💾 Filtered VCF saved to                             : {output_vcf}")
    log_print("\n🎉 bcftools filtering completed.")

    # ============================================================
    # Write log file
    # ============================================================
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    print("\t\t\t🧪 Variant filtering summary:")
    print("\t\t\t🔹 Variants are RETAINED if they pass the following active filters:")
    if pval_cutoff is not None:
        print(f"\t\t\t   • P-value evidence: LP ≥ {pval_cutoff}")
    if maf_cutoff is not None:
        print(f"\t\t\t   • Minor allele frequency range: {maf_cutoff} ≤ AF ≤ {1 - maf_cutoff}")
    if info_cutoff is not None:
        print(f"\t\t\t   • Imputation quality: INFO (SI) ≥ {info_cutoff}")
    if allelefreq_diff_cutoff is not None:
        print(f"\t\t\t   • External AF concordance: |AF − {external_af_name}| ≤ {allelefreq_diff_cutoff}")

    print(f"\t\t\t   • Allowed variant types: {'SNPs and indels' if include_indels else 'SNPs only'}")

    print("\n\t\t\t🔻 Variants are REMOVED if they meet any of the following conditions:")
    if exclude_palindromic:
        print(f"\t\t\t   • Palindromic SNPs with ambiguous frequency ({palindromic_af_lower} ≤ AF ≤ {palindromic_af_upper})")
    else:
        print("\t\t\t   • Palindromic SNPs: NOT removed")

    if remove_mhc:
        print(f"\t\t\t   • Located in MHC region (chr{mhc_chrom}:{mhc_start:,}-{mhc_end:,})")
    else:
        print("\t\t\t   • MHC region variants: NOT removed")

    print("")
    if extra_options is not None:
        print(f"\t\t\t📊 Variants in the input file                                     : {extra_options.get('total_variant_infile', 0):,}")
        print(f"\t\t\t📊 Variants successfully read by harmonisation module             : {extra_options.get('total_variant_read', 0):,}")
        print(f"\t\t\t📊 Variants USED for VCF creation                                 : {extra_options.get('total_variant_in_vcf_input', 0):,}")
    
    print(f"\t\t\t📊 Variants IN the VCF FILE BEFORE filtering                      : {pre_variants if pre_variants is not None else 0:,}")
    print(f"\t\t\t📊 Variants IN the VCF FILE AFTER filtering                       : {post_variants if post_variants is not None else 0:,}")

    if pre_variants is not None and post_variants is not None:
        print(f"\t\t\t🧮 Variants REMOVED during filtering                              : {pre_variants - post_variants:,}")
    
    print("")

    # --- CRITICAL FIX: Safe access for extra_options ---
    if extra_options is not None:
        infile = extra_options.get("total_variant_infile", 0)
        read = extra_options.get("total_variant_read", 0)
        vcf_input = extra_options.get("total_variant_in_vcf_input", 0)

        # 1️⃣ Input → Read
        if read < infile:
            diff = infile - read
            pct = (diff / infile * 100) if infile else 0
            msg = f"{BOLD}⚠️ Harmonisation did not read all variants: {diff:,} variants missing ({pct:.2f}%).{RESET}"
            log_print(msg)
            print(msg)

        # 2️⃣ Read → VCF input
        if vcf_input < read:
            diff = read - vcf_input
            pct = (diff / read * 100) if read else 0
            log_print(
                f"{BOLD}\t\t⚠️ Fewer variants used for VCF creation: "
                f"{diff:,} removed after harmonisation ({pct:.2f}%).{RESET}"
            )
            indent = "\t" * 4
            print(
                f"\n{indent}{BOLD}⚠️ ⚠️ VARIANT LOSS DURING HARMONISATION WARNING{RESET}\n"
                f"{indent}{'-'*60}\n"
                f"{indent}• Variants read from sumstats file: {read:,}\n"
                f"{indent}• Variants used for VCF creation  : {vcf_input:,}\n"
                f"{indent}• Removed after harmonisation     : {diff:,} ({pct:.2f}%)\n" 
                f"\n"
                f"{indent}Possible reasons variants in sumstat file may be :\n"
                f"{indent}    ├─ Missing in INFO file\n"
                f"{indent}    ├─ Missing in EAF file\n"
                f"{indent}Or Removed during QC:\n"
                f"{indent}         • Null BETA / SE\n"
                f"{indent}         • BETA = 0\n"
                f"{indent}         • SE ≤ 0\n"
                f"{indent}{'-'*60}\n" 
            )
        
        # 3️⃣ VCF input → Pre-filter
        if pre_variants is not None and pre_variants < vcf_input:
            diff = vcf_input - pre_variants
            pct = (diff / vcf_input * 100) if vcf_input else 0
            log_print(
                f"{BOLD}⚠️ Not all variants carried into VCF creation: "
                f"{diff:,} dropped during VCF generation ({pct:.2f}%).{RESET}"
            )

    print("")
    return {"filtered_vcf": output_vcf}
    

# def filter_gwas_vcf_bcftools(
#     vcf_path: str,
#     output_folder: str,
#     output_prefix: str,
#     pval_cutoff: Optional[float] = None,
#     maf_cutoff: Optional[float] = None,
#     allelefreq_diff_cutoff: Optional[float] = None,
#     info_cutoff: Optional[float] = None,
#     external_af_name: str = "EUR",
#     include_indels: bool = True,
#     exclude_palindromic: bool = False,
#     palindromic_af_lower: float = 0.4,
#     palindromic_af_upper: float = 0.6,
#     remove_mhc: bool = False,
#     mhc_chrom: str = "6",
#     mhc_start: int = 25000000,
#     mhc_end: int = 34000000,
#     threads: int = 5,
#     max_mem: str = "5G",
#     extra_options: Optional[Dict[str, Any]] = None,
# ) -> str:
#     BOLD = "\033[1m"
#     RESET = "\033[0m"
    
#     vcf_name = os.path.basename(vcf_path)
#     genomeversion = None
#     if "GRCh38" in vcf_name:
#         genomeversion = "GRCh38"
#     elif "GRCh37" in vcf_name:
#         genomeversion = "GRCh37"
#     output_prefix = (
#         f"{output_prefix}_{genomeversion}"
#         if genomeversion
#         else output_prefix)
    
#     # ============================================================
#     # Create log directory + log buffer
#     # ============================================================
#     step_dir = Path(output_folder)
#     step_dir.mkdir(parents=True, exist_ok=True)

#     log_dir = step_dir / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{output_prefix}_filter_gwas_vcf_bcftools.log"

#     log_buffer = io.StringIO()
#     def log_print(*args):
#         log_buffer.write(" ".join(str(a) for a in args) + "\n")

#     # ============================================================
#     # Header: Start
#     # ============================================================
#     if not os.path.exists(vcf_path):
#         msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
#         log_print(msg)
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())
#         raise FileNotFoundError(msg)

#     os.makedirs(output_folder, exist_ok=True)
#     output_vcf = os.path.join(output_folder, f"{output_prefix}_filtered.vcf.gz")

#     # ============================================================
#     # Count variants before filtering
#     # ============================================================
#     pre_cmd = f"bcftools view --threads {threads} -H {vcf_path} | wc -l"
#     pre = run_cmd(pre_cmd)
#     try:
#         pre_variants = int(pre.stdout.strip())
#     except Exception:
#         pre_variants = None
#     if extra_options !=None:
#         log_print(f"\n\t\t\t📊 Variants in the input file                        : {extra_options['total_variant_infile']}")
#         log_print(f"\n\t\t\t📊 Variants successfully read by harmonisation module: {extra_options['total_variant_read']}")
#         log_print(f"\n\t\t\t📊 Variants USED for VCF creation                    : {extra_options['total_variant_in_vcf_input']}")
#         log_print(f"\n\t\t\t📊 Variants BEFORE filtering                         : {pre_variants}")
#         log_print("")
#     # ============================================================
#     # Build inclusion expression
#     # ============================================================
#     include_parts = []

#     if pval_cutoff is not None:
#         include_parts.append(f"(FORMAT/LP >= {pval_cutoff})")

#     if maf_cutoff is not None:
#         include_parts.append(
#             f"(FORMAT/AF >= {maf_cutoff} & FORMAT/AF <= {1 - maf_cutoff})"
#         )

#     if info_cutoff is not None:
#         include_parts.append(f"(FORMAT/SI >= {info_cutoff})")

#     if allelefreq_diff_cutoff is not None:
#         include_parts.append(
#             f"(abs(INFO/AF - INFO/{external_af_name}) <= {allelefreq_diff_cutoff})"
#         )

#     include_expr = " & ".join(include_parts) if include_parts else "1"

#     log_print("🔧 Include expression (bcftools -i):")
#     log_print(include_expr)
#     log_print("")

#     # ============================================================
#     # Build exclusion expression
#     # ============================================================
#     palindromic_logic = (
#         '((REF=="A" & ALT=="T") | (REF=="T" & ALT=="A") | '
#         '(REF=="C" & ALT=="G") | (REF=="G" & ALT=="C"))'
#     )

#     exclude_expr = None
#     if exclude_palindromic:
#         exclude_expr = (
#             f"({palindromic_logic} & "
#             f"(FORMAT/AF >= {palindromic_af_lower} & FORMAT/AF <= {palindromic_af_upper}))"
#         )

#     log_print("🔧 Exclude expression (bcftools -e):")
#     log_print(str(exclude_expr))
#     log_print("")

#     # ============================================================
#     # Build bcftools pipeline
#     # ============================================================
#     cmd_parts = []

#     # 1) include expression
#     cmd_parts.append(
#         f"bcftools view --threads {threads} -i '{include_expr}' {vcf_path}"
#     )

#     # 2) remove indels if required
#     if not include_indels:
#         cmd_parts.append(f"bcftools view --threads {threads} --types snps")

#     # 3) exclude palindromic SNPs if required
#     if exclude_expr:
#         cmd_parts.append(f"bcftools view --threads {threads} -e '{exclude_expr}'")

#     # 4) remove MHC region if required
#     if remove_mhc:
#         mhc_bed = os.path.join(output_folder, f"{output_prefix}_mhc_exclude.bed")
#         with open(mhc_bed, "w") as f:
#             f.write(f"{mhc_chrom}\t{mhc_start}\t{mhc_end}\n")
#         cmd_parts.append(f"bcftools view --threads {threads} -T ^{mhc_bed}")

#     # 5) sort + compress
#     cmd_parts.append(
#         f"bcftools sort --temp-dir {output_folder} --max-mem {max_mem}"
#     )
#     cmd_parts.append("bgzip -c")

#     # final pipeline
#     pipeline = (
#         " | ".join(cmd_parts)
#         + f" > {output_vcf} && tabix -f -p vcf {output_vcf}"
#     )

#     log_print("🚀 Full bcftools pipeline:")
#     log_print(pipeline)
#     log_print("")

#     # ============================================================
#     # Run pipeline
#     # ============================================================
#     result = run_cmd(pipeline)
#     if result.returncode != 0:
#         log_print("❌ Error during bcftools pipeline execution.")
#         log_print("STDERR:")
#         log_print(result.stderr)
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())
#         raise RuntimeError("bcftools pipeline failed")

#     # ============================================================
#     # Count variants after filtering
#     # ============================================================
#     post_cmd = f"bcftools view --threads {threads} -H {output_vcf} | wc -l"
#     post = run_cmd(post_cmd)
#     try:
#         post_variants = int(post.stdout.strip())
#     except Exception:
#         post_variants = None

#     log_print("📊 Variants filtering command used:")
#     log_print(post_cmd)
#     log_print(f"✅ Variants AFTER filtering                          : {post_variants}")
#     log_print(f"💾 Filtered VCF saved to                             : {output_vcf}")
#     log_print("\n🎉 bcftools filtering completed.")

#     # ============================================================
#     # Write log file
#     # ============================================================
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     print("\t\t\t🧪 Variant filtering summary:")
#     # ----------------------------
#     # Inclusion (retention)
#     # ----------------------------
#     print("\t\t\t🔹 Variants are RETAINED if they pass the following active filters:")
#     if pval_cutoff is not None:
#         print(f"\t\t\t   • P-value evidence: LP ≥ {pval_cutoff}")
#     if maf_cutoff is not None:
#         print(f"\t\t\t   • Minor allele frequency range: {maf_cutoff} ≤ AF ≤ {1 - maf_cutoff}")

#     if info_cutoff is not None:
#         print(f"\t\t\t   • Imputation quality: INFO (SI) ≥ {info_cutoff}")

#     if allelefreq_diff_cutoff is not None:
#         print(
#             f"\t\t\t   • External AF concordance: "
#             f"|AF − {external_af_name}| ≤ {allelefreq_diff_cutoff}"
#         )

#     print(
#         f"\t\t\t   • Allowed variant types: "
#         f"{'SNPs and indels' if include_indels else 'SNPs only'}"
#     )

#     # ----------------------------
#     # Exclusions
#     # ----------------------------
#     print("\n\t\t\t🔻 Variants are REMOVED if they meet any of the following conditions:")

#     if exclude_palindromic:
#         print(
#             f"\t\t\t   • Palindromic SNPs with ambiguous frequency "
#             f"({palindromic_af_lower} ≤ AF ≤ {palindromic_af_upper})"
#         )
#     else:
#         print("\t\t\t   • Palindromic SNPs: NOT removed")

#     if remove_mhc:
#         print(
#             f"\t\t\t   • Located in MHC region "
#             f"(chr{mhc_chrom}:{mhc_start:,}-{mhc_end:,})"
#         )
#     else:
#         print("\t\t\t   • MHC region variants: NOT removed")

#     # ----------------------------
#     # Counts
#     # ----------------------------
#     print("")
#     if extra_options !=None:
#         print(f"\t\t\t📊 Variants in the input file                                     : {extra_options['total_variant_infile']:,}")
#         print(f"\t\t\t📊 Variants successfully read by harmonisation module             : {extra_options['total_variant_read']:,}")
#         print(f"\t\t\t📊 Variants USED for VCF creation                                 : {extra_options['total_variant_in_vcf_input']:,}")
#     print(f"\t\t\t📊 Variants IN the VCF FILE BEFORE filtering                      : {pre_variants:,}")
#     print(f"\t\t\t📊 Variants IN the VCF FILE AFTER filtering                       : {post_variants:,}")

#     if pre_variants is not None and post_variants is not None:
#         print(f"\t\t\t🧮 Variants REMOVED during filtering                              : {pre_variants - post_variants}")
    
#     print("")

#     infile = extra_options["total_variant_infile"]
#     read = extra_options["total_variant_read"]
#     vcf_input = extra_options["total_variant_in_vcf_input"]

#     # 1️⃣ Input → Read
#     if read < infile:
#         diff = infile - read
#         pct = (diff / infile * 100) if infile else 0
#         log_print(
#             f"{BOLD}⚠️ Harmonisation did not read all variants: "
#             f"{diff:,} variants missing ({pct:.2f}%).{RESET}"
#         )
#         print(
#             f"{BOLD}⚠️ Harmonisation did not read all variants: "
#             f"{diff:,} variants missing ({pct:.2f}%).{RESET}"
#         )
#     # 2️⃣ Read → VCF input
#     if vcf_input < read:
#         diff = read - vcf_input
#         pct = (diff / read * 100) if read else 0
#         log_print(
#             f"{BOLD}\t\t⚠️ Fewer variants used for VCF creation: "
#             f"{diff:,} removed after harmonisation ({pct:.2f}%).{RESET}"
#             f"{read:,} variants were read from the sumstat file by the harmonisation module; "
#             f"only {vcf_input:,} variants were used as input for VCF creation by gwas2vcf.{RESET}"
#         )
#         log_print(
#             "   Possible reasons: variants missing in INFO/EAF file or removed during QC."
#         )
#         log_print(
#             "   See: 00_harmonised_sumstat/qc_summary/*_gwas2vcf_summary.tsv"
#         )
#         indent = "\t" * 4

#         print(
#             f"\n{indent}{BOLD}⚠️ ⚠️ VARIANT LOSS DURING HARMONISATION WARNING{RESET}\n"
#             f"{indent}{'-'*60}\n"
#             f"{indent}• Variants read from sumstats file: {read:,}\n"
#             f"{indent}• Variants used for VCF creation  : {vcf_input:,}\n"
#             f"{indent}• Removed after harmonisation     : {diff:,} ({pct:.2f}%)\n" 
#             f"\n"
#             f"{indent}Possible reasons variants in sumstat file may be :\n"
#             f"{indent}    ├─ Missing in INFO file\n"
#             f"{indent}    ├─ Missing in EAF file\n"
#             f"{indent}Or Removed during QC:\n"
#             f"{indent}         • Null BETA / SE\n"
#             f"{indent}         • BETA = 0\n"
#             f"{indent}         • SE ≤ 0\n"
#             f"{indent}{'-'*60}\n" 
#         )
    
#     # 3️⃣ VCF input → Pre-filter
#     if pre_variants < vcf_input:
#         diff = vcf_input - pre_variants
#         pct = (diff / vcf_input * 100) if vcf_input else 0
#         log_print(
#             f"{BOLD}⚠️ Not all variants carried into VCF creation: "
#             f"{diff:,} dropped during VCF generation ({pct:.2f}%).{RESET}"
#         )

#     print("")
#     return {
#         "filtered_vcf": output_vcf
#     }
