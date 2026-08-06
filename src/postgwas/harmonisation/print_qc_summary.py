def print_qc_summary(final_qc_df, columns=None, warn_threshold=0.05):
    """
    Print GWAS QC summary in a human-readable, biologically meaningful format.
    Warnings are shown ONLY when thresholds are exceeded.
    """

    PREFIX = "\t\t"

    if isinstance(columns, str):
        columns = [columns]

    if columns is None:
        columns = ["raw_variant_count"]

    # Convert long → index-based
    final_qc_df = final_qc_df.set_index("index")

    QC_LABELS = {
        "num_records": "Total variant records in vcf file",
        "num_snps": "Total SNPs",
        "ts_tv_ratio": "Transition / Transversion ratio",
        "variants_with_missing_extrnal_af": "Variants missing external allele frequency",
        "variants_with_EAF_diff_gt_cutoff": "Variants with large EAF discrepancy & missing external allele frequency",
        "total_variant_with_missing_eaf": "Total variant with missing EAF",
        "total_variant_with_invalid_beta_se": "Total variant with invalid beta, SE, or Z-score",
    }

    available_cols = [c for c in columns if c in final_qc_df.columns]
    if not available_cols:
        print(f"{PREFIX}❌ No valid QC columns found.")
        return

    for col in available_cols:
        print(f"\n{PREFIX}📊 GWAS QC Summary — {col}")
        print(f"{PREFIX}" + "-" * (26 + len(col)))

        # Use num_records as denominator
        if "num_records" not in final_qc_df.index:
            print(f"{PREFIX}❌ Missing num_records — cannot compute QC ratios.")
            continue

        total_variants = final_qc_df.loc["num_records", col]

        def frac(val):
            return (val / total_variants) if total_variants and val is not None else 0.0

        # --------------------------------------------------
        # Core QC metrics
        # --------------------------------------------------
        for key, label in QC_LABELS.items():
            clean_key = key.strip()
            if clean_key not in final_qc_df.index:
                continue
            val = final_qc_df.loc[clean_key, col]

            if key == "ts_tv_ratio":
                print(f"{PREFIX}🔬 {label:45s}: {val:.2f}")
            else:
                print(f"{PREFIX}🧬 {label:45s}: {int(val):,}")

        # --------------------------------------------------
        # QC Warnings (only meaningful ones)
        # --------------------------------------------------
        warning_lines = []

        # Missing external AF
        if "variants_with_missing_extrnal_af" in final_qc_df.index:
            ratio = frac(final_qc_df.loc["variants_with_missing_extrnal_af", col])
            if ratio > warn_threshold:
                warning_lines.append(
                    f"{PREFIX}❗ {QC_LABELS['variants_with_missing_extrnal_af']}: "
                    f"{ratio*100:.2f}% "
                    f"({int(final_qc_df.loc['variants_with_missing_extrnal_af', col]):,} variants)"
                )
         
        # Variants with missing external AF
        if "total_variant_with_missing_eaf" in final_qc_df.index:
            ratio = frac(final_qc_df.loc["total_variant_with_missing_eaf", col])
            if ratio > warn_threshold:
                warning_lines.append(
                    f"{PREFIX}❗ {QC_LABELS['total_variant_with_missing_eaf']}: "
                    f"{ratio*100:.2f}% "
                    f"({int(final_qc_df.loc['total_variant_with_missing_eaf', col]):,} variants)"
                )

        # Large EAF discordance (allele flip warning)
        YELLOW = "\033[93m"
        RESET = "\033[0m"
        if "variants_with_EAF_diff_gt_cutoff" in final_qc_df.index:
            ratio = frac(final_qc_df.loc["variants_with_EAF_diff_gt_cutoff", col])
            if ratio > warn_threshold:
                warning_lines.append(
                    f"{PREFIX}❗ {QC_LABELS['variants_with_EAF_diff_gt_cutoff']}: "
                    f"{ratio*100:.2f}% "
                    f"({int(final_qc_df.loc['variants_with_EAF_diff_gt_cutoff', col]):,} variants)\n"
                    f"{PREFIX}   → {YELLOW}Possible allele flips, strand mismatches, "
                    f"REF/ALT mismatches (e.g. indels), or population mismatch{RESET}"
                )
        # --------------------------------------------------
        # Print warnings only if present
        # --------------------------------------------------
        if warning_lines:
            print(f"\n{PREFIX}⚠️ QC Warnings")
            print(f"{PREFIX}" + "-" * 13)
            for line in warning_lines:
                print(line)
    
