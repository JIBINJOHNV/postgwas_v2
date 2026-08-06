import polars as pl
import os
import subprocess

# --- Constants & Thresholds ---
LEAD_P = 5e-8        # P-value threshold for signal detection
GWAS_P = 0.05        # P-value for candidate SNPs in the locus
R2_CLUMP = 0.6       # Strict LD for independent signals
R2_LEAD = 0.1        # Loose LD for grouping Lead SNPs
MERGE_DIST = 250_000 # 250kb distance for physical merging

# --- Function 1: LD Retrieval (System Tabix) ---
def get_ld_partners(chrom, pos, rsid, ld_file_path, r2_threshold):
    """
    Returns a dictionary of {partner_id: r2_value} for SNPs in LD with rsid.
    """
    # Initialize with the SNP itself at r2=1.0
    partners_r2 = {rsid: 1.0}
    start, end = pos - 1000000, pos + 1000000
    region = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ['tabix', ld_file_path, region], 
            capture_output=True, text=True, check=True
        )
        if not result.stdout.strip():
            return partners_r2
        ld_df = pl.read_csv(
            result.stdout.encode(), 
            separator="\t", 
            has_header=False,
            new_columns=["chr_a", "pos_a", "snp_a", "chr_b", "pos_b", "snp_b", "r2"]
        )
        # Filter for rows involving the rsid and meeting the r2 threshold
        # We check both snp_a and snp_b because PLINK only stores the pair once
        filtered = ld_df.filter(
            (pl.col("r2") >= r2_threshold) & 
            ((pl.col("snp_a") == rsid) | (pl.col("snp_b") == rsid))
        )
        if not filtered.is_empty():
            # If rsid is snp_a, partner is snp_b. If rsid is snp_b, partner is snp_a.
            partners_df = filtered.select([
                pl.when(pl.col("snp_a") == rsid)
                .then(pl.col("snp_b"))
                .otherwise(pl.col("snp_a"))
                .alias("partner_id"),
                pl.col("r2")
            ])
            # Convert to dictionary for fast lookup
            # result: {'partner_snp_id': 0.85, ...}
            partners_r2.update(dict(zip(partners_df["partner_id"], partners_df["r2"])))
    except subprocess.CalledProcessError:
        pass
    return partners_r2

# --- Function 2: First Clumping (Independent Significant SNPs) ---
def find_ind_sig_snps(chr_df, ld_path):
    """
    Identifies Ind. sig. SNPs (P <= LEAD_P, r2 >= 0.6).
    Ensures each candidate SNP is assigned to ONLY the strongest Ind. Sig SNP.
    """
    # 1. This pool is for picking the NEXT Lead SNP
    remaining_leads = chr_df.filter(pl.col("pcol") <= LEAD_P).sort("pcol")
    # 2. This pool is for extracting Candidate SNPs (will be depleted)
    available_pool = chr_df.clone()
    ind_sig_clumps = [] 
    while not remaining_leads.is_empty():
        # Pick the strongest remaining SNP as the Ind. sig. SNP
        top_snp_row = remaining_leads.row(0, named=True)
        ind_sig_id = top_snp_row['uniq_id']
        # Get LD partners from reference panel
        partners_r2_map = get_ld_partners(
            top_snp_row['chrcol'], top_snp_row['poscol'], ind_sig_id, ld_path, R2_CLUMP
        )
        n_ref_count = len(partners_r2_map)
        candidate_ids = list(partners_r2_map.keys())
        # --- CRITICAL FIX: Extract only from SNPs not already assigned ---
        candidate_snps_df = available_pool.filter(pl.col("uniq_id").is_in(candidate_ids))
        if not candidate_snps_df.is_empty():
            # Map R2 values
            r2_map_df = pl.DataFrame({
                "uniq_id": list(partners_r2_map.keys()),
                "r2_with_IndSig": list(partners_r2_map.values())
            })
            candidate_snps_df = candidate_snps_df.join(r2_map_df, on="uniq_id", how="left")
            # Tag with the current Ind. sig. SNP ID
            candidate_snps_df = candidate_snps_df.with_columns([
                pl.lit(ind_sig_id).alias("ind_sig_SNP_id"),
                pl.lit(n_ref_count).alias("n_refsnps")
            ])
            ind_sig_clumps.append(candidate_snps_df)
            # --- UPDATE POOLS: Remove assigned SNPs from BOTH pools ---
            assigned_ids = candidate_snps_df["uniq_id"].to_list()
            # 1. Remove from the candidate pool so they can't be in another cluster
            available_pool = available_pool.filter(~pl.col("uniq_id").is_in(assigned_ids))
            # 2. Remove from the lead pool so they can't be picked as a new lead
            remaining_leads = remaining_leads.filter(~pl.col("uniq_id").is_in(assigned_ids))
            print(f"Ind. sig. SNP: {ind_sig_id} | New Candidates: {candidate_snps_df.height} | Ref SNPs: {n_ref_count}")
        else:
            # If the top lead itself was already consumed as a candidate elsewhere
            remaining_leads = remaining_leads.filter(pl.col("uniq_id") != ind_sig_id)
    return pl.concat(ind_sig_clumps) if ind_sig_clumps else pl.DataFrame()


# --- Function 3: Second Clumping (Lead SNPs) ---
def find_lead_snps(ind_sig_df, ld_path):
    """
    Identifies Lead SNPs by clumping Ind. sig. SNPs at r2 < 0.1.
    Input: long-format DF from find_ind_sig_snps.
    """
    if ind_sig_df.is_empty():
        return pl.DataFrame()
    # 1. Extract only the unique Ind. sig. SNPs (where the SNP is its own lead)
    ind_sig_heads = ind_sig_df.filter(
        pl.col("uniq_id") == pl.col("ind_sig_SNP_id")
    ).sort("pcol")
    lead_snp_clusters = []
    remaining_ind_sig = ind_sig_heads
    while not remaining_ind_sig.is_empty():
        # Top Ind. sig. SNP becomes a Lead SNP
        lead_row = remaining_ind_sig.row(0, named=True)
        lead_id = lead_row['uniq_id']
        # Find other Ind. sig. SNPs in LD (r2 >= 0.1) with this Lead SNP
        partners_r2_map = get_ld_partners(
            lead_row['chrcol'], lead_row['poscol'], lead_id, ld_path, 0.1 # r2 < 0.1 for independence
        )
        partner_ids = list(partners_r2_map.keys())
        # Identify which Ind. sig. SNPs are being clumped into this Lead SNP
        clumped_ind_sigs = ind_sig_heads.filter(pl.col("uniq_id").is_in(partner_ids))
        # Map R2 and tag with Lead SNP ID
        r2_map_df = pl.DataFrame({
            "uniq_id": partner_ids,
            "r2_with_Lead": list(partners_r2_map.values())
        })
        clumped_ind_sigs = clumped_ind_sigs.join(r2_map_df, on="uniq_id", how="left")
        clumped_ind_sigs = clumped_ind_sigs.with_columns(
            pl.lit(lead_id).alias("lead_SNP_id")
        )
        lead_snp_clusters.append(clumped_ind_sigs)
        # Remove processed Ind. sig. SNPs from the pool
        remaining_ind_sig = remaining_ind_sig.filter(~pl.col("uniq_id").is_in(partner_ids))
        print(f"Lead SNP: {lead_id} | Ind. sig. SNPs merged: {clumped_ind_sigs.height}")
    return pl.concat(lead_snp_clusters) if lead_snp_clusters else pl.DataFrame()

# --- Function 4: Boundary Definition ---
def define_genomic_risk_loci(ind_sig_df, leads_df, global_locus_start=1):
    if leads_df.is_empty() or ind_sig_df.is_empty():
        # Return empty DFs and the unchanged counter
        return pl.DataFrame(), pl.DataFrame(), pl.DataFrame(), pl.DataFrame(), global_locus_start
    # --- 1. IndSig Clusters (Aggregating Member Data) ---
    ind_sig_clusters = ind_sig_df.group_by("ind_sig_SNP_id").agg([
        pl.col("chrcol").first().alias("chr"),
        pl.col("poscol").min().alias("start"),
        pl.col("poscol").max().alias("end"),
        pl.col("n_refsnps").first().alias("n_refsnps"),
        pl.col("uniq_id").count().alias("n_IndSig_Members"),
        pl.col("pcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("min_p"),
        pl.col("becol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("beta"),
        pl.col("secol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("se"),
        pl.col("rsIDcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("rsID"),
        pl.col("eacol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("risk_allele"),
        pl.col("neacol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("non_risk_allele"),
        pl.col("eafcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("EAF"),
        pl.col("pcol").filter(pl.col("pcol") < 5e-8).count().alias("n_p_5e_8"),
        pl.col("pcol").filter(pl.col("pcol") < 5e-5).count().alias("n_p_5e_5"),
        pl.col("pcol").filter(pl.col("pcol") < 0.05).count().alias("n_p_0_05"),
        pl.concat_str([
            pl.col("uniq_id"), pl.lit("("),
            pl.col("r2_with_IndSig").round(4).cast(pl.Utf8), pl.lit(")")
        ]).filter(pl.col("uniq_id") != pl.col("ind_sig_SNP_id")).alias("temp_list")
    ]).with_columns([
        pl.concat_str([
            pl.lit("{"), pl.col("ind_sig_SNP_id"), pl.lit(": "),
            pl.col("temp_list").list.join("; "), pl.lit("}")
        ]).alias("IndSig_LD_Group")
    ]).drop("temp_list").sort(["chr", "start"])
    # --- 2. Initial Clusters (Group IndSigs into Lead SNPs) ---
    initial_clusters = ind_sig_clusters.join(
        leads_df.select(["uniq_id", "lead_SNP_id", "r2_with_Lead"]),
        left_on="ind_sig_SNP_id", right_on="uniq_id", how="inner"
    ).group_by("lead_SNP_id").agg([
        pl.col("chr").first(),
        pl.col("start").min(),
        pl.col("end").max(),
        pl.col("min_p").min(),
        pl.col("beta").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_beta"),
        pl.col("se").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_se"),
        pl.col("rsID").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_rsID"),
        pl.col("risk_allele").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_allele"),
        pl.col("non_risk_allele").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_non_risk_allele"),
        pl.col("EAF").filter(pl.col("min_p") == pl.col("min_p").min()).first().alias("lead_EAF"),
        pl.col("ind_sig_SNP_id").alias("ind_sig_list"),
        pl.col("IndSig_LD_Group").alias("IndSig_LD_Groups_list"),
        pl.col("n_refsnps").sum().alias("n_refsnps"),
        pl.col("n_IndSig_Members").sum().alias("total_Markers"),
        pl.col("n_p_5e_8").sum().alias("n_p_5e_8_total"),
        pl.col("n_p_5e_5").sum().alias("n_p_5e_5_total"),
        pl.col("n_p_0_05").sum().alias("n_p_0_05_total"),
        pl.concat_str([
            pl.col("ind_sig_SNP_id"), pl.lit("("),
            pl.col("r2_with_Lead").round(4).cast(pl.Utf8), pl.lit(")")
        ]).filter(pl.col("ind_sig_SNP_id") != pl.col("lead_SNP_id")).alias("temp_list")
    ]).with_columns([
        pl.concat_str([
            pl.lit("{"), pl.col("lead_SNP_id"), pl.lit(": "),
            pl.col("temp_list").list.join("; "), pl.lit("}")
        ]).alias("Lead_LD_Group_str")
    ]).drop("temp_list").sort(["chr", "start"])
    # --- 3. Physical Merging (250kb Window) ---
    merged_loci = []
    locus_idx = global_locus_start 
    if not initial_clusters.is_empty():
        rows = initial_clusters.iter_rows(named=True)
        current = next(rows)
        current["all_ind_sigs"] = set(current["ind_sig_list"])
        current["all_ind_sig_groups"] = set(current["IndSig_LD_Groups_list"])
        current["all_leads"] = {current["lead_SNP_id"]}
        current["all_lead_groups"] = {current["Lead_LD_Group_str"]}
        current["merged_by_dist"] = False
        for next_c in rows:
            if (next_c['chr'] == current['chr'] and next_c['start'] <= current['end'] + 250_000):
                current['start'] = min(current['start'], next_c['start'])
                current['end'] = max(current['end'], next_c['end'])
                current["all_ind_sigs"].update(next_c["ind_sig_list"])
                current["all_ind_sig_groups"].update(next_c["IndSig_LD_Groups_list"])
                current["all_leads"].add(next_c["lead_SNP_id"])
                current["all_lead_groups"].add(next_c["Lead_LD_Group_str"])
                current["n_refsnps"] += next_c["n_refsnps"]
                current["merged_by_dist"] = True
                # Update Statistics
                current["n_p_5e_8_total"] += next_c["n_p_5e_8_total"]
                current["n_p_5e_5_total"] += next_c["n_p_5e_5_total"]
                current["n_p_0_05_total"] += next_c["n_p_0_05_total"]
                current["total_Markers"] += next_c["total_Markers"]
                if next_c['min_p'] < current['min_p']:
                    current.update({
                        'lead_SNP_id': next_c['lead_SNP_id'], 'min_p': next_c['min_p'],
                        'lead_beta': next_c['lead_beta'], 'lead_se': next_c['lead_se'],
                        'lead_rsID': next_c['lead_rsID'], 'lead_allele': next_c['lead_allele'],
                        'lead_non_risk_allele': next_c['lead_non_risk_allele'],
                        'lead_EAF': next_c['lead_EAF']
                    })
            else:
                current["Genomic_locus"] = locus_idx
                merged_loci.append(current)
                locus_idx += 1
                current = next_c
                current["all_ind_sigs"] = set(current["ind_sig_list"])
                current["all_ind_sig_groups"] = set(current["IndSig_LD_Groups_list"])
                current["all_leads"] = {current["lead_SNP_id"]}
                current["all_lead_groups"] = {current["Lead_LD_Group_str"]}
                current["merged_by_dist"] = False
        current["Genomic_locus"] = locus_idx
        merged_loci.append(current)
        locus_idx += 1
    # --- 4. Final Formatting ---
    final_rows = []
    l_map = []
    for m in merged_loci:
        length_kb = (m["end"] - m["start"]) / 1000
        final_rows.append({
            "Genomic_locus": m["Genomic_locus"],
            "Locus_Length_kb": round(length_kb, 2),
            "Merged_by_Distance": 1 if m["merged_by_dist"] else 0,
            "Lead_uniqID": m["lead_SNP_id"], 
            "Lead_rsID": m["lead_rsID"],
            "Lead_Risk_Allele": m["lead_allele"],
            "Lead_NonRisk_Allele": m["lead_non_risk_allele"],
            "Lead_EAF": m["lead_EAF"],
            "chr": m["chr"], "start": m["start"], "end": m["end"],
            "P-value": f"{m['min_p']:.4e}",
            "Beta": round(m["lead_beta"], 4) if m["lead_beta"] is not None else None,
            "SE": round(m["lead_se"], 4) if m["lead_se"] is not None else None,
            "nIndSigSNPs": len(m["all_ind_sigs"]),
            "ind_sig_list": "; ".join(sorted(m["all_ind_sigs"])),
            "nLeadSNPs": len(m["all_leads"]),
            "LeadSNPs": "; ".join(sorted(m["all_leads"])),
            "n_refsnps": m["n_refsnps"],
            "total_Markers": m["total_Markers"],
            "n_p_5e_8": m["n_p_5e_8_total"],
            "n_p_5e_5": m["n_p_5e_5_total"],
            "n_p_0_05": m["n_p_0_05_total"],
            "IndSig_LD_Groups": " | ".join(sorted(m["all_ind_sig_groups"])),
            "Lead_LD_Groups": " | ".join(sorted(m["all_lead_groups"]))
        })
        for lead_id in m["all_leads"]:
            l_map.append({"lead_SNP_id": lead_id})
    # Prepare return DFs
    final_df = pl.DataFrame(final_rows)
    id_map_df = pl.DataFrame(l_map)
    full_hierarchy = ind_sig_df.join(
        leads_df.select(["uniq_id", "lead_SNP_id", "r2_with_Lead"]),
        left_on="ind_sig_SNP_id", right_on="uniq_id", how="inner"
    )
    if not id_map_df.is_empty():
        full_hierarchy = full_hierarchy.join(id_map_df, on="lead_SNP_id", how="left")
    return final_df, full_hierarchy, ind_sig_clusters, initial_clusters, locus_idx


# --- Main Pipeline ---
def run_pipeline(sumstat_path, ld_folder):
    # Load and prepare initial data
    gwas = pl.read_csv(sumstat_path, separator='\t').with_columns(
        pl.concat_str([
            pl.col("chrcol").cast(pl.Utf8), 
            pl.col("poscol").cast(pl.Utf8), 
            pl.col("neacol"), 
            pl.col("eacol")
        ], separator="_").alias("uniq_id")
    )
    # Filter for testing or specific chromosomes
    #gwas = gwas.filter(pl.col("chrcol").is_in([1, 2, 3]))
    res = {"summary": [], "hierarchy": [], "ind_sig": [], "lead_un": []}
    for chrom in gwas["chrcol"].unique().sort().to_list():
        print(f"--- Processing Chromosome {chrom} ---")
        ld_path = os.path.join(ld_folder, f"EUR_chr{chrom}.ld.gz")
        if not os.path.exists(ld_path): 
            print(f"Warning: LD file not found for chromosome {chrom}")
            continue
        chr_df = gwas.filter(pl.col("chrcol") == chrom)
        # Step 1: Clumping for Independent Sig SNPs (r2 >= 0.6)
        ind_sig = find_ind_sig_snps(chr_df, ld_path)
        if ind_sig.is_empty(): continue
        # Step 2: Clumping for Lead SNPs (r2 >= 0.1)
        leads = find_lead_snps(ind_sig, ld_path)
        if leads.is_empty(): continue
        # Step 3: Define Loci and Merge by 250kb
        current_locus_start = 1
        summ, hier, is_c, l_un, current_locus_start = define_genomic_risk_loci(
            ind_sig, leads, global_locus_start=current_locus_start
        )
        if not summ.is_empty():
            # Create a mapping to add the final Genomic Locus index back to the hierarchy
            #l_map = summ.select([pl.col("uniqID").alias("lead_SNP_id"), pl.col("Genomic_locus")])
            l_map = summ.select([pl.col("Lead_uniqID").alias("lead_SNP_id"), pl.col("Genomic_locus")])
            res["summary"].append(summ)
            # Join hierarchy with the locus map so every SNP knows its final Locus ID
            res["hierarchy"].append(hier.join(l_map, on="lead_SNP_id", how="left"))
            res["ind_sig"].append(is_c)
            res["lead_un"].append(l_un)
    # --- SAVE WITH AUTOMATIC LIST FLATTENING ---
    if res["summary"]:
        file_map = {
            "summary": "GenomicRiskLoci_Summary.txt",
            "hierarchy": "GenomicRiskLoci_Hierarchy.txt",
            "ind_sig": "IndSig_Clusters_Boundaries.txt",
            "lead_un": "Lead_Clusters_LD_Only.txt"
        }
        for key, filename in file_map.items():
            df_final = pl.concat(res[key])
            # --- THE FIX: Convert List columns to Strings instead of dropping them ---
            cols_to_write = []
            for col_name in df_final.columns:
                if df_final[col_name].dtype.is_nested():
                    # Turn [SNP1, SNP2] into "SNP1; SNP2"
                    cols_to_write.append(pl.col(col_name).list.join("; "))
                else:
                    cols_to_write.append(pl.col(col_name))
            df_final = df_final.select(cols_to_write)
            df_final.write_csv(filename, separator='\t')
        print("\nAll 4 files saved successfully with all member SNPs preserved.")

        

# --- Execution ---
if __name__ == "__main__":
    sumstat_path = "AllOA_formatted.tsv.gz"
    ld_folder = "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/"
    run_pipeline(sumstat_path, ld_folder)
    


LEAD_P=0.000007

sumstat_path="AllOA_formatted.tsv"
ld_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/"
