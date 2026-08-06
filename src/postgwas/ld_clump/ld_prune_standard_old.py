import polars as pl
import os
import math
import subprocess
from pathlib import Path
import textwrap
from concurrent.futures import ProcessPoolExecutor, as_completed
from io import StringIO



def vcf_to_standered_ldclump(sumstat_vcf: str, output_folder: str, sample_name: str, bcftools_path="bcftools", n_threads=4):
    """
    Converts VCF to TSV using specified bcftools path with multithreading support 
    and detailed command logging.
    """
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{sample_name}_formatted.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_standered_ldclump_conversion.log"
    # Constructing the bash command
    cmd = textwrap.dedent(f"""
        {{
            printf "chrcol\\tposcol\\tneacol\\teacol\\trsIDcol\\tpcol\\tbecol\\tsecol\\teafcol\\n"
            {bcftools_path} view --threads {n_threads} --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            {bcftools_path} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t[%LP]\\t[%ES]\\t[%SE]\\t[%AF]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                raw_p = exp(-$6 * log(10))
                $6 = sprintf("%.6g", raw_p)
                if ($7 == "" || $7 == ".") $7 = "NA"
                if ($8 == "" || $8 == ".") $8 = "NA"
                print
            }}'
        }} > "{output_file}"
    """)
    try:
        with open(log_file, "w") as lf:
            # Writing metadata and the exact command to the log
            lf.write(f"--- VCF Conversion Log ---\n")
            lf.write(f"Sample: {sample_name}\n")
            lf.write(f"Threads: {n_threads}\n")
            lf.write(f"Command used:\n{cmd}\n")
            lf.write(f"--- Execution Output ---\n")
            lf.flush() # Ensure the command is written before the process starts
            # Execute the command
            subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stdout=lf, stderr=lf)
            lf.write(f"\n[STAMP] Conversion completed successfully.\n")
    except subprocess.CalledProcessError:
        raise RuntimeError(f"VCF conversion failed for {sample_name}. Check log for details: {log_file}")
    return str(output_file)



def get_ld_partners(chrom, pos, rsid, ld_file_path, r2_threshold):
    partners_r2 = {rsid: 1.0}
    start = max(0, pos - 1_000_000)
    end   = pos + 1_000_000
    region = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ["tabix", ld_file_path, region],
            capture_output=True,
            text=True,
            check=True
        )
        if not result.stdout.strip():
            return partners_r2
        # 🔥 Explicit dtypes prevent schema inference crashes
        ld_df = pl.read_csv(
            StringIO(result.stdout),
            separator="\t",
            has_header=False,
            new_columns=[
                "chr_a", "pos_a", "snp_a",
                "chr_b", "pos_b", "snp_b", "r2"
            ],
            dtypes={
                "chr_a": pl.Utf8,
                "pos_a": pl.Int64,
                "snp_a": pl.Utf8,
                "chr_b": pl.Utf8,
                "pos_b": pl.Int64,
                "snp_b": pl.Utf8,
                "r2": pl.Float64,  # 🔥 critical fix
            },
        )
        filtered = ld_df.filter(
            (pl.col("r2") >= r2_threshold) &
            (
                (pl.col("snp_a") == rsid) |
                (pl.col("snp_b") == rsid)
            )
        )
        if not filtered.is_empty():
            partners_df = filtered.select([
                pl.when(pl.col("snp_a") == rsid)
                  .then(pl.col("snp_b"))
                  .otherwise(pl.col("snp_a"))
                  .alias("partner_id"),
                pl.col("r2")
            ])
            partners_r2.update(
                dict(
                    zip(
                        partners_df["partner_id"].to_list(),
                        partners_df["r2"].to_list()
                    )
                )
            )
    except subprocess.CalledProcessError as e:
        print(f"[LD ERROR] Tabix failed for {region}: {e}")
    except Exception as e:
        print(f"[LD ERROR] Unexpected parsing error in region {region}: {e}")
    return partners_r2



def find_ind_sig_snps(chr_df, ld_path, lead_p_threshold, r2_clump_threshold, n_threads=1):
    remaining_leads = chr_df.filter(pl.col("pcol") <= lead_p_threshold).sort("pcol")
    available_pool = chr_df.clone()
    ind_sig_clumps = []
    chrom = chr_df.select(pl.col("chrcol").first()).item()
    print(f"[CHR {chrom}] Independent variant detection started.")
    while not remaining_leads.is_empty():
        top_snp = remaining_leads.row(0, named=True)
        sid = top_snp['uniq_id']
        partners_map = get_ld_partners(top_snp['chrcol'], top_snp['poscol'], sid, ld_path, r2_clump_threshold)
        candidate_ids = list(partners_map.keys())
        members = available_pool.filter(pl.col("uniq_id").is_in(candidate_ids))
        if not members.is_empty():
            r2_df = pl.DataFrame({"uniq_id": list(partners_map.keys()), "r2_with_IndSig": list(partners_map.values())})
            members = members.join(r2_df, on="uniq_id", how="left", coalesce=True).with_columns([
                pl.lit(sid).alias("ind_sig_SNP_id"),
                pl.lit(len(partners_map)).alias("n_refsnps")
            ])
            ind_sig_clumps.append(members)
            assigned_ids = members["uniq_id"].to_list()
            available_pool = available_pool.filter(~pl.col("uniq_id").is_in(assigned_ids))
            remaining_leads = remaining_leads.filter(~pl.col("uniq_id").is_in(assigned_ids))
        else:
            remaining_leads = remaining_leads.filter(pl.col("uniq_id") != sid)
        print(f"[CHR {chrom}] Variants remaining: {remaining_leads.height}")
    print(f"[CHR {chrom}] Independent variant detection completed.")
    return pl.concat(ind_sig_clumps) if ind_sig_clumps else pl.DataFrame()




def find_lead_snps(ind_sig_df, ld_path, r2_lead_threshold, n_threads=1):
    if ind_sig_df.is_empty(): return pl.DataFrame()
    ind_sig_heads = ind_sig_df.filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).sort("pcol")
    lead_clusters = []
    remaining = ind_sig_heads
    while not remaining.is_empty():
        top_row = remaining.row(0, named=True)
        lid = top_row['uniq_id']
        partners_map = get_ld_partners(top_row['chrcol'], top_row['poscol'], lid, ld_path, r2_lead_threshold)
        partner_ids = list(partners_map.keys())
        members = ind_sig_heads.filter(pl.col("uniq_id").is_in(partner_ids))
        r2_df = pl.DataFrame({"uniq_id": partner_ids, "r2_with_Lead": list(partners_map.values())})
        members = members.join(r2_df, on="uniq_id", how="left", coalesce=True).with_columns(pl.lit(lid).alias("lead_SNP_id"))
        lead_clusters.append(members)
        remaining = remaining.filter(~pl.col("uniq_id").is_in(partner_ids))
    return pl.concat(lead_clusters) if lead_clusters else pl.DataFrame()



def define_genomic_risk_loci(ind_sig_df, leads_df, merge_dist, global_locus_start=1):
    if leads_df.is_empty():
        return pl.DataFrame(), pl.DataFrame(), pl.DataFrame(), pl.DataFrame(), global_locus_start
    # ------------------------------------------------------------------
    # 1) Aggregate Independent Significant SNP Clusters
    # ------------------------------------------------------------------
    is_clusters = (
        ind_sig_df
        .group_by("ind_sig_SNP_id")
        .agg([
            pl.col("chrcol").first().alias("chr"),
            pl.col("poscol").min().alias("start"),
            pl.col("poscol").max().alias("end"),
            pl.col("poscol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("pos"),
            pl.col("n_refsnps").first(),
            pl.col("uniq_id").count().alias("n_members"),
            pl.col("pcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("p"),
            pl.col("becol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("beta"),
            pl.col("secol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("se"),
            pl.col("rsIDcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("rsID"),
            pl.col("eacol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("ea"),
            pl.col("neacol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("nea"),
            pl.col("eafcol").filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id")).first().alias("eaf"),
            pl.col("pcol").filter(pl.col("pcol") < 5e-8).count().alias("n_5e_8"),
            pl.col("pcol").filter(pl.col("pcol") < 5e-5).count().alias("n_5e_5"),
            pl.col("pcol").filter(pl.col("pcol") < 0.05).count().alias("n_0_05"),
            pl.concat_str([
                pl.col("uniq_id"),
                pl.lit("("),
                pl.col("r2_with_IndSig").round(3).cast(pl.Utf8),
                pl.lit(")")
            ])
            .filter(pl.col("uniq_id") != pl.col("ind_sig_SNP_id"))
            .alias("ld_list")
        ])
        .with_columns([
            pl.concat_str([
                pl.lit("{"),
                pl.col("ind_sig_SNP_id"),
                pl.lit(": "),
                pl.col("ld_list").list.join("; "),
                pl.lit("}")
            ]).alias("IndSig_Group")
        ])
        .sort(["chr", "start"])
    )
    # ------------------------------------------------------------------
    # 2) Safe mapping: ind_sig_SNP_id -> lead_SNP_id (NO right_on join)
    # ------------------------------------------------------------------
    lead_map = (
        leads_df
        .select([
            pl.col("uniq_id").alias("ind_sig_SNP_id"),
            pl.col("lead_SNP_id"),
            pl.col("r2_with_Lead"),
        ])
        .unique(subset=["ind_sig_SNP_id", "lead_SNP_id"])
    )
    # ------------------------------------------------------------------
    # 3) Initial merging into Lead Clusters (same logic, safe join)
    # ------------------------------------------------------------------
    initial = (
        is_clusters
        .join(lead_map, on="ind_sig_SNP_id", how="left", coalesce=False)
    )
    initial = (
        initial
        .group_by("lead_SNP_id")
        .agg([
            pl.col("chr").first(),
            pl.col("start").min(),
            pl.col("end").max(),
            pl.col("p").min(),
            pl.col("pos").filter(pl.col("p") == pl.col("p").min()).first().alias("l_pos"),
            pl.col("beta").filter(pl.col("p") == pl.col("p").min()).first().alias("l_beta"),
            pl.col("se").filter(pl.col("p") == pl.col("p").min()).first().alias("l_se"),
            pl.col("rsID").filter(pl.col("p") == pl.col("p").min()).first().alias("l_rsid"),
            pl.col("ea").filter(pl.col("p") == pl.col("p").min()).first().alias("l_ea"),
            pl.col("nea").filter(pl.col("p") == pl.col("p").min()).first().alias("l_nea"),
            pl.col("eaf").filter(pl.col("p") == pl.col("p").min()).first().alias("l_eaf"),
            pl.col("ind_sig_SNP_id").alias("is_list"),
            pl.col("IndSig_Group").alias("is_groups"),
            pl.col("n_refsnps").sum().alias("sum_refsnps"),
            pl.col("n_members").sum().alias("sum_members"),
            pl.col("n_5e_8").sum().alias("sum_5e_8"),
            pl.col("n_5e_5").sum().alias("sum_5e_5"),
            pl.col("n_0_05").sum().alias("sum_0_05"),
            pl.concat_str([
                pl.col("ind_sig_SNP_id"),
                pl.lit("("),
                pl.col("r2_with_Lead").round(3).cast(pl.Utf8),
                pl.lit(")")
            ])
            .filter(pl.col("ind_sig_SNP_id") != pl.col("lead_SNP_id"))
            .alias("l_ld_list")
        ])
        .with_columns([
            pl.concat_str([
                pl.lit("{"),
                pl.col("lead_SNP_id"),
                pl.lit(": "),
                pl.col("l_ld_list").list.join("; "),
                pl.lit("}")
            ]).alias("Lead_Group")
        ])
        .sort(["chr", "start"])
    )
    # ------------------------------------------------------------------
    # 4) Physical merging based on MERGE_DIST (your original logic)
    # ------------------------------------------------------------------
    merged = []
    idx = global_locus_start
    if not initial.is_empty():
        rows = initial.iter_rows(named=True)
        curr = next(rows)
        curr.update({
            "all_is": set(curr["is_list"]),
            "all_is_g": set(curr["is_groups"]),
            "all_l": {curr["lead_SNP_id"]},
            "all_l_g": {curr["Lead_Group"]},
            "m_dist": False
        })
        for nxt in rows:
            if nxt["chr"] == curr["chr"] and nxt["start"] <= curr["end"] + merge_dist:
                curr["start"] = min(curr["start"], nxt["start"])
                curr["end"]   = max(curr["end"],   nxt["end"])
                curr["all_is"].update(nxt["is_list"])
                curr["all_is_g"].update(nxt["is_groups"])
                curr["all_l"].add(nxt["lead_SNP_id"])
                curr["all_l_g"].add(nxt["Lead_Group"])
                curr["sum_refsnps"] += nxt["sum_refsnps"]
                curr["sum_members"] += nxt["sum_members"]
                curr["sum_5e_8"]     += nxt["sum_5e_8"]
                curr["sum_5e_5"]     += nxt["sum_5e_5"]
                curr["sum_0_05"]     += nxt["sum_0_05"]
                curr["m_dist"] = True
                # update lead info if more significant lead is found
                if nxt["p"] < curr["p"]:
                    curr.update({
                        "lead_SNP_id": nxt["lead_SNP_id"],
                        "p": nxt["p"],
                        "l_pos": nxt["l_pos"],
                        "l_beta": nxt["l_beta"],
                        "l_se": nxt["l_se"],
                        "l_rsid": nxt["l_rsid"],
                        "l_ea": nxt["l_ea"],
                        "l_nea": nxt["l_nea"],
                        "l_eaf": nxt["l_eaf"],
                    })
            else:
                curr["Genomic_locus"] = idx
                merged.append(curr)
                idx += 1
                curr = nxt
                curr.update({
                    "all_is": set(curr["is_list"]),
                    "all_is_g": set(curr["is_groups"]),
                    "all_l": {curr["lead_SNP_id"]},
                    "all_l_g": {curr["Lead_Group"]},
                    "m_dist": False
                })
        curr["Genomic_locus"] = idx
        merged.append(curr)
        idx += 1
    # ------------------------------------------------------------------
    # 5) Build summary table + locus mapping (your original fields)
    # ------------------------------------------------------------------
    summary_rows = []
    l_map = []
    for m in merged:
        summary_rows.append({
            "Genomic_locus": m["Genomic_locus"],
            "Locus_Length_kb": round((m["end"] - m["start"]) / 1000, 2),
            "Merged_by_Distance": 1 if m["m_dist"] else 0,
            "Lead_uniqID": m["lead_SNP_id"],
            "Lead_rsID": m["l_rsid"],
            "CHR": m["chr"],
            "POS": m["l_pos"],
            "START": m["start"],
            "END": m["end"],
            "ea": m["l_ea"],
            "nea": m["l_nea"],
            "eaf": m["l_eaf"],
            "LP": -math.log10(m["p"]) if m["p"] > 0 else 300.0,
            "P_value": m["p"],
            "Beta": m["l_beta"],
            "SE": m["l_se"],
            "n_refsnps": m["sum_refsnps"],
            "n_members": m["sum_members"],
            "n_5e_8": m["sum_5e_8"],
            "n_5e_5": m["sum_5e_5"],
            "n_0_05": m["sum_0_05"],
            "nIndSigSNPs": len(m["all_is"]),
            "nLeadSNPs": len(m["all_l"]),
            "IndSig_LD_Groups": " | ".join(sorted(m["all_is_g"])),
            "Lead_LD_Groups": " | ".join(sorted(m["all_l_g"])),
        })
        for lid in m["all_l"]:
            l_map.append({"lead_SNP_id": lid, "Genomic_locus": m["Genomic_locus"]})
    summary_df = pl.DataFrame(summary_rows) if summary_rows else pl.DataFrame()
    # ------------------------------------------------------------------
    # 6) Hierarchy join (safe; no uniq_id_right ever created)
    # ------------------------------------------------------------------
    hier = ind_sig_df.join(lead_map, on="ind_sig_SNP_id", how="left", coalesce=False)
    if l_map:
        hier = hier.join(pl.DataFrame(l_map), on="lead_SNP_id", how="left", coalesce=False)
    return summary_df, hier, is_clusters, initial, idx




def process_chromosome(chrom, gwas_subset, ld_folder, pop, lead_p, r2_clump, r2_lead, merge_dist, n_threads=1):
    """
    Worker function to process a single chromosome.
    Now accepts n_threads to pass down to sub-functions.
    """
    ld_path = os.path.join(ld_folder, f"{pop}_chr{chrom}.ld.gz")
    if not os.path.exists(ld_path):
        return None
    # Step A: Find Independent Significant SNPs (Passing n_threads)
    ind_sig = find_ind_sig_snps(gwas_subset, ld_path, lead_p, r2_clump)
    if ind_sig.is_empty(): 
        return None
    # Step B: Find Lead SNPs (Passing n_threads)
    leads = find_lead_snps(ind_sig, ld_path, r2_lead, n_threads)
    if leads.is_empty(): 
        return None
    # Step C: Define boundaries
    summ, hier, is_c, l_un, _ = define_genomic_risk_loci(ind_sig, leads, merge_dist, 1)
    return {"summ": summ, "hier": hier, "is_c": is_c, "l_un": l_un, "chrom": chrom}



from concurrent.futures import ThreadPoolExecutor, as_completed

def ld_clump_standard(vcf_path, ld_folder, sample_id, outdir, pop="EUR", 
                 lead_p=5e-8, r2_clump=0.6, r2_lead=0.1, merge_dist=250000, 
                 bcftools_bin="bcftools", n_threads=5):
    from postgwas.utils.main import auto_detect_workers, safe_thread_count
    # 1. Conversion
    tsv_path = vcf_to_standered_ldclump(vcf_path, outdir, sample_id, bcftools_bin, n_threads=n_threads)
    # 2. Load Data
    gwas = pl.read_csv(tsv_path, separator='\t').with_columns(
        pl.concat_str([pl.col("chrcol").cast(pl.Utf8), pl.col("poscol").cast(pl.Utf8), 
                       pl.col("neacol"), pl.col("eacol")], separator="_").alias("uniq_id")
    )
    # --- BALANCED SCHEDULING ---
    chrom_stats = gwas.group_by("chrcol").count().sort("count", descending=True)
    balanced_chroms = []
    big_to_small = chrom_stats["chrcol"].to_list()
    while big_to_small:
        balanced_chroms.append(big_to_small.pop(0))
        if big_to_small:
            balanced_chroms.append(big_to_small.pop(-1))
    # In Docker/Threads, we can usually afford more workers than processes 
    # as long as the external tools (tabix) don't max out CPU.
    #n_workers = n_threads 
    # Use your functions for memory safety
    gwas_mem_gb = gwas.estimated_size() / (1024**3)
    dynamic_gb_per_thread = max(1.0, gwas_mem_gb * 1.1) 
    n_workers = safe_thread_count(auto_detect_workers(), gb_per_thread=dynamic_gb_per_thread)
# ---------------------------
    res_list = []
    print(f"[*] Docker-safe execution: Using ThreadPoolExecutor with {n_workers} workers.")
    # 3. Concurrent Processing using Threads
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                process_chromosome, chrom, gwas.filter(pl.col("chrcol") == chrom), 
                ld_folder, pop, lead_p, r2_clump, r2_lead, merge_dist# internal tool threads set to 1
            ): chrom for chrom in balanced_chroms
        }
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    res_list.append(result)
            except Exception as e:
                print(f"[!] Thread Error on chromosome: {e}")
    # 4. Final Aggregation
    if res_list:
        # ... (rest of your aggregation logic remains the same)
        def chr_key(x):
            c = str(x['chrom']).lower().replace('chr', '')
            return int(c) if c.isdigit() else 99
        res_list.sort(key=chr_key)
        final_res = {"summ": [], "hier": [], "is_c": [], "l_un": []}
        locus_counter = 1
        for r in res_list:
            unique_count = r["summ"]["Genomic_locus"].unique().len()
            r["summ"] = r["summ"].with_columns(pl.col("Genomic_locus") + (locus_counter - 1))
            r["hier"] = r["hier"].with_columns(pl.col("Genomic_locus") + (locus_counter - 1))
            for key in final_res.keys():
                final_res[key].append(r[key])
            locus_counter += unique_count
        for key, suffix in {"summ": "GenomicRiskLoci_Summary.txt", 
                            "hier": "GenomicRiskLoci_Hierarchy.txt", 
                            "is_c": "IndSig_Clusters_Boundaries.txt", 
                            "l_un": "Lead_Clusters_LD_Only.txt"}.items():
            df = pl.concat(final_res[key])
            # Ensure list types are converted to strings before CSV write
            df = df.select([
                pl.col(c).map_elements(lambda x: "; ".join(map(str, x)) if isinstance(x, list) else x, return_dtype=pl.Utf8) 
                if df[c].dtype.is_nested() else pl.col(c) 
                for c in df.columns
            ])
            df.write_csv(Path(outdir) / f"{sample_id}_{suffix}", separator='\t')
        return {"ldpruned_sig_file": str(Path(outdir) / f"{sample_id}_GenomicRiskLoci_Summary.txt")}
    return None
