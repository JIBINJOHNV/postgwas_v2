import os
import glob
import pandas as pd
import re
import shutil


def parse_cred_header(file_path):
    """
    Reads the first line to find the Posterior Probability.
    Example: # Post-Pr(# of causal SNPs is 2) = 0.241969
    """
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline()
        match = re.search(r'=\s*([0-9\.]+)', first_line)
        if match:
            return float(match.group(1))
        return 0.0
    except Exception:
        return 0.0

def reformat_snp_id(snp_id):
    """
    Transforms SNP ID from 'Chr_Pos_Ref_Alt' to 'Chr:Pos:Ref_Alt'.
    Replaces only the first 2 underscores with colons.
    """
    if isinstance(snp_id, str):
        # Turns "5_57744788_T_C" into "5:57744788:T_C"
        return snp_id.replace('_', ':', 2)
    return snp_id


def select_and_copy_best_models(input_dir, output_dir):
    """
    Iterates through loci folders, finds the .cred file with the highest
    posterior probability, and copies it to the output_dir.
    Returns a list of dictionaries for the summary CSV.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    all_files_tracking = []
    # Get loci folders
    loci_folders = [f for f in glob.glob(os.path.join(input_dir, "*")) if os.path.isdir(f)]
    print(f"Found {len(loci_folders)} loci folders.")
    for folder in loci_folders:
        locus_name = os.path.basename(folder)
        cred_files = glob.glob(os.path.join(folder, "*.cred*"))
        if not cred_files:
            continue
        # --- Sub-step: Find winner for this locus ---
        temp_locus_data = []
        max_prob = -1.0
        best_file_path = None
        for c_file in cred_files:
            prob = parse_cred_header(c_file)
            fname = os.path.basename(c_file)
            if prob > max_prob:
                max_prob = prob
                best_file_path = c_file
            # Extract K
            try:
                k_val = fname.split('.cred')[-1]
            except:
                k_val = "unknown"
            temp_locus_data.append({
                "Locus": locus_name,
                "File_Name": fname,
                "K_Causal": k_val,
                "Posterior_Prob": prob,
                "Full_Path": c_file
            })
        # --- Sub-step: Process and Copy ---
        for record in temp_locus_data:
            path_ref = record.pop("Full_Path")
            is_winner = (path_ref == best_file_path)
            record["Selected"] = is_winner
            all_files_tracking.append(record)
            if is_winner and path_ref:
                try:
                    dest_path = os.path.join(output_dir, record["File_Name"])
                    shutil.copy2(path_ref, dest_path)
                except Exception as e:
                    print(f"Error copying {path_ref}: {e}")
    return all_files_tracking

def save_summary_csv(tracking_data, output_dir):
    """
    Saves the tracking data to a CSV file in the output directory.
    """
    if tracking_data:
        summary_df = pd.DataFrame(tracking_data)
        summary_df = summary_df.sort_values(by=["Locus", "Posterior_Prob"], ascending=[True, False])
        
        summary_path = os.path.join(output_dir, "finemap_all_models_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"\t\t\t Summary CSV saved to: {summary_path}")
    else:
        print("No data to summarize.")


def filter_metadata(lines, k_val):
    """
    Parses header lines to extract only the stats relevant to the specific 
    credible set (k).
    
    Logic: The lines usually look like: Label V1 NA V2 NA V3 NA
    Indices (split): 0 (Label) 1 2 3 4 5 6
    
    For k=1: We want indices 1, 2
    For k=2: We want indices 3, 4
    
    Formula for value index: 1 + (k-1)*2
    """
    filtered_lines = []
    # Convert 1-based k string to 0-based integer
    target_idx = int(k_val) - 1 
    
    for line in lines:
        # 1. Keep the global Probability line as is (it applies to the whole model)
        if "Post-Pr" in line:
            filtered_lines.append(line)
            continue
            
        parts = line.split()
        
        # 2. Calculate indices for the specific k
        val_idx = 1 + (target_idx * 2)
        na_idx = val_idx + 1
        
        # 3. Check if line is long enough to have data for this k
        if len(parts) > na_idx:
            # Construct new line: Label + Value + Filler(NA)
            # Example: #log10bf 11.503 NA
            new_line = f"{parts[0]} {parts[val_idx]} {parts[na_idx]}\n"
            filtered_lines.append(new_line)
        else:
            # If the line doesn't match the expected "pairs" format, just keep it safe
            filtered_lines.append(line)
            
    return filtered_lines

def create_flames_file(input_dir, output_dir):
    """
    Reads .cred files, splits multiple credible sets (L1, L2...), 
    applies your reformat_snp_id logic, filters header metadata correctly,
    and saves standardized outputs.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cred_files = [f for f in glob.glob(os.path.join(input_dir, "*.cred*")) if os.path.isfile(f)]
    
    if not cred_files:
        print(f"No .cred files found in {input_dir}")
        return

    for file_path in cred_files:
        filename = os.path.basename(file_path)
        try:
            # --- B. Extract Metadata ---
            metadata_lines = []
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        metadata_lines.append(line)
                    else:
                        break
            
            # --- C. Read Data ---
            df = pd.read_csv(file_path, comment='#', delim_whitespace=True)

            # --- D. Find and Loop through all 'cred' columns ---
            cred_cols = [c for c in df.columns if c.startswith("cred")]
            
            if not cred_cols:
                continue

            for cred_col in cred_cols:
                # Extract k (e.g. '1', '2', '3')
                k_val = cred_col.replace("cred", "") 
                prob_col = f"prob{k_val}"

                if prob_col not in df.columns:
                    print(f"Warning: {cred_col} exists but {prob_col} missing in {filename}")
                    continue

                subset_df = df[["index", cred_col, prob_col]].copy()
                
                # Drop rows where this specific credible set has no data (NAs)
                subset_df = subset_df.dropna()

                if subset_df.empty:
                    continue

                # 1. Apply YOUR reformat logic (Chr:Pos:Ref_Alt)
                subset_df[cred_col] = subset_df[cred_col].apply(reformat_snp_id)
                
                # 2. Standardize column names for downstream
                subset_df.columns = ['index', 'cred1', 'prob1']

                # 3. Filter Metadata for this specific k
                specific_header = filter_metadata(metadata_lines, k_val)

                # 4. Write Output
                out_name = f"{filename}_L{k_val}.txt"
                out_path = os.path.join(output_dir, out_name)

                with open(out_path, 'w') as f_out:
                    f_out.writelines(specific_header)
                    subset_df.to_csv(f_out, sep=' ', index=False)
                    
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    print("Formatting complete.")



def process_finemap_output(raw_dir, inter_dir, final_dir):
    """
    Master workflow function.
    1. Selects best models -> Copies to intermediate dir.
    2. Saves summary CSV.
    3. Reformats copied files -> Saves to final directory.
    """
    print("\t\t\t Selecting Best Models & Copying")
    tracking_data = select_and_copy_best_models(raw_dir, inter_dir)
    
    print("\t\t\t Saving Summary Statistics")
    save_summary_csv(tracking_data, inter_dir)
    
    print("\t\t\t Creating Formatted FLAMES Input Files")
    create_flames_file(inter_dir, final_dir)
    print("\n\t\t\t Workflow Finished Successfully.")





# # ================= CONFIGURATION =================
# # 1. Where your original subfolders are (loci_data)
# RAW_DATA_DIR = "loci_data"

# # 2. Where to save the COPY of the best raw .cred files
# INTERMEDIATE_DIR = "finemap_best_results"

# # 3. Where to save the final reformatted/renamed files for FLAMES
# FINAL_OUTPUT_DIR = "flames_input"
# # =================================================


# if __name__ == "__main__":
#     # Execute the full workflow
#     run_flames_workflow(RAW_DATA_DIR, INTERMEDIATE_DIR, FINAL_OUTPUT_DIR)

