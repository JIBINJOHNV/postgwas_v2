import argparse
import multiprocessing
import sys
import json
import polars as pl
from rich_argparse import RichHelpFormatter 
from postgwas.utils.main import validate_path
import textwrap
from postgwas.clis.common_cli import get_defaultresourse_parser
from postgwas.harmonisation.io import read_config
from postgwas.harmonisation.main import run_harmonisation_pipeline


def get_harmonisation_parser():
    """Returns the ArgumentParser for the PostGWAS harmonisation step, structured with a group."""
    
    # add_help=False is VITAL for inheritance
    parser = argparse.ArgumentParser(
        add_help=False, 
        formatter_class=RichHelpFormatter
    ) 

    # 💥 CRITICAL FIX: DEFINE ARGUMENT GROUP HERE 💥
    harmonisation_group = parser.add_argument_group('HARMONISATION Arguments')
    
    # --- UPDATED --CONFIG HELP MESSAGE ---
    config_help_text =textwrap.dedent(
        '''[bold bright_red]Required[/bold bright_red]: Configuration CSV file with full path. 
        It should contain the following columns:\n
            sumstat_file gwas_outputname chr_col pos_col snp_id_col ea_col oa_col eaf_col beta_or_col se_col\n
            imp_z_col pval_col ncontrol_col ncase_col ncontrol ncase imp_info_col delimiter\n
            infofile infocolumn eaffile eafcolumn liftover chr_pos_col resourse_folder output_folder
            '''
    )

    # All arguments now use the 'harmonisation_group' object
    harmonisation_group.add_argument(
        "--config",
        type=validate_path(must_exist=True, must_be_file=True),
        metavar='', 
        help=config_help_text 
    )

    harmonisation_group.add_argument(
        "--defaults",
        type=validate_path(must_exist=True, must_be_file=True),
        metavar='', 
        help="[bold bright_red]Required[/bold bright_red]: YAML file with full path"
    )
    
    return parser

# --- 2. EXECUTABLE LOGIC (The "Pipeline Step") ---
def run_harmonisation(args):
    harmonised_vcfs_all = {}
    failed_datasets = []
    success_count = 0

    # -------------------------
    # Load config
    # -------------------------
    try:
        cfg_list = read_config(args.config)
        if not cfg_list:
            raise ValueError("read_config returned empty output.")
    except Exception as e:
        raise RuntimeError(f"❌ ERROR: Failed to load user config file '{args.config}'. Reason: {e}")

    # -------------------------
    # Core runner (reuse logic)
    # -------------------------
    def run_one(user_cfg):
        sample_name = user_cfg.get("gwas_outputname", "UNKNOWN")
        input_file = user_cfg.get("sumstat_file", "UNKNOWN")

        indent1 = "\t" * 2
        indent = "\t" * 3

        print("\n" + "="*80)
        print(f"{indent1}[DATASET] {sample_name}")
        print(f"{indent1}[INPUT FILE] {input_file}")
        print(f"{indent1}[DEFAULT CONFIG] {args.defaults}")
        print(f"{indent1}[THREADS] {args.nthreads}")
        print(indent1 + "-"*70)

        print(f"{indent}[USER CONFIG - FULL]")
        config_str = json.dumps(user_cfg, indent=4)
        config_str = "\n".join(indent + line for line in config_str.splitlines())
        print(config_str)

        print(indent + "="*70, flush=True)

        result = run_harmonisation_pipeline(
            sample_column_dict=user_cfg,
            default_cfg=args.defaults,
            nthreads=args.nthreads
        )

        return sample_name, input_file, result

    # -------------------------
    # First pass
    # -------------------------
    retry_queue = []

    for user_cfg in cfg_list:
        sample_name = user_cfg.get("gwas_outputname", "UNKNOWN")
        input_file = user_cfg.get("sumstat_file", "UNKNOWN")
        try:
            sample_name, input_file, result = run_one(user_cfg)

            if result:
                harmonised_vcfs_all[sample_name] = result
                success_count += 1
                print(f"[SUCCESS] Completed: {sample_name}", flush=True)
            else:
                raise ValueError("No output returned")

        except Exception as e:
            print("\n" + "!"*70)
            print(f"[ERROR] Harmonisation failed")
            print(f"[DATASET] {sample_name}")
            print(f"[INPUT FILE] {input_file}")
            print(f"[REASON] {str(e)}")
            print("!"*70 + "\n", flush=True)

            failed_datasets.append({
                "dataset": sample_name,
                "input": input_file,
                "reason": str(e)
            })

            retry_queue.append(user_cfg)

    # -------------------------
    # Retry failed (1 pass)
    # -------------------------
    if retry_queue:
        print("\n🔁 RETRYING FAILED DATASETS...\n")

    for user_cfg in retry_queue:
        sample_name = user_cfg.get("gwas_outputname", "UNKNOWN")
        input_file = user_cfg.get("sumstat_file", "UNKNOWN")

        try:
            _, _, result = run_one(user_cfg)

            if result:
                harmonised_vcfs_all[sample_name] = result
                success_count += 1
                print(f"[RECOVERED] {sample_name} succeeded on retry")

                # remove from failed list
                failed_datasets = [f for f in failed_datasets if f["dataset"] != sample_name]

        except Exception as e:
            print(f"[STILL FAILED] {sample_name}: {e}")

    # -------------------------
    # Save failed datasets
    # -------------------------
    # if failed_datasets:
    #     fail_df = pl.DataFrame(failed_datasets)

    #     fail_file = f"{args.run_name}_failed_harmonisation.tsv"
    #     fail_df.write_csv(fail_file, separator="\t")

    #     print(f"\n📁 Failed dataset list saved → {fail_file}")

    # -------------------------
    # Final Summary
    # -------------------------
    total = len(cfg_list)
    failed_count = len(failed_datasets)

    print("\n" + "="*70)
    print("📊 HARMONISATION SUMMARY")
    print("="*70)
    print(f"Total datasets : {total}")
    print(f"Successful     : {success_count}")
    print(f"Failed         : {failed_count}")

    if failed_datasets:
        print("\n❌ Failed datasets:")
        for f in failed_datasets:
            print(f"   • {f['dataset']} → {f['reason']}")

    print("="*70 + "\n")

    # -------------------------
    # Final Check
    # -------------------------
    if success_count == 0:
        raise RuntimeError("❌ All harmonisation attempts failed.")

    return harmonised_vcfs_all

# --- 3. STANDALONE ENTRY POINT (The "CLI") ---
def main():
    """
    Entry point only used when running: python cli.py
    """
    parser = argparse.ArgumentParser(
        description="Run gwas sumstat harmonisation",
        parents=[get_defaultresourse_parser(), get_harmonisation_parser()],
        formatter_class=RichHelpFormatter
    )

    args = parser.parse_args()

    # -------- FIX: If no arguments → show help instead of running -----
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    try:
        df = run_harmonisation(args)
        print("--- Harmonisation Complete ---")
        print(df)

    except Exception as e:
        print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()