import argparse
import sys
import yaml
import multiprocessing

from postgwas.harmonisation.main import run_harmonisation_pipeline
from postgwas.harmonisation.io import read_sumstats, read_config


def main():
    # ---------------------------------------------------------
    # CPU auto-detect: default threads = (cores - 2), minimum = 1
    # ---------------------------------------------------------
    total_cores = multiprocessing.cpu_count()
    auto_threads = max(1, total_cores - 2)

    parser = argparse.ArgumentParser(
        description="Run PostGWAS harmonisation + gwas2vcf pipeline"
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to the user configuration CSV"
    )

    parser.add_argument(
        "--defaults",
        required=True,
        help="Path to the default harmonisation YAML file"
    )

    parser.add_argument(
        "--nthreads",
        type=int,
        default=auto_threads,
        help=f"Number of threads to use (default: auto-detect = {auto_threads})"
    )

    args = parser.parse_args()
    print(f"🧵 Using {args.nthreads} threads (auto-detected {total_cores} cores).")

    # -----------------------------------------------------
    # Load USER config
    # -----------------------------------------------------
    try:
        cfg_list = read_config(args.config)
        if not cfg_list or len(cfg_list) == 0:
            raise ValueError("read_config returned empty output.")
        user_cfg = cfg_list[0]
    except Exception as e:
        print(f"❌ ERROR: Failed to load user config file '{args.config}'.")
        print(f"   Reason: {e}")
        sys.exit(1)


    # -----------------------------------------------------
    # Run pipeline
    # -----------------------------------------------------
    postgwas_qc_df=run_harmonisation_pipeline(
        sample_column_dict=user_cfg,
        default_cfg=args.defaults,
        nthreads=args.nthreads  # <<< pass to pipeline
    )
    
    print(postgwas_qc_df)

if __name__ == "__main__":
    main()
