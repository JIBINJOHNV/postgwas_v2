#!/usr/bin/env python3
import argparse, sys, subprocess
from pathlib import Path
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(
        prog="assoc-plot",
        description="Python wrapper to run GWAS association plots using assoc_plot.R (same folder).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # -------------------------------
    # Inputs
    # -------------------------------
    parser.add_argument("--genome")
    parser.add_argument("--cytoband")
    parser.add_argument("--vcf")
    parser.add_argument("--tbx")
    parser.add_argument("--pheno")

    # -------------------------------
    # Logical Flags
    # -------------------------------
    parser.add_argument("--as", dest="allelic_shift", action="store_true")
    parser.add_argument("--csq", action="store_true")

    # -------------------------------
    # Numeric options
    # -------------------------------
    parser.add_argument("--nauto", type=int, default=22)
    parser.add_argument("--min-af", type=float, default=0.0)
    parser.add_argument("--min-lp", type=int, default=2)
    parser.add_argument("--loglog-pval", type=int, default=10)
    parser.add_argument("--cyto-ratio", type=int, default=25)
    parser.add_argument("--max-height", type=int)
    parser.add_argument("--spacing", type=int, default=10)

    # -------------------------------
    # Outputs
    # -------------------------------
    parser.add_argument("--pdf")
    parser.add_argument("--png")

    # -------------------------------
    # Plot sizing
    # -------------------------------
    parser.add_argument("--width", type=float, default=7.0)
    parser.add_argument("--height", type=float)
    parser.add_argument("--fontsize", type=int, default=12)

    # No args → help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # -------------------------------
    # Require --genome or --cytoband
    # -------------------------------
    if args.genome is None and args.cytoband is None:
        parser.error("Either --genome or --cytoband is required.")

    # -------------------------------
    # Locate assoc_plot.R
    # -------------------------------
    script_dir = Path(__file__).resolve().parent
    R_SCRIPT = script_dir / "assoc_plot.R"

    if not R_SCRIPT.exists():
        sys.exit(f"❌ ERROR: assoc_plot.R not found in {script_dir}")

    # -------------------------------
    # Build command
    # -------------------------------
    cmd = ["Rscript", str(R_SCRIPT)]

    def add_option(flag, value):
        if value is not None:
            cmd.append(f"{flag}={value}")

    add_option("--genome", args.genome)
    add_option("--cytoband", args.cytoband)
    add_option("--vcf", args.vcf)
    add_option("--tbx", args.tbx)
    add_option("--pheno", args.pheno)
    add_option("--nauto", args.nauto)
    add_option("--min-af", args.min_af)
    add_option("--min-lp", args.min_lp)
    add_option("--loglog-pval", args.loglog_pval)
    add_option("--cyto-ratio", args.cyto_ratio)
    add_option("--max-height", args.max_height)
    add_option("--spacing", args.spacing)
    add_option("--pdf", args.pdf)
    add_option("--png", args.png)
    add_option("--width", args.width)
    add_option("--height", args.height)
    add_option("--fontsize", args.fontsize)

    if args.allelic_shift:
        cmd.append("--as")
    if args.csq:
        cmd.append("--csq")

    # ----------------------------------------------------
    # LOG PATH using user --pdf or --png + phenotype
    # ----------------------------------------------------
    if args.pdf:
        base = Path(args.pdf)
    elif args.png:
        base = Path(args.png)
    else:
        parser.error("Either --pdf or --png is required.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pheno_tag = args.pheno if args.pheno else "NO_PHENO"

    # Build final log file path
    log_file = base.parent / f"{base.stem}.{pheno_tag}.assocplot_{timestamp}.log"

    # ----------------------------------------------------
    # Execute Rscript and capture logs
    # ----------------------------------------------------
    with open(log_file, "w") as log:
        log.write("Executed command:\n")
        log.write(" ".join(cmd) + "\n\n")
        log.write("---- Rscript output ----\n\n")

        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=log,
            text=True
        )

    # ----------------------------------------------------
    # Report results
    # ----------------------------------------------------
    if process.returncode != 0:
        print(f"❌ assoc_plot.R failed. See log:\n   {log_file}")
        sys.exit(process.returncode)

    print(f"✅ Plot generated successfully.\n📄 Log saved at:\n   {log_file}")


if __name__ == "__main__":
    main()
