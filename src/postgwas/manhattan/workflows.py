#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from datetime import datetime
import subprocess
from rich_argparse import RichHelpFormatter





# =====================================================================
#   RUN MODES
# =====================================================================
import sys
import subprocess
from pathlib import Path
from datetime import datetime

def run_assoc_plot_direct(args):
    """Executes assoc_plot.R in direct mode."""
    script_dir = Path(__file__).resolve().parent
    R_SCRIPT = script_dir / "assoc_plot.R"

    if not R_SCRIPT.exists():
        raise FileNotFoundError(f"❌ assoc_plot.R not found in {script_dir}")

    # Build command
    cmd = ["Rscript", str(R_SCRIPT)]

    def add(flag, value):
        if value is not None:
            cmd.append(f"{flag}={value}")

    # 1. Resolve Phenotype (Auto-detect if missing)
    pheno = args.pheno
    if not pheno:
        try:
            # Query VCF for sample ID if not provided
            pheno = subprocess.check_output(
                ["bcftools", "query", "-l", args.vcf], text=True
            ).strip().splitlines()[0] # Take the first sample if multiple
        except Exception:
            pheno = "unknown_sample"

    # 2. Resolve Output Path (Fixing the Logic Bug)
    # If args.png is set, use it. Otherwise, use args.pdf OR the default generated PDF path.
    output_file = None
    
    if args.png:
        output_file = Path(args.png)
        add("--png", str(output_file))
    else:
        # Calculate default PDF path if args.pdf is None
        pdf_path = args.pdf or f"{args.outdir}/{args.sample_id}_manhattanplots.pdf"
        output_file = Path(pdf_path)
        add("--pdf", str(output_file))

    # 3. Add Flags
    add("--genome", args.genome_version)
    #add("--cytoband", args.cytoband)
    add("--vcf", args.vcf)
    add("--pheno", pheno)
    add("--nauto", args.nauto)
    add("--min-af", args.min_af)
    add("--min-lp", args.min_lp)
    add("--loglog-pval", args.loglog_pval)
    add("--cyto-ratio", args.cyto_ratio)
    add("--max-height", args.max_height)
    add("--spacing", args.spacing)

    add("--width", args.width)
    add("--height", args.height)
    add("--fontsize", args.fontsize)

    if args.allelic_shift:
        cmd.append("--as")
    if args.csq:
        cmd.append("--csq")

    # 4. Setup Logging
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Use the RESOLVED pheno variable, not args.pheno (which might be None)
    pheno_tag = pheno if pheno else "NO_PHENO"
    
    log_file = output_file.parent / f"{output_file.stem}.{pheno_tag}.assocplot_{timestamp}.log"

    print(f"        ℹ️Generating Manhattan plots: {output_file}")

    # 5. Run Execution
    with open(log_file, "w") as log:
        log.write("COMMAND:\n" + " ".join(cmd) + "\n\n")
        log.write("------ Rscript output ------\n\n")
        log.flush() # Ensure header is written before Rscript starts

        try:
            process = subprocess.run(
                cmd,
                stdout=log,
                stderr=log,
                text=True,
                check=True # This raises CalledProcessError if return code != 0
            )
        except subprocess.CalledProcessError:
            print(f"\n❌ assoc_plot failed. See log: {log_file}")
            raise RuntimeError("Rscript execution failed.")

    return {"plot": str(output_file), "log": str(log_file)}