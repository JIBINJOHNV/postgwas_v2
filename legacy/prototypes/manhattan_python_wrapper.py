from pathlib import Path
import subprocess
import shlex
import sys

def run_assoc_plot_r(
    genome=None, cytoband=None, nauto=22, vcf=None, pheno=None,
    as_flag=False, csq=False, tbx=None, region=None,
    min_af=0.0, min_lp=None, loglog_pval=10, cyto_ratio=25,
    max_height=None, spacing=20, pdf=None, png=None,
    width=7.0, height=None, fontsize=12, rscript="assoc_plot.R"
):

    # validation
    if not genome and not cytoband:
        raise ValueError("Either --genome or --cytoband is required")

    if genome and cytoband:
        raise ValueError("Cannot use both --genome and --cytoband")

    if not vcf and not tbx:
        raise ValueError("Either --vcf or --tbx is required")

    if pdf and png:
        raise ValueError("Cannot provide both --pdf and --png")

    if not pdf and not png:
        raise ValueError("Must provide either --pdf or --png output")

    # log file prefix
    prefix = Path(pdf if pdf else png).with_suffix("")
    log_file = Path(str(prefix) + ".log")

    # Build R command
    cmd = ["Rscript", rscript]

    if genome:   cmd += ["--genome", genome]
    if cytoband: cmd += ["--cytoband", cytoband]
    cmd += ["--nauto", str(nauto)]

    if vcf: cmd += ["--vcf", vcf]
    if tbx: cmd += ["--tbx", tbx]
    if pheno: cmd += ["--pheno", pheno]
    if as_flag: cmd += ["--as"]
    if csq: cmd += ["--csq"]
    if region: cmd += ["--region", region]

    cmd += ["--min-af", str(min_af)]
    if min_lp is not None:
        cmd += ["--min-lp", str(min_lp)]

    cmd += [
        "--loglog-pval", str(loglog_pval),
        "--cyto-ratio", str(cyto_ratio),
        "--spacing", str(spacing)
    ]

    if max_height is not None:
        cmd += ["--max-height", str(max_height)]

    if pdf: cmd += ["--pdf", pdf]
    if png: cmd += ["--png", png]

    cmd += ["--width", str(width)]
    if height is not None:
        cmd += ["--height", str(height)]

    cmd += ["--fontsize", str(fontsize)]

    final_cmd = " ".join(shlex.quote(x) for x in cmd)

    # run silently to log file
    with open(log_file, "w") as lf:
        process = subprocess.Popen(final_cmd, shell=True, stdout=lf, stderr=lf)
        process.wait()

    if process.returncode != 0:
        raise RuntimeError(f"assoc_plot.R failed. See log: {log_file}")

    return log_file
