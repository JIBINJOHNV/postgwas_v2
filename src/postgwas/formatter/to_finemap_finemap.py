import subprocess
from pathlib import Path
import polars as pl


from pathlib import Path
import subprocess
import textwrap


def vcf_to_finemap_finemap(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert a single summary-statistics VCF file into a FINEMAP-compatible
    tab-delimited text file.

    Parameters
    ----------
    sumstat_vcf : str
        Path to the input VCF (e.g., "PGC3_SCZ_european_chr1.vcf.gz")
    output_folder : str
        Folder where the converted text file will be written
    sample_name : str
        Prefix for output naming (e.g., "PGC3_SCZ_european")

    Output
    ------
    <output_folder>/<sample_name>_finemap.tsv
    """

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output file
    output_file = output_dir / f"{sample_name}_finemap_finemap.tsv"

    # Log file (prefix + function name)
    log_file = output_dir / f"{sample_name}_vcf_to_finemap_finemap.log"

    # Build bcftools → query → sed pipeline
    cmd = textwrap.dedent(f"""
    {{
        printf "rsid\\tchromosome\\tposition\\tallele1\\tallele2\\tmaf\\tbeta\\tse\\tNEF\\n"
        bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
        bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%ALT\\t%REF\\t[%AF]\\t[%ES]\\t[%SE]\\t[%NEF]\\n' | \
        sed 's|:|_|g'
    }} > "{output_file}"
    """)

    # ---------------------------------------------------------
    # Run command + log everything
    # ---------------------------------------------------------
    try:
        with open(log_file, "w") as lf:
            lf.write("### vcf_to_finemap command\n\n")
            lf.write(cmd + "\n")
            lf.write("\n### Command output (stdout / stderr)\n\n")

            subprocess.run(
                cmd,
                shell=True,
                check=True,
                executable="/bin/bash",
                stdout=lf,
                stderr=lf,
            )

    except subprocess.CalledProcessError as e:
        with open(log_file, "a") as lf:
            lf.write("\n❌ Conversion failed\n")
            lf.write(str(e) + "\n")

        raise RuntimeError(
            f"vcf_to_finemap failed for {vcf_path}. "
            f"See log file: {log_file}"
        )

    return {
        "finemap_input": str(output_file),
        "log_file": str(log_file),
    }
