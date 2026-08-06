from pathlib import Path
import subprocess
import textwrap


def vcf_to_standered_ldclump(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert summary-statistics VCF into tab-delimited format with columns:
    """
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{sample_name}_formatted.tsv"
    log_file = output_dir / f"{sample_name}_vcf_conversion.log"
    cmd = textwrap.dedent(f"""
        {{
            printf "chrcol\\tposcol\\tneacol\\teacol\\trsIDcol\\tpcol\\tbecol\\tsecol\\teafcol\\n"
            bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t[%LP]\\t[%ES]\\t[%SE]\\t[%AF]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                # Convert -log10(P) to P
                raw_p = exp(-$6 * log(10))
                $6 = sprintf("%.6g", raw_p)

                # Replace missing ES with NA
                if ($7 == "" || $7 == ".") $7 = "NA"

                # Replace missing SE with NA
                if ($8 == "" || $8 == ".") $8 = "NA"

                print
            }}'
        }} > "{output_file}"
    """)
    try:
        with open(log_file, "w") as lf:
            lf.write("### vcf_to_ldsc command\n\n")
            lf.write(cmd + "\n")
            lf.write("\n### Command output (stdout / stderr)\n\n")
            subprocess.run(
                cmd,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=lf,
                stderr=lf,
            )
    except subprocess.CalledProcessError as e:
        with open(log_file, "a") as lf:
            lf.write("\n❌ Conversion failed\n")
            lf.write(str(e) + "\n")
        raise RuntimeError(
            f"vcf_to_ldsc failed for {vcf_path}. "
            f"See log file: {log_file}"
        )
    return {
        "output_file": str(output_file),
        "log_file": str(log_file),
    }

