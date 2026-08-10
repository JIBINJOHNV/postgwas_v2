import shlex
from pathlib import Path
from typing import Dict, Tuple

from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.processes import run_checked_command


def run_gwas2vcf(
    gwas_outputname: str,
    chromosome: str,
    grch_version: str,
    output_folder: str,
    fasta: str,
    main_script_path: str,
    output_layout: Dict[str, str],
    python_bin: str,
    logger=None,
) -> Tuple[str, int]:
    """Run the bundled adapter without dbSNP and preserve its transcript.

    Record-level dbSNP assignment deliberately happens only after normalization
    in :func:`annotate_and_liftover_vcf`, where bcftools matches exact alleles.
    """

    # -------------------------------------------------------
    # Setup paths
    # -------------------------------------------------------
    outdir = Path(output_folder)
    outdir.mkdir(parents=True, exist_ok=True)
    values = {
        "dataset_id": gwas_outputname,
        "chromosome": chromosome,
        "build": grch_version,
    }
    input_tsv = configured_output_path(outdir, output_layout["adapter_input"], **values)
    output_vcf = configured_output_path(
        outdir, output_layout["adapter_output_vcf"], **values,
    )
    json_dict = configured_output_path(outdir, output_layout["adapter_mapping"], **values)
    transcript = configured_output_path(
        outdir, output_layout["adapter_transcript"], **values,
    )

    # -------------------------------------------------------
    # Input validation
    # -------------------------------------------------------
    for file_path, label in [
        (fasta, "FASTA"),
        (main_script_path, "Main script"),
    ]:
        if not Path(file_path).exists():
            message = "Missing %s file: %s" % (label, file_path)
            raise FileNotFoundError(message)

    # -------------------------------------------------------
    # Build command
    # -------------------------------------------------------
    cmd = [
        resolve_executable(
            python_bin, "Python executable", error_type=RuntimeError,
        ),
        main_script_path,
        "--data",
        str(input_tsv),
        "--ref",
        str(fasta),
        "--out",
        str(output_vcf),
        "--id",
        str(gwas_outputname),
        "--json",
        str(json_dict),
    ]

    command_str = " ".join(shlex.quote(part) for part in cmd)

    run_checked_command(
        cmd,
        "GWAS-to-VCF conversion for chromosome %s" % chromosome,
        logger=logger,
        error_type=RuntimeError,
        stdout_path=transcript,
        stderr_to_stdout=True,
        stdout_header="tool=gwas2vcf\ncommand=%s\n\n" % command_str,
    )
    if logger is not None:
        logger.record("OUTPUT", "gwas2vcf_transcript", path=str(transcript))
    return command_str, 0
