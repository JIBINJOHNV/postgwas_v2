import shlex
from pathlib import Path
from typing import Dict, Tuple

from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.processes import run_checked_command

from .vcf_processing import get_vcf_variant_count


def run_gwas2vcf(
    gwas_outputname: str,
    chromosome: str,
    grch_version: str,
    output_folder: str,
    fasta: str,
    main_script_path: str,
    output_layout: Dict[str, str],
    python_bin: str,
    bcftools_bin: str,
    expected_variants: int,
    logger=None,
) -> Tuple[str, int]:
    """Run the bundled adapter and prove that every exported row reached VCF.

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
    raw_vcf = configured_output_path(
        outdir, output_layout["chromosome_raw_vcf"], **values,
    )
    # Protocol invariant: this runner does not request the adapter's optional
    # CSI mode, so pysam.tabix_index creates the corresponding TBI index.
    raw_vcf_index = Path(str(raw_vcf) + ".tbi")
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
    if isinstance(expected_variants, bool) or not isinstance(expected_variants, int):
        raise ValueError("expected_variants must be a positive integer")
    expected_count = expected_variants
    if expected_count <= 0:
        raise ValueError("expected_variants must be a positive integer")

    # Remove only the exact, module-owned outputs for this dataset/chromosome.
    # Otherwise a stale VCF from an earlier run could satisfy output validation
    # after the current adapter invocation failed to create a new one.
    for stale_output in (output_vcf, raw_vcf, raw_vcf_index):
        try:
            stale_output.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            raise RuntimeError(
                "Could not remove stale GWAS-to-VCF output before chromosome "
                f"{chromosome}: {stale_output}: {exc}"
            ) from exc

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
        expected_outputs=[raw_vcf, raw_vcf_index],
    )
    if logger is not None:
        logger.record("OUTPUT", "gwas2vcf_transcript", path=str(transcript))
        logger.record("OUTPUT", "gwas2vcf_raw_vcf", path=str(raw_vcf))

    resolved_bcftools = resolve_executable(
        bcftools_bin, "bcftools executable", error_type=RuntimeError,
    )
    observed_count = get_vcf_variant_count(
        str(raw_vcf), resolved_bcftools, logger=logger,
    )
    if observed_count is None:
        raise RuntimeError(
            "Chromosome %s GWAS-to-VCF record accounting failed because the "
            "indexed raw VCF record count could not be read: %s"
            % (chromosome, raw_vcf)
        )
    if observed_count != expected_count:
        raise RuntimeError(
            f"Chromosome {chromosome} GWAS-to-VCF record accounting failed: "
            f"{expected_count:,} harmonised row(s) were exported, but the raw VCF "
            f"contains {observed_count:,} record(s). The adapter may have skipped "
            "variants because of reference-allele, FASTA, parsing, or writing "
            "errors; downstream annotation was not started."
        )
    if logger is not None:
        logger.record(
            "VALIDATION",
            "gwas2vcf_record_accounting",
            chromosome=str(chromosome),
            exported_rows=expected_count,
            raw_vcf_records=observed_count,
            balanced=True,
            path=str(raw_vcf),
        )
    return command_str, 0
