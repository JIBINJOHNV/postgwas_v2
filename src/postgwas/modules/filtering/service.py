
from postgwas.modules.filtering.sumstat_filter import filter_gwas_vcf_bcftools
from pathlib import Path
from postgwas.core.execution.runtime import require_executable
from postgwas.cli.compute import memory_limit, resolve_compute_args


def run_sumstat_filter_direct(args,ctx=None, config=None):
    resolve_compute_args(args)
    # Your logic — YOU said do not enforce checks here
    config = config or getattr(args, "module_config", None)
    configured_output_directory = config.output_directory if config else None
    output_directory = Path(configured_output_directory or args.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    require_executable("bcftools")
    require_executable("tabix")

    # -------------------------------------------------
    # Validation
    # -------------------------------------------------
    vcf = config.inputs.vcf if config and config.inputs.vcf else args.vcf
    dataset_id = (
        config.inputs.dataset_id
        if config and config.inputs.dataset_id
        else args.dataset_id
    )
    if not vcf:
        raise ValueError("run_sumstat_filter_direct: --vcf arguments are required.")

    if not dataset_id:
        raise ValueError("run_sumstat_filter_direct: --dataset-id is required for execution.")

    if not args.output_directory:
        raise ValueError(
            "run_sumstat_filter_direct: --output-directory is required for execution."
        )

    outputs=filter_gwas_vcf_bcftools(
        vcf_path=vcf,
        output_folder=str(output_directory),
        output_prefix=dataset_id,
        pval_cutoff=(config.minimum_neglog10_p if config else args.minimum_neglog10_p),
        maf_cutoff=(config.maf_min if config else args.minimum_maf),
        allelefreq_diff_cutoff=(config.frequency_difference_max if config else args.maximum_af_difference),
        info_cutoff=(config.info_min if config else args.minimum_info),
        info_missing=(config.missing_info_action if config else args.missing_info_action),
        info_min=(config.info_min if config else args.minimum_info),
        info_max=(config.info_max if config else args.maximum_info),
        external_af_name=(config.reference_population_tag if config else args.reference_af_column),
        include_indels=(config.include_indels if config else args.include_indels),
        exclude_palindromic=(config.remove_palindromic if config else args.remove_palindromic),
        palindromic_af_lower=(config.palindromic_lower if config else args.palindromic_af_lower),
        palindromic_af_upper=(config.palindromic_upper if config else args.palindromic_af_upper),
        remove_mhc=(config.remove_mhc if config else args.remove_mhc),
        mhc_chrom=(config.mhc.chromosome if config else args.mhc_chrom),
        mhc_start=(config.mhc.start if config else args.mhc_start),
        mhc_end=(config.mhc.end if config else args.mhc_end),
        threads=args.threads,
        max_mem=memory_limit(args.memory_gb),

    )
    # -------------------------------------------------
    # Pipeline mode: register outputs
    # -------------------------------------------------
    if ctx is not None:
        ctx["sumstat_filter"] = outputs
    return outputs
