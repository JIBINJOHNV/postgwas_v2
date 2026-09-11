"""Scientific provenance written into merged harmonisation VCF headers."""

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.vcf_provenance import (
    build_harmonisation_vcf_provenance,
)


def _vcf_config():
    return load_configuration().modules.harmonisation.vcf_processing.model_dump()


def _resources(tmp_path):
    resource_root = tmp_path / "resources"
    return {
        chromosome: {
            "default_eaf_file": str(
                resource_root / ("GRCh37_1000G_freq_chr%s.tsv.gz" % chromosome)
            ),
            "default_eaf_reference_column": "EUR",
            "default_comparison_af_file": str(
                resource_root / ("GRCh37_1000G_chr%s.vcf.gz" % chromosome)
            ),
            "user_eaf_file": str(
                resource_root / ("study_eaf_chr%s.tsv.gz" % chromosome)
            ),
            "user_info_file": str(
                resource_root / ("study_info_chr%s.tsv.gz" % chromosome)
            ),
            "dbsnp_file": str(
                resource_root / ("dbsnp_chr%s.vcf.gz" % chromosome)
            ),
            "genome_fasta_file": str(resource_root / "GRCh37.fa"),
            "target_fasta": str(resource_root / "GRCh38.fa"),
            "annot_path": str(resource_root / "genes_GRCh37.gff3.gz"),
            "chain_file": str(resource_root / "GRCh37_to_GRCh38.chain.gz"),
        }
        for chromosome in ("1", "2", "X")
    }


def test_external_sources_and_neglog_pvalue_have_complete_provenance(tmp_path):
    resource_directory = tmp_path / "resources"
    output_directory = tmp_path / "output" / "study" / "harmonisation"
    manifest = output_directory / "study_run_manifest.json"
    policies = default_policies().with_overrides({
        "sample_size.trait_type": "binary",
    })
    values = build_harmonisation_vcf_provenance(
        sample={
            "gwas_outputname": "study",
            "sumstat_file": str(tmp_path / "study.tsv.gz"),
            "output_root": str(tmp_path / "output"),
            "beta_or_col": "OR",
            "se_col": "SE",
            "imp_z_col": None,
            "pval_col": "NEG_LOG10_P",
            "eaf_col": None,
            "eaffile": str(resource_directory / "study_eaf_chr{chromosome}.tsv.gz"),
            "eafcolumn": "EAF",
            "info_source": "external",
            "imp_info_col": None,
            "infofile": str(resource_directory / "study_info_chr{chromosome}.tsv.gz"),
            "infocolumn": "INFO",
            "fixed_info": None,
            "ncase_col": "N_CASE",
            "ncontrol_col": "N_CONTROL",
            "ncase": None,
            "ncontrol": None,
        },
        study_decisions={
            "effect_type": "odds_ratio",
            "effect_type_source": "detector",
            "se_scale": "log_odds",
            "se_scale_source": "dataset_z_pvalue_cross_check",
            "pvalue_type": "neglog10",
            "pvalue_type_source": "detector",
            "info_score_type": "mach_rsq",
            "info_score_type_source": "automatic_full_dataset",
            "frequency_type": "effect_allele_frequency",
            "eaf_is_maf_source": "external_effect_allele_frequency",
            "strand": "forward",
        },
        resource_preflight={
            "require_default_eaf": True,
            "resources_by_chromosome": _resources(tmp_path),
        },
        completed_chromosomes=["1", "2"],
        policies=policies,
        resource_directory=str(resource_directory),
        output_directory=str(output_directory),
        run_manifest=str(manifest),
        vcf_config=_vcf_config(),
    )

    assert values["resource_directory"] == str(resource_directory.resolve())
    assert values["output_directory"] == str((tmp_path / "output").resolve())
    assert values["dataset_output_directory"] == str(output_directory.resolve())
    assert values["run_manifest"] == str(manifest.resolve())
    assert "GRCh37_1000G_freq_chr1.tsv.gz" in values["resolved_resource_files"]
    assert "GRCh37_1000G_freq_chr2.tsv.gz" in values["strand_reference_files"]
    assert values["af_reference_template"] == str(
        (resource_directory / "study_eaf_chr{chromosome}.tsv.gz").resolve()
    )
    assert "study_eaf_chr1.tsv.gz" in values["af_reference_files"]
    assert values["info_reference_template"] == str(
        (resource_directory / "study_info_chr{chromosome}.tsv.gz").resolve()
    )
    assert "study_info_chr2.tsv.gz" in values["info_reference_files"]
    assert values["effect_formula"] == "BETA = ln(OR)"
    assert values["se_formula"] == (
        "No SE-scale conversion; remaining missing |Z| = "
        "-ndtri_exp(ln(P)-ln(2)); SE = abs(BETA/Z)"
    )
    assert "null cells configured for recovery" in values["se_source"]
    assert values["info_interpretation"].startswith("External reference proxy")
    assert "MaCH Rsq scale" in values["info_interpretation"]
    assert values["pvalue_input_scale"] == "-log10(P)"
    assert values["pvalue_harmonisation_formula"] == "P = 10^(-input)"
    assert values["pvalue_vcf_formula"] == "LP = -log10(P)"
    assert values["sample_size_formula"] == "Neff = 4 / (1/Ncase + 1/Ncontrol)"
    assert set(values) | {
        "vcf_created_at", "input_genome_build", "output_genome_build", "liftover",
    } == set(_vcf_config()["provenance"]["headers"])


def test_raw_pvalue_fixed_info_and_raw_or_se_are_described(tmp_path):
    resources = _resources(tmp_path)
    for chromosome_resources in resources.values():
        chromosome_resources["user_eaf_file"] = None
        chromosome_resources["user_info_file"] = None
    values = build_harmonisation_vcf_provenance(
        sample={
            "gwas_outputname": "study",
            "sumstat_file": str(tmp_path / "study.tsv.gz"),
            "output_root": str(tmp_path / "output"),
            "beta_or_col": "OR",
            "se_col": "OR_SE",
            "imp_z_col": "Z",
            "pval_col": "P",
            "eaf_col": "EAF",
            "eaffile": None,
            "eafcolumn": None,
            "info_source": "fixed_cli",
            "imp_info_col": None,
            "infofile": None,
            "infocolumn": None,
            "fixed_info": 0.99,
            "ncase_col": None,
            "ncontrol_col": "N",
            "ncase": None,
            "ncontrol": None,
        },
        study_decisions={
            "effect_type": "odds_ratio",
            "effect_type_source": "sample_sheet",
            "se_scale": "as_given",
            "se_scale_source": "policy",
            "pvalue_type": "raw",
            "pvalue_type_source": "sample_sheet",
            "info_score_type": "standard_info",
            "info_score_type_source": "automatic_full_dataset",
            "frequency_type": "effect_allele_frequency",
            "eaf_is_maf_source": "study_level_statistic",
            "strand": "not_evaluated",
        },
        resource_preflight={
            "require_default_eaf": False,
            "resources_by_chromosome": resources,
        },
        completed_chromosomes=["1", "2"],
        policies=default_policies().with_overrides({
            "sample_size.trait_type": "quantitative",
        }),
        resource_directory=str(tmp_path / "resources"),
        output_directory=str(tmp_path / "output"),
        run_manifest=str(tmp_path / "output" / "manifest.json"),
        vcf_config=_vcf_config(),
    )

    assert values["se_input_scale"] == "raw odds-ratio scale"
    assert values["se_formula"] == (
        "SE_logOR = SE_OR / OR; missing SE = abs(BETA / Z) after signed "
        "BETA/Z agreement; "
        "remaining missing |Z| = -ndtri_exp(ln(P)-ln(2)); SE = abs(BETA/Z)"
    )
    assert values["pvalue_harmonisation_formula"] == "No harmonisation scale conversion"
    assert values["info_source"] == "user-assigned fixed value"
    assert values["info_fixed_value"] == "0.99"
    assert values["strand_reference_column"] == "not_provided"
    assert values["strand_reference_files"] == "not_provided"
    assert values["sample_size_formula"] == "Neff = supplied total sample size"


def test_multi_value_info_aggregation_preserves_source_column_provenance(tmp_path):
    resources = _resources(tmp_path)
    for chromosome_resources in resources.values():
        chromosome_resources["user_eaf_file"] = None
        chromosome_resources["user_info_file"] = None

    values = build_harmonisation_vcf_provenance(
        sample={
            "gwas_outputname": "study",
            "sumstat_file": str(tmp_path / "study.tsv.gz"),
            "output_root": str(tmp_path / "output"),
            "beta_or_col": "BETA",
            "se_col": "SE",
            "imp_z_col": "Z",
            "pval_col": "P",
            "eaf_col": "EAF",
            "eaffile": None,
            "eafcolumn": None,
            "info_source": "internal",
            "imp_info_col": "__postgwas_multi_value_info",
            "info_multi_value_source_column": "INFO",
            "info_multi_value_aggregation": "median",
            "infofile": None,
            "infocolumn": None,
            "fixed_info": None,
            "ncase_col": None,
            "ncontrol_col": "N",
            "ncase": None,
            "ncontrol": None,
        },
        study_decisions={
            "effect_type": "beta",
            "effect_type_source": "sample_sheet",
            "se_scale": None,
            "se_scale_source": "not_applicable_to_beta",
            "pvalue_type": "raw",
            "pvalue_type_source": "sample_sheet",
            "info_score_type": "standard_info",
            "info_score_type_source": "automatic_full_dataset",
            "frequency_type": "effect_allele_frequency",
            "eaf_is_maf_source": "study_level_statistic",
            "strand": "forward",
        },
        resource_preflight={
            "require_default_eaf": False,
            "resources_by_chromosome": resources,
        },
        completed_chromosomes=["1"],
        policies=default_policies().with_overrides({
            "sample_size.trait_type": "quantitative",
        }),
        resource_directory=str(tmp_path / "resources"),
        output_directory=str(tmp_path / "output"),
        run_manifest=str(tmp_path / "output" / "manifest.json"),
        vcf_config=_vcf_config(),
    )

    assert values["info_input_column"] == "INFO"
    assert "unweighted row-wise median" in values["info_interpretation"]
    assert "working column __postgwas_multi_value_info" in (
        values["info_interpretation"]
    )
    assert "sample-size weights were not inferred" in (
        values["info_interpretation"]
    )


def test_z_only_vcf_provenance_matches_configured_reconstruction_method(tmp_path):
    resources = _resources(tmp_path)
    for chromosome_resources in resources.values():
        chromosome_resources["user_eaf_file"] = None
        chromosome_resources["user_info_file"] = None

    def build(policies, chromosomes=("1", "2")):
        return build_harmonisation_vcf_provenance(
            sample={
                "gwas_outputname": "study",
                "sumstat_file": str(tmp_path / "study.tsv.gz"),
                "output_root": str(tmp_path / "output"),
                "beta_or_col": None,
                "se_col": None,
                "imp_z_col": "Z",
                "pval_col": "P",
                "eaf_col": "EAF",
                "eaffile": None,
                "eafcolumn": None,
                "info_source": "fixed_cli",
                "imp_info_col": None,
                "infofile": None,
                "infocolumn": None,
                "fixed_info": 1.0,
                "ncase_col": None,
                "ncontrol_col": "N",
                "ncase": None,
                "ncontrol": None,
            },
            study_decisions={
                "effect_type": None,
                "effect_type_source": "no_effect_column",
                "se_scale": None,
                "se_scale_source": "no_standard_error_column",
                "pvalue_type": "raw",
                "pvalue_type_source": "sample_sheet",
                "frequency_type": "effect_allele_frequency",
                "eaf_is_maf_source": "study_level_statistic",
                "strand": "not_evaluated",
            },
            resource_preflight={
                "require_default_eaf": False,
                "resources_by_chromosome": resources,
            },
            completed_chromosomes=list(chromosomes),
            policies=policies,
            resource_directory=str(tmp_path / "resources"),
            output_directory=str(tmp_path / "output"),
            run_manifest=str(tmp_path / "output" / "manifest.json"),
            vcf_config=_vcf_config(),
        )

    default_values = build(default_policies())
    assert default_values["effect_formula"] == (
        "BETA = Z / sqrt(2*EAF*(1-EAF)*Neff)"
    )
    assert default_values["se_formula"] == (
        "SE = 1 / sqrt(2*EAF*(1-EAF)*Neff)"
    )
    assert "effect_from_z.method=metal_large_n" in (
        default_values["effect_harmonisation"]
    )

    zhu_values = build(default_policies().with_overrides({
        "effect_from_z.method": "zhu_2016",
        "effect_from_z.phenotype_standard_deviation": 2.0,
    }))
    assert zhu_values["effect_formula"] == (
        "BETA = phenotype_SD * Z / sqrt(2*EAF*(1-EAF)*(Neff+Z^2)); "
        "phenotype_SD=2.0"
    )
    assert zhu_values["se_formula"] == (
        "SE = phenotype_SD / sqrt(2*EAF*(1-EAF)*(Neff+Z^2)); "
        "phenotype_SD=2.0"
    )
    assert zhu_values["effect_output_scale"] == (
        "additive effect approximation in configured phenotype units"
    )

    x_override_values = build(
        default_policies().with_overrides({
            "effect_from_z.x_chromosome_z_only_action": (
                "allow_autosomal_assumption"
            ),
        }),
        chromosomes=("1", "X"),
    )
    assert "x_chromosome_z_only_action=allow_autosomal_assumption" in (
        x_override_values["effect_harmonisation"]
    )
    assert "x_chromosome_z_only_action=allow_autosomal_assumption" in (
        x_override_values["se_harmonisation"]
    )
