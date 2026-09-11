"""Build dataset-level scientific provenance for merged GWAS-VCF headers.

This module only translates decisions and paths that the harmonisation run has
already resolved.  It never reads a variant table or resource file, so adding
the metadata does not add a dataset scan.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping, Sequence

from postgwas.core.values import optional_text
from postgwas.modules.harmonisation.effect_from_z import (
    describe_effect_from_z_reconstruction,
)


def _package_version() -> str:
    try:
        return version("postgwas")
    except PackageNotFoundError:
        return "unknown"


def _decision_source(value: Any) -> str:
    return {
        "detector": "automatically detected",
        "sample_sheet": "declared in sample sheet",
        "policy": "explicit run configuration",
        "dataset_z_pvalue_cross_check": "dataset-level Z/P-value cross-check",
        "not_applicable_to_beta": "not applicable to beta effects",
        "no_effect_column": "effect derived from Z",
        "no_standard_error_column": "standard error not supplied",
        "no_study_frequency_column": "external frequency source",
        "study_level_statistic": "initial study-level frequency distribution",
        "external_effect_allele_frequency": (
            "external source declared as effect allele frequency"
        ),
        "chromosome_reference_consensus": (
            "consistent chromosome-level reference comparisons"
        ),
        "chromosome_reference_unresolved": (
            "chromosome-level reference comparison was unresolved"
        ),
        "automatic_full_dataset": "automatically resolved from the complete source",
        "explicit_configuration": "explicit run configuration",
        "fixed_info": "fixed INFO value",
        "frequency_decision_missing": "frequency decision not available",
    }.get(str(value), str(value) if value is not None else "not available")


def _resolved_path(value: Any, label: str) -> str:
    text = optional_text(value)
    if text is None:
        raise ValueError("VCF provenance requires %s." % label)
    return str(Path(text).expanduser().resolve())


def _ordered_resource_paths(
    resources_by_chromosome: Mapping[str, Mapping[str, Any]],
    chromosomes: Sequence[str],
    key: str,
) -> list[str]:
    paths: list[str] = []
    seen: set[str] = set()
    for chromosome in chromosomes:
        value = optional_text(
            (resources_by_chromosome.get(str(chromosome)) or {}).get(key)
        )
        if value is not None:
            resolved = _resolved_path(value, "the resolved %s path" % key)
            if resolved not in seen:
                paths.append(resolved)
                seen.add(resolved)
    return paths


def _joined_paths(paths: Sequence[str], separator: str, missing: str) -> str:
    return separator.join(str(path) for path in paths) if paths else missing


def _effect_provenance(
    sample: Mapping[str, Any],
    decisions: Mapping[str, Any],
    reconstruction: Mapping[str, Any],
    missing: str,
) -> dict[str, str]:
    effect_column = optional_text(sample.get("beta_or_col"))
    se_column = optional_text(sample.get("se_col"))
    z_column = optional_text(sample.get("imp_z_col"))
    effect_type = optional_text(decisions.get("effect_type"))

    if effect_column is not None and effect_type not in ("beta", "odds_ratio"):
        raise ValueError(
            "VCF provenance requires a validated beta/odds-ratio decision when "
            "an effect column is supplied; got %r." % effect_type
        )

    if effect_column is not None and effect_type == "odds_ratio":
        effect_source = "study column"
        effect_harmonisation = "Converted supplied odds ratios to log-odds BETA"
        effect_formula = "BETA = ln(OR)"
        output_scale = "log-odds BETA"
    elif effect_column is not None:
        effect_source = "study column"
        effect_harmonisation = (
            "Supplied BETA retained on its declared scale and oriented to VCF ALT"
        )
        effect_formula = "No effect-scale conversion"
        output_scale = "BETA on the supplied study scale"
    elif z_column is not None and se_column is not None:
        effect_source = "derived from supplied Z and SE"
        effect_harmonisation = "Derived BETA from supplied Z and SE"
        effect_formula = "BETA = Z * SE"
        output_scale = "BETA on the supplied SE scale"
    else:
        effect_source = "derived from Z, EAF and effective sample size"
        effect_harmonisation = "Derived BETA with effect_from_z.method=%s" % (
            reconstruction["method"]
        )
        effect_formula = str(reconstruction["beta_formula"])
        output_scale = str(reconstruction["output_scale"])

    return {
        "effect_input_column": effect_column or missing,
        "effect_input_type": effect_type or missing,
        "effect_type_decision": _decision_source(
            decisions.get("effect_type_source")
        ),
        "effect_source": effect_source,
        "effect_harmonisation": effect_harmonisation,
        "effect_formula": effect_formula,
        "effect_orientation": "Relative to the final VCF ALT allele",
        "effect_output_scale": output_scale,
    }


def _se_provenance(
    sample: Mapping[str, Any],
    decisions: Mapping[str, Any],
    *,
    reconstruction: Mapping[str, Any],
    pvalue_tail: int,
    derive_partial_missing_se: bool,
    missing: str,
) -> dict[str, str]:
    effect_column = optional_text(sample.get("beta_or_col"))
    se_column = optional_text(sample.get("se_col"))
    z_column = optional_text(sample.get("imp_z_col"))
    effect_type = optional_text(decisions.get("effect_type"))
    se_scale = optional_text(decisions.get("se_scale"))

    if se_column is not None:
        source = "study column"
        if effect_type == "odds_ratio" and se_scale not in ("log_odds", "as_given"):
            raise ValueError(
                "VCF provenance requires the validated OR standard-error scale; "
                "expected 'log_odds' or 'as_given', got %r." % se_scale
            )
        if effect_type == "odds_ratio" and se_scale == "as_given":
            input_scale = "raw odds-ratio scale"
            harmonisation = "Converted raw odds-ratio SE to log-odds SE"
            formula = "SE_logOR = SE_OR / OR"
        else:
            input_scale = (
                "log-odds" if effect_type == "odds_ratio" else "BETA scale"
            )
            harmonisation = "Supplied SE retained without scale conversion"
            formula = "No SE-scale conversion"
        if z_column is not None:
            source += "; null cells configured for recovery from harmonised BETA and Z"
            harmonisation += (
                "; finite supplied values retained and null cells recovered where "
                "BETA and signed Z were usable and directionally consistent"
            )
            formula += (
                "; missing SE = abs(BETA / Z) after signed BETA/Z agreement"
            )
        if derive_partial_missing_se:
            source += (
                "; remaining null cells configured for recovery from raw P-value"
                if z_column is not None
                else "; null cells configured for recovery from harmonised BETA "
                "and raw P-value"
            )
            harmonisation += (
                "; remaining null cells recovered where BETA and P were usable"
            )
            formula += (
                "; remaining missing |Z| = -ndtri_exp(ln(P)-ln(%d)); "
                "SE = abs(BETA/Z)"
            ) % pvalue_tail
    elif effect_column is not None and z_column is not None:
        source = "derived from harmonised BETA and supplied Z"
        input_scale = missing
        harmonisation = (
            "Derived positive SE from directionally consistent BETA and signed Z"
        )
        formula = "SE = abs(BETA / Z) after signed BETA/Z agreement"
    elif effect_column is not None:
        source = "derived from harmonised BETA and raw P-value"
        input_scale = missing
        harmonisation = (
            "Derived SE from BETA and the configured %d-sided normal P-value"
            % pvalue_tail
        )
        formula = "|Z| = -ndtri_exp(ln(P)-ln(tail)); SE = abs(BETA/Z)"
    else:
        source = "derived from Z, EAF and effective sample size"
        input_scale = missing
        harmonisation = "Derived SE with effect_from_z.method=%s" % (
            reconstruction["method"]
        )
        formula = str(reconstruction["se_formula"])

    return {
        "se_input_column": se_column or missing,
        "se_source": source,
        "se_input_scale": input_scale,
        "se_scale_decision": _decision_source(decisions.get("se_scale_source")),
        "se_harmonisation": harmonisation,
        "se_formula": formula,
    }


def _sample_size_provenance(
    sample: Mapping[str, Any], trait_type: str, missing: str,
) -> dict[str, str]:
    case_column = optional_text(sample.get("ncase_col"))
    control_column = optional_text(sample.get("ncontrol_col"))
    case_value = optional_text(sample.get("ncase"))
    control_value = optional_text(sample.get("ncontrol"))
    has_cases = case_column is not None or case_value is not None
    has_controls = control_column is not None or control_value is not None

    if trait_type in ("auto", "binary") and has_cases and has_controls:
        source = "case and control counts supplied by column or fixed value"
        formula = "Neff = 4 / (1/Ncase + 1/Ncontrol)"
        trait_label = (
            "case_control (inferred from case and control inputs)"
            if trait_type == "auto" else "case_control"
        )
    elif has_controls:
        source = "single total sample-size input"
        formula = "Neff = supplied total sample size"
        trait_label = trait_type
    else:
        source = "not available"
        formula = "Neff not computed"
        trait_label = trait_type

    return {
        "trait_type": trait_label,
        "case_count_column": case_column or missing,
        "case_count_value": case_value or missing,
        "control_count_column": control_column or missing,
        "control_count_value": control_value or missing,
        "sample_size_source": source,
        "sample_size_formula": formula,
    }


def build_harmonisation_vcf_provenance(
    *,
    sample: Mapping[str, Any],
    study_decisions: Mapping[str, Any],
    resource_preflight: Mapping[str, Any],
    completed_chromosomes: Sequence[str],
    policies: Any,
    resource_directory: str,
    output_directory: str,
    run_manifest: str,
    vcf_config: Mapping[str, Any],
) -> dict[str, str]:
    """Return configured-header logical values for one harmonised dataset."""
    provenance_config = dict(vcf_config["provenance"])
    missing = str(provenance_config["missing_value"])
    separator = str(provenance_config["path_separator"])
    output_fields = dict(provenance_config["output_fields"])
    resources = dict(resource_preflight.get("resources_by_chromosome") or {})
    chromosomes = [str(chromosome) for chromosome in completed_chromosomes]
    if not chromosomes:
        raise ValueError(
            "VCF provenance cannot be built without a completed chromosome."
        )
    absent_maps = [
        chromosome for chromosome in chromosomes if chromosome not in resources
    ]
    if absent_maps:
        raise ValueError(
            "VCF provenance is missing the preflighted resource map for "
            "chromosome(s): %s." % ", ".join(absent_maps)
        )

    def paths(key: str) -> str:
        return _joined_paths(
            _ordered_resource_paths(resources, chromosomes, key),
            separator,
            missing,
        )

    all_paths: list[str] = []
    seen_paths: set[str] = set()
    resource_keys = [
        "default_comparison_af_file", "user_eaf_file", "user_info_file",
        "dbsnp_file", "genome_fasta_file", "target_fasta", "annot_path",
        "chain_file",
    ]
    if bool(resource_preflight.get("require_default_eaf")):
        resource_keys.insert(0, "default_eaf_file")
    for key in resource_keys:
        for path in _ordered_resource_paths(resources, chromosomes, key):
            if path not in seen_paths:
                all_paths.append(path)
                seen_paths.add(path)

    af_external = optional_text(sample.get("eaffile")) is not None
    af_internal = optional_text(sample.get("eaf_col")) is not None
    if af_external == af_internal:
        raise ValueError(
            "VCF provenance requires exactly one validated EAF source: a study "
            "column or an external file and column."
        )
    if af_external and optional_text(sample.get("eafcolumn")) is None:
        raise ValueError(
            "VCF provenance received an external EAF file without its validated "
            "value column."
        )
    af_column = (
        optional_text(sample.get("eafcolumn"))
        if af_external else optional_text(sample.get("eaf_col"))
    )
    af_source = "external user-provided reference" if af_external else "study column"

    info_source = optional_text(sample.get("info_source")) or missing
    if info_source == "internal":
        if optional_text(sample.get("imp_info_col")) is None:
            raise ValueError(
                "VCF provenance received info_source='internal' without the "
                "validated study INFO column."
            )
        info_source_label = "study column"
        working_info_column = optional_text(sample.get("imp_info_col")) or missing
        info_column = (
            optional_text(sample.get("info_multi_value_source_column"))
            or working_info_column
        )
        info_aggregation = optional_text(
            sample.get("info_multi_value_aggregation")
        )
        if info_aggregation is not None:
            info_interpretation = (
                "Study-provided imputation-quality values separated by %r; "
                "reduced by unweighted row-wise %s into working column %s; "
                "cohort sample-size weights were not inferred"
                % (
                    str(policies.get("info.multi_value_delimiter")),
                    info_aggregation,
                    working_info_column,
                )
            )
        else:
            info_interpretation = "Study-provided imputation-quality score"
        info_matching = "Not applicable"
        info_swap = "Not applicable"
    elif info_source == "external":
        if (
            optional_text(sample.get("infofile")) is None
            or optional_text(sample.get("infocolumn")) is None
        ):
            raise ValueError(
                "VCF provenance received info_source='external' without both "
                "the reference file and value column."
            )
        info_source_label = "external user-provided reference"
        info_column = optional_text(sample.get("infocolumn")) or missing
        info_interpretation = (
            "External reference proxy; not a study-measured imputation-quality score"
        )
        info_matching = (
            "Matched by chromosome, position and direct or swapped alleles"
        )
        info_swap = "INFO value retained unchanged for an allele-swapped match"
    elif info_source == "fixed_cli":
        if sample.get("fixed_info") is None:
            raise ValueError(
                "VCF provenance received info_source='fixed_cli' without the "
                "validated fixed INFO value."
            )
        info_source_label = "user-assigned fixed value"
        info_column = missing
        info_interpretation = (
            "User-assigned constant; not a variant-level measured imputation-quality score"
        )
        info_matching = "Not applicable"
        info_swap = "Not applicable"
    else:
        raise ValueError(
            "VCF provenance requires a validated INFO source; expected "
            "'internal', 'external', or 'fixed_cli', got %r." % info_source
        )
    info_score_type = optional_text(study_decisions.get("info_score_type"))
    if info_score_type == "standard_info":
        info_scale_interpretation = (
            "standard INFO scale; values above %g and no higher than %g were "
            "treated as rounding tolerance and rescaled to %g"
            % (
                float(policies.get("info.clip_max")),
                float(policies.get("info.clip_tolerance")),
                float(policies.get("info.clip_max")),
            )
        )
    elif info_score_type == "mach_rsq":
        info_scale_interpretation = (
            "MaCH Rsq scale; values through the configured maximum %g were "
            "retained without standard-INFO rescaling"
            % float(policies.get("info.mach_rsq_max"))
        )
    else:
        info_scale_interpretation = "INFO score type not available"
    info_interpretation = "%s; %s; decision %s" % (
        info_interpretation,
        info_scale_interpretation,
        _decision_source(study_decisions.get("info_score_type_source")),
    )

    pvalue_type = optional_text(study_decisions.get("pvalue_type"))
    if pvalue_type == "neglog10":
        pvalue_scale = "-log10(P)"
        pvalue_harmonisation = (
            "Converted supplied -log10(P) values to raw P where representable; "
            "preserved exact text and ln(P) below Float64"
        )
        pvalue_formula = "P = 10^(-input)"
    elif pvalue_type == "raw":
        pvalue_scale = "raw P-value"
        pvalue_harmonisation = (
            "Input was already raw P; preserved its source token and ln(P)"
        )
        pvalue_formula = "No harmonisation scale conversion"
    else:
        raise ValueError(
            "VCF provenance requires the validated study-level p-value type; "
            "expected 'raw' or 'neglog10', got %r." % pvalue_type
        )

    trait_type = str(policies.get("sample_size.trait_type"))
    reconstruction = describe_effect_from_z_reconstruction(
        str(policies.get("effect_from_z.method")),
        policies.get("effect_from_z.phenotype_standard_deviation"),
    )
    pvalue_output_field = str(output_fields["pvalue"])
    pvalue_output_tag = pvalue_output_field.split("/", 1)[1]
    values = {
        "postgwas_version": _package_version(),
        "dataset_id": str(sample.get("gwas_outputname") or missing),
        "input_file": _resolved_path(
            sample.get("sumstat_file"), "the summary-statistics input path"
        ),
        "resource_directory": _resolved_path(
            resource_directory, "the resource directory"
        ),
        "output_directory": _resolved_path(
            sample.get("output_root"), "the run output directory"
        ),
        "dataset_output_directory": _resolved_path(
            output_directory, "the dataset output directory"
        ),
        "run_manifest": _resolved_path(run_manifest, "the run-manifest path"),
        "resolved_resource_files": _joined_paths(all_paths, separator, missing),
        "source_fasta_files": paths("genome_fasta_file"),
        "target_fasta_files": paths("target_fasta"),
        "liftover_chain_files": paths("chain_file"),
        "annotation_files": paths("annot_path"),
        "dbsnp_reference_files": paths("dbsnp_file"),
        "strand_consensus": str(study_decisions.get("strand") or missing),
        "strand_reference_column": str(
            next(
                (
                    (resources.get(chromosome) or {}).get(
                        "default_eaf_reference_column"
                    )
                    for chromosome in chromosomes
                    if (resources.get(chromosome) or {}).get(
                        "default_eaf_reference_column"
                    ) is not None
                ),
                missing,
            )
            if bool(resource_preflight.get("require_default_eaf")) else missing
        ),
        "strand_reference_files": (
            paths("default_eaf_file")
            if bool(resource_preflight.get("require_default_eaf")) else missing
        ),
        "af_source": af_source,
        "af_input_column": af_column or missing,
        "af_input_type": str(
            study_decisions.get("frequency_type") or missing
        ),
        "af_type_decision": _decision_source(
            study_decisions.get("eaf_is_maf_source")
        ),
        "af_reference_template": (
            _resolved_path(sample.get("eaffile"), "the external EAF template")
            if af_external else missing
        ),
        "af_reference_files": paths("user_eaf_file") if af_external else missing,
        "af_harmonisation": "Frequency oriented to the final VCF ALT allele",
        "af_swapped_formula": "AF = 1 - input AF for an allele-swapped match",
        "af_output_field": str(output_fields["af"]),
        "info_source": info_source_label,
        "info_input_column": info_column,
        "info_reference_template": (
            _resolved_path(sample.get("infofile"), "the external INFO template")
            if info_source == "external" else missing
        ),
        "info_reference_files": (
            paths("user_info_file") if info_source == "external" else missing
        ),
        "info_fixed_value": str(
            sample.get("fixed_info")
            if sample.get("fixed_info") is not None else missing
        ),
        "info_matching": info_matching,
        "info_swapped_action": info_swap,
        "info_interpretation": info_interpretation,
        "info_output_field": str(output_fields["info"]),
        "se_output_field": str(output_fields["se"]),
        "z_input_column": str(sample.get("imp_z_col") or missing),
        "z_source": (
            "study column"
            if optional_text(sample.get("imp_z_col")) is not None
            else "derived from harmonised BETA and SE"
        ),
        "z_harmonisation": (
            "Supplied Z retained and oriented to the final VCF ALT allele"
            if optional_text(sample.get("imp_z_col")) is not None
            else "Calculated from harmonised BETA and SE"
        ),
        "z_formula": (
            "No Z calculation; supplied Z retained"
            if optional_text(sample.get("imp_z_col")) is not None
            else "Z = BETA / SE"
        ),
        "z_output_field": str(output_fields["z"]),
        "pvalue_input_column": str(sample.get("pval_col") or missing),
        "pvalue_input_scale": pvalue_scale,
        "pvalue_scale_decision": _decision_source(
            study_decisions.get("pvalue_type_source")
        ),
        "pvalue_harmonisation": pvalue_harmonisation,
        "pvalue_harmonisation_formula": pvalue_formula,
        "pvalue_vcf_export": (
            "Converted harmonised raw P-values to %s" % pvalue_output_field
        ),
        "pvalue_vcf_formula": "%s = -log10(P)" % pvalue_output_tag,
        "pvalue_output_field": pvalue_output_field,
        "sample_size_output_field": str(output_fields["sample_size"]),
        "population_af_reference_files": paths("default_comparison_af_file"),
        "population_af_fields": separator.join(
            str(field) for field in vcf_config["external_frequency_columns"]
            if str(field).startswith("INFO/")
        ),
        "vcf_status": "created and structurally validated",
    }
    values.update(
        _effect_provenance(sample, study_decisions, reconstruction, missing)
    )
    values["effect_output_field"] = str(output_fields["effect"])
    values.update(_se_provenance(
        sample,
        study_decisions,
        reconstruction=reconstruction,
        pvalue_tail=int(policies.get("pvalue.se_tail")),
        derive_partial_missing_se=bool(
            policies.get("pvalue.derive_partial_missing_se")
        ),
        missing=missing,
    ))
    if (
        "X" in chromosomes
        and optional_text(sample.get("beta_or_col")) is None
        and optional_text(sample.get("se_col")) is None
        and str(
            policies.get("effect_from_z.x_chromosome_z_only_action")
        ) == "allow_autosomal_assumption"
    ):
        x_assumption = (
            "chromosome X explicitly used the diploid autosomal genotype-"
            "variance assumption under effect_from_z."
            "x_chromosome_z_only_action=allow_autosomal_assumption"
        )
        values["effect_harmonisation"] += "; " + x_assumption
        values["se_harmonisation"] += "; " + x_assumption
    values.update(_sample_size_provenance(sample, trait_type, missing))
    return values


__all__ = ["build_harmonisation_vcf_provenance"]
