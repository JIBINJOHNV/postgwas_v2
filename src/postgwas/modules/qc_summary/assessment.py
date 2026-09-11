"""Fast QC assessment of an input GWAS-VCF without creating a filtered VCF.

The VCF is extracted once with ``bcftools query``. Polars evaluates every
configured QC rule independently against the raw records, combines all active
failure masks, and calculates the final virtual QC-passed metrics in two
streaming aggregation passes over the extracted temporary table.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import polars as pl

from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.core.dataframes import collect_streaming, streaming_collect_options
from postgwas.core.io.artifacts import publish_artifact_set
from postgwas.core.io.reports import (
    write_delimited_report,
    write_html_report,
)
from postgwas.core.paths import configured_output_path
from postgwas.core.polars_runtime import (
    current_polars_thread_runtime,
    run_in_bounded_polars_process,
)
from postgwas.core.ui import StageProgress
from postgwas.core.variant_qc import (
    VariantQCPolicy,
    variant_qc_rule_display_groups,
)
from postgwas.core.vcf import (
    VCF_TAG,
    declared_vcf_metadata_values,
    declared_vcf_tags,
    extract_vcf_table,
    read_vcf_header,
    required_vcf_query_tags,
    validate_vcf_header_contract,
)
from postgwas.modules.qc_summary.reporting import (
    QC_SUMMARY_COLUMNS,
    qc_summary_csv_records,
    render_qc_summary_html,
)


class VcfAssessmentError(RuntimeError):
    """The input VCF could not be converted to a trustworthy QC assessment."""


_QC_PROVENANCE_FIELDS = (
    "postgwas_version",
    "dataset_id",
    "vcf_status",
    "input_genome_build",
    "output_genome_build",
    "liftover",
    "af_source",
    "af_input_column",
    "af_input_type",
    "af_type_decision",
    "af_harmonisation",
    "af_output_field",
    "info_source",
    "info_input_column",
    "info_interpretation",
    "info_output_field",
)

# Build history belongs to harmonisation. QC still reads the harmonised output
# build to reject a contradictory VCF header, but its reports expose only the
# build of the GWAS-VCF that was actually assessed.
_QC_BUILD_HISTORY_PROVENANCE_FIELDS = frozenset({
    "input_genome_build",
    "output_genome_build",
    "liftover",
})


def _reported_scientific_provenance(
    provenance_header_names: Sequence[tuple[str, str]],
    scientific_provenance: Mapping[str, str],
) -> dict[str, Any]:
    """Return scientific field provenance without harmonisation build history."""
    expected = {
        key: header_name
        for key, header_name in provenance_header_names
        if key not in _QC_BUILD_HISTORY_PROVENANCE_FIELDS
    }
    values = {
        key: value
        for key, value in scientific_provenance.items()
        if key in expected
    }
    missing = [key for key in expected if key not in values]
    if not expected:
        status = "not_checked"
    elif not values:
        status = "not_declared"
    elif missing:
        status = "partial"
    else:
        status = "available"
    return {
        "status": status,
        "expected_fields": list(expected),
        "missing_fields": missing,
        "values": values,
    }


@dataclass(frozen=True)
class VcfAssessmentHeaderValidation:
    """Evidence that every configured QC query field is declared by the VCF."""

    vcf_path: str
    genome_build: str
    contigs: tuple[str, ...]
    genome_build_header: str
    supported_genome_builds: tuple[str, ...]
    required_fields: tuple[str, ...]
    declared_info_field_count: int
    declared_format_field_count: int
    provenance_header_names: tuple[tuple[str, str], ...] = ()
    scientific_provenance: dict[str, str] = field(default_factory=dict)

    def as_report(self) -> dict[str, Any]:
        return {
            "status": "passed",
            "genome_build": self.genome_build,
            "genome_build_source": "vcf_header",
            "declared_contig_count": len(self.contigs),
            "required_field_count": len(self.required_fields),
            "required_fields": list(self.required_fields),
            "declared_info_field_count": self.declared_info_field_count,
            "declared_format_field_count": self.declared_format_field_count,
        }

    def provenance_report(self) -> dict[str, Any]:
        """Return optional PostGWAS scientific metadata without enforcing origin."""
        return _reported_scientific_provenance(
            self.provenance_header_names,
            self.scientific_provenance,
        )


_TABLE_COLUMNS = (
    "CHROM", "POS", "REF", "ALT",
    "INFO_AF", "INFO_EXTERNAL_AF", "FORMAT_AF", "FORMAT_SI", "FORMAT_LP",
    "FORMAT_NEF",
)


def extract_vcf_assessment_table(
    vcf_path: str | Path,
    table_path: str | Path,
    dataset_id: str,
    external_af_name: str,
    vcf_fields: dict[str, str],
    table_delimiter: str,
    io_buffer_bytes: int,
    bcftools_bin: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    header_validation: VcfAssessmentHeaderValidation | None = None,
    provenance_headers: Mapping[str, str] | None = None,
    logger=None,
) -> str:
    """Extract the ten fields needed by all QC calculations in one VCF pass."""
    if not VCF_TAG.fullmatch(str(external_af_name)):
        raise VcfAssessmentError(
            "The comparison allele-frequency column must be a valid VCF INFO tag "
            "beginning with a letter: %r" % external_af_name
        )
    vcf = Path(vcf_path).expanduser().resolve()
    table = Path(table_path)
    if not vcf.is_file() or vcf.stat().st_size <= 0:
        raise VcfAssessmentError(
            "Input GWAS-VCF does not exist or is empty: %s" % vcf
        )
    columns = _resolved_vcf_fields(vcf_fields, external_af_name)
    header_validation = _validated_header_evidence(
        header_validation,
        vcf=vcf,
        external_af_name=external_af_name,
        vcf_fields=vcf_fields,
        columns=columns,
        bcftools_bin=bcftools_bin,
        genome_build_header=genome_build_header,
        supported_genome_builds=supported_genome_builds,
        provenance_headers=provenance_headers,
        logger=logger,
    )
    if logger is not None:
        logger.record(
            "INPUT", "vcf_qc_extraction",
            vcf=str(vcf), external_af="INFO/%s" % external_af_name,
            validated_fields=len(header_validation.required_fields),
        )
    return extract_vcf_table(
        vcf, table, dataset_id, columns, bcftools_bin,
        delimiter=table_delimiter,
        io_buffer_bytes=io_buffer_bytes,
        logger=logger,
        error_type=VcfAssessmentError,
        purpose="Extracting VCF fields for QC assessment",
    )


def _resolved_vcf_fields(
    configured: dict[str, str],
    external_af_name: str,
) -> dict[str, str]:
    """Map configured bcftools queries onto the stable internal QC schema."""
    keys = (
        "chromosome", "position", "reference_allele", "alternate_allele",
        "study_info_af", "external_info_af", "study_format_af",
        "imputation_format", "log_pvalue_format", "effective_sample_size_format",
    )
    missing = [key for key in keys if not str(configured.get(key, "")).strip()]
    if missing:
        raise VcfAssessmentError(
            "QC VCF field mapping is missing: %s" % ", ".join(missing)
        )
    values = [
        str(configured[key]).format(external_af=external_af_name)
        for key in keys
    ]
    return dict(zip(_TABLE_COLUMNS, values))


def validate_vcf_assessment_header(
    *,
    vcf_path: str | Path,
    external_af_name: str,
    vcf_fields: dict[str, str],
    bcftools_bin: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    header: str | None = None,
    provenance_headers: Mapping[str, str] | None = None,
    logger=None,
) -> VcfAssessmentHeaderValidation:
    """Infer the build and validate every configured INFO/FORMAT field."""
    if not VCF_TAG.fullmatch(str(external_af_name)):
        raise VcfAssessmentError(
            "The comparison allele-frequency column must be a valid VCF INFO tag "
            "beginning with a letter: %r" % external_af_name
        )
    vcf = Path(vcf_path).expanduser().resolve()
    if not vcf.is_file() or vcf.stat().st_size <= 0:
        raise VcfAssessmentError(
            "Input GWAS-VCF does not exist or is empty: %s" % vcf
        )
    columns = _resolved_vcf_fields(vcf_fields, external_af_name)
    required_fields = tuple(required_vcf_query_tags(columns.values()))
    if header is None:
        header = read_vcf_header(
            vcf,
            bcftools_bin,
            logger=logger,
            error_type=VcfAssessmentError,
        )
    elif not str(header).strip():
        raise VcfAssessmentError("VCF header evidence is empty for %s" % vcf)
    declared = {
        "INFO": declared_vcf_tags(header, "INFO"),
        "FORMAT": declared_vcf_tags(header, "FORMAT"),
    }
    provenance_names = _selected_provenance_headers(provenance_headers)
    try:
        metadata = declared_vcf_metadata_values(
            header, provenance_names.values(),
        )
    except ValueError as exc:
        raise VcfAssessmentError(
            "VCF scientific provenance metadata is malformed: %s" % exc
        ) from exc
    scientific_provenance = {
        logical_name: metadata[header_name]
        for logical_name, header_name in provenance_names.items()
        if header_name in metadata
    }
    missing_fields = []
    for field in required_fields:
        category, tag = field.split("/", 1)
        if tag not in declared[category]:
            missing_fields.append(field)

    provenance_report = _reported_scientific_provenance(
        tuple(provenance_names.items()),
        scientific_provenance,
    )
    log_fields = {
        "vcf": str(vcf),
        "required_field_count": len(required_fields),
        "required_fields": ",".join(required_fields),
        "declared_info_field_count": len(declared["INFO"]),
        "declared_format_field_count": len(declared["FORMAT"]),
        "scientific_provenance_status": provenance_report["status"],
    }
    try:
        genome_build, contigs = validate_vcf_header_contract(
            header=header,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            required_fields=(),
        )
    except ValueError as exc:
        if logger is not None:
            logger.record(
                "VALIDATION",
                "vcf_qc_header_contract",
                status="FAILED",
                failure="genome_build",
                **log_fields,
            )
        raise VcfAssessmentError(
            "QC assessment cannot infer a supported genome build from the VCF "
            "header: %s" % exc
        ) from exc
    log_fields.update({
        "genome_build": genome_build,
        "genome_build_source": "vcf_header",
        "declared_contig_count": len(contigs),
    })
    provenance_output_build = scientific_provenance.get("output_genome_build")
    if (
        provenance_output_build is not None
        and provenance_output_build != genome_build
    ):
        if logger is not None:
            logger.record(
                "VALIDATION",
                "vcf_qc_header_contract",
                status="FAILED",
                failure="provenance_output_build",
                provenance_output_build=provenance_output_build,
                **log_fields,
            )
        raise VcfAssessmentError(
            "QC assessment cannot start because the assessed GWAS-VCF build "
            "%s conflicts with its PostGWAS harmonisation output-build "
            "provenance %s. Provide a VCF with consistent ##genome_build and "
            "##%s metadata."
            % (
                genome_build,
                provenance_output_build,
                provenance_names["output_genome_build"],
            )
        )
    if missing_fields:
        if logger is not None:
            logger.record(
                "VALIDATION",
                "vcf_qc_header_contract",
                status="FAILED",
                missing_fields=",".join(missing_fields),
                **log_fields,
            )
        raise VcfAssessmentError(
            "QC assessment cannot start because the VCF header does not declare "
            "required field(s): %s. Provide a VCF containing these fields or "
            "correct modules.qc_summary.vcf_fields and "
            "modules.qc_summary.reference_af_column."
            % ", ".join(missing_fields)
        )
    if logger is not None:
        logger.record(
            "VALIDATION",
            "vcf_qc_header_contract",
            status="PASSED",
            **log_fields,
        )
    return VcfAssessmentHeaderValidation(
        vcf_path=str(vcf),
        genome_build=genome_build,
        contigs=tuple(contigs),
        genome_build_header=str(genome_build_header),
        supported_genome_builds=tuple(
            str(build) for build in supported_genome_builds
        ),
        required_fields=required_fields,
        declared_info_field_count=len(declared["INFO"]),
        declared_format_field_count=len(declared["FORMAT"]),
        provenance_header_names=tuple(provenance_names.items()),
        scientific_provenance=scientific_provenance,
    )


def _selected_provenance_headers(
    provenance_headers: Mapping[str, str] | None,
) -> dict[str, str]:
    """Select report-relevant names from the configured harmonisation contract."""
    if provenance_headers is None:
        return {}
    missing = [key for key in _QC_PROVENANCE_FIELDS if key not in provenance_headers]
    if missing:
        raise VcfAssessmentError(
            "Configured PostGWAS provenance mapping is missing QC field(s): %s"
            % ", ".join(missing)
        )
    return {
        key: str(provenance_headers[key]) for key in _QC_PROVENANCE_FIELDS
    }


def _validated_header_evidence(
    validation: VcfAssessmentHeaderValidation | None,
    *,
    vcf: Path,
    external_af_name: str,
    vcf_fields: dict[str, str],
    columns: dict[str, str],
    bcftools_bin: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    provenance_headers: Mapping[str, str] | None = None,
    logger=None,
) -> VcfAssessmentHeaderValidation:
    """Reuse one header read only when it proves this exact extraction contract."""
    if validation is None:
        return validate_vcf_assessment_header(
            vcf_path=vcf,
            external_af_name=external_af_name,
            vcf_fields=vcf_fields,
            bcftools_bin=bcftools_bin,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            provenance_headers=provenance_headers,
            logger=logger,
        )
    expected_fields = tuple(required_vcf_query_tags(columns.values()))
    expected_provenance_names = tuple(
        _selected_provenance_headers(provenance_headers).items()
    )
    if (
        validation.vcf_path != str(vcf)
        or validation.genome_build_header != str(genome_build_header)
        or validation.supported_genome_builds
        != tuple(str(build) for build in supported_genome_builds)
        or validation.required_fields != expected_fields
        or validation.provenance_header_names != expected_provenance_names
    ):
        raise VcfAssessmentError(
            "Internal QC header-validation evidence does not match the VCF or "
            "configured extraction fields."
        )
    return validation


def _field_label(query: str) -> str:
    """Render a bcftools query expression as a readable VCF field name."""
    query = query.strip()
    if query.startswith("[") and query.endswith("]"):
        return "FORMAT/" + query[1:-1].lstrip("%")
    return query.lstrip("%")


def resolve_qc_field_labels(
    vcf_fields: dict[str, str],
    external_af_name: str,
) -> dict[str, str]:
    """Resolve configured bcftools queries to user-facing VCF field names."""
    resolved = _resolved_vcf_fields(vcf_fields, external_af_name)
    return {
        "study_info_af": _field_label(resolved["INFO_AF"]),
        "external_info_af": _field_label(resolved["INFO_EXTERNAL_AF"]),
        "study_format_af": _field_label(resolved["FORMAT_AF"]),
        "imputation_format": _field_label(resolved["FORMAT_SI"]),
        "log_pvalue_format": _field_label(resolved["FORMAT_LP"]),
        "effective_sample_size_format": _field_label(resolved["FORMAT_NEF"]),
    }


def _text(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.String, strict=False)


def _number(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.Float64, strict=False)


def _missing(column: str) -> pl.Expr:
    return _text(column).is_null()


def _present(column: str) -> pl.Expr:
    return ~_missing(column)


def _usable_sample_size() -> pl.Expr:
    sample_size = _number("FORMAT_NEF")
    return _present("FORMAT_NEF") & sample_size.is_finite() & (sample_size > 0.0)


def _safe(mask: pl.Expr) -> pl.Expr:
    return mask.fill_null(False)


def _missing_outcome(policy_value: str) -> str:
    return (
        "exclude from the virtual subset"
        if policy_value == "remove"
        else "retain in the virtual subset"
    )


def _count(mask: pl.Expr, name: str) -> pl.Expr:
    return _safe(mask).cast(pl.UInt64).sum().alias(name)


def qc_assessment_runtime(expected_threads: int | None = None) -> dict[str, Any]:
    """Return durable execution metadata for QC aggregation provenance."""
    _, collection_api = streaming_collect_options(
        pl.LazyFrame.collect,
        error_type=VcfAssessmentError,
    )
    runtime = {
        "aggregation_strategy": "two_pass_streaming",
        "temporary_table_scans": 2,
        "streaming_collection_api": collection_api,
        "polars_version": getattr(pl, "__version__", "unknown"),
    }
    if expected_threads is not None:
        runtime.update(current_polars_thread_runtime(expected_threads))
    return runtime


def _resolved_qc_policy(
    configuration: QCSummaryConfig,
    genome_build: str,
) -> VariantQCPolicy:
    """Map the QC schema onto the policy shared with physical filtering."""
    rules = configuration.rules
    mhc_region = (
        configuration.mhc_region(genome_build) if rules.remove_mhc else None
    )
    return VariantQCPolicy(
        minimum_neglog10_p=rules.minimum_neglog10_p,
        missing_pvalue_action=rules.missing_pvalue_action,
        maf_min=rules.maf_min,
        missing_af_action=rules.missing_af_action,
        info_min=rules.info_min,
        info_max=rules.info_max,
        missing_info_action=rules.missing_info_action,
        maximum_af_difference=rules.maximum_af_difference,
        include_indels=rules.include_indels,
        remove_palindromic=rules.remove_palindromic,
        palindromic_lower=rules.palindromic_lower,
        palindromic_upper=rules.palindromic_upper,
        remove_mhc=rules.remove_mhc,
        mhc_chromosome=(None if mhc_region is None else str(mhc_region.chromosome)),
        mhc_start=(None if mhc_region is None else int(mhc_region.start)),
        mhc_end=(None if mhc_region is None else int(mhc_region.end)),
    )


def _all_rows() -> pl.Expr:
    """A true expression with the table's row cardinality, including null CHROM."""
    chromosome = pl.col("CHROM")
    return chromosome.is_null() | chromosome.is_not_null()


def _variant_masks() -> dict[str, pl.Expr]:
    ref = _text("REF").str.to_uppercase()
    alt = _text("ALT").str.to_uppercase()
    bases = ["A", "C", "G", "T"]
    snp = (
        (ref.str.len_chars() == 1)
        & (alt.str.len_chars() == 1)
        & ref.is_in(bases)
        & alt.is_in(bases)
    )
    transition = snp & (
        ((ref == "A") & (alt == "G"))
        | ((ref == "G") & (alt == "A"))
        | ((ref == "C") & (alt == "T"))
        | ((ref == "T") & (alt == "C"))
    )
    palindromic = snp & (
        ((ref == "A") & (alt == "T"))
        | ((ref == "T") & (alt == "A"))
        | ((ref == "C") & (alt == "G"))
        | ((ref == "G") & (alt == "C"))
    )
    return {
        "snp": _safe(snp),
        "transition": _safe(transition),
        "transversion": _safe(snp & ~transition),
        "palindromic": _safe(palindromic),
    }


def build_qc_rules(
    configuration: QCSummaryConfig,
    genome_build: str,
    external_af_name: str,
    field_labels: dict[str, str],
    policy: VariantQCPolicy | None = None,
) -> list[dict[str, Any]]:
    """Build independent rules from the resolved QC configuration."""
    if not VCF_TAG.fullmatch(str(external_af_name)):
        raise VcfAssessmentError(
            "Invalid comparison allele-frequency INFO tag: %r" % external_af_name
        )
    masks = _variant_masks()
    af = _number("FORMAT_AF")
    study_info_af = _number("INFO_AF")
    external_af = _number("INFO_EXTERNAL_AF")
    info = _number("FORMAT_SI")
    lp = _number("FORMAT_LP")
    rules: list[dict[str, Any]] = []
    format_af_label = field_labels["study_format_af"]
    format_info_label = field_labels["imputation_format"]
    format_lp_label = field_labels["log_pvalue_format"]
    study_info_af_label = field_labels["study_info_af"]
    external_info_af_label = field_labels["external_info_af"]
    study_info_af_tag = study_info_af_label.split("/", 1)[-1]
    external_info_af_tag = external_info_af_label.split("/", 1)[-1]
    resolved_fields = _resolved_vcf_fields(
        configuration.vcf_fields.model_dump(), external_af_name,
    )
    rule_display_groups = variant_qc_rule_display_groups({
        "chromosome": _field_label(resolved_fields["CHROM"]),
        "position": _field_label(resolved_fields["POS"]),
        "reference_allele": _field_label(resolved_fields["REF"]),
        "alternate_allele": _field_label(resolved_fields["ALT"]),
        "study_af": format_af_label,
        "imputation_quality": format_info_label,
        "log_pvalue": format_lp_label,
        "study_info_af": study_info_af_label,
        "external_info_af": external_info_af_label,
    })

    policy = policy or _resolved_qc_policy(configuration, genome_build)
    lp_cutoff = policy.minimum_neglog10_p
    if lp_cutoff is not None:
        lp_missing = policy.missing_pvalue_action
        missing_lp = _missing("FORMAT_LP")
        below = _present("FORMAT_LP") & (lp < float(lp_cutoff))
        rules.append({
            "key": "significance",
            "label": "Association significance",
            "category": "Analysis scope",
            "purpose": (
                "Restricts the virtual subset to associations meeting the "
                "configured significance threshold."
            ),
            "criterion": "%s ≥ %.6g; missing values %s"
            % (format_lp_label, float(lp_cutoff), _missing_outcome(lp_missing)),
            **rule_display_groups["significance"],
            "plan_items": [
                {
                    "label": "Missing %s" % format_lp_label,
                    "removes": lp_missing == "remove",
                },
                {
                    "label": "P-value evidence below %s %.6g"
                    % (format_lp_label, float(lp_cutoff)),
                    "removes": True,
                },
            ],
            "failure": below | (missing_lp if lp_missing == "remove" else pl.lit(False)),
            "details": [
                {"key": "missing", "label": "Missing %s" % format_lp_label, "mask": missing_lp,
                 "removes": lp_missing == "remove"},
                {"key": "below_minimum", "label": "Below the LP threshold", "mask": below, "removes": True},
            ],
        })

    maf = float(policy.maf_min)
    af_missing_action = policy.missing_af_action
    missing_format_af = _missing("FORMAT_AF")
    below_af = _present("FORMAT_AF") & (af < maf)
    above_af = _present("FORMAT_AF") & (af > 1.0 - maf)
    outside_maf = below_af | above_af
    rules.append({
        "key": "allele_frequency",
        "label": "Study allele-frequency range",
        "category": "Frequency evidence",
        "purpose": (
            "Checks whether the study ALT-allele frequency implies a minor "
            "allele frequency at or above the configured minimum."
        ),
        "criterion": "%.6g ≤ %s ≤ %.6g; missing values %s"
        % (maf, format_af_label, 1.0 - maf, _missing_outcome(af_missing_action)),
        **rule_display_groups["allele_frequency"],
        "plan_items": [
            {
                "label": "Missing %s" % format_af_label,
                "removes": af_missing_action == "remove",
            },
            {
                "label": "Minor allele frequency outside %.6g ≤ AF ≤ %.6g"
                % (maf, 1.0 - maf),
                "removes": True,
            },
        ],
        "failure": outside_maf | (
            missing_format_af if af_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"key": "missing", "label": "Missing %s" % format_af_label, "mask": missing_format_af,
             "removes": af_missing_action == "remove"},
            {"key": "below_minimum", "label": "%s below %.6g" % (format_af_label, maf),
             "mask": below_af, "removes": True},
            {"key": "above_maximum", "label": "%s above %.6g" % (format_af_label, 1.0 - maf),
             "mask": above_af, "removes": True},
        ],
    })

    info_min = float(policy.info_min)
    info_max = float(policy.info_max)
    info_missing_action = policy.missing_info_action
    missing_info = _missing("FORMAT_SI")
    below_info = _present("FORMAT_SI") & (info < info_min)
    above_info = _present("FORMAT_SI") & (info > info_max)
    outside_info = below_info | above_info
    rules.append({
        "key": "imputation_quality",
        "label": "Imputation quality score",
        "category": "Score quality and completeness",
        "purpose": (
            "Applies the configured imputation-quality score range and "
            "missing-value policy. "
            "The VCF provenance determines whether this is study-measured or "
            "an external proxy."
        ),
        "criterion": "%.6g ≤ %s ≤ %.6g; missing values %s"
        % (info_min, format_info_label, info_max, _missing_outcome(info_missing_action)),
        **rule_display_groups["imputation_quality"],
        "plan_items": [
            {
                "label": "Missing %s" % format_info_label,
                "removes": info_missing_action == "remove",
            },
            {
                "label": "Imputation quality below %s %.6g"
                % (format_info_label, info_min),
                "removes": True,
            },
            {
                "label": "Imputation quality above %s %.6g"
                % (format_info_label, info_max),
                "removes": True,
            },
        ],
        "failure": outside_info | (
            missing_info if info_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"key": "missing", "label": "Missing %s" % format_info_label, "mask": missing_info,
             "removes": info_missing_action == "remove"},
            {"key": "below_minimum", "label": "%s below %.6g" % (format_info_label, info_min),
             "mask": below_info, "removes": True},
            {"key": "above_maximum", "label": "%s above %.6g" % (format_info_label, info_max),
             "mask": above_info, "removes": True},
        ],
    })

    difference_cutoff = float(policy.maximum_af_difference)
    missing_study_info_af = _missing("INFO_AF")
    missing_external_af = _missing("INFO_EXTERNAL_AF")
    discordant = (
        _present("INFO_AF")
        & _present("INFO_EXTERNAL_AF")
        & ((study_info_af - external_af).abs() > difference_cutoff)
    )
    missing_for_comparison = missing_study_info_af | missing_external_af
    rules.append({
        "key": "external_af_concordance",
        "label": "Reference-frequency concordance",
        "category": "Reference concordance",
        "purpose": (
            "Checks agreement between study and configured reference "
            "ALT-allele frequencies on the same allele orientation."
        ),
        "criterion": "|%s − %s| ≤ %.6g; missing values %s"
        % (
            study_info_af_label, external_info_af_label, difference_cutoff,
            _missing_outcome(af_missing_action),
        ),
        **rule_display_groups["external_af_concordance"],
        "plan_items": [
            {
                "label": "Missing %s" % study_info_af_label,
                "removes": af_missing_action == "remove",
            },
            {
                "label": "Missing %s" % external_info_af_label,
                "removes": af_missing_action == "remove",
            },
            {
                "label": "External AF absolute difference |%s − %s| > %.6g"
                % (
                    study_info_af_tag,
                    external_info_af_tag,
                    difference_cutoff,
                ),
                "removes": True,
            },
        ],
        "failure": discordant | (
            missing_for_comparison
            if af_missing_action == "remove" else pl.lit(False)
        ),
        "details": [
            {"key": "missing_study", "label": "Missing %s" % study_info_af_label, "mask": missing_study_info_af,
             "removes": af_missing_action == "remove"},
            {"key": "missing_reference", "label": "Missing %s" % external_info_af_label,
             "mask": missing_external_af, "removes": af_missing_action == "remove"},
            {"key": "difference_above_maximum", "label": "Absolute frequency difference above %.6g" % difference_cutoff,
             "mask": discordant, "removes": True},
        ],
    })

    if not policy.include_indels:
        non_snp = ~masks["snp"]
        rules.append({
            "key": "variant_type",
            "label": "Variant-type scope",
            "category": "Analysis scope",
            "purpose": (
                "Restricts the virtual subset to single-nucleotide variants; "
                "matching non-SNP records are not inherently defective."
            ),
            "criterion": "single-nucleotide variants only",
            **rule_display_groups["variant_type"],
            "plan_items": [{
                "label": "Indels and other non-SNP variants",
                "removes": True,
            }],
            "failure": non_snp,
            "details": [
                {"key": "non_snp", "label": "Indels and other non-SNP variants", "mask": non_snp,
                 "removes": True},
            ],
        })

    if policy.remove_palindromic:
        lower = float(policy.palindromic_lower)
        upper = float(policy.palindromic_upper)
        ambiguous = (
            masks["palindromic"]
            & _present("FORMAT_AF")
            & (af >= lower)
            & (af <= upper)
        )
        rules.append({
            "key": "palindromic_variants",
            "label": "Frequency-ambiguous palindromic variants",
            "category": "Allele compatibility",
            "purpose": (
                "Conservatively excludes palindromic SNPs whose frequency "
                "does not reliably resolve strand orientation."
            ),
            "criterion": "palindromic SNP with %.6g ≤ %s ≤ %.6g"
            % (lower, format_af_label, upper),
            **rule_display_groups["palindromic_variants"],
            "plan_items": [{
                "label": (
                    "Palindromic SNPs with ambiguous frequency "
                    "(%.6g ≤ AF ≤ %.6g)" % (lower, upper)
                ),
                "removes": True,
            }],
            "failure": ambiguous,
            "details": [
                {"key": "frequency_ambiguous", "label": "Frequency-ambiguous palindromic SNPs", "mask": ambiguous,
                 "removes": True},
            ],
        })

    if policy.remove_mhc:
        configured_chromosome = str(policy.mhc_chromosome)
        chromosome = configured_chromosome.lower()
        if chromosome.startswith("chr"):
            chromosome = chromosome[3:]
        start = int(policy.mhc_start)
        end = int(policy.mhc_end)
        normalized_chromosome = _text("CHROM").str.to_lowercase().str.replace(
            r"^chr", ""
        )
        in_mhc = (
            (normalized_chromosome == chromosome)
            & (pl.col("POS").cast(pl.Int64, strict=False) >= start)
            & (pl.col("POS").cast(pl.Int64, strict=False) <= end)
        )
        rules.append({
            "key": "mhc_region",
            "label": "MHC-region scope",
            "category": "Analysis scope",
            "purpose": (
                "Excludes the configured MHC interval as an analysis-policy "
                "choice; matching variants are not automatically low quality."
            ),
            "criterion": "%s:%s-%s inclusive" % (
                configured_chromosome, format(start, ","), format(end, ",")
            ),
            **rule_display_groups["mhc_region"],
            "plan_items": [{
                "label": "Variants in the MHC region (%s:%s-%s)"
                % (
                    configured_chromosome,
                    format(start, ","),
                    format(end, ","),
                ),
                "removes": True,
            }],
            "failure": in_mhc,
            "details": [
                {"key": "inside_region", "label": "Variants inside the configured MHC region", "mask": in_mhc,
                 "removes": True},
            ],
        })

    rule_keys = tuple(rule["key"] for rule in rules)
    if rule_keys != policy.active_rule_keys():
        raise VcfAssessmentError(
            "QC rule construction drifted from the shared scientific-policy "
            "contract: %r != %r" % (rule_keys, policy.active_rule_keys())
        )
    return rules


def _qc_decision_plan(rules: Sequence[Mapping[str, Any]]) -> list[dict[str, str]]:
    """Flatten active QC triggers in the same order as their evaluated rules."""
    plan = []
    for rule in rules:
        items = list(rule.get("plan_items") or ())
        if not items:
            raise VcfAssessmentError(
                "Active QC rule %s has no decision-plan description."
                % rule.get("key", "unknown")
            )
        for item in items:
            plan.append({
                "rule_key": str(rule["key"]),
                "action": "EXCLUDE" if item.get("removes") else "KEEP",
                "label": str(item["label"]),
            })
    return plan


def inactive_qc_rules(policy: VariantQCPolicy) -> list[dict[str, str]]:
    """Explain configurable rules that do not contribute to the combined mask."""
    inactive = []
    if policy.minimum_neglog10_p is None:
        inactive.append({
            "key": "significance",
            "label": "Association significance",
            "category": "Analysis scope",
            "reason": "No minimum −log10(P) threshold is configured.",
        })
    if policy.include_indels:
        inactive.append({
            "key": "variant_type",
            "label": "SNP-only scope",
            "category": "Analysis scope",
            "reason": "Indels and other non-SNP variants are configured for inclusion.",
        })
    if not policy.remove_palindromic:
        inactive.append({
            "key": "palindromic_variants",
            "label": "Frequency-ambiguous palindromic variants",
            "category": "Allele compatibility",
            "reason": "Palindromic-variant exclusion is disabled.",
        })
    if not policy.remove_mhc:
        inactive.append({
            "key": "mhc_region",
            "label": "MHC-region scope",
            "category": "Analysis scope",
            "reason": "MHC-region exclusion is disabled.",
        })
    return inactive


def resolve_qc_rule_contract(
    configuration: QCSummaryConfig,
    genome_build: str,
    external_af_name: str,
) -> dict[str, Any]:
    """Resolve and validate the exact rule metadata used by the QC backend."""
    field_labels = resolve_qc_field_labels(
        configuration.vcf_fields.model_dump(), external_af_name,
    )
    policy = _resolved_qc_policy(configuration, genome_build)
    rules = build_qc_rules(
        configuration,
        genome_build,
        external_af_name,
        field_labels,
        policy,
    )
    return {
        "variant_qc_policy": policy.as_dict(),
        "field_labels": field_labels,
        "active_rule_count": len(rules),
        "decision_plan": _qc_decision_plan(rules),
        "rules": [
            {
                "key": rule["key"],
                "label": rule["label"],
                "category": rule["category"],
                "display_group": rule["display_group"],
                "display_group_kind": rule["display_group_kind"],
                "purpose": rule["purpose"],
                "criterion": rule["criterion"],
                "details": [
                    {
                        "key": detail["key"],
                        "label": detail["label"],
                        "decision": (
                            "exclude_from_virtual_subset"
                            if detail.get("removes")
                            else "retain_by_configured_policy"
                        ),
                    }
                    for detail in rule["details"]
                ],
            }
            for rule in rules
        ],
        "inactive_rules": inactive_qc_rules(policy),
    }


def validate_qc_rule_contract(
    assessment: Mapping[str, Any],
    rule_contract: Mapping[str, Any],
) -> None:
    """Prove that assessed rule metadata matches the pre-conversion contract."""
    scientific_rule_keys = (
        "key", "label", "category", "purpose", "criterion", "details",
    )

    def scientific_rules(source: Sequence[Mapping[str, Any]]) -> list[dict]:
        projected = [
            {key: rule[key] for key in scientific_rule_keys}
            for rule in source
        ]
        for rule in projected:
            rule["details"] = [
                {
                    key: detail[key]
                    for key in ("key", "label", "decision")
                }
                for detail in rule["details"]
            ]
        return projected

    observed_contract = {
        "variant_qc_policy": assessment["variant_qc_policy"],
        "field_labels": assessment["field_labels"],
        "active_rule_count": assessment["active_rule_count"],
        "rules": scientific_rules(assessment["rules"]),
        "inactive_rules": assessment["inactive_rules"],
    }
    # The numbered decision plan is derived presentation metadata.  Compare the
    # scientific contract itself so completed assessments written before that
    # display metadata was introduced remain valid resume inputs.
    expected_contract = dict(rule_contract)
    expected_contract["rules"] = scientific_rules(rule_contract["rules"])
    expected_contract.pop("decision_plan", None)
    if observed_contract != expected_contract:
        raise VcfAssessmentError(
            "The QC rules recorded by the assessment do not match the rules "
            "resolved before VCF conversion."
        )


def _metric_expressions(
    prefix: str,
    subset: pl.Expr,
    sample_size_statistics: dict[str, Any],
    external_difference_cutoff: float,
    sample_size_sd_multiplier: float,
    sample_size_reference_value: float | None,
    sample_size_minimum_threshold: float | None,
) -> list[pl.Expr]:
    masks = _variant_masks()
    study_af_missing = _missing("INFO_AF")
    external_af_missing = _missing("INFO_EXTERNAL_AF")
    comparable = ~study_af_missing & ~external_af_missing
    discordant = comparable & (
        (_number("INFO_AF") - _number("INFO_EXTERNAL_AF")).abs()
        > external_difference_cutoff
    )
    sample_size = _number("FORMAT_NEF")
    usable_sample_size = _usable_sample_size()
    stage_sample_size = _safe(subset) & usable_sample_size
    mean = pl.lit(
        sample_size_statistics["__%s_sample_size_mean" % prefix],
        dtype=pl.Float64,
    )
    standard_deviation = pl.lit(
        sample_size_statistics["__%s_sample_size_sd" % prefix],
        dtype=pl.Float64,
    )
    outlier_threshold = (
        mean + sample_size_sd_multiplier * standard_deviation
    )
    low_reference = pl.lit(sample_size_reference_value, dtype=pl.Float64)
    minimum_threshold = pl.lit(
        sample_size_minimum_threshold,
        dtype=pl.Float64,
    )
    return [
        _count(subset, prefix + "__num_records"),
        _count(subset & masks["snp"], prefix + "__num_snps"),
        _count(subset & ~masks["snp"], prefix + "__num_non_snps"),
        _count(subset & masks["transition"], prefix + "__transitions"),
        _count(subset & masks["transversion"], prefix + "__transversions"),
        _count(subset & _missing("FORMAT_AF"), prefix + "__format_af_missing"),
        _count(subset & _missing("FORMAT_SI"), prefix + "__format_si_missing"),
        _count(subset & study_af_missing, prefix + "__study_af_missing"),
        _count(subset & external_af_missing, prefix + "__external_af_missing"),
        _count(subset & comparable, prefix + "__af_comparable"),
        _count(subset & discordant, prefix + "__af_difference_above_cutoff"),
        _count(stage_sample_size, prefix + "__effective_sample_size_available"),
        _count(
            _safe(subset) & ~usable_sample_size,
            prefix + "__effective_sample_size_missing_or_invalid",
        ),
        pl.when(stage_sample_size).then(sample_size).min().alias(
            prefix + "__effective_sample_size_minimum"
        ),
        pl.when(stage_sample_size).then(sample_size).max().alias(
            prefix + "__effective_sample_size_maximum"
        ),
        mean.alias(prefix + "__effective_sample_size_mean"),
        standard_deviation.alias(
            prefix + "__effective_sample_size_standard_deviation"
        ),
        outlier_threshold.alias(
            prefix + "__effective_sample_size_outlier_threshold"
        ),
        _count(
            stage_sample_size & (sample_size > outlier_threshold),
            prefix + "__effective_sample_size_above_outlier_threshold",
        ),
        low_reference.alias(
            prefix + "__effective_sample_size_reference_quantile_value"
        ),
        minimum_threshold.alias(
            prefix + "__effective_sample_size_minimum_threshold"
        ),
        _count(
            stage_sample_size & (sample_size < minimum_threshold),
            prefix + "__effective_sample_size_below_minimum_threshold",
        ),
    ]


def _sample_size_statistics(
    frame: pl.LazyFrame,
    stage_subsets: dict[str, pl.Expr],
    reference_quantile: float,
) -> dict[str, Any]:
    """Collect stage Neff scalars and the raw reference quantile."""
    sample_size = _number("FORMAT_NEF")
    usable = _usable_sample_size()
    expressions = []
    for prefix, subset in stage_subsets.items():
        values = pl.when(_safe(subset) & usable).then(sample_size)
        expressions.extend([
            values.mean().alias("__%s_sample_size_mean" % prefix),
            values.std(ddof=1).alias("__%s_sample_size_sd" % prefix),
        ])
    raw_values = pl.when(usable).then(sample_size)
    expressions.append(
        raw_values.quantile(
            reference_quantile,
            interpolation="linear",
        ).alias("__raw_sample_size_reference_quantile")
    )
    return collect_streaming(
        frame.select(expressions),
        error_type=VcfAssessmentError,
    ).to_dicts()[0]


def assess_variant_table(
    table_path: str | Path,
    configuration: QCSummaryConfig,
    genome_build: str,
    external_af_name: str,
    vcf_fields: dict[str, str],
    *,
    delimiter: str,
    null_values: list[str],
    expected_threads: int | None = None,
    logger=None,
) -> dict[str, Any]:
    """Calculate raw-rule and combined QC metrics in two streaming table scans."""
    table = Path(table_path)
    schema = dict((column, pl.String) for column in _TABLE_COLUMNS)
    frame = pl.scan_csv(
        table,
        separator=delimiter,
        null_values=null_values,
        schema_overrides=schema,
        infer_schema_length=0,
    ).with_columns(pl.col("POS").cast(pl.Int64, strict=False))
    field_labels = resolve_qc_field_labels(vcf_fields, external_af_name)
    policy = _resolved_qc_policy(configuration, genome_build)
    rules = build_qc_rules(
        configuration, genome_build, external_af_name, field_labels, policy,
    )
    difference_cutoff = float(policy.maximum_af_difference)
    sample_size_sd_multiplier = float(
        configuration.rules.sample_size_outlier_standard_deviations
    )
    sample_size_reference_quantile = float(
        configuration.rules.sample_size_reference_quantile
    )
    sample_size_minimum_fraction = float(
        configuration.rules.sample_size_minimum_fraction_of_reference
    )
    rule_selections = []
    rule_failures = []
    combined_failure = pl.lit(False)
    failure_count = pl.lit(0, dtype=pl.UInt16)
    for index, rule in enumerate(rules, 1):
        failure = _safe(rule["failure"])
        rule_failures.append(failure)
        for detail_index, detail in enumerate(rule["details"], 1):
            detail_mask = _safe(detail["mask"])
            rule_selections.extend([
                _count(
                    detail_mask,
                    "rule_%02d__detail_%02d__matched_raw" % (index, detail_index),
                ),
            ])
        combined_failure = combined_failure | failure
        failure_count = failure_count + failure.cast(pl.UInt16)

    for index, failure in enumerate(rule_failures, 1):
        rule_selections.extend([
            _count(failure, "rule_%02d__failed_raw" % index),
            _count(
                failure & (failure_count == 1),
                "rule_%02d__unique_only_raw" % index,
            ),
            _count(
                failure & (failure_count > 1),
                "rule_%02d__overlap_raw" % index,
            ),
        ])

    combined_failure_column = "__combined_qc_failure"
    frame = frame.with_columns(
        _safe(combined_failure).alias(combined_failure_column)
    )
    qc_passed_subset = ~pl.col(combined_failure_column)
    stage_subsets = {"raw": _all_rows(), "qc_passed": qc_passed_subset}
    aggregation = qc_assessment_runtime(expected_threads)
    if logger is not None:
        logger.record("PARAM", "qc_summary_aggregation", **aggregation)
    sample_size_statistics = _sample_size_statistics(
        frame,
        stage_subsets,
        sample_size_reference_quantile,
    )
    sample_size_reference_value = sample_size_statistics[
        "__raw_sample_size_reference_quantile"
    ]
    sample_size_minimum_threshold = (
        None
        if sample_size_reference_value is None
        else float(sample_size_reference_value) * sample_size_minimum_fraction
    )
    selections = _metric_expressions(
        "raw", stage_subsets["raw"], sample_size_statistics,
        difference_cutoff, sample_size_sd_multiplier,
        sample_size_reference_value, sample_size_minimum_threshold,
    )
    selections.extend(rule_selections)
    selections.extend([
        _count(failure_count > 1, "assessment__overlap_variants"),
        (failure_count.cast(pl.Int64) - 1)
        .clip(lower_bound=0)
        .sum()
        .alias("assessment__extra_rule_matches"),
    ])
    selections.extend(_metric_expressions(
        "qc_passed", stage_subsets["qc_passed"], sample_size_statistics,
        difference_cutoff, sample_size_sd_multiplier,
        sample_size_reference_value, sample_size_minimum_threshold,
    ))
    values = collect_streaming(
        frame.select(selections),
        error_type=VcfAssessmentError,
    ).to_dicts()[0]

    def stage_metrics(prefix: str) -> dict[str, Any]:
        integer_metrics = {
            "num_records", "num_snps", "num_non_snps", "transitions",
            "transversions", "format_af_missing", "format_si_missing",
            "study_af_missing", "external_af_missing", "af_comparable",
            "af_difference_above_cutoff", "effective_sample_size_available",
            "effective_sample_size_missing_or_invalid",
            "effective_sample_size_above_outlier_threshold",
            "effective_sample_size_below_minimum_threshold",
        }
        metrics = {}
        for key, value in values.items():
            if not key.startswith(prefix + "__"):
                continue
            metric = key.split("__", 1)[1]
            metrics[metric] = (
                int(value or 0) if metric in integer_metrics
                else (None if value is None else float(value))
            )
        transversions = metrics.get("transversions", 0)
        metrics["ts_tv_ratio"] = (
            metrics.get("transitions", 0) / transversions
            if transversions else None
        )
        usable_sample_sizes = metrics.get("effective_sample_size_available", 0)
        low_sample_sizes = metrics.get(
            "effective_sample_size_below_minimum_threshold", 0,
        )
        metrics["effective_sample_size_below_minimum_threshold_fraction"] = (
            low_sample_sizes / usable_sample_sizes
            if usable_sample_sizes else None
        )
        return metrics

    rule_results = []
    for index, rule in enumerate(rules, 1):
        prefix = "rule_%02d__" % index
        details = []
        for detail_index, detail in enumerate(rule["details"], 1):
            detail_prefix = prefix + "detail_%02d__" % detail_index
            details.append({
                "key": detail["key"],
                "label": detail["label"],
                "decision": (
                    "exclude_from_virtual_subset"
                    if detail.get("removes")
                    else "retain_by_configured_policy"
                ),
                "matched_raw": int(values.get(detail_prefix + "matched_raw") or 0),
            })
        rule_results.append({
            "number": index,
            "key": rule["key"],
            "label": rule["label"],
            "category": rule["category"],
            "display_group": rule["display_group"],
            "display_group_kind": rule["display_group_kind"],
            "purpose": rule["purpose"],
            "criterion": rule["criterion"],
            "decision": "exclude_from_virtual_subset",
            "failed_raw": int(values.get(prefix + "failed_raw") or 0),
            "unique_only_raw": int(
                values.get(prefix + "unique_only_raw") or 0
            ),
            "overlap_raw": int(values.get(prefix + "overlap_raw") or 0),
            "details": details,
        })

    raw = stage_metrics("raw")
    qc_passed = stage_metrics("qc_passed")
    rule_match_total = sum(rule["failed_raw"] for rule in rule_results)
    excluded_total = raw["num_records"] - qc_passed["num_records"]
    extra_rule_matches = int(values.get("assessment__extra_rule_matches") or 0)
    for rule in rule_results:
        rule["failed_fraction_raw"] = (
            rule["failed_raw"] / raw["num_records"]
            if raw["num_records"] else 0.0
        )
        if rule["failed_raw"] != rule["unique_only_raw"] + rule["overlap_raw"]:
            raise VcfAssessmentError(
                "QC rule overlap accounting failed for %s: %s matches do not "
                "equal %s unique-only plus %s overlapping variants."
                % (
                    rule["key"], rule["failed_raw"], rule["unique_only_raw"],
                    rule["overlap_raw"],
                )
            )
        for detail in rule["details"]:
            detail["matched_fraction_raw"] = (
                detail["matched_raw"] / raw["num_records"]
                if raw["num_records"] else 0.0
            )
    unique_rule_only_total = sum(
        rule["unique_only_raw"] for rule in rule_results
    )
    overlapping_rule_match_total = sum(
        rule["overlap_raw"] for rule in rule_results
    )
    overlap_variants = int(values.get("assessment__overlap_variants") or 0)
    if excluded_total != unique_rule_only_total + overlap_variants:
        raise VcfAssessmentError(
            "QC attribution accounting failed: %s excluded variants do not "
            "equal %s unique-to-one-rule plus %s multi-rule variants."
            % (excluded_total, unique_rule_only_total, overlap_variants)
        )
    if overlapping_rule_match_total != overlap_variants + extra_rule_matches:
        raise VcfAssessmentError(
            "QC overlap attribution failed: %s overlapping rule matches do not "
            "equal %s multi-rule variants plus %s additional matches."
            % (overlapping_rule_match_total, overlap_variants, extra_rule_matches)
        )
    if raw["num_records"] != excluded_total + qc_passed["num_records"]:
        raise VcfAssessmentError(
            "QC accounting failed: %s raw records do not equal %s excluded plus "
            "%s retained. No assessment report was accepted."
            % (raw["num_records"], excluded_total, qc_passed["num_records"])
        )
    if rule_match_total != excluded_total + extra_rule_matches:
        raise VcfAssessmentError(
            "QC rule accounting failed: %s raw-rule matches do not equal %s unique "
            "failed variants plus %s overlapping matches."
            % (rule_match_total, excluded_total, extra_rule_matches)
        )
    # Aggregate match multiplicity is retained only for these integrity checks.
    # Public reports use distinct variant counts, which cannot be mistaken for
    # additional excluded variants.
    return {
        "definition": (
            "all active configured QC conditions are applied together to every "
            "input GWAS-VCF record; a variant belongs to the virtual QC-passed "
            "subset only when it passes every condition. The input VCF is "
            "unchanged and no filtered VCF is created"
        ),
        "raw": raw,
        "rules": rule_results,
        "decision_plan": _qc_decision_plan(rules),
        "inactive_rules": inactive_qc_rules(policy),
        "active_rule_count": len(rule_results),
        "qc_passed": qc_passed,
        "excluded_total": excluded_total,
        "unique_rule_only_total": unique_rule_only_total,
        "accounting_balanced": True,
        "retained_fraction": (
            qc_passed["num_records"] / raw["num_records"]
            if raw["num_records"] else 0.0
        ),
        "overlap_variants": overlap_variants,
        "external_af_name": external_af_name,
        "variant_qc_policy": policy.as_dict(),
        "field_labels": field_labels,
        "af_difference_cutoff": difference_cutoff,
        "sample_size_outlier_standard_deviations": sample_size_sd_multiplier,
        "sample_size_reference_quantile": sample_size_reference_quantile,
        "sample_size_minimum_fraction_of_reference": (
            sample_size_minimum_fraction
        ),
        "aggregation": aggregation,
    }


@contextmanager
def qc_assessment_stage(
    progress: StageProgress,
    logger,
    number: int,
    total: int,
    title: str,
):
    """Record one truthful QC stage without nesting the caller's logger steps."""
    started = time.monotonic()
    if logger is not None:
        logger.record(
            "STEP", "qc_summary_stage", number=number, total=total, name=title,
        )
    try:
        with progress.step(number, total, title):
            yield
    except BaseException as exc:
        if logger is not None:
            logger.record(
                "STATUS",
                "qc_summary_stage",
                number=number,
                total=total,
                name=title,
                status="FAILED",
                duration_seconds=time.monotonic() - started,
                exception_type=type(exc).__name__,
            )
        raise
    else:
        if logger is not None:
            logger.record(
                "STATUS",
                "qc_summary_stage",
                number=number,
                total=total,
                name=title,
                status="COMPLETED",
                duration_seconds=time.monotonic() - started,
            )


def _metric_report_records(
    assessment: dict[str, Any], null_output: str,
) -> list[dict[str, Any]]:
    records = []
    header_validation = assessment["vcf_header_validation"]
    for metric in (
        "status", "genome_build", "genome_build_source",
        "declared_contig_count", "required_field_count", "required_fields",
        "declared_info_field_count", "declared_format_field_count",
    ):
        value = header_validation[metric]
        if isinstance(value, list):
            value = ",".join(value)
        records.append({
            "stage": "vcf_header_validation",
            "metric": metric,
            "value": null_output if value is None else value,
        })
    provenance = assessment.get("vcf_provenance") or {}
    records.append({
        "stage": "vcf_scientific_provenance",
        "metric": "status",
        "value": provenance.get("status", "not_checked"),
    })
    for key, value in (provenance.get("values") or {}).items():
        records.append({
            "stage": "vcf_scientific_provenance",
            "metric": key,
            "value": null_output if value is None else value,
        })
    for stage in ("raw", "qc_passed"):
        for metric, value in assessment[stage].items():
            records.append({
                "stage": stage,
                "metric": metric,
                "value": null_output if value is None else value,
            })
    for metric, value in assessment["aggregation"].items():
        records.append({
            "stage": "aggregation",
            "metric": metric,
            "value": null_output if value is None else value,
        })
    for metric in (
        "active_rule_count", "excluded_total", "unique_rule_only_total",
        "accounting_balanced", "retained_fraction", "overlap_variants",
        "sample_size_outlier_standard_deviations",
        "sample_size_reference_quantile",
        "sample_size_minimum_fraction_of_reference", "definition",
    ):
        value = assessment[metric]
        records.append({
            "stage": "assessment",
            "metric": metric,
            "value": null_output if value is None else value,
        })
    return records


def _rule_report_records(assessment: dict[str, Any]) -> list[dict[str, Any]]:
    records = []
    for rule in assessment["rules"]:
        records.append({
            "record_type": "rule",
            "rule_number": rule["number"],
            "category": rule["category"],
            "rule": rule["label"],
            "purpose": rule["purpose"],
            "criterion": rule["criterion"],
            "decision": rule["decision"],
            "variants_matching_in_raw_vcf": rule["failed_raw"],
            "fraction_matching_in_raw_vcf": rule["failed_fraction_raw"],
            "variants_unique_to_rule": rule["unique_only_raw"],
            "variants_overlapping_other_rules": rule["overlap_raw"],
        })
        for detail in rule["details"]:
            records.append({
                "record_type": "detail",
                "rule_number": rule["number"],
                "category": rule["category"],
                "rule": rule["label"],
                "purpose": rule["purpose"],
                "criterion": rule["criterion"],
                "detail_key": detail["key"],
                "detail": detail["label"],
                "decision": detail["decision"],
                "variants_matching_in_raw_vcf": detail["matched_raw"],
                "fraction_matching_in_raw_vcf": detail[
                    "matched_fraction_raw"
                ],
            })
    for inactive in assessment.get("inactive_rules") or ():
        records.append({
            "record_type": "inactive_rule",
            "category": inactive["category"],
            "rule": inactive["label"],
            "purpose": inactive["reason"],
            "decision": "inactive",
        })
    return records


def _write_assessment_reports(
    assessment: dict[str, Any],
    report_paths: dict[str, Path],
    *,
    delimiter: str,
    null_output: str,
    progress: StageProgress,
    logger=None,
    stage_number_offset: int = 0,
) -> dict[str, str]:
    for path in report_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    staged_paths = {}
    for name, destination in report_paths.items():
        handle = tempfile.NamedTemporaryFile(
            prefix=".%s." % destination.name,
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        )
        handle.close()
        staged_paths[name] = Path(handle.name)

    try:
        with qc_assessment_stage(
            progress,
            logger,
            stage_number_offset + 3,
            stage_number_offset + 4,
            "Generate and validate the QC reports",
        ):
            write_delimited_report(
                _metric_report_records(assessment, null_output),
                staged_paths["summary"],
                fieldnames=("stage", "metric", "value"),
                delimiter=delimiter,
                null_value=null_output,
            )
            write_delimited_report(
                _rule_report_records(assessment),
                staged_paths["rules"],
                fieldnames=(
                    "record_type", "rule_number", "category", "rule",
                    "purpose", "criterion", "detail_key", "detail", "decision",
                    "variants_matching_in_raw_vcf",
                    "fraction_matching_in_raw_vcf",
                    "variants_unique_to_rule",
                    "variants_overlapping_other_rules",
                ),
                delimiter=delimiter,
                null_value=null_output,
            )
            staged_paths["json"].write_text(
                json.dumps(assessment, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            write_delimited_report(
                qc_summary_csv_records(assessment),
                staged_paths["csv"],
                fieldnames=QC_SUMMARY_COLUMNS,
                delimiter=",",
                null_value=null_output,
            )
            write_html_report(
                render_qc_summary_html(assessment), staged_paths["html"],
            )
            invalid = [
                str(path) for path in staged_paths.values()
                if not path.is_file() or path.stat().st_size <= 0
            ]
            if invalid:
                raise VcfAssessmentError(
                    "QC report rendering produced empty artifact(s): %s"
                    % ", ".join(invalid)
                )
            try:
                json.loads(staged_paths["json"].read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise VcfAssessmentError(
                    "Rendered QC assessment JSON is invalid: %s" % exc
                ) from exc

        with qc_assessment_stage(
            progress,
            logger,
            stage_number_offset + 4,
            stage_number_offset + 4,
            "Publish the validated QC report set",
        ):
            try:
                publish_artifact_set(tuple(
                    (staged_paths[name], destination)
                    for name, destination in report_paths.items()
                ))
            except RuntimeError as exc:
                raise VcfAssessmentError(str(exc)) from exc
    finally:
        for path in staged_paths.values():
            path.unlink(missing_ok=True)
    return {name: str(path) for name, path in report_paths.items()}


def run_vcf_qc_assessment(
    *,
    vcf_path: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    external_af_name: str,
    configuration: QCSummaryConfig,
    bcftools_bin: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    threads: int,
    header_validation: VcfAssessmentHeaderValidation | None = None,
    provenance_headers: Mapping[str, str] | None = None,
    logger=None,
    stage_progress: StageProgress | None = None,
    stage_number_offset: int = 0,
    rule_contract: Mapping[str, Any] | None = None,
    metrics_completed: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Extract once, assess raw rules, persist summaries, and remove the TSV."""
    progress = stage_progress or StageProgress(
        "GWAS-VCF QC assessment progress", enabled=False,
    )
    vcf_fields = configuration.vcf_fields.model_dump()
    resolved_vcf = Path(vcf_path).expanduser().resolve()
    columns = _resolved_vcf_fields(vcf_fields, external_af_name)
    header_validation = _validated_header_evidence(
        header_validation,
        vcf=resolved_vcf,
        external_af_name=external_af_name,
        vcf_fields=vcf_fields,
        columns=columns,
        bcftools_bin=bcftools_bin,
        genome_build_header=genome_build_header,
        supported_genome_builds=supported_genome_builds,
        provenance_headers=provenance_headers,
        logger=logger,
    )
    genome_build = header_validation.genome_build
    output_layout = configuration.output_layout
    table_configuration = configuration.table
    report_paths = {
        "summary": configured_output_path(
            output_directory, output_layout.metric_report,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "rules": configured_output_path(
            output_directory, output_layout.rule_report,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "json": configured_output_path(
            output_directory, output_layout.assessment_json,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "csv": configured_output_path(
            output_directory, output_layout.summary_csv,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
        "html": configured_output_path(
            output_directory, output_layout.html_report,
            error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
        ),
    }
    temporary_prefix = configured_output_path(
        output_directory, output_layout.temporary_table_prefix,
        error_type=VcfAssessmentError, dataset_id=dataset_id, build=genome_build,
    )
    temporary_prefix.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=temporary_prefix.name,
        suffix=table_configuration.temporary_suffix,
        dir=temporary_prefix.parent,
        delete=False,
    )
    table_path = Path(handle.name)
    handle.close()
    try:
        with qc_assessment_stage(
            progress,
            logger,
            stage_number_offset + 1,
            stage_number_offset + 4,
            "Convert required GWAS-VCF fields to a temporary TSV",
        ):
            extract_vcf_assessment_table(
                vcf_path=vcf_path,
                table_path=table_path,
                dataset_id=dataset_id,
                external_af_name=external_af_name,
                vcf_fields=vcf_fields,
                table_delimiter=table_configuration.delimiter,
                io_buffer_bytes=table_configuration.io_buffer_bytes,
                bcftools_bin=bcftools_bin,
                genome_build_header=genome_build_header,
                supported_genome_builds=supported_genome_builds,
                header_validation=header_validation,
                provenance_headers=provenance_headers,
                logger=logger,
            )
        with qc_assessment_stage(
            progress,
            logger,
            stage_number_offset + 2,
            stage_number_offset + 4,
            "Evaluate the configured QC rules and calculate metrics",
        ):
            assessment = run_in_bounded_polars_process(
                assess_variant_table,
                table_path,
                configuration,
                genome_build,
                external_af_name,
                vcf_fields,
                threads=threads,
                call_kwargs={
                    "delimiter": table_configuration.delimiter,
                    "null_values": table_configuration.null_values,
                    "expected_threads": threads,
                },
            )
        if logger is not None:
            logger.record(
                "PARAM", "qc_summary_aggregation", **assessment["aggregation"]
            )
        assessment["raw_vcf"] = str(resolved_vcf)
        assessment["dataset_id"] = str(dataset_id)
        assessment["genome_build"] = str(genome_build)
        assessment["vcf_header_validation"] = header_validation.as_report()
        assessment["vcf_provenance"] = header_validation.provenance_report()
        assessment["reports"] = {
            name: str(path) for name, path in report_paths.items()
        }
        if rule_contract is not None:
            validate_qc_rule_contract(assessment, rule_contract)
        if metrics_completed is not None:
            metrics_completed(assessment)
        _write_assessment_reports(
            assessment,
            report_paths,
            delimiter=table_configuration.delimiter,
            null_output=table_configuration.null_output,
            progress=progress,
            logger=logger,
            stage_number_offset=stage_number_offset,
        )
        if logger is not None:
            provenance = assessment["vcf_provenance"]
            logger.record(
                "INPUT",
                "vcf_scientific_provenance",
                status=provenance["status"],
                missing_fields=",".join(provenance["missing_fields"]),
                **provenance["values"],
            )
            for rule in assessment["rules"]:
                logger.record(
                    "RESULT",
                    "qc_summary_rule",
                    number=rule["number"],
                    key=rule["key"],
                    category=rule["category"],
                    criterion=rule["criterion"],
                    decision=rule["decision"],
                    variants_matching=rule["failed_raw"],
                    fraction_matching=rule["failed_fraction_raw"],
                    variants_unique_to_rule=rule["unique_only_raw"],
                    variants_overlapping_other_rules=rule["overlap_raw"],
                )
            for inactive in assessment["inactive_rules"]:
                logger.record(
                    "RESULT",
                    "qc_summary_inactive_rule",
                    key=inactive["key"],
                    category=inactive["category"],
                    reason=inactive["reason"],
                )
            logger.record(
                "RESULT", "vcf_qc_assessment",
                raw=assessment["raw"]["num_records"],
                qc_passed=assessment["qc_passed"]["num_records"],
                excluded=assessment["excluded_total"],
                variants_unique_to_one_rule=assessment[
                    "unique_rule_only_total"
                ],
                variants_matching_multiple_rules=assessment[
                    "overlap_variants"
                ],
                sample_size_outlier_standard_deviations=assessment[
                    "sample_size_outlier_standard_deviations"
                ],
                raw_sample_size_outliers=assessment["raw"][
                    "effective_sample_size_above_outlier_threshold"
                ],
                qc_passed_sample_size_outliers=assessment["qc_passed"][
                    "effective_sample_size_above_outlier_threshold"
                ],
                sample_size_reference_quantile=assessment[
                    "sample_size_reference_quantile"
                ],
                sample_size_reference_value=assessment["raw"][
                    "effective_sample_size_reference_quantile_value"
                ],
                sample_size_minimum_fraction_of_reference=assessment[
                    "sample_size_minimum_fraction_of_reference"
                ],
                sample_size_minimum_threshold=assessment["raw"][
                    "effective_sample_size_minimum_threshold"
                ],
                raw_sample_size_below_minimum_threshold=assessment["raw"][
                    "effective_sample_size_below_minimum_threshold"
                ],
                raw_sample_size_below_minimum_threshold_fraction=assessment[
                    "raw"
                ]["effective_sample_size_below_minimum_threshold_fraction"],
                qc_passed_sample_size_below_minimum_threshold=assessment[
                    "qc_passed"
                ]["effective_sample_size_below_minimum_threshold"],
                qc_passed_sample_size_below_minimum_threshold_fraction=(
                    assessment["qc_passed"][
                        "effective_sample_size_below_minimum_threshold_fraction"
                    ]
                ),
                vcf_header_contract=assessment["vcf_header_validation"]["status"],
                required_vcf_fields=assessment["vcf_header_validation"][
                    "required_field_count"
                ],
                aggregation_strategy=assessment["aggregation"][
                    "aggregation_strategy"
                ],
                temporary_table_scans=assessment["aggregation"][
                    "temporary_table_scans"
                ],
                streaming_collection_api=assessment["aggregation"][
                    "streaming_collection_api"
                ],
                polars_version=assessment["aggregation"]["polars_version"],
                vcf_output="input_unchanged_no_filtered_vcf",
            )
        return assessment
    finally:
        try:
            table_path.unlink()
        except OSError:
            pass
