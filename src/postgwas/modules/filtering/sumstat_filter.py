"""Filter a harmonised GWAS-VCF with a validated bcftools configuration."""

import csv
import io
import json
import os
import re
import shlex
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import polars as pl

from postgwas.core.dataframes import collect_streaming
from postgwas.core.execution.runtime import run_cmd
from postgwas.core.input_validation import (
    current_validation_recorder, current_validation_session, record_file_validation,
)
from postgwas.core.io.artifacts import publish_artifact_set
from postgwas.core.paths import configured_output_path, validate_filename_component
from postgwas.core.polars_runtime import (
    current_polars_thread_runtime,
    run_in_bounded_polars_process,
)
from postgwas.core.ui import StageProgress, style_screen_block
from postgwas.core.variant_qc import (
    VariantQCPolicy,
    variant_qc_rule_display_groups,
)
from postgwas.core.vcf import (
    VCF_TAG,
    count_indexed_vcf_records,
    declared_vcf_metadata_values,
    declared_vcf_tags,
    read_vcf_header,
    required_vcf_field_presence,
    validate_postgwas_vcf_provenance,
    validate_vcf_header_contract,
)
from postgwas.modules.filtering.reporting import (
    build_filtering_summary,
    render_filtering_plan,
    render_input_vcf_validation,
    render_filtering_summary,
    write_filtering_html_report,
    write_filtering_summary_csv,
)


MISSING_CHOICES = ("keep", "remove")
EMPTY_EXPRESSION_CHOICES = ("match_all", "match_none")


def _is_missing(value) -> bool:
    """Return whether a diagnostic count is unavailable or NaN-like."""
    if value is None:
        return True
    try:
        return value != value
    except Exception:
        return True


def _fmt_count(value, missing: str = "unavailable") -> str:
    """Thousands-separated count, or a word - never ``<NA>`` and never a fake 0."""
    if _is_missing(value):
        return missing
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def _filtering_summary_lines(
    missing_checks: List[Dict[str, Any]],
    condition_checks: List[Dict[str, Any]],
    variants_before: Optional[int],
    variants_after: Optional[int],
    terminal_label_width: int,
    reason_statistics: Optional[Dict[str, Any]] = None,
    reason_report: Optional[str] = None,
    data_flow: Optional[Dict[str, Any]] = None,
    soft_vcf: Optional[str] = None,
    soft_variants: Optional[int] = None,
    dataset_id: str = "Not provided",
    genome_build: str = "Not provided",
    summary_csv: Optional[str] = None,
    html_report: Optional[str] = None,
    input_validation: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """Compatibility wrapper around the canonical filtering summary renderer."""
    statistics = dict(reason_statistics or {})
    removed = (
        None
        if variants_before is None or variants_after is None
        else variants_before - variants_after
    )
    if "reconciled" not in statistics:
        statistics["reconciled"] = (
            removed is not None
            and statistics.get("primary_removed_total") is not None
            and int(statistics["primary_removed_total"]) == int(removed)
        )
    summary = build_filtering_summary(
        dataset_id=dataset_id,
        genome_build=genome_build,
        missing_checks=missing_checks,
        condition_checks=condition_checks,
        reason_statistics=statistics,
        variants_before=variants_before,
        variants_after=variants_after,
        soft_filter_variants=soft_variants,
        soft_filter_enabled=soft_vcf is not None,
        input_validation=input_validation,
        data_flow=data_flow,
        outputs={
            "soft_filtered_vcf": soft_vcf,
            "reason_summary": reason_report,
            "summary_csv": summary_csv,
            "html_report": html_report,
        },
        resolved_configuration={},
        runtime_seconds=0.0,
        display_missing_counts=bool(missing_checks),
    )
    return render_filtering_summary(
        summary, label_width=terminal_label_width,
    ).splitlines()


def _summarize_filter_tags_in_process(
    tag_file: Path,
    checks: List[Dict[str, Any]],
    expected_threads: int | None = None,
) -> Dict[str, Any]:
    """Use one lazy Polars aggregation to count raw, overlapping and primary reasons."""
    aggregation = (
        current_polars_thread_runtime(expected_threads)
        if expected_threads is not None else {
            "requested_threads": None,
            "polars_thread_pool_size": int(pl.thread_pool_size()),
            "thread_budget_enforced": False,
            "thread_environment_variable": None,
        }
    )
    aggregation["polars_version"] = getattr(pl, "__version__", "unknown")
    if not checks or tag_file.stat().st_size == 0:
        return {
            "checks": [
                dict(check, count=0, primary_count=0) for check in checks
            ],
            "primary_removed_total": 0,
            "overlap_variants": 0,
            "extra_rule_matches": 0,
            "aggregation": aggregation,
        }
    scan = pl.scan_csv(
        tag_file,
        has_header=False,
        separator="\t",
        schema={"filter_tags": pl.String},
    )
    tags = pl.col("filter_tags").fill_null("")
    masks = {
        check["tag"]: tags.str.contains(
            ";%s;" % check["tag"],
            literal=True,
        )
        for check in checks
    }
    selections = []
    for check in checks:
        selections.append(masks[check["tag"]].sum().alias(check["tag"] + "__observed"))

    removal_checks = [check for check in checks if check.get("removes")]
    seen = pl.lit(False)
    for check in removal_checks:
        mask = masks[check["tag"]]
        selections.append((mask & ~seen).sum().alias(check["tag"] + "__primary"))
        seen = seen | mask

    if removal_checks:
        failure_count = pl.sum_horizontal([
            masks[check["tag"]].cast(pl.UInt16) for check in removal_checks
        ])
        selections.extend([
            seen.sum().alias("__primary_removed_total"),
            (failure_count > 1).sum().alias("__overlap_variants"),
            (failure_count.cast(pl.Int32) - 1)
            .clip(lower_bound=0)
            .sum()
            .alias("__extra_rule_matches"),
        ])
    else:
        selections.extend([
            pl.lit(0).alias("__primary_removed_total"),
            pl.lit(0).alias("__overlap_variants"),
            pl.lit(0).alias("__extra_rule_matches"),
        ])

    values = collect_streaming(scan.select(selections)).to_dicts()[0]
    rows = []
    for check in checks:
        row = dict(check)
        row["count"] = int(values.get(check["tag"] + "__observed") or 0)
        row["primary_count"] = (
            int(values.get(check["tag"] + "__primary") or 0)
            if check.get("removes") else 0
        )
        rows.append(row)
    return {
        "checks": rows,
        "primary_removed_total": int(values.get("__primary_removed_total") or 0),
        "overlap_variants": int(values.get("__overlap_variants") or 0),
        "extra_rule_matches": int(values.get("__extra_rule_matches") or 0),
        "aggregation": aggregation,
    }


def _summarize_filter_tags(
    tag_file: Path,
    checks: List[Dict[str, Any]],
    *,
    threads: int | None = None,
) -> Dict[str, Any]:
    """Aggregate in-process for unit use or in a pre-import bounded worker."""
    if threads is None:
        return _summarize_filter_tags_in_process(tag_file, checks)
    return run_in_bounded_polars_process(
        _summarize_filter_tags_in_process,
        tag_file,
        checks,
        threads=threads,
        call_kwargs={"expected_threads": threads},
    )


def _collect_filter_reason_statistics(
    vcf_path: str,
    checks: List[Dict[str, Any]],
    output_folder: str,
    output_prefix: str,
    log_warn,
    *,
    bcftools_bin: str,
    bash_bin: str,
    tabix_bin: str,
    threads: int,
    max_mem: str,
    sort_output: bool,
    soft_vcf_path: Path | None,
) -> Optional[Dict[str, Any]]:
    """Tag removal reasons, optionally retain that stream, and aggregate counts."""
    if not checks and soft_vcf_path is None:
        return {
            "checks": [], "primary_removed_total": 0,
            "overlap_variants": 0, "extra_rule_matches": 0,
            "aggregation": {
                "requested_threads": int(threads),
                "polars_thread_pool_size": None,
                "thread_budget_enforced": None,
                "thread_environment_variable": None,
                "status": "not_required_no_active_reason_checks",
            },
        }

    handle = tempfile.NamedTemporaryFile(
        prefix=".%s_filter_reason_tags_" % output_prefix,
        suffix=".tsv",
        dir=output_folder,
        delete=False,
    )
    tag_file = Path(handle.name)
    handle.close()
    quoted_bcftools = shlex.quote(bcftools_bin)
    quoted_tabix = shlex.quote(tabix_bin)
    stages = [
        "%s view -Ou %s" % (quoted_bcftools, shlex.quote(str(vcf_path)))
    ]
    for check in checks:
        stages.append(
            "%s filter -Ou -m + -s %s -e %s"
            % (quoted_bcftools, check["tag"], shlex.quote(check["expr"]))
        )
    if soft_vcf_path is None:
        stages.append("%s query -f ';%%FILTER;\\n'" % quoted_bcftools)
        inner = "set -euo pipefail; %s > %s" % (
            " | ".join(stages), shlex.quote(str(tag_file)),
        )
    else:
        if sort_output:
            stages.append(
                "%s sort --temp-dir %s --max-mem %s"
                % (
                    quoted_bcftools,
                    shlex.quote(str(soft_vcf_path.parent)),
                    shlex.quote(str(max_mem)),
                )
            )
        stages.append("%s view -Oz --threads %d" % (quoted_bcftools, threads))
        quoted_soft_vcf = shlex.quote(str(soft_vcf_path))
        soft_pipeline = "%s > %s" % (" | ".join(stages), quoted_soft_vcf)
        query = "%s query -f ';%%FILTER;\\n' %s > %s" % (
            quoted_bcftools,
            quoted_soft_vcf,
            shlex.quote(str(tag_file)),
        )
        inner = "set -euo pipefail; %s && %s -f -p vcf %s" % (
            soft_pipeline,
            quoted_tabix,
            quoted_soft_vcf,
        )
        if checks:
            inner += " && %s" % query
    try:
        result = run_cmd(
            "%s -c %s" % (shlex.quote(bash_bin), shlex.quote(inner))
        )
        if result.returncode != 0:
            log_warn(
                "Filter-reason tagging failed (exit %s): %s"
                % (result.returncode, (result.stderr or "").strip())
            )
            return None
        if soft_vcf_path is not None:
            soft_index = Path(str(soft_vcf_path) + ".tbi")
            if (
                not soft_vcf_path.is_file()
                or soft_vcf_path.stat().st_size <= 0
                or not soft_index.is_file()
                or soft_index.stat().st_size <= 0
            ):
                log_warn(
                    "Soft-filter tagging did not create a non-empty VCF and index."
                )
                return None
        if not checks:
            return {
                "checks": [], "primary_removed_total": 0,
                "overlap_variants": 0, "extra_rule_matches": 0,
                "aggregation": {
                    "requested_threads": int(threads),
                    "polars_thread_pool_size": None,
                    "thread_budget_enforced": None,
                    "thread_environment_variable": None,
                    "status": "not_required_no_active_reason_checks",
                },
            }
        return _summarize_filter_tags(tag_file, checks, threads=threads)
    except Exception as exc:
        log_warn(
            "Filter-reason tagging failed: %s: %s" % (type(exc).__name__, exc)
        )
        return None
    finally:
        try:
            tag_file.unlink()
        except OSError:
            pass


def _write_filter_reason_report(
    statistics: Optional[Dict[str, Any]],
    report_path: str | Path,
) -> Optional[str]:
    """Persist the exact reason accounting as a reloadable TSV."""
    if statistics is None:
        return None
    path = Path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        encoding="utf-8",
        prefix=".%s." % path.name,
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(handle.name)
    try:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "order", "filter_id", "category", "reason", "action",
                "variants_matching_reason", "variants_removed_for_this_reason",
                "bcftools_expression", "status",
            ],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for order, check in enumerate(statistics.get("checks") or [], 1):
            writer.writerow({
                "order": order,
                "filter_id": (
                    check.get("tag")
                    if check.get("removes", check.get("action") == "remove")
                    else None
                ),
                "category": check.get("category"),
                "reason": check.get("label"),
                "action": check.get("action", "remove"),
                "variants_matching_reason": check.get("count"),
                "variants_removed_for_this_reason": check.get("primary_count"),
                "bcftools_expression": check.get("expr"),
            })
        for reason, observed, primary, status in (
            ("Variants failing multiple active rules",
             statistics.get("overlap_variants"), None, None),
            ("Rule failures beyond the first",
             statistics.get("extra_rule_matches"), None, None),
            ("All primary removal assignments",
             None, statistics.get("primary_removed_total"), None),
            ("Removed from final VCF",
             None, statistics.get("actual_removed"), None),
            ("Reason reconciliation",
             None, None, "PASS" if statistics.get("reconciled") else "FAILED"),
        ):
            writer.writerow({
                "category": "summary",
                "reason": reason,
                "variants_matching_reason": observed,
                "variants_removed_for_this_reason": primary,
                "status": status,
            })
        handle.flush()
        os.fsync(handle.fileno())
        handle.close()
        os.replace(temporary_path, path)
    except BaseException:
        handle.close()
        temporary_path.unlink(missing_ok=True)
        raise
    return str(path)


def _filtering_field_provenance(
    *,
    header: str,
    provenance_headers: Mapping[str, str],
    provenance_missing_value: str,
    pvalue_field: str | None,
    study_af_field: str | None,
    study_info_af_field: str | None,
    external_info_af_field: str | None,
    imputation_field: str | None,
    af_uses: List[str],
) -> tuple[Dict[str, str], List[Dict[str, str]]]:
    """Summarize configured harmonisation provenance for active VCF fields.

    Detailed field provenance may be incomplete even after the separate
    PostGWAS-origin contract passes. Actual header names remain
    configuration-driven.
    """
    logical_names = ["postgwas_version", "vcf_status"]
    if pvalue_field is not None:
        logical_names.append("pvalue_harmonisation")
    if study_af_field is not None:
        logical_names.append("af_source")
    if study_info_af_field is not None:
        logical_names.append("population_af_fields")
    if imputation_field is not None:
        logical_names.append("info_source")
    missing_mapping = [
        name for name in logical_names if name not in provenance_headers
    ]
    if missing_mapping:
        raise ValueError(
            "Configured PostGWAS provenance mapping is missing filtering "
            "field(s): %s" % ", ".join(missing_mapping)
        )
    selected_headers = {
        name: str(provenance_headers[name]) for name in logical_names
    }
    metadata = declared_vcf_metadata_values(
        header, selected_headers.values(),
    )
    values = {
        name: metadata.get(header_name)
        for name, header_name in selected_headers.items()
    }

    def available(value: str | None) -> bool:
        return bool(value and value != provenance_missing_value)

    available_count = sum(available(value) for value in values.values())
    if available_count == len(values):
        provenance_status = "AVAILABLE"
        detail = "PostGWAS %s; VCF status %s" % (
            values["postgwas_version"], values["vcf_status"],
        )
    elif available_count:
        provenance_status = "PARTIAL"
        detail = "%d/%d relevant PostGWAS metadata values available" % (
            available_count, len(values),
        )
    else:
        provenance_status = "UNAVAILABLE"
        detail = "PostGWAS harmonisation provenance headers not declared"

    summaries: List[Dict[str, str]] = []
    if pvalue_field is not None:
        pvalue_source = values.get("pvalue_harmonisation")
        summaries.append({
            "kind": "analysis",
            "label": "P-value evidence",
            "value": "%s · %s" % (
                pvalue_field,
                pvalue_source if available(pvalue_source) else "source unavailable",
            ),
        })
    af_fields = []
    if study_af_field is not None:
        af_fields.append(
            "%s for %s" % (study_af_field, " and ".join(af_uses))
        )
    if study_info_af_field is not None and external_info_af_field is not None:
        af_fields.append(
            "%s versus %s for external AF concordance"
            % (study_info_af_field, external_info_af_field)
        )
    if af_fields:
        summaries.append({
            "kind": "genetic",
            "label": "Allele-frequency fields",
            "value": "; ".join(af_fields),
        })
    if imputation_field is not None:
        info_source = values.get("info_source")
        summaries.append({
            "kind": "analysis",
            "label": "Imputation quality score",
            "value": "%s · %s" % (
                imputation_field,
                info_source if available(info_source) else "source unavailable",
            ),
        })
    return {
        "status": provenance_status,
        "detail": detail,
    }, summaries


def filter_gwas_vcf_bcftools(
    vcf_path: str,
    output_folder: str,
    output_prefix: str,
    *,
    genome_build_header: str,
    supported_genome_builds: List[str],
    required_provenance_headers: Mapping[str, str],
    provenance_headers: Mapping[str, str],
    provenance_missing_value: str,
    pval_cutoff: Optional[float],
    maf_cutoff: Optional[float],
    allelefreq_diff_cutoff: Optional[float],
    info_min: Optional[float],
    info_max: Optional[float],
    info_missing: str,
    external_af_name: str,
    include_indels: bool,
    exclude_palindromic: bool,
    palindromic_af_lower: float,
    palindromic_af_upper: float,
    remove_mhc: bool,
    mhc_regions: Mapping[str, Mapping[str, Any]],
    mhc_region_override: Mapping[str, Any],
    threads: int,
    max_mem: str,
    lp_missing: str,
    af_missing: str,
    empty_expression: str,
    sort_output: bool,
    write_soft_filter_vcf: bool,
    filter_reason_ids: Mapping[str, str],
    vcf_fields: Mapping[str, str],
    output_layout: Mapping[str, str],
    bcftools_bin: str,
    tabix_bin: str,
    bash_bin: str,
    resolved_configuration: Mapping[str, Any],
    report_missing_counts: bool,
    terminal_label_width: int,
    stage_progress: StageProgress | None = None,
) -> Dict[str, Any]:
    """Filter a GWAS-VCF and report, honestly, how many variants each rule cost.

    All scientific policies, VCF fields, output paths, and executable paths are
    supplied from the schema-validated canonical filtering configuration.
    """
    progress = stage_progress or StageProgress(
        "Summary-statistics filtering progress", enabled=False,
    )
    total_stages = 5
    if info_missing not in MISSING_CHOICES:
        raise ValueError("❌ info_missing must be either 'keep' or 'remove'")
    if af_missing not in MISSING_CHOICES:
        raise ValueError("❌ af_missing must be either 'keep' or 'remove'")
    if lp_missing not in MISSING_CHOICES:
        raise ValueError("❌ lp_missing must be either 'keep' or 'remove'")
    if empty_expression not in EMPTY_EXPRESSION_CHOICES:
        raise ValueError("❌ empty_expression must be either 'match_all' or 'match_none'")
    if threads < 1:
        raise ValueError("❌ threads must be at least 1")
    output_prefix = validate_filename_component(output_prefix, "dataset_id")
    if not VCF_TAG.fullmatch(external_af_name):
        raise ValueError("❌ external_af_name must be a valid VCF tag")
    configured_reason_ids = list(filter_reason_ids.values())
    if (
        not configured_reason_ids
        or any(
            not VCF_TAG.fullmatch(filter_id) or filter_id == "PASS"
            for filter_id in configured_reason_ids
        )
        or len(configured_reason_ids) != len(set(configured_reason_ids))
    ):
        raise ValueError("❌ Configured filtering reason IDs are invalid or duplicated")

    started_at = time.monotonic()
    fixed_pattern = re.compile(r"^[A-Z][A-Z0-9_]*$")
    format_pattern = re.compile(r"^FORMAT/[A-Za-z][A-Za-z0-9_.-]*$")
    info_pattern = re.compile(r"^INFO/[A-Za-z][A-Za-z0-9_.-]*$")
    chromosome_field = str(vcf_fields["chromosome"])
    position_field = str(vcf_fields["position"])
    reference_field = str(vcf_fields["reference_allele"])
    alternate_field = str(vcf_fields["alternate_allele"])
    variant_type_field = str(vcf_fields["variant_type"])
    study_af_format = str(vcf_fields["study_af_format"])
    imputation_quality_format = str(vcf_fields["imputation_quality_format"])
    log_pvalue_format = str(vcf_fields["log_pvalue_format"])
    study_af_info = str(vcf_fields["study_af_info"])
    external_af_info = str(vcf_fields["external_af_info"]).format(
        reference_population_tag=external_af_name
    )
    rule_display_groups = variant_qc_rule_display_groups({
        "chromosome": chromosome_field,
        "position": position_field,
        "reference_allele": reference_field,
        "alternate_allele": alternate_field,
        "variant_type": variant_type_field,
        "study_af": study_af_format,
        "imputation_quality": imputation_quality_format,
        "log_pvalue": log_pvalue_format,
        "study_info_af": study_af_info,
        "external_info_af": external_af_info,
    })
    if any(
        not fixed_pattern.fullmatch(field)
        for field in (
            chromosome_field,
            position_field,
            reference_field,
            alternate_field,
            variant_type_field,
        )
    ) or any(
        not format_pattern.fullmatch(field)
        for field in (study_af_format, imputation_quality_format, log_pvalue_format)
    ) or any(
        not info_pattern.fullmatch(field)
        for field in (study_af_info, external_af_info)
    ):
        raise ValueError("❌ Configured filtering VCF fields are invalid")

    step_dir = Path(output_folder).expanduser().resolve()
    step_dir.mkdir(parents=True, exist_ok=True)
    log_file = configured_output_path(
        step_dir, output_layout["preflight_log"], dataset_id=output_prefix
    )
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_buffer = io.StringIO()

    def log_print(*args):
        text = " ".join(str(a) for a in args)
        log_buffer.write(text + "\n")

    def log_warn(*args):
        text = " ".join(str(a) for a in args)
        log_buffer.write(text + "\n")

    def write_log():
        handle = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".%s." % log_file.name,
            suffix=".tmp",
            dir=log_file.parent,
            delete=False,
        )
        temporary_log = Path(handle.name)
        try:
            handle.write(log_buffer.getvalue())
            handle.flush()
            os.fsync(handle.fileno())
            handle.close()
            os.replace(temporary_log, log_file)
        except BaseException:
            handle.close()
            temporary_log.unlink(missing_ok=True)
            raise

    quoted_bcftools = shlex.quote(bcftools_bin)
    quoted_bash = shlex.quote(bash_bin)

    def fast_count(vcf: str) -> Optional[int]:
        """Number of VCF records, using an index when one is available."""
        if current_validation_session() is not None:
            # Pipeline entry and generated VCFs require validated indexes.
            # Do not bypass an identity/index failure with the direct-mode
            # unindexed-input fallback below.
            return count_indexed_vcf_records(vcf, bcftools_bin)
        quoted = shlex.quote(str(vcf))
        problems = []
        try:
            res = run_cmd(
                "%s index -n %s" % (quoted_bcftools, quoted), check=False,
            )
            if res.returncode == 0:
                return int(res.stdout.strip())
            problems.append(
                "`bcftools index -n` exit %s: %s"
                % (res.returncode, res.stderr.strip())
            )
        except Exception as exc:
            problems.append(
                "`bcftools index -n` raised %s: %s"
                % (type(exc).__name__, exc)
            )

        fallback = (
            "set -o pipefail; %s view --threads %d -H %s | wc -l"
            % (quoted_bcftools, threads, quoted)
        )
        try:
            res = run_cmd(
                "%s -c %s" % (quoted_bash, shlex.quote(fallback)),
                check=False,
            )
            if res.returncode == 0:
                return int(res.stdout.strip())
            problems.append(
                "`bcftools view -H | wc -l` exit %s: %s"
                % (res.returncode, res.stderr.strip())
            )
        except Exception as exc:
            problems.append(
                "`bcftools view -H | wc -l` raised %s: %s"
                % (type(exc).__name__, exc)
            )

        log_warn(
            "⚠️ Could not count the variants in %s, so every count derived from "
            "it is reported as 'unavailable' rather than 0. Reasons: %s"
            % (vcf, "; ".join(problems))
        )
        return None

    # ============================================================
    # VALIDATION
    # ============================================================
    log_print(
        "Resolved configuration:\n%s"
        % json.dumps(resolved_configuration, indent=2, sort_keys=True)
    )
    with progress.step(1, total_stages, "Validate input VCF and count variants"):
        input_vcf = Path(vcf_path).expanduser().resolve()
        if not input_vcf.is_file() or input_vcf.stat().st_size <= 0:
            msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
            log_print(msg)
            log_print("❌ STATUS: FAILED")
            write_log()
            raise FileNotFoundError(msg)
        genome_build = None
        contigs = None
        pre_variants = None
        required_field_checks: List[Dict[str, Any]] = []
        required_fields_evaluated = False
        field_provenance = {
            "status": "UNAVAILABLE",
            "detail": "PostGWAS harmonisation provenance was not evaluated",
        }
        field_summaries: List[Dict[str, str]] = []
        try:
            header = read_vcf_header(input_vcf, bcftools_bin, error_type=ValueError)
            origin = validate_postgwas_vcf_provenance(
                header,
                required_provenance_headers,
                vcf_path=input_vcf,
                error_type=ValueError,
            )
            if current_validation_recorder() is None:
                log_print("PostGWAS VCF provenance: %s" % json.dumps(origin, sort_keys=True))
            required_field_reasons: Dict[str, List[str]] = {}

            def require_field(field: str, reason: str) -> None:
                reasons = required_field_reasons.setdefault(field, [])
                if reason not in reasons:
                    reasons.append(reason)

            if pval_cutoff is not None:
                require_field(log_pvalue_format, "minimum -log10(P) filter")
            if maf_cutoff is not None:
                require_field(study_af_format, "minor-allele-frequency filter")
            if exclude_palindromic:
                require_field(
                    study_af_format, "palindromic-ambiguity filter",
                )
            if info_min is not None or info_max is not None:
                require_field(
                    imputation_quality_format, "imputation-quality filter",
                )
            if allelefreq_diff_cutoff is not None:
                require_field(
                    study_af_info, "study/reference AF-difference filter",
                )
                require_field(
                    external_af_info, "study/reference AF-difference filter",
                )
            required_fields = list(required_field_reasons)
            field_presence = required_vcf_field_presence(
                header, required_fields,
            )
            required_field_checks = [
                {
                    "field": field,
                    "declared": field_presence[field],
                    "required_by": required_field_reasons[field],
                }
                for field in required_fields
            ]
            required_fields_evaluated = True
            af_uses = []
            if maf_cutoff is not None:
                af_uses.append("MAF")
            if exclude_palindromic:
                af_uses.append("palindromic ambiguity")
            field_provenance, field_summaries = _filtering_field_provenance(
                header=header,
                provenance_headers=provenance_headers,
                provenance_missing_value=provenance_missing_value,
                pvalue_field=(
                    log_pvalue_format if pval_cutoff is not None else None
                ),
                study_af_field=study_af_format if af_uses else None,
                study_info_af_field=(
                    study_af_info
                    if allelefreq_diff_cutoff is not None else None
                ),
                external_info_af_field=(
                    external_af_info
                    if allelefreq_diff_cutoff is not None else None
                ),
                imputation_field=(
                    imputation_quality_format
                    if info_min is not None or info_max is not None else None
                ),
                af_uses=af_uses,
            )
            genome_build, contigs = validate_vcf_header_contract(
                header=header,
                genome_build_header=genome_build_header,
                supported_genome_builds=supported_genome_builds,
                required_fields=required_fields,
            )
            colliding_filter_ids = sorted(
                declared_vcf_tags(header, "FILTER") & set(configured_reason_ids)
            )
            if colliding_filter_ids:
                raise ValueError(
                    "Input VCF already declares configured PostGWAS FILTER IDs: %s. "
                    "Use an unfiltered harmonisation VCF or configure non-colliding "
                    "modules.filtering.filter_reason_ids values."
                    % ", ".join(colliding_filter_ids)
                )
            mhc_chrom = None
            configured_mhc_chrom = None
            mhc_start = None
            mhc_end = None
            if remove_mhc:
                mhc_region = mhc_regions.get(genome_build)
                if mhc_region is None:
                    raise ValueError(
                        "MHC removal is enabled, but modules.filtering.mhc_regions "
                        "does not define %s" % genome_build
                    )
                try:
                    configured_mhc_chrom = str(
                        mhc_region_override.get("chromosome")
                        or mhc_region["chromosome"]
                    )
                    mhc_chrom = configured_mhc_chrom
                    mhc_start = int(
                        mhc_region["start"]
                        if mhc_region_override.get("start") is None
                        else mhc_region_override["start"]
                    )
                    mhc_end = int(
                        mhc_region["end"]
                        if mhc_region_override.get("end") is None
                        else mhc_region_override["end"]
                    )
                except (KeyError, TypeError, ValueError) as exc:
                    raise ValueError(
                        "The configured MHC region for %s is incomplete or invalid"
                        % genome_build
                    ) from exc
                if mhc_end <= mhc_start:
                    raise ValueError(
                        "The effective MHC end must be greater than its start for %s"
                        % genome_build
                    )
                mhc_chrom = _match_contig_naming(
                    str(mhc_chrom), contigs, log_warn,
                )
            policy = VariantQCPolicy(
                minimum_neglog10_p=pval_cutoff,
                missing_pvalue_action=lp_missing,
                maf_min=maf_cutoff,
                missing_af_action=af_missing,
                info_min=info_min,
                info_max=info_max,
                missing_info_action=info_missing,
                maximum_af_difference=allelefreq_diff_cutoff,
                include_indels=include_indels,
                remove_palindromic=exclude_palindromic,
                palindromic_lower=palindromic_af_lower,
                palindromic_upper=palindromic_af_upper,
                remove_mhc=remove_mhc,
                mhc_chromosome=configured_mhc_chrom,
                mhc_start=mhc_start,
                mhc_end=mhc_end,
            )
            pre_variants = fast_count(str(input_vcf))
            if pre_variants is None:
                raise RuntimeError("Input VCF variant counting failed")
            input_validation = {
                "status": "PASS",
                "input_vcf": str(input_vcf),
                "genome_build": genome_build,
                "declared_contigs": len(contigs),
                "total_variants": pre_variants,
                "required_fields": required_field_checks,
                "required_fields_evaluated": True,
                "field_provenance": field_provenance,
                "field_summaries": field_summaries,
                "error": None,
            }
            record_file_validation(
                input_vcf, "Filtering input VCF",
                checks=("configured header contract", "non-colliding filter identifiers", "record count"),
                metrics={"declared_genome_build": genome_build, "contigs": list(contigs),
                         "total_variants": pre_variants, "required_fields": required_field_checks},
            )
        except BaseException as exc:
            log_print(
                "❌ Input VCF preflight failed: %s: %s"
                % (type(exc).__name__, exc)
            )
            failed_validation = {
                "status": "FAILED",
                "input_vcf": str(input_vcf),
                "genome_build": genome_build,
                "declared_contigs": (
                    len(contigs) if contigs is not None else None
                ),
                "total_variants": pre_variants,
                "required_fields": required_field_checks,
                "required_fields_evaluated": required_fields_evaluated,
                "field_provenance": field_provenance,
                "field_summaries": field_summaries,
                "error": "%s: %s" % (type(exc).__name__, exc),
            }
            validation_text = render_input_vcf_validation(
                failed_validation, label_width=terminal_label_width,
            )
            for line in validation_text.splitlines():
                log_print(line)
            log_print("❌ STATUS: FAILED")
            write_log()
            progress.print_block(style_screen_block(validation_text))
            raise
    validation_text = render_input_vcf_validation(
        input_validation, label_width=terminal_label_width,
    )
    for line in validation_text.splitlines():
        log_print(line)
    progress.print_block(style_screen_block(validation_text))
    if current_validation_recorder() is not None:
        log_print("VCF preflight passed; file details are combined in the validation audit.")
    else:
        log_print(
            "✅ VCF preflight passed: genome_build=%s; variants=%s; active_fields=%s"
            % (
                genome_build,
                pre_variants,
                ", ".join(required_fields) or "none",
            )
        )
    if remove_mhc:
        override_active = any(
            mhc_region_override.get(key) is not None
            for key in ("chromosome", "start", "end")
        )
        log_print(
            "✅ MHC region selected: %s:%d-%d (%s; %s)"
            % (
                mhc_chrom,
                mhc_start,
                mhc_end,
                genome_build,
                "custom override" if override_active else "build-specific default",
            )
        )
    resolved_configuration = dict(resolved_configuration)
    resolved_configuration["variant_qc_policy"] = policy.as_dict()
    log_print(
        "Shared variant-QC policy:\n%s"
        % json.dumps(policy.as_dict(), indent=2, sort_keys=True)
    )
    path_values = {"dataset_id": output_prefix, "genome_build": genome_build}
    log_file = configured_output_path(
        step_dir, output_layout["log_file"], **path_values
    )
    output_vcf = configured_output_path(
        step_dir, output_layout["filtered_vcf"], **path_values
    )
    soft_output_vcf = configured_output_path(
        step_dir, output_layout["soft_filtered_vcf"], **path_values
    )
    reason_report_path = configured_output_path(
        step_dir, output_layout["reason_summary"], **path_values
    )
    summary_csv_path = configured_output_path(
        step_dir, output_layout["summary_csv"], **path_values
    )
    html_report_path = configured_output_path(
        step_dir, output_layout["html_report"], **path_values
    )
    mhc_bed = configured_output_path(
        step_dir, output_layout["mhc_exclusion_bed"], **path_values
    )
    for path in (
        log_file,
        output_vcf,
        soft_output_vcf,
        reason_report_path,
        summary_csv_path,
        html_report_path,
        mhc_bed,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
    for label, executable in (
        ("bcftools", bcftools_bin),
        ("tabix", tabix_bin),
        ("bash", bash_bin),
    ):
        version = run_cmd([executable, "--version"], check=False, shell=False)
        version_text = (version.stdout or version.stderr or "").splitlines()
        log_print(
            "🔧 %s version: %s"
            % (label, version_text[0].strip() if version_text else "unavailable")
        )
    log_print("🔧 Polars version: %s" % getattr(pl, "__version__", "unavailable"))
    quoted_input = shlex.quote(str(input_vcf))
    quoted_tabix = shlex.quote(tabix_bin)

    def count_matching(vcf: str, expr: str) -> Optional[int]:
        """Count records matching a bcftools expression.

        ``set -o pipefail`` is what makes this trustworthy: piping into ``wc -l``
        otherwise hides a bcftools abort behind wc's exit status 0.
        """
        # `query -f '.\n'` counts one line per record without formatting the
        # record itself, which is markedly cheaper than `view -H` on a big file.
        inner = "set -o pipefail; %s query -i %s -f '.\\n' %s | wc -l" % (
            quoted_bcftools, shlex.quote(expr), shlex.quote(str(vcf)))
        try:
            res = run_cmd(
                "%s -c %s" % (quoted_bash, shlex.quote(inner)),
                check=False,
            )
        except Exception as exc:
            log_warn(f"⚠️ Could not count `{expr}`: {type(exc).__name__}: {exc}")
            return None
        if res.returncode != 0:
            log_warn(
                "⚠️ Could not count `%s` (exit %s): %s"
                % (expr, res.returncode, (res.stderr or "").strip())
            )
            return None
        try:
            return int(res.stdout.strip())
        except (TypeError, ValueError):
            log_warn(f"⚠️ Unexpected output counting `{expr}`: {res.stdout!r}")
            return None

    # ============================================================
    # P-VALUE / LP THRESHOLD
    # The configured log-P field is -log10(p): bigger means more significant.
    # ============================================================
    resolved_lp = float(pval_cutoff) if pval_cutoff is not None else None

    # ============================================================
    # INCLUDE LOGIC
    # bcftools treats a MISSING value as failing an -i test, so every
    # `*_missing` option below is the difference between a documented removal
    # and a silent one.
    # ============================================================
    include_parts: List[str] = []
    missing_checks: List[Dict[str, Any]] = []
    condition_checks: List[Dict[str, Any]] = []

    if resolved_lp is not None:
        lp_expr = f"({log_pvalue_format} >= {resolved_lp})"
        if lp_missing == "keep":
            lp_expr = f"({lp_expr} | ({log_pvalue_format} == '.'))"
        missing_checks.append({
            "field": log_pvalue_format,
            "expr": f"{log_pvalue_format} == '.'",
            "action": lp_missing, "key": "lp_missing", "priority": 10,
            "reason_key": "missing_pvalue",
            **rule_display_groups["significance"],
        })
        include_parts.append(lp_expr)
        condition_checks.append({
            "label": "P-value evidence below LP %.6g" % resolved_lp,
            "expr": (
                f"({log_pvalue_format} != '.' & "
                f"{log_pvalue_format} < {resolved_lp})"
            ),
            "priority": 11,
            "reason_key": "pvalue_below_threshold",
            "policy_key": "significance",
            **rule_display_groups["significance"],
        })

    if maf_cutoff is not None:
        maf_expr = (
            f"({study_af_format} >= {maf_cutoff} & "
            f"{study_af_format} <= {1 - maf_cutoff})"
        )
        if af_missing == "keep":
            maf_expr = f"({maf_expr} | ({study_af_format} == '.'))"
        missing_checks.append({
            "field": study_af_format,
            "expr": f"{study_af_format} == '.'",
            "action": af_missing, "key": "af_missing", "priority": 20,
            "reason_key": "missing_study_af",
            **rule_display_groups["allele_frequency"],
        })
        include_parts.append(maf_expr)
        condition_checks.append({
            "label": "Minor allele frequency outside %.6g ≤ AF ≤ %.6g"
            % (maf_cutoff, 1 - maf_cutoff),
            "expr": (
                f"({study_af_format} != '.' & "
                f"({study_af_format} < {maf_cutoff} | "
                f"{study_af_format} > {1 - maf_cutoff}))"
            ),
            "priority": 21,
            "reason_key": "maf_outside_range",
            "policy_key": "allele_frequency",
            **rule_display_groups["allele_frequency"],
        })

    info_expr = None
    base_expr = None
    info_failure_expr = None
    if info_min is not None or info_max is not None:
        if info_min is not None and info_max is not None:
            if info_min > info_max:
                raise ValueError("❌ Invalid INFO range: info_min > info_max")
            base_expr = (
                f"({imputation_quality_format} >= {info_min} & "
                f"{imputation_quality_format} <= {info_max})"
            )
            info_failure_expr = (
                f"({imputation_quality_format} != '.' & "
                f"({imputation_quality_format} < {info_min} | "
                f"{imputation_quality_format} > {info_max}))"
            )
        elif info_min is not None:
            base_expr = f"({imputation_quality_format} >= {info_min})"
            info_failure_expr = (
                f"({imputation_quality_format} != '.' & "
                f"{imputation_quality_format} < {info_min})"
            )
        else:
            base_expr = f"({imputation_quality_format} <= {info_max})"
            info_failure_expr = (
                f"({imputation_quality_format} != '.' & "
                f"{imputation_quality_format} > {info_max})"
            )

    if base_expr is not None:
        if info_missing == "keep":
            # Inside the filter, '.' must be in single quotes
            info_expr = f"({base_expr} | ({imputation_quality_format} == '.'))"
        else:
            info_expr = f"({base_expr} & ({imputation_quality_format} != '.'))"
        missing_checks.append({
            "field": imputation_quality_format,
            "expr": f"{imputation_quality_format} == '.'",
            "action": info_missing, "key": "info_missing", "priority": 30,
            "reason_key": "missing_imputation_quality",
            **rule_display_groups["imputation_quality"],
        })

    if info_expr is not None:
        include_parts.append(info_expr)
        if info_min is not None and info_max is not None:
            info_label = "Imputation quality outside %.6g ≤ INFO (SI) ≤ %.6g" % (
                info_min, info_max,
            )
        elif info_min is not None:
            info_label = "Imputation quality below INFO (SI) %.6g" % info_min
        else:
            info_label = "Imputation quality above INFO (SI) %.6g" % info_max
        condition_checks.append({
            "label": info_label,
            "expr": info_failure_expr,
            "priority": 31,
            "reason_key": "imputation_quality_outside_range",
            "policy_key": "imputation_quality",
            **rule_display_groups["imputation_quality"],
        })

    if allelefreq_diff_cutoff is not None:
        af_diff_expr = (
            f"(abs({study_af_info} - {external_af_info}) "
            f"<= {allelefreq_diff_cutoff})"
        )
        if af_missing == "keep":
            af_diff_expr = (
                f"({af_diff_expr} | ({study_af_info} == '.') | "
                f"({external_af_info} == '.'))"
            )
        missing_checks.extend([{
            "field": study_af_info,
            "expr": f"{study_af_info} == '.'",
            "action": af_missing,
            "key": "af_missing",
            "priority": 40,
            "reason_key": "missing_study_info_af",
            **rule_display_groups["external_af_concordance"],
        }, {
            "field": external_af_info,
            "expr": f"{external_af_info} == '.'",
            "action": af_missing,
            "key": "af_missing",
            "priority": 41,
            "reason_key": "missing_external_af",
            **rule_display_groups["external_af_concordance"],
        }])
        include_parts.append(af_diff_expr)
        condition_checks.append({
            "label": "External AF absolute difference |AF − %s| > %.6g"
            % (external_af_name, allelefreq_diff_cutoff),
            "expr": (
                f"({study_af_info} != '.' & {external_af_info} != '.' & "
                f"abs({study_af_info} - {external_af_info}) "
                f"> {allelefreq_diff_cutoff})"
            ),
            "priority": 42,
            "reason_key": "frequency_difference",
            "policy_key": "external_af_concordance",
            **rule_display_groups["external_af_concordance"],
        })

    # CRITICAL: `bcftools view -i "1"` matches NOTHING (verified against bcftools
    # 1.23: 6 records in, 0 out, exit 0, empty stderr).  "match everything" means
    # omitting -i altogether, not passing a truthy-looking placeholder.
    include_expr = " & ".join(include_parts) if include_parts else None
    if include_expr is None and empty_expression == "match_none":
        include_expr = "1"
        condition_checks.append({
            "label": "Empty active-filter policy configured to remove every variant",
            "expr": f"({position_field} >= 1)",
            "priority": 1,
            "reason_key": "empty_expression",
            **rule_display_groups["filtering_policy"],
        })

    log_print("🔧 Include expression (bcftools -i):")
    if include_expr is None:
        log_print("(none - no include filter is active, so -i is omitted and every variant is kept)")
    else:
        log_print(include_expr)
        if not include_parts:
            log_warn(
                "⚠️ empty_expression='match_none': no include filter is active, and the "
                "placeholder expression '1' makes bcftools match NO records. The output "
                "VCF will be empty."
            )
    log_print("")

    # ============================================================
    # EXCLUDE LOGIC
    # ============================================================
    palindromic_logic = (
        f"(({reference_field}=='A' & {alternate_field}=='T') | "
        f"({reference_field}=='T' & {alternate_field}=='A') | "
        f"({reference_field}=='C' & {alternate_field}=='G') | "
        f"({reference_field}=='G' & {alternate_field}=='C'))"
    )
    exclude_parts: List[str] = []
    if exclude_palindromic:
        palindromic_expr = (
            f"({palindromic_logic} & "
            f"({study_af_format} >= {palindromic_af_lower} & "
            f"{study_af_format} <= {palindromic_af_upper}))"
        )
        exclude_parts.append(palindromic_expr)
        condition_checks.append({
            "label": "Palindromic SNPs with ambiguous frequency (%.6g ≤ AF ≤ %.6g)"
            % (palindromic_af_lower, palindromic_af_upper),
            "expr": (
                "(%s & %s != '.' & %s >= %s & %s <= %s)"
                % (
                    palindromic_logic,
                    study_af_format,
                    study_af_format,
                    palindromic_af_lower,
                    study_af_format,
                    palindromic_af_upper,
                )
            ),
            "priority": 60,
            "reason_key": "palindromic",
            "policy_key": "palindromic_variants",
            **rule_display_groups["palindromic_variants"],
        })
    if not include_indels:
        non_snp_expr = f"({variant_type_field} != 'snp')"
        exclude_parts.append(non_snp_expr)
        condition_checks.append({
            "label": "Indels and other non-SNP variants",
            "expr": non_snp_expr,
            "priority": 50,
            "reason_key": "non_snp",
            "policy_key": "variant_type",
            **rule_display_groups["variant_type"],
        })
    if remove_mhc:
        condition_checks.append({
            "label": "Variants in the MHC region (%s:%s-%s)"
            % (mhc_chrom, f"{int(mhc_start):,}", f"{int(mhc_end):,}"),
            "expr": (
                f"({chromosome_field} == '{mhc_chrom}' & "
                f"{position_field} >= {int(mhc_start)} & "
                f"{position_field} <= {int(mhc_end)})"
            ),
            "priority": 70,
            "reason_key": "mhc",
            "policy_key": "mhc_region",
            **rule_display_groups["mhc_region"],
        })
    constructed_policy_keys = {
        check["policy_key"] for check in condition_checks
        if check.get("policy_key") is not None
    }
    expected_policy_keys = set(policy.active_rule_keys())
    if constructed_policy_keys != expected_policy_keys:
        raise RuntimeError(
            "Filtering expressions drifted from the shared scientific-policy "
            "contract: %r != %r"
            % (sorted(constructed_policy_keys), sorted(expected_policy_keys))
        )
    exclude_expr = " | ".join(exclude_parts) if exclude_parts else None

    # The display, primary attribution and report all use the same scientific
    # execution order. This makes "first failed rule" deterministic and visible.
    missing_checks.sort(key=lambda item: item.get("priority", 999))
    condition_checks.sort(key=lambda item: item.get("priority", 999))

    reason_checks = []
    for check in missing_checks:
        check.update({
            "category": "missing_value",
            "label": "Missing %s" % check["field"],
            "removes": check["action"] == "remove",
            "tag": filter_reason_ids[check["reason_key"]],
        })
        reason_checks.append(check)
    for check in condition_checks:
        check.update({
            "category": "filter_condition", "action": "remove", "removes": True,
            "tag": filter_reason_ids[check["reason_key"]],
        })
        reason_checks.append(check)
    reason_checks.sort(key=lambda item: item.get("priority", 999))
    plan_text = render_filtering_plan(reason_checks)
    for line in plan_text.splitlines():
        log_print(line)
    progress.print_block(style_screen_block(plan_text))

    temporary_soft_vcf = None
    temporary_soft_index = None
    if write_soft_filter_vcf:
        soft_handle = tempfile.NamedTemporaryFile(
            prefix=".%s." % soft_output_vcf.name,
            suffix=".vcf.gz",
            dir=soft_output_vcf.parent,
            delete=False,
        )
        temporary_soft_vcf = Path(soft_handle.name)
        soft_handle.close()
        temporary_soft_index = Path(str(temporary_soft_vcf) + ".tbi")

    audit_checks = (
        [check for check in reason_checks if check.get("removes")]
        if write_soft_filter_vcf else reason_checks
    )
    with progress.step(2, total_stages, "Audit active filter-rule failures"):
        reason_statistics = _collect_filter_reason_statistics(
            vcf_path=str(input_vcf),
            checks=audit_checks,
            output_folder=str(step_dir),
            output_prefix=output_prefix,
            log_warn=log_warn,
            bcftools_bin=bcftools_bin,
            bash_bin=bash_bin,
            tabix_bin=tabix_bin,
            threads=threads,
            max_mem=max_mem,
            sort_output=sort_output,
            soft_vcf_path=temporary_soft_vcf,
        )
        if reason_statistics is None:
            if temporary_soft_vcf is not None:
                temporary_soft_vcf.unlink(missing_ok=True)
            if temporary_soft_index is not None:
                temporary_soft_index.unlink(missing_ok=True)
            log_print(
                "❌ Filter-reason tagging failed; no filtering outputs were published."
            )
            log_print("❌ STATUS: FAILED")
            write_log()
            raise RuntimeError("Filter-reason tagging failed")

        by_tag = {
            check["tag"]: check for check in reason_statistics["checks"]
        }
        for check in reason_checks:
            measured = by_tag.get(check["tag"])
            if measured is not None:
                check["count"] = measured["count"]
                check["primary_count"] = measured["primary_count"]
            else:
                check["count"] = count_matching(str(input_vcf), check["expr"])
                check["primary_count"] = 0
                if check["count"] is None:
                    if temporary_soft_vcf is not None:
                        temporary_soft_vcf.unlink(missing_ok=True)
                    if temporary_soft_index is not None:
                        temporary_soft_index.unlink(missing_ok=True)
                    log_print(
                        "❌ A non-removing diagnostic count failed; no filtering "
                        "outputs were published."
                    )
                    log_print("❌ STATUS: FAILED")
                    write_log()
                    raise RuntimeError("Filter-reason diagnostic counting failed")
        reason_statistics["checks"] = [dict(check) for check in reason_checks]
        thread_runtime = reason_statistics["aggregation"]
        log_print(
            "Thread control: requested=%d; bcftools_compression=%d; "
            "polars_pool=%s; polars_budget_enforced=%s"
            % (
                threads,
                threads,
                thread_runtime.get("polars_thread_pool_size"),
                thread_runtime.get("thread_budget_enforced"),
            )
        )

    missing_value_counts = dict(
        (check["field"], check.get("count")) for check in missing_checks
    )
    missing_value_drops = dict(
        (check["field"], check.get("primary_count"))
        for check in missing_checks if check["action"] == "remove"
    )

    log_print("🔧 Exclude expression (bcftools -e):")
    log_print(str(exclude_expr))
    log_print("")

    # ============================================================
    # PIPELINE
    # ============================================================
    cmd_parts = []
    temporary_bed = None
    # speed edit: no --threads on intermediate uncompressed stream step
    cmd1 = "%s view -Ou" % quoted_bcftools
    if include_expr is not None:
        cmd1 += " -i %s" % shlex.quote(include_expr)
    cmd1 += f" {quoted_input}"
    cmd_parts.append(cmd1)
    # keep -e separate from -i
    if exclude_expr:
        # speed edit: no --threads on intermediate uncompressed stream step
        cmd_parts.append(
            "%s view -Ou -e %s"
            % (quoted_bcftools, shlex.quote(exclude_expr))
        )
    if remove_mhc:
        # BED is 0-based half-open; mhc_start is a 1-based inclusive coordinate,
        # so the start field must be one less or the first base of the region
        # survives the exclusion.
        bed_start = max(0, int(mhc_start) - 1)
        bed_handle = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".%s." % mhc_bed.name,
            suffix=".bed",
            dir=mhc_bed.parent,
            delete=False,
        )
        temporary_bed = Path(bed_handle.name)
        try:
            bed_handle.write(f"{mhc_chrom}\t{bed_start}\t{int(mhc_end)}\n")
            bed_handle.flush()
            os.fsync(bed_handle.fileno())
            bed_handle.close()
        except BaseException:
            bed_handle.close()
            temporary_bed.unlink(missing_ok=True)
            raise
        # speed edit: no --threads on intermediate uncompressed stream step
        cmd_parts.append(
            "%s view -Ou -T %s"
            % (quoted_bcftools, shlex.quote("^" + str(temporary_bed)))
        )
    # ============================================================
    # SPEED OPTIMIZATION
    # - keep -i and -e separate
    # - skip expensive sort unless configured
    # ============================================================
    if sort_output:
        cmd_parts.append(
            "%s sort --temp-dir %s --max-mem %s"
            % (
                quoted_bcftools,
                shlex.quote(str(step_dir)),
                shlex.quote(str(max_mem)),
            )
        )
    # keep threads here where they help most: final bgzip compression
    cmd_parts.append("%s view -Oz --threads %d" % (quoted_bcftools, threads))
    pipeline_core = " | ".join(cmd_parts)

    temporary_handle = tempfile.NamedTemporaryFile(
        prefix=".%s." % output_vcf.name,
        suffix=".vcf.gz",
        dir=output_vcf.parent,
        delete=False,
    )
    temporary_vcf = Path(temporary_handle.name)
    temporary_handle.close()
    temporary_index = Path(str(temporary_vcf) + ".tbi")

    # Without pipefail only the exit status of the LAST stage is seen, so an
    # abort in the middle of the pipe produced a silently TRUNCATED VCF whose
    # reduced count was then reported as "variants removed by filtering".
    inner = (
        f"set -euo pipefail; {pipeline_core} > {shlex.quote(str(temporary_vcf))} "
        f"&& {quoted_tabix} -f -p vcf {shlex.quote(str(temporary_vcf))}"
    )
    pipeline = "%s -c %s" % (quoted_bash, shlex.quote(inner))
    log_print("🚀 Full bcftools pipeline:")
    log_print(inner)
    log_print("")
    with progress.step(3, total_stages, "Create and index hard-filtered VCF"):
        try:
            run_cmd(pipeline)
            if (
                not temporary_vcf.is_file()
                or temporary_vcf.stat().st_size <= 0
                or not temporary_index.is_file()
                or temporary_index.stat().st_size <= 0
            ):
                raise RuntimeError(
                    "bcftools/tabix completed without a non-empty temporary VCF and index"
                )
        except BaseException as exc:
            log_print(
                "❌ Error during bcftools pipeline execution: %s: %s"
                % (type(exc).__name__, exc)
            )
            log_print("❌ STATUS: FAILED")
            write_log()
            temporary_vcf.unlink(missing_ok=True)
            temporary_index.unlink(missing_ok=True)
            if temporary_soft_vcf is not None:
                temporary_soft_vcf.unlink(missing_ok=True)
            if temporary_soft_index is not None:
                temporary_soft_index.unlink(missing_ok=True)
            if temporary_bed is not None:
                temporary_bed.unlink(missing_ok=True)
            raise RuntimeError("bcftools filtering pipeline failed") from exc
        post_variants = fast_count(str(temporary_vcf))
        log_print(
            "📊 Variants counted with: bcftools index -n "
            "(fallback: bcftools view -H | wc -l)"
        )
        log_print(
            "✅ Variants AFTER filtering                          : %s"
            % _fmt_count(post_variants)
        )

        if pre_variants and post_variants == 0:
            log_warn(
                "⚠️ Every variant was removed by filtering. Check the include "
                "expression above: a filter on a tag the VCF does not carry "
                "removes everything."
            )

    with progress.step(4, total_stages, "Reconcile filter-removal accounting"):
        actual_removed = (
            None if post_variants is None else pre_variants - post_variants
        )
        soft_variants = (
            fast_count(str(temporary_soft_vcf))
            if temporary_soft_vcf is not None else None
        )
        reason_statistics["actual_removed"] = actual_removed
        reason_statistics["reconciled"] = (
            actual_removed is not None
            and reason_statistics.get("primary_removed_total") == actual_removed
        )
        soft_reconciled = (
            temporary_soft_vcf is None
            or (
                soft_variants is not None
                and pre_variants == soft_variants
            )
        )
        if not reason_statistics["reconciled"] or not soft_reconciled:
            log_print(
                "❌ Filter validation failed: hard VCF removed %s variants; reason "
                "tags account for %s; soft VCF retained %s of %s input variants."
                % (
                    _fmt_count(actual_removed),
                    _fmt_count(reason_statistics.get("primary_removed_total")),
                    _fmt_count(soft_variants),
                    _fmt_count(pre_variants),
                )
            )
            log_print("❌ No filtering outputs were published.")
            log_print("❌ STATUS: FAILED")
            temporary_vcf.unlink(missing_ok=True)
            temporary_index.unlink(missing_ok=True)
            if temporary_soft_vcf is not None:
                temporary_soft_vcf.unlink(missing_ok=True)
            if temporary_soft_index is not None:
                temporary_soft_index.unlink(missing_ok=True)
            if temporary_bed is not None:
                temporary_bed.unlink(missing_ok=True)
            write_log()
            raise RuntimeError("Filter-reason reconciliation failed")
    def temporary_report_path(path: Path, suffix: str) -> Path:
        handle = tempfile.NamedTemporaryFile(
            prefix=".%s." % path.name,
            suffix=suffix,
            dir=path.parent,
            delete=False,
        )
        temporary_path = Path(handle.name)
        handle.close()
        temporary_path.unlink()
        return temporary_path

    temporary_report = temporary_report_path(reason_report_path, ".tsv")
    temporary_summary_csv = temporary_report_path(summary_csv_path, ".csv")
    temporary_html_report = temporary_report_path(html_report_path, ".html")
    output_index = Path(str(output_vcf) + ".tbi")
    soft_output_index = (
        Path(str(soft_output_vcf) + ".tbi") if write_soft_filter_vcf else None
    )
    report_outputs = {
        "input_vcf": str(input_vcf),
        "filtered_vcf": str(output_vcf),
        "filtered_vcf_index": str(output_index),
        "soft_filtered_vcf": (
            str(soft_output_vcf) if write_soft_filter_vcf else None
        ),
        "soft_filtered_vcf_index": (
            str(soft_output_index) if soft_output_index is not None else None
        ),
        "reason_summary": str(reason_report_path),
        "summary_csv": str(summary_csv_path),
        "html_report": str(html_report_path),
        "filter_log": str(log_file),
        "mhc_exclusion_bed": str(mhc_bed) if remove_mhc else None,
    }
    summary = build_filtering_summary(
        dataset_id=output_prefix,
        genome_build=genome_build,
        missing_checks=missing_checks,
        condition_checks=condition_checks,
        reason_statistics=reason_statistics,
        variants_before=pre_variants,
        variants_after=post_variants,
        soft_filter_variants=soft_variants,
        soft_filter_enabled=write_soft_filter_vcf,
        input_validation=input_validation,
        data_flow=None,
        outputs=report_outputs,
        resolved_configuration=resolved_configuration,
        runtime_seconds=time.monotonic() - started_at,
        display_missing_counts=report_missing_counts,
    )
    summary_text = render_filtering_summary(
        summary, label_width=terminal_label_width,
    )
    summary_lines = summary_text.splitlines()
    with progress.step(5, total_stages, "Save and publish filtering reports"):
        try:
            _write_filter_reason_report(reason_statistics, temporary_report)
            write_filtering_summary_csv(summary, temporary_summary_csv)
            write_filtering_html_report(summary, temporary_html_report)
            artifacts = [
                (temporary_vcf, output_vcf),
                (temporary_index, output_index),
            ]
            if temporary_soft_vcf is not None and temporary_soft_index is not None:
                artifacts.extend((
                    (temporary_soft_vcf, soft_output_vcf),
                    (temporary_soft_index, soft_output_index),
                ))
            if temporary_bed is not None:
                artifacts.append((temporary_bed, mhc_bed))
            artifacts.extend((
                (temporary_report, reason_report_path),
                (temporary_summary_csv, summary_csv_path),
                (temporary_html_report, html_report_path),
            ))
            publish_artifact_set(artifacts)
        except BaseException as exc:
            temporary_vcf.unlink(missing_ok=True)
            temporary_index.unlink(missing_ok=True)
            if temporary_soft_vcf is not None:
                temporary_soft_vcf.unlink(missing_ok=True)
            if temporary_soft_index is not None:
                temporary_soft_index.unlink(missing_ok=True)
            temporary_report.unlink(missing_ok=True)
            temporary_summary_csv.unlink(missing_ok=True)
            temporary_html_report.unlink(missing_ok=True)
            if temporary_bed is not None:
                temporary_bed.unlink(missing_ok=True)
            log_print(
                "❌ Filtering output finalisation failed: %s: %s"
                % (type(exc).__name__, exc)
            )
            log_print("❌ STATUS: FAILED")
            write_log()
            raise
    reason_report = str(reason_report_path)
    if reason_report:
        log_print("Detailed filter reason report path: %s" % reason_report)
    log_print(f"💾 Filtered VCF saved to                             : {output_vcf}")
    if write_soft_filter_vcf:
        log_print(
            "💾 Soft-filter audit VCF saved to                   : %s"
            % soft_output_vcf
        )
    log_print("Summary CSV path: %s" % summary_csv_path)
    log_print("Detailed HTML report path: %s" % html_report_path)
    for line in summary_lines:
        log_print(line)
    log_print("⏱️ Filtering runtime seconds                         : %.3f" % (
        time.monotonic() - started_at
    ))
    log_print("✅ STATUS: COMPLETED")
    log_print("\n🎉 bcftools filtering completed.")
    write_log()

    # The terminal shows observed consequences, not a second copy of the
    # configuration. Emit the complete block while any enclosing progress is
    # paused so a live refresh cannot overwrite wrapped summary content.
    progress.print_block(style_screen_block(summary_text))
    return {
        "genome_build": genome_build,
        "filtered_vcf": str(output_vcf),
        "filtered_vcf_index": str(Path(str(output_vcf) + ".tbi")),
        "soft_filtered_vcf": (
            str(soft_output_vcf) if write_soft_filter_vcf else None
        ),
        "soft_filtered_vcf_index": (
            str(Path(str(soft_output_vcf) + ".tbi"))
            if write_soft_filter_vcf else None
        ),
        "filter_log": str(log_file),
        "mhc_exclusion_bed": str(mhc_bed) if remove_mhc else None,
        "variants_before": pre_variants,
        "variants_after": post_variants,
        "missing_value_drops": missing_value_drops,
        "missing_value_counts": missing_value_counts,
        "filter_condition_counts": {
            check["label"]: check.get("count") for check in condition_checks
        },
        "filter_primary_reason_counts": {
            check["label"]: check.get("primary_count") for check in reason_checks
            if check.get("removes")
        },
        "filter_overlap_variants": (
            reason_statistics.get("overlap_variants")
            if reason_statistics is not None else None
        ),
        "filter_reason_reconciled": (
            reason_statistics.get("reconciled")
            if reason_statistics is not None else None
        ),
        "thread_runtime": reason_statistics.get("aggregation"),
        "variant_qc_policy": policy.as_dict(),
        "filter_reason_summary": reason_report,
        "filter_summary_csv": str(summary_csv_path),
        "filter_html_report": str(html_report_path),
        "filter_summary": summary,
        "input_validation": input_validation,
    }


def _match_contig_naming(
    mhc_chrom: str,
    contigs: List[str],
    log_warn,
) -> str:
    """Return the MHC contig name as the VCF actually spells it.

    ``bcftools view -T ^bed`` silently excludes nothing when the BED names a
    contig the VCF does not contain, so `6` against a `chr6` VCF (or the
    reverse) turns MHC removal into a no-op with no error anywhere.
    """
    if not contigs:
        raise ValueError(
            "MHC removal requires ##contig declarations in the input VCF header; "
            "none were found"
        )
    if mhc_chrom in contigs:
        return mhc_chrom

    alternative = (
        mhc_chrom[3:] if mhc_chrom.lower().startswith("chr") else "chr" + mhc_chrom
    )
    if alternative in contigs:
        log_warn(
            "⚠️ MHC contig '%s' is not in the VCF, but '%s' is — using '%s' so the "
            "MHC exclusion actually applies." % (mhc_chrom, alternative, alternative)
        )
        return alternative

    raise ValueError(
        "MHC contig %r is absent from the input VCF header (contigs seen: %s%s). "
        "Set the matching modules.filtering.mhc_regions.<build>.chromosome value."
        % (
            mhc_chrom,
            ", ".join(contigs[:10]),
            ", ..." if len(contigs) > 10 else "",
        )
    )
