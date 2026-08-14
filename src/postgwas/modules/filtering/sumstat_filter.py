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

from postgwas.core.execution.runtime import run_cmd
from postgwas.core.paths import configured_output_path, validate_filename_component
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.vcf import (
    VCF_TAG,
    read_vcf_header,
    validate_vcf_header_contract,
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
) -> List[str]:
    """Build a detailed report with raw and mutually exclusive reason counts."""
    reason_statistics = reason_statistics or {}
    data_flow = data_flow or {}
    lines = [
        screen_line("analysis", "Variant filtering summary", indent=4),
        screen_field(
            "info", "Attribution method",
            "a variant can fail more than one rule. To avoid double-counting, each removed "
            "variant is assigned to the first rule it fails; the assignment order is saved in the report",
            indent=8, label_width=terminal_label_width,
        ),
        "",
        screen_line("count", "Variant flow", indent=6),
    ]
    flow_fields = [
        ("Input summary statistics", "total_variant_infile"),
        ("Read by harmonisation", "total_variant_read"),
        ("Invalid coordinates removed", "total_variant_removed_null_coords"),
        ("Non-standard alleles removed", "total_variant_removed_non_standard_alleles"),
        ("Remaining for harmonisation", "total_variant_remaining_for_harmonisation"),
        ("Missing effect frequency", "total_variant_with_missing_eaf"),
        ("Invalid effect statistics", "total_variant_with_invalid_beta_se"),
        ("Used for VCF creation", "total_variant_in_vcf_input"),
    ]
    for label, key in flow_fields:
        if key in data_flow:
            lines.append(screen_field(
                "loss" if "removed" in label.lower() or "invalid" in label.lower()
                else "count",
                label, _fmt_count(data_flow.get(key)),
                indent=8, label_width=terminal_label_width,
            ))
    lines.extend([
        screen_field(
            "count", "VCF before filtering", _fmt_count(variants_before),
            indent=8, label_width=terminal_label_width,
        ),
        screen_field(
            "count", "VCF after filtering", _fmt_count(variants_after),
            indent=8, label_width=terminal_label_width,
        ),
        "",
        screen_line("genetic", "Missing values used by active filters", indent=6),
    ])
    if missing_checks:
        for check in missing_checks:
            lines.extend(screen_field(
                "loss" if check.get("action") == "remove" else "info",
                check["field"],
                "missing %s; removed for this reason %s; action %s" % (
                    _fmt_count(check.get("count")),
                    _fmt_count(check.get("primary_count")),
                    check.get("action"),
                ),
                indent=8, label_width=terminal_label_width,
            ).splitlines())
    else:
        lines.append(screen_line(
            "info", "No active filter depends on a nullable value", indent=8,
        ))

    lines.extend([
        "",
        screen_line("analysis", "Active filtering conditions", indent=6),
    ])
    if condition_checks:
        for check in condition_checks:
            lines.extend(screen_field(
                "loss", check["label"],
                "matched %s; removed for this reason %s" % (
                    _fmt_count(check.get("count")),
                    _fmt_count(check.get("primary_count")),
                ),
                indent=8, label_width=terminal_label_width,
            ).splitlines())
    else:
        lines.append(screen_line(
            "info", "No variant-level filtering condition is active", indent=8,
        ))

    lines.extend([
        "",
        screen_line("count", "Exact removal accounting", indent=6),
        screen_field(
            "analysis", "Multiple failed rules",
            _fmt_count(reason_statistics.get("overlap_variants")),
            indent=8, label_width=terminal_label_width,
        ),
        screen_field(
            "analysis", "Overlapping rule matches",
            _fmt_count(reason_statistics.get("extra_rule_matches")),
            indent=8, label_width=terminal_label_width,
        ),
    ])
    removed = None
    if variants_before is not None and variants_after is not None:
        removed = variants_before - variants_after
    lines.append(screen_field(
        "loss", "Removed by all filters", _fmt_count(removed),
        indent=8, label_width=terminal_label_width,
    ))

    attributed = reason_statistics.get("primary_removed_total")
    reconciled = (
        removed is not None and attributed is not None and int(removed) == int(attributed)
    )
    if removed is not None and attributed is not None:
        lines.append(screen_field(
            "success" if reconciled else "error",
            "Removal count check",
            "%s assigned to reasons %s %s removed in total" % (
                _fmt_count(attributed), "=" if reconciled else "!=",
                _fmt_count(removed),
            ),
            indent=8, label_width=terminal_label_width,
        ))
    else:
        lines.append(screen_field(
            "warning", "Removal count check", "unavailable",
            indent=8, label_width=terminal_label_width,
        ))
    if reason_report:
        lines.append(screen_field(
            "info", "Detailed reason report", os.path.basename(reason_report),
            indent=8, label_width=terminal_label_width,
        ))
    return lines


def _summarize_filter_tags(tag_file: Path, checks: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Use one lazy Polars aggregation to count raw, overlapping and primary reasons."""
    scan = pl.scan_csv(
        tag_file,
        has_header=False,
        separator="\t",
        schema={"filter_tags": pl.String},
    )
    tags = pl.col("filter_tags").fill_null("").str.split(";")
    masks = dict((check["tag"], tags.list.contains(check["tag"])) for check in checks)
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

    values = scan.select(selections).collect().to_dicts()[0]
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
    }


def _collect_filter_reason_statistics(
    vcf_path: str,
    checks: List[Dict[str, Any]],
    output_folder: str,
    output_prefix: str,
    log_warn,
    *,
    bcftools_bin: str,
    bash_bin: str,
) -> Optional[Dict[str, Any]]:
    """Tag all filter reasons in one bcftools stream, then aggregate with Polars."""
    if not checks:
        return {
            "checks": [], "primary_removed_total": 0,
            "overlap_variants": 0, "extra_rule_matches": 0,
        }

    for index, check in enumerate(checks, 1):
        check["tag"] = "PGWAS_FILTER_%02d" % index

    handle = tempfile.NamedTemporaryFile(
        prefix=".%s_filter_reason_tags_" % output_prefix,
        suffix=".tsv",
        dir=output_folder,
        delete=False,
    )
    tag_file = Path(handle.name)
    handle.close()
    quoted_bcftools = shlex.quote(bcftools_bin)
    stages = [
        "%s view -Ou %s" % (quoted_bcftools, shlex.quote(str(vcf_path)))
    ]
    for check in checks:
        stages.append(
            "%s filter -Ou -m + -s %s -e %s"
            % (quoted_bcftools, check["tag"], shlex.quote(check["expr"]))
        )
    stages.append("%s query -f '%%FILTER\\n'" % quoted_bcftools)
    inner = "set -euo pipefail; %s > %s" % (
        " | ".join(stages), shlex.quote(str(tag_file)),
    )
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
        return _summarize_filter_tags(tag_file, checks)
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
                "order", "category", "reason", "action",
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
            ("Overlapping rule matches beyond the first",
             statistics.get("extra_rule_matches"), None, None),
            ("All mutually exclusive primary reasons",
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


def _read_vcf_header(vcf_path: Path, bcftools_bin: str) -> str:
    """Compatibility wrapper around the shared VCF contract reader."""
    return read_vcf_header(vcf_path, bcftools_bin, error_type=ValueError)


def _validate_vcf_contract(
    *,
    header: str,
    genome_build_header: str,
    supported_genome_builds: List[str],
    required_fields: List[str],
) -> tuple[str, List[str]]:
    """Compatibility wrapper around the shared VCF contract validator."""
    return validate_vcf_header_contract(
        header=header,
        genome_build_header=genome_build_header,
        supported_genome_builds=supported_genome_builds,
        required_fields=required_fields,
    )


def filter_gwas_vcf_bcftools(
    vcf_path: str,
    output_folder: str,
    output_prefix: str,
    *,
    genome_build_header: str,
    supported_genome_builds: List[str],
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
    threads: int,
    max_mem: str,
    lp_missing: str,
    af_missing: str,
    empty_expression: str,
    sort_output: bool,
    vcf_fields: Mapping[str, str],
    output_layout: Mapping[str, str],
    bcftools_bin: str,
    tabix_bin: str,
    bash_bin: str,
    resolved_configuration: Mapping[str, Any],
    report_missing_counts: bool,
    terminal_label_width: int,
) -> Dict[str, Any]:
    """Filter a GWAS-VCF and report, honestly, how many variants each rule cost.

    All scientific policies, VCF fields, output paths, and executable paths are
    supplied from the schema-validated canonical filtering configuration.
    """
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

    # ============================================================
    # VALIDATION
    # ============================================================
    log_print(
        "Resolved configuration:\n%s"
        % json.dumps(resolved_configuration, indent=2, sort_keys=True)
    )
    input_vcf = Path(vcf_path).expanduser().resolve()
    if not input_vcf.is_file() or input_vcf.stat().st_size <= 0:
        msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
        log_print(msg)
        log_print("❌ STATUS: FAILED")
        write_log()
        raise FileNotFoundError(msg)
    try:
        header = _read_vcf_header(input_vcf, bcftools_bin)
        required_fields = []
        if pval_cutoff is not None:
            required_fields.append(log_pvalue_format)
        if maf_cutoff is not None or exclude_palindromic:
            required_fields.append(study_af_format)
        if info_min is not None or info_max is not None:
            required_fields.append(imputation_quality_format)
        if allelefreq_diff_cutoff is not None:
            required_fields.extend((study_af_info, external_af_info))
        genome_build, contigs = _validate_vcf_contract(
            header=header,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            required_fields=required_fields,
        )
        mhc_chrom = None
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
                mhc_chrom = str(mhc_region["chromosome"])
                mhc_start = int(mhc_region["start"])
                mhc_end = int(mhc_region["end"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    "The configured MHC region for %s is incomplete or invalid"
                    % genome_build
                ) from exc
            mhc_chrom = _match_contig_naming(
                str(mhc_chrom), contigs, log_warn,
            )
    except BaseException as exc:
        log_print("❌ VCF validation failed: %s: %s" % (type(exc).__name__, exc))
        log_print("❌ STATUS: FAILED")
        write_log()
        raise
    log_print(
        "✅ VCF contract validated: genome_build=%s; active_fields=%s"
        % (genome_build, ", ".join(required_fields) or "none")
    )
    if remove_mhc:
        log_print(
            "✅ Build-specific MHC region selected: %s:%d-%d (%s)"
            % (mhc_chrom, mhc_start, mhc_end, genome_build)
        )
    path_values = {"dataset_id": output_prefix, "genome_build": genome_build}
    log_file = configured_output_path(
        step_dir, output_layout["log_file"], **path_values
    )
    output_vcf = configured_output_path(
        step_dir, output_layout["filtered_vcf"], **path_values
    )
    reason_report_path = configured_output_path(
        step_dir, output_layout["reason_summary"], **path_values
    )
    mhc_bed = configured_output_path(
        step_dir, output_layout["mhc_exclusion_bed"], **path_values
    )
    for path in (log_file, output_vcf, reason_report_path, mhc_bed):
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
    quoted_input = shlex.quote(str(input_vcf))
    quoted_bcftools = shlex.quote(bcftools_bin)
    quoted_bash = shlex.quote(bash_bin)
    quoted_tabix = shlex.quote(tabix_bin)

    # ============================================================
    # FAST COUNT
    # ============================================================
    def fast_count(vcf: str) -> Optional[int]:
        """Number of records in an indexed VCF, or None (never a silent 0)."""
        quoted = shlex.quote(str(vcf))
        problems = []
        try:
            res = run_cmd(
                f"{quoted_bcftools} index -n {quoted}", check=False
            )
            if res.returncode == 0:
                return int(res.stdout.strip())
            problems.append(f"`bcftools index -n` exit {res.returncode}: {res.stderr.strip()}")
        except Exception as exc:  # run_cmd itself blew up
            problems.append(f"`bcftools index -n` raised {type(exc).__name__}: {exc}")

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
            problems.append(f"`bcftools view -H | wc -l` exit {res.returncode}: {res.stderr.strip()}")
        except Exception as exc:
            problems.append(f"`bcftools view -H | wc -l` raised {type(exc).__name__}: {exc}")

        log_warn(
            "⚠️ Could not count the variants in %s, so every count derived from it is "
            "reported as 'unavailable' rather than 0. Reasons: %s" % (vcf, "; ".join(problems))
        )
        return None

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

    pre_variants = fast_count(str(input_vcf))
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
        })
        include_parts.append(lp_expr)
        condition_checks.append({
            "label": "P-value evidence below LP %.6g" % resolved_lp,
            "expr": (
                f"({log_pvalue_format} != '.' & "
                f"{log_pvalue_format} < {resolved_lp})"
            ),
            "priority": 11,
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
            "label": info_label, "expr": info_failure_expr, "priority": 31,
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
        }, {
            "field": external_af_info,
            "expr": f"{external_af_info} == '.'",
            "action": af_missing,
            "key": "af_missing",
            "priority": 41,
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
        })

    # CRITICAL: `bcftools view -i "1"` matches NOTHING (verified against bcftools
    # 1.23: 6 records in, 0 out, exit 0, empty stderr).  "match everything" means
    # omitting -i altogether, not passing a truthy-looking placeholder.
    include_expr = " & ".join(include_parts) if include_parts else None
    if include_expr is None and empty_expression == "match_none":
        include_expr = "1"

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
    exclude_expr = None
    if exclude_palindromic:
        exclude_expr = (
            f"({palindromic_logic} & "
            f"({study_af_format} >= {palindromic_af_lower} & "
            f"{study_af_format} <= {palindromic_af_upper}))"
        )
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
        })
    if not include_indels:
        condition_checks.append({
            "label": "Indels and other non-SNP variants",
            "expr": f"({variant_type_field} != 'snp')",
            "priority": 50,
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
        })

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
        })
        reason_checks.append(check)
    for check in condition_checks:
        check.update({
            "category": "filter_condition", "action": "remove", "removes": True,
        })
        reason_checks.append(check)
    reason_checks.sort(key=lambda item: item.get("priority", 999))

    reason_statistics = _collect_filter_reason_statistics(
        vcf_path=str(input_vcf),
        checks=reason_checks,
        output_folder=str(step_dir),
        output_prefix=output_prefix,
        log_warn=log_warn,
        bcftools_bin=bcftools_bin,
        bash_bin=bash_bin,
    )
    if reason_statistics is None:
        # The exact single-pass audit failed, but retain the previous independent
        # diagnostics rather than fabricating zeroes.
        for check in reason_checks:
            check["count"] = count_matching(vcf_path, check["expr"])
            check["primary_count"] = None
    else:
        by_tag = dict((check["tag"], check) for check in reason_statistics["checks"])
        for check in reason_checks:
            measured = by_tag[check["tag"]]
            check["count"] = measured["count"]
            check["primary_count"] = measured["primary_count"]

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
    if not include_indels:
        cmd1 += " --types snps"
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
        if temporary_bed is not None:
            temporary_bed.unlink(missing_ok=True)
        raise RuntimeError("bcftools filtering pipeline failed") from exc
    # ============================================================
    # Count variants after filtering
    # ============================================================
    post_variants = fast_count(str(temporary_vcf))
    log_print("📊 Variants counted with: bcftools index -n (fallback: bcftools view -H | wc -l)")
    log_print(f"✅ Variants AFTER filtering                          : {_fmt_count(post_variants)}")

    if pre_variants and post_variants == 0:
        log_warn(
            "⚠️ Every variant was removed by filtering. Check the include expression above: "
            "a filter on a tag the VCF does not carry removes everything."
        )

    actual_removed = (
        pre_variants - post_variants
        if pre_variants is not None and post_variants is not None else None
    )
    if reason_statistics is not None:
        reason_statistics["actual_removed"] = actual_removed
        reason_statistics["reconciled"] = (
            actual_removed is not None
            and reason_statistics.get("primary_removed_total") == actual_removed
        )
        if actual_removed is not None and not reason_statistics["reconciled"]:
            log_warn(
                "Filter-reason reconciliation failed: %s variants were removed from "
                "the final VCF, but the mutually exclusive reasons account for %s."
                % (
                    _fmt_count(actual_removed),
                    _fmt_count(reason_statistics.get("primary_removed_total")),
                )
            )
    report_statistics = reason_statistics or {
        "checks": reason_checks,
        "primary_removed_total": None,
        "overlap_variants": None,
        "extra_rule_matches": None,
        "actual_removed": actual_removed,
        "reconciled": False,
    }
    report_handle = tempfile.NamedTemporaryFile(
        prefix=".%s." % reason_report_path.name,
        suffix=".tsv",
        dir=reason_report_path.parent,
        delete=False,
    )
    temporary_report = Path(report_handle.name)
    report_handle.close()
    temporary_report.unlink()
    summary_lines = _filtering_summary_lines(
        missing_checks if report_missing_counts else [],
        condition_checks,
        pre_variants,
        post_variants,
        terminal_label_width,
        reason_statistics=reason_statistics,
        reason_report=str(reason_report_path),
    )
    try:
        _write_filter_reason_report(
            report_statistics, temporary_report,
        )
        output_index = Path(str(output_vcf) + ".tbi")
        os.replace(temporary_vcf, output_vcf)
        os.replace(temporary_index, output_index)
        if temporary_bed is not None:
            os.replace(temporary_bed, mhc_bed)
        os.replace(temporary_report, reason_report_path)
        reason_report = str(reason_report_path)
    except BaseException as exc:
        temporary_vcf.unlink(missing_ok=True)
        temporary_index.unlink(missing_ok=True)
        temporary_report.unlink(missing_ok=True)
        if temporary_bed is not None:
            temporary_bed.unlink(missing_ok=True)
        log_print(
            "❌ Filtering output finalisation failed: %s: %s"
            % (type(exc).__name__, exc)
        )
        log_print("❌ STATUS: FAILED")
        write_log()
        raise
    if reason_report:
        log_print("Detailed filter reason report path: %s" % reason_report)
    log_print(f"💾 Filtered VCF saved to                             : {output_vcf}")
    for line in summary_lines:
        log_print(line)
    log_print("⏱️ Filtering runtime seconds                         : %.3f" % (
        time.monotonic() - started_at
    ))
    log_print("✅ STATUS: COMPLETED")
    log_print("\n🎉 bcftools filtering completed.")
    write_log()

    # The terminal shows observed consequences, not a second copy of the
    # configuration. The same lines are retained in the filtering log.
    for line in summary_lines:
        print(line)
    print("")
    return {
        "genome_build": genome_build,
        "filtered_vcf": str(output_vcf),
        "filtered_vcf_index": str(Path(str(output_vcf) + ".tbi")),
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
        "filter_reason_summary": reason_report,
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
