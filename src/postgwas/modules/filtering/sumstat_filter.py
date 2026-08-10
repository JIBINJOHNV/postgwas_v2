"""
Filter a harmonised GWAS-VCF with bcftools.

Every threshold that used to be hardcoded at the call site is now a policy key
(see ``harmonisation/policies.py``).  Passing ``policies=None`` - which is what
every existing call site does - reproduces the previous behaviour exactly,
because each historical default below is the value the parameter had before.

Policy keys read by this module::

    filter.lp_cutoff              null   (-log10 p threshold; see `pval_cutoff` below)
    filter.maf_cutoff             0.01
    filter.af_diff_cutoff         0.2
    filter.af_missing             remove | keep
    filter.info_cutoff            0.7
    filter.info_max               1.05
    filter.info_missing           remove | keep
    filter.include_indels         false
    filter.exclude_palindromic    true
    filter.palindromic_af_lower   0.4
    filter.palindromic_af_upper   0.6
    filter.remove_mhc             true
    filter.mhc_chrom              '6'
    filter.mhc_start              25000000
    filter.mhc_end                34000000
    filter.empty_expression       match_all | match_none
    execution.threads_per_chromosome  5
    execution.max_mem_per_job         '12G'
"""

import csv
import io
import math
import os
import re
import shlex
from pathlib import Path
from typing import Any, Dict, List, Optional

import polars as pl

from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.execution.runtime import run_cmd

try:  # pandas is a hard dependency of the pipeline, but this module can live without it
    import pandas as pd
except ImportError:  # pragma: no cover - defensive
    pd = None


# ``_UNSET`` distinguishes "the caller said nothing" from "the caller explicitly
# asked for None/False".  Only then can a policy value fill a parameter in
# without ever overriding something the caller actually passed.
_UNSET = object()

#: Historical hardcoded defaults, kept so that a call with policies=None and no
#: explicit argument behaves exactly as it did before the policy registry existed.
HISTORICAL_DEFAULTS = {
    "pval_cutoff": None,
    "maf_cutoff": None,
    "allelefreq_diff_cutoff": None,
    "info_cutoff": None,
    "info_max": None,
    "info_missing": "keep",
    "include_indels": True,
    "exclude_palindromic": False,
    "palindromic_af_lower": 0.4,
    "palindromic_af_upper": 0.6,
    "remove_mhc": False,
    "mhc_chrom": "6",
    "mhc_start": 25000000,
    "mhc_end": 34000000,
    "threads": 5,
    "max_mem": "5G",
    "af_missing": "remove",
    "lp_missing": "remove",
    "empty_expression": "match_all",
}

MISSING_CHOICES = ("keep", "remove")
EMPTY_EXPRESSION_CHOICES = ("match_all", "match_none")


def _resolve(value, policies, key: Optional[str], historical):
    """Explicit argument wins, then the policy, then the historical hardcoded value."""
    if value is not _UNSET:
        return value
    if policies is None or key is None:
        return historical
    return policies.get(key)


def _is_missing(value) -> bool:
    """True for None, NaN and pandas' NA.

    ``main.py`` puts ``pd.NA`` *inside* ``extra_options``, so ``dict.get(k, 0)``
    never falls back and ``pd.NA < 5`` raises "boolean value of NA is ambiguous"
    - after the filtered VCF has already been written.
    """
    if value is None:
        return True
    if pd is not None:
        try:
            result = pd.isna(value)
        except (TypeError, ValueError):
            return False
        try:
            return bool(result)
        except (TypeError, ValueError):  # array-like
            return False
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
            indent=6, label_width=22,
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
                indent=8, label_width=31,
            ))
    lines.extend([
        screen_field(
            "count", "VCF before filtering", _fmt_count(variants_before),
            indent=8, label_width=31,
        ),
        screen_field(
            "count", "VCF after filtering", _fmt_count(variants_after),
            indent=8, label_width=31,
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
                indent=8, label_width=20,
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
        label_width = max(
            [56] + [len(str(check["label"])) for check in condition_checks]
        )
        for check in condition_checks:
            lines.extend(screen_field(
                "loss", check["label"],
                "matched %s; removed for this reason %s" % (
                    _fmt_count(check.get("count")),
                    _fmt_count(check.get("primary_count")),
                ),
                indent=8, label_width=label_width,
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
            indent=8, label_width=28,
        ),
        screen_field(
            "analysis", "Overlapping rule matches",
            _fmt_count(reason_statistics.get("extra_rule_matches")),
            indent=8, label_width=28,
        ),
    ])
    removed = None
    if variants_before is not None and variants_after is not None:
        removed = variants_before - variants_after
    lines.append(screen_field(
        "loss", "Removed by all filters", _fmt_count(removed),
        indent=8, label_width=28,
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
            indent=8, label_width=28,
        ))
    else:
        lines.append(screen_field(
            "warning", "Removal count check", "unavailable",
            indent=8, label_width=28,
        ))
    if reason_report:
        lines.append(screen_field(
            "info", "Detailed reason report", os.path.basename(reason_report),
            indent=8, label_width=28,
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
) -> Optional[Dict[str, Any]]:
    """Tag all filter reasons in one bcftools stream, then aggregate with Polars."""
    if not checks:
        return {
            "checks": [], "primary_removed_total": 0,
            "overlap_variants": 0, "extra_rule_matches": 0,
        }

    for index, check in enumerate(checks, 1):
        check["tag"] = "PGWAS_FILTER_%02d" % index

    tag_file = Path(output_folder) / (
        ".%s_filter_reason_tags_%d.tsv" % (output_prefix, os.getpid())
    )
    stages = ["bcftools view -Ou %s" % shlex.quote(str(vcf_path))]
    for check in checks:
        stages.append(
            "bcftools filter -Ou -m + -s %s -e %s"
            % (check["tag"], shlex.quote(check["expr"]))
        )
    stages.append("bcftools query -f '%FILTER\\n'")
    inner = "set -euo pipefail; %s > %s" % (
        " | ".join(stages), shlex.quote(str(tag_file)),
    )
    try:
        result = run_cmd("bash -c %s" % shlex.quote(inner))
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
    output_folder: str,
    output_prefix: str,
) -> Optional[str]:
    """Persist the exact reason accounting as a reloadable TSV."""
    if statistics is None:
        return None
    path = Path(output_folder) / "qc_summary" / (
        "%s_filter_reason_summary.tsv" % output_prefix
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
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
    return str(path)


def filter_gwas_vcf_bcftools(
    vcf_path: str,
    output_folder: str,
    output_prefix: str,
    # Every parameter defaulted to `_UNSET` below keeps its historical default
    # value - see HISTORICAL_DEFAULTS - and only takes its value from `policies`
    # when the caller passed nothing at all. An explicitly passed value, even
    # None or False, always wins.
    pval_cutoff: Optional[float] = _UNSET,
    maf_cutoff: Optional[float] = _UNSET,
    allelefreq_diff_cutoff: Optional[float] = _UNSET,
    info_cutoff: Optional[float] = _UNSET,
    info_min: Optional[float] = None,
    info_max: Optional[float] = _UNSET,
    info_missing: str = _UNSET,
    external_af_name: str = "EUR",
    include_indels: bool = _UNSET,
    exclude_palindromic: bool = _UNSET,
    palindromic_af_lower: float = _UNSET,
    palindromic_af_upper: float = _UNSET,
    remove_mhc: bool = _UNSET,
    mhc_chrom: str = _UNSET,
    mhc_start: int = _UNSET,
    mhc_end: int = _UNSET,
    threads: int = _UNSET,
    max_mem: str = _UNSET,
    extra_options: Optional[Dict[str, Any]] = None,
    policies=None,
    logger=None,
    lp_cutoff: Optional[float] = None,
    max_pvalue: Optional[float] = None,
    lp_missing: Optional[str] = None,
    af_missing: Optional[str] = None,
    empty_expression: Optional[str] = None,
    report_missing_counts: bool = True,
) -> Dict[str, str]:
    """Filter a GWAS-VCF and report, honestly, how many variants each rule cost.

    New keyword arguments (all optional, all defaulting to the previous behaviour):

    ``policies``        a ``harmonisation.policies.Policies`` object; fills in any
                        threshold the caller did not pass.
    ``logger``          a ``core.pipeline_logging.PipelineLogger``; every
                        line that goes to the log file also goes to it, and the
                        filtering itself is recorded as one QC action with a
                        before and an after count.
    ``lp_cutoff``       threshold on ``FORMAT/LP``, which is -log10(p).
    ``max_pvalue``      a *raw* p-value, converted to ``lp_cutoff`` with -log10.
    ``lp_missing``      ``keep``/``remove`` for variants whose LP is missing.
    ``af_missing``      ``keep``/``remove`` for variants whose AF is missing.
    ``empty_expression````match_all``/``match_none`` when no include filter is active.

    ``pval_cutoff`` is kept for compatibility and still means "minimum
    ``FORMAT/LP``", i.e. a -log10 p value, NOT a raw p-value.  Passing 5e-8
    therefore filters nothing; a warning now says so.  Use ``max_pvalue=5e-8``.
    """
    # ------------------------------------------------------------------
    # Policy resolution (explicit argument > policy > historical default)
    # ------------------------------------------------------------------
    pval_cutoff = _resolve(pval_cutoff, policies, None, HISTORICAL_DEFAULTS["pval_cutoff"])
    maf_cutoff = _resolve(maf_cutoff, policies, "filter.maf_cutoff", HISTORICAL_DEFAULTS["maf_cutoff"])
    allelefreq_diff_cutoff = _resolve(
        allelefreq_diff_cutoff, policies, "filter.af_diff_cutoff",
        HISTORICAL_DEFAULTS["allelefreq_diff_cutoff"])
    info_cutoff = _resolve(info_cutoff, policies, "filter.info_cutoff", HISTORICAL_DEFAULTS["info_cutoff"])
    info_max = _resolve(info_max, policies, "filter.info_max", HISTORICAL_DEFAULTS["info_max"])
    info_missing = _resolve(info_missing, policies, "filter.info_missing", HISTORICAL_DEFAULTS["info_missing"])
    include_indels = _resolve(include_indels, policies, "filter.include_indels",
                              HISTORICAL_DEFAULTS["include_indels"])
    exclude_palindromic = _resolve(exclude_palindromic, policies, "filter.exclude_palindromic",
                                   HISTORICAL_DEFAULTS["exclude_palindromic"])
    palindromic_af_lower = _resolve(palindromic_af_lower, policies, "filter.palindromic_af_lower",
                                    HISTORICAL_DEFAULTS["palindromic_af_lower"])
    palindromic_af_upper = _resolve(palindromic_af_upper, policies, "filter.palindromic_af_upper",
                                    HISTORICAL_DEFAULTS["palindromic_af_upper"])
    remove_mhc = _resolve(remove_mhc, policies, "filter.remove_mhc", HISTORICAL_DEFAULTS["remove_mhc"])
    mhc_chrom = _resolve(mhc_chrom, policies, "filter.mhc_chrom", HISTORICAL_DEFAULTS["mhc_chrom"])
    mhc_start = _resolve(mhc_start, policies, "filter.mhc_start", HISTORICAL_DEFAULTS["mhc_start"])
    mhc_end = _resolve(mhc_end, policies, "filter.mhc_end", HISTORICAL_DEFAULTS["mhc_end"])
    threads = _resolve(threads, policies, "execution.threads_per_chromosome",
                       HISTORICAL_DEFAULTS["threads"])
    max_mem = _resolve(max_mem, policies, "execution.max_mem_per_job", HISTORICAL_DEFAULTS["max_mem"])

    af_missing = _resolve(
        _UNSET if af_missing is None else af_missing,
        policies, "filter.af_missing", HISTORICAL_DEFAULTS["af_missing"])
    empty_expression = _resolve(
        _UNSET if empty_expression is None else empty_expression,
        policies, "filter.empty_expression", HISTORICAL_DEFAULTS["empty_expression"])
    # NOTE: there is no `filter.lp_missing` key in the frozen policy registry, so
    # this one is a plain keyword argument.  Its default matches today's silent
    # behaviour (bcftools drops a record whose LP is missing).
    lp_missing = HISTORICAL_DEFAULTS["lp_missing"] if lp_missing is None else lp_missing

    if info_missing not in MISSING_CHOICES:
        raise ValueError("❌ info_missing must be either 'keep' or 'remove'")
    if af_missing not in MISSING_CHOICES:
        raise ValueError("❌ af_missing must be either 'keep' or 'remove'")
    if lp_missing not in MISSING_CHOICES:
        raise ValueError("❌ lp_missing must be either 'keep' or 'remove'")
    if empty_expression not in EMPTY_EXPRESSION_CHOICES:
        raise ValueError("❌ empty_expression must be either 'match_all' or 'match_none'")

    vcf_name = os.path.basename(vcf_path)
    genomeversion = "GRCh38" if "GRCh38" in vcf_name else "GRCh37" if "GRCh37" in vcf_name else None
    output_prefix = f"{output_prefix}_{genomeversion}" if genomeversion else output_prefix
    step_dir = Path(output_folder)
    step_dir.mkdir(parents=True, exist_ok=True)
    log_dir = step_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{output_prefix}_filter_gwas_vcf_bcftools.log"
    log_buffer = io.StringIO()

    def log_print(*args):
        text = " ".join(str(a) for a in args)
        log_buffer.write(text + "\n")
        if logger is not None:
            for line in text.splitlines():
                if line.strip():
                    logger.info(line.strip())

    def log_warn(*args):
        text = " ".join(str(a) for a in args)
        log_buffer.write(text + "\n")
        if logger is not None:
            for line in text.splitlines():
                if line.strip():
                    logger.warn(line.strip())

    def write_log():
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

    # ============================================================
    # VALIDATION
    # ============================================================
    if not os.path.exists(vcf_path):
        msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
        log_print(msg)
        write_log()
        raise FileNotFoundError(msg)
    os.makedirs(output_folder, exist_ok=True)
    output_vcf = os.path.join(output_folder, f"{output_prefix}_filtered.vcf.gz")

    quoted_input = shlex.quote(str(vcf_path))

    # ============================================================
    # FAST COUNT
    # ============================================================
    def fast_count(vcf: str) -> Optional[int]:
        """Number of records in an indexed VCF, or None (never a silent 0)."""
        quoted = shlex.quote(str(vcf))
        problems = []
        try:
            res = run_cmd(f"bcftools index -n {quoted}")
            if res.returncode == 0:
                return int(res.stdout.strip())
            problems.append(f"`bcftools index -n` exit {res.returncode}: {res.stderr.strip()}")
        except Exception as exc:  # run_cmd itself blew up
            problems.append(f"`bcftools index -n` raised {type(exc).__name__}: {exc}")

        fallback = "set -o pipefail; bcftools view --threads %d -H %s | wc -l" % (threads, quoted)
        try:
            res = run_cmd("bash -c %s" % shlex.quote(fallback))
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
        inner = "set -o pipefail; bcftools query -i %s -f '.\\n' %s | wc -l" % (
            shlex.quote(expr), shlex.quote(str(vcf)))
        try:
            res = run_cmd("bash -c %s" % shlex.quote(inner))
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

    pre_variants = fast_count(vcf_path)
    if extra_options is not None:
        log_print(f"\n\t\t\t📊 Variants in the input file                        : {_fmt_count(extra_options.get('total_variant_infile'), 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants successfully read by harmonisation module: {_fmt_count(extra_options.get('total_variant_read'), 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants USED for VCF creation                    : {_fmt_count(extra_options.get('total_variant_in_vcf_input'), 'N/A')}")
        log_print(f"\n\t\t\t📊 Variants BEFORE filtering                         : {_fmt_count(pre_variants)}")
        log_print("")

    # ============================================================
    # P-VALUE / LP THRESHOLD
    # `FORMAT/LP` is -log10(p): bigger means MORE significant.
    # ============================================================
    resolved_lp = None
    lp_source = None
    if lp_cutoff is not None:
        resolved_lp, lp_source = float(lp_cutoff), "lp_cutoff"
    elif max_pvalue is not None:
        if not (0 < float(max_pvalue) <= 1):
            raise ValueError("❌ max_pvalue must be a raw p-value in (0, 1], got %r" % (max_pvalue,))
        resolved_lp, lp_source = -math.log10(float(max_pvalue)), "max_pvalue"
        log_print(
            "🔧 max_pvalue %g converted to a -log10 threshold: FORMAT/LP >= %.6f"
            % (float(max_pvalue), resolved_lp)
        )
    elif pval_cutoff is not None:
        resolved_lp, lp_source = float(pval_cutoff), "pval_cutoff"
        if 0 < resolved_lp < 1:
            log_warn(
                "⚠️ pval_cutoff=%r is compared against FORMAT/LP, which is -log10(p), "
                "so this keeps every variant with p <= %g - i.e. it filters nothing. "
                "Use max_pvalue=%r (or lp_cutoff=%.4f) if you meant a raw p-value."
                % (pval_cutoff, 10 ** (-resolved_lp), pval_cutoff, -math.log10(resolved_lp))
            )
    elif policies is not None:
        policy_lp = policies.get("filter.lp_cutoff")
        if policy_lp is not None:
            resolved_lp, lp_source = float(policy_lp), "filter.lp_cutoff"

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
        lp_expr = f"(FORMAT/LP >= {resolved_lp})"
        if lp_missing == "keep":
            lp_expr = f"({lp_expr} | (FORMAT/LP == '.'))"
        missing_checks.append({
            "field": "FORMAT/LP", "expr": "FORMAT/LP == '.'",
            "action": lp_missing, "key": "lp_missing", "priority": 10,
        })
        include_parts.append(lp_expr)
        condition_checks.append({
            "label": "P-value evidence below LP %.6g" % resolved_lp,
            "expr": f"(FORMAT/LP != '.' & FORMAT/LP < {resolved_lp})",
            "priority": 11,
        })

    if maf_cutoff is not None:
        maf_expr = f"(FORMAT/AF >= {maf_cutoff} & FORMAT/AF <= {1 - maf_cutoff})"
        if af_missing == "keep":
            maf_expr = f"({maf_expr} | (FORMAT/AF == '.'))"
        missing_checks.append({
            "field": "FORMAT/AF", "expr": "FORMAT/AF == '.'",
            "action": af_missing, "key": "af_missing", "priority": 20,
        })
        include_parts.append(maf_expr)
        condition_checks.append({
            "label": "Minor allele frequency outside %.6g ≤ AF ≤ %.6g"
            % (maf_cutoff, 1 - maf_cutoff),
            "expr": (
                f"(FORMAT/AF != '.' & "
                f"(FORMAT/AF < {maf_cutoff} | FORMAT/AF > {1 - maf_cutoff}))"
            ),
            "priority": 21,
        })

    info_expr = None
    base_expr = None
    info_failure_expr = None
    if info_cutoff is not None:
        if info_max is not None:
            if info_cutoff > info_max:
                log_print("⚠️ info_cutoff > info_max → ignoring info_max")
                base_expr = f"(FORMAT/SI >= {info_cutoff})"
                info_failure_expr = f"(FORMAT/SI != '.' & FORMAT/SI < {info_cutoff})"
            else:
                base_expr = f"(FORMAT/SI >= {info_cutoff} & FORMAT/SI <= {info_max})"
                info_failure_expr = (
                    f"(FORMAT/SI != '.' & "
                    f"(FORMAT/SI < {info_cutoff} | FORMAT/SI > {info_max}))"
                )
        else:
            base_expr = f"(FORMAT/SI >= {info_cutoff})"
            info_failure_expr = f"(FORMAT/SI != '.' & FORMAT/SI < {info_cutoff})"
        if info_min is not None:
            log_print("⚠️ info_cutoff provided → info_min ignored")
    elif info_min is not None or info_max is not None:
        if info_min is not None and info_max is not None:
            if info_min > info_max:
                raise ValueError("❌ Invalid INFO range: info_min > info_max")
            base_expr = f"(FORMAT/SI >= {info_min} & FORMAT/SI <= {info_max})"
            info_failure_expr = (
                f"(FORMAT/SI != '.' & "
                f"(FORMAT/SI < {info_min} | FORMAT/SI > {info_max}))"
            )
        elif info_min is not None:
            base_expr = f"(FORMAT/SI >= {info_min})"
            info_failure_expr = f"(FORMAT/SI != '.' & FORMAT/SI < {info_min})"
        else:
            base_expr = f"(FORMAT/SI <= {info_max})"
            info_failure_expr = f"(FORMAT/SI != '.' & FORMAT/SI > {info_max})"

    if base_expr is not None:
        if info_missing == "keep":
            # Inside the filter, '.' must be in single quotes
            info_expr = f"({base_expr} | (FORMAT/SI == '.'))"
        else:
            info_expr = f"({base_expr} & (FORMAT/SI != '.'))"
        missing_checks.append({
            "field": "FORMAT/SI", "expr": "FORMAT/SI == '.'",
            "action": info_missing, "key": "info_missing", "priority": 30,
        })

    if info_expr is not None:
        include_parts.append(info_expr)
        if info_cutoff is not None and info_max is not None and info_cutoff <= info_max:
            info_label = "Imputation quality outside %.6g ≤ INFO (SI) ≤ %.6g" % (
                info_cutoff, info_max,
            )
        elif info_cutoff is not None:
            info_label = "Imputation quality below INFO (SI) %.6g" % info_cutoff
        elif info_min is not None and info_max is not None:
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
        af_diff_expr = f"(abs(INFO/AF - INFO/{external_af_name}) <= {allelefreq_diff_cutoff})"
        if af_missing == "keep":
            af_diff_expr = (
                f"({af_diff_expr} | (INFO/AF == '.') | (INFO/{external_af_name} == '.'))"
            )
        missing_checks.extend([{
            "field": "INFO/AF",
            "expr": "INFO/AF == '.'",
            "action": af_missing,
            "key": "af_missing",
            "priority": 40,
        }, {
            "field": f"INFO/{external_af_name}",
            "expr": f"INFO/{external_af_name} == '.'",
            "action": af_missing,
            "key": "af_missing",
            "priority": 41,
        }])
        include_parts.append(af_diff_expr)
        condition_checks.append({
            "label": "External AF absolute difference |AF − %s| > %.6g"
            % (external_af_name, allelefreq_diff_cutoff),
            "expr": (
                f"(INFO/AF != '.' & INFO/{external_af_name} != '.' & "
                f"abs(INFO/AF - INFO/{external_af_name}) > {allelefreq_diff_cutoff})"
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
    "((REF=='A' & ALT=='T') | (REF=='T' & ALT=='A') | "
    "(REF=='C' & ALT=='G') | (REF=='G' & ALT=='C'))"
    )
    exclude_expr = None
    if exclude_palindromic:
        exclude_expr = (
            f"({palindromic_logic} & "
            f"(FORMAT/AF >= {palindromic_af_lower} & FORMAT/AF <= {palindromic_af_upper}))"
        )
        condition_checks.append({
            "label": "Palindromic SNPs with ambiguous frequency (%.6g ≤ AF ≤ %.6g)"
            % (palindromic_af_lower, palindromic_af_upper),
            "expr": (
                "(%s & FORMAT/AF != '.' & FORMAT/AF >= %s & FORMAT/AF <= %s)"
                % (palindromic_logic, palindromic_af_lower, palindromic_af_upper)
            ),
            "priority": 60,
        })
    if not include_indels:
        condition_checks.append({
            "label": "Indels and other non-SNP variants",
            "expr": "(TYPE != 'snp')",
            "priority": 50,
        })
    if remove_mhc:
        mhc_chrom = _match_contig_naming(
            str(mhc_chrom), vcf_path, threads, log_print, log_warn,
        )
        condition_checks.append({
            "label": "Variants in the MHC region (%s:%s-%s)"
            % (mhc_chrom, f"{int(mhc_start):,}", f"{int(mhc_end):,}"),
            "expr": (
                f"(CHROM == '{mhc_chrom}' & "
                f"POS >= {int(mhc_start)} & POS <= {int(mhc_end)})"
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
        vcf_path=vcf_path,
        checks=reason_checks,
        output_folder=output_folder,
        output_prefix=output_prefix,
        log_warn=log_warn,
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
    # speed edit: no --threads on intermediate uncompressed stream step
    cmd1 = "bcftools view -Ou"
    if include_expr is not None:
        cmd1 += f' -i "{include_expr}"'
    if not include_indels:
        cmd1 += " --types snps"
    cmd1 += f" {quoted_input}"
    cmd_parts.append(cmd1)
    # keep -e separate from -i
    if exclude_expr:
        # speed edit: no --threads on intermediate uncompressed stream step
        cmd_parts.append(f'bcftools view -Ou -e "{exclude_expr}"')
    if remove_mhc:
        mhc_bed = os.path.join(output_folder, f"{output_prefix}_mhc_exclude.bed")
        # BED is 0-based half-open; mhc_start is a 1-based inclusive coordinate,
        # so the start field must be one less or the first base of the region
        # survives the exclusion.
        bed_start = max(0, int(mhc_start) - 1)
        with open(mhc_bed, "w") as f:
            f.write(f"{mhc_chrom}\t{bed_start}\t{int(mhc_end)}\n")
        # speed edit: no --threads on intermediate uncompressed stream step
        cmd_parts.append(f"bcftools view -Ou -T ^{shlex.quote(mhc_bed)}")
    # ============================================================
    # SPEED OPTIMIZATION
    # - keep -i and -e separate
    # - skip expensive sort unless explicitly requested
    #   use extra_options={"force_sort": True} if needed
    # ============================================================
    need_sort = extra_options.get("force_sort", False) if extra_options else False
    if need_sort:
        cmd_parts.append(
            f"bcftools sort --temp-dir {shlex.quote(str(output_folder))} --max-mem {shlex.quote(str(max_mem))}"
        )
    # keep threads here where they help most: final bgzip compression
    cmd_parts.append(f"bcftools view -Oz --threads {threads}")
    pipeline_core = " | ".join(cmd_parts)

    # Without pipefail only the exit status of the LAST stage is seen, so an
    # abort in the middle of the pipe produced a silently TRUNCATED VCF whose
    # reduced count was then reported as "variants removed by filtering".
    inner = (
        f"set -euo pipefail; {pipeline_core} > {shlex.quote(output_vcf)} "
        f"&& tabix -f -p vcf {shlex.quote(output_vcf)}"
    )
    pipeline = "bash -c %s" % shlex.quote(inner)
    log_print("🚀 Full bcftools pipeline:")
    log_print(inner)
    log_print("")
    result = run_cmd(pipeline)
    if result.returncode != 0:
        log_print("❌ Error during bcftools pipeline execution.")
        log_print("STDERR:")
        log_print(result.stderr)
        write_log()
        raise RuntimeError("bcftools pipeline failed")
    # ============================================================
    # Count variants after filtering
    # ============================================================
    post_variants = fast_count(output_vcf)
    log_print("📊 Variants counted with: bcftools index -n (fallback: bcftools view -H | wc -l)")
    log_print(f"✅ Variants AFTER filtering                          : {_fmt_count(post_variants)}")
    log_print(f"💾 Filtered VCF saved to                             : {output_vcf}")
    log_print("\n🎉 bcftools filtering completed.")

    # Rule 5: a QC action is never recorded without a before and an after count.
    if logger is not None and pre_variants is not None and post_variants is not None:
        logger.qc(
            "vcf filtering",
            "bcftools removed the variants that failed the active quality filters.",
            pre_variants,
            post_variants,
            step="04 filter_vcf",
            warn=post_variants == 0 and pre_variants > 0,
        )

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
    reason_report = _write_filter_reason_report(
        reason_statistics, output_folder, output_prefix,
    )
    if reason_report:
        log_print("Detailed filter reason report path: %s" % reason_report)

    # The terminal shows observed consequences, not a second copy of the
    # configuration. The same lines are retained in the filtering log.
    for line in _filtering_summary_lines(
        missing_checks if report_missing_counts else [],
        condition_checks,
        pre_variants,
        post_variants,
        reason_statistics=reason_statistics,
        reason_report=reason_report,
        data_flow=extra_options,
    ):
        print(line)
        log_print(line)
    print("")
    # --- Safe access for extra_options: main.py stores pd.NA IN the dict, so
    # --- dict.get(key, 0) never falls back and every comparison below has to be
    # --- guarded, or `pd.NA < 5` raises after the filtered VCF was written.
    if extra_options is not None:
        infile = extra_options.get("total_variant_infile")
        read = extra_options.get("total_variant_read")
        vcf_input = extra_options.get("total_variant_in_vcf_input")
        # 1️⃣ Input → Read
        if not _is_missing(read) and not _is_missing(infile) and read < infile:
            diff = infile - read
            pct = (diff / infile * 100) if infile else 0
            msg = "Harmonisation did not read all variants: %s missing (%.2f%%)." % (
                _fmt_count(diff), pct,
            )
            log_print(msg)
            print(screen_field(
                "warning", "Input loss", msg,
                indent=4, label_width=18,
            ))
        # 2️⃣ Read → VCF input
        if not _is_missing(vcf_input) and not _is_missing(read) and vcf_input < read:
            diff = read - vcf_input
            pct = (diff / read * 100) if read else 0
            message = "Fewer variants were used for VCF creation: %s removed after harmonisation (%.2f%%)." % (
                _fmt_count(diff), pct,
            )
            log_print(message)
            print("\n" + "\n".join([
                screen_line("warning", "Variant loss during harmonisation", indent=4),
                screen_field(
                    "count", "Summary-stat rows", _fmt_count(read),
                    indent=6, label_width=24,
                ),
                screen_field(
                    "count", "Used for VCF", _fmt_count(vcf_input),
                    indent=6, label_width=24,
                ),
                screen_field(
                    "loss", "Removed", "%s (%.2f%%)" % (_fmt_count(diff), pct),
                    indent=6, label_width=24,
                ),
                screen_field(
                    "info", "Detailed reasons",
                    "see the harmonisation reject-reason report and chromosome logs",
                    indent=6, label_width=24,
                ),
            ]))
        # 3️⃣ VCF input → Pre-filter
        if (pre_variants is not None and not _is_missing(vcf_input)
                and pre_variants < vcf_input):
            diff = vcf_input - pre_variants
            pct = (diff / vcf_input * 100) if vcf_input else 0
            log_print(
                "Not all variants reached the VCF: %s dropped during VCF generation (%.2f%%)."
                % (_fmt_count(diff), pct)
            )
    print("")
    write_log()
    return {
        "filtered_vcf": output_vcf,
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


def _match_contig_naming(mhc_chrom: str, vcf_path: str, threads: int, log_print, log_warn) -> str:
    """Return the MHC contig name as the VCF actually spells it.

    ``bcftools view -T ^bed`` silently excludes nothing when the BED names a
    contig the VCF does not contain, so `6` against a `chr6` VCF (or the
    reverse) turns MHC removal into a no-op with no error anywhere.
    """
    contigs = _vcf_contigs(vcf_path, threads)
    if not contigs:
        log_warn(
            "⚠️ Could not read the contig names from %s, so the MHC contig name "
            "'%s' could not be checked against the VCF." % (vcf_path, mhc_chrom)
        )
        return mhc_chrom
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

    log_warn(
        "⚠️ MHC contig '%s' is not present in %s (contigs seen: %s%s). The MHC "
        "exclusion will remove NOTHING. Set filter.mhc_chrom to match the VCF."
        % (mhc_chrom, vcf_path, ", ".join(contigs[:10]),
           ", ..." if len(contigs) > 10 else "")
    )
    return mhc_chrom


def _vcf_contigs(vcf_path: str, threads: int = 1) -> List[str]:
    """Contig names in a VCF: the ones with data if indexed, else the header's."""
    quoted = shlex.quote(str(vcf_path))
    names: List[str] = []
    try:
        res = run_cmd(f"bcftools index -s {quoted}")
        if res.returncode == 0 and (res.stdout or "").strip():
            for line in res.stdout.splitlines():
                parts = line.split("\t")
                if parts and parts[0].strip():
                    names.append(parts[0].strip())
    except Exception:
        names = []
    if names:
        return names
    try:
        res = run_cmd(f"bcftools view -h {quoted}")
        if res.returncode == 0:
            names = re.findall(r"##contig=<ID=([^,>]+)", res.stdout or "")
    except Exception:
        names = []
    return names
