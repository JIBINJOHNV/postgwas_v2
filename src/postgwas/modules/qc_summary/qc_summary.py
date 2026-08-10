"""
QC summary for a harmonised GWAS-VCF, built on ``bcftools stats``.

Policy keys read by this module::

    filter.af_diff_cutoff             0.2   (backs `allelefreq_diff_cutoff`)
    execution.threads_per_chromosome  5     (backs `threads`)
    qc.report_singleton_stats         true  (recommended: false - see below)

NOTE: the plan asks for ``qc.af_diff_cutoff`` and ``qc.threads``.  Neither key
exists in the frozen policy registry, so the two equivalent registered keys above
are used instead; both carry the same default as the function argument they back,
so nothing changes for a caller that passes neither.
"""

import pandas as pd
import re
import subprocess
from typing import Optional, Tuple

from postgwas.core.processes import run_checked_command


# ``_UNSET`` distinguishes "the caller said nothing" from an explicit value, so a
# policy can fill a parameter in without ever overriding the caller.
_UNSET = object()

HISTORICAL_DEFAULTS = {
    "allelefreq_diff_cutoff": 0.2,
    "threads": 5,
    "report_singleton_stats": True,
}


class MissingVcfTagError(RuntimeError):
    """A tag the QC metrics are computed from is not declared in the VCF header."""


class BcftoolsCommandError(RuntimeError):
    """A bcftools command used for a QC metric failed."""


def _resolve(value, policies, key: Optional[str], historical):
    """Explicit argument wins, then the policy, then the historical hardcoded value."""
    if value is not _UNSET:
        return value
    if policies is None or key is None:
        return historical
    return policies.get(key)


def vcf_header_tags(vcf_path: str, bcftools_bin: str = "bcftools") -> Tuple[set, set]:
    """Return ``(info_ids, format_ids)`` declared in a VCF header."""
    header = run_checked_command(
        [bcftools_bin, "view", "-h", str(vcf_path)],
        "Reading the VCF header of %s" % vcf_path,
        error_type=BcftoolsCommandError,
    )
    return (
        set(re.findall(r"##INFO=<ID=([^,>]+)", header)),
        set(re.findall(r"##FORMAT=<ID=([^,>]+)", header)),
    )


def count_matching(vcf_path: str, expression: Optional[str], threads: int = 5,
                   bcftools_bin: str = "bcftools") -> int:
    """Count the records matching a bcftools expression.

    ``bcftools view ... | wc -l`` under ``subprocess.check_output`` only ever sees
    wc's exit status, which is 0 even when bcftools aborted on the first record -
    so a wrong tag name produced the metric 0 and a clean bill of health on a file
    that was never actually checked.  ``set -o pipefail`` plus a checked exit
    status is what makes the number trustworthy.
    """
    inner = "set -o pipefail; %s view --threads %d --no-header" % (bcftools_bin, threads)
    if expression:
        inner += " -i '%s'" % expression.replace("'", "'\\''")
    inner += " '%s' | wc -l" % str(vcf_path).replace("'", "'\\''")
    out = run_checked_command(
        ["bash", "-c", inner],
        "Counting records where %s" % (expression or "TRUE"),
        error_type=BcftoolsCommandError,
    )
    return int(out.strip())


def run_bcftools_stats(
    vcf_path: str,
    external_af_name: str = "EUR",
    allelefreq_diff_cutoff: float = _UNSET,
    threads: int = _UNSET,
    bcftools_bin: str = "bcftools",
    policies=None,
    logger=None,
    on_missing_tag: str = "fail",
):
    """Run the three bcftools operations behind the AF-concordance QC metrics.

    The three allele-frequency categories are now **disjoint**, so they can be
    added up and reported as percentages of one denominator:

    ``external_af_missing``          the external panel has no frequency
    ``study_af_missing``             the panel has one, the study does not
    ``af_comparable``                both frequencies are present
    ``af_difference_above_cutoff``   comparable frequencies differ too much

    Previously the discordance count used ``-e`` on an expression that is *undefined*
    when either frequency is missing; ``-e`` only excludes where the expression is
    TRUE, so undefined records survived and were counted. The old discordance metric was
    therefore the union of all three categories and the missing-external count was a
    strict subset of it - two overlapping counts printed as separate percentages.

    ``on_missing_tag`` is ``fail`` (raise ``MissingVcfTagError``) or ``warn``
    (report the metrics as ``None``).
    """
    allelefreq_diff_cutoff = _resolve(
        allelefreq_diff_cutoff, policies, "filter.af_diff_cutoff",
        HISTORICAL_DEFAULTS["allelefreq_diff_cutoff"])
    threads = _resolve(threads, policies, "execution.threads_per_chromosome",
                       HISTORICAL_DEFAULTS["threads"])
    if on_missing_tag not in ("fail", "warn"):
        raise ValueError("on_missing_tag must be 'fail' or 'warn', got %r" % (on_missing_tag,))

    results = {
        "vcf_path": str(vcf_path),
        "external_af_name": external_af_name,
        "af_diff_cutoff": allelefreq_diff_cutoff,
    }

    def note(message, level="info"):
        if logger is not None:
            if level == "warn":
                logger.warn(message)
            else:
                logger.info(message)

    # ----------------------------------------------------------------------
    # 0️⃣ The tags the metrics are computed from must exist FIRST
    # ----------------------------------------------------------------------
    info_ids, format_ids = vcf_header_tags(vcf_path, bcftools_bin=bcftools_bin)
    missing_tags = []
    if external_af_name not in info_ids:
        missing_tags.append("INFO/%s" % external_af_name)
    if "AF" not in format_ids:
        missing_tags.append("FORMAT/AF")

    if missing_tags:
        message = (
            "%s is/are not declared in the header of %s, so the allele-frequency "
            "concordance metrics cannot be computed. Tags present: INFO=[%s] FORMAT=[%s]. "
            "Check external_af_name."
            % (", ".join(missing_tags), vcf_path,
               ", ".join(sorted(info_ids)), ", ".join(sorted(format_ids)))
        )
        if on_missing_tag == "fail":
            note(message, level="warn")
            raise MissingVcfTagError(message)
        note(message, level="warn")
        results.update({
            "total_records": None,
            "external_af_missing": None,
            "study_af_missing": None,
            "af_comparable": None,
            "af_difference_above_cutoff": None,
        })
    else:
        ext = external_af_name
        results["total_records"] = count_matching(
            vcf_path, None, threads=threads, bcftools_bin=bcftools_bin)

        # 1️⃣ external frequency absent
        results["external_af_missing"] = count_matching(
            vcf_path, 'INFO/%s=="."' % ext, threads=threads, bcftools_bin=bcftools_bin)

        # 2️⃣ external present, study frequency absent
        results["study_af_missing"] = count_matching(
            vcf_path, 'INFO/%s!="." && FORMAT/AF=="."' % ext,
            threads=threads, bcftools_bin=bcftools_bin)

        # 3️⃣ both present and genuinely discordant
        results["af_difference_above_cutoff"] = count_matching(
            vcf_path,
            'INFO/%s!="." && FORMAT/AF!="." && abs(FORMAT/AF - INFO/%s) > %s'
            % (ext, ext, allelefreq_diff_cutoff),
            threads=threads, bcftools_bin=bcftools_bin)

        total = results["total_records"]
        comparable = None
        if total is not None:
            comparable = total - results["external_af_missing"] - results["study_af_missing"]
        results["af_comparable"] = comparable
        note(
            "Allele-frequency concordance against INFO/%s at cutoff %s: %s comparable, "
            "%s discordant, %s missing external AF, %s missing study AF (of %s records)."
            % (ext, allelefreq_diff_cutoff, comparable,
               results["af_difference_above_cutoff"],
               results["external_af_missing"], results["study_af_missing"], total)
        )

    # ----------------------------------------------------------------------
    # 4️⃣ Run bcftools stats
    # ----------------------------------------------------------------------
    stats_file = f"{vcf_path}.stats"
    with open(stats_file, "w") as handle:
        result = subprocess.run(
            [bcftools_bin, "stats", str(vcf_path)],
            stdout=handle, stderr=subprocess.PIPE, text=True,
        )
    if result.returncode != 0:
        raise BcftoolsCommandError(
            "`%s stats %s` failed (exit status %d). bcftools said: %s"
            % (bcftools_bin, vcf_path, result.returncode, (result.stderr or "").strip())
        )
    results["stats_file"] = stats_file

    return results


def _to_numeric_column(series: pd.Series) -> pd.Series:
    """Coerce a column to numbers, keeping it as text when it is not numeric at all.

    ``pd.to_numeric(errors="ignore")`` warns on pandas 2.3, is removed in pandas 3,
    and is all-or-nothing per column: one unparseable cell left the whole column as
    ``object``, and ``print_qc_summary``'s ``f"{val:.2f}"`` then raised ValueError.
    ``errors="coerce"`` fixes that, but would also blank a genuinely textual column
    such as the substitution ``type`` (``A>C``), so a column that coerces to *all*
    NaN keeps its original text.
    """
    coerced = pd.to_numeric(series, errors="coerce")
    if coerced.notna().any() or series.isna().all():
        return coerced
    return series


def parse_bcftools_stats(stats_file: str):
    """Parse bcftools stats output into multiple pandas DataFrames."""
    sections = {
        "SN": [],
        "TSTV": [],
        "SiS": [],
        "AF": [],
        "QUAL": [],
        "IDD": [],
        "ST": [],
        "DP": []
    }

    with open(stats_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            tag = parts[0]
            if tag in sections:
                sections[tag].append(parts[1:])

    dfs = {}

    # SN
    if sections["SN"]:
        df_sn = pd.DataFrame(sections["SN"], columns=["id", "key", "value"])
        df_sn["value"] = _to_numeric_column(df_sn["value"])
        dfs["SN"] = df_sn

    # All other sections (USE EXACT HEADERS FROM YOUR ORIGINAL CODE)
    header_map = {
        "TSTV": ["id", "ts", "tv", "ts/tv", "ts_1st", "tv_1st", "ts/tv_1st"],
        "SiS": ["id", "allele_count", "num_snps", "num_transitions",
                "num_transversions", "num_indels",
                "repeat_consistent", "repeat_inconsistent", "not_applicable"],
        "AF": ["id", "allele_freq", "num_snps", "num_ts", "num_tv",
               "num_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"],
        "QUAL": ["id", "quality", "num_snps", "ts_1st", "tv_1st", "num_indels"],
        "IDD": ["id", "length", "num_sites", "num_genotypes", "mean_vaf"],
        "ST": ["id", "type", "count"],
        "DP": ["id", "depth_bin", "num_genotypes", "frac_genotypes",
               "num_sites", "frac_sites"]
    }

    for sec, header in header_map.items():
        rows = sections.get(sec, [])
        if rows:
            df = pd.DataFrame(rows, columns=header)
            for c in header:
                df[c] = _to_numeric_column(df[c])
            dfs[sec] = df

    return dfs


def _select_row(df: pd.DataFrame, column: str, value, logger=None, what: str = ""):
    """Return the row whose `column` equals `value`, or the first row as a fallback."""
    if column in df.columns:
        match = df[df[column] == value]
        if not match.empty:
            return match.iloc[0]
    if logger is not None:
        logger.warn(
            "bcftools stats has no %s row with %s == %r; falling back to the first row, "
            "which may not be the bin these metrics claim to describe."
            % (what or "matching", column, value)
        )
    return df.iloc[0]


def bcftools_essential_summary(dfs: dict, policies=None, report_singleton_stats=_UNSET,
                               logger=None) -> pd.DataFrame:
    """Combine key bcftools-stats metrics into a single flat DataFrame.

    ``qc.report_singleton_stats`` (default ``True`` == the previous behaviour,
    recommended ``False``) controls the five ``singleton_*`` columns.  A GWAS-VCF
    carries no genotypes, so bcftools files every variant into the singleton and
    AF=0 bins and the values describe nothing about the study; ``singleton_count``
    in particular is the allele-count BIN label, which is always 1, not a count.
    """
    report_singleton_stats = _resolve(
        report_singleton_stats, policies, "qc.report_singleton_stats",
        HISTORICAL_DEFAULTS["report_singleton_stats"])

    summary = {}

    # SN
    if "SN" in dfs:
        sn = dfs["SN"].set_index("key")["value"]
        summary.update({
            "num_samples": sn.get("number of samples:", None),
            "num_records": sn.get("number of records:", None),
            "num_snps": sn.get("number of SNPs:", None),
            "num_indels": sn.get("number of indels:", None),
            "num_mnps": sn.get("number of MNPs:", None),
            "num_others": sn.get("number of others:", None),
            "num_multiallelic": sn.get("number of multiallelic sites:", None),
            "num_multiallelic_snp": sn.get("number of multiallelic SNP sites:", None),
        })

    # TSTV
    if "TSTV" in dfs:
        tstv = dfs["TSTV"].iloc[0]
        summary.update({
            "ts": tstv["ts"],
            "tv": tstv["tv"],
            "ts_tv_ratio": tstv["ts/tv"]
        })

    # SiS — singleton bin. Selected by its allele_count label rather than by
    # position, so the column names describe the row actually taken.
    if "SiS" in dfs and report_singleton_stats:
        sis = _select_row(dfs["SiS"], "allele_count", 1, logger=logger, what="singleton (AC=1)")
        summary.update({
            "singleton_count": sis["allele_count"],
            "singleton_snps": sis["num_snps"],
            "singleton_ts": sis["num_transitions"],
            "singleton_tv": sis["num_transversions"],
            "singleton_indels": sis["num_indels"]
        })

    # AF — the AF=0 bin, selected explicitly instead of by position.
    if "AF" in dfs:
        af = _select_row(dfs["AF"], "allele_freq", 0, logger=logger, what="AF=0")
        summary.update({
            "af_0_snps": af["num_snps"],
            "af_0_ts": af["num_ts"],
            "af_0_tv": af["num_tv"],
            "af_0_indels": af["num_indels"]
        })

    return pd.DataFrame([summary])
