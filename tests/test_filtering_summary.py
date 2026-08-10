"""User-facing filtering summary and exact-reason accounting contracts."""

import csv

from postgwas.modules.filtering.sumstat_filter import (
    _filtering_summary_lines,
    _summarize_filter_tags,
    _write_filter_reason_report,
)


def test_filtering_summary_separates_missing_values_from_failed_conditions():
    lines = _filtering_summary_lines(
        [
            {
                "field": "FORMAT/AF", "count": 0,
                "primary_count": 0, "action": "remove",
            },
            {
                "field": "FORMAT/SI", "count": 4,
                "primary_count": 4, "action": "remove",
            },
            {
                "field": "INFO/EUR", "count": 1028,
                "primary_count": 1028, "action": "remove",
            },
        ],
        [
            {
                "label": "Minor allele frequency outside 0.01 ≤ AF ≤ 0.99",
                "count": 250, "primary_count": 250,
            },
            {
                "label": "Imputation quality outside 0.7 ≤ INFO (SI) ≤ 1.05",
                "count": 1200, "primary_count": 1196,
            },
            {
                "label": "External AF absolute difference |AF − EUR| > 0.2",
                "count": 800, "primary_count": 788,
            },
            {
                "label": "Palindromic SNPs with ambiguous frequency (0.4 ≤ AF ≤ 0.6)",
                "count": 410, "primary_count": 400,
            },
            {
                "label": "Variants in the MHC region (6:25,000,000-34,000,000)",
                "count": 600, "primary_count": 594,
            },
            {
                "label": "Indels and other non-SNP variants",
                "count": 90, "primary_count": 88,
            },
        ],
        variants_before=999999,
        variants_after=979918,
        reason_statistics={
            "overlap_variants": 34,
            "extra_rule_matches": 38,
            "primary_removed_total": 20081,
        },
        reason_report="/analysis/qc_summary/study_filter_reason_summary.tsv",
        data_flow={
            "total_variant_infile": 1000500,
            "total_variant_read": 1000000,
            "total_variant_removed_null_coords": 1,
            "total_variant_removed_non_standard_alleles": 0,
            "total_variant_remaining_for_harmonisation": 999999,
            "total_variant_with_missing_eaf": 0,
            "total_variant_with_invalid_beta_se": 0,
            "total_variant_in_vcf_input": 999999,
        },
    )
    text = "\n".join(lines)

    assert "🔬  Variant filtering summary" in text
    assert "🧮  Variant flow" in text
    assert "Input summary statistics" in text and "1,000,500" in text
    assert "Used for VCF creation" in text and "999,999" in text
    assert "🧬  Missing values used by active filters" in text
    assert "FORMAT/AF" in text and "missing 0; removed for this reason 0" in text
    assert "FORMAT/SI" in text and "missing 4; removed for this reason 4" in text
    assert "INFO/EUR" in text and "missing 1,028; removed for this reason 1,028" in text
    assert "🔬  Active filtering conditions" in text
    assert "matched 800; removed for this reason 788" in text
    assert "Multiple failed rules" in text and "34" in text
    assert "Overlapping rule matches" in text and "38" in text
    assert "VCF before filtering" in text and "999,999" in text
    assert "VCF after filtering" in text and "979,918" in text
    assert "Removed by all filters" in text and "20,081" in text
    assert "20,081 assigned to reasons = 20,081 removed in total" in text
    assert "study_filter_reason_summary.tsv" in text
    assert "/analysis/qc_summary" not in text


def test_filtering_summary_omits_inactive_filters():
    text = "\n".join(
        _filtering_summary_lines([], [], variants_before=10, variants_after=10)
    )

    assert "No variant-level filtering condition is active" in text
    assert "MHC" not in text
    assert "Palindromic" not in text


def test_filter_tag_summary_separates_observed_overlap_from_primary_reason(tmp_path):
    tag_file = tmp_path / "filter_tags.tsv"
    tag_file.write_text("PASS\nR1\nR1;R2\nR2;R3\n", encoding="utf-8")
    checks = [
        {"tag": "R1", "label": "Reason one", "removes": True},
        {"tag": "R2", "label": "Reason two", "removes": True},
        {"tag": "R3", "label": "Observed only", "removes": False},
    ]

    summary = _summarize_filter_tags(tag_file, checks)
    by_tag = {row["tag"]: row for row in summary["checks"]}

    assert by_tag["R1"]["count"] == 2
    assert by_tag["R1"]["primary_count"] == 2
    assert by_tag["R2"]["count"] == 2
    assert by_tag["R2"]["primary_count"] == 1
    assert by_tag["R3"]["count"] == 1
    assert by_tag["R3"]["primary_count"] == 0
    assert summary["primary_removed_total"] == 3
    assert summary["overlap_variants"] == 1
    assert summary["extra_rule_matches"] == 1


def test_filter_reason_report_is_reloadable_and_records_reconciliation(tmp_path):
    statistics = {
        "checks": [{
            "tag": "R1", "category": "filter_condition", "label": "Reason one",
            "action": "remove", "count": 7, "primary_count": 6,
            "expr": "FORMAT/AF < 0.01",
        }],
        "overlap_variants": 1,
        "extra_rule_matches": 1,
        "primary_removed_total": 6,
        "actual_removed": 6,
        "reconciled": True,
    }

    report = _write_filter_reason_report(statistics, str(tmp_path), "study")

    with open(report, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["reason"] == "Reason one"
    assert rows[0]["variants_matching_reason"] == "7"
    assert rows[0]["variants_removed_for_this_reason"] == "6"
    assert rows[-1]["reason"] == "Reason reconciliation"
    assert rows[-1]["status"] == "PASS"
