"""User-facing filtering summary and exact-reason accounting contracts."""

import csv
from unittest.mock import patch

from rich.cells import cell_len

from postgwas.core.dataframes import collect_streaming
from postgwas.core.ui.screen import SYMBOLS
from postgwas.modules.filtering.reporting import (
    SUMMARY_COLUMNS,
    build_filtering_summary,
    render_filtering_html,
    render_filtering_plan,
    render_input_vcf_validation,
    write_filtering_summary_csv,
)
from postgwas.modules.filtering.sumstat_filter import (
    _filtering_summary_lines,
    _summarize_filter_tags,
    _write_filter_reason_report,
)


def test_input_validation_card_is_concise_and_records_field_provenance():
    text = render_input_vcf_validation(
        {
            "status": "PASS",
            "input_vcf": "/analysis/study_GRCh37_merged.vcf.gz",
            "genome_build": "GRCh37",
            "declared_contigs": 23,
            "total_variants": 999999,
            "required_fields_evaluated": True,
            "required_fields": [
                {
                    "field": "FORMAT/AF",
                    "declared": True,
                    "required_by": ["minor-allele-frequency filter"],
                },
                {
                    "field": "FORMAT/SI",
                    "declared": True,
                    "required_by": ["imputation-quality filter"],
                },
            ],
            "field_provenance": {
                "status": "AVAILABLE",
                "detail": "PostGWAS 2.0; VCF status complete",
            },
            "field_summaries": [
                {
                    "kind": "genetic",
                    "label": "Allele-frequency fields",
                    "value": "FORMAT/AF for MAF",
                },
                {
                    "kind": "analysis",
                    "label": "Imputation quality score",
                    "value": "FORMAT/SI · external user-provided reference",
                },
            ],
        },
        label_width=42,
    )

    assert "✅  Input VCF validation" in text
    assert "study_GRCh37_merged.vcf.gz" in text
    assert "Genome build" in text and "GRCh37" in text
    assert "Declared contigs" in text and "23" in text
    assert "Total variants" in text and "999,999" in text
    assert "PASSED · 2/2 required INFO/FORMAT fields declared" in text
    assert "Field provenance" in text and "AVAILABLE" in text
    assert "PostGWAS 2.0; VCF status complete" in text
    assert "Allele-frequency fields" in text and "FORMAT/AF for MAF" in text
    assert "FORMAT/SI · external user-provided reference" in text
    assert "required by minor-allele-frequency filter" not in text


def test_filtering_plan_preserves_the_first_failure_audit_order():
    text = render_filtering_plan([
        {
            "label": "Missing FORMAT/AF",
            "action": "remove",
            "removes": True,
        },
        {
            "label": "Minor allele frequency outside range",
            "action": "remove",
            "removes": True,
        },
        {
            "label": "Missing FORMAT/SI",
            "action": "keep",
            "removes": False,
        },
    ])

    assert "Filtering plan · applied in this order" in text
    assert text.index("1. EXCLUDE · Missing FORMAT/AF") < text.index(
        "2. EXCLUDE · Minor allele frequency outside range"
    )
    assert text.index(
        "2. EXCLUDE · Minor allele frequency outside range"
    ) < text.index("3. KEEP · Missing FORMAT/SI")
    plan_rows = [
        line for line in text.splitlines()
        if "EXCLUDE ·" in line or "KEEP ·" in line
    ]
    assert all(line.startswith("        ") for line in plan_rows)
    assert all(
        not line.lstrip().startswith(tuple(SYMBOLS.values()))
        for line in plan_rows
    )
    assert "REMOVE ·" not in text
    assert "recorded removal reason" not in text


def test_filtering_summary_merges_missing_values_with_filter_conditions():
    af_group = {
        "display_group": "Allele frequency · FORMAT/AF",
        "display_group_kind": "genetic",
    }
    info_group = {
        "display_group": "Imputation quality · FORMAT/SI",
        "display_group_kind": "analysis",
    }
    concordance_group = {
        "display_group": (
            "Study/reference frequency concordance · INFO/AF and INFO/EUR"
        ),
        "display_group_kind": "genetic",
    }
    lines = _filtering_summary_lines(
        [
            {
                "field": "FORMAT/AF", "count": 0,
                "primary_count": 0, "action": "remove", "priority": 20,
                **af_group,
            },
            {
                "field": "FORMAT/SI", "count": 4,
                "primary_count": 4, "action": "remove", "priority": 30,
                **info_group,
            },
            {
                "field": "INFO/AF", "count": 0,
                "primary_count": 0, "action": "remove", "priority": 40,
                **concordance_group,
            },
            {
                "field": "INFO/EUR", "count": 1028,
                "primary_count": 1028, "action": "remove", "priority": 41,
                **concordance_group,
            },
        ],
        [
            {
                "label": "Minor allele frequency outside 0.01 ≤ AF ≤ 0.99",
                "count": 250, "primary_count": 250, "priority": 21,
                **af_group,
            },
            {
                "label": "Imputation quality outside 0.7 ≤ INFO (SI) ≤ 1.05",
                "count": 1200, "primary_count": 1196, "priority": 31,
                **info_group,
            },
            {
                "label": "External AF absolute difference |AF − EUR| > 0.2",
                "count": 800, "primary_count": 788, "priority": 42,
                **concordance_group,
            },
            {
                "label": "Palindromic SNPs with ambiguous frequency (0.4 ≤ AF ≤ 0.6)",
                "count": 410, "primary_count": 400, "priority": 60,
                "display_group": (
                    "Palindromic allele ambiguity · REF/ALT and FORMAT/AF"
                ),
                "display_group_kind": "genetic",
            },
            {
                "label": "Variants in the MHC region (6:25,000,000-34,000,000)",
                "count": 600, "primary_count": 594, "priority": 70,
                "display_group": "Genomic region · CHROM/POS",
                "display_group_kind": "genetic",
            },
            {
                "label": "Indels and other non-SNP variants",
                "count": 90, "primary_count": 88, "priority": 50,
                "display_group": "Variant type · TYPE",
                "display_group_kind": "genetic",
            },
        ],
        variants_before=999999,
        variants_after=979918,
        terminal_label_width=42,
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

    assert "Input VCF validation" not in text
    assert "🔬  Variant filtering summary" in text
    assert "🧮  Upstream harmonisation flow" in text
    assert "Input summary statistics" in text and "1,000,500" in text
    assert "Used for VCF creation" in text and "999,999" in text
    assert "Missing values used by active filters" not in text
    assert "🔬  Filtering results by rule" in text
    assert (
        "🧬  Allele frequency · FORMAT/AF\n\n"
        "            🔻  Missing FORMAT/AF\n"
        "                🔹  0 variants failed this rule and were removed:\n"
        "                    • 0 failed no earlier EXCLUDE rule\n"
        "                    • 0 also failed one or more earlier EXCLUDE rules"
    ) in text
    assert (
        "🔬  Imputation quality · FORMAT/SI\n\n"
        "            🔻  Missing FORMAT/SI\n"
        "                🔹  4 variants failed this rule and were removed:\n"
        "                    • 4 failed no earlier EXCLUDE rule\n"
        "                    • 0 also failed one or more earlier EXCLUDE rules"
    ) in text
    assert (
        "Missing INFO/EUR\n"
        "                🔹  1,028 variants failed this rule and were removed:\n"
        "                    • 1,028 failed no earlier EXCLUDE rule\n"
        "                    • 0 also failed one or more earlier EXCLUDE rules"
    ) in text
    assert (
        "800 variants failed this rule and were removed:\n"
        "                    • 788 failed no earlier EXCLUDE rule\n"
        "                    • 12 also failed one or more earlier EXCLUDE rules"
    ) in text
    assert (
        "Variants in the MHC region (6:25,000,000-34,000,000)\n"
        "                🔹  600 variants failed this rule and were removed:\n"
        "                    • 594 failed no earlier EXCLUDE rule\n"
        "                    • 6 also failed one or more earlier EXCLUDE rules"
    ) in text
    groups = (
        "Allele frequency · FORMAT/AF",
        "Imputation quality · FORMAT/SI",
        "Study/reference frequency concordance · INFO/AF and INFO/EUR",
        "Variant type · TYPE",
        "Palindromic allele ambiguity · REF/ALT and FORMAT/AF",
        "Genomic region · CHROM/POS",
    )
    assert all(text.count(group) == 1 for group in groups)
    assert [text.index(group) for group in groups] == sorted(
        text.index(group) for group in groups
    )
    assert "Variants failing multiple rules" in text and "34" in text
    assert "Rule failures beyond the first" in text and "38" in text
    assert "VCF before filtering" in text and "999,999" in text
    assert "VCF after filtering" in text and "979,918" in text
    assert "Removed by all filters" in text and "20,081" in text
    assert "20,081 assigned to reasons = 20,081 removed in total" in text
    assert "study_filter_reason_summary.tsv" in text
    assert "/analysis/qc_summary" not in text

    aligned_lines = [
        line
        for line in text.splitlines()
        if " : " in line
    ]
    assert len(aligned_lines) > 10
    assert len({cell_len(line.split(" : ", 1)[0]) for line in aligned_lines}) == 1
    assert "\n\n      🧮  Upstream harmonisation flow\n\n" in text
    assert "\n\n      🔬  Filtering results by rule\n\n" in text
    assert "\n\n      🧮  Exact removal accounting\n\n" in text
    assert "\n\n      🧮  Final filtering outcome\n\n" in text
    assert text.index("Saved reports") < text.index("Final filtering outcome")
    assert text.rstrip().endswith("97.99%")


def test_filtering_summary_omits_inactive_filters():
    text = "\n".join(
        _filtering_summary_lines(
            [], [], variants_before=10, variants_after=10,
            terminal_label_width=42,
        )
    )

    assert "No filtering rule is active" in text
    assert "Upstream harmonisation flow" not in text
    assert "Final filtering outcome" in text
    assert text.rstrip().endswith("100.00%")
    assert "MHC" not in text
    assert "Palindromic" not in text


def test_filtering_summary_merges_keep_missing_rule_with_other_results():
    group = {
        "display_group": "Imputation quality · FORMAT/SI",
        "display_group_kind": "analysis",
    }
    text = "\n".join(
        _filtering_summary_lines(
            [{
                "field": "FORMAT/SI",
                "label": "Missing FORMAT/SI",
                "category": "missing_value",
                "count": 12,
                "primary_count": 0,
                "action": "keep",
                "removes": False,
                "priority": 10,
                **group,
            }],
            [{
                "label": "Imputation quality below INFO (SI) 0.7",
                "count": 5,
                "primary_count": 5,
                "priority": 20,
                **group,
            }],
            variants_before=20,
            variants_after=15,
            terminal_label_width=42,
            reason_statistics={"primary_removed_total": 5},
        )
    )

    assert text.count("Filtering results by rule") == 1
    assert text.count("Imputation quality · FORMAT/SI") == 1
    assert "Missing values used by active filters" not in text
    assert text.index("Missing FORMAT/SI") < text.index(
        "Imputation quality below INFO (SI) 0.7"
    )
    assert (
        "12 variants matched this KEEP rule; this rule did not exclude any "
        "variants"
    ) in text


def test_filter_tag_summary_separates_observed_overlap_from_primary_reason(tmp_path):
    tag_file = tmp_path / "filter_tags.tsv"
    tag_file.write_text(
        ";PASS;\n;R1;\n;R1;R2;\n;R2;R3;\n;R10;\n",
        encoding="utf-8",
    )
    checks = [
        {"tag": "R1", "label": "Reason one", "removes": True},
        {"tag": "R2", "label": "Reason two", "removes": True},
        {"tag": "R3", "label": "Observed only", "removes": False},
        {"tag": "R10", "label": "Reason ten", "removes": True},
    ]

    with patch(
        "postgwas.modules.filtering.sumstat_filter.collect_streaming",
        wraps=collect_streaming,
    ) as collect:
        summary = _summarize_filter_tags(tag_file, checks)
    by_tag = {row["tag"]: row for row in summary["checks"]}

    assert collect.call_count == 1
    assert by_tag["R1"]["count"] == 2
    assert by_tag["R1"]["primary_count"] == 2
    assert by_tag["R2"]["count"] == 2
    assert by_tag["R2"]["primary_count"] == 1
    assert by_tag["R3"]["count"] == 1
    assert by_tag["R3"]["primary_count"] == 0
    assert by_tag["R10"]["count"] == 1
    assert by_tag["R10"]["primary_count"] == 1
    assert summary["primary_removed_total"] == 4
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

    report = _write_filter_reason_report(
        statistics, tmp_path / "study_filter_reason_summary.tsv"
    )

    with open(report, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["reason"] == "Reason one"
    assert rows[0]["filter_id"] == "R1"
    assert rows[0]["variants_matching_reason"] == "7"
    assert rows[0]["variants_removed_for_this_reason"] == "6"
    reasons = {row["reason"] for row in rows}
    assert "Rule failures beyond the first" in reasons
    assert "All primary removal assignments" in reasons
    assert rows[-1]["reason"] == "Reason reconciliation"
    assert rows[-1]["status"] == "PASS"


def test_filtering_csv_and_html_reuse_the_reconciled_summary(tmp_path):
    statistics = {
        "checks": [
            {
                "tag": "PGWAS_MAF",
                "category": "filter_condition",
                "label": "MAF outside range <unsafe>",
                "action": "remove",
                "removes": True,
                "count": 9,
                "primary_count": 8,
                "expr": "FORMAT/AF < 0.01",
                "priority": 20,
            },
            {
                "tag": "PGWAS_MISSING_SI",
                "category": "missing_value",
                "label": "Missing FORMAT/SI",
                "field": "FORMAT/SI",
                "action": "keep",
                "removes": False,
                "count": 2,
                "primary_count": 0,
                "expr": "FORMAT/SI == '.'",
                "priority": 30,
            },
        ],
        "overlap_variants": 1,
        "extra_rule_matches": 1,
        "primary_removed_total": 8,
        "actual_removed": 8,
        "reconciled": True,
    }
    outputs = {
        "filtered_vcf": str(tmp_path / "study_filtered.vcf.gz"),
        "reason_summary": str(tmp_path / "study_reasons.tsv"),
        "summary_csv": str(tmp_path / "study_summary.csv"),
        "html_report": str(tmp_path / "study_report.html"),
        "filter_log": str(tmp_path / "study.log"),
    }
    summary = build_filtering_summary(
        dataset_id="study<&",
        genome_build="GRCh37",
        missing_checks=(),
        condition_checks=(),
        reason_statistics=statistics,
        variants_before=100,
        variants_after=92,
        soft_filter_variants=None,
        soft_filter_enabled=False,
        input_validation={
            "status": "PASS",
            "input_vcf": str(tmp_path / "study.vcf.gz"),
            "genome_build": "GRCh37",
            "declared_contigs": 23,
            "total_variants": 100,
            "required_fields_evaluated": True,
            "required_fields": [
                {
                    "field": "FORMAT/AF",
                    "declared": True,
                    "required_by": ["minor-allele-frequency filter"],
                },
                {
                    "field": "INFO/EUR",
                    "declared": True,
                    "required_by": ["study/reference AF-difference filter"],
                },
            ],
            "field_provenance": {
                "status": "AVAILABLE",
                "detail": "PostGWAS 2.0; VCF status complete",
            },
            "field_summaries": [
                {
                    "kind": "genetic",
                    "label": "Allele-frequency fields",
                    "value": "FORMAT/AF for MAF; INFO/AF versus INFO/EUR",
                },
            ],
        },
        data_flow=None,
        outputs=outputs,
        resolved_configuration={
            "filtering": {"maf_min": 0.01, "remove_mhc": True},
            "execution": {"threads": 2, "memory_gb": 4},
            "executables": {"bcftools": "/tools/bcftools"},
        },
        runtime_seconds=1.25,
        display_missing_counts=True,
    )

    destination = write_filtering_summary_csv(
        summary, tmp_path / "study_summary.csv",
    )
    with destination.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert tuple(reader.fieldnames or ()) == SUMMARY_COLUMNS
    assert rows[0]["record_type"] == "overall"
    assert rows[0]["input_vcf"].endswith("study.vcf.gz")
    assert rows[0]["input_validation_status"] == "PASS"
    assert rows[0]["declared_contigs"] == "23"
    assert rows[0]["required_declared_fields_status"] == "PASS"
    assert rows[0]["required_declared_fields_present"] == "2"
    assert rows[0]["required_declared_fields_total"] == "2"
    assert rows[0]["variants_before"] == "100"
    assert rows[0]["variants_after"] == "92"
    assert rows[0]["variants_removed"] == "8"
    assert rows[0]["reconciliation_status"] == "PASS"
    assert rows[1]["filter_id"] == "PGWAS_MAF"
    assert rows[1]["variants_matching_reason"] == "9"
    assert rows[1]["variants_removed_for_this_reason"] == "8"
    assert rows[2]["filter_id"] == ""
    assert rows[2]["action"] == "keep"

    html = render_filtering_html(summary)
    assert html.startswith("<!doctype html>")
    assert "Variant filtering report" in html
    assert "Removal reconciliation:" in html and "PASS" in html
    assert "PGWAS_MAF" in html
    assert "8 primary removal assignments and 8 hard removals" in html
    assert "study&lt;&amp;" in html
    assert "MAF outside range &lt;unsafe&gt;" in html
    assert "study<&" not in html
    assert "MAF outside range <unsafe>" not in html
    assert "Evidence reuse:" in html
    assert "Input VCF validation" in html
    assert "Contract status:" in html and "PASS" in html
    assert "Header contract" in html
    assert "PASSED · 2/2 required INFO/FORMAT fields declared" in html
    assert "Field provenance" in html
    assert "PostGWAS 2.0; VCF status complete" in html
    assert "Allele-frequency fields" in html
    assert "FORMAT/AF for MAF; INFO/AF versus INFO/EUR" in html
    assert "Required field declarations" in html
    assert "FORMAT/AF" in html and "INFO/EUR" in html
    assert "minor-allele-frequency filter" in html
    assert "does not mean every record contains a value" in html
    assert "How to read rule counts" in html
    assert "1. Variants failing this rule" in html
    assert "2. Primary removals assigned here" in html
    assert "3. Overlapping earlier removal rules" in html
    assert (
        "overlapping earlier removal rules = variants failing this rule − "
        "primary removals assigned here"
    ) in html
    assert "Worked example from this run" in html
    assert "Variants failing this rule</dt><dd>9" in html
    assert "Primary removals assigned here</dt><dd>8" in html
    assert "Overlapping earlier removal rules</dt><dd>1" in html
    assert "2 observed; not a removal rule" in html
    assert "Not applicable" in html
    assert "every applied EXCLUDE or KEEP policy" in html
    assert "<td>EXCLUDE</td>" in html
    assert "<td>KEEP</td>" in html
