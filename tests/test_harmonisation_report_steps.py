"""Workflow cards interpret persisted evidence, never scientific inputs."""

from copy import deepcopy

import pytest

from postgwas.modules.harmonisation.report_steps import (
    chromosome_steps,
    dataset_steps,
    post_merge_steps,
)


@pytest.mark.parametrize("builder,total", [
    (dataset_steps, 8), (chromosome_steps, 16), (post_merge_steps, 5),
])
def test_missing_evidence_keeps_every_stage_without_inventing_success(builder, total):
    cards = builder({"status": "OK"})
    assert [card["number"] for card in cards] == list(range(1, total + 1))
    assert {card["status"] for card in cards} == {"Not recorded"}
    assert all(card["purpose"] and card["summary"] for card in cards)
    assert all(card["metrics"] == {} for card in cards)


def test_chromosome_order_and_recorded_zero_are_preserved():
    cards = chromosome_steps({"steps": [
        {"number": 3, "total": 16, "status": "ok", "rows_in": 100,
         "rows_out": 100, "removed": 0, "elapsed": 0.01},
    ], "stage_qc": {"beta_or_oddsratio_qc": {
        "effect_col": "OR_study", "effect_type": "odds_ratio",
        "conversion": "OR_to_Beta_log_transform_applied", "se_rescaled_by_or": False,
    }}})
    assert cards[2]["title"] == "Effect-scale conversion"
    assert cards[3]["title"] == "Allele alignment and effect frequency"
    assert cards[2]["status"] == "Completed"
    assert cards[2]["metrics"]["Rows removed at this step"] == 0
    assert cards[2]["summary"][0] == "Effect input column: OR_study."
    assert "Transformation: odds ratios converted to BETA using ln(OR)." in cards[2]["summary"]
    assert "SE converted from OR scale: no." in cards[2]["summary"]
    assert cards[3]["metrics"] == {}


def test_failed_step_is_not_masked_by_successful_nested_evidence():
    cards = chromosome_steps({"status": "failed", "steps": [
        {"number": 4, "total": 16, "status": "failed", "error": "Reference unavailable"},
    ], "stage_qc": {"eaf_qc": {"status": "success"}}})
    assert cards[3]["status"] == "Failed"
    assert cards[3]["summary"][0] == "Reference unavailable"
    assert cards[4]["status"] == "Not recorded"


def test_skipped_effect_conversion_is_not_needed():
    cards = chromosome_steps({"steps": [{"number": 3, "total": 16, "status": "ok"}],
                              "stage_qc": {"beta_or_oddsratio_qc": {
        "status": "skipped_no_input_effect_column", "effect_type": "beta",
    }}})
    assert cards[2]["status"] == "Not needed"


def test_strand_and_eaf_rejections_are_not_added_twice():
    payload = {"variants_removed": 51, "strand_orientation": {
        "palindromic_frequency_conflict": 3, "final_variants": 197084,
    }, "final_variants": 197036}
    card = chromosome_steps({"steps": [
        {"number": 4, "total": 16, "status": "ok", "removed": 51},
    ], "stage_qc": {"eaf_qc": payload}})[3]
    assert card["metrics"]["Rows removed at this step"] == 51
    assert card["details"]["Observed evidence"] == payload
    assert "Total removed within this combined stage: 51." in card["summary"]


def test_liftover_uses_actual_records_not_logged_harmonised_rows():
    cards = chromosome_steps({"steps": [
        {"number": 16, "total": 16, "status": "ok", "rows_in": 100,
         "rows_out": 100, "removed": 0, "extra": {"liftover_accounting": {
             "input": 99, "rejected": 1, "swap_excluded": 2, "final": 96,
         }}},
    ]})
    assert "Rows leaving step" not in cards[15]["metrics"]
    assert cards[15]["metrics"]["Target-build VCF records"] == 96
    assert cards[15]["metrics"]["Liftover plugin rejections"] == 1
    assert cards[15]["metrics"]["Successfully lifted variants excluded by swap policy"] == 2


def test_dataset_evidence_is_not_a_fabricated_completion_record():
    cards = dataset_steps({"status": "OK", "sumstat_file": "/input/a_b.gz", "dataset": {
        "genome_build": {"inferred_build": "GRCh37", "testable_variants": 80},
        "study_decisions": {"strand": "forward", "strand_consensus": {"forward": 70}},
    }})
    assert cards[4]["status"] == "Evidence recorded"
    assert "Build-testable study variants: 80." in cards[4]["summary"]
    assert cards[1]["status"] == "Not recorded"
    assert "Study input: /input/a_b.gz." in cards[2]["summary"]


def test_postmerge_missing_failed_disabled_and_virtual_qc_are_distinct():
    cards = post_merge_steps({"merged_vcfs": {"merge_status": "PARTIAL"},
                              "population_frequency_qc": {"status": "disabled"},
                              "qc_assessment": {"raw_variants": 100, "qc_passed_variants": 95}})
    assert cards[0]["status"] == "Partial"
    assert cards[1]["status"] == "Not needed"
    assert cards[2]["status"] == "Not recorded"
    assert cards[3]["status"] == "Evidence recorded"
    assert "does not itself filter" in cards[3]["purpose"]
    assert "Virtual QC-passed records: 95." in cards[3]["summary"]
    failed = post_merge_steps({"population_frequency_qc": {"status": "failed_warning"}})[1]
    assert failed["status"] == "Failed — continued by policy"


def test_recorded_policies_are_reused_without_defaults_and_inputs_are_immutable():
    summary = {"stage_qc": {"eaf_qc": {"strand_orientation": {"actions": {"forward": 3}}}},
               "qc_actions": [{"step": "04 strand_orientation", "removed": 2}]}
    policies = {"eaf.external_min_match_fraction": 0.73, "strand.mode": "auto",
                "final_check.require": ["beta", "se"], "other.setting": 4}
    original = deepcopy((summary, policies))
    cards = chromosome_steps(summary, policies)
    assert cards[3]["policies"] == {"eaf.external_min_match_fraction": 0.73, "strand.mode": "auto"}
    assert cards[12]["policies"] == {"final_check.require": ["beta", "se"]}
    cards[3]["details"]["Observed evidence"]["strand_orientation"]["actions"]["forward"] = 999
    cards[12]["policies"]["final_check.require"].append("pval")
    assert (summary, policies) == original


def test_step_totals_separate_study_and_postmerge_numbering():
    manifest = {"steps": [
        {"number": 1, "total": 8, "status": "ok"},
        {"number": 1, "total": 5, "status": "failed", "error": "Merge failed"},
    ]}
    assert dataset_steps(manifest)[0]["status"] == "Completed"
    assert post_merge_steps(manifest)[0]["status"] == "Failed"


def test_input_stage_does_not_relabel_downstream_removals_as_input_qc():
    card = dataset_steps({"pre_vcf_summary": {
        "total_variant_read": 100, "total_variant_remaining_for_harmonisation": 95,
        "total_variant_removed_chromosome_harmonisation": 10,
        "total_variant_in_vcf_input": 85,
    }})[2]
    assert card["metrics"] == {"Variants read": 100, "Ready for chromosome analysis": 95}
    assert "total_variant_in_vcf_input" not in card["details"]["Observed evidence"]["Input accounting"]


def test_frequency_groups_do_not_describe_warn_only_disagreements_as_changes():
    card = chromosome_steps({"stage_qc": {"eaf_qc": {
        "strand_orientation": {"actions": {"forward_swapped": 197084},
                               "reference_file": "/reference/AF_ref.tsv.gz",
                               "palindromic_frequency_conflict": 3},
        "strand_non_palindromic_af_comparable": 161866,
        "strand_non_palindromic_af_discordant": 200,
        "strand_af_discordance_action": "warn", "variants_removed": 51,
    }}, "qc_actions": [{"step": "04 strand_orientation", "changed": 200}]})[3]
    strand, frequency = card["groups"]
    assert strand["metrics"]["Alleles swapped"] == 197084
    assert strand["metrics"]["Reference file"] == "/reference/AF_ref.tsv.gz"
    assert "Alleles unchanged" not in strand["metrics"]
    assert frequency["metrics"]["Non-palindromic discordant"] == 200
    assert frequency["metrics"]["Combined strand/EAF removals"] == 51
    assert "retained with warning; no frequencies changed by this check" in frequency["note"]
    assert "already include strand losses" in frequency["note"]


def test_vcf_groups_reuse_successive_counts_and_distinguish_liftover_losses():
    card = chromosome_steps({"steps": [{"number": 16, "total": 16, "status": "ok",
        "extra": {"liftover_accounting": {"input": 100, "rejected": 2, "swap_excluded": 3, "final": 95}}}],
        "qc_actions": [{"step": "16 bcftools_annotate_liftover", "check": "NORM", "before": 101, "after": 100},
                       {"step": "16 bcftools_annotate_liftover", "check": "ID", "after": 100}]})[15]
    normal, annotation, lifted = card["groups"]
    assert normal["metrics"]["Normalized VCF records"] == 100
    assert annotation["metrics"] == {"After ID annotation": 100}
    assert lifted["metrics"]["Plugin rejections"] == 2
    assert lifted["metrics"]["Successfully lifted variants excluded by swap policy"] == 3
    assert "not failed liftover" in lifted["note"]


def test_sample_size_outcome_is_readable_without_renaming_columns():
    card = chromosome_steps({"stage_qc": {"sample_size_qc": {
        "Neff_status": "fallback_ncontrol_only:N_study",
        "ncontrol_source": "column:N_study",
    }}})[4]
    assert "Sample-size handling: control-count input used as total sample size; column N_study." in card["summary"]
    assert "Control or total-count source: column:N_study." in card["summary"]
