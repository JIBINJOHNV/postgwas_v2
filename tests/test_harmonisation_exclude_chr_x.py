"""Explicit chromosome-X exclusions precede scientific input QC."""

import polars as pl
import pytest

from postgwas.modules.harmonisation.coordinates import harmonise_coordinates_and_alleles
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import RejectCollector, SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.summary_statistics_io import read_summary_statistics


MAPPING = {
    "gwas_outputname": "study", "chr_col": "CHR", "pos_col": "POS",
    "ea_col": "EA", "oa_col": "OA", "snp_id_col": "SNP",
    "beta_or_col": "BETA", "pval_col": "P",
}
LAYOUT = {
    "rejected_directory": "rejected", "input_reject": "rejected/{dataset_id}_input.tsv",
    "duplicates": "{dataset_id}_duplicates.tsv",
}


def _read(tmp_path, text, *, exclude=True, mapping=None, overrides=None):
    path = tmp_path / "study.tsv"
    path.write_text(text)
    policies = default_policies().with_overrides({
        "chromosome.exclude_x": exclude, **(overrides or {}),
    })
    return read_summary_statistics(
        str(path), str(tmp_path), dict(mapping or MAPPING), LAYOUT, "\t",
        policies=policies, input_line_count=len(text.splitlines()) - 1,
    )


def test_explicit_exclusion_precedes_missingness_position_allele_and_duplicate_qc(tmp_path):
    text = (
        "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n"
        "1\t100\tA\tG\ta\t0.1\t0.01\n"
        "chrX\tbad\tN\tG\tx\tNA\tNA\n"
        "23\t100\tA\tG\tx\t0.1\t0.01\n"
        "X\t100\tA\tG\tx\t0.2\t0.02\n"
        "2\t200\tC\tT\tb\t-0.2\t0.02\n"
    )
    result = _read(tmp_path, text)
    retained, source = result[:2]
    assert retained["SNP"].to_list() == ["a", "b"]
    assert retained["BETA"].to_list() == [0.1, -0.2]
    assert source["CHR"].to_list() == ["1", "chrX", "23", "X", "2"]
    assert result[5:10] == (0, 0, 0, 0, 0)
    report = result[-1]
    assert report["input_rows"] == 2
    assert report["chromosome_exclusion"] == {
        "enabled": True, "excluded_chromosomes": ["X"], "excluded_rows": 3,
        "input_rows": 5, "retained_rows": 2,
    }
    rejected = pl.read_csv(tmp_path / "rejected/study_input.tsv.gz", separator="\t")
    assert rejected["SNP"].to_list() == ["x", "x", "x"]
    assert rejected["reject_reason"].to_list() == ["chromosome_excluded_by_user"] * 3
    assert rejected["POS"].cast(pl.String).to_list() == ["bad", "100", "100"]
    assert 5 == retained.height + rejected.height


@pytest.mark.parametrize("combined", [False, True])
def test_all_packaged_x_aliases_excluded_and_autosome_order_preserved(tmp_path, combined):
    labels = ["2", "X", "chrX", "23", "25", "XY", "PAR1", "PAR2", "chr023", "1"]
    header = "COORD" if combined else "CHR\tPOS"
    rows = [
        (label + ":bad" if label != "1" and label != "2" else label + ":100")
        if combined else label + "\t100"
        for label in labels
    ]
    text = header + "\tEA\tOA\tSNP\tBETA\tP\n" + "".join(
        row + "\tA\tG\tv%d\t0.1\t0.01\n" % index
        for index, row in enumerate(rows)
    )
    mapping = {**MAPPING, "chr_col": None, "pos_col": None, "chr_pos_col": "COORD"} if combined else MAPPING
    result = _read(tmp_path, text, mapping=mapping)
    assert result[0]["SNP"].to_list() == ["v0", "v9"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 8


def test_excluded_combined_x_with_empty_position_cannot_break_coordinate_split(tmp_path):
    text = "COORD\tEA\tOA\tSNP\tBETA\tP\n1:100\tA\tG\ta\t0.1\t0.01\nX:\tN\tN\tx\tNA\tNA\n"
    mapping = {**MAPPING, "chr_col": None, "pos_col": None, "chr_pos_col": "COORD"}
    result = _read(tmp_path, text, mapping=mapping)
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1


def test_explicit_chromosome_mapping_is_authoritative_for_early_exclusion(tmp_path):
    text = "CHR\tCOORD\tEA\tOA\tSNP\tBETA\tP\n1\t1:100\tA\tG\ta\t0.1\t0.01\nX\t2:bad\tN\tN\tx\tNA\tNA\n"
    mapping = {**MAPPING, "pos_col": None, "chr_pos_col": "COORD"}
    result = _read(tmp_path, text, mapping=mapping)
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1


def test_x_is_included_by_default_and_exclusion_without_x_preserves_values(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n1\t100\tA\tG\ta\t0.1\t0.01\nX\t200\tC\tT\tx\t-0.2\t0.02\n"
    result = _read(tmp_path, text, exclude=False)
    assert result[0]["SNP"].to_list() == ["a", "x"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 0
    only_autosome = tmp_path / "autosome"
    only_autosome.mkdir()
    result = _read(only_autosome, text.split("X\t200")[0])
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 0


def test_x_only_failure_flushes_original_exclusion_provenance(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\nX\tbad\tN\tN\tx\tNA\tNA\n"
    with pytest.raises(ValueError, match="No variants remain.*--exclude-chr-x"):
        _read(tmp_path, text)
    rejected = pl.read_csv(tmp_path / "rejected/study_input.tsv.gz", separator="\t")
    assert rejected["reject_reason"].to_list() == ["chromosome_excluded_by_user"]
    assert rejected["POS"].to_list() == ["bad"]


def test_disabled_rejection_files_still_balance_explicit_exclusion(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n1\t100\tA\tG\ta\t0.1\t0.01\nX\tbad\tN\tN\tx\tNA\tNA\n"
    result = _read(tmp_path, text, overrides={"rejects.enabled": False})
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1
    assert not (tmp_path / "rejected").exists()


def test_direct_coordinate_call_uses_same_exclusion_before_invalid_position(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["X", "1"], "POS": ["bad", "100"],
        "EA": ["N", "A"], "OA": ["N", "G"],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = RejectCollector(source_snapshot=frame, logger=None,
        out_path=str(tmp_path / "rejected.tsv"), delimiter="\t", compress=False)
    retained, _, qc = harmonise_coordinates_and_alleles(
        "All_Chrs", frame, dict(MAPPING),
        policies=default_policies().with_overrides({"chromosome.exclude_x": True}),
        rejects=collector, return_qc_info=True,
    )
    assert retained["CHR"].to_list() == ["1"]
    assert qc["removed_by_reason"]["chromosome_excluded_by_user"] == 1
    assert qc["removed_by_reason"]["invalid_position"] == 0
    assert collector.frame()["POS"].to_list() == ["bad"]


def test_missing_chromosomes_remain_subject_to_existing_invalid_coordinate_policy(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n1\t100\tA\tG\ta\t0.1\t0.01\nNA\t100\tA\tG\tunknown\t0.1\t0.01\n"
    result = _read(tmp_path, text)
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 0


def test_custom_alias_mapping_is_respected(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n1\t100\tA\tG\ta\t0.1\t0.01\nCUSTOMX\tbad\tN\tN\tx\tNA\tNA\n"
    result = _read(tmp_path, text, overrides={"chromosome.rename_map": {"CUSTOMX": "X"}})
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1


def test_numeric_inference_cannot_reject_values_only_on_excluded_x(tmp_path):
    text = "CHR\tPOS\tEA\tOA\tSNP\tBETA\tP\n1\t100\tA\tG\ta\t0.1\t0.01\nX\tbad\tN\tN\tx\tbad\tNA\n"
    result = _read(tmp_path, text, overrides={"input.schema_inference_rows": 1})
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[0]["BETA"].to_list() == [0.1]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1


def test_reference_scope_uses_configured_aliases_and_preserves_retained_order():
    from postgwas.modules.harmonisation.shared.variant_columns import retained_chromosome_scope
    policies = default_policies().with_overrides({"chromosome.exclude_x": True})
    assert retained_chromosome_scope(["2", "chrX", "PAR1", "1", "MT"], policies) == ["2", "1", "MT"]
    assert retained_chromosome_scope(["X", "23"], policies) == []
    assert retained_chromosome_scope([], policies) == []


def test_combined_scope_respects_configured_prefix_normalization(tmp_path):
    text = "COORD\tEA\tOA\tSNP\tBETA\tP\n1:100\tA\tG\ta\t0.1\t0.01\nchrX:200\tA\tG\tprefix\t0.1\t0.01\nX:bad\tN\tN\tx\tNA\tNA\n"
    mapping = {**MAPPING, "chr_col": None, "pos_col": None, "chr_pos_col": "COORD"}
    result = _read(tmp_path, text, mapping=mapping,
                   overrides={"chromosome.strip_chr_prefix": False})
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1
    assert result[6] == 1


def test_shared_external_reference_x_cannot_affect_duplicate_or_frequency_screen(tmp_path):
    from postgwas.modules.harmonisation.allele_frequency import merge_external_allele_frequencies
    reference = tmp_path / "reference.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t0.2\nX\t200\tG\tA\t0.8\nX\t200\tG\tA\t0.9\n")
    original = reference.read_bytes()
    frame = pl.DataFrame({"CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"]})
    policies = default_policies().with_overrides({
        "chromosome.exclude_x": True, "external_reference.non_identical_duplicate_action": "fail",
    })
    result, column, stats = merge_external_allele_frequencies(
        frame, str(reference), "EUR",
        {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"},
        {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"},
        policies=policies, chromosome="1",
    )
    assert result[column].to_list() == [0.2]
    assert stats["external_rows_loaded"] == 1
    assert stats["raw_frequency_screen"]["usable"] == 1
    assert stats["raw_frequency_screen"]["maf_like"] is True
    assert reference.read_bytes() == original


@pytest.mark.parametrize("alias,separator", [("SEX_X", ":"), ("SEX-X", "_"), ("SEX_X", "-")])
def test_combined_custom_x_aliases_use_existing_coordinate_separator(tmp_path, alias, separator):
    text = ("COORD\tEA\tOA\tSNP\tBETA\tP\n"
            + "1" + separator + "100\tA\tG\ta\t0.1\t0.01\n"
            + alias + separator + "bad\tN\tN\tx\tNA\tNA\n")
    mapping = {**MAPPING, "chr_col": None, "pos_col": None, "chr_pos_col": "COORD"}
    result = _read(tmp_path, text, mapping=mapping,
                   overrides={"chromosome.rename_map": {alias: "X"}})
    assert result[0]["SNP"].to_list() == ["a"]
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 1
    assert result[6] == 0


def test_concordance_partition_expression_resolves_custom_combined_x_aliases():
    from types import SimpleNamespace
    from postgwas.modules.harmonisation.concordance.analysis import input_chromosome_expression
    row = SimpleNamespace(chromosome_column=None, position_column=None, chromosome_position_column="CP")
    policies = default_policies().with_overrides({
        "chromosome.exclude_x": True, "chromosome.rename_map": {"SEX_X": "X", "SEX-X": "X"},
    })
    frame = pl.DataFrame({"CP": ["1:100", "SEX_X:200", "SEX-X_200"]})
    assert frame.select(input_chromosome_expression(row, policies)).to_series().to_list() == ["1", "X", "X"]


def test_unknown_combined_label_is_not_guessed_as_x(tmp_path):
    text = "COORD\tEA\tOA\tSNP\tBETA\tP\n1:100\tA\tG\ta\t0.1\t0.01\nX_1:200\tA\tG\tunknown\t0.1\t0.01\n"
    mapping = {**MAPPING, "chr_col": None, "pos_col": None, "chr_pos_col": "COORD"}
    result = _read(tmp_path, text, mapping=mapping)
    assert result[-1]["chromosome_exclusion"]["excluded_rows"] == 0
    assert result[6] == 1
    assert result[0]["SNP"].to_list() == ["a"]


def test_reference_scope_records_original_and_retained_counts_in_shared_audit(tmp_path):
    from postgwas.core.input_validation import InputValidationSession
    from postgwas.modules.harmonisation.shared.variant_columns import read_reference_variant_table
    source = tmp_path / "reference.tsv"
    source.write_text("CHR\tPOS\tALT\tREF\n1\t100\tA\tG\nX\tbad\tN\tN\n")
    policies = default_policies().with_overrides({"chromosome.exclude_x": True})
    with InputValidationSession(pipeline=False) as audit:
        frame, _ = read_reference_variant_table(source,
            {"chr": "CHR", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"},
            policies)
    assert frame["CHR"].to_list() == ["1"]
    scoped = [record for record in audit.records
              if "reference_chromosome_excluded_rows" in record.metrics]
    assert len(scoped) == 1
    assert scoped[0].metrics["reference_rows_before_chromosome_exclusion"] == 2
    assert scoped[0].metrics["reference_chromosome_excluded_rows"] == 1
    assert scoped[0].metrics["reference_rows_after_chromosome_exclusion"] == 1
