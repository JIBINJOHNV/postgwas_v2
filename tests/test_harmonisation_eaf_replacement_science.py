"""Conditional replacement uses the same reference evidence and keeps raw MAF."""
from unittest.mock import patch

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import allele_frequency as af, strand
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN
from postgwas.modules.harmonisation.chromosome_partition import ChromosomePartitions
from postgwas.modules.harmonisation.shared.runtime import NullStepContext
from postgwas.modules.harmonisation.study_properties import _decide_eaf_is_maf


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
COLUMNS = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}


def fixture(tmp_path, *, genuine_eaf=False):
    reference = tmp_path / "reference.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t0.8\n1\t200\tA\tC\t0.7\n1\t300\tA\tAC\t0.6\n2\t100\tC\tT\t0.2\n2\t200\tC\tA\t0.3\n")
    frame = pl.DataFrame({"CHR": ["1", "1", "1", "2", "2"], "POS": [100, 200, 300, 100, 200],
        "EA": ["A", "C", "AC", "C", "C"], "OA": ["G", "A", "A", "T", "A"],
        "AF": [0.8, 0.7, 0.6, 0.8, 0.7] if genuine_eaf else [0.2, 0.3, 0.4, 0.2, 0.3]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2})
    return reference, frame, policies


def prepare(frame, reference, policies, **kwargs):
    decisions = {"strand": "forward", "effect_type": "beta"}
    _decide_eaf_is_maf(frame, COLUMNS, policies, NullStepContext(), decisions, {})
    return af.validate_study_allele_frequency(frame, COLUMNS, decisions,
        {chrom: {"default_eaf_file": str(reference), "default_eaf_reference_column": "EUR"}
         for chrom in frame["CHR"].unique()}, MAPPING, policies=policies, **kwargs)


def approve(prepared, reference, policies, decision):
    """Exercise the real dataset gate before the existing adoption regressions."""
    paths, rows, maps = {}, {}, {}
    replacement = decision["eaf_replacement"]["file"]
    decision["eaf_replacement"]["resolved_files"] = {
        chrom: replacement for chrom in prepared["CHR"].unique()}
    for chrom in prepared["CHR"].unique():
        path = reference.parent / ("approved_" + chrom + ".parquet")
        subset = prepared.filter(pl.col("CHR") == chrom)
        subset.write_parquet(path)
        paths[chrom], rows[chrom] = str(path), subset.height
        maps[chrom] = {"default_eaf_file": str(reference), "default_eaf_reference_column": "EUR",
            "user_eaf_file": replacement if replacement != str(reference) else None}
    return af.validate_replacement_allele_frequency(ChromosomePartitions(paths, rows), COLUMNS,
        decision, maps, MAPPING, policies=policies)


def test_supported_maf_selects_replacement_with_one_pooled_match_and_keeps_indel(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    with patch.object(strand, "_load_reference", wraps=strand._load_reference) as reads:
        decision, prepared = prepare(frame, reference, policies,
            eaf_replacement_file=str(reference), eaf_replacement_column="EUR")
    assert reads.call_count == 1
    assert decision["eaf_selected_source"] == "replacement_reference"
    assert decision["eaf_input_interpretation"] == "maf"
    assert decision["eaf_reference_validation"]["status"] == "replacement_selected"
    assert decision["eaf_reference_validation"]["comparable_variants"] == 5
    assert decision["eaf_is_maf"] is True
    assert prepared.select(frame.columns).equals(frame)


def test_supported_maf_without_replacement_names_both_options(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    with pytest.raises(af.StudyAlleleFrequencyError, match="(?s)MAF.*eaf-replacement-file.*eaf-replacement-column"):
        prepare(frame, reference, policies)


def test_minor_ea_genuine_eaf_retained_even_with_replacement_configured(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    # List the minor allele as EA; this is genuine EAF, not invalid MAF.
    frame = frame.with_columns(pl.col("OA").alias("EA"), pl.col("EA").alias("OA"))
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file="/unused/not/read.tsv", eaf_replacement_column="EUR")
    assert decision["eaf_is_maf"] is False
    assert decision["eaf_selected_source"] == "study_supplied"
    assert decision["eaf_reference_validation"]["status"] == "passed"
    assert prepared.select(frame.columns).equals(frame)


def test_same_reference_replacement_is_fresh_column_without_reclassification_or_rematch(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(reference), eaf_replacement_column="EUR")
    decision["eaf_replacement"]["resolved_files"] = {"1": str(reference), "2": str(reference)}
    decision = approve(prepared, reference, policies, decision)
    columns = dict(COLUMNS, output_folder=str(tmp_path / "out"), gwas_outputname="study")
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    work = prepared.filter(pl.col("CHR") == "1")
    with patch.object(strand, "_load_reference", side_effect=AssertionError("reference reread")), \
         patch.object(af, "confirm_eaf_with_reference", side_effect=AssertionError("replacement classified")):
        result, qc, selected = af.harmonise_allele_frequency("1", work, columns, layout, "\t",
            external_eaf_colmap=MAPPING, default_eaf_file=str(reference), default_eaf_column="EUR",
            default_eaf_colmap=MAPPING, study_decision=decision, policies=policies)
    assert result["AF"].to_list() == pytest.approx([0.2, 0.3, 0.4])
    assert selected["eaf_col"] not in ("AF", strand.REFERENCE_AF_COLUMN)
    assert result[selected["eaf_col"]].to_list() == pytest.approx([0.8, 0.7, 0.6])
    assert qc["eaf_provenance"] == "reference_imputed"
    assert qc["eaf_selected_source"] == "replacement_reference"


def test_split_multiallelic_site_is_not_semantic_diagnostic_evidence(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    with reference.open("a") as handle:
        handle.write("1\t100\tG\tC\t0.1\n")
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(reference), eaf_replacement_column="EUR")
    assert decision["eaf_reference_validation"]["comparable_variants"] == 4
    assert decision["eaf_reference_validation"]["diagnostic_excluded_non_biallelic_rows"] == 1
    assert prepared.select(frame.columns).equals(frame)


def test_constant_frequency_does_not_activate_replacement(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    frame = frame.with_columns(pl.lit(.2).alias("AF"))
    with pytest.raises(af.StudyAlleleFrequencyError, match="inconclusive|correlation_unavailable"):
        prepare(frame, reference, policies, eaf_replacement_file=str(reference), eaf_replacement_column="EUR")


def test_declared_replacement_does_not_screen_its_minor_alt_values(tmp_path):
    reference, frame, policies = fixture(tmp_path)
    replacement = tmp_path / "replacement.tsv"
    replacement.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t0.2\n1\t200\tA\tC\t0.3\n1\t300\tA\tAC\t0.4\n2\t100\tC\tT\t0.2\n2\t200\tC\tA\t0.3\n")
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(replacement), eaf_replacement_column="EUR")
    decision["eaf_replacement"]["resolved_files"] = {"1": str(replacement)}
    columns = dict(COLUMNS, output_folder=str(tmp_path / "out"), gwas_outputname="study")
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    with patch("postgwas.modules.harmonisation.shared.allele_join.frequency_maf_screen",
               side_effect=AssertionError("replacement distribution screened")):
        decision = approve(prepared, reference, policies, decision)
        result, qc, selected = af.harmonise_allele_frequency("1", prepared.filter(pl.col("CHR") == "1"),
            columns, layout, "\t", eaffile=str(replacement), external_eaf_colmap=MAPPING,
            default_eaf_file=str(reference), default_eaf_column="EUR", default_eaf_colmap=MAPPING,
            study_decision=decision, policies=policies)
    assert result[selected["eaf_col"]].to_list() == pytest.approx([0.2, 0.3, 0.4])
    assert result["AF"].to_list() == pytest.approx([0.2, 0.3, 0.4])
    assert "external_raw_af_maf_like" not in qc


def test_declared_replacement_ref_frequency_uses_all_known_alt_alleles(tmp_path):
    replacement = tmp_path / "replacement.tsv"
    replacement.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tA\tG\t0.2\n1\t100\tA\tC\t0.1\n")
    frame = pl.DataFrame({"CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"]})
    result, selected, stats = af.merge_external_allele_frequencies(
        frame, str(replacement), "EUR", MAPPING,
        {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"},
        declared_alt_frequency=True,
    )
    assert result[selected].to_list() == pytest.approx([0.7])
    assert stats["reference_frequency_reconstructed_rows"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 0


def test_diagnostic_pairs_require_site_eligibility_evidence():
    candidates = pl.DataFrame({"row": [0], "AF": [.2], "__strand_af": [.8],
        "__strand_swapped": [False], "__strand_palindromic": [False], "__strand_candidates": [1]})
    with pytest.raises(strand.StrandOrientationError, match="biallelic.*eligibility"):
        strand._study_reference_frequency_pairs(candidates, "row", "AF")


@pytest.mark.parametrize("chromosome,expected_beta,expected_af", [
    ("1", [.5, .5, .5], [.8, .7, .6]), ("2", [-.5, -.5], [.2, .3]),
])
def test_replacement_composes_with_swap_complement_and_beta_z_orientation(tmp_path, chromosome, expected_beta, expected_af):
    reference, frame, policies = fixture(tmp_path)
    frame = frame.with_columns(pl.lit(.5).alias("BETA"), pl.lit(2.).alias("Z"),
        pl.when(pl.col(SOURCE_INPUT_ROW_COLUMN) == 0).then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"),
        pl.when(pl.col(SOURCE_INPUT_ROW_COLUMN) == 0).then(pl.lit("C")).otherwise(pl.col("OA")).alias("OA"))
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(reference), eaf_replacement_column="EUR")
    decision["eaf_replacement"]["resolved_files"] = {chromosome: str(reference)}
    decision = approve(prepared, reference, policies, decision)
    columns = dict(COLUMNS, beta_or_col="BETA", imp_z_col="Z", output_folder=str(tmp_path / "out"), gwas_outputname="study")
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    result, qc, selected = af.harmonise_allele_frequency(chromosome, prepared.filter(pl.col("CHR") == chromosome),
        columns, layout, "\t", external_eaf_colmap=MAPPING, default_eaf_file=str(reference),
        default_eaf_column="EUR", default_eaf_colmap=MAPPING, study_decision=decision, policies=policies)
    assert result["BETA"].to_list() == pytest.approx(expected_beta)
    assert result["Z"].to_list() == pytest.approx([2., 2., 2.] if chromosome == "1" else [-2., -2.])
    assert result[selected["eaf_col"]].to_list() == pytest.approx(expected_af)
    assert result["AF"].to_list() == pytest.approx([.2, .3, .4] if chromosome == "1" else [.2, .3])


@pytest.mark.parametrize("discarded", [None, float("nan"), float("inf"), -1.])
def test_selected_reference_qc_retains_invalid_original_maf_as_provenance(tmp_path, discarded):
    reference, frame, policies = fixture(tmp_path)
    frame = frame.with_columns(pl.Series("AF", [.2, .3, discarded, .2, .3]))
    policies = policies.with_overrides({"eaf.out_of_range": "fail"})
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(reference), eaf_replacement_column="EUR")
    decision["eaf_replacement"]["resolved_files"] = {"1": str(reference)}
    decision = approve(prepared, reference, policies, decision)
    columns = dict(COLUMNS, output_folder=str(tmp_path / "out"), gwas_outputname="study")
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    result, qc, selected = af.harmonise_allele_frequency("1", prepared.filter(pl.col("CHR") == "1"),
        columns, layout, "\t", external_eaf_colmap=MAPPING, default_eaf_file=str(reference),
        default_eaf_column="EUR", default_eaf_colmap=MAPPING, study_decision=decision, policies=policies)
    assert result["AF"].equals(frame.filter(pl.col("CHR") == "1")["AF"])
    assert result[selected["eaf_col"]].to_list() == pytest.approx([.8, .7, .6])
    assert qc["eaf_replacement_rows"] == 3


@pytest.mark.parametrize("endpoint_action,expected_rows", [("reject", 3), ("null", 5)])
def test_replacement_count_reports_final_nonmissing_values_after_endpoint_policy(tmp_path, endpoint_action, expected_rows):
    reference = tmp_path / "comparison.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t1\n1\t200\tG\tA\t1\n1\t300\tG\tA\t0.6\n1\t400\tG\tA\t0.65\n1\t500\tG\tA\t0.75\n")
    replacement = tmp_path / "replacement.tsv"
    replacement.write_text("CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t0\n1\t200\tG\tA\t1\n1\t300\tG\tA\t0.6\n1\t400\tG\tA\tNA\n1\t500\tG\tA\t0.75\n")
    frame = pl.DataFrame({"CHR": ["1"] * 5, "POS": [100, 200, 300, 400, 500],
        "EA": ["A"] * 5, "OA": ["G"] * 5, "AF": [0., 0., .4, .35, .25]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2, "eaf.degenerate": endpoint_action})
    decision, prepared = prepare(frame, reference, policies,
        eaf_replacement_file=str(replacement), eaf_replacement_column="EUR")
    decision["eaf_replacement"]["resolved_files"] = {"1": str(replacement)}
    decision = approve(prepared, reference, policies, decision)
    columns = dict(COLUMNS, output_folder=str(tmp_path / "out"), gwas_outputname="study")
    layout = load_configuration().modules.harmonisation.output_layout.model_dump()
    result, qc, selected = af.harmonise_allele_frequency("1", prepared, columns, layout, "\t",
        eaffile=str(replacement), external_eaf_colmap=MAPPING, default_eaf_file=str(reference),
        default_eaf_column="EUR", default_eaf_colmap=MAPPING, study_decision=decision, policies=policies)
    assert result.height == expected_rows
    assert result[selected["eaf_col"]].drop_nulls().to_list() == pytest.approx([.6, .75])
    assert qc["eaf_replacement_rows"] == 2
