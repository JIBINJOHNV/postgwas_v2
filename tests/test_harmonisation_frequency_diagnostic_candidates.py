"""Scoped AF diagnostics may include resolved multiallelic DNA allele pairs."""

import polars as pl
import pytest

from postgwas.modules.harmonisation import strand
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN


def _candidates():
    return pl.DataFrame({
        "row": [0, 1, 2, 3, 3, 4, 5, 6],
        "__strand_ref": ["A", "A", "A", "A", "A", "A", "A", "AT"],
        "__strand_alt": ["G", "C", "T", "G", "C", "*", "A", "A"],
        "__strand_biallelic": [False, True, True, False, False, False, False, False],
        "__strand_palindromic": [False, False, True, False, False, False, False, False],
    })


def test_default_diagnostic_selection_keeps_biallelic_requirement():
    assert strand.frequency_diagnostic_candidates(_candidates(), "row")["row"].to_list() == [1]


def test_scoped_diagnostics_drop_only_multiplicity_not_orientation_or_dna_guards():
    selected = strand.frequency_diagnostic_candidates(_candidates(), "row", require_biallelic=False)

    assert selected["row"].to_list() == [0, 1, 6]
    assert selected["__strand_biallelic"].to_list() == [False, True, False]


def test_scoped_diagnostics_preserve_precomputed_candidate_ambiguity():
    candidates = _candidates().filter(pl.col("row") == 0).with_columns(pl.lit(2).alias("__strand_candidates"))

    assert strand.frequency_diagnostic_candidates(candidates, "row", require_biallelic=False).is_empty()


@pytest.mark.parametrize("require_biallelic", [True, False])
def test_structural_membership_metadata_remains_required(require_biallelic):
    with pytest.raises(strand.StrandOrientationError, match="biallelic.*eligibility"):
        strand.frequency_diagnostic_candidates(_candidates().drop("__strand_biallelic"), "row",
                                               require_biallelic=require_biallelic)


def test_scoped_frequency_pairs_keep_original_af_and_requested_effect_direction():
    candidates = pl.DataFrame({
        "row": [0, 1, 2, 3, 4], "AF": [.12, .23, .3, None, .4],
        "__strand_ref": ["A"] * 5, "__strand_alt": ["G"] * 5,
        "__strand_af": [.8, .8, .3, .4, float("nan")],
        "__strand_ref_af": [.1, .6, .7, .6, .5],
        "__strand_swapped": [False, True, False, False, False],
        "__strand_biallelic": [False, False, True, True, True],
        "__strand_palindromic": [False] * 5,
    })
    before = candidates.clone()

    selected = strand._study_reference_frequency_pairs(candidates, "row", "AF", require_biallelic=False)

    assert selected["__study_frequency"].to_list() == [.12, .23, .3]
    assert selected["__expected_effect_frequency"].to_list() == pytest.approx([.8, .6, .3])
    assert candidates.equals(before)
    assert strand._study_reference_frequency_pairs(candidates, "row", "AF")["__study_frequency"].to_list() == [.3]


@pytest.mark.parametrize("extra,ea,oa,frequency,ref_frequency", [
    (.1, "G", "A", .3, .6), (.1, "A", "G", .6, .6),
    (0., "G", "A", .3, .7), (0., "A", "G", .7, .7),
])
def test_preparation_scopes_diagnostics_without_changing_packed_reference_contract(
    tmp_path, extra, ea, oa, frequency, ref_frequency,
):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n1\t10\tA\tG\t.3\n1\t10\tA\tC\t%s\n" % extra)
    frame = pl.DataFrame({"CHR": ["1"], "POS": [10], "EA": [ea], "OA": [oa], "AF": [frequency]})
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
    args = (frame, columns, "AF", str(reference), "AF", mapping, default_policies())
    old_pairs, old_packed, old_metadata = strand.prepare_study_reference_matches(*args)

    pairs, packed, metadata = strand.prepare_study_reference_matches(*args, diagnostic_require_biallelic=False)

    assert old_pairs.is_empty()
    assert pairs["__study_frequency"].to_list() == [frequency]
    assert pairs["__expected_effect_frequency"].to_list() == pytest.approx([frequency])
    assert packed.equals(old_packed)
    assert metadata["contract"] == old_metadata["contract"]
    assert old_metadata["diagnostic_statistics"]["excluded_non_biallelic_rows"] == 1
    assert metadata["diagnostic_statistics"]["excluded_non_biallelic_rows"] == 0
    candidates = strand.unpack_prepared_reference_matches(packed, SOURCE_INPUT_ROW_COLUMN)
    assert candidates["__strand_biallelic"].to_list() == [False]
    assert candidates["__strand_ref_af"].to_list() == pytest.approx([ref_frequency])
    assert set(metadata["duplicate_statistics"]) == {
        "reference_duplicate_groups", "reference_exact_duplicate_groups", "reference_non_identical_duplicate_groups",
        "reference_duplicate_groups_discarded", "reference_exact_duplicate_rows", "reference_non_identical_duplicate_rows",
        "raw_frequency_screen",
    }


@pytest.mark.parametrize("extra", [None, .9])
def test_incomplete_or_invalid_complete_site_does_not_guess_reference_frequency(tmp_path, extra):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n1\t10\tA\tG\t.3\n1\t10\tA\tC\t%s\n"
                         % ("NA" if extra is None else extra))
    frame = pl.DataFrame({"CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"], "AF": [.7]})
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}

    pairs, packed, _ = strand.prepare_study_reference_matches(frame, columns, "AF", str(reference), "AF", mapping,
        default_policies(), diagnostic_require_biallelic=False)

    assert pairs.is_empty()
    candidates = strand.unpack_prepared_reference_matches(packed, SOURCE_INPUT_ROW_COLUMN)
    assert candidates["__strand_ref_af"].to_list() == [None]


@pytest.mark.parametrize("role", ["chr_col", "pos_col", "ea_col", "oa_col", "eaf_col"])
def test_new_private_ref_frequency_cannot_shadow_a_mapped_input_column(tmp_path, role):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n1\t10\tA\tG\t.3\n1\t10\tA\tC\t.1\n")
    frame = pl.DataFrame({"CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"], "AF": [.2]})
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    frame = frame.rename({columns[role]: "__strand_ref_af"}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns[role] = "__strand_ref_af"
    before = frame.clone()
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}

    pairs, packed, _ = strand.prepare_study_reference_matches(frame, columns, columns["eaf_col"], str(reference),
        "AF", mapping, default_policies(), diagnostic_require_biallelic=False)

    assert pairs["__study_frequency"].to_list() == [.2]
    assert pairs["__expected_effect_frequency"].to_list() == pytest.approx([.6])
    assert frame.equals(before)
    candidates = strand.unpack_prepared_reference_matches(packed, SOURCE_INPUT_ROW_COLUMN)
    assert candidates["__strand_ref_af"].to_list() == pytest.approx([.6])


@pytest.mark.parametrize("acceptance", ["native_not_required", "native_passed", "external", "replacement"])
def test_recorded_dataset_frequency_acceptance_skips_unused_worker_panel_screen(tmp_path, monkeypatch, acceptance):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n1\t10\tA\tG\t.2\n1\t20\tA\tG\t.8\n")
    frame = pl.DataFrame({"CHR": ["1", "1"], "POS": [10, 20], "EA": ["G", "G"],
                          "OA": ["A", "A"], "AF": [.2, .8]})
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
    decision = {"eaf_is_maf": False, "strand": "forward"}
    if acceptance.startswith("native_"):
        decision["eaf_reference_validation"] = {"status": acceptance.removeprefix("native_")}
    else:
        decision["eaf_" + acceptance + "_validation"] = {"status": "passed", "contract": {"approved": True}}
        columns["eaf_col"] = None
    monkeypatch.setattr(strand, "frequency_maf_screen",
                        lambda *a, **k: pytest.fail("worker repeated a raw-panel frequency screen"))

    result, qc, _ = strand.harmonise_strand_orientation("1", frame, columns, str(reference), "AF", mapping,
                                                       study_decision=decision, policies=default_policies())

    assert result.select(frame.columns).equals(frame)
    assert qc["reference_matches_reused"] is False
    assert qc["actions"] == {"forward": 2}


def test_scoped_diagnostic_statistics_count_unsupported_dna_in_fixed_input_rows(tmp_path):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n1\t10\tA\tG\t.2\n1\t20\tA\t*\t.3\n")
    frame = pl.DataFrame({"CHR": ["1", "1"], "POS": [10, 20], "EA": ["G", "*"],
                          "OA": ["A", "A"], "AF": [.2, .3]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}

    pairs, packed, metadata = strand.prepare_study_reference_matches(frame, columns, "AF", str(reference), "AF", mapping,
        default_policies(), diagnostic_require_biallelic=False)

    assert pairs.height == 1
    assert packed.height == 2
    assert metadata["diagnostic_statistics"]["excluded_non_biallelic_rows"] == 0
    assert metadata["diagnostic_statistics"]["excluded_orientation_rows"] == 1
    assert metadata["diagnostic_statistics"]["excluded_counts"] == {
        "unmatched_rows": 0, "palindromic_rows": 0, "ambiguous_rows": 0,
        "unsupported_allele_rows": 1, "non_biallelic_rows": 0,
    }


@pytest.mark.parametrize("require_biallelic,expected_excluded,expected_multi,expected_eligible", [
    (False, 4, 0, 2), (True, 5, 1, 1),
])
def test_diagnostic_exclusion_reasons_are_disjoint_and_reconcile_to_source_rows(
    tmp_path, require_biallelic, expected_excluded, expected_multi, expected_eligible,
):
    reference = tmp_path / "panel.tsv"
    reference.write_text("CHROM\tPOS\tREF\tALT\tAF\n"
        "1\t10\tA\tG\t.2\n1\t20\tA\t*\t.3\n1\t30\tA\tT\t.4\n"
        "1\t40\tA\tG\t.2\n1\t40\tT\tC\t.4\n1\t60\tA\tG\t.3\n1\t60\tA\tC\t.1\n")
    frame = pl.DataFrame({"CHR": ["1"] * 6, "POS": [10, 20, 30, 40, 50, 60],
                          "EA": ["G", "*", "T", "G", "C", "G"], "OA": ["A"] * 6,
                          "AF": [.2, .3, .4, .2, .2, .3]}).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF"}
    mapping = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}

    pairs, packed, metadata = strand.prepare_study_reference_matches(frame, columns, "AF", str(reference), "AF", mapping,
        default_policies(), diagnostic_require_biallelic=require_biallelic)

    statistics = metadata["diagnostic_statistics"]
    assert statistics["excluded_orientation_rows"] == expected_excluded
    assert statistics["excluded_counts"] == {
        "unmatched_rows": 1, "palindromic_rows": 1, "ambiguous_rows": 1,
        "unsupported_allele_rows": 1, "non_biallelic_rows": expected_multi,
    }
    assert sum(statistics["excluded_counts"].values()) == expected_excluded
    assert pairs.height == expected_eligible
    assert packed.height == frame.height == expected_excluded + expected_eligible
    assert statistics["excluded_non_biallelic_rows"] == (3 if require_biallelic else 0)


def _native_multiallelic_swap_fixture(tmp_path, *, action="reject", replacement=False):
    from postgwas.modules.harmonisation import allele_frequency as af
    from postgwas.modules.harmonisation.shared.runtime import NullStepContext
    from postgwas.modules.harmonisation.study_properties import _decide_eaf_is_maf
    from test_harmonisation_external_af_preflight import fixture as external_fixture

    raw = [(.05 if replacement else .55) + i % 20 / 100 for i in range(1200)]
    def add_population_alternates(source, default, external):
        for path in (default, external):
            with path.open("a") as handle:
                handle.writelines("%s\t%d\tA\tG\t.1\n" % (1 + i % 2, 100 + i) for i in range(1200))
        return source
    parts = external_fixture(tmp_path, study_swap=True, values=raw, default_values=raw,
        modify=add_population_alternates, policies={"strand.ambiguous_action": action})
    values = [(.15 + i % 20 / 100) if replacement else (.35 - i % 20 / 100) for i in range(1200)]
    source = parts[0].with_columns(pl.Series("INPUT_AF", values))
    columns = dict(parts[1], eaf_col="INPUT_AF")
    decisions = {"strand": "forward", "effect_type": "beta"}
    _decide_eaf_is_maf(source, columns, parts[4], NullStepContext(), decisions, {})
    selected, prepared = af.validate_study_allele_frequency(source, columns, decisions, parts[3],
        {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}, policies=parts[4],
        eaf_replacement_file=parts[3]["1"]["user_eaf_file"] if replacement else None,
        eaf_replacement_column="AF" if replacement else None)
    return parts, raw, source, columns, selected, prepared


@pytest.mark.parametrize("action", ["reject", "fail"])
def test_accepted_native_multiallelic_ref_eaf_is_not_inverted_as_an_alt_frequency(tmp_path, monkeypatch, action):
    from postgwas.config import load_configuration
    from postgwas.modules.harmonisation import allele_frequency as af
    from postgwas.modules.harmonisation.rejects import RejectCollector
    from postgwas.modules.harmonisation.shared.runtime import NullStepContext
    from test_harmonisation_external_af_preflight import MAPPING

    parts, raw, source, columns, selected, prepared = _native_multiallelic_swap_fixture(tmp_path, action=action)
    evidence = selected["eaf_reference_validation"]
    assert evidence["status"] == "passed"
    assert evidence["comparable_variants"] == evidence["study_variant_count"] == 1200
    assert evidence["mean_eaf_error"] == pytest.approx(0.)
    assert prepared["INPUT_AF"].equals(source["INPUT_AF"])
    rejects = RejectCollector(source, None, str(tmp_path / "rejected.tsv"), "\t", compress=False)
    context = NullStepContext(rows_in=600)
    monkeypatch.setattr(strand, "frequency_maf_screen", lambda *a, **k: pytest.fail("worker repeated AF screening"))
    arguments = dict(chromosome="1", df=prepared.filter(pl.col("CHR") == "1"), sample_column_dict=dict(columns),
        output_layout=dict(load_configuration().modules.harmonisation.output_layout.root), output_delimiter="\t",
        external_eaf_colmap=MAPPING, default_eaf_file=parts[3]["1"]["default_eaf_file"],
        default_eaf_column="AF", default_eaf_colmap=MAPPING, policies=parts[4], study_decision=selected,
        rejects=rejects, step_context=context)
    if action == "fail":
        with pytest.raises(strand.StrandOrientationError, match="(?i)native.*multiallelic.*frequency"):
            af.harmonise_allele_frequency(**arguments)
    else:
        # Every chromosome row has this unsupported native swap, so the
        # existing empty-chromosome gate stops after row rejection.
        with pytest.raises(af.AlleleFrequencyError, match="No variants remain"):
            af.harmonise_allele_frequency(**arguments)
    reason = "native_multiallelic_eaf_swap_unsupported"
    if action == "reject":
        assert context.extra["strand_orientation"][reason] == 600
    assert rejects.counts() == {reason: 600}
    rejected = rejects.frame()
    assert rejected["reject_reason"].unique().to_list() == [reason]
    assert rejected["INPUT_AF"].cast(pl.Float64).to_list() == pytest.approx(source.filter(pl.col("CHR") == "1")["INPUT_AF"].to_list())
    assert rejected["BETA"].cast(pl.Float64).to_list() == [.2] * 600
    assert prepared["INPUT_AF"].equals(source["INPUT_AF"])


@pytest.mark.parametrize("role", ["external", "replacement"])
def test_approved_reference_frequency_sources_keep_complete_site_multiallelic_swaps(tmp_path, monkeypatch, role):
    from postgwas.config import load_configuration
    from postgwas.modules.harmonisation import allele_frequency as af
    from test_harmonisation_external_af_preflight import MAPPING, adopt, validate

    parts, raw, source, columns, selected, prepared = _native_multiallelic_swap_fixture(
        tmp_path, replacement=role == "replacement")
    if role == "external":
        selected = validate(parts)
        monkeypatch.setattr(strand, "frequency_maf_screen", lambda *a, **k: pytest.fail("worker repeated AF screening"))
        result, qc, selected_columns = adopt(parts, selected, "1")
    else:
        selected["eaf_replacement"]["resolved_files"] = {c: r["user_eaf_file"] for c, r in parts[3].items()}
        for c, path in parts[5].paths.items():
            prepared.filter(pl.col("CHR") == c).write_parquet(path)
        selected = af.validate_replacement_allele_frequency(parts[5], columns, selected, parts[3], MAPPING, policies=parts[4])
        monkeypatch.setattr(strand, "frequency_maf_screen", lambda *a, **k: pytest.fail("worker repeated AF screening"))
        result, qc, selected_columns = af.harmonise_allele_frequency("1", pl.read_parquet(parts[5].paths["1"]),
            dict(columns), dict(load_configuration().modules.harmonisation.output_layout.root), "\t",
            eaffile=parts[3]["1"]["user_eaf_file"], external_eaf_colmap=MAPPING,
            default_eaf_file=parts[3]["1"]["default_eaf_file"], default_eaf_column="AF", default_eaf_colmap=MAPPING,
            policies=parts[4], study_decision=selected)
        assert result["INPUT_AF"].to_list() == source.filter(pl.col("CHR") == "1")["INPUT_AF"].to_list()
    assert result.height == 600
    assert result[selected_columns["eaf_col"]].to_list() == pytest.approx(raw[::2])
    assert qc["strand_orientation"].get("native_multiallelic_eaf_swap_unsupported", 0) == 0
