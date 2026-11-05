"""The selected AF source must pass a row-bound whole-study gate before adoption."""
from copy import deepcopy
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.chromosome_partition import ChromosomePartitions
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import SOURCE_INPUT_ROW_COLUMN


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}


def test_prepared_work_fingerprint_preserves_each_callers_error_type(tmp_path):
    from postgwas.core.input_validation import InputValidationSession
    from postgwas.modules.harmonisation import service
    from postgwas.modules.harmonisation.shared import runtime
    path = tmp_path / "work.parquet"
    with InputValidationSession(pipeline=False):
        for error_type in (service.PipelineError, af.StudyAlleleFrequencyError):
            with pytest.raises(error_type, match="work.parquet"):
                runtime.prepared_work_fingerprint(path, error_type=error_type)
        path.write_bytes(b"validated work content")
        first = runtime.prepared_work_fingerprint(path, error_type=service.PipelineError)
        assert runtime.prepared_work_fingerprint(path, error_type=af.StudyAlleleFrequencyError) == first
        path.unlink()
        with pytest.raises(af.StudyAlleleFrequencyError, match="work.parquet"):
            runtime.prepared_work_fingerprint(path, error_type=af.StudyAlleleFrequencyError)


def fixture(tmp_path, *, same=False, chromosomes=("1",), replacement_values=None, modify=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    original_b = [.8, .7, .6, .1, .2, .3] * 2
    alleles = [("C", "A"), ("A", "C"), ("G", "T"), ("T", "G")] * 3
    chrom = [chromosomes[i % len(chromosomes)] for i in range(12)]
    reference_b = [1 - b if i % 4 in (1, 3) else b for i, b in enumerate(original_b)]
    source = pl.DataFrame({"CHR": chrom, "POS": list(range(100, 112)),
        "EA": [x[0] for x in alleles], "OA": [x[1] for x in alleles],
        "MAF": [min(b, 1-b) for b in original_b], "BETA": [.2] * 12, "Z": [2.] * 12,
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN)
    columns = {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": "MAF", "beta_or_col": "BETA", "imp_z_col": "Z",
        "gwas_outputname": "study", "output_folder": str(tmp_path / "out")}
    reference = pl.DataFrame({"CHROM": chrom, "POS": list(range(100,112)),
        "REF": ["A"] * 12, "ALT": ["C"] * 12, "AF": reference_b})
    default = tmp_path / "default.tsv"
    reference.write_csv(default, separator="\t")
    replacement = default if same else tmp_path / "replacement.tsv"
    if not same:
        reference.with_columns(pl.Series("AF", reference_b if replacement_values is None else replacement_values)).write_csv(replacement, separator="\t")
    if modify:
        source = modify(source, default, replacement)
    policies = default_policies().with_overrides({"eaf.maf_reference_min_overlap": 2})
    maps = {c: {"default_eaf_file": str(default), "default_eaf_reference_column": "AF",
        "user_eaf_file": None if same else str(replacement)} for c in chromosomes}
    decisions, work = af.validate_study_allele_frequency(source, columns,
        {"eaf_column": "MAF", "eaf_is_maf": True, "effect_type": "beta", "strand": "forward"},
        maps, MAPPING, policies=policies, eaf_replacement_file=str(replacement), eaf_replacement_column="AF")
    decisions["eaf_replacement"]["resolved_files"] = {c: str(replacement) for c in chromosomes}
    paths, rows = {}, {}
    for c in chromosomes:
        path = tmp_path / ("work_" + c + ".parquet")
        subset = work.filter(pl.col("CHR") == c)
        subset.write_parquet(path)
        paths[c], rows[c] = str(path), subset.height
    return source, columns, decisions, maps, policies, ChromosomePartitions(paths, rows), original_b


def validate(parts):
    _, columns, decisions, maps, policies, partitions, _ = parts
    return af.validate_replacement_allele_frequency(partitions, columns, decisions, maps, MAPPING, policies=policies)


def adopt(parts, decisions, frame=None, *, mapping=MAPPING, policies=None):
    _, columns, _, maps, default_pol, partitions, _ = parts
    c = next(iter(partitions.paths))
    return af.harmonise_allele_frequency(c, pl.read_parquet(partitions.paths[c]) if frame is None else frame,
        dict(columns), dict(load_configuration().modules.harmonisation.output_layout.root), "\t",
        eaffile=maps[c]["user_eaf_file"], external_eaf_colmap=mapping,
        default_eaf_file=maps[c]["default_eaf_file"], default_eaf_column="AF", default_eaf_colmap=MAPPING,
        policies=default_pol if policies is None else policies, study_decision=decisions)


def test_nonselected_replacement_gate_never_reads_partitions_or_candidate(monkeypatch):
    monkeypatch.setattr(pl, "read_parquet", lambda *a, **k: pytest.fail("candidate was opened"))
    original = {"eaf_selected_source": "study_supplied"}
    assert af.validate_replacement_allele_frequency(None, {}, original, {}, {}) == original


@pytest.mark.parametrize("same", [False, True])
@pytest.mark.parametrize("chromosomes", [("1",), ("1", "2")])
def test_gate_pools_original_ea_orientation_without_mutation_or_repeated_merge(tmp_path, monkeypatch, same, chromosomes):
    parts = fixture(tmp_path, same=same, chromosomes=chromosomes)
    source, columns, original, maps, policies, partitions, original_b = parts
    before = {p: Path(p).read_bytes() for p in partitions.paths.values()}
    calls = []
    merge = af.merge_external_allele_frequencies
    def observed(*args, **kwargs):
        calls.append(kwargs)
        return merge(*args, **kwargs)
    monkeypatch.setattr(af, "merge_external_allele_frequencies", observed)
    result = validate(parts)
    evidence = result["eaf_replacement_validation"]
    assert evidence["status"] == "passed"
    assert evidence["study_variant_count"] == evidence["comparable_variants"] == 12
    assert evidence["major_count"] == 6
    assert evidence["source"]["same_source_evidence_reused"] is same
    assert "evidence_is_independent_of_default" not in evidence["source"]
    assert evidence["correlation"] == pytest.approx(1.)
    assert evidence["mean_absolute_difference"] == pytest.approx(0.)
    assert result["eaf_reference_validation"] == original["eaf_reference_validation"]
    assert "eaf_replacement_validation" not in original
    assert len(calls) == (0 if same else len(chromosomes))
    assert all(call["min_match_fraction"] == 0.0 and call["declared_alt_frequency"] for call in calls)
    assert {p: Path(p).read_bytes() for p in partitions.paths.values()} == before
    assert set(tmp_path.iterdir()) == {Path(p) for p in partitions.paths.values()} | {tmp_path / "default.tsv"} | (set() if same else {tmp_path / "replacement.tsv"})


def test_failed_gate_carries_original_and_new_evidence(tmp_path):
    parts = fixture(tmp_path, replacement_values=[.25] * 12)
    with pytest.raises(af.StudyAlleleFrequencyError, match="replacement.*inconclusive") as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    assert decisions["eaf_reference_validation"] == parts[2]["eaf_reference_validation"]
    assert decisions["eaf_replacement_validation"]["status"] == "inconclusive"
    assert decisions["eaf_replacement_validation"]["failed_requirements"]


def test_provisional_or_status_only_approval_cannot_adopt(tmp_path):
    parts = fixture(tmp_path)
    for approval in (None, {"status": "passed"}):
        decisions = deepcopy(parts[2])
        if approval:
            decisions["eaf_replacement_validation"] = approval
        with pytest.raises(af.AlleleFrequencyError, match="replacement.*validation"):
            adopt(parts, decisions)


def test_passed_adoption_accepts_validated_subset_and_preserves_original_maf(tmp_path):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"]).slice(2, 7)
    result, qc, columns = adopt(parts, decisions, frame)
    assert result["MAF"].to_list() == frame["MAF"].to_list()
    assert result[columns["eaf_col"]].to_list() == pytest.approx([
        1-b if i % 4 in (1,3) else b for i,b in enumerate(parts[6])][2:9])
    assert qc["eaf_replacement_validation"]["status"] == "passed"


@pytest.mark.parametrize("change", ["frequency", "row_id", "alleles", "policy", "mapping", "source", "partition"])
def test_adoption_rejects_stale_or_changed_validation_contract(tmp_path, change):
    parts = fixture(tmp_path)
    decisions = validate(parts)
    frame = pl.read_parquet(parts[5].paths["1"])
    kwargs = {}
    if change == "frequency":
        frame = frame.with_columns(pl.lit(.123).alias("MAF"))
    elif change == "row_id":
        frame = frame.with_columns((pl.col(SOURCE_INPUT_ROW_COLUMN) + 100).alias(SOURCE_INPUT_ROW_COLUMN))
    elif change == "alleles":
        frame = frame.with_columns(pl.lit("T").alias("EA"))
    elif change == "policy":
        kwargs["policies"] = parts[4].with_overrides({"eaf.study_reference_min_correlation": .99})
    elif change == "mapping":
        kwargs["mapping"] = dict(MAPPING, a1="REF", a2="ALT")
    elif change == "source":
        path = tmp_path / "replacement.tsv"
        path.write_text(path.read_text().replace("0.8", "0.81"))
    else:
        frame.with_columns(pl.lit(.123).alias("MAF")).write_parquet(parts[5].paths["1"])
    with pytest.raises(af.AlleleFrequencyError, match="replacement.*validation|replacement.*changed"):
        adopt(parts, decisions, frame, **kwargs)


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -1.])
def test_unchanged_invalid_original_values_remain_provenance(tmp_path, value):
    def modify(source, default, replacement):
        return source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit(value, dtype=pl.Float64)).otherwise(pl.col("MAF")).alias("MAF"))
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["comparable_variants"] == 11
    frame = pl.read_parquet(parts[5].paths["1"]).reverse()
    result, _, columns = adopt(parts, decisions, frame)
    assert result["MAF"].equals(frame["MAF"])
    assert result.height == 12 and result[columns["eaf_col"]].null_count() == 0


@pytest.mark.parametrize("extra", [.1, None, float("nan"), .9])
@pytest.mark.parametrize("swapped", [False, True])
def test_multiallelic_diagnostics_exclude_rows_without_changing_attachment(tmp_path, extra, swapped):
    def modify(source, default, replacement):
        ref = pl.read_csv(replacement, separator="\t", infer_schema_length=0)
        if swapped:
            ref = ref.with_columns(pl.when(pl.col("POS") == "100").then(pl.lit("C")).otherwise(pl.col("REF")).alias("REF"),
                pl.when(pl.col("POS") == "100").then(pl.lit("A")).otherwise(pl.col("ALT")).alias("ALT"),
                pl.when(pl.col("POS") == "100").then(pl.lit("0.2")).otherwise(pl.col("AF")).alias("AF"))
        ref.write_csv(replacement, separator="\t")
        with replacement.open("a") as handle:
            handle.write("1\t100\t%s\tG\t%s\n" % ("C" if swapped else "A", "NA" if extra is None else extra))
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    evidence = decisions["eaf_replacement_validation"]
    assert evidence["study_variant_count"] == 12 and evidence["comparable_variants"] == 11
    assert evidence["excluded_counts"]["replacement_non_biallelic_or_unresolved_rows"] == 1
    result, _, columns = adopt(parts, decisions)
    expected = (None if extra is None or extra != extra or extra == .9 else .7) if swapped else .8
    if expected is None:
        assert result[columns["eaf_col"]][0] is None
    else:
        assert result[columns["eaf_col"]][0] == pytest.approx(expected)
    assert result.height == 12 and result["MAF"].equals(parts[0]["MAF"])


def test_global_half_coverage_does_not_override_worker_eighty_percent_gate(tmp_path):
    def modify(source, default, replacement):
        pl.read_csv(replacement, separator="\t").head(6).write_csv(replacement, separator="\t")
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["observed_overlap_fraction"] == .5
    with pytest.raises(af.AlleleFrequencyError, match="external_min_match_fraction"):
        adopt(parts, decisions)


def test_palindromic_and_unmatched_original_orientation_remain_in_fixed_n(tmp_path):
    def modify(source, default, replacement):
        # Two rows have no diagnostic orientation; other rows still support R2.
        source = source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"))
        reference = pl.read_csv(default, separator="\t").filter(pl.col("POS") != 101)
        reference.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, modify=modify)
    decisions = validate(parts)
    evidence = decisions["eaf_replacement_validation"]
    assert evidence["study_variant_count"] == 12 and evidence["comparable_variants"] == 10
    assert evidence["excluded_counts"]["original_orientation_ineligible_rows"] == 2
    assert pl.read_parquet(parts[5].paths["1"]).height == 12


def test_indel_diagnostic_keys_share_default_reference_canonicalization(tmp_path):
    def modify(source, default, replacement):
        for path in (default, replacement):
            ref = pl.read_csv(path, separator="\t", infer_schema_length=0)
            ref.with_columns(pl.when(pl.col("POS") == "100").then(pl.lit("aac")).otherwise(pl.col("ALT")).alias("ALT"),
                pl.when(pl.col("POS") == "100").then(pl.lit("aa")).otherwise(pl.col("REF")).alias("REF")).write_csv(path, separator="\t")
        return source.with_columns(pl.when(pl.col("POS") == 100).then(pl.lit(" AAC ")).otherwise(pl.col("EA")).alias("EA"),
            pl.when(pl.col("POS") == 100).then(pl.lit(" aa ")).otherwise(pl.col("OA")).alias("OA"))
    parts = fixture(tmp_path, modify=modify)
    result = validate(parts)
    assert result["eaf_replacement_validation"]["comparable_variants"] == 12
    assert parts[0]["EA"][0] == " AAC "


def test_gate_refuses_source_changed_during_numeric_replay(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    original = af.validate_replacement_frequency_pairs
    def mutate(*args, **kwargs):
        result = original(*args, **kwargs)
        path = tmp_path / "replacement.tsv"
        path.write_text(path.read_text().replace("0.8", "0.81"))
        return result
    monkeypatch.setattr(af, "validate_replacement_frequency_pairs", mutate)
    with pytest.raises(af.StudyAlleleFrequencyError, match="changed") as caught:
        validate(parts)
    evidence = caught.value.study_decisions["eaf_replacement_validation"]
    assert evidence["status"] == "invalid" and "contract" not in evidence
    assert evidence["comparable_variants"] == 12


def test_distinct_partition_directories_can_share_a_basename(tmp_path):
    parts = list(fixture(tmp_path, chromosomes=("1", "2")))
    paths = {}
    for chromosome, original in parts[5].paths.items():
        destination = tmp_path / chromosome / "work.parquet"
        destination.parent.mkdir()
        Path(original).rename(destination)
        paths[chromosome] = str(destination)
    parts[5] = ChromosomePartitions(paths, parts[5].rows)
    assert validate(parts)["eaf_replacement_validation"]["status"] == "passed"
    assert list((tmp_path / "1").iterdir()) == [Path(paths["1"])]


def test_inconclusive_evidence_cannot_be_approved_by_changing_status(tmp_path):
    parts = fixture(tmp_path, replacement_values=[.25] * 12)
    with pytest.raises(af.StudyAlleleFrequencyError) as caught:
        validate(parts)
    decisions = caught.value.study_decisions
    decisions["eaf_replacement_validation"]["status"] = "passed"
    with pytest.raises(af.AlleleFrequencyError, match="replacement validation"):
        adopt(parts, decisions)


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_default_source_alias_reuses_prepared_af_without_source_read(tmp_path, monkeypatch, alias_kind):
    alias = tmp_path / "alias.tsv"
    def modify(source, default, replacement):
        getattr(alias, "symlink_to" if alias_kind == "symlink" else "hardlink_to")(default)
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    parts[2]["eaf_replacement"].update(file=str(alias), resolved_files={"1": str(alias)})
    monkeypatch.setattr(af, "read_reference_variant_table", lambda *a, **k: pytest.fail("identical source reread"))
    decisions = validate(parts)
    assert decisions["eaf_replacement_validation"]["source"]["reused_default_chromosomes"] == ["1"]
    assert adopt(parts, decisions)[0].height == 12


@pytest.mark.parametrize("route", ["direct", "activated", "staged"])
@pytest.mark.parametrize("when", ["after_measurement", "after_approval", "during_strand"])
def test_retargeted_alias_never_adopts_unvalidated_frequencies(tmp_path, monkeypatch, route, when):
    from postgwas.modules.harmonisation import service
    parts = list(fixture(tmp_path, same=route != "staged"))
    alias = tmp_path / "alias_GRCh37_1.tsv"
    original = tmp_path / ("replacement.tsv" if route == "staged" else "default.tsv")
    alias.symlink_to(original)
    alternate = tmp_path / "unvalidated.tsv"
    pl.read_csv(original, separator="\t").with_columns(pl.lit(.25).alias("AF")).write_csv(alternate, separator="\t")
    if route == "staged":
        parts[3]["1"]["user_eaf_file"] = str(alias)
    else:
        parts[2]["eaf_replacement"].update(file=str(alias), resolved_files={"1": str(alias)})
    if route == "activated":
        parts[2]["eaf_replacement"]["file"] = str(tmp_path / "alias_{build}_{chromosome}.tsv")
        parts[3], _ = service._activate_eaf_replacement_resources(
            parts[3], parts[2], grch_version="GRCh37", default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING, executables={"bcftools": "must-not-run"},
            vcf_config=load_configuration().modules.harmonisation.vcf_processing.model_dump(), policies=parts[4])
    def retarget():
        alias.unlink()
        alias.symlink_to(alternate)
    if when == "after_measurement":
        measure = af.validate_replacement_frequency_pairs
        def measured(*args, **kwargs):
            result = measure(*args, **kwargs)
            retarget()
            return result
        monkeypatch.setattr(af, "validate_replacement_frequency_pairs", measured)
        with pytest.raises(af.StudyAlleleFrequencyError, match="changed") as caught:
            validate(parts)
        evidence = caught.value.study_decisions["eaf_replacement_validation"]
        assert evidence["status"] == "invalid" and "contract" not in evidence
        assert evidence["comparable_variants"] == 12 and evidence["correlation"] == pytest.approx(1.)
        return
    decisions = validate(parts)
    if when == "after_approval":
        retarget()
        monkeypatch.setattr(af, "harmonise_strand_orientation", lambda *a, **k: pytest.fail("strand ran after alias retarget"))
        with pytest.raises(af.AlleleFrequencyError, match="changed"):
            adopt(parts, decisions)
        return
    orient = af.harmonise_strand_orientation
    def oriented(*args, **kwargs):
        result = orient(*args, **kwargs)
        retarget()
        return result
    monkeypatch.setattr(af, "harmonise_strand_orientation", oriented)
    # Once verified, attachment must use the pinned resolved source, never the
    # new target of a mutable alias passed between orientation and attachment.
    result, _, selected = adopt(parts, decisions)
    expected = pl.read_csv(original, separator="\t")["AF"].to_list()
    assert result[selected["eaf_col"]].to_list() == pytest.approx(expected)


def test_changed_frequency_column_same_file_runs_shared_merge(tmp_path, monkeypatch):
    def modify(source, default, replacement):
        pl.read_csv(default, separator="\t").with_columns(pl.col("AF").alias("AF2")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    parts[2]["eaf_replacement"]["column"] = "AF2"
    calls = []
    original = af.merge_external_allele_frequencies
    def merge(*args, **kwargs):
        calls.append(args[2])
        return original(*args, **kwargs)
    monkeypatch.setattr(af, "merge_external_allele_frequencies", merge)
    assert validate(parts)["eaf_replacement_validation"]["status"] == "passed"
    assert calls == ["AF2"]


def test_changed_allele_mapping_same_file_runs_shared_merge(tmp_path, monkeypatch):
    def modify(source, default, replacement):
        pl.read_csv(default, separator="\t").with_columns(pl.col("ALT").alias("ALT2")).write_csv(default, separator="\t")
        return source
    parts = fixture(tmp_path, same=True, modify=modify)
    calls = []
    original = af.merge_external_allele_frequencies
    def merge(*args, **kwargs):
        calls.append(args[3])
        return original(*args, **kwargs)
    mapping = dict(MAPPING, a1="ALT2")
    monkeypatch.setattr(af, "merge_external_allele_frequencies", merge)
    decisions = af.validate_replacement_allele_frequency(parts[5], parts[1], parts[2], parts[3], mapping, policies=parts[4])
    assert decisions["eaf_replacement_validation"]["status"] == "passed"
    assert calls == [mapping]


def test_unreadable_candidate_preserves_io_failure_and_original_evidence(tmp_path, monkeypatch):
    parts = fixture(tmp_path)
    def unreadable(*args, **kwargs):
        raise PermissionError("selected replacement cannot be read")
    monkeypatch.setattr(af, "read_reference_variant_table", unreadable)
    with pytest.raises(PermissionError, match="selected replacement cannot be read") as caught:
        validate(parts)
    assert caught.value.study_decisions["eaf_reference_validation"] == parts[2]["eaf_reference_validation"]
    assert "eaf_replacement_validation" not in caught.value.study_decisions


def test_chromosome_without_diagnostic_candidates_does_not_shrink_fixed_n(tmp_path):
    def modify(source, default, replacement):
        for path in (default, replacement):
            pl.read_csv(path, separator="\t", infer_schema_length=0).with_columns(
                pl.when(pl.col("CHROM") == "2").then(pl.lit("T")).otherwise(pl.col("ALT")).alias("ALT")
            ).write_csv(path, separator="\t")
        return source.with_columns(
            pl.when(pl.col("CHR") == "2").then(pl.lit("T")).otherwise(pl.col("EA")).alias("EA"),
            pl.when(pl.col("CHR") == "2").then(pl.lit("A")).otherwise(pl.col("OA")).alias("OA"))
    parts = fixture(tmp_path, chromosomes=("1", "2"), modify=modify)
    evidence = validate(parts)["eaf_replacement_validation"]
    assert evidence["comparable_variants"] == 6 and evidence["study_variant_count"] == 12
    assert evidence["observed_overlap_fraction"] == .5


def test_chromosome_template_matches_one_global_whole_file_decision(tmp_path):
    parts = fixture(tmp_path, chromosomes=("1", "2"))
    whole = validate(parts)["eaf_replacement_validation"]
    reference = pl.read_csv(tmp_path / "replacement.tsv", separator="\t", infer_schema_length=0)
    for chromosome in parts[5].paths:
        path = tmp_path / ("replacement_" + chromosome + ".tsv")
        reference.filter(pl.col("CHROM") == chromosome).write_csv(path, separator="\t")
        parts[2]["eaf_replacement"]["resolved_files"][chromosome] = str(path)
        parts[3][chromosome]["user_eaf_file"] = str(path)
    parts[2]["eaf_replacement"]["file"] = str(tmp_path / "replacement_{chromosome}.tsv")
    templated = validate(parts)["eaf_replacement_validation"]
    for key in ("status", "comparable_variants", "major_count", "observed_overlap_fraction", "correlation", "mean_absolute_difference"):
        assert templated[key] == whole[key]


@pytest.mark.parametrize("problem", ["missing", "malformed", "column", "mapping", "wrong_build"])
def test_invalid_selected_sources_never_supply_approval(tmp_path, problem):
    parts = fixture(tmp_path)
    path = tmp_path / "replacement.tsv"
    if problem == "missing":
        path.unlink()
    elif problem == "malformed":
        path.write_text("not a valid table\n")
    elif problem == "column":
        parts[2]["eaf_replacement"]["column"] = "ABSENT_AF"
    elif problem == "mapping":
        with pytest.raises(af.AlleleFrequencyError, match="mapping.*missing|missing.*mapping"):
            af.validate_replacement_allele_frequency(parts[5], parts[1], parts[2], parts[3], {}, policies=parts[4])
        return
    else:
        pl.read_csv(path, separator="\t").with_columns((pl.col("POS") + 100000).alias("POS")).write_csv(path, separator="\t")
    with pytest.raises(af.AlleleFrequencyError) as caught:
        validate(parts)
    if problem == "wrong_build":
        assert caught.value.study_decisions["eaf_replacement_validation"]["comparable_variants"] == 0
    else:
        assert "eaf_replacement_validation" not in caught.value.study_decisions
