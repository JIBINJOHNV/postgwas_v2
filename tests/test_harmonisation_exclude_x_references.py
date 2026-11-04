"""Excluded X rows and references cannot change retained concordance results."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.concordance.analysis import (
    INPUT_SOURCE_ROW_COLUMN,
    PARTITION_COLUMN,
    VCF_SOURCE_ROW_COLUMN,
    compare_input_to_vcf,
    prepare_input_table,
)
from postgwas.modules.harmonisation.concordance.service import (
    _compare_staged_partitions,
)
from postgwas.modules.harmonisation.external_reference_staging import (
    stage_shared_external_reference_files,
)
from postgwas.modules.harmonisation.imputation_quality import (
    resolve_dataset_info_score_type,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.resource_preflight import (
    ResourcePreflightError,
    _external_frequency_screen,
    validate_harmonisation_resource_maps,
)
from postgwas.modules.harmonisation.resource_paths import (
    validate_external_resource_spec,
)


MAPPING = {
    "chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab",
}
VCF_CONFIG = {
    "csq_prefixes_by_build": {
        build: {"vcf": None, "gff": None, "fasta": None}
        for build in ("GRCh37", "GRCh38")
    },
    "external_frequency_columns": ["CHROM", "POS", "REF", "ALT", "INFO/EUR"],
}


def _row(**updates):
    values = {
        "dataset_id": "study", "input_file": Path("input.tsv"), "delimiter": "tab",
        "chromosome_column": "CHR", "position_column": "BP",
        "chromosome_position_column": None, "effect_allele_column": "EA",
        "other_allele_column": "OA", "effect_allele_frequency_column": "EAF",
        "effect_column": "BETA", "standard_error_column": "SE", "z_score_column": "Z",
        "p_value_column": "P", "effect_type": "beta", "p_value_type": "raw",
        "external_eaf_file": None, "external_eaf_column": None,
    }
    values.update(updates)
    return SimpleNamespace(**values)


def _resource_maps(tmp_path):
    """Write the two-chromosome resource bundle used by structural preflight."""
    def write(name, text):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    annotation = write(
        "genes.gff3",
        "##gff-version 3\n1\tsource\tgene\t1\t10\t.\t+\t.\tID=gene1\n"
        "2\tsource\tgene\t1\t10\t.\t+\t.\tID=gene2\n",
    )
    chain = write(
        "source_to_target.chain",
        "chain 1 1 100 + 0 100 1 100 + 0 100 1\n100\n\n"
        "chain 1 2 100 + 0 100 2 100 + 0 100 2\n100\n",
    )
    maps = {}
    for chromosome in ("1", "2"):
        reference = write("reference_chr%s.vcf.gz" % chromosome, "reference\n")
        write("reference_chr%s.vcf.gz.tbi" % chromosome, "index\n")
        fasta = write("reference_chr%s.fa" % chromosome, ">%s\nA\n" % chromosome)
        write("reference_chr%s.fa.fai" % chromosome, "%s\t100\t3\t100\t101\n" % chromosome)
        eaf = write(
            "default_eaf_chr%s.tsv" % chromosome,
            "CHROM\tPOS\tREF\tALT\tEUR\n%s\t10\tA\tG\t0.2\n" % chromosome,
        )
        maps[chromosome] = {
            "user_eaf_file": None, "user_eaf_column": None,
            "user_info_file": None, "user_info_column": None,
            "default_eaf_file": eaf, "default_eaf_reference_column": "EUR",
            "default_comparison_af_file": reference, "default_comparison_af_column": "EUR",
            "dbsnp_file": reference, "genome_fasta_file": fasta, "target_fasta": fasta,
            "annot_path": annotation, "chain_file": chain,
        }
    return maps


def _mock_bcftools(arguments, _purpose, **_kwargs):
    if arguments[1:3] == ["index", "--stats"]:
        chromosome = Path(arguments[-1]).name.split("chr", 1)[1].split(".", 1)[0]
        return "%s\t100\t1\n" % chromosome
    if arguments[1:3] == ["view", "--header-only"]:
        return "##INFO=<ID=EUR,Number=A,Type=Float,Description=frequency>"
    raise AssertionError("unexpected command: %r" % arguments)


def _policies(**overrides):
    return default_policies().with_overrides({
        "chromosome.exclude_x": True,
        **overrides,
    })


def _validate(maps, policies):
    return validate_harmonisation_resource_maps(
        maps, grch_version="GRCh37", bcftools="bcftools",
        vcf_config=VCF_CONFIG, default_eaf_colmap=MAPPING,
        external_eaf_colmap=MAPPING, external_info_colmap=MAPPING,
        policies=policies, require_default_eaf=True,
    )


def _input(chromosomes):
    size = len(chromosomes)
    return pl.DataFrame({
        "CHR": chromosomes,
        "BP": list(range(101, 101 + size)),
        "EA": ["G"] * size, "OA": ["A"] * size,
        "BETA": [0.2] * size, "SE": [0.1] * size,
        "Z": [2.0] * size, "EAF": [0.2] * size,
        "P": [0.01] * size,
    })


def _vcf(chromosomes):
    size = len(chromosomes)
    return pl.DataFrame({
        "CHROM": chromosomes,
        "POS": list(range(101, 101 + size)),
        "ID": ["rs%s" % number for number in range(size)],
        "ALT": ["G"] * size, "REF": ["A"] * size,
        "ES": [0.2] * size, "SE": [0.1] * size,
        "EZ": [2.0] * size, "AF": [0.2] * size,
        "LP": [2.0] * size,
    })


def _compare(input_frame, vcf_frame, *, policies, row=None, **kwargs):
    return compare_input_to_vcf(
        input_frame, vcf_frame, row or _row(), effect_type="beta",
        p_value_type="raw", eaf_is_maf=False,
        settings=load_configuration().modules.harmonisation.concordance_validation,
        policies=policies, **kwargs,
    )


def test_missing_x_resources_are_ignored_but_retained_resources_are_checked(
    tmp_path, monkeypatch,
):
    maps = _resource_maps(tmp_path)
    maps["X"] = {key: str(tmp_path / "missing_x" / key) for key in maps["1"]}
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        _mock_bcftools,
    )

    summary = _validate(maps, _policies())

    assert summary["status"] == "passed"
    assert summary["chromosomes"] == ["1", "2"]
    assert set(summary["resources_by_chromosome"]) == {"1", "2"}
    assert all("missing_x" not in Path(entry["path"]).parts for entry in summary["file_identities"])

    Path(maps["2"]["default_eaf_file"]).unlink()
    with pytest.raises(ResourcePreflightError) as caught:
        _validate(maps, _policies())
    assert caught.value.issues[0]["chromosomes"] == ["2"]
    assert "default allele-frequency table" in str(caught.value)


def test_default_resource_preflight_still_requires_x(tmp_path, monkeypatch):
    maps = _resource_maps(tmp_path)
    maps["X"] = {key: str(tmp_path / "missing_x" / key) for key in maps["1"]}
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.resource_preflight.run_checked_command",
        _mock_bcftools,
    )

    with pytest.raises(ResourcePreflightError) as caught:
        _validate(maps, default_policies())

    assert all(issue["chromosomes"] == ["X"] for issue in caught.value.issues)


def test_shared_external_frequency_screen_uses_only_retained_chromosomes(tmp_path):
    source = tmp_path / "shared_af.tsv"
    source.write_text(
        "CHROM\tPOS\tALT\tREF\tAF\n1\t101\tG\tA\t0.2\nX\t102\tG\tA\t0.8\n",
        encoding="utf-8",
    )

    screened = _external_frequency_screen(source, "AF", MAPPING, _policies())

    assert screened["rows"] == 1
    assert screened["usable"] == 1
    assert screened["maf_like"] is True


@pytest.mark.parametrize("exclude_x,expected_count", [(False, 8), (True, 1)])
def test_concordance_scope_normalizes_x_aliases_and_preserves_source_ids(
    exclude_x, expected_count,
):
    chromosomes = ["X", "chrX", "23", "25", "XY", "PAR1", "PAR2", "1"]
    source = _input(chromosomes)
    original = source.clone()

    result = _compare(
        source, _vcf(chromosomes),
        policies=_policies(**{"chromosome.exclude_x": exclude_x}),
    )

    assert result.status == "PASS"
    assert result.summary["input_rows"] == expected_count
    assert result.summary["vcf_records"] == expected_count
    assert result.summary["matched_variants"] == expected_count
    assert result.metric_summary["effect"]["checked"] == expected_count
    assert result.input_only.is_empty()
    assert result.vcf_only.is_empty()
    assert source.equals(original)
    if exclude_x:
        assert result.matched["input_row"].to_list() == [8]
        assert result.matched["vcf_row"].to_list() == [8]
        assert result.matched["input_beta"].to_list() == [0.2]


def test_excluded_x_does_not_enter_pvalue_validation_or_coordinate_agreement():
    source = _input(["X", "1"]).with_columns(
        pl.Series("CP", ["X:999", "1:102"]),
        pl.Series("P", [2.0, 0.01]),
    )
    prepared = prepare_input_table(
        source,
        _row(chromosome_column=None, chromosome_position_column="CP"),
        effect_type="beta", p_value_type="raw", eaf_is_maf=False,
        settings=load_configuration().modules.harmonisation.concordance_validation,
        policies=_policies(**{"pvalue.out_of_range": "fail"}),
    )

    assert prepared["input_row"].to_list() == [2]
    assert prepared["input_chrom"].to_list() == ["1"]
    assert prepared["input_pos"].to_list() == [102]


def test_combined_x_without_position_is_excluded_before_scientific_validation():
    source = _input(["X", "1"]).drop("CHR", "BP").with_columns(
        pl.Series("CP", ["chrX:", "1:102"]),
        pl.Series("P", [2.0, 0.01]),
    )
    prepared = prepare_input_table(
        source,
        _row(chromosome_column=None, position_column=None, chromosome_position_column="CP"),
        effect_type="beta", p_value_type="raw", eaf_is_maf=False,
        settings=load_configuration().modules.harmonisation.concordance_validation,
        policies=_policies(**{"pvalue.out_of_range": "fail"}),
    )

    assert prepared["input_row"].to_list() == [2]
    assert prepared["input_pos"].to_list() == [102]


def test_excluded_x_archived_strand_errors_do_not_affect_retained_concordance():
    actions = pl.DataFrame({
        "CHROM": ["X", "1"], "POS": [101, 102],
        "REF": ["A", "A"], "ALT": ["G", "G"],
        "strand_action": ["invalid_x_action", "forward"],
    })

    result = _compare(
        _input(["X", "1"]), _vcf(["X", "1"]),
        policies=_policies(), strand_action_frame=actions,
    )

    assert result.status == "PASS"
    assert result.summary["matched_variants"] == 1


def test_excluded_x_reference_duplicates_do_not_affect_retained_concordance():
    external = pl.DataFrame({
        "CHROM": ["X", "X", "1"], "POS": [101, 101, 102],
        "REF": ["A", "A", "A"], "ALT": ["G", "G", "G"],
        "EAF": [0.1, 0.9, 0.2],
    })

    result = _compare(
        _input(["X", "1"]), _vcf(["X", "1"]),
        row=_row(external_eaf_column="EAF"),
        policies=_policies(**{"external_reference.non_identical_duplicate_action": "fail"}),
        external_eaf_frame=external,
        external_eaf_mapping=load_configuration().modules.harmonisation.external_eaf_mapping,
    )

    assert result.status == "PASS"
    assert result.metric_summary["allele_frequency"]["concordant"] == 1


def test_staged_concordance_never_resolves_an_excluded_x_reference(tmp_path):
    source = _input(["X", "1"]).with_row_index(INPUT_SOURCE_ROW_COLUMN, offset=1)
    source = source.with_columns(pl.col("CHR").alias(PARTITION_COLUMN))
    vcf = _vcf(["X", "1"]).with_row_index(VCF_SOURCE_ROW_COLUMN, offset=1)
    vcf = vcf.with_columns(pl.col("CHROM").alias(PARTITION_COLUMN))
    input_path, vcf_path = tmp_path / "input.parquet", tmp_path / "vcf.parquet"
    source.write_parquet(input_path)
    vcf.write_parquet(vcf_path)
    (tmp_path / "af_1.tsv").write_text(
        "CHROM\tPOS\tREF\tALT\tEAF\n1\t102\tA\tG\t0.2\n", encoding="utf-8",
    )
    config = load_configuration().modules.harmonisation
    row = _row(
        effect_allele_frequency_column=None,
        external_eaf_file=str(tmp_path / "af_{chromosome}.tsv"),
        external_eaf_column="EAF",
    )

    result = _compare_staged_partitions(
        input_path=input_path, vcf_path=vcf_path, duplicate_path=None,
        external_eaf_path=None, workspace=tmp_path, row=row,
        effect_type="beta", p_value_type="raw", eaf_is_maf=False,
        settings=config.concordance_validation, policies=_policies(),
        external_eaf_mapping=config.external_eaf_mapping, external_eaf_build="GRCh37",
        log=SimpleNamespace(info=lambda *_args: None),
    )

    assert result.status == "PASS"
    assert result.summary["input_rows"] == 1
    assert result.summary["vcf_records"] == 1
    assert result.matched.collect()["input_row"].to_list() == [2]
    assert not (tmp_path / "af_X.tsv").exists()


def test_checkpoint_template_dependencies_exclude_only_the_declared_x_token(tmp_path):
    for name in ("panel_X_build37.tsv", "panel_1_build37.tsv", "panel_1_buildX.tsv"):
        (tmp_path / name).write_text("reference\n", encoding="utf-8")
    template = str(tmp_path / "panel_{chromosome}_build{build}.tsv")

    files = validate_external_resource_spec(template, excluded_chromosomes=("X",))

    assert {path.name for path in files} == {"panel_1_build37.tsv", "panel_1_buildX.tsv"}
    assert len(validate_external_resource_spec(template)) == 3
    assert validate_external_resource_spec(
        tmp_path / "panel_X_build37.tsv", excluded_chromosomes=("X",),
    ) == ((tmp_path / "panel_X_build37.tsv").resolve(),)


def test_checkpoint_excluded_only_template_is_not_a_missing_individual_reference(tmp_path):
    (tmp_path / "panel_X.tsv").write_text("reference\n", encoding="utf-8")

    assert validate_external_resource_spec(
        str(tmp_path / "panel_{chromosome}.tsv"), excluded_chromosomes=("X",),
    ) == ()


def test_streamed_shared_reference_ignores_x_bad_values_with_positive_inference(tmp_path):
    source = tmp_path / "shared_reference.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\tEAF\tINFO\n"
        "chr1\t101\tA\tG\t0.2\t0.8\n"
        "chrX\tbad_position\tA\tG\tbad_frequency\tbad_info\n"
        "chr2\t102\tC\tT\t0.3\t0.9\n",
        encoding="utf-8",
    )
    original = source.read_bytes()
    config = load_configuration().modules.harmonisation
    settings = config.external_reference_staging.model_dump()
    settings["batch_rows"] = 1
    resource_maps = {
        chromosome: {"user_eaf_file": str(source), "user_info_file": str(source)}
        for chromosome in ("1", "2")
    }

    staged, summary = stage_shared_external_reference_files(
        resource_maps, chromosomes=("1", "2"), dataset_id="study",
        output_directory=tmp_path / "staged", output_layout=dict(config.output_layout.root),
        user_eaf_specification=str(source), user_eaf_column="EAF",
        external_eaf_mapping=MAPPING, user_info_specification=str(source),
        user_info_column="INFO", external_info_mapping=MAPPING, settings=settings,
        policies=_policies(**{"input.schema_inference_rows": 1}),
    )

    assert summary["source_files_read"] == 1
    assert summary["resources"][0]["rows_read"] == 3
    assert summary["resources"][0]["rows_by_chromosome"] == {"1": 1, "2": 1}
    assert summary["resources"][0]["rows_not_in_observed_chromosomes"] == 1
    for chromosome, position, eaf, info in (("1", 101, 0.2, 0.8), ("2", 102, 0.3, 0.9)):
        assert staged[chromosome]["user_eaf_file"] == staged[chromosome]["user_info_file"]
        partition = pl.read_parquet(staged[chromosome]["user_eaf_file"])
        assert partition["CHROM"].to_list() == ["chr" + chromosome]
        assert partition["POS"].cast(pl.Int64).to_list() == [position]
        assert partition["EAF"].cast(pl.Float64).to_list() == [eaf]
        assert partition["INFO"].cast(pl.Float64).to_list() == [info]
    assert source.read_bytes() == original


def test_streamed_external_info_ignores_x_bad_score_with_positive_inference(tmp_path):
    source = tmp_path / "shared_info.tsv"
    source.write_text(
        "CHROM\tINFO\nchr1\t0.8\nchrX\tbad_info\nchr2\t0.9\n",
        encoding="utf-8",
    )
    original = source.read_bytes()
    settings = load_configuration().modules.harmonisation.external_reference_staging.model_dump()
    settings["batch_rows"] = 1

    decision = resolve_dataset_info_score_type(
        df=None, sample_column_dict={"info_source": "external"},
        resource_maps={
            "1": {"user_info_file": str(source)},
            "2": {"user_info_file": str(source)},
        },
        chromosomes=("1", "2"), external_info_column="INFO",
        external_info_colmap=MAPPING, external_reference_staging=settings,
        policies=_policies(**{"input.schema_inference_rows": 1}),
    )

    assert decision["resolved_type"] == "standard_info"
    assert decision["source_rows"] == 2
    assert decision["finite_values"] == 2
    assert decision["unusable_values"] == 0
    assert decision["minimum"] == 0.8
    assert decision["maximum"] == 0.9
    assert source.read_bytes() == original
