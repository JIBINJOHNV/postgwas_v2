"""Activate replacement resources only after the default-panel decision."""

from copy import deepcopy
from pathlib import Path
import shutil
import sys

import polars as pl
import pysam
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import service
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.resource_preflight import ResourcePreflightError
from postgwas.core.pipeline_logging import PipelineLogger


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}


def _table(path, chromosome, column="EUR"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"CHROM\tPOS\tREF\tALT\t{column}\n{chromosome}\t101\tA\tC\t0.2\n"
    )
    return str(path)


def _maps(tmp_path):
    return {
        chrom: {
            "default_eaf_file": _table(tmp_path / f"default_chr{chrom}.tsv", chrom),
            "default_eaf_reference_column": "EUR",
            "user_eaf_file": None,
            "user_eaf_column": None,
        }
        for chrom in ("1", "2")
    }


def _decision(template, column):
    return {
        "eaf_selected_source": "replacement_reference",
        "eaf_input_interpretation": "maf",
        "eaf_reference_validation": {"status": "replacement_selected"},
        "eaf_replacement": {"file": str(template), "column": column},
    }


def _activate(maps, decisions):
    return service._activate_eaf_replacement_resources(
        maps, decisions, grch_version="GRCh37",
        default_eaf_colmap=MAPPING, external_eaf_colmap=MAPPING,
        executables={"bcftools": "must-not-be-run"},
        vcf_config=load_configuration().modules.harmonisation.vcf_processing.model_dump(),
        policies=default_policies(),
    )


def test_retained_internal_eaf_never_opens_replacement(tmp_path):
    maps = _maps(tmp_path)
    original = deepcopy(maps)
    selected, summary = _activate(maps, {"eaf_selected_source": "study"})
    assert selected == original
    assert summary is None


def test_replacement_resolves_inferred_build_and_retains_original_paths(tmp_path):
    maps = _maps(tmp_path)
    original = deepcopy(maps)
    template = tmp_path / "{build}" / "replacement_chr{chromosome}.tsv"
    expected = {
        chrom: _table(tmp_path / "GRCh37" / f"replacement_chr{chrom}.tsv", chrom, "PANEL_AF")
        for chrom in maps
    }
    decisions = _decision(template, "PANEL_AF")
    selected, summary = _activate(maps, decisions)
    assert maps == original
    assert decisions["eaf_replacement"]["resolved_files"] == expected
    assert {chrom: row["user_eaf_file"] for chrom, row in selected.items()} == expected
    assert all(row["user_eaf_column"] == "PANEL_AF" for row in selected.values())
    assert {row["path"] for row in summary["file_identities"]} == set(expected.values())
    assert summary["frequency_vcf_headers_checked"] == 0
    assert summary["table_headers_checked"] == 2


def test_same_default_source_does_not_enter_external_staging_role(tmp_path):
    maps = _maps(tmp_path)
    template = tmp_path / "default_chr{chromosome}.tsv"
    decisions = _decision(template, "EUR")
    selected, summary = _activate(maps, decisions)
    assert all(row["user_eaf_file"] is None for row in selected.values())
    assert summary["table_headers_checked"] == 0
    assert summary["reused_default_chromosomes"] == ["1", "2"]
    assert decisions["eaf_replacement"]["resolved_files"] == {
        chrom: row["default_eaf_file"] for chrom, row in maps.items()
    }


def test_active_replacement_reports_every_missing_chromosome(tmp_path):
    maps = _maps(tmp_path)
    decisions = _decision(tmp_path / "missing_chr{chromosome}.tsv", "AF")
    with pytest.raises(ResourcePreflightError) as caught:
        _activate(maps, decisions)
    assert {c for issue in caught.value.issues for c in issue["chromosomes"]} == {"1", "2"}
    assert all(row["user_eaf_file"] is None for row in maps.values())
    assert "resolved_files" not in decisions["eaf_replacement"]


def test_active_replacement_rejects_missing_selected_frequency_column(tmp_path):
    maps = _maps(tmp_path)
    path = _table(tmp_path / "replacement.tsv", "1", "WRONG")
    with pytest.raises(ResourcePreflightError, match="SELECTED_AF"):
        _activate(maps, _decision(path, "SELECTED_AF"))


def test_replacement_cannot_activate_without_internal_maf_decision(tmp_path):
    decisions = _decision(tmp_path / "unused.tsv", "AF")
    decisions["eaf_input_interpretation"] = "unresolved"
    with pytest.raises(service.PipelineError, match="MAF"):
        _activate(_maps(tmp_path), decisions)


@pytest.mark.parametrize("reference_mode", ["distinct", "same", "partial", "inconclusive"], ids=["distinct-staged", "same-reused", "partly-reused", "inconclusive-no-workers"])
def test_real_dataset_preparation_selects_then_stages_replacement(tmp_path, monkeypatch, reference_mode):
    same_reference = reference_mode in ("same", "partial")
    module = load_configuration().modules.harmonisation
    count = 2400
    chromosome = ["1"] * 1200 + ["2"] * 1200
    position = [1001 + i % 1200 for i in range(count)]
    frequency = [(.1, .2, .3, .7, .8, .9)[i % 6] for i in range(count)]
    study = tmp_path / "study.tsv"
    pl.DataFrame({
        "CHR": chromosome, "POS": position, "EA": ["C"] * count, "OA": ["A"] * count,
        "AF": [min(p, 1-p) for p in frequency], "BETA": [.2] * count,
        "SE": [.1] * count, "Z": [2.] * count, "P": [.0455002639] * count,
        "INFO": [.9] * count,
    }).write_csv(study, separator="\t")
    default = tmp_path / "default.tsv"
    replacement = tmp_path / "replacement.tsv"
    reference = pl.DataFrame({
        "CHROM": chromosome, "POS": position, "REF": ["A"] * count,
        "ALT": ["C"] * count, "EUR": frequency,
    })
    reference.write_csv(default, separator="\t")
    reference.with_columns((pl.lit(.25) if reference_mode == "inconclusive" else pl.col("EUR")).alias("PANEL_AF")).write_csv(replacement, separator="\t")
    selected_file = default if same_reference else replacement
    selected_column = "EUR" if same_reference else "PANEL_AF"
    resource_maps = {
        chrom: {"default_eaf_file": str(default), "default_eaf_reference_column": "EUR",
                "user_eaf_file": None, "user_eaf_column": None, "user_info_file": None}
        for chrom in ("1", "2")
    }
    if reference_mode == "partial":
        second_default = tmp_path / "second_default.tsv"
        reference.write_csv(second_default, separator="\t")
        resource_maps["2"]["default_eaf_file"] = str(second_default)
    summary = dict.fromkeys(("unique_resource_files_checked", "indexed_vcfs_checked",
        "frequency_vcf_headers_checked", "fasta_indexes_checked", "table_headers_checked",
        "structured_text_files_checked"), 0)
    summary.update(chromosomes=["1", "2"], file_identities=[], resources_by_chromosome=resource_maps)
    # Only build/strand inference and unrelated native resource inventory are
    # replaced. Reader, property detection, AF comparison, resource activation,
    # real partition/staging and logging all remain production implementations.
    monkeypatch.setattr(service, "infer_genome_build", lambda *_a, **_k: {
        "inferred_build": "GRCh37", "input_variants": count,
        "testable_variants": count, "untestable_variants": 0,
    })
    monkeypatch.setattr(service, "resolve_strand_consensus", lambda *_a, **_k: {
        "strand": "forward", "forward": count, "reverse": 0, "ambiguous": 0,
        "dominant_fraction": 1., "reason": "known-forward fixture",
    })
    monkeypatch.setattr(service, "preflight_harmonisation_resources", lambda *_a, **_k: (resource_maps, summary))
    mapping = {
        "gwas_outputname": "study", "sumstat_file": str(study),
        "output_folder": str(tmp_path / "out"), "chr_col": "CHR", "pos_col": "POS",
        "ea_col": "EA", "oa_col": "OA", "eaf_col": "AF", "eaffile": None,
        "eafcolumn": None, "beta_or_col": "BETA", "se_col": "SE", "imp_z_col": "Z",
        "pval_col": "P", "declared_effect_type": "beta", "declared_pvalue_type": "raw",
        "ncontrol": "1000", "imp_info_col": "INFO", "info_source": "internal",
    }
    with PipelineLogger("study", "dataset", str(tmp_path / "logs")) as logger:
        arguments = dict(
            sumstat_file=str(study), sample_column_dict=mapping,
            output_dir_path=tmp_path / "out", resource_folder=str(tmp_path),
            default_eaf_reference_source="fixture", default_comparison_af_file="fixture",
            dbsnp="fixture", resource_layout={}, output_layout=dict(module.output_layout.root),
            gwas2vcf_input=module.gwas2vcf_input.model_dump(), executables={"bcftools": "unused"},
            vcf_config=module.vcf_processing.model_dump(), user_eaf_file=None,
            user_eaf_column=None, user_info_file=None, user_info_column=None,
            default_eaf_reference_column="EUR", default_comparison_af_column="EUR",
            build_reference_files={}, build_check_colmap={}, default_eaf_colmap=MAPPING,
            external_eaf_colmap=MAPPING, external_info_colmap={},
            external_reference_staging=module.external_reference_staging.model_dump(),
            policies=default_policies(), logger=logger, input_line_count=count,
            eaf_replacement_file=str(selected_file), eaf_replacement_column=selected_column,
        )
        if reference_mode == "inconclusive":
            def forbidden(*args, **kwargs):
                pytest.fail("Chromosome workers must not run after an inconclusive replacement gate")
            monkeypatch.setattr(service, "_run_one_round", forbidden)
            arguments["output_dir"] = str(arguments.pop("output_dir_path"))
            arguments.update(max_workers=1, gwas2vcf_main_script_path="unused")
            with pytest.raises(service.StudyAlleleFrequencyError, match="replacement.*inconclusive") as caught:
                service.harmonise_chromosomes(**arguments)
            assert caught.value.study_decisions["eaf_replacement_validation"]["status"] == "inconclusive"
            assert logger.summary()["steps"][-1]["status"] == "failed"
            log = Path(logger.log_path).read_text()
            assert "Whole-study replacement AF validation" in log and "correlation_unavailable" in log
            assert not list((tmp_path / "out").rglob("*.vcf*"))
            return
        result = service._prepare_dataset_for_chromosome_processing(**arguments)
    assert result.study_decisions["eaf_selected_source"] == "replacement_reference"
    assert result.study_decisions["eaf_replacement_validation"]["status"] == "passed"
    assert result.dataset_study_decisions["eaf_replacement_validation"] == result.study_decisions["eaf_replacement_validation"]
    assert result.study_decisions["eaf_replacement"]["resolved_files"] == {
        "1": str(selected_file), "2": str(selected_file),
    }
    assert result.resource_preflight["external_reference_staging"]["source_files_read"] == (0 if reference_mode == "same" else 1)
    for chrom, resources in result.resource_maps.items():
        if reference_mode == "same":
            assert resources["user_eaf_file"] is None
        else:
            staged = pl.read_parquet(resources["user_eaf_file"])
            assert staged.height == 1200
            if reference_mode == "distinct":
                assert staged["PANEL_AF"].cast(pl.Float64).to_list() == frequency[:1200]
            assert staged["CHROM"].unique().to_list() == [chrom]
        assert pl.read_parquet(result.chromosome_files[chrom])["AF"].max() <= .5

    # Continue through the real selected-AF stage, numerical gates, exporter
    # and native adapter process. Known build/strand above are the only
    # synthetic upstream evidence; no population-comparison source is changed.
    from postgwas.modules.harmonisation.allele_frequency import harmonise_allele_frequency
    from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
    from postgwas.modules.harmonisation.effect_validation import final_completeness_check, validate_effect_statistics
    from postgwas.modules.harmonisation.gwas2vcf_export import export_gwas2vcf_input
    from postgwas.modules.harmonisation.gwas2vcf_runner import run_gwas2vcf
    from postgwas.modules.harmonisation.sample_size import harmonise_sample_sizes
    from postgwas.modules.harmonisation.variant_identifiers import harmonise_variant_identifiers
    from postgwas.modules.harmonisation.adapters.gwas2vcf import gwas as adapter_gwas

    frame = pl.read_parquet(result.chromosome_files["1"])
    original_af = frame["AF"].clone()
    columns = dict(result.sample_columns)
    policies = default_policies()
    frame, _effect, columns = harmonise_effect_estimates(
        "1", frame, columns, policies=policies, decision=result.study_decisions,
    )
    frame, qc, columns = harmonise_allele_frequency(
        "1", frame, columns, dict(module.output_layout.root), "\t",
        eaffile=result.resource_maps["1"]["user_eaf_file"], external_eaf_colmap=MAPPING,
        default_eaf_file=str(default), default_eaf_column="EUR", default_eaf_colmap=MAPPING,
        policies=policies, study_decision=result.study_decisions,
    )
    expected_af = frequency[:1200]
    assert frame[columns["eaf_col"]].to_list() == pytest.approx(expected_af)
    assert frame["AF"].to_list() == original_af.to_list()
    assert qc["eaf_provenance"] == "reference_imputed"
    frame, _sizes, columns = harmonise_sample_sizes("1", frame, columns, policies=policies)
    frame, columns = harmonise_variant_identifiers("1", frame, columns, policies=policies)
    frame, _validation, columns = validate_effect_statistics("1", frame, columns, policies=policies)
    frame, _complete, columns = final_completeness_check("1", frame, columns, policies=policies)
    assert frame.height == 1200
    output = tmp_path / "out"
    export_gwas2vcf_input(
        frame, columns, str(output), "study", "1", "GRCh37",
        dict(module.output_layout.root), module.gwas2vcf_input.model_dump(),
    )
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">1\n" + "A" * 3000 + "\n")
    pysam.faidx(str(fasta))
    bcftools = shutil.which("bcftools")
    assert bcftools is not None, "Native indexed-VCF verification needs the installed bcftools."
    run_gwas2vcf(
        "study", "1", "GRCh37", str(output), str(fasta),
        str(Path(adapter_gwas.__file__).with_name("main.py")),
        dict(module.output_layout.root), sys.executable, bcftools, 1200,
    )
    vcf_path = output / module.output_layout.root["chromosome_raw_vcf"].format(
        dataset_id="study", chromosome="1", build="GRCh37",
    )
    assert Path(str(vcf_path) + ".tbi").is_file()
    with pysam.VariantFile(vcf_path) as vcf:
        records = list(vcf.fetch("1"))
        assert len(records) == 1200
        for index, record in enumerate(records):
            sample = record.samples["study"]
            assert "GT" not in sample
            assert sample["AF"] == pytest.approx((expected_af[index],))
            assert sample["ES"] == pytest.approx((.2,))
            assert sample["SE"] == pytest.approx((.1,))
            assert sample["EZ"] == pytest.approx((2.,))
